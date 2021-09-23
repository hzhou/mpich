/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include <mpidimpl.h>
#include "ofi_impl.h"
#include "ofi_events.h"

#define MPIDI_OFI_MR_KEY_PREFIX_SHIFT 63

int MPIDI_OFI_retry_progress(void)
{
    /* We do not call progress on hooks form netmod level
     * because it is not reentrant safe.
     */
    return MPID_Progress_test(NULL);
}

typedef struct MPIDI_OFI_mr_key_allocator_t {
    uint64_t chunk_size;
    uint64_t num_ints;
    uint64_t last_free_mr_key;
    uint64_t *bitmask;
} MPIDI_OFI_mr_key_allocator_t;

static MPIDI_OFI_mr_key_allocator_t mr_key_allocator;

int MPIDI_OFI_mr_key_allocator_init(void)
{
    int mpi_errno = MPI_SUCCESS;

    mr_key_allocator.chunk_size = 128;
    mr_key_allocator.num_ints = mr_key_allocator.chunk_size;
    mr_key_allocator.last_free_mr_key = 0;
    mr_key_allocator.bitmask = MPL_malloc(sizeof(uint64_t) * mr_key_allocator.num_ints,
                                          MPL_MEM_RMA);
    MPIR_ERR_CHKANDSTMT(mr_key_allocator.bitmask == NULL, mpi_errno,
                        MPI_ERR_NO_MEM, goto fn_fail, "**nomem");
    memset(mr_key_allocator.bitmask, 0xFF, sizeof(uint64_t) * mr_key_allocator.num_ints);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

#define MPIDI_OFI_INDEX_CALC(val,nval,shift,mask) \
    if ((val & mask) == 0) {                              \
        val >>= shift##ULL;                               \
        nval += shift;                                    \
    }

/* when key_type is MPIDI_OFI_LOCAL_MR_KEY, the input requested_key is ignored
 * and can be passed as MPIDI_OFI_INVALID_MR_KEY because mr key allocator will
 * decide which key to use. When key_type is MPIDI_OFI_COLL_MR_KEY, user should
 * pass a collectively unique key as requested_key and mr key allocator will mark
 * coll bit of the key and return to user.
 * we use highest bit of key to distinguish coll (user-specific key) and local
 * (auto-generated) 64-bits key; since the highest bit is reserved for key type,
 * a valid key has maximal 63 bits. */
uint64_t MPIDI_OFI_mr_key_alloc(int key_type, uint64_t requested_key)
{
    uint64_t ret_key = MPIDI_OFI_INVALID_MR_KEY;

    switch (key_type) {
        case MPIDI_OFI_LOCAL_MR_KEY:
            {
                uint64_t i;
                for (i = mr_key_allocator.last_free_mr_key; i < mr_key_allocator.num_ints; i++) {
                    if (mr_key_allocator.bitmask[i]) {
                        register uint64_t val, nval;
                        val = mr_key_allocator.bitmask[i];
                        nval = 2;
                        MPIDI_OFI_INDEX_CALC(val, nval, 32, 0xFFFFFFFFULL);
                        MPIDI_OFI_INDEX_CALC(val, nval, 16, 0xFFFFULL);
                        MPIDI_OFI_INDEX_CALC(val, nval, 8, 0xFFULL);
                        MPIDI_OFI_INDEX_CALC(val, nval, 4, 0xFULL);
                        MPIDI_OFI_INDEX_CALC(val, nval, 2, 0x3ULL);
                        nval -= val & 0x1ULL;
                        mr_key_allocator.bitmask[i] &= ~(0x1ULL << (nval - 1));
                        mr_key_allocator.last_free_mr_key = i;
                        ret_key = i * sizeof(uint64_t) * 8 + (nval - 1);
                        /* assert local key does not exceed its range */
                        MPIR_Assert((ret_key & (1ULL << MPIDI_OFI_MR_KEY_PREFIX_SHIFT)) == 0);
                        break;
                    }
                    if (i == mr_key_allocator.num_ints - 1) {
                        mr_key_allocator.num_ints += mr_key_allocator.chunk_size;
                        mr_key_allocator.bitmask = MPL_realloc(mr_key_allocator.bitmask,
                                                               sizeof(uint64_t) *
                                                               mr_key_allocator.num_ints,
                                                               MPL_MEM_RMA);
                        MPIR_Assert(mr_key_allocator.bitmask);
                        memset(&mr_key_allocator.bitmask[i + 1], 0xFF,
                               sizeof(uint64_t) * mr_key_allocator.chunk_size);
                    }
                }
                break;
            }

        case MPIDI_OFI_COLL_MR_KEY:
            {
                MPIR_Assert(requested_key != MPIDI_OFI_INVALID_MR_KEY);
                ret_key = requested_key | (1ULL << MPIDI_OFI_MR_KEY_PREFIX_SHIFT);
                break;
            }

        default:
            {
                MPIR_Assert(0);
            }
    }

    return ret_key;
}

void MPIDI_OFI_mr_key_free(int key_type, uint64_t alloc_key)
{

    switch (key_type) {
        case MPIDI_OFI_LOCAL_MR_KEY:
            {
                uint64_t int_index, bitpos, numbits;

                numbits = sizeof(uint64_t) * 8;
                int_index = alloc_key / numbits;
                bitpos = alloc_key % numbits;
                mr_key_allocator.last_free_mr_key =
                    MPL_MIN(int_index, mr_key_allocator.last_free_mr_key);
                mr_key_allocator.bitmask[int_index] |= (0x1ULL << bitpos);
                break;
            }

        case MPIDI_OFI_COLL_MR_KEY:
            {
                MPIR_Assert(alloc_key != MPIDI_OFI_INVALID_MR_KEY);
                break;
            }

        default:
            {
                MPIR_Assert(0);
            }
    }
}

void MPIDI_OFI_mr_key_allocator_destroy(void)
{
    MPL_free(mr_key_allocator.bitmask);
}

int MPIDI_OFI_huge_ack(int rank, MPIR_Comm * comm, MPI_Request handle, int vni_src, int vni_dst)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_OFI_send_control_t ctrl;
    ctrl.type = MPIDI_OFI_CTRL_HUGE_ACK;
    ctrl.u.huge_ack.ackreq = handle;

    /* note: this is receiver ack sender */
    mpi_errno = MPIDI_NM_am_send_hdr(rank, comm, MPIDI_OFI_INTERNAL_HANDLER_CONTROL,
                                     &ctrl, sizeof(ctrl), vni_dst, vni_src);
    return mpi_errno;
}

int MPIDI_OFI_huge_probe_msgsize(int rank, MPIR_Comm * comm, MPI_Request handle, int vni_src,
                                 int vni_dst, MPI_Aint * msgsize)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_OFI_send_control_t ctrl;
    MPIDI_OFI_huge_probe_reply_t reply;

    reply.done = false;
    ctrl.type = MPIDI_OFI_CTRL_HUGE_PROBE;
    ctrl.u.huge_probe.ackreq = handle;
    ctrl.u.huge_probe.reply_ptr = &reply;
    ctrl.u.huge_probe.dst_rank = comm->rank;
    ctrl.u.huge_probe.vni_src = vni_src;
    ctrl.u.huge_probe.vni_dst = vni_dst;

    /* note: this is receiver querying sender */
    mpi_errno = MPIDI_NM_am_send_hdr(rank, comm,
                                     MPIDI_OFI_INTERNAL_HANDLER_CONTROL,
                                     &ctrl, sizeof(ctrl), vni_dst, vni_src);
    MPIR_ERR_CHECK(mpi_errno);

    MPIDI_OFI_PROGRESS_WHILE(!reply.done, vni_dst);

    *msgsize = reply.msgsize;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_OFI_control_handler(void *am_hdr, void *data, MPI_Aint data_sz,
                              uint32_t attr, MPIR_Request ** req)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_OFI_send_control_t *ctrlsend = (MPIDI_OFI_send_control_t *) am_hdr;

    if (attr & MPIDIG_AM_ATTR__IS_ASYNC) {
        *req = NULL;
    }

    switch (ctrlsend->type) {
        case MPIDI_OFI_CTRL_HUGE_ACK:{
                MPIR_Request *sreq;
                MPIR_Request_get_ptr(ctrlsend->u.huge_ack.ackreq, sreq);
                MPIR_Assert(sreq);
                mpi_errno = MPIDI_OFI_dispatch_function(NULL, sreq);
            }
            break;

        case MPIDI_OFI_CTRL_HUGE_PROBE:{
                MPIR_Request *sreq;
                MPIR_Request_get_ptr(ctrlsend->u.huge_probe.ackreq, sreq);
                MPIR_Assert(sreq);

                MPIDI_OFI_huge_info_t *info = MPIDI_OFI_REQUEST(sreq, util.inject_buf);

                MPIDI_OFI_send_control_t reply;
                reply.type = MPIDI_OFI_CTRL_HUGE_PROBE_REPLY;
                reply.u.huge_probe_reply.reply_ptr = ctrlsend->u.huge_probe.reply_ptr;
                reply.u.huge_probe_reply.msgsize = info->msgsize;

                int dst_rank = ctrlsend->u.huge_probe.dst_rank;
                int vni_src = ctrlsend->u.huge_probe.vni_src;
                int vni_dst = ctrlsend->u.huge_probe.vni_dst;
                /* note: this is sender replying receiver's probe */
                mpi_errno = MPIDI_NM_am_send_hdr_reply(sreq->comm, dst_rank,
                                                       MPIDI_OFI_INTERNAL_HANDLER_CONTROL,
                                                       &reply, sizeof(reply), vni_src, vni_dst);
            }
            break;

        case MPIDI_OFI_CTRL_HUGE_PROBE_REPLY:{
                MPIDI_OFI_huge_probe_reply_t *reply_ptr = ctrlsend->u.huge_probe_reply.reply_ptr;
                reply_ptr->msgsize = ctrlsend->u.huge_probe_reply.msgsize;
                reply_ptr->done = true;
            }
            break;

        default:
            fprintf(stderr, "Bad control type: 0x%08x  %d\n", ctrlsend->type, ctrlsend->type);
            MPIR_Assert(0);
    }

  fn_exit:
    return mpi_errno;
}


/* MPI Datatype Processing for RMA */
#define isS_INT(x) ((x)==MPI_INTEGER ||                                \
    (x) == MPI_INT32_T || (x) == MPI_INTEGER4 ||       \
                     (x) == MPI_INT)
#define isUS_INT(x) ((x) == MPI_UINT32_T || (x) == MPI_UNSIGNED)
#define isS_SHORT(x) ((x) == MPI_SHORT || (x) == MPI_INT16_T ||        \
                       (x) == MPI_INTEGER2)
#define isUS_SHORT(x) ((x) == MPI_UNSIGNED_SHORT || (x) == MPI_UINT16_T)
#define isS_CHAR(x) ((x) == MPI_SIGNED_CHAR || (x) == MPI_INT8_T ||    \
                      (x) == MPI_INTEGER1 || (x) == MPI_CHAR)
#define isUS_CHAR(x) ((x) == MPI_BYTE ||                               \
                       (x) == MPI_UNSIGNED_CHAR || (x) == MPI_UINT8_T)
#define isS_LONG(x) ((x) == MPI_LONG || (x) == MPI_AINT)
#define isUS_LONG(x) ((x) == MPI_UNSIGNED_LONG)
#define isS_LONG_LONG(x) ((x) == MPI_INT64_T || (x) == MPI_OFFSET ||   \
    (x) == MPI_INTEGER8 || (x) == MPI_LONG_LONG || \
                           (x) == MPI_LONG_LONG_INT || (x) == MPI_COUNT)
#define isUS_LONG_LONG(x) ((x) == MPI_UINT64_T || (x) == MPI_UNSIGNED_LONG_LONG)
#define isFLOAT(x) ((x) == MPI_FLOAT || (x) == MPI_REAL)
#define isDOUBLE(x) ((x) == MPI_DOUBLE || (x) == MPI_DOUBLE_PRECISION)
#define isLONG_DOUBLE(x) ((x) == MPI_LONG_DOUBLE)
#define isLOC_TYPE(x) ((x) == MPI_2REAL || (x) == MPI_2DOUBLE_PRECISION || \
    (x) == MPI_2INTEGER || (x) == MPI_FLOAT_INT ||  \
    (x) == MPI_DOUBLE_INT || (x) == MPI_LONG_INT || \
    (x) == MPI_2INT || (x) == MPI_SHORT_INT ||      \
                        (x) == MPI_LONG_DOUBLE_INT)
#define isBOOL(x) ((x) == MPI_C_BOOL)
#define isLOGICAL(x) ((x) == MPI_LOGICAL)
#define isSINGLE_COMPLEX(x) ((x) == MPI_COMPLEX || (x) == MPI_C_FLOAT_COMPLEX)
#define isDOUBLE_COMPLEX(x) ((x) == MPI_DOUBLE_COMPLEX || (x) == MPI_COMPLEX8 || \
                              (x) == MPI_C_DOUBLE_COMPLEX)

static bool check_mpi_acc_valid(MPI_Datatype dtype, MPI_Op op)
{
    bool valid_flag = false;

    /* Check if the <datatype, op> is supported by MPICH. Note that MPICH
     * supports more combinations than that specified in standard (see definition
     * of these checking routines for extended support). */

    /* Predefined reduce operation, NO_OP, REPLACE */
    if (op != MPI_OP_NULL) {
        int mpi_errno;
        mpi_errno = (*MPIR_OP_HDL_TO_DTYPE_FN(op)) (dtype);
        if (mpi_errno == MPI_SUCCESS)
            valid_flag = TRUE;
    } else {
        /* Compare and swap */
        if (MPIR_Type_is_rma_atomic(dtype))
            valid_flag = true;
    }

    return valid_flag;
}

static int mpi_to_ofi(MPI_Datatype dt, enum fi_datatype *fi_dt, MPI_Op op, enum fi_op *fi_op)
{
    *fi_dt = FI_DATATYPE_LAST;
    *fi_op = FI_ATOMIC_OP_LAST;

    if (isS_INT(dt))
        *fi_dt = FI_INT32;
    else if (isUS_INT(dt))
        *fi_dt = FI_UINT32;
    else if (isFLOAT(dt))
        *fi_dt = FI_FLOAT;
    else if (isDOUBLE(dt))
        *fi_dt = FI_DOUBLE;
    else if (isLONG_DOUBLE(dt))
        *fi_dt = FI_LONG_DOUBLE;
    else if (isS_CHAR(dt))
        *fi_dt = FI_INT8;
    else if (isUS_CHAR(dt))
        *fi_dt = FI_UINT8;
    else if (isS_SHORT(dt))
        *fi_dt = FI_INT16;
    else if (isUS_SHORT(dt))
        *fi_dt = FI_UINT16;
    else if (isS_LONG(dt))
        *fi_dt = FI_INT64;
    else if (isUS_LONG(dt))
        *fi_dt = FI_UINT64;
    else if (isS_LONG_LONG(dt))
        *fi_dt = FI_INT64;
    else if (isUS_LONG_LONG(dt))
        *fi_dt = FI_UINT64;
    else if (isSINGLE_COMPLEX(dt))
        *fi_dt = FI_FLOAT_COMPLEX;
    else if (isDOUBLE_COMPLEX(dt))
        *fi_dt = FI_DOUBLE_COMPLEX;
    else if (isLOC_TYPE(dt))
        *fi_dt = FI_DATATYPE_LAST;
    else if (isLOGICAL(dt))
        *fi_dt = FI_UINT32;
    else if (isBOOL(dt))
        *fi_dt = FI_UINT8;

    if (*fi_dt == FI_DATATYPE_LAST)
        goto fn_fail;

    *fi_op = FI_ATOMIC_OP_LAST;

    switch (op) {
        case MPI_SUM:
            *fi_op = FI_SUM;
            goto fn_exit;

        case MPI_PROD:
            *fi_op = FI_PROD;
            goto fn_exit;

        case MPI_MAX:
            *fi_op = FI_MAX;
            goto fn_exit;

        case MPI_MIN:
            *fi_op = FI_MIN;
            goto fn_exit;

        case MPI_BAND:
            *fi_op = FI_BAND;
            goto fn_exit;

        case MPI_BOR:
            *fi_op = FI_BOR;
            goto fn_exit;

        case MPI_BXOR:
            *fi_op = FI_BXOR;
            goto fn_exit;
            break;

        case MPI_LAND:
            if (isLONG_DOUBLE(dt))
                goto fn_fail;

            *fi_op = FI_LAND;
            goto fn_exit;

        case MPI_LOR:
            if (isLONG_DOUBLE(dt))
                goto fn_fail;

            *fi_op = FI_LOR;
            goto fn_exit;

        case MPI_LXOR:
            if (isLONG_DOUBLE(dt))
                goto fn_fail;

            *fi_op = FI_LXOR;
            goto fn_exit;

        case MPI_REPLACE:{
                *fi_op = FI_ATOMIC_WRITE;
                goto fn_exit;
            }

        case MPI_NO_OP:{
                *fi_op = FI_ATOMIC_READ;
                goto fn_exit;
            }

        case MPI_OP_NULL:{
                *fi_op = FI_CSWAP;
                goto fn_exit;
            }

        default:
            goto fn_fail;
    }

  fn_exit:
    return MPI_SUCCESS;
  fn_fail:
    return -1;
}

#define _TBL MPIDI_OFI_global.win_op_table[i][j]
#define CHECK_ATOMIC(fcn,field1,field2)            \
  atomic_count = 0;                                \
  ret = fcn(ep, fi_dt, fi_op, &atomic_count);      \
  if (ret == 0 && atomic_count != 0) {             \
    _TBL.field1 = 1;                               \
    _TBL.field2 = atomic_count;                    \
  }

static void create_dt_map(struct fid_ep *ep)
{
    int i, j;
    size_t dtsize[FI_DATATYPE_LAST];
    dtsize[FI_INT8] = sizeof(int8_t);
    dtsize[FI_UINT8] = sizeof(uint8_t);
    dtsize[FI_INT16] = sizeof(int16_t);
    dtsize[FI_UINT16] = sizeof(uint16_t);
    dtsize[FI_INT32] = sizeof(int32_t);
    dtsize[FI_UINT32] = sizeof(uint32_t);
    dtsize[FI_INT64] = sizeof(int64_t);
    dtsize[FI_UINT64] = sizeof(uint64_t);
    dtsize[FI_FLOAT] = sizeof(float);
    dtsize[FI_DOUBLE] = sizeof(double);
    dtsize[FI_FLOAT_COMPLEX] = sizeof(float complex);
    dtsize[FI_DOUBLE_COMPLEX] = sizeof(double complex);
    dtsize[FI_LONG_DOUBLE] = sizeof(long double);
    dtsize[FI_LONG_DOUBLE_COMPLEX] = sizeof(long double complex);

    /* when atomics are disabled and atomics capability are not
     * enabled call of fi_atomic*** may crash */
    MPIR_Assert(MPIDI_OFI_ENABLE_ATOMICS);

    memset(MPIDI_OFI_global.win_op_table, 0, sizeof(MPIDI_OFI_global.win_op_table));

    for (i = 0; i < MPIR_DATATYPE_N_PREDEFINED; i++) {
        MPI_Datatype dt = MPIR_Datatype_predefined_get_type(i);

        /* MPICH sets predefined datatype handles to MPI_DATATYPE_NULL if they are not
         * supported on the target platform. Skip it. */
        if (dt == MPI_DATATYPE_NULL)
            continue;

        for (j = 0; j < MPIDIG_ACCU_NUM_OP; j++) {
            MPI_Op op = MPIDIU_win_acc_get_op(j);
            enum fi_datatype fi_dt = (enum fi_datatype) -1;
            enum fi_op fi_op = (enum fi_op) -1;

            mpi_to_ofi(dt, &fi_dt, op, &fi_op);
            MPIR_Assert(fi_dt != (enum fi_datatype) -1);
            MPIR_Assert(fi_op != (enum fi_op) -1);
            _TBL.dt = fi_dt;
            _TBL.op = fi_op;
            _TBL.atomic_valid = 0;
            _TBL.max_atomic_count = 0;
            _TBL.max_fetch_atomic_count = 0;
            _TBL.max_compare_atomic_count = 0;
            _TBL.mpi_acc_valid = check_mpi_acc_valid(dt, op);
            ssize_t ret;
            size_t atomic_count;

            if (fi_dt != FI_DATATYPE_LAST && fi_op != FI_ATOMIC_OP_LAST) {
                CHECK_ATOMIC(fi_atomicvalid, atomic_valid, max_atomic_count);
                CHECK_ATOMIC(fi_fetch_atomicvalid, fetch_atomic_valid, max_fetch_atomic_count);
                CHECK_ATOMIC(fi_compare_atomicvalid, compare_atomic_valid,
                             max_compare_atomic_count);
                _TBL.dtsize = dtsize[fi_dt];
            }
        }
    }
}

void MPIDI_OFI_index_datatypes(struct fid_ep *ep)
{
    /* do not generate map when atomics are not enabled */
    if (MPIDI_OFI_ENABLE_ATOMICS) {
        create_dt_map(ep);
    }
}
