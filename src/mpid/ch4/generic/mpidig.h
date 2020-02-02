/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef MPIDIG_H_INCLUDED
#define MPIDIG_H_INCLUDED

#include <mpidimpl.h>

#define MPIDI_AM_HANDLERS_MAX (64)

typedef int (*MPIDIG_am_target_cmpl_cb) (MPIR_Request * req);
typedef int (*MPIDIG_am_origin_cb) (MPIR_Request * req);

/* Callback function setup by handler register function */
/* for short cases, output arguments are NULL */
typedef int (*MPIDIG_am_target_msg_cb)
 (int handler_id, void *am_hdr, void **data,    /* data should be iovs if *is_contig is false */
  size_t * data_sz, int is_local,       /* SHM or NM directly specifies locality */
  int *is_contig, MPIDIG_am_target_cmpl_cb * target_cmpl_cb,    /* completion handler */
  MPIR_Request ** req);         /* if allocated, need pointer to completion function */

typedef struct MPIDIG_global_t {
    MPIDIG_am_target_msg_cb target_msg_cbs[MPIDI_AM_HANDLERS_MAX];
    MPIDIG_am_origin_cb origin_cbs[MPIDI_AM_HANDLERS_MAX];
} MPIDIG_global_t;
extern MPIDIG_global_t MPIDIG_global;

int MPIDIG_am_reg_cb(int handler_id,
                     MPIDIG_am_origin_cb origin_cb, MPIDIG_am_target_msg_cb target_msg_cb);
int MPIDIG_init(void);
void MPIDIG_finalize(void);

int MPIDIG_RMA_Init_targetcb_pvars(void);

int MPIDIG_am_put_init(void);
int MPIDIG_do_put(const void *origin_addr, int origin_count,
                  MPI_Datatype origin_datatype, int target_rank,
                  MPI_Aint target_disp, int target_count,
                  MPI_Datatype target_datatype, MPIR_Win * win, MPIR_Request ** sreq_ptr);

int MPIDIG_am_get_init(void);
int MPIDIG_do_get(void *origin_addr, int origin_count, MPI_Datatype origin_datatype,
                  int target_rank, MPI_Aint target_disp, int target_count,
                  MPI_Datatype target_datatype, MPIR_Win * win, MPIR_Request ** sreq_ptr);

int MPIDIG_am_cswap_init(void);
int MPIDIG_do_cswap(const void *origin_addr, const void *compare_addr,
                    void *result_addr, MPI_Datatype datatype,
                    int target_rank, MPI_Aint target_disp, MPIR_Win * win);

int MPIDIG_am_acc_init(void);
int MPIDIG_do_accumulate(const void *origin_addr, int origin_count,
                         MPI_Datatype origin_datatype, int target_rank,
                         MPI_Aint target_disp, int target_count,
                         MPI_Datatype target_datatype,
                         MPI_Op op, MPIR_Win * win, MPIR_Request ** sreq_ptr);

int MPIDIG_am_get_acc_init(void);
int MPIDIG_do_get_accumulate(const void *origin_addr, int origin_count,
                             MPI_Datatype origin_datatype, void *result_addr, int result_count,
                             MPI_Datatype result_datatype, int target_rank, MPI_Aint target_disp,
                             int target_count, MPI_Datatype target_datatype,
                             MPI_Op op, MPIR_Win * win, MPIR_Request ** sreq_ptr);

int MPIDIG_am_win_sync_init(void);

int MPIDIG_am_comm_abort_init(void);
int MPIDIG_comm_abort(MPIR_Comm * comm, int exit_code);

#ifdef MPIDI_CH4_DIRECT_NETMOD
#define MPIDIG_AM_SEND(is_local, rank, comm, ID, data, count, datatype, req) \
    mpi_errno = MPIDI_NM_am_isend(rank, comm, ID, &am_hdr, sizeof(am_hdr), \
                                  data, count, datatype, req)

#define MPIDIG_AM_SENDV(is_local, rank, comm, ID, iov, iov_len, data, count, datatype, req) \
    mpi_errno = MPIDI_NM_am_isendv(rank, comm, ID, iov, iov_len, \
                                   data, count, datatype, req)

#define MPIDIG_AM_SEND_HDR(is_local, rank, comm, ID) \
    mpi_errno = MPIDI_NM_am_send_hdr(rank, comm, ID, &am_hdr, sizeof(am_hdr))

#else
#define MPIDIG_AM_SEND(is_local, rank, comm, ID, data, count, datatype, req) \
    do { \
        if (is_local) { \
            mpi_errno = MPIDI_SHM_am_isend(rank, comm, ID, &am_hdr, sizeof(am_hdr), \
                                           data, count, datatype, req); \
        } else { \
            mpi_errno = MPIDI_NM_am_isend(rank, comm, ID, &am_hdr, sizeof(am_hdr), \
                                          data, count, datatype, req); \
        } \
    } while (0)

#define MPIDIG_AM_SENDV(is_local, rank, comm, ID, iov, iov_len, data, count, datatype, req) \
    do { \
        if (is_local) { \
            mpi_errno = MPIDI_SHM_am_isendv(rank, comm, ID, iov, iov_len, \
                                            data, count, datatype, req); \
        } else { \
            mpi_errno = MPIDI_NM_am_isendv(rank, comm, ID, iov, iov_len, \
                                           data, count, datatype, req); \
        } \
    } while (0)

#define MPIDIG_AM_SEND_HDR(is_local, rank, comm, ID) \
    do { \
        if (is_local) { \
            mpi_errno = MPIDI_SHM_am_send_hdr(rank, comm, ID, &am_hdr, sizeof(am_hdr)); \
        } else { \
            mpi_errno = MPIDI_NM_am_send_hdr(rank, comm, ID, &am_hdr, sizeof(am_hdr)); \
        } \
    } while (0)
#endif

#endif /* MPIDIG_H_INCLUDED */
