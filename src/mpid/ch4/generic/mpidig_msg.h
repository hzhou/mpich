/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2020 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */

#ifndef MPIDIG_MSG_H_INCLUDED
#define MPIDIG_MSG_H_INCLUDED

#include "ch4_impl.h"

static inline int handle_unexp_cmpl(MPIR_Request * rreq);
static inline int recv_target_cmpl_cb(MPIR_Request * rreq);

#ifndef MPIDI_CH4_DIRECT_NETMOD
static inline int anysrc_try_cancel(MPIR_Request * rreq)
{
    MPIR_Request *anysource_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq);
    if (!MPIR_STATUS_GET_CANCEL_BIT(anysource_partner->status)) {
        int mpi_errno = MPID_Cancel_recv(anysource_partner);
        MPIR_Assert(mpi_errno == MPI_SUCCESS);
        return MPIR_STATUS_GET_CANCEL_BIT(anysource_partner->status);
    }
    return 1;
}

static inline void anysrc_free_partner(MPIR_Request * rreq)
{
    MPIR_Request *anysource_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq);
    MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq) = NULL;
    MPIDI_REQUEST_ANYSOURCE_PARTNER(anysource_partner) = NULL;
    anysource_partner->status = rreq->status;

    MPIR_Request_free(anysource_partner);
}
#endif

static inline MPIR_Request *MPIDIGI_match_posted(int rank, int tag,
                                                 MPIR_Context_id_t context_id, MPIR_Comm * comm)
{
    MPIR_Request *rreq;
    /* MPIDI_CS_ENTER(); */
#ifdef MPIDI_CH4_DIRECT_NETMOD
    rreq = MPIDIG_dequeue_posted(rank, tag, context_id, &MPIDIG_COMM(comm, posted_list));
#else /* MPIDI_CH4_DIRECT_NETMOD */
    while (TRUE) {
        rreq = MPIDIG_dequeue_posted(rank, tag, context_id, &MPIDIG_COMM(comm, posted_list));

        if (rreq && MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq)) {
            if (!anysrc_try_cancel(rreq)) {
                MPIR_Comm_release(comm);        /* -1 for posted_list */
                MPIR_Datatype_release_if_not_builtin(MPIDIG_REQUEST(rreq, datatype));
                continue;
            }
        }
        break;
    }
#endif /* MPIDI_CH4_DIRECT_NETMOD */
    /* MPIDI_CS_EXIT(); */
    return rreq;
}

static inline int MPIDIGI_enqueue_unexp(int rank, int tag,
                                        MPIR_Context_id_t context_id, MPIR_Comm * comm,
                                        int data_sz, int is_local, MPIR_Request ** p_req)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_Request *rreq;
    rreq = MPIDIG_request_create(MPIR_REQUEST_KIND__RECV, 2);
    MPIR_ERR_CHKANDSTMT(rreq == NULL, mpi_errno, MPIX_ERR_NOREQ, goto fn_fail, "**nomemreq");
    MPIDIG_REQUEST(rreq, datatype) = MPI_BYTE;
    if (data_sz) {
        MPIDIG_REQUEST(rreq, buffer) = (char *) MPL_malloc(data_sz, MPL_MEM_BUFFER);
        MPIDIG_REQUEST(rreq, count) = data_sz;
    } else {
        MPIDIG_REQUEST(rreq, buffer) = NULL;
        MPIDIG_REQUEST(rreq, count) = 0;
    }
    MPIDIG_REQUEST(rreq, rank) = rank;
    MPIDIG_REQUEST(rreq, tag) = tag;
    MPIDIG_REQUEST(rreq, context_id) = context_id;
    MPIDIG_REQUEST(rreq, req->status) |= MPIDIG_REQ_BUSY;
    MPIDIG_REQUEST(rreq, req->status) |= MPIDIG_REQ_UNEXPECTED;
#ifndef MPIDI_CH4_DIRECT_NETMOD
    MPIDI_REQUEST(rreq, is_local) = is_local;
#endif
    MPID_THREAD_CS_ENTER(VCI, MPIDIU_THREAD_MPIDIG_GLOBAL_MUTEX);
    if (comm) {
        MPIR_Comm_add_ref(comm);
        MPIDIG_enqueue_unexp(rreq, &MPIDIG_COMM(comm, unexp_list));
    } else {
        MPIDIG_enqueue_unexp(rreq, MPIDIG_context_id_to_uelist(context_id));
    }
    MPID_THREAD_CS_EXIT(VCI, MPIDIU_THREAD_MPIDIG_GLOBAL_MUTEX);

    *p_req = rreq;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

/* Got matched `rreq`, parse recv datatype, update status */
static inline int do_send_target(void **data, size_t * p_data_sz, int *is_contig,
                                 MPIDIG_am_target_cmpl_cb * target_cmpl_cb, MPIR_Request * rreq)
{
    int dt_contig;
    MPI_Aint dt_true_lb, num_iov;
    MPIR_Datatype *dt_ptr;
    size_t data_sz;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDIG_DO_SEND_TARGET);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDIG_DO_SEND_TARGET);

    *target_cmpl_cb = recv_target_cmpl_cb;
    MPIDIG_REQUEST(rreq, req->seq_no) = MPL_atomic_fetch_add_uint64(&MPIDI_global.nxt_seq_no, 1);

    if (p_data_sz == NULL || 0 == MPIDIG_REQUEST(rreq, count))
        return MPI_SUCCESS;

    MPIDI_Datatype_get_info(MPIDIG_REQUEST(rreq, count), MPIDIG_REQUEST(rreq, datatype), dt_contig,
                            data_sz, dt_ptr, dt_true_lb);
    *is_contig = dt_contig;

    if (dt_contig) {
        *p_data_sz = data_sz;
        *data = (char *) MPIDIG_REQUEST(rreq, buffer) + dt_true_lb;
    } else {
        if (*p_data_sz > data_sz) {
            rreq->status.MPI_ERROR = MPI_ERR_TRUNCATE;
            *p_data_sz = data_sz;
        }

        MPIR_Typerep_iov_len(MPIDIG_REQUEST(rreq, buffer), MPIDIG_REQUEST(rreq, count),
                             MPIDIG_REQUEST(rreq, datatype), 0, data_sz, &num_iov);

        MPIR_Assert(num_iov > 0);
        MPIDIG_REQUEST(rreq, req->iov) =
            (struct iovec *) MPL_malloc(num_iov * sizeof(struct iovec), MPL_MEM_BUFFER);
        MPIR_Assert(MPIDIG_REQUEST(rreq, req->iov));

        int actual_iov_len;
        MPI_Aint actual_iov_bytes;
        MPIR_Typerep_to_iov(MPIDIG_REQUEST(rreq, buffer), MPIDIG_REQUEST(rreq, count),
                            MPIDIG_REQUEST(rreq, datatype), 0, MPIDIG_REQUEST(rreq, req->iov),
                            (int) num_iov, *p_data_sz, &actual_iov_len, &actual_iov_bytes);

        if (actual_iov_bytes != (MPI_Aint) * p_data_sz) {
            rreq->status.MPI_ERROR = MPI_ERR_TYPE;
        }
        *data = MPIDIG_REQUEST(rreq, req->iov);
        *p_data_sz = actual_iov_len;
        MPIDIG_REQUEST(rreq, req->status) |= MPIDIG_REQ_RCV_NON_CONTIG;
    }

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDIG_DO_SEND_TARGET);
    return MPI_SUCCESS;
}

static inline void MPIDIG_recv_copy(void *in_data, int in_data_sz,
                                    void *p_data, int data_sz, int is_contig, MPIR_Request * rreq)
{
    if (!p_data || !data_sz) {
        MPIR_STATUS_SET_COUNT(rreq->status, data_sz);
    } else {
        if (is_contig) {
            if (in_data_sz > data_sz) {
                rreq->status.MPI_ERROR = MPIR_Err_create_code(rreq->status.MPI_ERROR,
                                                              MPIR_ERR_RECOVERABLE, __func__,
                                                              __LINE__, MPI_ERR_TRUNCATE,
                                                              "**truncate",
                                                              "**truncate %d %d %d %d",
                                                              rreq->status.MPI_SOURCE,
                                                              rreq->status.MPI_TAG, data_sz,
                                                              in_data_sz);
            }

            data_sz = MPL_MIN(data_sz, in_data_sz);
            MPIR_Memcpy(p_data, in_data, data_sz);
            MPIR_STATUS_SET_COUNT(rreq->status, data_sz);
        } else {
            int done = 0;
            int rem = in_data_sz;
            struct iovec *iov = (struct iovec *) p_data;
            int iov_len = data_sz;

            for (int i = 0; i < iov_len && rem > 0; i++) {
                int curr_len = MPL_MIN(rem, iov[i].iov_len);
                MPIR_Memcpy(iov[i].iov_base, (char *) in_data + done, curr_len);
                rem -= curr_len;
                done += curr_len;
            }

            if (rem) {
                rreq->status.MPI_ERROR = MPIR_Err_create_code(rreq->status.MPI_ERROR,
                                                              MPIR_ERR_RECOVERABLE, __func__,
                                                              __LINE__, MPI_ERR_TRUNCATE,
                                                              "**truncate",
                                                              "**truncate %d %d %d %d",
                                                              rreq->status.MPI_SOURCE,
                                                              rreq->status.MPI_TAG, data_sz,
                                                              in_data_sz);
            }

            MPIR_STATUS_SET_COUNT(rreq->status, done);
        }
    }
}

static inline int MPIDIG_match_msg(int rank, int tag, MPIR_Context_id_t context_id,
                                   int error_bits, size_t data_sz, int is_local,
                                   int *got_match, MPIR_Request ** req)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Request *rreq = NULL;
    MPIR_Comm *root_comm;

    root_comm = MPIDIG_context_id_to_comm(context_id);
    if (root_comm) {
        rreq = MPIDIGI_match_posted(rank, tag, context_id, root_comm);
    }

    if (rreq == NULL) {
        mpi_errno =
            MPIDIGI_enqueue_unexp(rank, tag, context_id, root_comm, data_sz, is_local, &rreq);
        MPIR_ERR_CHECK(mpi_errno);
        *got_match = 0;
    } else {
        /* rreq != NULL <=> root_comm != NULL */
        MPIR_Assert(root_comm);
        /* Decrement the refcnt when popping a request out from posted_list */
        MPIR_Comm_release(root_comm);
        MPIDIG_REQUEST(rreq, rank) = rank;
        MPIDIG_REQUEST(rreq, tag) = tag;
        MPIDIG_REQUEST(rreq, context_id) = context_id;
        *got_match = 1;
    }

    *req = rreq;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDIG_handle_send(int rank, int tag, MPIR_Context_id_t context_id,
                                     int error_bits, void **data, size_t * p_data_sz,
                                     int is_local, int *is_contig,
                                     MPIDIG_am_target_cmpl_cb * target_cmpl_cb, MPIR_Request ** req)
{
    int mpi_errno = MPI_SUCCESS;

    size_t in_data_sz = p_data_sz ? *p_data_sz : 0;

    int got_match;
    MPIR_Request *rreq;
    mpi_errno = MPIDIG_match_msg(rank, tag, context_id, error_bits, in_data_sz, is_local,
                                 &got_match, &rreq);
    MPIR_ERR_CHECK(mpi_errno);

    rreq->status.MPI_ERROR = error_bits;
    MPIDIG_REQUEST(rreq, req->status) |= MPIDIG_REQ_IN_PROGRESS;

    mpi_errno = do_send_target(data, p_data_sz, is_contig, target_cmpl_cb, rreq);
    MPIR_ERR_CHECK(mpi_errno);

    /* TODO::
     * void *in_data = *data;
     * MPIDIG_recv_copy(in_data, in_data_sz, *data, *p_data_sz, is_contig, rreq);
     * recv_target_cmpl_cb(rreq);
     */

#ifndef MPIDI_CH4_DIRECT_NETMOD
    if (unlikely(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq))) {
        anysrc_free_partner(rreq);
    }
#endif /* MPIDI_CH4_DIRECT_NETMOD */

    *req = rreq;
  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDIG_reply_ssend(MPIR_Request * rreq);
/* This function is called when a receive has completed on the receiver side. The input is the
 * request that has been completed. */
static inline int recv_target_cmpl_cb(MPIR_Request * rreq)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDIG_RECV_TARGET_CMPL_CB);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDIG_RECV_TARGET_CMPL_CB);

    /* Check if this request is supposed to complete next or if it should be delayed. */
    if (!MPIDIG_check_cmpl_order(rreq, recv_target_cmpl_cb))
        return mpi_errno;

    /* If the request contained noncontiguous data, free the iov array that described it. */
    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_RCV_NON_CONTIG) {
        MPL_free(MPIDIG_REQUEST(rreq, req->iov));
    }

    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_UNEXPECTED) {
        mpi_errno = handle_unexp_cmpl(rreq);
        MPIR_ERR_CHECK(mpi_errno);
        goto fn_exit;
    }

    rreq->status.MPI_SOURCE = MPIDIG_REQUEST(rreq, rank);
    rreq->status.MPI_TAG = MPIDIG_REQUEST(rreq, tag);

    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_PEER_SSEND) {
        mpi_errno = MPIDIG_reply_ssend(rreq);
        MPIR_ERR_CHECK(mpi_errno);
    }
#ifndef MPIDI_CH4_DIRECT_NETMOD
    if (unlikely(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq))) {
        int continue_matching = 1;
        if (MPIDI_REQUEST(rreq, is_local)) {
            MPIDI_anysource_matched(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq), MPIDI_SHM,
                                    &continue_matching);
        } else {
            MPIDI_anysource_matched(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq),
                                    MPIDI_NETMOD, &continue_matching);
        }

        MPIR_Request_free(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq));
        MPIDI_REQUEST_ANYSOURCE_PARTNER(MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq)) = NULL;
        MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq) = NULL;
        /* In case of anysource recieve, a request object generated by
         * the shmmod is always returned to the user. Thus, the user-level
         * request cleanup functions (e.g. MPI_Wait) will only free the
         * request object from shmmod. Therefore, netmod should decrement
         * its request object by itself, otherwise it will leak.
         *
         * This logic follows what the OFI netmod is already doing (see
         * MPIDI_OFI_recv_event.) */
        MPIR_Request_free(rreq);
    }
#endif

    MPIR_Datatype_release_if_not_builtin(MPIDIG_REQUEST(rreq, datatype));
    if ((MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_LONG_RTS) &&
        MPIDIG_REQUEST(rreq, req->rreq.match_req) != NULL) {
        /* This block is executed only when the receive is enqueued (trylock/handoff) &&
         * receive was matched with an unexpected long RTS message.
         * `rreq` is the unexpected message received and `sigreq` is the message
         * that came from CH4 (e.g. MPIDI_recv_safe) */
        MPIR_Request *sigreq = MPIDIG_REQUEST(rreq, req->rreq.match_req);
        sigreq->status = rreq->status;
        MPIR_Request_add_ref(sigreq);
        MPID_Request_complete(sigreq);
        /* Free the unexpected request on behalf of the user */
        MPIR_Request_free(rreq);
    }
    MPID_Request_complete(rreq);
  fn_exit:
    MPIDIG_progress_compl_list();
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDIG_RECV_TARGET_CMPL_CB);
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDIG_handle_unexp_mrecv(MPIR_Request * rreq);

static inline int handle_unexp_cmpl(MPIR_Request * rreq)
{
    int mpi_errno = MPI_SUCCESS, in_use;
    MPIR_Comm *root_comm;
    MPIR_Request *match_req = NULL;
    size_t nbytes;
    int dt_contig;
    MPI_Aint dt_true_lb;
    MPIR_Datatype *dt_ptr;
    size_t dt_sz;

#ifndef MPIDI_CH4_DIRECT_NETMOD
    MPIR_Request *anysource_partner = NULL;
#endif

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDIG_HANDLE_UNEXP_CMPL);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDIG_HANDLE_UNEXP_CMPL);

    /* Check if this message has already been claimed by mprobe. */
    /* MPIDI_CS_ENTER(); */
    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_UNEXP_DQUED) {
        /* This request has been claimed by mprobe */
        if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_UNEXP_CLAIMED) {
            /* mrecv has been already called */
            MPIDIG_handle_unexp_mrecv(rreq);
        } else {
            /* mrecv has not been called yet -- just take out the busy flag so that
             * mrecv in future knows this request is ready */
            MPIDIG_REQUEST(rreq, req->status) &= ~MPIDIG_REQ_BUSY;
        }
        /* MPIDI_CS_EXIT(); */
        goto fn_exit;
    }
    /* MPIDI_CS_EXIT(); */

    root_comm = MPIDIG_context_id_to_comm(MPIDIG_REQUEST(rreq, context_id));

    /* If this request was previously matched, but not handled */
    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_MATCHED) {
        match_req = (MPIR_Request *) MPIDIG_REQUEST(rreq, req->rreq.match_req);

#ifndef MPIDI_CH4_DIRECT_NETMOD
        if (unlikely(match_req && MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req))) {
            anysource_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req);
            if (!MPIR_STATUS_GET_CANCEL_BIT(anysource_partner->status)) {
                mpi_errno = MPID_Cancel_recv(anysource_partner);
                if (mpi_errno != MPI_SUCCESS) {
                    goto fn_fail;
                }
                /* What should we do if the anysource partner request is not canceled? */
                MPIR_Assertp(MPIR_STATUS_GET_CANCEL_BIT(anysource_partner->status));
            }
            MPIR_Request_free(MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req));
            MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req) = NULL;
            MPIDI_REQUEST_ANYSOURCE_PARTNER(anysource_partner) = NULL;
        }
#endif /* MPIDI_CH4_DIRECT_NETMOD */

    } else {
        /* If this message hasn't been matched yet, look for it in the posted queue. */
        /* MPIDI_CS_ENTER(); */
        if (root_comm) {
#ifdef MPIDI_CH4_DIRECT_NETMOD
            match_req =
                MPIDIG_dequeue_posted(MPIDIG_REQUEST(rreq, rank),
                                      MPIDIG_REQUEST(rreq, tag),
                                      MPIDIG_REQUEST(rreq, context_id),
                                      &MPIDIG_COMM(root_comm, posted_list));
#else /* MPIDI_CH4_DIRECT_NETMOD */
            int continue_matching = 1;
            while (continue_matching) {
                match_req =
                    MPIDIG_dequeue_posted(MPIDIG_REQUEST(rreq, rank),
                                          MPIDIG_REQUEST(rreq, tag),
                                          MPIDIG_REQUEST(rreq, context_id),
                                          &MPIDIG_COMM(root_comm, posted_list));

                if (match_req && MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req)) {
                    anysource_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req);

                    mpi_errno = MPIDI_anysource_matched(anysource_partner,
                                                        MPIDI_REQUEST(rreq, is_local) ?
                                                        MPIDI_SHM : MPIDI_NETMOD,
                                                        &continue_matching);

                    MPIR_ERR_CHECK(mpi_errno);

                    MPIR_Request_free(MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req));
                    MPIDI_REQUEST_ANYSOURCE_PARTNER(match_req) = NULL;
                    MPIDI_REQUEST_ANYSOURCE_PARTNER(anysource_partner) = NULL;
                }

                break;
            }
#endif /* MPIDI_CH4_DIRECT_NETMOD */
        }

        /* If we found a matching request, remove it from the unexpected queue and clean things up
         * before we move the data around. */
        if (match_req) {
            MPIDIG_delete_unexp(rreq, &MPIDIG_COMM(root_comm, unexp_list));
            /* Decrement the counter twice, one for posted_list and the other for unexp_list */
            MPIR_Comm_release(root_comm);
            MPIR_Comm_release(root_comm);
        }
        /* MPIDI_CS_EXIT(); */
    }

    /* If we didn't match the request, unmark the busy bit and skip the data movement below. */
    if (!match_req) {
        MPIDIG_REQUEST(rreq, req->status) &= ~MPIDIG_REQ_BUSY;
        goto fn_exit;
    }

    match_req->status.MPI_SOURCE = MPIDIG_REQUEST(rreq, rank);
    match_req->status.MPI_TAG = MPIDIG_REQUEST(rreq, tag);

    /* Figure out how much data needs to be moved. */
    MPIDI_Datatype_get_info(MPIDIG_REQUEST(match_req, count),
                            MPIDIG_REQUEST(match_req, datatype),
                            dt_contig, dt_sz, dt_ptr, dt_true_lb);
    MPIR_Datatype_get_size_macro(MPIDIG_REQUEST(match_req, datatype), dt_sz);

    /* Make sure this request has the right amount of data in it. */
    if (MPIDIG_REQUEST(rreq, count) > dt_sz * MPIDIG_REQUEST(match_req, count)) {
        rreq->status.MPI_ERROR = MPI_ERR_TRUNCATE;
        nbytes = dt_sz * MPIDIG_REQUEST(match_req, count);
    } else {
        rreq->status.MPI_ERROR = MPI_SUCCESS;
        nbytes = MPIDIG_REQUEST(rreq, count);   /* incoming message is always count of bytes. */
    }

    MPIR_STATUS_SET_COUNT(match_req->status, nbytes);
    MPIDIG_REQUEST(rreq, count) = dt_sz > 0 ? nbytes / dt_sz : 0;

    /* Perform the data copy (using the datatype engine if necessary for non-contig transfers) */
    if (!dt_contig) {
        MPI_Aint actual_unpack_bytes;
        mpi_errno = MPIR_Typerep_unpack(MPIDIG_REQUEST(rreq, buffer), nbytes,
                                        MPIDIG_REQUEST(match_req, buffer),
                                        MPIDIG_REQUEST(match_req, count),
                                        MPIDIG_REQUEST(match_req, datatype), 0,
                                        &actual_unpack_bytes);
        MPIR_ERR_CHECK(mpi_errno);

        if (actual_unpack_bytes != (MPI_Aint) nbytes) {
            mpi_errno = MPIR_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                             __FUNCTION__, __LINE__,
                                             MPI_ERR_TYPE, "**dtypemismatch", 0);
            match_req->status.MPI_ERROR = mpi_errno;
        }
    } else {
        MPIR_Memcpy((char *) MPIDIG_REQUEST(match_req, buffer) + dt_true_lb,
                    MPIDIG_REQUEST(rreq, buffer), nbytes);
    }

    /* Now that the unexpected message has been completed, unset the status bit. */
    MPIDIG_REQUEST(rreq, req->status) &= ~MPIDIG_REQ_UNEXPECTED;

    /* If this is a synchronous send, send the reply back to the sender to unlock them. */
    if (MPIDIG_REQUEST(rreq, req->status) & MPIDIG_REQ_PEER_SSEND) {
        mpi_errno = MPIDIG_reply_ssend(rreq);
        MPIR_ERR_CHECK(mpi_errno);
    }
#ifndef MPIDI_CH4_DIRECT_NETMOD
    if (unlikely(anysource_partner)) {
        anysource_partner->status = match_req->status;
    }
#endif /* MPIDI_CH4_DIRECT_NETMOD */

    MPIR_Datatype_release_if_not_builtin(MPIDIG_REQUEST(match_req, datatype));
    MPL_free(MPIDIG_REQUEST(rreq, buffer));
    MPIR_Object_release_ref(rreq, &in_use);
    MPID_Request_complete(rreq);
    MPID_Request_complete(match_req);
  fn_exit:
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDIG_HANDLE_UNEXP_CMPL);
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

#endif /* MPIDIG_MSG_H_INCLUDED */
