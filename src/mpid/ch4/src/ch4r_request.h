/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef CH4R_REQUEST_H_INCLUDED
#define CH4R_REQUEST_H_INCLUDED

#include "ch4_types.h"
#include "ch4r_buf.h"

static inline MPIR_Request *MPIDIG_request_create(MPIR_Request_kind_t kind, int ref_count)
{
    MPIR_Request *req;
    int i;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDIG_REQUEST_CREATE);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDIG_REQUEST_CREATE);

    req = MPIR_Request_create(kind, 0);
    if (req == NULL)
        goto fn_fail;

    /* as long as ref_count is a constant, any compiler should be able
     * to unroll the below loop.  when threading is not enabled, the
     * compiler should be able to combine the below individual
     * increments to a single increment of "ref_count - 1".
     *
     * FIXME: when threading is enabled, the ref_count increase is an
     * atomic operation, so it might be more inefficient.  we should
     * use a new API to increase the ref_count value instead of the
     * for loop. */
    for (i = 0; i < ref_count - 1; i++)
        MPIR_Request_add_ref(req);

    MPIDI_NM_am_request_init(req);
#ifndef MPIDI_CH4_DIRECT_NETMOD
    MPIDI_SHM_am_request_init(req);
#endif

    MPL_COMPILE_TIME_ASSERT(sizeof(MPIDIG_req_ext_t) <= MPIDIU_BUF_POOL_SZ);
    MPIDIG_REQUEST(req, req) = (MPIDIG_req_ext_t *) MPIDIU_get_buf(MPIDI_global.buf_pool);
    MPIR_Assert(MPIDIG_REQUEST(req, req));
    MPIDIG_REQUEST(req, req->status) = 0;

  fn_exit:
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDIG_REQUEST_CREATE);
    return req;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX MPIR_Request *MPIDIG_request_init(MPIR_Request * req,
                                                           MPIR_Request_kind_t kind)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDIG_REQUEST_INIT);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDIG_REQUEST_INIT);

    MPIR_Assert(req != NULL);
    /* Increment the refcount by one to account for the MPIDIG layer */
    MPIR_Request_add_ref(req);

    MPIDI_NM_am_request_init(req);
#ifndef MPIDI_CH4_DIRECT_NETMOD
    MPIDI_SHM_am_request_init(req);
#endif

    MPL_COMPILE_TIME_ASSERT(sizeof(MPIDIG_req_ext_t) <= MPIDIU_BUF_POOL_SZ);
    MPIDIG_REQUEST(req, req) = (MPIDIG_req_ext_t *) MPIDIU_get_buf(MPIDI_global.buf_pool);
    MPIR_Assert(MPIDIG_REQUEST(req, req));
    MPIDIG_REQUEST(req, req->status) = 0;

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDIG_REQUEST_INIT);
    return req;
}

MPL_STATIC_INLINE_PREFIX void MPIDIG_request_finalize(MPIR_Request * req)
{
    if (MPIDIG_REQUEST(req, req)) {
        MPIDIU_release_buf(MPIDIG_REQUEST(req, req));
        MPIDIG_REQUEST(req, req) = NULL;
        MPIDI_NM_am_request_finalize(req);
#ifndef MPIDI_CH4_DIRECT_NETMOD
        MPIDI_SHM_am_request_finalize(req);
#endif
    }
}

/* This function should be called any time an anysource request is matched so
 * the upper layer will have a chance to arbitrate who wins the race between
 * the netmod and the shmod. This will cancel the request of the other side and
 * take care of copying any relevant data. */
static inline int MPIDI_anysource_matched(MPIR_Request * rreq, int caller, int *continue_matching)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_ANYSOURCE_MATCHED);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_ANYSOURCE_MATCHED);

    MPIR_Assert(MPIDI_NETMOD == caller || MPIDI_SHM == caller);

    if (MPIDI_NETMOD == caller) {
#ifndef MPIDI_CH4_DIRECT_NETMOD
        int c;
        /* MPIDI_SHM_mpi_cancel_recv may complete the request, but we do not
         * want it to be visible until we copy the request status below.
         * This is important when an asynchronous progress is on, because
         * a user thread may immediately see the (incomplete) status.
         * Bump the request refcount by one to prevent request completion inside
         * MPIDI_SHM_mpi_cancel_recv. */
        MPIR_cc_incr(rreq->cc_ptr, &c);
        mpi_errno = MPIDI_SHM_mpi_cancel_recv(rreq);

        /* If the netmod is cancelling the request, then shared memory will
         * just copy the status from the shared memory side because the netmod
         * will always win the race condition here. */
        if (MPIR_STATUS_GET_CANCEL_BIT(rreq->status)) {
            /* If the request is cancelled, copy the status object from the
             * partner request */
            rreq->status = MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq)->status;
        }
        /* Now it's safe to complete the request */
        MPID_Request_complete(rreq);
#endif
        *continue_matching = 0;
    } else if (MPIDI_SHM == caller) {
        mpi_errno = MPIDI_NM_mpi_cancel_recv(rreq);

        /* If the netmod has already matched this request, shared memory will
         * lose and should stop matching this request */
        *continue_matching = !MPIR_STATUS_GET_CANCEL_BIT(rreq->status);
    }

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_ANYSOURCE_MATCHED);
    return mpi_errno;
}

#ifndef MPIDI_CH4_DIRECT_NETMOD
/* mpi-layer request handle is connected to one of the request (the one with shm),
 * so we can't simply free one at matching time. Both need be completed and freed
 * at the same time.
 *
 * At matching time, we cancel one of them (but not freed).
 * At completion time, we complete and free both.
 */
static inline int MPIDI_anysrc_try_cancel(MPIR_Request * rreq, int *is_cancelled)
{
    int mpi_errno = MPI_SUCCESS;

    *is_cancelled = 1;  /* overwrite if cancel fails */

    MPIR_Request *anysrc_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq);
    if (unlikely(anysrc_partner)) {
        if (!MPIR_STATUS_GET_CANCEL_BIT(anysrc_partner->status)) {
            if (MPIDI_REQUEST(rreq, is_local)) {
                /* SHM, cancel NM partner */
                mpi_errno = MPIDI_NM_mpi_cancel_recv(anysrc_partner);
                MPIR_ERR_CHECK(mpi_errno);
                if (!MPIR_STATUS_GET_CANCEL_BIT(anysrc_partner->status)) {
                    /* failed, cancel SHM rreq instead
                     * NOTE: comm and datatype will be defreferenced at caller site
                     */
                    MPIR_STATUS_SET_CANCEL_BIT(rreq->status, TRUE);
                    *is_cancelled = 0;
                } else {
                    /* NM partner is not user-visible, free it now, which means
                     * explicit SHM code can skip MPIDI_anysrc_free_partner
                     */
                    MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq) = NULL;
                    MPIDI_REQUEST_ANYSOURCE_PARTNER(anysrc_partner) = NULL;
                    /* cancel freed it once, freed once more on behalf of mpi-layer */
                    MPIR_Request_free(anysrc_partner);
                }
            } else {
                /* NM, cancel SHM partner */
                int c;
                /* prevent free, we'll complete it separately */
                MPIR_cc_incr(anysrc_partner->cc_ptr, &c);
                mpi_errno = MPIDI_SHM_mpi_cancel_recv(anysrc_partner);
                MPIR_ERR_CHECK(mpi_errno);

                /* NOTE: We assume the cancel is always successful.
                 *       This is because SHM won't match unless NM cancels.
                 *       We only reach here when NM is still active.
                 */
            }
        }
    }

  fn_exit:
    return MPI_SUCCESS;
  fn_fail:
    goto fn_exit;
}

static inline void MPIDI_anysrc_free_partner(MPIR_Request * rreq)
{
    MPIR_Request *anysrc_partner = MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq);
    if (unlikely(anysrc_partner)) {
        if (!MPIDI_REQUEST(rreq, is_local)) {
            /* NM, complete and free SHM partner */
            MPIDI_REQUEST_ANYSOURCE_PARTNER(rreq) = NULL;
            MPIDI_REQUEST_ANYSOURCE_PARTNER(anysrc_partner) = NULL;
            MPID_Request_complete(anysrc_partner);
            /* copy status to SHM partner */
            anysrc_partner->status = rreq->status;
            /* free NM request on behalf of mpi-layer
             * final free for NM request happen at call-site MPID_Request_complete
             * final free for SHM partner happen at mpi-layer
             */
            MPIR_Request_free(rreq);
        } else {
            /* SHM, NM partner should already been freed (this branch can't happen) */
            MPIR_Assert(0);
        }
    }
}
#endif /* MPIDI_CH4_DIRECT_NETMOD */

#endif /* CH4R_REQUEST_H_INCLUDED */
