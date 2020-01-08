/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2018 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef CH4_STARTALL_H_INCLUDED
#define CH4_STARTALL_H_INCLUDED

#include "ch4_impl.h"

MPL_STATIC_INLINE_PREFIX int MPID_Startall(int count, MPIR_Request * requests[])
{
    int mpi_errno = MPI_SUCCESS, i;

    for (i = 0; i < count; i++) {
        MPIR_Request *const preq = requests[i];
        /* continue if the source/dest is MPI_PROC_NULL */
        if (MPIDIG_REQUEST(preq, rank) == MPI_PROC_NULL)
            continue;

        switch (MPIDI_PREQUEST(preq, p_type)) {
            case MPIDI_PTYPE_RECV:
                mpi_errno = MPID_Irecv(MPIDI_PREQUEST(preq, buffer), MPIDI_PREQUEST(preq, count),
                                       MPIDI_PREQUEST(preq, datatype), MPIDI_PREQUEST(preq, rank),
                                       MPIDI_PREQUEST(preq, tag), preq->comm,
                                       MPIDI_prequest_get_context_offset(preq),
                                       &preq->u.persist.real_request);
                break;

            case MPIDI_PTYPE_SEND:
                mpi_errno = MPID_Isend(MPIDI_PREQUEST(preq, buffer), MPIDI_PREQUEST(preq, count),
                                       MPIDI_PREQUEST(preq, datatype), MPIDI_PREQUEST(preq, rank),
                                       MPIDI_PREQUEST(preq, tag), preq->comm,
                                       MPIDI_prequest_get_context_offset(preq),
                                       &preq->u.persist.real_request);
                break;

            case MPIDI_PTYPE_SSEND:
                mpi_errno = MPID_Issend(MPIDI_PREQUEST(preq, buffer), MPIDI_PREQUEST(preq, count),
                                        MPIDI_PREQUEST(preq, datatype), MPIDI_PREQUEST(preq, rank),
                                        MPIDI_PREQUEST(preq, tag), preq->comm,
                                        MPIDI_prequest_get_context_offset(preq),
                                        &preq->u.persist.real_request);
                break;

            case MPIDI_PTYPE_BSEND:{
                    MPI_Request sreq_handle;
                    mpi_errno =
                        MPIR_Ibsend_impl(MPIDI_PREQUEST(preq, buffer), MPIDI_PREQUEST(preq, count),
                                         MPIDI_PREQUEST(preq, datatype), MPIDI_PREQUEST(preq, rank),
                                         MPIDI_PREQUEST(preq, tag), preq->comm, &sreq_handle);
                    if (mpi_errno == MPI_SUCCESS)
                        MPIR_Request_get_ptr(sreq_handle, preq->u.persist.real_request);

                    break;
                }

            default:
                mpi_errno = MPIR_Err_create_code(MPI_SUCCESS, MPIR_ERR_FATAL, __FUNCTION__,
                                                 __LINE__, MPI_ERR_INTERN, "**ch4|badreqtype",
                                                 "**ch4|badreqtype %d", MPIDI_PREQUEST(preq,
                                                                                       p_type));
        }

        if (mpi_errno == MPI_SUCCESS) {
            preq->status.MPI_ERROR = MPI_SUCCESS;

            if (MPIDI_PREQUEST(preq, p_type) == MPIDI_PTYPE_BSEND) {
                preq->cc_ptr = &preq->cc;
                MPID_Request_set_completed(preq);
            } else
                preq->cc_ptr = &preq->u.persist.real_request->cc;
        } else {
            preq->u.persist.real_request = NULL;
            preq->status.MPI_ERROR = mpi_errno;
            preq->cc_ptr = &preq->cc;
            MPID_Request_set_completed(preq);
        }
    }

    return mpi_errno;
}

#endif /* CH4_STARTALL_H_INCLUDED */
