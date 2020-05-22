/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef UCX_AM_RECV_H_INCLUDED
#define UCX_AM_RECV_H_INCLUDED

#include "ucx_impl.h"

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_mpi_recv(void *buf,
                                               MPI_Aint count,
                                               MPI_Datatype datatype,
                                               int rank,
                                               int tag,
                                               MPIR_Comm * comm,
                                               int context_offset, MPIDI_av_entry_t * addr,
                                               MPI_Status * status, MPIR_Request ** request)
{
    int ret;


    ret = MPIDIG_mpi_recv(buf, count, datatype, rank, tag, comm, context_offset, status, request);

    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_mpi_imrecv(void *buf,
                                                 MPI_Aint count, MPI_Datatype datatype,
                                                 MPIR_Request * message)
{
    int ret;


    ret = MPIDIG_mpi_imrecv(buf, count, datatype, message);

    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_mpi_irecv(void *buf,
                                                MPI_Aint count,
                                                MPI_Datatype datatype,
                                                int rank,
                                                int tag,
                                                MPIR_Comm * comm, int context_offset,
                                                MPIDI_av_entry_t * addr, MPIR_Request ** request)
{
    int ret;


    ret = MPIDIG_mpi_irecv(buf, count, datatype, rank, tag, comm, context_offset, request);

    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_mpi_cancel_recv(MPIR_Request * rreq)
{
    int ret;


    ret = MPIDIG_mpi_cancel_recv(rreq);

    return ret;
}

#endif /* UCX_AM_RECV_H_INCLUDED */
