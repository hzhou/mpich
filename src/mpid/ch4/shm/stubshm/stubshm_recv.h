/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef STUBSHM_RECV_H_INCLUDED
#define STUBSHM_RECV_H_INCLUDED

#include "stubshm_impl.h"

static inline int MPIDI_STUBSHM_mpi_recv(void *buf,
                                         MPI_Aint count,
                                         MPI_Datatype datatype,
                                         int rank,
                                         int tag,
                                         MPIR_Comm * comm,
                                         int context_offset, MPI_Status * status,
                                         MPIR_Request ** request)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_mpi_imrecv(void *buf,
                                           MPI_Aint count,
                                           MPI_Datatype datatype,
                                           MPIR_Request * message, MPIR_Request ** rreqp)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_mpi_irecv(void *buf,
                                          MPI_Aint count,
                                          MPI_Datatype datatype,
                                          int rank,
                                          int tag,
                                          MPIR_Comm * comm, int context_offset,
                                          MPIR_Request ** request)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_mpi_cancel_recv(MPIR_Request * rreq)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

#endif /* STUBSHM_RECV_H_INCLUDED */
