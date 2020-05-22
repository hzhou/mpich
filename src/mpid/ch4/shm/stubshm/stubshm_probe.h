/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef STUBSHM_PROBE_H_INCLUDED
#define STUBSHM_PROBE_H_INCLUDED

#include "stubshm_impl.h"


static inline int MPIDI_STUBSHM_mpi_improbe(int source,
                                            int tag,
                                            MPIR_Comm * comm,
                                            int context_offset,
                                            int *flag, MPIR_Request ** message, MPI_Status * status)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_mpi_iprobe(int source,
                                           int tag,
                                           MPIR_Comm * comm,
                                           int context_offset, int *flag, MPI_Status * status)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

#endif /* STUBSHM_PROBE_H_INCLUDED */
