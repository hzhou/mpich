/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "stubshm_impl.h"

int MPIDI_STUBSHM_mpi_comm_commit_pre_hook(MPIR_Comm * comm)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_mpi_comm_commit_post_hook(MPIR_Comm * comm)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_mpi_comm_free_hook(MPIR_Comm * comm)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}
