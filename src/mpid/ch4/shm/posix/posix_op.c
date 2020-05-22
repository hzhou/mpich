/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "posix_noinline.h"

int MPIDI_POSIX_mpi_op_commit_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}

int MPIDI_POSIX_mpi_op_free_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}
