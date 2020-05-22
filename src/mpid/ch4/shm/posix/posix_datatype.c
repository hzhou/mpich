/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "posix_noinline.h"

int MPIDI_POSIX_mpi_type_commit_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}

int MPIDI_POSIX_mpi_type_free_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}
