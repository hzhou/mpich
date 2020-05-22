/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "posix_noinline.h"
#include "release_gather.h"

int MPIDI_POSIX_mpi_comm_commit_pre_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

    /* Release_gather primitives based collective algorithm works for Intra-comms only */
    if (comm->comm_kind == MPIR_COMM_KIND__INTRACOMM) {
        mpi_errno = MPIDI_POSIX_mpi_release_gather_comm_init_null(comm);
        MPIR_ERR_CHECK(mpi_errno);
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_POSIX_mpi_comm_commit_post_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;


    /* prune selection tree */
    mpi_errno = MPIR_Csel_prune(MPIDI_global.shm.posix.csel_root, comm,
                                &MPIDI_POSIX_COMM(comm, csel_comm));
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_POSIX_mpi_comm_free_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

    if (MPIDI_POSIX_COMM(comm, csel_comm)) {
        mpi_errno = MPIR_Csel_free(MPIDI_POSIX_COMM(comm, csel_comm));
        MPIR_ERR_CHECK(mpi_errno);
    }

    /* Release_gather primitives based collective algorithm works for Intra-comms only */
    if (comm->comm_kind == MPIR_COMM_KIND__INTRACOMM) {
        mpi_errno = MPIDI_POSIX_mpi_release_gather_comm_free(comm);
        MPIR_ERR_CHECK(mpi_errno);
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}
