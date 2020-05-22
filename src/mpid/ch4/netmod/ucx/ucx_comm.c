/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "ucx_impl.h"
#include "mpidu_bc.h"
#ifdef HAVE_LIBHCOLL
#include "../../common/hcoll/hcoll.h"
#endif

int MPIDI_UCX_mpi_comm_commit_pre_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

#if defined HAVE_LIBHCOLL
    hcoll_comm_create(comm, NULL);
#endif

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_UCX_mpi_comm_commit_post_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_UCX_mpi_comm_free_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

#ifdef HAVE_LIBHCOLL
    hcoll_comm_destroy(comm, NULL);
#endif
  fn_exit:
    return mpi_errno;
}
