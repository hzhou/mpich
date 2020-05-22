/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "ofi_impl.h"
#include "mpidu_bc.h"
#include "ofi_noinline.h"

int MPIDI_OFI_mpi_comm_commit_pre_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

    MPIDIU_map_create(&MPIDI_OFI_COMM(comm).huge_send_counters, MPL_MEM_COMM);
    MPIDIU_map_create(&MPIDI_OFI_COMM(comm).huge_recv_counters, MPL_MEM_COMM);

    /* no connection for non-dynamic or non-root-rank of intercomm */
    MPIDI_OFI_COMM(comm).conn_id = -1;

    /* eagain defaults to off */
    if (comm->hints[MPIR_COMM_HINT_EAGAIN] == 0) {
        comm->hints[MPIR_COMM_HINT_EAGAIN] = FALSE;
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_OFI_mpi_comm_commit_post_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIDI_OFI_mpi_comm_free_hook(MPIR_Comm * comm)
{
    int mpi_errno = MPI_SUCCESS;

    MPIDIU_map_destroy(MPIDI_OFI_COMM(comm).huge_send_counters);
    MPIDIU_map_destroy(MPIDI_OFI_COMM(comm).huge_recv_counters);

    return mpi_errno;
}
