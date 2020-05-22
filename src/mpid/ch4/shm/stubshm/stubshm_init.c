/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "stubshm_impl.h"

int MPIDI_STUBSHM_mpi_init_hook(int rank, int size)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_vci_attr(int vci)
{
    int ret = 0;


    MPIR_Assert(0);

    return ret;
}

int MPIDI_STUBSHM_mpi_finalize_hook(void)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

void *MPIDI_STUBSHM_mpi_alloc_mem(size_t size, MPIR_Info * info_ptr)
{

    MPIR_Assert(0);

    return NULL;
}

int MPIDI_STUBSHM_mpi_free_mem(void *ptr)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_node_id(MPIR_Comm * comm, int rank, int *id_p)
{

    *id_p = 0;

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_max_node_id(MPIR_Comm * comm, int *max_id_p)
{

    *max_id_p = 0;

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_local_upids(MPIR_Comm * comm, size_t ** local_upid_size, char **local_upids)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_upids_to_lupids(int size, size_t * remote_upid_size, char *remote_upids,
                                  int **remote_lupids)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_create_intercomm_from_lpids(MPIR_Comm * newcomm_ptr, int size, const int lpids[])
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_mpi_type_commit_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_type_free_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_op_commit_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_op_free_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;

    return mpi_errno;
}
