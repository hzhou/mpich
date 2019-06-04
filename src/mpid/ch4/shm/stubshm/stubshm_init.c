/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */

#include "stubshm_impl.h"

int MPIDI_STUBSHM_mpi_init_hook(int rank, int size, int *n_vcis_provided)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_INIT_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_INIT_HOOK);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_INIT_HOOK);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_vci_attr(int vci)
{
    int ret = 0;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_QUERY_VCI);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_QUERY_VCI);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_QUERY_VCI);
    return ret;
}

int MPIDI_STUBSHM_mpi_finalize_hook(void)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_FINALIZE_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_FINALIZE_HOOK);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_FINALIZE_HOOK);
    return MPI_SUCCESS;
}

void *MPIDI_STUBSHM_mpi_alloc_mem(MPI_Aint size, MPIR_Info * info_ptr)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_ALLOC_MEM);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_ALLOC_MEM);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_ALLOC_MEM);
    return NULL;
}

int MPIDI_STUBSHM_mpi_free_mem(void *ptr)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_FREE_MEM);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_FREE_MEM);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_FREE_MEM);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_node_id(MPIR_Comm * comm, int rank, int *id_p)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_GET_NODE_ID);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_GET_NODE_ID);

    *id_p = 0;

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_GET_NODE_ID);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_max_node_id(MPIR_Comm * comm, int *max_id_p)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_GET_MAX_NODE_ID);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_GET_MAX_NODE_ID);

    *max_id_p = 0;

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_GET_MAX_NODE_ID);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_get_local_upids(MPIR_Comm * comm, MPI_Aint ** local_upid_size, char **local_upids)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_GET_LOCAL_UPIDS);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_GET_LOCAL_UPIDS);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_GET_LOCAL_UPIDS);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_upids_to_lupids(int size, MPI_Aint * remote_upid_size, char *remote_upids,
                                  int **remote_lupids)
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_UPIDS_TO_LUPIDS);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_UPIDS_TO_LUPIDS);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_UPIDS_TO_LUPIDS);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_create_intercomm_from_lpids(MPIR_Comm * newcomm_ptr, int size, const int lpids[])
{
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_CREATE_INTERCOMM_FROM_LPIDS);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_CREATE_INTERCOMM_FROM_LPIDS);

    MPIR_Assert(0);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_CREATE_INTERCOMM_FROM_LPIDS);
    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_mpi_type_commit_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_CREATE_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_CREATE_HOOK);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_CREATE_HOOK);
    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_type_free_hook(MPIR_Datatype * type)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_FREE_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_FREE_HOOK);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_TYPE_FREE_HOOK);
    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_op_commit_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_OP_CREATE_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_OP_CREATE_HOOK);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_OP_CREATE_HOOK);
    return mpi_errno;
}

int MPIDI_STUBSHM_mpi_op_free_hook(MPIR_Op * op)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_STUBSHM_MPI_OP_FREE_HOOK);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_STUBSHM_MPI_OP_FREE_HOOK);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_STUBSHM_MPI_OP_FREE_HOOK);
    return mpi_errno;
}
