/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef STUBSHM_PROC_H_INCLUDED
#define STUBSHM_PROC_H_INCLUDED

#include "stubshm_impl.h"

static inline int MPIDI_STUBSHM_rank_is_local(int rank, MPIR_Comm * comm)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_av_is_local(MPIDI_av_entry_t * av)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

static inline int MPIDI_STUBSHM_comm_get_lpid(MPIR_Comm * comm_ptr,
                                              int idx, int *lpid_ptr, bool is_remote)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

#endif /* STUBSHM_PROC_H_INCLUDED */
