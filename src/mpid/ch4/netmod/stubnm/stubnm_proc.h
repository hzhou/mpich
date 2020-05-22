/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef STUBNM_PROC_H_INCLUDED
#define STUBNM_PROC_H_INCLUDED

#include "stubnm_impl.h"

static inline int MPIDI_NM_rank_is_local(int rank, MPIR_Comm * comm)
{
    int ret;

    MPIR_Assert(0);
    ret = 0;

    return ret;
}

static inline int MPIDI_NM_av_is_local(MPIDI_av_entry_t * av)
{
    int ret;

    MPIR_Assert(0);
    ret = 0;

    return ret;
}

static inline int MPIDI_NM_comm_get_lpid(MPIR_Comm * comm_ptr, int idx, int *lpid_ptr,
                                         bool is_remote)
{
    MPIR_Assert(0);
    return MPI_SUCCESS;
}

#endif /* STUBNM_PROC_H_INCLUDED */
