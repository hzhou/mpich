/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "stubshm_impl.h"

int MPIDI_STUBSHM_do_progress_recv(int blocking, int *completion_count)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_do_progress_send(int blocking, int *completion_count)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress(int vci, int blocking)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_test(void)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_poke(void)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

void MPIDI_STUBSHM_progress_start(MPID_Progress_state * state)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

void MPIDI_STUBSHM_progress_end(MPID_Progress_state * state)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_wait(MPID_Progress_state * state)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_register(int (*progress_fn) (int *))
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_deregister(int id)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_activate(int id)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}

int MPIDI_STUBSHM_progress_deactivate(int id)
{

    MPIR_Assert(0);

    return MPI_SUCCESS;
}
