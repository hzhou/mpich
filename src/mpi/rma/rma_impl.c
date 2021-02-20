/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"

int MPIR_Win_shared_query_impl(MPIR_Win * win_ptr, int rank, MPI_Aint * size, MPI_Aint * disp_unit,
                               void *baseptr)
{
    return MPID_Win_shared_query(win_ptr, rank, size, disp_unit, baseptr);
}
