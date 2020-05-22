/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "ofi_impl.h"
#include "ofi_noinline.h"

int MPIDI_OFI_mpi_op_free_hook(MPIR_Op * op_p)
{

    return 0;
}

int MPIDI_OFI_mpi_op_commit_hook(MPIR_Op * op_p)
{

    return 0;
}
