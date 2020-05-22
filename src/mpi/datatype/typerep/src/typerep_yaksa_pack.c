/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"
#include "yaksa.h"
#include "typerep_internal.h"

int MPIR_Typerep_pack(const void *inbuf, MPI_Aint incount, MPI_Datatype datatype,
                      MPI_Aint inoffset, void *outbuf, MPI_Aint max_pack_bytes,
                      MPI_Aint * actual_pack_bytes)
{

    int mpi_errno = MPI_SUCCESS;
    int rc;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(datatype);

    yaksa_request_t request;
    uintptr_t real_pack_bytes;
    rc = yaksa_ipack(inbuf, incount, type, inoffset, outbuf, max_pack_bytes, &real_pack_bytes,
                     NULL, &request);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

    rc = yaksa_request_wait(request);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

    *actual_pack_bytes = (MPI_Aint) real_pack_bytes;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_unpack(const void *inbuf, MPI_Aint insize, void *outbuf, MPI_Aint outcount,
                        MPI_Datatype datatype, MPI_Aint outoffset, MPI_Aint * actual_unpack_bytes)
{

    int mpi_errno = MPI_SUCCESS;
    int rc;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(datatype);

    uintptr_t size;
    rc = yaksa_type_get_size(type, &size);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

    uintptr_t real_insize = MPL_MIN(insize, size * outcount);

    yaksa_request_t request;
    uintptr_t real_unpack_bytes;
    rc = yaksa_iunpack(inbuf, real_insize, outbuf, outcount, type, outoffset, &real_unpack_bytes,
                       NULL, &request);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

    rc = yaksa_request_wait(request);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

    *actual_unpack_bytes = (MPI_Aint) real_unpack_bytes;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}
