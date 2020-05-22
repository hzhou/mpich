/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"
#include "typerep_internal.h"
#include "yaksa.h"

int MPIR_Typerep_create_vector(int count, int blocklength, int stride, MPI_Datatype oldtype,
                               void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_vector(count, blocklength, stride, type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_hvector(int count, int blocklength, MPI_Aint stride, MPI_Datatype oldtype,
                                void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_hvector(count, blocklength, stride, type,
                                       (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_contig(int count, MPI_Datatype oldtype, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_contig(count, type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_dup(MPI_Datatype oldtype, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_dup(type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_indexed_block(int count, int blocklength, const int *array_of_displacements,
                                      MPI_Datatype oldtype, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_indexed_block(count, blocklength, array_of_displacements,
                                             type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_hindexed_block(int count, int blocklength,
                                       const MPI_Aint * array_of_displacements,
                                       MPI_Datatype oldtype, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_hindexed_block(count, blocklength, array_of_displacements,
                                              type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_indexed(int count, const int *array_of_blocklengths,
                                const int *array_of_displacements, MPI_Datatype oldtype,
                                void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_indexed(count, array_of_blocklengths, array_of_displacements,
                                       type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_hindexed(int count, const int *array_of_blocklengths,
                                 const MPI_Aint * array_of_displacements, MPI_Datatype oldtype,
                                 void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_hindexed(count, array_of_blocklengths, array_of_displacements,
                                        type, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_resized(MPI_Datatype oldtype, MPI_Aint lb, MPI_Aint extent, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;

    yaksa_type_t type = MPII_Typerep_get_yaksa_type(oldtype);

    int rc = yaksa_type_create_resized(type, lb, extent, (yaksa_type_t *) typerep);
    MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_struct(int count, const int *array_of_blocklengths,
                               const MPI_Aint * array_of_displacements,
                               const MPI_Datatype * array_of_types, void **typerep)
{

    int mpi_errno = MPI_SUCCESS;
    yaksa_type_t *array_of_yaksa_types = (yaksa_type_t *) MPL_malloc(count * sizeof(yaksa_type_t),
                                                                     MPL_MEM_DATATYPE);
    int *array_of_real_blocklengths = (int *) MPL_malloc(count * sizeof(int), MPL_MEM_DATATYPE);
    intptr_t *array_of_real_displacements = (intptr_t *) MPL_malloc(count * sizeof(intptr_t),
                                                                    MPL_MEM_DATATYPE);


    /* MPI_LB and MPI_UB were in the MPI-1 standard and were
     * deprecated in MPI-2.  We should probably stop supporting them
     * in the future.  For now, we can simply convert them to a
     * resized type. */
    int real_count = 0;
    for (int i = 0; i < count; i++) {
        if (array_of_types[i] != MPI_LB && array_of_types[i] != MPI_UB) {
            array_of_yaksa_types[real_count] = MPII_Typerep_get_yaksa_type(array_of_types[i]);
            array_of_real_blocklengths[real_count] = array_of_blocklengths[i];
            array_of_real_displacements[real_count] = (intptr_t) array_of_displacements[i];
            real_count++;
        }
    }

    if (count == real_count) {
        int rc = yaksa_type_create_struct(count, array_of_blocklengths, array_of_real_displacements,
                                          array_of_yaksa_types, (yaksa_type_t *) typerep);
        MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");
    } else {
        yaksa_type_t tmp;

        int rc = yaksa_type_create_struct(real_count, array_of_real_blocklengths,
                                          array_of_real_displacements,
                                          array_of_yaksa_types, &tmp);
        MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

        intptr_t lb, ub;
        uintptr_t extent;
        rc = yaksa_type_get_extent(tmp, &lb, &extent);
        MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");
        ub = lb + extent;

        /* see if the user set MPI_LB/UB to override these values */
        for (int i = 0; i < count; i++) {
            if (array_of_types[i] == MPI_LB)
                lb = array_of_displacements[i];
            if (array_of_types[i] == MPI_UB)
                ub = array_of_displacements[i];
        }

        rc = yaksa_type_create_resized(tmp, lb, ub - lb, (yaksa_type_t *) typerep);
        MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");

        rc = yaksa_type_free(tmp);
        MPIR_ERR_CHKANDJUMP(rc, mpi_errno, MPI_ERR_INTERN, "**yaksa");
    }

    MPL_free(array_of_real_displacements);
    MPL_free(array_of_real_blocklengths);
    MPL_free(array_of_yaksa_types);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

int MPIR_Typerep_create_subarray(int ndims, const int *array_of_sizes, const int *array_of_subsizes,
                                 const int *array_of_starts, int order,
                                 MPI_Datatype oldtype, void **typerep)
{

    /* MPICH breaks down subarrays into smaller types, so we don't
     * need to use yaksa subarray types */
    return MPI_SUCCESS;
}

int MPIR_Typerep_create_darray(int size, int rank, int ndims, const int *array_of_gsizes,
                               const int *array_of_distribs, const int *array_of_dargs,
                               const int *array_of_psizes, int order, MPI_Datatype oldtype,
                               void **typerep)
{

    return MPI_SUCCESS;
}
