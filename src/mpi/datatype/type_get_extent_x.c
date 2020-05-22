/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Type_get_extent_x */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Type_get_extent_x = PMPI_Type_get_extent_x
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Type_get_extent_x  MPI_Type_get_extent_x
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Type_get_extent_x as PMPI_Type_get_extent_x
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count * lb, MPI_Count * extent)
    __attribute__ ((weak, alias("PMPI_Type_get_extent_x")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Type_get_extent_x
#define MPI_Type_get_extent_x PMPI_Type_get_extent_x

/* any non-MPI functions go here, especially non-static ones */

void MPIR_Type_get_extent_x_impl(MPI_Datatype datatype, MPI_Count * lb, MPI_Count * extent)
{
    MPIR_Datatype *datatype_ptr = NULL;

    MPIR_Datatype_get_ptr(datatype, datatype_ptr);

    if (HANDLE_IS_BUILTIN(datatype)) {
        *lb = 0;
        *extent = MPIR_Datatype_get_basic_size(datatype);
    } else {
        *lb = datatype_ptr->lb;
        *extent = datatype_ptr->extent; /* derived, should be same as ub - lb */
    }
}

#endif /* MPICH_MPI_FROM_PMPI */

/*@
MPI_Type_get_extent_x - Get the lower bound and extent as MPI_Count values
 for a Datatype

Input Parameters:
. datatype - datatype (handle)

Output Parameters:
+ lb - lower bound of datatype (integer)
- extent - extent of datatype (integer)

.N ThreadSafe

.N Fortran

.N Errors
@*/
int MPI_Type_get_extent_x(MPI_Datatype datatype, MPI_Count * lb, MPI_Count * extent)
{
    int mpi_errno = MPI_SUCCESS;

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

    /* Validate parameters, especially handles needing to be converted */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_DATATYPE(datatype, "datatype", mpi_errno);

            /* TODO more checks may be appropriate */
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* Convert MPI object handles to object pointers */

    /* Validate parameters and objects (post conversion) */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            if (!HANDLE_IS_BUILTIN(datatype)) {
                MPIR_Datatype *datatype_ptr = NULL;
                MPIR_Datatype_get_ptr(datatype, datatype_ptr);
                MPIR_Datatype_valid_ptr(datatype_ptr, mpi_errno);
            }

            /* TODO more checks may be appropriate (counts, in_place, buffer aliasing, etc) */
            if (mpi_errno != MPI_SUCCESS)
                goto fn_fail;
            MPIR_ERRTEST_ARGNULL(extent, "extent", mpi_errno);
            MPIR_ERRTEST_ARGNULL(lb, "lb", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ...  */

    MPIR_Type_get_extent_x_impl(datatype, lb, extent);

    /* ... end of body of routine ... */

#ifdef HAVE_ERROR_CHECKING
  fn_exit:
#endif
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

#ifdef HAVE_ERROR_CHECKING
  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_type_get_extent_x", "**mpi_type_get_extent_x %D %p %p",
                                 datatype, lb, extent);
    }
    mpi_errno = MPIR_Err_return_comm(NULL, __func__, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
#endif
}
