/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2011 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_T_category_get_pvars */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_T_category_get_pvars = PMPI_T_category_get_pvars
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_T_category_get_pvars  MPI_T_category_get_pvars
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_T_category_get_pvars as PMPI_T_category_get_pvars
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_T_category_get_pvars(int cat_index, int len, int indices[])
    __attribute__ ((weak, alias("PMPI_T_category_get_pvars")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_T_category_get_pvars
#define MPI_T_category_get_pvars PMPI_T_category_get_pvars

/* any non-MPI functions go here, especially non-static ones */

int MPIR_T_category_get_pvars_impl(int cat_index, int len, int indices[])
{
    int mpi_errno = MPI_SUCCESS;
    cat_table_entry_t *cat;
    int i, num_pvars, count;

    cat = (cat_table_entry_t *) utarray_eltptr(cat_table, cat_index);
    num_pvars = utarray_len(cat->pvar_indices);
    count = len < num_pvars ? len : num_pvars;

    for (i = 0; i < count; i++) {
        indices[i] = *(int *) utarray_eltptr(cat->pvar_indices, i);
    }

    return mpi_errno;
}

#endif /* MPICH_MPI_FROM_PMPI */

/*@
MPI_T_category_get_pvars - Get performance variables in a category

Input Parameters:
+ cat_index - index of the category to be queried, in the range [0,N-1] (integer)
- len - the length of the indices array (integer)

Output Parameters:
. indices - an integer array of size len, indicating performance variable indices (array of integers)

.N ThreadSafe

.N Errors
.N MPI_SUCCESS
.N MPI_T_ERR_NOT_INITIALIZED
.N MPI_T_ERR_INVALID_INDEX
@*/
int MPI_T_category_get_pvars(int cat_index, int len, int indices[])
{
    int mpi_errno = MPI_SUCCESS;


    MPIR_ERRTEST_MPIT_INITIALIZED(mpi_errno);
    MPIR_T_THREAD_CS_ENTER();


    /* Validate parameters */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_CAT_INDEX(cat_index, mpi_errno);
            if (len != 0)
                MPIR_ERRTEST_ARGNULL(indices, "indices", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif /* HAVE_ERROR_CHECKING */

    /* ... body of routine ...  */

    if (len == 0)
        goto fn_exit;

    mpi_errno = MPIR_T_category_get_pvars_impl(cat_index, len, indices);
    MPIR_ERR_CHECK(mpi_errno);

    /* ... end of body of routine ... */

  fn_exit:

    MPIR_T_THREAD_CS_EXIT();
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#ifdef HAVE_ERROR_CHECKING
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_t_category_get_pvars",
                                 "**mpi_t_category_get_pvars %d %d %p", cat_index, len, indices);
    }
#endif
    mpi_errno = MPIR_Err_return_comm(NULL, __func__, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}
