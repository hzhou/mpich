/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Free_mem */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Free_mem = PMPI_Free_mem
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Free_mem  MPI_Free_mem
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Free_mem as PMPI_Free_mem
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Free_mem(void *base) __attribute__ ((weak, alias("PMPI_Free_mem")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Free_mem
#define MPI_Free_mem PMPI_Free_mem

#endif

/*@
   MPI_Free_mem - Free memory allocated with MPI_Alloc_mem

Input Parameters:
.  base - initial address of memory segment allocated by 'MPI_ALLOC_MEM'
       (choice)

.N ThreadSafe

.N Fortran

.N Errors
.N MPI_SUCCESS
@*/
int MPI_Free_mem(void *base)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_ERRTEST_INITIALIZED_ORDIE();

    MPID_THREAD_CS_ENTER(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);

    /* ... body of routine ...  */

    if (base == NULL)
        goto fn_exit;

    mpi_errno = MPID_Free_mem(base);
    if (mpi_errno)
        goto fn_fail;

    /* ... end of body of routine ... */

  fn_exit:
    MPID_THREAD_CS_EXIT(GLOBAL, MPIR_THREAD_GLOBAL_ALLFUNC_MUTEX);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#ifdef HAVE_ERROR_CHECKING
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_free_mem", "**mpi_free_mem %p", base);
    }
#endif
    /* MPI_Free_mem must invoke the error handler on MPI_COMM_WORLD if there
     * is an error */
    mpi_errno = MPIR_Err_return_comm(NULL, __func__, mpi_errno);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}
