/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Session_finalize */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Session_finalize = PMPI_Session_finalize
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Session_finalize  MPI_Session_finalize
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Session_finalize as PMPI_Session_finalize
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Session_finalize(MPI_Session * session)
    __attribute__ ((weak, alias("PMPI_Session_finalize")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Session_finalize
#define MPI_Session_finalize PMPI_Session_finalize
#endif

/*@
   MPI_Session_finalize - Create a new MPI session.

Input Parameters:
+  info - info object to specify thread support level and MPI implementation specific resources (handle)
-  errhandler - error handler to invoke in the event that an error is encountered during this function call (handle)

Output Parameters:
. session - new session (handle)

.N ThreadSafe

.N Fortran

.N Errors
.N MPI_SUCCESS

@*/
int MPI_Session_finalize(MPI_Session * session)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Session *session_ptr;
    MPIR_FUNC_TERSE_STATE_DECL(MPID_STATE_MPI_SESSION_FINALIZE);
    MPIR_FUNC_TERSE_ENTER(MPID_STATE_MPI_SESSION_FINALIZE);

    /* Validate parameters, especially handles needing to be converted */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_ARGNULL(session, "session", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif
    MPIR_Session_get_ptr(*session, session_ptr);

    /* Validate parameters and objects (post conversion) */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_Session_valid_ptr(session_ptr, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif

    /* ... body of routine ...  */

    mpi_errno = MPIR_Session_finalize_impl(session_ptr);
    if (mpi_errno) {
        goto fn_fail;
    }
    *session = MPI_SESSION_NULL;

    /* ... end of body of routine ... */

    MPIR_FUNC_TERSE_EXIT(MPID_STATE_MPI_SESSION_FINALIZE);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#ifdef HAVE_ERROR_REPORTING
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_session_finalize", "**mpi_session_finalize %p", session);
    }
#endif
    mpi_errno = MPIR_Err_return_comm(0, __func__, mpi_errno);
    return mpi_errno;
    /* --END ERROR HANDLING-- */
}
