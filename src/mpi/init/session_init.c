/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpiimpl.h"

/* -- Begin Profiling Symbol Block for routine MPI_Session_init */
#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_Session_init = PMPI_Session_init
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_Session_init  MPI_Session_init
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_Session_init as PMPI_Session_init
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_Session_init(int *argc, char ***argv) __attribute__ ((weak, alias("PMPI_Session_init")));
#endif
/* -- End Profiling Symbol Block */

/* Define MPICH_MPI_FROM_PMPI if weak symbols are not supported to build
   the MPI routines */
#ifndef MPICH_MPI_FROM_PMPI
#undef MPI_Session_init
#define MPI_Session_init PMPI_Session_init
#endif

/*@
   MPI_Session_init - Create a new MPI session.

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
int MPI_Session_init(MPI_Info info, MPI_Errhandler errhandler, MPI_Session * session)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Info *info_ptr;
    MPIR_Errhandler *errhan_ptr;
    MPIR_FUNC_TERSE_STATE_DECL(MPID_STATE_MPI_SESSION_INIT);
    MPIR_FUNC_TERSE_ENTER(MPID_STATE_MPI_SESSION_INIT);

    /* Validate parameters, especially handles needing to be converted */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            MPIR_ERRTEST_INFO_OR_NULL(info, mpi_errno);
            MPIR_ERRTEST_ARGNULL(session, "session", mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif
    MPIR_Info_get_ptr(info, info_ptr);
    MPIR_Errhandler_get_ptr(errhandler, errhan_ptr);

    /* Validate parameters and objects (post conversion) */
#ifdef HAVE_ERROR_CHECKING
    {
        MPID_BEGIN_ERROR_CHECKS;
        {
            if (!HANDLE_IS_BUILTIN(errhandler)) {
                MPIR_Errhandler_valid_ptr(errhan_ptr, mpi_errno);
            }
            MPIR_ERRTEST_ERRHANDLER(errhandler, mpi_errno);
        }
        MPID_END_ERROR_CHECKS;
    }
#endif

    /* ... body of routine ...  */

    mpi_errno = MPIR_Session_init_impl(info_ptr, errhan_ptr, session);
    if (mpi_errno) {
        goto fn_fail;
    }

    /* ... end of body of routine ... */

    MPIR_FUNC_TERSE_EXIT(MPID_STATE_MPI_SESSION_INIT);
    return mpi_errno;

  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
#ifdef HAVE_ERROR_REPORTING
    {
        mpi_errno =
            MPIR_Err_create_code(mpi_errno, MPIR_ERR_RECOVERABLE, __func__, __LINE__, MPI_ERR_OTHER,
                                 "**mpi_session_init", "**mpi_session_init %I %C %p", info,
                                 errhandler, session);
    }
#endif
    mpi_errno = MPIR_Err_return_comm(0, __func__, mpi_errno);
    return mpi_errno;
    /* --END ERROR HANDLING-- */
}
