/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpioimpl.h"

#ifdef HAVE_WEAK_SYMBOLS

#if defined(HAVE_PRAGMA_WEAK)
#pragma weak MPI_File_open = PMPI_File_open
#elif defined(HAVE_PRAGMA_HP_SEC_DEF)
#pragma _HP_SECONDARY_DEF PMPI_File_open MPI_File_open
#elif defined(HAVE_PRAGMA_CRI_DUP)
#pragma _CRI duplicate MPI_File_open as PMPI_File_open
/* end of weak pragmas */
#elif defined(HAVE_WEAK_ATTRIBUTE)
int MPI_File_open(MPI_Comm comm, const char *filename, int amode, MPI_Info info, MPI_File * fh)
    __attribute__ ((weak, alias("PMPI_File_open")));
#endif

/* Include mapping from MPI->PMPI */
#define MPIO_BUILD_PROFILING
#include "mpioprof.h"
#endif

/* for user-definde reduce operator */
#include "adio_extern.h"

/* */
#if !defined(HAVE_WEAK_SYMBOLS) && !defined(MPIO_BUILD_PROFILING)
void *dummy_refs_MPI_File_open[] = {
    (void *) ADIO_ImmediateOpen,
    (void *) MPIU_datatype_full_size,
};
#endif

extern int ADIO_Init_keyval;

/*@
    MPI_File_open - Opens a file

Input Parameters:
. comm - communicator (handle)
. filename - name of file to open (string)
. amode - file access mode (integer)
. info - info object (handle)

Output Parameters:
. fh - file handle (handle)

.N fortran
@*/
int MPI_File_open(MPI_Comm comm, ROMIO_CONST char *filename, int amode,
                  MPI_Info info, MPI_File * fh)
{
    int error_code = MPI_SUCCESS, flag;
    static char myname[] = "MPI_FILE_OPEN";
#ifdef MPI_hpux
    int fl_xmpi;

    HPMP_IO_OPEN_START(fl_xmpi, comm);
#endif /* MPI_hpux */

    ROMIO_THREAD_CS_ENTER();

    /* --BEGIN ERROR HANDLING-- */
    MPIO_CHECK_COMM(comm, myname, error_code);
    MPIO_CHECK_INFO_ALL(info, error_code, comm);
    /* --END ERROR HANDLING-- */

    error_code = MPI_Comm_test_inter(comm, &flag);
    /* --BEGIN ERROR HANDLING-- */
    if (error_code || flag) {
        error_code = MPIO_Err_create_code(error_code, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__, MPI_ERR_COMM, "**commnotintra", 0);
        goto fn_fail;
    }

    if (((amode & MPI_MODE_RDONLY) ? 1 : 0) + ((amode & MPI_MODE_RDWR) ? 1 : 0) +
        ((amode & MPI_MODE_WRONLY) ? 1 : 0) != 1) {
        error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__, MPI_ERR_AMODE, "**fileamodeone", 0);
        goto fn_fail;
    }

    if ((amode & MPI_MODE_RDONLY) && ((amode & MPI_MODE_CREATE) || (amode & MPI_MODE_EXCL))) {
        error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__, MPI_ERR_AMODE, "**fileamoderead", 0);
        goto fn_fail;
    }

    if ((amode & MPI_MODE_RDWR) && (amode & MPI_MODE_SEQUENTIAL)) {
        error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__, MPI_ERR_AMODE, "**fileamodeseq", 0);
        goto fn_fail;
    }

    error_code = MPIR_File_open_impl(comm, filename, amode, info, fh);
    if (error_code != MPI_SUCCESS)
        goto fn_fail;

#ifdef MPI_hpux
    HPMP_IO_OPEN_END(fl_xmpi, *fh, comm);
#endif /* MPI_hpux */

  fn_exit:
    ROMIO_THREAD_CS_EXIT();
    return error_code;
  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
    error_code = MPIO_Err_return_file(MPI_FILE_NULL, error_code);
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}
