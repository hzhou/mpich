/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpioimpl.h"

int MPIR_File_open_impl(MPI_Comm comm, ROMIO_CONST char *filename, int amode,
                        MPI_Info info, MPI_File * fh)
{
    int error_code = MPI_SUCCESS;
    MPI_Comm dupcomm = MPI_COMM_NULL;

    MPI_Comm_dup(comm, &dupcomm);

    /* Check if ADIO has been initialized. If not, initialize it */
    MPIR_MPIOInit(&error_code);
    if (error_code != MPI_SUCCESS)
        goto fn_fail;

    /* check if amode is the same on all processes: at first glance, one might try
     * to use a built-in operator like MPI_BAND, but we need every mpi process to
     * agree the amode was not the same.  Consider process A with
     * MPI_MODE_CREATE|MPI_MODE_RDWR, and B with MPI_MODE_RDWR:  MPI_BAND yields
     * MPI_MODE_RDWR.  A determines amodes are different, but B proceeds having not
     * detected an error */
    int tmp_amode = 0;
    MPI_Allreduce(&amode, &tmp_amode, 1, MPI_INT, ADIO_same_amode, dupcomm);

    if (tmp_amode == ADIO_AMODE_NOMATCH) {
        error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__, MPI_ERR_NOT_SAME, "**fileamodediff", 0);
        goto fn_fail;
    }

    int file_system = -1, known_fstype;
    ADIOI_Fns *fsops;
    /* resolve file system type from file name; this is a collective call */
    known_fstype = ADIO_ResolveFileType(dupcomm, filename, &file_system, &fsops, &error_code);
    if (error_code != MPI_SUCCESS) {
        /* ADIO_ResolveFileType() will print as informative a message as it
         * possibly can or call MPIO_Err_setmsg.  We just need to propagate
         * the error up.
         */
        goto fn_fail;
    }

    if (known_fstype) {
        /* filename contains a known file system type prefix, such as "ufs:".
         * strip off prefix if there is one, but only skip prefixes
         * if they are greater than length one to allow for windows
         * drive specifications (e.g. c:\...)
         */
        char *tmp = strchr(filename, ':');
        if (tmp > filename + 1) {
            filename = tmp + 1;
        }
    }

    /* use default values for disp, etype, filetype */
    *fh = ADIO_Open(comm, dupcomm, filename, file_system, fsops, amode, 0,
                    MPI_BYTE, MPI_BYTE, info, ADIO_PERM_NULL, &error_code);

    if (error_code != MPI_SUCCESS) {
        goto fn_fail;
    }

    /* if MPI_MODE_SEQUENTIAL requested, file systems cannot do explicit offset
     * or independent file pointer accesses, leaving not much else aside from
     * shared file pointer accesses. */
    if (!ADIO_Feature((*fh), ADIO_SHARED_FP) && (amode & MPI_MODE_SEQUENTIAL)) {
        error_code = MPIO_Err_create_code(MPI_SUCCESS, MPIR_ERR_RECOVERABLE,
                                          myname, __LINE__,
                                          MPI_ERR_UNSUPPORTED_OPERATION, "**iosequnsupported", 0);
        ADIO_Close(*fh, &error_code);
        goto fn_fail;
    }

    /* determine name of file that will hold the shared file pointer */
    /* can't support shared file pointers on a file system that doesn't
     * support file locking. */
    if ((error_code == MPI_SUCCESS) && ADIO_Feature((*fh), ADIO_SHARED_FP)) {
        int rank;
        MPI_Comm_rank(dupcomm, &rank);
        ADIOI_Shfp_fname(*fh, rank, &error_code);
        if (error_code != MPI_SUCCESS)
            goto fn_fail;

        /* if MPI_MODE_APPEND, set the shared file pointer to end of file.
         * indiv. file pointer already set to end of file in ADIO_Open.
         * Here file view is just bytes. */
        if ((*fh)->access_mode & MPI_MODE_APPEND) {
            if (rank == (*fh)->hints->ranklist[0])      /* only one person need set the sharedfp */
                ADIO_Set_shared_fp(*fh, (*fh)->fp_ind, &error_code);
            MPI_Barrier(dupcomm);
        }
    }

  fn_exit:
    return error_code;
  fn_fail:
    if (dupcomm != MPI_COMM_NULL)
        MPI_Comm_free(&dupcomm);
    goto fn_exit;
}
