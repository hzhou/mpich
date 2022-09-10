/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "ad_testfs.h"
#include "adioi.h"
#include "adio_extern.h"

/* ADIOI_TESTFS_SeekIndividual()
 *
 * Implements SEEK_SET only (and doesn't test for whence type); all
 * other types of whence must be converted before calling this.
 *
 * Returns an absolute offset in bytes.  The offset passed into the call is in
 * terms of the etype relative to the filetype, so some calculations are
 * necessary.
 */
ADIO_Offset ADIOI_TESTFS_SeekIndividual(ADIO_File fd, ADIO_Offset offset,
                                        int whence, int *error_code)
{
    int myrank, nprocs;

    ADIO_Offset off;
    ADIOI_Flatlist_node *flat_file;
    int i;
    ADIO_Offset abs_off_in_filetype = 0, sum;
    int filetype_is_contig;
    MPI_Count filetype_size;
    MPI_Aint etype_size, lb, filetype_extent;

    *error_code = MPI_SUCCESS;

    MPI_Comm_size(fd->comm, &nprocs);
    MPI_Comm_rank(fd->comm, &myrank);
    FPRINTF(stdout, "[%d/%d] ADIOI_TESTFS_SeekIndividual called on %s\n",
            myrank, nprocs, fd->filename);

    ADIOI_Datatype_iscontig(fd->filetype, &filetype_is_contig);
    etype_size = (MPI_Aint) fd->etype_size;

    if (filetype_is_contig)
        off = fd->disp + etype_size * offset;
    else {
        flat_file = ADIOI_Flatten_and_find(fd->filetype);

        MPI_Type_get_extent(fd->filetype, &lb, &filetype_extent);
        MPI_Type_size_x(fd->filetype, &filetype_size);
        if (!filetype_size) {
            *error_code = MPI_SUCCESS;
            return 0;
        }

        MPI_Count n_etypes_in_filetype = filetype_size / etype_size;
        MPI_Count n_filetypes = offset / n_etypes_in_filetype;
        MPI_Count etype_in_filetype = offset % n_etypes_in_filetype;
        MPI_Count size_in_filetype = etype_in_filetype * etype_size;

        sum = 0;
        for (i = 0; i < flat_file->count; i++) {
            sum += flat_file->blocklens[i];

            if (sum > size_in_filetype) {
                abs_off_in_filetype = flat_file->indices[i] +
                    size_in_filetype - (sum - flat_file->blocklens[i]);
                break;
            }
        }

        /* abs. offset in bytes in the file */
        off = fd->disp + (ADIO_Offset) n_filetypes *(ADIO_Offset) filetype_extent +
            abs_off_in_filetype;
    }

    fd->fp_ind = off;

    return off;
}
