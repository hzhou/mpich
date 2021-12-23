/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpitest.h"
#include <string.h>

#ifdef MULTI_TESTS
#define run pt2pt_sendrecv2
int run(const char *arg);
#endif

int run(const char *arg)
{
    int i, j, errs = 0;
    int rank, size;
    MPI_Datatype newtype;
    char *buf = NULL;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size < 2) {
        fprintf(stderr, "comm size must be > 1\n");
        errs++;
        goto fn_exit;
    }

    buf = malloc(64 * 129);
    if (buf == NULL) {
        fprintf(stderr, "error allocating buffer\n");
        errs++;
        goto fn_exit;
    }

    for (i = 8; i < 64; i += 4) {
        MPI_Type_vector(i, 128, 129, MPI_CHAR, &newtype);

        MPI_Type_commit(&newtype);
        memset(buf, 0, 64 * 129);

        if (rank == 0) {
            /* init buffer */
            for (j = 0; j < i; j++) {
                int k;
                for (k = 0; k < 129; k++) {
                    buf[129 * j + k] = (char) j;
                }
            }

            /* send */
            MPI_Send(buf, 1, newtype, 1, i, MPI_COMM_WORLD);
        } else if (rank == 1) {
            /* recv */
            MPI_Recv(buf, 1, newtype, 0, i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            /* check buffer */
            for (j = 0; j < i; j++) {
                int k;
                for (k = 0; k < 129; k++) {
                    if (k < 128 && buf[129 * j + k] != (char) j) {
                        if (errs < 10) {
                            printf("(i=%d, pos=%d) should be %d but is %d\n",
                                   i, 129 * j + k, j, (int) buf[129 * j + k]);
                        }
                        errs++;
                    } else if (k == 128 && buf[129 * j + k] != (char) 0) {
                        if (errs < 10) {
                            printf("(i=%d, pos=%d) should be %d but is %d\n",
                                   i, 129 * j + k, 0, (int) buf[129 * j + k]);
                        }
                        errs++;
                    }
                }
            }
        }

        MPI_Type_free(&newtype);
    }

    if (rank == 0) {
        int recv_errs = 0;

        MPI_Recv(&recv_errs, 1, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (recv_errs) {
            printf("%d errors reported from receiver\n", recv_errs);
        }
    } else if (rank == 1) {
        MPI_Send(&errs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

  fn_exit:
    free(buf);
    return errs;
}
