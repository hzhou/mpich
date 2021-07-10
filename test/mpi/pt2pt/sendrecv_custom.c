/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpitest.h"
#include "dtpools.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main(int argc, char *argv[])
{
    int errs = 0;
    int ret;

    MTest_Init(&argc, &argv);

    MPI_Comm comm = MPI_COMM_WORLD;
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int src = 0;
    int dst = 1;
    int tag = 0;

    MPI_Aint maxbufsize = MTestDefaultMaxBufferSize();

    int count = 100;
    int seed = 1;
    DTP_pool_s dtp;
    ret = DTP_pool_create("MPI_INT", count, seed, &dtp);

    for (int iter = 0; iter < 3; iter++) {
        DTP_obj_s send_obj, recv_obj;
        if (iter == 0) {
            DTP_obj_create(dtp, &send_obj, maxbufsize);
            DTP_obj_create(dtp, &recv_obj, maxbufsize);
        } else if (iter == 1) {
            if (rank == src) {
                DTP_obj_create_idx(dtp, &send_obj, maxbufsize, 0);
            }
            if (rank == dst) {
                DTP_obj_create_idx(dtp, &recv_obj, maxbufsize, 1);
            }
        } else if (iter == 2) {
            if (rank == src) {
                DTP_obj_create_custom(dtp, &send_obj,
                                      "0: structure [numblks 2, blklen (5,5), displs (8,208)]"
                                      "1: blkhindx [numblks 2, blklen 1, displs (8,40)]"
                                      "2: resized [lb 8, extent 8]");
            }
            if (rank == dst) {
                DTP_obj_create_custom(dtp, &recv_obj,
                                      "0: blkhindx [numblks 2, blklen 1, displs (8,40)]"
                                      "1: resized [lb 8, extent 8]");
            }
        } else {
            assert(0);
        }

        if (rank == src) {
            MTestPrintfMsg(0, "iter=%d, Send obj: count=%ld, %s\n", iter, send_obj.DTP_type_count,
                           DTP_obj_get_description(send_obj));
            void *send_buf = malloc(send_obj.DTP_bufsize);
            DTP_obj_buf_init(send_obj, send_buf, 1, 2, count);
            ret = MPI_Send((char *) send_buf + send_obj.DTP_buf_offset,
                           send_obj.DTP_type_count, send_obj.DTP_datatype, dst, tag, comm);
            free(send_buf);
        } else if (rank == dst) {
            MTestPrintfMsg(0, "iter=%d, Recv obj: count=%ld, %s\n", iter, recv_obj.DTP_type_count,
                           DTP_obj_get_description(recv_obj));
            void *recv_buf = malloc(recv_obj.DTP_bufsize);
            ret = MPI_Recv((char *) recv_buf + recv_obj.DTP_buf_offset,
                           recv_obj.DTP_type_count,
                           recv_obj.DTP_datatype, src, tag, comm, MPI_STATUS_IGNORE);
            ret = DTP_obj_buf_check(recv_obj, recv_buf, 1, 2, count);
            if (ret != DTP_SUCCESS) {
                errs++;
                printf("Send datatype: count = %ld, %s\n", send_obj.DTP_type_count,
                       DTP_obj_get_description(send_obj));
                printf("Recv datatype: count = %ld, %s\n", recv_obj.DTP_type_count,
                       DTP_obj_get_description(recv_obj));
            }
            free(recv_buf);
        }
        if (iter == 0 || rank == src) {
            DTP_obj_free(send_obj);
        }
        if (iter == 0 || rank == dst) {
            DTP_obj_free(recv_obj);
        }
    }

    DTP_pool_free(dtp);
    MTest_Finalize(errs);
    return MTestReturnValue(errs);
}
