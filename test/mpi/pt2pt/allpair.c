/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "mpitest.h"

/*
 * This test is based on the allpair.f test from the MPICH-1 test
 * (test/pt2pt/allpair.f), which in turn was inspired by a bug report from
 * fsset@corelli.lerc.nasa.gov (Scott Townsend)
 * This is the C equivallent version of f77/pt2pt/sendf.f and friends.
 */

#define TEST_SIZE 2000

int errs;
int rank, size;
int next, prev;
int tag, count;
float send_buf[TEST_SIZE];
float recv_buf[TEST_SIZE];

void init_global(MPI_Comm comm);
void test_pair_send(MPI_Comm comm);
void test_pair_psend(MPI_Comm comm);
void test_pair_prsend(MPI_Comm comm);
void init_test_data(float *buf, int count);
void clear_test_data(float *buf, int count);
void msg_check(float *buf, int source, int tag, int count, MPI_Status status, int max);

int main(int argc, char *argv[])
{
    MTest_Init(&argc, &argv);

    MPI_Comm comm;
#if 1
    while (MTestGetIntracommGeneral(&comm, 2, 0)) {
        if (comm == MPI_COMM_NULL) {
            continue;
        }
#else
    comm = MPI_COMM_WORLD;
    if (1) {
#endif
        init_global(comm);
        test_pair_prsend(comm);
        MTestFreeComm(&comm);
    }
    MTest_Finalize(errs);
    return errs;
}

void init_global(MPI_Comm comm)
{
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    printf("Testing comm %d / %d ...\n", rank, size);

    next = (rank == size - 1) ? 0 : rank + 1;
    prev = (rank == 0) ? size - 1 : rank - 1;

    tag = 3123;
    count = TEST_SIZE / 5;
}

void test_pair_send(MPI_Comm comm)
{
    MPI_Status statuses[2];

    clear_test_data(recv_buf, TEST_SIZE);

    if (rank == 0) {
        init_test_data(send_buf, TEST_SIZE);
        MPI_Send(send_buf, count, MPI_FLOAT, next, tag, comm);
        MPI_Recv(recv_buf, count, MPI_FLOAT, next, tag, comm, &statuses[1]);
        msg_check(recv_buf, next, tag, count, statuses[1], TEST_SIZE);
    } else if (prev == 0) {
        MPI_Recv(recv_buf, count, MPI_FLOAT, prev, tag, comm, &statuses[1]);
        msg_check(recv_buf, next, tag, count, statuses[1], TEST_SIZE);
        MPI_Send(recv_buf, count, MPI_FLOAT, next, tag, comm);
    }
}

void test_pair_psend(MPI_Comm comm)
{
    MPI_Status statuses[2];

    clear_test_data(recv_buf, TEST_SIZE);

    MPI_Request preqs[2];
    MPI_Recv_init(recv_buf, count, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &preqs[1]);
    MPI_Send_init(send_buf, count, MPI_FLOAT, next, tag, comm, &preqs[0]);

    if (rank == 0) {
        init_test_data(send_buf, TEST_SIZE);
        MPI_Startall(2, preqs);
        MPI_Waitall(2, preqs, statuses);
        msg_check(recv_buf, next, tag, count, statuses[1], TEST_SIZE);
    } else if (prev == 0) {
        MPI_Start(&preqs[1]);
        MPI_Wait(&preqs[1], &statuses[1]);
        msg_check(recv_buf, next, tag, count, statuses[1], TEST_SIZE);
        for (int i = 0; i < count; i++) {
            send_buf[i] = recv_buf[i];
        }
        MPI_Start(&preqs[0]);
        MPI_Wait(&preqs[0], &statuses[0]);
    }

    MPI_Request_free(&preqs[0]);
    MPI_Request_free(&preqs[1]);
}

void test_pair_prsend(MPI_Comm comm)
{
    MPI_Status statuses[2];

    clear_test_data(recv_buf, TEST_SIZE);

    MPI_Request preqs[2];
    MPI_Recv_init(recv_buf, count, MPI_FLOAT, MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &preqs[1]);
    MPI_Rsend_init(send_buf, count, MPI_FLOAT, next, tag, comm, &preqs[0]);

    if (rank == 0) {
        init_test_data(send_buf, TEST_SIZE);
        MPI_Recv(MPI_BOTTOM, 0, MPI_INT, next, tag, comm, MPI_STATUS_IGNORE);
        MPI_Startall(2, preqs);
        int index = -1;
        while (index != 2) {
            int outcount, indices[2];
            MPI_Waitsome(2, preqs, &outcount, indices, statuses);
            for (int i = 0; i < outcount; i++) {
                if (indices[i] == 1) {
                    msg_check(recv_buf, next, tag, count, statuses[i], TEST_SIZE);
                    index = 2;
                }
            }
        }
    } else if (prev == 0) {
        MPI_Send(MPI_BOTTOM, 0, MPI_INT, prev, tag, comm);
        MPI_Start(&preqs[1]);
        int flag = 0;
        while (!flag) {
            MPI_Test(&preqs[1], &flag, &statuses[1]);
        }
        msg_check(recv_buf, next, tag, count, statuses[1], TEST_SIZE);
        for (int i = 0; i < count; i++) {
            send_buf[i] = recv_buf[i];
        }
        MPI_Start(&preqs[0]);
        MPI_Wait(&preqs[0], &statuses[0]);
    }

    MPI_Request_free(&preqs[0]);
    MPI_Request_free(&preqs[1]);
}

/* ---------------------------------------- */
void init_test_data(float *buf, int count)
{
    for (int i = 0; i < count; i++) {
        buf[i] = (float) i;
    }
}

void clear_test_data(float *buf, int count)
{
    for (int i = 0; i < count; i++) {
        buf[i] = 0;
    }
}

void verify_test_data(float *buf, int count, int max)
{
    for (int i = 0; i < count; i++) {
        if (buf[i] != (float) i) {
            if (errs < 10) {
                printf("buffer[%d] got %f, expect %f\n", i, buf[i], (float) i);
            }
            errs++;
        }
    }
    for (int i = count; i < max; i++) {
        if (buf[i] != 0) {
            if (errs < 10) {
                printf("buffer[%d] got %f, expect 0\n", i, buf[i]);
            }
            errs++;
        }
    }
}

void msg_check(float *buf, int source, int tag, int count, MPI_Status status, int max)
{
    if (status.MPI_SOURCE != source) {
        printf("Unexpected source: %d\n", status.MPI_SOURCE);
        errs++;
    }

    if (status.MPI_TAG != tag) {
        printf("Unexpected tag: %d\n", status.MPI_TAG);
        errs++;
    }

    int recv_count;
    MPI_Get_count(&status, MPI_FLOAT, &recv_count);
    if (recv_count != count) {
        printf("Unexpected count: %d\n", recv_count);
        errs++;
    }

    verify_test_data(buf, count, max);
}
