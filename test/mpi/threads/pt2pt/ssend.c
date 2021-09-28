/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include <stdio.h>
#include <mpi.h>
#include "mpitest.h"
#include "mpithreadtest.h"

#define MAX_COUNT 1024 * 1600

#define NUM_THREADS        4
#define NUM_MSG_PER_THREAD 10
#define NUM_CHECK          4

#define TOTAL  NUM_THREADS * NUM_MSG_PER_THREAD

int buf[TOTAL][MAX_COUNT];
MPI_Request reqs[TOTAL];
int counts[NUM_THREADS];

MPI_Comm comm = MPI_COMM_WORLD;
int tag = 1;

static MTEST_THREAD_RETURN_TYPE do_ssend(void *arg)
{
    int id = (long) arg;
    int base = id * NUM_MSG_PER_THREAD;

    for (int i = 0; i < NUM_MSG_PER_THREAD; i++) {
        buf[base + i][0] = id;
        int count = 1;
        if (i % 2 == 0) {
            count = MAX_COUNT;
        }
        MPI_Issend(buf[base + i], count, MPI_INT, 1, tag, comm, &reqs[base + i]);
    }
    return NULL;
}

static MTEST_THREAD_RETURN_TYPE do_recv(void *arg)
{
    int id = (long) arg;
    int base = id * NUM_MSG_PER_THREAD;

    for (int i = 0; i < NUM_CHECK; i++) {
        MPI_Irecv(buf[base + i], MAX_COUNT, MPI_INT, 0, tag, comm, &reqs[base + i]);
    }
    MPI_Waitall(NUM_CHECK, reqs + base, MPI_STATUSES_IGNORE);
    return NULL;
}

int main(int argc, char *argv[])
{
    int errs = 0;

    int provided;
    MTest_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (provided != MPI_THREAD_MULTIPLE) {
        printf("MPI_THREAD_MULTIPLE not supported by the MPI implementation\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    if (size != 2) {
        printf("This test require 2 processes\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    if (rank == 0) {
        for (int i = 0; i < NUM_THREADS; i++) {
            /* Issend NUM_MSG_PER_THREAD * NUM_THREADS messages */
            MTest_Start_thread(do_ssend, (void *) (long) i);
        }
        MTest_Join_threads();

        for (int i = 0; i < NUM_CHECK * NUM_THREADS; i++) {
            int id, indx;
            MPI_Waitany(TOTAL, reqs, &indx, MPI_STATUS_IGNORE);
            id = indx / NUM_MSG_PER_THREAD;
            printf(" - %d - %d send complete\n", id, indx);
            counts[id]++;
        }

        MPI_Send(counts, NUM_THREADS, MPI_INT, 1, tag + 1, comm);

        MPI_Barrier(comm);

        MPI_Waitall(TOTAL, reqs, MPI_STATUSES_IGNORE);
    } else {
#if 0
        for (int i = 0; i < NUM_THREADS; i++) {
            /* Receive NUM_CHECK * NUM_THREADS messages */
            MTest_Start_thread(do_recv, (void *) (long) i);
        }
        MTest_Join_threads();

        for (int j = 0; j < NUM_THREADS; j++) {
            for (int i = 0; i < NUM_CHECK; i++) {
                int id = buf[j * NUM_MSG_PER_THREAD + i][0];
                counts[id]++;
            }
        }
#else
        for (int i = 0; i < NUM_CHECK * NUM_THREADS; i++) {
            MPI_Irecv(buf[i], MAX_COUNT, MPI_INT, 0, tag, comm, &reqs[i]);
        }
        MPI_Waitall(NUM_CHECK * NUM_THREADS, reqs, MPI_STATUSES_IGNORE);
        for (int i = 0; i < NUM_CHECK * NUM_THREADS; i++) {
            int id = buf[i][0];
            counts[id]++;
        }
#endif
        int recv_counts[TOTAL];
        MPI_Recv(recv_counts, NUM_THREADS, MPI_INT, 0, tag + 1, comm, MPI_STATUS_IGNORE);
        for (int i = 0; i < NUM_THREADS; i++) {
            if (counts[i] != recv_counts[i]) {
                errs++;
            }
            printf("From thread %d, received %d messages, sender reported %d ssend completed\n", i,
                   counts[i], recv_counts[i]);
        }

        MPI_Barrier(comm);

        for (int i = NUM_CHECK * NUM_THREADS; i < TOTAL; i++) {
            MPI_Recv(&buf[i], MAX_COUNT, MPI_INT, 0, tag, comm, MPI_STATUS_IGNORE);
        }
    }

    MTest_Finalize(errs);
    return 0;
}
