/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef POSIX_COLL_H_INCLUDED
#define POSIX_COLL_H_INCLUDED

#include "posix_impl.h"
#include "ch4_impl.h"
#include "posix_coll_release_gather.h"
#include "posix_csel_container.h"


/*
=== BEGIN_MPI_T_CVAR_INFO_BLOCK ===

cvars:
    - name        : MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM
      category    : COLLECTIVE
      type        : enum
      default     : auto
      class       : none
      verbosity   : MPI_T_VERBOSITY_USER_BASIC
      scope       : MPI_T_SCOPE_ALL_EQ
      description : |-
        Variable to select algorithm for intra-node bcast
        mpir           - Fallback to MPIR collectives
        release_gather - Force shm optimized algo using release, gather primitives
        auto - Internal algorithm selection (can be overridden with MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE)

    - name        : MPIR_CVAR_REDUCE_POSIX_INTRA_ALGORITHM
      category    : COLLECTIVE
      type        : enum
      default     : auto
      class       : none
      verbosity   : MPI_T_VERBOSITY_USER_BASIC
      scope       : MPI_T_SCOPE_ALL_EQ
      description : |-
        Variable to select algorithm for intra-node reduce
        mpir           - Fallback to MPIR collectives
        release_gather - Force shm optimized algo using release, gather primitives
        auto - Internal algorithm selection (can be overridden with MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE)

    - name        : MPIR_CVAR_ALLREDUCE_POSIX_INTRA_ALGORITHM
      category    : COLLECTIVE
      type        : enum
      default     : auto
      class       : none
      verbosity   : MPI_T_VERBOSITY_USER_BASIC
      scope       : MPI_T_SCOPE_ALL_EQ
      description : |-
        Variable to select algorithm for intra-node allreduce
        mpir           - Fallback to MPIR collectives
        release_gather - Force shm optimized algo using release, gather primitives
        auto - Internal algorithm selection (can be overridden with MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE)

    - name        : MPIR_CVAR_BARRIER_POSIX_INTRA_ALGORITHM
      category    : COLLECTIVE
      type        : enum
      group       : MPIR_CVAR_GROUP_COLL_ALGO
      default     : auto
      class       : none
      verbosity   : MPI_T_VERBOSITY_USER_BASIC
      scope       : MPI_T_SCOPE_ALL_EQ
      description : |-
        Variable to select algorithm for intra-node barrier
        mpir           - Fallback to MPIR collectives
        release_gather - Force shm optimized algo using release, gather primitives
        auto - Internal algorithm selection (can be overridden with MPIR_CVAR_CH4_POSIX_COLL_SELECTION_TUNING_JSON_FILE)

=== END_MPI_T_CVAR_INFO_BLOCK ===
*/


static inline int MPIDI_POSIX_mpi_barrier(MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Csel_coll_sig_s coll_sig = {
        .coll_type = MPIR_CSEL_COLL_TYPE__BARRIER,
        .comm_ptr = comm,
    };
    MPIDI_POSIX_csel_container_s *cnt;


    switch (MPIR_CVAR_BARRIER_POSIX_INTRA_ALGORITHM) {
        case MPIR_CVAR_BARRIER_POSIX_INTRA_ALGORITHM_release_gather:
            MPII_COLLECTIVE_FALLBACK_CHECK(comm->rank, !MPIR_IS_THREADED, mpi_errno,
                                           "Barrier release_gather cannot be applied.\n");
            mpi_errno = MPIDI_POSIX_mpi_barrier_release_gather(comm, errflag);
            break;

        case MPIR_CVAR_BARRIER_POSIX_INTRA_ALGORITHM_mpir:
            goto fallback;

        case MPIR_CVAR_BARRIER_POSIX_INTRA_ALGORITHM_auto:
            cnt = MPIR_Csel_search(MPIDI_POSIX_COMM(comm, csel_comm), coll_sig);
            if (cnt == NULL)
                goto fallback;

            switch (cnt->id) {
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIDI_POSIX_mpi_barrier_release_gather:
                    mpi_errno =
                        MPIDI_POSIX_mpi_barrier_release_gather(comm, errflag);
                    break;
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIR_Barrier_impl:
                    goto fallback;
                default:
                    MPIR_Assert(0);
            }
            break;

        default:
            MPIR_Assert(0);
    }

    MPIR_ERR_CHECK(mpi_errno);
    goto fn_exit;

  fallback:
    mpi_errno = MPIR_Barrier_impl(comm, errflag);
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_bcast(void *buffer, int count, MPI_Datatype datatype,
                                        int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Csel_coll_sig_s coll_sig = {
        .coll_type = MPIR_CSEL_COLL_TYPE__BCAST,
        .comm_ptr = comm,
        .u.bcast.buffer = buffer,
        .u.bcast.count = count,
        .u.bcast.datatype = datatype,
        .u.bcast.root = root,
    };
    MPIDI_POSIX_csel_container_s *cnt;


    switch (MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM) {
        case MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM_release_gather:
            MPII_COLLECTIVE_FALLBACK_CHECK(comm->rank, !MPIR_IS_THREADED, mpi_errno,
                                           "Bcast release_gather cannot be applied.\n");
            mpi_errno =
                MPIDI_POSIX_mpi_bcast_release_gather(buffer, count, datatype, root, comm, errflag);
            break;

        case MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM_mpir:
            goto fallback;

        case MPIR_CVAR_BCAST_POSIX_INTRA_ALGORITHM_auto:
            cnt = MPIR_Csel_search(MPIDI_POSIX_COMM(comm, csel_comm), coll_sig);
            if (cnt == NULL)
                goto fallback;

            switch (cnt->id) {
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIDI_POSIX_mpi_bcast_release_gather:
                    mpi_errno =
                        MPIDI_POSIX_mpi_bcast_release_gather(buffer, count, datatype, root, comm,
                                                             errflag);
                    break;
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIR_Bcast_impl:
                    goto fallback;
                default:
                    MPIR_Assert(0);
            }
            break;

        default:
            MPIR_Assert(0);
    }

    MPIR_ERR_CHECK(mpi_errno);
    goto fn_exit;

  fallback:
    mpi_errno = MPIR_Bcast_impl(buffer, count, datatype, root, comm, errflag);
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_allreduce(const void *sendbuf, void *recvbuf, int count,
                                            MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                            MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Csel_coll_sig_s coll_sig = {
        .coll_type = MPIR_CSEL_COLL_TYPE__ALLREDUCE,
        .comm_ptr = comm,
        .u.allreduce.sendbuf = sendbuf,
        .u.allreduce.recvbuf = recvbuf,
        .u.allreduce.count = count,
        .u.allreduce.datatype = datatype,
        .u.allreduce.op = op,
    };
    MPIDI_POSIX_csel_container_s *cnt;


    switch (MPIR_CVAR_ALLREDUCE_POSIX_INTRA_ALGORITHM) {
        case MPIR_CVAR_ALLREDUCE_POSIX_INTRA_ALGORITHM_release_gather:
            MPII_COLLECTIVE_FALLBACK_CHECK(comm->rank, !MPIR_IS_THREADED &&
                                           MPIR_Op_is_commutative(op), mpi_errno,
                                           "Allreduce release_gather cannot be applied.\n");
            mpi_errno =
                MPIDI_POSIX_mpi_allreduce_release_gather(sendbuf, recvbuf, count, datatype, op,
                                                         comm, errflag);
            break;

        case MPIR_CVAR_ALLREDUCE_POSIX_INTRA_ALGORITHM_mpir:
            goto fallback;

        case MPIR_CVAR_ALLREDUCE_POSIX_INTRA_ALGORITHM_auto:
            cnt = MPIR_Csel_search(MPIDI_POSIX_COMM(comm, csel_comm), coll_sig);
            if (cnt == NULL)
                goto fallback;

            switch (cnt->id) {
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIDI_POSIX_mpi_allreduce_release_gather:
                    mpi_errno =
                        MPIDI_POSIX_mpi_allreduce_release_gather(sendbuf, recvbuf, count, datatype,
                                                                 op, comm, errflag);
                    break;

                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIR_Allreduce_impl:
                    goto fallback;

                default:
                    MPIR_Assert(0);
            }
            break;

        default:
            MPIR_Assert(0);
    }

    MPIR_ERR_CHECK(mpi_errno);
    goto fn_exit;

  fallback:
    mpi_errno = MPIR_Allreduce_impl(sendbuf, recvbuf, count, datatype, op, comm, errflag);
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_allgather(const void *sendbuf, int sendcount,
                                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                            MPI_Datatype recvtype, MPIR_Comm * comm,
                                            MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Allgather_impl(sendbuf, sendcount, sendtype,
                                    recvbuf, recvcount, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_allgatherv(const void *sendbuf, int sendcount,
                                             MPI_Datatype sendtype, void *recvbuf,
                                             const int *recvcounts, const int *displs,
                                             MPI_Datatype recvtype, MPIR_Comm * comm,
                                             MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Allgatherv_impl(sendbuf, sendcount, sendtype,
                                     recvbuf, recvcounts, displs, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                         void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                         int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Gather_impl(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                                 recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, const int *recvcounts, const int *displs,
                                          MPI_Datatype recvtype, int root, MPIR_Comm * comm,
                                          MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Gatherv_impl(sendbuf, sendcount, sendtype, recvbuf, recvcounts,
                                  displs, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;

}

static inline int MPIDI_POSIX_mpi_scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                          int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Scatter_impl(sendbuf, sendcount, sendtype, recvbuf, recvcount,
                                  recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_scatterv(const void *sendbuf, const int *sendcounts,
                                           const int *displs, MPI_Datatype sendtype,
                                           void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                           int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Scatterv_impl(sendbuf, sendcounts, displs, sendtype,
                                   recvbuf, recvcount, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_alltoall(const void *sendbuf, int sendcount,
                                           MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                           MPI_Datatype recvtype, MPIR_Comm * comm,
                                           MPIR_Errflag_t * errflag)
{
    int mpi_errno;


    mpi_errno = MPIR_Alltoall_impl(sendbuf, sendcount, sendtype,
                                   recvbuf, recvcount, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_alltoallv(const void *sendbuf, const int *sendcounts,
                                            const int *sdispls, MPI_Datatype sendtype,
                                            void *recvbuf, const int *recvcounts,
                                            const int *rdispls, MPI_Datatype recvtype,
                                            MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno;


    mpi_errno = MPIR_Alltoallv_impl(sendbuf, sendcounts, sdispls,
                                    sendtype, recvbuf, recvcounts,
                                    rdispls, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_alltoallw(const void *sendbuf, const int sendcounts[],
                                            const int sdispls[], const MPI_Datatype sendtypes[],
                                            void *recvbuf, const int recvcounts[],
                                            const int rdispls[], const MPI_Datatype recvtypes[],
                                            MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno;


    mpi_errno = MPIR_Alltoallw_impl(sendbuf, sendcounts, sdispls,
                                    sendtypes, recvbuf, recvcounts,
                                    rdispls, recvtypes, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_reduce(const void *sendbuf, void *recvbuf, int count,
                                         MPI_Datatype datatype, MPI_Op op, int root,
                                         MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_Csel_coll_sig_s coll_sig = {
        .coll_type = MPIR_CSEL_COLL_TYPE__REDUCE,
        .comm_ptr = comm,
        .u.reduce.sendbuf = sendbuf,
        .u.reduce.recvbuf = recvbuf,
        .u.reduce.count = count,
        .u.reduce.datatype = datatype,
        .u.reduce.op = op,
        .u.reduce.root = root,
    };
    MPIDI_POSIX_csel_container_s *cnt;


    switch (MPIR_CVAR_REDUCE_POSIX_INTRA_ALGORITHM) {
        case MPIR_CVAR_REDUCE_POSIX_INTRA_ALGORITHM_release_gather:
            MPII_COLLECTIVE_FALLBACK_CHECK(comm->rank, !MPIR_IS_THREADED &&
                                           MPIR_Op_is_commutative(op), mpi_errno,
                                           "Reduce release_gather cannot be applied.\n");
            mpi_errno =
                MPIDI_POSIX_mpi_reduce_release_gather(sendbuf, recvbuf, count, datatype, op, root,
                                                      comm, errflag);
            break;

        case MPIR_CVAR_REDUCE_POSIX_INTRA_ALGORITHM_mpir:
            goto fallback;

        case MPIR_CVAR_REDUCE_POSIX_INTRA_ALGORITHM_auto:
            cnt = MPIR_Csel_search(MPIDI_POSIX_COMM(comm, csel_comm), coll_sig);
            if (cnt == NULL)
                goto fallback;

            switch (cnt->id) {
                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIDI_POSIX_mpi_reduce_release_gather:
                    mpi_errno =
                        MPIDI_POSIX_mpi_reduce_release_gather(sendbuf, recvbuf, count, datatype, op,
                                                              root, comm, errflag);
                    break;

                case MPIDI_POSIX_CSEL_CONTAINER_TYPE__ALGORITHM__MPIR_Reduce_impl:
                    goto fallback;

                default:
                    MPIR_Assert(0);
            }
            break;

        default:
            MPIR_Assert(0);
    }

    MPIR_ERR_CHECK(mpi_errno);
    goto fn_exit;

  fallback:
    mpi_errno = MPIR_Reduce_impl(sendbuf, recvbuf, count, datatype, op, root, comm, errflag);
    MPIR_ERR_CHECK(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_reduce_scatter(const void *sendbuf, void *recvbuf,
                                                 const int recvcounts[], MPI_Datatype datatype,
                                                 MPI_Op op, MPIR_Comm * comm,
                                                 MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Reduce_scatter_impl(sendbuf, recvbuf, recvcounts, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_reduce_scatter_block(const void *sendbuf, void *recvbuf,
                                                       int recvcount, MPI_Datatype datatype,
                                                       MPI_Op op, MPIR_Comm * comm,
                                                       MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    MPIR_Reduce_scatter_block_impl(sendbuf, recvbuf, recvcount, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_scan(const void *sendbuf, void *recvbuf, int count,
                                       MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                       MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Scan_impl(sendbuf, recvbuf, count, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_exscan(const void *sendbuf, void *recvbuf, int count,
                                         MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                         MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Exscan_impl(sendbuf, recvbuf, count, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_POSIX_mpi_neighbor_allgather(const void *sendbuf, int sendcount,
                                                     MPI_Datatype sendtype, void *recvbuf,
                                                     int recvcount, MPI_Datatype recvtype,
                                                     MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno =
        MPIR_Neighbor_allgather(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_neighbor_allgatherv(const void *sendbuf, int sendcount,
                                                      MPI_Datatype sendtype, void *recvbuf,
                                                      const int recvcounts[], const int displs[],
                                                      MPI_Datatype recvtype, MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_allgatherv(sendbuf, sendcount, sendtype,
                                         recvbuf, recvcounts, displs, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_neighbor_alltoall(const void *sendbuf, int sendcount,
                                                    MPI_Datatype sendtype, void *recvbuf,
                                                    int recvcount, MPI_Datatype recvtype,
                                                    MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoall(sendbuf, sendcount, sendtype,
                                       recvbuf, recvcount, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_neighbor_alltoallv(const void *sendbuf, const int sendcounts[],
                                                     const int sdispls[], MPI_Datatype sendtype,
                                                     void *recvbuf, const int recvcounts[],
                                                     const int rdispls[], MPI_Datatype recvtype,
                                                     MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
                                        recvbuf, recvcounts, rdispls, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_neighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                                                     const MPI_Aint sdispls[],
                                                     const MPI_Datatype sendtypes[], void *recvbuf,
                                                     const int recvcounts[],
                                                     const MPI_Aint rdispls[],
                                                     const MPI_Datatype recvtypes[],
                                                     MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
                                        recvbuf, recvcounts, rdispls, recvtypes, comm);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ineighbor_allgather(const void *sendbuf, int sendcount,
                                                      MPI_Datatype sendtype, void *recvbuf,
                                                      int recvcount, MPI_Datatype recvtype,
                                                      MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_allgather(sendbuf, sendcount, sendtype,
                                         recvbuf, recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ineighbor_allgatherv(const void *sendbuf, int sendcount,
                                                       MPI_Datatype sendtype, void *recvbuf,
                                                       const int recvcounts[], const int displs[],
                                                       MPI_Datatype recvtype, MPIR_Comm * comm,
                                                       MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_allgatherv(sendbuf, sendcount, sendtype,
                                          recvbuf, recvcounts, displs, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ineighbor_alltoall(const void *sendbuf, int sendcount,
                                                     MPI_Datatype sendtype, void *recvbuf,
                                                     int recvcount, MPI_Datatype recvtype,
                                                     MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoall(sendbuf, sendcount, sendtype,
                                        recvbuf, recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ineighbor_alltoallv(const void *sendbuf, const int sendcounts[],
                                                      const int sdispls[], MPI_Datatype sendtype,
                                                      void *recvbuf, const int recvcounts[],
                                                      const int rdispls[], MPI_Datatype recvtype,
                                                      MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoallv(sendbuf, sendcounts, sdispls, sendtype,
                                         recvbuf, recvcounts, rdispls, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                                                      const MPI_Aint sdispls[],
                                                      const MPI_Datatype sendtypes[], void *recvbuf,
                                                      const int recvcounts[],
                                                      const MPI_Aint rdispls[],
                                                      const MPI_Datatype recvtypes[],
                                                      MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoallw(sendbuf, sendcounts, sdispls, sendtypes,
                                         recvbuf, recvcounts, rdispls, recvtypes, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ibarrier(MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ibarrier(comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ibcast(void *buffer, int count, MPI_Datatype datatype,
                                         int root, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ibcast(buffer, count, datatype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iallgather(const void *sendbuf, int sendcount,
                                             MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                             MPI_Datatype recvtype, MPIR_Comm * comm,
                                             MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallgather(sendbuf, sendcount, sendtype, recvbuf,
                                recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iallgatherv(const void *sendbuf, int sendcount,
                                              MPI_Datatype sendtype, void *recvbuf,
                                              const int *recvcounts, const int *displs,
                                              MPI_Datatype recvtype, MPIR_Comm * comm,
                                              MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallgatherv(sendbuf, sendcount, sendtype,
                                 recvbuf, recvcounts, displs, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ialltoall(const void *sendbuf, int sendcount,
                                            MPI_Datatype sendtype, void *recvbuf, int recvcount,
                                            MPI_Datatype recvtype, MPIR_Comm * comm,
                                            MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoall(sendbuf, sendcount, sendtype, recvbuf,
                               recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ialltoallv(const void *sendbuf, const int *sendcounts,
                                             const int *sdispls, MPI_Datatype sendtype,
                                             void *recvbuf, const int *recvcounts,
                                             const int *rdispls, MPI_Datatype recvtype,
                                             MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoallv(sendbuf, sendcounts, sdispls,
                                sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ialltoallw(const void *sendbuf, const int *sendcounts,
                                             const int *sdispls, const MPI_Datatype sendtypes[],
                                             void *recvbuf, const int *recvcounts,
                                             const int *rdispls, const MPI_Datatype recvtypes[],
                                             MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoallw(sendbuf, sendcounts, sdispls,
                                sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iexscan(const void *sendbuf, void *recvbuf, int count,
                                          MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                          MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iexscan(sendbuf, recvbuf, count, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                          int root, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Igather(sendbuf, sendcount, sendtype, recvbuf,
                             recvcount, recvtype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_igatherv(const void *sendbuf, int sendcount,
                                           MPI_Datatype sendtype, void *recvbuf,
                                           const int *recvcounts, const int *displs,
                                           MPI_Datatype recvtype, int root, MPIR_Comm * comm,
                                           MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Igatherv(sendbuf, sendcount, sendtype,
                              recvbuf, recvcounts, displs, recvtype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ireduce_scatter_block(const void *sendbuf, void *recvbuf,
                                                        int recvcount, MPI_Datatype datatype,
                                                        MPI_Op op, MPIR_Comm * comm,
                                                        MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce_scatter_block(sendbuf, recvbuf, recvcount, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ireduce_scatter(const void *sendbuf, void *recvbuf,
                                                  const int recvcounts[], MPI_Datatype datatype,
                                                  MPI_Op op, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce_scatter(sendbuf, recvbuf, recvcounts, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_ireduce(const void *sendbuf, void *recvbuf, int count,
                                          MPI_Datatype datatype, MPI_Op op, int root,
                                          MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce(sendbuf, recvbuf, count, datatype, op, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iallreduce(const void *sendbuf, void *recvbuf, int count,
                                             MPI_Datatype datatype, MPI_Op op,
                                             MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallreduce(sendbuf, recvbuf, count, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iscan(const void *sendbuf, void *recvbuf, int count,
                                        MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                        MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscan(sendbuf, recvbuf, count, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iscatter(const void *sendbuf, int sendcount,
                                           MPI_Datatype sendtype, void *recvbuf,
                                           int recvcount, MPI_Datatype recvtype,
                                           int root, MPIR_Comm * comm, MPIR_Request ** request)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscatter(sendbuf, sendcount, sendtype, recvbuf,
                              recvcount, recvtype, root, comm, request);

    return mpi_errno;
}

static inline int MPIDI_POSIX_mpi_iscatterv(const void *sendbuf, const int *sendcounts,
                                            const int *displs, MPI_Datatype sendtype,
                                            void *recvbuf, int recvcount,
                                            MPI_Datatype recvtype, int root,
                                            MPIR_Comm * comm, MPIR_Request ** request)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscatterv(sendbuf, sendcounts, displs, sendtype,
                               recvbuf, recvcount, recvtype, root, comm, request);

    return mpi_errno;
}

#endif /* POSIX_COLL_H_INCLUDED */
