/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef OFI_COLL_H_INCLUDED
#define OFI_COLL_H_INCLUDED

#include "ofi_impl.h"
#include "ch4_csel_container.h"

static inline int MPIDI_NM_mpi_barrier(MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Barrier_impl(comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_bcast(void *buffer, int count, MPI_Datatype datatype,
                                     int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Bcast_impl(buffer, count, datatype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_allreduce(const void *sendbuf, void *recvbuf, int count,
                                         MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                         MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Allreduce_impl(sendbuf, recvbuf, count, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_allgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                         void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                         MPIR_Comm * comm, MPIR_Errflag_t * errflag)
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

static inline int MPIDI_NM_mpi_allgatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, const int *recvcounts, const int *displs,
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

static inline int MPIDI_NM_mpi_gather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                      void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                      int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Gather_impl(sendbuf, sendcount, sendtype, recvbuf,
                                 recvcount, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_gatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                       void *recvbuf, const int *recvcounts, const int *displs,
                                       MPI_Datatype recvtype, int root, MPIR_Comm * comm,
                                       MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Gatherv_impl(sendbuf, sendcount, sendtype, recvbuf,
                                  recvcounts, displs, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;

}

static inline int MPIDI_NM_mpi_scatter(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                       void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                       int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Scatter_impl(sendbuf, sendcount, sendtype, recvbuf,
                                  recvcount, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_scatterv(const void *sendbuf, const int *sendcounts,
                                        const int *displs, MPI_Datatype sendtype,
                                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                        int root, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Scatterv_impl(sendbuf, sendcounts, displs, sendtype, recvbuf,
                                   recvcount, recvtype, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_alltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                        void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                        MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Alltoall_impl(sendbuf, sendcount, sendtype,
                                   recvbuf, recvcount, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_alltoallv(const void *sendbuf, const int *sendcounts,
                                         const int *sdispls, MPI_Datatype sendtype,
                                         void *recvbuf, const int *recvcounts,
                                         const int *rdispls, MPI_Datatype recvtype,
                                         MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Alltoallv_impl(sendbuf, sendcounts, sdispls,
                                    sendtype, recvbuf, recvcounts,
                                    rdispls, recvtype, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_alltoallw(const void *sendbuf, const int sendcounts[],
                                         const int sdispls[], const MPI_Datatype sendtypes[],
                                         void *recvbuf, const int recvcounts[],
                                         const int rdispls[], const MPI_Datatype recvtypes[],
                                         MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Alltoallw_impl(sendbuf, sendcounts, sdispls,
                                    sendtypes, recvbuf, recvcounts,
                                    rdispls, recvtypes, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_reduce(const void *sendbuf, void *recvbuf, int count,
                                      MPI_Datatype datatype, MPI_Op op, int root,
                                      MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Reduce_impl(sendbuf, recvbuf, count, datatype, op, root, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_reduce_scatter(const void *sendbuf, void *recvbuf,
                                              const int recvcounts[], MPI_Datatype datatype,
                                              MPI_Op op, MPIR_Comm * comm, MPIR_Errflag_t * errflag)
{
    int mpi_errno = MPI_SUCCESS;


    mpi_errno = MPIR_Reduce_scatter_impl(sendbuf, recvbuf, recvcounts, datatype, op, comm, errflag);

    MPIR_ERR_CHECK(mpi_errno);


  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

static inline int MPIDI_NM_mpi_reduce_scatter_block(const void *sendbuf, void *recvbuf,
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

static inline int MPIDI_NM_mpi_scan(const void *sendbuf, void *recvbuf, int count,
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

static inline int MPIDI_NM_mpi_exscan(const void *sendbuf, void *recvbuf, int count,
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

static inline int MPIDI_NM_mpi_neighbor_allgather(const void *sendbuf, int sendcount,
                                                  MPI_Datatype sendtype, void *recvbuf,
                                                  int recvcount, MPI_Datatype recvtype,
                                                  MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno =
        MPIR_Neighbor_allgather_impl(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype,
                                     comm);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_neighbor_allgatherv(const void *sendbuf, int sendcount,
                                                   MPI_Datatype sendtype, void *recvbuf,
                                                   const int recvcounts[], const int displs[],
                                                   MPI_Datatype recvtype, MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_allgatherv_impl(sendbuf, sendcount, sendtype,
                                              recvbuf, recvcounts, displs, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_neighbor_alltoall(const void *sendbuf, int sendcount,
                                                 MPI_Datatype sendtype, void *recvbuf,
                                                 int recvcount, MPI_Datatype recvtype,
                                                 MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoall_impl(sendbuf, sendcount, sendtype,
                                            recvbuf, recvcount, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_neighbor_alltoallv(const void *sendbuf, const int sendcounts[],
                                                  const int sdispls[], MPI_Datatype sendtype,
                                                  void *recvbuf, const int recvcounts[],
                                                  const int rdispls[], MPI_Datatype recvtype,
                                                  MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoallv_impl(sendbuf, sendcounts, sdispls, sendtype,
                                             recvbuf, recvcounts, rdispls, recvtype, comm);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_neighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                                                  const MPI_Aint sdispls[],
                                                  const MPI_Datatype sendtypes[], void *recvbuf,
                                                  const int recvcounts[], const MPI_Aint rdispls[],
                                                  const MPI_Datatype recvtypes[], MPIR_Comm * comm)
{
    int mpi_errno;

    mpi_errno = MPIR_Neighbor_alltoallw_impl(sendbuf, sendcounts, sdispls, sendtypes,
                                             recvbuf, recvcounts, rdispls, recvtypes, comm);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ineighbor_allgather(const void *sendbuf, int sendcount,
                                                   MPI_Datatype sendtype, void *recvbuf,
                                                   int recvcount, MPI_Datatype recvtype,
                                                   MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_allgather_impl(sendbuf, sendcount, sendtype,
                                              recvbuf, recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ineighbor_allgatherv(const void *sendbuf, int sendcount,
                                                    MPI_Datatype sendtype, void *recvbuf,
                                                    const int recvcounts[], const int displs[],
                                                    MPI_Datatype recvtype, MPIR_Comm * comm,
                                                    MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_allgatherv_impl(sendbuf, sendcount, sendtype,
                                               recvbuf, recvcounts, displs, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ineighbor_alltoall(const void *sendbuf, int sendcount,
                                                  MPI_Datatype sendtype, void *recvbuf,
                                                  int recvcount, MPI_Datatype recvtype,
                                                  MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoall_impl(sendbuf, sendcount, sendtype,
                                             recvbuf, recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ineighbor_alltoallv(const void *sendbuf, const int sendcounts[],
                                                   const int sdispls[], MPI_Datatype sendtype,
                                                   void *recvbuf, const int recvcounts[],
                                                   const int rdispls[], MPI_Datatype recvtype,
                                                   MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoallv_impl(sendbuf, sendcounts, sdispls, sendtype,
                                              recvbuf, recvcounts, rdispls, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ineighbor_alltoallw(const void *sendbuf, const int sendcounts[],
                                                   const MPI_Aint sdispls[],
                                                   const MPI_Datatype sendtypes[], void *recvbuf,
                                                   const int recvcounts[], const MPI_Aint rdispls[],
                                                   const MPI_Datatype recvtypes[],
                                                   MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ineighbor_alltoallw_impl(sendbuf, sendcounts, sdispls, sendtypes,
                                              recvbuf, recvcounts, rdispls, recvtypes, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ibarrier(MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ibarrier_impl(comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ibcast(void *buffer, int count, MPI_Datatype datatype,
                                      int root, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ibcast_impl(buffer, count, datatype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iallgather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                          void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                          MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallgather_impl(sendbuf, sendcount, sendtype, recvbuf,
                                     recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iallgatherv(const void *sendbuf, int sendcount,
                                           MPI_Datatype sendtype, void *recvbuf,
                                           const int *recvcounts, const int *displs,
                                           MPI_Datatype recvtype, MPIR_Comm * comm,
                                           MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallgatherv_impl(sendbuf, sendcount, sendtype,
                                      recvbuf, recvcounts, displs, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iallreduce(const void *sendbuf, void *recvbuf, int count,
                                          MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                          MPIR_Request ** request)
{
    int mpi_errno;

    mpi_errno = MPIR_Iallreduce_impl(sendbuf, recvbuf, count, datatype, op, comm, request);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ialltoall(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                         void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                         MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoall_impl(sendbuf, sendcount, sendtype, recvbuf,
                                    recvcount, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ialltoallv(const void *sendbuf, const int *sendcounts,
                                          const int *sdispls, MPI_Datatype sendtype,
                                          void *recvbuf, const int *recvcounts,
                                          const int *rdispls, MPI_Datatype recvtype,
                                          MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoallv_impl(sendbuf, sendcounts, sdispls,
                                     sendtype, recvbuf, recvcounts, rdispls, recvtype, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ialltoallw(const void *sendbuf, const int *sendcounts,
                                          const int *sdispls, const MPI_Datatype sendtypes[],
                                          void *recvbuf, const int *recvcounts,
                                          const int *rdispls, const MPI_Datatype recvtypes[],
                                          MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ialltoallw_impl(sendbuf, sendcounts, sdispls,
                                     sendtypes, recvbuf, recvcounts, rdispls, recvtypes, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iexscan(const void *sendbuf, void *recvbuf, int count,
                                       MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                       MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iexscan_impl(sendbuf, recvbuf, count, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_igather(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                       void *recvbuf, int recvcount, MPI_Datatype recvtype,
                                       int root, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Igather_impl(sendbuf, sendcount, sendtype, recvbuf,
                                  recvcount, recvtype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
                                        void *recvbuf, const int *recvcounts, const int *displs,
                                        MPI_Datatype recvtype, int root, MPIR_Comm * comm,
                                        MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Igatherv_impl(sendbuf, sendcount, sendtype,
                                   recvbuf, recvcounts, displs, recvtype, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ireduce_scatter_block(const void *sendbuf, void *recvbuf,
                                                     int recvcount, MPI_Datatype datatype,
                                                     MPI_Op op, MPIR_Comm * comm,
                                                     MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce_scatter_block_impl(sendbuf, recvbuf, recvcount,
                                                datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ireduce_scatter(const void *sendbuf, void *recvbuf,
                                               const int recvcounts[], MPI_Datatype datatype,
                                               MPI_Op op, MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce_scatter_impl(sendbuf, recvbuf, recvcounts, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_ireduce(const void *sendbuf, void *recvbuf, int count,
                                       MPI_Datatype datatype, MPI_Op op, int root,
                                       MPIR_Comm * comm, MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Ireduce_impl(sendbuf, recvbuf, count, datatype, op, root, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iscan(const void *sendbuf, void *recvbuf, int count,
                                     MPI_Datatype datatype, MPI_Op op, MPIR_Comm * comm,
                                     MPIR_Request ** req)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscan_impl(sendbuf, recvbuf, count, datatype, op, comm, req);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iscatter(const void *sendbuf, int sendcount,
                                        MPI_Datatype sendtype, void *recvbuf,
                                        int recvcount, MPI_Datatype recvtype,
                                        int root, MPIR_Comm * comm, MPIR_Request ** request)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscatter_impl(sendbuf, sendcount, sendtype, recvbuf,
                                   recvcount, recvtype, root, comm, request);

    return mpi_errno;
}

static inline int MPIDI_NM_mpi_iscatterv(const void *sendbuf, const int *sendcounts,
                                         const int *displs, MPI_Datatype sendtype,
                                         void *recvbuf, int recvcount,
                                         MPI_Datatype recvtype, int root,
                                         MPIR_Comm * comm, MPIR_Request ** request)
{
    int mpi_errno;

    mpi_errno = MPIR_Iscatterv_impl(sendbuf, sendcounts, displs, sendtype,
                                    recvbuf, recvcount, recvtype, root, comm, request);

    return mpi_errno;
}

#endif /* OFI_COLL_H_INCLUDED */
