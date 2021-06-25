/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#ifndef UCX_AM_H_INCLUDED
#define UCX_AM_H_INCLUDED

#include "ucx_impl.h"

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_am_isend(int rank,
                                               MPIR_Comm * comm,
                                               int handler_id,
                                               const void *am_hdr,
                                               MPI_Aint am_hdr_sz,
                                               const void *data,
                                               MPI_Aint count, MPI_Datatype datatype,
                                               MPIR_Request * sreq)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_UCX_ucp_request_t *ucp_request;
    ucp_ep_h ep;
    MPIDI_UCX_am_header_t ucx_hdr;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_NM_AM_ISEND);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_NM_AM_ISEND);

    ep = MPIDI_UCX_COMM_TO_EP(comm, rank, 0, 0);

    int dt_contig;
    size_t data_sz;
    MPI_Aint dt_true_lb;
    MPIR_Datatype *dt_ptr;
    MPIDI_Datatype_get_info(count, datatype, dt_contig, data_sz, dt_ptr, dt_true_lb);

    /* initialize our portion of the hdr */
    ucx_hdr.handler_id = handler_id;
    ucx_hdr.data_sz = data_sz;

#ifdef HAVE_UCP_AM_NBX
    size_t header_size = sizeof(ucx_hdr) + am_hdr_sz;
    void *send_buf, *header, *data_ptr;
    /* note: since we are not copying large contig gpu data, it is less useful
     * to use MPIR_gpu_malloc_host */
    if (dt_contig) {
        /* only need copy headers */
        send_buf = MPL_malloc(header_size, MPL_MEM_OTHER);
        MPIR_Assert(send_buf);
        header = send_buf;

        MPIR_Memcpy(header, &ucx_hdr, sizeof(ucx_hdr));
        MPIR_Memcpy((char *) header + sizeof(ucx_hdr), am_hdr, am_hdr_sz);

        data_ptr = (char *) data + dt_true_lb;
    } else {
        /* need copy headers and pack data */
        send_buf = MPL_malloc(header_size + data_sz, MPL_MEM_OTHER);
        MPIR_Assert(send_buf);
        header = send_buf;
        data_ptr = (char *) send_buf + header_size;

        MPIR_Memcpy(header, &ucx_hdr, sizeof(ucx_hdr));
        MPIR_Memcpy((char *) header + sizeof(ucx_hdr), am_hdr, am_hdr_sz);

        MPI_Aint actual_pack_bytes;
        mpi_errno = MPIR_Typerep_pack(data, count, datatype, 0, data_ptr, data_sz,
                                      &actual_pack_bytes);
        MPIR_ERR_CHECK(mpi_errno);
        MPIR_Assert(actual_pack_bytes == data_sz);
    }
    ucp_request_param_t param = {
        .op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK,
        .cb.send = &MPIDI_UCX_am_isend_callback_nbx,
    };
    ucp_request = (MPIDI_UCX_ucp_request_t *) ucp_am_send_nbx(ep, MPIDI_UCX_AM_NBX_HANDLER_ID,
                                                              header, header_size,
                                                              data_ptr, data_sz, &param);
    MPIDI_UCX_CHK_REQUEST(ucp_request);
    /* if send is done, free all resources and complete the request */
    if (ucp_request == NULL) {
        MPL_free(send_buf);
        MPIDIG_global.origin_cbs[handler_id] (sreq);
        goto fn_exit;
    }

    MPIDI_UCX_AM_SEND_REQUEST(sreq, pack_buffer) = send_buf;
    MPIDI_UCX_AM_SEND_REQUEST(sreq, handler_id) = handler_id;
    ucp_request->req = sreq;
    ucp_request_release(ucp_request);

#else /* !HAVE_UCP_AM_NBX */

    MPL_pointer_attr_t attr;
    MPIR_GPU_query_pointer_attr(data, &attr);
    if (attr.type == MPL_GPU_POINTER_DEV) {
        /* Force packing of GPU buffer in host memory */
        dt_contig = 0;
    }

    char *send_buf;
    void *send_buf_p;
    size_t total_sz;
    ucp_datatype_t dt;
    if (dt_contig) {
        MPIR_gpu_malloc_host((void **) &send_buf, sizeof(ucx_hdr) + am_hdr_sz);
        MPIR_Memcpy(send_buf, &ucx_hdr, sizeof(ucx_hdr));
        MPIR_Memcpy(send_buf + sizeof(ucx_hdr), am_hdr, am_hdr_sz);

        ucp_dt_iov_t *iov = MPIDI_UCX_AM_SEND_REQUEST(sreq, iov);
        iov[0].buffer = send_buf;
        iov[0].length = sizeof(ucx_hdr) + am_hdr_sz;
        iov[1].buffer = (char *) data + dt_true_lb;
        iov[1].length = data_sz;

        send_buf_p = iov;
        dt = ucp_dt_make_iov();
        total_sz = 2;
    } else {
        MPIR_gpu_malloc_host((void **) &send_buf, data_sz + am_hdr_sz + sizeof(ucx_hdr));
        MPIR_Memcpy(send_buf, &ucx_hdr, sizeof(ucx_hdr));
        MPIR_Memcpy(send_buf + sizeof(ucx_hdr), am_hdr, am_hdr_sz);

        MPI_Aint actual_pack_bytes;
        mpi_errno =
            MPIR_Typerep_pack(data, count, datatype, 0, send_buf + am_hdr_sz + sizeof(ucx_hdr),
                              data_sz, &actual_pack_bytes);
        MPIR_ERR_CHECK(mpi_errno);
        MPIR_Assert(actual_pack_bytes == data_sz);

        send_buf_p = send_buf;
        dt = ucp_dt_make_contig(1);
        total_sz = data_sz + am_hdr_sz + sizeof(ucx_hdr);
    }

    ucp_request = (MPIDI_UCX_ucp_request_t *) ucp_am_send_nb(ep, MPIDI_UCX_AM_HANDLER_ID,
                                                             send_buf_p, total_sz, dt,
                                                             &MPIDI_UCX_am_isend_callback, 0);
    MPIDI_UCX_CHK_REQUEST(ucp_request);

    /* send is done. free all resources and complete the request */
    if (ucp_request == NULL) {
        MPIR_gpu_free_host(send_buf);
        MPIDIG_global.origin_cbs[handler_id] (sreq);
        goto fn_exit;
    }

    /* set the ch4r request inside the UCP request */
    MPIDI_UCX_AM_SEND_REQUEST(sreq, pack_buffer) = send_buf;
    MPIDI_UCX_AM_SEND_REQUEST(sreq, handler_id) = handler_id;
    ucp_request->req = sreq;
    ucp_request_release(ucp_request);

#endif /* HAVE_UCP_AM_NBX */

  fn_exit:
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_NM_AM_ISEND);
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_am_isend_reply(MPIR_Comm * comm,
                                                     int src_rank,
                                                     int handler_id,
                                                     const void *am_hdr,
                                                     MPI_Aint am_hdr_sz,
                                                     const void *data, MPI_Aint count,
                                                     MPI_Datatype datatype, MPIR_Request * sreq)
{
    int mpi_errno = MPI_SUCCESS;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_NM_AM_ISEND_REPLY);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_NM_AM_ISEND_REPLY);

    mpi_errno = MPIDI_NM_am_isend(src_rank, comm, handler_id,
                                  am_hdr, am_hdr_sz, data, count, datatype, sreq);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_NM_AM_ISEND_REPLY);
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX MPI_Aint MPIDI_NM_am_hdr_max_sz(void)
{
    return (MPIDI_UCX_MAX_AM_EAGER_SZ - sizeof(MPIDI_UCX_am_header_t));
}

MPL_STATIC_INLINE_PREFIX MPI_Aint MPIDI_NM_am_eager_limit(void)
{
    return (MPIDI_UCX_MAX_AM_EAGER_SZ - sizeof(MPIDI_UCX_am_header_t));
}

MPL_STATIC_INLINE_PREFIX MPI_Aint MPIDI_NM_am_eager_buf_limit(void)
{
    return MPIDI_UCX_MAX_AM_EAGER_SZ;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_am_send_hdr(int rank,
                                                  MPIR_Comm * comm,
                                                  int handler_id, const void *am_hdr,
                                                  MPI_Aint am_hdr_sz)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_UCX_ucp_request_t *ucp_request;
    ucp_ep_h ep;
    char *send_buf;
    MPIDI_UCX_am_header_t ucx_hdr;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_NM_AM_SEND_HDR);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_NM_AM_SEND_HDR);

    ep = MPIDI_UCX_COMM_TO_EP(comm, rank, 0, 0);

    /* initialize our portion of the hdr */
    ucx_hdr.handler_id = handler_id;
    ucx_hdr.data_sz = 0;

    /* just pack and send for now */
    send_buf = MPL_malloc(am_hdr_sz + sizeof(ucx_hdr), MPL_MEM_BUFFER);
    MPIR_Memcpy(send_buf, &ucx_hdr, sizeof(ucx_hdr));
    MPIR_Memcpy(send_buf + sizeof(ucx_hdr), am_hdr, am_hdr_sz);

    ucp_request = (MPIDI_UCX_ucp_request_t *) ucp_am_send_nb(ep, MPIDI_UCX_AM_HANDLER_ID, send_buf,
                                                             am_hdr_sz + sizeof(ucx_hdr),
                                                             ucp_dt_make_contig(1),
                                                             &MPIDI_UCX_am_send_callback, 0);
    MPIDI_UCX_CHK_REQUEST(ucp_request);

    if (ucp_request == NULL) {
        /* inject is done */
        MPL_free(send_buf);
    } else {
        ucp_request->buf = send_buf;
    }

  fn_exit:
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_NM_AM_SEND_HDR);
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_am_send_hdr_reply(MPIR_Comm * comm,
                                                        int rank,
                                                        int handler_id, const void *am_hdr,
                                                        MPI_Aint am_hdr_sz)
{
    int mpi_errno = MPI_SUCCESS;
    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_NM_AM_SEND_HDR_REPLY);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_NM_AM_SEND_HDR_REPLY);

    mpi_errno = MPIDI_NM_am_send_hdr(rank, comm, handler_id, am_hdr, am_hdr_sz);

    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_NM_AM_SEND_HDR_REPLY);
    return mpi_errno;
}

MPL_STATIC_INLINE_PREFIX bool MPIDI_NM_am_check_eager(MPI_Aint am_hdr_sz, MPI_Aint data_sz,
                                                      const void *data, MPI_Aint count,
                                                      MPI_Datatype datatype, MPIR_Request * sreq)
{
#ifdef HAVE_UCP_AM_NBX
    /* ucx will handle rndv */
    return true;
#else
    return (am_hdr_sz + data_sz) <= (MPIDI_UCX_MAX_AM_EAGER_SZ - sizeof(MPIDI_UCX_am_header_t));
#endif
}

#ifdef HAVE_UCP_AM_NBX
MPL_STATIC_INLINE_PREFIX void MPIDI_UCX_do_am_recv(MPIR_Request * rreq)
{
    void *recv_buf;
    bool is_contig;
    MPI_Aint data_sz, in_data_sz;

    MPIDIG_get_recv_buffer(&recv_buf, &data_sz, &is_contig, &in_data_sz, rreq);
    if (!is_contig || in_data_sz > data_sz) {
        /* non-contig datatype, need receive into pack buffer */
        /* ucx will error out if buffer size is less than the promised data size,
         * also use a pack buffer in this case */
        recv_buf = MPL_malloc(in_data_sz, MPL_MEM_OTHER);
        MPIR_Assert(recv_buf);
        MPIDI_UCX_AM_RECV_REQUEST(rreq, pack_buffer) = recv_buf;
    } else {
        MPIDI_UCX_AM_RECV_REQUEST(rreq, pack_buffer) = NULL;
    }

    MPIDI_UCX_ucp_request_t *ucp_request;
    size_t received_length;
    ucp_request_param_t param = {
        .op_attr_mask = UCP_OP_ATTR_FIELD_CALLBACK | UCP_OP_ATTR_FIELD_RECV_INFO,
        .cb.recv_am = &MPIDI_UCX_am_recv_callback_nbx,
        .recv_info.length = &received_length,
    };
    void *data_desc = MPIDI_UCX_AM_RECV_REQUEST(rreq, data_desc);
    /* note: use in_data_sz to match promised data size */
    ucp_request = ucp_am_recv_data_nbx(MPIDI_UCX_global.ctx[0].worker,
                                       data_desc, recv_buf, in_data_sz, &param);
    if (ucp_request == NULL) {
        /* completed immediately */
        MPIDI_UCX_ucp_request_t tmp_ucp_request;
        tmp_ucp_request.req = rreq;
        MPIDI_UCX_am_recv_callback_nbx(&tmp_ucp_request, UCS_OK, received_length, NULL);
    } else {
        ucp_request->req = rreq;
    }
}
#endif

MPL_STATIC_INLINE_PREFIX int MPIDI_NM_am_recv(MPIR_Request * rreq)
{
    int ret = MPI_SUCCESS;

    MPIR_FUNC_VERBOSE_STATE_DECL(MPID_STATE_MPIDI_NM_AM_RECV);
    MPIR_FUNC_VERBOSE_ENTER(MPID_STATE_MPIDI_NM_AM_RECV);
#ifdef HAVE_UCP_AM_NBX
    MPIDI_UCX_do_am_recv(rreq);
#endif
    MPIR_FUNC_VERBOSE_EXIT(MPID_STATE_MPIDI_NM_AM_RECV);
    return ret;
}

#endif /* UCX_AM_H_INCLUDED */
