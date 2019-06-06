/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2016 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Mellanox Technologies Ltd.
 *  Copyright (C) Mellanox Technologies Ltd. 2016. ALL RIGHTS RESERVED
 */

#include "mpidimpl.h"
#include "ucx_impl.h"
#include <ucp/api/ucp.h>
#ifdef HAVE_LIBHCOLL
#include "../../../common/hcoll/hcoll.h"
#endif

struct pack_state {
    void *buffer;
    MPI_Aint count;
    MPI_Datatype datatype;
};

static void *start_pack(void *context, const void *buffer, MPI_Aint count);
static void *start_unpack(void *context, void *buffer, MPI_Aint count);
static MPI_Aint pack_size(void *state);
static MPI_Aint pack(void *state, MPI_Aint offset, void *dest, MPI_Aint max_length);
static ucs_status_t unpack(void *state, MPI_Aint offset, const void *src, MPI_Aint count);
static void finish_pack(void *state);

static void *start_pack(void *context, const void *buffer, MPI_Aint count)
{
    struct pack_state *state;

    /* Todo: Add error handling */
    state = MPL_malloc(sizeof(struct pack_state), MPL_MEM_DATATYPE);

    state->buffer = (void *) buffer;
    state->count = count;
    state->datatype = *((MPI_Datatype *) context);

    return (void *) state;
}

static void *start_unpack(void *context, void *buffer, MPI_Aint count)
{
    struct pack_state *state;

    /* Todo: Add error handling */
    state = MPL_malloc(sizeof(struct pack_state), MPL_MEM_DATATYPE);

    state->buffer = buffer;
    state->count = count;
    state->datatype = *((MPI_Datatype *) context);

    return (void *) state;
}

static MPI_Aint pack_size(void *state)
{
    struct pack_state *pack_state = (struct pack_state *) state;
    MPI_Aint packsize;

    MPIR_Pack_size_impl(pack_state->count, pack_state->datatype, &packsize);

    return (MPI_Aint) packsize;
}

static MPI_Aint pack(void *state, MPI_Aint offset, void *dest, MPI_Aint max_length)
{
    struct pack_state *pack_state = (struct pack_state *) state;
    MPI_Aint actual_pack_bytes;

    MPIR_Typerep_pack(pack_state->buffer, pack_state->count, pack_state->datatype, offset,
                      dest, max_length, &actual_pack_bytes);

    return actual_pack_bytes;
}

static ucs_status_t unpack(void *state, MPI_Aint offset, const void *src, MPI_Aint count)
{
    struct pack_state *pack_state = (struct pack_state *) state;
    MPI_Aint max_unpack_bytes;
    MPI_Aint actual_unpack_bytes;
    MPI_Aint packsize;

    MPIR_Pack_size_impl(pack_state->count, pack_state->datatype, &packsize);
    max_unpack_bytes = MPL_MIN(packsize, count);

    MPIR_Typerep_unpack(src, max_unpack_bytes, pack_state->buffer, pack_state->count,
                        pack_state->datatype, offset, &actual_unpack_bytes);
    if (unlikely(actual_unpack_bytes != max_unpack_bytes)) {
        return UCS_ERR_MESSAGE_TRUNCATED;
    }

    return UCS_OK;
}

static void finish_pack(void *state)
{
    MPIR_Datatype *dt_ptr;
    struct pack_state *pack_state = (struct pack_state *) state;
    MPIR_Datatype_get_ptr(pack_state->datatype, dt_ptr);
    MPIR_Datatype_ptr_release(dt_ptr);
    MPL_free(pack_state);
}

int MPIDI_UCX_mpi_type_free_hook(MPIR_Datatype * datatype_p)
{
    if (datatype_p->is_committed && (int) datatype_p->dev.netmod.ucx.ucp_datatype >= 0) {
        ucp_dt_destroy(datatype_p->dev.netmod.ucx.ucp_datatype);
        datatype_p->dev.netmod.ucx.ucp_datatype = -1;
    }
#ifdef HAVE_LIBHCOLL
    hcoll_type_free_hook(datatype_p);
#endif

    return 0;
}

int MPIDI_UCX_mpi_type_commit_hook(MPIR_Datatype * datatype_p)
{
    ucp_datatype_t ucp_datatype;
    ucs_status_t status;
    int is_contig;

    datatype_p->dev.netmod.ucx.ucp_datatype = -1;
    MPIR_Datatype_is_contig(datatype_p->handle, &is_contig);

    if (!is_contig) {

        ucp_generic_dt_ops_t dt_ops = {
            .start_pack = start_pack,
            .start_unpack = start_unpack,
            .packed_size = pack_size,
            .pack = pack,
            .unpack = unpack,
            .finish = finish_pack
        };
        status = ucp_dt_create_generic(&dt_ops, datatype_p, &ucp_datatype);
        MPIR_Assertp(status == UCS_OK);
        datatype_p->dev.netmod.ucx.ucp_datatype = ucp_datatype;

    }
#ifdef HAVE_LIBHCOLL
    hcoll_type_commit_hook(datatype_p);
#endif

    return 0;
}
