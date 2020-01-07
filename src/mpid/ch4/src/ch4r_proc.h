/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef CH4R_PROC_H_INCLUDED
#define CH4R_PROC_H_INCLUDED

#include "ch4_types.h"

int MPIDIU_get_n_avts(void);
int MPIDIU_get_max_n_avts(void);
int MPIDIU_get_avt_size(int avtid);
int MPIDIU_alloc_globals_for_avtid(int avtid);
int MPIDIU_free_globals_for_avtid(int avtid);
int MPIDIU_get_next_avtid(int *avtid);
int MPIDIU_free_avtid(int avtid);
int MPIDIU_new_avt(int size, int *avtid);
int MPIDIU_free_avt(int avtid);
int MPIDIU_avt_add_ref(int avtid);
int MPIDIU_avt_release_ref(int avtid);
int MPIDIU_avt_init(void);
int MPIDIU_avt_destroy(void);
int MPIDIU_get_node_id(MPIR_Comm * comm, int rank, int *id_p);
int MPIDIU_get_max_node_id(MPIR_Comm * comm, int *max_id_p);
int MPIDIU_build_nodemap(int myrank, MPIR_Comm * comm, int sz, int *out_nodemap, int *sz_out);
int MPIDIU_build_nodemap_avtid(int myrank, MPIR_Comm * comm, int sz, int avtid);

MPL_STATIC_INLINE_PREFIX int MPIDIU_comm_rank_to_pid(MPIR_Comm * comm, int rank, int *idx,
                                                     int *avtid)
{



    *avtid = 0;

    switch (MPIDI_COMM(comm, map).mode) {
        case MPIDI_RANK_MAP_DIRECT:
            *avtid = MPIDI_COMM(comm, map).avtid;
            *idx = rank;
            break;
        case MPIDI_RANK_MAP_DIRECT_INTRA:
            *idx = rank;
            break;
        case MPIDI_RANK_MAP_OFFSET:
            *avtid = MPIDI_COMM(comm, map).avtid;
            *idx = rank + MPIDI_COMM(comm, map).reg.offset;
            break;
        case MPIDI_RANK_MAP_OFFSET_INTRA:
            *idx = rank + MPIDI_COMM(comm, map).reg.offset;
            break;
        case MPIDI_RANK_MAP_STRIDE:
            *avtid = MPIDI_COMM(comm, map).avtid;
            *idx = MPIDI_CALC_STRIDE_SIMPLE(rank, MPIDI_COMM(comm, map).reg.stride.stride,
                                            MPIDI_COMM(comm, map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_STRIDE_INTRA:
            *idx = MPIDI_CALC_STRIDE_SIMPLE(rank, MPIDI_COMM(comm, map).reg.stride.stride,
                                            MPIDI_COMM(comm, map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_STRIDE_BLOCK:
            *avtid = MPIDI_COMM(comm, map).avtid;
            *idx = MPIDI_CALC_STRIDE(rank, MPIDI_COMM(comm, map).reg.stride.stride,
                                     MPIDI_COMM(comm, map).reg.stride.blocksize,
                                     MPIDI_COMM(comm, map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_STRIDE_BLOCK_INTRA:
            *idx = MPIDI_CALC_STRIDE(rank, MPIDI_COMM(comm, map).reg.stride.stride,
                                     MPIDI_COMM(comm, map).reg.stride.blocksize,
                                     MPIDI_COMM(comm, map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_LUT:
            *avtid = MPIDI_COMM(comm, map).avtid;
            *idx = MPIDI_COMM(comm, map).irreg.lut.lpid[rank];
            break;
        case MPIDI_RANK_MAP_LUT_INTRA:
            *idx = MPIDI_COMM(comm, map).irreg.lut.lpid[rank];
            break;
        case MPIDI_RANK_MAP_MLUT:
            *idx = MPIDI_COMM(comm, map).irreg.mlut.gpid[rank].lpid;
            *avtid = MPIDI_COMM(comm, map).irreg.mlut.gpid[rank].avtid;
            break;
        case MPIDI_RANK_MAP_NONE:
            MPIR_Assert(0);
            break;
    }
    MPL_DBG_MSG_FMT(MPIDI_CH4_DBG_MAP, VERBOSE,
                    (MPL_DBG_FDEST, " comm_to_pid: rank=%d, avtid=%d idx=%d", rank, *avtid, *idx));

    return *idx;
}

MPL_STATIC_INLINE_PREFIX MPIDI_av_entry_t *MPIDIU_comm_rank_to_av(MPIR_Comm * comm, int rank)
{
    MPIDI_av_entry_t *ret = NULL;



    switch (MPIDI_COMM(comm, map).mode) {
        case MPIDI_RANK_MAP_DIRECT:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).avtid]->table[rank];
            break;
        case MPIDI_RANK_MAP_DIRECT_INTRA:
            ret = &MPIDI_av_table0->table[rank];
            break;
        case MPIDI_RANK_MAP_OFFSET:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).avtid]
                ->table[rank + MPIDI_COMM(comm, map).reg.offset];
            break;
        case MPIDI_RANK_MAP_OFFSET_INTRA:
            ret = &MPIDI_av_table0->table[rank + MPIDI_COMM(comm, map).reg.offset];
            break;
        case MPIDI_RANK_MAP_STRIDE:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).avtid]
                ->table[MPIDI_CALC_STRIDE_SIMPLE(rank,
                                                 MPIDI_COMM(comm, map).reg.stride.stride,
                                                 MPIDI_COMM(comm, map).reg.stride.offset)];
            break;
        case MPIDI_RANK_MAP_STRIDE_INTRA:
            ret = &MPIDI_av_table0->table[MPIDI_CALC_STRIDE_SIMPLE(rank,
                                                                   MPIDI_COMM(comm,
                                                                              map).reg.
                                                                   stride.stride, MPIDI_COMM(comm,
                                                                                             map).
                                                                   reg.stride.offset)];
            break;
        case MPIDI_RANK_MAP_STRIDE_BLOCK:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).avtid]
                ->table[MPIDI_CALC_STRIDE(rank,
                                          MPIDI_COMM(comm, map).reg.stride.stride,
                                          MPIDI_COMM(comm, map).reg.stride.blocksize,
                                          MPIDI_COMM(comm, map).reg.stride.offset)];
            break;
        case MPIDI_RANK_MAP_STRIDE_BLOCK_INTRA:
            ret = &MPIDI_av_table0->table[MPIDI_CALC_STRIDE(rank,
                                                            MPIDI_COMM(comm, map).reg.stride.stride,
                                                            MPIDI_COMM(comm,
                                                                       map).reg.stride.blocksize,
                                                            MPIDI_COMM(comm,
                                                                       map).reg.stride.offset)];
            break;
        case MPIDI_RANK_MAP_LUT:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).avtid]
                ->table[MPIDI_COMM(comm, map).irreg.lut.lpid[rank]];
            break;
        case MPIDI_RANK_MAP_LUT_INTRA:
            ret = &MPIDI_av_table0->table[MPIDI_COMM(comm, map).irreg.lut.lpid[rank]];
            break;
        case MPIDI_RANK_MAP_MLUT:
            ret = &MPIDI_av_table[MPIDI_COMM(comm, map).irreg.mlut.gpid[rank].avtid]
                ->table[MPIDI_COMM(comm, map).irreg.mlut.gpid[rank].lpid];
            break;
        case MPIDI_RANK_MAP_NONE:
            MPIR_Assert(0);
            break;
    }

    MPL_DBG_MSG_FMT(MPIDI_CH4_DBG_MAP, VERBOSE,
                    (MPL_DBG_FDEST, " comm_to_av_addr: rank=%d, av addr=%p", rank, (void *) ret));

    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDIU_comm_rank_to_pid_local(MPIR_Comm * comm, int rank, int *idx,
                                                           int *avtid)
{



    *avtid = MPIDI_COMM(comm, local_map).avtid;
    switch (MPIDI_COMM(comm, local_map).mode) {
        case MPIDI_RANK_MAP_DIRECT:
        case MPIDI_RANK_MAP_DIRECT_INTRA:
            *idx = rank;
            break;
        case MPIDI_RANK_MAP_OFFSET:
        case MPIDI_RANK_MAP_OFFSET_INTRA:
            *idx = rank + MPIDI_COMM(comm, local_map).reg.offset;
            break;
        case MPIDI_RANK_MAP_STRIDE:
        case MPIDI_RANK_MAP_STRIDE_INTRA:
            *idx = MPIDI_CALC_STRIDE_SIMPLE(rank, MPIDI_COMM(comm, local_map).reg.stride.stride,
                                            MPIDI_COMM(comm, local_map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_STRIDE_BLOCK:
        case MPIDI_RANK_MAP_STRIDE_BLOCK_INTRA:
            *idx = MPIDI_CALC_STRIDE(rank, MPIDI_COMM(comm, local_map).reg.stride.stride,
                                     MPIDI_COMM(comm, local_map).reg.stride.blocksize,
                                     MPIDI_COMM(comm, local_map).reg.stride.offset);
            break;
        case MPIDI_RANK_MAP_LUT:
        case MPIDI_RANK_MAP_LUT_INTRA:
            *idx = MPIDI_COMM(comm, local_map).irreg.lut.lpid[rank];
            break;
        case MPIDI_RANK_MAP_MLUT:
            *idx = MPIDI_COMM(comm, local_map).irreg.mlut.gpid[rank].lpid;
            *avtid = MPIDI_COMM(comm, local_map).irreg.mlut.gpid[rank].avtid;
            break;
        case MPIDI_RANK_MAP_NONE:
            MPIR_Assert(0);
            break;
    }
    MPL_DBG_MSG_FMT(MPIDI_CH4_DBG_MAP, VERBOSE,
                    (MPL_DBG_FDEST, " comm_to_pid_local: rank=%d, avtid=%d idx=%d",
                     rank, *avtid, *idx));

    return *idx;
}

MPL_STATIC_INLINE_PREFIX int MPIDIU_rank_is_local(int rank, MPIR_Comm * comm)
{
    int ret = 0;



#ifdef MPIDI_BUILD_CH4_LOCALITY_INFO
    ret = MPIDIU_comm_rank_to_av(comm, rank)->is_local;
    MPL_DBG_MSG_FMT(MPIDI_CH4_DBG_MAP, VERBOSE,
                    (MPL_DBG_FDEST, " is_local=%d, rank=%d", ret, rank));
#endif


    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDIU_av_is_local(MPIDI_av_entry_t * av)
{
    int ret = 0;



#ifdef MPIDI_BUILD_CH4_LOCALITY_INFO
    ret = av->is_local;
    MPL_DBG_MSG_FMT(MPIDI_CH4_DBG_MAP, VERBOSE,
                    (MPL_DBG_FDEST, " is_local=%d, av=%p", ret, (void *) av));
#endif


    return ret;
}

MPL_STATIC_INLINE_PREFIX int MPIDIU_rank_to_lpid(int rank, MPIR_Comm * comm)
{
    int ret;



    int avtid = 0, lpid = 0;
    MPIDIU_comm_rank_to_pid(comm, rank, &lpid, &avtid);
    if (avtid == 0) {
        ret = lpid;
    } else {
        ret = -1;
    }


    return ret;
}

#endif /* CH4R_PROC_H_INCLUDED */
