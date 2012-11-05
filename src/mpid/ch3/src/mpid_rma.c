/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2001 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "mpidimpl.h"
#include "mpidrma.h"


MPIU_THREADSAFE_INIT_DECL(initRMAoptions);
#ifdef USE_MPIU_INSTR
MPIU_INSTR_DURATION_DECL(wincreate_allgather);
MPIU_INSTR_DURATION_DECL(winfree_rs);
MPIU_INSTR_DURATION_DECL(winfree_complete);
MPIU_INSTR_DURATION_DECL(rmaqueue_alloc);
MPIU_INSTR_DURATION_DECL(rmaqueue_set);
extern void MPIDI_CH3_RMA_InitInstr(void);
#endif


static int win_init(MPI_Aint size, int disp_unit, int create_flavor, int model,
                    MPID_Comm *comm_ptr, MPID_Win **win_ptr);


#define MPID_WIN_FTABLE_SET_DEFAULTS(win_ptr)                   \
    do {                                                        \
        /* Get ptr to RMAFns, which is embedded in MPID_Win */  \
        MPID_RMAFns *ftable         = &(*(win_ptr))->RMAFns;    \
        ftable->Win_free            = MPIDI_Win_free;           \
        ftable->Win_attach          = MPIDI_Win_attach;         \
        ftable->Win_detach          = MPIDI_Win_detach;         \
        ftable->Win_shared_query    = MPIDI_Win_shared_query;   \
                                                                \
        ftable->Win_set_info        = MPIDI_Win_set_info;       \
        ftable->Win_get_info        = MPIDI_Win_get_info;       \
                                                                \
        ftable->Put                 = MPIDI_Put;                \
        ftable->Get                 = MPIDI_Get;                \
        ftable->Accumulate          = MPIDI_Accumulate;         \
        ftable->Get_accumulate      = MPIDI_Get_accumulate;     \
        ftable->Fetch_and_op        = MPIDI_Fetch_and_op;       \
        ftable->Compare_and_swap    = MPIDI_Compare_and_swap;   \
                                                                \
        ftable->Rput                = MPIDI_Rput;               \
        ftable->Rget                = MPIDI_Rget;               \
        ftable->Raccumulate         = MPIDI_Raccumulate;        \
        ftable->Rget_accumulate     = MPIDI_Rget_accumulate;    \
                                                                \
        ftable->Win_fence           = MPIDI_Win_fence;          \
        ftable->Win_post            = MPIDI_Win_post;           \
        ftable->Win_start           = MPIDI_Win_start;          \
        ftable->Win_complete        = MPIDI_Win_complete;       \
        ftable->Win_wait            = MPIDI_Win_wait;           \
        ftable->Win_test            = MPIDI_Win_test;           \
                                                                \
        ftable->Win_lock            = MPIDI_Win_lock;           \
        ftable->Win_unlock          = MPIDI_Win_unlock;         \
        ftable->Win_lock_all        = MPIDI_Win_lock_all;       \
        ftable->Win_unlock_all      = MPIDI_Win_unlock_all;     \
                                                                \
        ftable->Win_flush           = MPIDI_Win_flush;          \
        ftable->Win_flush_all       = MPIDI_Win_flush_all;      \
        ftable->Win_flush_local     = MPIDI_Win_flush_local;    \
        ftable->Win_flush_local_all = MPIDI_Win_flush_local_all;\
        ftable->Win_sync            = MPIDI_Win_sync;           \
    } while (0)


#undef FUNCNAME
#define FUNCNAME MPID_Win_create
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_create(void *base, MPI_Aint size, int disp_unit, MPID_Info *info, 
                    MPID_Comm *comm_ptr, MPID_Win **win_ptr)
{
    int mpi_errno=MPI_SUCCESS;

    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_CREATE);
    
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPID_WIN_CREATE);

    mpi_errno = win_init(size, disp_unit, MPI_WIN_FLAVOR_CREATE, MPI_WIN_SEPARATE, comm_ptr, win_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    (*win_ptr)->base = base;

    mpi_errno = MPIDI_CH3U_Win_fns.create(base, size, disp_unit, info, comm_ptr, win_ptr); 
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

 fn_fail:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_CREATE);
    return mpi_errno;
}


#undef FUNCNAME
#define FUNCNAME MPID_Win_allocate
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_allocate(MPI_Aint size, int disp_unit, MPID_Info *info, 
                    MPID_Comm *comm_ptr, void *baseptr, MPID_Win **win_ptr)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_ALLOCATE);
    
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPID_WIN_ALLOCATE);

    mpi_errno = win_init(size, disp_unit, MPI_WIN_FLAVOR_ALLOCATE, MPI_WIN_SEPARATE, comm_ptr, win_ptr);
    if (mpi_errno != MPI_SUCCESS) { MPIU_ERR_POP(mpi_errno); }

    mpi_errno = MPIDI_CH3U_Win_fns.allocate(size, disp_unit, info, comm_ptr, baseptr, win_ptr);
    if (mpi_errno != MPI_SUCCESS) { MPIU_ERR_POP(mpi_errno); }

 fn_fail:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_ALLOCATE);
    return mpi_errno;
}


#undef FUNCNAME
#define FUNCNAME MPID_Win_create_dynamic
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_create_dynamic(MPID_Info *info, MPID_Comm *comm_ptr,
                            MPID_Win **win_ptr)
{
    int mpi_errno=MPI_SUCCESS;

    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_CREATE_DYNAMIC);

    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPID_WIN_CREATE_DYNAMIC);

    mpi_errno = win_init(0 /* spec defines size to be 0 */,
                         1 /* spec defines disp_unit to be 1 */,
                         MPI_WIN_FLAVOR_DYNAMIC, MPI_WIN_SEPARATE,
                         comm_ptr, win_ptr);

    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    (*win_ptr)->base = MPI_BOTTOM;

    mpi_errno = MPIDI_CH3U_Win_fns.create_dynamic(info, comm_ptr, win_ptr);
    if (mpi_errno != MPI_SUCCESS) { MPIU_ERR_POP(mpi_errno); }

 fn_fail:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_CREATE_DYNAMIC);
    return mpi_errno;
}


/* The memory allocation functions */
#undef FUNCNAME
#define FUNCNAME MPID_Alloc_mem
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
void *MPID_Alloc_mem( size_t size, MPID_Info *info_ptr )
{
    void *ap=NULL;
    MPIDI_STATE_DECL(MPID_STATE_MPID_ALLOC_MEM);

    MPIDI_FUNC_ENTER(MPID_STATE_MPID_ALLOC_MEM);

    ap = MPIDI_Alloc_mem(size, info_ptr);
    
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_ALLOC_MEM);
    return ap;
}


#undef FUNCNAME
#define FUNCNAME MPID_Free_mem
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Free_mem( void *ptr )
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_STATE_DECL(MPID_STATE_MPID_FREE_MEM);

    MPIDI_FUNC_ENTER(MPID_STATE_MPID_FREE_MEM);

    mpi_errno = MPIDI_Free_mem(ptr);
    if (mpi_errno != MPI_SUCCESS) {
	MPIU_ERR_POP(mpi_errno);
    }
        
 fn_fail:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_FREE_MEM);
    return mpi_errno;
}


#undef FUNCNAME
#define FUNCNAME MPID_Win_allocate_shared
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_allocate_shared(MPI_Aint size, int disp_unit, MPID_Info *info, MPID_Comm *comm_ptr,
                             void **base_ptr, MPID_Win **win_ptr)
{
    int mpi_errno=MPI_SUCCESS;

    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_ALLOCATE_SHARED);
    
    MPIDI_RMA_FUNC_ENTER(MPID_STATE_MPID_WIN_ALLOCATE_SHARED);

    mpi_errno = win_init(size, disp_unit, MPI_WIN_FLAVOR_SHARED, MPI_WIN_UNIFIED, comm_ptr, win_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    mpi_errno = MPIDI_CH3U_Win_fns.allocate_shared(size, disp_unit, info, comm_ptr, base_ptr, win_ptr);
    if (mpi_errno != MPI_SUCCESS) MPIU_ERR_POP(mpi_errno);

 fn_fail:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_ALLOCATE_SHARED);
    return mpi_errno;
}


#undef FUNCNAME
#define FUNCNAME win_init
#undef FCNAME
#define FCNAME MPIU_QUOTE(FUNCNAME)
static int win_init(MPI_Aint size, int disp_unit, int create_flavor, int model,
                          MPID_Comm *comm_ptr, MPID_Win **win_ptr)
{
    int mpi_errno = MPI_SUCCESS;
    int i;
    MPID_Comm *win_comm_ptr;
    MPIU_CHKPMEM_DECL(1);
    MPIDI_STATE_DECL(MPID_STATE_WIN_INIT);

    MPIDI_FUNC_ENTER(MPID_STATE_WIN_INIT);

    if(initRMAoptions) {
        MPIU_THREADSAFE_INIT_BLOCK_BEGIN(initRMAoptions);
#ifdef USE_MPIU_INSTR
        /* Define all instrumentation handles used in the CH3 RMA here*/
        MPIU_INSTR_DURATION_INIT(wincreate_allgather,0,"WIN_CREATE:Allgather");
        MPIU_INSTR_DURATION_INIT(winfree_rs,0,"WIN_FREE:ReduceScatterBlock");
        MPIU_INSTR_DURATION_INIT(winfree_complete,0,"WIN_FREE:Complete");
        MPIU_INSTR_DURATION_INIT(rmaqueue_alloc,0,"Allocate RMA Queue element");
        MPIU_INSTR_DURATION_INIT(rmaqueue_set,0,"Set fields in RMA Queue element");
        MPIDI_CH3_RMA_InitInstr();
#endif

        MPIU_THREADSAFE_INIT_CLEAR(initRMAoptions);
        MPIU_THREADSAFE_INIT_BLOCK_END(initRMAoptions);
    }

    *win_ptr = (MPID_Win *)MPIU_Handle_obj_alloc( &MPID_Win_mem );
    MPIU_ERR_CHKANDJUMP1(!(*win_ptr),mpi_errno,MPI_ERR_OTHER,"**nomem",
                         "**nomem %s","MPID_Win_mem");

    mpi_errno = MPIR_Comm_dup_impl(comm_ptr, &win_comm_ptr);
    if (mpi_errno) MPIU_ERR_POP(mpi_errno);

    MPIU_Object_set_ref(*win_ptr, 1);

    (*win_ptr)->fence_cnt           = 0;
    /* (*win_ptr)->errhandler is set by upper level; */
    /* (*win_ptr)->base is set by caller; */
    (*win_ptr)->size                = size;
    (*win_ptr)->disp_unit           = disp_unit;
    (*win_ptr)->create_flavor       = create_flavor;
    (*win_ptr)->model               = model;
    (*win_ptr)->attributes          = NULL;
    (*win_ptr)->start_group_ptr     = NULL;
    (*win_ptr)->start_assert        = 0;
    (*win_ptr)->comm_ptr            = win_comm_ptr;
    (*win_ptr)->myrank              = comm_ptr->rank;

    (*win_ptr)->my_counter          = 0;
    /* (*win_ptr)->base_addrs[] is set by caller; */
    /* (*win_ptr)->sizes[] is set by caller; */
    /* (*win_ptr)->disp_units[] is set by caller; */
    /* (*win_ptr)->all_win_handles[] is set by caller; */
    (*win_ptr)->current_lock_type   = MPID_LOCK_NONE;
    (*win_ptr)->shared_lock_ref_cnt = 0;
    (*win_ptr)->lock_queue          = NULL;
    (*win_ptr)->pt_rma_puts_accs    = NULL;
    (*win_ptr)->my_pt_rma_puts_accs = 0;

    /* Initialize the passive target lock state */
    MPIU_CHKPMEM_MALLOC((*win_ptr)->targets, struct MPIDI_Win_target_state_s *,
                        sizeof(struct MPIDI_Win_target_state_s)*MPIR_Comm_size(win_comm_ptr),
                        mpi_errno, "RMA target states array");

    for (i = 0; i < MPIR_Comm_size(win_comm_ptr); i++) {
        (*win_ptr)->targets[i].rma_ops_list = NULL;
        (*win_ptr)->targets[i].remote_lock_state = MPIDI_CH3_WIN_LOCK_NONE;
    }

    /* Initialize the info flags */
    (*win_ptr)->info_args.no_locks            = 0;
    (*win_ptr)->info_args.accumulate_ordering = MPIDI_ACC_ORDER_RAR | MPIDI_ACC_ORDER_RAW |
                                                MPIDI_ACC_ORDER_WAR | MPIDI_ACC_ORDER_WAW;
    (*win_ptr)->info_args.accumulate_ops      = MPIDI_ACC_OPS_SAME_OP_NO_OP;
    (*win_ptr)->info_args.same_size           = 0;
    (*win_ptr)->info_args.alloc_shared_noncontig = 0;

    MPID_WIN_FTABLE_SET_DEFAULTS(win_ptr);

fn_exit:
    MPIDI_FUNC_EXIT(MPID_STATE_WIN_INIT);
    return mpi_errno;
fn_fail:
    MPIU_CHKPMEM_REAP();
    goto fn_exit;
}


#undef FUNCNAME
#define FUNCNAME MPID_Win_set_info
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_set_info(MPID_Win *win, MPID_Info *info)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_SET_INFO);

    MPIDI_FUNC_ENTER(MPID_STATE_MPID_WIN_SET_INFO);

    mpi_errno = win->RMAFns.Win_set_info(win, info);
    if (mpi_errno != MPI_SUCCESS) { MPIU_ERR_POP(mpi_errno); }

 fn_exit:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_SET_INFO);
    return mpi_errno;
 fn_fail:
    goto fn_exit;
}


#undef FUNCNAME
#define FUNCNAME MPID_Win_get_info
#undef FCNAME
#define FCNAME MPIDI_QUOTE(FUNCNAME)
int MPID_Win_get_info(MPID_Win *win, MPID_Info **info_used)
{
    int mpi_errno = MPI_SUCCESS;
    MPIDI_STATE_DECL(MPID_STATE_MPID_WIN_GET_INFO);

    MPIDI_FUNC_ENTER(MPID_STATE_MPID_WIN_GET_INFO);

    mpi_errno = win->RMAFns.Win_get_info(win, info_used);
    if (mpi_errno != MPI_SUCCESS) { MPIU_ERR_POP(mpi_errno); }

 fn_exit:
    MPIDI_FUNC_EXIT(MPID_STATE_MPID_WIN_GET_INFO);
    return mpi_errno;
 fn_fail:
    goto fn_exit;
}
