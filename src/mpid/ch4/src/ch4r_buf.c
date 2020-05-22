/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpidimpl.h"
#include "ch4r_buf.h"

MPIDIU_buf_pool_t *MPIDIU_create_buf_pool(int num, int size)
{
    MPIDIU_buf_pool_t *buf_pool;

    buf_pool = MPIDIU_create_buf_pool_internal(num, size, NULL);

    return buf_pool;
}

void MPIDIU_destroy_buf_pool(MPIDIU_buf_pool_t * pool)
{
    int ret;


    if (pool->next)
        MPIDIU_destroy_buf_pool(pool->next);

    MPID_Thread_mutex_destroy(&pool->lock, &ret);
    MPL_free(pool->memory_region);
    MPL_free(pool);

}
