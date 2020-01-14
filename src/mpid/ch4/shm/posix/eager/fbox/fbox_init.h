/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2006 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2017 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef POSIX_EAGER_FBOX_INIT_H_INCLUDED
#define POSIX_EAGER_FBOX_INIT_H_INCLUDED

#include "fbox_types.h"
#include "fbox_impl.h"

/*
=== BEGIN_MPI_T_CVAR_INFO_BLOCK ===

cvars:
    - name        : MPIR_CVAR_CH4_POSIX_EAGER_FBOX_POLL_CACHE_SIZE
      category    : CH4
      type        : int
      default     : 3
      class       : device
      verbosity   : MPI_T_VERBOSITY_USER_BASIC
      scope       : MPI_T_SCOPE_ALL_EQ
      description : >-
        The size of the array to store expected receives to speed up fastbox polling.

=== END_MPI_T_CVAR_INFO_BLOCK ===
*/

#define MPIDI_POSIX_MAILBOX_INDEX(sender, receiver) ((MPIDI_POSIX_global.num_local) * (sender) + (receiver))

extern MPIDI_POSIX_eager_fbox_control_t MPIDI_POSIX_eager_fbox_control_global;

MPL_STATIC_INLINE_PREFIX int MPIDI_POSIX_eager_init(int rank, int size)
{
    int mpi_errno = MPI_SUCCESS;
    int i;
    MPIDI_POSIX_fastbox_t *fastboxes_p = NULL;


#ifdef MPL_USE_DBG_LOGGING
    MPIDI_CH4_SHM_POSIX_FBOX_GENERAL = MPL_dbg_class_alloc("SHM_POSIX_FBOX", "shm_posix_fbox");
#endif /* MPL_USE_DBG_LOGGING */

    MPIR_CHKPMEM_DECL(3);

    MPIR_CHKPMEM_MALLOC(MPIDI_POSIX_eager_fbox_control_global.first_poll_local_ranks,
                        int16_t *,
                        sizeof(*MPIDI_POSIX_eager_fbox_control_global.first_poll_local_ranks) *
                        (MPIR_CVAR_CH4_POSIX_EAGER_FBOX_POLL_CACHE_SIZE + 1), mpi_errno,
                        "first_poll_local_ranks", MPL_MEM_SHM);

    /* -1 means we aren't looking for anything in particular. */
    for (i = 0; i < MPIR_CVAR_CH4_POSIX_EAGER_FBOX_POLL_CACHE_SIZE; i++) {
        MPIDI_POSIX_eager_fbox_control_global.first_poll_local_ranks[i] = -1;
    }

    /* The final entry in the cache array is for tracking the fallback place to start looking for
     * messages if all other entries in the cache are empty. */
    MPIDI_POSIX_eager_fbox_control_global.first_poll_local_ranks[i] = 0;

    /* Create region with one fastbox for every pair of local processes. */
    size_t len =
        MPIDI_POSIX_global.num_local * MPIDI_POSIX_global.num_local * sizeof(MPIDI_POSIX_fastbox_t);

    /* Actually allocate the segment and assign regions to the pointers */
    mpi_errno = MPIDU_Init_shm_alloc(len, &MPIDI_POSIX_eager_fbox_control_global.shm_ptr);
    MPIR_ERR_CHECK(mpi_errno);

    fastboxes_p = (MPIDI_POSIX_fastbox_t *) MPIDI_POSIX_eager_fbox_control_global.shm_ptr;

    /* Allocate table of pointers to fastboxes */
    MPIR_CHKPMEM_MALLOC(MPIDI_POSIX_eager_fbox_control_global.mailboxes.in,
                        MPIDI_POSIX_fastbox_t **,
                        MPIDI_POSIX_global.num_local * sizeof(MPIDI_POSIX_fastbox_t *), mpi_errno,
                        "fastboxes", MPL_MEM_SHM);
    MPIR_CHKPMEM_MALLOC(MPIDI_POSIX_eager_fbox_control_global.mailboxes.out,
                        MPIDI_POSIX_fastbox_t **,
                        MPIDI_POSIX_global.num_local * sizeof(MPIDI_POSIX_fastbox_t *), mpi_errno,
                        "fastboxes", MPL_MEM_SHM);

    /* Fill in fbox tables */
    for (i = 0; i < MPIDI_POSIX_global.num_local; i++) {
        MPIDI_POSIX_eager_fbox_control_global.mailboxes.in[i] =
            &fastboxes_p[MPIDI_POSIX_MAILBOX_INDEX(i, MPIDI_POSIX_global.my_local_rank)];
        MPIDI_POSIX_eager_fbox_control_global.mailboxes.out[i] =
            &fastboxes_p[MPIDI_POSIX_MAILBOX_INDEX(MPIDI_POSIX_global.my_local_rank, i)];

        memset(MPIDI_POSIX_eager_fbox_control_global.mailboxes.in[i], 0,
               sizeof(MPIDI_POSIX_fastbox_t));
    }

    mpi_errno = MPIDU_Init_shm_barrier();
    MPIR_ERR_CHECK(mpi_errno);

    MPIR_CHKPMEM_COMMIT();

  fn_exit:
    return mpi_errno;
  fn_fail:
    /* --BEGIN ERROR HANDLING-- */
    MPIR_CHKPMEM_REAP();
    goto fn_exit;
    /* --END ERROR HANDLING-- */
}

MPL_STATIC_INLINE_PREFIX int MPIDI_POSIX_eager_finalize()
{
    int mpi_errno = MPI_SUCCESS;


    MPL_free(MPIDI_POSIX_eager_fbox_control_global.mailboxes.in);
    MPL_free(MPIDI_POSIX_eager_fbox_control_global.mailboxes.out);
    MPL_free(MPIDI_POSIX_eager_fbox_control_global.first_poll_local_ranks);

    mpi_errno = MPIDU_Init_shm_free(MPIDI_POSIX_eager_fbox_control_global.shm_ptr);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

#endif /* POSIX_EAGER_FBOX_INIT_H_INCLUDED */
