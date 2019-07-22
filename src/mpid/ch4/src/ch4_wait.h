/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/*
 *  (C) 2019 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 *
 *  Portions of this code were written by Intel Corporation.
 *  Copyright (C) 2011-2016 Intel Corporation.  Intel provides this material
 *  to Argonne National Laboratory subject to Software Grant and Corporate
 *  Contributor License Agreement dated February 8, 2012.
 */
#ifndef CH4_WAIT_H_INCLUDED
#define CH4_WAIT_H_INCLUDED

/* MPID_{Wait,Waitall,Waitany,Waitsome,Test,Testall,Testany,Testsome}
 *   Default calls MPIR_Xxx_impl(...)
 *   Depend on configure flag, they may be replaced with MDTA versions
 *   or PER VCI version.
 */
#ifdef MPICH_THREAD_USE_MDTA

MPL_STATIC_INLINE_PREFIX int MPID_Waitall(int count, MPIR_Request * request_ptrs[],
                                          MPI_Status array_of_statuses[], int request_properties)
{
    const int single_threaded = !MPIR_ThreadInfo.isThreaded;
    int mpi_errno = MPI_SUCCESS;
    int i;
    MPIR_Thread_sync_t *sync = NULL;

    if (likely(single_threaded))
        return MPIR_Waitall_impl(count, request_ptrs, array_of_statuses, request_properties);

    MPIR_Thread_sync_alloc(&sync, count);

    /* Fix up number of pending requests and attach the sync. */
    for (i = 0; i < count; i++) {
        if (request_ptrs[i] == NULL || MPIR_Request_is_complete(request_ptrs[i])) {
            MPIR_Thread_sync_signal(sync, 0);
        } else {
            MPIR_Request_attach_sync(request_ptrs[i], sync);
        }
    }

    /* Wait on the synchronization object. */
    MPIR_Thread_sync_wait(sync);

    /* Either being signaled, or become a server, so we poll from now. */
    mpi_errno = MPIR_Waitall_impl(count, request_ptrs, array_of_statuses, request_properties);

    MPIR_Thread_sync_free(sync);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Wait(MPIR_Request * request_ptr, MPI_Status * status)
{
    const int single_threaded = !MPIR_ThreadInfo.isThreaded;
    int mpi_errno = MPI_SUCCESS;
    MPIR_Thread_sync_t *sync = NULL;

    if (likely(single_threaded))
        return MPIR_Wait_impl(request_ptr, status);

    if (request_ptr == NULL || MPIR_Request_is_complete(request_ptr))
        goto fn_exit;

    /* The request cannot be completed immediately, wait on a sync. */
    MPIR_Thread_sync_alloc(&sync, 1);
    MPIR_Request_attach_sync(request_ptr, sync);
    MPIR_Thread_sync_wait(sync);

    /* Either being signaled, or become a server, so we poll from now. */
    mpi_errno = MPIR_Wait_impl(request_ptr, status);

    MPIR_Thread_sync_free(sync);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

/* default */
MPL_STATIC_INLINE_PREFIX int MPID_Test(MPIR_Request * request_ptr, int *flag, MPI_Status * status)
{
    return MPIR_Test_impl(request_ptr, flag, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testall(int count, MPIR_Request * request_ptrs[],
                                          int *flag, MPI_Status array_of_statuses[],
                                          int requests_property)
{
    return MPIR_Testall_impl(count, request_ptrs, flag, array_of_statuses, requests_property);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, int *flag, MPI_Status * status)
{
    return MPIR_Testany_impl(count, request_ptrs, indx, flag, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    return MPIR_Testsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, MPI_Status * status)
{
    return MPIR_Waitany_impl(count, request_ptrs, indx, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    return MPIR_Waitsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
}

#elif MPICH_THREAD_GRANULARITY == MPICH_THREAD_GRANULARITY__VCI

MPL_STATIC_INLINE_PREFIX int MPID_Waitall(int count, MPIR_Request * request_ptrs[],
                                          MPI_Status array_of_statuses[], int request_properties)
{
    int mpi_errno = MPI_SUCCESS;
    MPID_Progress_state progress_state;
    int i;

    if (request_properties & MPIR_REQUESTS_PROPERTY__NO_NULL) {
        MPID_Progress_start(&progress_state);
        for (i = 0; i < count; i++) {
            while (!MPIR_Request_is_complete(request_ptrs[i])) {
                mpi_errno = MPID_Progress_wait_req(request_ptrs[i]);
                /* must check and handle the error, can't guard with HAVE_ERROR_CHECKING, but it's
                 * OK for the error case to be slower */
                if (unlikely(mpi_errno)) {
                    /* --BEGIN ERROR HANDLING-- */
                    MPID_Progress_end(&progress_state);
                    MPIR_ERR_POP(mpi_errno);
                    /* --END ERROR HANDLING-- */
                }
            }
        }
        MPID_Progress_end(&progress_state);
    } else {
        return MPIR_Waitall_impl(count, request_ptrs, array_of_statuses, request_properties);
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Wait(MPIR_Request * request_ptr, MPI_Status * status)
{
    int mpi_errno = MPI_SUCCESS;
    MPID_Progress_state progress_state;
    if (request_ptr == NULL)
        goto fn_exit;

    MPID_Progress_start(&progress_state);
    while (!MPIR_Request_is_complete(request_ptr)) {
        mpi_errno = MPID_Progress_wait_req(request_ptr);
        if (mpi_errno) {
            /* --BEGIN ERROR HANDLING-- */
            MPID_Progress_end(&progress_state);
            MPIR_ERR_POP(mpi_errno);
            /* --END ERROR HANDLING-- */
        }

        if (unlikely(MPIR_Request_is_anysrc_mismatched(request_ptr))) {
            mpi_errno = MPIR_Request_handle_proc_failed(request_ptr);
            goto fn_fail;
        }
    }
    MPID_Progress_end(&progress_state);

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Test(MPIR_Request * request_ptr, int *flag, MPI_Status * status)
{
    int mpi_errno = MPI_SUCCESS;

    MPIDI_vci_progress_counter++;
    if (MPIDI_vci_progress_counter > MPIDI_MAX_VCI_PROGRESS_ATTEMPTS) {
        return MPIR_Test_impl(request_ptr, flag, status);
    }

    int vci = MPIDI_REQUEST(request_ptr, vci);
    mpi_errno = MPIDI_Progress_test(MPIDI_PROGRESS_ALL, MPIDI_PROGRESS_TYPE__VCI, vci);
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);

    if (MPIR_Request_has_poll_fn(request_ptr)) {
        mpi_errno = MPIR_Grequest_poll(request_ptr, status);
        if (mpi_errno)
            MPIR_ERR_POP(mpi_errno);
    }

    if (MPIR_Request_is_complete(request_ptr))
        *flag = TRUE;
    else
        *flag = FALSE;

  fn_exit:
    return mpi_errno;

  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Testall(int count, MPIR_Request * request_ptrs[],
                                          int *flag, MPI_Status array_of_statuses[],
                                          int requests_property)
{
    return MPIR_Testall_impl(count, request_ptrs, flag, array_of_statuses, requests_property);
}

MPL_STATIC_INLINE_PREFIX int requests_has_single_vci(int count, MPIR_Request * request_ptrs[])
{
    int vci = -1;
    int i;
    for (i = 0; i < count; i++) {
        if (request_ptrs[i]) {
            int t_vci = MPIDI_REQUEST(request_ptrs[i], vci);
            if (vci >= 0) {
                if (vci != t_vci)
                    return -1;
            }
            else {
                vci = t_vci;
            }
        }
    }
    return vci;
}

MPL_STATIC_INLINE_PREFIX int MPID_Testany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, int *flag, MPI_Status * status)
{
    int vci = requests_has_single_vci(count, request_ptrs);
    if (vci < 0) {
        return MPIR_Testany_impl(count, request_ptrs, indx, flag, status);
    }

    MPIDI_vci_progress_counter++;
    if (MPIDI_vci_progress_counter > MPIDI_MAX_VCI_PROGRESS_ATTEMPTS) {
        return MPIR_Testany_impl(count, request_ptrs, indx, flag, status);
    }

    int i;
    int n_inactive = 0;
    int mpi_errno = MPI_SUCCESS;

    mpi_errno = MPIDI_Progress_test(MPIDI_PROGRESS_ALL, MPIDI_PROGRESS_TYPE__VCI, vci);
    /* --BEGIN ERROR HANDLING-- */
    if (mpi_errno)
        MPIR_ERR_POP(mpi_errno);
    /* --END ERROR HANDLING-- */

    for (i = 0; i < count; i++) {
        if ((i + 1) % MPIR_CVAR_REQUEST_POLL_FREQ == 0) {
            mpi_errno = MPIDI_Progress_test(MPIDI_PROGRESS_ALL, MPIDI_PROGRESS_TYPE__VCI, vci);
            if (mpi_errno)
                MPIR_ERR_POP(mpi_errno);
        }

        if (request_ptrs[i] != NULL && MPIR_Request_has_poll_fn(request_ptrs[i])) {
            mpi_errno = MPIR_Grequest_poll(request_ptrs[i], status);
            if (mpi_errno != MPI_SUCCESS)
                goto fn_fail;
        }
        if (!MPIR_Request_is_active(request_ptrs[i])) {
            n_inactive += 1;
        } else if (MPIR_Request_is_complete(request_ptrs[i])) {
            *flag = TRUE;
            *indx = i;
            goto fn_exit;
        }
    }

    if (n_inactive == count) {
        *flag = TRUE;
        *indx = MPI_UNDEFINED;
    }

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Testsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    int i;
    int n_inactive;
    int mpi_errno = MPI_SUCCESS;

    int vci = requests_has_single_vci(incount, request_ptrs);
    if (vci < 0) {
        return MPIR_Testsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
    }

    MPIDI_vci_progress_counter++;
    if (MPIDI_vci_progress_counter > MPIDI_MAX_VCI_PROGRESS_ATTEMPTS) {
        return MPIR_Testsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
    }

    mpi_errno = MPIDI_Progress_test(MPIDI_PROGRESS_ALL, MPIDI_PROGRESS_TYPE__VCI, vci);
    /* --BEGIN ERROR HANDLING-- */
    if (mpi_errno != MPI_SUCCESS)
        MPIR_ERR_POP(mpi_errno);
    /* --END ERROR HANDLING-- */

    n_inactive = 0;
    *outcount = 0;

    for (i = 0; i < incount; i++) {
        if ((i + 1) % MPIR_CVAR_REQUEST_POLL_FREQ == 0) {
            mpi_errno = MPIDI_Progress_test(MPIDI_PROGRESS_ALL, MPIDI_PROGRESS_TYPE__VCI, vci);
            if (mpi_errno)
                MPIR_ERR_POP(mpi_errno);
        }

        if (request_ptrs[i] != NULL && MPIR_Request_has_poll_fn(request_ptrs[i])) {
            mpi_errno = MPIR_Grequest_poll(request_ptrs[i], &array_of_statuses[i]);
            if (mpi_errno != MPI_SUCCESS)
                goto fn_fail;
        }
        if (!MPIR_Request_is_active(request_ptrs[i])) {
            n_inactive += 1;
        } else if (MPIR_Request_is_complete(request_ptrs[i])) {
            array_of_indices[*outcount] = i;
            *outcount += 1;
        }
    }

    if (n_inactive == incount)
        *outcount = MPI_UNDEFINED;

  fn_exit:
    return mpi_errno;
  fn_fail:
    goto fn_exit;
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, MPI_Status * status)
{
    return MPIR_Waitany_impl(count, request_ptrs, indx, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    return MPIR_Waitsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
}

#else

MPL_STATIC_INLINE_PREFIX int MPID_Test(MPIR_Request * request_ptr, int *flag, MPI_Status * status)
{
    return MPIR_Test_impl(request_ptr, flag, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testall(int count, MPIR_Request * request_ptrs[],
                                          int *flag, MPI_Status array_of_statuses[],
                                          int requests_property)
{
    return MPIR_Testall_impl(count, request_ptrs, flag, array_of_statuses, requests_property);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, int *flag, MPI_Status * status)
{
    return MPIR_Testany_impl(count, request_ptrs, indx, flag, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Testsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    return MPIR_Testsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
}

MPL_STATIC_INLINE_PREFIX int MPID_Wait(MPIR_Request * request_ptr, MPI_Status * status)
{
    return MPIR_Wait_impl(request_ptr, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitall(int count, MPIR_Request * request_ptrs[],
                                          MPI_Status array_of_statuses[], int request_properties)
{
    return MPIR_Waitall_impl(count, request_ptrs, array_of_statuses, request_properties);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitany(int count, MPIR_Request * request_ptrs[],
                                          int *indx, MPI_Status * status)
{
    return MPIR_Waitany_impl(count, request_ptrs, indx, status);
}

MPL_STATIC_INLINE_PREFIX int MPID_Waitsome(int incount, MPIR_Request * request_ptrs[],
                                           int *outcount, int array_of_indices[],
                                           MPI_Status array_of_statuses[])
{
    return MPIR_Waitsome_impl(incount, request_ptrs, outcount, array_of_indices, array_of_statuses);
}

#endif

#endif  /* CH4_WAIT_H_INCLUDED */
