/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "hydra_utils.h"

HYD_Status HYDU_Create_process(char **client_arg, int *in, int *out, int *err, int *pid)
{
    int inpipe[2], outpipe[2], errpipe[2], tpid;
    HYD_Status status = HYD_SUCCESS;

    HYDU_FUNC_ENTER();

    if (in != NULL) {
        if (pipe(inpipe) < 0) {
            HYDU_Error_printf("pipe error (errno: %d)\n", errno);
            status = HYD_SOCK_ERROR;
            goto fn_fail;
        }
    }

    if (pipe(outpipe) < 0) {
        HYDU_Error_printf("pipe error (errno: %d)\n", errno);
        status = HYD_SOCK_ERROR;
        goto fn_fail;
    }

    if (pipe(errpipe) < 0) {
        HYDU_Error_printf("pipe error (errno: %d)\n", errno);
        status = HYD_SOCK_ERROR;
        goto fn_fail;
    }

    /* Fork off the process */
    tpid = fork();
    if (tpid == 0) {    /* Child process */
        close(outpipe[0]);
        close(1);
        if (dup2(outpipe[1], 1) < 0) {
            HYDU_Error_printf("dup2 error (errno: %d)\n", errno);
            status = HYD_SOCK_ERROR;
            goto fn_fail;
        }

        close(errpipe[0]);
        close(2);
        if (dup2(errpipe[1], 2) < 0) {
            HYDU_Error_printf("dup2 error (errno: %d)\n", errno);
            status = HYD_SOCK_ERROR;
            goto fn_fail;
        }

        close(inpipe[1]);
        close(0);
        if (in != NULL) {
            if (dup2(inpipe[0], 0) < 0) {
                HYDU_Error_printf("dup2 error (errno: %d)\n", errno);
                status = HYD_SOCK_ERROR;
                goto fn_fail;
            }
        }

        if (execvp(client_arg[0], client_arg) < 0) {
            HYDU_Error_printf("execvp error\n");
            status = HYD_INTERNAL_ERROR;
            goto fn_fail;
        }
    }
    else {      /* Parent process */
        close(outpipe[1]);
        close(errpipe[1]);
        if (in)
            *in = inpipe[1];
        if (out)
            *out = outpipe[0];
        if (err)
            *err = errpipe[0];
    }

    if (pid)
        *pid = tpid;

  fn_exit:
    HYDU_FUNC_EXIT();
    return status;

  fn_fail:
    goto fn_exit;
}
