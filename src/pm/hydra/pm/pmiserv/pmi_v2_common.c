/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*
 *  (C) 2008 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include "hydra.h"
#include "hydra_utils.h"
#include "bsci.h"
#include "pmiserv_pmi.h"
#include "pmi_v2_common.h"

HYD_status HYD_pmcd_pmi_v2_queue_req(int fd, int pid, int pgid, enum type type, char *args[],
                                     struct HYD_pmcd_pmi_v2_reqs **pending_reqs)
{
    struct HYD_pmcd_pmi_v2_reqs *req, *tmp;
    HYD_status status = HYD_SUCCESS;

    HYDU_MALLOC(req, struct HYD_pmcd_pmi_v2_reqs *, sizeof(struct HYD_pmcd_pmi_v2_reqs),
                status);
    req->fd = fd;
    req->pid = pid;
    req->pgid = pgid;
    req->type = type;
    req->next = NULL;

    status = HYDU_strdup_list(args, &req->args);
    HYDU_ERR_POP(status, "unable to dup args\n");

    if (*pending_reqs == NULL)
        *pending_reqs = req;
    else {
        for (tmp = *pending_reqs; tmp->next; tmp = tmp->next);
        tmp->next = req;
    }

  fn_exit:
    return status;

  fn_fail:
    goto fn_exit;
}

void HYD_pmcd_pmi_v2_print_req_list(struct HYD_pmcd_pmi_v2_reqs *pending_reqs)
{
    struct HYD_pmcd_pmi_v2_reqs *req;

    if (pending_reqs)
        HYDU_dump_noprefix(stdout, "(");
    for (req = pending_reqs; req; req = req->next)
        HYDU_dump_noprefix(stdout, "%s ",
                           (req->type == NODE_ATTR_GET) ? "NODE_ATTR_GET" : "KVS_GET");
    if (pending_reqs)
        HYDU_dump_noprefix(stdout, ")\n");
}
