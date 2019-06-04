/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/* vim: set ft=c.mpich : */
/*
 *  (C) 2018 by Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef MPIDU_BC_H_INCLUDED
#define MPIDU_BC_H_INCLUDED

#include "build_nodemap.h"

int MPIDU_bc_table_create(int rank, int size, int *nodemap, void *bc, int bc_len, int same_len,
                          int roots_only, void **bc_table, MPI_Aint ** bc_indices);
int MPIDU_bc_table_destroy(void *bc_table);
int MPIDU_bc_allgather(MPIR_Comm * comm, int *nodemap, void *bc, int bc_len, int same_len,
                       void **bc_table, MPI_Aint ** bc_indices);

#endif /* MPIDU_BC_H_INCLUDED */
