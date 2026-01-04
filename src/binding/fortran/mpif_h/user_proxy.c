/*
 * Copyright (C) by Argonne National Laboratory
 *     See COPYRIGHT in top-level directory
 */

#include "mpi_fortimpl.h"
#include <assert.h>

#include "uthash.h"

/* ---- attr -----*/
struct F77_attr_state {
    int keyval;
    MPI_Copy_function copy_fn;
    MPI_Delete_function delete_fn;
    void *extra_state;
    UT_hash_handle hh;
};

static struct F77_attr_state *keyval_hash = NULL;

static int F77_attr_copy_proxy(MPI_Comm comm, int keyval, void *context,
                               void *value_in, void *value_out, int *flag)
{
    struct F77_attr_state *p = (struct F77_attr_state *) context;

    MPI_Fint ierr = 0;
    MPI_Fint fhandle = MPI_Comm_toint(comm);
    MPI_Fint fkeyval = (MPI_Fint) keyval;
    MPI_Aint fvalue = (MPI_Aint) value_in;
    MPI_Aint *fextra = (MPI_Aint *) p->extra_state;
    MPI_Aint fnew = 0;
    MPI_Fint fflag = 0;

    p->copy_fn(&fhandle, &fkeyval, fextra, &fvalue, &fnew, &fflag, &ierr);
    *flag = MPII_FROM_FLOG(fflag);
    *value_out = (void *) fnew;
}

static int F77_attr_delete_proxy(MPI_Comm comm, int keyval, void *value, void *context)
{
    struct F77_attr_state *p = (struct F77_attr_state *) context;

    MPI_Fint ierr = 0;
    MPI_Fint fhandle = MPI_Comm_toint(comm);
    MPI_Fint fkeyval = (MPI_Fint) keyval;
    MPI_Aint fvalue = (MPI_Aint) value;
    MPI_Aint *fextra = (MPI_Aint *) p->extra_state;

    p->delete_fn(&fhandle, &fkeyval, &fvalue, fextra, &ierr);
    return (int) ierr;
}

static int F77_Comm_attr_copy_proxy(MPI_Comm comm, int keyval, void *context,
                                    void *value_in, void *value_out, int *flag)
{
    int handle = MPI_Comm_toint(comm);
    return F77_attr_copy_proxy(handle, keyval, context, value_in, value_out, flag);
}

static int F77_Comm_attr_delete_proxy(MPI_Comm comm, int keyval, void *value, void *context)
{
    int handle = MPI_Comm_toint(comm);
    return F77_attr_delete_proxy(handle, keyval, context, value);
}

static int F77_Win_attr_copy_proxy(MPI_Win win, int keyval, void *context,
                                   void *value_in, void *value_out, int *flag)
{
    int handle = MPI_Win_toint(win);
    return F77_attr_copy_proxy(handle, keyval, context, value_in, value_out, flag);
}

static int F77_Win_attr_delete_proxy(MPI_Win win, int keyval, void *value, void *context)
{
    int handle = MPI_Win_toint(win);
    return F77_attr_delete_proxy(handle, keyval, context, value);
}

static int F77_Type_attr_copy_proxy(MPI_Datatype datatype, int keyval, void *context,
                                    void *value_in, void *value_out, int *flag)
{
    int handle = MPI_Type_toint(datatype);
    return F77_attr_copy_proxy(handle, keyval, context, value_in, value_out, flag);
}

static int F77_Type_attr_delete_proxy(MPI_Datatype datatype, int keyval, void *value, void *context)
{
    int handle = MPI_Type_toint(datatype);
    return F77_attr_delete_proxy(handle, keyval, context, value);
}

int MPII_Keyval_create(MPI_Copy_function copy_fn, MPI_Delete_function delete_fn, int *keyval_out,
                       void *extra_state, enum F77_handle_type type)
{
    struct F77_attr_state *p = MPL_malloc(sizeof(struct F77_attr_state), MPL_MEM_OTHER);
    p->copy_fn = copy_fn;
    p->delete_fn = delete_fn;
    p->extra_state = extra_state;

    int mpi_errno;
    switch (type) {
        case F77_COMM:
            mpi_errno = MPI_Comm_create_keyval(F77_Comm_attr_copy_proxy, F77_Comm_attr_delete_proxy,
                                               keyval_out, p);
            break;
        case F77_WIN:
            mpi_errno = MPI_Win_create_keyval(F77_Win_attr_copy_proxy, F77_Win_attr_delete_proxy,
                                              keyval_out, p);
            break;
        case F77_DATATYPE:
            mpi_errno = MPI_Type_create_keyval(F77_Type_attr_copy_proxy, F77_Type_attr_delete_proxy,
                                               keyval_out, p);
            break;
        default:
            assert(0);
    }

    if (mpi_errno == MPI_SUCCESS) {
        p->keyval = *keyval_out;
        HASH_ADD_INT(keyval_hash, keyval, p, MPL_MEM_OTHER);
    } else {
        MPL_free(p);
    }

    return mpi_errno;
}

int MPII_Keyval_free(int keyval)
{
    struct F77_attr_state *p;
    HASH_FIND_INT(keyval_hash, &keyval, p);
    if (p) {
        MPL_free(p);
    } else {
        /* TODO: error handling */
        assert(p);
    }
    return MPI_Keyval_free(keyval);
}

/* When keyval is builtin (e.g. MPI_TAG_UB), C returns a pointer to an internal int instead of
 * directly returning the int value casted as a (void *), which is the behavior in fortran when
 * the keyval is not builtin.
 * It's a very stupid divergence of behavior, but now we have to deal with it in the binding.
 */
bool MPII_Keyval_is_builtin(int keyval)
{
    /* There is no API for checking builtin handles, so we are *hacking* it.
     * Yeah, we know this fragile and ugly */
    if (keyval > 0 && keyval < 4096) {
        /* in MPI ABI, they are between 5xx (comm) or 6xx (win) */
        return true;
    } else if ((keyval & 0xffff0000) == 0x64400000) {
        return true;
    } else {
        return false;
    }
}

/* ---- user op ----------------- */
struct F77_op_state {
    F77_OpFunction *opfn;
};

static void F77_op_proxy(void *invec, void *inoutvec, MPI_Count len, MPI_Datatype datatype,
                         void *extra_state)
{
    struct F77_op_state *p = extra_state;
    MPI_Fint len_i = (MPI_Fint) len;
    MPI_Fint datatype_i = MPI_Type_c2f(datatype);

    p->opfn(invec, inoutvec, &len_i, &datatype_i);
}

static void F77_op_free(void *extra_state)
{
    free(extra_state);
}

int MPII_op_create(F77_OpFunction * opfn, MPI_Fint commute, MPI_Fint * op)
{
    struct F77_op_state *p = malloc(sizeof(struct F77_op_state));
    p->opfn = opfn;

    MPI_Op op_i;
    int ret = MPIX_Op_create_x(F77_op_proxy, F77_op_free, commute, p, &op_i);
    if (ret == MPI_SUCCESS) {
        *op = MPI_Op_c2f(op_i);
    } else {
        free(p);
    }

    return ret;
}

/* ---- user errhandler  ----------------- */
struct F77_errhan_state {
    F77_ErrFunction *err_fn;
};

static void F77_comm_errhan_proxy(MPI_Comm comm, int error_code, void *extra_state)
{
    struct F77_errhan_state *p = extra_state;
    MPI_Fint comm_i = MPI_Comm_c2f(comm);
    MPI_Fint error_code_i = error_code;

    p->err_fn(&comm_i, &error_code_i);
}

static void F77_win_errhan_proxy(MPI_Win win, int error_code, void *extra_state)
{
    struct F77_errhan_state *p = extra_state;
    MPI_Fint win_i = MPI_Win_c2f(win);
    MPI_Fint error_code_i = error_code;

    p->err_fn(&win_i, &error_code_i);
}

static void F77_file_errhan_proxy(MPI_File file, int error_code, void *extra_state)
{
    struct F77_errhan_state *p = extra_state;
    MPI_Fint file_i = MPI_File_c2f(file);
    MPI_Fint error_code_i = error_code;

    p->err_fn(&file_i, &error_code_i);
}

static void F77_session_errhan_proxy(MPI_Session session, int error_code, void *extra_state)
{
    struct F77_errhan_state *p = extra_state;
    MPI_Fint session_i = MPI_Session_c2f(session);
    MPI_Fint error_code_i = error_code;

    p->err_fn(&session_i, &error_code_i);
}

static void F77_errhan_free(void *extra_state)
{
    free(extra_state);
}

int MPII_errhan_create(F77_ErrFunction * err_fn, MPI_Fint * errhandler, enum F77_handle_type type)
{
    struct F77_errhan_state *p = malloc(sizeof(struct F77_errhan_state));
    p->err_fn = err_fn;

    MPI_Errhandler errhandler_i;
    int ret = MPI_SUCCESS;
    switch (type) {
        case F77_COMM:
            ret = MPIX_Comm_create_errhandler_x(F77_comm_errhan_proxy, F77_errhan_free,
                                                p, &errhandler_i);
            break;
        case F77_WIN:
            ret = MPIX_Win_create_errhandler_x(F77_win_errhan_proxy, F77_errhan_free,
                                               p, &errhandler_i);
            break;
        case F77_FILE:
            ret = MPIX_File_create_errhandler_x(F77_file_errhan_proxy, F77_errhan_free,
                                                p, &errhandler_i);
            break;
        case F77_SESSION:
            ret = MPIX_Session_create_errhandler_x(F77_session_errhan_proxy, F77_errhan_free,
                                                   p, &errhandler_i);
            break;
        default:
            assert(0);
    }
    if (ret == MPI_SUCCESS) {
        *errhandler = MPI_Errhandler_c2f(errhandler_i);
    } else {
        free(p);
    }

    return ret;
}

/* For use by MPI_F08 */
int MPII_Comm_create_errhandler(F77_ErrFunction * err_fn, MPI_Fint * errhandler)
{
    struct F77_errhan_state *p = malloc(sizeof(struct F77_errhan_state));
    p->err_fn = err_fn;

    MPI_Errhandler errhandler_i;
    int ret = MPI_SUCCESS;

    ret = MPIX_Comm_create_errhandler_x(F77_comm_errhan_proxy, F77_errhan_free, p, &errhandler_i);
    if (ret == MPI_SUCCESS) {
        *errhandler = MPI_Errhandler_c2f(errhandler_i);
    } else {
        free(p);
    }

    return ret;
}

int MPII_File_create_errhandler(F77_ErrFunction * err_fn, MPI_Fint * errhandler)
{
    struct F77_errhan_state *p = malloc(sizeof(struct F77_errhan_state));
    p->err_fn = err_fn;

    MPI_Errhandler errhandler_i;
    int ret = MPI_SUCCESS;

    ret = MPIX_File_create_errhandler_x(F77_file_errhan_proxy, F77_errhan_free, p, &errhandler_i);
    if (ret == MPI_SUCCESS) {
        *errhandler = MPI_Errhandler_c2f(errhandler_i);
    } else {
        free(p);
    }

    return ret;
}

int MPII_Win_create_errhandler(F77_ErrFunction * err_fn, MPI_Fint * errhandler)
{
    struct F77_errhan_state *p = malloc(sizeof(struct F77_errhan_state));
    p->err_fn = err_fn;

    MPI_Errhandler errhandler_i;
    int ret = MPI_SUCCESS;

    ret = MPIX_Win_create_errhandler_x(F77_win_errhan_proxy, F77_errhan_free, p, &errhandler_i);
    if (ret == MPI_SUCCESS) {
        *errhandler = MPI_Errhandler_c2f(errhandler_i);
    } else {
        free(p);
    }

    return ret;
}

int MPII_Session_create_errhandler(F77_ErrFunction * err_fn, MPI_Fint * errhandler)
{
    struct F77_errhan_state *p = malloc(sizeof(struct F77_errhan_state));
    p->err_fn = err_fn;

    MPI_Errhandler errhandler_i;
    int ret = MPI_SUCCESS;

    ret = MPIX_Session_create_errhandler_x(F77_session_errhan_proxy, F77_errhan_free,
                                           p, &errhandler_i);
    if (ret == MPI_SUCCESS) {
        *errhandler = MPI_Errhandler_c2f(errhandler_i);
    } else {
        free(p);
    }

    return ret;
}
