/*
*  Copyright (C) by Argonne National Laboratory.
*      See COPYRIGHT in top-level directory.
*/

/** Rationale:
 *    An open addressing hash implementation specialized for string keys and string
 *    values. All strings are stored in a managed string pool in a insertion only
 *    manner. Deletion is achieved by simply overwrite the hash item with sentinel
 *    value/NULL. Deletion and reassignement will result in expired strings still
 *    being held in memory, which may be a concern if the application requires
 *    frequent deletion and modification and memory consumption is of concern. All
 *    memory will be released upon final free.
 *
 *    Compared to uthash.h, the usage interface is cleaner and more intuitive --
 *    hash_new, hash_set, hash_get, hash_has, and hash_free. In comparison, uthash
 *    requires extra data structure, does not manage string memory, and the code is
 *    more complex. This implementation is faster than uthash as well.
 */

#include "mplconfig.h"
#include "mpl.h"
#include <stdlib.h>
#include <assert.h>

struct strpool {
    int i_str;
    int i_pool;
    int n_str;
    size_t pool_size;
    int *pn_str;
    unsigned char *pc_pool;
};

struct MPL_hash {
    int n_size;
    int n_count;
    int n_val_size;
    char *p_exist;
    int *p_key;
    struct strpool pool;
    int *p_val;
};

/* Internal functions */
static int f_strhash_lookup(struct MPL_hash *hash, unsigned char *pc_key, int keylen);
static int f_strhash_lookup_add(struct MPL_hash *hash, unsigned char *pc_key, int keylen);
static int f_addto_strpool(struct strpool *pool, char *s, int n);
static void f_strhash_resize(struct MPL_hash *hash, int n_size);
static void f_resize_strpool_n(struct strpool *pool, int n_size, int n_avg);

MPL_hash_t MPL_hash_new()
{
    struct MPL_hash *hv;
    hv = (struct MPL_hash *) calloc(1, sizeof(struct MPL_hash));
    hv->n_val_size = sizeof(int);
    return (MPL_hash_t) hv;
}

bool MPL_hash_has(MPL_hash_t hv, char *s_key, int keylen)
{
    struct MPL_hash *h = hv;
    int k = f_strhash_lookup(h, (unsigned char *) s_key, keylen);
    return h->p_exist[k];
}

void MPL_hash_set(MPL_hash_t hv, char *s_key, int keylen, char *s_value, int valuelen)
{
    struct MPL_hash *h = hv;
    int k = f_strhash_lookup_add(h, (unsigned char *) s_key, keylen);
    h->p_val[k] = f_addto_strpool(&h->pool, s_value, valuelen);
}

void MPL_hash_get(MPL_hash_t hv, char *s_key, int keylen, const char **value_out, int *valuelen_out)
{
    struct MPL_hash *h = hv;
    int k = f_strhash_lookup(h, (unsigned char *) s_key, keylen);
    if (h->p_exist[k]) {
        struct strpool *pool = &h->pool;
        int strnum = h->p_val[k];
        *value_out = (char *) ((pool->pc_pool + pool->pn_str[strnum]));
        /* note: we allow valuelen_out to be NULL (if user don't care) */
        if (valuelen_out) {
            *valuelen_out = pool->pn_str[strnum + 1] - pool->pn_str[strnum] - 1;
        }
    } else {
        *value_out = NULL;
        if (valuelen_out) {
            *valuelen_out = -1;
        }
    }
}

void MPL_hash_free(MPL_hash_t hv)
{
    struct MPL_hash *h = hv;
    free(h->p_key);
    free(h->p_exist);
    free(h->p_val);
    free(h->pool.pn_str);
    free(h->pool.pc_pool);
    free(h);
}

int MPL_hash_get_keynum(MPL_hash_t hv, char *key, int keylen)
{
    struct MPL_hash *h = hv;
    int k = f_strhash_lookup(h, (unsigned char *) key, keylen);
    if (h->p_exist[k]) {
        return h->p_key[k];
    } else {
        return -1;
    }
}

void MPL_hash_get_key(MPL_hash_t hv, int keynum, char **key_out, int *keylen_out)
{
    struct MPL_hash *h = hv;
    struct strpool *pool = &h->pool;
    if (keynum >= pool->i_str) {
        *key_out = NULL;
        *keylen_out = 0;
    } else {
        *key_out = (char *) (pool->pc_pool + pool->pn_str[keynum]);
        *keylen_out = pool->pn_str[keynum + 1] - pool->pn_str[keynum] - 1;
    }
}

int f_strhash_lookup(struct MPL_hash *hash, unsigned char *pc_key, int keylen)
{
    unsigned int tu_h;
    int k;

    if (hash->n_size == 0) {
        return 0;
    }
    /* FNV-1a hash */
    tu_h = 2166136261u;
    for (int i = 0; i < keylen; i++) {
        tu_h = tu_h ^ pc_key[i];
        tu_h = tu_h * 16777619;
    }
    k = (int) (tu_h % hash->n_size);

    struct strpool *pool = &hash->pool;
    while (1) {
        if (!hash->p_exist[k]) {
            return k;
        } else {
            int key = hash->p_key[k];
            int strlen = pool->pn_str[key + 1] - pool->pn_str[key] - 1;
            unsigned char *s = pool->pc_pool + pool->pn_str[key];
            if (strlen == keylen && memcmp(s, pc_key, keylen) == 0) {
                return k;
            }
        }

        if (k == 0) {
            k = hash->n_size - 1;
        } else {
            k--;
        }
    }
}

int f_strhash_lookup_add(struct MPL_hash *hash, unsigned char *pc_key, int keylen)
{
    int k;

    if (hash->n_size == 0 || hash->n_count + 1 >= hash->n_size ||
        (hash->n_size > 20 && hash->n_count >= hash->n_size * 85 / 100)) {
        f_strhash_resize(hash, 0);
    }
    k = f_strhash_lookup(hash, pc_key, keylen);
    if (!hash->p_exist[k]) {
        hash->p_key[k] = f_addto_strpool(&hash->pool, (char *) pc_key, keylen);
        hash->p_exist[k] = 1;
        hash->n_count++;
    }
    return k;
}

int f_addto_strpool(struct strpool *pool, char *s, int n)
{
    int n_size;
    int tn_avg_str_size;

    if (pool->i_str + 2 >= pool->n_str) {
        n_size = pool->i_str * 3 / 2 + 10;
        f_resize_strpool_n(pool, n_size, 0);
    }

    if (pool->i_pool + n >= pool->pool_size) {
        tn_avg_str_size = 6;
        if (pool->i_str > 0) {
            tn_avg_str_size = pool->i_pool / pool->i_str + 1;
        }
        pool->pool_size = tn_avg_str_size * pool->n_str + n;
        pool->pc_pool = realloc(pool->pc_pool, pool->pool_size);
    }
    memcpy(pool->pc_pool + pool->i_pool, s, n);
    pool->i_str++;
    pool->i_pool += n;
    assert(pool->i_str > 0);
    assert(pool->i_pool > 0);
    pool->pc_pool[pool->i_pool] = '\0';
    pool->i_pool += 1;
    pool->pn_str[pool->i_str] = pool->i_pool;

    return pool->i_str - 1;
}

void f_strhash_resize(struct MPL_hash *hash, int n_size)
{
    char *tp_exist = hash->p_exist;
    int *tp_key = hash->p_key;
    int tn_old_size;
    void *tp_val = hash->p_val;
    int k;

    if (n_size == 0) {
        if (hash->n_size <= 0) {
            n_size = 10;
        } else {
            n_size = hash->n_size * 5 / 3;
        }
    } else if (n_size <= hash->n_size) {
        return;
    }

    tn_old_size = hash->n_size;

    hash->n_size = n_size;
    hash->p_key = (int *) calloc(n_size, sizeof(int));
    hash->p_exist = (char *) calloc(n_size, sizeof(char));
    if (hash->n_val_size > 0) {
        hash->p_val = (void *) calloc((n_size * hash->n_val_size), sizeof(void));
    }

    struct strpool *pool = &hash->pool;
    if (tn_old_size > 0) {
        for (int i = 0; i < tn_old_size; i++) {
            if (tp_exist[i]) {

                k = f_strhash_lookup(hash, (pool->pc_pool + pool->pn_str[tp_key[i]]),
                                     (pool->pn_str[tp_key[i] + 1] - pool->pn_str[tp_key[i]] - 1));
                hash->p_exist[k] = 1;
                hash->p_key[k] = tp_key[i];
                if (hash->n_val_size > 0) {
                    memcpy(hash->p_val + k * hash->n_val_size, tp_val + i * hash->n_val_size,
                           hash->n_val_size);
                }
            }
        }
        free(tp_exist);
        free(tp_key);
        if (hash->n_val_size > 0) {
            free(tp_val);
        }
    }
}

void f_resize_strpool_n(struct strpool *pool, int n_size, int n_avg)
{
    if (n_size <= pool->n_str) {
        return;
    }
    pool->n_str = n_size;
    pool->pn_str = realloc(pool->pn_str, pool->n_str * sizeof(int));

    if (n_avg == 0) {
        n_avg = 6;
        if (pool->i_str > 0) {
            n_avg = pool->i_pool / pool->i_str + 1;
        }
    }
    if (n_avg * pool->n_str > pool->pool_size) {
        pool->pool_size = n_avg * pool->n_str;
        pool->pc_pool = realloc(pool->pc_pool, pool->pool_size);
        pool->pc_pool[pool->pool_size - 1] = 0;
    }

    if (pool->i_str == 0) {
        pool->pn_str[0] = 0;
    }
}
