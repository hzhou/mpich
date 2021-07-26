#ifndef MPL_HASH_H_INCLUDED
#define MPL_HASH_H_INCLUDED

#include <string.h>

typedef void *MPL_hash_t;

MPL_hash_t MPL_hash_new();
bool MPL_hash_has(MPL_hash_t hv, char *key, int keylen);
void MPL_hash_set(MPL_hash_t hv, char *key, int keylen, char *value, int valuelen);
void MPL_hash_get(MPL_hash_t hv, char *key, int keylen, const char **value_out, int *valuelen_out);
void MPL_hash_free(MPL_hash_t hv);

/* Following routines are provided so we can use the internal strpool
 * for string management. It converts string into integers.
 */
int MPL_hash_get_keynum(MPL_hash_t hv, char *key, int keylen);
void MPL_hash_get_key(MPL_hash_t hv, int keynum, char **key_out, int *keylen_out);

/* Convenience wrappers when values are C string (null terminated) */
static inline void MPL_hash_set_str(MPL_hash_t hv, char *key, char *value)
{
    MPL_hash_set(hv, key, strlen(key), value, strlen(value));
}

static inline const char *MPL_hash_get_str(MPL_hash_t hv, char *key)
{
    const char *value;
    MPL_hash_get(hv, key, strlen(key), &value, NULL);
    return value;
}

#endif /* MPL_HASH_H_INCLUDED */
