#ifndef ETALE_BITSET_H
#define ETALE_BITSET_H

#include <stddef.h>
#include <stdint.h>

/*
 * Dense bitset indexed by the integers 0, ..., bit_count - 1.
 *
 * Storage is ceil(bit_count / 64) uint64_t words.  All bits are initially 0.
 */
typedef struct
{
    uint64_t *words;
    size_t bit_count;
    size_t word_count;
} bitset_t;

/*
 * Fast unchecked membership test for hot loops.  The caller must guarantee
 * index < set->bit_count.
 */
static inline int bitset_contains_unchecked(const bitset_t *set, size_t index)
{
    return (set->words[index >> 6] &
            (UINT64_C(1) << (index & 63u))) != 0;
}

/* Allocate a zero-initialized bitset. Return 1 on success, 0 on failure. */
int bitset_init(bitset_t *set, size_t bit_count);

/* Set or clear one bit. Return 1 on success, 0 if the index is out of range. */
int bitset_add(bitset_t *set, size_t index);
int bitset_remove(bitset_t *set, size_t index);

/* Return nonzero precisely when index is in range and its bit is set. */
int bitset_contains(const bitset_t *set, size_t index);

/* Set every bit to zero. */
void bitset_clear(bitset_t *set);

/* Release all memory owned by the bitset. */
void bitset_destroy(bitset_t *set);

#endif /* ETALE_BITSET_H */
