#include "etale_bitset.h"

#include <stdlib.h>
#include <string.h>

#define BITSET_WORD_BITS ((size_t)64)

int bitset_init(bitset_t *set, size_t bit_count)
{
    size_t word_count;

    if (set == NULL)
        return 0;

    set->words = NULL;
    set->bit_count = 0;
    set->word_count = 0;

    if (bit_count == 0)
        return 1;

    word_count = (bit_count + BITSET_WORD_BITS - 1) / BITSET_WORD_BITS;
    set->words = calloc(word_count, sizeof(*set->words));
    if (set->words == NULL)
        return 0;

    set->bit_count = bit_count;
    set->word_count = word_count;
    return 1;
}

int bitset_add(bitset_t *set, size_t index)
{
    if (set == NULL || index >= set->bit_count)
        return 0;

    set->words[index / BITSET_WORD_BITS] |=
        UINT64_C(1) << (index % BITSET_WORD_BITS);
    return 1;
}

int bitset_remove(bitset_t *set, size_t index)
{
    if (set == NULL || index >= set->bit_count)
        return 0;

    set->words[index / BITSET_WORD_BITS] &=
        ~(UINT64_C(1) << (index % BITSET_WORD_BITS));
    return 1;
}

int bitset_contains(const bitset_t *set, size_t index)
{
    if (set == NULL || index >= set->bit_count)
        return 0;

    return bitset_contains_unchecked(set, index);
}

void bitset_clear(bitset_t *set)
{
    if (set == NULL || set->words == NULL)
        return;

    memset(set->words, 0, set->word_count * sizeof(*set->words));
}

void bitset_destroy(bitset_t *set)
{
    if (set == NULL)
        return;

    free(set->words);
    set->words = NULL;
    set->bit_count = 0;
    set->word_count = 0;
}
