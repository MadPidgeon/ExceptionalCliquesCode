/*
 * Compute the clique number of the graph attached to the finite etale algebra
 *
 *     A = F_{p^{d_1}} x ... x F_{p^{d_r}}
 *
 * where x and y are adjacent when N_{A/F_p}(x-y) is +1 or -1.
 *
 * Every clique containing an edge is normalized to contain 0 and 1.  The
 * remaining normalized vertices are the exceptional units
 *
 *     E = { x in A^*: N(x)=+-1 and N(1-x)=+-1 }.
 *
 * Norm-+-1 units are indexed directly by a mixed-radix coordinate system
 * of size |N^{-1}({+-1})|.  Primitive generators in all field factors are
 * normalized to have the same norm zeta in F_p^*, so H is described by the
 * single congruence sum_i e_i=0 mod (p-1)/2.  A smallest-degree factor is used
 * as the dependent mixed-radix pivot.  Equal-degree factors share the same
 * PARI finite-field model, so their permutations are literal exponent swaps
 * followed by re-encoding.
 *
 * The computation has four main phases:
 *   1. construct finite-field models and index H directly;
 *   2. enumerate exceptional units and classify symmetry orbits;
 *   3. solve one canonical maximum-clique problem per orbit;
 *   4. verify that every orbit completed successfully before reporting the
 *      final clique number.
 *
 */

#define _POSIX_C_SOURCE 200809L

#include <pari/pari.h>
#include <cliquer.h>
#include "etale_bitset.h"

#include <dirent.h>
#include <errno.h>
#include <inttypes.h>
#include <limits.h>
#include <signal.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>

#define PROGRAM_NAME "EtaleUpperBound"
#define DEFAULT_V4_THRESHOLD 0
/* Decimal digits used for every printed elapsed time.  Override at build
 * time, e.g. `make TIME_PRECISION=6`, without changing the source. */
#ifndef TIME_PRECISION
#define TIME_PRECISION 9
#endif
#if TIME_PRECISION < 0 || TIME_PRECISION > 12
#error "TIME_PRECISION must be between 0 and 12"
#endif
#define UNSET_CODE UINT64_MAX
#define UNSET_EXP UINT64_MAX

#define USAGE \
"Usage: " PROGRAM_NAME " p d1 [d2 ...] [options]\n" \
"Options:\n" \
"  -cpus n           number of worker processes (default 1)\n" \
"  -pari MB          PARI stack size in MB (default 100)\n" \
"  -v4-threshold n   use V4 when the canonical local graph has >= n vertices\n" \
"                    (default 0: disabled)\n"

typedef uint64_t h_code_t;
typedef uint64_t factor_exp_t;

#if defined(__SIZEOF_INT128__)
__extension__ typedef unsigned __int128 wide_uint_t;
#endif

typedef struct
{
    int degree;
    uint64_t q;
    uint64_t order;

    /* Mixed-radix data for the direct H=N^{-1}({+-1}) indexing. */
    uint64_t h_radix;
    uint64_t h_stride;

    int model_index;
} factor_t;

typedef struct
{
    int degree;
    uint64_t q;
    uint64_t order;
    int prime_field;

    GEN field_generator;
    GEN primitive_generator;
    GEN one;
    GEN order_data;

    factor_exp_t *one_minus_exp;
} field_model_t;

typedef struct
{
    int left;
    int right;
} swap_generator_t;

typedef struct
{
    int size;
    int *elements;
} orbit_t;

typedef struct
{
    int size;
    int elements[4];
} v4_orbit_t;

typedef struct
{
    int candidate_count;
    int clique_size;
    int skipped;
    double elapsed;
} v4_iteration_stats_t;

typedef struct
{
    int used_v4;
    int retained_suborbits;
    v4_iteration_stats_t *iterations;
    int iteration_count;
} local_stats_t;

/* Data returned by a completed local orbit computation. */
typedef struct
{
    int clique_number;
    int vertex_count;
    double elapsed;
} orbit_result_t;

/* Global read-only data used after fork(). */
static int p;
static int factor_count;
static int total_dimension;
static int exceptional_count;
static int orbit_count;
static int v4_threshold = DEFAULT_V4_THRESHOLD;
static int pivot_factor;

static factor_t *factors;
static uint64_t norm_half_order;
static h_code_t total_h_elements;
static uint64_t base_field_generator;
static int *base_field_log;
static bitset_t exceptional_codes = {0};
static factor_exp_t *exceptional_exponents;
static int *exceptional_id_by_code;
static int *one_minus_id;
static orbit_t *orbits;

static char computation_path[MAXPATHLEN];
static char local_path[MAXPATHLEN];
static char error_path[MAXPATHLEN];

/* Idempotent cleanup, also used by the parent after its final fork. */
static void free_global_data(void);

/* -------------------------------------------------------------------------
 * General utilities
 * ------------------------------------------------------------------------- */

static void checked_snprintf(char *destination, size_t destination_size,
                             const char *format, ...)
{
    va_list arguments;
    int written;

    va_start(arguments, format);
    written = vsnprintf(destination, destination_size, format, arguments);
    va_end(arguments);

    if (written < 0 || (size_t)written >= destination_size)
    {
        fprintf(stderr, "path or formatted output is too long\n");
        exit(EXIT_FAILURE);
    }
}

static void append_text(char *destination, size_t destination_size,
                        const char *format, ...)
{
    size_t used = strlen(destination);
    va_list arguments;
    int written;

    if (used >= destination_size)
    {
        fprintf(stderr, "path is too long\n");
        exit(EXIT_FAILURE);
    }

    va_start(arguments, format);
    written = vsnprintf(destination + used, destination_size - used,
                        format, arguments);
    va_end(arguments);

    if (written < 0 || (size_t)written >= destination_size - used)
    {
        fprintf(stderr, "path is too long\n");
        exit(EXIT_FAILURE);
    }
}

/* Create a directory if needed, but reject other filesystem errors. */
static void ensure_directory(const char *path)
{
    if (mkdir(path, 0700) != 0 && errno != EEXIST)
    {
        perror(path);
        exit(EXIT_FAILURE);
    }
}

static int has_suffix(const char *text, const char *suffix)
{
    size_t text_length = strlen(text);
    size_t suffix_length = strlen(suffix);

    return text_length >= suffix_length &&
           strcmp(text + text_length - suffix_length, suffix) == 0;
}

/*
 * A run reuses the directory determined by p and the factor degrees.  Remove
 * only files produced by earlier local orbit computations; unrelated files in
 * the directory are left untouched.
 */
static void clear_previous_local_results(void)
{
    DIR *directory = opendir(local_path);
    struct dirent *entry;

    if (directory == NULL)
    {
        perror("opendir local output directory");
        exit(EXIT_FAILURE);
    }

    while ((entry = readdir(directory)) != NULL)
    {
        char path[MAXPATHLEN];

        if (strncmp(entry->d_name, "local_", 6) != 0 ||
            !has_suffix(entry->d_name, ".out"))
            continue;

        checked_snprintf(path, sizeof(path), "%s/%s",
                         local_path, entry->d_name);
        if (unlink(path) != 0 && errno != ENOENT)
        {
            perror("unlink old local result");
            closedir(directory);
            exit(EXIT_FAILURE);
        }
    }

    if (closedir(directory) != 0)
    {
        perror("closedir local output directory");
        exit(EXIT_FAILURE);
    }
}

/*
 * Keep monotonic timestamps as integer seconds/nanoseconds and convert only the
 * elapsed interval to double.  Subtracting two large absolute double-valued
 * timestamps can lose precision on long-running systems.
 */
static struct timespec monotonic_now(void)
{
    struct timespec timestamp;

    if (clock_gettime(CLOCK_MONOTONIC, &timestamp) != 0)
    {
        perror("clock_gettime");
        exit(EXIT_FAILURE);
    }
    return timestamp;
}

static double elapsed_seconds(struct timespec start)
{
    struct timespec end = monotonic_now();
    time_t seconds = end.tv_sec - start.tv_sec;
    long nanoseconds = end.tv_nsec - start.tv_nsec;

    if (nanoseconds < 0)
    {
        --seconds;
        nanoseconds += 1000000000L;
    }
    return (double)seconds + (double)nanoseconds / 1000000000.0;
}

/* Flush progress data immediately so it is visible during long searches. */
static void flush_progress(FILE *stream, int outer_orbit)
{
    if (fflush(stream) != 0)
    {
        fprintf(stderr, "could not flush local progress [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }
}

/*
 * A V4 subproblem uses one fixed-width line while it is running.  When Cliquer
 * returns we overwrite that same line and seek back to EOF, so live output
 * never contains a separate "running" and "finished" line for one suborbit.
 */
#define V4_PROGRESS_LINE_WIDTH 160

static long append_v4_progress_line(FILE *stream, int outer_orbit,
                                    const char *format, ...)
{
    char text[256];
    va_list arguments;
    long position = ftell(stream);

    if (position < 0)
    {
        perror("ftell local progress");
        exit(EXIT_FAILURE);
    }
    va_start(arguments, format);
    vsnprintf(text, sizeof(text), format, arguments);
    va_end(arguments);
    fprintf(stream, "%-*.*s\n", V4_PROGRESS_LINE_WIDTH,
            V4_PROGRESS_LINE_WIDTH, text);
    flush_progress(stream, outer_orbit);
    return position;
}

static void replace_v4_progress_line(FILE *stream, long position,
                                     int outer_orbit, const char *format, ...)
{
    char text[256];
    va_list arguments;

    va_start(arguments, format);
    vsnprintf(text, sizeof(text), format, arguments);
    va_end(arguments);
    if (fseek(stream, position, SEEK_SET) != 0)
    {
        perror("fseek local progress");
        exit(EXIT_FAILURE);
    }
    fprintf(stream, "%-*.*s\n", V4_PROGRESS_LINE_WIDTH,
            V4_PROGRESS_LINE_WIDTH, text);
    flush_progress(stream, outer_orbit);
    if (fseek(stream, 0, SEEK_END) != 0)
    {
        perror("fseek local progress end");
        exit(EXIT_FAILURE);
    }
}

static void write_preparation_status(FILE *stream, const char *stage)
{
    if (fflush(stream) != 0 ||
        ftruncate(fileno(stream), 0) != 0 ||
        fseek(stream, 0, SEEK_SET) != 0)
    {
        fprintf(stderr, "could not update preparation status\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stream, "Status: preparing\nStage: %s\n", stage);
    if (fflush(stream) != 0)
    {
        fprintf(stderr, "could not flush preparation status\n");
        exit(EXIT_FAILURE);
    }
}

static int parse_int_at_least(const char *text, const char *name, int minimum)
{
    char *end = NULL;
    long value;

    errno = 0;
    value = strtol(text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' ||
        value < minimum || value > INT_MAX)
    {
        fprintf(stderr, "invalid %s: %s\n", name, text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static int parse_positive_int(const char *text, const char *name)
{
    return parse_int_at_least(text, name, 1);
}

static int parse_nonnegative_int(const char *text, const char *name)
{
    return parse_int_at_least(text, name, 0);
}

static int is_prime_int(int n)
{
    int divisor;

    if (n < 2)
        return 0;
    if ((n & 1) == 0)
        return n == 2;

    for (divisor = 3; (int64_t)divisor * divisor <= n; divisor += 2)
    {
        if (n % divisor == 0)
            return 0;
    }
    return 1;
}

static uint64_t checked_u64_power(uint64_t base, int exponent)
{
    uint64_t result = 1;
    uint64_t factor_value = base;

    while (exponent > 0)
    {
        if (exponent & 1)
        {
            if (factor_value != 0 && result > UINT64_MAX / factor_value)
            {
                fprintf(stderr, "field cardinality exceeds 64-bit range\n");
                exit(EXIT_FAILURE);
            }
            result *= factor_value;
        }

        exponent >>= 1;
        if (exponent > 0)
        {
            if (factor_value != 0 &&
                factor_value > UINT64_MAX / factor_value)
            {
                fprintf(stderr, "field cardinality exceeds 64-bit range\n");
                exit(EXIT_FAILURE);
            }
            factor_value *= factor_value;
        }
    }

    return result;
}

/* -------------------------------------------------------------------------
 * Integer arithmetic used by the H-coordinate representation
 * ------------------------------------------------------------------------- */

static uint64_t add_mod_u64(uint64_t left, uint64_t right,
                            uint64_t modulus)
{
    left %= modulus;
    right %= modulus;
    if (left >= modulus - right)
        return left - (modulus - right);
    return left + right;
}

static uint64_t subtract_mod_u64(uint64_t left, uint64_t right,
                                 uint64_t modulus)
{
    left %= modulus;
    right %= modulus;
    if (left >= right)
        return left - right;
    return modulus - (right - left);
}

static uint64_t multiply_mod_u64(uint64_t left, uint64_t right,
                                 uint64_t modulus)
{
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((wide_uint_t)(left % modulus) *
                       (wide_uint_t)(right % modulus)) %
                      (wide_uint_t)modulus);
#else
    uint64_t result = 0;

    left %= modulus;
    while (right > 0)
    {
        if (right & 1)
            result = add_mod_u64(result, left, modulus);
        right >>= 1;
        if (right > 0)
            left = add_mod_u64(left, left, modulus);
    }
    return result;
#endif
}

static uint64_t modular_power(uint64_t base, uint64_t exponent,
                              uint64_t modulus)
{
    uint64_t result = 1;

    base %= modulus;
    while (exponent > 0)
    {
        if (exponent & 1)
            result = multiply_mod_u64(result, base, modulus);
        exponent >>= 1;
        if (exponent > 0)
            base = multiply_mod_u64(base, base, modulus);
    }
    return result;
}

static uint64_t primitive_root_prime(int prime)
{
    uint64_t phi = (uint64_t)prime - 1;
    uint64_t remaining = phi;
    uint64_t prime_divisors[64];
    int divisor_count = 0;
    uint64_t divisor;
    uint64_t candidate;
    int i;

    for (divisor = 2; divisor * divisor <= remaining; ++divisor)
    {
        if (remaining % divisor != 0)
            continue;

        prime_divisors[divisor_count++] = divisor;
        while (remaining % divisor == 0)
            remaining /= divisor;
    }
    if (remaining > 1)
        prime_divisors[divisor_count++] = remaining;

    for (candidate = 2; candidate < (uint64_t)prime; ++candidate)
    {
        int primitive = 1;

        for (i = 0; i < divisor_count; ++i)
        {
            if (modular_power(candidate, phi / prime_divisors[i],
                              (uint64_t)prime) == 1)
            {
                primitive = 0;
                break;
            }
        }
        if (primitive)
            return candidate;
    }

    fprintf(stderr, "could not find a primitive root modulo %d\n", prime);
    exit(EXIT_FAILURE);
}

static uint64_t inverse_mod_u64(uint64_t value, uint64_t modulus)
{
    int64_t old_r;
    int64_t r;
    int64_t old_t = 0;
    int64_t t = 1;

    if (modulus == 1)
        return 0;

    old_r = (int64_t)modulus;
    r = (int64_t)(value % modulus);

    while (r != 0)
    {
        int64_t quotient = old_r / r;
        int64_t next_r = old_r - quotient * r;
        int64_t next_t = old_t - quotient * t;

        old_r = r;
        r = next_r;
        old_t = t;
        t = next_t;
    }

    if (old_r != 1)
    {
        fprintf(stderr, "noninvertible modular value\n");
        exit(EXIT_FAILURE);
    }

    old_t %= (int64_t)modulus;
    if (old_t < 0)
        old_t += (int64_t)modulus;
    return (uint64_t)old_t;
}


static uint64_t gcd_u64(uint64_t a, uint64_t b)
{
    while (b != 0)
    {
        uint64_t remainder = a % b;
        a = b;
        b = remainder;
    }
    return a;
}

/*
 * PARI returns an arbitrary primitive generator g of F_{p^d}^*.  Its norm is
 * another primitive element of F_p^*, say N(g)=zeta^c, where zeta is the
 * primitive root fixed by initialize_base_field_log().  Replace g by g^u with
 *
 *     c u = 1 mod (p-1)
 *
 * and gcd(u,p^d-1)=1.  The first condition gives N(g^u)=zeta; the second
 * preserves primitivity.  A suitable lift u of c^{-1} modulo p-1 always
 * exists, and the short search below finds one directly.
 */
static uint64_t normalized_generator_exponent(uint64_t norm_log,
                                               uint64_t order)
{
    uint64_t base_order = (uint64_t)p - 1;
    uint64_t u = inverse_mod_u64(norm_log, base_order);

    while (u < order)
    {
        if (gcd_u64(u, order) == 1)
            return u;

        if (u > UINT64_MAX - base_order)
            break;
        u += base_order;
    }

    fprintf(stderr, "could not normalize a primitive field generator\n");
    exit(EXIT_FAILURE);
}

static void initialize_base_field_log(void)
{
    uint64_t exponent;
    uint64_t value = 1;

    base_field_generator = primitive_root_prime(p);
    base_field_log = malloc((size_t)p * sizeof(*base_field_log));
    if (base_field_log == NULL)
    {
        fprintf(stderr, "could not allocate base-field logarithm table\n");
        exit(EXIT_FAILURE);
    }

    for (exponent = 0; exponent < (uint64_t)p; ++exponent)
        base_field_log[exponent] = -1;

    for (exponent = 0; exponent < (uint64_t)p - 1; ++exponent)
    {
        base_field_log[value] = (int)exponent;
        value = value * base_field_generator % (uint64_t)p;
    }
}

/* -------------------------------------------------------------------------
 * Direct indexing of H = N^{-1}({+1,-1})
 * ------------------------------------------------------------------------- */

static factor_exp_t *element_exponents(int exceptional_id)
{
    return exceptional_exponents +
           (size_t)exceptional_id * (size_t)factor_count;
}

/*
 * Encode an exponent tuple in H=N^{-1}({+-1}).  All field generators have
 * been normalized to have norm zeta, so membership is simply
 *
 *     e_0 + ... + e_{r-1} = 0 mod h,   h=(p-1)/2.
 *
 * One factor is chosen as the mixed-radix pivot.  Its exponent is reconstructed
 * from the other exponents and a reduced digit y of radix n_pivot/h:
 *
 *     e_pivot = h*y - sum_{i != pivot} e_i  mod n_pivot.
 *
 * The pivot is chosen among the factors of smallest degree.  This keeps the
 * only dependent coordinate and its modulus as small as possible, while the
 * larger field exponents remain direct mixed-radix digits.
 */
static h_code_t encode_h_exponents_known_member(
    const factor_exp_t *exponents)
{
    uint64_t pivot_order = factors[pivot_factor].order;
    uint64_t pivot_sum = 0;
    h_code_t code;
    int i;

    for (i = 0; i < factor_count; ++i)
        pivot_sum = add_mod_u64(pivot_sum, exponents[i], pivot_order);

    if (pivot_sum % norm_half_order != 0)
    {
        fprintf(stderr, "internal H-coordinate divisibility failure\n");
        exit(EXIT_FAILURE);
    }

    code = (pivot_sum / norm_half_order) *
           factors[pivot_factor].h_stride;
    for (i = 0; i < factor_count; ++i)
    {
        if (i != pivot_factor)
            code += exponents[i] * factors[i].h_stride;
    }

    if (code >= total_h_elements)
    {
        fprintf(stderr, "internal H-coordinate overflow\n");
        exit(EXIT_FAILURE);
    }
    return code;
}

static int encode_h_exponents(const factor_exp_t *exponents,
                              h_code_t *code_out)
{
    uint64_t residue = 0;
    int i;

    if (factor_count < 1 || code_out == NULL)
        return 0;

    if (norm_half_order > 1)
    {
        for (i = 0; i < factor_count; ++i)
            residue = add_mod_u64(residue, exponents[i], norm_half_order);
        if (residue != 0)
            return 0;
    }

    *code_out = encode_h_exponents_known_member(exponents);
    return 1;
}

static h_code_t encode_h_exponents_or_die(const factor_exp_t *exponents,
                                           const char *context)
{
#ifndef NDEBUG
    h_code_t code;

    if (!encode_h_exponents(exponents, &code))
    {
        fprintf(stderr, "%s left the norm-+-1 group\n", context);
        exit(EXIT_FAILURE);
    }
    return code;
#else
    (void)context;
    return encode_h_exponents_known_member(exponents);
#endif
}

static void decode_h_code(h_code_t code, factor_exp_t *exponents)
{
    uint64_t pivot_value;
    int i;

    if (code >= total_h_elements)
    {
        fprintf(stderr, "invalid H code\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < factor_count; ++i)
    {
        if (i != pivot_factor)
            exponents[i] =
                (code / factors[i].h_stride) % factors[i].h_radix;
    }

    pivot_value =
        norm_half_order *
        ((code / factors[pivot_factor].h_stride) %
         factors[pivot_factor].h_radix);

    for (i = 0; i < factor_count; ++i)
    {
        if (i != pivot_factor)
            pivot_value = subtract_mod_u64(
                pivot_value, exponents[i], factors[pivot_factor].order);
    }

    exponents[pivot_factor] = pivot_value;
}

static h_code_t quotient_h_code(int numerator_id, int denominator_id)
{
    const factor_exp_t *numerator = element_exponents(numerator_id);
    const factor_exp_t *denominator = element_exponents(denominator_id);
    uint64_t pivot_order = factors[pivot_factor].order;
    uint64_t pivot_sum = 0;
    h_code_t code = 0;
    int i;

    for (i = 0; i < factor_count; ++i)
    {
        factor_exp_t value =
            numerator[i] >= denominator[i]
                ? numerator[i] - denominator[i]
                : factors[i].order - (denominator[i] - numerator[i]);

        pivot_sum = add_mod_u64(pivot_sum, value, pivot_order);
        if (i != pivot_factor)
            code += value * factors[i].h_stride;
    }

#ifndef NDEBUG
    if (pivot_sum % norm_half_order != 0)
    {
        fprintf(stderr, "quotient left the norm-+-1 group\n");
        exit(EXIT_FAILURE);
    }
#endif
    return code + (pivot_sum / norm_half_order) *
                      factors[pivot_factor].h_stride;
}

static h_code_t inverse_h_code(int exceptional_id)
{
    const factor_exp_t *exponents = element_exponents(exceptional_id);
    uint64_t pivot_order = factors[pivot_factor].order;
    uint64_t pivot_sum = 0;
    h_code_t code = 0;
    int i;

    for (i = 0; i < factor_count; ++i)
    {
        factor_exp_t value =
            exponents[i] == 0 ? 0 : factors[i].order - exponents[i];

        pivot_sum = add_mod_u64(pivot_sum, value, pivot_order);
        if (i != pivot_factor)
            code += value * factors[i].h_stride;
    }

#ifndef NDEBUG
    if (pivot_sum % norm_half_order != 0)
    {
        fprintf(stderr, "inversion left the norm-+-1 group\n");
        exit(EXIT_FAILURE);
    }
#endif
    return code + (pivot_sum / norm_half_order) *
                      factors[pivot_factor].h_stride;
}

/* Return the H-code of -a*b/c, componentwise in the unit group. */
static h_code_t minus_times_product_quotient_h_code(int a_id, int b_id,
                                                     int c_id)
{
    const factor_exp_t *a = element_exponents(a_id);
    const factor_exp_t *b = element_exponents(b_id);
    const factor_exp_t *c = element_exponents(c_id);
    uint64_t pivot_order = factors[pivot_factor].order;
    uint64_t pivot_sum = 0;
    h_code_t code = 0;
    int i;

    for (i = 0; i < factor_count; ++i)
    {
        factor_exp_t value = factors[i].order / 2;

        value = add_mod_u64(value, a[i], factors[i].order);
        value = add_mod_u64(value, b[i], factors[i].order);
        value = subtract_mod_u64(value, c[i], factors[i].order);

        pivot_sum = add_mod_u64(pivot_sum, value, pivot_order);
        if (i != pivot_factor)
            code += value * factors[i].h_stride;
    }

#ifndef NDEBUG
    if (pivot_sum % norm_half_order != 0)
    {
        fprintf(stderr, "Möbius map left the norm-+-1 group\n");
        exit(EXIT_FAILURE);
    }
#endif
    return code + (pivot_sum / norm_half_order) *
                      factors[pivot_factor].h_stride;
}

static int exceptional_id_from_code(h_code_t code, const char *context)
{
    int id;

    if (code >= total_h_elements)
    {
        fprintf(stderr, "%s produced an invalid H code\n", context);
        exit(EXIT_FAILURE);
    }

    id = exceptional_id_by_code[(size_t)code];
    if (id < 0)
    {
        fprintf(stderr,
                "%s left the exceptional set (H-code=%" PRIu64 ")\n",
                context, code);
        exit(EXIT_FAILURE);
    }
    return id;
}

/* x-y is a norm-+-1 unit iff x/y-1 is exceptional. */
static inline int is_edge(int left_id, int right_id)
{
    h_code_t ratio_code = quotient_h_code(left_id, right_id);
    return bitset_contains_unchecked(&exceptional_codes,
                                     (size_t)ratio_code);
}

static void sort_small_int_array(int *values, int count)
{
    int i;

    for (i = 1; i < count; ++i)
    {
        int value = values[i];
        int j = i;

        while (j > 0 && values[j - 1] > value)
        {
            values[j] = values[j - 1];
            --j;
        }
        values[j] = value;
    }
}

static int compare_ints(const void *left, const void *right)
{
    int a = *(const int *)left;
    int b = *(const int *)right;
    return (a > b) - (a < b);
}

static int compare_v4_orbits(const void *left, const void *right)
{
    const v4_orbit_t *a = left;
    const v4_orbit_t *b = right;
    return (a->elements[0] > b->elements[0]) -
           (a->elements[0] < b->elements[0]);
}

/* The V4 images of x for the fixed quadruple {0,1,o,infinity}. */
/* -------------------------------------------------------------------------
 * Graph construction and clique search
 * ------------------------------------------------------------------------- */

static int v4_orbit_images(int base_id, int x_id, int images[4])
{
    int quotient_id;
    int one_minus_x;
    int one_minus_quotient;
    h_code_t code;
    int count;
    int i;

    images[0] = x_id;

    code = quotient_h_code(base_id, x_id);
    quotient_id = exceptional_id_from_code(code, "V4: o/x");
    images[1] = quotient_id;

    one_minus_x = one_minus_id[x_id];
    one_minus_quotient = one_minus_id[quotient_id];
    if (one_minus_x < 0 || one_minus_quotient < 0)
    {
        fprintf(stderr, "missing one-minus image in V4 computation\n");
        exit(EXIT_FAILURE);
    }

    /* (x-o)/(x-1) = -x(1-o/x)/(1-x). */
    code = minus_times_product_quotient_h_code(
        x_id, one_minus_quotient, one_minus_x);
    images[2] = exceptional_id_from_code(code, "V4: (x-o)/(x-1)");

    /* o(x-1)/(x-o) = -(o/x)(1-x)/(1-o/x). */
    code = minus_times_product_quotient_h_code(
        quotient_id, one_minus_x, one_minus_quotient);
    images[3] = exceptional_id_from_code(code, "V4: o(x-1)/(x-o)");

    sort_small_int_array(images, 4);

    count = 0;
    for (i = 0; i < 4; ++i)
    {
        if (i == 0 || images[i] != images[i - 1])
            images[count++] = images[i];
    }
    return count;
}

/*
 * Materialize the graph induced by the supplied exceptional-unit IDs and ask
 * Cliquer for its exact maximum clique.  The algebraic edge predicate stays
 * outside Cliquer, so no finite-field objects are retained in the graph.
 */
static int exact_internal_clique(const int *vertices, int vertex_count,
                                 int outer_orbit, int suborbit)
{
    graph_t *graph;
    clique_options options;
    int i;
    int j;
    int maximum;

    if (vertex_count == 0)
        return 0;

    graph = graph_new(vertex_count);
    if (graph == NULL)
    {
        fprintf(stderr, "graph_new failed [orbit %d, suborbit %d]\n",
                outer_orbit, suborbit);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < vertex_count; ++i)
    {
        for (j = i + 1; j < vertex_count; ++j)
        {
            if (is_edge(vertices[i], vertices[j]))
                GRAPH_ADD_EDGE(graph, i, j);
        }
    }

    options = *clique_default_options;
    options.time_function = NULL;
    options.output = NULL;
    maximum = clique_unweighted_max_weight(graph, &options);
    graph_free(graph);

    if (maximum <= 0)
    {
        fprintf(stderr,
                "Cliquer maximum search failed [orbit %d, suborbit %d]\n",
                outer_orbit, suborbit);
        exit(EXIT_FAILURE);
    }
    return maximum;
}

static int local_clique_direct(int outer_orbit,
                               const int *canonical_candidates,
                               int canonical_count)
{
    if (canonical_count == 0)
        return 3;

    return 3 + exact_internal_clique(canonical_candidates, canonical_count,
                                     outer_orbit, -1);
}

/*
 * Quotient the canonical candidate set by the Klein-four action attached to
 * {0,1,o,infinity}.  The full V4 orbit of a candidate need not stay inside the
 * canonical set: some images may belong to outer symmetry orbits that were
 * processed earlier.  For duplicate suppression we only need the nonempty
 * intersection
 *
 *     (V4 . x) intersect canonical_candidates.
 *
 * Every candidate in that intersection represents the same choice of the
 * fourth normalized point, so exactly one of them is used as a representative.
 * Images outside the candidate set are simply ignored.
 */
static int local_clique_v4(int outer_orbit, int base_id,
                           int *canonical_candidates,
                           int canonical_count, local_stats_t *stats,
                           FILE *progress_stream)
{
    /*
     * This function takes ownership of canonical_candidates.  Once the V4
     * classes have been materialized, that array is no longer needed and is
     * released before any maximum-clique subproblem is solved.
     */
    bitset_t remaining_candidates = {0};
    v4_orbit_t *suborbits;
    int suborbit_capacity;
    int suborbit_count = 0;
    int partitioned_vertex_count = 0;
    int *subproblem_vertices;
    int local_maximum = 3;
    int i;
    int k;

    if (!bitset_init(&remaining_candidates, (size_t)exceptional_count))
    {
        fprintf(stderr, "could not allocate local V4 bitset [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < canonical_count; ++i)
        bitset_add(&remaining_candidates, (size_t)canonical_candidates[i]);

    suborbit_capacity = (canonical_count + 3) / 4;
    if (suborbit_capacity < 1)
        suborbit_capacity = 1;

    suborbits = malloc((size_t)suborbit_capacity * sizeof(*suborbits));
    if (suborbits == NULL)
    {
        fprintf(stderr, "could not allocate local V4 orbits [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < canonical_count; ++i)
    {
        int x_id = canonical_candidates[i];
        int images[4];
        int intersection[4];
        int image_count;
        int intersection_count = 0;
        int j;

        if (!bitset_contains_unchecked(&remaining_candidates, (size_t)x_id))
            continue;

        image_count = v4_orbit_images(base_id, x_id, images);

        /*
         * v4_orbit_images returns the distinct full V4 orbit in sorted order.
         * remaining_candidates contains precisely the canonical vertices not
         * assigned to a previous V4 class.  Since full V4 orbits are disjoint,
         * filtering through it gives the desired intersection and removing the
         * images simultaneously marks the class as processed.
         */
        for (j = 0; j < image_count; ++j)
        {
            int image = images[j];

            if (!bitset_contains_unchecked(&remaining_candidates,
                                              (size_t)image))
                continue;

            intersection[intersection_count++] = image;
            bitset_remove(&remaining_candidates, (size_t)image);
        }

        /*
         * x itself is a V4 image and belongs to the candidate set, so the
         * intersection cannot be empty.  Treat failure as an internal error.
         */
        if (intersection_count == 0)
        {
            fprintf(stderr,
                    "empty V4 intersection [orbit %d, candidate %d]\n",
                    outer_orbit + 1, x_id);
            exit(EXIT_FAILURE);
        }

        if (suborbit_count == suborbit_capacity)
        {
            v4_orbit_t *grown;
            suborbit_capacity *= 2;
            grown = realloc(suborbits,
                            (size_t)suborbit_capacity * sizeof(*suborbits));
            if (grown == NULL)
            {
                fprintf(stderr,
                        "could not grow local V4 orbit list [orbit %d]\n",
                        outer_orbit + 1);
                exit(EXIT_FAILURE);
            }
            suborbits = grown;
        }

        suborbits[suborbit_count].size = intersection_count;
        memcpy(suborbits[suborbit_count].elements, intersection,
               (size_t)intersection_count * sizeof(*intersection));
        partitioned_vertex_count += intersection_count;
        ++suborbit_count;
    }

    if (partitioned_vertex_count != canonical_count)
    {
        fprintf(stderr,
                "V4 candidate partition mismatch [orbit %d]: %d != %d\n",
                outer_orbit + 1, partitioned_vertex_count, canonical_count);
        exit(EXIT_FAILURE);
    }

    /* Classification is finished: neither object is needed by Cliquer. */
    bitset_destroy(&remaining_candidates);
    free(canonical_candidates);

    qsort(suborbits, (size_t)suborbit_count, sizeof(*suborbits),
          compare_v4_orbits);

    stats->retained_suborbits = suborbit_count;
    stats->iteration_count = suborbit_count;

    /*
     * Orbit classification is complete before the first clique search starts,
     * so expose the number of V4 suborbits immediately.  This makes a local
     * file useful for monitoring even when the first Cliquer call is long.
     */
    if (progress_stream != NULL)
    {
        fprintf(progress_stream, "Suborbits: %d\n\n", suborbit_count);
        flush_progress(progress_stream, outer_orbit);
    }

    if (suborbit_count > 0)
    {
        stats->iterations = calloc((size_t)suborbit_count,
                                   sizeof(*stats->iterations));
        if (stats->iterations == NULL)
        {
            fprintf(stderr,
                    "could not allocate V4 statistics [orbit %d]\n",
                    outer_orbit + 1);
            exit(EXIT_FAILURE);
        }
    }

    if (suborbit_count == 0)
    {
        free(suborbits);
        return 3;
    }

    subproblem_vertices = malloc((size_t)canonical_count *
                                 sizeof(*subproblem_vertices));
    if (subproblem_vertices == NULL)
    {
        fprintf(stderr,
                "could not allocate V4 subproblem [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    for (k = 0; k < suborbit_count; ++k)
    {
        int representative = suborbits[k].elements[0];
        int subproblem_count = 0;
        int internal_maximum;
        int candidate_maximum;
        int t;
        int j;
        struct timespec iteration_start = monotonic_now();
        v4_iteration_stats_t *iteration = &stats->iterations[k];

        /*
         * Only suborbits at or after k are considered, so every clique is
         * assigned to the first V4 suborbit that it meets.  Within those
         * suborbits, retain exactly the vertices adjacent to the fixed
         * representative.
         */
        for (t = k; t < suborbit_count; ++t)
        {
            for (j = 0; j < suborbits[t].size; ++j)
            {
                int candidate = suborbits[t].elements[j];
                if (is_edge(representative, candidate))
                    subproblem_vertices[subproblem_count++] = candidate;
            }
        }

        iteration->candidate_count = subproblem_count;

        /* This local cardinality bound cannot improve the best clique found. */
        if (4 + subproblem_count <= local_maximum)
        {
            iteration->skipped = 1;
            iteration->elapsed = elapsed_seconds(iteration_start);
            if (progress_stream != NULL)
                append_v4_progress_line(
                    progress_stream, outer_orbit,
                    "Suborbit %d/%d: Vertices=%d  skipped  time=%.*f s",
                    k + 1, suborbit_count, subproblem_count,
                    TIME_PRECISION, iteration->elapsed);
            continue;
        }

        /*
         * Show the known subgraph size before entering Cliquer, then overwrite
         * this same physical line with the result when the search returns.
         */
        long progress_position = -1;
        if (progress_stream != NULL)
            progress_position = append_v4_progress_line(
                progress_stream, outer_orbit,
                "Suborbit %d/%d: Vertices=%d  running",
                k + 1, suborbit_count, subproblem_count);

        internal_maximum = exact_internal_clique(
            subproblem_vertices, subproblem_count, outer_orbit, k);
        candidate_maximum = 4 + internal_maximum;
        iteration->clique_size = candidate_maximum;
        iteration->elapsed = elapsed_seconds(iteration_start);

        if (candidate_maximum > local_maximum)
            local_maximum = candidate_maximum;

        if (progress_stream != NULL)
            replace_v4_progress_line(
                progress_stream, progress_position, outer_orbit,
                "Suborbit %d/%d: Vertices=%d  clique=%d  time=%.*f s",
                k + 1, suborbit_count, subproblem_count, candidate_maximum,
                TIME_PRECISION, iteration->elapsed);
    }

    free(subproblem_vertices);
    free(suborbits);
    return local_maximum;
}
/*
 * Compute the local clique number attached to one outer symmetry orbit.
 *
 * Once the canonical graph is known, the worker creates local_N.out and flushes
 * Orbit, Method and Vertices before any Cliquer call.  In V4 mode it also writes
 * the suborbit count and updates one fixed line per subproblem in place.
 *
 * After all local searches finish, the worker rewrites the file in its compact
 * final format.  The parent then appends the global Best value and
 * Status: complete, so interrupted or missing orbit computations remain visible
 * on disk while completed files keep the stable output format documented below.
 */
static int local_clique(int outer_orbit, orbit_result_t *result_out)
{
    FILE *stream;
    char path[MAXPATHLEN];
    int *canonical_candidates = NULL;
    int canonical_count = 0;
    int local_clique_number;
    int base_id;
    int use_v4;
    int i;
    int j;
    struct timespec start = monotonic_now();
    double elapsed;
    local_stats_t stats = {0};

    /*
     * The local graph assigned to outer orbit i contains only vertices from
     * orbit i and later orbits.  A vertex x is therefore retained precisely
     * when x is exceptional, x is adjacent to the chosen base o, and the
     * outer symmetry orbit of x has not already been processed.
     *
     * This canonical restriction is essential: it assigns every normalized
     * clique to the first outer orbit that represents it and makes the local
     * graph sizes decrease as the outer-orbit loop advances.
     */
    {
        size_t candidate_capacity = 0;

        for (i = outer_orbit; i < orbit_count; ++i)
            candidate_capacity += (size_t)orbits[i].size;

        canonical_candidates = malloc(candidate_capacity *
                                      sizeof(*canonical_candidates));
        if (canonical_candidates == NULL && candidate_capacity != 0)
        {
            fprintf(stderr, "candidate allocation failed [orbit %d]\n",
                    outer_orbit + 1);
            exit(EXIT_FAILURE);
        }
    }

    base_id = orbits[outer_orbit].elements[0];

    checked_snprintf(path, sizeof(path), "%s/local_%d.out",
                     local_path, outer_orbit + 1);

    /*
     * Every exceptional candidate already satisfies N(x)=N(1-x)=+-1.
     * is_edge(base_id, candidate) adds the remaining condition
     * N(o-x)=+-1.  Starting at outer_orbit excludes all earlier symmetry
     * orbits from the actual local graph.
     */
    for (i = outer_orbit; i < orbit_count; ++i)
    {
        for (j = 0; j < orbits[i].size; ++j)
        {
            int candidate = orbits[i].elements[j];

            if (is_edge(base_id, candidate))
                canonical_candidates[canonical_count++] = candidate;
        }
    }

    /*
     * The V4 threshold is applied to the canonical graph just constructed.
     * Hence sufficiently late outer orbits automatically return to the direct
     * Cliquer path as the number of remaining candidates decreases.
     */
    use_v4 = v4_threshold > 0 && canonical_count >= v4_threshold;

    /*
     * At this point the local graph and the algorithm choice are both known,
     * while no maximum-clique search has started yet.  Rewrite the progress
     * file with every piece of information already available and flush it so
     * users can inspect long-running jobs in real time.
     */
    stream = fopen(path, "w");
    if (stream == NULL)
    {
        fprintf(stderr, "fopen failed [orbit %d]\n", outer_orbit + 1);
        exit(EXIT_FAILURE);
    }
    fprintf(stream, "Status: running\n");
    fprintf(stream, "Orbit: %d/%d\n", outer_orbit + 1, orbit_count);
    fprintf(stream, "Method: %s\n", use_v4 ? "V4" : "direct");
    fprintf(stream, "Vertices: %d\n", canonical_count);
    flush_progress(stream, outer_orbit);

    if (use_v4)
    {
        stats.used_v4 = 1;
        local_clique_number = local_clique_v4(
            outer_orbit, base_id,
            canonical_candidates, canonical_count, &stats, stream);
        canonical_candidates = NULL;
    }
    else
    {
        local_clique_number = local_clique_direct(
            outer_orbit, canonical_candidates, canonical_count);
    }

    if (fclose(stream) != 0)
    {
        fprintf(stderr, "fclose failed [orbit %d]\n", outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    elapsed = elapsed_seconds(start);

    stream = fopen(path, "w");
    if (stream == NULL)
    {
        fprintf(stderr, "fopen failed [orbit %d]\n", outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    fprintf(stream, "Clique: %d\n", local_clique_number);
    fprintf(stream, "Time: %.*f s\n", TIME_PRECISION, elapsed);
    fprintf(stream, "Orbit: %d/%d\n", outer_orbit + 1, orbit_count);
    fprintf(stream, "Method: %s\n", stats.used_v4 ? "V4" : "direct");
    fprintf(stream, "Vertices: %d\n", canonical_count);

    if (stats.used_v4)
    {
        fprintf(stream, "Suborbits: %d\n\n", stats.retained_suborbits);

        for (i = 0; i < stats.iteration_count; ++i)
        {
            const v4_iteration_stats_t *iteration = &stats.iterations[i];

            fprintf(stream,
                    "Suborbit %d/%d: Vertices=%d  ",
                    i + 1, stats.iteration_count,
                    iteration->candidate_count);

            if (iteration->skipped)
                fprintf(stream, "skipped  ");
            else
                fprintf(stream, "clique=%d  ", iteration->clique_size);

            fprintf(stream, "time=%.*f s\n", TIME_PRECISION, iteration->elapsed);
        }
    }

    free(stats.iterations);
    free(canonical_candidates);

    if (fclose(stream) != 0)
    {
        fprintf(stderr, "fclose failed [orbit %d]\n", outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    if (result_out != NULL)
    {
        result_out->clique_number = local_clique_number;
        result_out->vertex_count = canonical_count;
        result_out->elapsed = elapsed;
    }
    return local_clique_number;
}

static void append_orbit_result(const char *result_path, int outer_orbit,
                                const orbit_result_t *result,
                                int best_clique)
{
    FILE *stream = fopen(result_path, "a");

    if (stream == NULL)
    {
        fprintf(stderr, "could not append result.out\n");
        exit(EXIT_FAILURE);
    }

    fprintf(stream,
            "%d/%d  Vertices: %d  Clique: %d  Best: %d  Time: %.*f s\n",
            outer_orbit + 1, orbit_count, result->vertex_count,
            result->clique_number, best_clique, TIME_PRECISION, result->elapsed);

    if (fclose(stream) != 0)
    {
        fprintf(stderr, "could not close result.out\n");
        exit(EXIT_FAILURE);
    }
}

/*
 * Read the worker-produced numerical fields.  The parent does not trust a
 * successful exit status alone: a truncated or malformed local file makes the
 * complete computation fail instead of silently omitting an orbit.
 */
static void read_local_result(int outer_orbit, orbit_result_t *result)
{
    char path[MAXPATHLEN];
    char line[256];
    FILE *stream;
    int got_clique = 0;
    int got_elapsed = 0;
    int got_vertices = 0;

    checked_snprintf(path, sizeof(path), "%s/local_%d.out",
                     local_path, outer_orbit + 1);
    stream = fopen(path, "r");
    if (stream == NULL)
    {
        fprintf(stderr, "could not read local result [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, sizeof(line), stream) != NULL)
    {
        int integer_value;
        double double_value;

        if (sscanf(line, "Clique: %d", &integer_value) == 1)
        {
            result->clique_number = integer_value;
            got_clique = 1;
        }
        else if (sscanf(line, "Time: %lf s", &double_value) == 1)
        {
            result->elapsed = double_value;
            got_elapsed = 1;
        }
        else if (sscanf(line, "Vertices: %d", &integer_value) == 1)
        {
            result->vertex_count = integer_value;
            got_vertices = 1;
        }
    }

    if (ferror(stream) || !got_clique || !got_elapsed || !got_vertices)
    {
        fprintf(stderr, "invalid local result [orbit %d]\n",
                outer_orbit + 1);
        fclose(stream);
        exit(EXIT_FAILURE);
    }

    if (fclose(stream) != 0)
    {
        fprintf(stderr, "could not close local result [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }
}

/*
 * The parent appends information that depends on completion order.  In a
 * parallel run the worker cannot know the current global best, and only the
 * parent can certify that the result was received and accounted for.
 */
static void finalize_local_result(int outer_orbit, int best_clique)
{
    char path[MAXPATHLEN];
    FILE *stream;

    checked_snprintf(path, sizeof(path), "%s/local_%d.out",
                     local_path, outer_orbit + 1);
    stream = fopen(path, "a");
    if (stream == NULL)
    {
        fprintf(stderr, "could not finalize local result [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    fprintf(stream, "Best: %d\n", best_clique);
    fprintf(stream, "Status: complete\n");

    if (fclose(stream) != 0)
    {
        fprintf(stderr, "could not close local result [orbit %d]\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }
}

static void terminate_workers(const pid_t *pid_table, int process_count,
                              pid_t except)
{
    int i;

    for (i = 0; i < process_count; ++i)
    {
        pid_t pid = pid_table[i];
        if (pid > 0 && pid != except)
            kill(pid, SIGTERM);
    }

    while (waitpid(-1, NULL, 0) > 0)
        ;
}

/*
 * Record one completed orbit exactly once.  The completion bitmap belongs to
 * the parent process and is therefore independent of the order in which
 * parallel workers finish.
 */
static void record_completed_orbit(const char *result_path, int outer_orbit,
                                   const orbit_result_t *result,
                                   unsigned char *completed,
                                   int *completed_count, int *best_clique)
{
    if (outer_orbit < 0 || outer_orbit >= orbit_count ||
        completed[outer_orbit])
    {
        fprintf(stderr, "duplicate or invalid orbit completion: %d\n",
                outer_orbit + 1);
        exit(EXIT_FAILURE);
    }

    if (result->clique_number > *best_clique)
        *best_clique = result->clique_number;

    completed[outer_orbit] = 1;
    ++*completed_count;

    finalize_local_result(outer_orbit, *best_clique);
    append_orbit_result(result_path, outer_orbit, result, *best_clique);
}

static int compute_all_orbits(const char *result_path, int cpus,
                              int *completed_count_out)
{
    int best_clique = exceptional_count > 0 ? 3 : 2;
    int completed_count = 0;
    unsigned char *completed = calloc((size_t)orbit_count, sizeof(*completed));
    int i;

    if (completed == NULL && orbit_count != 0)
    {
        fprintf(stderr, "could not allocate orbit-completion table\n");
        exit(EXIT_FAILURE);
    }

    if (cpus <= 1)
    {
        for (i = 0; i < orbit_count; ++i)
        {
            orbit_result_t result = {0};

            local_clique(i, &result);

            /* Later local graphs start at orbit i+1, so orbit i is now dead. */
            free(orbits[i].elements);
            orbits[i].elements = NULL;

            record_completed_orbit(result_path, i, &result, completed,
                                   &completed_count, &best_clique);
        }
    }
    else
    {
        int active = 0;
        int status;
        pid_t *pid_table = calloc((size_t)orbit_count, sizeof(*pid_table));

        if (pid_table == NULL && orbit_count != 0)
        {
            fprintf(stderr, "could not allocate process table\n");
            free(completed);
            exit(EXIT_FAILURE);
        }

        /*
         * Launch orbit jobs in canonical order, but consume them in completion
         * order.  pid_table maps each child back to its orbit number.
         */
        for (i = 0; i < orbit_count; ++i)
        {
            pid_t pid = fork();

            if (pid == 0)
            {
                if (freopen(error_path, "a", stderr) == NULL)
                    _exit(EXIT_FAILURE);
                local_clique(i, NULL);
                _exit(EXIT_SUCCESS);
            }
            if (pid < 0)
            {
                fprintf(stderr, "fork failed\n");
                terminate_workers(pid_table, i, -1);
                free(pid_table);
                free(completed);
                exit(EXIT_FAILURE);
            }

            pid_table[i] = pid;
            ++active;

            /*
             * The child inherited orbit i before this free, while every future
             * child starts at a later orbit and will never inspect it.
             */
            free(orbits[i].elements);
            orbits[i].elements = NULL;

            /*
             * Once the last worker has been forked, the parent no longer needs
             * any graph/algebra data.  Children keep their private inherited
             * address spaces, so releasing the parent's copies is safe and
             * lets the OS reclaim pages as workers finish.
             */
            if (i == orbit_count - 1)
                free_global_data();

            while (active >= cpus ||
                   (i == orbit_count - 1 && active > 0))
            {
                pid_t finished = wait(&status);
                int finished_orbit = -1;
                int j;
                orbit_result_t result = {0};

                if (finished < 0)
                {
                    fprintf(stderr, "wait failed\n");
                    terminate_workers(pid_table, i + 1, -1);
                    free(pid_table);
                    free(completed);
                    exit(EXIT_FAILURE);
                }

                if (!WIFEXITED(status) ||
                    WEXITSTATUS(status) != EXIT_SUCCESS)
                {
                    fprintf(stderr,
                            "subprocess %ld failed; terminating workers\n",
                            (long)finished);
                    fflush(stderr);
                    terminate_workers(pid_table, i + 1, finished);
                    free(pid_table);
                    free(completed);
                    exit(EXIT_FAILURE);
                }

                --active;
                for (j = 0; j <= i; ++j)
                {
                    if (pid_table[j] == finished)
                    {
                        finished_orbit = j;
                        pid_table[j] = 0;
                        break;
                    }
                }

                if (finished_orbit < 0)
                {
                    fprintf(stderr, "completed process not found\n");
                    terminate_workers(pid_table, i + 1, finished);
                    free(pid_table);
                    free(completed);
                    exit(EXIT_FAILURE);
                }

                read_local_result(finished_orbit, &result);
                record_completed_orbit(result_path, finished_orbit, &result,
                                       completed, &completed_count,
                                       &best_clique);
            }
        }

        free(pid_table);
    }

    /* A successful run must account for every orbit exactly once. */
    for (i = 0; i < orbit_count; ++i)
    {
        if (!completed[i])
        {
            fprintf(stderr, "orbit %d/%d was not completed\n",
                    i + 1, orbit_count);
            free(completed);
            exit(EXIT_FAILURE);
        }
    }

    free(completed);
    *completed_count_out = completed_count;
    return best_clique;
}

/* -------------------------------------------------------------------------
 * Finite-field preparation
 * ------------------------------------------------------------------------- */

static int find_model(const field_model_t *models, int model_count, int degree)
{
    int i;
    for (i = 0; i < model_count; ++i)
    {
        if (models[i].degree == degree)
            return i;
    }
    return -1;
}

/*
 * Build one reusable model for each distinct factor degree.  Equal-degree
 * factors deliberately share the same generator and logarithm table, which
 * makes their algebra automorphisms literal permutations of exponent slots.
 * The full e -> log_g(1-g^e) table is skipped for the specialized one-factor
 * path, where only subgroup elements are visited.
 */
static void initialize_field_model(field_model_t *model, int degree,
                                   FILE *status_stream,
                                   int build_full_one_minus_table)
{
    uint64_t norm_base;
    uint64_t norm_log;
    uint64_t generator_power;
    uint64_t value;
    uint64_t exponent;
    uint64_t progress_step;
    char stage[128];

    memset(model, 0, sizeof(*model));
    model->degree = degree;
    model->q = checked_u64_power((uint64_t)p, degree);
    if (model->q < 2)
    {
        fprintf(stderr, "invalid field cardinality\n");
        exit(EXIT_FAILURE);
    }
    model->order = model->q - 1;

    if (build_full_one_minus_table)
    {
        if (model->order > SIZE_MAX / sizeof(*model->one_minus_exp))
        {
            fprintf(stderr, "field factor is too large for this build\n");
            exit(EXIT_FAILURE);
        }

        model->one_minus_exp = malloc((size_t)model->order *
                                      sizeof(*model->one_minus_exp));
        if (model->one_minus_exp == NULL)
        {
            fprintf(stderr,
                    "could not allocate local table for degree %d\n",
                    degree);
            exit(EXIT_FAILURE);
        }

        for (exponent = 0; exponent < model->order; ++exponent)
            model->one_minus_exp[exponent] = UNSET_EXP;
    }

    if (degree == 1)
    {
        model->prime_field = 1;

        if (build_full_one_minus_table)
        {
            value = 1;
            for (exponent = 0; exponent < model->order; ++exponent)
            {
                uint64_t one_minus =
                    (1 + (uint64_t)p - value) % (uint64_t)p;

                if (one_minus != 0)
                    model->one_minus_exp[exponent] =
                        (factor_exp_t)base_field_log[one_minus];
                value = value * base_field_generator % (uint64_t)p;
            }
        }
        return;
    }

    model->field_generator = ffgen(utoipos(model->q), -1);
    model->primitive_generator =
        FF_primroot(model->field_generator, &model->order_data);
    model->one = FF_1(model->field_generator);

    norm_base = (uint64_t)itos(FF_norm(model->primitive_generator));
    if (norm_base == 0 || norm_base >= (uint64_t)p ||
        base_field_log[norm_base] < 0)
    {
        fprintf(stderr, "invalid primitive-generator norm in degree %d\n",
                degree);
        exit(EXIT_FAILURE);
    }

    /*
     * Normalize every field model to the same base-field generator zeta.
     * This removes all per-factor norm coefficients from the H coordinates:
     * for x=(g_i^{e_i}), N(x)=zeta^{sum e_i}.
     */
    norm_log = (uint64_t)base_field_log[norm_base];
    generator_power = normalized_generator_exponent(norm_log, model->order);
    if (generator_power != 1)
        model->primitive_generator =
            FF_pow(model->primitive_generator, utoipos(generator_power));

#ifndef NDEBUG
    norm_base = (uint64_t)itos(FF_norm(model->primitive_generator));
    if (norm_base != base_field_generator)
    {
        fprintf(stderr,
                "normalized generator has wrong norm in degree %d\n",
                degree);
        exit(EXIT_FAILURE);
    }
#endif

    if (!build_full_one_minus_table)
        return;

    /*
     * Precompute e -> log_g(1-g^e) once for this field model.  The
     * involution 1-(1-x)=x halves the number of discrete logarithms.
     */
    if (model->order > 1)
    {
        GEN x = model->primitive_generator;

        progress_step = model->order / 20;
        if (progress_step < 1)
            progress_step = 1;

        for (exponent = 1; exponent < model->order; ++exponent)
        {
            pari_sp av = avma;
            GEN next;

            if (status_stream != NULL &&
                (exponent == 1 || exponent % progress_step == 0))
            {
                checked_snprintf(stage, sizeof(stage),
                                 "field degree %d table: %" PRIu64
                                 "/%" PRIu64,
                                 degree, exponent, model->order - 1);
                write_preparation_status(status_stream, stage);
            }

            if (model->one_minus_exp[exponent] == UNSET_EXP)
            {
                GEN y = FF_sub(model->one, x);
                GEN logarithm = FF_log(y, model->primitive_generator,
                                       model->order_data);
                factor_exp_t mate = (factor_exp_t)itou(logarithm);

                mate %= model->order;
                if (mate >= model->order || mate == 0)
                {
                    fprintf(stderr,
                            "invalid local one-minus logarithm in degree %d\n",
                            degree);
                    exit(EXIT_FAILURE);
                }
                if (model->one_minus_exp[mate] != UNSET_EXP &&
                    model->one_minus_exp[mate] != exponent)
                {
                    fprintf(stderr,
                            "inconsistent local one-minus involution in degree %d\n",
                            degree);
                    exit(EXIT_FAILURE);
                }

                model->one_minus_exp[exponent] = mate;
                model->one_minus_exp[mate] = exponent;
            }

            next = FF_mul(x, model->primitive_generator);
            x = gerepilecopy(av, next);
        }
    }

    for (exponent = 1; exponent < model->order; ++exponent)
    {
        factor_exp_t mate = model->one_minus_exp[exponent];

        if (mate == UNSET_EXP || mate >= model->order ||
            model->one_minus_exp[mate] != exponent)
        {
            fprintf(stderr,
                    "incomplete local one-minus table in degree %d\n",
                    degree);
            exit(EXIT_FAILURE);
        }
    }
}

static void free_field_models(field_model_t *models, int model_count)
{
    int i;
    for (i = 0; i < model_count; ++i)
        free(models[i].one_minus_exp);
    free(models);
}

static int norm_is_plus_or_minus_one(int norm)
{
    return norm == 1 || norm == p - 1;
}

static void mark_exceptional(h_code_t code, int *count)
{
    if (bitset_contains_unchecked(&exceptional_codes, (size_t)code))
        return;

    if (*count == INT_MAX)
    {
        fprintf(stderr, "too many exceptional units for Cliquer\n");
        exit(EXIT_FAILURE);
    }

    bitset_add(&exceptional_codes, (size_t)code);
    ++*count;
}

static void record_exceptional_pair(h_code_t *one_minus_code,
                                    h_code_t code, h_code_t mate,
                                    int *count)
{
    if (code >= total_h_elements || mate >= total_h_elements)
    {
        fprintf(stderr, "invalid one-minus H code\n");
        exit(EXIT_FAILURE);
    }

    if (one_minus_code[(size_t)code] != UNSET_CODE &&
        one_minus_code[(size_t)code] != mate)
    {
        fprintf(stderr, "inconsistent one-minus image\n");
        exit(EXIT_FAILURE);
    }
    if (one_minus_code[(size_t)mate] != UNSET_CODE &&
        one_minus_code[(size_t)mate] != code)
    {
        fprintf(stderr, "inconsistent one-minus involution\n");
        exit(EXIT_FAILURE);
    }

    one_minus_code[(size_t)code] = mate;
    one_minus_code[(size_t)mate] = code;
    mark_exceptional(code, count);
    mark_exceptional(mate, count);
}

/*
 * Specialized preparation for A=F_{p^d}.  In this case
 *
 *     H=N^{-1}({+-1})=<g^((p-1)/2)>,
 *
 * so enumerate H directly and compute a discrete logarithm only when
 * 1-x also lies in H.  This avoids constructing the full q-1 entry local
 * table e -> log_g(1-g^e), which is useful only in the multifactor case.
 */
static int enumerate_single_field_exceptional(
    const field_model_t *model, h_code_t *one_minus_code,
    FILE *status_stream)
{
    uint64_t progress_step = total_h_elements / 20;
    h_code_t code;
    int count = 0;
    char stage[160];

    if (progress_step < 1)
        progress_step = 1;

    if (model->prime_field)
    {
        uint64_t subgroup_generator = modular_power(
            base_field_generator, norm_half_order, (uint64_t)p);
        uint64_t x = 1;

        for (code = 0; code < total_h_elements; ++code)
        {
            if (status_stream != NULL &&
                (code == 0 || code % progress_step == 0))
            {
                checked_snprintf(stage, sizeof(stage),
                                 "single-field subgroup: %" PRIu64
                                 "/%" PRIu64,
                                 code, total_h_elements);
                write_preparation_status(status_stream, stage);
            }

            if (one_minus_code[(size_t)code] == UNSET_CODE)
            {
                uint64_t one_minus =
                    (1 + (uint64_t)p - x) % (uint64_t)p;

                if (norm_is_plus_or_minus_one((int)one_minus))
                {
                    int mate_exponent = base_field_log[one_minus];
                    h_code_t mate;

                    if (mate_exponent < 0 ||
                        (uint64_t)mate_exponent % norm_half_order != 0)
                    {
                        fprintf(stderr,
                                "invalid prime-field subgroup logarithm\n");
                        exit(EXIT_FAILURE);
                    }
                    mate = (h_code_t)mate_exponent / norm_half_order;
                    record_exceptional_pair(one_minus_code, code, mate,
                                            &count);
                }
            }

            x = x * subgroup_generator % (uint64_t)p;
        }

        if (x != 1)
        {
            fprintf(stderr, "prime-field subgroup order mismatch\n");
            exit(EXIT_FAILURE);
        }
        return count;
    }
    else
    {
        GEN subgroup_generator = FF_pow(
            model->primitive_generator, utoipos(norm_half_order));
        GEN subgroup_order = utoipos(total_h_elements);
        GEN subgroup_order_data = mkvec2(
            subgroup_order, factorint(subgroup_order, 0));
        GEN x = model->one;

        for (code = 0; code < total_h_elements; ++code)
        {
            pari_sp av = avma;
            GEN next;

            if (status_stream != NULL &&
                (code == 0 || code % progress_step == 0))
            {
                checked_snprintf(stage, sizeof(stage),
                                 "single-field subgroup: %" PRIu64
                                 "/%" PRIu64,
                                 code, total_h_elements);
                write_preparation_status(status_stream, stage);
            }

            if (one_minus_code[(size_t)code] == UNSET_CODE)
            {
                GEN one_minus = FF_sub(model->one, x);
                int norm = itos(FF_norm(one_minus));

                if (norm_is_plus_or_minus_one(norm))
                {
                    GEN logarithm = FF_log(one_minus, subgroup_generator,
                                           subgroup_order_data);
                    h_code_t mate = (h_code_t)itou(logarithm);

                    mate %= total_h_elements;
                    record_exceptional_pair(one_minus_code, code, mate,
                                            &count);
                }
            }

            next = FF_mul(x, subgroup_generator);
            x = gerepilecopy(av, next);
        }

        if (!FF_equal1(x))
        {
            fprintf(stderr, "single-field subgroup order mismatch\n");
            exit(EXIT_FAILURE);
        }
        return count;
    }
}

static h_code_t compute_one_minus_code(const field_model_t *models,
                                          const factor_exp_t *exponents,
                                          factor_exp_t *mate_exponents,
                                          int *is_exceptional)
{
    h_code_t mate_code;
    int i;

    *is_exceptional = 0;

    for (i = 0; i < factor_count; ++i)
    {
        const field_model_t *model = &models[factors[i].model_index];
        factor_exp_t mate = model->one_minus_exp[(size_t)exponents[i]];

        if (mate == UNSET_EXP)
            return 0;
        mate_exponents[i] = mate;
    }

    if (!encode_h_exponents(mate_exponents, &mate_code))
        return 0;

    *is_exceptional = 1;
    return mate_code;
}

/* -------------------------------------------------------------------------
 * Symmetry-orbit classification
 * ------------------------------------------------------------------------- */

static swap_generator_t *build_swap_generators(int *swap_count_out)
{
    swap_generator_t *swaps;
    int capacity = factor_count > 1 ? factor_count - 1 : 1;
    int count = 0;
    int i;

    swaps = malloc((size_t)capacity * sizeof(*swaps));
    if (swaps == NULL)
    {
        fprintf(stderr, "could not allocate factor-permutation generators\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < factor_count; ++i)
    {
        int previous = -1;
        int j;

        for (j = i - 1; j >= 0; --j)
        {
            if (factors[j].degree == factors[i].degree)
            {
                previous = j;
                break;
            }
        }

        if (previous >= 0)
        {
            swaps[count].left = previous;
            swaps[count].right = i;
            ++count;
        }
    }

    *swap_count_out = count;
    return swaps;
}

/*
 * In exponent coordinates the local Frobenius is diagonal:
 * g_i^{e_i} -> g_i^{p e_i}, hence e_i -> p*e_i mod (p^{d_i}-1).
 * Since p=1 mod (p-1), it preserves the norm congruence; re-encoding only
 * updates the mixed-radix pivot digit when necessary.
 */
static h_code_t frobenius_code(int exceptional_id, int factor_index)
{
    const factor_exp_t *exponents = element_exponents(exceptional_id);
    factor_exp_t values[factor_count];

    memcpy(values, exponents, (size_t)factor_count * sizeof(*values));
    values[factor_index] = multiply_mod_u64(
        (uint64_t)p, values[factor_index], factors[factor_index].order);

    return encode_h_exponents_or_die(values, "local Frobenius");
}

/*
 * Equal-degree factors share the same field model and normalized generator.
 * Their algebra permutation is therefore literally an exponent-slot swap,
 * followed by re-encoding (important when one of the slots is the pivot).
 */
static h_code_t swap_code(int exceptional_id, int left, int right)
{
    const factor_exp_t *exponents = element_exponents(exceptional_id);
    factor_exp_t values[factor_count];
    factor_exp_t temporary;

    memcpy(values, exponents, (size_t)factor_count * sizeof(*values));
    temporary = values[left];
    values[left] = values[right];
    values[right] = temporary;

    return encode_h_exponents_or_die(values, "equal-factor permutation");
}

/*
 * Classify exceptional units under generators of the full symmetry group:
 * x -> 1-x, inversion, one Frobenius generator per non-prime factor, and
 * adjacent swaps among equal-degree factors.  Breadth-first search produces
 * disjoint orbit lists and records the canonical orbit index of every vertex.
 */
static void classify_outer_orbits(const swap_generator_t *swaps,
                                  int swap_count)
{
    int galois_count = 0;
    int *galois_factors;
    int generator_count;
    int *generator_images;
    bitset_t unclassified = {0};
    int *orbit_buffer;
    size_t orbit_capacity = 128;
    int i;
    int g;

    if (exceptional_count == 0)
    {
        orbit_count = 0;
        orbits = NULL;
        return;
    }

    galois_factors = malloc((size_t)factor_count * sizeof(*galois_factors));
    if (galois_factors == NULL)
    {
        fprintf(stderr, "could not allocate Galois generator list\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < factor_count; ++i)
    {
        if (factors[i].degree > 1)
            galois_factors[galois_count++] = i;
    }

    generator_count = 2 + galois_count + swap_count;
    if ((size_t)generator_count > SIZE_MAX / (size_t)exceptional_count ||
        (size_t)generator_count * (size_t)exceptional_count >
            SIZE_MAX / sizeof(*generator_images))
    {
        fprintf(stderr, "automorphism table is too large\n");
        exit(EXIT_FAILURE);
    }

    generator_images = malloc((size_t)generator_count *
                              (size_t)exceptional_count *
                              sizeof(*generator_images));
    if (generator_images == NULL)
    {
        fprintf(stderr, "could not allocate automorphism table\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i < exceptional_count; ++i)
    {
        int offset = 0;
        h_code_t code;

        generator_images[(size_t)offset++ * exceptional_count + i] =
            one_minus_id[i];

        code = inverse_h_code(i);
        generator_images[(size_t)offset++ * exceptional_count + i] =
            exceptional_id_from_code(code, "inversion");

        for (g = 0; g < galois_count; ++g)
        {
            code = frobenius_code(i, galois_factors[g]);
            generator_images[(size_t)offset++ * exceptional_count + i] =
                exceptional_id_from_code(code, "local Frobenius");
        }

        for (g = 0; g < swap_count; ++g)
        {
            code = swap_code(i, swaps[g].left, swaps[g].right);
            generator_images[(size_t)offset++ * exceptional_count + i] =
                exceptional_id_from_code(code, "equal-factor permutation");
        }
    }

    if (!bitset_init(&unclassified, (size_t)exceptional_count))
    {
        fprintf(stderr, "could not allocate orbit-classification bitset\n");
        exit(EXIT_FAILURE);
    }
    for (i = 0; i < exceptional_count; ++i)
        bitset_add(&unclassified, (size_t)i);

    orbit_buffer = malloc((size_t)exceptional_count * sizeof(*orbit_buffer));
    orbits = malloc(orbit_capacity * sizeof(*orbits));
    if ((orbit_buffer == NULL || orbits == NULL) && exceptional_count != 0)
    {
        fprintf(stderr, "orbit allocation failed\n");
        exit(EXIT_FAILURE);
    }

    orbit_count = 0;
    for (i = 0; i < exceptional_count; ++i)
    {
        int head;
        int size;
        int *elements;

        if (!bitset_contains_unchecked(&unclassified, (size_t)i))
            continue;

        head = 0;
        size = 0;
        orbit_buffer[size++] = i;
        bitset_remove(&unclassified, (size_t)i);

        while (head < size)
        {
            int current = orbit_buffer[head++];

            for (g = 0; g < generator_count; ++g)
            {
                int image = generator_images[
                    (size_t)g * exceptional_count + current];

                if (bitset_contains_unchecked(&unclassified,
                                              (size_t)image))
                {
                    bitset_remove(&unclassified, (size_t)image);
                    orbit_buffer[size++] = image;
                }
            }
        }

        qsort(orbit_buffer, (size_t)size, sizeof(*orbit_buffer), compare_ints);

        if ((size_t)orbit_count == orbit_capacity)
        {
            orbit_t *grown;
            orbit_capacity *= 2;
            grown = realloc(orbits, orbit_capacity * sizeof(*orbits));
            if (grown == NULL)
            {
                fprintf(stderr, "orbit-list realloc failed\n");
                exit(EXIT_FAILURE);
            }
            orbits = grown;
        }

        elements = malloc((size_t)size * sizeof(*elements));
        if (elements == NULL)
        {
            fprintf(stderr, "orbit allocation failed\n");
            exit(EXIT_FAILURE);
        }
        memcpy(elements, orbit_buffer, (size_t)size * sizeof(*elements));

        orbits[orbit_count].size = size;
        orbits[orbit_count].elements = elements;
        ++orbit_count;
    }

    free(orbit_buffer);
    bitset_destroy(&unclassified);
    free(generator_images);
    free(galois_factors);
}

/* -------------------------------------------------------------------------
 * Program setup and cleanup
 * ------------------------------------------------------------------------- */

static void free_global_data(void)
{
    int i;

    bitset_destroy(&exceptional_codes);

    free(exceptional_exponents);
    exceptional_exponents = NULL;
    free(exceptional_id_by_code);
    exceptional_id_by_code = NULL;
    free(one_minus_id);
    one_minus_id = NULL;

    if (orbits != NULL)
    {
        for (i = 0; i < orbit_count; ++i)
            free(orbits[i].elements);
    }
    free(orbits);
    orbits = NULL;

    free(base_field_log);
    base_field_log = NULL;
    free(factors);
    factors = NULL;
}

int main(int argc, char *argv[])
{
    int degree_argument_end;
    int cpus = 1;
    int pari_stack_mb = 100;
    int option_index;
    int model_count = 0;
    field_model_t *models;
    swap_generator_t *swaps;
    int swap_count;
    FILE *result_stream;
    char result_path[MAXPATHLEN];
    struct timespec total_start;
    struct timespec preparation_start;
    double preparation_time;
    double total_time;
    int i;

    if (argc == 2 &&
        (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0))
    {
        fputs(USAGE, stdout);
        return 0;
    }

    if (argc < 3)
    {
        fputs(USAGE, stderr);
        return 1;
    }

    /* Parse the prime, factor degrees, and optional runtime parameters. */
    p = parse_positive_int(argv[1], "prime");
    if (p < 3 || (p & 1) == 0 || !is_prime_int(p))
    {
        fprintf(stderr, "p must be an odd prime\n");
        return 1;
    }

    degree_argument_end = 2;
    while (degree_argument_end < argc &&
           argv[degree_argument_end][0] != '-')
        ++degree_argument_end;

    factor_count = degree_argument_end - 2;
    if (factor_count < 1)
    {
        fputs(USAGE, stderr);
        return 1;
    }

    factors = calloc((size_t)factor_count, sizeof(*factors));
    if (factors == NULL)
    {
        fprintf(stderr, "could not allocate factor list\n");
        return 1;
    }

    total_dimension = 0;
    for (i = 0; i < factor_count; ++i)
    {
        factors[i].degree = parse_positive_int(argv[i + 2], "factor degree");
        if (total_dimension > INT_MAX - factors[i].degree)
        {
            fprintf(stderr, "total dimension is too large\n");
            return 1;
        }
        total_dimension += factors[i].degree;
    }

    if ((argc - degree_argument_end) % 2 != 0)
    {
        fputs(USAGE, stderr);
        return 1;
    }

    option_index = degree_argument_end;
    while (option_index < argc)
    {
        if (strcmp(argv[option_index], "-cpus") == 0)
            cpus = parse_positive_int(argv[option_index + 1], "CPU count");
        else if (strcmp(argv[option_index], "-pari") == 0)
            pari_stack_mb = parse_positive_int(argv[option_index + 1],
                                               "PARI stack size");
        else if (strcmp(argv[option_index], "-v4-threshold") == 0)
            v4_threshold = parse_nonnegative_int(argv[option_index + 1],
                                                 "V4 threshold");
        else
        {
            fputs(USAGE, stderr);
            return 1;
        }
        option_index += 2;
    }

    /* Create a deterministic output directory before initializing PARI. */
    if (getcwd(computation_path, sizeof(computation_path)) == NULL)
    {
        perror("getcwd");
        return 1;
    }
    append_text(computation_path, sizeof(computation_path), "/p=%d_d=", p);
    for (i = 0; i < factor_count; ++i)
        append_text(computation_path, sizeof(computation_path),
                    "%s%d", i == 0 ? "" : "x", factors[i].degree);

    ensure_directory(computation_path);

    checked_snprintf(local_path, sizeof(local_path), "%s/local",
                     computation_path);
    ensure_directory(local_path);
    clear_previous_local_results();

    checked_snprintf(error_path, sizeof(error_path), "%s/stderr",
                     computation_path);
    if (freopen(error_path, "w", stderr) == NULL)
    {
        perror("freopen");
        return 1;
    }

    checked_snprintf(result_path, sizeof(result_path), "%s/result.out",
                     computation_path);
    result_stream = fopen(result_path, "w");
    if (result_stream == NULL)
    {
        fprintf(stderr, "could not create result.out\n");
        return 1;
    }
    write_preparation_status(result_stream, "initializing PARI");

    total_start = monotonic_now();
    preparation_start = total_start;

    if ((size_t)pari_stack_mb > SIZE_MAX / 1000000u)
    {
        fprintf(stderr, "PARI stack size is too large\n");
        return 1;
    }
    pari_init_opts((size_t)pari_stack_mb * 1000000u,
                   2, INIT_JMPm | INIT_DFTm);

    models = calloc((size_t)factor_count, sizeof(*models));
    if (models == NULL)
    {
        fprintf(stderr, "could not allocate field models\n");
        return 1;
    }

    write_preparation_status(result_stream, "initializing field factors");

    initialize_base_field_log();

    /* Initialize one field model for every distinct factor degree. */
    for (i = 0; i < factor_count; ++i)
    {
        int model_index = find_model(models, model_count, factors[i].degree);

        if (model_index < 0)
        {
            char stage[128];

            checked_snprintf(stage, sizeof(stage),
                             "initializing field degree %d",
                             factors[i].degree);
            write_preparation_status(result_stream, stage);
            model_index = model_count++;
            initialize_field_model(&models[model_index],
                                   factors[i].degree, result_stream,
                                   factor_count > 1);
        }

        factors[i].model_index = model_index;
        factors[i].q = models[model_index].q;
        factors[i].order = models[model_index].order;
    }

    /*
     * All primitive field generators have norm zeta, so for
     * x=(g_i^{e_i}) the condition N(x)=+-1 is simply
     *
     *     sum_i e_i = 0 mod h,   h=(p-1)/2.
     *
     * Use a smallest-degree factor as the mixed-radix pivot.  Its exponent is
     * the only dependent coordinate, so this keeps the reconstructed modulus
     * and the rapidly varying pivot digit as small as possible; exponents in
     * larger factors remain direct digits.  In case of a tie, keep the first
     * occurrence in the user's factor order.
     */
    norm_half_order = ((uint64_t)p - 1) / 2;
    pivot_factor = 0;
    for (i = 1; i < factor_count; ++i)
    {
        if (factors[i].degree < factors[pivot_factor].degree)
            pivot_factor = i;
    }

    if (norm_half_order == 0 ||
        factors[pivot_factor].order % norm_half_order != 0)
    {
        fprintf(stderr, "invalid norm-+-1 pivot radix\n");
        return 1;
    }

    /*
     * Give the pivot digit stride 1, then append the remaining factors in
     * their original order.  The radices are n_pivot/h for the pivot and n_i
     * for every other factor, so their product is exactly |H|.
     */
    total_h_elements = 1;
    factors[pivot_factor].h_radix =
        factors[pivot_factor].order / norm_half_order;
    factors[pivot_factor].h_stride = total_h_elements;
    total_h_elements *= factors[pivot_factor].h_radix;

    for (i = 0; i < factor_count; ++i)
    {
        uint64_t radix;

        if (i == pivot_factor)
            continue;

        radix = factors[i].order;
        factors[i].h_radix = radix;
        factors[i].h_stride = total_h_elements;

        if (radix != 0 && total_h_elements > UINT64_MAX / radix)
        {
            fprintf(stderr,
                    "norm-+-1 group exceeds 64-bit mixed-radix range\n");
            return 1;
        }
        total_h_elements *= radix;
    }

    if (total_h_elements > SIZE_MAX)
    {
        fprintf(stderr,
                "norm-+-1 group exceeds addressable memory range\n");
        return 1;
    }

    write_preparation_status(result_stream, "enumerating exceptional units");

    if (!bitset_init(&exceptional_codes, (size_t)total_h_elements))
    {
        fprintf(stderr, "could not allocate exceptional-unit bitset\n");
        return 1;
    }

    /* Enumerate exceptional units and build constant-time lookup tables. */
    {
        h_code_t *one_minus_code;
        h_code_t *exceptional_code;
        factor_exp_t *exponents;
        factor_exp_t *mate_exponents;
        h_code_t code;
        int count = 0;

        if ((size_t)total_h_elements >
            SIZE_MAX / sizeof(*one_minus_code))
        {
            fprintf(stderr, "one-minus table is too large\n");
            return 1;
        }

        one_minus_code = malloc((size_t)total_h_elements *
                                sizeof(*one_minus_code));
        exponents = factor_count > 1
                        ? malloc((size_t)factor_count * sizeof(*exponents))
                        : NULL;
        mate_exponents = factor_count > 1
                             ? malloc((size_t)factor_count *
                                      sizeof(*mate_exponents))
                             : NULL;
        if (one_minus_code == NULL ||
            (factor_count > 1 &&
             (exponents == NULL || mate_exponents == NULL)))
        {
            fprintf(stderr, "exceptional enumeration allocation failed\n");
            return 1;
        }

        for (code = 0; code < total_h_elements; ++code)
            one_minus_code[(size_t)code] = UNSET_CODE;

        if (factor_count == 1)
        {
            count = enumerate_single_field_exceptional(
                &models[factors[0].model_index], one_minus_code,
                result_stream);
        }
        else
        {
            for (code = 0; code < total_h_elements; ++code)
            {
                int is_exceptional;
                h_code_t mate;

                if (one_minus_code[(size_t)code] != UNSET_CODE)
                    continue;

                decode_h_code(code, exponents);
                mate = compute_one_minus_code(
                    models, exponents, mate_exponents, &is_exceptional);

                if (!is_exceptional)
                    continue;

                record_exceptional_pair(one_minus_code, code, mate,
                                        &count);
            }
        }

        exceptional_count = count;

        if ((size_t)exceptional_count >
                SIZE_MAX / (size_t)factor_count ||
            (size_t)exceptional_count * (size_t)factor_count >
                SIZE_MAX / sizeof(*exceptional_exponents) ||
            (size_t)total_h_elements >
                SIZE_MAX / sizeof(*exceptional_id_by_code))
        {
            fprintf(stderr, "exceptional tables are too large\n");
            return 1;
        }

        exceptional_code = malloc((size_t)exceptional_count *
                                  sizeof(*exceptional_code));
        exceptional_exponents = malloc((size_t)exceptional_count *
                                       (size_t)factor_count *
                                       sizeof(*exceptional_exponents));
        exceptional_id_by_code = malloc((size_t)total_h_elements *
                                        sizeof(*exceptional_id_by_code));
        one_minus_id = malloc((size_t)exceptional_count *
                              sizeof(*one_minus_id));
        if ((exceptional_code == NULL || exceptional_exponents == NULL ||
             one_minus_id == NULL) && exceptional_count != 0)
        {
            fprintf(stderr, "exceptional data allocation failed\n");
            return 1;
        }
        if (exceptional_id_by_code == NULL && total_h_elements != 0)
        {
            fprintf(stderr, "H-code-to-id allocation failed\n");
            return 1;
        }

        for (code = 0; code < total_h_elements; ++code)
            exceptional_id_by_code[(size_t)code] = -1;

        count = 0;
        for (code = 0; code < total_h_elements; ++code)
        {
            if (!bitset_contains_unchecked(&exceptional_codes, (size_t)code))
                continue;

            exceptional_code[count] = code;
            decode_h_code(code,
                        exceptional_exponents +
                        (size_t)count * (size_t)factor_count);
            exceptional_id_by_code[(size_t)code] = count;
            ++count;
        }

        for (i = 0; i < exceptional_count; ++i)
        {
            h_code_t mate = one_minus_code[(size_t)exceptional_code[i]];
            int mate_id;

            if (mate == UNSET_CODE || mate >= total_h_elements)
            {
                fprintf(stderr, "missing one-minus H code\n");
                return 1;
            }
            mate_id = exceptional_id_by_code[(size_t)mate];
            if (mate_id < 0)
            {
                fprintf(stderr, "one-minus image is not exceptional\n");
                return 1;
            }
            one_minus_id[i] = mate_id;
        }

        free(exceptional_code);
        free(mate_exponents);
        free(exponents);
        free(one_minus_code);
    }

    /*
     * No later phase performs finite-field arithmetic.  Release the base-field
     * logarithm table, the per-degree 1-x tables, and the PARI stack before
     * allocating the symmetry-generator table used by orbit classification.
     * This avoids overlapping the two largest preparation-time memory peaks.
     */
    free(base_field_log);
    base_field_log = NULL;
    free_field_models(models, model_count);
    models = NULL;
    pari_close();

    write_preparation_status(result_stream, "classifying symmetry orbits");

    swaps = build_swap_generators(&swap_count);
    classify_outer_orbits(swaps, swap_count);
    free(swaps);

    /*
     * Direct clique searches need only exceptional_codes, exceptional_exponents
     * and the orbit lists.  The code-to-ID and 1-x tables are retained only
     * when V4 may be requested, because V4 transformations use them.
     */
    if (v4_threshold == 0)
    {
        free(exceptional_id_by_code);
        exceptional_id_by_code = NULL;
        free(one_minus_id);
        one_minus_id = NULL;
    }

    preparation_time = elapsed_seconds(preparation_start);

    if (fclose(result_stream) != 0)
    {
        fprintf(stderr, "could not close preparation status file\n");
        return 1;
    }
    result_stream = fopen(result_path, "w");
    if (result_stream == NULL)
    {
        fprintf(stderr, "could not rewrite result.out\n");
        return 1;
    }

    fprintf(result_stream, "Program: %s\n", PROGRAM_NAME);
    fprintf(result_stream, "Factors:");
    for (i = 0; i < factor_count; ++i)
        fprintf(result_stream, " %d", factors[i].degree);
    fputc('\n', result_stream);
    fprintf(result_stream, "Dimension: %d\n", total_dimension);
    fprintf(result_stream, "CPUs: %d\n", cpus);
    fprintf(result_stream, "PARI stack: %d MB\n", pari_stack_mb);
    fprintf(result_stream, "Preparation mode: %s\n",
            factor_count == 1 ? "single-field subgroup"
                              : "multifactor local tables");
    fprintf(result_stream, "Norm +-1 units: %" PRIu64 "\n",
            total_h_elements);
    fprintf(result_stream, "Exceptional units: %d\n", exceptional_count);
    fprintf(result_stream, "Orbits: %d\n", orbit_count);
    fprintf(result_stream, "Preparation: %.*f s\n",
            TIME_PRECISION, preparation_time);
    fprintf(result_stream, "V4 threshold: %d%s\n", v4_threshold,
            v4_threshold == 0 ? " (disabled)" : "");
    fprintf(result_stream, "Status: computing\n\n");
    if (fclose(result_stream) != 0)
    {
        fprintf(stderr, "could not close result header\n");
        return 1;
    }

    {
        int completed_orbits = 0;
        int clique_number = compute_all_orbits(result_path, cpus,
                                                &completed_orbits);

        total_time = elapsed_seconds(total_start);
        result_stream = fopen(result_path, "a");
        if (result_stream == NULL)
        {
            fprintf(stderr, "could not finalize result.out\n");
            return 1;
        }
        fprintf(result_stream,
                "\nOrbit computations: %d/%d successfully completed\n",
                completed_orbits, orbit_count);
        fprintf(result_stream, "Clique number: %d\n", clique_number);
        fprintf(result_stream, "Total time: %.*f s\n",
                TIME_PRECISION, total_time);
        fprintf(result_stream, "Status: complete\n");
        if (fclose(result_stream) != 0)
        {
            fprintf(stderr, "could not close final result.out\n");
            return 1;
        }
    }

    free_global_data();
    return 0;
}
