/* Modified by Oleksandr Kulkov with AI assistance, 2026-10-31:
   preserve near-collinear vertices in reentrant 2-d convex hulls.
   Qhull license: COPYING.txt; original source: https://www.qhull.org/ */
#ifndef qhDEFexactsum
#define qhDEFexactsum
#if defined(_MSC_VER) && _MSC_VER < 1600
/* Visual Studio 2005/2008 do not provide stdint.h. */
typedef unsigned __int32 qh_uint32;
typedef unsigned __int64 qh_uint64;
#else
#include <stdint.h>
typedef uint32_t qh_uint32;
typedef uint64_t qh_uint64;
#endif
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/* Other floating formats continue to use the original Qhull path. */
#define HULL_BINARY64 (sizeof(double) == 8 && FLT_RADIX == 2 && \
                       DBL_MANT_DIG == 53 && DBL_MAX_EXP == 1024)

/* Preserve the platform hypot implementation without requiring C99 headers. */
#if defined(_MSC_VER) && _MSC_VER < 1800
#define hull_hypot _hypot
#else
#if !defined(__cplusplus) && (!defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L)
extern double hypot(double, double);
#endif
#define hull_hypot hypot
#endif
#define hull_isfinite(x) ((x) <= DBL_MAX && (x) >= -DBL_MAX)
/* The sign arguments are nonzero determinants, never signed zero or NaN. */
#define hull_copysign(x,y) ((y) < 0 ? -(x) : (x))
#define HULL_TRUE_MIN (DBL_MIN * DBL_EPSILON)

/* The arguments used here are non-negative finite magnitudes. */
#define hull_fmax(a,b) ((a) > (b) ? (a) : (b))

/* A binary64 product has exponents from -2148 through 2047. The
 * additional words allow summing up to INT_MAX products without overflow.
 * Positive and negative terms are accumulated separately before subtraction.
 */
#define HULL_WORDS 136
/* [first,last) bounds occupied words; last == 0 denotes an empty sum. */
typedef struct { qh_uint32 pos[HULL_WORDS], neg[HULL_WORDS]; int first, last; } hull_sum;
static int hull_add_word(qh_uint32 *a, int i, qh_uint64 v) {
    while (v) {
        v += a[i];
        a[i++] = (qh_uint32)v;
        v >>= 32;
    }
    return i;
}
static qh_uint64 hull_decode(double x, int *exponent, int *negative) {
    qh_uint64 bits, m;
    int e;
    memcpy(&bits, &x, sizeof(bits));
    *negative = (int)(bits >> 63);
    e = (int)((bits >> 52) & 2047);
    m = bits & (((qh_uint64)0xfffff << 32) | 0xffffffffU);
    *exponent = e ? e - 1075 : -1074;
    return e ? m | (((qh_uint64)1) << 52) : m;
}
static void hull_product(hull_sum *s, double x, double y, int minus) {
    int ex, ey, sx, sy, i, j, shift;
    qh_uint32 product[4] = {0};
    qh_uint32 *a;
    qh_uint64 mx = hull_decode(x, &ex, &sx);
    qh_uint64 my = hull_decode(y, &ey, &sy);
    if (!mx || !my) return;
    for (i = 0; i < 2; ++i) {
        qh_uint64 carry = 0;
        for (j = 0; j < 2; ++j) {
            qh_uint64 v = (qh_uint64)(qh_uint32)(mx >> (32*i)) *
                         (qh_uint32)(my >> (32*j)) + product[i+j] + carry;
            product[i+j] = (qh_uint32)v;
            carry = v >> 32;
        }
        product[i+2] = (qh_uint32)carry;
    }
    shift = ex + ey + 2148;
    if (!s->last || shift/32 < s->first) s->first=shift/32;
    a = (sx ^ sy ^ minus) ? s->neg : s->pos;
    for (i = 0; i < 4; ++i) {
        int end=hull_add_word(a, shift/32+i, (qh_uint64)product[i] << (shift%32));
        if (end>s->last) s->last=end;
    }
}
static int hull_sum_sign(const hull_sum *s) {
    int i;
    for (i = s->last - 1; i >= s->first; --i) {
        if (s->pos[i] != s->neg[i]) return s->pos[i] > s->neg[i] ? 1 : -1;
    }
    return 0;
}
static double hull_sum_divide(const hull_sum *s, int scale, double divisor) {
    int sign = hull_sum_sign(s);
    const qh_uint32 *a = sign > 0 ? s->pos : s->neg;
    const qh_uint32 *b = sign > 0 ? s->neg : s->pos;
    qh_uint32 d[HULL_WORDS];
    qh_uint64 borrow = 0;
    int top = s->first, bottom, i;
    double value = 0;
    if (!sign) return 0;
    for (i = s->first; i < s->last; ++i) {
        qh_uint64 v = (qh_uint64)b[i] + borrow;
        d[i] = (qh_uint32)((qh_uint64)a[i] - v);
        borrow = (qh_uint64)a[i] < v;
        if (d[i]) top = i;
    }
    bottom = top > s->first+2 ? top - 2 : s->first;
    for (i = top; i >= bottom; --i) value = ldexp(value, 32) + d[i];
    return sign * ldexp(value/divisor, 32*bottom - 2148 + scale);
}

static double hull_sum_value(const hull_sum *s, int scale) {
    return hull_sum_divide(s, scale, 1);
}
/* Bound the edge components for safe norm evaluation. Return the binary
 * scale exponent; ordinary edges need no rescaling. */
static int hull_edge(const coordT *a, const coordT *b, double *dx, double *dy) {
    int exponent, extra=0;
    *dx=b[0]-a[0];
    *dy=b[1]-a[1];
    if (!hull_isfinite(*dx) || !hull_isfinite(*dy)) {
        *dx=b[0]/2-a[0]/2;
        *dy=b[1]/2-a[1]/2;
        extra=1;
    }
    if (!extra && hull_fmax(fabs(*dx),fabs(*dy)) >= 1e-100 &&
        hull_fmax(fabs(*dx),fabs(*dy)) <= 1e100)
        return 0;
    frexp(hull_fmax(fabs(*dx),fabs(*dy)),&exponent);
    *dx=ldexp(*dx,-exponent);
    *dy=ldexp(*dy,-exponent);
    return exponent+extra;
}
#endif
