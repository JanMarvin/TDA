/* tda_rtfmt.c -- printing through the formats makefmt()/makenfmt() build
   at run time ("%10.4f ", "%012.4e\t", "(%g) ", "%6d ").  The conversion
   goes through a literal format and the width and flags are applied per
   C99, so -Wformat=2 has nothing to flag; tools/check_rtfmt.c holds the
   output byte-identical to the C library's.  Only the first conversion is
   performed, the rest is copied as text. */

#include "tda.h"
#include "t_gen.h"
#include <limits.h>

typedef struct {
    const char *pre;
    const char *pre_end;
    const char *post;
    int minus, plus, space, hash, zero;
    int width;
    int prec;
    int lmod;
    char conv;
} RtSpec;

enum { RT_LMOD_NONE, RT_LMOD_H, RT_LMOD_HH, RT_LMOD_OTHER };

static int rt_digits(const char **pp)
{
    const char *p = *pp;
    long v = 0;

    while (*p >= '0' && *p <= '9') {
        if (v < INT_MAX / 10)
            v = v * 10 + (*p - '0');
        ++p;
    }
    *pp = p;
    return (int)v;
}

static void rt_parse(const char *fmt, RtSpec *s)
{
    const char *p = fmt;
    const char *spec;

    memset(s, 0, sizeof(*s));
    s->pre = fmt;
    s->width = -1;
    s->prec = -1;

    while (*p) {
        if (p[0] == '%' && p[1] == '%') {
            p += 2;
            continue;
        }
        if (p[0] == '%')
            break;
        ++p;
    }
    s->pre_end = p;
    if (!*p) {
        s->post = p;
        return;
    }

    spec = p++;
    for (;; ++p) {
        if (*p == '-')
            s->minus = 1;
        else if (*p == '+')
            s->plus = 1;
        else if (*p == ' ')
            s->space = 1;
        else if (*p == '#')
            s->hash = 1;
        else if (*p == '0')
            s->zero = 1;
        else
            break;
    }
    if (*p == '*')
        ++p;
    else if (*p >= '0' && *p <= '9')
        s->width = rt_digits(&p);
    if (*p == '.') {
        ++p;
        if (*p == '*')
            ++p;
        else
            s->prec = rt_digits(&p);
    }
    if (*p == 'h') {
        ++p;
        s->lmod = RT_LMOD_H;
        if (*p == 'h') {
            ++p;
            s->lmod = RT_LMOD_HH;
        }
    }
    else
        while (*p && strchr("lLqjzt", *p)) {
            ++p;
            s->lmod = RT_LMOD_OTHER;
        }

    if (*p && strchr("fFeEgGaAdiouxXcs", *p)) {
        s->conv = *p;
        s->post = p + 1;
    }
    else
        s->post = spec;
}

static int rt_is_dbl(char c)
{
    return c != 0 && strchr("fFeEgGaA", c) != NULL;
}

static int rt_dbl(char *b, size_t n, char conv, int hash, int prec, double x)
{
    switch (conv) {
        case 'F': return hash ? snprintf(b, n, "%#.*F", prec, x) : snprintf(b, n, "%.*F", prec, x);
        case 'e': return hash ? snprintf(b, n, "%#.*e", prec, x) : snprintf(b, n, "%.*e", prec, x);
        case 'E': return hash ? snprintf(b, n, "%#.*E", prec, x) : snprintf(b, n, "%.*E", prec, x);
        case 'g': return hash ? snprintf(b, n, "%#.*g", prec, x) : snprintf(b, n, "%.*g", prec, x);
        case 'G': return hash ? snprintf(b, n, "%#.*G", prec, x) : snprintf(b, n, "%.*G", prec, x);
        case 'a': return hash ? snprintf(b, n, "%#.*a", prec, x) : snprintf(b, n, "%.*a", prec, x);
        case 'A': return hash ? snprintf(b, n, "%#.*A", prec, x) : snprintf(b, n, "%.*A", prec, x);
        default:  return hash ? snprintf(b, n, "%#.*f", prec, x) : snprintf(b, n, "%.*f", prec, x);
    }
}

static int rt_int(char *b, size_t n, const RtSpec *s, int x)
{
    unsigned int u = (unsigned int)x;

    if (s->lmod == RT_LMOD_H) {
        x = (short)x;
        u = (unsigned short)u;
    }
    else if (s->lmod == RT_LMOD_HH) {
        x = (signed char)x;
        u = (unsigned char)u;
    }
    switch (s->conv) {
        case 'o': return s->hash ? snprintf(b, n, "%#.*o", s->prec, u) : snprintf(b, n, "%.*o", s->prec, u);
        case 'u': return snprintf(b, n, "%.*u", s->prec, u);
        case 'x': return s->hash ? snprintf(b, n, "%#.*x", s->prec, u) : snprintf(b, n, "%.*x", s->prec, u);
        case 'X': return s->hash ? snprintf(b, n, "%#.*X", s->prec, u) : snprintf(b, n, "%.*X", s->prec, u);
        case 'c': return snprintf(b, n, "%c", (unsigned char)u);
        default:  return snprintf(b, n, "%.*d", s->prec, x);
    }
}

static void rt_put(char *out, size_t cap, size_t *pos, char c)
{
    if (*pos + 1 < cap)
        out[*pos] = c;
    ++*pos;
}

static void rt_text(char *out, size_t cap, size_t *pos, const char *p, const char *end)
{
    while (p < end && *p) {
        if (p[0] == '%' && p[1] == '%')
            ++p;
        rt_put(out, cap, pos, *p++);
    }
}

static int rt_render(char *out, size_t cap, const char *fmt, int isdbl, double dx, int ix)
{
    RtSpec s;
    char small[512];
    char *cv = small;
    const char *lead;
    size_t pos = 0, len, k, pad = 0;
    int n, sign, zero;

    rt_parse(fmt, &s);
    rt_text(out, cap, &pos, s.pre, s.pre_end);

    if (s.conv) {
        if (isdbl && !rt_is_dbl(s.conv)) {
            s.conv = 'f';
            s.prec = 0;
        }
        else if (!isdbl && rt_is_dbl(s.conv))
            dx = (double)ix;
        else if (!isdbl && s.conv == 's') {
            s.conv = 'd';
            s.prec = -1;
        }

        sign = rt_is_dbl(s.conv) || s.conv == 'd' || s.conv == 'i';
        lead = sign ? (s.plus ? "+" : s.space ? " " : "") : "";
        n = rt_is_dbl(s.conv) ? rt_dbl(small, sizeof(small), s.conv, s.hash, s.prec, dx)
                              : rt_int(small, sizeof(small), &s, ix);
        if (n < 0)
            n = 0;
        if ((size_t)n >= sizeof(small)) {
            cv = (char *)malloc((size_t)n + 1);
            if (cv == NULL) {
                cv = small;
                n = (int)sizeof(small) - 1;
            }
            else if (rt_is_dbl(s.conv))
                rt_dbl(cv, (size_t)n + 1, s.conv, s.hash, s.prec, dx);
            else
                rt_int(cv, (size_t)n + 1, &s, ix);
        }
        if (cv[0] == '-')
            lead = "";

        len = strlen(lead) + (size_t)n;
        if (s.width > 0 && (size_t)s.width > len)
            pad = (size_t)s.width - len;

        zero = s.zero && !s.minus && s.conv != 'c';
        if (rt_is_dbl(s.conv))
            zero = zero && isfinite(dx);
        else if (s.prec >= 0)
            zero = 0;

        if (pad && !s.minus && !zero)
            for (k = 0; k < pad; ++k)
                rt_put(out, cap, &pos, ' ');
        rt_text(out, cap, &pos, lead, lead + strlen(lead));
        k = 0;
        if (zero) {
            if (cv[0] == '-' || cv[0] == '+' || cv[0] == ' ')
                rt_put(out, cap, &pos, cv[k++]);
            if (cv[k] == '0' && (cv[k + 1] == 'x' || cv[k + 1] == 'X')) {
                rt_put(out, cap, &pos, cv[k++]);
                rt_put(out, cap, &pos, cv[k++]);
            }
            for (; pad; --pad)
                rt_put(out, cap, &pos, '0');
        }
        for (; k < (size_t)n; ++k)
            rt_put(out, cap, &pos, cv[k]);
        if (s.minus)
            for (; pad; --pad)
                rt_put(out, cap, &pos, ' ');
        if (cv != small)
            free(cv);
    }

    rt_text(out, cap, &pos, s.post, s.post + strlen(s.post));
    if (cap > 0)
        out[pos < cap ? pos : cap - 1] = '\0';
    return pos > (size_t)INT_MAX ? INT_MAX : (int)pos;
}

int rt_snprintf_d(char *buf, size_t sz, const char *fmt, double x)
{
    return rt_render(buf, sz, fmt, 1, x, 0);
}

int rt_snprintf_i(char *buf, size_t sz, const char *fmt, int x)
{
    return rt_render(buf, sz, fmt, 0, 0.0, x);
}

static char *rt_alloc_render(char *small, size_t sz, const char *fmt, int isdbl, double dx, int ix)
{
    char *big;
    int n = rt_render(small, sz, fmt, isdbl, dx, ix);

    if ((size_t)n < sz)
        return small;
    big = (char *)malloc((size_t)n + 1);
    if (big == NULL)
        return small;
    rt_render(big, (size_t)n + 1, fmt, isdbl, dx, ix);
    return big;
}

int rt_fprintf_d(TDAContext *ctx, FILE *fd, const char *fmt, double x)
{
    char small[512];
    char *t = rt_alloc_render(small, sizeof(small), fmt, 1, x, 0);
    int r;

#ifdef TDA_R_PACKAGE
    r = tda_fprintf_rt(ctx, fd, t, fmt, x);
#else
    (void)ctx;
    r = fprintf(fd, "%s", t);
#endif
    if (t != small)
        free(t);
    return r;
}

int rt_fprintf_i(TDAContext *ctx, FILE *fd, const char *fmt, int x)
{
    char small[512];
    char *t = rt_alloc_render(small, sizeof(small), fmt, 0, 0.0, x);
    int r;

#ifdef TDA_R_PACKAGE
    r = tda_fprintf_rt(ctx, fd, t, fmt, x);
#else
    (void)ctx;
    r = fprintf(fd, "%s", t);
#endif
    if (t != small)
        free(t);
    return r;
}

void rt_printf1_d(TDAContext *ctx, const char *fmt, double x)
{
    char small[512];
    char *t = rt_alloc_render(small, sizeof(small), fmt, 1, x, 0);

    printf1_rt(ctx, fmt, t, x);
    if (t != small)
        free(t);
}

void rt_printf1_i(TDAContext *ctx, const char *fmt, int x)
{
    char small[512];
    char *t = rt_alloc_render(small, sizeof(small), fmt, 0, 0.0, x);

    printf1_rt(ctx, fmt, t, x);
    if (t != small)
        free(t);
}

void rt_printf2_d(TDAContext *ctx, const char *fmt, double x)
{
    char small[512];
    char *t = rt_alloc_render(small, sizeof(small), fmt, 1, x, 0);

    printf2(ctx, "%s", t);
    if (t != small)
        free(t);
}
