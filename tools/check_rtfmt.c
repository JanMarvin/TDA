/* Checks tda_rtfmt.c against the C library: every format shape makefmt()
   and makenfmt() can build, plus random flag/width/precision combinations,
   rendered by rt_snprintf_*() and by snprintf() itself.  The oracle side
   needs a run-time format, so this file alone is built without
   -Wformat-nonliteral:

       sh tools/check_rtfmt.sh
*/
#include "../src/tda_rtfmt.c"
#include <float.h>

void printf1_rt(TDAContext *ctx, const char *fmt, const char *text, ...)
{
    (void)ctx; (void)fmt; (void)text;
}

void printf2(TDAContext *ctx, const char *fmt, ...)
{
    (void)ctx; (void)fmt;
}

static unsigned long rng = 12345UL;

static unsigned int rnd(unsigned int n)
{
    rng = rng * 1103515245UL + 12345UL;
    return (unsigned int)((rng >> 8) % n);
}

static long nfail = 0, ncase = 0;

static void cmp_d(const char *fmt, double x)
{
    char a[2048], b[2048], c[8], d[8];
    int na, nb, nc, nd;

    na = rt_snprintf_d(a, sizeof(a), fmt, x);
    nb = snprintf(b, sizeof(b), fmt, x);
    nc = rt_snprintf_d(c, sizeof(c), fmt, x);
    nd = snprintf(d, sizeof(d), fmt, x);
    ++ncase;
    if (na != nb || nc != nd || memcmp(a, b, (size_t)na + 1) != 0 ||
        memcmp(c, d, (size_t)(nc < 7 ? nc : 7) + 1) != 0) {
        if (++nfail <= 20)
            printf("FAIL d fmt=\"%s\" x=%.17g: rt=\"%s\"(%d) libc=\"%s\"(%d)\n",
                   fmt, x, a, na, b, nb);
    }
}

static void cmp_i(const char *fmt, int x)
{
    char a[2048], b[2048], c[8], d[8];
    int na, nb, nc, nd;

    na = rt_snprintf_i(a, sizeof(a), fmt, x);
    nb = snprintf(b, sizeof(b), fmt, x);
    nc = rt_snprintf_i(c, sizeof(c), fmt, x);
    nd = snprintf(d, sizeof(d), fmt, x);
    ++ncase;
    if (na != nb || nc != nd || memcmp(a, b, (size_t)na + 1) != 0 ||
        memcmp(c, d, (size_t)(nc < 7 ? nc : 7) + 1) != 0) {
        if (++nfail <= 20)
            printf("FAIL i fmt=\"%s\" x=%d: rt=\"%s\"(%d) libc=\"%s\"(%d)\n",
                   fmt, x, a, na, b, nb);
    }
}

static const double dv[] = {
    0.0, -0.0, 1.0, -1.0, 0.5, -0.5, 0.05, 123.456, -123.456, 1e-5, 9.99995,
    0.000123456789, 1e15, -1e15, 1e-300, 1e300, 5e-324, DBL_MAX, -DBL_MAX,
    DBL_MIN, 2.5, 3.5, -2.5, 1234567.891, 0.1 + 0.2
};
static const int iv[] = { 0, 1, -1, 7, -7, 42, 99999, -99999, 123456789, INT_MAX, INT_MIN };

static double special(int k)
{
    volatile double z = 0.0;
    return k == 0 ? 1.0 / z : k == 1 ? -1.0 / z : z / z;
}

static void all_d(const char *fmt)
{
    size_t i;
    int k;

    for (i = 0; i < sizeof(dv) / sizeof(dv[0]); ++i)
        cmp_d(fmt, dv[i]);
    for (k = 0; k < 3; ++k)
        cmp_d(fmt, special(k));
    for (k = 0; k < 20; ++k)
        cmp_d(fmt, ((double)rnd(2000000) - 1e6) * pow(10.0, (double)rnd(40) - 20.0));
}

static void all_i(const char *fmt)
{
    size_t i;
    int k;

    for (i = 0; i < sizeof(iv) / sizeof(iv[0]); ++i)
        cmp_i(fmt, iv[i]);
    for (k = 0; k < 10; ++k)
        cmp_i(fmt, (int)rnd(2000000) - 1000000);
}

int main(void)
{
    static const char seps[] = { ' ', '\t', ',', ';', '\n', '\0' };
    static const char dconv[] = "fFeEgGaA";
    static const char iconv[] = "diouxXc";
    char fmt[64], flags[8];
    size_t s;
    int w, p, k, n;

    for (s = 0; s < sizeof(seps); ++s) {
        for (w = 1; w < 100; w += (w < 25 ? 1 : 7))
            for (p = 0; p <= w && p < 100; p += (p < 12 ? 1 : 9)) {
                snprintf(fmt, sizeof(fmt), "%%%d.%df%c", w, p, seps[s]);  all_d(fmt);
                snprintf(fmt, sizeof(fmt), "%%0%d.%df%c", w, p, seps[s]); all_d(fmt);
                snprintf(fmt, sizeof(fmt), "(%%%d.%df)%c", w, p, seps[s]); all_d(fmt);
                snprintf(fmt, sizeof(fmt), "%%%d.%de%c", w + 8, p, seps[s]);  all_d(fmt);
                snprintf(fmt, sizeof(fmt), "%%0%d.%de%c", w + 8, p, seps[s]); all_d(fmt);
                snprintf(fmt, sizeof(fmt), "(%%%d.%de)%c", w + 8, p, seps[s]); all_d(fmt);
            }
        snprintf(fmt, sizeof(fmt), "%%lg%c", seps[s]);  all_d(fmt);
        snprintf(fmt, sizeof(fmt), "%%0lg%c", seps[s]); all_d(fmt);
        snprintf(fmt, sizeof(fmt), "(%%g)%c", seps[s]); all_d(fmt);
        for (w = 1; w < 100; ++w) {
            snprintf(fmt, sizeof(fmt), "%%%dd%c", w, seps[s]);
            all_i(fmt);
        }
    }

    for (k = 0; k < 20000; ++k) {
        const char *fl = "-+ #0";
        size_t nf = 0;
        int isd = (int)rnd(2);
        char conv = isd ? dconv[rnd(8)] : iconv[rnd(7)];

        for (n = 0; n < 5; ++n)
            if (rnd(3) == 0)
                flags[nf++] = fl[n];
        flags[nf] = '\0';
        if (conv == 'c' || conv == 'd' || conv == 'i' || conv == 'u')
            for (n = 0; flags[n]; ++n)
                if (flags[n] == '#' || (conv == 'c' && flags[n] == '0'))
                    flags[n] = '-';
        w = (int)rnd(40) - 5;
        p = (int)rnd(30) - 5;
        n = snprintf(fmt, sizeof(fmt), "%s%%%s", rnd(4) ? "" : "a%%b", flags);
        if (w > 0)
            n += snprintf(fmt + n, sizeof(fmt) - (size_t)n, "%d", w);
        if (p >= 0 && conv != 'c')
            n += snprintf(fmt + n, sizeof(fmt) - (size_t)n, ".%d", p);
        else if (p == -1 && conv != 'c')
            n += snprintf(fmt + n, sizeof(fmt) - (size_t)n, ".");
        if (isd && rnd(3) == 0)
            n += snprintf(fmt + n, sizeof(fmt) - (size_t)n, "l");
        if (!isd && conv != 'c' && rnd(4) == 0)
            n += snprintf(fmt + n, sizeof(fmt) - (size_t)n, rnd(2) ? "h" : "hh");
        snprintf(fmt + n, sizeof(fmt) - (size_t)n, "%c%s", conv, rnd(3) ? " " : "%%x");
        if (isd)
            all_d(fmt);
        else
            all_i(fmt);
    }

    printf("rtfmt: %ld cases, %ld failures\n", ncase, nfail);
    return nfail != 0;
}
