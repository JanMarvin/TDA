/* tda_out.c                                                                */
/*                                                                          */
/*  TDA's console primitives.  Outside the R package they are printf,       */
/*  vprintf and fprintf(stderr, ...).                                       */
/*                                                                          */
/*  Inside the package they must not name printf, puts, stdout or stderr -- */
/*  R CMD check refuses a shared object that does -- but they must still    */
/*  write to descriptors 1 and 2.  That is not a detail: tda_run() captures */
/*  a run by pointing those descriptors at a file and reading it back, so   */
/*  TDA's output is the result the R side parses, not messages for the      */
/*  user, and it has to land in that file.  Sending it to Rprintf() instead */
/*  would reach R's console rather than the capture, and would only look    */
/*  correct under Rscript, where R's console happens to be descriptor 1.    */
/*                                                                          */
/*  Nothing escapes to the user's terminal either way: tda_run() redirects  */
/*  both descriptors for the whole of every run.                            */

#include "tda.h"

#ifdef TDA_R_PACKAGE

#ifdef TDA_R_PACKAGE
/* Console output line counter, so the print tap (t_gen.c) can tag each
   captured number with the output line it appears on.  Every console
   write goes through one of the two functions below, including the ones
   that bypass printf1(), which is why the count lives here rather than
   in the tap.  Newlines are counted in the format string: a "%s" whose
   own argument contains one would be missed, and TDA has none. */
static long tda_out_lines = 0;
static long tda_out_line0 = 0;   /* line the current write started on */

/* Counting newlines in the format string is not enough: TDA echoes the
   command file with printf1("%s", line), so the newline arrives in the
   ARGUMENT, and the count drifts against the real output exactly where
   that happens.  The text is formatted a second time into a scratch
   buffer and the newlines counted there, which is what the reader
   actually sees. */
static void tda_count_lines(const char *fmt, va_list ap) TDA_PRINTF(1, 0);

static void tda_count_lines(const char *fmt, va_list ap)
{
    char buf[1024];
    char *big = NULL;
    const char *q;
    int n;
    va_list aq;

    if (fmt == NULL)
        return;
    tda_out_line0 = tda_out_lines;
    va_copy(aq, ap);
    n = vsnprintf(buf, sizeof(buf), fmt, aq);
    va_end(aq);
    if (n < 0)
        return;
    if ((size_t)n >= sizeof(buf)) {
        big = (char *)malloc((size_t)n + 1);
        if (big != NULL) {
            va_copy(aq, ap);
            vsnprintf(big, (size_t)n + 1, fmt, aq);
            va_end(aq);
        }
    }
    for (q = big != NULL ? big : buf; *q; ++q)
        if (*q == '\n')
            tda_out_lines++;
    free(big);
}

/* The line the write currently being counted began on -- the print tap
   needs this, not the line it ended on, and cannot work it out itself
   once a newline can arrive inside an argument. */
long tda_out_line_start(void)
{
    return tda_out_line0;
}

void tda_out_line_reset(void)
{
    tda_out_lines = 0;
    tda_out_line0 = 0;
}

#else
#define tda_count_lines(fmt, ap) ((void)0)
#endif

/* Under R, the console is a pair of growing buffers in memory: what TDA
   prints on stdout and stderr accumulates here and is handed to R as
   text at the end of the run.  Nothing is written to a file and read
   back; the descriptors 1 and 2 are not touched.  A vsnprintf that
   overruns the stack scratch is formatted again into the buffer itself. */
typedef struct { char *s; size_t len, cap; } tda_cbuf;
static tda_cbuf tda_obuf, tda_ebuf;

/* The console as a FILE*: a few writers (prn_mdes, the idf table in
   t_imat.c) take a stream and are handed TDA_CONSOLE.  Under R that is
   a sentinel -- a real stream opened on the null device, never written
   to -- which tda_fprintf() recognises and turns into a console write
   into the buffer above.  fd 1 is not involved. */
FILE *tda_console(void)
{
    static FILE *fp = NULL;

    if (fp == NULL) {
#ifdef _WIN32
        fp = fopen("NUL", "w");
#else
        fp = fopen("/dev/null", "w");
#endif
    }
    return fp;
}

static void cbuf_reserve(tda_cbuf *b, size_t extra)
{
    size_t need = b->len + extra + 1;
    char *n;

    if (need <= b->cap)
        return;
    while (b->cap < need)
        b->cap = b->cap ? b->cap * 2 : 65536;
    n = (char *)realloc(b->s, b->cap);
    if (n == NULL)
        return;
    b->s = n;
}

static void cbuf_vappend(tda_cbuf *b, const char *fmt, va_list ap)
{
    char scratch[1024];
    va_list aq;
    int n;

    va_copy(aq, ap);
    n = vsnprintf(scratch, sizeof scratch, fmt, aq);
    va_end(aq);
    if (n < 0)
        return;
    cbuf_reserve(b, (size_t)n);
    if (b->s == NULL || b->len + (size_t)n + 1 > b->cap)
        return;
    if ((size_t)n < sizeof scratch)
        memcpy(b->s + b->len, scratch, (size_t)n + 1);
    else {
        va_copy(aq, ap);
        vsnprintf(b->s + b->len, (size_t)n + 1, fmt, aq);
        va_end(aq);
    }
    b->len += (size_t)n;
}

void tda_console_reset(void)
{
    tda_obuf.len = 0;
    tda_ebuf.len = 0;
    if (tda_obuf.s) tda_obuf.s[0] = '\0';
    if (tda_ebuf.s) tda_ebuf.s[0] = '\0';
}

const char *tda_console_text(int which, size_t *len)
{
    tda_cbuf *b = which == 2 ? &tda_ebuf : &tda_obuf;

    if (len) *len = b->len;
    return b->s ? b->s : "";
}

void tda_out(const char *fmt, ...)
{
    va_list ap;

    va_start(ap,fmt);
    tda_count_lines(fmt,ap);
    cbuf_vappend(&tda_obuf, fmt, ap);
    va_end(ap);
}

void tda_vout(const char *fmt, va_list ap)
{
    tda_count_lines(fmt,ap);
    cbuf_vappend(&tda_obuf, fmt, ap);
}

void tda_err(const char *fmt, ...)
{
    va_list ap;

    va_start(ap,fmt);
    cbuf_vappend(&tda_ebuf, fmt, ap);
    va_end(ap);
}

void tda_vout_err(const char *fmt, va_list ap)
{
    cbuf_vappend(&tda_ebuf, fmt, ap);
}

void tda_out_flush(void)
{
}

void tda_err_flush(void)
{
}

#else

void tda_out(const char *fmt, ...)
{
    va_list ap;

    va_start(ap,fmt);
    vprintf(fmt,ap);
    va_end(ap);
}

void tda_vout(const char *fmt, va_list ap)
{
    vprintf(fmt,ap);
}

void tda_err(const char *fmt, ...)
{
    va_list ap;

    va_start(ap,fmt);
    vfprintf(stderr,fmt,ap);
    va_end(ap);
}

void tda_vout_err(const char *fmt, va_list ap)
{
    vfprintf(stderr,fmt,ap);
}

void tda_out_flush(void)
{
    fflush(stdout);
}

void tda_err_flush(void)
{
    fflush(stderr);
}

#endif /* TDA_R_PACKAGE */
