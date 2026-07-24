/* tda_export.c -- direct handout of numeric results from TDA to R.
   Package-only (lives beside tda_r.c, not in the standalone build).

   Producers inside the TDA sources call tda_export_mat() at the point
   where the same numbers are written as text; the text output stays
   byte-identical and every existing parser keeps working.  The list is
   drained into an R object by C_tda_run (tda_r.c) after tda_main()
   returns, and freed with the context.  See CONTRIBUTING.md for
   the migration plan this belongs to. */

#ifdef TDA_R_PACKAGE

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include "tda_context.h"
#include "tda_rhooks.h"

typedef struct TDAExport {
    struct TDAExport *next;
    char name[64];
    int type;               /* 0: numeric matrix, 1: string vector */
    int nrow;
    int ncol;               /* 1 for strings */
    double *v;              /* type 0: column-major copy, nrow * ncol */
    char **s;               /* type 1: nrow copied strings */
} TDAExport;

static TDAExport *node_new(TDAContext *ctx, const char *name)
{
    TDAExport *e, *p;
    int dup;

    e = (TDAExport *)malloc(sizeof(TDAExport));
    if (e == NULL)
        return NULL;
    memset(e, 0, sizeof(*e));
    strncpy(e->name, name, sizeof(e->name) - 1);
    if (ctx->RExports == NULL) {
        ctx->RExports = e;
        return e;
    }
    /* Count only nodes with the SAME name, or that name plus a ".N"
       suffix.  strncmp() on the prefix also matched siblings like
       "coeff.names", so a second "coeff" table was numbered ".3" and
       .exkey_n()'s "key.2" lookup found nothing at all. */
    dup = 0;
    for (p = (TDAExport *)ctx->RExports; ; p = p->next) {
        size_t l = strlen(e->name);
        if (strncmp(p->name, e->name, l) == 0 &&
            (p->name[l] == '\0' ||
             (p->name[l] == '.' && p->name[l + 1] >= '0' &&
              p->name[l + 1] <= '9')))
            dup++;
        if (p->next == NULL)
            break;
    }
    if (dup > 0) {
        size_t l = strlen(e->name);
        if (l < sizeof(e->name) - 4)
            snprintf(e->name + l, sizeof(e->name) - l, ".%d", dup + 1);
    }
    p->next = e;
    return e;
}

void tda_export_mat(TDAContext *ctx, const char *name,
                    const double *v, int nrow, int ncol)
{
    TDAExport *e;
    size_t n;

    if (ctx == NULL || name == NULL || v == NULL || nrow < 1 || ncol < 1)
        return;
    n = (size_t)nrow * (size_t)ncol;
    e = node_new(ctx, name);
    if (e == NULL)
        return;                  /* exporting is best-effort by design */
    e->v = (double *)malloc(n * sizeof(double));
    if (e->v == NULL)
        return;                  /* an empty node is drained as 0x0 */
    memcpy(e->v, v, n * sizeof(double));
    e->nrow = nrow;
    e->ncol = ncol;
}

void tda_export_strings(TDAContext *ctx, const char *name,
                        const char **v, int n)
{
    TDAExport *e;
    int i;

    if (ctx == NULL || name == NULL || v == NULL || n < 1)
        return;
    e = node_new(ctx, name);
    if (e == NULL)
        return;
    e->type = 1;
    e->ncol = 1;
    e->s = (char **)calloc((size_t)n, sizeof(char *));
    if (e->s == NULL)
        return;
    for (i = 0; i < n; ++i)
        e->s[i] = v[i] ? strdup(v[i]) : strdup("");
    e->nrow = n;
}

/* Per-row string staging, the label-side twin of tda_export_row(): one
   call per printed label, flushed into a string-vector export where the
   table closes. */

typedef struct TDAExportSPend {
    struct TDAExportSPend *next;
    char name[64];
    int n;
    int cap;
    char **s;
} TDAExportSPend;

static TDAExportSPend *spend_find(TDAContext *ctx, const char *name)
{
    TDAExportSPend *p;

    for (p = (TDAExportSPend *)ctx->RExportsSPend; p != NULL; p = p->next)
        if (strcmp(p->name, name) == 0)
            return p;
    return NULL;
}

void tda_export_str_row(TDAContext *ctx, const char *name, const char *v)
{
    TDAExportSPend *p;

    if (ctx == NULL || name == NULL || v == NULL)
        return;
    p = spend_find(ctx, name);
    if (p == NULL) {
        p = (TDAExportSPend *)malloc(sizeof(TDAExportSPend));
        if (p == NULL)
            return;
        memset(p, 0, sizeof(*p));
        strncpy(p->name, name, sizeof(p->name) - 1);
        p->cap = 32;
        p->s = (char **)malloc((size_t)p->cap * sizeof(char *));
        if (p->s == NULL) {
            free(p);
            return;
        }
        p->next = (TDAExportSPend *)ctx->RExportsSPend;
        ctx->RExportsSPend = p;
    }
    if (p->n == p->cap) {
        char **nv = (char **)realloc(p->s, (size_t)p->cap * 2u *
                                     sizeof(char *));
        if (nv == NULL)
            return;
        p->s = nv;
        p->cap *= 2;
    }
    p->s[p->n] = strdup(v);
    if (p->s[p->n] != NULL)
        p->n++;
}

void tda_export_str_flush(TDAContext *ctx, const char *name)
{
    TDAExportSPend *p, *q;
    int i;

    if (ctx == NULL || name == NULL)
        return;
    p = spend_find(ctx, name);
    if (p == NULL)
        return;
    if (p->n > 0)
        tda_export_strings(ctx, name, (const char **)p->s, p->n);
    if (ctx->RExportsSPend == p)
        ctx->RExportsSPend = p->next;
    else
        for (q = (TDAExportSPend *)ctx->RExportsSPend; q; q = q->next)
            if (q->next == p) {
                q->next = p->next;
                break;
            }
    for (i = 0; i < p->n; ++i)
        free(p->s[i]);
    free(p->s);
    free(p);
}

/* Row accumulator: producers that print a table one row at a time call
   tda_export_row() beside each row's fprintf and tda_export_flush()
   where the table closes; the flush turns the pending rows into one
   column-major matrix under the given name (so repeated tables in one
   run get .2, .3 through tda_export_mat's own suffixing). */

typedef struct TDAExportPend {
    struct TDAExportPend *next;
    char name[64];
    int ncol;
    int nrow;
    int cap;
    double *rows;            /* row-major while accumulating */
} TDAExportPend;

static TDAExportPend *pend_find(TDAContext *ctx, const char *name)
{
    TDAExportPend *p;

    for (p = (TDAExportPend *)ctx->RExportsPend; p != NULL; p = p->next)
        if (strcmp(p->name, name) == 0)
            return p;
    return NULL;
}

void tda_export_row(TDAContext *ctx, const char *name,
                    const double *row, int ncol)
{
    TDAExportPend *p;

    if (ctx == NULL || name == NULL || row == NULL || ncol < 1)
        return;
    p = pend_find(ctx, name);
    if (p == NULL) {
        p = (TDAExportPend *)malloc(sizeof(TDAExportPend));
        if (p == NULL)
            return;
        strncpy(p->name, name, sizeof(p->name) - 1);
        p->name[sizeof(p->name) - 1] = '\0';
        p->ncol = ncol;
        p->nrow = 0;
        p->cap = 64;
        p->rows = (double *)malloc((size_t)p->cap * (size_t)ncol *
                                   sizeof(double));
        if (p->rows == NULL) {
            free(p);
            return;
        }
        p->next = (TDAExportPend *)ctx->RExportsPend;
        ctx->RExportsPend = p;
    }
    if (ncol != p->ncol)
        return;                  /* shape drift: drop the row, not the run */
    if (p->nrow == p->cap) {
        double *nv = (double *)realloc(p->rows,
                                       (size_t)p->cap * 2u *
                                       (size_t)p->ncol * sizeof(double));
        if (nv == NULL)
            return;
        p->rows = nv;
        p->cap *= 2;
    }
    memcpy(p->rows + (size_t)p->nrow * (size_t)p->ncol, row,
           (size_t)p->ncol * sizeof(double));
    p->nrow++;
}

void tda_export_flush(TDAContext *ctx, const char *name)
{
    TDAExportPend *p, *q;
    double *cm;
    int i, j;

    if (ctx == NULL || name == NULL)
        return;
    p = pend_find(ctx, name);
    if (p == NULL || p->nrow == 0)
        goto unlink;
    cm = (double *)malloc((size_t)p->nrow * (size_t)p->ncol *
                          sizeof(double));
    if (cm != NULL) {
        for (i = 0; i < p->nrow; ++i)
            for (j = 0; j < p->ncol; ++j)
                cm[(size_t)j * (size_t)p->nrow + (size_t)i] =
                    p->rows[(size_t)i * (size_t)p->ncol + (size_t)j];
        tda_export_mat(ctx, name, cm, p->nrow, p->ncol);
        free(cm);
    }
unlink:
    if (p != NULL) {
        if (ctx->RExportsPend == p)
            ctx->RExportsPend = p->next;
        else
            for (q = (TDAExportPend *)ctx->RExportsPend; q; q = q->next)
                if (q->next == p) {
                    q->next = p->next;
                    break;
                }
        free(p->rows);
        free(p);
    }
}

/* Cell accumulator, for rows whose width is only known as they are
   written: the seqmd data matrix appends period dummies, then
   time-independent covariates, then a variable number of
   event-specific columns, each from its own branch.  Producers call
   tda_export_cell() beside each value's fprintf and
   tda_export_endrow() where the row's newline goes; the committed row
   then goes through tda_export_row() like any other. */

typedef struct TDAExportCPend {
    struct TDAExportCPend *next;
    char name[64];
    int n;
    int cap;
    double *v;
} TDAExportCPend;

static TDAExportCPend *cpend_find(TDAContext *ctx, const char *name)
{
    TDAExportCPend *p;

    for (p = (TDAExportCPend *)ctx->RExportsCPend; p != NULL; p = p->next)
        if (strcmp(p->name, name) == 0)
            return p;
    return NULL;
}

void tda_export_cell(TDAContext *ctx, const char *name, double v)
{
    TDAExportCPend *p;

    if (ctx == NULL || name == NULL)
        return;
    p = cpend_find(ctx, name);
    if (p == NULL) {
        p = (TDAExportCPend *)malloc(sizeof(TDAExportCPend));
        if (p == NULL)
            return;
        memset(p, 0, sizeof(*p));
        strncpy(p->name, name, sizeof(p->name) - 1);
        p->cap = 64;
        p->v = (double *)malloc((size_t)p->cap * sizeof(double));
        if (p->v == NULL) {
            free(p);
            return;
        }
        p->next = (TDAExportCPend *)ctx->RExportsCPend;
        ctx->RExportsCPend = p;
    }
    if (p->n == p->cap) {
        double *nv = (double *)realloc(p->v, (size_t)p->cap * 2u *
                                       sizeof(double));
        if (nv == NULL)
            return;
        p->v = nv;
        p->cap *= 2;
    }
    p->v[p->n++] = v;
}

void tda_export_endrow(TDAContext *ctx, const char *name)
{
    TDAExportCPend *p;

    if (ctx == NULL || name == NULL)
        return;
    p = cpend_find(ctx, name);
    if (p == NULL)
        return;
    if (p->n > 0)
        tda_export_row(ctx, name, p->v, p->n);
    p->n = 0;
}

static void cpend_free(TDAContext *ctx)
{
    TDAExportCPend *p, *q;

    for (p = (TDAExportCPend *)ctx->RExportsCPend; p != NULL; p = q) {
        q = p->next;
        free(p->v);
        free(p);
    }
    ctx->RExportsCPend = NULL;
}

/* Producers whose table has no natural closing point (prn_sfmt's
   label/value pairs run to the end of whatever command is printing)
   leave rows pending; the bridge calls this once after tda_main()
   returns so nothing accumulated is lost. */

void tda_export_flush_all(TDAContext *ctx)
{
    if (ctx == NULL)
        return;
    while (ctx->RExportsPend != NULL)
        tda_export_flush(ctx, ((TDAExportPend *)ctx->RExportsPend)->name);
    while (ctx->RExportsSPend != NULL)
        tda_export_str_flush(ctx,
                             ((TDAExportSPend *)ctx->RExportsSPend)->name);
    cpend_free(ctx);
}

/* File-writer tap.  TDA writes its data files with fprintf() to a
   stream held on the context (PMFd, PMF1d, PSFd, ...), ~1900 numeric
   sites across the sources.  tda.h redirects those calls here under the
   package build, so the numbers reach R without a producer beside each
   one; the write itself is unchanged.  Values are staged per stream,
   under "file.<stream>.values", as (line, value) -- the same shape the
   console tap uses -- with line numbers counted per stream. */

#define TDA_FDTAP_MAX 12

static const char *fdtap_stream_name(int slot)
{
    static const char *nm[TDA_FDTAP_MAX] = {
        "PMFd", "PMF1d", "PSFd", "PMTabFd", "PMResFd", "PMProtFd",
        "PMPCFd", "PMPPFd", "PMTDAFd", "PMDVARFd", "PRTFd", "other"
    };
    return (slot >= 0 && slot < TDA_FDTAP_MAX) ? nm[slot] : "other";
}

static const char *fdtap_name(TDAContext *ctx, FILE *fd, int *slot)
{
    static const char *nm[TDA_FDTAP_MAX] = {
        "PMFd", "PMF1d", "PSFd", "PMTabFd", "PMResFd", "PMProtFd",
        "PMPCFd", "PMPPFd", "PMTDAFd", "PMDVARFd", "PRTFd", "other"
    };
    FILE *fds[TDA_FDTAP_MAX];
    int i;

    fds[0] = ctx->PMFd;      fds[1] = ctx->PMF1d;   fds[2] = ctx->PSFd;
    fds[3] = ctx->PMTabFd;   fds[4] = ctx->PMResFd; fds[5] = ctx->PMProtFd;
    fds[6] = ctx->PMPCFd;    fds[7] = ctx->PMPPFd;  fds[8] = ctx->PMTDAFd;
    fds[9] = NULL;           fds[10] = NULL;        fds[11] = NULL;
    for (i = 0; i < 9; ++i)
        if (fd != NULL && fds[i] == fd) {
            *slot = i;
            return nm[i];
        }
    *slot = TDA_FDTAP_MAX - 1;
    return nm[TDA_FDTAP_MAX - 1];
}

/* The file each stream is writing, as TDA has it on the context, so R
   can tell which tapped stream a "df=" or "pres=" file is without
   opening the file: handed over once per stream and run, as the string
   export "file.<stream>.name". */
static const char *fdtap_filename(TDAContext *ctx, int slot)
{
    switch (slot) {
        case 0: return ctx->PMFdName;
        case 1: return ctx->PMF1dName;
        case 3: return ctx->PMTabFName;
        case 4: return ctx->PMResFName;
        case 5: return ctx->PMProtFName;
        case 6: return ctx->PMPCFName;
        case 7: return ctx->PMPPFName;
        case 8: return ctx->PMTDAFName;
        default: return NULL;
    }
}

static int fdtap_slot(TDAContext *ctx, FILE *fd, char *key, size_t keysz)
{
    const char *nm, *fn;
    int slot;

    if (ctx == NULL || fd == NULL)
        return -1;
    nm = fdtap_name(ctx, fd, &slot);
    if (slot < 0 || slot >= TDA_FDTAP_MAX)
        return -1;
    snprintf(key, keysz, "file.%s.values", nm);
    if (!ctx->RExportFdNamed[slot]) {
        fn = fdtap_filename(ctx, slot);
        if (fn != NULL && fn[0]) {
            char nkey[64];
            snprintf(nkey, sizeof nkey, "file.%s.name", nm);
            tda_export_str_row(ctx, nkey, fn);
            tda_export_str_flush(ctx, nkey);
        }
        ctx->RExportFdNamed[slot] = 1;
    }
    return slot;
}

static void fdtap_values(TDAContext *ctx, int slot, const char *key,
                         const char *fmt, va_list ap)
{
    const char *p;
    double v, row[2];
    int lmod, has, line;

    line = ctx->RExportFdLine[slot] + 1;

    /* walk the arguments exactly as printf would, so the sequence stays
       in step; the numeric ones are staged, the rest only consumed */
    for (p = fmt; *p; ++p) {
        if (*p != '%')
            continue;
        ++p;
        if (*p == '%')
            continue;
        while (*p && strchr("-+ #0'", *p)) ++p;
        if (*p == '*') { (void)va_arg(ap, int); ++p; }
        else while (*p >= '0' && *p <= '9') ++p;
        if (*p == '.') {
            ++p;
            if (*p == '*') { (void)va_arg(ap, int); ++p; }
            else while (*p >= '0' && *p <= '9') ++p;
        }
        lmod = 0;
        while (*p && strchr("hlLqjzt", *p)) {
            if (*p == 'l') lmod++;
            else if (*p == 'L') lmod = 3;
            else if (*p == 'h') lmod = -1;
            ++p;
        }
        if (!*p) break;
        has = 1; v = 0.0;
        switch (*p) {
            case 'd': case 'i':
                if (lmod >= 2)      v = (double)va_arg(ap, long long);
                else if (lmod == 1) v = (double)va_arg(ap, long);
                else                v = (double)va_arg(ap, int);
                break;
            case 'u': case 'o': case 'x': case 'X':
                if (lmod >= 2)      v = (double)va_arg(ap, unsigned long long);
                else if (lmod == 1) v = (double)va_arg(ap, unsigned long);
                else                v = (double)va_arg(ap, unsigned int);
                break;
            case 'e': case 'E': case 'f': case 'F':
            case 'g': case 'G': case 'a': case 'A':
                if (lmod == 3)      v = (double)va_arg(ap, long double);
                else                v = va_arg(ap, double);
                break;
            case 'c': (void)va_arg(ap, int);   has = 0; break;
            case 's': (void)va_arg(ap, char *); has = 0; break;
            case 'p': (void)va_arg(ap, void *); has = 0; break;
            case 'n': (void)va_arg(ap, int *);  has = 0; break;
            default:  has = 0; break;
        }
        if (has) {
            row[0] = (double)line;
            row[1] = v;
            tda_export_row(ctx, key, row, 2);
        }
    }
}

/* Numbers TDA writes as literals in the format itself -- the "1 3 4"
   object header of a spatial data file, a fixed "0" -- have no argument
   to tap; they are integer constants in the C, exact by definition, so
   they are staged from the format text.  Only a format without any
   conversion is scanned, and only whole numbers standing on their own. */
static void fdtap_literals(TDAContext *ctx, int slot, const char *key,
                           const char *fmt)
{
    const char *p;
    double row[2];
    int line, neg;
    long v;

    if (fmt == NULL || strchr(fmt, '%') != NULL)
        return;
    line = ctx->RExportFdLine[slot] + 1;
    for (p = fmt; *p; ) {
        if (*p == '-' && p[1] >= '0' && p[1] <= '9' &&
            (p == fmt || p[-1] == ' ' || p[-1] == '\t')) {
            neg = 1; ++p;
        } else if (*p >= '0' && *p <= '9' &&
                   (p == fmt || p[-1] == ' ' || p[-1] == '\t')) {
            neg = 0;
        } else {
            ++p;
            continue;
        }
        v = 0;
        while (*p >= '0' && *p <= '9') { v = v * 10 + (*p - '0'); ++p; }
        if (*p == '\0' || *p == ' ' || *p == '\t' || *p == '\n') {
            row[0] = (double)line;
            row[1] = neg ? -(double)v : (double)v;
            tda_export_row(ctx, key, row, 2);
        } else {
            while (*p && *p != ' ' && *p != '\t' && *p != '\n') ++p;
        }
    }
}

/* Every line TDA writes to a tapped stream is handed over as text as
   well: "file.<stream>.lines", one string per line, so an output whose
   content is text -- gdot's Graphviz source, pcyc's cycle notation,
   rplz's place names, a protocol -- reaches R without the file.  The
   "#" comment lines above a data table also go to
   "file.<stream>.comments", which is what the table readers use for
   their column labels. */
static void fdtap_comment(TDAContext *ctx, int slot, const char *text, int n)
{
    const char *q;
    char *buf = ctx->RExportFdCom[slot];
    int *len = &ctx->RExportFdComLen[slot];
    char key[64];

    for (q = text; n > 0 && *q; ++q) {
        if (ctx->RExportFdAtBOL[slot]) {
            ctx->RExportFdAtBOL[slot] = 0;
            *len = 0;
            ctx->RExportFdIsCom[slot] = (*q == '#');
        }
        if (*q == '\n') {
            buf[*len] = '\0';
            snprintf(key, sizeof key, "file.%s.lines", fdtap_stream_name(slot));
            tda_export_str_row(ctx, key, buf);
            if (ctx->RExportFdIsCom[slot]) {
                snprintf(key, sizeof key, "file.%s.comments",
                         fdtap_stream_name(slot));
                tda_export_str_row(ctx, key, buf + 1);
            }
            ctx->RExportFdAtBOL[slot] = 1;
            ctx->RExportFdIsCom[slot] = 0;
            *len = 0;
        }
        else if (*len < (int)sizeof(ctx->RExportFdCom[slot]) - 1) {
            buf[(*len)++] = *q;
        }
    }
}

/* Is the fragment about to be written part of a "#" comment line? */
static int fdtap_on_comment(TDAContext *ctx, int slot, const char *text)
{
    if (ctx->RExportFdAtBOL[slot])
        return text != NULL && text[0] == '#';
    return ctx->RExportFdIsCom[slot];
}

static void fdtap_text(TDAContext *ctx, FILE *fd, int slot, const char *text, int n)
{
    const char *q;

    fdtap_comment(ctx, slot, text, n);
    for (q = text; n > 0 && *q; ++q)
        if (*q == '\n')
            ctx->RExportFdLine[slot]++;
    /* The dtda description stream carries the COLUMN NAMES of whatever
       data file the command just wrote.  A name is text, and the export
       channel carries typed string vectors perfectly well -- coeff.names
       and the rest are exactly that -- so the names are handed over AS
       a string vector rather than as raw text for R to pick apart.
       Every such line has the shape "  NAME <w>[fmt] = cK, # comment",
       so the leading token of a fragment containing "= c" is the name.
       Reading our own writer's fixed format at the point of writing is
       a producer, not a parse of a file after the fact. */
    if (fd == ctx->PMTDAFd && n > 0) {
        const char *t = text;
        char nm1[64];
        size_t k = 0;

        while (*t == ' ' || *t == '\t' || *t == '\n')
            ++t;
        while (*t && *t != ' ' && *t != '\t' && *t != '<' && *t != ',' &&
               k < sizeof(nm1) - 1)
            nm1[k++] = *t++;
        nm1[k] = '\0';
        if (strstr(text, "= c") != NULL) {
            /* A name starts with a letter.  qreg writes the name and the
               "<w>[fmt] = cK," part in SEPARATE calls, so this fragment's
               own leading token is "=" -- the name is the one remembered
               from the previous fragment. */
            const char *use = NULL;
            if (k > 0 && ((nm1[0] >= 'A' && nm1[0] <= 'Z') ||
                          (nm1[0] >= 'a' && nm1[0] <= 'z')))
                use = nm1;
            else if (ctx->RExportDtdaPend[0])
                use = ctx->RExportDtdaPend;
            if (use != NULL)
                tda_export_str_row(ctx, "dtda.names", use);
            ctx->RExportDtdaPend[0] = '\0';
        }
        else if (k > 0 && ((nm1[0] >= 'A' && nm1[0] <= 'Z') ||
                           (nm1[0] >= 'a' && nm1[0] <= 'z'))) {
            /* a bare token: remember it in case its "= c" half follows */
            strncpy(ctx->RExportDtdaPend, nm1,
                    sizeof(ctx->RExportDtdaPend) - 1);
            ctx->RExportDtdaPend[sizeof(ctx->RExportDtdaPend) - 1] = '\0';
        }
    }
}

int tda_fprintf(TDAContext *ctx, FILE *fd, const char *fmt, ...)
{
    va_list ap, aq;
    char key[80];
    char buf[1024];
    char *big = NULL;
    int r, n, slot;

    if (fd == tda_console()) {
        /* a console write through a FILE*: into the buffered console,
           line-counted and value-tapped like any printf1() */
        va_start(ap, fmt);
        tda_vout(fmt, ap);
        va_end(ap);
        return 0;
    }
    va_start(ap, fmt);
    r = vfprintf(fd, fmt, ap);
    va_end(ap);

    if (fmt == NULL || (slot = fdtap_slot(ctx, fd, key, sizeof(key))) < 0)
        return r;

    /* newlines are counted on the formatted text, not the format: a
       "%s" argument carrying one would otherwise desynchronise the
       line numbers, exactly as it did for the console tap */
    va_start(aq, fmt);
    n = vsnprintf(buf, sizeof(buf), fmt, aq);
    va_end(aq);
    if (n > 0 && (size_t)n >= sizeof(buf)) {
        big = (char *)malloc((size_t)n + 1);
        if (big != NULL) {
            va_start(aq, fmt);
            vsnprintf(big, (size_t)n + 1, fmt, aq);
            va_end(aq);
        }
    }
    /* a number inside a "#" comment line ("# c3 : Vx") is text, not a
       value of the table: the fragment is judged before its values are
       staged, and skipped when it is on a comment line */
    if (!fdtap_on_comment(ctx, slot, big != NULL ? big : buf)) {
        va_start(ap, fmt);
        fdtap_values(ctx, slot, key, fmt, ap);
        va_end(ap);
        fdtap_literals(ctx, slot, key, fmt);
    }
    fdtap_text(ctx, fd, slot, big != NULL ? big : buf, n);
    free(big);
    return r;
}

/* tda_fprintf() with the text already formatted (tda_rtfmt.c); the
   arguments after fmt are fmt's own and only feed the tap. */
int tda_fprintf_rt(TDAContext *ctx, FILE *fd, const char *text, const char *fmt, ...)
{
    va_list ap;
    char key[80];
    size_t len = strlen(text);
    int r, slot;

    if (fd == tda_console()) {
        tda_out("%s", text);
        return (int)len;
    }
    r = fwrite(text, 1, len, fd) == len ? (int)len : -1;

    if (fmt == NULL || (slot = fdtap_slot(ctx, fd, key, sizeof(key))) < 0)
        return r;

    va_start(ap, fmt);
    if (!fdtap_on_comment(ctx, slot, text))
        fdtap_values(ctx, slot, key, fmt, ap);
    va_end(ap);
    fdtap_text(ctx, fd, slot, text, len > (size_t)INT_MAX ? INT_MAX : (int)len);
    return r;
}

/* Read-side accessors, so the .Call bridge never has to know the list
   node layout.  Index i is 0-based, in production order. */
static TDAExport *export_at(TDAContext *ctx, int i)
{
    TDAExport *e;

    if (ctx == NULL)
        return NULL;
    for (e = (TDAExport *)ctx->RExports; e != NULL && i > 0; e = e->next)
        i--;
    return i == 0 ? e : NULL;
}

int tda_export_count(TDAContext *ctx)
{
    TDAExport *e;
    int n = 0;

    if (ctx == NULL)
        return 0;
    for (e = (TDAExport *)ctx->RExports; e != NULL; e = e->next)
        n++;
    return n;
}

const char *tda_export_name(TDAContext *ctx, int i)
{
    TDAExport *e = export_at(ctx, i);
    return e ? e->name : "";
}

int tda_export_nrow(TDAContext *ctx, int i)
{
    TDAExport *e = export_at(ctx, i);
    return e ? e->nrow : 0;
}

int tda_export_ncol(TDAContext *ctx, int i)
{
    TDAExport *e = export_at(ctx, i);
    return e ? e->ncol : 0;
}

const double *tda_export_data(TDAContext *ctx, int i)
{
    TDAExport *e = export_at(ctx, i);
    return e ? e->v : NULL;
}

int tda_export_type(TDAContext *ctx, int i)
{
    TDAExport *e = export_at(ctx, i);
    return e ? e->type : 0;
}

const char *tda_export_str(TDAContext *ctx, int i, int j)
{
    TDAExport *e = export_at(ctx, i);

    if (e == NULL || e->type != 1 || e->s == NULL || j < 0 || j >= e->nrow)
        return "";
    return e->s[j] ? e->s[j] : "";
}

void tda_export_free(TDAContext *ctx)
{
    TDAExport *e, *n;

    if (ctx == NULL)
        return;
    for (e = (TDAExport *)ctx->RExports; e != NULL; e = n) {
        n = e->next;
        free(e->v);
        if (e->s != NULL) {
            int i;
            for (i = 0; i < e->nrow; ++i)
                free(e->s[i]);
            free(e->s);
        }
        free(e);
    }
    ctx->RExports = NULL;
    {
        TDAExportPend *p, *pn;
        for (p = (TDAExportPend *)ctx->RExportsPend; p != NULL; p = pn) {
            pn = p->next;
            free(p->rows);
            free(p);
        }
        ctx->RExportsPend = NULL;
    }
    {
        TDAExportSPend *p, *pn;
        int i;
        for (p = (TDAExportSPend *)ctx->RExportsSPend; p != NULL; p = pn) {
            pn = p->next;
            for (i = 0; i < p->n; ++i)
                free(p->s[i]);
            free(p->s);
            free(p);
        }
        ctx->RExportsSPend = NULL;
    }
}

#else

/* the standalone build never compiles this file, but keep it valid C
   if it ever does */
typedef int tda_export_not_used;

#endif /* TDA_R_PACKAGE */
