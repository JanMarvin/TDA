#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>

#include "tda_context.h"
#include "tda_rhooks.h"

void tda_rdf_set(SEXP df);
void tda_rdf_clear(void);

int tda_main(TDAContext *ctx, int argc, char *argv[]);

static jmp_buf tda_r_jb;
static int tda_r_jb_active = 0;
static int tda_r_status = 0;

void tda_r_exit(int status)
{
    if (tda_r_jb_active) {
        tda_r_status = status;
        longjmp(tda_r_jb, 1);
    }
    Rf_error("TDA called exit(%d) outside of a run", status);
}


/* One R character vector per console stream: the buffered text split at
   its newlines, the way readLines() would have split a file. */
static SEXP console_lines(int which)
{
    size_t len, i, start, n;
    const char *t = tda_console_text(which, &len);
    SEXP v;

    n = 0;
    for (i = 0; i < len; ++i)
        if (t[i] == '\n') n++;
    if (len > 0 && t[len - 1] != '\n') n++;
    v = PROTECT(Rf_allocVector(STRSXP, (R_xlen_t)n));
    start = 0; n = 0;
    for (i = 0; i <= len; ++i) {
        if (i == len || t[i] == '\n') {
            if (i == len && start == len)
                break;
            SET_STRING_ELT(v, (R_xlen_t)n++,
                           Rf_mkCharLen(t + start, (int)(i - start)));
            start = i + 1;
        }
    }
    UNPROTECT(1);
    return v;
}

SEXP C_tda_run(SEXP args, SEXP data, SEXP commands)
{
    int argc, i, rc;
    char **argv;
    TDAContext *ctx;
    char *cftext = NULL;

    if (TYPEOF(args) != STRSXP)
        Rf_error("args must be a character vector");

    argc = (int) LENGTH(args) + 1;
    argv = (char **)(void *) R_alloc((size_t) argc, sizeof(char *));
    argv[0] = "tda";
    for (i = 1; i < argc; ++i)
        argv[i] = (char *) CHAR(STRING_ELT(args, i - 1));

    ctx = tda_context_new();
    if (!ctx)
        Rf_error("cannot allocate TDA context");
    tda_context_init(ctx);
    tda_reset_globals();
    tda_interrupt_reset();

    tda_out_line_reset();
    memset(ctx->RExportFdLine, 0, sizeof(ctx->RExportFdLine));
    memset(ctx->RExportFdNamed, 0, sizeof(ctx->RExportFdNamed));
    memset(ctx->RExportFdComLen, 0, sizeof(ctx->RExportFdComLen));
    memset(ctx->RExportFdIsCom, 0, sizeof(ctx->RExportFdIsCom));
    for (i = 0; i < 12; ++i) ctx->RExportFdAtBOL[i] = 1;
    ctx->RExportPleRows = 0;
    ctx->RExportDtdaPend[0] = '\0';
    tda_console_reset();

    /* the commands, one string per line, joined into the text the
       parser reads in place of a command file */
    if (TYPEOF(commands) == STRSXP && LENGTH(commands) > 0) {
        size_t tot = 0, off = 0;
        for (i = 0; i < LENGTH(commands); ++i)
            tot += strlen(CHAR(STRING_ELT(commands, i))) + 1;
        cftext = (char *)R_alloc(tot + 1, 1);
        for (i = 0; i < LENGTH(commands); ++i) {
            const char *l = CHAR(STRING_ELT(commands, i));
            size_t n = strlen(l);
            memcpy(cftext + off, l, n);
            off += n;
            cftext[off++] = '\n';
        }
        cftext[off] = '\0';
    }
    tda_cf_set(cftext);

    /*  The data frame, if there is one, stays reachable from the caller for
        the whole run, so rdataframe() can read it without protecting it
        again.  It is cleared afterwards so a later run cannot see a stale
        pointer. */
    tda_rdf_set(data == R_NilValue ? NULL : data);

    tda_r_status = 0;
    tda_r_jb_active = 1;
    if (setjmp(tda_r_jb) == 0)
        rc = tda_main(ctx, argc, argv);
    else
        rc = tda_r_status;
    tda_r_jb_active = 0;

    tda_rdf_clear();
    tda_cf_set(NULL);

    /* Drain the direct exports (tda_export.c) into an R list BEFORE the
       context -- and with it the export list -- is freed.  The result
       becomes list(rc = <int>, exports = <named list of matrices>,
       errors = <int>, the count of errors TDA reported);
       tda_run() on the R side unpacks it, and everything downstream of
       the printed output is untouched (CONTRIBUTING.md). */
    {
        SEXP exps, nms, res2;
        int nexp, k, errors;

        errors = ctx->ErrCnt;

        tda_export_flush_all(ctx);
        nexp = tda_export_count(ctx);
        exps = PROTECT(Rf_allocVector(VECSXP, nexp));
        nms = PROTECT(Rf_allocVector(STRSXP, nexp));
        for (k = 0; k < nexp; ++k) {
            int nr = tda_export_nrow(ctx, k);
            int nc = tda_export_ncol(ctx, k);
            SEXP m;
            if (tda_export_type(ctx, k) == 1) {
                int j;
                m = PROTECT(Rf_allocVector(STRSXP, nr));
                for (j = 0; j < nr; ++j)
                    SET_STRING_ELT(m, j,
                                   Rf_mkChar(tda_export_str(ctx, k, j)));
            }
            else {
                m = PROTECT(Rf_allocMatrix(REALSXP, nr, nc));
                if (nr > 0 && nc > 0 && tda_export_data(ctx, k) != NULL)
                    memcpy(REAL(m), tda_export_data(ctx, k),
                           (size_t) nr * (size_t) nc * sizeof(double));
            }
            SET_VECTOR_ELT(exps, k, m);
            SET_STRING_ELT(nms, k, Rf_mkChar(tda_export_name(ctx, k)));
            UNPROTECT(1);
        }
        Rf_setAttrib(exps, R_NamesSymbol, nms);

        tda_context_free(ctx);

        res2 = PROTECT(Rf_allocVector(VECSXP, 5));
        SET_VECTOR_ELT(res2, 0, Rf_ScalarInteger(rc));
        SET_VECTOR_ELT(res2, 1, exps);
        SET_VECTOR_ELT(res2, 2, Rf_ScalarInteger(errors));
        SET_VECTOR_ELT(res2, 3, console_lines(1));
        SET_VECTOR_ELT(res2, 4, console_lines(2));
        {
            SEXP rn = PROTECT(Rf_allocVector(STRSXP, 5));
            SET_STRING_ELT(rn, 0, Rf_mkChar("rc"));
            SET_STRING_ELT(rn, 1, Rf_mkChar("exports"));
            SET_STRING_ELT(rn, 2, Rf_mkChar("errors"));
            SET_STRING_ELT(rn, 3, Rf_mkChar("output"));
            SET_STRING_ELT(rn, 4, Rf_mkChar("stderr"));
            Rf_setAttrib(res2, R_NamesSymbol, rn);
            UNPROTECT(1);
        }
        UNPROTECT(3);
        return res2;
    }
}
