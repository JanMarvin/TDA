/* tda_rdf.c                                                                */
/*                                                                          */
/*  rdataframe -- build TDA's internal data matrix from an R data frame     */
/*  that is already in memory, instead of writing it out as text and        */
/*  parsing it back.                                                        */
/*                                                                          */
/*  This command exists only in the R package.  The standalone program has  */
/*  no way to be handed an R object, so t_cmd.c dispatches it under         */
/*  TDA_R_PACKAGE and nothing else in the command language changes.  It     */
/*  follows the same path as rstata() and rspss1(): save_var() for each     */
/*  column, then alloc_vdat() and put_data() for the values, so every       */
/*  model, plot and table downstream reads it through get_data() without    */
/*  knowing where the numbers came from.                                    */
/*                                                                          */
/*  The values are copied rather than borrowed.  Pointing VDPtr[] straight  */
/*  at REAL(x) would save the copy, but TDA owns and reallocates that       */
/*  memory and R's garbage collector owns the vector, so the two would have */
/*  to agree about lifetime; a single pass with no parsing and no           */
/*  formatting is already the whole of the win over a text file.            */
/*                                                                          */
/*  Missing values follow the convention rstata() and rspss1() use: NA is   */
/*  replaced by ctx->PMMSYS, which msys= sets and which defaults to -5.     */
/*                                                                          */
/*  Character and factor columns are refused, and the reason is one small   */
/*  addition away: it is the right-hand side of a variable definition that  */
/*  sets the type, and save_var() in t_var.c has branches for =str(n,m),    */
/*  =spss(len) and =stata(len) that set VTyp = 1.  "=c%d" always means a    */
/*  number, so a string column gets its storage but stays typed numeric.    */
/*  An =rdataframe(len) branch beside the others, and put_str() in place of */
/*  put_data(), is all it needs -- put_str() takes the string itself as its */
/*  buffer and pads to the width, so no record buffer has to be built.      */

#include "tda.h"

/*  tda.h does not pull these in; each module keeps its own header, and the
    functions this command builds on live in three of them. */
#include "t_gen.h"      /* printf1                                          */
#include "t_var.h"      /* get_nidx, save_var, prn_var, alloc_vdat          */
#include "t_gdat.h"     /* put_data                                         */

#ifdef TDA_R_PACKAGE

#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

/*  The frame the R side parks here for the duration of one run.  It is
    protected on the R side, so nothing is allocated or freed here. */
static SEXP tda_rdf_data = NULL;

void tda_rdf_set(SEXP df)
{
    tda_rdf_data = df;
}

void tda_rdf_clear(void)
{
    tda_rdf_data = NULL;
}

/*  A column's print format.  TDA needs a width and a number of decimals for
    its own output; integers get none, doubles enough to survive the round
    trip through any file TDA writes later. */

/*  The <n> in a variable definition is the storage code get_data() switches
    on, not a placeholder: 5 is an int array and 8 a double array, and 0 is a
    bitmap, which is what an integer column silently became on the first
    attempt.  Doubles are stored as doubles; everything else fits in an int. */

static void rdf_fmt(SEXP col, int *n, int *w, int *d)
{
    if (TYPEOF(col) == INTSXP || TYPEOF(col) == LGLSXP) {
        *n = 5;
        *w = 12;
        *d = 0;
    }
    else {
        *n = 8;
        *w = 24;
        *d = 16;
    }
}

int rdataframe(TDAContext *ctx)
{
    int err, idxn, nvar, noc, i, j, n, w, d, nmiss;
    double x;
    SEXP names, col;
    char vd[VNLMax + 64];

    err = -1;

    if (tda_rdf_data == NULL) {
        printf1(ctx, "Error: no data frame was passed to this run.\n");
        return(-1);
    }
    if (TYPEOF(tda_rdf_data) != VECSXP) {
        printf1(ctx, "Error: rdataframe needs a list of columns.\n");
        return(-1);
    }

    nvar = (int)Rf_length(tda_rdf_data);
    if (nvar < 1) {
        printf1(ctx, "Error: the data frame has no columns.\n");
        return(-1);
    }
    names = Rf_getAttrib(tda_rdf_data,R_NamesSymbol);
    if (TYPEOF(names) != STRSXP || Rf_length(names) != nvar) {
        printf1(ctx, "Error: the data frame has no column names.\n");
        return(-1);
    }
    noc = (int)Rf_length(VECTOR_ELT(tda_rdf_data,0));
    for (j = 0; j < nvar; ++j) {
        if ((int)Rf_length(VECTOR_ELT(tda_rdf_data,j)) != noc) {
            printf1(ctx, "Error: column %d has a different length.\n",j + 1);
            return(-1);
        }
    }
    if (noc < 1) {
        printf1(ctx, "Error: the data frame has no rows.\n");
        return(-1);
    }

    /*  register one variable per column */

    idxn = get_nidx(ctx);

    for (j = 0; j < nvar; ++j) {
        col = VECTOR_ELT(tda_rdf_data,j);
        if (TYPEOF(col) != INTSXP && TYPEOF(col) != REALSXP &&
            TYPEOF(col) != LGLSXP) {
            printf1(ctx,
                "Error: column %s is not numeric; character and factor\n"
                "columns need a string branch in save_var() -- see the note\n"
                "at the top of this file.\n",
                CHAR(STRING_ELT(names,j)));
            return(-1);
        }
        rdf_fmt(col,&n,&w,&d);
        if ((int)strlen(CHAR(STRING_ELT(names,j))) > VNLMax) {
            printf1(ctx, "Error: variable name is too long: %s\n",
                    CHAR(STRING_ELT(names,j)));
            return(-1);
        }
        snprintf(vd,sizeof(vd),"%s<%d>[%d.%d]=c%d",
                 CHAR(STRING_ELT(names,j)),n,w,d,j + 1);
        if (save_var(ctx, vd,0)) {
            printf1(ctx, "\nError: can't save variable definitions.\n");
            return(-1);
        }
    }
    prn_var(ctx, idxn);
    printf1(ctx, "\n");

    /*  allocate the data matrix and fill it */

    ctx->NOCMaxA = noc;
    printf1(ctx, "Reading a data frame to create internal data matrix.\n");
    printf1(ctx, "Maximum number of cases: %d\n",noc);

    if (alloc_vdat(ctx, idxn,1)) {
        printf1(ctx, "Error: insufficient memory for data matrix (%d cases).\n",
                noc);
        return(-1);
    }

    nmiss = 0;
    for (j = 0; j < nvar; ++j) {
        col = VECTOR_ELT(tda_rdf_data,j);
        for (i = 0; i < noc; ++i) {
            if (TYPEOF(col) == REALSXP) {
                x = REAL(col)[i];
                if (ISNA(x) || ISNAN(x)) {
                    x = ctx->PMMSYS;
                    nmiss++;
                }
            }
            else {
                int v = INTEGER(col)[i];

                if (v == NA_INTEGER) {
                    x = ctx->PMMSYS;
                    nmiss++;
                }
                else
                    x = (double)v;
            }
            put_data(ctx, x,idxn + j,i);
        }
    }

    ctx->NOCDM = ctx->NOC = noc;
    ctx->DMDef = 1;
    /* what nvar sets once a data matrix exists (t_gdat.c): the sum of
       case weights, NOC until a cwt command says otherwise.  Left at
       zero, every statistic normalised by it -- the Breslow and
       Tarone-Ware weights of ple's csf, for one -- divides by zero. */
    ctx->WIVar = -1;
    ctx->WSum = (double)ctx->NOC;
    ctx->WSumS = 0.0;
    ctx->WNorm = 1.0;
    ctx->WNormFlag = 0;

    printf1(ctx, "Created a data matrix with %d variables and %d cases.\n",
            ctx->NVAR,ctx->NOC);
    if (nmiss > 0)
        printf1(ctx, "%d missing values were substituted by: %g\n",
                nmiss,ctx->PMMSYS);

    err = 0;
    return(err);
}

#endif /* TDA_R_PACKAGE */
