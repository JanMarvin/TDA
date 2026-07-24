#include <stdlib.h>
#include "tda_context.h"

TDAContext *tda_context_new(void) {
    return calloc(1, sizeof(TDAContext));
}

/* The teardown functions live across the t_*.c files; declaring them
   here keeps tda_context.c free of the per-file headers.  Signatures
   verified against the definitions (a first draft invented them and
   clear_var(ctx, 0) cleared VARIABLE 0 -- intermittent corruption
   until audited). */
void clear_dm(TDAContext *ctx);
void clear_avar(TDAContext *ctx, int idx);
void mat_free(TDAContext *ctx);
void gdd_free(TDAContext *ctx, int opt);
void sdnvar_close(TDAContext *ctx);
void edat_off(TDAContext *ctx, int opt);
void seq_afree(TDAContext *ctx, int opt);

void tda_context_free(TDAContext *ctx) {
    int i;

    if (ctx == NULL)
        return;

    /* the teardown reuses the clear command's paths, which report
       what they remove via printf1; in the R package that printing
       happens after the run's output redirection has closed and must
       stay off (SILENTFlg = 3 mutes printf1 and printf2 alike) */
    ctx->SILENTFlg = 3;

    if (ctx->VIFirst >= 0)
        clear_avar(ctx, ctx->VIFirst);   /* all variables, like clear */
    clear_dm(ctx);
    mat_free(ctx);
    /* the eight alloc_mat top arrays; mat_free releases per-matrix
       storage but is also a mid-run "clear all" (t_cmd.c), so the
       arrays themselves can only go here */
    if (ctx->MatName != NULL && ctx->MatDef != NULL &&
        ctx->MatLocDef != NULL)
        for (i = 0; i < ctx->MaxMat; ++i) {
            free(ctx->MatName[i]);
            free(ctx->MatDef[i]);
            free(ctx->MatLocDef[i]);
        }
    free(ctx->MatAlloc);  free(ctx->MatLoc);  free(ctx->MatRow);
    free(ctx->MatCol);    free(ctx->MatLocDef);
    free(ctx->MatName);   free(ctx->MatDef);  free(ctx->MatVal);
    gdd_free(ctx, 1);
    sdnvar_close(ctx);
    edat_off(ctx, 1);
    seq_afree(ctx, 1);
    /* the expression stack: alloc_est's own free form trips its
       internal sanity check when called on this state, so release
       the three arrays directly */
    /* the per-command constraint flags p_con() builds (t_con.c); read
       by prn1_coeff() while the fit prints, so they outlive p_con and
       can only be released here (valgrind: glm6.cf, lsreg2.cf) */
    free(ctx->ParFixed);  ctx->ParFixed = NULL;
    free(ctx->ESTyp);  ctx->ESTyp = NULL;
    free(ctx->ESIdx);  ctx->ESIdx = NULL;
    free(ctx->ESVal);  ctx->ESVal = NULL;
    ctx->ESMaxLen = 0;


    /* the same teardown the clear command performs (its callees all
       tolerate an empty state), plus the long-lived arrays clear
       keeps around for reuse -- without these, every run leaked the
       whole variable system (about 390 KB per run) */
    if (ctx->VPFmtS != NULL)
        for (i = 0; i < ctx->MaxNV; ++i)
            free(ctx->VPFmtS[i]);        /* allocated by alloc_vmax;
                                            clear_avar frees per-var
                                            data, not these slots */
    free(ctx->VAlloc);  free(ctx->VTyp);   free(ctx->VTypA);
    free(ctx->VSLen);   free(ctx->VPFmt1); free(ctx->VPFmt2);
    free(ctx->VPFmtS);  free(ctx->VNxt);   free(ctx->VESCnt);
    free(ctx->VESTyp);  free(ctx->AVVAL);  free(ctx->VESVal);
    /* the five remaining alloc_vmax top arrays -- their per-variable
       contents are released by clear_avar above, the pointer arrays
       themselves were not freed anywhere */
    free(ctx->VDef);    free(ctx->VName);  free(ctx->VLabel);
    free(ctx->VDPtr);   free(ctx->VStrN);
    free(ctx->CmdBuf);

#ifdef TDA_R_PACKAGE
    tda_export_free(ctx);
#endif
    free(ctx);
}
