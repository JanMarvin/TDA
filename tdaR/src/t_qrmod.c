/****************************************************************************/
/*  t_qrmod                                                                 */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-97 Goetz Rohwer. All rights reserved.           */
/*                                                                          */
/*  This file is part of TDA.                                               */
/*                                                                          */
/*  TDA is free software; you can redistribute it and/or modify             */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  TDA is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program (it should be in a file named COPYING),         */
/*  if not, write to the Free Software Foundation, Inc.,                    */
/*  675 Mass Ave, Cambridge, MA 02139, USA.                                 */
/*                                                                          */

#include "tda.h"
#include "t_gen.h"   
#include "t_pgen.h"   
#include "t_parm.h"   
#include "t_gdat.h"   
#include "t_freq.h"   
#include "t_ml.h"   
#include "t_gf.h"   
#include "t_var.h"   
#include "t_cdf.h"   
#include "t_alloc.h"   
#include "t_sort.h"   
#include "t_con.h"   
#include "t_min.h"   
#include "t_fnqr.h"   
#include "tda_context.h"

/*  functions in t_qrmod.c */

int qreg(TDAContext *ctx);
int prn_qrmod(TDAContext *ctx, int mod);
int qr_check(TDAContext *ctx);
void qr_nparm(TDAContext *ctx);
int qr_alloc(TDAContext *ctx, int typ,int opt);
void qr_sval(TDAContext *ctx);
void prn_qcoeff(TDAContext *ctx, int mod);
void prn_qc1(TDAContext *ctx, int k);
void qr_options(TDAContext *ctx, int mod);
void qr_resid(TDAContext *ctx, int mod);
void qr_resida(TDAContext *ctx, int k, int iv);
int qr_svar(TDAContext *ctx, int opt);
int qr_depcat(TDAContext *ctx, int opt);
int pdmc_init(TDAContext *ctx, int opt);
void qr_pdata(TDAContext *ctx);
void qr_prob1(TDAContext *ctx, int mod,int icase,int wave,double *par,double *prob);
void qr_prob2(TDAContext *ctx, int mod,int icase,int wave,double *par,double *prob);
void qr_log3(TDAContext *ctx, int icase,int wave,double *par,double *prob);
void qr_corr(TDAContext *ctx);
int rmod(TDAContext *ctx);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */





/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */


/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/*  qreg()          Quantal response models.                                */
/*                                                                          */
/*                  qreg(                                                   */
/*                      m=...,      model number, def. 1                    */
/*                                   1 = binary logit                       */
/*                                   2 = binary probit                      */
/*                                   3 = ordinal logit                      */
/*                                   4 = oridnal probit                     */
/*                                   5 = multinomial logit                  */
/*                      nw=...,     number of waves, def. 1                 */
/*                      ni=...,      1 if without intercept, def. 0         */
/*                      maxcat=..., max number of categories, def. 1000     */
/*                      pmin                                                */
/*                      dsv=...,    starting values                         */
/*                      ppar=...,   print parameter                         */
/*                      pcov=...,   cov matrix                              */ 
/*                      tfmt=...,   print format parameter, def. 10.4       */
/*                      mfmt=...,   print format pcov, def. 10.6            */
/*                      con=...,    constraints                             */
/*                  ) = varlist;                                            */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int qreg(TDAContext *ctx)
{
    int err,nv;
           
    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Quantal response models. Current memory: %d bytes.\n",ctx->MemReq);

    set_mldef(ctx);                /* set defaults for ML estimation */
    ctx->PDMCFlg = 0;                      

    if (parm(ctx, ctx->CmdBuf + 4,4,1))   /* get parameters */
        goto QRFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMMaxCatFlg == 0)
        ctx->PMMaxCat = MaxYC;

    ctx->NWave = ctx->PMNW;

    newline(ctx);
    if ((ctx->PMM == 5 || ctx->PMM == 6) && ctx->PMNQ < 1) {
        printf1(ctx, "Error: multinomial logit and probit models needs nq parameter.\n");
        goto QRFin;
    }
    if ((ctx->PMM == 7 || ctx->PMM == 8) && ctx->NWave < 2) {
        printf1(ctx, "Error: this model needs at least two waves.\n");
        goto QRFin;
    }
    if (ctx->PMM == 6 && ctx->PMNQ > DMVMax) {
        printf1(ctx, "Error: max number of categories is %d.\n",DMVMax);
        goto QRFin;
    }
    if (ctx->PMM == 8 && ctx->NWave > DMVMax) {
        printf1(ctx, "Error: max number of waves is %d.\n",DMVMax);
        goto QRFin;
    }
    nv = check_nvar(ctx, 1);         /* check variables */
    if (nv == 0)                /* PMNV is number of variables */       
        goto QRFin;        

    if (prn_qrmod(ctx, ctx->PMM))         /* print type of model */
        goto QRFin;
    newline(ctx);
       
    prn_nwvar(ctx, 1);               /* print variables */
          
    if (ctx->NPZ > 0) {
        if (ctx->PMM != 5 && ctx->PMM != 6) {
            printf1(ctx, "Error: z variables only with m=5 or m=6.\n");
            goto QRFin;
        }
    }
    if (ctx->PMM == 7)   /* always without intercept */
        ctx->PMNI = 1;

    ctx->NPX1 = ctx->NPX;
    if (ctx->PMNI) {
        printf1(ctx, "Model without intercept.\n");
        ctx->PMNI = 1;
        ctx->IntFlg = 0;
    }
    else {
        ctx->PMNI = 0;
        ctx->IntFlg = 1;
        ctx->NPX1++;
        if (ctx->NPX + ctx->NPZ == 0)  
            printf1(ctx, "Model without independent variables.\n");
    }
    if (ctx->NPX1 + ctx->NPZ == 0) {
        printf1(ctx, "Error: no model parameters.\n");
        goto QRFin;
    }
    if (qr_svar(ctx, 1))         /* save variables */
        goto QRFin;

    ctx->NCV = pdmc_init(ctx, 1);     /* init missing value flags */
    if (ctx->NCV < 2) {
        if (ctx->NCV >= 0)
            printf1(ctx, "Error: insufficient number of cases.\n");
        goto QRFin;
    }
    newline(ctx);

    if (qr_depcat(ctx, 1))       /* distribution of categories */
        goto QRFin;

    if (ctx->NPZ > 0 && ctx->NYCat != ctx->PMNQ) {
        printf1(ctx, "Error: nq=%d inconsistent with data.\n",ctx->PMNQ);
        goto QRFin;
    }
    if (ctx->WIVar >= 0) {
        prn_cwt(ctx);          /* case weight information */
        newline(ctx);
    }
    if (qr_check(ctx))         /* check model requirements */
        goto QRFin;

    qr_nparm(ctx);             /* set NParm = number of parameters */

    set_mlopt(ctx);            /* adjust ML options */         

    if (ml_init(ctx, 1,0,0))     /* init ML estimation */
        goto QRFin;

    if (ctx->QRTyp == QRPROB3 || ctx->QRTyp == QRPROB4) {
        if (ctx->PMNHP < 2)
            ctx->PMNHP = 2;
        else if (ctx->PMNHP > 20)
            ctx->PMNHP = 20;
        else if (ctx->PMNHP == 11)
            ctx->PMNHP = 10;

        if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
            ctx->PMOPT = 1;

        printf1(ctx, "Control of integration (nhp): %d",ctx->PMNHP);
        if (ctx->PMNHP >= 12)
            printf1(ctx, " [eps=%lg]",ctx->PMEPS);
        printf1(ctx, "\nType of parameterization: %d\n",ctx->PMOPT);
    }

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto QRFin;
  
    if (ctx->DSVFlg == 0)
        qr_sval(ctx);          /* default starting values */

    prot_init(ctx, 2,ctx->QRTyp,0);   /* init protocol file */

    if (ctx->PMNConS > 0) {      /* process constraints */
        newline(ctx);
        if (con_proc(ctx, ctx->CmdBuf))
            goto QRFin;
    }
    if (ctx->NCV < ctx->NParm1) {
        printf1(ctx, "Error: insufficient number of cases for %d parameters.\n",ctx->NParm1);
        goto QRFin;
    }

    if (qr_alloc(ctx, ctx->QRTyp,1))  /* additional storage */
        goto QRFin;

    ctx->QR_INIT = 1;            /* flag first call of fn_logit(), fn_probit() */

    if (ffmin(ctx, ctx->QRTyp,0,1,1))   /* minimization */
        goto QRFin;
           
    prn_mlres(ctx, 0);           /* print results of min algorithm */

    if (ctx->QRTyp == QRLOG4)  
        printf1(ctx, "Number of cases used for likelihood: %d\n\n",ctx->NFLUsed);

    prn_qcoeff(ctx, ctx->QRTyp);      /* print estimated parameters etc. */   
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "qreg.est");
    tda_export_str_flush(ctx, "qreg.est.names");
#endif
           
    if (ctx->LConv >= 0)  
        prn_ml1res(ctx, 0);      /* additional results */
      
    if (ctx->PMF1Def)            /* print data to output file */
        qr_pdata(ctx);       

    qr_options(ctx, ctx->QRTyp);      /* additional options */

    if ((ctx->QRTyp == QRPROB3 || ctx->QRTyp == QRPROB4) && ctx->PMPPFDef)  /* write correlation matrix */
        qr_corr(ctx);       

    err = 0;

QRFin:
    ctx->QR_INIT = 0; 
    con_free(ctx);
    ml_init(ctx, 0,0,0);      
    qr_alloc(ctx, 0,0);
    qr_depcat(ctx, 0);
    pdmc_init(ctx, 0);
    qr_svar(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_qrmod(mod)      Print name of quantal response model and set QRTyp. */
/*                                                                          */
/*              1 : QRLOG1      binary logit                                */
/*              2 : QRPROB1     binary probit                               */
/*              3 : QRLOG2      ordinal logit                               */
/*              4 : QRPROB2     ordinal probit                              */
/*              5 : QRLOG3      multinomial logit                           */
/*              6 : QRPROB3     multinomial probit                          */
/*              7 : QRLOG4      conditional logit                           */
/*              8 : QRPROB4     simultaneous binary probit                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int prn_qrmod(TDAContext *ctx, int mod)
{
    printf1(ctx, "Model: ");
    switch (mod) {
        case  1:        printf1(ctx, "binary logit.\n");
                        ctx->QRTyp = QRLOG1;
                        break;
        case  2:        printf1(ctx, "binary probit.\n");
                        ctx->QRTyp = QRPROB1;
                        break;
        case  3:        printf1(ctx, "ordinal logit.\n");
                        ctx->PMNI = 1;
                        ctx->QRTyp = QRLOG2;
                        break;
        case  4:        printf1(ctx, "ordinal probit.\n");
                        ctx->PMNI = 1;
                        ctx->QRTyp = QRPROB2;
                        break;
        case  5:        printf1(ctx, "multinomial logit.\n");
                        ctx->QRTyp = QRLOG3;
                        break;
        case  6:        printf1(ctx, "multivariate probit.\n");
                        ctx->QRTyp = QRPROB3;
                        break;
        case  7:        printf1(ctx, "conditional logit.\n");
                        ctx->QRTyp = QRLOG4;
                        break;
        case  8:        printf1(ctx, "simultaneous binary probit.\n");
                        ctx->QRTyp = QRPROB4;
                        break;

        default:        ++ctx->ErrCnt;
                        printf1(ctx, "%d -- not defined.\n",mod);
                        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  qr_check()      check model requirements, model is QRTyp                */
/*                  and set PMDOPT according to derivatives.                */
/*  Return 0 if OK, -1 if error.                                            */

int qr_check(TDAContext *ctx)
{
    if (ctx->QRTyp == QRLOG1   || ctx->QRTyp == QRPROB1 || ctx->QRTyp == QRLOG4 || ctx->QRTyp == QRPROB4) {
        if (ctx->NYCat != 2) {
            printf1(ctx, "Error: model requires two categories in dependent variable.\n");
            return(-1);
        }
    }
    else if (ctx->QRTyp == QRLOG2 || ctx->QRTyp == QRPROB2) {
        if (ctx->NYCat < 3) {
            printf1(ctx, "Error: model requires three or more categories in dependent variable.\n");
            return(-1);
        }
    }
    if (ctx->QRTyp == QRPROB3 || ctx->QRTyp == QRPROB4) { /* only function values */
        ctx->PMDOPT = 0;
        if (ctx->MINA < 7 || ctx->MINA > 8)
            ctx->MINA = 8;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  qr_nparm()      set number of model parameters for QRTyp                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

void qr_nparm(TDAContext *ctx)
{
    ctx->NParm = 0;       

    switch (ctx->QRTyp) {

        case QRLOG1: 
        case QRLOG4:
        case QRPROB1:   ctx->NParm = ctx->NPX1;
                        break;
        case QRLOG2:
        case QRPROB2:   ctx->NParm = ctx->NPX + ctx->NYCat - 1;
                        break;

        case QRLOG3:    ctx->NParm = (ctx->NYCat - 1) * ctx->NPX1 + ctx->NPZ;
                        break;

        case QRPROB3:   ctx->NParm = (ctx->NYCat - 1) * ctx->NPX1 + ctx->NPZ + (ctx->NYCat * (ctx->NYCat - 1)) / 2;
                        break;

        case QRPROB4:   ctx->NParm = ctx->NWave * ctx->NPX1 + (ctx->NWave * (ctx->NWave - 1)) / 2;
                        break;

        default:        printf1(ctx, "Error: model %d not available.\n",ctx->QRTyp);
                        return;       
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_alloc(typ,opt)   If opt != 0 allocate additional storage, depending  */
/*                      on typ. Otherwise free previously allocated         */
/*                      memory.                                             */
/*                      Return 0 if OK, otherwise -1.                       */

int qr_alloc(TDAContext *ctx, int typ,int opt)
{
    int err,n;

    if (opt == 0) {
        err = 0;
        goto QRALFin;
    }
    err = -1;

    n = ctx->NYCat;
    if (n < ctx->NWave)
        n = ctx->NWave;

    if (!(ctx->QRTmp = (double *)calloc((size_t)(n + 1),sizeof(double))))  
        goto QRALFin;  
    memrq(ctx, n + 1,sizeof(double));
    ctx->QRTmpA = n + 1;

    if (typ == QRPROB3 || typ == QRPROB4)
        n *= n;                

    if (!(ctx->QRProb = (double *)calloc((size_t)(n + 1),sizeof(double))))  
        goto QRALFin;  
    memrq(ctx, n + 1,sizeof(double));
    ctx->QRProbA = n + 1;

    if (typ == QRLOG4) {       /* conditional logit */

        n = ctx->NWave + 1;

        if (!(ctx->L4_ISET = (int *)calloc((size_t)(n),sizeof(int))))  
            goto QRALFin;  
        memrq(ctx, n,sizeof(int));
        ctx->L4_ISETA = n;
       
        if (!(ctx->L4_IPTR = (int *)calloc((size_t)(n),sizeof(int))))  
            goto QRALFin;  
        memrq(ctx, n,sizeof(int));
        ctx->L4_IPTRA = n;
         
        if (!(ctx->L4_XB = (double *)calloc((size_t)(n),sizeof(double))))  
            goto QRALFin;  
        memrq(ctx, n,sizeof(double));
        ctx->L4_XBA = n;
    }
    return(0);

QRALFin:
    if (err)
        p_err(ctx, -2,1);

    if (ctx->QRProbA > 0) {
        free((char *)ctx->QRProb);
        memrq(ctx, -ctx->QRProbA,sizeof(double));
        ctx->QRProbA = 0;        
    }
    if (ctx->QRTmpA > 0) {
        free((char *)ctx->QRTmp);
        memrq(ctx, -ctx->QRTmpA,sizeof(double));
        ctx->QRTmpA = 0;        
    }
    if (ctx->L4_ISETA > 0) {
        free((char *)ctx->L4_ISET);
        memrq(ctx, -ctx->L4_ISETA,sizeof(int));
        ctx->L4_ISETA = 0;        
    }
    if (ctx->L4_IPTRA > 0) {
        free((char *)ctx->L4_IPTR);
        memrq(ctx, -ctx->L4_IPTRA,sizeof(int));
        ctx->L4_IPTRA = 0;        
    }
    if (ctx->L4_XBA > 0) {
        free((char *)ctx->L4_XB);
        memrq(ctx, -ctx->L4_XBA,sizeof(double));
        ctx->L4_XBA = 0;        
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  qr_sval()       Calculate  default starting values for model QRTyp      */

void qr_sval(TDAContext *ctx) 
{
    register int j;
    double tmp,tmp1;

    for (j = 1; j <= ctx->NParm; ++j)
        ctx->Par[j] = 0.0;

    if (ctx->QRTyp == QRLOG2 || ctx->QRTyp == QRPROB2) {
        tmp = tmp1 = 0.0;
        for (j = 0; j < ctx->NYCat; ++j)  
            tmp1 += ctx->PYFreq[j];
        for (j = 1; j < ctx->NYCat; ++j) {
            tmp += ctx->PYFreq[j - 1];
            if (tmp > 0.0 && tmp1 > tmp)
                ctx->Par[j] = rlog(ctx, (tmp1 - tmp) / tmp);           
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  Function: prn_qcoeff(mod)                                               */
/*      Print estimated coefficients for quantal response model mod.        */

#ifdef TDA_R_PACKAGE
/* stage one "Cat|Term|Variable" label for qreg.est.names, mirroring the
   printf beside each call -- never re-deriving the walk */
static void qr_explab(TDAContext *ctx, const char *cat, const char *term,
                      const char *var)
{
    char b[160];

    snprintf(b, sizeof(b), "%s|%s|%s", cat, term, var);
    tda_export_str_row(ctx, "qreg.est.names", b);
}
static void qr_explab_i(TDAContext *ctx, int cat, const char *term,
                        const char *var)
{
    char c[16];

    snprintf(c, sizeof(c), "%d", cat);
    qr_explab(ctx, c, term, var);
}
#endif

void prn_qcoeff(TDAContext *ctx, int mod)
{
    register int j,k,l,li;
    int len,nvlen,wave;

    nvlen = ctx->PMNVLEN;
    if (nvlen < 10)
        nvlen = 10;

    len = 23;
    printf1(ctx, "Idx Cat Term   ");
    printf1(ctx, "Variable ");
    prnchar(ctx, ' ',nvlen - 8,0);
    prnchar(ctx, ' ',ctx->PMTFmt1 - 5,0); printf1(ctx, "Coeff"); 
    prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Error"); 
    prnchar(ctx, ' ',ctx->PMTFmt1 - 6,0); printf1(ctx, "C/Error  Signif");
       
    k = 1;
    wave = 1;

NXT:
    if (ctx->NPX1) {
        for (j = 1; j < ctx->NYCat; ++j) {

            printf1(ctx, "\n");
            prnchar(ctx, '-',len + 3 * (ctx->PMTFmt1 + 1) + nvlen,0);

            for (l = 1; l <= ctx->NParm; ++l) {

                if (mod == QRLOG2 || mod == QRPROB2) {
                    printf1(ctx, "\n%3d   - ",k);
                    if (l < ctx->NYCat) {
                        printf1(ctx, " I %-3d Alpha %-3d ",l,l); 
                        prnchar(ctx, ' ',nvlen - 9,0);
#ifdef TDA_R_PACKAGE
                        {
                            char t[16], v[32];
                            snprintf(t, sizeof(t), "I %d", l);
                            snprintf(v, sizeof(v), "Alpha %d", l);
                            qr_explab(ctx, "-", t, v);
                        }
#endif
                    }
                    else {
                        li = l - ctx->NYCat;
                        printf1(ctx, " X     %s ",ctx->VName[ctx->PXVar[li][0]]);
                        prnchar(ctx, ' ',nvlen - (int)strlen(ctx->VName[ctx->PXVar[li][0]]),0);
#ifdef TDA_R_PACKAGE
                        qr_explab(ctx, "-", "X",
                                  ctx->VName[ctx->PXVar[li][0]]);
#endif
                    }
                }
                else {
                    printf1(ctx, "\n%3d %3d ",k,ctx->YCat[j]);
                    if (ctx->IntFlg && l == 1) {
                        if (mod == QRPROB4)
                            printf1(ctx, " W%2d",wave);
                        else
                            printf1(ctx, " I  ");
                        printf1(ctx, "   Intercept ");
                        prnchar(ctx, ' ',nvlen - 9,0);
#ifdef TDA_R_PACKAGE
                        {
                            char t[16];
                            if (mod == QRPROB4)
                                snprintf(t, sizeof(t), "W %d", wave);
                            else
                                snprintf(t, sizeof(t), "I");
                            qr_explab_i(ctx, ctx->YCat[j], t, "Intercept");
                        }
#endif
                    }
                    else {
                        li = l - ctx->IntFlg - 1;
                        if (mod == QRPROB4)
                            printf1(ctx, " W%2d",wave);
                        else
                            printf1(ctx, " X  ");
                        printf1(ctx, "   %s ",ctx->VName[ctx->PXVar[li][0]]);
                        prnchar(ctx, ' ',nvlen - (int)strlen(ctx->VName[ctx->PXVar[li][0]]),0);
#ifdef TDA_R_PACKAGE
                        {
                            char t[16];
                            if (mod == QRPROB4)
                                snprintf(t, sizeof(t), "W %d", wave);
                            else
                                snprintf(t, sizeof(t), "X");
                            qr_explab_i(ctx, ctx->YCat[j], t,
                                        ctx->VName[ctx->PXVar[li][0]]);
                        }
#endif
                    }
                }
                prn_qc1(ctx, k);
                k++;
                if (l >= ctx->NPX1 && mod != QRLOG2 && mod != QRPROB2)    
                    break;
            }
            if (mod == QRLOG2 || mod == QRPROB2)
                break;
        }
    }
   
    if (ctx->NPZ && (mod == QRLOG3 || mod == QRPROB3)) {     /* generic variables */ 
        printf1(ctx, "\n");
        prnchar(ctx, '-',len + 3 * (ctx->PMTFmt1 + 1) + nvlen,0);

        for (j = 0; j < ctx->NPZ; ++j) {
            printf1(ctx, "\n%3d   - ",k);
            printf1(ctx, " Z%-4d %s ",j + 1,ctx->VName[ctx->PZVar[j][0]]);
            prnchar(ctx, ' ',nvlen - (int)strlen(ctx->VName[ctx->PZVar[j][0]]),0);
#ifdef TDA_R_PACKAGE
            {
                char t[16];
                snprintf(t, sizeof(t), "Z%d", j + 1);
                qr_explab(ctx, "-", t, ctx->VName[ctx->PZVar[j][0]]);
            }
#endif

            prn_qc1(ctx, k);
            k++;
        }
    }
    if (mod == QRPROB4 && wave++ < ctx->NWave)
        goto NXT;

    if (mod == QRPROB3 || mod == QRPROB4) {
        printf1(ctx, "\n");
        prnchar(ctx, '-',len + 3 * (ctx->PMTFmt1 + 1) + nvlen,0);
        for (j = 2; j <= ctx->NYCat; ++j) {
            for (l = 1; l < j; ++l) {
                printf1(ctx, "\n%3d   - ",k);
                printf1(ctx, " S     Sigma%2d,%2d ",j,l);
                prnchar(ctx, ' ',nvlen - 10,0);
#ifdef TDA_R_PACKAGE
                {
                    char v[32];
                    snprintf(v, sizeof(v), "Sigma %d, %d", j, l);
                    qr_explab(ctx, "-", "S", v);
                }
#endif
                prn_qc1(ctx, k);
                k++;
            }
        }
    }
    newline(ctx);
    newline(ctx);
}

void prn_qc1(TDAContext *ctx, int k)     /* print coefficients */
{
    double tmp,tmp1,tmp2;
#ifdef TDA_R_PACKAGE
    double erow[5];
    erow[0] = (double)k;
    erow[1] = ctx->Par[k];
    erow[2] = erow[3] = erow[4] = (double)NAN;
#endif

    rt_printf1_d(ctx, ctx->PMTFmtS,ctx->Par[k]);
    if (ctx->LConv >= 0 && ctx->Diag[k] >= 0.0) {
        if (ctx->Diag[k] > ctx->EPSI1)
            tmp = sqrt(ctx->Diag[k]);
        else
            tmp = 0.0;

#ifdef TDA_R_PACKAGE
        erow[2] = tmp;
#endif
        rt_printf1_d(ctx, ctx->PMTFmtS,tmp);        
    }
    else {
        tmp = 0.0;
        prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0);
        printf1(ctx, "--- ");
    }
    if (ctx->LConv >= 0 && tmp > 0.0) {
        tmp1 = ctx->Par[k] / tmp;
        rt_printf1_d(ctx, ctx->PMTFmtS,tmp1);          
        tmp2 = 2.0 * cdnf(ctx, fabs(tmp1)) - 1.0;
        printf1(ctx, " %6.4lf",tmp2);          
#ifdef TDA_R_PACKAGE
        erow[3] = tmp1;
        erow[4] = tmp2;
#endif
    }
    else {
        prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0);
        printf1(ctx, "---     ---");         
    }
#ifdef TDA_R_PACKAGE
    /* the same row the table above prints, staged for the qreg.est
       export; --- prints as NaN (CONTRIBUTING.md) */
    tda_export_row(ctx, "qreg.est", erow, 5);
#endif
}

/* ------------------------------------------------------------------------ */
/*  qr_options(mod)     process additional options for model mod.           */

void qr_options(TDAContext *ctx, int mod)
{
    int i;

    for (i = 0; i < ctx->PMResN; ++i) {
        switch (ctx->PMRes[i]) {
            case  1:                        /* stand. residuals */
                if (ctx->NWave != 1)
                    break;
                if (mod == QRLOG1 || mod == QRPROB1)
                    qr_resid(ctx, mod);
                break;

            default:
                break;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_resid()                                                              */
/*      Print standardized coefficients, only for binary and multinomial    */
/*      logit models with cross-sectional data (NWave = 1).                 */
       
void qr_resid(TDAContext *ctx, int mod)
{
    (void)mod;        /* unused: the signature is shared */
    register int j,k,l,i,len,nvlen;

    nvlen = ctx->PMNVLEN;
    if (nvlen < 9)
        nvlen = 9;
        
    printf1(ctx, "\nStandardized coefficients.\n\n");

    len = 15;
    printf1(ctx, "Idx Cat Term   Variable ");
    prnchar(ctx, ' ',nvlen - 8,0);

    prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0); printf1(ctx, "     Coeff"); 
    prnchar(ctx, ' ',ctx->PMTFmt1 -  5,0); printf1(ctx, "Exp(C)"); 
    prnchar(ctx, ' ',ctx->PMTFmt1 -  8,0); printf1(ctx, "Exp(C*SD)");
    prnchar(ctx, ' ',ctx->PMTFmt1 -  7,0); printf1(ctx, "Std.Dev.");
    k = 1;

    if (ctx->NPX1) {
        for (j = 1; j < ctx->NYCat; ++j) {

            newline(ctx);           
            prnchar(ctx, '-',len + 4 * (ctx->PMTFmt1 + 1) + nvlen,0);

            for (l = 1; l <= ctx->NParm; ++l) {

                printf1(ctx, "\n%3d %3d ",k,ctx->YCat[j]);
                if (ctx->IntFlg && l == 1) {
                    printf1(ctx, " I     Intercept ");
                    prnchar(ctx, ' ',nvlen - 9,0);
                    qr_resida(ctx, k,-1);
                }
                else {
                    i = l - ctx->IntFlg - 1;
                    printf1(ctx, " X     %s ",ctx->VName[ctx->PXVar[i][0]]);           
                    prnchar(ctx, ' ',nvlen - (int)strlen(ctx->VName[ctx->PXVar[i][0]]),0);
                    qr_resida(ctx, k,ctx->PXVar[i][0]);
                }
                k++;
                if (l >= ctx->NPX1)
                    break;
            }
        }
    }
    /***
    if (NPZ && mod == QRLOG3) {        generic variables   
        printf1(ctx, "\n");
        prnchar(ctx, '-',len + 3 * (PMTFmt1 + 1) + nvlen,0);
        kk = k;
        for (j = 0; j < NPZ; ++j) {
            for (jj = 0; jj < PZN[j]; ++jj) {
                i = PZVar[j][jj];
                printf1(ctx, "\n%3d %3d ",kk++,YCat[jj]);
                printf1(ctx, " Z%-4d %s ",PZNum[j],VName[i]);
                prnchar(ctx, ' ',nvlen - strlen(VName[i]),0);
                qr_resida(ctx, k,i);
            }
            k++;
        }
    }
    **/
    newline(ctx);
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "qreg.standardized");
#endif
}

void qr_resida(TDAContext *ctx, int k, int iv)     /* print standardized coefficients */
{
    register int i;
    double wt,sdev,sum,ws,tmp;

    sdev = 0.0;

    if (iv >= 0) {
        wt = 1.0;  
        sum = ws = 0.0;

        for (i = 0; i < ctx->NOC; ++i) {
                     
            if (ctx->PDMCFlg && ctx->PDMCTi[i] == 0)  
                continue;

            if (!ctx->PDMCFlg || ctx->PDMCIdx[i * ctx->NWave]) {

                if (ctx->WIVar >= 0)
                    wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
                ws += wt;
                tmp = get_data(ctx, iv,i);
                sum += tmp * wt;
                sdev += tmp * tmp * wt;
            }
        }
        /*  tdaR 6.4q: this block stood inside the loop above and ran
            after every case -- subtracting the running mean and taking
            the square root, so the next iteration added tmp*tmp*wt to an
            already rooted value.  Every covariate's Std.Dev. came out 0
            and Exp(C*SD) 1 (the manual's qr1 box shows that output).
            Moved after the loop: for qr1.dat with cwt = Weight the
            weighted SD of Log10Dose is 0.27370532608834974 (the sum of
            the case weights, 74, is the divisor's n).  The Intercept row
            keeps 0: iv < 0 there, and a constant has no SD.  */
        if (ws > 0.0) {
            tmp = sum / ws;
            sdev -= tmp * tmp * ws;
            if (sdev > 0.0 && ws > 1.0)
                sdev = sqrt(sdev / (ws - 1.0));
            else
                sdev = 0.0;
        }
        else
            sdev = 0.0;
    }
    rt_printf1_d(ctx, ctx->PMTFmtS,ctx->Par[k]);
    rt_printf1_d(ctx, ctx->PMTFmtS,rexp(ctx, ctx->Par[k]));
    rt_printf1_d(ctx, ctx->PMTFmtS,rexp(ctx, ctx->Par[k] * sdev));
    rt_printf1_d(ctx, ctx->PMTFmtS,sdev);
#ifdef TDA_R_PACKAGE
    {
        double erow[4];
        erow[0] = ctx->Par[k];
        erow[1] = rexp(ctx, ctx->Par[k]);
        erow[2] = rexp(ctx, ctx->Par[k] * sdev);
        erow[3] = sdev;
        tda_export_row(ctx, "qreg.standardized", erow, 4);
    }
#endif
}

   
/* ------------------------------------------------------------------------ */
/*  qr_svar(opt)    If opt = 1, set var structure for qr models.            */
/*                  if opt = 0, free previously allocated memory.           */
/*                                                                          */
/*                  NWave               number of waves (= PMNW)            */
/*                  short PYVar[w]      indices of dependent variables      */
/*                                      w = 0,...,NWave - 1.                */
/*                  short PXVar[i][w]   indices of independent variables    */  
/*                                      i = 0,...,NPX-1, w=0,...,NWave-1    */
/*                                                                          */
/*                  short PZVar[i][w]   indices of z variables              */  
/*                                      i = 0,...,NPZ-1                     */
/*                                      w = 0,...,NWave - 1                 */
/*                                      j = 0,...,PMNQ - 1                  */
/*                                      PZVar[i][w * PMNQ + j]              */
/*                                      later: NYCat = PMNQ                 */
/*                                                                          */
/*                  Return 0 if successful, otherwise -1.                   */

int qr_svar(TDAContext *ctx, int opt)
{
    register int i,j,k,kk;
    int err;

    err = -1;
    if (opt == 0) {
        err = 0;
        goto QRSVFree;
    }
    if (!(ctx->PYVar = (short *)calloc((size_t)(ctx->NWave),sizeof(short))))  
        goto QRSVFree;
    memrq(ctx, ctx->NWave,sizeof(short));
    ctx->PYVarA = ctx->NWave;

    for (i = 0; i < ctx->NWave; ++i)     /* save indices of dep variable */
        ctx->PYVar[i] = ctx->PMVIdx[i];

    if (ctx->NPX > 0) {                  /* if there are independent variables */

        if (!(ctx->PXVar = (short **)calloc((size_t)(ctx->NPX),sizeof(short *))))  
            goto QRSVFree;
        memrq(ctx, ctx->NPX,sizeof(short *));
        ctx->PXVarA = ctx->NPX;

        for (i = 0; i < ctx->NPX; ++i) {
            if (!(ctx->PXVar[i] = (short *)calloc((size_t)(ctx->NWave),sizeof(short))))  
                goto QRSVFree;
            memrq(ctx, ctx->NWave,sizeof(short));
            ctx->PXNVarA++;           
        }
        k = ctx->NWave;
        for (i = 0; i < ctx->NPX; ++i) {
            for (j = 0; j < ctx->NWave; ++j)
                ctx->PXVar[i][j] = ctx->PMVIdx[k++];
        }
    }

    if (ctx->NPZ > 0) {                  /* if there are z variables */

        if (!(ctx->PZVar = (short **)calloc((size_t)(ctx->NPZ),sizeof(short *))))  
            goto QRSVFree;
        memrq(ctx, ctx->NPZ,sizeof(short *));
        ctx->PZVarA = ctx->NPZ;
        ctx->PZNVarAL = ctx->NWave * ctx->PMNQ;
        for (i = 0; i < ctx->NPZ; ++i) {
            if (!(ctx->PZVar[i] = (short *)calloc((size_t)(ctx->PZNVarAL),sizeof(short))))  
                goto QRSVFree;
            memrq(ctx, ctx->PZNVarAL,sizeof(short));
            ctx->PZNVarA++;           
        }
        kk = 0;
        for (i = 0; i < ctx->NPZ; ++i) {
            for (k = 0; k < ctx->NWave; ++k) {
                for (j = 0; j < ctx->PMNQ; ++j) 
                    ctx->PZVar[i][k * ctx->PMNQ + j] = ctx->PMZIdx[kk++];
            }
        }
    }
    return(0);

QRSVFree:
    if (ctx->PYVarA) {
        free((char *)ctx->PYVar);
        memrq(ctx, -ctx->NWave,sizeof(short));
        ctx->PYVarA = 0;       
    }
    if (ctx->PXVarA) {
        for (i = 0; i < ctx->PXNVarA; ++i) {
            free((char *)ctx->PXVar[i]);
            memrq(ctx, -ctx->NWave,sizeof(short));
        }
        ctx->PXNVarA = 0;       
        free((char *)ctx->PXVar);
        memrq(ctx, -ctx->PXVarA,sizeof(short *));
        ctx->PXVarA = 0;       
    }
    if (ctx->PZVarA) {
        for (i = 0; i < ctx->PZNVarA; ++i) {
            free((char *)ctx->PZVar[i]);
            memrq(ctx, -ctx->PZNVarAL,sizeof(short));
        }
        ctx->PZNVarAL = ctx->PZNVarA = 0;       
        free((char *)ctx->PZVar);
        memrq(ctx, -ctx->PZVarA,sizeof(short *));
        ctx->PZVarA = 0;       
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  Function: qr_depcat(opt)                                                */
/*                                                                          */
/*      If opt = 1 calculate and print categories of dependent variable.    */
/*      Also print weights of the categories.                               */
/*      If opt = 0 free previously allocated memory.                        */
/*                                                                          */
/*      NYCat    =   Number of categories                                   */
/*      YCat[i]  =   Array with values of categories                        */
/*                                                                          */
/*      The categories are mapped to internal cat numbers: 0,...,NYCat-1.   */
/*      Mapping is through: YCatI[j], so that: if j is the value of a       */
/*      category (in the input data), then YCatI[j] is the internal cat     */
/*      number. Therefore: If YCat[i] = j, then YCatI[j] = i, i = 0,NYCat-1 */
/*                                                                          */
/*      Generally, the access is as follows:                                */
/*                                                                          */
/*      There are NWave dependent variables with associated variable        */
/*      numbers: PYVar[l], l = 0,NWave-1. Then, for any l (0,NWave-1):      */
/*      PYVar[l] is the internal variable number associated with            */
/*      the l.th dependent variable. The value for the record icase is      */
/*      found by val = get_data(PYVar[l],icase).                            */
/*                                                                          */
/*      Then, the associated internal cat number is found by                */
/*      cat = YCatI[val].                                                   */
/*                                                                          */
/*      The global array PYFreq is build with the weighted frequencies      */
/*      of the categories of the dependent variable, pooled over waves.     */
/*                                                                          */
/*      Note: If panel missing codes are defined, these observations        */
/*      are excluded from the frequency distribution.                       */
/*                                                                          */
/*      Return 0 if successful, -1 if insuff memory, -2 other error.        */

int qr_depcat(TDAContext *ctx, int opt)
{
    register int i,j,k,l,kk;
    int err,m,yca,mycat,yfreqa,yptra,jfya,wyfreqa;
    int *yc,*yfreq,*yptr,*jfy;
    double tmp,wt,*wyfreq[MaxW],d1,d2;

    if (opt == 0) {          /* just free memory */
        err = 0;
        goto QRDEPC1Fin;
    }
    wyfreqa = yca = yfreqa = yptra = jfya = 0;
    err = -1;

    printf1(ctx, "Categories of dependent variable.\n");
    printf1(ctx, "Maximum number of categories: %d\n\n",ctx->PMMaxCat);

    if (!(yc = (int *)calloc ((size_t)(ctx->PMMaxCat + 1),sizeof(int))))  
        goto QRDEPCFin;
    memrq(ctx, ctx->PMMaxCat + 1,sizeof(int));
    yca = 1;

    if (!(yfreq = (int *)calloc ((size_t)(ctx->PMMaxCat + 1),sizeof(int))))  
        goto QRDEPCFin;
    memrq(ctx, ctx->PMMaxCat + 1,sizeof(int));
    yfreqa = 1;

    if (!(yptr = (int *)calloc ((size_t)(ctx->PMMaxCat + 1),sizeof(int))))  
        goto QRDEPCFin;
    memrq(ctx, ctx->PMMaxCat + 1,sizeof(int));
    yptra = 1;

    if (!(jfy = (int *)calloc ((size_t)(ctx->NWave * ctx->NOC + 1),sizeof(int))))  
        goto QRDEPCFin;
    memrq(ctx, ctx->NWave * ctx->NOC + 1,sizeof(int));
    jfya = 1;

    k = 0;
    for (j = 0; j < ctx->NWave; ++j) {
        l = ctx->PYVar[j];
        for (i = 0; i < ctx->NOC; ++i) { 

            if (!ctx->PDMCFlg || ctx->PDMCIdx[i * ctx->NWave + j]) {
                m = (int)get_data(ctx, l,i);
                if (m < 0) {
                    printf1(ctx, "Error: dependent variable in wave %d has negative categories.\n",j + 1);
                    err = -2;
                    goto QRDEPCFin;
                }
                jfy[++k] = m;
            }
        }
    }

    ctx->NYCat = cfreq(ctx, 1,k,jfy,ctx->PMMaxCat,yc,yfreq,yptr,0,&d1,&d2);
    
    if (ctx->NYCat <= 0) {
        printf1(ctx, "Can't create distribution of categories of dependent variables.\n");
        printf1(ctx, "Possibly too many categories, current maximum is %d.\n",ctx->PMMaxCat);
        err = -2;
        goto QRDEPCFin;
    }
    if (ctx->NYCat < 2) {
        printf1(ctx, "Error: dependent variable has less than two categories.\n");
        err = -2;
        goto QRDEPCFin;
    }

    /*  Now we have NYCat categories in the array yc. Build another array
        of categories, YCat, where the categories are sorted in ascending
        order: YCat[j], j = 0, NYCat - 1. */        
                                                    
    if (!(ctx->YCat = (int *)calloc ((size_t)(ctx->NYCat),sizeof(int))))   
        goto QRDEPCFin; 
    memrq(ctx, ctx->NYCat,sizeof(int));
    ctx->YCatA = ctx->NYCat;

    if (!(ctx->PYFreq = (double *)calloc ((size_t)(ctx->NYCat),sizeof(double))))  
        goto QRDEPCFin; 
    memrq(ctx, ctx->NYCat,sizeof(double));
    ctx->PYFreqA = ctx->NYCat;

    mycat = 0;
    for (i = 1; i <= ctx->NYCat; ++i) {
        j = yptr[i];
        ctx->YCat[i - 1] = m = yc[j];
        if (mycat < m)
            mycat = m;
    }

    /*  build inverse array YCatI to map the categories to internal
        category numbers: 0,..., NYCat - 1. */

    if (!(ctx->YCatI = (int *)calloc ((size_t)(mycat + 1),sizeof(int))))  
        goto QRDEPCFin; 
    memrq(ctx, mycat + 1,sizeof(int));
    ctx->YCatIA = mycat + 1;

    for (i = 0; i < ctx->NYCat; ++i) {
        j = ctx->YCat[i];
        ctx->YCatI[j] = i;
    }

    /*  Build wyfreq[i], weighted frequencies, separately for each wave i */
    /*  PYFreq contains the same but pooled over all waves. */

    k = 0;
    wt = 1.0;            
    k = 0;
    for (i = 0; i < ctx->NWave; ++i) {

        if (!(wyfreq[i] = (double *)calloc ((size_t)(ctx->NYCat + 1),sizeof(double))))  
            goto QRDEPCFin;
        memrq(ctx, ctx->NYCat + 1,sizeof(double));
        wyfreqa++;

        for (j = 0; j < ctx->NOC; ++j) { 

            if (!ctx->PDMCFlg || ctx->PDMCIdx[j * ctx->NWave + i]) {

                kk = jfy[++k];
                l = ctx->YCatI[kk];

                if (ctx->WIVar >= 0)                      /* get weights */
                    wt = get_data(ctx, ctx->WIVar,j) * ctx->WNorm;
       
                wyfreq[i][l] += wt;
                ctx->PYFreq[l] += wt;
            }
        }
    }
    printf1(ctx, "Index      ");
    for (j = 0; j < ctx->NYCat; ++j)  
        printf1(ctx, "  %8d",j);
       
    printf1(ctx, "   (Weighted)\nCategory   ");

    for (j = 0; j < ctx->NYCat; ++j)  
        printf1(ctx, "  %8d",ctx->YCat[j]);
       
    printf1(ctx, "  Observations\n");
                 
    prnchar(ctx, '-',25 + ctx->NYCat * 10,0);
    for (i = 0; i < ctx->NWave; ++i) {
        wt = 0.0;
        for (j = 0; j < ctx->NYCat; ++j)  
            wt += wyfreq[i][j];


        printf1(ctx, "\nWave %-3d N ",i + 1);
        for (j = 0; j < ctx->NYCat; ++j) { 
            printf1(ctx, "%10.2lf",wyfreq[i][j]);
#ifdef TDA_R_PACKAGE
            {
                double erow[4];
                erow[0] = (double)(i + 1);
                erow[1] = (double)ctx->YCat[j];
                erow[2] = wyfreq[i][j];
                erow[3] = wt > 0.0 ? 100.0 * wyfreq[i][j] / wt : 0.0;
                tda_export_row(ctx, "qreg.categories", erow, 4);
            }
#endif
        }
             
        printf1(ctx, " %13.2lf\n         Pct ",wt);
        for (j = 0; j < ctx->NYCat; ++j) {
            if (wt > 0.0)
                tmp = 100.0 * wyfreq[i][j] / wt;
            else
                tmp = 0.0;

            printf1(ctx, "%8.2lf  ",tmp);
        }
    }
    newline(ctx);                     
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "qreg.categories");
#endif
    err = 0;
    newline(ctx);

QRDEPCFin:
    if (yca) {
        free((char *)yc);
        memrq(ctx, -ctx->PMMaxCat - 1,sizeof(int));
    }
    if (yfreqa) {
        free((char *)yfreq);
        memrq(ctx, -ctx->PMMaxCat - 1,sizeof(int));
    }
    if (yptra) {
        free((char *)yptr);
        memrq(ctx, -ctx->PMMaxCat - 1,sizeof(int));
    }
    if (jfya) {
        free((char *)jfy);
        memrq(ctx, -ctx->NWave * ctx->NOC - 1,sizeof(int));
    }
    for (i = 0; i < wyfreqa; ++i) {
        free((char *)wyfreq[i]);
        memrq(ctx, -ctx->NYCat - 1,sizeof(double));
    }

QRDEPC1Fin:
    if (err || opt == 0) {
        if (ctx->YCatA) {
            free((char *)ctx->YCat);
            memrq(ctx, -ctx->YCatA,sizeof(int));
            ctx->YCatA = 0;
        }
        if (ctx->YCatIA) {
            free((char *)ctx->YCatI);
            memrq(ctx, -ctx->YCatIA,sizeof(int));
            ctx->YCatIA = 0;
        }
        if (ctx->PYFreqA) {
            free((char *)ctx->PYFreq);
            memrq(ctx, -ctx->PYFreqA,sizeof(double));
            ctx->PYFreqA = 0;
        }
    }
    if (err == -1)
        p_err(ctx, -2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  pdmc_init   If opt = 1 init data structure for panel miss. values.      */
/*              If opt = 0 free previously allocated memory.                */
/*                                                                          */
/*              PDMCIdx of length NOC * NWave is set to 1 if an observation */
/*                      is available                                        */
/*              PDMCTi  of length NOC is the number of waves an individual  */
/*                      has participated.                                   */
/*                                                                          */
/*              A case is missing if the dependent variable has a negative  */
/*              value.                                                      */
/*              PMPMin = min number of waves for participation (pmin)       */
/*              PDMCTi is only set to a nonzero value if an individual has  */
/*              participated at least PMPMin waves with valid values.       */
/*              Otherwise PDMCTi is zero and also PDMCIdx is zero for       */
/*              all waves.                                                  */
/*                                                                          */
/*              If there is at least one missing value then PMMCFlg = 1,    */  
/*              otherwise PDMCFlg = 0.                                      */
/*                                                                          */
/*              Return  number of valid cases.                              */
/*                     -1  if insufficient memory                           */

int pdmc_init(TDAContext *ctx, int opt)
{
    register int i,j;
    int err,n,nn,ni;

    if (opt == 0) {
        err = 0;
        goto PDMCFin;
    }
    printf1(ctx, "Checking available data (pmin=%d)\n",ctx->PMPMin);

    err = -1;
    ctx->PDMCFlg = 0;         

    if (!(ctx->PDMCIdx = (char *) calloc((size_t)(ctx->NOC) * (size_t)(ctx->NWave) + 1,sizeof(char)))) {
        p_err(ctx, -2,1);
        goto PDMCFin;
    }
    ctx->PDMCIdxA = ctx->NOC * ctx->NWave + 1;
    memrq(ctx, ctx->PDMCIdxA,sizeof(char));

    if (!(ctx->PDMCTi = (short *) calloc((size_t)(ctx->NOC + 1),sizeof(short)))) {
        p_err(ctx, -2,1);
        goto PDMCFin;
    }
    ctx->PDMCTiA = ctx->NOC + 1;
    memrq(ctx, ctx->PDMCTiA,sizeof(short));

    nn = ni = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        ctx->PDMCTi[i] = 0;
        n = 0;
        for (j = 0; j < ctx->NWave; ++j) {
            if (get_data(ctx, ctx->PMVIdx[j],i) >= 0) {
                ctx->PDMCIdx[i * ctx->NWave + j] = 1;
                n++;
            }
        }
        if (n >= ctx->PMPMin) {
            ctx->PDMCTi[i] = (short)(n);
            ni++;
            nn += n;
        }
        else
            ctx->PDMCFlg = 1;
    }
    printf1(ctx, "Number of cases with valid data: %d\n",ni);
    if (ctx->NWave > 1)
        printf1(ctx, "Number of valid observations: %d\n",nn);

    return(ni);

PDMCFin:
    if (err || opt == 0) {
        if (ctx->PDMCIdxA) {
            free(ctx->PDMCIdx);
            memrq(ctx, -ctx->PDMCIdxA,sizeof(char));
            ctx->PDMCIdxA = 0;                 
        } 
        if (ctx->PDMCTiA) {
            free((char *)ctx->PDMCTi);
            memrq(ctx, -ctx->PDMCTiA,sizeof(short));
            ctx->PDMCTiA = 0;                 
        } 
        ctx->PDMCFlg = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  qr_pdata()      Print data to PMF1d.                                    */
/*                  Use PMFmtS.                                             */
/*                                                                          */
/*                  Return 0 if successful, otherwise -1.                   */

void qr_pdata(TDAContext *ctx)
{
    register int i,j,k,l;
    int r,n,nrec,pflag = 0;
    double wt,prob;

    if (ctx->PMF1Def == 0)
        return;

    nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) { 
 
        if (ctx->PDMCFlg && ctx->PDMCTi[i] == 0)  
            continue;

        for (j = 0; j < ctx->NWave; ++j) {

            if (!ctx->PDMCFlg || ctx->PDMCIdx[i * ctx->NWave + j]) {

                fprintf(ctx->PMF1d,"%6d %2d ",i + 1,j + 1);
#ifdef TDA_R_PACKAGE
                /* the prediction row, cell by cell: its width depends on
                   the covariate count, the z-variable block, whether
                   case weights were given, and how many response
                   categories the model has */
                tda_export_cell(ctx, "qreg.predictions", (double)(i + 1));
                tda_export_cell(ctx, "qreg.predictions", (double)(j + 1));
#endif

                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,get_data(ctx, ctx->PYVar[j],i));
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "qreg.predictions",
                                get_data(ctx, ctx->PYVar[j],i));
#endif
                for (k = 0; k < ctx->NPX; ++k) {
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,get_data(ctx, ctx->PXVar[k][j],i));
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "qreg.predictions",
                                    get_data(ctx, ctx->PXVar[k][j],i));
#endif
                }
   
                for (k = 0; k < ctx->NPZ; ++k) {
                    for (l = 0; l < ctx->NYCat; ++l) {
                        n = ctx->PZVar[k][j * ctx->NYCat + l];
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,get_data(ctx, n,i));
#ifdef TDA_R_PACKAGE
                        tda_export_cell(ctx, "qreg.predictions",
                                        get_data(ctx, n,i));
#endif
                    }
                }
   
                if (ctx->WIVar >= 0) {   
                    wt = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,wt);
#ifdef TDA_R_PACKAGE
                    tda_export_cell(ctx, "qreg.predictions", wt);
#endif
                }
                pflag = 1;

                switch (ctx->QRTyp) {
                    case QRLOG1:
                    case QRPROB1:   qr_prob1(ctx, ctx->QRTyp,i,j,ctx->Par,ctx->QRTmp);
                                    break;
                    case QRLOG2:
                    case QRPROB2:   qr_prob2(ctx, ctx->QRTyp,i,j,ctx->Par,ctx->QRTmp);
                                    break;

                    case QRLOG3:    qr_log3(ctx, i,j,ctx->Par,ctx->QRTmp);
                                    break;

                    case QRPROB3:   for (k = -1; k < ctx->NYCat; ++k) {
                                        r = fn_prob3(ctx, i,j,k,ctx->Par,&prob);
                                        if (r)
                                            prob = -1.0;
                                        ctx->QRTmp[k + 1] = prob;
                                    }
                                    break;
                    default:        pflag = 0;
                                    break;
                }
                if (pflag) {
                    /* braces: the export call must sit INSIDE this loop */
                    for (k = 0; k <= ctx->NYCat; ++k) {
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->QRTmp[k]);
#ifdef TDA_R_PACKAGE
                        tda_export_cell(ctx, "qreg.predictions",
                                        ctx->QRTmp[k]);
#endif
                    }
                }
                fprintf(ctx->PMF1d,"\n");
#ifdef TDA_R_PACKAGE
                tda_export_endrow(ctx, "qreg.predictions");
#endif
                nrec++;
            }
        }
    }
    printf1(ctx, "Data (%d records) written to: %s\n",nrec,ctx->PMF1dName);

    if (ctx->PMTDAFDef) {
        fprintf(ctx->PMTDAFd,"# data written by qreg command.\n");
        fprintf(ctx->PMTDAFd,"nvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMF1dName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nrec);
        k = 0;
        fprintf(ctx->PMTDAFd,"  CaseID");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 6,0);
        fprintf(ctx->PMTDAFd," [6.0] = c%-2d,\n",++k);

        fprintf(ctx->PMTDAFd,"  Wave");
        fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 4,0);
        fprintf(ctx->PMTDAFd," [2.0] = c%-2d,\n",++k);

        for (j = -1; j < ctx->NPX; ++j) {
            if (j < 0)
                i = ctx->PYVar[0];
            else if (j < ctx->NPX)
                i = ctx->PXVar[j][0];

            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[i]);
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[i]),0);
            fprintf(ctx->PMTDAFd," [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++k);
        }
        for (j = 0; j < ctx->NPZ; ++j) {
            for (l = 0; l < ctx->NYCat; ++l) {
                i = ctx->PZVar[j][l];
                fprintf(ctx->PMTDAFd,"  %s",ctx->VName[i]);
                fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[i]),0);
                fprintf(ctx->PMTDAFd," [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++k);
            }
        }
        if (ctx->WIVar >= 0) {
            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[ctx->WIVar]);
            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[ctx->WIVar]),0);
            fprintf(ctx->PMTDAFd," [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++k);
        }
        if (pflag) {
            for (j = 0; j <= ctx->NYCat; ++j) {
                if (j == 0)
                    fprintf(ctx->PMTDAFd,"  PROB  ");
                else
                    fprintf(ctx->PMTDAFd,"  PROB%-2d",j - 1);
                fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - 6,0);
                fprintf(ctx->PMTDAFd," [%d.%d] = c%-2d,\n",ctx->PMFmt1,ctx->PMFmt2,++k);
            }
        }
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  qr_prob1(mod,icase,wave,par,prob)                                       */
/*                                                                          */
/*  mod is QRLOG1 or QRPROB1. Return probability for icase and wave:        */
/*  prob[0] = prob for actual choice                                        */
/*  prob[1] = prob for Y = 0                                                */
/*  prob[2] = prob for Y = 1                                                */

void qr_prob1(TDAContext *ctx, int mod,int icase,int wave,double *par,double *prob)
{
    register int j,k;
    double xb,tmp;

    tmp = 0.0;
    if (ctx->IntFlg)
        xb = par[1];
    else
        xb = 0.0;
    for (j = ctx->IntFlg; j < ctx->NParm; ++j) {
        k = j - ctx->IntFlg;
        xb += par[j + 1] * get_data(ctx, ctx->PXVar[k][wave],icase);
    }
    if (mod == QRLOG1) {            /* binary logit */
        tmp = rexp(ctx, xb);
        tmp /= (tmp + 1.0);
    }
    else if (mod == QRPROB1)        /* binary probit */
        tmp = cdnf(ctx, xb);

    prob[2] = tmp;
    prob[1] = 1.0 - tmp;

    k = (int) ctx->PYVar[wave];         /* dependent variable */

    j = ctx->YCatI[(int)get_data(ctx, k,icase)];   /* internal category */
    if (j == 0)  
        prob[0] = prob[1]; 
    else
        prob[0] = prob[2];
}

/* ------------------------------------------------------------------------ */
/*  qr_prob2(mod,icase,wave,par,prob)                                       */
/*                                                                          */
/*  mod is QRLOG2 or QRPROB2. Return probabilitiesfor icase and wave        */
/*  in prob[].                                                              */

void qr_prob2(TDAContext *ctx, int mod,int icase,int wave,double *par,double *prob)
{
    register int i,j,k;
    double xb,tmp;

    tmp = 0.0;

    k = (int) ctx->PYVar[wave];                /* dependent variable */
    k = ctx->YCatI[(int)get_data(ctx, k,icase)];    /* internal category */
    
    j = 0;
    xb = 0.0;
    for (i = ctx->NYCat; i <= ctx->NParm; ++i)   
        xb += par[i] * get_data(ctx, ctx->PXVar[j++][wave],icase);
    
    for (k = 1; k < ctx->NYCat; ++k) {
        if (mod == QRLOG2) {
            tmp = rexp(ctx, par[k] + xb);
            prob[k] = tmp / (1.0 + tmp);
        }
        else
            prob[k] = cdnf(ctx, par[k] + xb);
    }
    prob[ctx->NYCat] = 0.0;

    k = (int) ctx->PYVar[wave];         /* dependent variable */

    j = ctx->YCatI[(int)get_data(ctx, k,icase)];   /* internal category */
    prob[0] = prob[j + 1];
}

/* ------------------------------------------------------------------------ */
/*  qr_log3(icase,wave,par,prob)                                            */
/*                                                                          */
/*  Return probabilities for icase and wave based on                        */
/*  multinomial logit model. Prob for actual choice in prob[0].             */
/*                                                                          */

void qr_log3(TDAContext *ctx, int icase,int wave,double *par,double *prob)
{
    register int i,j,l,u,k1,k2;
    double tmp,d;  

    d = 1.0;        /* d is the inclusive value */
    u = 1;
    for (j = 2; j <= ctx->NYCat; ++j) {   

        if (ctx->IntFlg)
            tmp = par[u++];
        else
            tmp = 0.0;

        for (i = 0; i < ctx->NPX; ++i)  
            tmp += par[u++] * get_data(ctx, ctx->PXVar[i][wave],icase);
           
        l = (ctx->NYCat - 1) * ctx->NPX1 + 1;
    
        for (i = 0; i < ctx->NPZ; ++i) {
            k1 = ctx->PZVar[i][wave * ctx->NYCat + j - 1];
            k2 = ctx->PZVar[i][wave * ctx->NYCat];
            tmp += par[l++] * (get_data(ctx, k1,icase) - get_data(ctx, k2,icase));
        }
        prob[j] = rexp(ctx, tmp);
        d += prob[j];
    }
    prob[1] = 1.0;
    for (i = 2; i <= ctx->NYCat; ++i) {   
        prob[i] /= d;
        prob[1] -= prob[i];
    }
    u = (int) ctx->PYVar[wave];         /* dependent variable */

    j = ctx->YCatI[(int)get_data(ctx, u,icase)];   /* internal category */
    prob[0] = prob[j + 1];
}

/* ------------------------------------------------------------------------ */
/*  qr_corr()   print correlation matrix for multinomial probit model.      */

void qr_corr(TDAContext *ctx)
{
    register int k,l,k1,l1,kk;
    int m,jptr,nn;
    double s,ss[DMVMax + 1];
         
    if (ctx->QRTyp == QRPROB3) {
        nn = ctx->NYCat;
        jptr = (nn - 1) * ctx->NPX1 + ctx->NPZ;    /* points to sigma in par */
    }
    else {
        nn = ctx->NWave;
        jptr = nn * ctx->NPX1; 
    }
    ctx->Par[0] = 1.0;
    for (k = 1; k <= nn; ++k) {
        for (l = 1; l <= nn; ++l) {
                    
            if (ctx->PMOPT == 3) {
                if (k == 1)  
                    k1 = 0;
                else  
                    k1 = jptr + ((k - 1) * (k - 2) / 2) + 1;
                if (l == 1)  
                    l1 = 0;
                else  
                    l1 = jptr + ((l - 1) * (l - 2) / 2) + 1;

                m = imin(ctx, k,l);
                s = 0.0;
                for (kk = 1; kk < m; ++kk)  
                    s += ctx->Par[l1++] * ctx->Par[k1++];
                if (k > m)
                    s += ctx->Par[k1];
                else if (l > m)
                    s += ctx->Par[l1];
                else
                    s += 1.0;
            }
            else {
                if (k == l)
                    s = 1.0;
                else {
                    if (l < k)      
                        s = ctx->Par[jptr + ((k - 1) * (k - 2)) / 2 + l];
                    else               
                        s = ctx->Par[jptr + ((l - 1) * (l - 2)) / 2 + k];

                    if (ctx->PMOPT == 2) {
                        s = rexp(ctx, s);
                        s = 2.0 * s / (1.0 + s) - 1.0;
                    }
                }
            }
            ctx->QRProb[(k - 1) * nn + l] = s;           
        }
    }
    for (k = 1; k <= nn; ++k) {
        s = ctx->QRProb[(k - 1) * nn + k];
        if (s <= 0.0) {
            printf1(ctx, "Warning: estimated correlation matrix not positive definite.\n");
            return;
        }
        ss[k] = sqrt(s);
    }
    for (k = 1; k <= nn; ++k) {
        fprintf(ctx->PMPPFd,"# ");
        for (l = 1; l <= nn; ++l) {
            s = ctx->QRProb[(k - 1) * nn + l] / (ss[k] * ss[l]);     
            rt_fprintf_d(ctx, ctx->PMPPFd,ctx->PMTFmtS,s);
        }
        fprintf(ctx->PMPPFd,"\n");
    }
    printf1(ctx, "Estimated correlation matrix written to: %s\n",ctx->PMPPFName);
}

/* ------------------------------------------------------------------------ */
/*  rmod()          Rasch model.                                            */
/*  ##                                                                      */
/*                  rmod(                                                   */
/*                      m=...,      model number, def. 1                    */
/*                                   1 = simple Rasch model                 */
/*                      maxcat=..., max number of categories, def. 1000     */
/*                      dsv=...,    starting values                         */
/*                      ppar=...,   print parameter                         */
/*                      pcov=...,   cov matrix                              */ 
/*                      tfmt=...,   print format parameter, def. 10.4       */
/*                      mfmt=...,   print format pcov, def. 10.6            */
/*                      df=...,     write combinations to output file       */
/*                  ) = varlist;                                            */
/*                                                                          */
/*                  Return 0 if OK, otherwise -1.                           */

int rmod(TDAContext *ctx)
{
    register int i,j,k;
    int err;
    double d1,d2;
           
    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    printf1(ctx, "Rasch models. Current memory: %d bytes.\n",ctx->MemReq);

    set_mldef(ctx);                /* set defaults for ML estimation */

    if (parm(ctx, ctx->CmdBuf + 4,4,1))   /* get parameters */
        goto RMFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMMaxCatFlg == 0)
        ctx->PMMaxCat = MaxYC;

    printf1(ctx, "Number of variables: %d\n",ctx->PMNV);

    /* make frequency distribution from data */
    /* max number of categories is: PMMaxCat */

    ctx->RMPV = ctx->PMNV;

    if (alloc_acr(ctx, ctx->NOC * ctx->RMPV + 1))          /* raw data */
        goto RMFin;
    if (alloc_acs(ctx, ctx->PMMaxCat * ctx->RMPV + 1))     /* data combinations */
        goto RMFin;
    if (alloc_aci(ctx, ctx->PMMaxCat + 1))            /* frequencies */
        goto RMFin;
    if (alloc_acj(ctx, ctx->PMMaxCat + 1))            /* bf pointer */
        goto RMFin;

    for (i = 0; i < ctx->NOC; ++i) {
        for (j = 0; j < ctx->RMPV; ++j) {
            k = (int)get_data(ctx, ctx->PMVIdx[j],i);
            if (k < 0 || k > 1) {
                printf1(ctx, "Error: all values must be 0 or 1.\n");
                goto RMFin;
            }
            ctx->AcR[i * ctx->RMPV + j + 1] = k;
        }
    }
    ctx->RMPN = cfreq(ctx, ctx->RMPV,ctx->NOC,ctx->AcR,ctx->PMMaxCat,ctx->AcS,ctx->AcI,ctx->AcJ,0,&d1,&d2);

    if (ctx->RMPN <= 0) {
        printf1(ctx, "Can't create distribution of categories of variables.\n");
        printf1(ctx, "Possibly too many categories, current maximum is %d.\n",ctx->PMMaxCat);
        goto RMFin;
    }
    printf1(ctx, "Number of patterns: %d\n\n",ctx->RMPN);

    /* save in RMPAT, RMPF */

    if (!(ctx->RMPAT = (short *)calloc((size_t)(ctx->RMPN) * (size_t)(ctx->RMPV) + 1,sizeof(short))))  
        goto RMFin;  
    ctx->RMPATA = ctx->RMPN * ctx->RMPV + 1;
    memrq(ctx, ctx->RMPATA,sizeof(short));

    if (!(ctx->RMPF = (int *)calloc((size_t)(ctx->RMPN + 1),sizeof(int))))  
        goto RMFin;  
    ctx->RMPFA = ctx->RMPN + 1;
    memrq(ctx, ctx->RMPFA,sizeof(int));

    for (i = 1; i <= ctx->RMPN; ++i) {
        k = ctx->AcJ[i];
        for (j = 1; j <= ctx->RMPV; ++j) {
            ctx->RMPAT[(i - 1) * ctx->RMPV + j] = (short)(ctx->AcS[(k - 1) * ctx->RMPV + j]);

            if (ctx->PMF1Def)  
                fprintf(ctx->PMF1d,"%d ",ctx->RMPAT[(i-1) * ctx->RMPV + j]);
        }
        ctx->RMPF[i] = ctx->AcI[k];

        if (ctx->PMF1Def)  
            fprintf(ctx->PMF1d,"%6d\n",ctx->RMPF[i]);
    }
    ctx->NParm = ctx->RMPV + ctx->RMPN;

    alloc_acr(ctx, 0);
    alloc_acs(ctx, 0);
    alloc_aci(ctx, 0);
    alloc_acj(ctx, 0);




    ctx->PMDOPT = 0;
    if (ctx->MINA == 5 || ctx->MINA == 6)
        ctx->MINA = 4;
    else if (ctx->MINA == 7 || ctx->MINA > 8)
        ctx->MINA = 8;

    set_mlopt(ctx);            /* adjust ML options */         
            
    if (ml_init(ctx, 1,0,0))     /* init ML estimation */
        goto RMFin;

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto RMFin;
  
/*  if (DSVFlg == 0)
        qr_sval(ctx);     */   /* default starting values */

    prot_init(ctx, 2,RMOD1,0);   /* init protocol file */

    /**
    if (RMPN < NParm1) {
        printf1(ctx, "Error: insufficient number of cases for %d parameters.\n",NParm1);
        goto RMFin;
    }
    **/

    if (ffmin(ctx, RMOD1,0,1,1))   /* minimization */
        goto RMFin;
    /* an unconditional goto RMFin here skipped everything below --
       result printing, the additional options, and the err = 0 that
       marks success -- so rmod estimated and then discarded the
       estimates.  Removed so the documented output is produced. */
    prn_mlres(ctx, 0);           /* print results of min algorithm */

#ifdef TDA_R_PACKAGE
    {   /* the estimated parameters: RMPV item parameters followed by
           RMPN pattern parameters, in Par[1..NParm] */
        int rmi;
        for (rmi = 1; rmi <= ctx->NParm; ++rmi)
            tda_export_cell(ctx, "rmod.par", ctx->Par[rmi]);
        tda_export_endrow(ctx, "rmod.par");
        tda_export_flush(ctx, "rmod.par");
    }
#endif

    if (ctx->QRTyp == QRLOG4)  
        printf1(ctx, "Number of cases used for likelihood: %d\n\n",ctx->NFLUsed);

    prn_qcoeff(ctx, ctx->QRTyp);      /* print estimated parameters etc. */   
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "qreg.est");
    tda_export_str_flush(ctx, "qreg.est.names");
#endif
           
    if (ctx->LConv >= 0)  
        prn_ml1res(ctx, 0);      /* additional results */
      
    if (ctx->PMF1Def)            /* print data to output file */
        qr_pdata(ctx);       

    qr_options(ctx, ctx->QRTyp);      /* additional options */

    if ((ctx->QRTyp == QRPROB3 || ctx->QRTyp == QRPROB4) && ctx->PMPPFDef)  /* write correlation matrix */
        qr_corr(ctx);       

    err = 0;

RMFin:
    if (ctx->RMPATA > 0) {
        free((char *)ctx->RMPAT);
        memrq(ctx, -ctx->RMPATA,sizeof(short));
        ctx->RMPATA = 0;        
    }
    if (ctx->RMPFA > 0) {
        free((char *)ctx->RMPF);
        memrq(ctx, -ctx->RMPFA,sizeof(int));
        ctx->RMPFA = 0;        
    }
    ml_init(ctx, 0,0,0);      
    p_clean(ctx);
    return(err);
}


