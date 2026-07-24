/****************************************************************************/
/*  t_con                                                                   */
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
#include "t_ml.h"
#include "t_alloc.h"
#include "t_lsei.h"
#include "tda_context.h"

/*  functions in t_con.c */

int p_con(TDAContext *ctx, char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni);
int get_con(TDAContext *ctx, char *pcmd);
int con_proc(TDAContext *ctx, char *pcmd);
int con_build(TDAContext *ctx);
void con_free(TDAContext *ctx);

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/*  p_con(pcmd,nx,nif,mw,nw,w,ne,ni)                                        */
/*                                                                          */
/*  get constraints: lsecon, lsicon from string pcmd.                       */
/*  nx = number of b parameters. if nif = 1 then without intercept.         */
/*  put lsecon constraints at begin, lsicon constraints at the end of       */
/*  data matrix w[] which has mw rows and nw columns.                       */
/*  return *ne = number of lsecon, *ni = number of lsicon constraints.      */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

int p_con(TDAContext *ctx, char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni)
{
    int nb,eflag,fflag,err,ir,ctyp,nterms,only_nb;
    register char *p;
    double tmp;

    p = pcmd;         

    *ne = 0;
    *ni = 0;
    err = 0;
    /* PATCH (tdaR, 2026-08): track, per lsecon=/lsicon=, whether it
       solely fixes one parameter (con=b3=1, not con=b1+b2=1) -- used
       by prn1_coeff() (t_pgen.c) to show "---" for that parameter's
       own standard error instead of whatever floating-point noise the
       constrained-optimisation covariance computation happened to
       leave there, rather than an absolute magnitude threshold, which
       cannot tell that noise apart from a genuinely tiny but real
       variance (examples/exam/lsreg3.cf's Y2 has a real variance of 6.3e-19 and is not constrained). Rank-deficient fits are caught
       separately, upstream (ranka != nx1-ranke forces df = 0, making
       this whole branch unreachable), so a constraint is the only
       other reason this value could be near zero. */
    /* rebuilt per command rather than allocated once: the previous
       command's flags must not leak into this one (an unconstrained fit
       would inherit its "---" markers), and a later fit with more
       parameters would index past an array sized for the earlier one.
       tda_context_free() releases the last one. */
    free((char *)ctx->ParFixed);
    ctx->ParFixed = NULL;
    if (ctx->NCONSTR > 0)
        /* nb below is 1-indexed and shifted up by one whenever nif == 0
           (with intercept), so it can reach nx + 1, one past the nx
           the nb > nx bounds check itself allows through before that
           shift -- sized nx + 2 here, not nx + 1, to safely hold that;
           confirmed the hard way, a real crash (malloc(): invalid
           size) on examples/exam/glm6.cf's own con=b9=1 with the
           smaller size. */
        ctx->ParFixed = (int *)calloc((size_t)(nx + 2),sizeof(int));
    while (*p) {

        ctyp = 0;
        if (!strncmp(p,"lsecon",6)) {     
            printf1(ctx, "LSECon: ");
            ctyp = 1;
            *ne += 1;
            ir = *ne;
            p += 6;
        }
        else if (!strncmp(p,"lsicon",6)) {     
            printf1(ctx, "LSICon: ");
            ctyp = 2;
            *ni += 1;
            ir = mw + 1 - *ni;
            p += 6;
        }
        if (ctyp) {
            if (*p++ != '=')
                goto PCONFin;

            err = -1;
            fflag = eflag = 0;
            nterms = 0;
            only_nb = 0;
            while (*p) {
                if (sscanf(p,"%lg",&tmp) != 1) {

                    if (*p == 'b')
                        tmp = 1.0;

                    else if (*p == '+' && *(p + 1) == 'b') {
                        tmp = 1.0;
                        p++;
                    }
                    else if (*p == '-' && *(p + 1) == 'b') {
                        tmp = -1.0;
                        p++;
                    }
                    else
                        goto PCONFin;
                }
                else
                    p = skip_dbl(ctx, p);

                printf1(ctx, "%g ",tmp);
                if (eflag) {
                    w[(ir - 1) * nw + nw] = tmp;
                    break;
                }
                if (*p == '*')
                    p++;

                if (sscanf(p,"b%d",&nb) != 1) 
                    goto PCONFin;
                p = skip_int(ctx, p + 1);

                printf1(ctx, "b%d ",nb);
                if (nb < nif || nb > nx)
                    fflag = 1; 
                else {
                    if (nif == 0)
                        nb++;
                    w[(ir - 1) * nw + nb] = tmp;
                    nterms++;
                    only_nb = nb;
                }
                if (*p != '+' && *p != '-' && *p != '=')
                    goto PCONFin;

                if (*p != '-')  
                    printf1(ctx, "%c ",*p);
                if (*p == '=') {
                    p++;
                    eflag = 1;
                }
            }
            if (eflag == 0 || (*p != ',' && *p != ')')) 
                goto PCONFin;
 
            if (fflag) {
                err = -2;
                goto PCONFin;
            }
            if (nterms == 1 && ctx->ParFixed)
                ctx->ParFixed[only_nb] = 1;
            newline(ctx);
        }
        p++;
    }
    err = 0;  

PCONFin:
    if (err == -1)   
        printf1(ctx, "\nSyntax error.\n");
    else if (err == -2)  
        printf1(ctx, "\nError in parameter index.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_con(pcmd)                                                           */
/*  ##                                                                      */
/*  Get constraints: con=... from string pcmd.                              */
/*  There should be PMNConS con=... strings in pcmd (set in parm()).        */
/*  It is assumed that there are NParm parameters, corresponding to         */
/*  b1,b2,... in the definition of constraints.                             */
/*                                                                          */
/*  Save constraints in ConR and ConD.                                      */
/*                                                                          */
/*  return 0 if OK, -1 if insufficient memory, -2 if syntax error.          */

int get_con(TDAContext *ctx, char *pcmd)
{
    int nb,eflag,fflag,err,first;
    register char *p,*pp;
    double tmp;

    err = -1;
    if (!(ctx->ConR = (double *)calloc((size_t)(ctx->NParm) * (size_t)(ctx->PMNConS) + 1,sizeof(double)))) {
        p_err(ctx, -2,1);
        goto GCONFin;
    }
    ctx->ConRA = ctx->NParm * ctx->PMNConS + 1;
    memrq(ctx, ctx->ConRA,sizeof(double));

    if (!(ctx->ConD = (double *)calloc((size_t)(ctx->PMNConS + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto GCONFin;
    }
    ctx->ConDA = ctx->PMNConS + 1;
    memrq(ctx, ctx->ConDA,sizeof(double));

    p = pcmd;         

    ctx->NCon = 0;           /* number of constraints */
    err = -2;
    while (*p) {

        if (!strncmp(p,"con=",4)) {     
            pp = p;
            first = 1;
            p += 4;

            if (++ctx->NCon > ctx->PMNConS)
                goto GCONFin;

            fflag = eflag = 0;
            while (*p) {
                if (sscanf(p,"%lg",&tmp) != 1) {
                    if (*p == 'b')
                        tmp = 1.0;
                    else if (*p == '+') {
                        tmp = 1.0;
                        p++;
                    }
                    else if (*p == '-') {
                        tmp = -1.0;
                        p++;
                    }
                    else
                        goto GCONFin;
                }
                else
                    p = skip_dbl(ctx, p);

                if (first) {
                    printf1(ctx, "Con: ");
                    first = 0;
                }
                printf1(ctx, "%g",fabs(tmp));
                if (eflag) {
                    ctx->ConD[ctx->NCon] = tmp;
                    break;
                }
                if (*p == '*')
                    p++;

                if (sscanf(p,"b%d",&nb) != 1)                          
                    goto GCONFin;

                if (nb < 1 || nb > ctx->NParm) {
                    fflag = 1;
                    goto GCONFin;
                }
                p = skip_int(ctx, p + 1);

                printf1(ctx, " * b%d ",nb);
                ctx->ConR[(ctx->NCon - 1) * ctx->NParm + nb] = tmp;

                if (*p != '+' && *p != '-' && *p != '=')
                    goto GCONFin;

                printf1(ctx, "%c ",*p);
                if (*p == '=') {
                    p++;
                    eflag = 1;
                }
            }
            if (eflag == 0 || (*p != ',' && *p != ')')) 
                goto GCONFin;
            newline(ctx);
        }
        p++;
    }
    err = 0;  

GCONFin:
    if (err)  
        con_free(ctx);
    if (err == -2) {
        p = skip_com(ctx, pp);
        *p = '\0';
        printf1(ctx, "\nError: %s\n",pp);
        if (fflag)  
            printf1(ctx, "Error in parameter index.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_proc()      process constraints for ML estimation.                  */
/*                  write into protocol file.                               */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

int con_proc(TDAContext *ctx, char *pcmd)
{
    register int i,j;
    int err,r,sflag;
    double rnorm,tmp;

    err = -1;
    printf1(ctx, "Checking constraints.\n");
    if (ctx->PMNConS >= ctx->NParm) {
        printf1(ctx, "Error: number of constraints should be less than number of parameters.\n");
        goto CONPROCFin;
    }
    if (get_con(ctx, pcmd))
        goto CONPROCFin;

    printf1(ctx, "\nNumber of constraints: %d\n",ctx->NCon);         
    if (ctx->NCon == 0)
        goto CONPROCFin;

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Constraints (R,D).\n");
        for (i = 1; i <= ctx->NCon; ++i) {
            for (j = 1; j <= ctx->NParm; ++j)
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->ConR[(i - 1) * ctx->NParm + j]);
            rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->ConD[i]);
            fprintf(ctx->PMProtFd,"\n");
        }
        fprintf(ctx->PMProtFd,"\n");
    }

    /* we need additional parameter vectors. ParS for (new) starting values,
       and ParC for expansion. */

    if (!(ctx->ParS = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto CONPROCFin;
    }
    ctx->ParSA = ctx->NParm + 1;
    memrq(ctx, ctx->ParSA,sizeof(double));

    if (!(ctx->ParC = (double *)calloc((size_t)(ctx->NParm + 1),sizeof(double)))) {
        p_err(ctx, -2,1);
        goto CONPROCFin;
    }
    ctx->ParCA = ctx->NParm + 1;
    memrq(ctx, ctx->ParCA,sizeof(double));

    /* check whether starting values fullfil constraints */

    sflag = 0;
    for (i = 1; i <= ctx->NCon; ++i) {
        tmp = 0.0;
        for (j = 1; j <= ctx->NParm; ++j)
            tmp += ctx->ConR[(i - 1) * ctx->NParm + j] * ctx->Par[j];
        if (fabs(tmp - ctx->ConD[i]) > ctx->EPSI1) {
            sflag = 1;
            break;
        }
    }
    if (sflag) {
        printf1(ctx, "Try to find consistent starting values.\n");         

        if (alloc_acu(ctx, ctx->ConRA))
            goto CONPROCFin;

        if (alloc_acv(ctx, ctx->NParm + 1))
            goto CONPROCFin;

        if (alloc_acn(ctx, ctx->NParm + 1))
            goto CONPROCFin;

        for (i = 1; i < ctx->ConRA; ++i)
            ctx->AcU[i] = ctx->ConR[i];

        for (i = 1; i <= ctx->NCon; ++i)
            ctx->AcV[i] = ctx->ConD[i];

        r = lhhfti(ctx, ctx->NCon,ctx->NParm,ctx->NParm,1,ctx->AcU,ctx->AcV,&rnorm,1,ctx->WrkD,ctx->WrkH,ctx->AcN,ctx->EPSI1);

        printf1(ctx, "Rank of constraints: %d\n",r);          
        if (r != ctx->NCon) {
            printf1(ctx, "Error: should be equal to number of constraints.\n");
            goto CONPROCFin;
        }
        if (rnorm > ctx->EPSI1) {
            printf1(ctx, "Can't find a solution of constraints.\n");
            goto CONPROCFin;
        }
        for (i = 1; i <= ctx->NParm; ++i)    /* save new starting values */
            ctx->ParS[i] = ctx->AcV[i];

        if (ctx->PMProtFDef) {
            fprintf(ctx->PMProtFd,"Solution of constraints (new starting values).\n");
            for (i = 1; i <= ctx->NParm; ++i)  
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,ctx->AcV[i]);
            fprintf(ctx->PMProtFd,"\n\n");
        }
        alloc_acu(ctx, 0);
        alloc_acv(ctx, 0);
        alloc_acn(ctx, 0);
    }
    else {
        for (i = 1; i <= ctx->NParm; ++i)    /* use original starting values */
            ctx->ParS[i] = ctx->Par[i];
    }
    ctx->NCon1 = ctx->NCon;               /* used for estimation */
    ctx->NParm1 = ctx->NParm - ctx->NCon;      /* dimensions of reduced parameter space */
    ctx->HSiz1 = ctx->NParm1 * (ctx->NParm1 - 1) / 2;  

    if (con_build(ctx))     
        goto CONPROCFin;

    /* set Par[], used for minimization, to zero */

    for (i = 1; i <= ctx->NParm1; ++i) 
        ctx->Par[i] = 0.0;

    err = 0;

CONPROCFin:
    if (err)
        con_free(ctx);

    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_build.   Build the projection matrix ConQ for function minimization */
/*  ##           in the reduced parameter space. ConQ is a NParm * NParm1   */
/*               matrix. Cf. McCormick 1983, p.265.                         */
/*                                                                          */
/*  Return 0 if OK, -1 if insufficient memory.                              */
/*                  -2 if less than NCon indep constraints.                 */

int con_build(TDAContext *ctx)
{
    register int i,j,k,l;
    int r,n,n1,err;

    err = -1;
    if (!(ctx->ConQ = (double *)calloc((size_t)(ctx->NParm) * (size_t)(ctx->NParm1) + 1,sizeof(double)))) {
        p_err(ctx, -2,1);
        goto BCONFin;
    }
    ctx->ConQA = ctx->NParm * ctx->NParm1 + 1;
    memrq(ctx, ctx->ConQA,sizeof(double));

    if (!(ctx->CIP = (int *)calloc((size_t)(ctx->NParm + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto BCONFin;
    }
    ctx->CIPA = ctx->NParm + 1;
    memrq(ctx, ctx->CIPA,sizeof(int));

    if (alloc_acu(ctx, ctx->ConRA))           /* used as copy of ConR */
        goto BCONFin;

    if (alloc_acw(ctx, ctx->NParm + 1))
        goto BCONFin;

    if (alloc_acn(ctx, ctx->NParm + 1))
        goto BCONFin;

    /*  Try to find NCon linear independent columns of ConR. */

    n = 0;  /* number of columns already found, indices in CIP    */
    k = 0;  /* CIP is used to record the column interchange. */

    for (j = 1; j <= ctx->NParm; ++j)
        ctx->CIP[j] = j;

    for (j = 1; j <= ctx->NParm; ++j) {
        n1 = n + 1;
        for (i = 1; i <= n; ++i) {
            for (l = 1; l <= ctx->NCon; ++l)
                ctx->AcU[(l - 1) * n1 + i] = ctx->ConR[(l - 1) * ctx->NParm + ctx->CIP[i]];
        }
        for (l = 1; l <= ctx->NCon; ++l)
            ctx->AcU[(l - 1) * n1 + n1] = ctx->ConR[(l - 1) * ctx->NParm + j];

        r = lhhfti(ctx, ctx->NCon,n1,n1,0,ctx->AcU,ctx->AcU,ctx->AcW,0,ctx->WrkD,ctx->WrkH,ctx->AcN,ctx->EPSI1);
        if (r == n1) {
            k = ctx->CIP[++n];
            ctx->CIP[n] = j;
            ctx->CIP[j] = k;
        }
        if (n == ctx->NCon)
            break;
    }

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Selected columns for projection matrix.\n ");
        for (i = 1; i <= ctx->NCon; ++i)  
            fprintf(ctx->PMProtFd,"%2d ",ctx->CIP[i]);
        fprintf(ctx->PMProtFd,"\n\n");
    }
    if (n != ctx->NCon) {
        err = -2;
        goto BCONFin;
    }
   
    /*  Now build the first NCon rows of the projection matrix ConQ.    */
    
    for (j = 1; j <= ctx->NCon; ++j) {
        for (i = 1; i <= ctx->NCon; ++i)
            ctx->AcU[(j - 1) * ctx->NCon + i] = ctx->ConR[(j - 1) * ctx->NParm + ctx->CIP[i]];

        for (i = 1; i <= ctx->NParm1; ++i)  
            ctx->ConQ[(j - 1) * ctx->NParm1 + i] =
                             -ctx->ConR[(j - 1) * ctx->NParm + ctx->CIP[ctx->NCon + i]];
    }
    r = lhhfti(ctx, ctx->NCon,ctx->NCon,ctx->NCon,ctx->NParm1,ctx->AcU,ctx->ConQ,ctx->AcW,1,ctx->WrkD,ctx->WrkH,ctx->AcN,ctx->EPSI1);

    if (r != ctx->NCon) {
        n = r;
        err = -2;
        goto BCONFin;
    }
    
    /*  The lower part of ConQ is an identity matrix. */
   
    j = 1;
    for (i = ctx->NCon + 1; i <= ctx->NParm; ++i) {
        ctx->ConQ[(i - 1) * ctx->NParm1 + j] = 1.0;
        j++;
    }
    for (i = 1; i <= ctx->NParm; ++i) {
        for (j = 1; j <= ctx->NParm1; ++j) {
            if (fabs(ctx->ConQ[(i - 1) * ctx->NParm1 + j]) <= ctx->EPSI1)
                ctx->ConQ[(i - 1) * ctx->NParm1 + j] = 0.0;
        }
    }
    
    /*  Build pointers to the interchanged rows of ConQ in CIP. */
   
    for (i = 1; i <= ctx->NParm; ++i)
        ctx->AcN[ctx->CIP[i]] = i;

    for (i = 1; i <= ctx->NParm; ++i)
        ctx->CIP[i] = ctx->AcN[i];

    if (ctx->PMProtFDef) {
        fprintf(ctx->PMProtFd,"Parameter interchange for projection.\n ");
        for (i = 1; i <= ctx->NParm; ++i)
            fprintf(ctx->PMProtFd,"%2d ",ctx->CIP[i]);

        fprintf(ctx->PMProtFd,"\n\nProjection matrix.\n ");
        for (i = 1; i <= ctx->NParm; ++i) {
            for (j = 1; j <= ctx->NParm1; ++j)  
                fprintf(ctx->PMProtFd," %g",ctx->ConQ[(i - 1) * ctx->NParm1 + j]);
            fprintf(ctx->PMProtFd,"\n ");
        }
        fprintf(ctx->PMProtFd,"\n");
    }
    err = 0;

BCONFin:
    if (err == -2)  
        printf1(ctx, "Error: only %d linear independent constraint(s).\n",n);
    alloc_acu(ctx, 0);
    alloc_acw(ctx, 0);
    alloc_acn(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_free()      free memory used for constraints.                       */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

void con_free(TDAContext *ctx)
{
    if (ctx->ParSA) {
        free((char *)ctx->ParS);
        memrq(ctx, -ctx->ParSA,sizeof(double));
        ctx->ParSA = 0;
    }
    if (ctx->ParCA) {
        free((char *)ctx->ParC);
        memrq(ctx, -ctx->ParCA,sizeof(double));
        ctx->ParCA = 0;
    }
    if (ctx->ConRA) {
        free((char *)ctx->ConR);
        memrq(ctx, -ctx->ConRA,sizeof(double));
        ctx->ConRA = 0;
    }
    if (ctx->ConDA) {
        free((char *)ctx->ConD);
        memrq(ctx, -ctx->ConDA,sizeof(double));
        ctx->ConDA = 0;
    }
    if (ctx->ConQA > 0) {
        free((char *)ctx->ConQ);
        memrq(ctx, -ctx->ConQA,sizeof(double));
        ctx->ConQA = 0;
    }
    if (ctx->CIPA > 0) {
        free((char *)ctx->CIP);
        memrq(ctx, -ctx->CIPA,sizeof(int));
        ctx->CIPA = 0;
    }
    ctx->NParm1 = ctx->NParm;
    ctx->HSiz1 = ctx->HSiz;
    ctx->NCon1 = ctx->NCon = 0;
}
