/****************************************************************************/
/*  t_gmin                                                                  */
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
#include "t_tmin.h"  
#include "t_var.h"  
#include "t_gdat.h"  
#include "t_eval.h"  
#include "t_eval1.h"  
#include "t_eval3.h"  
#include "t_lsei.h"  
#include "t_cdf.h"  
#include "t_min.h"  
#include "t_con.h"  
#include "t_edat.h"  
#include "t_mat.h"  
#include "t_matf.h"  
#include "tda_context.h"
#include "tda_compat.h"

/*  functions in t_gmin.c */

int f_min(TDAContext *ctx, int typ);
void gm_coeff(TDAContext *ctx, int n,double *x,int cov,double *diag);
int get_func(TDAContext *ctx, char *s,int opt,int *n3,int iflag,char *msg);
char *get_level(TDAContext *ctx, char *s,int *err); 
void prn_fexp(TDAContext *ctx, char *s,int opt); 
int get_farg(TDAContext *ctx, char *s);         
int fcomp(const void *, const void *, void *);
void prn_feval(TDAContext *ctx);
int get_flval(TDAContext *ctx, double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d);
int get_filval(TDAContext *ctx, int n,double *x,int deriv,int *ni);
int get_fival(TDAContext *ctx, int n,double *x,int deriv,int aflag);
int prn_fsval(TDAContext *ctx, int n,double *x,int gmina,double *lb,double *ub,int opt);
void gm_pres(TDAContext *ctx);          
int get_ifval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu,int deriv, double *gl,double *gu);
int get_ifival(TDAContext *ctx, int n,double *x,int deriv);
int get_icval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu);          

/* ------------------------------------------------------------------------ */
/*  f_min(typ)      User-defined models/functions.                          */
/*                  typ = 0 : ML                                            */
/*                        1 : general functions                             */
/*                        2 : nonlinear regression                          */
/*                        3 : frml (user-def rate models)                   */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int f_min(TDAContext *ctx, int typ)
{
    int err,n,cov,typ1,iflag;

    err = -1;
    iflag = 0;
    if (check_cmd(ctx, 2))
        return(-1);

    n = 4;
    typ1 = typ;

    if (typ == 0) {
        printf1(ctx, "ML estimation of user-defined model. ");
        n--;
    }
    else if (typ == 1)
        printf1(ctx, "Function minimization. ");
    else if (typ == 2)
        printf1(ctx, "Nonlinear regression. ");
    else if (typ == 3) {
        printf1(ctx, "ML estimation of user-defined rate model. ");
        typ1 = 0;
    }
    else
        return(-1);

    printf1(ctx, "Current memory: %d bytes.\n",ctx->MemReq);

    if (typ == 3 && ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(-1);
    }

    set_mldef(ctx);                /* set defaults for ML estimation */
                             
    if (parm(ctx, ctx->CmdBuf + n,12,1))   /* get parameters */
        goto FMLFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 12,4);
           
    newline(ctx);

    if (ctx->FNArgN < 1) {
        printf1(ctx, "Error: function should contain at least one argument.\n");
        goto FMLFin;
    }
    if (typ != 3 && ctx->PMFTYP5 > 0) {
        p_err(ctx, -40,1);
        goto FMLFin;
    }
    if (typ == 3)
        ctx->FVFlg = 2;

    if (ctx->FVFlg) {
        prn_feval(ctx);
        if (typ == 3) {
            ctx->FVFlg = 2;
            printf1(ctx, "Using episode data");
            if (ctx->NSP > 0)  
                printf1(ctx, " (%d episode splits)",ctx->NSplits);
            printf1(ctx, ".\n");
        }
        newline(ctx);
    }
    if (ctx->PMDOPT < 0)             /* default with derivatives */
        ctx->PMDOPT = 2;

    newline(ctx);
    set_mlopt(ctx);                /* adjust ML options */         
    ctx->NParm = ctx->FNArgN;             /* number of parameters */
           
    /* covariance matrix only calculated if function refers to data matrix
       variables and there are degrees of freedom. */

    /**************
    if (CCTyp > 2)
        CCTyp = 2;     
    *************/

    cov = 1;
    if (ctx->FVFlg == 0 || ctx->NParm >= ctx->NOC)
        ctx->CCTyp = cov = 0;

    if (ml_init(ctx, 1,typ1,0))  /* init optimization */
        goto FMLFin;

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto FMLFin;

    /* print parameters and starting values */

    if (prn_fsval(ctx, ctx->NParm,ctx->Par + 1,0,ctx->ParLB + 1,ctx->ParUB + 1,0))
        goto FMLFin;

    if (fnd_alloc(ctx, 1,2,ctx->FNArgN,ctx->FNPN,iflag))   /* add memory for derivatives */
        goto FMLFin;
               
    prot_init(ctx, 3 + typ,UMOD,0);          /* init protocol file */

    if (ctx->PMNConS > 0) {                  /* process constraints */
        newline(ctx);
        if (con_proc(ctx, ctx->CmdBuf))
            goto FMLFin;
    }

    if (ffmin(ctx, UMOD,typ1,1,1))       /* minimization */
        goto FMLFin;                /* insuff memory or fn error */
                   
    prn_mlres(ctx, typ1);                /* print info about min algorithm */

    if (ctx->LConv < 0)
        cov = 0;
   
    gm_coeff(ctx, ctx->NParm,ctx->Par,cov,ctx->Diag);   /* print coefficients */
             
    if (ctx->LConv >= 0) {               /* if convergence */
        newline(ctx);
        prn_ml1res(ctx, typ1);           /* additional ML results */

        if (ctx->PMResFDef && ctx->FVFlg)     /* print residuals */
            gm_pres(ctx);
    }
    err = 0;

FMLFin:
    if (ctx->FVFlg)
        ctx->FVFlg = 1;
           
    con_free(ctx);
    ml_init(ctx, 0,0,0);     
    fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gm_coeff(n,x,cov,diag)                                                  */
/*                                                                          */
/*  Print estimated coefficients x[i], i = 1,...,n.                         */
/*  If cov = 1  diag[] is main diagonal of cov matrix.                      */
/*  If cov = 2  diag[] is full cov matrix.                                  */
   
void gm_coeff(TDAContext *ctx, int n,double *x,int cov,double *diag)
{
    register int i,k,l;
    double tmp,tmp1;

    l = 9;
    if (l < ctx->FNMLen)
        l = ctx->FNMLen;

    printf1(ctx, "Idx  Parameter ");
    prnchar(ctx, ' ',l - 9,0);
    prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Value");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 4,0); printf1(ctx, "Error");
    prnchar(ctx, ' ',ctx->PMTFmt1 - 6,0); printf1(ctx, "Value/E  Signif\n");
    prnchar(ctx, '-',14 + l + 3 * (ctx->PMTFmt1 + 1),1);
    
    for (i = 1; i <= ctx->FNArgN; ++i) {
        k = ctx->FNArgSP[i - 1];
        printf1(ctx, "%3d  %s  ",i,ctx->FNArgDef[k]);
        prnchar(ctx, ' ',l - (int)strlen(ctx->FNArgDef[k]),0);
    
        rt_printf1_d(ctx, ctx->PMTFmtS,x[i]); 

        if (cov) {
            if (cov == 1)
                tmp = diag[i];
            else
                tmp = diag[(i - 1) * n + i];


            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                rt_printf1_d(ctx, ctx->PMTFmtS,tmp);
                if (tmp > 0.0) {
                    tmp1 = x[i] / tmp;
                    tmp = 2.0 * cdnf(ctx, fabs(tmp1)) - 1.0;
                    rt_printf1_d(ctx, ctx->PMTFmtS,tmp1);
                    printf1(ctx, "%7.4lf",tmp);
#ifdef TDA_R_PACKAGE
                    /*  fml/freg/frml have their OWN estimate printer --
                        not prn1_coeff -- so the coeff export never
                        covered them.  Value and Error reached R at full
                        precision only because a ppar= file is written at
                        24.16; Signif is printed at %7.4lf and came back
                        as 1.0000 for anything significant.  */
                    {
                        double erow[4];
                        erow[0] = x[i];
                        erow[1] = tmp;
                        erow[2] = tmp1;
                        erow[3] = 2.0 * cdnf(ctx, fabs(tmp1)) - 1.0;
                        erow[1] = sqrt((cov == 1) ? diag[i]
                                       : diag[(i - 1) * n + i]);
                        tda_export_row(ctx, "fml.est", erow, 4);
                        tda_export_str_row(ctx, "fml.est.names",
                                           ctx->FNArgDef[k]);
                    }
#endif
                }
                else {
                    prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0);
                    printf1(ctx, "---     ---"); 
                }
            }
            else {
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "--- "); 
                prnchar(ctx, ' ',ctx->PMTFmt1 - 3,0); printf1(ctx, "---     ---"); 
            }
        }
        newline(ctx);           
    }
}

/* ------------------------------------------------------------------------ */
/*  get_func(s,opt,n3,iflag,msg)                                            */
/*                                                                          */
/*                      If opt != 0 get definition of function in s.        */
/*                      If opt == 0 free previously allocated memory.       */
/*                      If iflag != 0 allocate also FNArgVal1[].            */
/*                                                                          */
/*                      FNFlg = 1 if valid function definition.             */
/*                      FCFlg = 1 if valid constraint function definition.  */
/*                      FVFlg = 1 if references to data matrix variables.   */
/*                      FNMLen = max length of parameter strings.           */
/*                                                                          */
/*                      FNArgN = number of arguments. Sorted in FNArgDef,   */
/*                      meaning that FNArgSP[k] is the k.th argument.       */
/*                                                                          */
/*                      n3 is set if function contains type 3 operators     */
/*                      or refers to type 5 variables.                      */
/*                                                                          */
/*                      Introduced optional levels, see get_level().        */  
/*                                                                          */
/*                      Return 0 if OK, -1 if error.                        */

int get_func(TDAContext *ctx, char *s,int opt,int *n3,int iflag,char *msg)
{
    register int i;
    int err,np,fin,l,ll,n,nl,nv,nc,n2,fflag,cflag;
    register char c,*p,*q;

    ctx->FNNLev = fflag = err = np = 0;
    if (opt == 0)
        goto GFFin;

    printf1(ctx, "\nFunction definition");
    if (msg != NULL)
        printf1(ctx, " (%s):\n",msg);
    else
        printf1(ctx, ":\n");

    if (!strncmp(s,"level",5)) {    /* get level */
        s =  get_level(ctx, s + 5,&err);
        if (err)
            goto GFFin;
    }
    err = -1;
    p = s;
      
    cflag = fin = 0;
    ctx->FNMLen = ctx->FVFlg = ctx->FNFlg = ctx->FCFlg = 0;
    *n3 = 0;

    while (*p) {

        s = p;
        if (!strncmp(p,"fn=",3)) {
            fin = 1;
            p += 2;
            *p = '\0';
        }
        else if (!strncmp(p,"con=",4)) {
            fin = -1;
            p += 3;
            *p = '\0';
        }
        else {
            l = 0;
            while (*p && ((*p >= 'a' && *p <= 'z') || isdigit((int)*p))) {
                l++;
                p++;
            }
            if (l == 0 || *p != '=' || !*(p + 1)) {
                q = skip_expr(ctx, s);
                if (!*q) {
                    fin = 2;
                    p = s;
                    goto GFCONT;
                }
                goto GFFin;
            }
            *p = '\0';

            if (v_search(ctx, s,&nl,&n) >= 0) {
                err = -3;
                goto GFFin;
            }
            if (np >= MaxFNP) {
                err = -4;
                goto GFFin;
            }   
            if (!(ctx->FNPDef[np] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                err = -2;
                goto GFFin;   
            }
            ctx->FNPLen[np] = (short)(l);
            memrq(ctx, l + 1,sizeof(char));
            strcpy(ctx->FNPDef[np],s);
            ctx->FNPN++;
        }
        *p++ = '=';
GFCONT:
        q = skip_expr(ctx, p);
        c = *q;
        *q = '\0';

        prn_fexp(ctx, s,fin);       /* print expression */

        if (fin >= 0) {
            if (get_farg(ctx, p)) { /* check expression for undefined arguments */
                err = -5;
                goto GFFin;
            }
        }

        if ((n = v_parse(ctx, p,iflag)) < 0 || ctx->ESCnt <= 0) {
            printf1(ctx, "Syntax error (%d) in expression: %s.\n",n,s);
            if (n < 0)
                prn_emsg1(ctx, n);
            err = -5;                       
            goto GFFin;
        }
        if (fin > 0) {
            if (!(ctx->FNETyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(ctx->FNEVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
                free((char *)ctx->FNETyp);
                err = -2;
                goto GFFin;   
            }
            ctx->FNECnt = ctx->ESCnt;
            memrq(ctx, ctx->FNECnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ctx->ESCnt; ++i) {
                ctx->FNETyp[i] = ctx->ESTyp[i];
                ctx->FNEVal[i] = ctx->ESVal[i];
            }

            /* check for unused arguments */

            /*******************************
            for (i = 0; i < FNArgN; ++i) {
                fflag = 0;
                for (j = 0; j < ESCnt; ++j) {
                    if (FNETyp[j] == FAOFFS + i) {
                        fflag = 1;
                        break;
                    }
                }
                if (fflag == 0) {
                    printf1(ctx, "Error: at least one argument does not occur in function.\n");
                    err = -5;
                    goto GFFin;
                }
            }
            *****************/
            fflag = 1;
        }
        else if (fin < 0) {
            if (fflag == 0) {
                printf1(ctx, "\nError: using constraints requires previous function definition.\n");
                goto GFFin;
            }
            if (!(ctx->FNECTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(ctx->FNECVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
                free((char *)ctx->FNECTyp);
                err = -2;
                goto GFFin;   
            }
            ctx->FNECCnt = ctx->ESCnt;
            memrq(ctx, ctx->FNECCnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ctx->ESCnt; ++i) {
                ctx->FNECTyp[i] = ctx->ESTyp[i];
                ctx->FNECVal[i] = ctx->ESVal[i];
            }
            cflag = 1;
        }
        else {
            if (!(ctx->FNPETyp[np] = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(ctx->FNPEVal[np] = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
                free((char *)ctx->FNPETyp[np]);
                err = -2;
                goto GFFin;   
            }
            ctx->FNPECnt[np] = ctx->ESCnt;
            memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ctx->ESCnt; ++i) {
                ctx->FNPETyp[np][i] = ctx->ESTyp[i];
                ctx->FNPEVal[np][i] = ctx->ESVal[i];
            }
            np++;
        }
        *q = c;
        if (*q != ',') {
            p = q;
            break;
        }
        p = q + 1;

    }
    if (*p) {        
        printf1(ctx, "Error: definition must end with fn or con expression.\n");
        err = -5;
        goto GFFin;
    }
    if (fin == 0) {
        printf1(ctx, "Error: need an fn=... expression.\n");
        err = -5;
        goto GFFin;
    }

    /* allocate memory for parameter values */

    if (ctx->FNArgN > 0) {
        if (!(ctx->FNArgVal = (double *)calloc((size_t)(ctx->FNArgN),sizeof(double)))) {
            err = -2;
            goto GFFin;   
        }
        ctx->FNArgValA = ctx->FNArgN;
        memrq(ctx, ctx->FNArgN,sizeof(double));

        if (iflag) {
            if (!(ctx->FNArgVal1 = (double *)calloc((size_t)(ctx->FNArgN),sizeof(double)))) {
                err = -2;
                goto GFFin;   
            }
            ctx->FNArgVal1A = ctx->FNArgN;
            memrq(ctx, ctx->FNArgN,sizeof(double));
        }
    }
   
    if (ctx->FNArgN > 0) {       /* sorting */

        if (!(ctx->FNArgSP = (short *)calloc((size_t)(ctx->FNArgN),sizeof(short)))) { 
            err = -2;     
            goto GFFin;    
        }
        memrq(ctx, ctx->FNArgN,sizeof(short));
        ctx->FNArgSPA = ctx->FNArgN;

        if (!(ctx->FNArgSPI = (short *)calloc((size_t)(ctx->FNArgN),sizeof(short)))) { 
            err = -2;     
            goto GFFin;    
        }
        memrq(ctx, ctx->FNArgN,sizeof(short));
        ctx->FNArgSPIA = ctx->FNArgN;

        for (i = 0; i < ctx->FNArgN; ++i)  
            ctx->FNArgSP[i] = (short)(i);

        tda_qsort_r((char *)ctx->FNArgSP,(size_t)(ctx->FNArgN),sizeof(short), fcomp, ctx);

        for (i = 0; i < ctx->FNArgN; ++i)  
            ctx->FNArgSPI[ctx->FNArgSP[i]] = (short)(i);
    }

    /* check for data matrix variables */

    n2 = nc = nv = 0;
    for (i = 0; i < ctx->FNPN; ++i) {

        check_expr(ctx, ctx->FNPECnt[i],ctx->FNPETyp[i],&n,&nl,&l,&ll,0);
        nv += n;
        nc += nl;
        n2 += l;
        *n3 += ll;
    }
    check_expr(ctx, ctx->FNECnt,ctx->FNETyp,&n,&nl,&l,&ll,0);
    nv += n;
    nc += nl;
    n2 += l;
    *n3 += ll;
    if (nc > 0) {
        printf1(ctx, "\nError: function may not contain ci references.\n");
        err = -5;  
        goto GFFin;
    }
    if (n2 > 0) {
        printf1(ctx, "\nError: function may not contain type 2 operators or type 4 variables.\n");
        err = -5;  
        goto GFFin;
    }
    if (nv > 0)
        ctx->FVFlg = 1;
    else if (ctx->FNNLev) {
        printf1(ctx, "\nError: multilevel function requires reference to data matrix.\n");
        err = -6; 
        goto GFFin;
    }
    ctx->FNFlg = 1;

    if (cflag) {        /* check constraint function */

        check_expr(ctx, ctx->FNECCnt,ctx->FNECTyp,&n,&nl,&l,&ll,0);

        if (n > 0 || nl > 0 || l > 0 || ll > 0) {
            printf1(ctx, "\nError: constraint function may not contain ");
            if (n > 0)  
                printf1(ctx, "data matrix variables.\n");
            else if (nl > 0)  
                printf1(ctx, "ci references.\n");
            else if (l > 0)  
                printf1(ctx, "type 2 operators.\n");
            else if (ll > 0)  
                printf1(ctx, "type 3 operators.\n");
            err = -5;  
            goto GFFin;
        }
        ctx->FCFlg = 1;
    }
    if (ctx->FNNLev > 0) {       /* check for correct recognition of levels */
        n = ctx->FNNLev;
        for (i = 0; i < ctx->FNPN; ++i) {
            if (sscanf(ctx->FNPDef[i],"fn%d",&nl) == 1 && nl == n) {
                ctx->FNPELev[i] = (short)(n);
                n--;
            }
        }
        if (n != 0) {
            printf1(ctx, "\nError: function definition inconsistent with levels.\n");
            printf1(ctx, "Need: fn%d",ctx->FNNLev);
            for (i = ctx->FNNLev - 1; i > 0; --i)
                printf1(ctx, ", fn%d",i);
            printf1(ctx, ", and fn expression.\n");   
            err = -6;
            goto GFFin;
        }
    }
    err = 0;

GFFin:
    if (err == -1)
        printf1(ctx, "Syntax error: %s\n",s);
    else if (err == -2)
        p_err(ctx, -2,1);
    else if (err == -3) {
        printf1(ctx, "Error: %s\n",s);
        printf1(ctx, "Parameter is reserved, already used, or not allowed.\n");
    }
    else if (err == -4)  
        printf1(ctx, "Error: exceeded max number of intermediate parameters.\n");

    if (err)
        err = -1;

    if (opt == 0 || err) {

        for (i = 0; i < ctx->FNPN; ++i) {
            if (ctx->FNPLen[i] > 0) {
                free((char *)ctx->FNPDef[i]);
                memrq(ctx, -ctx->FNPLen[i] - 1,sizeof(char));
                ctx->FNPLen[i] = 0;
            }
            if (ctx->FNPECnt[i] > 0) {
                free((char *)ctx->FNPETyp[i]);
                free((char *)ctx->FNPEVal[i]);
                memrq(ctx, -ctx->FNPECnt[i],sizeof(int) + sizeof(double));
                ctx->FNPECnt[i] = 0;
            }
        }
        for (i = 0; i < ctx->FNArgN; ++i) {
            if (ctx->FNArgLen[i] > 0) {
                free((char *)ctx->FNArgDef[i]);
                memrq(ctx, -ctx->FNArgLen[i] - 1,sizeof(char));
                ctx->FNArgLen[i] = 0;
            }
        }
        if (ctx->FNECnt > 0) {
            free((char *)ctx->FNETyp);
            free((char *)ctx->FNEVal);
            memrq(ctx, -ctx->FNECnt,sizeof(int) + sizeof(double));
            ctx->FNECnt = 0;
        }
        if (ctx->FNECCnt > 0) {
            free((char *)ctx->FNECTyp);
            free((char *)ctx->FNECVal);
            memrq(ctx, -ctx->FNECCnt,sizeof(int) + sizeof(double));
            ctx->FNECCnt = 0;
        }
        if (ctx->FNArgSPA > 0) {
            free((char *)ctx->FNArgSP);
            memrq(ctx, -ctx->FNArgSPA,sizeof(short));
            ctx->FNArgSPA = 0;         
        }
        if (ctx->FNArgSPIA > 0) {
            free((char *)ctx->FNArgSPI);
            memrq(ctx, -ctx->FNArgSPIA,sizeof(short));
            ctx->FNArgSPIA = 0;         
        }
        if (ctx->FNArgValA > 0) {
            free((char *)ctx->FNArgVal);
            memrq(ctx, -ctx->FNArgValA,sizeof(double));
            ctx->FNArgValA = 0;         
        }
        if (ctx->FNArgVal1A > 0) {
            free((char *)ctx->FNArgVal1);
            memrq(ctx, -ctx->FNArgVal1A,sizeof(double));
            ctx->FNArgVal1A = 0;         
        }
        if (ctx->FNLVarA > 0) {
            free((char *)ctx->FNLVar);
            memrq(ctx, -ctx->FNLVarA,sizeof(int));
            ctx->FNLVarA = 0;         
        }
        ctx->FNNLev = ctx->FNPN = ctx->FNArgN = ctx->FNFlg = ctx->FCFlg = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_level(s,err)    scan level expression.                              */
/*                                                                          */
/*                      level(1+VName,2+VName,...)                          */
/*                      there might be * instead of +                       */
/*                                                                          */
/*                      FNNLev = number of level                            */
/*                      FNLVar = index of variable                          */
/*                                                                          */
/*  Return err = 0 if OK, -1 if syntax error, -2 if insufficient memory,    */
/*         err = -9 if undefined variable.                                  */
/*                                                                          */
/*  Return pointer to next character in string s.                           */
   
char *get_level(TDAContext *ctx, char *s,int *err)  
{
    register int i;
    int n;
    register char *p,*q;
    char vname[VNLMax + 1];

    *err = -1;
    p = s; 
    if (*p++ != '(')
        return(s);
    q = p;
    ctx->FNNLev = 1;
    while (*q && *q != ')') {
        if (*q++ == ',')
            ctx->FNNLev++; 
    }
    if (*q++ != ')' || *q++ != ',')
        return(s);

    if (!(ctx->FNLVar = (int *)calloc((size_t)(ctx->FNNLev),sizeof(int)))) {
        *err = -2;
        return(s);
    }
    memrq(ctx, ctx->FNNLev,sizeof(int));
    ctx->FNLVarA = ctx->FNNLev;

    for (i = 1; i <= ctx->FNNLev; ++i) {

        if (sscanf(p,"%d",&n) != 1 || n != i)
            return(s);

        p = skip_int(ctx, p);

        if (*p++ != '=')
            return(s);

        n = get_vidx1(ctx, p,vname);
        if (n < 0) {
            printf1(ctx, "Undefined variable name in level expression.\n");
            *err = -9;
            return(s);
        }
        ctx->FNLVar[i - 1] = n;
        p += strlen(vname);
        if (*p != ',' && *p != ')')
            return(s);
        p++;
    }
    newline(ctx);
    for (i = 0; i < ctx->FNNLev; ++i)  
        printf1(ctx, "Level%3d defined by: %s\n",i + 1,ctx->VName[ctx->FNLVar[i]]);
    newline(ctx);
    *err = 0;
    return(q);
}

/* ------------------------------------------------------------------------ */
/*  prn_fexp(s,opt)     Print function expression.                          */

void prn_fexp(TDAContext *ctx, char *s,int opt)  
{
    register int l;
    register char *p;

    if (opt == 2)
        printf1(ctx, "fn      = ");
        
    l = 8;
    p = s;
    while (*p) {
        if (*p == '=') {
            while (l-- > 0)
                printf1(ctx, " ");
            printf1(ctx, "= ");
        }
        else
            printf1(ctx, "%c",*p);
        if (l > 0)
            l--;
        p++;
    }
    printf1(ctx, "\n");
}

/* ------------------------------------------------------------------------ */
/*  get_farg(s)     Check expression s for new arguments. Save arguments    */
/*                  in FNArgLen[] and FNArgDef[]. Count FNArgN.             */
/*                  Max number of arguments is MaxP.                        */
/*                  Return 0 if OK, -1 if error (exceeded MaxP, or if       */
/*                  insufficient memory).                                   */

int get_farg(TDAContext *ctx, char *s)          
{
    register int l;
    register char c,*p;
    int n,nl,len; 

    while (*s) {

        if (check_vname(ctx, s)) {       /* check for variable */
            if (v_search1(ctx, s,&len) >= 0)
                s += len;
        }         
        if (!*s)
            break;
        p = s;
        l = 0;

        if (*p >= 'a' && *p <= 'z') {

            while (*p && ((*p >= 'a' && *p <= 'z') || isdigit((int)*p))) {
                l++;
                p++;
            }
            if (*p == '(')  
                l = 0;

            if (l > 0) {
                c = *p;
                *p = '\0';
                if (v_search(ctx, s,&nl,&n) < 0) {
                    if (ctx->FNArgN >= MaxP) {
                        printf1(ctx, "Error: exceeded max number of arguments.\n");
                        return(-1);
                    }
                    if (!(ctx->FNArgDef[ctx->FNArgN] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                        p_err(ctx, -2,1);
                        return(-1);
                    }
                    ctx->FNArgLen[ctx->FNArgN] = (short)(l);
                    memrq(ctx, l + 1,sizeof(char));
                    strcpy(ctx->FNArgDef[ctx->FNArgN],s);
                    ctx->FNArgN++;
                    if (ctx->FNMLen < l)
                        ctx->FNMLen = l;
                }
                *p = c;
            }
            else
                p++;
        }
        else
            p++;
        s = p;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fcomp()    compare function                                             */

int fcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;     
    char *p1,*p2;

    p1 = ctx->FNArgDef[*(short *)arg1];
    p2 = ctx->FNArgDef[*(short *)arg2];
    return(strcmp(p1,p2));

}

/* ------------------------------------------------------------------------ */
/*  prn_feval()     print info about function evaluation.                   */

void prn_feval(TDAContext *ctx)
{
    register int i,ix,n;
    double tmp,tmp1;

    if (ctx->FVFlg) {
        printf1(ctx, "Function evaluation: ");
        if (ctx->FNNLev == 0)
            printf1(ctx, "sum over %d data matrix cases.\n",ctx->NOC);
        else {
            printf1(ctx, "according to %d level(s) in data matrix.\n",ctx->FNNLev);
            n = 0;
            ix = ctx->FNLVar[0];
            tmp1 = get_data(ctx, ix,0) - 1.0;
            for (i = 0; i < ctx->NOC; ++i) {
                tmp = get_data(ctx, ix,i);
                if (tmp != tmp1) {
                    n++;
                    tmp1 = tmp;
                }
            }
            printf1(ctx, "Number of level 1 units: %d\n",n);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_flval(f,n,x,deriv,opt,g,h,d)                                        */
/*  ##                                                                      */
/*                      get value from function expression.                 */  
/*                      if FVFlg == 1 sum over all data matrix cases.       */
/*                      if FVFlg == 2 use episode data.                     */
/*                                                                          */
/*                      Function arguments: x[i], i = 0,...,n-1.            */
/*                      Return value in f.                                  */
/*                                                                          */
/*  If deriv >= 1  calculate gradient.      FNGrad[0][...]                  */
/*  If deriv == 2  calculate hessian.       FNHess[0][...]                  */
/*  If CCovFlg ==1 use outer product of gradients instead of hessian        */
/*                 (only if reference to data matrix)                       */
/*  If CGradFlg == 1 calculate gradients in MatVal[MPGradIdx]               */
/*                                                                          */
/*  opt = 0 : do not use g,h,d                                              */
/*  opt = 1 : save gradient in g[1,...,n], full hessian in h[1,...],        */
/*            and do not use d.                                             */
/*  opt = 2 : save gradient in g[1,...,n], lower triangle of hessian in     */
/*            h[1,...] and diagonal of hessian in d[1,...,n].               */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_flval(TDAContext *ctx, double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d)
{
    register int i,j,k,l;
    int err,nn,ii,nii,deriv1,sn,org,des,spl,nspl;
    double tmp,ts,tf;

    ctx->NOCUsed = err = 0;
    for (i = 0; i < n; ++i)         /* set values for evaluation */
        ctx->FNArgVal[i] = x[i];

    if (ctx->FVFlg == 0) {               /* no reference to data matrix, evaluate
                                       expression */
        err = get_fival(ctx, 0,f,deriv,0);

        if (err == 0 && deriv && opt) {

            for (i = 1; i <= n; ++i) {              /* save gradient */
                g[i] = ctx->FNGrad[0][i - 1];

                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                    ctx->MatVal[ctx->MPGradIdx][i] = g[i];        
            }

            if (deriv > 1) {                        /* save hessian */
                if (opt == 1) {
                    k = 0;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = ctx->FNHess[0][k++];
                            h[(i - 1) * n + j] = tmp;
                            if (j < i)
                                h[(j - 1) * n + i] = tmp;
                        }
                    }
                }
                else if (opt == 2) {

                    k = 0;
                    l = 1;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = ctx->FNHess[0][k++];
                            if (j < i)
                                h[l++] = tmp;
                            else
                                d[i] = tmp;
                        }
                    }
                }
            }
        }
        ctx->NOCUsed = 1;
        return(err);
    }

    /*  reference to data matrix. Sum over all level 1 units */

    *f = 0.0;
    if (deriv && opt) {
        for (i = 1; i <= n; ++i)
            g[i] = 0.0;
        if (deriv > 1) {
            if (opt == 1) {
                nn = n * n;
                for (i = 1; i <= nn; ++i)
                    h[i] = 0.0;
            }
            else if (opt == 2) {
                for (i = 1; i <= n; ++i)
                    d[i] = 0.0;
                nn = n * (n - 1) / 2;
                for (i = 1; i <= nn; ++i)
                    h[i] = 0.0;
            }   
        }
    }
    deriv1 = deriv;
    if (deriv1 == 2 && ctx->CCovFlg != 0)
        deriv1 = 1;

    if (ctx->FVFlg == 2)
        get_spell(ctx, 1,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    else
        ii = 0;

    while (1) {
    
        if (ctx->FVFlg == 2) {
            if (get_spell(ctx, 0,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl) == 0) 
                break;
        }
        err = get_filval(ctx, ii,&tmp,deriv1,&nii);

        if (err)
            break;

        *f += tmp;

        if (deriv && opt) {

            for (i = 1; i <= n; ++i) {          /* save gradient */
                g[i] += ctx->FNGrad[0][i - 1];

                if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow)
                    ctx->MatVal[ctx->MPGradIdx][ctx->NOCUsed * ctx->MPGradCol + i] = ctx->FNGrad[0][i - 1];
            }
            if (deriv > 1) {                    /* save hessian */
                if (opt == 1) {
                    k = 0;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            if (ctx->CCovFlg == 0)
                                tmp = ctx->FNHess[0][k++];
                            else
                                tmp = ctx->FNGrad[0][i - 1] * ctx->FNGrad[0][j - 1];

                            h[(i - 1) * n + j] += tmp;
                            if (j < i)
                                h[(j - 1) * n + i] += tmp;
                        }
                    }
                }
                else if (opt == 2) {

                    k = 0;
                    l = 1;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {

                            if (ctx->CCovFlg == 0)
                                tmp = ctx->FNHess[0][k++];
                            else
                                tmp = ctx->FNGrad[0][i - 1] * ctx->FNGrad[0][j - 1];
                            if (j < i)
                                h[l++] += tmp;
                            else
                                d[i] += tmp;
                        }
                    }
                }
            }
        }
        if (ctx->CGradFlg && ctx->NOCUsed < ctx->MPGradRow && ctx->PM2NV > 0)
            mp_putvar(ctx, ctx->NOCUsed,ctx->MPGradCol,ctx->MatVal[ctx->MPGradIdx],ctx->PM2NV,ctx->PM2VIdx,ii);

        ctx->NOCUsed++;
        if (ctx->FVFlg != 2) {
            ii = nii;
            if (ii >= ctx->NOC)
                break;
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_filval(n,x,deriv,ni)    get value from multilevel function for      */  
/*                              level 1 unit beginning with data            */
/*  ##                          matrix case n.                              */
/*                              Return value in x, return number of next    */
/*                              data matrix case in ni.                     */
/*                                                                          */  
/*  If deriv >= 1  Gradient is in FNGrad[0][...]                            */
/*  If deriv == 2  Hessian is  in FNHess[0][...]                            */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_filval(TDAContext *ctx, int n,double *x,int deriv,int *ni)
{
    register int i,j,l;
    int err,last;
    int first[100];
    double ltmp,lcase[100];

    err = 0;
    for (l = 0; l < ctx->FNNLev; ++l) {
        first[l] = 1;
        lcase[l] = get_data(ctx, ctx->FNLVar[l],n);
    }
    l = ctx->FNNLev - 1;

FICONT1:
    for (i = 0; i < ctx->FNPN; ++i) {

        err = v_eval1(ctx, n,ctx->FNPECnt[i],ctx->FNPETyp[i],ctx->FNPEVal[i],ctx->ESIdx,x,0,deriv,0,0,0);
        if (err)  
            return(err);

        if (l < 0 || ctx->FNPELev[i] != l + 1) {  
            ctx->FNPVal[i] = *x;
            if (deriv > 0) {        /* copy gradient and hessian */

                for (j = 0; j < ctx->FNGradL; ++j)
                    ctx->FNPGrad[i][j] = ctx->FNGrad[0][j];

                if (deriv > 1) {    
                    for (j = 0; j < ctx->FNHessL; ++j)
                        ctx->FNPHess[i][j] = ctx->FNHess[0][j];
                }
                ctx->FNPEATyp[i] = (short)(ctx->FNEVATyp);
            }
        }
        else {
            if (first[l]) {
                ctx->FNPVal[i] = *x;
                if (deriv > 0) {        /* copy gradient and hessian */

                    for (j = 0; j < ctx->FNGradL; ++j)
                        ctx->FNPGrad[i][j] = ctx->FNGrad[0][j];

                    if (deriv > 1) {    
                        for (j = 0; j < ctx->FNHessL; ++j)
                            ctx->FNPHess[i][j] = ctx->FNHess[0][j];
                    }
                    ctx->FNPEATyp[i] = (short)(ctx->FNEVATyp);
                }
                first[l] = 0;
            }
            else {   
                ctx->FNPVal[i] += *x;
                if (deriv > 0) {        /* copy gradient and hessian */

                    for (j = 0; j < ctx->FNGradL; ++j)
                        ctx->FNPGrad[i][j] += ctx->FNGrad[0][j];

                    if (deriv > 1) {    
                        for (j = 0; j < ctx->FNHessL; ++j)
                            ctx->FNPHess[i][j] += ctx->FNHess[0][j];
                    }
                    ctx->FNPEATyp[i] = (short)(ctx->FNEVATyp);
                }
            }
            last = 0;
            if (n + 1 >= ctx->NOC)
                last = 1;   
            else {          
                for (j = l; j >= 0; --j) {
                    ltmp = get_data(ctx, ctx->FNLVar[j],n + 1);
                    if (ltmp != lcase[j]) {
                        first[l] = 1;
                        last = 1;
                        break;
                    }
                }
            }
            if (last == 0) {
                n++;
                for (j = 0; j < ctx->FNNLev; ++j) {
                    lcase[j] = get_data(ctx, ctx->FNLVar[j],n);
                    if (j > l)
                        first[j] = 1;
                }
                l = ctx->FNNLev - 1;
                goto FICONT1;
            }
            l--;
        }
    }

    *ni = n + 1;

    /* evaluate final function expression */

    err = v_eval1(ctx, n,ctx->FNECnt,ctx->FNETyp,ctx->FNEVal,ctx->ESIdx,x,0,deriv,0,0,0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_fival(n,x,deriv,aflag)                                              */
/*                                                                          */
/*  Get value from function expression.                                     */  
/*  if FVFlg != 0 use data matrix case n.                                   */
/*  if aflag != 0 use AVVAL in call of v_eval1().                           */
/*  Return value in x.                                                      */
/*                                                                          */  
/*  If deriv >= 1  Gradient is in FNGrad[0][...]                            */
/*  If deriv == 2  Hessian is  in FNHess[0][...]                            */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_fival(TDAContext *ctx, int n,double *x,int deriv,int aflag)
{
    register int i,j;
    int err = 0;

    /* first evaluate intermediate parameters */

    for (i = 0; i < ctx->FNPN; ++i) {
        if (ctx->FNPECnt[i] > 0) {

            err = v_eval1(ctx, n,ctx->FNPECnt[i],ctx->FNPETyp[i],ctx->FNPEVal[i],ctx->ESIdx,x,aflag,deriv,0,0,0);
            if (err)  
                return(err);
            ctx->FNPVal[i] = *x;

            if (deriv > 0) {        /* copy gradient and hessian */

                for (j = 0; j < ctx->FNGradL; ++j)
                    ctx->FNPGrad[i][j] = ctx->FNGrad[0][j];

                if (deriv > 1) {    
                    for (j = 0; j < ctx->FNHessL; ++j)
                        ctx->FNPHess[i][j] = ctx->FNHess[0][j];
                }
                ctx->FNPEATyp[i] = (short)(ctx->FNEVATyp);
            }
        }
    }

    /* then evaluate function */

    err = v_eval1(ctx, n,ctx->FNECnt,ctx->FNETyp,ctx->FNEVal,ctx->ESIdx,x,aflag,deriv,0,0,0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_fsval(n,x,gmina,lb,ub,opt)                                          */
/*                                                                          */
/*  Print parameters and starting values: x[i], i = 0,...,n-1               */                    
/*  If gmina != 0 print also bounding box. If opt != 0 do not print         */
/*  parameters.                                                             */  
/*                                                                          */
/*  Return 0 if OK, -1 if error in bounds.                                  */

int prn_fsval(TDAContext *ctx, int n,double *x,int gmina,double *lb,double *ub,int opt)
{
    register int i,j;
    int l,err;

    err = 0;
    l = 10;
    if (l < ctx->FNMLen)
        l = ctx->FNMLen;

    printf1(ctx, "Idx  Parameter ");
    if (opt == 0) {
        prnchar(ctx, ' ',l - 9,0);
        printf1(ctx, " Starting value");
    }
    if (gmina)
        printf1(ctx, "      Lower bound     Upper bound");
    newline(ctx);
    for (i = 0; i < n; ++i) {
        j = ctx->FNArgSP[i];
        printf1(ctx, "%3d  %s  ",i + 1,ctx->FNArgDef[j]);
        prnchar(ctx, ' ',l - (int)strlen(ctx->FNArgDef[j]),0);
        if (opt == 0)
            printf1(ctx, "%15.8e ",x[i]);
        if (gmina) {
            printf1(ctx, "%15.8e ",lb[i]);
            printf1(ctx, "%15.8e ",ub[i]);
            if (opt == 0) {
                if (lb[i] > x[i] + ctx->EPSI || ub[i] < x[i] - ctx->EPSI)
                    err = -1;
            }
            else {
                if (lb[i] >= ub[i] + ctx->EPSI)
                    err = -1;
            }
        }
        printf1(ctx, "\n");
    }
    newline(ctx);
    if (err)  
        printf1(ctx, "Error in bounding box.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gm_pres()   print residuals to PMResFd.                                 */

void gm_pres(TDAContext *ctx)           
{
    register int i,j,k;
    int err;
    double tmp;

    if (ctx->PMResFDef == 0)
        return;

    for (i = 0; i < ctx->FNArgN; ++i)         /* set values for evaluation */
        ctx->FNArgVal[i] = ctx->Par[i + 1];
  
    for (i = 0; i < ctx->NOC; ++i) {
        fprintf(ctx->PMResFd,"%6d ",i + 1);
        err = get_fival(ctx, i,&tmp,0,0);
        if (err) {
            fprintf(ctx->PMResFd,"error\n");
            continue;
        }
        rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMFmtS,tmp);
        for (j = 0; j < ctx->PMNV; ++j) {   /* add variables */
            k = ctx->PMVIdx[j];
            tmp = get_data(ctx, k,i);
            rt_fprintf_d(ctx, ctx->PMResFd,ctx->VPFmtS[k],tmp);
        }
        fprintf(ctx->PMResFd,"\n");
    }
    printf1(ctx, "Residuals written to: %s\n",ctx->PMResFName);
}

/* ------------------------------------------------------------------------ */
/*  get_ifval(fl,fu,n,xl,xu,deriv,gl,gu)                                    */
/*  ##                                                                      */
/*  Calculates an inclusion function. Function arguments are                */
/*  [xl[i],xu[i]] for i = 0,...,n - 1. Return interval in [fl,fu].          */
/*                                                                          */
/*  If deriv = 1 calculate also gradiend bounds in gl[i],gu[i], i=0,n-1     */
/*                                                                          */
/*  if FVFlg == 1 sum over all data matrix cases.                           */
/*  if FVFlg == 2 use episode data.                                         */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_ifval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu,int deriv, double *gl,double *gu)
{
    register int i;
    int err,ii,sn,org,des,spl,nspl;
    double ts,tf,xx[2];

    err = 0;
    for (i = 0; i < n; ++i) {       /* set values for evaluation */
        ctx->FNArgVal[i] = xl[i];
        ctx->FNArgVal1[i] = xu[i];
    }
    *fl = *fu = 0.0;

    if (ctx->FVFlg == 0) {

        err = get_ifival(ctx, 0,xx,deriv);
        if (err == 0) {
            *fl = xx[0];
            *fu = xx[1];

            if (deriv) {
                for (i = 0; i < n; ++i) {              /* save gradient */
                    gl[i] = ctx->FNGrad[0][i];
                    gu[i] = ctx->FNGrad[1][i];
                }
            }
        }
    }
    else {                      /* sum over data matrix cases */
                                /* or use episode data */
        if (deriv) {
            for (i = 0; i < n; ++i)            
                gu[i] = gl[i] = 0.0;           
        }
        if (ctx->FVFlg == 2)
            get_spell(ctx, 1,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl);
        else
            ii = 0;

        while (1) {
    
            if (ctx->FVFlg == 2) {
                if (get_spell(ctx, 0,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl) == 0) 
                    break;
            }
            err = get_ifival(ctx, ii,xx,deriv);
            if (err)
                break;

            *fl += xx[0];
            *fu += xx[1];

            if (deriv) {
                for (i = 0; i < n; ++i) {              /* save gradient */
                    gl[i] += ctx->FNGrad[0][i];
                    gu[i] += ctx->FNGrad[1][i];
                }
            }
    
            if (ctx->FVFlg == 1) {
                if (++ii >= ctx->NOC)
                    break;
            }
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_ifival(n,x,deriv)                                                   */
/*                                                                          */
/*  Get value from function expression, using intervals. If FVFlg != 0      */
/*  use data matrix case n. If deriv = 1 calculate also gradient range.     */
/*                                                                          */
/*  Return lower bound: x[0]                                                */
/*         upper bound: x[1]                                                */
/*                                                                          */  
/*  If deriv = 1 then Gradient lower bound in FNGrad[0][...]                */
/*                             upper bound in FNGrad[1][...]                */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_ifival(TDAContext *ctx, int n,double *x,int deriv)
{
    register int i,j;
    int err = 0;

    /* first evaluate intermediate parameters */

    for (i = 0; i < ctx->FNPN; ++i) {
        if (ctx->FNPECnt[i] > 0) {
            err = v_eval1(ctx, n,ctx->FNPECnt[i],ctx->FNPETyp[i],ctx->FNPEVal[i],ctx->ESIdx,x,0,deriv,0,0,1);
            if (err)  
                return(err);
            ctx->FNPVal[i] = x[0];
            ctx->FNPVal1[i] = x[1];

            if (deriv) {        /* copy gradient */

                for (j = 0; j < ctx->FNGradL; ++j) {
                    ctx->FNPGrad[i][j] = ctx->FNGrad[0][j];
                    ctx->FNPGrad1[i][j] = ctx->FNGrad[1][j];
                }
                ctx->FNPEATyp[i] = (short)(ctx->FNEVATyp);
            }
        }
    }
   
    /* then evaluate function */
           
    err = v_eval1(ctx, n,ctx->FNECnt,ctx->FNETyp,ctx->FNEVal,ctx->ESIdx,x,0,deriv,0,0,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_icval(fl,fu,n,xl,xu)                                                */
/*  ##                                                                      */
/*  Calculates an inclusion function for a constraint function.             */  
/*  Function arguments are [ xl[i], xu[i] ], i = 0,...,n - 1.               */
/*  Return interval in [fl,fu].                                             */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_icval(TDAContext *ctx, double *fl,double *fu,int n,double *xl,double *xu)           
{
    register int i;
    int err;
    double xx[2];

    for (i = 0; i < n; ++i) {       /* set values for evaluation */
        ctx->FNArgVal[i] = xl[i];
        ctx->FNArgVal1[i] = xu[i];
    }

    /* first evaluate intermediate parameters */

    for (i = 0; i < ctx->FNPN; ++i) {
        if (ctx->FNPECnt[i] > 0) {
            err = v_eval1(ctx, 0,ctx->FNPECnt[i],ctx->FNPETyp[i],ctx->FNPEVal[i],ctx->ESIdx,xx,0,0,0,0,1);
            if (err)  
                return(err);
            ctx->FNPVal[i] =  xx[0];
            ctx->FNPVal1[i] = xx[1];
        }
    }
   
    /* then evaluate constraint function */
           
    err = v_eval1(ctx, 0,ctx->FNECCnt,ctx->FNECTyp,ctx->FNECVal,ctx->ESIdx,xx,0,0,0,0,1);
    if (err == 0) {
        *fl = xx[0];
        *fu = xx[1];
    }
    return(err);
}

