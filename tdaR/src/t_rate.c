/****************************************************************************/
/*  t_rate                                                                  */
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
#include "t_var.h"
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_edat.h"
#include "t_ml.h"
#include "t_gf.h"
#include "t_min.h"
#include "t_cdf.h"
#include "t_prate.h"
#include "t_fnrta.h"
#include "t_fnrtb.h"
#include "t_con.h"
#include "t_int.h"
#include "tda_context.h"

/*  functions in t_rate.c */

int rate(TDAContext *ctx);
int prn_rmod(TDAContext *ctx, int n);
int alloc_vdef(TDAContext *ctx, int mode);
int get_rvar(TDAContext *ctx, char *pcmd);
int save_vdef(TDAContext *ctx, char *s);
int check_mod(TDAContext *ctx);
int build_pvec(TDAContext *ctx);
void exp_null(TDAContext *ctx);
void ml_rstart(TDAContext *ctx);
void prn_rcoeff(TDAContext *ctx, int mode);
#ifdef TDA_R_PACKAGE
static void rt_explab(TDAContext *ctx, int sn, int org, int des, int mt,
                      int tp, const char *name);
#endif
void prn_rc(TDAContext *ctx, int k,int mode);
void prn_resid(TDAContext *ctx, int mod);
void prn_resid1(TDAContext *ctx, int icase, int org, int des, double ts, double tf, double rate, double surv, double resid, double wt);
int prn_tper(TDAContext *ctx);
int rt_alloc(TDAContext *ctx, int opt);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */





/* ------------------------------------------------------------------------ */
/*  rate()          Transition rate models. Command in CmdBuf.              */
/*                  Return 0 if OK, otherwise -1.                           */

int rate(TDAContext *ctx)
{
    int err;

    err = -1;         
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Transition rate models. Current memory: %d bytes.\n",ctx->MemReq);

    if (ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(0);
    }
    set_mldef(ctx);                /* set defaults for ML estimation */

    if (parm(ctx, ctx->CmdBuf + 4,3,1))   /* get parameters */
        goto RTFin;

    ctx->MixTyp = 0;
    if (ctx->PMMIX)
        ctx->MixTyp = 1;

    /* check option for derivatives with algorithms 7 and 8 */

    if (ctx->MINA == 7 || ctx->MINA == 8) {
        if (ctx->PMDOPT < 0 || ctx->PMDOPT > 2)
            ctx->PMDOPT = 2;
    }

    if (prn_rmod(ctx, ctx->PMRHSI))       /* check/print type of model */
        goto RTFin;

    if (ctx->RMTyp == MEXP1 || ctx->RMTyp == MEXP2) {        /* time periods */
        if (ctx->PMNTP == 0) {
            p_err(ctx, -17,1);
            goto RTFin;
        }
        if (fabs(ctx->PMTP[0]) > ctx->EPSI) {      /* must begin with zero */
            p_err(ctx, -33,1);
            goto RTFin;
        }
        if (prn_tper(ctx))
            goto RTFin;
    }
    if (alloc_vdef(ctx, 1)) {        /* allocate memory for model terms */
        p_err(ctx, -2,1);
        goto RTFin;
    }
    if (ctx->PMNXA) {                /* get covariable definitions */
        if (get_rvar(ctx, ctx->CmdBuf))
            goto RTFin;
    }
    if (ctx->NDef == 0)  
        printf1(ctx, "Model without covariates.\n\n");
       
    if (check_mod(ctx))            /* check model definition */
        goto RTFin;

    if (build_pvec(ctx))           /* create list of model parameters */
        goto RTFin;

    set_mlopt(ctx);                /* adjust ML options */         

    /* allocate additional memory */

    if (rt_alloc(ctx, 1))
        goto RTFin;

    /* for transition rate model we currently have only the standard form
       of covariance matrix calculation. */

    ctx->CCTyp = 2;
    if (ml_init(ctx, 1,0,0))             /* init ML estimation */
        goto RTFin;

    exp_null(ctx);         /* calculate likelihood of exp null model */
   
    printf1(ctx, "Log-likelihood of exponential null model: %lg\n",ctx->FMax0);
#ifdef TDA_R_PACKAGE
    tda_export_mat(ctx, "rate.logLik.null", &ctx->FMax0, 1, 1);
#endif
    if (ctx->DScalFlg == 0)
        printf1(ctx, "Changed scaling factor for log-likelihood: %lg\n",ctx->DScal);

    if (get_dsv(ctx, ctx->NParm,ctx->Par,0,ctx->ParLB,ctx->ParUB,1))   /* try to get starting values */
        goto RTFin;
       
    if (ctx->DSVFlg == 0)
        ml_rstart(ctx);                /* default starting values */

    prot_init(ctx, 1,ctx->RMTyp,0);           /* init protocol file */

    prn_cwt(ctx);                      /* case weight information */
       
    if (ctx->PMNConS > 0) {              /* process constraints */
        newline(ctx);
        if (con_proc(ctx, ctx->CmdBuf))
            goto RTFin;
    }

    if (ctx->RMTyp == MCOX && ctx->ESortFlg == 0) {        /* sort data */
        if (sort_ed(ctx, 1))  
            goto RTFin;
    }

    /* -------------------------------------------------------------------- */
    /*  call ffmin() for maximization of log likelihood.                    */

    if (ffmin(ctx, ctx->RMTyp,0,1,1)) /* minimization */
        goto RTFin;         /* insuff memory or fn error */
       
    prn_mlres(ctx, 0);           /* print info about min algorithm */
    prn_rcoeff(ctx, 0);          /* print parameter estimates */

    if (ctx->LConv >= 0) {       /* if convergence */

        prn_ml1res(ctx, 0);      /* additional ML results */

        if (ctx->PMRRFlg)        /* print relative risks  */
            prn_rcoeff(ctx, 1);

        if (ctx->PMResFDef)          /* generalized residuals */
            prn_resid(ctx, ctx->RMTyp);

        if (ctx->PMPRN > 0)          /* prate option */
            prate(ctx, ctx->CmdBuf,ctx->RMTyp);

        if (ctx->RMTyp == MCOX && ctx->PMNTP > 0)     /* goodness of fit for Cox models */
            test_cox(ctx);
    }
    err = 0;

RTFin:
    if (ctx->ESortFlg)
        sort_ed(ctx, 0);     
    rt_alloc(ctx, 0);
    con_free(ctx);
    ml_init(ctx, 0,0,0);     
    alloc_vdef(ctx, 0);      
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_rmod(n)     Print name of rate model n and set RMTyp. If not        */
/*                  available error message and RMTyp = 0.                  */
/*                  Return 0 if OK, -1 if error                             */

int prn_rmod(TDAContext *ctx, int n)
{
    printf1(ctx, "\nModel: ");

    switch (n) {

        case  MCOX:     printf1(ctx, "Cox (partial likelihood)\n"); 
                        if (ctx->PM1NV > 0)  
                            printf1(ctx, "Stratified with %d groups.\n",ctx->PM1NV);
                        break;
        case  MEXP:     printf1(ctx, "Exponential\n");         
                        break;
        case  MEXP1:    printf1(ctx, "Exponential with time-periods\n"); 
                        break;
        case  MEXP2:    printf1(ctx, "Exponential with period-specific effects\n");          
                        break;
        case  MPOL:     printf1(ctx, "Polynomial (I) with degree %d\n",ctx->PMDEG);        
                        break;
        case  MPOL1:    printf1(ctx, "Polynomial (II) with degree %d\n",ctx->PMDEG);         
                        if (ctx->PMDEG > 0) {
                            printf1(ctx, "Numerical integration with method %d. Relative error: %lg (%lg).\n",
                                                    ctx->NINTMETH,ctx->NINTRERR,ctx->NINTAERR);
                        }
                        break;
        case  MWEI:     printf1(ctx, "Weibull\n");
                        break;
        case  MLL:      printf1(ctx, "Log-logistic (I)\n");
                        break;
        case  MLL2:     printf1(ctx, "Log-logistic (II)\n");        
                        break;
        case  MLN:      printf1(ctx, "Log-normal\n");        
                        break;
        case  MIG:      printf1(ctx, "Inverse Gaussian\n");         
                        break;
        case  MSIC:     printf1(ctx, "Sickle\n");         
                        break;
       
        case  MGAM:     printf1(ctx, "Generalized Gamma (kgam=%lg).\n",ctx->PMKGam);          
                        break;
    
        case  MGM:      printf1(ctx, "Gompertz-Makeham\n");           
                        break;
  
        case  DLR:      printf1(ctx, "Discrete time logistic regression, degree %d.\n",ctx->PMDEG);          
                        break;
        case  CLL:      printf1(ctx, "Discrete time complementary log-log, degree %d.\n",ctx->PMDEG);         
                        break;

        default:        ++ctx->ErrCnt;
                        printf1(ctx, "%d not defined.\n\n",n);
                        return(-1);
    }
                 
    if (ctx->MixTyp == 1)  
        printf1(ctx, "Gamma mixture.\n");
    ctx->RMTyp = n;
    newline(ctx);             
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_vdef(mode)    If mode != 0 allocate                               */
/*                      VTerm, VOrg, VDes, VSnum for a maximum of MaxVP     */
/*                      entries and set VDAlloc = 1; otherwise free the     */
/*                      memory and set VDAlloc = 0.                         */
/*                                                                          */
/*                      In addition, allocate/free:                         */
/*                      VNVar with number of variables,                     */  
/*                      PIdx with variable numbers for model parameters     */
/*                      PIdxM with corresponding model terms.               */
/*                      PIdxPtr with pointers to trans-specific parameters  */
/*                                                                          */
/*                      TranTVar[i] = 1 for VTyp 5 variables                */
/*                      i = 0,...,NTran1 - 1.                               */
/*                                                                          */
/*                      Return 0 if successful, -1 if insufficient memory.  */

int alloc_vdef(TDAContext *ctx, int mode)
{
    register int i;
    int err = -1;

    if (mode == 0) {
        if (ctx->VDAlloc == 0)
            return(0);

        for (i = 0; i < ctx->NDef; ++i) {
            if (ctx->VNVar[i] > 0) {
                free((char *)ctx->VNIVar[i]);
                memrq(ctx, -ctx->VNVar[i],sizeof(short));
                ctx->VNVar[i] = 0;
            }
        }
        ctx->NDef = err = 0;
        ctx->XMFlg[1] = ctx->XMFlg[2] = ctx->XMFlg[3] = ctx->XMFlg[4] = 0;
        goto VDA10;
    }

    if (!(ctx->TranTVar = (short *)calloc((size_t)(ctx->NTran1),sizeof(short))))   
        goto VDA00;           
    memrq(ctx, ctx->NTran1,sizeof(short));

    if (!(ctx->VTNum = (short *)calloc((size_t)(ctx->NTran1),sizeof(short))))   
        goto VDA0;           
    memrq(ctx, ctx->NTran1,sizeof(short));

    if (!(ctx->VTerm = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA1;           
    memrq(ctx, MaxVP,sizeof(short));
              
    if (!(ctx->VOrg = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA2;           
    memrq(ctx, MaxVP,sizeof(short));
              
    if (!(ctx->VDes = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA3;           
    memrq(ctx, MaxVP,sizeof(short));
              
    if (!(ctx->VSNum = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA4;           
    memrq(ctx, MaxVP,sizeof(short));
              
    if (!(ctx->VNVar = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA5;           
    memrq(ctx, MaxVP,sizeof(short));
              
    if (!(ctx->VNIVar = (short **)calloc(MaxVP,sizeof(short *))))   
        goto VDA6;           
    memrq(ctx, MaxVP,sizeof(short *));
              
    if (!(ctx->PIdx = (short *)calloc(MaxP + 2,sizeof(short))))   
        goto VDA7;           
    memrq(ctx, MaxP + 2,sizeof(short));
              
    if (!(ctx->PIdxM = (char *)calloc(MaxP + 2,sizeof(char))))   
        goto VDA8;           
    memrq(ctx, MaxP + 2,sizeof(char));
              
    if (!(ctx->PIdxPtr = (short *)calloc((size_t)(ctx->NTran1),sizeof(short))))   
        goto VDA9;           
    memrq(ctx, ctx->NTran1,sizeof(short));
              
    ctx->VDAlloc = 1;
    return(0);

VDA10:
    free((char *)ctx->PIdxPtr);
    memrq(ctx, -ctx->NTran1,sizeof(short));
VDA9:
    free((char *)ctx->PIdxM);
    memrq(ctx, -MaxP - 2,sizeof(char));
VDA8:
    free((char *)ctx->PIdx);
    memrq(ctx, -MaxP - 2,sizeof(short));
VDA7:
    free((char *)ctx->VNIVar);
    memrq(ctx, -MaxVP,sizeof(short *));
VDA6:
    free((char *)ctx->VNVar);
    memrq(ctx, -MaxVP,sizeof(short));
VDA5:
    free((char *)ctx->VSNum);
    memrq(ctx, -MaxVP,sizeof(short));
VDA4:
    free((char *)ctx->VDes);
    memrq(ctx, -MaxVP,sizeof(short));
VDA3:
    free((char *)ctx->VOrg);
    memrq(ctx, -MaxVP,sizeof(short));
VDA2:
    free((char *)ctx->VTerm);
    memrq(ctx, -MaxVP,sizeof(short));
VDA1:
    free((char *)ctx->VTNum);
    memrq(ctx, -ctx->NTran1,sizeof(short));
VDA0:
    free((char *)ctx->TranTVar);
    memrq(ctx, -ctx->NTran1,sizeof(short));
VDA00:
    ctx->VDAlloc = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_rvar(pcmd)  get covariates for transition rate model.               */
/*  ##              Return 0 if OK, otherwise -1.                           */

int get_rvar(TDAContext *ctx, char *pcmd)
{
    int err;
    register char c,*p,*q;

    err = 0;
    p = pcmd;
    while (*p) {
        if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) || 
            !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3)) { 

            q = skip_xa(ctx, p);
            c = *q;
            *q = '\0';
            err = save_vdef(ctx, p);
            if (err)
                break;
            *q = c;
            p = q;
        }
        p++;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  save_vdef(s)    s is a string of the form xa(org,des)=variable_list,    */
/*                  or xa(org,des,sn) = variable_list; and analogously for  */
/*                  b,c,d model terms. This function saves the definition   */
/*                  in:                                                     */
/*                                                                          */
/*                  VTerm[NDef] = model term (a=1,b=2,c=3,d=4)              */
/*                  VSNum[NDef] = spell number (or 1 if not defined)        */  
/*                  VOrg[NDef]  = origin state                              */
/*                  VDes[NDef]  = destination state                         */  
/*                  VNVar[NDef] = number of variables                       */
/*                  VNIVar[NDef][i] = indices of variables                  */
/*                                                                          */
/*                  XMFlg[1...4] is an array with flags for model terms.    */
/*                                                                          */  
/*                  If NDef >= MaxVP then error!                            */
/*                                                                          */
/*                  Return 0 if successful, -1 if error.                    */

int save_vdef(TDAContext *ctx, char *s)
{
    register int i;
    int n,sn,org,des,nb;
    register char *p;

    if (ctx->NDef + 3 >= MaxVP) {
        printf1(ctx, "Exceeded maximum number of vector definitions in model specification.\n");
        return(-1);
    }
    p = s;
    if (*p++ != 'x')
        return(-1);

    i = 0;
    switch (*p) {
        case 'a':   i = 1; break;
        case 'b':   i = 2; break;
        case 'c':   i = 3; break;
        case 'd':   i = 4; break;
        default:    goto SVD1;
    }
    ctx->XMFlg[i] = 1;
    ctx->VTerm[ctx->NDef] = (short)(i);

    if (*++p != '(')
        goto SVD1;

    if (sscanf(p,"(%d,%d,%d)",&org,&des,&sn) == 3) ;
    else if (sscanf(p,"(%d,%d)",&org,&des) == 2)
        sn = 1;
    else
        goto SVD1;

    n = 0;
    for (i = 0; i < ctx->NTran1; ++i) {  /* check if sn,org defined */
        if (des == ctx->DesTran1[i] && org == ctx->OrgTran1[i] && sn == ctx->SnTran1[i]) {
            n = 1;
            break;
        }
    }
    if (n == 0) {
        printf1(ctx, "Error: %s\nNo transition with (org,des,sn) = (%d,%d,%d).\n",s,org,des,sn);
        return(-1);
    }
    ctx->VOrg[ctx->NDef] = (short)(org);
    ctx->VDes[ctx->NDef] = (short)(des);
    ctx->VSNum[ctx->NDef] = (short)(sn);

    while (*p++ && *p != ')');
    if (*p++ != ')' || *p++ != '=')
        goto SVD1;

    if (*p == '1') {            /* only constant */
        ctx->VNVar[ctx->NDef++] = 0;
        return(0);
    }
    p = get_nvia(ctx, p,&n,0,&nb);
    if (n <= 0) {
        printf1(ctx, "Error: %s\nDefinition contains undefined variables.\n",s);
        return(-1);
    }
    if (nb) {
        p_err(ctx, -42,1);
        return(-1);
    }
    ctx->VNVar[ctx->NDef] = (short)(n);

    if (!(ctx->VNIVar[ctx->NDef] = (short *)calloc((size_t)(n),sizeof(short)))) { 
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, n,sizeof(short));
    for (i = 0; i < n; ++i)
        ctx->VNIVar[ctx->NDef][i] = ctx->VLVIdx[i]; 

    ctx->NDef++;
    alloc_vl(ctx, 0);
    return(0);

SVD1:
    printf1(ctx, "Syntax error: %s\n",s);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  check_mod       Check rate model specification.                         */
/*                  Return 0 if OK, -1 if error.                            */

int check_mod(TDAContext *ctx)
{
    switch (ctx->RMTyp) {
        case MCOX:  if (!ctx->XMFlg[1])
                        goto CMErr1;
                    if (ctx->XMFlg[2])
                        goto CMErr4;
                    if (ctx->XMFlg[3])
                        goto CMErr5;
                    break;

        case MEXP:  if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2])
                        goto CMErr4;
                    if (ctx->XMFlg[3]) 
                        goto CMErr5;
                    break;
    
        case MEXP1: 
        case MEXP2: if (ctx->PMNTP == 0)   
                        goto CMErr9; 
                    if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2]) 
                        goto CMErr4;
                    if (ctx->XMFlg[3])
                        goto CMErr5;
                    break;
   
        case MPOL:  if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2]) 
                        goto CMErr4;
                    if (ctx->XMFlg[3]) 
                        goto CMErr5;
                    if (ctx->PMDEG)
                        ctx->XMFlg[2] = 1;
                    else
                        ctx->XMFlg[2] = 0;
                    break;

        case MPOL1: if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2])
                        goto CMErr4;
                    if (ctx->XMFlg[3]) 
                        goto CMErr5;
                    if (ctx->PMDEG)
                        ctx->XMFlg[2] = 1;
                    else
                        ctx->XMFlg[2] = 0;
                    break;
    
        case MWEI:  
        case MLL:
        case MLL2:  
        case MLN:   
        case MIG:   
        case MSIC:  
        case MGAM:  if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (!ctx->XMFlg[2])
                        ctx->XMFlg[2] = 1;
                    if (ctx->RMTyp != MLL2 && ctx->RMTyp != MLN) {
                        if (ctx->XMFlg[3])
                            goto CMErr5;
                    }
                    else if (ctx->RMTyp != MLN && !ctx->XMFlg[3])
                        ctx->XMFlg[3] = 1;

                    if (ctx->RMTyp == MLN && !ctx->XMFlg[3] && ctx->PMNI != 0)
                        ctx->XMFlg[3] = 1;

                    break;
         
        case MGM:   if (!ctx->XMFlg[1] && !ctx->XMFlg[2] && !ctx->XMFlg[3])
                        ctx->XMFlg[2] = 1;
                    if (ctx->XMFlg[2] && !ctx->XMFlg[3])
                        ctx->XMFlg[3] = 1;
                    else if (!ctx->XMFlg[2] && ctx->XMFlg[3])
                        ctx->XMFlg[2] = 1;
                    /***
                    if (!XMFlg[1] && !XMFlg[2] && !XMFlg[3])
                        XMFlg[1] = 1;
                    ***/
                    if (!ctx->XMFlg[1] && ctx->PMNI != 0)
                        ctx->XMFlg[1] = 1;
                    break;
         
        case DLR:   if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2])
                        goto CMErr4;
                    if (ctx->XMFlg[3]) 
                        goto CMErr5;
                    if (ctx->PMDEG)
                        ctx->XMFlg[2] = 1;
                    else
                        ctx->XMFlg[2] = 0;
                    ctx->DisFlg = 1;
                    break;

        case CLL:   if (!ctx->XMFlg[1])
                        ctx->XMFlg[1] = 1;
                    if (ctx->XMFlg[2]) 
                        goto CMErr4;
                    if (ctx->XMFlg[3]) 
                        goto CMErr4;
                    if (ctx->PMDEG)
                        ctx->XMFlg[2] = 1;
                    else
                        ctx->XMFlg[2] = 0;
                    ctx->DisFlg = 1;
                    break;

        default:    return(-1);
    }
    if (ctx->MixTyp == 1) {
        if (ctx->DisFlg || ctx->RMTyp == MEXP1  || ctx->RMTyp == MPOL   || ctx->RMTyp == MPOL1 ||
                      ctx->RMTyp == MLL2   || ctx->RMTyp == MGAM   ||
                      ctx->RMTyp == MIG    || ctx->RMTyp == MGM    || ctx->RMTyp == MEXP2 ||
                      ctx->RMTyp == MCOX   || (ctx->RMTyp == MLN   && ctx->XMFlg[3])) {

            goto CMErr7;
        }
        if (!ctx->XMFlg[4])  
            ctx->XMFlg[4] = 1;
    }
    else if (ctx->XMFlg[4])  
        goto CMErr8;
         
    return(0);

CMErr1:     printf1(ctx, "Error: model 1 (partial likelihood) requires covariates.\n");
            return(-1);
CMErr4:     printf1(ctx, "Error: B term variables not possible.\n");
            return(-1);
CMErr5:     printf1(ctx, "Error: C term variables not possible,\n");
            return(-1);
CMErr7:     printf1(ctx, "Error: Gamma mixture not possible.\n");
            return(-1);
CMErr8:     printf1(ctx, "Error: D term variables possible only with mixtures.\n");
            return(-1);
CMErr9:     printf1(ctx, "Error: this model needs time periods.\n");
            return(-1);
}

/* ------------------------------------------------------------------------ */
/*  build_pvec. Calculate NParm, the number of model parameters. And        */
/*              create references to model parameters:                      */
/*                                                                          */
/*      PIdx(i)         i = 1,NParm. Contains the number of the variable    */
/*                      associated with the i.th parameter.                 */
/*      PIdxM(i)        i = 1,NParm. Contains the number of the model term  */
/*                      associated with the i.th parameter.                 */
/*      PIdxPtr[j]      pointer to begin of parameters for transition j.    */
/*      TranTVar[j]     set to 1 if there are type 5 varibales for          */
/*                      transition j.                                       */
/*      VTNum[j]        number of covariates for transition j.              */
/*      VTMax           maximum of VTNum[j]                                 */
/*                                                                          */
/*      Return 0 if OK, -1 if error.                                        */

int build_pvec(TDAContext *ctx)
{
    register int ii,i,j,k,kk,l;
    int nv,sn,org,des,kkn,first,npx;

    ctx->VTMax = 0;
    kkn = 1;

    if (ctx->RMTyp == MEXP2)
        kkn = ctx->PMNTP;
      
    ctx->NParm = 0;          /* number of model parameters */

    for (ii = 0; ii < ctx->NTran1; ++ii) {       /* for all transitions */

        npx = ctx->NParm;

        ctx->TranTVar[ii] = 0;           /* set for VTyp 5 variables */
        sn = ctx->SnTran1[ii];
        org = ctx->OrgTran1[ii];
        des = ctx->DesTran1[ii];
        first = 1;
        nv = 0;

        for (k = 1; k <= ctx->MTerm; ++k) {      

            if (ctx->XMFlg[k]) {        

                if (first) {
                    ctx->PIdxPtr[ii] = (short)(ctx->NParm + 1);
                    first = 0;
                }   

                /* first insert a constant, skip Cox models */

                if (ctx->RMTyp != MCOX) {    
                    if (++ctx->NParm > MaxP) 
                        goto BPErr1;

                    ctx->PIdx[ctx->NParm] = 0;
                    ctx->PIdxM[ctx->NParm] = (char)(k);
                }

                /* add period-specific constants */
         
                if (ctx->RMTyp == MEXP1 || ctx->RMTyp == MEXP2) {
                    for (i = 1; i < ctx->PMNTP; ++i) {
                        if (++ctx->NParm > MaxP)  
                            goto BPErr1;

                        ctx->PIdx[ctx->NParm] = (short)(-i);
                        ctx->PIdxM[ctx->NParm] = (char)(k);
                    }
                }
                if (k == 2 && (ctx->RMTyp == MPOL || ctx->RMTyp == MPOL1 ||
                               ctx->RMTyp == DLR  || ctx->RMTyp == CLL)) {   

                    for (i = 1; i < ctx->PMDEG; ++i) {
                        if (++ctx->NParm > MaxP)  
                            goto BPErr1;

                        ctx->PIdx[ctx->NParm] = (short)(-i);
                        ctx->PIdxM[ctx->NParm] = (char)(k);
                    }
                }
     
                /* now insert indices for the user defined covariates */

                for (i = 0; i < ctx->NDef; ++i) {

                    if ((int)ctx->VTerm[i] == k && (int)ctx->VOrg[i] == org 
                                           && (int)ctx->VDes[i] == des 
                                           && (int)ctx->VSNum[i] == sn) {

                        for (j = 0; j < ctx->VNVar[i]; ++j) {
                    
                            /* repeat only for MEXP2 model */           
  
                            for (kk = 1; kk <= kkn; ++kk) {

                                if (++ctx->NParm > MaxP)  
                                    goto BPErr1;
     
                                l = ctx->VNIVar[i][j];
                                nv++;

                                if (ctx->RMTyp == MCOX && !j)
                                    ctx->PIdx[ctx->NParm] = (short)(-l - 1);
                                else
                                    ctx->PIdx[ctx->NParm] = (short)(l + 1);

                                ctx->PIdxM[ctx->NParm] = (char)(k);
                                if (ctx->VTyp[l] == 5)  /* VTyp 5 variable */
                                    ctx->TranTVar[ii] = 1;
                            }   
                        }
                    }
                }
            }
        }
        ctx->VTNum[ii] = (short)(nv);
        if (ctx->VTMax < nv)
            ctx->VTMax = nv;

        if (ctx->NParm == npx)           
            goto BPErr2;

    }
    if (ctx->RMTyp == MCOX)
        ctx->PIdx[ctx->NParm + 1] = -1;
    else
        ctx->PIdx[ctx->NParm + 1] = 0;

    ctx->PIdxM[ctx->NParm + 1] = 0;
    if (ctx->NParm <= 0) {
        printf1(ctx, "Error: number of model parameters is %d.\n",ctx->NParm);
        return(-1);
    }

    /***********************************************************
    printf1(ctx, "\n  PIdx : ");
    for (i = 1; i <= NParm + 1 ; ++i) printf1(ctx, " %2d",PIdx[i]);
    printf1(ctx, "\n  PIdxM: ");
    for (i = 1; i <= NParm + 1; ++i) printf1(ctx, " %2d",PIdxM[i]);
    printf1(ctx, "\nPIdxPtr\n");
    for (i = 0; i < NTran1; ++i)
        printf1(ctx, "%2d ",PIdxPtr[i]);
    printf1(ctx, "\n");
    ***********************************************************/
    return(0);

BPErr1: 
    printf1(ctx, "Error: exceeded max number of parameters (%d)\n",MaxP);
    return(-1);
BPErr2: 
    printf1(ctx, "Error: need at least one parameter (covariate) for each transition.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  exp_null.   Calculate FMax0, the log-likelihood of an exponential       */
/*              null model. Also calculate the scaling factor DScal.        */
/*              DScal is only calculated if DScalFlg = 0.                   */

void exp_null(TDAContext *ctx)
{
    register int i;
    double w;

    ctx->FMax0 = 0.0;
    for (i = 0; i < ctx->NTran1; ++i) {
        if ((double)(ctx->TWFreq[i]) > 0.0 && (double)(ctx->TWDur[i]) > 0.0)  
            ctx->FMax0 += (double)ctx->TWFreq[i] * (rlog(ctx, (double)((ctx->TWFreq[i]))/(double)(ctx->TWDur[i])) - 1.0);
    }
    if (ctx->DScalFlg == 0) {
        ctx->DScal = -1.0;                                                     
        w = fabs(ctx->FMax0);
        while (w < 0.5) {
            w *= 10.0;
            ctx->DScal *= 10.0;
        }
        while (w > 5.0) {
            w /= 10.0;
            ctx->DScal /= 10.0;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ml_rstart.      Calculation of default starting values.                 */

void ml_rstart(TDAContext *ctx) 
{
    register int i,j,k;
    double tmp,freq,dur,rate,lrate;

    printf1(ctx, "Using default starting values.\n");

    for (j = 1; j <= ctx->NParm; ++j)
        ctx->Par[j] = 0.0;

    if (ctx->RMTyp == MCOX)        /* Cox model has zero starting values */
        return;
            
    j = 1;
    for (i = 0; i < ctx->NTran1; ++i) {

        freq = (double)(ctx->TWFreq[i]);
        dur = (double)(ctx->TWDur[i]);
        if (freq > 0.0 && dur > 0.0) {
            rate = freq / dur; 
            lrate = rlog(ctx, rate); 
        }
        else
            rate = lrate = 0.0;

        for (k = 1; k <= ctx->MTerm; ++k) {

            if (ctx->XMFlg[k]) {

                tmp = 0.0;
                if (k == 1) {
                    if (ctx->RMTyp == MLN || ctx->RMTyp == MGAM)
                        tmp = -lrate;
                    else if (ctx->RMTyp == MSIC)
                        tmp = 2.0 * lrate;
                    else if (ctx->RMTyp == MIG)
                        tmp = rate;
                    else if (ctx->RMTyp != MLL2)
                        tmp = lrate;
                }
                else if (k == 2) {
                    if (ctx->RMTyp == MGM)
                        tmp = lrate;
                    else if (ctx->RMTyp == MSIC)
                        tmp = -lrate;
                    else if (ctx->RMTyp == MIG)
                        tmp = rlog(ctx, sqrt(rate));
                }
                else if (k == 3) {
                    if (ctx->RMTyp == MGM) 
                        tmp = -rate;
                    else if (ctx->RMTyp == MLL2)
                        tmp = lrate;
                }
                ctx->Par[j] = tmp;
                if (ctx->RMTyp == MEXP1 || ctx->RMTyp == MEXP2) {
                    while (ctx->PIdx[++j] < 0)  
                        ctx->Par[j] = tmp;
                    j--;
                }
                while (ctx->PIdx[++j])  
                    ;
            }
        }
    }
}

#ifdef TDA_R_PACKAGE
/* one "sn|org|des|MT|name" label per printed row lead, staged beside
   each branch's own printf (rate.est.names / rate.rrisk share the row
   order via prn_rc) */
static void rt_explab(TDAContext *ctx, int sn, int org, int des, int mt,
                      int tp, const char *name)
{
    char b[160];

    if (ctx->RMTyp == MEXP2)
        snprintf(b, sizeof(b), "%d|%d|%d|%c|P%d|%s", sn, org, des,
                 (char)('A' + mt - 1), tp, name);
    else
        snprintf(b, sizeof(b), "%d|%d|%d|%c|%s", sn, org, des,
                 (char)('A' + mt - 1), name);
    tda_export_str_row(ctx, "rate.est.names", b);
}
#endif

/* ------------------------------------------------------------------------ */
/*  prn_rcoeff(mode)                                                        */
/*  ##  Print estimated parameters, standard errors, and test statistics    */
/*      of rate model estimation. Print format is PMTFmt.                   */
/*      If mode == 0 print parameter estimates.                             */  
/*      If mode == 1 print relative risks.                                  */

void prn_rcoeff(TDAContext *ctx, int mode)
{
    register int i,j,k,sn,idx;             
    int lflag,tp,len,m,org,des;

    if (mode == 1) {
        printf1(ctx, "Estimated relative risks");
        switch (ctx->RMTyp) {
            case MCOX:  
            case MEXP:  
            case MEXP1: 
            case MEXP2:   
            /*  case MPOL:  */
            case MPOL1: 
            case MGM:   
            case MWEI:  
            case MSIC:  
            /*  case MLL:   */
            case MLL2:     
            /*  case MLN:   */
            /*  case MIG:   */
            /*  case MGAM:  */
            /*  case DLR:   */
            /*  case CLL:   */
                        break;
            default:    printf1(ctx, ": not available for this model.\n\n");
                        return;
        }
        printf1(ctx, ".\n\n");
    }
    len = 0;
    printf1(ctx, "Idx SN Org Des MT ");          
    if (ctx->RMTyp == MEXP2) {
        printf1(ctx, " P ");         
        len = 3;
    }
    printf1(ctx, "Variable ");        
    prnchar(ctx, ' ',ctx->VNameLen - 8,0);

    lflag = 0;      /* check for var labels */

    for (i = 1; i <= ctx->NParm; ++i) {
        j = ctx->PIdx[i];   
        if (j > 0 || (j < 0 && ctx->RMTyp == MCOX)) {
            if (j < 0)
                j = -j;
            if (ctx->VLabel[j - 1] != NULL) {
                lflag = 1;
                break;
            }
        }
    }
    if (lflag) {
        printf1(ctx, "Label ");
        prnchar(ctx, ' ',ctx->VLabelLen - 5,0);
    }
    if (mode == 1) {
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0); printf1(ctx, "    R.Risk");        
        len += 18;
        m = 1;
    }
    else {
        prnchar(ctx, ' ',ctx->PMTFmt1 - 10,0); printf1(ctx, "     Coeff"); 
        prnchar(ctx, ' ',ctx->PMTFmt1 -  9,0); printf1(ctx, "     Error"); 
        prnchar(ctx, ' ',ctx->PMTFmt1 -  9,0); printf1(ctx, "   C/Error  Signif");
        len += 26;
        m = 3;
    }
    if (lflag)
        len += ctx->VLabelLen + 1;

    k = 1;

    for (i = 0; i < ctx->NTran1; ++i) {

        sn = ctx->SnTran1[i];
        org = ctx->OrgTran1[i];
        des = ctx->DesTran1[i];

        printf1(ctx, "\n");           
        prnchar(ctx, '-',len + m * (ctx->PMTFmt1 + 1) + ctx->VNameLen,0);

        for (j = 1; j <= ctx->MTerm; ++j) {

            if (ctx->XMFlg[j]) {

                tp = 1;  /* time period */

                while (1) {

                    /* relative risks only for selected model terms */

                    if (mode == 1) {
                        if (ctx->RMTyp == MPOL1 || ctx->RMTyp == MWEI || ctx->RMTyp == MSIC) {
                            if (j != 1)
                                goto PRCNXT;
                        }
                        else if (ctx->RMTyp == MGM) {
                            if (j != 2)
                                goto PRCNXT;
                        }
                        else if (ctx->RMTyp == MLL2) {
                            if (j != 3)
                                goto PRCNXT;
                        }
                    }
                    idx = ctx->PIdx[k];
                    if (idx || ctx->RMTyp != MCOX) {
                        printf1(ctx, "\n%3d %2d %3d %3d  %c ",k,sn,org,des,(char)('A'+j-1));

                        if (idx > 0 || ctx->RMTyp == MCOX) {

                            if (ctx->RMTyp == MEXP2)  
                                printf1(ctx, "%2d ",tp);

                            if (idx < 0 && ctx->RMTyp == MCOX)
                                idx = -idx;
                            prn_vname(ctx, idx - 1);
#ifdef TDA_R_PACKAGE
                            rt_explab(ctx, sn, org, des, j, tp,
                                      ctx->VName[idx - 1]);
#endif
                            if (lflag)
                                prn_vlabel(ctx, idx - 1);
                        }
                        else {
                            if (ctx->RMTyp == MEXP2)  
                                printf1(ctx, "%2d ",tp);
                                                 
                            if (ctx->RMTyp == MEXP1 || ctx->RMTyp == MEXP2) {
                                printf1(ctx, "Period-%-2d",-idx + 1);          
                                prnchar(ctx, ' ',ctx->VNameLen - 8,0);
#ifdef TDA_R_PACKAGE
                                {
                                    char pb[24];
                                    snprintf(pb, sizeof(pb), "Period-%d",
                                             -idx + 1);
                                    rt_explab(ctx, sn, org, des, j, tp, pb);
                                }
#endif
                            }
                            else if (j == 2 &&
                                (ctx->RMTyp == MPOL || ctx->RMTyp == MPOL1 ||
                                 ctx->RMTyp == DLR  || ctx->RMTyp == CLL)) {  

                                printf1(ctx, "Beta-%-2d",-idx + 1);         
                                prnchar(ctx, ' ',ctx->VNameLen - 6,0);
#ifdef TDA_R_PACKAGE
                                {
                                    char pb[24];
                                    snprintf(pb, sizeof(pb), "Beta-%d",
                                             -idx + 1);
                                    rt_explab(ctx, sn, org, des, j, tp, pb);
                                }
#endif
                            }
                            else {
                                printf1(ctx, "Constant");          
                                prnchar(ctx, ' ',ctx->VNameLen - 7,0);
#ifdef TDA_R_PACKAGE
                                rt_explab(ctx, sn, org, des, j, tp,
                                          "Constant");
#endif
                            }
                            if (lflag)
                                prnchar(ctx, ' ',ctx->VLabelLen + 1,0);
                        }
                        prn_rc(ctx, k,mode);  /* print coefficients */
                    }
PRCNXT:           
                    k++;
                    if (!ctx->PIdx[k] || (ctx->RMTyp == MCOX && ctx->PIdx[k] < 0))
                        break;
       
                    if (ctx->RMTyp == MEXP2) {
                        if (++tp > ctx->PMNTP)
                            tp = 1;
                    }
                }
            }
        }
    }
    printf1(ctx, "\n\n");
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "rate.est");
    tda_export_flush(ctx, "rate.rrisk");
    tda_export_str_flush(ctx, "rate.est.names");
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_rc                                                                  */
/*           Print coefficients, used by prn_rcoeff.                        */
/*           If mode == 1 print relative risks                              */

void prn_rc(TDAContext *ctx, int k,int mode)
{
    register int j;
    double tmp,tmp1,tmp2;
#ifdef TDA_R_PACKAGE
    double erow[5];
    erow[0] = (double)k;
    erow[2] = erow[3] = erow[4] = (double)NAN;
#endif

    if (mode == 1) {
        tmp = 1.0;
        if (ctx->RMTyp == MWEI) {
            j = k;
            while (ctx->PIdx[++j]) ;
            tmp = rexp(ctx, ctx->Par[j]);
        }
        rt_printf1_d(ctx, ctx->PMTFmtS,rexp(ctx, tmp * ctx->Par[k]));
#ifdef TDA_R_PACKAGE
        erow[1] = rexp(ctx, tmp * ctx->Par[k]);
        tda_export_row(ctx, "rate.rrisk", erow, 2);
#endif
        return;
    }
    else  
        rt_printf1_d(ctx, ctx->PMTFmtS,ctx->Par[k]);
#ifdef TDA_R_PACKAGE
    erow[1] = ctx->Par[k];
#endif
       
    if (ctx->LConv >= 0 && ctx->Diag[k] > 0.0) {
        tmp = sqrt(ctx->Diag[k]);
        rt_printf1_d(ctx, ctx->PMTFmtS,tmp);        
#ifdef TDA_R_PACKAGE
        erow[2] = tmp;
#endif
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
    /* the same row the rate table prints; --- as NaN */
    tda_export_row(ctx, "rate.est", erow, 5);
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_resid(mod)                                                          */
/*      Print residuals for model mod to output file PMResFd.               */

void prn_resid(TDAContext *ctx, int mod)
{
    if (ctx->NTran1 > 1 || ctx->MEFlg || ctx->MixTyp || mod == MCOX || mod == DLR || mod == CLL) {    
        printf1(ctx, "Cannot calculate generalized residuals; pres option ignored.\n");
        return;
    }
    printf1(ctx, "Calculating generalized residuals; will be written to: %s\n",ctx->PMResFName);

    fprintf(ctx->PMResFd,"# Generalized residuals.\n");
    fprintf(ctx->PMResFd,"# ");
    fprnchar(ctx, ctx->PMResFd,' ',16 + 4 * ctx->PMMFmt1 - 8,0); fprintf(ctx->PMResFd,"Survivor\n# ");
    fprintf(ctx->PMResFd,"Case Org Des ");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 2,0); fprintf(ctx->PMResFd,"TS");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 1,0); fprintf(ctx->PMResFd,"TF");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 3,0); fprintf(ctx->PMResFd,"Rate");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 7,0); fprintf(ctx->PMResFd,"Function ");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 8,0); fprintf(ctx->PMResFd,"Residual ");
    fprnchar(ctx, ctx->PMResFd,' ',ctx->PMMFmt1 - 6,0); fprintf(ctx->PMResFd,"Weight\n# ");
    fprnchar(ctx, ctx->PMResFd,'-',12 + 6 * (ctx->PMMFmt1 + 1),1); 

    ctx->LFunc = ctx->LGrad = ctx->LSec = 0;

    switch (mod) {

        case  MEXP:  fn_exp(ctx, 1);     /* exponential model */
                     break;
    
        case  MEXP1: fn_exp1(ctx, 1);    /* exp. model with time periods */
                     break;
     
        case  MEXP2: fn_exp2(ctx, 1);    /* exp. model with time periods, II */
                     break;     
     
        case  MPOL:  fn_pol(ctx, 1);     /* polynomial, I */
                     break;     

        case  MPOL1: fn_pol1(ctx, 1);    /* polynomial, II */
                     break;     
       
        case  MGM:   fn_gm(ctx, 1);      /* Gompertz model */
                     break;

        case  MWEI:  fn_wei(ctx, 1);     /* Weibull model */
                     break;

        case  MSIC:  fn_sic(ctx, 1);     /* sickle model */
                     break;
   
        case  MLL:   fn_ll(ctx, 1);      /* log-logistic model, I */
                     break;

        case  MLL2:  fn_ll2(ctx, 1);     /* log-logistic model, II */     
                     break;
     
        case  MLN:   fn_ln(ctx, 1);      /* log-normal model */
                     break;

        case  MIG:   fn_ig(ctx, 1);      /* inverse gaussian model */
                     break;

        case  MGAM:  fn_gam(ctx, 1);     /* gamma */
                     break;

        default:
        fprintf(ctx->PMResFd,"  Calculation in current version of TDA not possible.\n");
    }
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  prn_resid1                                                              */
/*      Print residuals to output file PMResFd.                            */

void prn_resid1(TDAContext *ctx, int icase, int org, int des, double ts, double tf, double rate, double surv, double resid, double wt)
{
    fprintf(ctx->PMResFd,"%6d ",icase);
    fprintf(ctx->PMResFd,"%3d %3d ",org,des);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,ts);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,tf);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,rate);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,surv);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,resid);
    rt_fprintf_d(ctx, ctx->PMResFd,ctx->PMMFmtS,wt);
    fprintf(ctx->PMResFd,"\n");
#ifdef TDA_R_PACKAGE
    /* the same nine values res.out gets, before its print format
       rounds them */
    {
        double erow[9];
        erow[0] = (double)icase;
        erow[1] = (double)org;
        erow[2] = (double)des;
        erow[3] = ts;
        erow[4] = tf;
        erow[5] = rate;
        erow[6] = surv;
        erow[7] = resid;
        erow[8] = wt;
        tda_export_row(ctx, "rate.residuals", erow, 9);
    }
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_tper()      Print table with time periods and number of ending      */
/*                  times and events.                                       */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int prn_tper(TDAContext *ctx)
{
    register int i;
    int err,icase,sn,org,des,spl,nspl;
    double ts,tf;

    err = -1;

    if (alloc_acn(ctx, ctx->PMNTP))
        goto TPERFin;

    if (alloc_aci(ctx, ctx->PMNTP))
        goto TPERFin;

    if (alloc_acj(ctx, ctx->PMNTP))
        goto TPERFin;
     
    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {
        for (i = ctx->PMNTP - 1; i >= 0; --i) {
            if (ts >= ctx->PMTP[i]) {
                ctx->AcN[i] += 1;
                break;
            }
        }
        for (i = ctx->PMNTP - 1; i >= 0; --i) {
            if (tf >= ctx->PMTP[i]) {
                ctx->AcI[i] += 1;
                if (des != org)
                    ctx->AcJ[i] += 1;
                break;
            }
        }
    }
    printf1(ctx, "    Time period        Starting times  Ending times    Events\n  ");         
    prnchar(ctx, '-',59,1);

    for (i = 0; i < ctx->PMNTP; ++i) {
        printf1(ctx, "%9.4lf - ",ctx->PMTP[i]);           
        if (i < ctx->PMNTP - 1)
            printf1(ctx, "%9.4lf",ctx->PMTP[i + 1]);
        else
            printf1(ctx, "         ");
        printf1(ctx, "  %14d  %12d %9d\n",ctx->AcN[i],ctx->AcI[i],ctx->AcJ[i]);
    }
    newline(ctx);
    err = 0;

TPERFin:
    alloc_acn(ctx, 0);       
    alloc_aci(ctx, 0);       
    alloc_acj(ctx, 0);       
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rt_alloc(opt)   if opt != 0 allocate additional work arrays for rate    */
/*                  models, otherwise free previously allocated memory.     */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int rt_alloc(TDAContext *ctx, int opt)
{
    int err = 0;

    if (opt == 0)
        goto RTAFin;

    err = -1;

    if (ctx->RMTyp == MPOL || ctx->RMTyp == MPOL1) {

        if (!(ctx->POLWrk1 = (double *)calloc((size_t)(2) * (size_t)(ctx->PMDEG) + 1,sizeof(double))))   
            goto RTAFin;         
        ctx->POLWrk1A = 2 * ctx->PMDEG + 1;
        memrq(ctx, ctx->POLWrk1A,sizeof(double));

        if (!(ctx->POLWrk2 = (double *)calloc((size_t)(ctx->PMDEG + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->POLWrk2A = ctx->PMDEG + 1;
        memrq(ctx, ctx->POLWrk2A,sizeof(double));

        if (!(ctx->BCBeta = (double *)calloc((size_t)(ctx->PMDEG + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->BCBetaA = ctx->PMDEG + 1;
        memrq(ctx, ctx->BCBetaA,sizeof(double));

        if (!(ctx->BCTG = (double *)calloc((size_t)(ctx->PMDEG + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->BCTGA = ctx->PMDEG + 1;
        memrq(ctx, ctx->BCTGA,sizeof(double));

    }
    if (ctx->RMTyp == DLR || ctx->RMTyp == CLL) {

        if (!(ctx->DRJKW = (double *)calloc((size_t)(ctx->NTran1 + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->DRJKWA = ctx->NTran1 + 1;
        memrq(ctx, ctx->DRJKWA,sizeof(double));
    }
    if (ctx->RMTyp == CLL) {

        if (!(ctx->DAJKW = (double *)calloc((size_t)(ctx->NTran1 + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->DAJKWA = ctx->NTran1 + 1;
        memrq(ctx, ctx->DAJKWA,sizeof(double));

        if (!(ctx->DBJKW = (double *)calloc((size_t)(ctx->NTran1 + 1),sizeof(double))))   
            goto RTAFin;         
        ctx->DBJKWA = ctx->NTran1 + 1;
        memrq(ctx, ctx->DBJKWA,sizeof(double));
    }

    return(0);

RTAFin:
    if (ctx->POLWrk1A) {
        free((char *)ctx->POLWrk1);
        memrq(ctx, -ctx->POLWrk1A,sizeof(double));
        ctx->POLWrk1A = 0;               
    }
    if (ctx->POLWrk2A) {
        free((char *)ctx->POLWrk2);
        memrq(ctx, -ctx->POLWrk2A,sizeof(double));
        ctx->POLWrk2A = 0;               
    }
    if (ctx->BCBetaA) {
        free((char *)ctx->BCBeta);
        memrq(ctx, -ctx->BCBetaA,sizeof(double));
        ctx->BCBetaA = 0;               
    }
    if (ctx->BCTGA) {
        free((char *)ctx->BCTG);
        memrq(ctx, -ctx->BCTGA,sizeof(double));
        ctx->BCTGA = 0;               
    }
    if (ctx->BCGamA) {
        free((char *)ctx->BCGam);
        memrq(ctx, -ctx->BCGamA,sizeof(double));
        ctx->BCGamA = 0;               
    }
    if (ctx->DRJKWA) {
        free((char *)ctx->DRJKW);
        memrq(ctx, -ctx->DRJKWA,sizeof(double));
        ctx->DRJKWA = 0;               
    }
    if (ctx->DAJKWA) {
        free((char *)ctx->DAJKW);
        memrq(ctx, -ctx->DAJKWA,sizeof(double));
        ctx->DAJKWA = 0;               
    }
    if (ctx->DBJKWA) {
        free((char *)ctx->DBJKW);
        memrq(ctx, -ctx->DBJKWA,sizeof(double));
        ctx->DBJKWA = 0;               
    }
    return(err);
} 




