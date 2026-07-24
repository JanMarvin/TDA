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

/*  functions in t_rate.c */

int rate(void);
int prn_rmod(int n);
int alloc_vdef(int mode);
int get_rvar(char *pcmd);
int save_vdef(char *s);
int check_mod(void);
int build_pvec(void);
void exp_null(void);
void ml_rstart(void);
void prn_rcoeff(int mode);
void prn_rc(int k,int mode);
void prn_resid(int mod);
void prn_resid1(int icase, int org, int des, double ts, double tf,
     double rate, double surv, double resid, double wt);
int prn_tper(void);
int rt_alloc(int opt);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int MTerm = 4;          /* max number of model terms                        */

int RMTyp = 0;          /* type number of rate model                        */
int MixTyp = 0;         /* set for mixture models                           */
int DisFlg = 0;         /* set for discrete time models                     */
int VDAlloc = 0;        /* 1 if model specification allocated               */
int XMFlg[5];           /* flags for model terms (a,b,c,d)                  */
int NDef = 0;           /* number of vector definitions                     */
short *VTerm;           /* model specification: term                        */
short *VOrg;            /* model specification: origin state                */
short *VDes;            /* model specification: destination state           */
short *VSNum;           /* model specification: spell number                */
short *VNVar;           /* model specification: number of variables         */
short **VNIVar;         /* pointer to variable numbers                      */

short *PIdx;            /* array with variables numbers of parameters       */
char *PIdxM;            /* array with corresponding model terms             */
short *PIdxPtr;         /* pointer to transition-specific parameters        */
short *TranTVar;        /* set for VTyp 5 variables                         */
short *VTNum;           /* max number of covariates in transition           */
int VTMax = 0;          /* max( VTNum[i] ).                                 */

double *POLWrk1;        /* work array                                       */
int POLWrk1A = 0;       /* if allocated                                     */
double *POLWrk2;        /* work array                                       */
int POLWrk2A = 0;       /* if allocated                                     */
double *BCBeta;         /* work array                                       */
int BCBetaA = 0;        /* if allocated                                     */
double *BCTG;           /* work array                                       */
int BCTGA = 0;          /* if allocated                                     */
double *BCGam;          /* work array                                       */
int BCGamA = 0;         /* if allocated                                     */
double *DRJKW;          /* work array                                       */
int DRJKWA = 0;         /* if allocated                                     */
double *DAJKW;          /* work array                                       */
int DAJKWA = 0;         /* if allocated                                     */
double *DBJKW;          /* work array                                       */
int DBJKWA = 0;         /* if allocated                                     */

/* ------------------------------------------------------------------------ */
/*  rate()          Transition rate models. Command in CmdBuf.              */
/*                  Return 0 if OK, otherwise -1.                           */

int rate(void)
{
    int err;

    err = -1;         
    if (check_cmd(0))
        return(-1);

    printf1("Transition rate models. Current memory: %d bytes.\n",MemReq);

    if (EDAvail == 0) {
        p_err(-15,1);
        return(0);
    }
    set_mldef();                /* set defaults for ML estimation */

    if (parm(CmdBuf + 4,3,1))   /* get parameters */
        goto RTFin;

    MixTyp = 0;
    if (PMMIX)
        MixTyp = 1;

    /* check option for derivatives with algorithms 7 and 8 */

    if (MINA == 7 || MINA == 8) {
        if (PMDOPT < 0 || PMDOPT > 2)
            PMDOPT = 2;
    }

    if (prn_rmod(PMRHSI))       /* check/print type of model */
        goto RTFin;

    if (RMTyp == MEXP1 || RMTyp == MEXP2) {        /* time periods */
        if (PMNTP == 0) {
            p_err(-17,1);
            goto RTFin;
        }
        if (fabs(PMTP[0]) > EPSI) {      /* must begin with zero */
            p_err(-33,1);
            goto RTFin;
        }
        if (prn_tper())
            goto RTFin;
    }
    if (alloc_vdef(1)) {        /* allocate memory for model terms */
        p_err(-2,1);
        goto RTFin;
    }
    if (PMNXA) {                /* get covariable definitions */
        if (get_rvar(CmdBuf))
            goto RTFin;
    }
    if (NDef == 0)  
        printf1("Model without covariates.\n\n");
       
    if (check_mod())            /* check model definition */
        goto RTFin;

    if (build_pvec())           /* create list of model parameters */
        goto RTFin;

    set_mlopt();                /* adjust ML options */         

    /* allocate additional memory */

    if (rt_alloc(1))
        goto RTFin;

    /* for transition rate model we currently have only the standard form
       of covariance matrix calculation. */

    CCTyp = 2;
    if (ml_init(1,0,0))             /* init ML estimation */
        goto RTFin;

    exp_null();         /* calculate likelihood of exp null model */
   
    printf1("Log-likelihood of exponential null model: %lg\n",FMax0);
    if (DScalFlg == 0)
        printf1("Changed scaling factor for log-likelihood: %lg\n",DScal);

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))   /* try to get starting values */
        goto RTFin;
       
    if (DSVFlg == 0)
        ml_rstart();                /* default starting values */

    prot_init(1,RMTyp,0);           /* init protocol file */

    prn_cwt();                      /* case weight information */
       
    if (PMNConS > 0) {              /* process constraints */
        newline();
        if (con_proc(CmdBuf))
            goto RTFin;
    }

    if (RMTyp == MCOX && ESortFlg == 0) {        /* sort data */
        if (sort_ed(1))  
            goto RTFin;
    }

    /* -------------------------------------------------------------------- */
    /*  call ffmin() for maximization of log likelihood.                    */

    if (ffmin(RMTyp,0,1,1)) /* minimization */
        goto RTFin;         /* insuff memory or fn error */
       
    prn_mlres(0);           /* print info about min algorithm */
    prn_rcoeff(0);          /* print parameter estimates */

    if (LConv >= 0) {       /* if convergence */

        prn_ml1res(0);      /* additional ML results */

        if (PMRRFlg)        /* print relative risks  */
            prn_rcoeff(1);

        if (PMResFDef)          /* generalized residuals */
            prn_resid(RMTyp);

        if (PMPRN > 0)          /* prate option */
            prate(CmdBuf,RMTyp);

        if (RMTyp == MCOX && PMNTP > 0)     /* goodness of fit for Cox models */
            test_cox();
    }
    err = 0;

RTFin:
    if (ESortFlg)
        sort_ed(0);     
    rt_alloc(0);
    con_free();
    ml_init(0,0,0);     
    alloc_vdef(0);      
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_rmod(n)     Print name of rate model n and set RMTyp. If not        */
/*                  available error message and RMTyp = 0.                  */
/*                  Return 0 if OK, -1 if error                             */

int prn_rmod(int n)
{
    printf1("\nModel: ");

    switch (n) {

        case  MCOX:     printf1("Cox (partial likelihood)\n"); 
                        if (PM1NV > 0)  
                            printf1("Stratified with %d groups.\n",PM1NV);
                        break;
        case  MEXP:     printf1("Exponential\n");         
                        break;
        case  MEXP1:    printf1("Exponential with time-periods\n"); 
                        break;
        case  MEXP2:    printf1("Exponential with period-specific effects\n");          
                        break;
        case  MPOL:     printf1("Polynomial (I) with degree %d\n",PMDEG);        
                        break;
        case  MPOL1:    printf1("Polynomial (II) with degree %d\n",PMDEG);         
                        if (PMDEG > 0) {
                            printf1("Numerical integration with method %d. Relative error: %lg (%lg).\n",
                                                    NINTMETH,NINTRERR,NINTAERR);
                        }
                        break;
        case  MWEI:     printf1("Weibull\n");
                        break;
        case  MLL:      printf1("Log-logistic (I)\n");
                        break;
        case  MLL2:     printf1("Log-logistic (II)\n");        
                        break;
        case  MLN:      printf1("Log-normal\n");        
                        break;
        case  MIG:      printf1("Inverse Gaussian\n");         
                        break;
        case  MSIC:     printf1("Sickle\n");         
                        break;
       
        case  MGAM:     printf1("Generalized Gamma (kgam=%lg).\n",PMKGam);          
                        break;
    
        case  MGM:      printf1("Gompertz-Makeham\n");           
                        break;
  
        case  DLR:      printf1("Discrete time logistic regression, degree %d.\n",PMDEG);          
                        break;
        case  CLL:      printf1("Discrete time complementary log-log, degree %d.\n",PMDEG);         
                        break;

        default:        printf1("%d not defined.\n\n",n);
                        return(-1);
    }
                 
    if (MixTyp == 1)  
        printf1("Gamma mixture.\n");
    RMTyp = n;
    newline();             
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

int alloc_vdef(int mode)
{
    register int i;
    int err = -1;

    if (mode == 0) {
        if (VDAlloc == 0)
            return(0);

        for (i = 0; i < NDef; ++i) {
            if (VNVar[i] > 0) {
                free((char *)VNIVar[i]);
                memrq(-VNVar[i],sizeof(short));
                VNVar[i] = 0;
            }
        }
        NDef = err = 0;
        XMFlg[1] = XMFlg[2] = XMFlg[3] = XMFlg[4] = 0;
        goto VDA10;
    }

    if (!(TranTVar = (short *)calloc(NTran1,sizeof(short))))   
        goto VDA00;           
    memrq(NTran1,sizeof(short));

    if (!(VTNum = (short *)calloc(NTran1,sizeof(short))))   
        goto VDA0;           
    memrq(NTran1,sizeof(short));

    if (!(VTerm = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA1;           
    memrq(MaxVP,sizeof(short));
              
    if (!(VOrg = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA2;           
    memrq(MaxVP,sizeof(short));
              
    if (!(VDes = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA3;           
    memrq(MaxVP,sizeof(short));
              
    if (!(VSNum = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA4;           
    memrq(MaxVP,sizeof(short));
              
    if (!(VNVar = (short *)calloc(MaxVP,sizeof(short))))   
        goto VDA5;           
    memrq(MaxVP,sizeof(short));
              
    if (!(VNIVar = (short **)calloc(MaxVP,sizeof(short *))))   
        goto VDA6;           
    memrq(MaxVP,sizeof(short *));
              
    if (!(PIdx = (short *)calloc(MaxP + 2,sizeof(short))))   
        goto VDA7;           
    memrq(MaxP + 2,sizeof(short));
              
    if (!(PIdxM = (char *)calloc(MaxP + 2,sizeof(char))))   
        goto VDA8;           
    memrq(MaxP + 2,sizeof(char));
              
    if (!(PIdxPtr = (short *)calloc(NTran1,sizeof(short))))   
        goto VDA9;           
    memrq(NTran1,sizeof(short));
              
    VDAlloc = 1;
    return(0);

VDA10:
    free((char *)PIdxPtr);
    memrq(-NTran1,sizeof(short));
VDA9:
    free((char *)PIdxM);
    memrq(-MaxP - 2,sizeof(char));
VDA8:
    free((char *)PIdx);
    memrq(-MaxP - 2,sizeof(short));
VDA7:
    free((char *)VNIVar);
    memrq(-MaxVP,sizeof(short *));
VDA6:
    free((char *)VNVar);
    memrq(-MaxVP,sizeof(short));
VDA5:
    free((char *)VSNum);
    memrq(-MaxVP,sizeof(short));
VDA4:
    free((char *)VDes);
    memrq(-MaxVP,sizeof(short));
VDA3:
    free((char *)VOrg);
    memrq(-MaxVP,sizeof(short));
VDA2:
    free((char *)VTerm);
    memrq(-MaxVP,sizeof(short));
VDA1:
    free((char *)VTNum);
    memrq(-NTran1,sizeof(short));
VDA0:
    free((char *)TranTVar);
    memrq(-NTran1,sizeof(short));
VDA00:
    VDAlloc = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_rvar(pcmd)  get covariates for transition rate model.               */
/*  ##              Return 0 if OK, otherwise -1.                           */

int get_rvar(char *pcmd)
{
    int err;
    register char c,*p,*q;

    err = 0;
    p = pcmd;
    while (*p) {
        if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) || 
            !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3)) { 

            q = skip_xa(p);
            c = *q;
            *q = '\0';
            err = save_vdef(p);
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

int save_vdef(char *s)
{
    register int i;
    int n,sn,org,des,nb;
    register char *p;

    if (NDef + 3 >= MaxVP) {
        printf1("Exceeded maximum number of vector definitions in model specification.\n");
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
    XMFlg[i] = 1;
    VTerm[NDef] = i;

    if (*++p != '(')
        goto SVD1;

    if (sscanf(p,"(%d,%d,%d)",&org,&des,&sn) == 3) ;
    else if (sscanf(p,"(%d,%d)",&org,&des) == 2)
        sn = 1;
    else
        goto SVD1;

    n = 0;
    for (i = 0; i < NTran1; ++i) {  /* check if sn,org defined */
        if (des == DesTran1[i] && org == OrgTran1[i] && sn == SnTran1[i]) {
            n = 1;
            break;
        }
    }
    if (n == 0) {
        printf1("Error: %s\nNo transition with (org,des,sn) = (%d,%d,%d).\n",s,org,des,sn);
        return(-1);
    }
    VOrg[NDef] = org;
    VDes[NDef] = des;
    VSNum[NDef] = sn;

    while (*p++ && *p != ')');
    if (*p++ != ')' || *p++ != '=')
        goto SVD1;

    if (*p == '1') {            /* only constant */
        VNVar[NDef++] = 0;
        return(0);
    }
    p = get_nvia(p,&n,0,&nb);
    if (n <= 0) {
        printf1("Error: %s\nDefinition contains undefined variables.\n",s);
        return(-1);
    }
    if (nb) {
        p_err(-42,1);
        return(-1);
    }
    VNVar[NDef] = n;

    if (!(VNIVar[NDef] = (short *)calloc(n,sizeof(short)))) { 
        p_err(-2,1);
        return(-1);
    }
    memrq(n,sizeof(short));
    for (i = 0; i < n; ++i)
        VNIVar[NDef][i] = VLVIdx[i]; 

    NDef++;
    alloc_vl(0);
    return(0);

SVD1:
    printf1("Syntax error: %s\n",s);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  check_mod       Check rate model specification.                         */
/*                  Return 0 if OK, -1 if error.                            */

int check_mod(void)
{
    switch (RMTyp) {
        case MCOX:  if (!XMFlg[1])
                        goto CMErr1;
                    if (XMFlg[2])
                        goto CMErr4;
                    if (XMFlg[3])
                        goto CMErr5;
                    break;

        case MEXP:  if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2])
                        goto CMErr4;
                    if (XMFlg[3]) 
                        goto CMErr5;
                    break;
    
        case MEXP1: 
        case MEXP2: if (PMNTP == 0)   
                        goto CMErr9; 
                    if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2]) 
                        goto CMErr4;
                    if (XMFlg[3])
                        goto CMErr5;
                    break;
   
        case MPOL:  if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2]) 
                        goto CMErr4;
                    if (XMFlg[3]) 
                        goto CMErr5;
                    if (PMDEG)
                        XMFlg[2] = 1;
                    else
                        XMFlg[2] = 0;
                    break;

        case MPOL1: if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2])
                        goto CMErr4;
                    if (XMFlg[3]) 
                        goto CMErr5;
                    if (PMDEG)
                        XMFlg[2] = 1;
                    else
                        XMFlg[2] = 0;
                    break;
    
        case MWEI:  
        case MLL:
        case MLL2:  
        case MLN:   
        case MIG:   
        case MSIC:  
        case MGAM:  if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (!XMFlg[2])
                        XMFlg[2] = 1;
                    if (RMTyp != MLL2 && RMTyp != MLN) {
                        if (XMFlg[3])
                            goto CMErr5;
                    }
                    else if (RMTyp != MLN && !XMFlg[3])
                        XMFlg[3] = 1;

                    if (RMTyp == MLN && !XMFlg[3] && PMNI != 0)
                        XMFlg[3] = 1;

                    break;
         
        case MGM:   if (!XMFlg[1] && !XMFlg[2] && !XMFlg[3])
                        XMFlg[2] = 1;
                    if (XMFlg[2] && !XMFlg[3])
                        XMFlg[3] = 1;
                    else if (!XMFlg[2] && XMFlg[3])
                        XMFlg[2] = 1;
                    /***
                    if (!XMFlg[1] && !XMFlg[2] && !XMFlg[3])
                        XMFlg[1] = 1;
                    ***/
                    if (!XMFlg[1] && PMNI != 0)
                        XMFlg[1] = 1;
                    break;
         
        case DLR:   if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2])
                        goto CMErr4;
                    if (XMFlg[3]) 
                        goto CMErr5;
                    if (PMDEG)
                        XMFlg[2] = 1;
                    else
                        XMFlg[2] = 0;
                    DisFlg = 1;
                    break;

        case CLL:   if (!XMFlg[1])
                        XMFlg[1] = 1;
                    if (XMFlg[2]) 
                        goto CMErr4;
                    if (XMFlg[3]) 
                        goto CMErr4;
                    if (PMDEG)
                        XMFlg[2] = 1;
                    else
                        XMFlg[2] = 0;
                    DisFlg = 1;
                    break;

        default:    return(-1);
    }
    if (MixTyp == 1) {
        if (DisFlg || RMTyp == MEXP1  || RMTyp == MPOL   || RMTyp == MPOL1 ||
                      RMTyp == MLL2   || RMTyp == MGAM   ||
                      RMTyp == MIG    || RMTyp == MGM    || RMTyp == MEXP2 ||
                      RMTyp == MCOX   || (RMTyp == MLN   && XMFlg[3])) {

            goto CMErr7;
        }
        if (!XMFlg[4])  
            XMFlg[4] = 1;
    }
    else if (XMFlg[4])  
        goto CMErr8;
         
    return(0);

CMErr1:     printf1("Error: model 1 (partial likelihood) requires covariates.\n");
            return(-1);
CMErr4:     printf1("Error: B term variables not possible.\n");
            return(-1);
CMErr5:     printf1("Error: C term variables not possible,\n");
            return(-1);
CMErr7:     printf1("Error: Gamma mixture not possible.\n");
            return(-1);
CMErr8:     printf1("Error: D term variables possible only with mixtures.\n");
            return(-1);
CMErr9:     printf1("Error: this model needs time periods.\n");
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

int build_pvec(void)
{
    register int ii,i,j,k,kk,l;
    int nv,sn,org,des,kkn,first,npx;

    VTMax = 0;
    kkn = 1;

    if (RMTyp == MEXP2)
        kkn = PMNTP;
      
    NParm = 0;          /* number of model parameters */

    for (ii = 0; ii < NTran1; ++ii) {       /* for all transitions */

        npx = NParm;

        TranTVar[ii] = 0;           /* set for VTyp 5 variables */
        sn = SnTran1[ii];
        org = OrgTran1[ii];
        des = DesTran1[ii];
        first = 1;
        nv = 0;

        for (k = 1; k <= MTerm; ++k) {      

            if (XMFlg[k]) {        

                if (first) {
                    PIdxPtr[ii] = NParm + 1;
                    first = 0;
                }   

                /* first insert a constant, skip Cox models */

                if (RMTyp != MCOX) {    
                    if (++NParm > MaxP) 
                        goto BPErr1;

                    PIdx[NParm] = 0;
                    PIdxM[NParm] = k;
                }

                /* add period-specific constants */
         
                if (RMTyp == MEXP1 || RMTyp == MEXP2) {
                    for (i = 1; i < PMNTP; ++i) {
                        if (++NParm > MaxP)  
                            goto BPErr1;

                        PIdx[NParm] = -i;
                        PIdxM[NParm] = k;
                    }
                }
                if (k == 2 && (RMTyp == MPOL || RMTyp == MPOL1 ||
                               RMTyp == DLR  || RMTyp == CLL)) {   

                    for (i = 1; i < PMDEG; ++i) {
                        if (++NParm > MaxP)  
                            goto BPErr1;

                        PIdx[NParm] = -i;
                        PIdxM[NParm] = k;
                    }
                }
     
                /* now insert indices for the user defined covariates */

                for (i = 0; i < NDef; ++i) {

                    if ((int)VTerm[i] == k && (int)VOrg[i] == org 
                                           && (int)VDes[i] == des 
                                           && (int)VSNum[i] == sn) {

                        for (j = 0; j < VNVar[i]; ++j) {
                    
                            /* repeat only for MEXP2 model */           
  
                            for (kk = 1; kk <= kkn; ++kk) {

                                if (++NParm > MaxP)  
                                    goto BPErr1;
     
                                l = VNIVar[i][j];
                                nv++;

                                if (RMTyp == MCOX && !j)
                                    PIdx[NParm] = -l - 1;
                                else
                                    PIdx[NParm] = l + 1;

                                PIdxM[NParm] = k;
                                if (VTyp[l] == 5)  /* VTyp 5 variable */
                                    TranTVar[ii] = 1;
                            }   
                        }
                    }
                }
            }
        }
        VTNum[ii] = nv;
        if (VTMax < nv)
            VTMax = nv;

        if (NParm == npx)           
            goto BPErr2;

    }
    if (RMTyp == MCOX)
        PIdx[NParm + 1] = -1;
    else
        PIdx[NParm + 1] = 0;

    PIdxM[NParm + 1] = 0;
    if (NParm <= 0) {
        printf1("Error: number of model parameters is %d.\n",NParm);
        return(-1);
    }

    /***********************************************************
    printf1("\n  PIdx : ");
    for (i = 1; i <= NParm + 1 ; ++i) printf1(" %2d",PIdx[i]);
    printf1("\n  PIdxM: ");
    for (i = 1; i <= NParm + 1; ++i) printf1(" %2d",PIdxM[i]);
    printf1("\nPIdxPtr\n");
    for (i = 0; i < NTran1; ++i)
        printf1("%2d ",PIdxPtr[i]);
    printf1("\n");
    ***********************************************************/
    return(0);

BPErr1: 
    printf1("Error: exceeded max number of parameters (%d)\n",MaxP);
    return(-1);
BPErr2: 
    printf1("Error: need at least one parameter (covariate) for each transition.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  exp_null.   Calculate FMax0, the log-likelihood of an exponential       */
/*              null model. Also calculate the scaling factor DScal.        */
/*              DScal is only calculated if DScalFlg = 0.                   */

void exp_null(void)
{
    register int i;
    double w;

    FMax0 = 0.0;
    for (i = 0; i < NTran1; ++i) {
        if (TWFreq[i] > 0.0 && TWDur[i] > 0.0)  
            FMax0 += (double)TWFreq[i] * (rlog((double)TWFreq[i]/TWDur[i]) - 1.0);
    }
    if (DScalFlg == 0) {
        DScal = -1.0;                                                     
        w = fabs(FMax0);
        while (w < 0.5) {
            w *= 10.0;
            DScal *= 10.0;
        }
        while (w > 5.0) {
            w /= 10.0;
            DScal /= 10.0;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ml_rstart.      Calculation of default starting values.                 */

void ml_rstart(void) 
{
    register int i,j,k;
    double tmp,freq,dur,rate,lrate;

    printf1("Using default starting values.\n");

    for (j = 1; j <= NParm; ++j)
        Par[j] = 0.0;

    if (RMTyp == MCOX)        /* Cox model has zero starting values */
        return;
            
    j = 1;
    for (i = 0; i < NTran1; ++i) {

        freq = TWFreq[i];
        dur  = TWDur[i];
        if (freq > 0.0 && dur > 0.0) {
            rate = freq / dur; 
            lrate = rlog(rate); 
        }
        else
            rate = lrate = 0.0;

        for (k = 1; k <= MTerm; ++k) {

            if (XMFlg[k]) {

                tmp = 0.0;
                if (k == 1) {
                    if (RMTyp == MLN || RMTyp == MGAM)
                        tmp = -lrate;
                    else if (RMTyp == MSIC)
                        tmp = 2.0 * lrate;
                    else if (RMTyp == MIG)
                        tmp = rate;
                    else if (RMTyp != MLL2)
                        tmp = lrate;
                }
                else if (k == 2) {
                    if (RMTyp == MGM)
                        tmp = lrate;
                    else if (RMTyp == MSIC)
                        tmp = -lrate;
                    else if (RMTyp == MIG)
                        tmp = rlog(sqrt(rate));
                }
                else if (k == 3) {
                    if (RMTyp == MGM) 
                        tmp = -rate;
                    else if (RMTyp == MLL2)
                        tmp = lrate;
                }
                Par[j] = tmp;
                if (RMTyp == MEXP1 || RMTyp == MEXP2) {
                    while (PIdx[++j] < 0)  
                        Par[j] = tmp;
                    j--;
                }
                while (PIdx[++j])  
                    ;
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_rcoeff(mode)                                                        */
/*  ##  Print estimated parameters, standard errors, and test statistics    */
/*      of rate model estimation. Print format is PMTFmt.                   */
/*      If mode == 0 print parameter estimates.                             */  
/*      If mode == 1 print relative risks.                                  */

void prn_rcoeff(int mode)
{
    register int i,j,k,sn,idx;             
    int lflag,tp,len,m,org,des;

    if (mode == 1) {
        printf1("Estimated relative risks");
        switch (RMTyp) {
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
            default:    printf1(": not available for this model.\n\n");
                        return;
        }
        printf1(".\n\n");
    }
    len = 0;
    printf1("Idx SN Org Des MT ");          
    if (RMTyp == MEXP2) {
        printf1(" P ");         
        len = 3;
    }
    printf1("Variable ");        
    prnchar(' ',VNameLen - 8,0);

    lflag = 0;      /* check for var labels */

    for (i = 1; i <= NParm; ++i) {
        j = PIdx[i];   
        if (j > 0 || (j < 0 && RMTyp == MCOX)) {
            if (j < 0)
                j = -j;
            if (VLabel[j - 1] != NULL) {
                lflag = 1;
                break;
            }
        }
    }
    if (lflag) {
        printf1("Label ");
        prnchar(' ',VLabelLen - 5,0);
    }
    if (mode == 1) {
        prnchar(' ',PMTFmt1 - 10,0); printf1("    R.Risk");        
        len += 18;
        m = 1;
    }
    else {
        prnchar(' ',PMTFmt1 - 10,0); printf1("     Coeff"); 
        prnchar(' ',PMTFmt1 -  9,0); printf1("     Error"); 
        prnchar(' ',PMTFmt1 -  9,0); printf1("   C/Error  Signif");
        len += 26;
        m = 3;
    }
    if (lflag)
        len += VLabelLen + 1;

    k = 1;

    for (i = 0; i < NTran1; ++i) {

        sn = SnTran1[i];
        org = OrgTran1[i];
        des = DesTran1[i];

        printf1("\n");           
        prnchar('-',len + m * (PMTFmt1 + 1) + VNameLen,0);

        for (j = 1; j <= MTerm; ++j) {

            if (XMFlg[j]) {

                tp = 1;  /* time period */

                while (1) {

                    /* relative risks only for selected model terms */

                    if (mode == 1) {
                        if (RMTyp == MPOL1 || RMTyp == MWEI || RMTyp == MSIC) {
                            if (j != 1)
                                goto PRCNXT;
                        }
                        else if (RMTyp == MGM) {
                            if (j != 2)
                                goto PRCNXT;
                        }
                        else if (RMTyp == MLL2) {
                            if (j != 3)
                                goto PRCNXT;
                        }
                    }
                    idx = PIdx[k];
                    if (idx || RMTyp != MCOX) {
                        printf1("\n%3d %2d %3d %3d  %c ",k,sn,org,des,(char)('A'+j-1));

                        if (idx > 0 || RMTyp == MCOX) {

                            if (RMTyp == MEXP2)  
                                printf1("%2d ",tp);

                            if (idx < 0 && RMTyp == MCOX)
                                idx = -idx;
                            prn_vname(idx - 1);
                            if (lflag)
                                prn_vlabel(idx - 1);
                        }
                        else {
                            if (RMTyp == MEXP2)  
                                printf1("%2d ",tp);
                                                 
                            if (RMTyp == MEXP1 || RMTyp == MEXP2) {
                                printf1("Period-%-2d",-idx + 1);          
                                prnchar(' ',VNameLen - 8,0);
                            }
                            else if (j == 2 &&
                                (RMTyp == MPOL || RMTyp == MPOL1 ||
                                 RMTyp == DLR  || RMTyp == CLL)) {  

                                printf1("Beta-%-2d",-idx + 1);         
                                prnchar(' ',VNameLen - 6,0);
                            }
                            else {
                                printf1("Constant");          
                                prnchar(' ',VNameLen - 7,0);
                            }
                            if (lflag)
                                prnchar(' ',VLabelLen + 1,0);
                        }
                        prn_rc(k,mode);  /* print coefficients */
                    }
PRCNXT:           
                    k++;
                    if (!PIdx[k] || (RMTyp == MCOX && PIdx[k] < 0))
                        break;
       
                    if (RMTyp == MEXP2) {
                        if (++tp > PMNTP)
                            tp = 1;
                    }
                }
            }
        }
    }
    printf1("\n\n");
}

/* ------------------------------------------------------------------------ */
/*  prn_rc                                                                  */
/*           Print coefficients, used by prn_rcoeff.                        */
/*           If mode == 1 print relative risks                              */

void prn_rc(int k,int mode)
{
    register int j;
    double tmp,tmp1,tmp2;

    if (mode == 1) {
        tmp = 1.0;
        if (RMTyp == MWEI) {
            j = k;
            while (PIdx[++j]) ;
            tmp = rexp(Par[j]);
        }
        printf1(PMTFmtS,rexp(tmp * Par[k]));
        return;
    }
    else  
        printf1(PMTFmtS,Par[k]);
       
    if (LConv >= 0 && Diag[k] > 0.0) {
        tmp = sqrt(Diag[k]);
        printf1(PMTFmtS,tmp);        
    }
    else {
        tmp = 0.0;
        prnchar(' ',PMTFmt1 - 3,0);
        printf1("--- ");
    }
    if (LConv >= 0 && tmp > 0.0) {
        tmp1 = Par[k] / tmp;
        printf1(PMTFmtS,tmp1);
        tmp2 = 2.0 * cdnf(fabs(tmp1)) - 1.0;
        printf1(" %6.4lf",tmp2);
    }
    else {
        prnchar(' ',PMTFmt1 - 3,0);
        printf1("---     ---");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_resid(mod)                                                          */
/*      Print residuals for model mod to output file PMResFd.               */

void prn_resid(int mod)
{
    if (NTran1 > 1 || MEFlg || MixTyp || mod == MCOX || mod == DLR || mod == CLL) {    
        printf1("Cannot calculate generalized residuals; pres option ignored.\n");
        return;
    }
    printf1("Calculating generalized residuals; will be written to: %s\n",PMResFName);

    fprintf(PMResFd,"# Generalized residuals.\n");
    fprintf(PMResFd,"# ");
    fprnchar(PMResFd,' ',16 + 4 * PMMFmt1 - 8,0); fprintf(PMResFd,"Survivor\n# ");
    fprintf(PMResFd,"Case Org Des ");
    fprnchar(PMResFd,' ',PMMFmt1 - 2,0); fprintf(PMResFd,"TS");
    fprnchar(PMResFd,' ',PMMFmt1 - 1,0); fprintf(PMResFd,"TF");
    fprnchar(PMResFd,' ',PMMFmt1 - 3,0); fprintf(PMResFd,"Rate");
    fprnchar(PMResFd,' ',PMMFmt1 - 7,0); fprintf(PMResFd,"Function ");
    fprnchar(PMResFd,' ',PMMFmt1 - 8,0); fprintf(PMResFd,"Residual ");
    fprnchar(PMResFd,' ',PMMFmt1 - 6,0); fprintf(PMResFd,"Weight\n# ");
    fprnchar(PMResFd,'-',12 + 6 * (PMMFmt1 + 1),1); 

    LFunc = LGrad = LSec = 0;

    switch (mod) {

        case  MEXP:  fn_exp(1);     /* exponential model */
                     break;
    
        case  MEXP1: fn_exp1(1);    /* exp. model with time periods */
                     break;
     
        case  MEXP2: fn_exp2(1);    /* exp. model with time periods, II */
                     break;     
     
        case  MPOL:  fn_pol(1);     /* polynomial, I */
                     break;     

        case  MPOL1: fn_pol1(1);    /* polynomial, II */
                     break;     
       
        case  MGM:   fn_gm(1);      /* Gompertz model */
                     break;

        case  MWEI:  fn_wei(1);     /* Weibull model */
                     break;

        case  MSIC:  fn_sic(1);     /* sickle model */
                     break;
   
        case  MLL:   fn_ll(1);      /* log-logistic model, I */
                     break;

        case  MLL2:  fn_ll2(1);     /* log-logistic model, II */     
                     break;
     
        case  MLN:   fn_ln(1);      /* log-normal model */
                     break;

        case  MIG:   fn_ig(1);      /* inverse gaussian model */
                     break;

        case  MGAM:  fn_gam(1);     /* gamma */
                     break;

        default:
        fprintf(PMResFd,"  Calculation in current version of TDA not possible.\n");
    }
    newline();
}

/* ------------------------------------------------------------------------ */
/*  prn_resid1                                                              */
/*      Print residuals to output file PMResFd.                            */

void prn_resid1(int icase, int org, int des, double ts, double tf,
    double rate, double surv, double resid, double wt)
{
    fprintf(PMResFd,"%6d ",icase);
    fprintf(PMResFd,"%3d %3d ",org,des);
    fprintf(PMResFd,PMMFmtS,ts);
    fprintf(PMResFd,PMMFmtS,tf);
    fprintf(PMResFd,PMMFmtS,rate);
    fprintf(PMResFd,PMMFmtS,surv);
    fprintf(PMResFd,PMMFmtS,resid);
    fprintf(PMResFd,PMMFmtS,wt);
    fprintf(PMResFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prn_tper()      Print table with time periods and number of ending      */
/*                  times and events.                                       */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int prn_tper(void)
{
    register int i;
    int err,icase,sn,org,des,spl,nspl;
    double ts,tf;

    err = -1;

    if (alloc_acn(PMNTP))
        goto TPERFin;

    if (alloc_aci(PMNTP))
        goto TPERFin;

    if (alloc_acj(PMNTP))
        goto TPERFin;
     
    get_spell(1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {
        for (i = PMNTP - 1; i >= 0; --i) {
            if (ts >= PMTP[i]) {
                AcN[i] += 1;
                break;
            }
        }
        for (i = PMNTP - 1; i >= 0; --i) {
            if (tf >= PMTP[i]) {
                AcI[i] += 1;
                if (des != org)
                    AcJ[i] += 1;
                break;
            }
        }
    }
    printf1("    Time period        Starting times  Ending times    Events\n  ");         
    prnchar('-',59,1);

    for (i = 0; i < PMNTP; ++i) {
        printf1("%9.4lf - ",PMTP[i]);           
        if (i < PMNTP - 1)
            printf1("%9.4lf",PMTP[i + 1]);
        else
            printf1("         ");
        printf1("  %14d  %12d %9d\n",AcN[i],AcI[i],AcJ[i]);
    }
    newline();
    err = 0;

TPERFin:
    alloc_acn(0);       
    alloc_aci(0);       
    alloc_acj(0);       
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rt_alloc(opt)   if opt != 0 allocate additional work arrays for rate    */
/*                  models, otherwise free previously allocated memory.     */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int rt_alloc(int opt)
{
    int err = 0;

    if (opt == 0)
        goto RTAFin;

    err = -1;

    if (RMTyp == MPOL || RMTyp == MPOL1) {

        if (!(POLWrk1 = (double *)calloc(2 * PMDEG + 1,sizeof(double))))   
            goto RTAFin;         
        POLWrk1A = 2 * PMDEG + 1;
        memrq(POLWrk1A,sizeof(double));

        if (!(POLWrk2 = (double *)calloc(PMDEG + 1,sizeof(double))))   
            goto RTAFin;         
        POLWrk2A = PMDEG + 1;
        memrq(POLWrk2A,sizeof(double));

        if (!(BCBeta = (double *)calloc(PMDEG + 1,sizeof(double))))   
            goto RTAFin;         
        BCBetaA = PMDEG + 1;
        memrq(BCBetaA,sizeof(double));

        if (!(BCTG = (double *)calloc(PMDEG + 1,sizeof(double))))   
            goto RTAFin;         
        BCTGA = PMDEG + 1;
        memrq(BCTGA,sizeof(double));

    }
    if (RMTyp == DLR || RMTyp == CLL) {

        if (!(DRJKW = (double *)calloc(NTran1 + 1,sizeof(double))))   
            goto RTAFin;         
        DRJKWA = NTran1 + 1;
        memrq(DRJKWA,sizeof(double));
    }
    if (RMTyp == CLL) {

        if (!(DAJKW = (double *)calloc(NTran1 + 1,sizeof(double))))   
            goto RTAFin;         
        DAJKWA = NTran1 + 1;
        memrq(DAJKWA,sizeof(double));

        if (!(DBJKW = (double *)calloc(NTran1 + 1,sizeof(double))))   
            goto RTAFin;         
        DBJKWA = NTran1 + 1;
        memrq(DBJKWA,sizeof(double));
    }

    return(0);

RTAFin:
    if (POLWrk1A) {
        free((char *)POLWrk1);
        memrq(-POLWrk1A,sizeof(double));
        POLWrk1A = 0;               
    }
    if (POLWrk2A) {
        free((char *)POLWrk2);
        memrq(-POLWrk2A,sizeof(double));
        POLWrk2A = 0;               
    }
    if (BCBetaA) {
        free((char *)BCBeta);
        memrq(-BCBetaA,sizeof(double));
        BCBetaA = 0;               
    }
    if (BCTGA) {
        free((char *)BCTG);
        memrq(-BCTGA,sizeof(double));
        BCTGA = 0;               
    }
    if (BCGamA) {
        free((char *)BCGam);
        memrq(-BCGamA,sizeof(double));
        BCGamA = 0;               
    }
    if (DRJKWA) {
        free((char *)DRJKW);
        memrq(-DRJKWA,sizeof(double));
        DRJKWA = 0;               
    }
    if (DAJKWA) {
        free((char *)DAJKW);
        memrq(-DAJKWA,sizeof(double));
        DAJKWA = 0;               
    }
    if (DBJKWA) {
        free((char *)DBJKW);
        memrq(-DBJKWA,sizeof(double));
        DBJKWA = 0;               
    }
    return(err);
} 




