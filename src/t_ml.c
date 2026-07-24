/****************************************************************************/
/*  t_ml                                                                    */
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
#include "t_gf.h"
#include "t_fnrta.h"
#include "t_fnrtb.h"
#include "t_fnrtc.h"
#include "t_fnqr.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_con.h"
#include "t_tmin.h"
#include "t_edat.h"
#include "t_rate.h"
#include "t_gmin.h"
#include "t_lin.h"
#include "t_eval.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_mds.h"
#include "t_cl.h"
#include "t_ireg.h"

/*  functions in t_ml.c */

void set_mldef(void);
void set_mlopt(void);
int get_dsv(int n,double *par,int gmina,double *lb,double *ub,int nl);
int ml_init(int opt,int typ,int gmina);
int syminv1(int n, double *d, double *h);
int syminv2(int n, double *x);
void fn(double *x, int mod, int typ, int scal, int opt,int *err);
int fn_umod(void);
void prjvec(double *x);
void prjhess(void);
int checkov(void);
void prn_mlres(int typ);
void prn_ml1res(int typ);
void prvec(char *txt, int n, double *x);
void prmat(char *txt,int n,int m,double *x);
void prhess(char *txt, int n);
void prval(char *txt, double x);
void prot_init(int typ,int mod,int gmina);
void prot_prni(char *s,int n,int *x);
void prot_prns(char *s,int n,short *x);
void prot_prnc(char *s,int n,char *x);

/*--------------------------------------------------------------------------*/
/*  global variables defined in t_ml.c                                      */

int MLAlloc = 0;            /* set if arrays for ML estimation allocated    */

double  OFMax =  3.3333333e33;  /* Values used for overflow and underflows  */
double  OFMin = -3.3333333e33;  /* checks in checkov().                     */

double FNTStat = -1.0;      /* Test statistic calculated by fn()            */
int FNTRank = 0;            /* Rank of hessian with FNTStat calculation     */
int GCTestFlg = 0;          /* Set if gctest command                        */
double GCTest = -1.0;       /* GC test statistic                            */
int GCRank = 0;             /* Rank of hessian with GC calculation          */
int LMTestFlg = 0;          /* Set if lmtest command                        */
double LMTest = -1.0;       /* LM test statistic                            */
int LMRank = 0;             /* Rank of hessian with LM calculation          */
int NWTest = 0;             /* Number of wtest commands                     */
short **WTest;              /* Storage for WTest[i][j]                      */
short *WTestN;              /* Number of entries in WTest                   */
/*--------------------------------------------------------------------------*/
int MINA = 5;               /* Type of minimization algorithm               */
int CCTyp = 2;              /* Type of covariance matrix calculation        */
int CCovFlg = 0;            /* Set during cov matrix calculation            */
int MCovFlg = 0;            /* Set during covariance matrix calculation     */
int CGradFlg = 0;           /* Set when calculating gradients into matrix   */

int ParTyp = 1;             /* Type of model parameterization               */
int LConv = 0;              /* Return Code of minimization algorithms       */
int NOCUsed = 0;            /* number of cases used                         */
int NumF = 0;               /* Number of function calls                     */
int NumFG = 0;              /* Number of function calls with gradient       */
int NumFH = 0;              /* Number of function calls with hessian        */

int NParm = 0;              /* Total number of model parameter              */
int NParm1 = 0;             /* Same in reduced space                        */
int NParmA = 0;             /* Used to alloc Par[]                          */
int HSiz = 0;               /* Size of lower hessian triangle               */
int HSiz1 = 0;              /* Same in reduced space                        */
int HSizA = 0;              /* Used to alloc Hess[]                         */
/*--------------------------------------------------------------------------*/
int Iter = 0;               /* Counter for iterations                       */
int MxIter = -1;            /* Max number of iterations                     */
int MxItFlg = 0;            /* Set if mxit command used                     */
int MxIt1  = 50;            /* Max number of iterations for linear search   */
int MxItR  = 100;           /* Max number of iterations for random search   */

int Crite = 1;              /* Convergence Criterion (1 - 6)                */
int CritFlg = 0;            /* Set if user defined criterion                */
double TOLG = 1.e-6;        /* Tolerance for gradient                       */
double TOLF = 1.e-12;       /* Tolerance for function value changes         */
double TOLP = 1.e-4;        /* Tolerance for parameter changes              */
double TOLV = 1.e-10;       /* Tolerance for Simplex algorithm              */
double TOLS = 1.e-4;        /* Tolerance for direct search algorithm        */
double TOLI = 1.e-3;        /* Tolerance for Monte Carlo integration        */
double TOLSG = 1.e-5;       /* Tolerance for scaled gradient                */
double TOLSP = 1.e-8;       /* Tolerance for scaled parameter change        */
double TOLSGF = 0.0;        /* Final scaled gradient                        */
double TOLSPF = 0.0;        /* Final scaled parameter change                */
double TOLA = 1.e-20;       /* A-convergence tolerance, used for rsearch    */

int CCheck = 5;             /* Check convergence interval in Simplex alg.   */
int STFlg = 0;              /* Set if function evaluation without           */
                            /* derivatives in step size search.             */

double AMue = 0.2;          /* Mue for Armijo conditions in fmin            */
double SMin = 0.0;          /* Minimal step size                            */
double SLen = 1.0;          /* Step length in direct search algorithm       */
double SRed = 0.5;          /* Step reduction in direct search algorithm    */
double SRLen = 1.0;         /* Step length for random search algorithm      */
double SRRed = 0.1;         /* Step length reduction for random search      */

double CValG = 0.0;         /* Value for convergence criterion 1            */
double CValF = 0.0;         /* Value for convergence criterion 2            */
double CValP = 0.0;         /* Value for convergence criterion 3            */
double CValV = 0.0;         /* Value for convergence criterion 4            */

int SVFlg = 0;              /* Set if singular value calculation for final  */
                            /* hessian requested                            */
/*--------------------------------------------------------------------------*/
int LLTN = 0;               /* Set for log likelihood sensitivity test      */
double LLTA = 0.01;
/*--------------------------------------------------------------------------*/
double FMax0 = 0.0;         /* Max log-likelihood (exponential null model)  */
double FMax  = 0.0;         /* Max log-likelihood (estimated model)         */
double FMax1 = 0.0;         /* Log-likelihood (at starting values)          */
double FMin  = 0.0;         /* Value of function after minimization         */
double FTmp  = 0.0;         /* Same, but temporary                          */

double *Par;                /* Parameter vector (NParm)                     */
int ParA = 0;               /* if allocated                                 */
double *Par1;               /* Parameter vector (NParm)                     */
int Par1A = 0;              /* if allocated                                 */
double *ParLB;              /* lower bound                                  */
int ParLBA = 0;
double *ParUB;              /* upper bound                                  */
int ParUBA = 0;

double *ParS;               /* Parameter vector with starting values        */
int ParSA = 0;              /* if allocated                                 */
double *ParC;               /* Parameter vector used with constraints       */
int ParCA = 0;              /* if allocated                                 */
double *TPar;               /* Temporary pointer for Par                    */

double *Grad;               /* Gradient vector                              */
int GradA = 0;              /* if allocated                                 */
double *GTmp;               /* Temporary for outer product                  */
int GTmpA = 0;              /* if allocated                                 */
double *Srch;               /* Search vektor                                */
int SrchA = 0;              /* if allocated                                 */
double *Diag;               /* Diagonal of Hessian                          */
int DiagA = 0;              /* if allocated                                 */
double *Diag1;              /* Diagonal of Hessian                          */
int Diag1A = 0;             /* if allocated                                 */
double *WrkD;               /* Working area with NParm dimension            */
                            /* used to handle constraints                   */
int WrkDA = 0;              /* if allocated                                 */
double *WrkE;               /* Working area with NParm dimension            */
int WrkEA = 0;              /* if allocated                                 */
double *WrkH;               /* Working area with HSiz dimension             */
                            /* used to handle constraints!                  */
int WrkHA = 0;              /* if allocated                                 */
double *Hess;               /* Strict lower triangle of Hessian (HSiz)      */
int HessA = 0;              /* if allocated                                 */
double *Hess1;              /* Strict lower triangle of Hessian (HSiz)      */
int Hess1A = 0;             /* if allocated                                 */
double *TSHess;             /* Full Hessian for MINA 7 and 8.               */
int TSHessA = 0;            /* if allocated                                 */
double *TSGrad;             /* Gradient for MINA 7 and 8.                   */
int TSGradA = 0;            /* if allocated                                 */

double DScal = 1.0;         /* Scaling factor                               */
int DScalFlg = 0;           /* Set if user defined scaling factor           */

int LFunc  = 0;             /* Flag: function value calculation             */
int LGrad  = 0;             /* Flag: calculation of gradient                */
int LSec   = 0;             /* Flag: calculation of Hessian                 */
int FN_First = 0;           /* set for first call of fn()                   */
int DSVFlg = 0;             /* set if externally supplied starting values   */

/*--------------------------------------------------------------------------*/
/*  set_mldef()     set default options for ML estimation.                  */
         
void set_mldef(void)
{
    MINA = 5;          /* type of min algorithm */
    MxIter = -1;        /* max number of iterations */
    MxItFlg = 0;        /* set if mxit command used */
    MxIt1 = -1;         /* max number of subiterations */
    MxItR = 100;        /* max number of random search iterations */
    Crite = 1;          /* type of convergence criterion */
    CritFlg = 0;        /* Set if user defined criterion */
    DScal = -1.0;       /* Default scaling factor */
    DScalFlg = 0;       /* Set to 1 if dscal command */
    TOLG = 1.e-6;       /* Tolerance for gradient */
    TOLF = 1.e-12;      /* Tolerance for function value changes */
    TOLP = 1.e-4;       /* Tolerance for parameter changes */
    TOLV = 1.e-10;      /* Tolerance for Simplex algorithm */
    TOLS = 1.e-4;       /* Tolerance for direct search algorithm */
    TOLSG = 1.e-5;      /* Tolerance for scaled gradient */
    TOLSP = 1.e-8;      /* Tolerance for scaled parameter change */

    AMue = 0.2;         /* Mue for Armijo conditions in fmin */
    SMin = 1.e-10;      /* Minimal step size */
    SLen = 1.0;         /* Step length in direct search algorithm */
    SRed = 0.5;         /* Step reduction in direct search algorithm */
    CCheck = 5;         /* Check convergence interval in Simplex alg. */
    STFlg = 0;          /* Set if function evaluation without derivatives */
    CCTyp = 2;          /* Type of covariance matrix calculation */
    SRLen = 1.0;        /* Step length for random search algorithm */
    SRRed = 0.1;        /* Step length reduction for random search */
}

/*--------------------------------------------------------------------------*/
/*  set_mlopt()     adjust ML options after the interpreation in parm().    */  
         
void set_mlopt(void)
{
    if (MINA == 7 || MINA == 8) {

        if (PMDOPT < 0 || PMDOPT > 2)
            PMDOPT = 2;

        TSDERIV = PMDOPT;
        if (MxIter < 0) {
            if (TSDERIV == 0)
                MxIter = 200;
            else if (TSDERIV == 1)
                MxIter = 100;
            else
                MxIter = 20;
        }
        Crite = 6;
    }
    if (MINA == 1)
        Crite = 5;
    else if (MINA == 2) 
        Crite = 4;
    else if (MINA == 3 || MINA == 4)
        Crite = 1;

    if (MxIter < 0) {
        if (MINA <= 4)
            MxIter = 100;
        else
            MxIter = 20;
    }
    if (MxIt1 < 0)
        MxIt1 = 50;
}

/* ------------------------------------------------------------------------ */
/*  get_dsv(n,par,gmina,lb,ub,nl)                                           */
/*                                                                          */
/*              Get starting values into par[i], i = 1,...,n.               */
/*              First check xp=... , then input file PMDSVName.             */
/*              If gmina != 0 put lower and upper bounds into lb[],ub[].    */
/*              Default is 0 +/- 1. If nl != 0 additional newline.          */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */
/*              If successful set DSVFlg = 1.                               */
/*                                                                          */

int get_dsv(int n,double *par,int gmina,double *lb,double *ub,int nl) 
{
    FILE *fd;
    register int i;
    int err,m,fflag;
    char buf[202],*p;
    double tmp;

    DSVFlg = 0;
    err = -1;
    fflag = m = 0;
    if (PMNX > 0) {
        m = imin(n,PMNX);
        for (i = 1; i <= m; ++i) {
            par[i] = PMXX[i - 1];
            if (gmina) {
                lb[i] = PMXA[i - 1];
                ub[i] = PMXB[i - 1];
            }
        }
        printf1("%d starting value(s) from xp parameter.\n",m);
        DSVFlg = 1;
    }
    else if (PMDSVFlg) {
                     
        printf1("Reading starting values from: %s\n",PMDSVName);
        if (!(fd = fopen(PMDSVName,OPEN_RD))) {   
            printf1("Error: can't open %s\n",PMDSVName);
            goto GETDSVFin;
        }
        fflag = 1;

        while (fgets(buf,200,fd)) {
               
            if (check_drec(buf)) {  /* if a data record */
                p = skip_b(buf);    /* skip blanks */
                if (sscanf(p,"%lg",&tmp) != 1)  
                    goto GETDSVFin;
                par[m + 1] = tmp;
                if (gmina) {
                    p = skip_dbl(p);
                    p = skip_b(p);
                    if (sscanf(p,"%lg",&tmp) == 1) {
                        lb[m + 1] = tmp;
                        p = skip_dbl(p);
                        p = skip_b(p);
                    }
                    else
                        lb[m + 1] = par[m + 1] - 1.0;

                    if (sscanf(p,"%lg",&tmp) == 1)  
                        ub[m + 1] = tmp;
                    else
                        ub[m + 1] = par[m + 1] + 1.0;
                }
                m++;
            }
            if (m >= n)
                break;
        }
        if (m > 0)  
            DSVFlg = 1;
    }
    if (DSVFlg) {
        if (m < n)  
            printf1("Warning: found %d, would need %d starting values.\n",m,n);
        if (nl)
            newline();
    }
    for (i = m + 1; i <= n; ++i) {
        par[i] = 0.0;       
        if (gmina) {
            lb[i] = -1.0;
            ub[i] =  1.0;
        }
    }
    err = 0;   

GETDSVFin:
    if (fflag)
        fclose(fd);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ml_init(opt,typ,gmina)                                                  */  
/*                                                                          */
/*  Initialization for function minimization.                               */
/*                                                                          */
/*  typ = 0 : ML                                                            */
/*        1 : general func minimization                                     */
/*        2 : nonlinear regression                                          */
/*                                                                          */
/*  If opt != 0 allocate required memory for model estimation with          */
/*  NParm parameters. If opt == 0 free previously allocated memory.         */
/*  If gmina != 0 allocate additional memory for parameter bounds.          */
/*                                                                          */
/*  Return 0 if OK, -1 if error (insufficient memory).                      */

int ml_init(int opt,int typ,int gmina)
{
    int err = 0;

    DSVFlg = 0;             /* flag for externally supplied starting values */
    clear_npflg();          /* clear numerical problem flags */
    NINTMUsed = -1;

    if (opt == 0) {
        if (MLAlloc == 0)
            return(0);
        goto MLIFin;
    }
    err = -1;

    if (typ == 0)
        printf1("Maximum likelihood estimation.");
    else if (typ == 1)
        printf1("Function minimization.");
    else
        printf1("Nonlinear regression.");

    printf1("\nAlgorithm %d: ",MINA);

    switch (MINA) {
        case  1:    printf1("Direct Search\n\n");
                    break;
        case  2:    printf1("Simplex\n\n");
                    break;
        case  3:    printf1("conjugate gradients\n\n");
                    break;
        case  4:    printf1("BFGS\n\n");
                    break;
        case  6:    printf1("Newton (II)\n\n");
                    break;
        case  7:    printf1("CES (quadratic model [dopt=%d])\n\n",TSDERIV);
                    break;
        case  8:    printf1("CES (tensor model [dopt=%d])\n\n",TSDERIV);
                    break;
        default:    printf1("Newton (I)\n\n"); /* MINA 5 is default */
                    break;
    }
    NParmA = NParm;
    HSiz = NParm * (NParm - 1) / 2;     /* size of hessian lower triangular */
    NParm1 = NParm;                     /* dimensions in reduced space */
    HSiz1 = HSiz;

    HSizA = HSiz;
    if (MINA == 3) {                  /* Hess is used as a working area */
        if (HSizA < 5 * NParm + 2)     /* in ffmin5(). */
            HSizA = 5 * NParm + 2;
    }
    else if (MINA == 4) {
        if (HSizA < NParm * (NParm + 7) / 2)
            HSizA = NParm * (NParm + 7) / 2;
    }

    /* print info about ML estimation */

    printf1("Number of model parameters: %d\n",NParm);
    if (CCTyp)
        printf1("Type of covariance matrix: %d\n",CCTyp);
    printf1("Maximum number of iterations: %d\n",MxIter);

    if (MINA == 1 || MINA == 2)  
        printf1("Starting value of step length: %lg\n",SLen);
       
    if (MINA == 1)   
        printf1("Step length reduction factor: %lg\n",SRed);
       
    printf1("Convergence criterion: %d\n",Crite);

    if (Crite == 1)
        printf1("Tolerance for norm of final gradient: %lg\n",TOLG);
    else if (Crite == 2)
        printf1("Tolerance for change in function value: %lg\n",TOLF);
    else if (Crite == 3)
        printf1("Tolerance for change in parameters: %lg\n",TOLP);
    else if (Crite == 4) 
        printf1("Tolerance for variance of function values: %lg\n",TOLV);
    else if (Crite == 5)  
        printf1("Minimum value of step length: %lg\n",TOLS);
    else if (Crite == 6) {
        printf1("Tolerance for scaled gradient: %lg\n",TOLSG); 
        printf1("Tolerance for scaled parameter change: %lg\n",TOLSP);
    }
    if (MINA == 5 || MINA == 6)  
        printf1("Mue of Armijo condition: %lg\n",AMue);
       
    if (MINA >= 3 && MINA <= 6)  
        printf1("Minimum of step size value: %lg\n",SMin);
       
    if (typ > 0 && DScal <= 0.0)
        DScal = 1.0;
   
    printf1("Scaling factor: %lg\n\n",DScal);

    /* set initial conditions */

    NumF = NumFG = NumFH = 0;

    /* memory allocation */

    if (!(Par = (double *)calloc(NParmA + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(NParmA + 1,sizeof(double));
    ParA = NParmA + 1;
              
    if (gmina) {
        if (!(ParLB = (double *)calloc(NParmA + 1,sizeof(double))))   
            goto MLIFin;         
        memrq(NParmA + 1,sizeof(double));
        ParLBA = NParmA + 1;

        if (!(ParUB = (double *)calloc(NParmA + 1,sizeof(double))))   
            goto MLIFin;         
        memrq(NParmA + 1,sizeof(double));
        ParUBA = NParmA + 1;
    }         
    if (!(WrkD = (double *)calloc(NParm + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(NParm + 1,sizeof(double));
    WrkDA = NParm + 1;

    if (!(WrkE = (double *)calloc(NParm + 1,sizeof(double))))   
        goto MLIFin;          
    memrq(NParm + 1,sizeof(double));
    WrkEA = NParm + 1;

    if (!(WrkH = (double *)calloc(HSiz + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(HSiz + 1,sizeof(double));
    WrkHA = HSiz + 1;

    if (!(Par1 = (double *)calloc(NParm + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(NParm + 1,sizeof(double));
    Par1A = NParm + 1;
              
    if (!(Grad = (double *)calloc(NParm + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(NParm + 1,sizeof(double));
    GradA = NParm + 1;
              
    if (!(Diag = (double *)calloc(NParmA + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(NParmA + 1,sizeof(double));
    DiagA = NParmA + 1;
              
    if (!(Hess = (double *)calloc(HSizA + 1,sizeof(double))))   
        goto MLIFin;         
    memrq(HSizA + 1,sizeof(double));
    HessA = HSizA + 1;
              
    if (CCTyp != 2 || MINA == 7 || MINA == 8) {
        if (!(GTmp = (double *)calloc(NParm + 1,sizeof(double))))   
            goto MLIFin;         
        memrq(NParm + 1,sizeof(double));
        GTmpA = NParm + 1;
    }         
    if (MINA >= 5 && MINA <= 8) {
        if (!(Srch = (double *)calloc(NParm + 1,sizeof(double))))   
            goto MLIFin;         
        memrq(NParm + 1,sizeof(double));
        SrchA = NParm + 1;
    }         
    if (MINA == 7 || MINA == 8) {
        if (!(TSHess = (double *)calloc(NParm * NParm + 1,sizeof(double))))   
            goto MLIFin;         
        TSHessA = NParm * NParm + 1;
        memrq(TSHessA,sizeof(double));
              
        if (!(TSGrad = (double *)calloc(NParm + 1,sizeof(double))))   
            goto MLIFin;         
        TSGradA = NParm + 1;
        memrq(TSGradA,sizeof(double));
    }         
    MLAlloc = 1;
    return(0);

MLIFin:
    if (SrchA) {
        free((char *)Srch);
        memrq(-SrchA,sizeof(double));
        SrchA = 0;
    }
    if (GTmpA) {
        free((char *)GTmp);
        memrq(-GTmpA,sizeof(double));
        GTmpA = 0;
    }
    if (HessA) {
        free((char *)Hess);
        memrq(-HessA,sizeof(double));
        HessA = 0;
    }
    if (DiagA) {
        free((char *)Diag);
        memrq(-DiagA,sizeof(double));
        DiagA = 0;
    }
    if (GradA) {
        free((char *)Grad);
        memrq(-GradA,sizeof(double));
        GradA = 0;
    }
    if (Par1A) {
        free((char *)Par1);
        memrq(-Par1A,sizeof(double));
        Par1A = 0;
    }
    if (WrkHA) {
        free((char *)WrkH);
        memrq(-WrkHA,sizeof(double));
        WrkHA = 0;
    }
    if (WrkEA) {
        free((char *)WrkE);
        memrq(-WrkEA,sizeof(double));
        WrkEA = 0;
    }
    if (WrkDA) {
        free((char *)WrkD);
        memrq(-WrkDA,sizeof(double));
        WrkDA = 0;
    }
    if (ParA) {
        free((char *)Par);
        memrq(-ParA,sizeof(double));
        ParA = 0;
    }
    if (ParLBA) {
        free((char *)ParLB);
        memrq(-ParLBA,sizeof(double));
        ParLBA = 0;
    }
    if (ParUBA) {
        free((char *)ParUB);
        memrq(-ParUBA,sizeof(double));
        ParUBA = 0;
    }
    if (TSHessA) {
        free((char *)TSHess);
        memrq(-TSHessA,sizeof(double));
        TSHessA = 0;
    }
    if (TSGradA) {
        free((char *)TSGrad);
        memrq(-TSGradA,sizeof(double));
        TSGradA = 0;
    }
    MLAlloc = 0;
    if (err)  
        p_err(-2,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  syminv1 Inversion of the hessian matrix with dimension                  */
/*          n = NParm or NParm1.                                            */
/*          The diagonal elements are assumed to be in d, the lower         */
/*          tiangular is assumed to be in h. The inversion is performed     */
/*          directly in the lower triangular. The function returns the      */
/*          rank of the matrix. For the determination of the pseudo rank    */
/*          the tolerance EPSI1 is used.                                    */
/*                                                                          */
/*  Return -1 if insufficient memory.                                       */

int syminv1(int n, double *d, double *h)
{
    register int i,j,k,kk,ii,jj;
    int aflag,r;
    double crit,tmp1,tmp2,tmp3,*wrk;

    if (!(wrk = (double *)calloc(n + 1,sizeof(double))))
        return(-1); 
    memrq(n + 1,sizeof(double));
                
    r = n;         
    for (i = 1; i <= n; ++i)            
        wrk[i] = d[i];

    for (k = 1; k <= n; ++k) {

        kk = (k - 1) * (k - 2) / 2;

        aflag = 0;                          /* flag for aliased parameters  */

        /* check for current sweep-column */

        tmp1 = d[k];
        if (tmp1 <= EPSI1 * wrk[k]) {
            aflag = 1;
        }
        else {
            for (i = 1; i < k; ++i) {
                tmp2 = d[i];
                if (tmp2 != 0.0) {
                    tmp3 = h[kk + i];
                    crit = 1.0 / (tmp2 + (tmp3 * tmp3) / tmp1);
                    if (crit < EPSI1 * wrk[i]) {
                        aflag = 1;
                        break;
                    }
                }
            }
        }
        if (!aflag) {

            /* collinearity check passed, so updating triangle */
 
            for (j = 1; j < k; ++j) {
                jj = (j - 1) * (j - 2) / 2;
                tmp2 = h[kk + j];
                tmp3 = tmp2 / tmp1;
                for (i = 1; i < j; ++i)  
                    h[jj + i] += h[kk + i] * tmp3;
                d[j] += tmp2 * tmp3;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                tmp2 = h[ii + k] / tmp1;
                for (j = k + 1; j <= i; ++j) {
                    jj = (j - 1) * (j - 2) / 2;
                    if (j < i)
                        h[ii + j] -= h[jj + k] * tmp2;
                    else
                        d[j] -= h[jj + k] * tmp2;
                }
            }
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                tmp2 = h[ii + k] / tmp1;
                for (j = 1; j < k; ++j)  
                    h[ii + j] -= h[kk + j] * tmp2;
            }
            for (j = 1; j < k; ++j)
                h[kk + j] /= tmp1;
     
            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                h[ii + k] /= -tmp1;
            }
            d[k] = 1.0 / tmp1;
        }
        else {
 
            /*  collinearity check not passed, so updating the k column     */
            /*  and k row with zeroes                                       */
 
            for (i = 1; i < k; ++i)  
                h[kk + i] = 0.0;
            d[k] = 0.0;

            for (i = k + 1; i <= n; ++i) {
                ii = (i - 1) * (i - 2) / 2;
                h[ii + k] = 0.0;
            }
            r--;
        }
    }
    memrq(-n - 1,sizeof(double));
    free((char *)wrk);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  syminv2 Inversion of a pos. def. symmetric matrix x with dimension n.   */
/*          It is assumed that the array x only contains the lower          */
/*          triangle, including the diagonal elements.                      */
/*                                                                          */
/*          The inversion is performed directly in the lower triangle.      */
/*          The function returns the rank of the matrix. For the            */
/*          determination of the pseudo rank tolerance EPSI1 is used.       */
/*                                                                          */
/*  Return -1 if insufficient memory.                                       */

int syminv2(int n, double *x)
{
    register int i,j,k,kk,ii,jj;
    int aflag,r;
    double crit,tmp1,tmp2,tmp3,*wrk;

    if (!(wrk = (double *)calloc(n + 1,sizeof(double))))
        return(-1); 
    memrq(n + 1,sizeof(double));
                
    r = n;         
    for (i = 1; i <= n; ++i)            
        wrk[i] = x[(i * i + i) / 2];

    for (k = 1; k <= n; ++k) {

        kk = ((k - 1) * k) / 2;

        aflag = 0;                          /* flag for aliased parameters  */

        /* check for current sweep-column */

        tmp1 = x[kk + k];

        if (tmp1 <= EPSI1 * wrk[k]) {
            aflag = 1;
        }
        else {
            for (i = 1; i < k; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + i];
                if (tmp2 != 0.0) {
                    tmp3 = x[kk + i];
                    crit = 1.0 / (tmp2 + (tmp3 * tmp3) / tmp1);
                    if (crit < EPSI1 * wrk[i]) {
                        aflag = 1;
                        break;
                    }
                }
            }
        }
        if (!aflag) {

            /* collinearity check passed, so updating triangle */
 
            for (j = 1; j < k; ++j) {
                jj = ((j - 1) * j) / 2;
                tmp2 = x[kk + j];
                tmp3 = tmp2 / tmp1;
                for (i = 1; i <= j; ++i)  
                    x[jj + i] += x[kk + i] * tmp3;
            }
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + k] / tmp1;
                for (j = k + 1; j <= i; ++j) {
                    jj = ((j - 1) * j) / 2;
                    x[ii + j] -= x[jj + k] * tmp2;
                }
            }
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                tmp2 = x[ii + k] / tmp1;
                for (j = 1; j < k; ++j)  
                    x[ii + j] -= x[kk + j] * tmp2;
            }
            for (j = 1; j < k; ++j)
                x[kk + j] /= tmp1;
     
            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                x[ii + k] /= -tmp1;
            }
            x[kk + k] = 1.0 / tmp1;
        }
        else {
 
            /*  collinearity check not passed, so updating the k column     */
            /*  and k row with zeroes                                       */
 
            for (i = 1; i <= k; ++i)  
                x[kk + i] = 0.0;

            for (i = k + 1; i <= n; ++i) {
                ii = ((i - 1) * i) / 2;
                x[ii + k] = 0.0;
            }
            r--;
        }
    }
    memrq(-n - 1,sizeof(double));
    free((char *)wrk);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  fn  Call of model-specific functions for calculation of loglikelihood,  */
/*      gradient and hessian. Parameters are assumed to be in x[], a        */
/*      vector of length NParm1. If constraints are defined the vector is   */
/*      expanded before calling the model-specific functions.               */
/*                                                                          */
/*      The function value is calculated in FTmp. The number of function    */
/*      calls is counted in the global variable NumF, calls with            */
/*      derivatives are counted in NumFG and NumFH.                         */
/*                                                                          */
/*      If scal != 0, the calculated function value, the gradient, and      */
/*      the hessian, are scaled with DScal before returning.                */
/*                                                                          */
/*      If typ = 0  only function value is calculated                       */
/*      If typ = 1  function value and gradient are calculated              */
/*      If typ = 2  function value, gradient, and hessian are calculated    */
/*                                                                          */
/*      If opt == 1, calculation of score test statistics based on log-     */
/*      likelihood in the unrestricted parameter space.                     */
/*      Results saved in FNTStat and FNTRank.                               */
/*                                                                          */
/*      If first call of fn(), and if GCTestFlg = 1, calculation and        */
/*      print of global chi square (score test statistic with actual        */
/*      starting values). Only done if typ = 2.                             */
/*      First call recognized by global variable FN_First.                  */
/*                                                                          */
/*      The actual number of cases used to evaluate the function is         */
/*      counted in NOCUsed.                                                 */
/*                                                                          */
/*      Return err = -1 if insufficient memory, 1 if error in function.     */

void fn(double *x, int mod, int typ, int scal, int opt,int *err)
{
    register int i,j;
    double tmp; /** ,tmp1, *d,*h **/  

    *err = 0;
    NumF++;
    LFunc = 1;
    NOCUsed = LGrad = LSec = 0;
    if (typ == 1) {
        LGrad = 1;
        NumFG++;
    }
    else if (typ == 2) {    /* LSec = 1 implies LGrad = 1 */
        LGrad = LSec = 1;
        NumFG++;
        NumFH++;
    }
    if (NCon1) {                             /* expand x */
        for (i = 1; i <= NParm; ++i) {
            tmp = 0.0;
            for (j = 1; j <= NParm1; ++j)
                tmp += ConQ[(CIP[i] - 1) * NParm1 + j] * x[j];
            ParC[i] = ParS[i] + tmp;
        }
        TPar = ParC;
    }
    else
        TPar = x;

    if (PMProtFDef >= 3) {
        fprintf(PMProtFd,"\nFunction evaluation %d (%d,%d). ",NumF,NumFG,NumFH);
        fprintf(PMProtFd,"Func %d Grad %d Sec %d\n",LFunc,LGrad,LSec);
    }
    if (LFunc)  
        FTmp = 0.0;                                              
    if (LGrad)
        dclear(NParm,Grad);                               
    if (LSec) {
        dclear(NParm,Diag);                                 
        dclear(HSiz,Hess);                                   
    }
    switch (mod) {

        /* quantal response models */

        case QRLOG1:    *err = fn_logit1();
                        break;
        case QRLOG2:    *err = fn_logit2();
                        break;
        case QRLOG3:    *err = fn_logit3();
                        break;
        case QRLOG4:    *err = fn_logit4();
                        break;

        case QRPROB1:   *err = fn_probit1();
                        break;
        case QRPROB2:   *err = fn_probit2();
                        break;
        case QRPROB3:   *err = fn_probit3();
                        break;
        case QRPROB4:   *err = fn_probit4();
                        break;

        case RMOD1:     *err = fn_rmod1();      
                        break;
                        
        /* transition rate models */

        case MCOX:      *err = fn_cox();      
                        break;
        case MEXP:      *err = fn_exp(0);   
                        break;
        case MEXP1:     *err = fn_exp1(0);  
                        break;
        case MEXP2:     *err = fn_exp2(0);  
                        break;
        case MPOL:      *err = fn_pol(0);  
                        break;
        case MPOL1:     *err = fn_pol1(0);
                        break;
        case MGM:       *err = fn_gm(0);   
                        break;
        case MWEI:      *err = fn_wei(0);
                        break;
        case MSIC:      *err = fn_sic(0);   
                        break;
        case MLL:       *err = fn_ll(0);
                        break;
        case MLL2:      *err = fn_ll2(0);
                        break;
        case MLN:       *err = fn_ln(0);   
                        break;
        case MIG:       *err = fn_ig(0);   
                        break;
        case MGAM:      *err = fn_gam(0);      
                        break;
        case DLR:   
        case CLL:       *err = fn_dis(mod);  
                        break;

        /* general functions */

        case CLMOD1:    *err = clu_fn();   /* clu */
                        break;

        /*******************************************
        case MDSMOD1:   *err = mds1_fn();   
                        break;
        case MDSMOD2:   *err = mds2_fn();   
                        break;
        ********************************************/
        case MDSMOD3:   *err = mds3_fn();   /* MDS */
                        break;
        case MDSMOD4:   *err = mds4_fn();   /* MDS */
                        break;

        case IVRMOD1:   *err = ivr1_fn();   /* interval regression */
                        break;

        case UMOD:      *err = fn_umod();   /* user-defined models */
                        break;

        default:        printf1("Error: fn %d not available\n",mod);
                        exit(0);
    }
    if (PMProtFDef >= 3) {
        prval("Function",FTmp);
        prvec("Parameter",NParm,TPar);
        if (LGrad)  
            prvec("Gradient",NParm,Grad);
        if (LSec)   
            prhess("Hessian",NParm);
    }   
    if (*err)
        return;
           
   
    /* calculate score test statistic if opt == 1, or if first call
       of fn() and GCTestFlg = 1. */
    /***
    if (opt == 1 || (GCTestFlg && FN_First)) {

        memrq(NParm + HSiz + 2,sizeof(double));
        if (!(d = (double *)calloc(NParm + 1,sizeof(double))) ||
            !(h = (double *)calloc(HSiz + 1,sizeof(double))))   
            p_err(-2,1);
         
        for (i = 1; i <= NParm; ++i)    
            d[i] = -Diag[i];
        for (i = 1; i <= HSiz; ++i)
            h[i] = -Hess[i];

        check memory    

        if ((FNTRank = syminv1(NParm,d,h)) != NParm) {   
            FNTStat = -1.0;
        }
        else {
            FNTStat = 0.0;
            for (i = 1; i <= NParm; ++i) {
                tmp1 = 0.0;
                for (j = 1; j <= NParm; ++j) {

                    if (i == j)
                        tmp = d[j];
                    else if (i < j)
                        tmp = h[(j - 1) * (j - 2) / 2 + i];
                    else
                        tmp = h[(i - 1) * (i - 2) / 2 + j];

                    tmp1 += Grad[j] * tmp;
                }
                FNTStat += tmp1 * Grad[i];
            }
        }
        memrq(-NParm - HSiz - 2,sizeof(double));
        free((char *)d);
        free((char *)h);

        if (GCTestFlg && FN_First) {
            GCTest = FNTStat;
            GCRank = FNTRank;
        }
    }
    *******/

    if (LFunc && FN_First)
        FMax1 = FTmp;

    if (LFunc && scal)
        FTmp *= DScal;                                      

    if (LGrad) {
        if (scal) {
            for (j = 1; j <= NParm; ++j)  
                Grad[j] *= DScal;                              
        }
        if (NCon1)   
            prjvec(Grad);
    }
    if (LSec) {
        if (scal) {
            for (j = 1; j <= NParm; ++j)   
                Diag[j] *= DScal;                                 
            for (j = 1; j <= HSiz; ++j)
                Hess[j] *= DScal;                                             
        }
        if (NCon1)  
            prjhess();
    }
    FN_First = 0;
}

/* ------------------------------------------------------------------------ */
/*  fn_umod()   calculate function, gradient, hessian for user-defined      */
/*              models. Use parameter values in TPar[] (i=1,NParm).         */
/*              if LFunc  return function value in FTmp.                    */
/*              if LGrad  return gradient in Grad[]                         */
/*              if LSec   return hessian in Hess[] and Diag[]               */
/*                                                                          */
/*              If CCovFlg = 1 return outer product of gradients instead    */
/*              of Hessian.                                                 */
/*                                                                          */
/*              Return 0 if OK, -1 if insuff memory, 1 if error in function */

int fn_umod(void)
{
    int i,err,deriv;

    NOCUsed = deriv = 0;
    if (LGrad)
        deriv = 1;
    if (LSec)  
        deriv = 2;

    err = get_flval(&FTmp,NParm,TPar + 1,deriv,2,Grad,Hess,Diag);

    if (err) {
        printf1("\nError in function evaluation.\n");
        prn_emsg2(err);
        printf1("Current function arguments:");
        for (i = 1; i <= NParm; ++i)
            printf1(" %lg",TPar[i]);
        printf1("\n");
        err = 1;
    }
    else {
        err = checkov();          /* check overflow */
        if (err) {
            printf1("\nError in function evaluation.\n");
            printf1("Numerical overflow.\n");
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prjvec(x)   Project the NParm vector x to NParm1 = NParm - NCon1 by     */
/*              multiplication with the transposed projection matrix ConQ.  */
/*              Use of WrkD as temporary storage.                           */

void prjvec(double *x)
{
    register int i,j;
    double tmp;

    for (i = 1; i <= NParm1; ++i) {
        tmp = 0.0;
        for (j = 1; j <= NParm; ++j)  
            tmp += ConQ[(CIP[j] - 1) * NParm1 + i] * x[j];
        WrkD[i] = tmp;
    }
    for (i = 1; i <= NParm1; ++i)  
        x[i] = WrkD[i];
}

/* ------------------------------------------------------------------------ */
/*  prjhess     Project the Hessian (Hess and Diag) into the NParm1 =       */
/*              NParm - NCon1 space.                                        */
/*              Use of WrkH as temporary storage.                           */

void prjhess(void)
{
    register int i,j,k,l;
    double tmp;

    for (i = 1; i <= NParm1; ++i) {
        for (j = 1; j <= i; ++j) {
            tmp = 0.0;
            for (k = 1; k <= NParm; ++k) {
                for (l = 1; l <= NParm; ++l) {
                    if (l == k)
                        tmp += Diag[k] * ConQ[(CIP[l] - 1) * NParm1 + i] * 
                                         ConQ[(CIP[k] - 1) * NParm1 + j];
                    else if (l < k)
                        tmp += Hess[(k - 1) * (k - 2) / 2 + l] *
                                         ConQ[(CIP[l] - 1) * NParm1 + i] * 
                                         ConQ[(CIP[k] - 1) * NParm1 + j];
                    else
                        tmp += Hess[(l - 1) * (l - 2) / 2 + k] *
                                         ConQ[(CIP[l] - 1) * NParm1 + i] * 
                                         ConQ[(CIP[k] - 1) * NParm1 + j];
                }
            }
            if (i == j)
                WrkD[i] = tmp;
            else
                WrkH[(i - 1) * (i - 2) / 2 + j] = tmp;
        }
    }
    for (i = 1; i <= NParm1; ++i)  
        Diag[i] = WrkD[i];
    for (i = 1; i <= HSiz1; ++i)
        Hess[i] = WrkH[i];
}

/* ------------------------------------------------------------------------ */
/*  checkov.  Check overflow and underflow. OFMax and OFMin are used as     */
/*            ceilings and floors.                                          */
/*            Return 0 if OK,  1 if error.                                  */

int checkov(void)
{
    register int j,k = 0;

    if (LFunc) {
        if (FTmp <= OFMin) {
            k = 1;
            FTmp  = OFMin;
        }
        if (FTmp >= OFMax) {
            k = 1;  
            FTmp  = OFMax;
        }
    }
    if (LGrad) {
        for (j = 1; j <= NParm; ++j) {
            if (Grad[j] <= OFMin) {
                k = 1;           
                Grad[j] =  OFMin;
            }
            if (Grad[j] >= OFMax) {
                k = 1;
                Grad[j]  = OFMax;
            }
        }
    }
    if (LSec) {
        for (j = 1; j <= NParm; ++j) {
            if (Diag[j] <= OFMin) {
                k = 1;  
                Diag[j] =  OFMin;
            }
            if (Diag[j] >= OFMax) {
                k = 1;   
                Diag[j]  = OFMax;
            }
        }
        for (j = 1; j <= HSiz; ++j) {
            if (Hess[j] <= OFMin) {
                k = 1;   
                Hess[j] =  OFMin;
            }
            if (Hess[j] >= OFMax) {
                k = 1;   
                Hess[j]  = OFMax;
            }
        }
    }
    NPFlgs[2] += k; 
    return(k); 
}

/* ------------------------------------------------------------------------ */
/*  prn_mlres(typ)      Print results of function minimization.             */
/*                      typ 0 : ML                                          */
/*                          1 : general function                            */
/*                          2 : nonlinear regression                        */

void prn_mlres(int typ)
{
    printf1("\nConvergence ");
    if (LConv) 
        printf1("not ");
    printf1("reached in %d iterations.\n",Iter);
    printf1("Number of function evaluations: %d (%d,%d)\n",NumF,NumFG,NumFH);
    if (NINTMUsed >= 0)
        ni_info();

    if (LConv) {
        printf1("\nProblem: ");

        switch (LConv) {
            case  1:    if (typ == 0)
                            printf1("cannot increase log likelihood.\n");
                        else
                            printf1("cannot decrease function.\n");
                        break;
            case  2:    printf1("reached max number of iterations.\n");
                        break;
            case -1:    printf1("final hessian (or outer product) not positive definite.\n");
                        break;
            case -2:    printf1("step size search failed.\n");
                        break;
            case -3:    printf1("cannot calculate eigenvalues and vectors of the hessian.\n");
                        break;
            case -4:    printf1("search vector is no descent direction.\n");
                        break;
            case -5:    printf1("exceeded max iterations in step size search.\n");
                        break;
            default:    printf1("LConv = %d\n",LConv);
                        break;
        }
    }
    if (typ == 0)
        printf1("\nMaximum of log likelihood: %lg",FMax);
    else
        printf1("\nMinimum of function: %lg",FMax);

    if (MINA == 2)  
        printf1("\nVariance of final function values: %lg",CValV);
       
    if (LConv >= 0) {

        if (CValG)
            printf1("\nNorm of final gradient vector: %lg",CValG);

        if (MINA == 5 || MINA == 6) {
            printf1("\nLast absolute change of function value: %lg",CValF);
            printf1("\nLast relative change in parameters: %lg",CValP);
        }
        else if (MINA == 7 || MINA == 8) {
            printf1("\nFinal scaled gradient: %lg",TOLSGF);
            printf1("\nFinal scaled parameter change: %lg",TOLSPF);
        }
    }
    printf1("\n\n");
    prn_npflg();           /* print numerical problems if any */
}

/* ------------------------------------------------------------------------ */
/*  prn_ml1res(typ)     Print additional results of func minimization.      */
/*                      typ 0 : ML                                          */
/*                          1 : general function                            */
/*                          2 : nonlinear regression                        */
/*                                                                          */
/*                      do only if LConv >= 0.                              */
/*                      if CCTyp > 0 and PMCovFDef print cov. matrix.       */

void prn_ml1res(int typ)
{
    if (typ == 0)
        printf1("Log likelihood");
    else
        printf1("Function");
    printf1(" (starting values): ");
    printf1(PMTFmtS,FMax1); 

    if (typ == 0)
        printf1("\nLog likelihood");
    else
        printf1("\nFunction");
    printf1(" (final estimates): "); 
    printf1(PMTFmtS,FMax); 
    newline();    
    newline();
    if (PMPPWFlg || PMCovWFlg) {
        if (PMPPWFlg)  
            printf1("Parameter estimates written to: %s\n",PMPPFName);
        if (PMCovWFlg)  
            printf1("Covariance matrix written to: %s\n",PMCovFName);
        newline();
    }
    mp_info();              /* print info about matrices */

    /***
    if (GCTestFlg)         score test with starting values   
        gctest();

    if (LMTestFlg)         LM test   
        lmtest();
      
    if (NWTest)            Wald tests   
        wtest();
    ****/
}

/* ------------------------------------------------------------------------ */
/*  prvec(txt,n,x)                                                          */
/*      Prints a vector with n elements. The string txt is printed first.   */
/*      Used for iteration protocol.                                        */

void prvec(char *txt, int n, double *x)
{
    register int i; 

    fprintf(PMProtFd,"%s\n",txt); 
    for (i = 1; i <= n; ++i)  
        fprintf(PMProtFd,PMPFmtS,x[i]);
    fprintf(PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prmat(txt,n,m,x)                                                        */
/*      Prints a matrix with n,m elements. The string txt is printed first. */
/*      Used for protocol file.                                             */

void prmat(char *txt,int n,int m,double *x)
{
    register int i,j;

    fprintf(PMProtFd,"%s\n",txt); 
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= m; ++j) 
            fprintf(PMProtFd,PMPFmtS,x[i * m + j]);
        fprintf(PMProtFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prhess(txt,n)       Print Hessian to protocol file.                     */

void prhess(char *txt, int n)
{
    register int i,j,k;

    fprintf(PMProtFd,"%s\n",txt); 
    k = 1;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j < i; ++j)
            fprintf(PMProtFd,PMPFmtS,Hess[k++]);
        fprintf(PMProtFd,PMPFmtS,Diag[i]);
        fprintf(PMProtFd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prval(txt,x)                                                            */
/*      Print a double value x. The string txt is printed first.            */
/*      Used for iteration protocol.                                        */

void prval(char *txt, double x)
{
    fprintf(PMProtFd,"%s\n",txt); 
    fprintf(PMProtFd,PMPFmtS,x);  
    fprintf(PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  prot_init(typ,mod,gmina)                                                */
/*                                                                          */
/*  Init protocol file                                                      */
/*                      typ 1 = transition rate models                      */
/*                          2 = quantal response models                     */
/*                          3 = user-def. log-likelihood                    */
/*                          4 = general function                            */
/*                          5 = nonlinear regression                        */
/*                          6 = user-def rate model                         */
/*                                                                          */
/*  If typ 1 add some information about transitions and parameter vector.   */
/*  If gmina != 0 print lower/upper parameter bounds.                       */
/*                                                                          */

void prot_init(int typ,int mod,int gmina)
{
    register int i,j;

    if (PMProtFDef == 0)
        return;

    printf1("Protocol will be written to: %s\n",PMProtFName);

    switch (typ) {
        case 1:     fprintf(PMProtFd,"Transition rate model: %d\n\n",mod);
                    break;
        case 2:     fprintf(PMProtFd,"Quantal response model: %d\n\n",mod);
                    break;
        case 3:     fprintf(PMProtFd,"User-defined log-likelihood.\n\n");
                    break;
        case 4:     fprintf(PMProtFd,"General function minimization.\n\n");
                    break;
        case 5:     fprintf(PMProtFd,"Nonlinear regression.\n\n");
                    break;
        case 6:     fprintf(PMProtFd,"User-defined rate model.\n\n");
                    break;
        default:    break;
    }
    fprintf(PMProtFd,"Number of model parameters: %d\n",NParm);
    fprintf(PMProtFd,"Starting values");
    if (gmina)
        fprintf(PMProtFd," and bounding box");
    fprintf(PMProtFd,".\n");

    for (i = 1; i <= NParm; ++i) {
        fprintf(PMProtFd,"%3d  ",i);
        fprintf(PMProtFd,PMPFmtS,Par[i]);
        if (gmina) {
            fprintf(PMProtFd,PMPFmtS,ParLB[i]);
            fprintf(PMProtFd,PMPFmtS,ParUB[i]);
        }
        fprintf(PMProtFd,"\n");
    }
    fprintf(PMProtFd,"\n");

    if (typ == 1) {     /* transition rate models */

        fprintf(PMProtFd,"Number of transitions: %d\n",NTran1);

        for (i = 1; i <= MaxSnn; ++i) {
            fprintf(PMProtFd,"TranPtr:  %2d  :",i);
            for (j = 0; j <= MaxOrg; ++j)
                fprintf(PMProtFd," %2d",TranPtr[(i - 1) * MaxOrg1 + j]);
            fprintf(PMProtFd,"\n");
        }
        prot_prni("SnTran1: ",NTran1,SnTran1);
        prot_prni("OrgTran1:",NTran1,OrgTran1);
        prot_prni("DesTran1:",NTran1,DesTran1);

        prot_prns("VTNum:   ",NTran1,VTNum);
        prot_prns("PIdxPtr: ",NTran1,PIdxPtr);
        prot_prns("PIdx:    ",NParm,PIdx + 1);
        prot_prnc("PIdxM:   ",NParm,PIdxM + 1);
        fprintf(PMProtFd,"\n");
    }
}

void prot_prni(char *s,int n,int *x)
{
    register int i;

    fprintf(PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(PMProtFd," %2d",x[i]);
    fprintf(PMProtFd,"\n");
}

void prot_prns(char *s,int n,short *x)
{
    register int i;

    fprintf(PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(PMProtFd," %2d",(int)x[i]);
    fprintf(PMProtFd,"\n");
}

void prot_prnc(char *s,int n,char *x)
{
    register int i;

    fprintf(PMProtFd,"%s",s);
    for (i = 0; i < n; ++i)
        fprintf(PMProtFd," %2d",(int)x[i]);
    fprintf(PMProtFd,"\n");
}




