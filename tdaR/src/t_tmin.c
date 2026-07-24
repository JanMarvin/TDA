/****************************************************************************/
/*  t_tmin                                                                  */
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
#include "t_gmin.h"  
#include "t_lin.h"  
#include "t_gf.h"
#include "tda_context.h"

/*  functions in t_tmin.c */

int ts_min(TDAContext *ctx, int n,double *x,double *typx,double *xpls,double *g,double *gpls, double *s,double *d,double *dn,double *e,double *wk1,double *wk2, double *h,int *pivot);

int lnsrch(TDAContext *ctx, int n,double *x,double f,double *g,double *p,double *xpls, double *fpls,int *mxtake,double stepmx,double *typx, double *w2);
void zhz(TDAContext *ctx, int nc,int n,double *y,double *h,double *u,double *t);
void solvew(TDAContext *ctx, int nc,int n,double *al,double *u,double *w,double *b);
void dstar(TDAContext *ctx, int nc,int n,double *u,double *s,double *w1,double *w2,double *w3, double sigma,double *al,double *d);
void mkmdl(TDAContext *ctx, int nc,int n,double f,double fp,double *g,double *gp,double *s, double *h,double *alpha,double *beta,double *sh,double *a);
void sigma(TDAContext *ctx, double *sgstar,double a,double b,double c,double d);
void roots(TDAContext *ctx, double *s1,double *s2,double *s3,double a,double b,double c,double d);
void sortrt(TDAContext *ctx, double *s1,double *s2,double *s3);
int fstofd(TDAContext *ctx, int nc,int m,int n,double *xpls,double *fpls,double *a,double *typx, double rnoise,double *fhat,int icase);
int sndofd(TDAContext *ctx, int nc,int n,double *xpls,double fpls,double *a,double *typx, double rnoise,double *stepsz,double *anbr);
void bakslv(TDAContext *ctx, int nc,int n,double *a,double *x,double *b);
void forslv(TDAContext *ctx, int nc,int n,double *a,double *x,double *b);
void choldr(TDAContext *ctx, int nc,int n,double *h,double *g,int *pivot,double *e, double *diag,double addmax);
void modchl(TDAContext *ctx, int ndim,int n,double *a,double *g,double tau1, double tau2,int *p,double *e);
void ts_init(TDAContext *ctx, int n,int ndim,double *a,int *phase1,double *delta,int *p,double *g, double *e,double *ming,double tau1,double *gamma,double *taugam);
void gersch(TDAContext *ctx, int ndim,int n,double *a,int j,double *g);
void fin2x2(TDAContext *ctx, int ndim,int n,double *a,double *e,double tau2,double *delta, double gamma);
void slvmdl(TDAContext *ctx, int nc,int n,double *h,double *u,double *t,double *e,double *diag, double *s,double *g,int *pivot,double *w1,double *w2,double *w3, double alpha,double beta,int *nomin);
int optstp(TDAContext *ctx, int n,double *xpls,double fpls,double *gpls,double *x,            int *icscmx,double fscale, int mxtake,double *rgx,double *rsx,int iretcd);
int ts_func(TDAContext *ctx, int n,double *x,double *f,int gflg,double *g,int hflg,double *h);
void ts_prn(TDAContext *ctx, int iter,double f,double gnorm,double pc,int n,double *x,double *g);
void ts_prn1(TDAContext *ctx, int n,double *h);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

                    /* 1 if tensor model                                    */

                    /* 1 if also analytical gradient                        */
                    /* 2 if also analytical hessian                         */


/*  variables to be used for calling ts_min(). */


/* ------------------------------------------------------------------------ */
/*  ts_min()                                                                */
/*                                                                          */
/*  Unconstrained function minimization with a tensor method.               */
/*                                                                          */
/*  Code adaped from: ACM algorithm 738. T. Chow, E. Eskow, R. Schnabel,    */
/*  A Software Package for Unconstrained Optimization Using Tensor Methods. */
/*  ACM Transactions on Mathematical Software 20 (1994), pp. 518 - 530.     */
/*                                                                          */
/*  ts_min(int n,double *x,double *typx,double *xpls,double *g,double *gpls,*/
/*        double *s,double *d,double *dn,double *e,double *wk1,double *wk2, */
/*        double *h,int *pivot)                                             */
/*                                                                          */
/*  n       ->  number of parameters                                        */
/*  x      <->  starting values, and final parameter vector                 */
/*  typx    ->  typical size of x, used for scaling                         */
/*  gpls   <-   final gradient                                              */
/*  h      <-   final hessian                                               */
/*                                                                          */
/*  There are two different methods.                                        */
/*  Method 0 : only Newton steps with a quadratic model,                    */
/*  Method 1 : Tensor model.                                                */
/*                                                                          */
/*  controlled by the global variables                                      */
/*                                                                          */
/*  TSMETH =           0 if quadratic (Newton) model                        */
/*                     1 if tensor model                                    */
/*                                                                          */
/*  TSDERIV =          0 if only function values                            */
/*                     1 if also analytical gradient                        */
/*                     2 if also analytical hessian                         */
/*                                                                          */
/*  The algorithm uses:                                                     */
/*  MxIter   =   max number of iterations                                   */  
/*  Iter     =   iteration counter                                          */  
/*  FMin     =   final function value                                       */
/*  TOLSG    =   tolerance for scaled gradient (= gradtl)                   */
/*  TOLSP    =   tolerance for scaled parameter change (= steptl)           */
/*  PMNDIGIT =  number of significant digits in func calculation.           */
/*                                                                          */
/*  Function evaluation is done by ts_func() and depends on:                */
/*  TSFNC = 0   then use standard models.                                   */
/*          1   then use user-defined functions.                            */
/*  TSMOD =     model number for TSFNC = 0.                                 */
/*                                                                          */
/*  Return:  1  convergence reached with gradient check                     */
/*           2  convergence reached with parameter change                   */
/*           3  no success in last line search                              */
/*           4  exceeded max number of iterations                           */
/*           5  five consecutive max steps in line search                   */
/*          -1  if error in function evaluation (only with TSFNC = 1)       */
 
int ts_min(TDAContext *ctx, int n,double *x,double *typx,double *xpls,double *g,double *gpls, double *s,double *d,double *dn,double *e,double *wk1,double *wk2, double *h,int *pivot)
{
    register int i;
    int err,gflg,hflg,icscmx,iretcd,itrmcd,nomin,mxtake;               
    double f,fn,fp,fscale,temp,fpls,rnf,gnorm;
    double addmax,tmp,rgx,rsx,alpha,beta,gd,stepmx;
    /* fstofd() indexes fpls from 1; the original passed &f - 1, which
       forms a pointer before the object.  A two-element array gives the
       same fpls[1] without leaving the object. */
    double f1[2];

    /* initialization */

    ctx->Iter = 0;                   /* iteration counter */
    ctx->NumF = ctx->NumFG = ctx->NumFH = 0;   /* function calls */

    iretcd = 2;                 /* for initial call of optstp(ctx) */
    icscmx = 0;                 /* accumulate consecutive max steps */
    mxtake = 0;                 /* indicates max step in line search */

    nomin = gflg = hflg = 0;
    addmax = 0.0;

    if (ctx->TSDERIV == 1)     /* if analytical gradient */
        gflg = 1;

    if (ctx->TSDERIV == 2)     /* if analytical hessian */
        hflg = gflg = 1;

    fscale = 1.0;
    for (i = 1; i <= n; ++i)
        typx[i] = 1.0;
                   
    /***************************
    temp = pow(EPSI,1.0 / 3.0);
    gradtl = temp;
    steptl = temp * temp;
    *************************/

    stepmx = 0.0;

    /* check input parameter */

    for (i = 1; i <= n; ++i) {                                          
        if (typx[i] < ctx->EPSI1)
            typx[i] = 1.0;                                
    }
    if (stepmx <= 0.0) {             
        tmp = 0.0;
        for (i = 1; i <= n; ++i)
            tmp += x[i] * x[i] / (typx[i] * typx[i]);
        tmp = sqrt(tmp);
        stepmx = dmax(ctx, 1000.0 * tmp,1000.0);
    }
    /****
    printf1(ctx, "stepmx=%16.8e\n",stepmx);
    ****/

    for (i = 1; i <= n; ++i)        /* scale x */
        x[i] /= typx[i];

    /* initial iteration */
    /* compute typical size of x */

    for (i = 1; i <= n; ++i)  
        dn[i] = 1.0 / typx[i];

    /*****************************
    tmp =  -log(EPSI) / log(10.0);
    printf1(ctx, "tmp... =%lg\n",tmp);

    rnf = dmax(ctx, pow(10.0,-tmp),EPSI);
    printf1(ctx, "rnf=...%lg\n",rnf);
    *******************************/

    rnf = dmax(ctx, pow(10.0,-ctx->PMNDIGIT),ctx->EPSI);                      
    dmax(ctx, 0.01,sqrt(rnf));

    /* unscale x and compute f and g */

    for (i = 1; i <= n; ++i)
        wk1[i] = x[i] * typx[i];
  
    err = ts_func(ctx, n,wk1,&f,gflg,g,hflg,h);
    if (err)
        return(-1);

    if (gflg == 0) {    /* finite difference gradient */

        f1[1] = f;
        err = fstofd(ctx, 1,1,n,wk1,f1,g,typx,rnf,wk2,1);
        if (err)
            return(-1);
    } 
    gnorm = twonrm(ctx, n,g);

    ts_prn(ctx, ctx->Iter,f,gnorm,0.0,n,wk1,g);
   
    /* TEST WHETHER THE INITIAL GUESS SATISFIES THE STOPPING CRITERIA */

    if (gnorm <= ctx->TOLSG) {       
        fpls = f;
        for (i = 1; i <= n; ++i) {
            xpls[i] = x[i];
            gpls[i] = g[i];
        }
        itrmcd = optstp(ctx, n,xpls,fpls,gpls,x,&icscmx,fscale,
                    mxtake,&rgx,&rsx,iretcd);
        goto L350;
    }

    /* iteration 1 */
  
    ctx->Iter++;         
 
    /* compute Hessian if not already in h[] */

    if (hflg == 0) {         
        /**********     not used
        if (gflg)           
            fstofd(ctx, n,n,n,wk1,g,h,typx,rnf,wk2,3);                   
        else
        ***********/
        err = sndofd(ctx, n,n,wk1,f,h,typx,rnf,wk2,d);
        if (err)
            return(-1); 
    }
    ts_prn1(ctx, n,h);

    for (i = 2; i <= n; ++i)  
        dcopy(ctx, i - 1,h + (i - 1) * n,1,h + i - 1,n);

    /* CHOLESKY DECOMPOSITION FOR H (H=LLT) */
   
    choldr(ctx, n,n,h,d,pivot,e,wk1,addmax);
 
    /* SOLVE FOR NEWTON STEP D */

    for (i = 1; i <= n; ++i)
        wk2[i] = -g[i];

    forslv(ctx, n,n,h,wk1,wk2);
    bakslv(ctx, n,n,h,d,wk1);
 
    /* APPLY LINESEARCH TO THE NEWTON STEP */
              
    iretcd = lnsrch(ctx, n,x,f,g,d,xpls,&fpls,&mxtake,stepmx,typx,wk1);
    if (iretcd < 0)
        return(-1);

    /* UPDATE G */
    /* CALL DCOPY(N,GPLS(1),1,GP(1),1)  */
    /* UNSCALE XPLS AND COMPUTE GPLS */

    for (i = 1; i <= n; ++i) 
        wk1[i] = xpls[i] * typx[i];

    if (gflg)                                       /* compute gradient */
        err = ts_func(ctx, n,wk1,&tmp,1,gpls,hflg,h);
    else {
        f1[1] = fpls;
        err = fstofd(ctx, 1,1,n,wk1,f1,gpls,typx,rnf,wk2,1);
    }

    if (err)
        return(-1);

    /* CHECK STOPPING CONDITIONS */
     
    itrmcd = optstp(ctx, n,xpls,fpls,gpls,x,&icscmx,fscale,
                    mxtake,&rgx,&rsx,iretcd);
 
    /* IF ITRMCD > 0 THEN STOPPING CONDITIONS SATISFIED */

    if (itrmcd > 0)             
        goto L350;

    /* UPDATE X,F AND S FOR TENSOR MODEL */

    fp = f;
    f = fpls;
    for (i = 1; i <= n; ++i) {
        temp = xpls[i];
        s[i] = x[i] - temp;
        x[i] = temp;
    }
    ts_prn(ctx, ctx->Iter,fpls,rgx,rsx,n,x,gpls);
      
    /* iteration > 1 */
   
    /* UNSCALE X AND COMPUTE H */

L200:
    for (i = 1; i <= n; ++i)
        wk1[i] = x[i] * typx[i];

    if (hflg == 0) {  /* calculate hessian if not already in h[] */
        /*************** not used
        if (gflg)           
            fstofd(ctx, n,n,n,wk1,g,h,typx,rnf,wk2,3);                   
        else
        *****************/
        err = sndofd(ctx, n,n,wk1,f,h,typx,rnf,wk2,d);
        if (err)
            return(-1);
    }
    ts_prn1(ctx, n,h);

    for (i = 2; i <= n; ++i)
        dcopy(ctx, i - 1,h + (i - 1) * n,1,h + i - 1,n);
 
    /* IF METHOD = 0 THEN USE NEWTON STEP ONLY */

    if (ctx->TSMETH == 0) { 

        /* CHOLESKY DECOMPOSITION FOR H */

        choldr(ctx, n,n,h,wk2,pivot,e,wk1,addmax);
 
        /* COMPUTE NEWTON STEP */

        for (i = 1; i <= n; ++i)
            wk1[i] = -gpls[i];

        forslv(ctx, n,n,h,wk2,wk1);
        bakslv(ctx, n,n,h,d,wk2);
 
        nomin = 1;      /* no tensor step */
        goto L300;
    }           

    /* FORM TENSOR MODEL */

    mkmdl(ctx, n,n,f,fp,gpls,g,s,h,&alpha,&beta,wk1,d);
 
    /* SOLVE TENSOR MODEL AND COMPUTE NEWTON STEP */
    /* ON INPUT : SH IS STORED IN WK1 */
    /*            A=(G-GPLS-SH-S*BETA/(6*STS)) IS STORED IN D */
    /* ON OUTPUT: NEWTON STEP IS STORED IN DN */
    /* TENSOR STEP IS STORED IN D */

    slvmdl(ctx, n,n,h,xpls,wk2,e,g,s,gpls,pivot,d,wk1,dn,alpha,beta,&nomin);
 
    /* IF TENSOR MODEL HAS NO MINIMIZER THEN USE NEWTON STEP */

    if (nomin) {     
        dcopy(ctx, n,dn,1,d,1);
        goto L300; 
    }             
 
    /* IF TENSOR STEP IS NOT IN DESCENT DIRECTION THEN USE NEWTON STEP */

    gd = ddot(ctx, n,gpls,1,d,1);
    if (gd > 0.0) {     
        dcopy(ctx, n,dn,1,d,1);
        nomin = 1;    
    }               
 
L300:
    ctx->Iter++;
    dcopy(ctx, n,gpls,1,g,1);
 
    /* APPLY LINESEARCH TO TENSOR (OR NEWTON) STEP */

    iretcd = lnsrch(ctx, n,x,f,g,d,xpls,&fpls,&mxtake,stepmx,typx,wk1);
    if (iretcd < 0)
        return(-1);
 
    if (nomin == 0) {      

        /* TENSOR STEP IS FOUND AND IN DESCENT DIRECTION, */
        /* APPLY LINESEARCH TO NEWTON STEP */
        /* NEW NEWTON POINT IN WK2 */

        iretcd = lnsrch(ctx, n,x,f,g,dn,wk2,&fn,&mxtake,stepmx,typx,wk1);
        if (iretcd < 0)
            return(-1);

        /* COMPARE TENSOR STEP TO NEWTON STEP */
        /* IF NEWTON STEP IS BETTER, SET NEXT ITERATE TO NEW NEWTON POINT */

        if (fn < fpls) {       
            fpls = fn;
            dcopy(ctx, n,dn,1,d,1);
            dcopy(ctx, n,wk2,1,xpls,1);
        }         
    }           

    /* current point in xpls[] */

    for (i = 1; i <= n; ++i)
        d[i] = xpls[i] - x[i];

    /* unscale xpls, and compute fpls and gpls */

    for (i = 1; i <= n; ++i)
        wk1[i] = xpls[i] * typx[i];

    /* new function value and gradient, if hflg then also new hessian */

    if (gflg)  
        err = ts_func(ctx, n,wk1,&fpls,1,gpls,hflg,h);
    else {
        f1[1] = fpls;
        err = fstofd(ctx, 1,1,n,wk1,f1,gpls,typx,rnf,wk2,1);
    }

    if (err)
        return(-1);
 
    /* CHECK STOPPING CONDITIONS */

    itrmcd = optstp(ctx, n,xpls,fpls,gpls,x,&icscmx,fscale,
                    mxtake,&rgx,&rsx,iretcd);

    ts_prn(ctx, ctx->Iter,fpls,rgx,rsx,n,x,gpls);
   
    if (itrmcd == 0)        /* not over yet */
        goto L500;
 
L350:
       
    /* TRANSFORM X BACK TO ORIGINAL SPACE */

    for (i = 1; i <= n; ++i)
        x[i] = xpls[i] * typx[i];
                   
    ctx->FMin = fpls;        /* final function value */

    /* update hessian if not already in h[] */

    if (hflg == 0) { 
        /*************** not used
        if (gflg)           
            fstofd(ctx, n,n,n,x,gpls,h,typx,rnf,wk2,3);                   
        else
        *****************/
        err = sndofd(ctx, n,n,x,fpls,h,typx,rnf,wk2,d);
        if (err)
            return(-1);
    }
    ts_prn1(ctx, n,h);

    return(itrmcd);

L500:   /* UPDATE INFORMATION AT THE CURRENT POINT */

    dcopy(ctx, n,xpls,1,x,1);

    for (i = 1; i <= n; ++i) 
        s[i] = -d[i];

    /* IF TOO MANY ITERATIONS THEN RETURN */

    if (ctx->Iter > ctx->MxIter) {           
        itrmcd = 4;
        goto L350;
    }
    fp = f;         /* update f */
    f = fpls;
    goto L200;      /* and perform next iteration */
}             

/* ------------------------------------------------------------------------ */
/*  lnsrch(n,x,f,g,p,xpls,fpls,mxtake,stepmx,typx,w2)                       */
/*                                                                          */
/*  FIND A NEXT NEWTON ITERATE BY LINE SEARCH.                              */
/*  THE ALPHA CONDITION ONLY LINE SEARCH                                    */
/*                                                                          */
/*  N            --> DIMENSION OF PROBLEM                                   */
/*  X(N)         --> OLD ITERATE:   X[K-1]                                  */
/*  F            --> FUNCTION VALUE AT OLD ITERATE, F(X)                    */
/*  G(N)         --> GRADIENT AT OLD ITERATE, G(X), OR APPROXIMATE          */
/*  P(N)         --> NON-ZERO NEWTON STEP                                   */
/*  XPLS(N)     <--  NEW ITERATE X[K]                                       */
/*  FPLS        <--  FUNCTION VALUE AT NEW ITERATE, F(XPLS)                 */
/*  MXTAKE      <--  BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED    */
/*  STEPMX       --> MAXIMUM ALLOWABLE STEP SIZE                            */
/*  TYPX(N)      --> DIAGONAL SCALING MATRIX FOR X (NOT IN UNCMIN)          */
/*  IPR          --> DEVICE TO WHICH TO SEND OUTPUT                         */
/*  W2           --> WORKING SPACE                                          */
/*                                                                          */
/*  Uses the following global variable:                                     */
/*  TOLSP   RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES                 */
/*          CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM                  */
/*                                                                          */
/*  Return: iretcd      = 0 if solution found                               */
/*                        1 if no new x vector found                        */
/*                       -1 if error in function evaluation                 */

int lnsrch(TDAContext *ctx, int n,double *x,double f,double *g,double *p,double *xpls, double *fpls,int *mxtake,double stepmx,double *typx, double *w2)
{
    register int i,k;
    int iretcd,err;
    double alpha,tmp,rln,sln,slp,scl,temp,temp1,temp2,almbmn,almbda,plmbda;
    double t1,t2,t3,a,b,disc,tlmbda,pfpls,ttmp;

    pfpls = plmbda = 0.0;
    *fpls = 0.0;
    *mxtake = 0;                                                       
    iretcd = 2;                                                       
    alpha = 1.e-4;

    tmp = 0.0;                                                        
    for (i = 1; i <= n; ++i)
        tmp += p[i] * p[i];

    sln = sqrt(tmp);        /* Newton length */                         

    if (sln <= stepmx)                                                 
        goto L10;

    /* NEWTON STEP LONGER THAN MAXIMUM ALLOWED */                            

    scl = stepmx / sln;                                                 
    dscal(ctx, n,scl,p,1);
    sln = stepmx;                                                     

L10:
    slp = ddot(ctx, n,g,1,p,1);
    rln = 0.0;                                                        
    for (i = 1; i <= n; ++i) {
        temp = 1.0;
        temp1 = fabs(x[i]);
        temp2 = dmax(ctx, temp1,temp);
        temp1 = fabs(p[i]);
        rln = dmax(ctx, rln,temp1 / temp2);
    }
    almbmn = ctx->TOLSP / rln;                                                
    almbda = 1.0;                                                       

    /* LOOP */                                                          
    /* CHECK IF NEW ITERATE SATISFACTORY.  GENERATE NEW LAMBDA IF NECESSARY. */
                                                                        
L100:
    if (iretcd < 2)                                                   
        return(iretcd);

    for (i = 1; i <= n; ++i)
        xpls[i] = x[i] + almbda * p[i];                                     

    for (k = 1; k <= n; ++k)
        w2[k] = xpls[k] * typx[k];

    err = ts_func(ctx, n,w2,fpls,0,&ttmp,0,&ttmp);
    if (err)
        return(-1); 

    if (*fpls > f + slp * alpha * almbda)                                   
        goto L130;

    /* IF(FPLS.LE. F+SLP*1.E-4*ALMBDA) THEN SOLUTION FOUND */            
                                                                        
    iretcd = 0;                                                     
    if (almbda == 1.0 && sln > 0.99 * stepmx)
        *mxtake = 1;        
    goto L100;                                                      
                                                                        
    /* SOLUTION NOT (YET) FOUND */                                           
                                                                        
L130:
    if (almbda >= almbmn)                                          
        goto L140;

    /* IF(ALMBDA .LT. ALMBMN) THEN */                                  
    /* NO SATISFACTORY XPLS FOUND SUFFICIENTLY DISTINCT FROM X */            
                                                                        
    iretcd = 1;                                                  
    goto L100;                                                    

    /* CALCULATE NEW LAMBDA */                                               
                                                                        
L140:
    if (almbda != 1.0)                                             
        goto L150;

    /* IF(ALMBDA.EQ.1.0) THEN FIRST BACKTRACK: QUADRATIC FIT */                                     
                                                                        
    tlmbda = -slp / (2.0 * (*fpls - f - slp));                             
    goto L170;                                                  

    /* ELSE ALL SUBSEQUENT BACKTRACKS: CUBIC FIT */                               
                                                                        
L150:
    t1 = *fpls - f - almbda * slp;                                       
    t2 = pfpls - f - plmbda * slp;                                      
    t3 = 1.0 / (almbda - plmbda);                                     
    a = t3 * (t1 / (almbda * almbda) - t2 / (plmbda * plmbda));             
    b = t3 * (t2 * almbda / (plmbda * plmbda) -                           
                       t1 * plmbda / (almbda * almbda));                         
    disc = b * b - 3.0 * a * slp;                                         
    if (disc <= b * b)                                            
        goto L160;

    /* IF(DISC.GT. B*B) THEN */                                    
    /* ONLY ONE POSITIVE CRITICAL POINT, MUST BE MINIMUM */                  
                                                                        
    tlmbda = (-b + fsign(ctx, 1.0,a) * sqrt(disc)) / (3.0 * a);             
    goto L165;                                                

    /* ELSE BOTH CRITICAL POINTS POSITIVE, FIRST IS MINIMUM */                    
                                                                        
L160:
    tlmbda = (-b - fsign(ctx, 1.0,a) * sqrt(disc)) / (3.0 * a);             

L165:
    if (tlmbda > 0.5 * almbda)
        tlmbda = 0.5 * almbda;              

L170:
    plmbda = almbda;                                                
    pfpls = *fpls;                                                   
    if (tlmbda >= almbda * 0.1)                                    
        goto L180;
    almbda *= 0.1;                                                
    goto L190;                                                  
                                                                        
L180:
    almbda = tlmbda;                                              
                                                                        
L190:
    goto L100;                                                        
}

/* ------------------------------------------------------------------------ */
/*  zhz(nc,n,y,h,u,t)                                                       */
/*                                                                          */  
/*  COMPUTE QTHQ(N,N) AND ZTHZ(N-1,N-1) = FIRST N-1 ROWS AND                */
/*  FIRST N-1 COLUMNS OF QTHQ                                               */
/*                                                                          */
/*  NC      --> Column DIMENSION OF MATRIX                                  */
/*  N       --> DIMENSION OF PROBLEM                                        */
/*  Y(N)    --> FIRST BASIS IN Q                                            */
/*  H(N,N) <--> ON INPUT : HESSIAN                                          */
/*              ON OUTPUT: QTHQ (ZTHZ)                                      */
/*  U(N)   <--  VECTOR TO FORM Q AND Z                                      */
/*  T(N)    --> WORKSPACE                                                   */
/*                                                                          */  

void zhz(TDAContext *ctx, int nc,int n,double *y,double *h,double *u,double *t)
{
    register int i,j;
    double ynorm,d,s,sgn;

    /* U=Y+SGN(Y(N))||Y||E(N) */

    if (y[n] != 0.0)       
        sgn = y[n] / fabs(y[n]);
    else   
        sgn = 1.0;
                
    ynorm = ddot(ctx, n,y,1,y,1);
    ynorm = sqrt(ynorm);
    u[n] = y[n] + sgn * ynorm;
    dcopy(ctx, n - 1,y,1,u,1);
 
    /* D=UTU/2 */
    d = ddot(ctx, n,u,1,u,1);
    d /= 2.0;
 
    /* T=2HU/UTU */

    for (i = 1; i <= n; ++i) {
        t[i] = 0.0;
        for (j = 1; j <= n; ++j)
            t[i] = h[(i - 1) * nc + j] * u[j] + t[i];
        t[i] /= d;  
    }

    /* S=4UHU/(UTU)**2 */

    s = ddot(ctx, n,u,1,t,1);
    s /= d; 
 
    /* COMPUTE QTHQ (ZTHZ) */

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j) {
            h[(i - 1) * nc + j] = 
            h[(i - 1) * nc + j] - u[i] * t[j] - t[i] * u[j] + u[i] * u[j] * s;
        }
    }
}           

/* ------------------------------------------------------------------------ */
/*  solvew(nc,n,al,u,w,b)                                                   */
/*                                                                          */
/*  SOLVE L*W=ZT*V                                                          */
/*                                                                          */
/*  NC            --> Column DIMENSION OF MATRIX                            */
/*  N             --> DIMENSION OF PROBLEM                                  */
/*  AL(N-1,N-1)   --> LOWER TRIAGULAR MATRIX                                */
/*  U(N)          --> VECTOR TO FORM Z                                      */
/*  W(N)          --> ON INPUT : VECTOR V IN SYSTEM OF LINEAR EQUATIONS     */
/*                    ON OUTPUT: SOLUTION OF SYSTEM OF LINEAR EQUATIONS     */
/*  B(N)          --> WORKSPACE TO STORE ZT*V                               */
/*                                                                          */

void solvew(TDAContext *ctx, int nc,int n,double *al,double *u,double *w,double *b)
{
    register int i,j;
    double d;

    /* FORM ZT*V (STORED IN B) */

    d = ddot(ctx, n,u,1,u,1);
    d /= 2.0; 

    for (i = 1; i < n; ++i) {
        b[i] = 0.0;
        for (j = 1; j <= n; ++j)
            b[i] += u[j] * u[i] * w[j] / d;
        b[i] = w[i] - b[i];
    }
 
    /* SOLVE LW=ZT*V */

    forslv(ctx, nc,n - 1,al,w,b);
}

/* ------------------------------------------------------------------------ */
/*  dstar(nc,n,u,s,w1,w2,w3,sigma,al,d)                                     */
/*                                                                          */
/*  COMPUTE TENSOR STEP D=SIGMA*S+ZT*T(SIGMA)                               */
/*                                                                          */
/*  NC          --> Column DIMENSION OF MATRIX                              */
/*  N           --> DIMENSION OF PROBLEM                                    */
/*  U(N)        --> VECTOR TO FORM Z                                        */
/*  S(N)        --> PREVIOUS STEP                                           */
/*  W1(N)       --> L**-1*ZT*A, WHERE A IS DESCRIBED IN SUBROUTINE SLVMDL   */
/*  W2(N)       --> L**-1*ZT*SH, WHERE H IS CURRENT HESSIAN                 */
/*  W3(N)       --> L**-1*ZT*G, WHERE G IS CURRENT GRADIENT                 */
/*  SIGMA       --> SOLUTION FOR REDUCED ONE VARIABLE MODEL                 */
/*  AL(N-1,N-1) --> LOWER TRIANGULAR MATRIX L                               */
/*  D(N)        --> TENSOR STEP                                             */
/*                                                                          */
      
void dstar(TDAContext *ctx, int nc,int n,double *u,double *s,double *w1,double *w2,double *w3, double sigma,double *al,double *d)
{
    register int i;
    double utu,utt,temp;

    if (n == 1)         
        d[1] = sigma * s[1];

    else {
 
        /* COMPUTE T(SIGMA)=-(ZTHZ)*ZT*(G+SIGMA*SH+SIGMA**2*A/2) (STORED IN D) */

        for (i = 1; i < n; ++i)  
            w2[i] = w3[i] + sigma * w2[i] + 0.5 * w1[i] * sigma * sigma;
	
        bakslv(ctx, nc,n - 1,al,d,w2);
	
        d[n] = 0.0; 
 
        /* COMPUTE TENSOR STEP D=SIGMA*S+ZT*T(SIGMA) */

        utu = ddot(ctx, n,u,1,u,1);
        utt = ddot(ctx, n,u,1,d,1);
        temp = utt / utu;
        for (i = 1; i <= n; ++i)
            d[i] = sigma * s[i] - (d[i] - 2.0 * u[i] * temp);
    }
}             

/* ------------------------------------------------------------------------ */
/*  mkmdl(nc,n,f,fp,g,gp,s,h,alpha,beta,sh,a)                               */
/*                                                                          */
/*  FORM TENSOR MODEL                                                       */
/*                                                                          */
/*  NC      --> Column DIMENSION OF MATRIX                                  */
/*  N       --> DIMENSION OF PROBLEM                                        */
/*  F       --> CURRENT FUNCTION VALUE                                      */
/*  FP      --> PREVIOUS FUNCTION VALUE                                     */
/*  G(N)    --> CURRENT GRADIENT                                            */
/*  GP(N)   --> PREVIOUS GRADIENT                                           */
/*  S(N)    --> STEP TO PREVIOUS POINT                                      */
/*  H(N,N)  --> HESSIAN                                                     */
/*  ALPHA  <--  SCALAR TO FORM 3RD ORDER TERM OF TENSOR MODEL               */
/*  BETA   <--  SCALAR TO FORM 4TH ORDER TERM OF TENSOR MODEL               */
/*  SH(N)  <--  SH                                                          */
/*  A(N)   <--  A=2*(GP-G-SH-S*BETA/(6*STS))                                */
/*                                                                          */
      
void mkmdl(TDAContext *ctx, int nc,int n,double f,double fp,double *g,double *gp,double *s, double *h,double *alpha,double *beta,double *sh,double *a)
{
    register int i,j;
    double gs,gps,shs,b1,b2,sts;

    /* COMPUTE SH */

    for (i = 1; i <= n; ++i) {
        sh[i] = 0.0;
        for (j = 1; j <= n; ++j)
            sh[i] += s[j] * h[(j - 1) * nc + i];
    }
    gs = ddot(ctx, n,g,1,s,1);
    gps = ddot(ctx, n,gp,1,s,1);
    shs = ddot(ctx, n,sh,1,s,1);
    b1 = gps - gs - shs;
    b2 = fp - f - gs - 0.5 * shs;
    *alpha = 24.0 * b2 - 6.0 * b1;
    *beta = 24.0 * b1 - 72.0 * b2;
 
    /* COMPUTE A */

    sts = ddot(ctx, n,s,1,s,1);
    for (i = 1; i <= n; ++i) 
        a[i] = 2.0 * (gp[i] - g[i] - sh[i] - s[i] * *beta / (6.0 * sts));
}

/* ------------------------------------------------------------------------ */
/*  sigma(sgstar,a,b,c,d)                                                   */
/*                                                                          */
/*  COMPUTE DESIRABLE ROOT OF REDUCED ONE VARIABLE EQUATION                 */
/*                                                                          */
/*  SGSTAR <--> DESIRABLE ROOT                                              */
/*  A       --> COEFFICIENT OF 3RD ORDER TERM                               */
/*  B       --> COEFFICIENT OF 2ND ORDER TERM                               */
/*  C       --> COEFFICIENT OF 1ST ORDER TERM                               */
/*  D       --> COEFFICIENT OF CONSTANT TERM                                */
/*                                                                          */

void sigma(TDAContext *ctx, double *sgstar,double a,double b,double c,double d)
{
    double s1,s2,s3;
 
    /* COMPUTE ALL THREE ROOTS */

    roots(ctx, &s1,&s2,&s3,a,b,c,d);
 
    /* SORT ROOTS */

    sortrt(ctx, &s1,&s2,&s3);
 
    /* CHOOSE DESIRABLE ROOT */

    if (a > 0.0) {       
        *sgstar = s3;
        if (s2 >= 0.0)         
            *sgstar = s1;
    }       
    else {
        *sgstar = s2;
        if ((s1 > 0.0) || (s3 < 0.0)) {      
            if (s1 > 0.0)      
                *sgstar = s1;
            else 
                *sgstar = s3;
            a = 0.0;
        }         
    }
}             

/* ------------------------------------------------------------------------ */
/*  roots(s1,s2,s3,a,b,c,d)                                                 */
/*                                                                          */
/*  COMPUTE ROOT(S) OF 3RD ORDER EQUATION                                   */
/*                                                                          */
/*  S1     <--  ROOT   (IF THREE ROOTS ARE                                  */
/*  S2     <--  ROOT    EQUAL, THEN S1=S2=S3)                               */
/*  S3     <--  ROOT                                                        */
/*  A       --> COEFFICIENT OF 3RD ORDER TERM                               */
/*  B       --> COEFFICIENT OF 2ND ORDER TERM                               */
/*  C       --> COEFFICIENT OF 1ST ORDER TERM                               */
/*  D       --> COEFFICIENT OF CONSTANT TERM                                */

void roots(TDAContext *ctx, double *s1,double *s2,double *s3,double a,double b,double c,double d)
{
    (void)ctx;        /* unused: the signature is shared */
    double a1,a2,a3,q,r,v,s,t,temp,theta;
    double pi = 3.141592653589793;

    a1 = b / a;
    a2 = c / a;
    a3 = d / a;
    q = (30.0 * a2 - a1 * a1) / 90.0;
    r = (90.0 * a1 * a2 - 27.0 * a3 - 2.0 * a1 * a1 * a1) / 54.0;  
    v= q * q * q + r * r;
    if (v > 0.0) {       
        s = r + sqrt(v);
        t = r - sqrt(v);
        if (t < 0.0)       
            t = -pow(-t,1.0 / 3.0);
        else
            t = pow(t,1.0 / 3.0);

        if (s < 0.0)         
            s = -pow(-s,1.0 / 3.0);
        else
            s = pow(s,1.0 / 3.0);

        *s1 = s + t - a1 / 3.0;
        *s3 = *s1;
        *s2 = *s1;
    }
    else {    
        temp = r / sqrt(-q * q * q);

        theta = acos(temp);
        theta /= 3.0;         
        temp = 2.0 * sqrt(-q);
        *s1 = temp * cos(theta) - a1 / 3.0; 
        *s2 = temp * cos(theta + pi * 20.0 / 30.0) - a1 / 30.0;
        *s3 = temp * cos(theta + pi * 40.0 / 30.0) - a1 / 30.0;
    }         
}             

/* ------------------------------------------------------------------------ */
/*  sortrt(s1,s2,s3)                                                        */
/*                                                                          */
/*  SORT ROOTS INTO ASCENDING ORDER                                         */
/*                                                                          */
/*  S1  <--> ROOT                                                           */
/*  S2  <--> ROOT                                                           */
/*  S3  <--> ROOT                                                           */
 
void sortrt(TDAContext *ctx, double *s1,double *s2,double *s3)
{
    (void)ctx;        /* unused: the signature is shared */
    double t;

    if (*s1 > *s2) {     
        t = *s1;
        *s1 = *s2;
        *s2 = t;
    }
    if (*s2 > *s3) {   
        t = *s2;
        *s2 = *s3;
        *s3 = t;
    }
    if (*s1 > *s2) {   
        t = *s1;
        *s1 = *s2;
        *s2 = t;
    }
}

/* ------------------------------------------------------------------------ */
/*  fstofd(nc,m,n,xpls,fpls,a,typx,rnoise,fhat,icase)                       */
/*                                                                          */
/*  FIND FIRST ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A" TO THE     */
/*  FIRST DERIVATIVE OF THE FUNCTION DEFINED BY THE SUBPROGRAM "FNAME"      */
/*  EVALUATED AT THE NEW ITERATE "XPLS".                                    */
/*                                                                          */
/*  FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE:                          */
/*  1) THE FIRST DERIVATIVE (GRADIENT) OF THE OPTIMIZATION FUNCTION "FCN"   */
/*                                                                          */
/*  2) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION         */
/*     IF NO ANALYTIC USER ROUTINE HAS BEEN SUPPLIED FOR THE HESSIAN BUT    */
/*     ONE HAS BEEN SUPPLIED FOR THE GRADIENT ("FCN") AND IF THE            */
/*     OPTIMIZATION FUNCTION IS INEXPENSIVE TO EVALUATE                     */
/*                                                                          */
/*  NOTE                                                                    */
/*  _M=1 (OPTIMIZATION) ALGORITHM ESTIMATES THE GRADIENT OF THE FUNCTION    */
/*       (FCN).   FCN(X) # F: R(N)-->R(1)                                   */
/*  _M=N (SYSTEMS) ALGORITHM ESTIMATES THE JACOBIAN OF THE FUNCTION         */
/*       FCN(X) # F: R(N)-->R(N).                                           */
/*  _M=N (OPTIMIZATION) ALGORITHM ESTIMATES THE HESSIAN OF THE OPTIMIZATIO  */
/*       FUNCTION, WHERE THE HESSIAN IS THE FIRST DERIVATIVE OF "FCN"       */
/*                                                                          */
/*  NC           --> Column DIMENSION OF MATRIX                             */
/*  M            --> NUMBER OF ROWS IN A                                    */
/*  N            --> NUMBER OF COLUMNS IN A; DIMENSION OF PROBLEM           */
/*  XPLS(N)      --> NEW ITERATE:  X[K]                                     */
/*  FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION                */
/*  FPLS(M)      --> _M=1 (OPTIMIZATION) FUNCTION VALUE AT NEW ITERATE:     */
/*                        FCN(XPLS)                                         */
/*                   _M=N (OPTIMIZATION) VALUE OF FIRST DERIVATIVE          */
/*                        (GRADIENT) GIVEN BY USER FUNCTION FCN             */
/*                   _M=N (SYSTEMS)  FUNCTION VALUE OF ASSOCIATED           */
/*                        MINIMIZATION FUNCTION                             */
/*  A(N,NC)     <--  FINITE DIFFERENCE APPROXIMATION (SEE NOTE).  ONLY      */
/*                   LOWER TRIANGULAR MATRIX AND DIAGONAL ARE RETURNED      */
/*  RNOISE       --> RELATIVE NOISE IN FCN [F(X)]                           */
/*  FHAT(M)      --> WORKSPACE                                              */
/*  ICASE        --> =1 OPTIMIZATION (GRADIENT)                             */
/*                   =2 SYSTEMS                                             */
/*                   =3 OPTIMIZATION (HESSIAN)                              */
/*  INTERNAL VARIABLES                                                      */
/*  STEPSZ - STEPSIZE IN THE J-TH VARIABLE DIRECTION                        */
/*                                                                          */
/*  Note: original code had a problem in this function concerning           */
/*  fhat array. This has been changed. However, this function is            */
/*  currently never used for calculation of hessian.                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error in function evaluation                      */
/*                                                                          */

int fstofd(TDAContext *ctx, int nc,int m,int n,double *xpls,double *fpls,double *a,double *typx, double rnoise,double *fhat,int icase)
{    
    (void)fhat;        /* unused: the signature is shared */
    register int i,j;
    int err;
    double xtmpj,stepsz,tmp,f;

    /* FIND J-TH COLUMN OF A */                                              
    /* EACH COLUMN IS DERIVATIVE OF F(FCN) WITH RESPECT TO XPLS(J) */        
                                                                        
    for (j = 1; j <= n; ++j) {
        xtmpj = xpls[j];                                                  
        stepsz = sqrt(rnoise) * dmax(ctx, fabs(xpls[j]),1.0);                       
        xpls[j] = xtmpj + stepsz;                                           

        err = ts_func(ctx, n,xpls,&f,0,&tmp,0,&tmp);
        if (err)
            return(-1);

        xpls[j] = xtmpj;                                                  
        for (i = 1; i <= m; ++i) {
            a[(i - 1) * nc + j] = (f - fpls[i]) / stepsz;                              
            a[(i - 1) * nc + j] *= typx[j];
        }
    }
    if (icase != 3)                                                    
        return(0);

    /* IF COMPUTING HESSIAN, A MUST BE SYMMETRIC */                          
                                                                        
    if (n == 1)                                                        
        return(0);

    for (j = 1; j < n; ++j) {
        for (i = j + 1; i <= m; ++i)
            a[(i - 1) * nc + j] =
                (a[(i - 1) * nc + j] + a[(j - 1) * nc + i]) / 2.0;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sndofd(nc,n,xpls,fpls,a,typx,rnoise,stepsz,anbr)                        */
/*                                                                          */
/*  FIND SECOND ORDER FORWARD FINITE DIFFERENCE APPROXIMATION "A"           */
/*  TO THE SECOND DERIVATIVE (HESSIAN) OF THE FUNCTION DEFINED BY THE SUBP  */
/*  "FCN" EVALUATED AT THE NEW ITERATE "XPLS"                               */
/*                                                                          */
/*  FOR OPTIMIZATION USE THIS ROUTINE TO ESTIMATE                           */
/*  1) THE SECOND DERIVATIVE (HESSIAN) OF THE OPTIMIZATION FUNCTION         */
/*     IF NO ANALYTICAL USER FUNCTION HAS BEEN SUPPLIED FOR EITHER          */
/*     THE GRADIENT OR THE HESSIAN AND IF THE OPTIMIZATION FUNCTION         */
/*     "FCN" IS INEXPENSIVE TO EVALUATE.                                    */ 
/*                                                                          */
/*  NC           --> Column DIMENSION OF MATRIX                             */
/*  N            --> DIMENSION OF PROBLEM                                   */
/*  XPLS(N)      --> NEW ITERATE:   X[K]                                    */
/*  FCN          --> NAME OF SUBROUTINE TO EVALUATE FUNCTION                */
/*  FPLS         --> FUNCTION VALUE AT NEW ITERATE, F(XPLS)                 */
/*  A(N,N)      <--  FINITE DIFFERENCE APPROXIMATION TO HESSIAN             */
/*                   ONLY LOWER TRIANGULAR MATRIX AND DIAGONAL              */
/*                   ARE RETURNED                                           */
/*  RNOISE       --> RELATIVE NOISE IN FNAME [F(X)]                         */
/*  STEPSZ(N)    --> WORKSPACE (STEPSIZE IN I-TH COMPONENT DIRECTION)       */
/*  ANBR(N)      --> WORKSPACE (NEIGHBOR IN I-TH DIRECTION)                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error in function evaluation.                     */
/*                                                                          */

int sndofd(TDAContext *ctx, int nc,int n,double *xpls,double fpls,double *a,double *typx, double rnoise,double *stepsz,double *anbr)
{    
    register int i,j;
    int err;
    double ov3,xtmpi,xtmpj,tmp,fhat;
                                                                        
    /* FIND I-TH STEPSIZE AND EVALUATE NEIGHBOR IN DIRECTION */              
    /* OF I-TH UNIT VECTOR. */                                               
                                                                        
    ov3 = 1.0 / 3.0;
    for (i = 1; i <= n; ++i) {
        xtmpi = xpls[i];                                                  
        stepsz[i] = pow(rnoise,ov3) * dmax(ctx, fabs(xpls[i]),1.0);                  
        xpls[i] = xtmpi + stepsz[i];                                        

        err = ts_func(ctx, n,xpls,anbr + i,0,&tmp,0,&tmp);
        if (err)
            return(-1);

        xpls[i] = xtmpi;                                                  
    }
                                                                        
    /* CALCULATE COLUMN I OF A */                                            
                                                                        
    for (i = 1; i <= n; ++i) {
        xtmpi = xpls[i];                                                  
        xpls[i] = xtmpi + 2.0 * stepsz[i];                                    

        err = ts_func(ctx, n,xpls,&fhat,0,&tmp,0,&tmp);
        if (err)
            return(-1);

        a[(i - 1) * nc + i] =
             ((fpls - anbr[i]) + (fhat - anbr[i])) / (stepsz[i] * stepsz[i]);   
        a[(i - 1) * nc + i] *= typx[i] * typx[i];

        /* CALCULATE SUB-DIAGONAL ELEMENTS OF COLUMN */                          

        if (i == n)                                                      
            goto L25;

        xpls[i] = xtmpi + stepsz[i];                                        

        for (j = i + 1; j <= n; ++j) {
            xtmpj = xpls[j];                                                
            xpls[j] = xtmpj + stepsz[j];                                      

            err = ts_func(ctx, n,xpls,&fhat,0,&tmp,0,&tmp);
            if (err)
                return(-1);

            a[(j - 1) * nc + i] =
                ((fpls - anbr[i]) + (fhat - anbr[j])) / (stepsz[i] * stepsz[j]); 

            a[(j - 1) * nc + i] *= typx[i] * typx[j];
            xpls[j] = xtmpj;                                                
        }
L25:  
        xpls[i] = xtmpi;                                                  
    }
    return(0);
}                                                                       

/* ------------------------------------------------------------------------ */
/*  bakslv(nc,n,a,x,b)                                                      */
/*                                                                          */
/*  SOLVE  AX=B  WHERE A IS UPPER TRIANGULAR MATRIX.                        */
/*  NOTE THAT A IS INPUT AS A LOWER TRIANGULAR MATRIX AND                   */
/*  THAT THIS ROUTINE TAKES ITS TRANSPOSE IMPLICITLY.                       */
/*                                                                          */
/*  NC           --> Column DIMENSION OF MATRIX                             */
/*  N            --> DIMENSION OF PROBLEM                                   */
/*  A(N,N)       --> LOWER TRIANGULAR MATRIX (PRESERVED)                    */
/*  X(N)        <--  SOLUTION VECTOR                                        */
/*  B(N)         --> RIGHT-HAND SIDE VECTOR                                 */
/*                                                                          */
/*  NOTE                                                                    */
/*  IF B IS NO LONGER REQUIRED BY CALLING ROUTINE,                          */
/*  THEN VECTORS B AND X MAY SHARE THE SAME STORAGE.                        */
/*                                                                          */

void bakslv(TDAContext *ctx, int nc,int n,double *a,double *x,double *b)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j,ip1;
    double sum;
                                                                        
    /* SOLVE (L-TRANSPOSE)X=B. (BACK SOLVE) */                               
                                                                        
    i = n;                                                             
    x[i] = b[i] / a[(i - 1) * nc + i];                                     
    if (n == 1)                                                        
        return;
L30:
    ip1 = i;                                                          
    i--;                                                              
    sum = 0.0;                                                         
    for (j = ip1; j <= n; ++j)
        sum += a[(j - 1) * nc + i] * x[j];                                            

    x[i] = (b[i] - sum) / a[(i - 1) * nc + i];                               
    if (i > 1)
        goto L30;
}                                                                       

/* ------------------------------------------------------------------------ */
/*  forslv(nc,n,a,x,b)                                                      */
/*                                                                          */
/*  SOLVE  AX=B  WHERE A  IS LOWER TRIANGULAR  MATRIX                       */
/*                                                                          */
/*  NC          ---> Column DIMENSION OF MATRIX                             */
/*  N           ---> DIMENSION OF PROBLEM                                   */
/*  A(N,N)      ---> LOWER TRIANGULAR MATRIX (PRESERVED)                    */
/*  X(N)       <---  SOLUTION VECTOR                                        */
/*  B(N)        ---> RIGHT-HAND SIDE VECTOR                                 */
/*                                                                          */
/*  NOTE                                                                    */
/*  THEN VECTORS B AND X MAY SHARE THE SAME STORAGE                         */
/*                                                                          */

void forslv(TDAContext *ctx, int nc,int n,double *a,double *x,double *b)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,j;
    double sum;
                                                                        
    /* SOLVE LX=B.  (FOREWARD  SOLVE) */                                     
                                                                        
    x[1] = b[1] / a[1];                                                   
    if (n == 1)                                                        
        return;

    for (i = 2; i <= n; ++i) {
        sum = 0.0;                                                        
        for (j = 1; j < i; ++j)
            sum += a[(i - 1) * nc + j] * x[j];                                          
        x[i] = (b[i] - sum) / a[(i - 1) * nc + i];                                         
    }
}

/* ------------------------------------------------------------------------ */
/*  choldr(nc,n,h,g,pivot,e,diag,addmax)                                    */
/*                                                                          */
/*  DRIVER FOR CHOLESKY DECOMPOSITION                                       */
/*                                                                          */
/*  NC          --> number of columns in h                                  */
/*  N           --> DIMENSION OF PROBLEM                                    */
/*  H(N,N)      --> MATRIX                                                  */
/*  G(N)        --> WORK SPACE                                              */
/*  EPS         --> MACHINE EPSILON                                         */
/*  PIVOT(N)    --> PIVOTING VECTOR                                         */
/*  E(N)        --> DIAGONAL MATRIX ADDED TO H FOR MAKING H P.D.            */
/*  DIAG(N)     --> DIAGONAL OF H                                           */
/*  ADDMAX      --> ADDMAX * I  IS ADDED TO H                               */
/*                                                                          */

void choldr(TDAContext *ctx, int nc,int n,double *h,double *g,int *pivot,double *e, double *diag,double addmax)
{
    register int i,j,k;
    int redo;
    double tau1,tau2,temp;

    redo = 0;      
 
    /*  SAVE DIAGONAL OF H  */

    for (i = 1; i <= n; ++i) 
        diag[i] = h[(i - 1) * nc + i];
 
    tau2 = tau1= pow(ctx->EPSI,1.0 / 3.0);

    modchl(ctx, nc,n,h,g,tau1,tau2,pivot,e);

    addmax = e[n];
    for (i = 1; i <= n; ++i) {
        if (pivot[i] != i)
            redo = 1;
    }
    if ((addmax > 0.0) || redo) {       

        /* H IS NOT P.D. */           
        /* H=H+UI */

        for (i = 2; i <= n; ++i) {
            for (j = 1; j < i; ++j)
                h[(i - 1) * nc + j] = h[(j - 1) * nc + i];
        }
        for (i = 1; i <= n; ++i) {
            pivot[i] = i;
            h[(i - 1) * nc + i] = diag[i] + addmax;
        }

        /* COMPUTE L */               

        for (j = 1; j <= n; ++j) {
 
            /* COMPUTE L(J,J) */

            temp = 0.0;
            if (j >  1) {    
                for (i = 1; i < j; ++i)
                    temp += h[(j - 1) * nc + i] * h[(j - 1) * nc + i];
            }
            h[(j - 1) * nc + j] -= temp;
            h[(j - 1) * nc + j] = sqrt(h[(j - 1) * nc + j]);
 
            /* COMPUTE L(I,J) */

            for (i = j + 1; i <= n; ++i) {
                temp = 0.0;
                if (j > 1) {       
                    for (k = 1; k < j; ++k)
                        temp += h[(i - 1) * nc + k] * h[(j - 1) * nc + k];
                }  
                h[(i - 1) * nc + j] = h[(j - 1) * nc + i] - temp;
                h[(i - 1) * nc + j] /= h[(j - 1) * nc + j];
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  modchl(ndmim,n,a,g,tau1,tau2,p,e)                                       */
/*                                                                          */
/*  PERFORM A MODIFIED CHOLESKY FACTORIZATION                               */
/*  OF THE FORM (PTRANSPOSE)AP  + E = L(LTRANSPOSE),                        */
/*  WHERE L IS STORED IN THE LOWER TRIANGLE OF THE                          */
/*  ORIGINAL MATRIX A.                                                      */
/*  THE FACTORIZATION HAS 2 PHASES:                                         */
/*  PHASE 1: PIVOT ON THE MAXIMUM DIAGONAL ELEMENT.                         */
/*           CHECK THAT THE NORMAL CHOLESKY UPDATE                          */
/*           WOULD RESULT IN A POSITIVE DIAGONAL                            */
/*           AT THE CURRENT ITERATION, AND                                  */
/*           IF SO, DO THE NORMAL CHOLESKY UPDATE,                          */
/*           OTHERWISE SWITCH TO PHASE 2.                                   */
/*  PHASE 2: PIVOT ON THE MINIMUM OF THE NEGATIVES                          */
/*           OF THE LOWER GERSCHGORIN BOUND                                 */
/*           ESTIMATES.                                                     */
/*           COMPUTE THE AMOUNT TO ADD TO THE                               */
/*           PIVOT ELEMENT AND ADD THIS                                     */
/*           TO THE PIVOT ELEMENT.                                          */
/*           DO THE CHOLESKY UPDATE.                                        */
/*           UPDATE THE ESTIMATES OF THE                                    */
/*           GERSCHGORIN BOUNDS.                                            */
/*                                                                          */
/*  INPUT   : NDIM    - LARGEST DIMENSION OF MATRIX THAT                    */
/*                      WILL BE USED (number of colums in a[] )             */
/*            N       - DIMENSION OF MATRIX A                               */
/*            A       - N*N SYMMETRIC MATRIX (ONLY LOWER TRIANGULAR         */
/*                      PORTION OF A, INCLUDING THE MAIN DIAGONAL, IS USED) */
/*            G       - N*1 WORK ARRAY                                      */
/*            MCHEPS  - MACHINE PRECISION                                   */
/*            TAU1    - TOLERANCE USED FOR DETERMINING WHEN TO SWITCH TO    */
/*                      PHASE 2                                             */
/*            TAU2    - TOLERANCE USED FOR DETERMINING THE MAXIMUM          */
/*                      CONDITION NUMBER OF THE FINAL 2X2 SUBMATRIX.        */
/*                                                                          */
/*  OUTPUT  : L     - STORED IN THE MATRIX A (IN LOWER TRIANGULAR           */
/*                    PORTION OF A, INCLUDING THE MAIN DIAGONAL)            */
/*            P     - A RECORD OF HOW THE ROWS AND COLUMNS                  */
/*                    OF THE MATRIX WERE PERMUTED WHILE                     */
/*                    PERFORMING THE DECOMPOSITION                          */
/*            E     - N*1 ARRAY, THE ITH ELEMENT IS THE                     */
/*                    AMOUNT ADDED TO THE DIAGONAL OF A                     */
/*                    AT THE ITH ITERATION                                  */
/*                                                                          */
/*  J              - CURRENT ITERATION NUMBER                               */
/*  IMING          - INDEX OF THE ROW WITH THE MIN. OF THE                  */
/*                   NEG. LOWER GERSCH. BOUNDS                              */
/*  IMAXD          - INDEX OF THE ROW WITH THE MAXIMUM DIAG.                */
/*                   ELEMENT                                                */
/*  I,ITEMP,JPL,K  - TEMPORARY INTEGER VARIABLES                            */
/*  DELTA          - AMOUNT TO ADD TO AJJ AT THE JTH ITERATION              */
/*  GAMMA          - THE MAXIMUM DIAGONAL ELEMENT OF THE ORIGINAL           */
/*                   MATRIX A.                                              */
/*  NORMJ          - THE 1 NORM OF A(COLJ), ROWS J+1 --> N.                 */
/*  MING           - THE MINIMUM OF THE NEG. LOWER GERSCH. BOUNDS           */
/*  MAXD           - THE MAXIMUM DIAGONAL ELEMENT                           */
/*  TAUGAM         - TAU1 * GAMMA                                           */
/*  PHASE1          - LOGICAL, TRUE IF IN PHASE1, OTHERWISE FALSE           */
/*  DELTA1,TEMP,JDMIN,TDMIN,TEMPJJ - TEMPORARY DOUBLE PRECISION VARS.       */
/*                                                                          */

void modchl(TDAContext *ctx, int ndim,int n,double *a,double *g,double tau1, double tau2,int *p,double *e)
{
    register int i,j,k,jp1; 
    int imaxd,itemp,phase1,iming;
    double delta,maxd,temp,tempjj,ming,normj,delta1,tmp,taugam,gamma;
    double jdmin,tdmin;

    ts_init(ctx, n,ndim,a,&phase1,&delta,p,g,e,&ming,tau1,&gamma,&taugam);

    /* CHECK FOR N=1 */

    if (n == 1) {         
        tmp = a[(1 - 1) * ndim + 1];
        delta = (tau2 * fabs(tmp)) - tmp;    
        if (delta > 0.0)
            e[1] = delta;
        if (tmp == 0.0)
            e[1] = tau2;
        a[(1 - 1) * n + 1] = sqrt(tmp + e[1]);
    }         
 
    for (j = 1; j < n; ++j) {
 
        /* PHASE 1 */
 
        if (phase1) {    
          
            /* FIND INDEX OF MAXIMUM DIAGONAL ELEMENT A(I,I) WHERE I>=J */
 
            maxd = a[(j - 1) * n + j];
            imaxd = j;
            for (i = j + 1; i <= n; ++i) {
                tmp = a[(i - 1) * ndim + i];
                if (maxd < tmp) {       
                    maxd = tmp;    
                    imaxd = i;
                }                   
            }

            /* PIVOT TO THE TOP THE ROW AND COLUMN WITH THE MAX DIAG */
 
            if (imaxd != j) {         
 
                /* SWAP ROW J WITH ROW OF MAX DIAG */
 
                for (i = 1; i < j; ++i) {
                    temp = a[(j - 1) * ndim + i];
                    a[(j - 1) * ndim + i] = a[(imaxd - 1) * ndim + i];
                    a[(imaxd - 1) * ndim + i] = temp;
                }
 
                /* SWAP COLJ AND ROW MAXDIAG BETWEEN J AND MAXDIAG */
 
                for (i = j + 1; i < imaxd; ++i) {
                    temp = a[(i - 1) * ndim + j];
                    a[(i - 1) * ndim + j] = a[(imaxd - 1) * ndim + i];
                    a[(imaxd - 1) * ndim + i] = temp;
                }
 
                /* SWAP COLUMN J WITH COLUMN OF MAX DIAG */
 
                for (i = imaxd + 1; i <= n; ++i) {
                    temp = a[(i - 1) * ndim + j];
                    a[(i - 1) * ndim + j] = a[(i - 1) * ndim + imaxd];
                    a[(i - 1) * ndim + imaxd] = temp;
                }
 
                /* SWAP DIAG ELEMENTS */
         
                temp = a[(j - 1) * ndim + j];
                a[(j - 1) * ndim + j] = a[(imaxd - 1) * ndim + imaxd];
                a[(imaxd - 1) * ndim + imaxd] = temp;
 
                /* SWAP ELEMENTS OF THE PERMUTATION VECTOR */
 
                itemp = p[j];
                p[j] = p[imaxd];
                p[imaxd] = itemp;
            }       

            /* CHECK TO SEE WHETHER THE NORMAL CHOLESKY UPDATE FOR THIS */
            /* ITERATION WOULD RESULT IN A POSITIVE DIAGONAL,  */
            /* AND IF NOT THEN SWITCH TO PHASE 2. */

            jp1 = j + 1;
            tempjj = a[(j - 1) * ndim + j];

            if (tempjj > 0.0) {         

                jdmin = a[(jp1 - 1) * ndim + jp1];
                for (i = jp1; i <= n; ++i) {
                    tmp = a[(i - 1) * ndim + j];
                    temp = tmp * tmp / tempjj;
                    tdmin = a[(i - 1) * ndim + i] - temp;
                    jdmin = dmin(ctx, jdmin,tdmin);
                }
                if (jdmin < taugam)
                    phase1 = 0;           
            }
            else 
                phase1 = 0;         

            if (phase1) {         
 
                /* DO THE NORMAL CHOLESKY UPDATE IF STILL IN PHASE 1 */
 
                a[(j - 1) * ndim + j] = sqrt(a[(j - 1) * ndim + j]);
                tempjj = a[(j - 1) * ndim + j];
                for (i = jp1; i <= n; ++i) 
                    a[(i - 1) * ndim + j] /= tempjj;

                for (i = jp1; i <= n; ++i) {
                    temp = a[(i - 1) * ndim + j];
                    for (k = jp1; k <= i; ++k)
                        a[(i - 1) * ndim + k] -= (temp * a[(k - 1) * ndim + j]);
                }
                if (j == n - 1)
                    a[(n - 1) * ndim + n] = sqrt(a[(n - 1) * ndim + n]);
            }
            else {

                /* CALCULATE THE NEGATIVES OF THE LOWER GERSCHGORIN BOUNDS */
 
                gersch(ctx, ndim,n,a,j,g);
            }
        }
 
        /* PHASE 2 */
 
        if (phase1 == 0) {         

            if (j != n - 1) {       
 
                /* FIND THE MINIMUM NEGATIVE GERSHGORIN BOUND */

                iming = j;
                ming = g[j];
                for (i = j + 1; i <= n; ++i) {
                    if (ming > g[i]) {       
                        ming = g[i];
                        iming = i;  
                    }     
                }
   
                /* PIVOT TO THE TOP THE ROW AND COLUMN WITH THE */
                /* MINIMUM NEGATIVE GERSCHGORIN BOUND */
 
                if (iming != j) {         
  
                    /* SWAP ROW J WITH ROW OF MIN GERSCH BOUND */
  
                    for (i = 1; i < j; ++i) {
                        temp = a[(j - 1) * ndim + i];
                        a[(j - 1) * ndim + i] = a[(iming - 1) * ndim + i];
                        a[(iming - 1) * ndim + i] = temp;
                    }
 
                    /* SWAP COLJ WITH ROW IMING FROM J TO IMING */
 
                    for (i = j + 1; i < iming; ++i) {
                        temp = a[(i - 1) * ndim + j];
                        a[(i - 1) * ndim + j] = a[(iming - 1) * ndim + i];
                        a[(iming - 1) * ndim + i] = temp;
                    }
  
                    /* SWAP COLUMN J WITH COLUMN OF MIN GERSCH BOUND */
          
                    for (i = iming + 1; i <= n; ++i) {
                        temp = a[(i - 1) * ndim + j];
                        a[(i - 1) * ndim + j] = a[(i - 1) * ndim + iming];
                        a[(i - 1) * ndim + iming] = temp;
                    }
 
                    /* SWAP DIAGONAL ELEMENTS  */
 
                    temp = a[(j - 1) * ndim + j];
                    a[(j - 1) * ndim + j] = a[(iming - 1) * ndim + iming];
                    a[(iming - 1) * ndim + iming] = temp;
  
                    /* SWAP ELEMENTS OF THE PERMUTATION VECTOR */
     
                    itemp = p[j];
                    p[j] = p[iming];
                    p[iming] = itemp;
 
                    /* SWAP ELEMENTS OF THE NEGATIVE GERSCHGORIN BOUNDS VECTOR  */
 
                    temp = g[j];
                    g[j] = g[iming];
                    g[iming] = temp;
                }            
 
                /* CALCULATE DELTA AND ADD TO THE DIAGONAL. */
                /* DELTA=MAX{0,-A(J,J) + MAX{NORMJ,TAUGAM},DELTA_PREVIOUS} */
                /* WHERE NORMJ=SUM OF |A(I,J)|,FOR I=1,N, */
                /* DELTA_PREVIOUS IS THE DELTA COMPUTED AT THE PREVIOUS ITERATION, */
                /* AND TAUGAM IS TAU1*GAMMA. */

                normj = 0.0;
                for (i = j + 1; i <= n; ++i) 
                    normj += fabs(a[(i - 1) * ndim + j]);

                temp = dmax(ctx, normj,taugam);
                delta1 = temp - a[(j - 1) * ndim + j];
                temp = 0.0;
                delta1 = dmax(ctx, temp,delta1);
                delta =  dmax(ctx, delta1,delta);    
                e[j] = delta;
                a[(j - 1) * ndim + j] += e[j];
  
                /* UPDATE THE GERSCHGORIN BOUND ESTIMATES */
                /* (NOTE: G(I) IS THE NEGATIVE OF THE */
                /* GERSCHGORIN LOWER BOUND.) */
               
                if (a[(j - 1) * ndim + j] != normj) {         
                    temp = (normj / a[(j - 1) * ndim + j]) - 1.0;
 
                    for (i = j + 1; i <= n; ++i)
                        g[i] += fabs(a[(i - 1) * ndim + j]) * temp;
                }                       
 
                /* DO THE CHOLESKY UPDATE */
 
                tempjj = a[(j - 1) * ndim + j] = sqrt(a[(j - 1) * ndim + j]);
                for (i = j + 1; i <= n; ++i)
                    a[(i - 1) * ndim + j] /= tempjj;

                for (i = j + 1; i <= n; ++i) {
                    temp = a[(i - 1) * ndim + j];
                    for (k = j + 1; k <= i; ++k)
                        a[(i - 1) * ndim + k] -= (temp * a[(k - 1) * n + j]);
                }
            }
            else  
                fin2x2(ctx, ndim,n,a,e,tau2,&delta,gamma);
        }           
    }
}

/* ------------------------------------------------------------------------ */
/*  ts_init(n,ndim,a,phase1,delta,p,g,e,ming,tau1,gamma,taugam)             */
/*                                                                          */
/*  SET UP FOR START OF CHOLESKY FACTORIZATION                              */  
/*                                                                          */
/*  INPUT :    N, NDIM, A, TAU1                                             */
/*                                                                          */
/*  OUTPUT :    PHASE1    - BOOLEAN VALUE SET TO TRUE IF IN PHASE ONE,      */
/*              OTHERWISE FALSE.                                            */
/*  DELTA     - AMOUNT TO ADD TO AJJ AT ITERATION J                         */
/*  P,G,E     - DESCRIBED ABOVE IN MODCHL                                   */
/*  MING      - THE MINIMUM NEGATIVE GERSCHGORIN BOUND                      */
/*  GAMMA     - THE MAXIMUM DIAGONAL ELEMENT OF A                           */
/*  TAUGAM    - TAU1 * GAMMA                                                */
/*                                                                          */ 

void ts_init(TDAContext *ctx, int n,int ndim,double *a,int *phase1,double *delta,int *p,double *g, double *e,double *ming,double tau1,double *gamma,double *taugam)
{
    register int i;
    double tmp;

    *phase1 = 1;       
    *delta = 0.0;
    *ming = 0.0;
    for (i = 1; i <= n; ++i) {
        p[i] = i;
        g[i] = 0.0;
        e[i] = 0.0;
    }
 
    /* FIND THE MAXIMUM MAGNITUDE OF THE DIAGONAL ELEMENTS. */
    /* IF ANY DIAGONAL ELEMENT IS NEGATIVE, THEN PHASE1 IS FALSE. */
 
    *gamma = 0.0;
    for (i = 1; i <= n; ++i) {
        tmp = a[(i - 1) * ndim + i];
        *gamma = dmax(ctx, *gamma,fabs(tmp));
        if (tmp < 0.0)
            *phase1 = 0;         
    }
    *taugam = tau1 * *gamma;

    /* IF NOT IN PHASE1, THEN CALCULATE THE INITIAL GERSCHGORIN BOUNDS */
    /* NEEDED FOR THE START OF PHASE2. */
 
    if (*phase1 == 0)
        gersch(ctx, ndim,n,a,1,g);
}  

/* ------------------------------------------------------------------------ */
/*  gersch(ndim,n,a,j,g)                                                    */
/*                                                                          */
/*  CALCULATE THE NEGATIVE OF THE GERSCHGORIN BOUNDS                        */
/*  CALLED ONCE AT THE START OF PHASE II.                                   */
/*                                                                          */
/*  INPUT   : NDIM, N, A, J                                                 */
/*                                                                          */
/*  OUTPUT  : G - AN N VECTOR CONTAINING THE NEGATIVES OF THE               */
/*            GERSCHGORIN BOUNDS.                                           */

void gersch(TDAContext *ctx, int ndim,int n,double *a,int j,double *g)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,k;
    double offrow;

    for (i = j; i <= n; ++i) {
        offrow = 0.0;
        for (k = j; k < i; ++k)
            offrow += fabs(a[(i - 1) * ndim + k]);

        for (k = i + 1; k <= n; ++k) 
            offrow += fabs(a[(k - 1) * ndim + i]);
        g[i] = offrow - a[(i - 1) * ndim + i];
    }
}

/* ------------------------------------------------------------------------ */
/*  fin2x2(ndim,n,a,e,tau2,delta,gamma)                                     */
/*                                                                          */
/*  HANDLES FINAL 2X2 SUBMATRIX IN PHASE II.                                */
/*  FINDS EIGENVALUES OF FINAL 2 BY 2 SUBMATRIX,                            */
/*  CALCULATES THE AMOUNT TO ADD TO THE DIAGONAL,                           */
/*  ADDS TO THE FINAL 2 DIAGONAL ELEMENTS,                                  */
/*  AND DOES THE FINAL UPDATE.                                              */
/*                                                                          */
/*  INPUT : NDIM, N, A, E, TAU2,                                            */
/*          DELTA - AMOUNT ADDED TO THE DIAGONAL IN THE                     */
/*                  PREVIOUS ITERATION                                      */
/*                                                                          */
/*  OUTPUT : A - MATRIX WITH COMPLETE L FACTOR IN THE LOWER TRIANGLE,       */
/*           E - N*1 VECTOR CONTAINING THE AMOUNT ADDED TO THE DIAGONAL     */
/*               AT EACH ITERATION,                                         */
/*           DELTA - AMOUNT ADDED TO DIAGONAL ELEMENTS N-1 AND N.           */
/*                                                                          */

void fin2x2(TDAContext *ctx, int ndim,int n,double *a,double *e,double tau2,double *delta, double gamma)
{                                                                            
    double t1,t2,t1a,t2a,t3,lmbd1,lmbd2,lmbdhi,lmbdlo,delta1,temp;
 
    /* FIND EIGENVALUES OF FINAL 2 BY 2 SUBMATRIX */
 
    t1 = a[(n - 2) * ndim + n - 1] + a[(n - 1) * ndim + n];
    t2 = a[(n - 2) * ndim + n - 1] - a[(n - 1) * ndim + n];

    t1a = fabs(t2);
    t2a = 2.0 * fabs(a[(n - 1) * ndim + n - 1]);

    if (t1a >= t2a) {       
        if (t1a > 0.0)
            t2a /= t1a;
        t3 = t1a * sqrt(1.0 + (t2a * t2a));
    }     
    else {
        t1a = t1a / t2a;
        t3 = t2a * sqrt(1.0 + (t1a * t1a));
    }       
    lmbd1 = (t1 - t3) / 2.0;
    lmbd2 = (t1 + t3) / 2.0;
    lmbdhi = dmax(ctx, lmbd1,lmbd2);
    lmbdlo = dmin(ctx, lmbd1,lmbd2);
 
    /* FIND DELTA SUCH THAT:
       1.  THE L2 CONDITION NUMBER OF THE FINAL 
           2X2 SUBMATRIX + DELTA*I <= TAU2 
       2. DELTA >= PREVIOUS DELTA,
       3. LMBDLO + DELTA >= TAU2 * GAMMA, 
          WHERE LMBDLO IS THE SMALLEST EIGENVALUE OF THE FINAL 
          2X2 SUBMATRIX */

    delta1 = (lmbdhi - lmbdlo) / (1.0 - tau2);
    delta1 = dmax(ctx, delta1,gamma);
    delta1 = tau2 * delta1 - lmbdlo;
    temp = 0.0;
    *delta = dmax(ctx, *delta,temp);
    *delta = dmax(ctx, *delta,delta1);

    if (*delta > 0.0) {         
        a[(n - 2) * ndim + n - 1] += *delta;
        a[(n - 1) * ndim + n] += *delta;
        e[n - 1] = *delta;
        e[n] = *delta;
    }         
 
    /* FINAL UPDATE */
 
    a[(n - 2) * ndim + n - 1] = sqrt(a[(n - 2) * ndim + n - 1]);
    a[(n - 1) * ndim + n - 1] /= a[(n - 2) * ndim + n - 1];
    a[(n - 1) * ndim + n] -= a[(n - 1) * ndim + n - 1] * a[(n - 1) * ndim + n - 1];
    a[(n - 1) * ndim + n] = sqrt(a[(n - 1) * ndim + n]);
}

/* ------------------------------------------------------------------------ */
/*  slvmdl(nc,n,h,u,t,e,diag,s,g,pivot,w1,w2,w3,alpha,beta,nomin)           */
/*                                                                          */
/*  COMPUTE TENSOR AND NEWTON STEPS                                         */
/*                                                                          */
/*  NC          --> Column DIMENSION OF MATRIX                              */
/*  N           --> DIMENSION OF PROBLEM                                    */
/*  H(N,N)      --> HESSIAN                                                 */
/*  U(N)        --> VECTOR TO FORM Q IN QR                                  */
/*  T(N)        --> WORKSPACE                                               */
/*  E(N)        --> DIAGONAL ADDED TO HESSIAN IN CHOLESKY DECOMPOSITION     */
/*  DIAG(N)     --> DIAGONAL OF HESSIAN                                     */
/*  S(N)        --> STEP TO PREVIOUS POINT (FOR TENSOR MODEL)               */
/*  G(N)        --> CURRENT GRADIENT                                        */
/*  PIVOT(N)    --> PIVOT VECTOR FOR CHOLESKY DECOMPOSITION                 */
/*  W1(N)      <--> ON INPUT: A=2*(GP-G-HS-S*BETA/(6*STS))                  */
/*                  ON OUTPUT: TENSOR STEP                                  */
/*  W2(N)       --> SH                                                      */
/*  W3(N)      <--  NEWTON STEP                                             */
/*  ALPHA       --> SCALAR FOR 3RD ORDER TERM OF TENSOR MODEL               */
/*  BETA        --> SCALAR FOR 4TH ORDER TERM OF TENSOR MODEL               */
/*  NOMIN      <--  =.TRUE. IF TENSOR MODEL HAS NO MINIMIZER                */

void slvmdl(TDAContext *ctx, int nc,int n,double *h,double *u,double *t,double *e,double *diag, double *s,double *g,int *pivot,double *w1,double *w2,double *w3, double alpha,double beta,int *nomin)
{
    register int i,j;
    double shs,w11,w12,w13,w22,w23,ca,cb,cc,sg,cd,r,r1,r2,temp,uu,ss;
    double addmax,tmp,sgstar;
 
    /* solve model */
 
    addmax = 0.0;
    *nomin = 0;       
 
    /* COMPUTE QTHQ(N,N), ZTHZ(N-1,N-1) = FIRST N-1 ROWS AND N-1 */
    /* COLUMNS OF QTHQ */

    if (n > 1) {       
        zhz(ctx, nc,n,s,h,u,t);
	
        /* IN CHOLESKY DECOMPOSITION WILL STORE H(1,1) ... H(N-1,N-1) */
        /* IN DIAG(1) ... DIAG(N-1), STORE H(N,N) IN DIAG(N) FIRST */

        diag[n] = h[(n - 1) * nc + n];
 
        /* COLESKY DECOMPOSITION FOR FIRST N-1 ROWS AND N-1 COLUMNS OF ZTHZ */
        /* ZTHZ(N-1,N-1)=LLT */

        choldr(ctx, nc,n - 1,h,t,pivot,e,diag,addmax);
    }           
 
    /*  ON INPUT: SH IS STORED IN W2 */

    shs = 0.0;
    for (i = 1; i <= n; ++i) {
        shs += w2[i] * s[i];
        w3[i] = g[i];
    }
 
    /* COMPUTE W1,W2,W3
       W1=L**-1*ZT*A
       W2=L**-1*ZT*SH
       W3=L**-1*ZT*G */
 
    if (n > 1) {       

        solvew(ctx, nc,n,h,u,w1,t);
        solvew(ctx, nc,n,h,u,w2,t);
        solvew(ctx, nc,n,h,u,w3,t);
    }
 
    /*  COMPUTE COEFFICIENTS CA,CB,CC AND CD FOR REDUCED ONE VARIABLE */
    /*  3RD ORDER EQUATION */

    w11 = 0.0;
    for (i = 1; i < n; ++i)
        w11 += w1[i] * w1[i];

    ca = beta / 6.0 - w11 / 2.0;
    w12 = 0.0;
    for (i = 1; i < n; ++i)
        w12 += w1[i] * w2[i];

    cb = alpha / 2.0 - 3.0 * w12 / 2.0;
    w13 = 0.0;
    for (i = 1; i < n; ++i)
        w13 += w1[i] * w3[i];

    w22 = 0.0;
    for (i = 1; i < n; ++i)
        w22 += w2[i] * w2[i];

    cc = shs - w22 - w13;
    sg = 0.0;
    for (i = 1; i <= n; ++i)
        sg += s[i] * g[i];

    w23 = 0.0;
    for (i = 1; i < n; ++i)
        w23 += w2[i] * w3[i];

    cd = sg - w23;
 
    /* COMPUTE DESIRABLE ROOT, SGSTAR, OF 3RD ORDER EQUATION */
  
    if (ca != 0.0) {         
        sigma(ctx, &sgstar,ca,cb,cc,cd);
        if (ca == 0.0) {      
            *nomin = 1;     
            goto L200;
        }         
    }
    else {

        /* 2ND ORDER ( CA=0 ) */

        if (cb != 0.0) {      
            r = cc * cc - 4.0 * cb * cd;
            if (r < 0.0) {     
                *nomin = 1;   
                goto L200;
            }
            else {
                r1 = (-cc + sqrt(r)) / (2.0 * cb);
                r2 = (-cc - sqrt(r)) / (2.0 * cb);
                if (r2 < r1) {       
                    temp = r1;
                    r1 = r2;
                    r2 = temp;
                }
                if (cb > 0.0)       
                    sgstar = r2;
                else
                    sgstar = r1;

                if(((r1 > 0.0) && (sgstar == r2)) ||  
                                        ((r2 < 0.0) && (sgstar == r1))) {       
                    *nomin = 1;   
                    goto L200;
                }        
            }         
        }
        else {

            /* 1ST ORDER (CA=0,CB=0) */

            if (cc > 0.0)         
                sgstar = -cd / cc;
            else {
                *nomin = 1;  
                goto L200;
            }           
        }         
    }           
 
    /* FIND TENSOR STEP, W1 (FUNCTION OF SGSTAR) */

    dstar(ctx, nc,n,u,s,w1,w2,w3,sgstar,h,w1);
 
    /* COMPUTE DN */

L200:
    uu = ss = 0.0;
    for (i = 1; i <= n; ++i) {
        uu += u[i] * u[i];
        ss += s[i] * s[i];
    }
    uu /= 2.0;
    ss = sqrt(ss);
      
    if (n == 1) {       
        choldr(ctx, nc,n,h,t,pivot,e,diag,addmax);
    }
    else {
 
        /* COMPUTE LAST ROW OF L(N,N) */

        for (i = 1; i < n; ++i) {
            temp = 0.0;
            if (i > 1) {       
                for (j = 1; j < i; ++j)
                    temp += h[(n - 1) * nc + j] * h[(i - 1) * nc + j];
            }
            h[(n - 1) * nc + i] = 
                    (h[(i - 1) * nc + n] - temp) / h[(i - 1) * nc + i];

        }
        temp = 0.0;
        for (i = 1; i < n; ++i)
            temp += h[(n - 1) * nc + i] * h[(n - 1) * nc + i];

        tmp = h[(n - 1) * nc + n] = diag[n] - temp + addmax;
        if (tmp > 0.0)         
            h[(n - 1) * nc + n] = sqrt(tmp);
           
        else {

            /* AFTER ADDING THE LAST COLUMN AND ROW */
            /* QTHQ IS NOT P.D., NEED TO REDO CHOLESKY DECOMPOSITION */

            for (i = 2; i <= n; ++i) {
                for (j = 1; j < i; ++j)
                    h[(i - 1) * nc + j] = h[(j - 1) * nc + i];
                h[(i - 1) * nc + i] = diag[i];
            }
            h[1] = diag[1];
            choldr(ctx, nc,n,h,t,pivot,e,diag,addmax);
        }         
    }       
 
    /* SOLVE QTHQ*QT*W3=-QT*G, WHERE W3 IS NEWTON STEP */
    /* W2=-QT*G */

    for (i = 1; i <= n; ++i) {
        w2[i] = 0.0;
        for (j = 1; j <= n; ++j)
            w2[i] += u[j] * u[i] * g[j] / uu;
        w2[i] -= g[i];
    }

    forslv(ctx, nc,n,h,w3,w2);
    bakslv(ctx, nc,n,h,w2,w3);               

    /* W2=QT*W3 => W3=Q*W2 --- NEWTON STEP */

    for (i = 1; i <= n; ++i) {
        w3[i] = 0.0;
        for (j = 1; j <= n; ++j)
            w3[i] += u[i] * u[j] * w2[j] / uu;
        w3[i] = w2[i] - w3[i];
    }
}           

/* ------------------------------------------------------------------------ */
/*  optstp(n,xpls,fpls,gpls,x,icscmx,fscale,mxtake,rgx,rsx,iretcd)          */
/*                                                                          */
/*  UNCONSTRAINED MINIMIZATION STOPPING CRITERIA                            */
/*                                                                          */
/*  FIND WHETHER THE ALGORITHM SHOULD TERMINATE, DUE TO ANY                 */
/*  OF THE FOLLOWING:                                                       */
/*  1) PROBLEM SOLVED WITHIN USER TOLERANCE                                 */
/*  2) CONVERGENCE WITHIN USER TOLERANCE                                    */
/*  3) ITERATION LIMIT REACHED                                              */
/*  4) DIVERGENCE OR TOO RESTRICTIVE MAXIMUM STEP (STEPMX) SUSPECTED        */
/*                                                                          */
/*  N            --> DIMENSION OF PROBLEM                                   */
/*  XPLS(N)      --> NEW ITERATE X[K]                                       */
/*  FPLS         --> FUNCTION VALUE AT NEW ITERATE F(XPLS)                  */
/*  GPLS(N)      --> GRADIENT AT NEW ITERATE, G(XPLS), OR APPROXIMATE       */
/*  X(N)         --> OLD ITERATE X[K-1]                                     */
/*  ICSCMX      <--> NUMBER CONSECUTIVE STEPS .GE. STEPMX                   */
/*                   [RETAIN VALUE BETWEEN SUCCESSIVE CALLS]                */
/*  FSCALE       --> ESTIMATE OF SCALE OF OBJECTIVE FUNCTION                */
/*  MXTAKE       --> BOOLEAN FLAG INDICATING STEP OF MAXIMUM LENGTH USED    */
/*                                                                          */
/*  Uses the following global variables:                                    */  
/*                                                                          */
/*  TOLSG   TOLERANCE AT WHICH RELATIVE GRADIENT CONSIDERED CLOSE           */
/*          ENOUGH TO ZERO TO TERMINATE ALGORITHM                           */
/*  TOLSP   RELATIVE STEP SIZE AT WHICH SUCCESSIVE ITERATES                 */
/*          CONSIDERED CLOSE ENOUGH TO TERMINATE ALGORITHM                  */
/*                                                                          */
/*  Input: iretcd = return code from line search.                           */
/*          0 if line search successful,                                    */
/*          1 if line search not successful.                                */
/*  Note: for initial check optstp() is called with iretcd = 2.             */
/*                                                                          */
/*  Output: rgx = scaled gradient used for convergence check.               */
/*          rsx = scaled parametr change for convergence check.             */
/*                                                                          */
/*  Return: itrmcd  = 1 convergence reached with gradient check             */
/*                    2 convergence reached with parameter change           */
/*                    3 no success in last line search                      */
/*                    4 exceeded max number of iterations                   */
/*                    5 five consecutive max steps in line search           */
/*                    0 otherwise                                           */

int optstp(TDAContext *ctx, int n,double *xpls,double fpls,double *gpls,double *x,            int *icscmx,double fscale,int mxtake,double *rgx,double *rsx,int iretcd)
{                                                                       
    register int i;
    double d,relgrd,relstp;

    *rgx = 0.0;
    *rsx = 0.0;
    ctx->TOLSGF = 0.0;       /* final scaled gradient */
    ctx->TOLSPF = 0.0;       /* final scaled parameter change */
                                                                        
    if (iretcd == 1)        /* last linesearch failed */               
        return(3);
                                                                        
    /* FIND DIRECTION IN WHICH RELATIVE GRADIENT MAXIMUM. */                 
    /* CHECK WHETHER WITHIN TOLERANCE */                                     
                                                                        
    d = dmax(ctx, fabs(fpls),fscale);                                          

    for (i = 1; i <= n; ++i) {
        relgrd = fabs(gpls[i]) * dmax(ctx, fabs(xpls[i]),1.0) / d;                     
        *rgx = dmax(ctx, *rgx,relgrd);                                            
    }
    ctx->TOLSGF = *rgx;
    if (ctx->PMProtFDef > 1)  
        prval(ctx, "\nScaled gradient",ctx->TOLSGF);

    if (ctx->Iter == 0) {                                                 
        if (*rgx <= ctx->TOLSG)
            return(1);                                                      
        return(0);                                                 
    }

    /* FIND DIRECTION IN WHICH RELATIVE STEPSIZE MAXIMUM */                  
    /* CHECK WHETHER WITHIN TOLERANCE. */                                    
                                                                        
    for (i = 1; i <= n; ++i) {
        relstp = fabs(xpls[i] - x[i]) / dmax(ctx, fabs(xpls[i]),1.0);                  
        *rsx = dmax(ctx, *rsx,relstp);                                            
    }
    ctx->TOLSPF = *rsx;
    if (ctx->PMProtFDef > 1)  
        prval(ctx, "Scaled parameter change",ctx->TOLSPF);

    if (*rgx <= ctx->TOLSG)
        return(1);                                                      

    if (*rsx <= ctx->TOLSP) 
        return(2);                                                      

    /* CHECK ITERATION LIMIT */                                              

    if (ctx->Iter >= ctx->MxIter)
        return(4);                                   
                                                                        
    if (mxtake) {            
        *icscmx += 1;
        if (*icscmx >= 5)
            return(5);                                                  
    }   
    else
        *icscmx = 0;        /* accumulate consecutive max steps */

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ts_func(n,x,f,gflg,g,hflg,h)                                            */
/*                                                                          */
/*  Calculation of function value for current parameter vector x[].         */
/*  If gflg return also gradient in g[].                                    */
/*  If hflg return also hessian in h[].                                     */
/*                                                                          */
/*  If TSFNC = 1 (user-defined functions) return error flag from t_eval().  */
/*                                                                          */

int ts_func(TDAContext *ctx, int n,double *x,double *f,int gflg,double *g,int hflg,double *h)
{
    register int i,j,k;
    int err,m,deriv;

    if (ctx->TSFNC == 0) {       /* standard models, use fn(ctx) */
        m = 0;
        if (gflg)
            m++;
        if (hflg)
            m++;

        fn(ctx, x,ctx->TSMOD,m,1,0,&err);
        if (err)
            return(err);

        *f = ctx->FTmp;

        if (gflg) {
            for (i = 1; i <= n; ++i)
                g[i] = ctx->Grad[i];
        }
        if (hflg) {
            k = 1;
            for (i = 1; i <= n; ++i) {
                h[(i - 1) * n + i] = ctx->Diag[i];
                for (j = 1; j < i; ++j) {
                    h[(j - 1) * n + i] =
                    h[(i - 1) * n + j] = ctx->Hess[k++];
                }
            }
        }
    }
    else {                  /* user-defined functions */
        ctx->NumF++;
        deriv = 0;
        if (gflg) {
            deriv = 1;
            ctx->NumFG++;
        }
        if (hflg) {
            deriv = 2;
            ctx->NumFH++;
        }
        m = get_flval(ctx, f,n,x + 1,deriv,1,g,h,h);  
        if (m) {
            *f = 0.0;
            return(m);
        }
    } 
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ts_prn(iter,f,gnorm,pc,n,x,g)                                           */  
/*                                                                          */

void ts_prn(TDAContext *ctx, int iter,double f,double gnorm,double pc,int n,double *x,double *g)
{
    if (iter == 0) {
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

        if (ctx->PMProtFDef == 1)  
            fprintf(ctx->PMProtFd,"\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");
    }
    if (ctx->SILENTFlg < 2) {
        printfe(ctx, "  %3d  %20.13e %17.10e ",iter,f,gnorm);
        printfe(ctx, "%11.4e  %6d (%d,%d)\n",pc,ctx->NumF,ctx->NumFG,ctx->NumFH);
    }
    if (ctx->PMProtFDef) {
        if (ctx->PMProtFDef > 1)  
            fprintf(ctx->PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
        fprintf(ctx->PMProtFd,"  %3d  %20.13e %17.10e ",iter,f,gnorm);
        fprintf(ctx->PMProtFd,"%11.4e  %6d (%d,%d)\n",pc,ctx->NumF,ctx->NumFG,ctx->NumFH);

        if (ctx->PMProtFDef > 1) {
            prvec(ctx, "\nParameter",n,x);
            prvec(ctx, "Gradient",n,g);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  ts_prn1(n,h)                                                            */  
/*                                                                          */

void ts_prn1(TDAContext *ctx, int n,double *h)
{
    register int i,j;

    if (ctx->PMProtFDef > 1) {
        fprintf(ctx->PMProtFd,"Hessian\n");
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j)
                rt_fprintf_d(ctx, ctx->PMProtFd,ctx->PMPFmtS,h[(i - 1) * n + j]);
            fprintf(ctx->PMProtFd,"\n");
        }
    }
}






