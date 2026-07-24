/****************************************************************************/
/*  t_min                                                                   */
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
#include "t_gf.h"
#include "t_ml.h"
#include "t_svd.h"
#include "t_rand.h"
#include "t_con.h"
#include "t_alloc.h"
#include "t_tmin.h"
#include "t_gdat.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"

/*  functions in t_min.c */

int ffmin(int mod,int typ,int cov,int msg);
int ffmin1(int mod);
int ffmin2(int mod);
int ffmin3(int mod);
int ffmin5(int mod);
void exphess(void);
int chinv(int opt);
int ffmin7(int mod);

/* ------------------------------------------------------------------------ */
/*  ffmin(mod,typ,cov,msg)      Function minimization                       */
/*                                                                          */
/*          typ 0 : ML                                                      */
/*              1 : gen function                                            */
/*              2 : nonlinear regression                                    */
/*                                                                          */
/*          If cov == 0 do not calculate covariance matrix.                 */
/*          If msg == 0 do not print messages.                              */
/*                                                                          */
/*          a) Call specified minimization algorithm. If constraints are    */
/*             defined, this is done in the restricted parameter space.     */
/*          b) If successful and cov != 0 additional call of fn_function to */
/*             calculate the covariance matrix (this depends on the user-   */
/*             specified type, given by CCTyp). If CCTyp = 0 the cov.       */
/*             matrix is not calculated.                                    */
/*                                                                          */
/*             If typ = 2, the cov matrix is scaled with                    */
/*                                                                          */
/*              2 * f / (NOC - NParm)                                       */
/*                                                                          */
/*          c) If constraints are defined, the final parameter vector and   */
/*             the covariance matrix are expanded to the orginal space.     */
/*          d) Print result to the iteration protocol.                      */
/*          e) If requested print covariance matrix to output file.         */
/*                                                                          */
/*          Result:                                                         */
/*              LConv depending on the success of the function.             */
/*              FMax  log likelihood at maximum                             */
/*              Par   parameter vector in the original parameter space      */
/*              Hess/Diag   covariance matrix in the original space.        */
/*                                                                          */
/*          LConv 1 cannot increase log likelihood                          */
/*                2 exceeded max number of iterations                       */
/*               -1 final hessian (or outer product) not positive definite  */
/*               -2 step size search failed                                 */
/*               -3 cannot calculate eigenvalues and vectors of the hessian */
/*               -4 search vector is no descent direction                   */
/*               -5 exceeded max iterations in step size search             */
/*               -6 insufficient memory.                                    */
/*               -7 error in function evaluation                            */
/*                                                                          */
/*      Note: This function also calculates:                                */
/*      a)  LM test statistics, if LMTestFlg = 1.                           */
/*                                                                          */  
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */                                
  
int ffmin(int mod,int typ,int cov,int msg)
{
    register int i,j,k,l,ij;
    int r,err,rr,n;
    double tmp,tmp1,tmp2,tmp3;

    err = -1;
    Iter = LConv = CGradFlg = CCovFlg = MCovFlg = NumF = NumFG = NumFH = 0;
    NINTMUsed = -1;
    clear_npflg();          /* clear numerical problem flags */

    if (MxIter < 1) {
        LConv = 2;
        return(0);
    }
    /* fflush(stdout); */

    FN_First = 1;   /* set for first call of fn() */


    /*  Call minimization algorithm specified by MINA.
        The resulting parameters are given in Par with dimension NParm1.
        So, if constraints are defined, Par is in the reduced space.
        The evaluated function value is returned in FMin. */

    if (PMProtFDef) {
        fprintf(PMProtFd,"Function minimization with algorithm %d",MINA);
        if (MINA == 7 || MINA == 8)
            fprintf(PMProtFd," [dopt=%d]",PMDOPT);
        fprintf(PMProtFd,".\n");
    }

    switch (MINA) {

        case 1:     /* Direct search */
                    r = ffmin1(mod);
                    break;
        case 2:     /* Simplex method */
                    r = ffmin2(mod);
                    break;
        case 3:     /* conjugate gradient or BFGS methods */
        case 4:     r = ffmin3(mod);
                    break;

        case 7:     /* CES (qudratic and tensor model) */
        case 8:     r = ffmin7(mod);
                    break;

        default:    /* Newton algorithms (MINA 5 and 6) */
                    r = ffmin5(mod);
                    break;
    }
    if (r) {            /* insufficient memory or error in function */
        if (r < 0) {
            err = -1;
            LConv = -6;
        }
        else {
            err = r;     
            LConv = -7;
        }
        goto FFMFin;
    }
    if (SILENTFlg < 2) {
        printfe("\n");
        fflushe();
    }   
    if (PMProtFDef) {
        switch (LConv) {
          case  0:  fprintf(PMProtFd,"\nConvergence reached with criterion %d.\n\n",Crite);
                    break;
          case  1:  fprintf(PMProtFd,"\nCannot increase log likelihood.\n\n");
                    break;
          case  2:  fprintf(PMProtFd,"\nExceeded max number of iterations.\n\n");
                    break;
          case -1:  fprintf(PMProtFd,"\nFinal hessian or outer product not positive definite.\n\n");
                    break;
          case -2:  fprintf(PMProtFd,"\nStep size search failed.\n\n");
                    break;
          case -3:  fprintf(PMProtFd,"\nCannot calculate eigenvalues and vectors of the hessian.\n\n");
                    break;
          case -4:  fprintf(PMProtFd,"\nSearch vector not in descent direction.\n\n");
                    break;
          case -5:  fprintf(PMProtFd,"\nExceeded max iterations in step size search.\n\n");
                    break;
        }
    }

    FMax = FMin / DScal;

    if (PMMPLogDef == 1)                /* create matrix for loglikelihood */
        mp_putlog(FMax);

    /*  If successful calculate covariance matrix. This depends on the
        type defined by CCTyp. All calculations are done in a parameter space
        which is reduced if constraints are present.
        Note: fn() is called without the scaling option.    

        Note: If cov = 0 or CCTyp = 0, cov matrix is not calculated */
       
    if (PMMPGradDef == 1) {
        if (cov == 0 || CCTyp == 0 || NOCUsed < 1)
            PMMPGradDef = -1;
        else if (mp_alloc(4,NOCUsed,NParm + PM2NV))                 
            PMMPGradDef = -1;
    }

    if (cov == 0 || CCTyp == 0)
        goto FFMCon;

    if (LConv >= 0 && CCTyp > 0) {

        if (PMProtFDef)  
            fprintf(PMProtFd,"Covariance matrix calculation with ccov=%d\n\n",CCTyp);
                 
        if (PMMPGradDef == 1) {     /* calculate gradients */
            CGradFlg = 1;
            fn(Par,mod,1,0,0,&rr);
            CGradFlg = 0;
            if (rr) {
                err = rr;
                goto FFMFin;
            }
        }

        MCovFlg = 1;    /* signals calculation of covariance matrix */

        /* calculate outer product */

        if (CCTyp == 1 || CCTyp == 3) {
    
            CCovFlg = 1;
            fn(Par,mod,2,0,0,&rr);
            CCovFlg = 0;
            if (rr) {
                err = rr;
                goto FFMFin;
            }

            if (CCTyp == 1) {
                if ((i = syminv1(NParm1,Diag,Hess)) != NParm1) { 
               
                    if (i < 0) {    /* insufficient memory */
                        LConv = -6;
                        err = -1;
                        goto FFMFin; 
                    }
                    LConv = -1;
                }
            }
            else {  /* save outer product in Hess1/Diag1 */

                if (!(Diag1 = (double *)calloc(NParm1 + 1,sizeof(double)))) { 
                    LConv = -6;
                    err = -1;
                    goto FFMFin; 
                }
                Diag1A = NParm1 + 1;
                memrq(Diag1A,sizeof(double));

                if (!(Hess1 = (double *)calloc(HSiz1 + 1,sizeof(double)))) { 
                    LConv = -6;
                    err = -1;
                    goto FFMFin; 
                }
                Hess1A = HSiz1 + 1;
                memrq(Hess1A,sizeof(double));

                for (i = 1; i <= NParm1; ++i)    
                    Diag1[i] = Diag[i];
                for (i = 1; i <= HSiz1; ++i)
                    Hess1[i] = Hess[i];
            }
        }

        /*  calculate hessian in Hess / Diag */

        if (CCTyp == 2 || CCTyp == 3) {   
   
            /* if algorithm 7 or 8 and we have only function values, then
               use the Hessian returned from ffmin7() */
       
            if ((MINA == 7 || MINA == 8) && PMDOPT == 0)
                goto FFMCONT;
       
            if (LMTestFlg && NCon1) { 
                fn(Par,mod,2,0,1,&rr);       /* LM test if LMTestFlg = 1 */
                LMTest = FNTStat;
                LMRank = FNTRank;
            }
            else
                fn(Par,mod,2,0,0,&rr);

            if (rr) {
                err = rr;
                goto FFMFin;
            }
FFMCONT:         
            if (typ == 0) {
                for (i = 1; i <= NParm1; ++i)    
                    Diag[i] *= -1.0;
                for (i = 1; i <= HSiz1; ++i)
                    Hess[i] *= -1.0;
            }
            if ((i = syminv1(NParm1,Diag,Hess)) != NParm1) {   
                if (i < 0) {
                    LConv = -6; /* insufficient memory */
                    err = -1;
                    goto FFMFin;
                }
                else  
                    LConv = -1;
 
            }

            if (CCTyp == 3) {   /* calculate cross product */

                if (LConv >= 0) {
            
                    for (i = 1; i <= NParm1; ++i)  
                        WrkD[i] = Diag[i];
                    for (i = 1; i <= HSiz1; ++i)
                        WrkH[i] = Hess[i];
                    ij = 1;
                    for (i = 1; i <= NParm1; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = 0.0;
                            for (k = 1; k <= NParm1; ++k) {
                                tmp1 = 0.0;
                                for (l = 1; l <= NParm1; ++l) {
                                    if (l == k)
                                        tmp3 = Diag1[k];
                                    else if (l < k)
                                        tmp3 = Hess1[(k - 1) * (k - 2) / 2 + l];
                                    else
                                        tmp3 = Hess1[(l - 1) * (l - 2) / 2 + k];
    
                                    if (j == l)
                                        tmp2 = WrkD[l];
                                    else if (j < l)
                                        tmp2 = WrkH[(l - 1) * (l - 2) / 2 + j];
                                    else
                                        tmp2 = WrkH[(j - 1) * (j - 2) / 2 + l];
 
                                    tmp1 += tmp3 * tmp2;
                                }
                                if (k == i)
                                    tmp2 = WrkD[i];
                                else if (k < i)
                                    tmp2 = WrkH[(i - 1) * (i - 2) / 2 + k];
                                else
                                    tmp2 = WrkH[(k - 1) * (k - 2) / 2 + i];
                            
                                tmp += tmp1 * tmp2;
                            }
                            if (i == j)
                                Diag[i] = tmp;
                            else
                                Hess[ij++] = tmp;
                        }
                    }
                }
            }
        }
        MCovFlg = 0;

        tmp = 0.0;                      /* calculate norm of gradient */
        for (i = 1; i <= NParm1; ++i)  
            tmp += Grad[i] * Grad[i];
        CValG = sqrt(tmp);
    }
        
    /*  If constraints are defined, expand the parameter vector and the
        covariance matrix. */

FFMCon:

    if (NCon1) {     

        for (i = 1; i <= NParm; ++i) {
            tmp = 0.0;
            for (j = 1; j <= NParm1; ++j)
                tmp += ConQ[(CIP[i] - 1) * NParm1 + j] * Par[j];
            WrkD[i] = tmp;
        }
        for (i = 1; i <= NParm; ++i)  
            Par[i] = ParS[i] + WrkD[i];

        if (CCTyp > 0)
            exphess();
    }

    if (PMMPParDef == 1)                /* create matrix for parameters */
        mp_putpar(NParm,Par);
   
    if (PMProtFDef) {
        prvec("Parameter vector",NParm,Par);

        if (LConv >= 0 && CCTyp > 0)
            prhess("Covariance matrix",NParm);

        fprintf(PMProtFd,"\nEnd of algorithm.\n\n");
    }
    if (cov == 0 || CCTyp == 0) {
        err = 0;
        goto FFMFin;
    }

    /* scale cov matrix for nonlinear regression */

    n = NOC - (NParm - NCon1);

    if (typ == 2 && CCTyp == 2 && LConv >= 0 && n > 0) {

        tmp = 2.0 * FMax / (double)n;

        for (j = 1; j <= NParm; ++j)
            Diag[j] *= tmp;

        for (j = 1; j <= HSiz; ++j)
            Hess[j] *= tmp;
    }
    if (PMPPFDef || PMCovFDef || PMMPCovDef) {
     
        if (PMMPCovDef == 1) {              /* allocate matrix for cov matrix */
            if (mp_alloc(3,NParm,NParm))                 
                PMMPCovDef = -1;
        }

        for (i = 1; i <= NParm; ++i) {

            if (PMPPFDef)  
                fprintf(PMPPFd,PMTFmtS,Par[i]);

            if (LConv >= 0) {

                if (PMPPFDef) {
                    tmp = Diag[i];
                    if (tmp > 0.0)
                        tmp = sqrt(tmp);
                    else
                        tmp = 0.0;
                    fprintf(PMPPFd,PMTFmtS,tmp);
                }
                if ((PMCovFDef || PMMPCovDef == 1) && CCTyp > 0) {        

                    for (j = 1; j <= NParm; ++j) {
                        if (j == i)
                            tmp = Diag[i];
                        else if (j < i)
                            tmp = Hess[(i - 1) * (i - 2) / 2 + j];
                        else
                            tmp = Hess[(j - 1) * (j - 2) / 2 + i];
                  
                        if (PMCovFDef)
                            fprintf(PMCovFd,PMMFmtS,tmp);

                        if (PMMPCovDef == 1)
                            MatVal[MPCovIdx][(i - 1) * NParm + j] = tmp;    
                    }
                    if (PMCovFDef)
                        fprintf(PMCovFd,"\n");
                }
            }
            if (PMPPFDef)
                fprintf(PMPPFd,"\n");
        }
        if (PMPPFDef)  
            PMPPWFlg = 1;
        if (LConv >= 0 && PMCovFDef && CCTyp > 0)
            PMCovWFlg = 1;
    }
    err = 0;

FFMFin:
    if (msg) {
        if (err < 0)
            printf1("\nError: insufficient memory for function minimization.\n");
        else if (err > 0) {
            printf1("\nError in function evaluation.\nCannot continue.\n");
            prn_npflg();
        }
    }
    if (Diag1A) {
        free((char *)Diag1);
        memrq(-Diag1A,sizeof(double));
        Diag1A = 0;
    }
    if (Hess1A) {
        free((char *)Hess1);
        memrq(-Hess1A,sizeof(double));
        Hess1A = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin1                                                                  */
/*      Function minimization with a direct search method.                  */
/*      References:                                                         */
/*      A. F. Kaupe, Direct Search, Algorithm 178, Communications of the    */
/*      ACM 6 (1963), 313                                                   */
/*      M. Bell, M.C. Pike, Comm. ACM 9 (1966), 684                         */
/*      R. de Vogelaere, Comm. ACM 11 (1968), 498                           */
/*      F.K. Tomlin, L.B. Smith, Comm. ACM 1969, 637                        */
/*      L.B. Smith, Comm. ACM 1969, 638                                     */
/*                                                                          */
/*      The current implementation corresponds to the proposals of Bell     */
/*      and Pike with modifications proposed by Tomlin & Smith and Smith.   */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */                                

int ffmin1(int mod)
{
    register int k;
    int err;
    double delta,fm,fm1,theta,tmp,*step;

    err = 0;

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    if (PMProtFDef == 1)
        fprintf(PMProtFd,"\n  Iter    Function Value       Step Length     Par Change   FCall\n");

    delta = SLen;
    if (!(step   = (double *) calloc(NParm + 1,sizeof(double))))
        return(-1); 
    memrq(NParm + 1,sizeof(double));
        
    for (k = 1; k <= NParm1; ++k) {
        tmp = fabs(Par[k]);
        if (tmp < EPSI1)
            step[k] = delta;
        else
            step[k] = delta * tmp;
    }
    fn(Par,mod,0,1,0,&err);
    if (err)
        goto FM1Fin;

    FMin = FTmp; 

FM1L1:  
    fm = FMin;
    for (k = 1; k <= NParm1; ++k)
        Par1[k] = Par[k];

    for (k = 1; k <= NParm1; ++k) {
        Par1[k] += step[k];
        fn(Par1,mod,0,1,0,&err);
        if (err)
            goto FM1Fin;

        fm1 = FTmp; 

        if (fm1 < fm) 
            fm = fm1;
        else {
            step[k] = -step[k];
            Par1[k] += 2.0 * step[k];

            fn(Par1,mod,0,1,0,&err);
            if (err)
                goto FM1Fin;

            fm1 = FTmp; 

            if (fm1 < fm)
                fm = fm1;
            else
                Par1[k] -= step[k];
        }
    }
    if (fm < FMin) {

FM1L2:
        for (k = 1; k <= NParm1; ++k) {
            if ((Par1[k] > Par[k] && step[k] <  0.0) ||
                        (Par1[k] <= Par[k] && step[k] >= 0.0))  
                step[k] = -step[k];
            theta = Par[k];
            Par[k] = Par1[k];
            Par1[k] = 2.0 * Par1[k] - theta;
        }
        FMin = fm;
        fn(Par1,mod,0,1,0,&err);
        if (err)
            goto FM1Fin;

        fm = fm1 = FTmp; 

        for (k = 1; k <= NParm1; ++k) {
            Par1[k] += step[k];

            fn(Par1,mod,0,1,0,&err);
            if (err)
                goto FM1Fin;

            fm1 = FTmp; 

            if (fm1 < fm) 
                fm = fm1;
            else {
                step[k] = -step[k];
                Par1[k] += 2.0 * step[k];

                fn(Par1,mod,0,1,0,&err);
                if (err)
                    goto FM1Fin;

                fm1 = FTmp; 
                if (fm1 < fm)
                    fm = fm1;
                else
                    Par1[k] -= step[k];
            }
        }
        if (fm >= FMin)
            goto FM1L1;

        for (k = 1; k <= NParm1; ++k) {
            if (fabs(Par1[k] - Par[k]) > 0.5 * fabs(step[k])) 
                goto FM1L2;
        }
    }
    Iter++;

    if (SILENTFlg < 2) {
        printfe("  %3d  %20.13e %17.10e ",Iter,FMin,delta);
        printfe("        --   %6d\n",NumF);
    }
    if (PMProtFDef) {
        if (PMProtFDef > 1)
            fprintf(PMProtFd,"\nITER      Function Value       Step Length     Par Change   FCall\n");
        fprintf(PMProtFd,"  %3d  %20.13e %17.10e ",Iter,FMin,delta);
        fprintf(PMProtFd,"        --   %6d\n",NumF);
    }

    if (delta > TOLS) {
        if (Iter >= MxIter) {
            LConv = 2;
            goto FM1Fin;
        }
        delta *= SRed;
        for (k = 1; k <= NParm1; ++k)
            step[k] *= SRed;
        goto FM1L1;
    }
    else  
        LConv = 0;

FM1Fin:
    if (err)
        err = 1;

    free((char *)step);
    memrq(-NParm - 1,sizeof(double));
    return(err);
}  
         
/* ------------------------------------------------------------------------ */
/*  ffmin2                                                                  */
/*      Function minimization with a Simplex algorithm.                     */
/*                                                                          */
/*      Ref.: R. O'Neill, Function Minimization using a Simplex Procedure,  */
/*      Algorithm AS 47, Applied Statistics 20 (1971), 338 - 345.           */
/*      Nelder and Mead, Computer Journal 7 (1965), 308-313.                */
/*                                                                          */
/*      Remarks:                                                            */
/*      1. J.M. Chambers, J.E. Ertel, Applied Statistics 23 (1974),250      */
/*      2. R. O'Neill, Applied Statistics 23 (1974), 252                    */
/*      3. P.R. Benyon, Applied Statistics 25 (1976), 97                    */
/*      4. I.D. Hill, Applied Statistics 27 (1978), 380                     */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */                                

int ffmin2(int mod)
{
    register int i,j,l;
    int err,rr,nn,ilo,ihi,ccount;
    int stepa,pstara,p2stara,pbara,ya,pa;
    double *step,*p,*pstar,*p2star,*pbar,*y;
    double dn,dnn,z,sum,summ,ylo,rcoeff,ystar,ecoeff,y2star;
    double ccoeff,del,ynewlo;
    double reduc = 0.001;
    
    stepa = pstara = p2stara = pbara = ya = pa = 0;

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value         Variance      Par Change   FCall\n");

    if (PMProtFDef == 1)
        fprintf(PMProtFd,"\n  Iter    Function Value         Variance      Par Change   FCall\n");

    err = -1;
    if (!(step = (double *) calloc(NParm + 1,sizeof(double))))  
        goto FM2End;
    memrq(NParm + 1,sizeof(double));
    stepa = NParm + 1;
   
    if (!(pstar = (double *) calloc(NParm + 1,sizeof(double))))  
        goto FM2End;
    memrq(NParm + 1,sizeof(double));
    pstara = NParm + 1;

    if (!(p2star = (double *) calloc(NParm + 1,sizeof(double))))  
        goto FM2End;
    memrq(NParm + 1,sizeof(double));
    p2stara = NParm + 1;

    if (!(pbar = (double *) calloc(NParm + 1,sizeof(double))))  
        goto FM2End;
    memrq(NParm + 1,sizeof(double));
    pbara = NParm + 1;

    if (!(y = (double *) calloc(NParm + 2,sizeof(double))))  
        goto FM2End;
    memrq(NParm + 2,sizeof(double));
    ya = NParm + 2;

    if (!(p = (double *) calloc((NParm + 1) * NParm + 1,sizeof(double))))  
        goto FM2End;
    pa = (NParm + 1) * NParm + 1;
    memrq(pa,sizeof(double));

    err = 0;

    for (i = 0; i <= NParm1; ++i) 
        step[i] = SLen;

    rcoeff = 1.0;
    ecoeff = 2.0;
    ccoeff = 0.5;
    ccount = CCheck;
    dn = NParm1;
    nn = NParm1 + 1;
    dnn = nn;
    del = 1.0;

FM2Res:
    for (i = 1; i <= NParm1; ++i)
        p[(i - 1) * nn + nn] = Par[i];

    fn(Par,mod,0,1,0,&rr);
    if (rr) {
        err = 1;
        goto FM2End;
    }
    z = FTmp; 

    y[nn] = z;
    sum = z;
    summ = z * z;
    for (j = 1; j <= NParm1; ++j) {
        Par[j] += step[j] * del;

        for (i = 1; i <= NParm1; ++i)
            p[(i - 1) * nn + j] = Par[i];

        fn(Par,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }
        z = FTmp; 

        y[j] = z;
        sum += z;
        summ += z * z;
        Par[j] -= step[j] * del;
    }

    while (1) {
        ylo = y[1];
        ynewlo = ylo;
        ilo = ihi = 1;

        for (i = 2; i <= nn; ++i) {
            if (y[i] < ylo) {
                ylo = y[i];
                ilo = i;
            }
            if (y[i] > ynewlo) {
                ynewlo = y[i];
                ihi = i;
            }
        }
        sum -= ynewlo;
        summ -= ynewlo * ynewlo;

        for (i = 1; i <= NParm1; ++i) {
            z = 0.0;
            for (j = 1; j <= nn; ++j)
                z += p[(i - 1) * nn + j];
     
            z -= p[(i - 1) * nn + ihi];
            pbar[i] = z / dn;
        }
        for (i = 1; i <= NParm1; ++i) 
           pstar[i] = (1.0 + rcoeff) * pbar[i] - rcoeff * p[(i - 1) * nn + ihi];
     
        fn(pstar,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        ystar = FTmp; 
        if (ystar >= ylo) 
            goto FM2L12;
     
        for (i = 1; i <= NParm1; ++i)  
            p2star[i] = ecoeff * pstar[i] + (1.0 - ecoeff) * pbar[i];

        fn(p2star,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }
        y2star = FTmp; 
        if (y2star >= ystar) {
            for (i = 1; i <= NParm1; ++i)
                p[(i - 1) * nn + ihi] = pstar[i];
            y[ihi] = ystar;
            sum += y[ihi];
            summ += y[ihi] * y[ihi];
            goto FM2L901;
        }
     
FM2L10:
        for (i = 1; i <= NParm1; ++i)
            p[(i - 1) * nn + ihi] = p2star[i];
        y[ihi] = y2star;
        sum += y[ihi];
        summ += y[ihi] * y[ihi];
        goto FM2L901;
     
FM2L12:
        l = 0;
        for (i = 1; i <= nn; ++i) {
            if (y[i] > ystar) 
                l++;
        }
        if (l > 1) {
            for (i = 1; i <= NParm1; ++i)
                p[(i - 1) * nn + ihi] = pstar[i];
            y[ihi] = ystar;
            sum += y[ihi];
            summ += y[ihi] * y[ihi];
        }
        else {
            if (l != 0) {
                for (i = 1; i <= NParm1; ++i)
                    p[(i - 1) * nn + ihi] = pstar[i];
                y[ihi] = ystar;
            }
            for (i = 1; i <= NParm1; ++i)
                p2star[i] = ccoeff * p[(i - 1) * nn + ihi] + (1.0 - ccoeff)
                                                                * pbar[i];
            fn(p2star,mod,0,1,0,&rr);
            if (rr) {
                err = 1;
                goto FM2End;
            }
            y2star = FTmp; 
            if (y2star <= y[ihi])
                goto FM2L10;
     
            sum = summ = 0.0;
            for (j = 1; j <= nn; ++j) {
                for (i = 1; i <= NParm1; ++i) {
                    p[(i - 1) * nn + j] = 
                           (p[(i - 1) * nn + j] + p[(i - 1) * nn + ilo]) * 0.5;
                    Par1[i] = p[(i - 1) * nn + j];
                }
                fn(Par1,mod,0,1,0,&rr);
                if (rr) {
                    err = 1;
                    goto FM2End;
                }
                y[j] = FTmp; 

                sum += y[j];
                summ += y[j] * y[j];
            }
        }
FM2L901:
        if (--ccount == 0) {
            Iter++;
            ccount = CCheck;
            CValV = fabs((summ - (sum * sum) / dnn) / dn);

            if (SILENTFlg < 2) {
                printfe("  %3d  %20.13e %17.10e ",Iter,FTmp,CValV);
                printfe("        --   %6d\n",NumF);
            }   
            if (PMProtFDef) {
                if (PMProtFDef > 1)
                    fprintf(PMProtFd,"\nITER      Function Value         Variance      Par Change   FCall\n");
                fprintf(PMProtFd,"  %3d  %20.13e %17.10e ",Iter,FTmp,CValV);
                fprintf(PMProtFd,"        --   %6d\n",NumF);
            }   
            if (Iter >= MxIter || CValV <= TOLV) 
                break;
        }
    }
    if (y[ihi] > y[ilo])
        ihi = ilo;

    for (i = 1; i <= NParm1; ++i)
        Par1[i] = p[(i - 1) * nn + ihi];

    FMin = ynewlo = y[ihi];

    if (Iter >= MxIter) {
        LConv = 2;
        goto FM2Fin;
    }
    for (i = 1; i <= NParm1; ++i) {     /* check final result */
        del = step[i] * reduc;
        Par1[i] += del;

        fn(Par1,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        z = FTmp; 
        if (ynewlo - z > TOLF)
            goto FM2Nxt;  
    
        Par1[i] = Par1[i] - del - del;

        fn(Par1,mod,0,1,0,&rr);
        if (rr) {
            err = 1;
            goto FM2End;
        }

        z = FTmp; 
        if (ynewlo - z > TOLF)
            goto FM2Nxt;  
    
        Par1[i] += del;
    }
    FMin = ynewlo;
    LConv = 0;

FM2Fin:
    for (i = 1; i <= NParm1; ++i)
        Par[i] = Par1[i];
    goto FM2End;

FM2Nxt:
    for (i = 1; i <= NParm1; ++i)
        Par[i] = Par1[i];
    del = reduc;
    if (SILENTFlg < 2)
        printfe("  Restart\n");
    if (PMProtFDef)
        fprintf(PMProtFd,"RESTART.\n");

    goto FM2Res;

FM2End:
    if (pa > 0) {
        free((char *)p);
        memrq(-pa,sizeof(double));
    }  
    if (ya > 0) {
        free((char *)y);
        memrq(-ya,sizeof(double));
    }
    if (pbara > 0) {
        free((char *)pbar);
        memrq(-pbara,sizeof(double));
    }   
    if (p2stara > 0) {
        free((char *)p2star);
        memrq(-p2stara,sizeof(double));
    }
    if (pstara > 0) {
        free((char *)pstar);
        memrq(-pstara,sizeof(double));
    }  
    if (stepa > 0) {
        free((char *)step);
        memrq(-stepa,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin3  Function minimization with conjugate gradient or BFGS           */
/*          Adapted and translated to C from: Algorithm 500 of Collected    */
/*          Algorithms from ACM., Trans. Math. Software 6 (1980), 618 - 22  */
/*          Written by D.F. Shanno and K.H. Phua.                           */
/*                                                                          */
/*          Note: most comments have been stripped. See the original        */
/*          Fortran source for a full documentation.                        */
/*                                                                          */
/*          Hess[] is used as a working area. It must have dimension of     */
/*          at least 5 * NParm + 2 if MINA = 3 (conjugate gradients)        */
/*          at least NParm * (NParm + 7) / 2 if MINA = 4 (BFGS).            */
/*                                                                          */
/*      Return:     0 if successful.                                        */
/*                 -1 if insufficient memory.                               */
/*                 +1 if error in function evaluation.                      */                                

int ffmin3(int mod)
{
    register int i,j,ii,ij;
    int nx,ng,nry,nrd,ncons,ncons1,ncons2,nrst,rsw,ifun,ncalls;
    int err,it1,nxpi,ngpi,nrdpi,nrypi,ngpj;
    double tmp,fm,fp,dal,dg1,xsq,gsq,alpha,at,ap,dp,dg,rtst,step,u1,u2,u3,u4;

    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    if (PMProtFDef == 1)
        fprintf(PMProtFd,"\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    alpha = 1.0;
    err = ifun = 0;
    nx = NParm1;
    ng = nx + NParm1;

    if (MINA == 3) {       /* conjugate gradient */
        nry = ng + NParm1;
        nrd = nry + NParm1;
        ncons = 5 * NParm1;
        ncons1 = ncons + 1;
        ncons2 = ncons + 2;
    }
    else
        ncons = 3 * NParm1;

FFL20:
    fn(Par,mod,1,1,0,&err);              /* initial function evaluation */
    if (err) 
        goto FM3Fin;

    FMin = FTmp;
    ifun++;

    nrst = NParm1;
    rsw = 1;
    xsq = dg1 = 0.0;
    for (i = 1; i <= NParm1; ++i) {
        Hess[i] = -Grad[i];
        xsq += Par[i] * Par[i];
        dg1 -= Grad[i] * Grad[i];
    }
    gsq = -dg1;

    dg = dg1;   /* inserted */

    CValG = 0.0;                     /* calculate norm of gradient */
    for (i = 1; i <= NParm1; ++i) {
        tmp = Grad[i];
        CValG += tmp * tmp;
    }
    CValG = sqrt(CValG);
    if (CValG <= TOLG)      /* check if initial point is minimizer */
        goto FM3Fin;

    /* entry point of major iteration loop */

    while (1) {

        if (++Iter > MxIter) {
            Iter--;
            LConv = 2;      /* exceeded max number of iterations (func calls) */
            goto FM3Fin;
        }
        if (SILENTFlg < 2) {
            printfe("  %3d  %20.13e %17.10e ",Iter,FMin,CValG);
            printfe("        --   %6d (%d,%d)\n",NumF,NumFG,NumFH);
        }
        if (PMProtFDef) {
            if (PMProtFDef > 1)
                fprintf(PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
            fprintf(PMProtFd,"  %3d  %20.13e %17.10e ",Iter,FMin,CValG);
            fprintf(PMProtFd,"        --   %6d (%d,%d)\n",NumF,NumFG,NumFH);

            if (PMProtFDef > 1) {
                prvec("\nParameter",NParm1,Par);
                prvec("Gradient",NParm1,Grad);
            } 
        }
        fm = FMin;
        ncalls = ifun;

        /* begin linear search */

        alpha *= dg / dg1;
    
        if (PMProtFDef > 1)
            fprintf(PMProtFd,"Begin line search: alpha = %g dg=%g dg1=%g\n",
                                               alpha,dg,dg1);
        if (nrst == 1 || MINA == 4)  
            alpha = 1.0;
        if (rsw)  
            alpha = 1.0 / sqrt(gsq);
        ap = 0.0;
        fp = fm;
        dg = dp = dg1;
    
        step = 0.0;
        for (i = 1; i <= NParm1; ++i) {
            step += Hess[i] * Hess[i];
            nxpi = nx + i;
            ngpi = ng + i;
            Hess[nxpi] = Par[i];
            Hess[ngpi] = Grad[i];
        }
        step = sqrt(step);
    
        /*  begin of linear search iterations */

        it1 = 0;
FFL80:
   
        if (PMProtFDef > 1)
            fprintf(PMProtFd,"Line search with alpha=%g step=%g\n",alpha,step);

        if (++it1 > MxIt1) {
            LConv = -5;     /* linear search failed */
            goto FM3Fin;
        }
        if (alpha * step <= SMin) {
            if (!rsw)
                goto FFL20;

            if (PMProtFDef > 1)  
                fprintf(PMProtFd,"Line search failed (smin=%g alpha=%g step=%g).\n",SMin,alpha,step);
                 
            LConv = -2;     /* linear search failed */
            goto FM3Fin;
        }
        for (i = 1; i <= NParm1; ++i) {
            nxpi = nx + i;
            Par[i] = Hess[nxpi] + alpha * Hess[i];
        }
        fn(Par,mod,1,1,0,&err);      /* new function evaluation */
        if (err) 
            goto FM3Fin;

        FMin = FTmp;
        ifun++;

        dal=0.0;
        for (i = 1; i <= NParm1; ++i) {
            dal += Grad[i] * Hess[i];
        }
        if (FMin > fm && dal < 0.0) {
            alpha /= 3.0;
            ap = 0.0;
            fp = fm;
            dp = dg;
            goto FFL80;
        }
        if (FMin > (fm + 0.0001 * alpha * dg) || fabs(dal / dg) > 0.9 ||
           ((ifun - ncalls) <= 1 && fabs(dal / dg) > EPSI1 && MINA == 3)) {

            u1 = dp + dal - 3.0 * (fp - FMin) / (ap - alpha);
            u2 = u1 * u1 - dp * dal;
            if (u2 <= 0.0)
                u2 = 0.0;
            else
                u2 = sqrt(u2);
            at = alpha - (alpha - ap) * (dal + u2 - u1) / (dal - dp + 2.0 * u2);
    
            if ((dal / dp) <= 0.0) {
                if (at < (1.01 * dmin(alpha,ap)) || at > (0.99 * dmax(alpha,ap)))
                    at = (alpha + ap) / 2.0;
            }
            else {
                if ((dal > 0.0 && 0.0 < at && at < (0.99 * dmin(alpha,ap))) ||
                    (dal <= 0.0 && at > (1.01 * dmax(alpha,ap)))) {
                    ;
                }
                else {
                    if (dal <= 0.0)  
                        at = 2.0 * dmax(alpha,ap);
                    else
                        at = dmin(alpha,ap) / 2.0;
                }
            }
            ap = alpha;
            fp = FMin;
            dp = dal;
            alpha = at;
            goto FFL80;
        }

        /*  line search was successful. test convergence */
  
        xsq = gsq = 0.0;
        for (i = 1; i <= NParm1; ++i) {
            gsq += Grad[i] * Grad[i];
            xsq += Par[i] * Par[i];
        }
        tmp = 0.0;                      /* calculate norm of gradient */
        for (i = 1; i <= NParm1; ++i)  
            tmp += Grad[i] * Grad[i];
        CValG = sqrt(tmp);
        if (CValG <= TOLG) {
            Iter++;
            goto FM3Prn;
        }
        for (i = 1; i <= NParm1; ++i) {
            Hess[i] *= alpha;
        }
        if (MINA == 3) {       /* conjugate gradients */
            rtst = 0.0;
            for (i = 1; i <= NParm1; ++i) {
                ngpi = ng + i;
                rtst += Grad[i] * Hess[ngpi];
            }
            if (fabs(rtst / gsq) > 0.2)
                nrst = NParm1;
    
            if (nrst == NParm1) {
                Hess[ncons + 1] = 0.0;
                Hess[ncons + 2] = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    ngpi = ng + i;
                    Hess[nrypi]  = Grad[i] - Hess[ngpi];
                    Hess[nrdpi]  = Hess[i];
                    Hess[ncons1] += Hess[nrypi] * Hess[nrypi];
                    Hess[ncons2] += Hess[i] * Hess[nrypi];
                }
            }
            u1 = u2 = 0.0;
            for (i = 1; i <= NParm1; ++i) {
                nrdpi = nrd + i;
                nrypi = nry + i;
                u1 = u1 - Hess[nrdpi] * Grad[i] / Hess[ncons1];
                u2 = u2 + Hess[nrdpi] * Grad[i] * 2.0 / Hess[ncons2] - 
                                               Hess[nrypi] * Grad[i] / Hess[ncons1];
            }
            u3 = Hess[ncons2] / Hess[ncons1];
            for (i = 1; i <= NParm1; ++i) {
                nxpi = nx + i;
                nrdpi = nrd + i;
                nrypi = nry + i;
                Hess[nxpi] = -u3 * Grad[i] - u1 * Hess[nrypi] - u2 * Hess[nrdpi];
            }
            if (nrst != NParm1) {
                u1 = u2 = u3 = u4 = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    u1 = u1 - (Grad[i] - Hess[ngpi]) * Hess[nrdpi] / Hess[ncons1];
                    u2 = u2 - (Grad[i] - Hess[ngpi]) * Hess[nrypi] / Hess[ncons1] +
                                      2.0 * Hess[nrdpi] * (Grad[i] - Hess[ngpi]) / Hess[ncons2];
                    u3 = u3 + Hess[i] * (Grad[i] - Hess[ngpi]);
                }
                step = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    ngpi = ng + i;
                    nrdpi = nrd + i;
                    nrypi = nry + i;
                    step = (Hess[ncons2] / Hess[ncons1]) * (Grad[i] - Hess[ngpi]) +
                                                  u1 * Hess[nrypi] + u2 * Hess[nrdpi];
                    u4 = u4 + step * (Grad[i] - Hess[ngpi]);
                    Hess[ngpi] = step;
                }
                u1 = u2 = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    u1 = u1 - Hess[i] * Grad[i] / u3;
                    ngpi = ng + i;
                    u2 = u2 + (1.0 + u4 / u3) * Hess[i] * Grad[i] / u3 - 
                                                         Hess[ngpi] * Grad[i] / u3;
                }
                for (i = 1; i <= NParm1; ++i) {
                    ngpi = ng + i;
                    nxpi = nx + i;
                    Hess[nxpi] = Hess[nxpi] - u1 * Hess[ngpi] - u2 * Hess[i];
                }
            }
            dg1 = 0.0;
            for (i = 1; i <= NParm1; ++i) {
                nxpi = nx + i;
                Hess[i] = Hess[nxpi];
                dg1 += Hess[i] * Grad[i];
            }
            if (dg1 > 0.0) {
                LConv = -4; /* search vector not a descent direction */
                goto FM3Fin;
            }
            if (nrst == NParm1)
                nrst = 0;
            nrst++;
            rsw = 0;
        }
        else {  /* BFGS */
            u1 = 0.0;
            for (i = 1; i <= NParm1; ++i) {
                ngpi = ng + i;
                Hess[ngpi] = Grad[i] - Hess[ngpi];
                u1 += Hess[i] * Hess[ngpi];
            }
            if (rsw) {
                u2 = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    ngpi = ng + i;
                    u2 += Hess[ngpi] * Hess[ngpi];
                }
                ij = 1;
                u3 = u1 / u2;
                for (i = 1; i <= NParm1; ++i) {
                    for (j = i; j <= NParm1; ++j) {
                        ncons1 = ncons + ij;
                        Hess[ncons1] = 0.0;
                        if (i == j)
                            Hess[ncons1] = u3;
                        ij++;
                    }
                    nxpi = nx + i;
                    ngpi = ng + i;
                    Hess[nxpi] = u3 * Hess[ngpi];
                }
                u2 *= u3;
            }
            else {
                u2 = 0.0;
                for (i = 1; i <= NParm1; ++i) {
                    u3 = 0.0;
                    ij = i;
                    if (i != 1) {
                        ii = i - 1;
                        for (j = 1; j <= ii; ++j) {
                            ngpj = ng + j;
                            ncons1 = ncons + ij;
                            u3 += Hess[ncons1] * Hess[ngpj];
                            ij = ij + NParm1 - j;
                        }
                    }
                    for (j = i; j <= NParm1; ++j) {
                        ncons1 = ncons + ij;
                        ngpj = ng + j;
                        u3 += Hess[ncons1] * Hess[ngpj];
                        ij++;
                    }
                    ngpi = ng + i;
                    u2 += u3 * Hess[ngpi];
                    nxpi = nx + i;
                    Hess[nxpi] = u3;
                }
            }
            u4 = 1.0 + u2 / u1;
            for (i = 1; i <= NParm1; ++i) {
                nxpi = nx + i;
                ngpi = ng + i;
                Hess[ngpi] = u4 * Hess[i] - Hess[nxpi];
            }
            ij = 1;
    
            for (i = 1; i <= NParm1; ++i) {
                nxpi = nx + i;
                u3 = Hess[i] / u1;
                u4 = Hess[nxpi] / u1;
                for (j = i; j <= NParm1; ++j) {
                    ncons1 = ncons + ij;
                    ngpj = ng + j;
                    Hess[ncons1] = Hess[ncons1] + u3 * Hess[ngpj] - u4 * Hess[j];
                    ij++;
                }
            }
            dg1 = 0.0;
            for (i = 1; i <= NParm1; ++i) {
                u3 = 0.0;
                ij = i;
                if (i != 1) {
                    ii = i - 1;
                    for (j = 1; j <= ii; ++j) {
                        ncons1 = ncons + ij;
                        u3 -= Hess[ncons1] * Grad[j];
                        ij = ij + NParm1 - j;
                    }
                }
                for (j = i; j <= NParm1; ++j) {
                    ncons1 = ncons + ij;
                    u3 -= Hess[ncons1] * Grad[j];
                    ij++;   
                }
                dg1 += u3 * Grad[i];
                Hess[i] = u3;
            }
            if (dg1 > 0.0) {
                LConv = -4; /* search vector not a descent direction */
                goto FM3Fin;
            }
            rsw = 0;  
        }
    }
FM3Prn:
    if (SILENTFlg < 2) {
        printfe("  %3d  %20.13e %17.10e ",Iter,FMin,CValG);
        printfe("        --   %6d (%d,%d)\n",NumF,NumFG,NumFH);
    }
    if (PMProtFDef) {
        if (PMProtFDef > 1)
            fprintf(PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
        fprintf(PMProtFd,"  %3d  %20.13e %17.10e ",Iter,FMin,CValG);
        fprintf(PMProtFd,"        --   %6d (%d,%d)\n",NumF,NumFG,NumFH);

        if (PMProtFDef > 1) {
            prvec("\nParameter",NParm1,Par);
            prvec("Gradient",NParm1,Grad);
        } 
    }

FM3Fin:
    if (err)
        err = 1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  ffmin5  Function minimization with a modified Newton method.            */
/*          If MINA == 5, modification is with gradient methods, if         */
/*          MINA == 6 use additional information of max eigenvector.        */
/*                                                                          */
/*          Return -1 if insufficient memory, otherwise 0.                  */

int ffmin5(int mod)
{
    register int i,j,k;
    int err,npflg,ssflg;
    double step,test,test1,test2,lval,tmp,emin,*ev,*evec;

    err = 0;
    lval = 0.0;
  
    if (SILENTFlg < 2)
        printfe("\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");

    if (PMProtFDef == 1)  
        fprintf(PMProtFd,"\n  Iter    Function Value     Norm of Gradient  Par Change   FCall\n");
      
    LConv = 2;      /* return status of this function */
    npflg = 0;      /* set if hessian not positive definite and we use the  */
                    /* modified algorithm.                                  */

                             /* first call of function evaluation, including   */
    fn(Par,mod,2,1,0,&err);  /* derivatives in Grad, Hess, Diag. The function  */
    if (err) 
        goto FM5Fin;

    FMin = FTmp;        /* value is returned in FTmp.                     */

    while (++Iter <= MxIter) {   /* Main loop over iterations */
        /******************
        if (FTmp <= 0.0)  
            NPFlgs[3] += 1;
        *********************/
        tmp = 0.0;                      /* calculate norm of gradient */
        for (j = 1; j <= NParm1; ++j)  
            tmp += Grad[j] * Grad[j];
        CValG = sqrt(tmp);

        if (Iter > 1) { 
            CValF = fabs(FMin - lval);
            if (lval)
                CValF /= lval;

            CValP = 0.0;
            for (j = 1; j <= NParm1; ++j) {
                tmp = fabs(Par[j] - Par1[j]);
                if (Par[j])
                    tmp /= fabs(Par[j]);
                if (CValP < tmp)
                    CValP = tmp;
            }
            if (PMProtFDef > 1) {
                prval("Change in function value",CValF);
                prval("Change in parameter estimates",CValP);
            }
        }
        if (SILENTFlg < 2)
            printfe("  %3d  %20.13e %17.10e ",Iter,FMin,CValG);

        if (PMProtFDef) {
            if (PMProtFDef > 1)  
                fprintf(PMProtFd,"\nITER      Function Value     Norm of Gradient  Par Change   FCall\n");
            fprintf(PMProtFd,"  %3d  %20.13e %17.10e ",Iter,FMin,CValG);
        }
        if (Iter == 1) {
            if (SILENTFlg < 2)
                printfe("         --  %6d (%d,%d)",NumF,NumFG,NumFH);
            if (PMProtFDef)
                fprintf(PMProtFd,"         --  %6d (%d,%d)\n",NumF,NumFG,NumFH);
        }
        else {
            if (SILENTFlg < 2)
                printfe("%11.4e  %6d (%d,%d)",CValP,NumF,NumFG,NumFH);
            if (PMProtFDef)
                fprintf(PMProtFd,"%11.4e  %6d (%d,%d)\n",CValP,NumF,NumFG,NumFH);
        }
        if (PMProtFDef > 1) {
            prvec("\nParameter",NParm1,Par);
            prvec("Gradient",NParm1,Grad);
            prhess("Hessian",NParm1);
        } 

        /*  Check for convergence. If Crite > 1 begin with second iteration */

        if (Crite == 1 || Iter > 1) { 

            switch (Crite) {
                case 1:     if (fabs(CValG) <= TOLG) LConv = 0;
                            break;
                case 2:     if (fabs(CValF) <= TOLF) LConv = 0;
                            break;
                case 3:     if (fabs(CValP) <= TOLP) LConv = 0;
                            break;
                default:    break;
            }
            if (!CritFlg && Iter > 1 && (fabs(CValG) <= TOLG || fabs(CValF) <= TOLF))
                LConv = 0;

            if (!LConv)                 /* convergence was reached.         */
                goto FM5Fin;

            /****************************************************************
            if (Iter > 1 && fabs(CValF) <= TOLF) {  check for minimum decrease 
                                                    of the function value.     
                LConv = 1;
                goto FM5Fin;
            }
            *****************************************************************/
        }
        lval = FMin;                        /* save actual function value   */
        for (j = 1; j <= NParm1; ++j)       /* and parameters.              */
            Par1[j] = Par[j];
        
        if (npflg) {     /* we need a copy of the hessian in this case */

            for (j = 1; j <= NParm1; ++j) {  
                evec[(j - 1) * NParm1 + j] = Diag[j];
                for (k = 1; k < j; ++k)  
                    evec[(j - 1) * NParm1 + k] = 
                    evec[(k - 1) * NParm1 + j] = Hess[(j - 1) * (j - 2) / 2 + k];
            }
        }

        /*  Try to calculate the Newton search direction by a Cholesky      */
        /*  decomposition of the hessian used to solve Hess * Srch = -Grad. */
        /*  This is done by chinv(). If the return is zero, is was success- */
        /*  full and Srch contains the Newton search vector. Otherwise the  */
        /*  hessian is not positive definite.                               */

        if (!chinv(1)) {

            /*  The hessian is positive definite. Srch contains now the     */
            /*  Newton search vector. If constraints, Srch is the search    */
            /*  vector in the reduced parameter space.                      */

            if (SILENTFlg < 2)
                printfe("\n");
            if (PMProtFDef > 1)
                prvec("Search vector",NParm1,Srch);
            ssflg = 1;
        }
        else {  /*  The hessian is not positive definite. */

            NPFlgs[4] += 1;
            if (SILENTFlg < 2)
                printfe("  (npd)\n");
            if (PMProtFDef > 1)
                fprintf(PMProtFd,"Hessian not positive definite.\n");
  
            /*  First part of the search vector is the negative gradient.   */

            for (j = 1; j <= NParm1; ++j)    
                Srch[j] = -Grad[j];

            if (MINA == 5) {

                if (PMProtFDef > 1)
                    prvec("Search vector",NParm1,Srch);
  
                ssflg = 1;
            }
            else {

                ssflg = 0;

                /*  The second part is the eigenvector of the hessian       */
                /*  associated with the minimal eigenvalue.                 */

                if (!npflg) {   /* we need some additional storage in this case */

                    if (!(ev = (double *)calloc(NParm1 + 1,sizeof(double))))  
                        return(-1);

                    if (!(evec = (double *)calloc(NParm1 * NParm1 + 1,sizeof(double)))) {
                        free((char *)ev);
                        return(-1);
                    }
                    memrq(NParm1 + NParm1 * NParm1 + 2,sizeof(double));
 
                    fn(Par,mod,2,1,0,&err);      
                    if (err) 
                        goto FM5Fin;

                    for (j = 1; j <= NParm1; ++j) {  
                        evec[(j - 1) * NParm1 + j] = Diag[j];
                        for (k = 1; k < j; ++k)  
                            evec[(j - 1) * NParm1 + k] = 
                            evec[(k - 1) * NParm1 + j] = Hess[(j - 1) * (j - 2) / 2 + k];
                    }
                    npflg = 1;
                }

                /*  Calculate eigenvalues and eigenvectors of the hessian */

                if (evecf(NParm1,ev,evec)) {
                    LConv = -3;
                    goto FM5Fin;
                }

                /*  determine the minimum eigenvalue emin, and the associated   */
                /*  eigenvector, stored again in evec.                          */

                k = 1;
                emin = ev[1];
                for (j = 2; j <= NParm1; ++j) {
                    if (emin > ev[j]) {
                        emin = ev[j];
                        k = j;
                    }
                }
                i = 1;
                for (j = 1; j <= NParm1; ++j) 
                    evec[i++] = evec[(j - 1) * NParm1 + k];

                if (PMProtFDef > 1) {
                    prval("Minimum eigenvalue",emin);
                    prvec("Associated eigenvector",NParm1,evec);
                    fprintf(PMProtFd,"Step size search.\n");
                }
  
                /*  Now select the sign so that evec is a nonascent direction   */
                /*  with negative curvature.                                    */

                test1 = 0.0;
                for (j = 1; j <= NParm1; ++j)  
                    test1 += evec[j] * Grad[j];

                if (test1 > 0.0) {
                    for (j = 1; j <= NParm1; ++j)
                        evec[j] = -evec[j];
                }

                /*  The step size search is one with the second-order Armijo    */
                /*  condition. Cf. McCormick 1983, p. 135.                      */

                test1 = test2 = 0.0;
                for (j = 1; j <= NParm1; ++j) {
                    test1 += Grad[j] * Srch[j];
                    test2 += evec[j] * evec[j];
                }
                test2 *= emin;
                test = AMue * (test1 + 0.5 * test2);
                step = 1.0;

                while (1) {
                    for (j = 1; j <= NParm1; ++j)
                        Par[j] = Par1[j] + step * Srch[j] + sqrt(step) * evec[j];


                    if (STFlg)          /* depending on STFlg the calculation   */
                        fn(Par,mod,0,1,0,&err);      /* includes derivatives or not.         */
                    else  
                        fn(Par,mod,2,1,0,&err);
                    if (err) 
                        goto FM5Fin;

                    if (PMProtFDef > 1)
                        prval("Function value",FTmp);
   
                    /*******************
                    if (FTmp <= 0.0)  
                        NPFlgs[3] += 1;
                    *********************/
                    if (FTmp - FMin <= test * step)
                        break;

                    /* check for minimum change in function value */
                    if (fabs(FTmp - FMin) <= TOLF) {   
                        break;
                        /*********** changed 
                        LConv = 1;  
                        goto FM5Fin;
                        ********************/
                    }
                    if ((step *= 0.5) < SMin) {    /* check for minimal step size  */
                        LConv = -2;
                        goto FM5Fin;
                    }
                }
            }
        }
        if (ssflg) {        /* standard line search */

            if (PMProtFDef > 1)
                fprintf(PMProtFd,"Step size search.\n");
  
            /*  Calculation of a suitable step size. We use the first-      */
            /*  order Armijo condition, cf. McCormick 1983, p. 134.         */

            test = 0.0;
            for (j = 1; j <= NParm1; ++j)       /* test criterion in the    */
                test += Grad[j] * Srch[j];      /* reduced space.           */

            test *= AMue;
            step = 1.0;

            while (1) {
                for (j = 1; j <= NParm1; ++j)
                    Par[j] = Par1[j] + step * Srch[j];

                if (STFlg)                    /* depending on STFlg the calculation   */
                    fn(Par,mod,0,1,0,&err);   /* includes derivatives or not.         */
                else  
                    fn(Par,mod,2,1,0,&err);
                if (err) 
                    goto FM5Fin;

                if (PMProtFDef > 1)
                    prval("Function value",FTmp);
  
                /******************
                if (FTmp <= 0.0)  
                    NPFlgs[3] += 1;
                *********************/
                if (FTmp - FMin <= test * step)
                    break;
   
                if (PMProtFDef > 1) {
                    prval("Change in function value",FTmp - FMin);
                    prval("Armijo test criterion",test * step);
                }
      
                /* check for minimum change in function value */
                if (fabs(FTmp - FMin) <= TOLF) {   
                    break;
                    /************* changed 
                    LConv = 1;  
                    goto FM5Fin;
                    *********************/
                }
                if ((step *= 0.5) < SMin) {    /* check for minimal step size  */
                    LConv = -2;
                    goto FM5Fin; 
                }
            }
        }
        if (STFlg) {                 /* in this case we need an additional function  */
            fn(Par,mod,2,1,0,&err);  /* evaluation including derivatives.            */
            if (err) 
                goto FM5Fin;
        }
        FMin = FTmp;    /* set new function value */

        if (PMProtFDef > 1)
            prval("Step size",step);
    }
    LConv = 2;                /* exceeded max number of iterations */
    Iter--;

FM5Fin:
    if (npflg) {
        memrq(NParm1 + NParm1 * NParm1 + 2,-sizeof(double));
        free((char *)ev);
        free((char *)evec);
    }
    if (err)
        err = 1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  exphess   Expand the Hessian matrix (given by Hess and Diag) to NParm   */
/*            by premultiplication with ConQ and by postmultiplication      */
/*            with the transpose of ConQ.                                   */
/*            Use of WrkD and WrkH as temporary storage.                    */

void exphess(void)
{
    register int i,j,k,l;
    double tmp;

    for (i = 1; i <= NParm; ++i) {
        for (j = 1; j <= i; ++j) {
            tmp = 0.0;
            for (k = 1; k <= NParm1; ++k) {
                for (l = 1; l <= NParm1; ++l) {
                    if (l == k)
                        tmp += Diag[k] * ConQ[(CIP[i] - 1) * NParm1 + l] * 
                                         ConQ[(CIP[j] - 1) * NParm1 + k];
                    else if (l < k)
                        tmp += Hess[(k - 1) * (k - 2) / 2 + l] *
                                         ConQ[(CIP[i] - 1) * NParm1 + l] * 
                                         ConQ[(CIP[j] - 1) * NParm1 + k];
                    else
                        tmp += Hess[(l - 1) * (l - 2) / 2 + k] *
                                         ConQ[(CIP[i] - 1) * NParm1 + l] * 
                                         ConQ[(CIP[j] - 1) * NParm1 + k];
                }
            }
            if (i == j)
                WrkD[i] = tmp;
            else
                WrkH[(i - 1) * (i - 2) / 2 + j] = tmp;
        }
    }
    for (i = 1; i <= NParm; ++i)  
        Diag[i] = WrkD[i];
    for (i = 1; i <= HSiz; ++i)
        Hess[i] = WrkH[i];
}

/* ------------------------------------------------------------------------ */
/*  chinv(opt)                                                              */
/*      If opt = 1, this function solves Hess * Srch = -Grad. First, a      */
/*      Cholesky decomposition of Hess is performed and the result is       */
/*      stored again in Hess (and Diag). If Hess is not positive definite,  */
/*      the function return with an positive index. Else the search vector  */
/*      is calculated by solving Hess * Srch = -Grad, and the function      */
/*      returns zero. If opt = 0, it is only checked if the hessian matrix  */
/*      is positive definite.                                               */
/*                                                                          */
/*      Return 0 if positive definite, else 1.                              */

int chinv(int opt)
{
    register int i,j,k,l,m;
    double tmp,tmp1;

    for (i = 1; i <= NParm1; ++i) {
        m = (i - 1) * (i - 2) / 2;
        tmp1 = Diag[i];
        for (j = 1; j < i; ++j) {
            l = (j - 1) * (j - 2) / 2;
            tmp = Hess[m + j];
            for (k = 1; k < j; ++k)  
                tmp -= Hess[m + k] * Hess[++l];
            tmp /= Diag[j];
            Hess[m + j] = tmp;
            tmp1 -= tmp * tmp;
        }
        if (tmp1 < EPSI1)
            return(i);
        Diag[i] = sqrt(tmp1);
    }
    if (!opt)
        return(0);

    /* now solution of Hess * Srch = -Grad */

    for (i = 1; i <= NParm1; ++i) {      /* forward substitution */
        k = (i - 1) * (i - 2) / 2;
        tmp = -Grad[i];
        for (j = 1; j < i; ++j)
            tmp -= Hess[++k] * WrkD[j];
        WrkD[i] = tmp / Diag[i];
    }
    for (i = NParm1; i >= 1; --i) {      /* backwards substitution */
        tmp = WrkD[i];
        for (j = i + 1; j <= NParm1; ++j) {
            k = (j - 1) * (j - 2) / 2;
            tmp -= Hess[k + i] * Srch[j];
        }
        Srch[i] = tmp / Diag[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ffmin7  Function minimization with CES (quadratic or tensor model).     */
/*                                                                          */
/*  Return 0 if OK, -1 if insufficient memory, 1 if error in function       */

int ffmin7(int mod)
{
    register int i,j,k;
    int err,r,wa,na,*acn;
    double *typx,*acw,*acy,*acu,*acv,*w;

    err = -1;
    na = wa = 0;

    if (!(w = (double *)calloc(5 * NParm1 + 15,sizeof(double))))   
        goto FFMIN7Fin;

    wa = 5 * NParm1 + 15;
    memrq(wa,sizeof(double));

    typx = w;
    acw  = typx + NParm1 + 1;
    acy  = acw  + NParm1 + 1;
    acu  = acy  + NParm1 + 1;
    acv  = acu  + NParm1 + 1;

    if (!(acn = (int *)calloc(NParm1 + 1,sizeof(int))))   
        goto FFMIN7Fin;

    na = NParm1 + 1;
    memrq(na,sizeof(int));

    for (i = 1; i <= NParm1; ++i)
        typx[i] = 1.0;

    TSFNC = 0;      /* use standard models */
    TSMOD = mod;
    if (MINA == 7)
        TSMETH = 0;
    else
        TSMETH = 1;

    r = ts_min(NParm1,Par,typx,Par1,GTmp,TSGrad,Srch,acw,acy,WrkE,acu,acv,
               TSHess,acn);

    if (r < 0) {            /* error in function evaluation */
        err = 1;
        goto FFMIN7Fin;
    }
    switch (r) {
        case 1:
        case 2:     LConv = 0;
                    break;
        case 3:
        case 5:     LConv = -2;
                    break;
        case 4:     LConv = 2;
                    break;
        default:    LConv = 1;
                    break;
    }

    /* if we only have function values then put gradient into Grad[] 
       and Hessian into Hess[] and Diag[] */

    if (PMDOPT == 0) {
        k = 1;
        for (i = 1; i <= NParm1; ++i) {
            for (j = 1; j < i; ++j)
                Hess[k++] = TSHess[(i - 1) * NParm1 + j] / DScal;
            Diag[i] = TSHess[(i - 1) * NParm1 + i] / DScal;
            Grad[i] = TSGrad[i] / DScal;
        }
    }
    err = 0;

FFMIN7Fin:

    if (wa > 0) {
        free((char *)w);
        memrq(-wa,sizeof(double));
    }
    if (na > 0) {
        free((char *)acn);
        memrq(-na,sizeof(int));
    }
    return(err);
}



