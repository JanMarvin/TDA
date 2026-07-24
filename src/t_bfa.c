/****************************************************************************/
/*  t_bfa                                                                   */
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
#include "t_gdat.h"
#include "t_alloc.h"
#include "t_var.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_freq.h"
#include "t_con.h"


/*  functions in t_bfa.c */

int bfa(void);
int bfa_a(int mmax);
int bfa_free(void);
void bfa_p(unsigned int x,int n,int m);
void bfa_ps(unsigned int x,unsigned int y,int n,int m,int opt);
void bfa_ps_p(unsigned int x,unsigned int y,int n,int m);
int bfa_len(unsigned int x,int n,int m);
int bfa_find(int m,int n,int maximp);
int bfa_elim(int m,int n,int nip);
int bfa_chart(int m,int n,int nip,int maximp);
int bfa_ess(int nip,int ny1); 
int bfa_red(int nip,int ny1,int *ni,int *nj);
int bfa_pm(int ni,int nj,int maximp);
int bfa_minc(int n,char *x, char *y);
int bfa_mince(int n,char *x, char *y);
int bfa_icheck(int n,char *x, char *y);
void bfa_copy(int n,char *x, char *y);
int bfa_fcheck(unsigned int x,int n,int m,int nei,int ni,char *sel);
void bfa_prn(int nei,int nf,int ni,int n,int m,int opt);
int bfa_lm(int ni,int nj);
int bfa_lm1(int ni,int nj);
int bfc(void);
void bfa_ps_n(unsigned int x,unsigned int y,int n,int m,int jv);
                                    

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define IMPMax 50000
#define BMMax 15            /* maximum number of arguments                  */ 

unsigned int *BMImplA;      /* table with implicants                        */
int BMImplA_A = 0;          /* allocated                                    */
unsigned int *BMImplA1; 
int BMImplA1_A = 0;
unsigned int *BMImplB;  
int BMImplB_A = 0;
unsigned int *BMImplB1; 
int BMImplB1_A = 0;
unsigned int *BMImplP;      /* prime implicants                             */
int BMImplP_A = 0;
unsigned int *BMImplP1; 
int BMImplP1_A = 0;
unsigned int *BMImplE;      /* essential prime implicants                   */
int BMImplE_A = 0;
unsigned int *BMImplE1; 
int BMImplE1_A = 0;

char *BMImplX;              /* indicates implicants                         */
int BMImplX_A = 0;      
char *BMImplC;              /* prime implicant chart                        */
int BMImplC_A = 0;      
char *BMImplS;              /* selections of implicants                     */
int BMImplS_A = 0;      
short *BMImplL;             /* index                                        */
int BMImplL_A = 0;      


/* ------------------------------------------------------------------------ */
/*  bfa     Boolean function analysis                                       */
/*                                                                          */
/*  bfa(                                                                    */
/*      opt=...,            treatment of undefined arguments, def. 1        */
/*                          1 treated as don't cares                        */
/*                          2 treated as logical 1                          */
/*                          3 treated as logical 0                          */
/*      alg=...,            algorithm for minimization, def. 0              */
/*                          0 not performed                                 */  
/*                          1 = Lawler's algorithm                          */
/*                          2 = Petrick's method                            */
/*                          3 = Lawler's algorithm II                       */
/*      ptab=...,           print table with input data                     */
/*      ptab1=...,          print table with minimized boolean functions    */
/*      prot=...,           protocol file with additional information       */
/*      max=...,            max number of implicants, def. 10000            */
/*                                                                          */
/*  ) = Y,X1,...,Xm;        variables for Boolean function                  */
/*                          (Y is the dependent variable)                   */
/*                                                                          */
/*  Boolean function: Y = f(X1,...,Xm). All variables are interepreted      */
/*  as binary variables. A value not equal to zero will be interpreted      */
/*  as a logical 1.                                                         */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int bfa(void)
{
    register int i,j,k,l,p;
    int err,maximp,m,n,ii,ny,ny0,ny1,ny01,d,ma,fnd,nip,ni,nj;
    int nei,nni,nf,mlen,nsel;
    unsigned int ui,ai,aj,ax,xx,am;

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Boolean function analysis. Current memory: %d bytes.\n\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto BFAFin;

    if (PMALG > 3)
        PMALG = 0;

    if (PMMax > 0)             
        maximp = PMMax;
    else
        maximp = IMPMax;

    m = PMNV - 1;
    if (m < 1) {
        printf1("Error: need at least two variables.\n");
        goto BFAFin;
    }
    else if (m > BMMax) {
        printf1("Error: maximum number of arguments is %d.\n",BMMax);
        goto BFAFin;
    }
    n = (int)pow(2.0,(double)m);

    printf1("Number of arguments: %d\n",m);
    printf1("Rows of truth table: %d\n",n);

    if (alloc_acn(n + 1))                 /* if Y=1 */
        goto BFAFin;
    if (alloc_acm(n + 1))                 /* if Y=0 */
        goto BFAFin;
                           
    for (i = 0; i < NOC; ++i) {
        ii = 0;
        k  = 1;
        for (j = PMNV - 1; j >= 1; --j) {
            if ((int)get_data(PMVIdx[j],i))
                ii += k;
            k *= 2;
        }
        if ((int)get_data(PMVIdx[0],i))
            AcN[ii] += 1;
        else   
            AcM[ii] += 1;
    }
    ny = ny0 = ny1 = ny01 = 0;
    for (i = 0; i < n; ++i) {
        if (AcM[i] > 0 && AcN[i] > 0)
            ny01++;
        else if (AcM[i] > 0)
            ny0++;
        else if (AcN[i] > 0)
            ny1++;
        else
            ny++;
    }
    printf1("- with y = 0 : %d\n",ny0);
    printf1("- with y = 1 : %d\n",ny1);
    printf1("- both       : %d\n",ny01);
    printf1("- undefined  : %d\n\n",ny);

    if (ny > 0) {
        printf1("Undefined truth table rows treated as ");
        if (PMOPT == 2) 
            printf1("logical 1.\n"); 
        else if (PMOPT == 3) 
            printf1("logical 0.\n"); 
        else {                  
            PMOPT = 1;
            printf1("don't cares.\n"); 
        }      
        if (PMOPT == 2 || PMOPT == 3) {
            for (i = 0; i < n; ++i) {
                if (AcM[i] == 0 && AcN[i] == 0) {
                    if (PMOPT == 2)  
                        AcN[i] += 1;
                    else               
                        AcM[i] += 1;
                }
            }
            if (PMOPT == 2)
                ny1 += ny;
            else
                ny0 += ny;
            ny = 0;
        }
    }

    /* -------------------------------------------------------------------- */
    /* If requested print the truth table to an output file.                */

    if (PMTabFDef) {            /* print table */ 
        for (i = 0; i < n; ++i) {
            fprintf(PMTabFd,"%10d  ",i);

            ui = n;
            for (j = 0; j < m; ++j) {
                ui = ui >> 1;
                if (ui & i)
                    fprintf(PMTabFd,"1 ");
                else
                    fprintf(PMTabFd,"0 ");
            }
            fprintf(PMTabFd,"%3d %3d\n",AcM[i],AcN[i]);
        }
        printf1("Table written to: %s\n",PMTabFName);
    }
    if (PMProtFDef)
        printf1("Protocol file: %s\n",PMProtFName);
    newline();

    if (ny1 < 1) {
        printf1("Number of true minterms is zero.\n");
        err = 0;
        goto BFAFin;
    }
    if (ny01 > 0) {
        printf1("Minimization requires unique function values.\n");
        err = 0;
        goto BFAFin;
    }

    if (bfa_a(maximp))              /* allocate memory */
        goto BFAFin;

    nei = nni = ni = nj = nf = 0;

    /* -------------------------------------------------------------------- */
    /* Find all prime implicants and store in BMImplP.                      */
    /* nip is the number of prime implicants.                               */

    nip = bfa_find(m,n,maximp);                                  
    if (nip < 0) {
        err = nip;
        goto BFAFin;
    }

    /* If table is incomplete and PMOPT = 1 eliminate redundant
       implicants */

    if (ny > 0 && PMOPT == 1)   
        nip = bfa_elim(m,n,nip);

    /* print prime implicants */

    printf1("Found %d prime implicants.\n",nip);
    for (i = 0; i < nip; ++i) {
        printf1("%5d : ",i);
        bfa_ps(BMImplP[i],BMImplP1[i],n,m,0);
        printf1("    ");
        bfa_ps(BMImplP[i],BMImplP1[i],n,m,1);
        newline();                 
    }
    newline();
        
    if (PMALG == 0) {   /* nothing more done if alg = 0 */
        err = 0;
        goto BFAFin;
    }

    /* -------------------------------------------------------------------- */
    /* Make prime implicant chart in BMImplC, size is nip x ny1             */

    ma = nip * ny1;
    ni = nip;
    nj = ny1;

    printf1("Allocating prime implicant chart (%d x %d = %d bytes).\n",nip,ny1,ma);
    if (!(BMImplC = (char *)calloc(ma,sizeof(char)))) { 
        p_err(-2,1);
        goto BFAFin;
    }
    BMImplC_A = ma;           
    memrq(ma,sizeof(char));

    if ((err = bfa_chart(m,n,nip,maximp)) < 0)  
        goto BFAFin;
        
    /* -------------------------------------------------------------------- */
    /* Find essential prime implicants. Mark them by a 1 in BMImplX.        */
    /* The number of essential prime implicants is nei.                     */

    nei = bfa_ess(nip,ny1);

    printf1("Found %d essential prime implicants.\n",nei);

    if (PMProtFDef) {
        fprintf(PMProtFd,"Prime implicants: %d. Essential: %d.\n",nip,nei);
        fprintf(PMProtFd,"Prime implicant chart (%d x %d).\n",nip,ny1);

        for (i = 0; i < nip; ++i) {
            fprintf(PMProtFd,"%5d : ",i);
            bfa_ps_p(BMImplP[i],BMImplP1[i],n,m);
            fprintf(PMProtFd," : %d : ",BMImplX[i]);
            for (j = 0; j < ny1; ++j)  
                fprintf(PMProtFd,"%d ",BMImplC[i * ny1 + j]);
            fprintf(PMProtFd,"\n");
        }
    }

    /* -------------------------------------------------------------------- */
    /* Store essential prime implicants in BMImplE.                         */

    if (nei > 0) {
        if (!(BMImplE = (unsigned int *)calloc(nei + 1,sizeof(int)))) { 
            p_err(-2,1);
            goto BFAFin;
        }
        BMImplE_A = nei + 1;        
        memrq(nei + 1,sizeof(int));

        if (!(BMImplE1 = (unsigned int *)calloc(nei + 1,sizeof(int)))) { 
            p_err(-2,1);
            goto BFAFin;
        }
        BMImplE1_A = nei + 1;           
        memrq(nei + 1,sizeof(int));

        j = 0;
        for (i = 0; i < nip; ++i) {
            if (BMImplX[i]) {
                BMImplE[j] = BMImplP[i];
                BMImplE1[j] = BMImplP1[i];
                j++;
            }
        }
    }

    /* -------------------------------------------------------------------- */
    /* Store non-essential prime implicants in BMImplP.                     */
    /* The number of non-essential prime implicants is nni.                 */

    nni = 0;
    for (i = 0; i < nip; ++i) {
        if (BMImplX[i] == 0) {
            BMImplP[nni] = BMImplP[i];
            BMImplP1[nni] = BMImplP1[i];
            nni++;
        }
    }
    if (nni <= 0)
        goto BFAFin1;

    /* -------------------------------------------------------------------- */
    /* If nei > 0 make reduced table with only the non-essential prime      */ 
    /* implicants.                                                          */
       
    if (nei > 0) {                                              
        if (bfa_red(nip,ny1,&ni,&nj))                                  
            goto BFAFin;
          
        printf1("Removing essential implicants; reduced chart: %d x %d\n",ni,nj);

        if (PMProtFDef) {
            fprintf(PMProtFd,"\nReduced chart (%d x %d).\n",ni,nj);
            if (ni > 0 && nj > 0) {
                for (i = 0; i < ni; ++i) {
                    fprintf(PMProtFd,"%5d : ",i);
                    bfa_ps_p(BMImplP[i],BMImplP1[i],n,m);
                    fprintf(PMProtFd," : ");

                    for (j = 0; j < nj; ++j)  
                        fprintf(PMProtFd,"%d ",BMImplC[i * nj + j]);
                    fprintf(PMProtFd,"\n");
                }
            }
        }
    }
    if (ni != nni) {  
        printf1("Error. Can't continue.\n");
        goto BFAFin;
    }
    if (ni <= 0 || nj <= 0)    
        goto BFAFin1;

    /* -------------------------------------------------------------------- */
    /* Find a minimal set of the non-essential implicants.                  */
    /* If alg == 0 use Lawler's method                                      */
    /*           1 use Petrick's method.                                    */  
    /*                                                                      */
    /* Store the covers in the rows of BMImplS; length in BMImplL.          */
    /* If alg = 1, use maximp as an upper limit for the number of covers.   */
    /* If alg = 0, an upper limit is determined in bfa_lm().                */
    /*                                                                      */
    /* nf is the final number of covers in PMImplS and PMImplL.             */

    if (PMALG == 1) {
        printf1("Find minimal covers with Lawler's method.\n");

        nf = bfa_lm(ni,nj);      /* Lawler's method */           
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    else if (PMALG == 3) {
        printf1("Find minimal covers with Lawler's method II.\n");

        nf = bfa_lm1(ni,nj);      /* Lawler's method */           
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    else if (PMALG == 2) {
        printf1("Find minimal covers with Petrick's method.\n");

        ma = ni * maximp + 1;
        if (!(BMImplS = (char *)calloc(ma,sizeof(char)))) { 
            p_err(-2,1);
            goto BFAFin;
        }
        BMImplS_A = ma;           
        memrq(ma,sizeof(char));

        ma = maximp + 1;
        if (!(BMImplL = (short *)calloc(ma,sizeof(short)))) { 
            p_err(-2,1);
            goto BFAFin;
        }
        BMImplL_A = ma;           
        memrq(ma,sizeof(short));

        nf = bfa_pm(ni,nj,maximp);      /* Petrick's method */           
        if (nf < 0) {
            err = nf;
            goto BFAFin;
        }
    }
    if (alloc_ack(ni + 1))                                    
        goto BFAFin;
         
    for (i = 0; i < ni; ++i)  
        AcK[i] = bfa_len(BMImplP1[i],n,m);

    mlen = INTMAX;
    for (i = 0; i < nf; ++i) {
        k = 0;
        for (j = 0; j < ni; ++j)
            k += BMImplS[i * ni + j] * AcK[j];
        BMImplL[i] = k;
        mlen = imin(k,mlen);
    }
    nsel = 0;
    for (i = 0; i < nf; ++i) {
        if (BMImplL[i] == mlen) {
            BMImplL[i] = 1;
            nsel++;
        }
        else
            BMImplL[i] = 0;
    }
    printf1("Number of covers for non-essential implicants: %d\n",nsel);

BFAFin1:    

    /* -------------------------------------------------------------------- */
    /* Print final covers.                                                  */  

    printf1("\nFinal selection of prime implicants.\n");
    bfa_prn(nei,nf,ni,n,m,0);
    newline();
    bfa_prn(nei,nf,ni,n,m,1);


    /* -------------------------------------------------------------------- */
    /* Final check of selections.                                           */  

    fnd = 0;
    for (i = 0; i < nf; ++i) {
        if (BMImplL[i]) {
            for (j = 0; j < n; ++j) {
                k = bfa_fcheck((unsigned int)j,n,m,nei,ni,BMImplS + i * ni);
                if ((AcN[j] > 0 && k == 0) || (AcM[j] > 0 && k == 1)) {
                    printf1("Final check - error in row %d - selection %d\n",j,i+1);
                    fnd++;
                }
            }
        }
    }
    printf1("\nFinal check: %d errors.\n",fnd);
    err = 0;

BFAFin:
    if (err == -2) {
        printf1("Exceeded maximum number of implicants (%d).\n",maximp);
        err = -1;
    }
    else if (err == -3) {
        printf1("Exceeded maximum levels in bfa_lm().\n");
        err = -1;
    }
    bfa_free();
    if (BMImplC_A > 0) {
        free((char *)BMImplC);
        memrq(-BMImplC_A,sizeof(char));
        BMImplC_A = 0;    
    }
    if (BMImplS_A > 0) {
        free((char *)BMImplS);
        memrq(-BMImplS_A,sizeof(char));
        BMImplS_A = 0;    
    }
    if (BMImplL_A > 0) {
        free((char *)BMImplL);
        memrq(-BMImplL_A,sizeof(short));
        BMImplL_A = 0;    
    }
    if (BMImplE_A > 0) {
        free((char *)BMImplE);
        memrq(-BMImplE_A,sizeof(int));
        BMImplE_A = 0;    
    }
    if (BMImplE1_A > 0) {
        free((char *)BMImplE1);
        memrq(-BMImplE1_A,sizeof(int));
        BMImplE1_A = 0;    
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  bfa_a   allocated memory. Return 0 if successful, else -1.              */

int bfa_a(int mmax)
{
    if (!(BMImplA = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplA_A = mmax;           
    memrq(mmax,sizeof(int));

    if (!(BMImplA1 = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplA1_A = mmax;           
    memrq(mmax,sizeof(int));

    if (!(BMImplB = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplB_A = mmax;           
    memrq(mmax,sizeof(int));

    if (!(BMImplB1 = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplB1_A = mmax;           
    memrq(mmax,sizeof(int));
  
    if (!(BMImplP = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplP_A = mmax;           
    memrq(mmax,sizeof(int));

    if (!(BMImplP1 = (unsigned int *)calloc(mmax,sizeof(int)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplP1_A = mmax;           
    memrq(mmax,sizeof(int));

    if (!(BMImplX = (char *)calloc(mmax,sizeof(char)))) { 
        p_err(-2,1);
        return(-1);;
    }
    BMImplX_A = mmax;           
    memrq(mmax,sizeof(char));

    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_free()  free allocated memory.                                      */

int bfa_free(void)
{
    if (BMImplA_A > 0) {
        free((char *)BMImplA);
        memrq(-BMImplA_A,sizeof(int));
        BMImplA_A = 0;    
    }
    if (BMImplA1_A > 0) {
        free((char *)BMImplA1);
        memrq(-BMImplA1_A,sizeof(int));
        BMImplA1_A = 0;    
    }
    if (BMImplB_A > 0) {
        free((char *)BMImplB);
        memrq(-BMImplB_A,sizeof(int));
        BMImplB_A = 0;    
    }
    if (BMImplB1_A > 0) {
        free((char *)BMImplB1);
        memrq(-BMImplB1_A,sizeof(int));
        BMImplB1_A = 0;    
    }
    if (BMImplP_A > 0) {
        free((char *)BMImplP);
        memrq(-BMImplP_A,sizeof(int));
        BMImplP_A = 0;    
    }
    if (BMImplP1_A > 0) {
        free((char *)BMImplP1);
        memrq(-BMImplP1_A,sizeof(int));
        BMImplP1_A = 0;    
    }
    if (BMImplX_A > 0) {
        free((char *)BMImplX);
        memrq(-BMImplX_A,sizeof(char));
        BMImplX_A = 0;    
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_p(x) Print x as a bit pattern.                                      */

void bfa_p(unsigned int x,int n,int m)
{
    register int j;
    unsigned int ui;

    ui = n;
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & x)
            printf1("1");
        else
            printf1("0");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps(x)   Print x as a bit pattern. If opt = 1 print variable names.  */

void bfa_ps(unsigned int x,unsigned int y,int n,int m,int opt)
{
    register int j,k,l;
    unsigned int ui;

    ui = n;
    if (opt == 0) {
        for (j = 0; j < m; ++j) {
            ui = ui >> 1;
            if (ui & (~y))  
                printf1("-");
            else if (ui & x)  
                printf1("1");
            else
                printf1("0");
        }
    }
    else {
        l = k = 0;
        for (j = 0; j < m; ++j) {
            ui = ui >> 1;
            if (ui & (~y))  
                ;                    
            else {
                l++;
                if (k)
                    printf1(".");
                printf1("%s",VName[PMVIdx[j + 1]]);
                k = 1;
                if (ui & x)  
                    ;
                else
                    printf1("'");
            }
        }
        if (l == 0)
            printf("[empty]");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps_p(x)   Print x as a bit pattern to the protocol file.            */

void bfa_ps_p(unsigned int x,unsigned int y,int n,int m)
{
    register int j;
    unsigned int ui;

    if (PMProtFDef == 0)  
        return;

    ui = n;
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & (~y))
            fprintf(PMProtFd,"-");
        else if (ui & x)
            fprintf(PMProtFd,"1");
        else
            fprintf(PMProtFd,"0");
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_len(x)  Return number of nonzero bits in x.                         */

int bfa_len(unsigned int x,int n,int m)
{
    register int j,r;
    unsigned int ui;

    r = 0;
    ui = n;
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & x)
            r++; 
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  bfa_find()  Find all implicants. The prime implicants are stored in     */
/*              BMImplP and BMImplP1. Return the number of implicants.      */
/*              Return -2 if exceeded storage for max number of implicants. */

int bfa_find(int m,int n,int maximp)
{
    register int i,j,k;
    int d,nia,nib,nip,fnd;
    unsigned int ai,aj,ax,xx,am;

    /* begin with implicants of order 0 */

    k = 0;
    for (i = 0; i < n; ++i) {
        if (AcM[i] == 0) {
            if (k >= maximp)  
                return(-2);
            BMImplA[k] = i;
            k++;
        }
    }
    nia = k;
    for (i = 0; i < nia; ++i)  
        BMImplA1[i] = 0xffffff;

    nip = 0;    /* number of prime implicants */

    for (d = 1; d <= m; ++d) {  /* d is order of implicant */

        fprintf(stderr,"Level %d of %d (nia=%d, nip=%d)\n",d,m,nia,nip);
        fflush(stderr);

        if (nia == 0)
            break;

        nib = 0;
        for (i = 1; i < nia; ++i) {
            ai = BMImplA[i];           
            ax = BMImplA1[i];
            ai = ai & ax;

            for (j = 0; j < i; ++j) {
                if (ax != BMImplA1[j])
                    continue;

                aj = BMImplA[j] & ax;      
                xx = ai ^ aj;
                fnd = 0;
                for (k = 0; k < m; ++k) {
                    if (xx & 0x01) {
                        if (xx == 0x01)  
                            fnd = 1;
                        break;
                    }
                    xx = xx >> 1;
                }  
                if (fnd) {
                    BMImplX[i] = 1;     /* mark i and j as used */
                    BMImplX[j] = 1;

                    /* check if previously found */

                    am = ax & (~(ai ^ aj));

                    for (k = 0; k < nib; ++k) {
                        if (BMImplB1[k] == am && (BMImplB[k] & am) == (ai & am)) {
                            fnd = 0;
                            break;
                        }
                    }
                    if (fnd) {
                        if (nib >= maximp)  
                           return(-2);
                        BMImplB[nib] = ai;
                        BMImplB1[nib] = am;
                        nib++;
                    }
                }
            }
        }

        /* save prime implicants  */

        for (i = 0; i < nia; ++i) {
            if (BMImplX[i] == 0) {
                if (nip >= maximp)  
                    return(-2);
                BMImplP[nip] = BMImplA[i];
                BMImplP1[nip] = BMImplA1[i];
                nip++;
            }
        }
        for (i = 0; i < nib; ++i) {
            BMImplA[i] = BMImplB[i];
            BMImplA1[i] = BMImplB1[i];
            BMImplX[i] = 0;        
        }
        nia = nib;
        nib = 0;
    }
    if (nia > 0) {
        for (i = 0; i < nia; ++i) {
            if (nip >= maximp)  
                return(-2);
            BMImplP[nip] = BMImplA[i];
            BMImplP1[nip] = BMImplA1[i];
            nip++;
        }
    }
    return(nip);
} 

/* ------------------------------------------------------------------------ */
/*  bfa_elim()  Eliminate redundant implicants from BMImplP, and BMImplP1.  */
/*              Return nip = new number of prime implicants.                */

int bfa_elim(int m,int n,int nip)
{
    register int i,j,k;
    unsigned int ax,ai;

    j = 0;
    for (i = 0; i < nip; ++i) {
        ai = BMImplP[i];           
        ax = BMImplP1[i];
        ai = ai & ax;
        for (k = 0; k < n; ++k) {
            if (AcN[k] > 0 && (k & ax) == ai) {
                BMImplP[j] = BMImplP[i];
                BMImplP1[j] = BMImplP1[i];
                j++;
                break; 
            }
        }
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  bfa_chart() Make prime implicant chart in BMImplC (nip x ny1)           */
/*              Return  0 if successful.                                    */
/*                     -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */

int bfa_chart(int m,int n,int nip,int maximp)
{
    register int i,j,k;
    int ma,nia; 
    unsigned int ai,ax,am;

    k = 0;
    for (i = 0; i < n; ++i) {
        if (AcM[i] == 0 && AcN[i] > 0) {
            if (k >= maximp)  
                return(-2);
            BMImplA[k] = i;
            k++;
        }
    }
    nia = k;
    for (i = 0; i < nip; ++i) {
        ai = BMImplP[i]; 
        ax = BMImplP1[i];
        am = ai & ax;
        for (j = 0; j < nia; ++j) {
            if (am == (BMImplA[j] & ax))
                BMImplC[i * nia + j] = 1;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_ess()   Find essential prime implicants form prime implicant chart. */
/*              Mark them by a 1 in BMImplX (i = 0,...,nip-1).              */
/*              Return number of essential implicants.                      */

int bfa_ess(int nip,int ny1)
{
    register int i,j,k,l;
    int nei;

    nei = 0;
    for (i = 0; i < nip; ++i)  
        BMImplX[i] = 0;

    for (j = 0; j < ny1; ++j) {
        k = 0;
        l = -1;
        for (i = 0; i < nip; ++i) {
            if (BMImplC[i * ny1 + j]) {
                k++;
                if (l < 0)
                    l = i;
            }
        }
        if (k == 1 && l >= 0)  
            BMImplX[l] = 1;
    }
    for (i = 0; i < nip; ++i) {
        if (BMImplX[i])
            nei++;
    }
    return(nei);
}

/* ------------------------------------------------------------------------ */
/*  bfa_red()   Reduce the prime implicant chart by deleting the essential  */
/*              implicants. It is assumed that BMImplX still marks the      */
/*              essential prime implicants.                                 */
/*                                                                          */
/*              Return 0 if successful (and ni and nj of reduced table)     */
/*                    -1 if insufficient memory.                            */

int bfa_red(int nip,int ny1,int *ni,int *nj)
{
    register int i,j,ii,jj;

    if (alloc_acj(ny1 + 1))   
        return(-1);  

    for (j = 0; j < ny1; ++j) {
        for (i = 0; i < nip; ++i) {
            if (BMImplX[i] && BMImplC[i * ny1 + j]) {
                AcJ[j] = 1;
                break;
            }
        }
    }
    *ni = *nj = 0;
    for (i = 0; i < nip; ++i) {
        if (BMImplX[i] == 0)
            *ni += 1;
    }
    for (j = 0; j < ny1; ++j) {
        if (AcJ[j] == 0)
            *nj += 1;
    }
    if (*ni == 0 || *nj == 0)
        return(0);

    ii = 0;
    for (i = 0; i < nip; ++i) {
        if (BMImplX[i] == 0) {
            jj = 0;
            for (j = 0; j < ny1; ++j) {
                if (AcJ[j] == 0) {
                    BMImplC[ii * *nj + jj] = BMImplC[i * ny1 + j];
                    jj++;
                }
            }
            ii++;
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_pm()    Enumeration (Petrick's method). The algorithm follows       */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */ 
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */  
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */

int bfa_pm(int ni,int nj,int maximp)
{
    register int i,j,k,l;
    int is,isel,ii,jj,fnd;
    char *xptr;

    if (alloc_acj(nj + 1))   
        return(-1);  
    if (alloc_ack(nj + 1))   
        return(-1);  

    for (i = 0; i < ni; ++i) {
        k = i * nj;
        for (j = 0; j < nj; ++j) {
            if (BMImplC[k + j])
                AcJ[j] += 1;
        }
    }
    if (sortdpi(nj,AcJ,AcK))
        return(-1);

    is = isel = 0;
    for (jj = 0; jj < nj; ++jj) {
        j = AcK[jj];
        /****************
        fprintf(stderr,"column %d of %d\n",jj,nj);
        fflush(stderr);
        *******************/
        if (isel == 0) {
            for (i = 0; i < ni; ++i) {
                if (BMImplC[i * nj + j]) {
                    if (isel >= maximp)  
                        return(-2);
                    for (k = 0; k < ni; ++k) {
                        if (k == i)
                            BMImplS[isel * ni + k] = 1;
                        else
                            BMImplS[isel * ni + k] = 0;
                    }
                    isel++;
                }
            }
            is = isel;
        }
        else {
            is = isel;
            for (i = 0; i < ni; ++i) {
                if (BMImplC[i * nj + j]) {
                    for (l = 0; l < isel; ++l) {
                        if (is >= maximp)  
                            return(-2);

                        for (k = 0; k < ni; ++k) {
                            if (k == i)
                                BMImplS[is * ni + k] = 1;
                            else
                                BMImplS[is * ni + k] = BMImplS[l * ni + k];
                        }
                        is++;
                    }
                }
            }
        }

        /* find minimal elements */

        if (isel == is)
            continue;
             
        for (l = isel; l < is; ++l) {
            BMImplL[l] = 1;
            xptr = BMImplS + l * ni;
            for (k = isel; k < is; ++k) {
                if (k != l && bfa_minc(ni,BMImplS + k * ni,xptr)) {
                    BMImplL[l] = 0;
                    break;
                }
            }
        }
        ii = 0;
        for (l = isel; l < is; ++l) {
            if (BMImplL[l]) {
                fnd = 0;
                for (k = 0; k < ii; ++k) {  /* does vector already exist? */
                    if (bfa_icheck(ni,BMImplS + k * ni,BMImplS + l * ni)) {
                        fnd = 1;
                        break;
                    }
                }
                if (fnd == 0) {
                    for (k = 0; k < ni; ++k) 
                        BMImplS[ii * ni + k] = BMImplS[l * ni + k];
                    ii++;
                }
            }   
        }
        isel = ii;
    }
    fnd = INTMAX;
    for (l = 0; l < isel; ++l) {
        ii = 0;
        for (k = 0; k < ni; ++k) {
            if (BMImplS[l * ni + k])
                ii++;
        } 
        BMImplL[l] = ii;
        fnd = imin(fnd,ii);
    }
    j = 0;
    for (i = 0; i < isel; ++i) {
        if (BMImplL[i] > fnd)
            continue;
        if (i > j) {
            bfa_copy(ni,BMImplS + j * ni,BMImplS + i * ni);  
            ii++;
        }
        j++; 
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  bfa_minc(n,x,y)     Return 1 if x is less than y, otherwise 0.          */

int bfa_minc(int n,char *x, char *y)
{
    register int i,m;

    m = 0;
    for (i = 0; i < n; ++i) {
        if (x[i] > y[i])  
            return(0);
        else if (x[i] < y[i])
            m++;
    }
    if (m)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_mince(n,x,y)    Return 1 if x is less or equal than y, otherwise 0. */

int bfa_mince(int n,char *x, char *y)
{
    register int i;

    for (i = 0; i < n; ++i) {
        if (x[i] > y[i])  
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  bfa_icheck(n,x,y)   Return 1 if x and y are identical, otherwise 0.     */

int bfa_icheck(int n,char *x, char *y)
{
    register int i;

    for (i = 0; i < n; ++i) {
        if (x[i] != y[i])  
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  bfa_copy(n,x,y)     Copy y into x.                                      */

void bfa_copy(int n,char *x, char *y)
{
    register int i;
    for (i = 0; i < n; ++i)   
        x[i] = y[i];  
}

/* ------------------------------------------------------------------------ */
/*  bfa_fcheck()    Check final selections. nei is the number of essential  */
/*                  implicants in BMImplE.                                  */
/*                                                                          */  

int bfa_fcheck(unsigned int x,int n,int m,int nei,int ni,char *sel)
{
    register int i;

    for (i = 0; i < nei; ++i) {
        if ((x & BMImplE1[i]) == (BMImplE[i] & BMImplE1[i]))  
            return(1);
    }
    for (i = 0; i < ni; ++i) {
        if (sel[i]) {
            if ((x & BMImplP1[i]) == (BMImplP[i] & BMImplP1[i]))  
                return(1);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_prn()   Print final selections.                                     */

void bfa_prn(int nei,int nf,int ni,int n,int m,int opt)
{
    register int i,j,k,l;

    if (nf == 0) {
        for (j = 0; j < nei; ++j) {
            if (j)  
                printf1(" + ");
            bfa_ps(BMImplE[j],BMImplE1[j],n,m,opt);
        }                  
        newline();
        return;
    }
    k = 0;
    for (i = 0; i < nf; ++i) {
        if (BMImplL[i]) {
            l = 0;
            printf1("%3d  ",++k);
            for (j = 0; j < nei; ++j) {
                if (l)  
                    printf1(" + ");
                bfa_ps(BMImplE[j],BMImplE1[j],n,m,opt);
                l = 1;
            }                  
            for (j = 0; j < ni; ++j) {
                if (BMImplS[i * ni + j]) {
                    if (l)            
                        printf1(" + ");
                    bfa_ps(BMImplP[j],BMImplP1[j],n,m,opt);
                    l = 1;
                }
            }
            newline();
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  bfa_lm()    Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */ 
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */  
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

int bfa_lm(int ni,int nj)
{
    register int i,j,k,l;
    int ii,ik,jk,cmin,nej,ne,nlev,maxlev,ma,nf,jmin;

int ncnt;


    maxlev = nj + 10;
    nf = 0;

    if (alloc_acl(ni * maxlev + 1))    
        return(-1);        
    if (alloc_aci(maxlev + 1))    
        return(-1);        
    if (alloc_acns(ni * nj + 1))    
        return(-1);        
    if (alloc_acms(ni + 1))    
        return(-1);        
    if (alloc_acr(ni + 1))    
        return(-1);        
    if (alloc_acs(ni + 1))    
        return(-1);        
    if (alloc_acc(ni + 1))    
        return(-1);        
    if (alloc_acd(ni + 1))    
        return(-1);        
    if (alloc_ace(ni * nj + 1))    
        return(-1);        
    if (alloc_acf(ni * nj + 1))    
        return(-1);        

    k = ni * nj;                /* use AcF as local copy of BMImplC */
    for (i = 0; i < k; ++i)
        AcF[i] = BMImplC[i];
    nej = nj;
    nlev = 0;

    while (nej > 0) {
printf("nej=%d\n",nej);            /* select column from AcF */
/**
        for (i = 0; i < ni; ++i) {
            printf("AcF  ");
            for (j = 0; j < nej; ++j)
                printf("%d ",AcF[i * nj + j]);
            newline();
        }
**/



/*********

        if (alloc_acj(nej + 1))   
            return(-1);  
        if (alloc_ack(nej + 1))   
            return(-1);  

        for (i = 0; i < ni; ++i) {
            k = i * nj;
            for (j = 0; j < nej; ++j) {
                if (AcF[k + j])
                    AcJ[j] += 1;
            }
        }
        if (sortdpi(nej,AcJ,AcK))
            return(-1);
        j = AcK[0];
**********/

        jmin = INTMAX;
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * nj + j])
                    k++;
            }
            jmin = imin(jmin,k);
        }
        printf("JMIN=%d\n",jmin);

        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * nj + j])
                    k++;
            }
            if (k == jmin)
                break;
        }

printf("Using j=%d nej=%d ni=%d nj=%d    \n",j,nej,ni,nj  );


        ik = 0;  
        for (i = 0; i < ni; ++i) {
            ii = i * nj;
            if (AcF[ii + j] == 0)
                continue;

            if (nlev >= maxlev)
                return(-3);

            AcL[nlev * ni + ik] = i;
            AcI[nlev] += 1;

            jk = 0;
            for (k = 0; k < nej; ++k) {
                if (AcF[ii + k] == 0) {
                    AcNS[ik * nj + jk] = k;
                    jk++;
                }
            }
            AcMS[ik] = jk;
            ik++;
        }
        nlev++;
        if (nej <= 1)
            break;

        ne = 0;
        for (i = 0; i < ik; ++i)
            AcR[i] = 0;

ncnt = 0;

        while (1) {

            ncnt++;
                                      /* create or'ed row */
            for (i = 0; i < ni; ++i) 
                AcC[i] = 0;                           
            for (k = 0; k < ik; ++k) {
                l = AcNS[k * nj + AcR[k]];
                for (i = 0; i < ni; ++i) 
                    AcC[i] = AcC[i] | AcF[i * nj + l];
            }
          
            /* --------------------------------------------- */

            cmin = 1;       /* assume AcC is minimal */
            for (i = 0; i < ik; ++i)
                AcS[i] = 0;

            while (1) {
                for (i = 0; i < ni; ++i) 
                    AcD[i] = 0;                           
                for (k = 0; k < ik; ++k) {
                    l = AcNS[k * nj + AcS[k]];
                    for (i = 0; i < ni; ++i) 
                        AcD[i] = AcD[i] | AcF[i * nj + l];
                }
                if (bfa_minc(ni,AcD,AcC)) {
                    cmin = 0;
                    break;
                }   
                for (k = 0; k < ik; ++k) {
                    AcS[k] += 1;
                    if (AcS[k] < AcMS[k])
                        break;
                    AcS[k] = 0;
                }
                if (k >= ik)
                    break;
            }

            /* --------------------------------------------- */

            if (cmin) {        /* save minimal column in AcE */
                for (i = 0; i < ni; ++i)
                    AcE[i * nj + ne] = AcC[i];
                ne++;
            }
            for (k = 0; k < ik; ++k) {
                AcR[k] += 1;
                if (AcR[k] < AcMS[k])
                    break;
                AcR[k] = 0;
            }
            if (k >= ik)
                break;
        }         
printf("ncnt=%d\n",ncnt);


                                     /* copy AcE to AcF */
        for (i = 0; i < ni; ++i) {
            for (j = 0; j < ne; ++j)  
                AcF[i * nj + j] = AcE[i * nj + j];
        }
        nej = ne;
    }

    printf("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        printf("AcI: %d ",AcI[i]);
        for (j = 0; j < AcI[i]; ++j)
            printf("%d ",AcL[i * ni + j]);        
        newline();
    }











    ma = 1;
    for (i = 0; i < nlev; ++i)
        ma *= AcI[i];

    if (!(BMImplS = (char *)calloc(ma * ni + 1,sizeof(char)))) { 
        p_err(-2,1);
        return(-1);  
    }
    BMImplS_A = ma * ni + 1;  
    memrq(ma * ni + 1,sizeof(char));

    if (!(BMImplL = (short *)calloc(ma + 1,sizeof(short)))) { 
        p_err(-2,1);
        return(-1);   
    }
    BMImplL_A = ma + 1;       
    memrq(ma + 1,sizeof(short));

    if (alloc_acr(nlev + 1))    
        return(-1);        

    while (1) {
        ii = 0;
        for (j = 0; j < nj; ++j) {
            ii = 0;
            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                ii += BMImplC[k * nj + j];
            }
            if (ii == 0)
                break;
        }
        if (ii) {
            if (nf >= ma)  
                return(-3);

            for (i = 0; i < ni; ++i)
                BMImplS[nf * ni + i] = 0;

            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                BMImplS[nf * ni + k] = 1;
            }
            nf++;
        }
        for (i = 0; i < nlev; ++i) {
            AcR[i] += 1;
            if (AcR[i] < AcI[i])
                break;
            AcR[i] = 0;
        }
        if (i >= nlev)
            break;
    }
    return(nf);
}


/* ------------------------------------------------------------------------ */
/*  bfa_lm1()   Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */ 
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */  
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

int bfa_lm2(int ni,int nj)
{
    register int i,j,k,l;
    int ii,ik,ll,jk,cmin,nej,ne,nlev,maxlev,ma,nf,nm,njmax;
    char x;

    njmax = 2 * nj;
    maxlev = nj + 10;
    nf = 0;

    if (alloc_acl(ni * maxlev + 1))    
        return(-1);        
    if (alloc_aci(maxlev + 1))    
        return(-1);        
    if (alloc_acns(ni * njmax + 1))    
        return(-1);        
    if (alloc_acs(ni + 1))    
        return(-1);        
    if (alloc_acr(ni + 1))    
        return(-1);        
    if (alloc_acc(ni + 1))    
        return(-1);        
    if (alloc_ace(ni * njmax + 1))    
        return(-1);        
    if (alloc_acf(ni * njmax + 1))    
        return(-1);        

                                /* use AcF as local copy of BMImplC */
    for (i = 0; i < ni; ++i) {
        for (j = 0; j < nj; ++j) 
            AcF[i * njmax + j] = BMImplC[i * nj + j];
    }
    nej = nj;
    ne = nlev = 0;

    while (nej > 0) {
printf("nej=%d\n",nej);            /* select column from AcF */
/**
        for (i = 0; i < ni; ++i) {
            printf("AcF  ");
            for (j = 0; j < nej; ++j)
                printf("%d ",AcF[i * njmax + j]);
            newline();
        }
**/
        ii = INTMAX;
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * njmax + j])
                    k++;
            }
            ii = imin(ii,k);
        }
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * njmax + j])
                    k++;
            }
            if (k == ii)
                break;
        }

printf("Using j=%d nej=%d ni=%d nj=%d ii=%d    \n",j,nej,ni,nj,ii);


        ik = 0;  
        for (i = 0; i < ni; ++i) {
            ii = i * njmax;
            if (AcF[ii + j] == 0)
                continue;

            if (nlev >= maxlev)     {
printf("nlev=%d maxlev=%d\n",nlev,maxlev);

                return(-3);
}
            AcL[nlev * ni + ik] = i;
            AcI[nlev] += 1;

            jk = 0;
            for (k = 0; k < nej; ++k) {
                if (AcF[ii + k] == 0) {
                    AcNS[ik * njmax + jk] = k;
                    jk++;
                }
            }
            AcS[ik] = jk;
            ik++;
        }
        printf("AcNS\n");
        for (i = 0; i < ik; ++i) {
            printf("AcS=%3d : ",AcS[i]);
            for (k = 0; k < AcS[i]; ++k)
                printf("%3d ",AcNS[i * njmax + k]);
            newline();
        }

        if (ik == 0)
            break;


        nlev++;
        if (nej <= 1)
            break;

        /* Create array with column numbers for or'ed columns */
        /* AcMS[nm x ik] for column numbers */
        /* AcD[nm] to index non-minimal elements */

        nm = 1;
        for (k = 0; k < ik; ++k) {
            if (AcS[k] > 0)
                nm *= AcS[k];
        }

printf("nm==========%d ik=%d AcS[0]=%d       \n",nm,ik,AcS[0]    );

        if (alloc_acms(nm * ik + 1))    
            return(-1);        
        if (alloc_acd(nm + 1))    
            return(-1);        

        for (i = 0; i < ik; ++i)
            AcR[i] = 0;

                   
        j = 0;
        while (1) {
            for (k = 0; k < ik; ++k) { 
                if (AcS[k])
                    AcMS[j * ik + k] = AcNS[k * njmax + AcR[k]];
                else
                    AcMS[j * ik + k] = -1;                        
            }
            j++;
            for (k = 0; k < ik; ++k) {
                AcR[k] += 1;
                if (AcR[k] < AcS[k])
                    break;
                AcR[k] = 0;
            }
            if (k >= ik)     
                break;
        }         
printf("hier ===%d nm=%d     \n",j,nm  );
/**
if (nm > 200000)
    return(-1);
**/
/*******  
        for (i = 0; i < nm; ++i) {
            for (k = 0; k < ik; ++k)
                printf("%4d ",AcMS[i * ik + k]);
            newline();
        }
**/      

        /* Find minimal elements. AcD[] = 1 if not minimal */

        /* first find ... */

        ne = 0;
        for (j = 0; j < nm; ++j) {
            k = 1;
            l = j * ik;
            ii = AcMS[l];
            for (i = 1; i < ik; ++i) {
                if (AcMS[++l] != ii) {
                    k = 0;
                    break;
                }
            }
            if (k) {        /* found minimal vector */

                AcD[j] = 2;

                for (i = 0; i < ni; ++i)  
                    AcC[i] = 0;

                for (k = 0; k < ik; ++k) {
                    if ((l = AcMS[j * ik + k]) >= 0) {
                        for (i = 0; i < ni; ++i)
                            AcC[i] |= AcF[i * njmax + l];
                    }
                }

                ii = 0;         
                for (k = 0; k < ne; ++k) {
                    l = 0;
                    for (i = 0; i < ni; ++i) {
                        if (AcE[k * ni + i] <= AcC[i] == 0)  
                            l++;
                    }
                    if (l == ni) {
                        ii = 1;
                        break;
                    }
                }
                if (ii)
                    printf("already found\n");

                if (ii == 0) {

                    if (ne >= njmax) {
                        printf1("Fatal error in bfa_lm() [nj=%d ne=%d].\n",nj,ne);
                        return(-1);
                    }
                    for (i = 0; i < ni; ++i)
                        AcE[ne * ni + i] = AcC[i];
                    ne++;
printf("ne=%d\n",ne);

                }

            }
        }
        nf = 0;

        for (j = 0; j < nm; ++j) {
            if (AcD[j])
                continue;


            for (i = 0; i < ni; ++i)  
                AcC[i] = 0;

            for (k = 0; k < ik; ++k) {
                if ((l = AcMS[j * ik + k]) >= 0) {
                    for (i = 0; i < ni; ++i)
                        AcC[i] |= AcF[i * njmax + l];
                }
            }
            for (i = 0; i < ne; ++i) {
                if (bfa_mince(ni,AcE + i * ni,AcC)) {         
                    AcD[j] = 1;
                    break;
                }
            }
            if (AcD[j] == 0) {
                nf++;
            }
        }
        printf("not found=%d\n",nf);

        if (nf > 0) {
            k = 0;
            for (j = 0; j < nm; ++j) {
                if (AcD[j] == 0) {
/*    printf("%4d: ",j);                            */
                    for (l = 0; l < ik; ++l) {
                       
/*  printf("%d ",AcMS[j * ik + l]);     */

                        AcMS[k * ik + l] = AcMS[j * ik + l];


                    }
/*  newline();  */


                    AcD[k] = 0;
                    k++;
                }
            }
            printf("kkkk=%d ne=%d      \n",k,ne);

            for (j = 0; j < nf; ++j) {
printf("> ne=%d\n",ne);
                for (i = 0; i < ni; ++i)  
                    AcC[i] = 0;

                for (k = 0; k < ik; ++k) {
                    if ((l = AcMS[j * ik + k]) >= 0) {
                        for (i = 0; i < ni; ++i)
                            AcC[i] |= AcF[i * njmax + l];
                    }
                }
                for (i = 0; i < ne; ++i) {
                    if (bfa_mince(ni,AcE + i * ni,AcC)) {         
                        AcD[j] = 1;
                        break;
                    }
                }
                if (AcD[j])
                    continue;

                for (ll = 0; ll < nf; ++ll) {
                    if (AcD[ll])   
                        continue;

                    ii = 0;                          
                    for (i = 0; i < ni; ++i) {
                        x = 0;
                        for (k = 0; k < ik; ++k) { 
                            if ((l = AcMS[ll * ik + k]) >= 0)
                                x |= AcF[i * njmax + l];
                        }
                        if (AcC[i]) {
                            if (x == 0)
                                ii++;
                        }
                        else if (x) {
                            ii = 0;
                            break;
                        }
                    }
                    if (ii) {  /* l is less than j */
                        AcD[j] = 1;
                        break;
                    }
                }
                if (AcD[j] == 0) {
                    if (ne >= njmax) {
                        printf1("Fatal error in bfa_lm() [nj=%d ne=%d].\n",nj,ne);
                        return(-1);
                    }
                    for (i = 0; i < ni; ++i)
                        AcE[ne * ni + i] = AcC[i];
                    ne++;
                }
            }
            printf("NEUES ne=%d\n",ne);
        }

        for (i = 0; i < ni; ++i) {
            for (j = 0; j < ne; ++j)  
                AcF[i * njmax + j] = AcE[j * ni + i];
        }
        nej = ne;


    }

    printf("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        printf("AcI: %d ",AcI[i]);
        for (j = 0; j < AcI[i]; ++j)
            printf("%d ",AcL[i * ni + j]);        
        newline();
    }



    ma = 1;
    for (i = 0; i < nlev; ++i)
        ma *= AcI[i];

    if (!(BMImplS = (char *)calloc(ma * ni + 1,sizeof(char)))) { 
        p_err(-2,1);
        return(-1);  
    }
    BMImplS_A = ma * ni + 1;  
    memrq(ma * ni + 1,sizeof(char));

    if (!(BMImplL = (short *)calloc(ma + 1,sizeof(short)))) { 
        p_err(-2,1);
        return(-1);   
    }
    BMImplL_A = ma + 1;       
    memrq(ma + 1,sizeof(short));

    if (alloc_acr(nlev + 1))    
        return(-1);        

    nf = 0;
    while (1) {
        ii = 0;
        for (j = 0; j < nj; ++j) {
            ii = 0;
            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                ii += BMImplC[k * nj + j];
            }
            if (ii == 0)
                break;
        }
        if (ii) {
            if (nf >= ma)  
                return(-3);
 
            for (i = 0; i < ni; ++i)
                BMImplS[nf * ni + i] = 0;

            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                BMImplS[nf * ni + k] = 1;
            }
            nf++;
        }
        for (i = 0; i < nlev; ++i) {
            AcR[i] += 1;
            if (AcR[i] < AcI[i])
                break;
            AcR[i] = 0;
        }
        if (i >= nlev)
            break;
    }
    return(nf);
}


/* ------------------------------------------------------------------------ */
/*  bfa_minc1(n,x,y)     Return 1 if x is less than y, otherwise 0.          */

int bfa_minc1(int n,char *x, char *y)
{
    register int i,m;

    printf("minc1 n=%d\n",n);
    printf("x ");
    for (i = 0; i < n; ++i)
        printf("%d ",x[i]);
    printf("\ny ");
    for (i = 0; i < n; ++i)
        printf("%d ",y[i]);
    newline();

    m = 0;
    for (i = 0; i < n; ++i) {
printf("%d (%d,%d)\n",i,x[i],y[i]);


        if (x[i] > y[i]) {
printf(" hier return 0\n");
            return(0);
        }
        else if (x[i] < y[i])
            m++;
    }
printf(" hier m=%d\n",m );
    if (m)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  bfa_lm1()   Find minimal covers. The algorithm is due to                */
/*              E.L. Lawler, Covering Problems: Duality Relations and a     */
/*              New Method of Solution, SIAM Journal on Applied             */ 
/*              Mathematics 14 (1966), 1115-32.                             */
/*                                                                          */
/*              Prime implicant chart is in BMImplC with ni rows and nj     */
/*              columns. The functions uses BMImplS with ni rows and        */  
/*              maximp columns for storage of selections.                   */
/*                                                                          */
/*              Return -1 if insufficient memory.                           */
/*                     -2 if exceeded storage for max number of implicants. */
/*                     -3 if maxlev (or ma) insufficient.                   */

int bfa_lm1(int ni,int nj)
{
    register int i,j,k,l;
    int ii,jj,ik,ll,jk,nej,ne,nlev,maxlev,ma,nf,nm,njmax,nec,maxnf;
    char x;

    njmax = 2 * nj;
    maxlev = nj + 10;
    maxnf = 200000;

    nf = 0;

    if (alloc_acl(ni * maxlev + 1))    
        return(-1);        
    if (alloc_aci(maxlev + 1))    
        return(-1);        
    if (alloc_ack(njmax + 1))    
        return(-1);        
    if (alloc_acns(ni * njmax + 1))    
        return(-1);        
    if (alloc_acs(ni + 1))    
        return(-1);        
    if (alloc_acr(ni + 1))    
        return(-1);        
    if (alloc_acc(ni + 1))    
        return(-1);        
    if (alloc_ace(ni * njmax + 1))    
        return(-1);        
    if (alloc_acf(ni * njmax + 1))    
        return(-1);        

                                /* use AcF as local copy of BMImplC */
    for (i = 0; i < ni; ++i) {
        for (j = 0; j < nj; ++j) 
            AcF[i * njmax + j] = BMImplC[i * nj + j];
    }
    nej = nj;
    ne = nlev = 0;

    while (nej > 0) {
printf("nej=%d ni=%d \n",nej,ni  );            /* select column from AcF */
          
        for (i = 0; i < ni; ++i) {
            printf("AcF  ");
            for (j = 0; j < nej; ++j)
                printf("%d ",AcF[i * njmax + j]);
            newline();
        }
       
        ii = INTMAX;
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * njmax + j])
                    k++;
            }
            ii = imin(ii,k);
        }
        for (j = 0; j < nej; ++j) {
            k = 0;
            for (i = 0; i < ni; ++i) {
                if (AcF[i * njmax + j])
                    k++;
            }
            if (k == ii)
                break;
        }

printf("Using j=%d nej=%d ni=%d nj=%d ii=%d    \n",j,nej,ni,nj,ii);


        ik = 0;  
        for (i = 0; i < ni; ++i) {
            ii = i * njmax;
            if (AcF[ii + j] == 0)
                continue;

            if (nlev >= maxlev)     {
printf("nlev=%d maxlev=%d\n",nlev,maxlev);

                return(-3);
}
            AcL[nlev * ni + ik] = i;
            AcI[nlev] += 1;

            jk = 0;
            for (k = 0; k < nej; ++k) {
                if (AcF[ii + k] == 0) {
                    AcNS[ik * njmax + jk] = k;
                    jk++;
                }
            }
            AcS[ik] = jk;
            ik++;
        }
        printf("AcNS\n");
        for (i = 0; i < ik; ++i) {
            printf("AcS=%3d : ",AcS[i]);
            for (k = 0; k < AcS[i]; ++k)
                printf("%3d ",AcNS[i * njmax + k]);
            newline();
        }
        j = 0;
        for (i = 0; i < ik; ++i) {
            ii = AcS[i];
            if (ii > 0) {
                for (k = 0; k < ii; ++k)
                    AcNS[j * njmax + k] = AcNS[i * njmax + k];
                j++;
            }
        }
/****************

        if (j < ik) {
            ik = j;

            printf("Reduced AcNS\n");
            for (i = 0; i < ik; ++i) {
                printf("AcS=%3d : ",AcS[i]);
                for (k = 0; k < AcS[i]; ++k)
                    printf("%3d ",AcNS[i * njmax + k]);
                newline();
            }
        }
        ik = j;
***********/


        if (ik == 0)
            break;


        nlev++;
        if (nej <= 1)
            break;

        /* Create array with column numbers for or'ed columns */
        /* AcMS[nm x ik] for column numbers */
        /* AcD[nm] to index non-minimal elements */

        /* first find the number of equal column numbers: nec and save */
        /* in AcK */

        nec = 0;
        nm = 0;
        for (i = 0; i < ik; ++i) {
            AcR[i] = 0;
            if (AcS[i])              
                nm = 1;
        }
        if (nm != 0) {
                     
            while (1) {
                l = -1;                           
                ii = 0;
                for (k = 0; k < ik; ++k) { 
                    if (AcS[k] > 0) {
                        if (l < 0) {
                            l = AcNS[k * njmax + AcR[k]];
                            ii = 1;
                        }
                        else if (l != AcNS[k * njmax + AcR[k]]) {
                            ii = 0;
                            break;
                        }
                    }
                }
                if (ii)  
                    AcK[nec++] = l;

                for (k = 0; k < ik; ++k) {
                    AcR[k] += 1;
                    if (AcR[k] < AcS[k])
                        break;
                    AcR[k] = 0;
                }
                if (k >= ik)     
                    break;
            }         

            printf("nec=%d\nAcK: ",nec );
            for (k = 0; k < nec; ++k)
                printf("%d ",AcK[k]);
            newline();

            /* reorganize AcS and AcNS */

            for (i = 0; i < ik; ++i) {
                jj = AcS[i];
                ll = 0;
                for (k = 0; k < jj; ++k) {
                    ii = AcNS[i * njmax + k];
                    nf = 1;
                    for (l = 0; l < nec; ++l) {
                        if (ii == AcK[l]) {
                            nf = 0;
                            break;
                        }
                    }
                    if (nf) { /* not found */
                        AcNS[i * njmax + ll] = ii;
                        ll++;
                    }
                }
                AcS[i] = ll;
            }

            printf("NEW AcNS\n");
            for (i = 0; i < ik; ++i) {
                printf("AcS=%3d : ",AcS[i]);
                for (k = 0; k < AcS[i]; ++k)
                    printf("%3d ",AcNS[i * njmax + k]);
                newline();
            }
            nm = 0;
            for (i = 0; i < ik; ++i) {
                if (AcS[i])
                    nm = 1;                      
            }
        }

printf("nm==========%d ik=%d \n",nm,ik);

        /* Save minimal vectors from equal column numbers */

        ne = 0;
        for (j = 0; j < nec; ++j) {
            k = AcK[j];
            if (ne >= njmax) {
                printf1("Fatal error in bfa_lm() [njmax=%d ne=%d].\n",njmax,ne);
                return(-1);
            }
            for (i = 0; i < ni; ++i)
                AcE[ne * ni + i] = AcF[i * njmax + k];
            ne++;
        }   


printf("actual ne=%d\n",ne);

        /* If nm > 0 find remaining minimal vectors */

        if (nm != 0) {

            nm = 1;
            for (i = 0; i < ik; ++i) {
                if (AcS[i])
                    nm *= AcS[i];                
            }

printf("Actual nm=%d\n",nm);


            if (alloc_acms(nm * ik + 1))    
                return(-1);        
            if (alloc_acd(nm + 1))    
                return(-1);        

            for (i = 0; i < ik; ++i)
                AcR[i] = 0;

            j = 0;
            while (1) {
                for (k = 0; k < ik; ++k) { 
                    if (AcS[k])
                        AcMS[j * ik + k] = AcNS[k * njmax + AcR[k]];
                    else
                        AcMS[j * ik + k] = -1;                        
                }
                j++;
                for (k = 0; k < ik; ++k) {
                    AcR[k] += 1;
                    if (AcR[k] < AcS[k])
                        break;
                    AcR[k] = 0;
                }
                if (k >= ik)     
                    break;
            }         
        printf("AcMS\n");
            for (i = 0; i < nm; ++i) {
                for (k = 0; k < ik; ++k)
                    printf("%4d ",AcMS[i * ik + k]);
                newline();
            }
           

            /* Find minimal elements. AcD[] = 1 if not minimal */

            for (j = 0; j < nm; ++j) {
                for (i = 0; i < ni; ++i)  
                    AcC[i] = 0;

                for (k = 0; k < ik; ++k) {
                    if ((l = AcMS[j * ik + k]) >= 0) {
                        for (i = 0; i < ni; ++i)
                            AcC[i] |= AcF[i * njmax + l];
                    }
                }
                for (i = 0; i < ne; ++i) {
                    if (bfa_mince(ni,AcE + i * ni,AcC)) {         
                        AcD[j] = 1;
                        break;
                    }
                }
                if (AcD[j])  
                    continue;

                for (ll = j + 1; ll < nm; ++ll) {
                    ii = 0;                          
                    for (i = 0; i < ni; ++i) {
                        x = 0;
                        for (k = 0; k < ik; ++k) { 
                            if ((l = AcMS[ll * ik + k]) >= 0)
                                x |= AcF[i * njmax + l];
                        }
                        if (AcC[i]) {
                            if (x == 0)
                                ii++;
                        }
                        else if (x) {
                            ii = 0;
                            break;
                        }
                    }
                    if (ii) {  /* l is less than j */
                        AcD[j] = 1;
                        break;
                    }
                }
                if (AcD[j] == 0) {      /* found new minimal vector */
                    if (ne >= njmax) {
                        printf1("Fatal error in bfa_lm() [njmax=%d ne=%d].\n",njmax,ne);
                        return(-1);
                    }
                    for (i = 0; i < ni; ++i)
                        AcE[ne * ni + i] = AcC[i];
                    ne++;
                }
            }
            printf("NEUES ne=%d\n",ne);
        }

        for (i = 0; i < ni; ++i) {
            for (j = 0; j < ne; ++j)  
                AcF[i * njmax + j] = AcE[j * ni + i];
        }
        nej = ne;


    }

    printf("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        printf("AcI: %d ",AcI[i]);
        for (j = 0; j < AcI[i]; ++j)
            printf("%d ",AcL[i * ni + j]);        
        newline();
    }

    /* Simplify levels */

    for (i = 1; i < nlev; ++i) {
        ii = AcI[i];
        for (j = 0; j < i; ++j) {
            jj = AcI[j];
            ll = 0; 
            for (k = 0; k < jj; ++k) {
                ma = AcL[j * ni + k];
                for (l = 0; l < ii; ++l) {
                    if (ma == AcL[i * ni + l]) {
                        ll++;
                        break;
                    }
                }
            }
            if (ll == jj) {
                for (k = 0; k < jj; ++k) {
                    ma = AcL[j * ni + k];
                    for (l = 0; l < ii; ++l) {
                        if (ma == AcL[i * ni + l]) {
                            AcL[i * ni + l] = -1;
                            break;
                        }
                    }
                }
            }
        }
        j = 0;
        for (l = 0; l < ii; ++l) {
            if ((k = AcL[i * ni + l]) >= 0)
                AcL[i * ni + j++] = k;
        }
        AcI[i] = j;
    }
    printf("LEVEL nlev=%d\n",nlev  );
    for (i = 0; i < nlev; ++i) {
        printf("AcI: %d ",AcI[i]);
        for (j = 0; j < AcI[i]; ++j)
            printf("%d ",AcL[i * ni + j]);        
        newline();
    }

          

    if (!(BMImplS = (char *)calloc(maxnf * ni + 1,sizeof(char)))) { 
        p_err(-2,1);
        return(-1);  
    }
    BMImplS_A = maxnf * ni + 1;  
    memrq(maxnf * ni + 1,sizeof(char));

    if (!(BMImplL = (short *)calloc(maxnf + 1,sizeof(short)))) { 
        p_err(-2,1);
        return(-1);   
    }
    BMImplL_A = maxnf + 1;       
    memrq(maxnf + 1,sizeof(short));

    if (alloc_acr(nlev + 1))    
        return(-1);        

    nf = 0;
    while (1) {
        ii = 0;
        for (j = 0; j < nj; ++j) {
            ii = 0;
            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                ii += BMImplC[k * nj + j];
            }
            if (ii == 0)
                break;
        }
        if (ii) {
            if (nf >= maxnf)  
                return(-3);
 
            for (i = 0; i < ni; ++i)
                BMImplS[nf * ni + i] = 0;

            for (i = 0; i < nlev; ++i) {
                k = AcL[i * ni + AcR[i]];
                BMImplS[nf * ni + k] = 1;
            }
            ii = 0;
            for (k = 0; k < nf; ++k) {
                ii = 1;
                for (i = 0; i < ni; ++i) {
                    if (BMImplS[nf * ni + i] != BMImplS[k * ni + i]) {
                        ii = 0;
                        break;
                    }
                }
                if (ii)     /* already found */
                    break;
            }
            if (ii == 0)
                nf++;
        }
        for (i = 0; i < nlev; ++i) {
            AcR[i] += 1;
            if (AcR[i] < AcI[i])
                break;
            AcR[i] = 0;
        }
        if (i >= nlev)
            break;
    }

    printf("final nf=%d\n",nf);


    return(nf);
}

/* ------------------------------------------------------------------------ */
/*  bfc     Context-dependencies in Boolean functions                       */
/*                                                                          */
/*  bfc(                                                                    */
/*      opt=...,            treatment of undefined arguments, def. 1        */
/*                          1 treated as don't cares                        */
/*                          2 treated as logical 1                          */
/*                          3 treated as logical 0                          */
/*      prot=...,           protocol file with additional information       */
/*      max=...,            max number of implicants, def. 10000            */
/*                                                                          */
/*  ) = Y,X1,...,Xm;        variables for Boolean function                  */
/*                          (Y is the dependent variable)                   */
/*                                                                          */
/*  Boolean function: Y = f(X1,...,Xm). All variables are interepreted      */
/*  as binary variables. A value not equal to zero will be interpreted      */
/*  as a logical 1.                                                         */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int bfc(void)
{
    register int i,j,k;
    int err,maximp,m,n,m1,jv,ny0,ny1,nyu,ii,nip;      

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("Context-dependencies in Boolean functions. Current memory: %d bytes.\n\n",MemReq);

    if (parm(CmdBuf + 3,4,1))     /* get parameters */
        goto BFCFin;

    if (PMMax > 0)             
        maximp = PMMax;
    else
        maximp = IMPMax;

    m = PMNV - 1;
    if (m < 2) {
        printf1("Error: need at least three variables.\n");
        goto BFCFin;
    }
    else if (m > BMMax) {
        printf1("Error: maximum number of arguments is %d.\n",BMMax);
        goto BFCFin;
    }
    m1 = m - 1;
    n = (int)pow(2.0,(double)m1);

    printf1("Number of arguments: %d\n",m);
    printf1("Rows of reduced truth table: %d\n",n);
    if (PMProtFDef)
        printf1("Protocol file: %s\n",PMProtFName);
    newline();

    if (bfa_a(maximp))              /* allocate memory */
        goto BFCFin;

    /* ---------------------------------------------------------------- */
    /* Do for all independent variables: ...                            */

    for (jv = 1; jv <= m; ++jv) {
        printf1("Variable: %s\n",VName[PMVIdx[jv]]);

        /* AcN: X = 0 and Y = 0 */
        /* AcM: X = 0 and Y = 1 */
        /* AcI: X = 1 and Y = 0 */
        /* AcJ: X = 1 and Y = 1 */
        /* AcR = 0 if not known
                 1 if positive effect
                 2 if negative effect
                 3 if no effect       
                 4 if inconsistent */

        if (alloc_acn(n + 1))                                  
            goto BFCFin;
        if (alloc_acm(n + 1))                                  
            goto BFCFin;
        if (alloc_aci(n + 1))                                  
            goto BFCFin;
        if (alloc_acj(n + 1))                                  
            goto BFCFin;
        if (alloc_acr(n + 1))                                  
            goto BFCFin;
               
        /* create reduced truth table */

        for (i = 0; i < NOC; ++i) {
            ii = 0;
            k  = 1;
            for (j = PMNV - 1; j >= 1; --j) {
                if (j == jv)
                    continue;

                if ((int)get_data(PMVIdx[j],i))
                    ii += k;
                k *= 2;
            }
            if ((int)get_data(PMVIdx[0],i)) {
                if ((int)get_data(PMVIdx[jv],i))   
                    AcJ[ii] += 1;
                else
                    AcM[ii] += 1;
            }
            else {
                if ((int)get_data(PMVIdx[jv],i))   
                    AcI[ii] += 1;
                else
                    AcN[ii] += 1;
            }
        }
        ii = 0;
           
        for (i = 0; i < n; ++i) {
            if (AcN[i] > 0 && AcM[i] > 0) {
                AcR[i] = 4;
                ii = 1;
            }
            else if (AcI[i] > 0 && AcJ[i] > 0) {
                AcR[i] = 4;
                ii = 1;
            }
            else if (AcN[i] > 0 && AcJ[i] > 0)
                AcR[i] = 1;
            else if (AcM[i] > 0 && AcI[i] > 0)
                AcR[i] = 2;
            else if (AcN[i] > 0 && AcI[i] > 0)
                AcR[i] = 3;
            else if (AcM[i] > 0 && AcJ[i] > 0)
                AcR[i] = 3;
        }
        if (PMProtFDef) {
            fprintf(PMProtFd,"Variable: %s\n",VName[PMVIdx[jv]]);
            for (i = 0; i < n; ++i) {
                fprintf(PMProtFd,"%4d ",i);
                fprintf(PMProtFd,"%4d ",AcN[i]);
                fprintf(PMProtFd,"%4d ",AcM[i]);
                fprintf(PMProtFd,"%4d ",AcI[i]);
                fprintf(PMProtFd,"%4d  ",AcJ[i]);
                if (AcR[i] == 0)
                    fprintf(PMProtFd,"(not known)\n");
                else if (AcR[i] == 1)
                    fprintf(PMProtFd,"(positive effect)\n");
                else if (AcR[i] == 2)
                    fprintf(PMProtFd,"(negative effect)\n");
                else if (AcR[i] == 3)
                    fprintf(PMProtFd,"(no effect)\n");
                else                   
                    fprintf(PMProtFd,"(inconsistent)\n");
            }
        }
        if (ii) {
            printf1("-- data are inconsistent.\n\n");
            continue;
        }

        for (ii = 0; ii <= 3; ++ii) {
            if (ii == 0)
                printf1("Not known: ");
            else if (ii == 1)
                printf1("Positive effect: ");
            else if (ii == 2)
                printf1("Negative effect: ");
            else            
                printf1("No effect: ");

            for (i = 0; i < n; ++i) {
                if (AcR[i] == ii)
                    AcM[i] = 0;
                else
                    AcM[i] = 1;
            }
            nip = bfa_find(m1,n,maximp);                                  
            if (nip < 0) {
                err = nip;
                goto BFCFin;
            }
            else if (nip == 0)
                printf1("none\n");
            else {
                for (i = 0; i < nip; ++i) {
                    if (i)
                         printf1(", ");
                    bfa_ps_n(BMImplP[i],BMImplP1[i],n,m1,jv-1);
                }
                newline();
            }
        }
        newline();
    }

    err = 0;

BFCFin:
    if (err == -2) {
        printf1("Exceeded maximum number of implicants (%d).\n",maximp);
        err = -1;
    }
    else if (err == -3) {
        printf1("Found inconsistent data (check with bfa).\n");
        err = -1;
    }
    bfa_free();
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  bfa_ps_n(x) Print variable names according to x, skip jv.               */

void bfa_ps_n(unsigned int x,unsigned int y,int n,int m,int jv)
{
    register int j,k,l,j1;
    unsigned int ui;

    ui = n;
    l = k = 0;
    for (j = 0; j < m; ++j) {
        ui = ui >> 1;
        if (ui & (~y))  
            ;                    
        else {
            l++;
            if (k)
                printf1(".");

            if (j < jv)
                j1 = j;
            else
                j1 = j + 1;

            printf1("%s",VName[PMVIdx[j1 + 1]]);
            k = 1;
            if (ui & x)  
                ;
            else
                printf1("'");
        }
    }
    if (l == 0)
        printf("always");
}


