/****************************************************************************/
/*  t_dens                                                                  */
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
#include "t_var.h"
#include "t_plot.h"
#include "t_alloc.h"
#include "t_cdf.h"

/*  functions in t_dens.c */

int kdens(int idx);
int kdensf(int n,float *x,int m,float a,float d,float sig,float *crit);
int kdcomp(const void *arg1,const void *arg2);

/* ------------------------------------------------------------------------ */
/*  kdens(idx)  Density estimation.                                         */
/*  ##          kdens(v=,sig=,x=,sc=,pl=) [=fname]                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int kdens(int idx)
{
    register int i;
    int err,n,m,mn,x_a;
    char *p;
    float xmin,xmax,a,b,d,crit;
    float *x;

    printfe("kdens not implemented.\n");
    gerr_exit(203);

    crit = 0.00001;

    x_a = 0;
    err = -1;
    if (check_cmd(2))
        return(-1);
/*  p = CmdDef[idx];        */
    p = CmdBuf;

    /** prnchar('-',LLEN,1); **/ 
    printf1("Density estimation: %s\n",p);

    if (parm(p + 5,1,0))    /* get parameters */
        goto KDENSFin;
        
    if (PMNV < 1) {
        printf1("Error: no variable.\n");
        goto KDENSFin;
    }
    if (PMSIG <= 0.0)   /* default smoothing factor */
        PMSIG = 0.5;

    if (PMSC <= 0.0)    /* default scaling factor */
        PMSC = 1.0;

    if (PMRXFlg) {      /* range */
        m = 1 + (int)ceil((double)((PMRXB - PMRXA) / PMRXD));
        a = PMRXA;
        d = PMRXD;
        mn = m;
    }
    else
        mn = m = 101;

    if (mn < NOC)
        mn = NOC;

    if (!(x = (float *)calloc(mn,sizeof(float)))) { 
        p_err(-2,1);
        goto KDENSFin;
    }
    memrq(mn,sizeof(float));
    x_a = 1;

    xmin = xmax = (float)get_data((int)PMVIdx[0],0);

    for (i = 1; i < NOC; ++i) {
        x[i] = (float)get_data((int)PMVIdx[0],i);
        if (xmin > x[i])
            xmin = x[i];
        if (xmax < x[i])
            xmax = x[i];
    }
    printf1("Variable: %s  Min: %g  Max: %g\n",VName[PMVIdx[0]],xmin,xmax);

    if (PMRXFlg == 0) {
        a = xmin;
        d = (xmax - xmin) / 100.0;
    }
    n = kdensf(NOC,x,m,a,d,PMSIG,&crit);
    if (n) {
        p_err(-2,1);
        goto KDENSFin;
    }
    printf1("Max value of density: %g\n",(double)crit);

    if (PMFDef) {          /* if output file defined */

        if (PMFmtF == 0) {
            PMFmt1 = 10;
            PMFmt2 = 4;
            makefmt(&PMFmt1,&PMFmt2,PMFmtS,0,SEPC,0);
        }
        fprintf(PMFd,"# Variable: %s  Sigma: %g ",VName[PMVIdx[0]],PMSIG);
        fprintf(PMFd,"Max value of density: %g\n",(double)crit);
        for (i = 0; i < m; ++i) {
            b = a + (double)i * d;
            fprintf(PMFd,"%5d ",i + 1);
            fprintf(PMFd,PMFmtS,(double)b);
            fprintf(PMFd,PMFmtS,(double)x[i]);
            fprintf(PMFd,"\n");
        }
    }  
    if (PMPLFlg && PSFFlg) {                  /* plot output required */

/*      fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);      */
/*      fprintf(PSFd,"\n%%#density\n");         */
        fprintf(PSFd,"gsave\n");

        ps_ltyp(PMLT);     /* set line type */
        ps_lwidth(PMLW);   /* and line width */

        ps_2dplot((double)a,(double)(PMSC * x[0]),0);    
         
        for (i = 1; i < m; ++i) {
            b = a + (double)i * d;
            ps_2dplot((double)b,(double)(PMSC * x[i]),1);    
        }
        fprintf(PSFd,"stroke\n");
        fprintf(PSFd,"grestore\n");
    }
    err = 0;

KDENSFin:
    if (x_a) {
        free((char *)x);
        memrq(-mn,sizeof(float));
    }
    p_clean();          /* close files and free memory in t_parm */
    newline();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  kdensf(n,x,m,a,d,sig,crit)                                              */
/*                                                                          */
/*  Density estimation with a normally distributed kernel.                  */
/*                                                                          */
/*  x[i] (i=0,n-1) are the input points. Return density.                    */
/*                                                                          */
/*  x[0] =  density at a                                                    */
/*  x[i] =  density at a + i * d     i = 0,1,...,m - 1                      */
/*                                                                          */
/*  The maximum value of densities will be returned in crit.                */
/*                                                                          */
/*  Return: 0 if successful, or -1 if insufficient memory.                  */

int kdensf(int n,float *x,int m,float a,float d,float sig,float *crit)
{
    register int i,j,k;
    int err,jflg,ptr_a;
    int *ptr;
    float max,t,y,tmp,critv,dmax;

    critv = *crit;
    ptr_a = 0;
    err = -1;

    if (alloc_acxf(NOC + 1))
        goto KDENSFin;

    if (!(ptr = (int *)calloc(n,sizeof(int))))  
        goto KDENSFin; 
    memrq(n,sizeof(int));
    ptr_a = 1;

    /* sort x */

    for (i = 0; i < n; ++i) {
        AcXF[i] = x[i];
        x[i] = 0.0;
        ptr[i] = i;
    }
    qsort((char *)ptr,n,sizeof(int),kdcomp);

    dmax = 0.0;
    max = AcXF[ptr[n - 1]]; 
    j = jflg = 0;
    t = a;

    for (k = 0; k < m; ++k) {

        y = 0.0;

        if (jflg >= 0) {
            jflg = 0;
            for (i = j; i < n; ++i) {
                tmp = (float)dnf((t - AcXF[ptr[i]]) / sig) / sig;
                if (tmp < critv) {
                    if (jflg)  
                        break;
                }
                else {
                    y += tmp;                   
                    if (jflg == 0) {
                        j = i;
                        jflg = 1;
                    }
                }
                /**
                if (jflg == 0 && t >= max)
                    jflg = -1;
                **/
            }
        }
        if (dmax < y)
            dmax = y;

        x[k] = y;
        t += d;
    }
    *crit = dmax;
    err = 0;

KDENSFin:
    if (ptr_a) {
        free((char *)ptr);
        memrq(-n,sizeof(int));
    }
    return(err);
}

/****************************************************************************/
/*  kdcomp()                                                                */
/*      Compare function used by kdensf().                                  */

int kdcomp(const void *arg1,const void *arg2)
{     
    float tmp;
    tmp = AcXF[*(int *)arg1] - AcXF[*(int *)arg2];
    if (tmp > 0.0)
        return(1);
    else if (tmp < 0.0)
        return(-1);
    return(0);
}


