/****************************************************************************/
/*  t_xplot                                                                 */
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
#include "t_plot.h"
#include "t_psinit.h"
#include "t_psf.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_sort.h"
#include "t_alloc.h"
#include "t_gf.h"
#include "t_lsei.h"
#include "t_l1reg.h"
#include "t_areg.h"
#include "t_eval.h"
#include "t_int.h"
#include "t_freq.h"
#include "t_graph.h"
#include "t_l1reg.h"
#include "t_gmin.h"
#include "t_gm.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_xplot.c                                                  */

int xplot(void);
int x_cpsf(double xmin,double xmax,double ymin,double ymax,int opt);
int alloc_scx(int n);
void x_psetup(double xa,double xb,double ya,double yb);
int x_getn(double x);
int x_getaxval(double xmin,double xmax,double *xa,double *xb,int opt);
void x_plaxis(int opt,int nx,double xa,double xb);
void x_scplot(int n,int ns,int opt,int clip,char *cmd);
int x_cchk(double xmin,double xmax,double ymin,double ymax);
int xcheck(void);
int xlog(void);
int xlog1(void);
int xlog_check(void);
void xlog_prn(int opt);
int xlogp(void);
int xconh(void);
int xconhp(int g,int lt,char *cmd);
int xreg(void);
int xplotf(void);
int xopen(void);
int xdelete(void);
int xdens(void);
int xfunc(void);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int SCN = 0;            /* number of data points                            */
double *SCXVar;         /* X variable for scatterplot                       */
double *SCYVar;         /* Y variable for scatterplot                       */
int *SCGVar;            /* group variable for scatterplot                   */
int SCXVarN = 0;        /* allocated for SCXVar                             */
int SCYVarN = 0;        /* allocated for SCYVar                             */
int SCGVarN = 0;        /* allocated for SCGVar                             */
double SCXMin = 0.0;    /* minimum of X variable                            */
double SCYMin = 0.0;    /* minimum of Y variable                            */
double SCXMax = 0.0;    /* maximum of X variable                            */
double SCYMax = 0.0;    /* maximum of Y variable                            */
int SCNGRP = 0;         /* number of groups                                 */
int SCGRP[MaxG];        /* array with group numbers                         */

short SYMLIST[17] = {11,13,12,10,9,1,2,16,3,6,8,17,14,4,5,7,15};

/*--------------------------------------------------------------------------*/
/*  xplot()     Scatterplot, optionally with grouping.                      */
/*              If there is no currently active PostScript file, the        */
/*              command creates a new one.                                  */
/*                                                                          */
/*              xplot(                                                      */
/*                  opt=...,        1 scatterplot, 2 line plot, def. 1      */
/*                  s=...,                                                  */
/*                  fs=...,         symbol size, def. 1.3                   */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  pxlen=...,      length of x axis, def. 120 mm           */
/*                  pylen=...,      length of y axis, def.  80 mm           */
/*                                                                          */
/*              ) = X,Y [,G];                                               */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xplot(void)
{
    register int i,j;
    int err,ix,iy,ig,n;

    ig = err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XPlot. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto XPLOTFin;

    PMOPT--;
    if (PMOPT != 1)
        PMOPT = 0;

    if (PMNV < 2 || PMNV > 3) {
        printf1("Error: need two or three variables.\n");
        goto XPLOTFin;
    }
    ix = PMVIdx[0];
    iy = PMVIdx[1];
    if (PMNV == 3)
        ig = PMVIdx[2];

    printf1("Scatterplot with: %s, %s\n",VName[ix],VName[iy]);
    if (ig >= 0)
        printf1("Grouping with: %s\n",VName[ig]);

    if (alloc_scx(NOC))
        goto XPLOTFin;

    SCNGRP = 1;                 /* number of groups */
    SCGRP[0] = 1;

    if (ig >= 0) {
        for (i = 0; i < NOC; ++i)  
            SCGVar[i] = (int)get_data(ig,i);       

        if (sorti(NOC,SCGVar,0))
            goto XPLOTFin;

        SCGRP[0] = SCGVar[0];
        for (i = 1; i < NOC; ++i) {
            if (SCGVar[i] != SCGVar[i - 1]) {
                if (SCNGRP < MaxG) 
                    SCGRP[SCNGRP++] = SCGVar[i];
                else {
                    printf1("Warning: exceeded max number of groups.\n");
                    break;
                }
            }
        }
        printf1("Number of groups: %d  [",SCNGRP);
        for (i = 0; i < SCNGRP; ++i)
            printf1(" %d",SCGRP[i]);
        printf1("]\n");
    }
    SCXMin = SCXMax = get_data(ix,0);
    SCYMin = SCYMax = get_data(iy,0);

    for (i = 0; i < NOC; ++i) {

        SCXVar[i] = get_data(ix,i);       
        SCYVar[i] = get_data(iy,i);       
        SCXMin = dmin(SCXMin,SCXVar[i]); 
        SCXMax = dmax(SCXMax,SCXVar[i]); 
        SCYMin = dmin(SCYMin,SCYVar[i]); 
        SCYMax = dmax(SCYMax,SCYVar[i]); 

        if (ig >= 0)
            SCGVar[i] = (int)get_data(ig,i);       
        else
            SCGVar[i] = 1;
    }
    SCN = NOC;

    printf1("X Min %lg  Max %lg\n",SCXMin,SCXMax);
    printf1("Y Min %lg  Max %lg\n",SCYMin,SCYMax);

    if (PSFFlg == 2) {
        xlog_prn(1);            /* using current PostScript file */

        if (SCXMin >= PA2[0] || SCXMax <= PA1[0] ||
            SCYMin >= PA2[1] || SCYMax <= PA1[1]) {  
            printf1("Empty intersection with current coordinates.\n");
            err = 0;
            goto XPLOTFin;
        }
    }
    else {
        if (x_cpsf(SCXMin,SCXMax,SCYMin,SCYMax,0))
            goto XPLOTFin;
    }

    /* make scatterplots for all groups */
   
    if (PMOPT == 0) {               /* scatterplot */

        if (PMFSFlg == 0)
            PMFS = 1.3;                 /* default size of symbols */

        if (SCNGRP == 1) {
            if (PMS >= 1 && PMS <= 17)
                n = PMS;
            else
                n = 1;
        }
        else
            n = 0;

        i = 0;
        for (j = 0; j < SCNGRP; ++j) {
            if (n)
                x_scplot(SCGRP[j],n,PMOPT,1,CmdBuf);          
            else {
                x_scplot(SCGRP[j],SYMLIST[i],PMOPT,1,CmdBuf);          
                if (++i >= 17)
                    i = 0;
            }
        }
    }
    else {                      /* line plot */

        if (PMLT < 1 || PMLT > 9)
            PMLT = 1;

        for (j = 0; j < SCNGRP; ++j) {
            x_scplot(SCGRP[j],PMLT,PMOPT,1,CmdBuf);          
            if (++PMLT > 9)
                PMLT = 1;  
        }
    }
    err = 0;

XPLOTFin:
    if (err)
        alloc_scx(0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  x_cpsf()    Create new PostScript file with standard coordinate system. */  
/*                                                                          */
/*  Return: 0 if OK, otherwise -1                                           */

int x_cpsf(double xmin,double xmax,double ymin,double ymax,int opt)
{
    int nx,ny;
    double xa,xb,ya,yb;

    if (x_cchk(xmin,xmax,ymin,ymax))
        return(-1);     
      
    nx = x_getaxval(xmin,xmax,&xa,&xb,opt);     /* get values for axes */
    ny = x_getaxval(ymin,ymax,&ya,&yb,opt);

    /* setup new Postscipt file */

    newline();
    strcpy(PSFName,"xplot.ps");

    if (!(PSFd = fopen(PSFName,OPEN_WR))) {
        printf1("Error: can't create: %s\n",PSFName);
        return(-1);    
    }
    ps_init();          /* write header to output file */
    PSFFlg = 1;
    printf1("New PostScript file: %s\n",PSFName);

    x_psetup(xa,xb,ya,yb);      /* set up coordinate system */

    x_plaxis(0,nx,xa,xb);
    x_plaxis(1,ny,ya,yb);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_scx(n)    If n > 0 allocate SCXVar,SCYVar,SCGVar, otherwise free. */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_scx(int n)
{
    if (SCXVarN > 0) {
        free((char *)SCXVar);
        memrq(-SCXVarN,sizeof(double));
        SCXVarN = 0;
    }
    if (SCYVarN > 0) {
        free((char *)SCYVar);
        memrq(-SCYVarN,sizeof(double));
        SCYVarN = 0;
    }
    if (SCGVarN > 0) {
        free((char *)SCGVar);
        memrq(-SCGVarN,sizeof(int));
        SCGVarN = 0;
    }
    SCN = 0;

    if (n > 0) {
        if (!(SCXVar = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        SCXVarN = n;
        if (!(SCYVar = (double *)calloc(n,sizeof(double)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(double));
        SCYVarN = n;
        if (!(SCGVar = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        SCGVarN = n;
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  x_psetup()  Set up coordinate system.                                   */

void x_psetup(double xa,double xb,double ya,double yb)
{
    PSLog[0] = PSLog[1] = 0;  

    XOrg = 150;         /* origin of PostScript figure */
    YOrg = 460;  
    SCALXFac = 1.0;     /* PostScript x scaling factor */
    SCALYFac = 1.0;     /* PostScript y scaling factor */
    ROTFac = 0.0;       /* PostScript rotation factor */

    PXLen = PMXLen;     /* phys length of X axis in mm */
    PYLen = PMYLen;     /* phys length of Y axis in mm */

    PA1[0] = xa;
    PA2[0] = xb;
    PA1[1] = ya;
    PA2[1] = yb;
         
    UXLen = PA2[0] - PA1[0];
    UYLen = PA2[1] - PA1[1];

    PSXLen = PtMM * PXLen;   /* phys. x axis length points */
    PSYLen = PtMM * PYLen;   /* phys. y axis length points */

    printf1("Size (width x height in mm): %g x %g\n",PXLen,PYLen);
    printf1("PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    XOrg,XOrg + (int)(PSXLen + 0.5),YOrg,YOrg + (int)(PSYLen + 0.5));   
    printf1("Origin (PostScript x and y coordinates): %d,%d\n",XOrg,YOrg);

    if (PSFFlg == 1) {      /* init bounding box parameters */
        BBLX = XOrg;
        BBLY = YOrg;
        BBUX = BBLX + (int)(PSXLen + 0.5);
        BBUY = BBLY + (int)(PSYLen + 0.5);
    }
    ps_inis(1);

    printf1("User coordinates X axis: %lg,%lg\n",PA1[0],PA2[0]);
    printf1("User coordinates Y axis: %lg,%lg\n",PA1[1],PA2[1]);

    PSFFlg = 2;     /* set to 2 if setup successful */
}

/*--------------------------------------------------------------------------*/
/*  x_getn()                                                                */
/*                                                                          */

int x_getn(double x) 
{
    int n; double tmp;

    tmp = pow(10.0,-10.0);
    for (n = -10; n <= 10; ++n) {
        if (tmp <= x && x <= 10.0 * tmp)   
            break;
        tmp *= 10.0;
    }
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  x_getaxval(xmin,xmax,xa,xb,opt)                                         */
/*                                                                          */
/*  calculate reasonable min and max values for axis. xmin and xmax are     */
/*  the min and max x values, respectively. Return xa and xb, and also      */
/*  n such that 10^-n is the unit.                                          */
/*                                                                          */
/*  If opt=0 check for boundaries.                                          */

int x_getaxval(double xmin,double xmax,double *xa,double *xb,int opt)
{
    int n,nn;
    double d;

    d = (xmax - xmin) / 5.0;
    if (d < EPSI1)
        d = EPSI1;

    /* printf1("xmin=%lg xmax=%lg d=%lg\n",xmin,xmax,d); */
    n = x_getn(d);
    d = pow(10.0,(double)n);

    /* printf1("nach n=%d d=%lg\n",n,d); */  

    if (xmin >= 0.0 && xmin <= EPSI1)
        nn = 0;
    else 
        nn = (int)floor(xmin / d);
    *xa = (double)nn * d;

    if (opt == 0) {
        if (fabs(*xa - xmin) <= EPSI1) 
            *xa -= d;
    }
    if (xmax <= 0.0 && xmax >= -EPSI1)
        nn = 0;
    else  
        nn = (int)ceil(xmax / d);
    *xb = (double)nn * d;

    if (opt == 0) {
        if (fabs(*xb - xmax) <= EPSI1) 
            *xb += d;
    }
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  x_plaxis(opt,nx,xa,xb)      plot axis for coordinate system.            */
/*                              opt 0 : x axis                              */
/*                              opt 1 : y axis                              */

void x_plaxis(int opt,int nx,double xa,double xb)         
{
    PMIC = 0.0;
    PMSC = pow(10.0,(double)nx);
    if ((xb - xa) / PMSC > 20) {    
        PMSC = pow(10.0,(double)(nx+1));
        PMIC = 10.0;
    }
    pl_axis(opt,1,1,1,0.2);
}

/*--------------------------------------------------------------------------*/
/*  x_scplot(n,ns,opt,clip,cmd)                                             */  
/*                                                                          */
/*                      make scatterplot for group n with symbol ns.        */
/*                      if opt = 0 scatterplot, ns is symbol type.          */
/*                      if opt = 1 then line plot, ns is line type.         */
/*                      if clip != 0 use clip option.                       */
/*                      cmd is the plot object command.                     */
/*                                                                          */

void x_scplot(int n,int ns,int opt,int clip,char *cmd)
{
    register int i;
    int first;
    double x,y;

    fprintf(PSFd,"\n%%#%d: %s [group %d]\n",++PSONUM,cmd,n);
    fprintf(PSFd,"gsave\n");
   
    if (clip)
        set_clip();

    ps_lwidth(PMLW);
    if (opt)
        ps_ltyp(ns);
    else
        ps_ltyp(1);

    if (opt == 0) {
        fprintf(PSFd,"gsave\n");
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);
    }
    first = 1;
    for (i = 0; i < SCN; ++i) {
        if (SCGVar[i] == n) {
            x = ps_2dx(SCXVar[i]);       
            y = ps_2dy(SCYVar[i]);       
            if (opt == 0)
                ps_sym(ns,x,y,PtMM * PMFS / 2.0);
            else if (first) {
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else              
                fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
        }
    }
    if (opt == 1)
        fprintf(PSFd,"stroke\n");
    else
        fprintf(PSFd,"grestore\n");
    fprintf(PSFd,"grestore\n");
}

/*--------------------------------------------------------------------------*/
/*  x_cchk      Check for possibility of coordinate system.                 */
/*              Return 0 if OK, -1 if not.                                  */

int x_cchk(double xmin,double xmax,double ymin,double ymax)
{
    if (fabs(xmin - xmax) < EPSI1 || fabs(ymin - ymax) < EPSI1) {
        printf1("\nError: can't create a coordinate system.\n");
        return(-1);    
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  xcheck      Check if basic data for xplot are available.                */
/*              Return 0 if OK, -1 if not.                                  */

int xcheck(void)
{
    if (SCN > 0)  
        return(0);
    printf1("Error: no current scatter plot.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  xlog        Print list of plot objects.                                 */
/*              xlog [=psfile];                                             */
/*              Return 0 if OK, -1 if not.                                  */

int xlog(void)
{
    int err, oflag;

    err = -1;
    oflag = 0;
    if (check_cmd(1))
        return(-1);

    printf1("List of plot objects. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,2,0))     /* get parameters */
        goto XLOGFin;

    if (PMFDef == 0) {
        if (xlog_check() != 0)
            goto XLOGFin;
        xlog_prn(1);

        fflush(PSFd);
        if (!(PMFd = fopen(PSFName,OPEN_RD))) {   
            printf1("Error: can't re-open this file.\n");
            goto XLOGFin;
        }
        oflag = 1;
    }
    else  
        printf1("Using PostScript file: %s\n",PMFdName);

    if (alloc_acc(RLMaxDef))
        goto XLOGFin;

    if (xlogp() == 0)            /* print list */
        printf1("No plot objects.\n");
    err = 0;

XLOGFin:
    if (oflag) {
        fclose(PMFd);
        PMFDef = 0;
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xlog1       Print current ps file and list of plot objects.             */
/*              Return 0 if OK, -1 if not.                                  */

int xlog1(void)
{
    int err;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (PSFFlg != 2) {
        printf1("No current PostScript file.\n");
        return(0);
    }
    xlog_prn(1);
    printf1("\nList of plot objects. ");
    fflush(PSFd);
    if (!(PMFd = fopen(PSFName,OPEN_RD))) {   
        printf1("Error: can't re-open this file.\n");
        goto XLOG1Fin;
    }
    if (alloc_acc(RLMaxDef))
        goto XLOG1Fin;

    if (xlogp() == 0)            /* print list */
        printf1("No plot objects.\n");
    fclose(PMFd);
    PMFDef = 0;
    err = 0;

XLOG1Fin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xlog_check  Check for current plot file. Return 0 if OK, -1 if error.   */

int xlog_check(void)
{
    if (PSFFlg != 2) {
        printf1("Error: no current plot file.\n");
        return(-1);   
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  xlog_prn    Print current plot file and, optionally, coordinates.       */

void xlog_prn(int opt)
{
    printf1("Using current plot file: %s\n",PSFName);
    if (opt)
        printf1("X axis: %lg,%lg  Y axis: %lg,%lg\n",PA1[0],PA2[0],PA1[1],PA2[1]);
}

/*--------------------------------------------------------------------------*/
/*  xlogp       Print list of plot objects.                                 */
/*              Return number of plot objects.                              */

int xlogp(void)
{
    int n,m,first;
    register char *p;

    m = 0; first = 1;
    fseek(PMFd,0,0);
    while (fgets(AcC,RLMaxDef,PMFd)) {
        p = AcC;
        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            p = skip_int(p + 2);
            if (first) {
                newline();
                first = 0;
            }
            printf1("%3d%s",n,p);
            m++;
        }
    }
    return(m);
}

/*--------------------------------------------------------------------------*/
/*  xconh       Convex hull. Uses the currently defined XData structure.    */
/*                                                                          */
/*              xconh(                                                      */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2mm                  */
/*              );                                                          */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xconh(void)
{
    register int i,j;
    int err,n;

    err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XConh. Current memory: %d bytes.\n",MemReq);
    if (xcheck())
        goto XCONHFin;

    if (parm(CmdBuf + 5,3,0))     /* get parameters */
        goto XCONHFin;

    printf1("Adding %d convex hull(s).\n",SCNGRP);

    j = PMLT;
    if (j < 1 || j > 9)
        j = 1;

    for (i = 0; i < SCNGRP; ++i) {
        n = xconhp(SCGRP[i],j,CmdBuf);
        if (n)  
            printf1("Error in group %d\n",SCGRP[i]);
        if (PMLTFlg == 0) {
            if (++j > 9)
                j = 1;
        }   
    }
    err = 0;
   
XCONHFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xconhp(g,lt,cmd)   Plot convex hull for group g with line type lt.      */
/*                     Use cmd for plot object description.                 */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xconhp(int g,int lt,char *cmd)
{
    register int i,j,ik,l;
    int err,nn,n;
    double x,y;

    err = -1;

    if (alloc_acxf(SCN + 1))
        goto XCONHPFin;
    if (alloc_acyf(SCN + 1))
        goto XCONHPFin;
    if (alloc_acn(SCN + 1))
        goto XCONHPFin;
    if (alloc_aci(SCN + 1))
        goto XCONHPFin;
    if (alloc_acj(SCN + 1))
        goto XCONHPFin;

    nn = 0;
    for (i = 0; i < SCN; ++i) {
        if (SCGVar[i] == g) {
            nn++;
            AcXF[nn] = (float)SCXVar[i];            
            AcYF[nn] = (float)SCYVar[i];            
            AcN[nn] = nn;
        }
    }
    if (nn < 1)
        goto XCONHPFin;

    n = g_chull(AcXF,AcYF,nn,AcN,AcJ,AcI);   
    if (n <= 0) {
        if (n < 0)
            p_err(-2,1);
        goto XCONHPFin;
    }

    fprintf(PSFd,"\n%%#%d: %s [group %d]\n",++PSONUM,cmd,g);
    fprintf(PSFd,"gsave\n");
    ps_lwidth(PMLW);
    ps_ltyp(lt);

    ik = AcI[1];
    l = 0;
    for (i = 1; i <= n; ++i) {
        j = AcJ[ik];
        x = (double)AcXF[j];
        y = (double)AcYF[j];
        ps_2dplot(x,y,l); 
        l = 1;
        ik = AcI[ik];
    }
    ik = AcI[1];
    j = AcJ[ik];
    ps_2dplot((double)AcXF[j],(double)AcYF[j],1);     /* close path */

    /********************************************
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");       
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");   
    }
    *********************************/
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
    err = 0;

XCONHPFin:
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xreg        Plot regression curve.                                      */
/*                                                                          */
/*              xreg(                                                       */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  sig=...,        sigma for lowess, def. 0.5              */
/*                                                                          */  
/*              ) = n;              1 : lowess                              */
/*                                  2 : least squares                       */
/*                                  3 : l1 norm regression                  */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xreg(void)
{
    register int i;
    int err,typ,n,ig,g,first,ra,re,r,nx;
    double x,y,xmin,xmax,rne,rnl,xx[4];

    err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XREG. Current memory: %d bytes.\n",MemReq);
    if (xcheck())
        goto XREGFin;

    if (parm(CmdBuf + 4,3,0))     /* get parameters */
        goto XREGFin;

    switch (PMRHSI) {
        case 2:         printf1("Least squares regression.\n");
                        typ = 1;
                        nx = 3;
                        break;
        case 3:         printf1("L1 norm regression.\n");
                        typ = 2;
                        nx = 2;
                        break;
        default:        printf1("Lowess.\n");
                        typ = 3;
                        nx = 1;
                        break;
    }
    if (PMSIG <= 0.0 || PMSIG >= 1.0)
        PMSIG = 0.5;

    if (alloc_acx(nx * SCN + 1))
        goto XREGFin;

    if (typ > 1) {
        if (alloc_acy(SCN + 1))
            goto XREGFin;
    }
    if (typ == 3) {
        if (alloc_acu(SCN + 1))
            goto XREGFin;
        if (alloc_acv(SCN + 1))
            goto XREGFin;
        if (alloc_acw(SCN + 1))
            goto XREGFin;
    }
    if (PMLT < 1 || PMLT > 9)
        PMLT = 1;

    for (ig = 0; ig < SCNGRP; ++ig) {
        xmin = xmax = 0.0;
        n = 0;
        g = SCGRP[ig];
        first = 1;
        for (i = 0; i < SCN; ++i) {
            if (SCGVar[i] == g) {
                x = SCXVar[i];            
                y = SCYVar[i];            
                if (first) {
                    xmax = xmin = x;
                    first = 0;
                }   
                else {
                    xmin = dmin(xmin,x);
                    xmax = dmax(xmax,x);
                }
                if (typ == 1) {
                    AcX[n * nx + 1] = 1.0;
                    AcX[n * nx + 2] = x;
                    AcX[n * nx + 3] = y;
                }
                else if (typ == 2) {
                    AcX[n * nx + 1] = 1.0;
                    AcX[n * nx + 2] = x;
                    AcY[n + 1] = y;
                }
                else if (typ == 3) {
                    AcX[n + 1] = x;
                    AcY[n + 1] = y;
                    AcW[n + 1] = 1.0;
                }
                n++;
            }
        }
        if (n < 2 || xmax < xmin + EPSI1)
            continue;

        printf1("Group %3d : ",g);

        switch (typ) {
            case 1:         /* least squares regression */

                r = lsei(AcX,0,n,0,2,xx,0,&rne,&rnl,&ra,&re);
                if (r)
                    printf1("error.\n");
                else  
                    printf1("rank=%d  alpha=%lg  beta=%lg  norm of residuals=%lg\n",ra,xx[1],xx[2],rnl);
                break;

            case 2:         /* l1 norm regression */

                r = l1regf(n,2,AcX,AcY,xx,&rnl,&ra);
                if (r)
                    printf1("error.\n");
                else  
                    printf1("rank=%d  alpha=%lg  beta=%lg  sum of residuals=%lg\n",ra,xx[1],xx[2],rnl);
                break;

            case 3:         /* lowess regression */

                if (sortd2(n,AcX + 1,AcY + 1))     /* sort simultaneously according to x */
                    r = -1;            
                else  
                    r = lowess(AcX,AcY,n,PMSIG,2,0.0,AcU,AcW,AcV);
                if (r)
                    printf1("error.\n");
                else
                    printf1("sigma=%lg\n",PMSIG);
                break;

            default: break;
        }
        if (r)
            continue;
   
        fprintf(PSFd,"\n%%#%d: %s [group %d]\n",++PSONUM,CmdBuf,g);
        fprintf(PSFd,"gsave\n");
        ps_lwidth(PMLW);
        ps_ltyp(PMLT);

        if (typ == 3) {
            ps_2dplot(AcX[1],AcU[1],0);
            for (i = 2; i <= n; ++i)  
                ps_2dplot(AcX[i],AcU[i],1);
        }
        else {
            x = xx[1] + xx[2] * xmin;
            ps_2dplot(xmin,x,0);
            x = xx[1] + xx[2] * xmax;
            ps_2dplot(xmax,x,1);
        }
        fprintf(PSFd,"stroke\n");
        fprintf(PSFd,"grestore\n");

        if (PMLTFlg == 0) {
            if (++PMLT > 9)
                PMLT = 1;
        }
    }
    err = 0;
   
XREGFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xplotf      Scatterplot or line plot that takes the values directly     */
/*              from a file.                                                */
/*                                                                          */
/*              xplotf(                                                     */
/*                  cn=...,         column number for variables, def. 1,2   */
/*                  gn=...,         column number for group variable,       */
/*                                  default no grouping.                    */
/*                  opt=...,        1 scatterplot, 2 line plot, def. 1      */
/*                  s=...,                                                  */
/*                  fs=...,         symbol size, def. 1.3                   */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  pxlen=...,      length of x axis, def. 120 mm           */
/*                  pylen=...,      length of y axis, def.  80 mm           */
/*                                                                          */
/*              ) = input data file;                                        */
/*                                                                          */
/*                                                                          */
/*  There are two options for reading the data.                             */
/*  a) gn not defined. Then                                                 */
/*                                                                          */
/*     cn=c1,c2,...,    c1 defines X values, each of the following          */
/*                      columns defines one group.                          */
/*                                                                          */
/*  b) gn is used. Then each block of values in this column defines one     */
/*     group. In this case, cn=c1,c2 must only define two columns, the      */
/*     first one for X, the second one for Y.                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int xplotf(void)
{
    register int i,j,k;
    int err,n,nn,cmax,cn,ng,cx,cy;
    double tmp,tmp1;              
    register char *p;

    err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XPLOTF. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,2,1))     /* get parameters */
        goto XSCFFin;

    PMOPT--;
    if (PMOPT != 1)
        PMOPT = 0;

    cmax = 1;
    if (PMOPT == 0)
        printf1("Scatterplot. ");
    else if (PMOPT == 1)
        printf1("Line plot. ");
    printf1("File: %s  ",PMFdName);
    printf1("Columns:");
    if (PMNCN > 0) {
        for (i = 0; i < PMNCN; ++i) {
            printf1(" %d",PMCN[i]);
            if (cmax < PMCN[i])
                cmax = PMCN[i];
        }
        newline();
        if (PMNCN < 2) {
            printf1("Error: need at least two columns.\n");
            goto XSCFFin;
        }
        cn = PMNCN;
        cx = PMCN[0];
        cy = PMCN[1];
    }
    else {
        printf1(" 1 2\n");
        cn = cmax = 2;
        cx = 1;
        cy = 2;
    }
    if (PMGN > 0) {
        printf1("Group definition by column %d\n",PMGN);
        if (cn > 2) {
            printf1("Error: cn must provide exactly two columns.\n");
            goto XSCFFin;
        }
        if (cmax < PMGN)
            cmax = PMGN;
    }

    if (alloc_acc(RLMaxDef + 1))
        goto XSCFFin;

    i = ng = n = 0;
    while (fgets(AcC,RLMaxDef,PMFd)) {
        i++;
        if (check_drec(AcC)) {      /* check for data records */
            n++;

            if (PMGN > 0) {
                p = AcC;
                for (j = 1; j <= PMGN; ++j) {
                    p = skip_sep(p);                    /* skip separator */
                    if (sscanf(p,"%lg",&tmp) != 1) {
                        printf1("Error: can't read c%d in record %d.\n",j,i);
                        goto XSCFFin;
                    }
                    if (j == PMGN) {
                        if (n == 1) {
                            ng = 1;
                            tmp1 = tmp;
                        }
                        else {
                            if (tmp1 != tmp) {
                                ng++;
                                tmp1 = tmp;
                            }
                        }
                    }
                    p = skip_dval(p);
                }
            }
        }
    }
    if (n < 1) {
        printf1("Error: file contains no data records.\n");
        goto XSCFFin;
    }
    if (PMGN == 0) {
        ng = cn - 1;
        nn = ng * n;
    }
    else
        nn = n;

    if (alloc_scx(nn))
        goto XSCFFin;

    if (alloc_actmp(cmax + 1))
        goto XSCFFin;

    fseek(PMFd,0,0);
    k = i = 0;
    while (fgets(AcC,RLMaxDef,PMFd)) {
        i++;
        if (check_drec(AcC)) {      /* check for data records */
            p = AcC;
            for (j = 1; j <= cmax; ++j) {
                p = skip_sep(p);                    /* skip separator */
                if (sscanf(p,"%lg",&tmp) != 1) {
                    printf1("Error: can't read c%d in record %d.\n",j,i);
                    goto XSCFFin;
                }
                AcTmp[j] = tmp;
                p = skip_dval(p);
            }
            if (cn == 2) {
                SCXVar[k] = AcTmp[cx];
                SCYVar[k] = AcTmp[cy];
                if (PMGN > 0)
                    SCGVar[k] = (int)AcTmp[PMGN];
                else
                    SCGVar[k] = 1;            
                k++;
            }
            else {
                for (j = 1; j < cn; ++j) {
                    SCXVar[(j - 1) * n + k] = AcTmp[cx];
                    SCYVar[(j - 1) * n + k] = AcTmp[PMCN[j]];
                    SCGVar[(j - 1) * n + k] = j;         
                }
                k++;
            }
        }
    }
    alloc_acc(0);

    SCNGRP = 1;
    k = SCGVar[0];
    SCGRP[0] = k;

    for (i = 1; i < nn; ++i) {
        if (SCGVar[i] != k) {
            if (SCNGRP >= MaxG) {
                printf1("Error: exceeded max number of groups.\n");
                goto XSCFFin;
            }
            k = SCGVar[i];          
            SCGRP[SCNGRP++] = k;
        }
    }
    SCN = nn;

    printf1("Number of groups: %d  [",SCNGRP);
    for (i = 0; i < SCNGRP; ++i)
        printf1(" %d",SCGRP[i]);
    printf1("]\n");
       
    SCXMin = SCXMax = SCXVar[0];      
    SCYMin = SCYMax = SCYVar[0];     

    for (i = 1; i < SCN; ++i) {
        SCXMin = dmin(SCXMin,SCXVar[i]); 
        SCXMax = dmax(SCXMax,SCXVar[i]); 
        SCYMin = dmin(SCYMin,SCYVar[i]); 
        SCYMax = dmax(SCYMax,SCYVar[i]); 
    }
    printf1("X Min %lg  Max %lg\n",SCXMin,SCXMax);
    printf1("Y Min %lg  Max %lg\n",SCYMin,SCYMax);

    if (PSFFlg == 2) {
        xlog_prn(1);            /* using current PostScript file */

        if (SCXMin >= PA2[0] || SCXMax <= PA1[0] ||
            SCYMin >= PA2[1] || SCYMax <= PA1[1]) {  
            printf1("Empty intersection with current coordinates.\n");
            err = 0;
            goto XSCFFin;
        }
    }
                /* new PostScript file */
    else {
        if (x_cpsf(SCXMin,SCXMax,SCYMin,SCYMax,PMOPT))
            goto XSCFFin;
    }

    /* make scatterplots for all groups */

    if (PMOPT == 0) {               /* scatterplot */

        if (PMFSFlg == 0)
            PMFS = 1.3;                 /* default size of symbols */
        if (SCNGRP == 1) {
            if (PMS >= 1 && PMS <= 17)
                n = PMS;
            else
                n = 1;
        }
        else
            n = 0;

        i = 0;
        for (j = 0; j < SCNGRP; ++j) {
            if (n)
                x_scplot(SCGRP[j],n,PMOPT,1,CmdBuf);          
            else {
                x_scplot(SCGRP[j],SYMLIST[i],PMOPT,1,CmdBuf);          
                if (++i >= 17)
                    i = 0;
            }
        }
    }
    else {                      /* line plot */

        if (PMLT < 1 || PMLT > 9)
            PMLT = 1;

        for (j = 0; j < SCNGRP; ++j) {
            x_scplot(SCGRP[j],PMLT,PMOPT,1,CmdBuf);          
            if (++PMLT > 9)
                PMLT = 1;  
        }
    }
    err = 0;

XSCFFin:
    if (err)
        alloc_scx(0);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xf          Plot function.                                              */
/*                                                                          */
/*              xf(                                                         */
/*                  rx-a(d)b,       optional range for x axis, def. 1-10    */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2                    */
/*                  pxlen=...,      length of x axis, def. 120 mm           */
/*                  pylen=...,      length of y axis, def.  80 mm           */
/*                  nc=...,         1 for no clipping, def. 0               */
/*              ) = function;                                               */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xfunc(void)
{
    register int i;
    int err,n,r,nflg;    
    double x,y,xmin,xmax,ymin,ymax;

    err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XF. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 2,9,1))    /* get parameters */
        goto PXFFin;

    if (PMLT < 1 || PMLT > 9)
        PMLT = 1;

    if (FNFlg == 0) {
        p_err(-1,1);
        goto PXFFin;
    }
    if (FNArgN != 1) {
        printf1("Error: there must be exactly one function argument.\n");
        goto PXFFin;
    }
    if (PMFTYP5) {
        p_err(-40,1);
        goto PXFFin;
    }
    if (FVFlg)
        printf1("Sum over %d data matrix cases.\n",NOC);

    if (PSFFlg == 2) {
        nflg = 0;
        PMRXA = PA1[0];
        PMRXB = PA2[0];
        PMRXD = (PMRXB - PMRXA) / 200.0;
        n = 200;
    }
    else {
        nflg = 1;
        if (PMRXFlg == 0) {
            PMRXA = 1.0;   
            PMRXB = 10.0;   
            PMRXD = (PMRXB - PMRXA) / 200.0;
            n = 200;
        }
        else
            n = (int)((PMRXB - PMRXA + EPSI1) / PMRXD);
    }
    xmin = PMRXA;
    xmax = PMRXB;

    if (alloc_acx(n + 1))
        goto PXFFin;
    if (alloc_acy(n + 1))
        goto PXFFin;

    NINTMUsed = -1;

    for (i = 0; i <= n; ++i) {

        FNArgVal[0] = AcX[i] = PMRXA + (double)i * PMRXD;

        r = get_flval(&y,FNArgN,FNArgVal,0,0,AcU,AcV,AcV);

        if (r) {     /* r from v_eval1() */  

            printf1("Can't evaluate function for x=%lg.\n",x);
            prn_emsg2(r);
            goto PXFFin;
        }
        AcY[i] = y;
        if (i == 0)
            ymin = ymax = y;
        else {
            ymin = dmin(ymin,y);
            ymax = dmax(ymax,y);
        }
    }
    if (NINTMUsed >= 0)  
        ni_info();

    printf1("\nX Min %lg  Max %lg\n",xmin,xmax);
    printf1("Y Min %lg  Max %lg\n",ymin,ymax);
   
    if (nflg) {         /* create new PostScript file */

        if (x_cpsf(xmin,xmax,ymin,ymax,1))
            goto PXFFin;
    }
    else
        printf1("\nUsing current PostScript file: %s\n",PSFName);

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);         
    fprintf(PSFd,"gsave\n");

    if (PMNC != 1)
        set_clip();

    ps_lwidth(PMLW);
    ps_ltyp(PMLT);

    for (i = 0; i <= n; ++i) {
        x = ps_2dx(AcX[i]);       
        y = ps_2dy(AcY[i]);       
        if (i == 0)
            fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
        else              
            fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
    }
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");

    err = 0;

PXFFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xopen       Open a new PostScript file.                                 */
/*              xopen = psfile;                                             */
/*              Return 0 if OK, -1 if not.                                  */

int xopen(void)
{
    register char *p;
    int err,fnd,n,len;

    err = -1;
    if (check_cmd(1))
        return(-1);

    printf1("Open new PostScript file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,2,1))     /* get parameters */
        goto XOPENFin;
  
    if (PSFFlg)  
        ps_close(0,0,0);
   
    if (alloc_acc(RLMaxDef))
        goto XOPENFin;

    len = fnd = 0;
    while (fgets(AcC,RLMaxDef,PMFd)) {

        len += strlen(AcC) + 1;

        if (fnd == 0 && !strncmp(AcC,"%%Creator: TDA",14))  
            fnd = 1;

        if (fnd == 1 && sscanf(AcC + 2,"BoundingBox: %d %d %d %d",&BBLX,&BBLY,&BBUX,&BBUY) == 4)  
            fnd = 2;

        if (fnd == 2 && sscanf(AcC,"%%#Parameter: %d %d %lg %lg %lg %lg %lg %lg %d %d %lg %lg %lg\n",
            &XOrg,&YOrg,&PXLen,&PYLen,&PA1[0],&PA1[1],&PA2[0],&PA2[1],&PSLog[0],&PSLog[1],
            &SCALXFac,&SCALYFac,&ROTFac) == 13) {
            fnd = 3;
        }
    }
    if (fnd != 3) {
        printf1("Can't read parameters from: %s\n",PMFdName);
        goto XOPENFin;
    }
    UXLen = PA2[0] - PA1[0];
    UYLen = PA2[1] - PA1[1];

    PSXLen = PtMM * PXLen;   /* phys. x axis length points */
    PSYLen = PtMM * PYLen;   /* phys. y axis length points */

    printf1("New PostScript file: %s\n\n",PMFdName);
    printf1("Size (width x height in mm): %g x %g\n",PXLen,PYLen);
    printf1("PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    XOrg,XOrg + (int)(PSXLen + 0.5),YOrg,YOrg + (int)(PSYLen + 0.5));   
    printf1("Origin (PostScript x and y coordinates): %d,%d\n",XOrg,YOrg);
    printf1("User coordinates X axis: %lg,%lg\n",PA1[0],PA2[0]);
    printf1("User coordinates Y axis: %lg,%lg\n",PA1[1],PA2[1]);
    printf1("BoundingBox: %d %d %d %d\n",BBLX,BBLY,BBUX,BBUY);       

    if (alloc_acd(len + 1))
        goto XOPENFin;

    PSONUM = 0;
    fseek(PMFd,0,0);
    p = AcD;
    while (fgets(p,RLMaxDef,PMFd)) {
        if (!strncmp(p,"showpage",8))  
            break;
        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            if (PSONUM < n)
                PSONUM = n;
        }
        p += strlen(p);
    }
    *p = '\0';

    strcpy(PSFName,PMFdName);
    err = -2;
    if (!(PSFd = fopen(PSFName,OPEN_WR)))    
        goto XOPENFin;

    fwrite(AcD,(int)(p - AcD),1,PSFd);   
    PSFFlg = 2;
    err = 0;

XOPENFin:
    if (err == -2)
        printf1("\nError: can't re-open: %s\n",PSFName);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xdelete     Delete plot objects from current PostScript file.           */
/*              xdelete = list of plot object numbers;                      */
/*              Return 0 if OK, -1 if not.                                  */

int xdelete(void)
{
    register int i;
    register char *p;
    int err,n,fnd,sflg,len,oflag;

    err = -1;
    oflag = 0;
    if (check_cmd(1))
        return(-1);

    printf1("Delete plot objects. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,7,1))     /* get parameters */
        goto XDELFin;

    if (xlog_check() != 0)
        goto XDELFin;

    printf1("Deleting plot object(s): %d",(int)PMRHSX[0]);
    for (i = 1; i < PMRHSN; ++i)
        printf1(",%d",(int)PMRHSX[i]);
    newline();

    if (alloc_acc(RLMaxDef))
        goto XDELFin;

    ps_close(0,0,1);    /* close current plot file */

    if (!(PMFd = fopen(PSFName,OPEN_RD))) {   
        printf1("Error: can't re-open %s.\n",PSFName);
        goto XDELFin;
    }
    oflag = 1;

    len = 0;
    while (fgets(AcC,RLMaxDef,PMFd))  
        len += strlen(AcC) + 1;

    if (alloc_acd(len + 1))
        goto XDELFin;

    PSONUM = 0;
    p = AcD;
    fseek(PMFd,0,0);
    sflg = 0;
    while (fgets(p,RLMaxDef,PMFd)) {

        if (!strncmp(p,"showpage",8))  
            break;

        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            fnd = -1;
            for (i = 0; i < PMRHSN; ++i) {
                if (n == (int)PMRHSX[i]) {
                    fnd = n;
                    sflg = 1;
                    PMRHSX[i] *= -1.0;
                    break;
                }
            }
            if (fnd < 0) {
                sflg = 0;
                if (PSONUM < n)
                    PSONUM = n;
            }
        }
        if (sflg == 0)  
            p += strlen(p);
    }
    *p = '\0';

    err = -2;
    if (!(PSFd = fopen(PSFName,OPEN_WR)))    
        goto XDELFin;

    fwrite(AcD,(int)(p - AcD),1,PSFd);   
    PSFFlg = 2;

    n = 0;
    for (i = 0; i < PMRHSN; ++i) {
        if (PMRHSX[i] > 0.0) {
            if (n == 0) 
                printf1("Ignored: %d",(int)PMRHSX[i]);
            else
                printf1(",%d",(int)PMRHSX[i]);
            n++;
        }
    }
    if (n)
        newline();
    err = 0;

XDELFin:
    if (err == -2)  
        printf1("\nError: can't re-open: %s\n",PSFName);
    if (oflag) {
        fclose(PMFd);
        PMFDef = 0;
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xdens       Plot of density, or frequency distribution.                 */
/*              This command always creates a new PostScript file.          */
/*                                                                          */
/*              xdens(                                                      */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  pxlen=...,      length of x axis, def. 120 mm           */
/*                  pylen=...,      length of y axis, def.  80 mm           */
/*                                                                          */
/*              ) = varlist;                                                */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xdens(void)
{
    register int i,j,k;
    int err,n,iflag,r,nmax;
    double tmp,xmin,xmax,ymin,ymax,d;

    err = -1;

    if (check_cmd(1))
        return(-1);

    printf1("XDens. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,4,1))     /* get parameters */
        goto XDENSFin;

    if (PMOPT != 2)
        PMOPT = 1;

    if (PMLT < 1 || PMLT > 9)
        PMLT = 1;

    if (PSFFlg)  
        ps_close(0,0,0);

    xmin = xmax = get_data(PMVIdx[0],0);

    iflag = 1;
    for (i = 0; i < NOC; ++i) {
        for (j = 0; j < PMNV; ++j) {
            tmp = get_data(PMVIdx[j],i);       
            if (iflag) {
                n = (int)tmp;
                if ((double)n != tmp)
                    iflag = 0;
            }
            xmin = dmin(xmin,tmp);
            xmax = dmax(xmax,tmp);
        }
    }
/*
printf1("iflag=%d minx=%lg %lg\n",iflag,xmin,xmax);
*/
    if (iflag) {            /* integer values, plot frequencies */

        printf1("Frequency plot of variable(s): %s",VName[PMVIdx[0]]);
        for (i = 1; i < PMNV; ++i)
            printf1(",%s",VName[PMVIdx[i]]);
        newline();

        if (PMLWFlg == 0)
            PMLW = 0.8;

        if (alloc_ack(NOC + 1))     /* value */
            goto XDENSFin;

        if (alloc_acn(NOC + 1))     /* freq */
            goto XDENSFin;

        if (alloc_aci(NOC + 1))     /* bf pointer */
            goto XDENSFin;

        nmax = 0;
        for (i = 0; i < PMNV; ++i) {    

            r = cfreq1(1,PMVIdx + i,NOC,AcK,AcN,AcI,-1,&tmp);

            if (r <= 0) {
                printf1("Error: can't create frequency distribution.\n");
                goto XDENSFin;
            }
            for (j = 1; j <= r; ++j) {
                k = AcI[j];
                nmax = imax(nmax,AcN[k]);
            }
        }
        ymin = 0.0;
        if (PMOPT == 1)
            ymax = (double)nmax / (double)NOC;
        else
            ymax = (double)nmax;

        if (x_cpsf(xmin - 1.0,xmax + 1.0,ymin,ymax,1))
            goto XDENSFin;

        d = 2.0 * UXLen * PMLW / PXLen;

        for (i = 0; i < PMNV; ++i) {    

            r = cfreq1(1,PMVIdx + i,NOC,AcK,AcN,AcI,-1,&tmp);

            if (r <= 0) {
                printf1("Error: can't create frequency distribution.\n");
                goto XDENSFin;
            }
            fprintf(PSFd,"\n%%#%d: %s [%s]\n",++PSONUM,CmdBuf,VName[PMVIdx[i]]);         
            fprintf(PSFd,"gsave\n");
            ps_lwidth(PMLW);
            ps_ltyp(PMLT);

            for (j = 1; j <= r; ++j) {
                k = AcI[j];
                tmp = (double)AcN[k];
                if (PMOPT == 1)
                    tmp /= (double)NOC;

                ps_2dplot((double)AcK[k] + d * (double)i,0.0,0);
                ps_2dplot((double)AcK[k] + d * (double)i,tmp,1);
                fprintf(PSFd,"stroke\n");
            }
            fprintf(PSFd,"grestore\n");
            if (++PMLT > 9)
                PMLT = 1;
        }
    }
    err = 0;

XDENSFin:
    p_clean();
    return(err);
}

