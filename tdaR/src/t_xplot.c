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
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_xplot.c                                                  */

int xplot(TDAContext *ctx);
int x_cpsf(TDAContext *ctx, double xmin,double xmax,double ymin,double ymax,int opt);
int alloc_scx(TDAContext *ctx, int n);
void x_psetup(TDAContext *ctx, double xa,double xb,double ya,double yb);
int x_getn(TDAContext *ctx, double x);
int x_getaxval(TDAContext *ctx, double xmin,double xmax,double *xa,double *xb,int opt);
void x_plaxis(TDAContext *ctx, int opt,int nx,double xa,double xb);
void x_scplot(TDAContext *ctx, int n,int ns,int opt,int clip,char *cmd);
int x_cchk(TDAContext *ctx, double xmin,double xmax,double ymin,double ymax);
int xcheck(TDAContext *ctx);
int xlog(TDAContext *ctx);
int xlog1(TDAContext *ctx);
int xlog_check(TDAContext *ctx);
void xlog_prn(TDAContext *ctx, int opt);
int xlogp(TDAContext *ctx);
int xconh(TDAContext *ctx);
int xconhp(TDAContext *ctx, int g,int lt,char *cmd);
int xreg(TDAContext *ctx);
int xplotf(TDAContext *ctx);
int xopen(TDAContext *ctx);
int xdelete(TDAContext *ctx);
int xdens(TDAContext *ctx);
int xfunc(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */



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

int xplot(TDAContext *ctx)
{
    register int i,j;
    int err,ix,iy,ig,n;

    ig = err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XPlot. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto XPLOTFin;

    ctx->PMOPT--;
    if (ctx->PMOPT != 1)
        ctx->PMOPT = 0;

    if (ctx->PMNV < 2 || ctx->PMNV > 3) {
        printf1(ctx, "Error: need two or three variables.\n");
        goto XPLOTFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    if (ctx->PMNV == 3)
        ig = ctx->PMVIdx[2];

    printf1(ctx, "Scatterplot with: %s, %s\n",ctx->VName[ix],ctx->VName[iy]);
    if (ig >= 0)
        printf1(ctx, "Grouping with: %s\n",ctx->VName[ig]);

    if (alloc_scx(ctx, ctx->NOC))
        goto XPLOTFin;

    ctx->SCNGRP = 1;                 /* number of groups */
    ctx->SCGRP[0] = 1;

    if (ig >= 0) {
        for (i = 0; i < ctx->NOC; ++i)  
            ctx->SCGVar[i] = (int)get_data(ctx, ig,i);       

        if (sorti(ctx, ctx->NOC,ctx->SCGVar,0))
            goto XPLOTFin;

        ctx->SCGRP[0] = ctx->SCGVar[0];
        for (i = 1; i < ctx->NOC; ++i) {
            if (ctx->SCGVar[i] != ctx->SCGVar[i - 1]) {
                if (ctx->SCNGRP < MaxG) 
                    ctx->SCGRP[ctx->SCNGRP++] = ctx->SCGVar[i];
                else {
                    printf1(ctx, "Warning: exceeded max number of groups.\n");
                    break;
                }
            }
        }
        printf1(ctx, "Number of groups: %d  [",ctx->SCNGRP);
        for (i = 0; i < ctx->SCNGRP; ++i)
            printf1(ctx, " %d",ctx->SCGRP[i]);
        printf1(ctx, "]\n");
    }
    ctx->SCXMin = ctx->SCXMax = get_data(ctx, ix,0);
    ctx->SCYMin = ctx->SCYMax = get_data(ctx, iy,0);

    for (i = 0; i < ctx->NOC; ++i) {

        ctx->SCXVar[i] = get_data(ctx, ix,i);       
        ctx->SCYVar[i] = get_data(ctx, iy,i);       
        ctx->SCXMin = dmin(ctx, ctx->SCXMin,ctx->SCXVar[i]); 
        ctx->SCXMax = dmax(ctx, ctx->SCXMax,ctx->SCXVar[i]); 
        ctx->SCYMin = dmin(ctx, ctx->SCYMin,ctx->SCYVar[i]); 
        ctx->SCYMax = dmax(ctx, ctx->SCYMax,ctx->SCYVar[i]); 

        if (ig >= 0)
            ctx->SCGVar[i] = (int)get_data(ctx, ig,i);       
        else
            ctx->SCGVar[i] = 1;
    }
    ctx->SCN = ctx->NOC;

    printf1(ctx, "X Min %lg  Max %lg\n",ctx->SCXMin,ctx->SCXMax);
    printf1(ctx, "Y Min %lg  Max %lg\n",ctx->SCYMin,ctx->SCYMax);

    if (ctx->PSFFlg == 2) {
        xlog_prn(ctx, 1);            /* using current PostScript file */

        if (ctx->SCXMin >= ctx->PA2[0] || ctx->SCXMax <= ctx->PA1[0] ||
            ctx->SCYMin >= ctx->PA2[1] || ctx->SCYMax <= ctx->PA1[1]) {  
            printf1(ctx, "Empty intersection with current coordinates.\n");
            err = 0;
            goto XPLOTFin;
        }
    }
    else {
        if (x_cpsf(ctx, ctx->SCXMin,ctx->SCXMax,ctx->SCYMin,ctx->SCYMax,0))
            goto XPLOTFin;
    }

    /* make scatterplots for all groups */
   
    if (ctx->PMOPT == 0) {               /* scatterplot */

        if (ctx->PMFSFlg == 0)
            ctx->PMFS = 1.3;                 /* default size of symbols */

        if (ctx->SCNGRP == 1) {
            if (ctx->PMS >= 1 && ctx->PMS <= 17)
                n = ctx->PMS;
            else
                n = 1;
        }
        else
            n = 0;

        i = 0;
        for (j = 0; j < ctx->SCNGRP; ++j) {
            if (n)
                x_scplot(ctx, ctx->SCGRP[j],n,ctx->PMOPT,1,ctx->CmdBuf);          
            else {
                x_scplot(ctx, ctx->SCGRP[j],ctx->SYMLIST[i],ctx->PMOPT,1,ctx->CmdBuf);          
                if (++i >= 17)
                    i = 0;
            }
        }
    }
    else {                      /* line plot */

        if (ctx->PMLT < 1 || ctx->PMLT > 9)
            ctx->PMLT = 1;

        for (j = 0; j < ctx->SCNGRP; ++j) {
            x_scplot(ctx, ctx->SCGRP[j],ctx->PMLT,ctx->PMOPT,1,ctx->CmdBuf);          
            if (++ctx->PMLT > 9)
                ctx->PMLT = 1;  
        }
    }
    err = 0;

XPLOTFin:
    if (err)
        alloc_scx(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  x_cpsf()    Create new PostScript file with standard coordinate system. */  
/*                                                                          */
/*  Return: 0 if OK, otherwise -1                                           */

int x_cpsf(TDAContext *ctx, double xmin,double xmax,double ymin,double ymax,int opt)
{
    int nx,ny;
    double xa,xb,ya,yb;

    if (x_cchk(ctx, xmin,xmax,ymin,ymax))
        return(-1);     
      
    nx = x_getaxval(ctx, xmin,xmax,&xa,&xb,opt);     /* get values for axes */
    ny = x_getaxval(ctx, ymin,ymax,&ya,&yb,opt);

    /* setup new Postscipt file */

    newline(ctx);
    strcpy(ctx->PSFName,"xplot.ps");

    if (!(ctx->PSFd = fopen(ctx->PSFName,OPEN_WR))) {
        printf1(ctx, "Error: can't create: %s\n",ctx->PSFName);
        return(-1);    
    }
    ps_init(ctx);          /* write header to output file */
    ctx->PSFFlg = 1;
    printf1(ctx, "New PostScript file: %s\n",ctx->PSFName);

    x_psetup(ctx, xa,xb,ya,yb);      /* set up coordinate system */

    x_plaxis(ctx, 0,nx,xa,xb);
    x_plaxis(ctx, 1,ny,ya,yb);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  alloc_scx(n)    If n > 0 allocate SCXVar,SCYVar,SCGVar, otherwise free. */  
/*                                                                          */
/*  Return: 0 if OK, otherwise insufficient memory.                         */

int alloc_scx(TDAContext *ctx, int n)
{
    if (ctx->SCXVarN > 0) {
        free((char *)ctx->SCXVar);
        memrq(ctx, -ctx->SCXVarN,sizeof(double));
        ctx->SCXVarN = 0;
    }
    if (ctx->SCYVarN > 0) {
        free((char *)ctx->SCYVar);
        memrq(ctx, -ctx->SCYVarN,sizeof(double));
        ctx->SCYVarN = 0;
    }
    if (ctx->SCGVarN > 0) {
        free((char *)ctx->SCGVar);
        memrq(ctx, -ctx->SCGVarN,sizeof(int));
        ctx->SCGVarN = 0;
    }
    ctx->SCN = 0;

    if (n > 0) {
        if (!(ctx->SCXVar = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->SCXVarN = n;
        if (!(ctx->SCYVar = (double *)calloc((size_t)(n),sizeof(double)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(double));
        ctx->SCYVarN = n;
        if (!(ctx->SCGVar = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->SCGVarN = n;
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  x_psetup()  Set up coordinate system.                                   */

void x_psetup(TDAContext *ctx, double xa,double xb,double ya,double yb)
{
    ctx->PSLog[0] = ctx->PSLog[1] = 0;  

    ctx->XOrg = 150;         /* origin of PostScript figure */
    ctx->YOrg = 460;  
    ctx->SCALXFac = 1.0;     /* PostScript x scaling factor */
    ctx->SCALYFac = 1.0;     /* PostScript y scaling factor */
    ctx->ROTFac = 0.0;       /* PostScript rotation factor */

    ctx->PXLen = ctx->PMXLen;     /* phys length of X axis in mm */
    ctx->PYLen = ctx->PMYLen;     /* phys length of Y axis in mm */

    ctx->PA1[0] = xa;
    ctx->PA2[0] = xb;
    ctx->PA1[1] = ya;
    ctx->PA2[1] = yb;
         
    ctx->UXLen = ctx->PA2[0] - ctx->PA1[0];
    ctx->UYLen = ctx->PA2[1] - ctx->PA1[1];

    ctx->PSXLen = ctx->PtMM * ctx->PXLen;   /* phys. x axis length points */
    ctx->PSYLen = ctx->PtMM * ctx->PYLen;   /* phys. y axis length points */

    printf1(ctx, "Size (width x height in mm): %g x %g\n",ctx->PXLen,ctx->PYLen);
    printf1(ctx, "PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    ctx->XOrg,ctx->XOrg + (int)(ctx->PSXLen + 0.5),ctx->YOrg,ctx->YOrg + (int)(ctx->PSYLen + 0.5));   
    printf1(ctx, "Origin (PostScript x and y coordinates): %d,%d\n",ctx->XOrg,ctx->YOrg);

    if (ctx->PSFFlg == 1) {      /* init bounding box parameters */
        ctx->BBLX = ctx->XOrg;
        ctx->BBLY = ctx->YOrg;
        ctx->BBUX = ctx->BBLX + (int)(ctx->PSXLen + 0.5);
        ctx->BBUY = ctx->BBLY + (int)(ctx->PSYLen + 0.5);
    }
    ps_inis(ctx, 1);

    printf1(ctx, "User coordinates X axis: %lg,%lg\n",ctx->PA1[0],ctx->PA2[0]);
    printf1(ctx, "User coordinates Y axis: %lg,%lg\n",ctx->PA1[1],ctx->PA2[1]);

    ctx->PSFFlg = 2;     /* set to 2 if setup successful */
}

/*--------------------------------------------------------------------------*/
/*  x_getn()                                                                */
/*                                                                          */

int x_getn(TDAContext *ctx, double x) 
{
    (void)ctx;        /* unused: the signature is shared */
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

int x_getaxval(TDAContext *ctx, double xmin,double xmax,double *xa,double *xb,int opt)
{
    int n,nn;
    double d;

    d = (xmax - xmin) / 5.0;
    if (d < ctx->EPSI1)
        d = ctx->EPSI1;

    /* printf1("xmin=%lg xmax=%lg d=%lg\n",xmin,xmax,d); */
    n = x_getn(ctx, d);
    d = pow(10.0,(double)n);

    /* printf1("nach n=%d d=%lg\n",n,d); */  

    if (xmin >= 0.0 && xmin <= ctx->EPSI1)
        nn = 0;
    else 
        nn = (int)floor(xmin / d);
    *xa = (double)nn * d;

    if (opt == 0) {
        if (fabs(*xa - xmin) <= ctx->EPSI1) 
            *xa -= d;
    }
    if (xmax <= 0.0 && xmax >= -ctx->EPSI1)
        nn = 0;
    else  
        nn = (int)ceil(xmax / d);
    *xb = (double)nn * d;

    if (opt == 0) {
        if (fabs(*xb - xmax) <= ctx->EPSI1) 
            *xb += d;
    }
    return(n);
}

/*--------------------------------------------------------------------------*/
/*  x_plaxis(opt,nx,xa,xb)      plot axis for coordinate system.            */
/*                              opt 0 : x axis                              */
/*                              opt 1 : y axis                              */

void x_plaxis(TDAContext *ctx, int opt,int nx,double xa,double xb)         
{
    ctx->PMIC = 0.0;
    ctx->PMSC = pow(10.0,(double)nx);
    if ((xb - xa) / ctx->PMSC > 20) {    
        ctx->PMSC = pow(10.0,(double)(nx+1));
        ctx->PMIC = 10.0;
    }
    pl_axis(ctx, opt,1,1,1,0.2);
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

void x_scplot(TDAContext *ctx, int n,int ns,int opt,int clip,char *cmd)
{
    register int i;
    int first;
    double x,y;

    fprintf(ctx->PSFd,"\n%%#%d: %s [group %d]\n",++ctx->PSONUM,cmd,n);
    fprintf(ctx->PSFd,"gsave\n");
   
    if (clip)
        set_clip(ctx);

    ps_lwidth(ctx, ctx->PMLW);
    if (opt)
        ps_ltyp(ctx, ns);
    else
        ps_ltyp(ctx, 1);

    if (opt == 0) {
        fprintf(ctx->PSFd,"gsave\n");
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);
    }
    first = 1;
    for (i = 0; i < ctx->SCN; ++i) {
        if (ctx->SCGVar[i] == n) {
            x = ps_2dx(ctx, ctx->SCXVar[i]);       
            y = ps_2dy(ctx, ctx->SCYVar[i]);       
            if (opt == 0)
                ps_sym(ctx, ns,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            else if (first) {
                fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else              
                fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
        }
    }
    if (opt == 1)
        fprintf(ctx->PSFd,"stroke\n");
    else
        fprintf(ctx->PSFd,"grestore\n");
    fprintf(ctx->PSFd,"grestore\n");
}

/*--------------------------------------------------------------------------*/
/*  x_cchk      Check for possibility of coordinate system.                 */
/*              Return 0 if OK, -1 if not.                                  */

int x_cchk(TDAContext *ctx, double xmin,double xmax,double ymin,double ymax)
{
    if (fabs(xmin - xmax) < ctx->EPSI1 || fabs(ymin - ymax) < ctx->EPSI1) {
        printf1(ctx, "\nError: can't create a coordinate system.\n");
        return(-1);    
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  xcheck      Check if basic data for xplot are available.                */
/*              Return 0 if OK, -1 if not.                                  */

int xcheck(TDAContext *ctx)
{
    if (ctx->SCN > 0)  
        return(0);
    printf1(ctx, "Error: no current scatter plot.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  xlog        Print list of plot objects.                                 */
/*              xlog [=psfile];                                             */
/*              Return 0 if OK, -1 if not.                                  */

int xlog(TDAContext *ctx)
{
    int err, oflag;

    err = -1;
    oflag = 0;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "List of plot objects. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,2,0))     /* get parameters */
        goto XLOGFin;

    if (ctx->PMFDef == 0) {
        if (xlog_check(ctx) != 0)
            goto XLOGFin;
        xlog_prn(ctx, 1);

        fflush(ctx->PSFd);
        if (!(ctx->PMFd = fopen(ctx->PSFName,OPEN_RD))) {   
            printf1(ctx, "Error: can't re-open this file.\n");
            goto XLOGFin;
        }
        oflag = 1;
    }
    else  
        printf1(ctx, "Using PostScript file: %s\n",ctx->PMFdName);

    if (alloc_acc(ctx, RLMaxDef))
        goto XLOGFin;

    if (xlogp(ctx) == 0)            /* print list */
        printf1(ctx, "No plot objects.\n");
    err = 0;

XLOGFin:
    if (oflag) {
        fclose(ctx->PMFd);
        ctx->PMFDef = 0;
    }
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xlog1       Print current ps file and list of plot objects.             */
/*              Return 0 if OK, -1 if not.                                  */

int xlog1(TDAContext *ctx)
{
    int err;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->PSFFlg != 2) {
        printf1(ctx, "No current PostScript file.\n");
        return(0);
    }
    xlog_prn(ctx, 1);
    printf1(ctx, "\nList of plot objects. ");
    fflush(ctx->PSFd);
    if (!(ctx->PMFd = fopen(ctx->PSFName,OPEN_RD))) {   
        printf1(ctx, "Error: can't re-open this file.\n");
        goto XLOG1Fin;
    }
    if (alloc_acc(ctx, RLMaxDef))
        goto XLOG1Fin;

    if (xlogp(ctx) == 0)            /* print list */
        printf1(ctx, "No plot objects.\n");
    fclose(ctx->PMFd);
    ctx->PMFDef = 0;
    err = 0;

XLOG1Fin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xlog_check  Check for current plot file. Return 0 if OK, -1 if error.   */

int xlog_check(TDAContext *ctx)
{
    if (ctx->PSFFlg != 2) {
        printf1(ctx, "Error: no current plot file.\n");
        return(-1);   
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  xlog_prn    Print current plot file and, optionally, coordinates.       */

void xlog_prn(TDAContext *ctx, int opt)
{
    printf1(ctx, "Using current plot file: %s\n",ctx->PSFName);
    if (opt)
        printf1(ctx, "X axis: %lg,%lg  Y axis: %lg,%lg\n",ctx->PA1[0],ctx->PA2[0],ctx->PA1[1],ctx->PA2[1]);
}

/*--------------------------------------------------------------------------*/
/*  xlogp       Print list of plot objects.                                 */
/*              Return number of plot objects.                              */

int xlogp(TDAContext *ctx)
{
    int n,m,first;
    register char *p;

    m = 0; first = 1;
    fseek(ctx->PMFd,0,0);
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {
        p = ctx->AcC;
        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            p = skip_int(ctx, p + 2);
            if (first) {
                newline(ctx);
                first = 0;
            }
            printf1(ctx, "%3d%s",n,p);
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

int xconh(TDAContext *ctx)
{
    register int i,j;
    int err,n;

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XConh. Current memory: %d bytes.\n",ctx->MemReq);
    if (xcheck(ctx))
        goto XCONHFin;

    if (parm(ctx, ctx->CmdBuf + 5,3,0))     /* get parameters */
        goto XCONHFin;

    printf1(ctx, "Adding %d convex hull(s).\n",ctx->SCNGRP);

    j = ctx->PMLT;
    if (j < 1 || j > 9)
        j = 1;

    for (i = 0; i < ctx->SCNGRP; ++i) {
        n = xconhp(ctx, ctx->SCGRP[i],j,ctx->CmdBuf);
        if (n)  
            printf1(ctx, "Error in group %d\n",ctx->SCGRP[i]);
        if (ctx->PMLTFlg == 0) {
            if (++j > 9)
                j = 1;
        }   
    }
    err = 0;
   
XCONHFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xconhp(g,lt,cmd)   Plot convex hull for group g with line type lt.      */
/*                     Use cmd for plot object description.                 */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int xconhp(TDAContext *ctx, int g,int lt,char *cmd)
{
    register int i,j,ik,l;
    int err,nn,n;
    double x,y;

    err = -1;

    if (alloc_acxf(ctx, ctx->SCN + 1))
        goto XCONHPFin;
    if (alloc_acyf(ctx, ctx->SCN + 1))
        goto XCONHPFin;
    if (alloc_acn(ctx, ctx->SCN + 1))
        goto XCONHPFin;
    if (alloc_aci(ctx, ctx->SCN + 1))
        goto XCONHPFin;
    if (alloc_acj(ctx, ctx->SCN + 1))
        goto XCONHPFin;

    nn = 0;
    for (i = 0; i < ctx->SCN; ++i) {
        if (ctx->SCGVar[i] == g) {
            nn++;
            ctx->AcXF[nn] = (float)ctx->SCXVar[i];            
            ctx->AcYF[nn] = (float)ctx->SCYVar[i];            
            ctx->AcN[nn] = nn;
        }
    }
    if (nn < 1)
        goto XCONHPFin;

    n = g_chull(ctx, ctx->AcXF,ctx->AcYF,nn,ctx->AcN,ctx->AcJ,ctx->AcI);   
    if (n <= 0) {
        if (n < 0)
            p_err(ctx, -2,1);
        goto XCONHPFin;
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s [group %d]\n",++ctx->PSONUM,cmd,g);
    fprintf(ctx->PSFd,"gsave\n");
    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, lt);

    ik = ctx->AcI[1];
    l = 0;
    for (i = 1; i <= n; ++i) {
        j = ctx->AcJ[ik];
        x = (double)ctx->AcXF[j];
        y = (double)ctx->AcYF[j];
        ps_2dplot(ctx, x,y,l); 
        l = 1;
        ik = ctx->AcI[ik];
    }
    ik = ctx->AcI[1];
    j = ctx->AcJ[ik];
    ps_2dplot(ctx, (double)ctx->AcXF[j],(double)ctx->AcYF[j],1);     /* close path */

    /********************************************
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");       
        ps_fill(ctx, PMGS);
        fprintf(PSFd,"grestore\n");   
    }
    *********************************/
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
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

int xreg(TDAContext *ctx)
{
    register int i;
    int err,typ,n,ig,g,first,ra,re,r,nx;
    double x,y,xmin,xmax,rne,rnl,xx[4];

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XREG. Current memory: %d bytes.\n",ctx->MemReq);
    if (xcheck(ctx))
        goto XREGFin;

    if (parm(ctx, ctx->CmdBuf + 4,3,0))     /* get parameters */
        goto XREGFin;

    switch (ctx->PMRHSI) {
        case 2:         printf1(ctx, "Least squares regression.\n");
                        typ = 1;
                        nx = 3;
                        break;
        case 3:         printf1(ctx, "L1 norm regression.\n");
                        typ = 2;
                        nx = 2;
                        break;
        default:        printf1(ctx, "Lowess.\n");
                        typ = 3;
                        nx = 1;
                        break;
    }
    if (ctx->PMSIG <= 0.0 || ctx->PMSIG >= 1.0)
        ctx->PMSIG = 0.5;

    if (alloc_acx(ctx, nx * ctx->SCN + 1))
        goto XREGFin;

    if (typ > 1) {
        if (alloc_acy(ctx, ctx->SCN + 1))
            goto XREGFin;
    }
    if (typ == 3) {
        if (alloc_acu(ctx, ctx->SCN + 1))
            goto XREGFin;
        if (alloc_acv(ctx, ctx->SCN + 1))
            goto XREGFin;
        if (alloc_acw(ctx, ctx->SCN + 1))
            goto XREGFin;
    }
    if (ctx->PMLT < 1 || ctx->PMLT > 9)
        ctx->PMLT = 1;

    for (ig = 0; ig < ctx->SCNGRP; ++ig) {
        xmin = xmax = 0.0;
        n = 0;
        g = ctx->SCGRP[ig];
        first = 1;
        for (i = 0; i < ctx->SCN; ++i) {
            if (ctx->SCGVar[i] == g) {
                x = ctx->SCXVar[i];            
                y = ctx->SCYVar[i];            
                if (first) {
                    xmax = xmin = x;
                    first = 0;
                }   
                else {
                    xmin = dmin(ctx, xmin,x);
                    xmax = dmax(ctx, xmax,x);
                }
                if (typ == 1) {
                    ctx->AcX[n * nx + 1] = 1.0;
                    ctx->AcX[n * nx + 2] = x;
                    ctx->AcX[n * nx + 3] = y;
                }
                else if (typ == 2) {
                    ctx->AcX[n * nx + 1] = 1.0;
                    ctx->AcX[n * nx + 2] = x;
                    ctx->AcY[n + 1] = y;
                }
                else if (typ == 3) {
                    ctx->AcX[n + 1] = x;
                    ctx->AcY[n + 1] = y;
                    ctx->AcW[n + 1] = 1.0;
                }
                n++;
            }
        }
        if (n < 2 || xmax < xmin + ctx->EPSI1)
            continue;

        printf1(ctx, "Group %3d : ",g);

        switch (typ) {
            case 1:         /* least squares regression */

                r = lsei(ctx, ctx->AcX,0,n,0,2,xx,0,&rne,&rnl,&ra,&re);
                if (r)
                    printf1(ctx, "error.\n");
                else  
                    printf1(ctx, "rank=%d  alpha=%lg  beta=%lg  norm of residuals=%lg\n",ra,xx[1],xx[2],rnl);
                break;

            case 2:         /* l1 norm regression */

                r = l1regf(ctx, n,2,ctx->AcX,ctx->AcY,xx,&rnl,&ra);
                if (r)
                    printf1(ctx, "error.\n");
                else  
                    printf1(ctx, "rank=%d  alpha=%lg  beta=%lg  sum of residuals=%lg\n",ra,xx[1],xx[2],rnl);
                break;

            case 3:         /* lowess regression */

                if (sortd2(ctx, n,ctx->AcX + 1,ctx->AcY + 1))     /* sort simultaneously according to x */
                    r = -1;            
                else  
                    r = lowess(ctx, ctx->AcX,ctx->AcY,n,ctx->PMSIG,2,0.0,ctx->AcU,ctx->AcW,ctx->AcV);
                if (r)
                    printf1(ctx, "error.\n");
                else
                    printf1(ctx, "sigma=%lg\n",ctx->PMSIG);
                break;

            default: break;
        }
        if (r)
            continue;
   
        fprintf(ctx->PSFd,"\n%%#%d: %s [group %d]\n",++ctx->PSONUM,ctx->CmdBuf,g);
        fprintf(ctx->PSFd,"gsave\n");
        ps_lwidth(ctx, ctx->PMLW);
        ps_ltyp(ctx, ctx->PMLT);

        if (typ == 3) {
            ps_2dplot(ctx, ctx->AcX[1],ctx->AcU[1],0);
            for (i = 2; i <= n; ++i)  
                ps_2dplot(ctx, ctx->AcX[i],ctx->AcU[i],1);
        }
        else {
            x = xx[1] + xx[2] * xmin;
            ps_2dplot(ctx, xmin,x,0);
            x = xx[1] + xx[2] * xmax;
            ps_2dplot(ctx, xmax,x,1);
        }
        fprintf(ctx->PSFd,"stroke\n");
        fprintf(ctx->PSFd,"grestore\n");

        if (ctx->PMLTFlg == 0) {
            if (++ctx->PMLT > 9)
                ctx->PMLT = 1;
        }
    }
    err = 0;
   
XREGFin:
    p_clean(ctx);
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

int xplotf(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,nn,cmax,cn,ng,cx,cy;
    double tmp,tmp1 = 0.0;              
    register char *p;

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XPLOTF. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,2,1))     /* get parameters */
        goto XSCFFin;

    ctx->PMOPT--;
    if (ctx->PMOPT != 1)
        ctx->PMOPT = 0;

    cmax = 1;
    if (ctx->PMOPT == 0)
        printf1(ctx, "Scatterplot. ");
    else if (ctx->PMOPT == 1)
        printf1(ctx, "Line plot. ");
    printf1(ctx, "File: %s  ",ctx->PMFdName);
    printf1(ctx, "Columns:");
    if (ctx->PMNCN > 0) {
        for (i = 0; i < ctx->PMNCN; ++i) {
            printf1(ctx, " %d",ctx->PMCN[i]);
            if (cmax < ctx->PMCN[i])
                cmax = ctx->PMCN[i];
        }
        newline(ctx);
        if (ctx->PMNCN < 2) {
            printf1(ctx, "Error: need at least two columns.\n");
            goto XSCFFin;
        }
        cn = ctx->PMNCN;
        cx = ctx->PMCN[0];
        cy = ctx->PMCN[1];
    }
    else {
        printf1(ctx, " 1 2\n");
        cn = cmax = 2;
        cx = 1;
        cy = 2;
    }
    if (ctx->PMGN > 0) {
        printf1(ctx, "Group definition by column %d\n",ctx->PMGN);
        if (cn > 2) {
            printf1(ctx, "Error: cn must provide exactly two columns.\n");
            goto XSCFFin;
        }
        if (cmax < ctx->PMGN)
            cmax = ctx->PMGN;
    }

    if (alloc_acc(ctx, RLMaxDef + 1))
        goto XSCFFin;

    i = ng = n = 0;
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {
        i++;
        if (check_drec(ctx, ctx->AcC)) {      /* check for data records */
            n++;

            if (ctx->PMGN > 0) {
                p = ctx->AcC;
                for (j = 1; j <= ctx->PMGN; ++j) {
                    p = skip_sep(ctx, p);                    /* skip separator */
                    if (sscanf(p,"%lg",&tmp) != 1) {
                        printf1(ctx, "Error: can't read c%d in record %d.\n",j,i);
                        goto XSCFFin;
                    }
                    if (j == ctx->PMGN) {
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
                    p = skip_dval(ctx, p);
                }
            }
        }
    }
    if (n < 1) {
        printf1(ctx, "Error: file contains no data records.\n");
        goto XSCFFin;
    }
    if (ctx->PMGN == 0) {
        ng = cn - 1;
        nn = ng * n;
    }
    else
        nn = n;

    if (alloc_scx(ctx, nn))
        goto XSCFFin;

    if (alloc_actmp(ctx, cmax + 1))
        goto XSCFFin;

    fseek(ctx->PMFd,0,0);
    k = i = 0;
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {
        i++;
        if (check_drec(ctx, ctx->AcC)) {      /* check for data records */
            p = ctx->AcC;
            for (j = 1; j <= cmax; ++j) {
                p = skip_sep(ctx, p);                    /* skip separator */
                if (sscanf(p,"%lg",&tmp) != 1) {
                    printf1(ctx, "Error: can't read c%d in record %d.\n",j,i);
                    goto XSCFFin;
                }
                ctx->AcTmp[j] = tmp;
                p = skip_dval(ctx, p);
            }
            if (cn == 2) {
                ctx->SCXVar[k] = ctx->AcTmp[cx];
                ctx->SCYVar[k] = ctx->AcTmp[cy];
                if (ctx->PMGN > 0)
                    ctx->SCGVar[k] = (int)ctx->AcTmp[ctx->PMGN];
                else
                    ctx->SCGVar[k] = 1;            
                k++;
            }
            else {
                for (j = 1; j < cn; ++j) {
                    ctx->SCXVar[(j - 1) * n + k] = ctx->AcTmp[cx];
                    ctx->SCYVar[(j - 1) * n + k] = ctx->AcTmp[ctx->PMCN[j]];
                    ctx->SCGVar[(j - 1) * n + k] = j;         
                }
                k++;
            }
        }
    }
    alloc_acc(ctx, 0);

    ctx->SCNGRP = 1;
    k = ctx->SCGVar[0];
    ctx->SCGRP[0] = k;

    for (i = 1; i < nn; ++i) {
        if (ctx->SCGVar[i] != k) {
            if (ctx->SCNGRP >= MaxG) {
                printf1(ctx, "Error: exceeded max number of groups.\n");
                goto XSCFFin;
            }
            k = ctx->SCGVar[i];          
            ctx->SCGRP[ctx->SCNGRP++] = k;
        }
    }
    ctx->SCN = nn;

    printf1(ctx, "Number of groups: %d  [",ctx->SCNGRP);
    for (i = 0; i < ctx->SCNGRP; ++i)
        printf1(ctx, " %d",ctx->SCGRP[i]);
    printf1(ctx, "]\n");
       
    ctx->SCXMin = ctx->SCXMax = ctx->SCXVar[0];      
    ctx->SCYMin = ctx->SCYMax = ctx->SCYVar[0];     

    for (i = 1; i < ctx->SCN; ++i) {
        ctx->SCXMin = dmin(ctx, ctx->SCXMin,ctx->SCXVar[i]); 
        ctx->SCXMax = dmax(ctx, ctx->SCXMax,ctx->SCXVar[i]); 
        ctx->SCYMin = dmin(ctx, ctx->SCYMin,ctx->SCYVar[i]); 
        ctx->SCYMax = dmax(ctx, ctx->SCYMax,ctx->SCYVar[i]); 
    }
    printf1(ctx, "X Min %lg  Max %lg\n",ctx->SCXMin,ctx->SCXMax);
    printf1(ctx, "Y Min %lg  Max %lg\n",ctx->SCYMin,ctx->SCYMax);

    if (ctx->PSFFlg == 2) {
        xlog_prn(ctx, 1);            /* using current PostScript file */

        if (ctx->SCXMin >= ctx->PA2[0] || ctx->SCXMax <= ctx->PA1[0] ||
            ctx->SCYMin >= ctx->PA2[1] || ctx->SCYMax <= ctx->PA1[1]) {  
            printf1(ctx, "Empty intersection with current coordinates.\n");
            err = 0;
            goto XSCFFin;
        }
    }
                /* new PostScript file */
    else {
        if (x_cpsf(ctx, ctx->SCXMin,ctx->SCXMax,ctx->SCYMin,ctx->SCYMax,ctx->PMOPT))
            goto XSCFFin;
    }

    /* make scatterplots for all groups */

    if (ctx->PMOPT == 0) {               /* scatterplot */

        if (ctx->PMFSFlg == 0)
            ctx->PMFS = 1.3;                 /* default size of symbols */
        if (ctx->SCNGRP == 1) {
            if (ctx->PMS >= 1 && ctx->PMS <= 17)
                n = ctx->PMS;
            else
                n = 1;
        }
        else
            n = 0;

        i = 0;
        for (j = 0; j < ctx->SCNGRP; ++j) {
            if (n)
                x_scplot(ctx, ctx->SCGRP[j],n,ctx->PMOPT,1,ctx->CmdBuf);          
            else {
                x_scplot(ctx, ctx->SCGRP[j],ctx->SYMLIST[i],ctx->PMOPT,1,ctx->CmdBuf);          
                if (++i >= 17)
                    i = 0;
            }
        }
    }
    else {                      /* line plot */

        if (ctx->PMLT < 1 || ctx->PMLT > 9)
            ctx->PMLT = 1;

        for (j = 0; j < ctx->SCNGRP; ++j) {
            x_scplot(ctx, ctx->SCGRP[j],ctx->PMLT,ctx->PMOPT,1,ctx->CmdBuf);          
            if (++ctx->PMLT > 9)
                ctx->PMLT = 1;  
        }
    }
    err = 0;

XSCFFin:
    if (err)
        alloc_scx(ctx, 0);
    p_clean(ctx);
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

int xfunc(TDAContext *ctx)
{
    register int i;
    int err,n,r,nflg;    
    double x,y,xmin,xmax,ymin = 0.0,ymax = 0.0;

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XF. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 2,9,1))    /* get parameters */
        goto PXFFin;

    if (ctx->PMLT < 1 || ctx->PMLT > 9)
        ctx->PMLT = 1;

    if (ctx->FNFlg == 0) {
        p_err(ctx, -1,1);
        goto PXFFin;
    }
    if (ctx->FNArgN != 1) {
        printf1(ctx, "Error: there must be exactly one function argument.\n");
        goto PXFFin;
    }
    if (ctx->PMFTYP5) {
        p_err(ctx, -40,1);
        goto PXFFin;
    }
    if (ctx->FVFlg)
        printf1(ctx, "Sum over %d data matrix cases.\n",ctx->NOC);

    if (ctx->PSFFlg == 2) {
        nflg = 0;
        ctx->PMRXA = (float)(ctx->PA1[0]);
        ctx->PMRXB = (float)(ctx->PA2[0]);
        ctx->PMRXD = (float)((double)((ctx->PMRXB - ctx->PMRXA)) / 200.0);
        n = 200;
    }
    else {
        nflg = 1;
        if (ctx->PMRXFlg == 0) {
            ctx->PMRXA = 1.0;   
            ctx->PMRXB = 10.0;   
            ctx->PMRXD = (float)((double)((ctx->PMRXB - ctx->PMRXA)) / 200.0);
            n = 200;
        }
        else
            n = (int)((double)(((((double)(ctx->PMRXB) - (double)(ctx->PMRXA) + (double)(ctx->EPSI1))))) / (double)(ctx->PMRXD));
    }
    xmin = (double)(ctx->PMRXA);
    xmax = (double)(ctx->PMRXB);

    if (alloc_acx(ctx, n + 1))
        goto PXFFin;
    if (alloc_acy(ctx, n + 1))
        goto PXFFin;

    ctx->NINTMUsed = -1;

    for (i = 0; i <= n; ++i) {

        ctx->FNArgVal[0] = ctx->AcX[i] = (double)(ctx->PMRXA) + (double)(i) * (double)(ctx->PMRXD);

        r = get_flval(ctx, &y,ctx->FNArgN,ctx->FNArgVal,0,0,ctx->AcU,ctx->AcV,ctx->AcV);

        if (r) {     /* r from v_eval1(ctx) */  

            printf1(ctx, "Can't evaluate function for x=%lg.\n",ctx->AcX[i]);
            prn_emsg2(ctx, r);
            goto PXFFin;
        }
        ctx->AcY[i] = y;
        if (i == 0)
            ymin = ymax = y;
        else {
            ymin = dmin(ctx, ymin,y);
            ymax = dmax(ctx, ymax,y);
        }
    }
    if (ctx->NINTMUsed >= 0)  
        ni_info(ctx);

    printf1(ctx, "\nX Min %lg  Max %lg\n",xmin,xmax);
    printf1(ctx, "Y Min %lg  Max %lg\n",ymin,ymax);
   
    if (nflg) {         /* create new PostScript file */

        if (x_cpsf(ctx, xmin,xmax,ymin,ymax,1))
            goto PXFFin;
    }
    else
        printf1(ctx, "\nUsing current PostScript file: %s\n",ctx->PSFName);

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);         
    fprintf(ctx->PSFd,"gsave\n");

    if (ctx->PMNC != 1)
        set_clip(ctx);

    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, ctx->PMLT);

    for (i = 0; i <= n; ++i) {
        x = ps_2dx(ctx, ctx->AcX[i]);       
        y = ps_2dy(ctx, ctx->AcY[i]);       
        if (i == 0)
            fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
        else              
            fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
    }
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");

    err = 0;

PXFFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xopen       Open a new PostScript file.                                 */
/*              xopen = psfile;                                             */
/*              Return 0 if OK, -1 if not.                                  */

int xopen(TDAContext *ctx)
{
    register char *p;
    int err,fnd,n,len;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Open new PostScript file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,2,1))     /* get parameters */
        goto XOPENFin;
  
    if (ctx->PSFFlg)  
        ps_close(ctx, 0,0,0);
   
    if (alloc_acc(ctx, RLMaxDef))
        goto XOPENFin;

    len = fnd = 0;
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd)) {

        len += (int)strlen(ctx->AcC) + 1;

        if (fnd == 0 && !strncmp(ctx->AcC,"%%Creator: TDA",14))  
            fnd = 1;

        if (fnd == 1 && sscanf(ctx->AcC + 2,"BoundingBox: %d %d %d %d",&ctx->BBLX,&ctx->BBLY,&ctx->BBUX,&ctx->BBUY) == 4)  
            fnd = 2;

        if (fnd == 2 && sscanf(ctx->AcC,"%%#Parameter: %d %d %lg %lg %lg %lg %lg %lg %d %d %lg %lg %lg\n",
            &ctx->XOrg,&ctx->YOrg,&ctx->PXLen,&ctx->PYLen,&ctx->PA1[0],&ctx->PA1[1],&ctx->PA2[0],&ctx->PA2[1],&ctx->PSLog[0],&ctx->PSLog[1],
            &ctx->SCALXFac,&ctx->SCALYFac,&ctx->ROTFac) == 13) {
            fnd = 3;
        }
    }
    if (fnd != 3) {
        printf1(ctx, "Can't read parameters from: %s\n",ctx->PMFdName);
        goto XOPENFin;
    }
    ctx->UXLen = ctx->PA2[0] - ctx->PA1[0];
    ctx->UYLen = ctx->PA2[1] - ctx->PA1[1];

    ctx->PSXLen = ctx->PtMM * ctx->PXLen;   /* phys. x axis length points */
    ctx->PSYLen = ctx->PtMM * ctx->PYLen;   /* phys. y axis length points */

    printf1(ctx, "New PostScript file: %s\n\n",ctx->PMFdName);
    printf1(ctx, "Size (width x height in mm): %g x %g\n",ctx->PXLen,ctx->PYLen);
    printf1(ctx, "PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    ctx->XOrg,ctx->XOrg + (int)(ctx->PSXLen + 0.5),ctx->YOrg,ctx->YOrg + (int)(ctx->PSYLen + 0.5));   
    printf1(ctx, "Origin (PostScript x and y coordinates): %d,%d\n",ctx->XOrg,ctx->YOrg);
    printf1(ctx, "User coordinates X axis: %lg,%lg\n",ctx->PA1[0],ctx->PA2[0]);
    printf1(ctx, "User coordinates Y axis: %lg,%lg\n",ctx->PA1[1],ctx->PA2[1]);
    printf1(ctx, "BoundingBox: %d %d %d %d\n",ctx->BBLX,ctx->BBLY,ctx->BBUX,ctx->BBUY);       

    if (alloc_acd(ctx, len + 1))
        goto XOPENFin;

    ctx->PSONUM = 0;
    fseek(ctx->PMFd,0,0);
    p = ctx->AcD;
    while (fgets(p,RLMaxDef,ctx->PMFd)) {
        if (!strncmp(p,"showpage",8))  
            break;
        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            if (ctx->PSONUM < n)
                ctx->PSONUM = n;
        }
        p += strlen(p);
    }
    *p = '\0';

    strcpy(ctx->PSFName,ctx->PMFdName);
    err = -2;
    if (!(ctx->PSFd = fopen(ctx->PSFName,OPEN_WR)))    
        goto XOPENFin;

    fwrite(ctx->AcD,(size_t)((int)(p - ctx->AcD)),1,ctx->PSFd);   
    ctx->PSFFlg = 2;
    err = 0;

XOPENFin:
    if (err == -2)
        printf1(ctx, "\nError: can't re-open: %s\n",ctx->PSFName);
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  xdelete     Delete plot objects from current PostScript file.           */
/*              xdelete = list of plot object numbers;                      */
/*              Return 0 if OK, -1 if not.                                  */

int xdelete(TDAContext *ctx)
{
    register int i;
    register char *p;
    int err,n,fnd,sflg,len,oflag;

    err = -1;
    oflag = 0;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Delete plot objects. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,7,1))     /* get parameters */
        goto XDELFin;

    if (xlog_check(ctx) != 0)
        goto XDELFin;

    printf1(ctx, "Deleting plot object(s): %d",(int)ctx->PMRHSX[0]);
    for (i = 1; i < ctx->PMRHSN; ++i)
        printf1(ctx, ",%d",(int)ctx->PMRHSX[i]);
    newline(ctx);

    if (alloc_acc(ctx, RLMaxDef))
        goto XDELFin;

    ps_close(ctx, 0,0,1);    /* close current plot file */

    if (!(ctx->PMFd = fopen(ctx->PSFName,OPEN_RD))) {   
        printf1(ctx, "Error: can't re-open %s.\n",ctx->PSFName);
        goto XDELFin;
    }
    oflag = 1;

    len = 0;
    while (fgets(ctx->AcC,RLMaxDef,ctx->PMFd))  
        len += (int)strlen(ctx->AcC) + 1;

    if (alloc_acd(ctx, len + 1))
        goto XDELFin;

    ctx->PSONUM = 0;
    p = ctx->AcD;
    fseek(ctx->PMFd,0,0);
    sflg = 0;
    while (fgets(p,RLMaxDef,ctx->PMFd)) {

        if (!strncmp(p,"showpage",8))  
            break;

        if (sscanf(p,"%%#%d",&n) == 1 && n >= 1) {
            fnd = -1;
            for (i = 0; i < ctx->PMRHSN; ++i) {
                if (n == (int)ctx->PMRHSX[i]) {
                    fnd = n;
                    sflg = 1;
                    ctx->PMRHSX[i] *= -1.0;
                    break;
                }
            }
            if (fnd < 0) {
                sflg = 0;
                if (ctx->PSONUM < n)
                    ctx->PSONUM = n;
            }
        }
        if (sflg == 0)  
            p += strlen(p);
    }
    *p = '\0';

    err = -2;
    if (!(ctx->PSFd = fopen(ctx->PSFName,OPEN_WR)))    
        goto XDELFin;

    fwrite(ctx->AcD,(size_t)((int)(p - ctx->AcD)),1,ctx->PSFd);   
    ctx->PSFFlg = 2;

    n = 0;
    for (i = 0; i < ctx->PMRHSN; ++i) {
        if (ctx->PMRHSX[i] > 0.0) {
            if (n == 0) 
                printf1(ctx, "Ignored: %d",(int)ctx->PMRHSX[i]);
            else
                printf1(ctx, ",%d",(int)ctx->PMRHSX[i]);
            n++;
        }
    }
    if (n)
        newline(ctx);
    err = 0;

XDELFin:
    if (err == -2)  
        printf1(ctx, "\nError: can't re-open: %s\n",ctx->PSFName);
    if (oflag) {
        fclose(ctx->PMFd);
        ctx->PMFDef = 0;
    }
    p_clean(ctx);
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

int xdens(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,iflag,r,nmax;
    double tmp,xmin,xmax,ymin,ymax,d;

    err = -1;

    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "XDens. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))     /* get parameters */
        goto XDENSFin;

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (ctx->PMLT < 1 || ctx->PMLT > 9)
        ctx->PMLT = 1;

    if (ctx->PSFFlg)  
        ps_close(ctx, 0,0,0);

    xmin = xmax = get_data(ctx, ctx->PMVIdx[0],0);

    iflag = 1;
    for (i = 0; i < ctx->NOC; ++i) {
        for (j = 0; j < ctx->PMNV; ++j) {
            tmp = get_data(ctx, ctx->PMVIdx[j],i);       
            if (iflag) {
                n = (int)tmp;
                if ((double)n != tmp)
                    iflag = 0;
            }
            xmin = dmin(ctx, xmin,tmp);
            xmax = dmax(ctx, xmax,tmp);
        }
    }
/*
printf1(ctx, "iflag=%d minx=%lg %lg\n",iflag,xmin,xmax);
*/
    if (iflag) {            /* integer values, plot frequencies */

        printf1(ctx, "Frequency plot of variable(s): %s",ctx->VName[ctx->PMVIdx[0]]);
        for (i = 1; i < ctx->PMNV; ++i)
            printf1(ctx, ",%s",ctx->VName[ctx->PMVIdx[i]]);
        newline(ctx);

        if (ctx->PMLWFlg == 0)
            ctx->PMLW = 0.8;

        if (alloc_ack(ctx, ctx->NOC + 1))     /* value */
            goto XDENSFin;

        if (alloc_acn(ctx, ctx->NOC + 1))     /* freq */
            goto XDENSFin;

        if (alloc_aci(ctx, ctx->NOC + 1))     /* bf pointer */
            goto XDENSFin;

        nmax = 0;
        for (i = 0; i < ctx->PMNV; ++i) {    

            r = cfreq1(ctx, 1,ctx->PMVIdx + i,ctx->NOC,ctx->AcK,ctx->AcN,ctx->AcI,-1,&tmp);

            if (r <= 0) {
                printf1(ctx, "Error: can't create frequency distribution.\n");
                goto XDENSFin;
            }
            for (j = 1; j <= r; ++j) {
                k = ctx->AcI[j];
                nmax = imax(ctx, nmax,ctx->AcN[k]);
            }
        }
        ymin = 0.0;
        if (ctx->PMOPT == 1)
            ymax = (double)nmax / (double)ctx->NOC;
        else
            ymax = (double)nmax;

        if (x_cpsf(ctx, xmin - 1.0,xmax + 1.0,ymin,ymax,1))
            goto XDENSFin;

        d = 2.0 * ctx->UXLen * ctx->PMLW / ctx->PXLen;

        for (i = 0; i < ctx->PMNV; ++i) {    

            r = cfreq1(ctx, 1,ctx->PMVIdx + i,ctx->NOC,ctx->AcK,ctx->AcN,ctx->AcI,-1,&tmp);

            if (r <= 0) {
                printf1(ctx, "Error: can't create frequency distribution.\n");
                goto XDENSFin;
            }
            fprintf(ctx->PSFd,"\n%%#%d: %s [%s]\n",++ctx->PSONUM,ctx->CmdBuf,ctx->VName[ctx->PMVIdx[i]]);         
            fprintf(ctx->PSFd,"gsave\n");
            ps_lwidth(ctx, ctx->PMLW);
            ps_ltyp(ctx, ctx->PMLT);

            for (j = 1; j <= r; ++j) {
                k = ctx->AcI[j];
                tmp = (double)ctx->AcN[k];
                if (ctx->PMOPT == 1)
                    tmp /= (double)ctx->NOC;

                ps_2dplot(ctx, (double)ctx->AcK[k] + d * (double)i,0.0,0);
                ps_2dplot(ctx, (double)ctx->AcK[k] + d * (double)i,tmp,1);
                fprintf(ctx->PSFd,"stroke\n");
            }
            fprintf(ctx->PSFd,"grestore\n");
            if (++ctx->PMLT > 9)
                ctx->PMLT = 1;
        }
    }
    err = 0;

XDENSFin:
    p_clean(ctx);
    return(err);
}

