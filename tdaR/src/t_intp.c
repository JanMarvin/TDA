/****************************************************************************/
/*  t_intp                                                                  */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-2002 Goetz Rohwer. All rights reserved.         */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gf.h"
#include "t_psf.h"
#include "t_plot.h"
#include "t_plot3.h"
#include "t_sd.h"
#include "tda_context.h"

/*  functions in t_intp.c */

int triang(TDAContext *ctx);
int intp(TDAContext *ctx);
int intp_grid(TDAContext *ctx, int alg,int n,double *x,double *y,double *z,int nu,int nv, double *u,double *v,double *w,int ncp);
void intp_dbox2(TDAContext *ctx, int n,double *x,double *y);
void intp_dbox3(TDAContext *ctx, int n,double *x,double *y,double *z);
void intp_info(TDAContext *ctx, int n,int nb,int na,int nt);

int intp_trmesh(TDAContext *ctx, int n,double *x,double *y,int *iadj,int *iend);
int intp_adnode(TDAContext *ctx, int kk,double *x,double *y,int *iadj,int *iend); 
void intp_shiftd(TDAContext *ctx, int nfrst,int nlast,int kk,int *iarr); 
void intp_bdyadd(TDAContext *ctx, int kk,int i1,int i2,int *iadj,int *iend);
void intp_intadd(TDAContext *ctx, int kk,int i1,int i2,int i3,int *iadj,int *iend);
void intp_swap(TDAContext *ctx, int nin1,int nin2,int nout1,int nout2,int *iadj,int *iend);
void intp_trfind(TDAContext *ctx, int nst,double px,double py,double *x,double *y,int *iadj, int *iend,int *i1,int *i2,int *i3);
int intp_swptst(TDAContext *ctx, int in1,int in2,int io1,int io2,double *x,double *y);            
int intp_index(TDAContext *ctx, int nvertx,int nabor,int *iadj,int *iend);            
int intp_left(TDAContext *ctx, double x1,double y1,double x2,double y2,double x0,double y0);
int intp_trmtst(TDAContext *ctx, int n,double *x,double *y,int *iadj,int *iend,double tol); 
void intp_trmtst_emsg(TDAContext *ctx, int r);
int intp_circum(TDAContext *ctx, double x1,double x2,double x3, double y1,double y2,double y3,double *cx,double *cy);
int intp_bnodes(TDAContext *ctx, int n,int *iadj,int *iend,int *na,int *nt,int *nodes);
int intp_coords(TDAContext *ctx, double x,double y,double x1,double x2,double x3, double y1,double y2,double y3,double *r);
int intp_intrc0(TDAContext *ctx, int n,double px,double py,double *x,double *y,double *z, int *iadj,int *iend,int *ist,double *pz);
int intp_intrc1(TDAContext *ctx, int n,double px,double py,double *x,double *y,double *z, int *iadj,int *iend,int *ist,double *pz);
void intp_rotate(TDAContext *ctx, int n,double c,double s,double *x,double *y);
void intp_givens(TDAContext *ctx, double *a,double *b,double *c,double *s);
void intp_setup(TDAContext *ctx, double xk,double yk,double zk,double xi,double yi,double zi, double s1,double s2,double r,double *row);
int intp_tval(TDAContext *ctx, double x,double y,double x1,double x2,double x3, double y1,double y2,double y3,double z1,double z2,double z3, double zx1,double zx2,double zx3,double zy1,double zy2,double zy3, double *w,double *wx,double *wy,int iflag);
int intp_gradl(TDAContext *ctx, int n,int k,double *x,double *y,double *z,int *iadj,int *iend, double *dx,double *dy);
int intp_getnp(TDAContext *ctx, double *x,double *y,int *iadj,int *iend,int l,int *npts, double *ds);

int ak_idsfft(TDAContext *ctx, int ncp,int ndp,double * xd,double *yd,double *zd, int nxi,int nyi,double * xi,double *yi,double *zi,int *iwk,double *wk);
int ak_idtang(TDAContext *ctx, int ndp,double *xd,double *yd,int *nt,int *ipt,int *nl,int *ipl, int *iwl,int *iwp,double *wk);
double ak_dsqf(TDAContext *ctx, double u1,double v1,double u2,double v2);
double ak_side(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3);
int ak_idxchg(TDAContext *ctx, double *x,double *y,int i1,int i2,int i3,int i4); 
void ak_idptip(TDAContext *ctx, double *xd,double *yd,double *zd,int nt,int *ipt,int nl, int *ipl,double *pdd,int iti,double xii,double yii,double *zii);
void ak_idpdrv(TDAContext *ctx, int ndp,double *xd,double *yd,double *zd,int ncp, int *ipc,double *pd);
void ak_idgrid(TDAContext *ctx, double *xd,double *yd,int nt,int *ipt,int nl,int *ipl, int nxi,int nyi,double *xi,double *yi,int *ngp,int *igp);
double ak_side1(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3);
double ak_spdt1(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3);
int ak_idcldp(TDAContext *ctx, int ndp,double *xd,double *yd,int ncp,int *ipc);

/* ------------------------------------------------------------------------ */
/*  triang.     Triangulation.                                              */
/*                                                                          */  
/*              triang(                                                     */
/*                  prn=...,    print option, def. 0 (nothing)              */
/*                  df=...,     output file                                 */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    print format, def. 10.4                     */
/*                                                                          */
/*                  plot=...,   plot option, def 0 (nothing)                */
/*                  zmin=...,   minimum  of z coordinate, def. 0.0          */
/*                  lt=...,     line type, def. 1 (solid line)              */
/*                  lw=...,     line width, def. 0.2 (mm)                   */
/*                  s=...,      number of marker symbol                     */  
/*                  fs=...,     size of marker symbol, def. 2 (mm)          */      
/*                  nc=...,     no clipping option, def. 0                  */            
/*                  gs=...,     grey scale value, def. 1 (white)            */
/*                              (only with plot = 2)                        */
/*              ) = X,Y [,Z];                                               */
/*                                                                          */
/*  The prn options require the definition of an output file with the df    */
/*  parameter. The following options are available:                         */
/*                                                                          */
/*  1   Print the triangulation structure as an edge list. The nodes are    */
/*      numbered by the case values of the points in TDA's data matrix.     */
/*  2   Like option 1 but in addition print the x and y values of the       */
/*      points.                                                             */
/*  3   Print the points on the boundary (convex hull) of the point set.    */
/*                                                                          */
/*  The plot options require a valid 2- or 3-dimensional PostScript         */
/*  coordinate defined with the psetup or psetup3 command. The following    */
/*  options are available:                                                  */
/*                                                                          */  
/*  1   Plot of all arcs that make up the mesh.                             */
/*  2   Plot of the convex hull. This option allows to use the gs parameter */
/*      to fill the region with a grey scale value.                         */
/*                                                                          */
/*      In both cases, if s specifies a valid symbol and fs > 0, the        */
/*      symbols are drawn at the location of all points (plot = 1) or       */
/*      at the location of the points on the boundary (plot = 2).           */
/*                                                                          */
/*      Note: For a 3d plot, the mesh is drawn with a z coordinate defined  */
/*      by the zmin parameter. If the right-hand side specifies 3 variables */
/*      then, at each node, an additional line is drawn from zmin to the    */
/*      corresponding z coordinate given by the third variable.             */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int triang(TDAContext *ctx)
{
    register int i,j,k,l; 
    int err,r,ix,iy,iz = 0,nrec,na,nt,nb;
    double x,y;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
    printf1(ctx, "Triangulation of a point set. Current memory: %d bytes.\n\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 6,4,1))     /* get parameters */
        goto TRIFin;

    if (ctx->PMNV < 2) {
        printf1(ctx, "Error: need at least two variables on right-hand side.\n");
        goto TRIFin;
    }
    if (ctx->PMOPT != 1)
        ctx->PMOPT = 1;

    if (ctx->NOC < 3) {
        printf1(ctx, "Error: need at least three data points.\n");
        goto TRIFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto TRIFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto TRIFin;

    if (ctx->PMNV >= 3) {
        iz = ctx->PMVIdx[2];
        if (alloc_acz(ctx, ctx->NOC + 1))
            goto TRIFin;
    }

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcX[i + 1] = get_data(ctx, ix,i);
        ctx->AcY[i + 1] = get_data(ctx, iy,i);
        if (ctx->PMNV >= 3)
            ctx->AcZ[i + 1] = get_data(ctx, iz,i);
    }
    if (ctx->PMFmtF == 0)  
        pmfmt(ctx, 10,4);

    if (ctx->PMNV == 2)              /* write bounding box of input data */
        intp_dbox2(ctx, ctx->NOC,ctx->AcX,ctx->AcY); 
    else
        intp_dbox3(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcZ); 
    newline(ctx);

    if (alloc_aci(ctx, 6 * ctx->NOC + 1))
        goto TRIFin;
    if (alloc_acj(ctx, ctx->NOC + 1))
        goto TRIFin;
    if (alloc_ack(ctx, ctx->NOC + 1))
        goto TRIFin;

    r = intp_trmesh(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcJ);
    if (r) {               
        printf1(ctx, "Error: all points are collinear.\n");
        goto TRIFin;
    }

    /* check for valid data structure */
   
    if ((r = intp_trmtst(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcJ,ctx->EPSI1)) > 0) {
        intp_trmtst_emsg(ctx, r);
        goto TRIFin;
    }
/**
    tda_out("IADJ: ");
    for (i = 1; i <= 6*NOC; ++i)
        tda_out("%d ",AcI[i]);

    tda_out("\nIEND: ");
    for (i = 1; i <= NOC; ++i)
        tda_out("%d ",AcJ[i]);

    newline(ctx);
***/


    nb = intp_bnodes(ctx, ctx->NOC,ctx->AcI,ctx->AcJ,&na,&nt,ctx->AcK);
    intp_info(ctx, ctx->NOC,nb,na,nt);
    newline(ctx);

    if (ctx->PMF1Def != 0 && ctx->PMPRNO >= 1 && ctx->PMPRNO <= 3) { /* print to output file */

        printf1(ctx, "Processing print option %d.\n",ctx->PMPRNO);
        nrec = 0;

        if (ctx->PMPRNO == 1 || ctx->PMPRNO == 2) {       /* edge list */

            for (i = 1; i <= ctx->NOC; ++i) {
                if (i == 1)
                    k = 1;
                else
                    k = ctx->AcJ[i - 1] + 1;

                l = ctx->AcJ[i];
                while (k <= l) {
                    j = ctx->AcI[k];
                    if (j > i) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j);

                        if (ctx->PMPRNO == 2) {
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[j]);
                            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[j]);
                        }
                        fprintf(ctx->PMF1d,"\n");
                        nrec++;
                    }
                    k++;
                }
            }
        }
        else if (ctx->PMPRNO == 3) {       /* boundary (convex hull) */

            for (i = 1; i <= nb; ++i) {
                j = ctx->AcK[i];
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[j]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[j]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
        }
        printf1(ctx, "%d records written to: %s\n\n",nrec,ctx->PMF1dName);
    }
    if (ctx->PSFFlg == 2 && ctx->PMPLOT >= 1 && ctx->PMPLOT <= 2) {

        printf1(ctx, "Processing plot option %d.\n",ctx->PMPLOT);
        if (ctx->PS3DFlg)  
            printf1(ctx, "3d plot (zmin = %g)\n",ctx->PMZMin);

        fprintf(ctx->PSFd,"\n%%#%d: triang (plot=%d)\n",++ctx->PSONUM,ctx->PMPLOT);
        fprintf(ctx->PSFd,"gsave\n");
        if (ctx->PMNC != 1)
            set_clip(ctx);
        ps_ltyp(ctx, ctx->PMLT);
        ps_lwidth(ctx, ctx->PMLW);

        if (ctx->PMPLOT == 1) {

            for (i = 1; i <= ctx->NOC; ++i) {
                if (i == 1)
                    k = 1;
                else
                    k = ctx->AcJ[i - 1] + 1;

                l = ctx->AcJ[i];
                while (k <= l) {
                    j = ctx->AcI[k];
                    if (j > i) {
                        if (ctx->PS3DFlg == 0) {
                            ps_2dplot(ctx, ctx->AcX[i],ctx->AcY[i],0);
                            ps_2dplot(ctx, ctx->AcX[j],ctx->AcY[j],1);
                        }
                        else {
                            ps_3dplot(ctx, ctx->AcX[i],ctx->AcY[i],ctx->PMZMin,0,0,ctx->PMNC);
                            ps_3dplot(ctx, ctx->AcX[j],ctx->AcY[j],ctx->PMZMin,1,0,ctx->PMNC);
                        }
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    k++;
                }
            }
        }
        else if (ctx->PMPLOT == 2) {
            k = 0;
            for (i = 1; i <= nb; ++i) {
                j = ctx->AcK[i];
                if (ctx->PS3DFlg == 0)  
                    ps_2dplot(ctx, ctx->AcX[j],ctx->AcY[j],k);
                else
                    ps_3dplot(ctx, ctx->AcX[j],ctx->AcY[j],ctx->PMZMin,k,0,ctx->PMNC);

                k = 1;
            }
            j = ctx->AcK[1];

            if (ctx->PS3DFlg == 0)  
                ps_2dplot(ctx, ctx->AcX[j],ctx->AcY[j],1);
            else
                ps_3dplot(ctx, ctx->AcX[j],ctx->AcY[j],ctx->PMZMin,1,0,ctx->PMNC);

            if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
                fprintf(ctx->PSFd,"gsave\n");
                ps_fill(ctx, ctx->PMGS);
                fprintf(ctx->PSFd,"grestore\n");
            }
            fprintf(ctx->PSFd,"stroke\n");
        }

        if (ctx->PS3DFlg && ctx->PMNV >= 3) {  /* plot heights */

            for (i = 1; i <= ctx->NOC; ++i) {
                if (ctx->PMPLOT == 1) {
                    ps_3dplot(ctx, ctx->AcX[i],ctx->AcY[i],ctx->PMZMin,0,0,ctx->PMNC);
                    ps_3dplot(ctx, ctx->AcX[i],ctx->AcY[i],ctx->AcZ[i],1,0,ctx->PMNC);
                }
                else {
                    if (i > nb)
                        break;
                    j = ctx->AcK[i];
                    ps_3dplot(ctx, ctx->AcX[j],ctx->AcY[j],ctx->PMZMin,0,0,ctx->PMNC);
                    ps_3dplot(ctx, ctx->AcX[j],ctx->AcY[j],ctx->AcZ[j],1,0,ctx->PMNC);
                }
                fprintf(ctx->PSFd,"stroke\n");
            }   
        }

        if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

            fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);
                                
            for (i = 1; i <= ctx->NOC; ++i) {
                if (ctx->PMPLOT == 1) {
                    if (ctx->PS3DFlg == 0) {
                        x = ps_2dx(ctx, ctx->AcX[i]);       
                        y = ps_2dy(ctx, ctx->AcY[i]);       
                    }
                    else {
                        ps_3dprj(ctx, ctx->AcX[i],ctx->AcY[i],ctx->PMZMin,&x,&y);
                        x = ps_2dx(ctx, x);
                        y = ps_2dy(ctx, y);
                    }
                }   
                else {
                    if (i > nb)
                        break;
                    j = ctx->AcK[i];
                    if (ctx->PS3DFlg == 0) {
                        x = ps_2dx(ctx, ctx->AcX[j]);       
                        y = ps_2dy(ctx, ctx->AcY[j]);       
                    }
                    else {
                        ps_3dprj(ctx, ctx->AcX[j],ctx->AcY[j],ctx->PMZMin,&x,&y);
                        x = ps_2dx(ctx, x);
                        y = ps_2dy(ctx, y);
                    }
                }
                if (ctx->PMNC)
                    upd_bbox(ctx, 1,x,y);                 /* update bounding box */
                ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            }   
        }
        fprintf(ctx->PSFd,"grestore\n");
        printf1(ctx, "Plot commands written to: %s\n",ctx->PSFName);
    }
    err = 0;

TRIFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  intp        Interpolation of a point set.                               */
/*                                                                          */  
/*              intp(                                                       */
/*                  alg=...,    algorithm, def. 1                           */
/*                  rx=...,     range of x values (required)                */
/*                  ry=...,     range of y values (required)                */
/*                  df=...,     output file (required)                      */
/*                  fmt=...,    print format, def. 10.4                     */
/*                                                                          */
/*              ) = X,Y,Z;                                                  */
/*                                                                          */
/*                                                                          */
/*  The command assumes a set of points given by x(i), y(i), z(i), 1=1,n    */
/*  defined by the variables X, Y and Z, and n equals NOC, the number of    */
/*  currently selected cases from the TDA data matrix.                      */
/*                                                                          */
/*  The command constructs a Thiessen-triangulation based on x(i) and y(i)  */
/*  using routines from CALGO 624 (Robert Renka).                           */
/*                                                                          */
/*  Then, for the rectangular x-y-grid defined by rx and ry, the command    */  
/*  interpolates values for z. The grid must be defined as follows:         */
/*                                                                          */
/*      rx = rxa (rxd) rxb,                                                 */
/*      ry = rya (ryd) ryb,                                                 */
/*                                                                          */
/*  This defines values                                                     */
/*                                                                          */
/*      x = rxa + i * rxd  (i = 0,1,2,... until x > rxb)                    */
/*      y = rya + j * ryd  (j = 0,1,2,... until y > ryb)                    */
/*                                                                          */  
/*  The following interpolation methods can be selected with the alg        */
/*  parameter:                                                              */
/*                                                                          */
/*  (1) ACM algorithm 624 (Renka). Linear interpolation.                    */
/*  (2) ACM algorithm 624 (Renka). C1 interpolation                         */
/*  (3) ACM algorithm 526 (Akima). For this option, ncp provides the        */
/*      number of additional points, must be >= 2 and <= 25.                */
/*                                                                          */
/*  The resulting values are written to the output file defined by the      */
/*  df parameter using the print format defined by fmt. Each record         */
/*  contains three entries:                                                 */
/*                                                                          */
/*  (1)  x value from grid definition.                                      */  
/*  (2)  y value from grid definition.                                      */
/*  (3)  the interpolated z value for the (x,y) point.                      */
/*                                                                          */
/*  Note: entry (4) is only present for opt = 1 and 2.                      */
/*                                                                          */
/*  Return 0 if OK, otherwise -1.                                           */

int intp(TDAContext *ctx)
{
    register int i,j;
    int err,r,ix,iy,iz,nrec,ncp,nu,nv,nuu,nvv;
    double u,v,w;  

    ncp = 3;

    err = -1;         
    if (check_cmd(ctx, 1))
        return(-1);
    printf1(ctx, "Interpolation of a point set. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 4,4,1))     /* get parameters */
        goto INTPFin;

    if (ctx->PMNV < 3) {
        printf1(ctx, "Error: need three variables on right-hand side.\n");
        goto INTPFin;
    }
    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;

    if (ctx->PMRXFlg == 0 || ctx->PMRYFlg == 0) {
        printf1(ctx, "Error: need valid rx and ry parameters.\n");
        goto INTPFin;
    }
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto INTPFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    iz = ctx->PMVIdx[2];

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto INTPFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto INTPFin;
    if (alloc_acz(ctx, ctx->NOC + 1))
        goto INTPFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcX[i + 1] = get_data(ctx, ix,i);
        ctx->AcY[i + 1] = get_data(ctx, iy,i);
        ctx->AcZ[i + 1] = get_data(ctx, iz,i);
    }
    if (ctx->PMFmtF == 0)  
        pmfmt(ctx, 10,4);

    printf1(ctx, "Algorithm: %d\n\n",ctx->PMALG);

    intp_dbox3(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcZ);     /* write bounding box of input data */

    printf1(ctx, "\nRequested range of interpolated points.\n");
    printf1(ctx, "x - from: %14.6f to %14.6f  increment: %14.6f\n",(double)(ctx->PMRXA),(double)(ctx->PMRXB),(double)(ctx->PMRXD));
    printf1(ctx, "y - from: %14.6f to %14.6f  increment: %14.6f\n\n",(double)(ctx->PMRYA),(double)(ctx->PMRYB),(double)(ctx->PMRYD));

    nu = (int)((ctx->PMRXB - ctx->PMRXA) / ctx->PMRXD) + 4;
    if (alloc_acu(ctx, nu + 1))
        goto INTPFin;

    nv = (int)((ctx->PMRYB - ctx->PMRYA) / ctx->PMRYD) + 4;
    if (alloc_acv(ctx, nv + 1))
        goto INTPFin;

    nuu = nvv = 0;
    for (u = (double)ctx->PMRXA; u <= (double)(ctx->PMRXB) + ctx->EPSI1; u += (double)(ctx->PMRXD)) {
        if (++nuu > nu)
            break;
        ctx->AcU[nuu] = u;
    }
    for (v = (double)ctx->PMRYA; v <= (double)(ctx->PMRYB) + ctx->EPSI1; v += (double)(ctx->PMRYD)) {
        if (++nvv > nv)
            break;
        ctx->AcV[nvv] = v;
    }
    if (alloc_acw(ctx, nuu * nvv + 2))
        goto INTPFin;

    r = intp_grid(ctx, ctx->PMALG,ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcZ,nuu,nvv,ctx->AcU,ctx->AcV,ctx->AcW,ncp);
    if (r)
        goto INTPFin;

    nrec = 0;
    for (i = 1; i <= nuu; ++i) {
        u = ctx->AcU[i];
        for (j = 1; j <= nvv; ++j) {
            v = ctx->AcV[j];
            w = ctx->AcW[(i - 1) * nvv + j];
                
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,u);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,v);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,w);
            fprintf(ctx->PMF1d,"\n");
            nrec++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

INTPFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  intp_grid(alg,n,x,y,z,nu,nv,u,v,w,ncp)                                  */
/*                                                                          */
/*  Get values for a rectangular grid by interpolation or smoothing.        */  
/*  Input parameters:                                                       */
/*                                                                          */
/*  alg = selected algorithm.                                               */
/*  (1) ACM algorithm 624 (Renka). Linear interpolation.                    */
/*  (2) ACM algorithm 624 (Renka). C1 interpolation                         */
/*  (3) ACM algorithm 526 (Akima). For this option, ncp provides the        */
/*      number of additional points, must be >= 2 and <= 25.                */
/*                                                                          */
/*                                                                          */
/*  x[i], y[i], z[i] for i = 1,...,n  are the data points, irregularly      */
/*  distributed in the x-y plane.                                           */
/*                                                                          */
/*  nu and nv provide the number of grid lines in the x and y direction,    */
/*  respectively. The coordinates are given by u[i] (i = 1,...,nu) and      */
/*  v[j] (j = 1,...,nv).                                                    */
/*                                                                          */
/*  w is a doubly indexed array w[i,j] = w[(i - 1) * nv + j] that will be   */
/*  used to save the estimated z-values at the grid points.                 */
/*                                                                          */
/*  Note: This function allocates and uses the global arrays:               */
/*  AcI, AcJ, AcK (for alg = 1 or 2) and AcI and AcTmp (for alg = 3).       */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int intp_grid(TDAContext *ctx, int alg,int n,double *x,double *y,double *z,int nu,int nv, double *u,double *v,double *w,int ncp)
{
    register int i,j;
    int r,na,nb,nt,ist;
    double xx,yy,zz;

    if (alg < 1 || alg > 3) {
        printf1(ctx, "Error: algorithm %d not available.\n",alg);
        return(-1);
    }
    if (alg == 1 || alg == 2) {
        if (n < 3) {
            printf1(ctx, "Error: need at least three data points.\n");
            return(-1);
        }
        if (alloc_aci(ctx, 6 * n + 1))
            return(-1);      
        if (alloc_acj(ctx, n + 1))
            return(-1);       
        if (alloc_ack(ctx, n + 1))
            return(-1);      

        r = intp_trmesh(ctx, n,x,y,ctx->AcI,ctx->AcJ);
        if (r) {               
            printf1(ctx, "Error: all points are collinear.\n");
            return(-1);   
        }

        /* check for valid data structure */
   
        if ((r = intp_trmtst(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcI,ctx->AcJ,ctx->EPSI1)) > 0) {
            intp_trmtst_emsg(ctx, r);
            return(-1);    
        }
        nb = intp_bnodes(ctx, ctx->NOC,ctx->AcI,ctx->AcJ,&na,&nt,ctx->AcK);
        intp_info(ctx, n,nb,na,nt);

        ist = 1;
        for (i = 1; i <= nu; ++i) {
            xx = u[i];
            for (j = 1; j <= nv; ++j) {
                yy = v[j];

                if (alg == 1)  
                    r = intp_intrc0(ctx, n,xx,yy,x,y,z,ctx->AcI,ctx->AcJ,&ist,&zz);
                else           
                    r = intp_intrc1(ctx, n,xx,yy,x,y,z,ctx->AcI,ctx->AcJ,&ist,&zz);
                
                w[(i - 1) * nv + j] = zz;
            }
        }
        alloc_aci(ctx, 0);
        alloc_acj(ctx, 0);
        alloc_ack(ctx, 0);
    }
    else if (alg == 3) {
        if (n < 4) {
            printf1(ctx, "Error: need at least four data points.\n");
            return(-1);
        }
        if (ncp < 2 || ncp > 25) {
            printf1(ctx, "Error in ncp parameter.\n");
            return(-1);
        }
        if (alloc_actmp(ctx, 5 * n + 10))
            return(-1);   

        if (alloc_aci(ctx, imax(ctx, 31,27 + ncp) * n + nu * nv + 10))
            return(-1);     

        r = ak_idsfft(ctx, ncp,n,x,y,z,nu,nv,u,v,w,ctx->AcI,ctx->AcTmp);
        if (r) {
            printf1(ctx, "Error: ak_idsfft (r=%d).\n",r);
            return(-1);     
        }
        alloc_aci(ctx, 0);
        alloc_actmp(ctx, 0);
    }
    return(0);
}

/* -##=-------------------------------------------------------------------- */
/*  intp_dbox2(n,x,y)  Write min and max to standard output. Assume arrays  */
/*                     x[i], y[i], i = 1,...,n.                             */

void intp_dbox2(TDAContext *ctx, int n,double *x,double *y)
{
    register int i;
    double xmin,xmax,ymin,ymax;

    xmin = xmax = x[1];
    ymin = ymax = y[1];
    for (i = 2; i <= n; ++i) {
        xmin = dmin(ctx, xmin,x[i]);
        xmax = dmax(ctx, xmax,x[i]);
        ymin = dmin(ctx, ymin,y[i]);
        ymax = dmax(ctx, ymax,y[i]);
    }
    printf1(ctx, "Bounding box of input data.\n");
    printf1(ctx, "x minimum: %14.6f  maximum: %14.6f\n",xmin,xmax);
    printf1(ctx, "y minimum: %14.6f  maximum: %14.6f\n",ymin,ymax);
}

/* -##=-------------------------------------------------------------------- */
/*  intp_dbox3(n,x,y,z) Write min and max to standard output. Assume arrays */
/*                      x[i], y[i], z[i],  i = 1,...,n.                     */

void intp_dbox3(TDAContext *ctx, int n,double *x,double *y,double *z)
{
    register int i;
    double xmin,xmax,ymin,ymax,zmin,zmax;

    xmin = xmax = x[1];
    ymin = ymax = y[1];
    zmin = zmax = z[1];
    for (i = 2; i <= n; ++i) {
        xmin = dmin(ctx, xmin,x[i]);
        xmax = dmax(ctx, xmax,x[i]);
        ymin = dmin(ctx, ymin,y[i]);
        ymax = dmax(ctx, ymax,y[i]);
        zmin = dmin(ctx, zmin,z[i]);
        zmax = dmax(ctx, zmax,z[i]);
    }
    printf1(ctx, "Bounding box of input data.\n");
    printf1(ctx, "x minimum: %14.6f  maximum: %14.6f\n",xmin,xmax);
    printf1(ctx, "y minimum: %14.6f  maximum: %14.6f\n",ymin,ymax);
    printf1(ctx, "z minimum: %14.6f  maximum: %14.6f\n",zmin,zmax);
}

/* -##=-------------------------------------------------------------------- */
/*  intp_info(n,x,y)    Write info about triangulation to standard output.  */

void intp_info(TDAContext *ctx, int n,int nb,int na,int nt)
{
    printf1(ctx, "Number of data points: %d\n",n);
    printf1(ctx, "Number of points on the boundary: %d\n",nb);
    printf1(ctx, "Number of arcs in the mesh: %d\n",na);
    printf1(ctx, "Number of triangles in the mesh: %d\n",nt);
}

/* -##=-------------------------------------------------------------------- */
/*  intp_trmesh(n,x,y,iadj,iend,ier)                                        */
/*                                                                          */
/*  THIS ROUTINE CREATES A THIESSEN TRIANGULATION OF N                      */
/*  ARBITRARILY SPACED POINTS IN THE PLANE REFERRED TO AS                   */
/*  NODES.  THE TRIANGULATION IS OPTIMAL IN THE SENSE THAT IT               */
/*  IS AS NEARLY EQUIANGULAR AS POSSIBLE.  TRMESH IS PART OF                */
/*  AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES SUBROUTINES                */
/*  TO REORDER THE NODES, ADD A NEW NODE, DELETE AN ARC, PLOT               */
/*  THE MESH, AND PRINT THE DATA STRUCTURE.                                 */
/*  UNLESS THE NODES ARE ALREADY ORDERED IN SOME REASONABLE                 */
/*  FASHION, THEY SHOULD BE REORDERED BY SUBROUTINE REORDR FOR              */
/*  INCREASED EFFICIENCY BEFORE CALLING TRMESH.                             */
/*                                                                          */
/*  INPUT PARAMETERS -     N - NUMBER OF NODES IN THE MESH.                 */
/*                             N .GE. 3.                                    */
/*                                                                          */
/*                       X,Y - N-VECTORS OF COORDINATES.                    */
/*                             (X(I),Y(I)) DEFINES NODE I.                  */
/*                                                                          */
/*                      IADJ - VECTOR OF LENGTH .GE. 6*N-9.                 */
/*                                                                          */
/*                      IEND - VECTOR OF LENGTH .GE. N.                     */
/*                                                                          */
/*  N, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.                            */
/*                                                                          */
/*  OUTPUT PARAMETERS - IADJ - ADJACENCY LISTS OF NEIGHBORS IN              */
/*                             COUNTERCLOCKWISE ORDER.  THE                 */
/*                             LIST FOR NODE I+1 FOLLOWS THAT               */
/*                             FOR NODE I WHERE X AND Y DEFINE              */
/*                             THE ORDER.  THE VALUE 0 DENOTES              */
/*                             THE BOUNDARY (OR A PSEUDO-NODE               */
/*                             AT INFINITY) AND IS ALWAYS THE               */
/*                             LAST NEIGHBOR OF A BOUNDARY                  */
/*                             NODE.  IADJ IS UNCHANGED IF IER              */
/*                             .NE. 0.                                      */
/*                                                                          */
/*                      IEND - POINTERS TO THE ENDS OF                      */
/*                             ADJACENCY LISTS (SETS OF                     */
/*                             NEIGHBORS) IN IADJ.  THE                     */
/*                             NEIGHBORS OF NODE 1 BEGIN IN                 */
/*                             iadj[1).  FOR K .GT. 1, THE                  */
/*                             NEIGHBORS OF NODE K BEGIN IN                 */
/*                             iadj[iend[K-1)+1) AND K HAS                  */
/*                             iend[K) - iend[K-1) NEIGHBORS                */
/*                             INCLUDING (POSSIBLY) THE                     */
/*                             BOUNDARY.  iadj[iend[K)) .EQ. 0              */
/*                             IFF NODE K IS ON THE BOUNDARY.               */
/*                             IEND IS UNCHANGED IF IER = 1.                */
/*                             IF IER = 2 IEND CONTAINS THE                 */
/*                             INDICES OF A SEQUENCE OF N                   */
/*                             NODES ORDERED FROM LEFT TO                   */
/*                             RIGHT WHERE LEFT AND RIGHT ARE               */
/*                             DEFINED BY ASSUMING NODE 1 IS                */
/*                             TO THE LEFT OF NODE 2.                       */
/*                                                                          */
/*  Return: 0 if successful, -1 if n < 3 or all points are collinear.       */
/*                                                                          */
/*  Functions called: intp_shiftd(), intp_adnode().                         */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_trmesh(TDAContext *ctx, int n,double *x,double *y,int *iadj,int *iend)
{
    register int i;
    int nn,k,km1,nl,nr,ind,indx,n0,itemp,km1d2,kmi,kmin;
    double xl,yl,xr,yr,dxr,dyr,xk,yk,dxk,dyk,cprod,sprod;
                                                                              
    nn = n;                                                                   
    if (nn < 3)
        return(-1);                                                         
                                                                              
    iend[1] = 1;                                                              
    iend[2] = 2;                                                              
    xl = x[1];                                                                
    yl = y[1];                                                                
    xr = x[2];                                                                
    yr = y[2];                                                                
    k = 2;                                                                    
                                                                              
L1:
    dxr = xr - xl;                                                            
    dyr = yr - yl;                                                            
                                                                              
L2:
    if (k == nn)
        return(-1);                                                 
    km1 = k;                                                                  
    k = km1 + 1;                                                              
    xk = x[k];                                                                
    yk = y[k];                                                                
    dxk = xk - xl;                                                            
    dyk = yk - yl;                                                            
    cprod = dxr * dyk - dxk * dyr;                                            
    if (cprod > 0.0)
        goto L6;                                               
    if (cprod < 0.0)
        goto L8;                                               

    sprod = dxr * dxk + dyr * dyk;                                            
    if (sprod > 0.0)
        goto L3;                                               
                                                                              
    intp_shiftd(ctx, 1,km1,1,iend);                                              

    iend[1] = k;                                                              
    xl = xk;                                                                  
    yl = yk;                                                                  
    goto L1;                                                                  
                                                                              
L3:
    for (ind = 2; ind <= km1; ++ind) {
        n0 = iend[ind];                                                       

        sprod = (xl - x[n0]) * (xk - x[n0]) + (yl - y[n0]) * (yk - y[n0]);    
        if (sprod >= 0.0)
            goto L5;                                                          
    }                                                                         
    iend[k] = k;                                                              
    xr = xk;                                                                  
    yr = yk;                                                                  
    goto L1;                                                                  
                                                                              
L5:
    intp_shiftd(ctx, ind,km1,1,iend);                                            
    iend[ind] = k;                                                            
    goto L2;                                                                  
                                                                              
L6:
    km1d2 = km1 / 2;                                                          
    for (i = 1; i <= km1d2; ++i) {
        kmi = k - i;                                                          

        itemp = iend[i];                                                      
        iend[i] = iend[kmi];                                                  
        iend[kmi] = itemp;                                                    

    }                                                                         
                                                                              
L8:
    nl = iend[1];                                                             
    nr = iend[km1];                                                           
                                                                              
    for (ind = 1; ind <= km1; ++ind) {
        n0 = iend[ind];                                                       
        indx = 4 * n0;                                                        
        if (n0 >= nl)
            indx--;                                                           
        if (n0 >= nr)
            indx--;                                                           
        iadj[indx] = 0;                                                       
        indx--;                                                               

        if (ind < km1)
            iadj[indx] = iend[ind + 1];                             
        if (ind < km1)
            indx--;                                                

        iadj[indx] = k;                                                       
        if (ind == 1)
            goto L9;                                                
        iadj[indx - 1] = iend[ind - 1];                                       
L9:     ;                                                                     
    }                                                                         
    indx = 5 * km1 - 1;                                                       
    iadj[indx] = 0;                                                           
    for (ind = 1; ind <= km1; ++ind) {
        indx--;                                                               
        iadj[indx] = iend[ind];                                               
    }                                                                         
    indx = 0;                                                                 
    for (ind = 1; ind <= km1; ++ind) {
        indx += 4;                                                            
        if (ind == nl || ind == nr)
            indx--;                                                     
        iend[ind] = indx;                                                     
    }
    indx += k;                                                                
    iend[k] = indx;                                                           
                                                                              
    if (k == nn)
        return(0);                                                          
    kmin = k + 1;                                                             

    for (k = kmin; k <= nn; ++k) {
        intp_adnode(ctx, k,x,y,iadj,iend);                                  
    }
    return(0);                                                                
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_adnode(kk,x,y,iadj,iend,ier)                                       */
/*                                                                          */
/*  THIS ROUTINE ADDS NODE KK TO A TRIANGULATION OF A SET                   */
/*  OF POINTS IN THE PLANE PRODUCING A NEW TRIANGULATION.  A                */
/*  SEQUENCE OF EDGE SWAPS IS THEN APPLIED TO THE MESH,                     */
/*  RESULTING IN AN OPTIMAL TRIANGULATION.  ADNODE IS PART                  */
/*  OF AN INTERPOLATION PACKAGE WHICH ALSO PROVIDES ROUTINES                */
/*  TO INITIALIZE THE DATA STRUCTURE, PLOT THE MESH, AND                    */
/*  DELETE ARCS.                                                            */
/*                                                                          */
/*  INPUT PARAMETERS -   KK - INDEX OF THE NODE TO BE ADDED                 */
/*                            TO THE MESH.  KK .GE. 4.                      */
/*                                                                          */
/*                      X,Y - VECTORS OF COORDINATES OF THE                 */
/*                            NODES IN THE MESH.  (X(I),Y(I))               */
/*                            DEFINES NODE I FOR I = 1,..,KK.               */
/*                                                                          */
/*                     IADJ - SET OF ADJACENCY LISTS OF NODES               */
/*                            1,..,KK-1.                                    */
/*                                                                          */
/*                     IEND - POINTERS TO THE ENDS OF                       */
/*                            ADJACENCY LISTS IN IADJ FOR                   */
/*                            EACH NODE IN THE MESH.                        */
/*                                                                          */
/*  IADJ AND IEND MAY BE CREATED BY TRMESH.                                 */
/*                                                                          */
/*  KK, X, AND Y ARE NOT ALTERED BY THIS ROUTINE.                           */
/*                                                                          */
/*  Output: iadj and iend - UPDATED WITH THE ADDITION                       */
/*                                  OF NODE KK AS THE LAST                  */
/*                                  ENTRY.                                  */
/*                                                                          */
/*  Return: ERROR INDICATOR                                                 */
/*                                  IER = 0 IF NO ERRORS                    */
/*                                          WERE ENCOUNTERED.               */
/*                                  IER = 1 IF ALL NODES                    */
/*                                          (INCLUDING KK) ARE              */
/*                                          COLLINEAR.                      */
/*                                                                          */
/*  Functions called: intp_index(), intp_swptst(), intp_trfind(),           */
/*  intp_swap(), intp_intadd(), intp_bdyadd().                              */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_adnode(TDAContext *ctx, int kk,double *x,double *y,int *iadj,int *iend)  
{
    int k,km1,i1,i2,i3,indkf,indkl,nabor1,io1,io2,in1,indk1,ind2f,ind21;
    double xk,yk;

    k = kk;                                                                   
    km1 = k - 1;                                                              
    xk = x[k];                                                                
    yk = y[k];                                                                

    intp_trfind(ctx, km1,xk,yk,x,y,iadj,iend,&i1,&i2,&i3);                           

    if (i1 == 0)
        return(1);
                                                                   
    if (i3 == 0)
        intp_bdyadd(ctx, k,i1,i2,iadj,iend);                          

    if (i3 != 0)
        intp_intadd(ctx, k,i1,i2,i3,iadj,iend);                       
                                                                              
    indkf = iend[km1] + 1;                                                    
    indkl = iend[k];                                                          
    nabor1 = iadj[indkf];                                                     
    io2 = nabor1;                                                             
    indk1 = indkf + 1;                                                        
    io1 = iadj[indk1];                                                        

L1:
    ind2f = 1;                                                                
    if (io2 != 1)
        ind2f = iend[io2 - 1] + 1;                                  
    ind21 = intp_index(ctx, io2,io1,iadj,iend);                                         
    if (ind2f == ind21)
        goto L2;                                            
    in1 = iadj[ind21 - 1];                                                    
    goto L3;                                                                  
                                                                              
L2:
    ind21 = iend[io2];                                                        
    in1 = iadj[ind21];                                                        
    if (in1 == 0)
        goto L4;                                                  
                                                                              
L3:
    if (intp_swptst(ctx, in1,k,io1,io2,x,y) == 0)
        goto L4;                           

    intp_swap(ctx, in1,k,io1,io2,iadj,iend);                                     

    io1 = in1;                                                                
    indk1--;                                                                  
    indkf--;                                                                  
    goto L1;                                                                  
                                                                              
L4:
    if (io1 == nabor1)
        return(0);                                                          
    io2 = io1;                                                                
    indk1++;                                                                  
    if (indk1 > indkl)
        indk1 = indkf;                                      
    io1 = iadj[indk1];                                                        
    if (io1 != 0)
        goto L1;                                                  
    return(0);                                                                
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_shiftd(nfrst,nlast,kk,iarr)                                        */
/*                                                                          */
/*  THIS ROUTINE SHIFTS A SET OF CONTIGUOUS ELEMENTS OF AN                  */
/*  INTEGER ARRAY KK POSITIONS DOWNWARD (UPWARD IF KK .LT. 0).              */
/*  THE LOOPS ARE UNROLLED IN ORDER TO INCREASE EFFICIENCY.                 */
/*                                                                          */
/*  INPUT PARAMETERS - NFRST,NLAST - BOUNDS ON THE PORTION OF               */
/*                                   IARR TO BE SHIFTED.  ALL               */
/*                                   ELEMENTS BETWEEN AND                   */
/*                                   INCLUDING THE BOUNDS ARE               */
/*                                   SHIFTED UNLESS NFRST .GT.              */
/*                                   NLAST, IN WHICH CASE NO                */
/*                                   SHIFT OCCURS.                          */
/*                                                                          */
/*                              KK - NUMBER OF POSITIONS EACH               */
/*                                   ELEMENT IS TO BE SHIFTED.              */
/*                                   IF KK .LT. 0 SHIFT UP.                 */
/*                                   IF KK .GT. 0 SHIFT DOWN.               */
/*                                                                          */
/*                            IARR - INTEGER ARRAY OF LENGTH                */
/*                                   .GE. NLAST + MAX(KK,0).                */
/*                                                                          */
/*  NFRST, NLAST, AND KK ARE NOT ALTERED BY THIS ROUTINE.                   */
/*                                                                          */
/*  Ouput: iarr - SHIFTED ARRAY.                                            */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_shiftd(TDAContext *ctx, int nfrst,int nlast,int kk,int *iarr)  
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    int inc,k,nf,nl,nlp1,ns,nsl,ibak,indx,iimax;

    inc = 5;

    k = kk;                                                                   
    nf = nfrst;                                                               
    nl = nlast;                                                               
    if (nf > nl || k == 0)
        return;                                                      
    nlp1 = nl + 1;                                                            
    ns = nlp1 - nf;                                                           
    nsl = inc * (ns / inc);                                                   
    if (k < 0)
        goto L4;                                                          
                                                                              
    if (nsl <= 0)
        goto L2;                                                  

    for (i = 1; i <= nsl; i += inc) {
        ibak = nlp1 - i;                                                      
        indx = ibak + k;                                                      
      
        iarr[indx] = iarr[ibak];                                              
        iarr[indx-1] = iarr[ibak-1];                                          
        iarr[indx-2] = iarr[ibak-2];                                          
        iarr[indx-3] = iarr[ibak-3];                                          
        iarr[indx-4] = iarr[ibak-4];                                          
    }
L2:
    ibak = nlp1 - nsl;                                                        
L3:
    if (ibak <= nf)
        return;                                                             

    ibak--;                                                                   
    indx = ibak + k;                                                                

    iarr[indx] = iarr[ibak];                                                  
    goto L3;                                                                  
                                                                              
L4:
    if (nsl <= 0)
        goto L6;                                                  
    iimax = nlp1 - inc;                                                       

    for (i = nf; i <= iimax; i += inc) {
        indx = i + k;                                                         
        iarr[indx] = iarr[i];                                                 
  
        iarr[indx + 1] = iarr[i + 1];                                             
        iarr[indx + 2] = iarr[i + 2];                                             
        iarr[indx + 3] = iarr[i + 3];                                             
        iarr[indx + 4] = iarr[i + 4];                                             
    }                                                                         

L6:
    i = nsl + nf;                                                             
L7:
    if (i > nl)
        return;                                                            
    indx = i + k;                                                             
    iarr[indx] = iarr[i];                                                     
    i++;                                                                      
    goto L7;                                                                  
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_bdyadd(kk,i1,i2,iadj,iend)                                         */
/*                                                                          */
/*  THIS ROUTINE ADDS A BOUNDARY NODE TO A TRIANGULATION OF                 */
/*  A SET OF KK-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND                 */
/*  ARE UPDATED WITH THE INSERTION OF NODE KK.                              */
/*                                                                          */
/*  INPUT PARAMETERS -   KK - INDEX OF AN EXTERIOR NODE TO BE               */
/*                            ADDED.  KK .GE. 4.                            */
/*                                                                          */
/*                       I1 - FIRST (RIGHTMOST AS VIEWED FROM               */
/*                            KK) BOUNDARY NODE IN THE MESH                 */
/*                            WHICH IS VISIBLE FROM KK - THE                */
/*                            LINE SEGMENT KK-I1 INTERSECTS                 */
/*                            NO ARCS.                                      */
/*                                                                          */
/*                       I2 - LAST (LEFTMOST) BOUNDARY NODE                 */
/*                            WHICH IS VISIBLE FROM KK.                     */
/*                                                                          */
/*                     IADJ - SET OF ADJACENCY LISTS OF NODES               */
/*                            IN THE MESH.                                  */
/*                                                                          */
/*                     IEND - POINTERS TO THE ENDS OF                       */
/*                            ADJACENCY LISTS IN IADJ FOR                   */
/*                            EACH NODE IN THE MESH.                        */
/*                                                                          */
/*  IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN                 */
/*  THE VERTICES I1 AND I2.  I1 AND I2 MAY BE DETERMINED BY                 */
/*  TRFIND.                                                                 */
/*                                                                          */
/*  KK, I1, AND I2 ARE NOT ALTERED BY THIS ROUTINE.                         */
/*                                                                          */
/*  Output: iadj and iend: UPDATED WITH THE ADDITION                        */
/*                                  OF NODE KK AS THE LAST                  */
/*                                  ENTRY.  NODE KK WILL BE                 */
/*                                  COnnECTED TO I1, I2, AND                */
/*                                  ALL BOUNDARY NODES BETWEEN              */
/*                                  THEM.  NO OPTIMIZATION OF               */
/*                                  THE MESH IS PERFORMED.                  */
/*                                                                          */
/*  Functions called: intp_shiftd()                                         */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_bdyadd(TDAContext *ctx, int kk,int i1,int i2,int *iadj,int *iend)
{
    register int i;
    int k,km1,nright,nleft,nf,nl,n1,n2,iimin,iimax,kend,next,indx;

    k = kk;                                                                   
    km1 = k - 1;                                                              
    nright = i1;                                                              
    nleft = i2;                                                               
                                                                              
    nl = iend[km1];                                                           
    n1 = 1;                                                                   
    if (nleft != 1)
        n1 = iend[nleft - 1] + 1;                                 
    n2 = iend[nright];                                                        
    nf = imax(ctx, n1,n2);                                                         
                                                                              
    intp_shiftd(ctx, nf,nl,2,iadj);                                               

    iadj[nf + 1] = k;                                                         
    iimin = imax(ctx, nright,nleft);                                               

    for (i = iimin; i <= km1; ++i)
        iend[i] += 2;
                                                                              
    kend = nl + 3;                                                            
    nl = nf - 1;                                                              
    nf = imin(ctx, n1,n2);                                                         
    intp_shiftd(ctx, nf,nl,1,iadj);                                               

    iadj[nf] = k;                                                             
    iimax = iimin - 1;                                                        
    iimin = imin(ctx, nright,nleft);                                               

    for (i = iimin; i <= iimax; ++i)
        iend[i] += 1;
                                                                              
    iadj[kend] = nright;                                                      
    indx = iend[nright] - 2;                                                  

L3:
    next = iadj[indx];                                                        
    if (next == nleft)
        goto L4;                                             
                                                                              
    kend++;                                                                   
    iadj[kend] = next;                                                        
    indx = iend[next];                                                        
    iadj[indx] = k;                                                           
    indx--;                                                                   
    goto L3;                                                                  
                                                                              
L4:
    iadj[kend + 1] = nleft;                                                   
    kend += 2;                                                                
    iadj[kend] = 0;                                                           
    iend[k] = kend;                                                           
    return;                                                                   
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_intadd(kk,i1,i2,i3,iadj,iend)                                      */
/*                                                                          */
/*  THIS ROUTINE ADDS AN INTERIOR NODE TO A TRIANGULATION                   */
/*  OF A SET OF kk-1 POINTS ON THE UNIT SPHERE.  IADJ AND IEND              */
/*  ARE UPDATED WITH THE INSERTION OF NODE kk IN THE TRIANGLE               */
/*  WHOSE VERTICES ARE I1, I2, AND I3.                                      */
/*                                                                          */
/*  INPUT PARAMETERS -        kk - INDEX OF NODE TO BE                      */
/*                                 INSERTED.  kk .GE. 4.                    */
/*                                                                          */
/*                      I1,I2,I3 - INDICES OF THE VERTICES OF               */
/*                                 A TRIANGLE CONTAINING NODE               */
/*                                 kk - IN COUNTERCLOCKWISE                 */
/*                                 ORDER.                                   */
/*                                                                          */
/*                          IADJ - SET OF ADJACENCY LISTS                   */
/*                                 OF NODES IN THE MESH.                    */
/*                                                                          */
/*                          IEND - POINTERS TO THE ENDS OF                  */
/*                                 ADJACENCY LISTS IN IADJ FOR              */
/*                                 EACH NODE IN THE MESH.                   */
/*                                                                          */
/*  IADJ AND IEND MAY BE CREATED BY TRMESH AND MUST CONTAIN                 */
/*  THE VERTICES I1, I2, AND I3.  I1,I2,I3 MAY BE DETERMINED                */
/*  BY TRFIND.                                                              */
/*                                                                          */
/*  kk, I1, I2, AND I3 ARE NOT ALTERED BY THIS ROUTINE.                     */
/*                                                                          */
/*  Output: iadj and iend: UPDATED WITH THE ADDITION                        */
/*                                  OF NODE kk AS THE LAST                  */
/*                                  ENTRY.  NODE kk WILL BE                 */
/*                                  COnnECTED TO NODES I1, I2,              */
/*                                  AND I3.  NO OPTIMIZATION                */
/*                                  OF THE MESH IS PERFORMED.               */
/*                                                                          */
/*  Functions called: intp_shiftd()                                         */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_intadd(TDAContext *ctx, int kk,int i1,int i2,int i3,int *iadj,int *iend)
{
    register int i;
    int k,km1,ip1,ip2,ip3,indx,nf,nl,n1,n2,iimin,iimax,itemp;
    int n[4],nft[4];
 
    k = kk;                                                                   
    n[1] = i1;                                                                
    n[2] = i2;                                                                
    n[3] = i3;                                                                
                                                                              
    for (i = 1; i <= 3; ++i) {

        n1 = n[i];                                                            
        indx = (i % 3) + 1;                                                  
        n2 = n[indx];                                                         
        indx = iend[n1] + 1;                                                  
                                                                              
L1:     indx--;                                                               
        if (iadj[indx] != n2)
            goto L1;                                        
        nft[i] = indx + 1;                                                    
    }
    ip1 = 1;                                                                  
    ip2 = 2;                                                                  
    ip3 = 3;                                                                  
    if (n[2] <= n[1])
        goto L3;                                            
    ip1 = 2;                                                                  
    ip2 = 1;                                                                  
L3:
    if (n[3] <= n[ip1])
        goto L4;                                          
    ip3 = ip1;                                                                
    ip1 = 3;                                                                  

L4:
    if (n[ip3] <= n[ip2])
        goto L5;                                       
    itemp = ip2;                                                              
    ip2 = ip3;                                                                
    ip3 = itemp;                                                              
                                                                              
L5:
    km1 = k - 1;                                                              
    nl = iend[km1];                                                           
    nf = nft[ip1];                                                            
    if (nf <= nl)
        intp_shiftd(ctx, nf,nl,3,iadj);                               

    iadj[nf + 2] = k;                                                         
    iimin = n[ip1];                                                           
    iimax = km1;                                                              

    for (i = iimin; i <= iimax; ++i) 
        iend[i] += 3;
                                                                              
    nl = nf - 1;                                                              
    nf = nft[ip2];                                                            
    intp_shiftd(ctx, nf,nl,2,iadj);                                               

    iadj[nf + 1] = k;                                                         
    iimax = iimin - 1;                                                        
    iimin = n[ip2];                                                           

    for (i = iimin; i <= iimax; ++i) 
        iend[i] += 2;
                                                                              
    nl = nf - 1;                                                              
    nf = nft[ip3];                                                            
    intp_shiftd(ctx, nf,nl,1,iadj);                                               

    iadj[nf] = k;                                                             
    iimax = iimin - 1;                                                        
    iimin = n[ip3];                                                           

    for (i = iimin; i <= iimax; ++i) 
        iend[i] += 1;
                                                                              
    indx = iend[km1];                                                         
    iend[k] = indx + 3;                                                       

    for (i = 1; i <= 3; ++i) {
        indx++;                                                               
        iadj[indx] = n[i];                                                    
    }
    return;                                                                   
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_swap(nin1,nin2,nout1,nout2,iadj,iend)                              */
/*                                                                          */
/*  THIS SUBROUTINE SWAPS THE DIAGONALS IN A CONVEX QUADRILATERAL.          */
/*                                                                          */
/*  INPUT PARAMETERS -  NIN1,NIN2,NOUT1,NOUT2 - NODAL INDICES               */
/*                             OF A PAIR OF ADJACENT TRIANGLES              */
/*                             WHICH FORM A CONVEX QUADRILAT-               */
/*                             ERAL.  NOUT1 AND NOUT2 ARE CON-              */
/*                             NECTED BY AN ARC WHICH IS TO BE              */
/*                             REPLACED BY THE ARC NIN1-NIN2.               */
/*                             (NIN1,NOUT1,NOUT2) MUST BE TRI-              */
/*                             ANGLE VERTICES IN COUNTERCLOCK-              */
/*                             WISE ORDER.                                  */
/*                                                                          */
/*  THE ABOVE PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.                   */
/*                                                                          */
/*                 iadj,IEND - TRIANGULATION DATA STRUCTURE                 */
/*                             (SEE SUBROUTINE TRMESH).                     */
/*                                                                          */
/*  Output: iadj, iend   updated with the arc replacement.                  */
/*                                                                          */
/*  Functions called: intp_index(), intp_shiftd()                           */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_swap(TDAContext *ctx, int nin1,int nin2,int nout1,int nout2,int *iadj,int *iend) 
{
    register int i;
    int in[3],io[3],ip1,ip2,j,k,nf,nl,iimin,iimax;

    in[1] = nin1;                                                             
    in[2] = nin2;                                                             
    io[1] = nout1;                                                            
    io[2] = nout2;                                                            
    ip1 = 1;                                                                  
                                                                              
    if (in[1] < in[2])
        goto L1;                                            
    in[1] = in[2];                                                            
    in[2] = nin1;                                                             
    ip1 = 2;                                                                  

L1:
    if (io[1] < io[2])
        goto L2;                                            
    io[1] = io[2];                                                            
    io[2] = nout1;                                                            
    ip1 = 3 - ip1;                                                            

L2:
    ip2 = 3 - ip1;                                                            
    if (io[2] < in[1])
        goto L8;                                            
    if (in[2] < io[1])
        goto L12;                                           
                                                                              
    for (j = 1; j <= 2; ++j) {

        k = 3 - j;                                                            
        if (in[j] > io[j])
            goto L4;                                          
                                                                              
        nf = 1 + intp_index(ctx, in[j],io[ip1],iadj,iend);                              
        nl = -1 + intp_index(ctx, io[j],io[k],iadj,iend);                               
        if (nf <= nl)
            intp_shiftd(ctx, nf,nl,1,iadj);                            

        iadj[nf] = in[k];                                                     
        iimin = in[j];                                                        
        iimax = io[j] - 1;                                                    

        for (i = iimin; i <= iimax; ++i)  
            iend[i] += 1;                                                     
        goto L6;                                                              
                                                                              
L4:     nf = 1 + intp_index(ctx, io[j],io[k],iadj,iend);                                
        nl = -1 + intp_index(ctx, in[j],io[ip2],iadj,iend);                             
        if (nf <= nl)
            intp_shiftd(ctx, nf,nl,-1,iadj);                           

        iadj[nl] = in[k];                                                     
        iimin = io[j];                                                        
        iimax = in[j] - 1;                                                    

        for (i = iimin; i <= iimax; ++i)
            iend[i] -= 1;
L6:
        ip1 = ip2;                                                            
        ip2 = 3 - ip1;                                                        
    }
    return;                                                                   
                                                                              
L8:
    nf = 1 + intp_index(ctx, io[1],io[2],iadj,iend);                                    
    nl = -1 + intp_index(ctx, io[2],io[1],iadj,iend);                                   
    if (nf <= nl)
        intp_shiftd(ctx, nf,nl,-1,iadj);                             

    iimin = io[1];                                                            
    iimax = io[2] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] -= 1;
                                                                              
    nf = nl + 2;                                                              
    nl = -1 + intp_index(ctx, in[1],io[ip2],iadj,iend);                                 
    if (nf <= nl)
        intp_shiftd(ctx, nf,nl,-2,iadj);                             

    iadj[nl - 1] = in[2];                                                     
    iimin = io[2];                                                            
    iimax = in[1] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] -= 2;
                                                                              
    nf = nl + 1;                                                              
    nl = -1 + intp_index(ctx, in[2],io[ip1],iadj,iend);                                 
    intp_shiftd(ctx, nf,nl,-1,iadj);                                             

    iadj[nl] = in[1];                                                         
    iimin = in[1];                                                            
    iimax = in[2] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] -= 1;
    return;                                                                   
                                                                              
L12:
    nf = 1 + intp_index(ctx, io[1],io[2],iadj,iend);                                    
    nl = -1 + intp_index(ctx, io[2],io[1],iadj,iend);                                   
    if (nf <= nl)
        intp_shiftd(ctx, nf,nl,1,iadj);                              

    iimin = io[1];                                                            
    iimax = io[2] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] += 1;
                                                                              
    nl = nf - 2;                                                              
    nf = 1 + intp_index(ctx, in[2],io[ip2],iadj,iend);                                  
    if (nf <= nl)
        intp_shiftd(ctx, nf,nl,2,iadj);                              

    iadj[nf+1] = in[1];                                                       
    iimin = in[2];                                                            
    iimax = io[1] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] += 2;
                                                                              
    nl = nf - 1;                                                              
    nf = 1 + intp_index(ctx, in[1],io[ip1],iadj,iend);                                  
    intp_shiftd(ctx, nf,nl,1,iadj);                                              

    iadj[nf] = in[2];                                                         
    iimin = in[1];                                                            
    iimax = in[2] - 1;                                                        

    for (i = iimin; i <= iimax; ++i)
        iend[i] += 1;                                                         
    return;                                                                   
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_trfind(nst,px,py,x,y,iadj,iend,i1,i2,i3)                           */
/*                                                                          */
/*  THIS ROUTINE LOCATES A POINT P IN A THIESSEN TRIANGU-                   */
/*  LATION, RETURNING THE VERTEX INDICES OF A TRIANGLE WHICH                */
/*  CONTAINS P.  TRFIND IS PART OF AN INTERPOLATION PACKAGE                 */
/*  WHICH PROVIDES SUBROUTINES FOR CREATING THE MESH.                       */
/*                                                                          */
/*  INPUT PARAMETERS -    NST - INDEX OF NODE AT WHICH TRFIND               */
/*                              BEGINS SEARCH.  SEARCH TIME                 */
/*                              DEPENDS ON THE PROXIMITY OF                 */
/*                              NST TO P.                                   */
/*                                                                          */
/*                      PX,PY - X AND Y-COORDINATES OF THE                  */
/*                              POINT TO BE LOCATED.                        */
/*                                                                          */
/*                        X,Y - VECTORS OF COORDINATES OF                   */
/*                              NODES IN THE MESH.  (X(I),Y(I))             */
/*                              DEFINES NODE I FOR I = 1,...,N              */
/*                              WHERE N .GE. 3.                             */
/*                                                                          */
/*                       iadj - SET OF ADJACENCY LISTS OF                   */
/*                              NODES IN THE MESH.                          */
/*                                                                          */
/*                       iend - POINTERS TO THE ENDS OF                     */
/*                              ADJACENCY LISTS IN iadj FOR                 */
/*                              EACH NODE IN THE MESH.                      */
/*                                                                          */
/*  iadj AND iend MAY BE CREATED BY TRMESH.                                 */
/*                                                                          */
/*  OUTPUT PARAMETERS - I1,I2,I3 - VERTEX INDICES IN COUNTER-               */
/*                               CLOCKWISE ORDER - VERTICES                 */
/*                               OF A TRIANGLE CONTAINING P                 */
/*                               IF P IS AN INTERIOR NODE.                  */
/*                               IF P IS OUTSIDE OF THE                     */
/*                               BOUNDARY OF THE MESH, I1                   */
/*                               AND I2 ARE THE FIRST (RIGHT                */
/*                               -MOST) AND LAST (LEFTMOST)                 */
/*                               NODES WHICH ARE VISIBLE                    */
/*                               FROM P, AND I3 = 0.  IF P                  */
/*                               AND ALL OF THE NODES LIE ON                */
/*                               A SINGLE LINE THEN I1 = I2                 */
/*                               = I3 = 0.                                  */
/*                                                                          */
/*  Functions called: intp_left()                                           */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_trfind(TDAContext *ctx, int nst,double px,double py,double *x,double *y,int *iadj, int *iend,int *i1,int *i2,int *i3)
{
    int n0,n1,n2,n3,n4,indx,ind,nf,nl,next;
    double xp,yp;

    xp = px;                                                                  
    yp = py;                                                                  
                                                                              
    n0 = imax(ctx, nst,1);                                                         
L1: indx = iend[n0];                                                          
    nl = iadj[indx];                                                          
    indx = 1;                                                                 
    if (n0 != 1)
        indx = iend[n0 - 1] + 1;                                     
    nf = iadj[indx];                                                          
    n1 = nf;                                                                  
    if (nl != 0)
        goto L3;                                                   
                                                                              
    ind = iend[n0] - 1;                                                       
    nl = iadj[ind];                                                           
    if (intp_left(ctx, x[n0],y[n0],x[nf],y[nf],xp,yp))
        goto L2;                       
                                                                              
    nl = n0;                                                                  
    goto L16;                                                                 

L2: if (intp_left(ctx, x[nl],y[nl],x[n0],y[n0],xp,yp))
        goto L4;                       
                                                                              
    *i1 = n0;                                                                 
    goto L18;                                                                 
                                                                              
L3: if (intp_left(ctx, x[n0],y[n0],x[n1],y[n1],xp,yp))
        goto L4;                       

    indx++;                                                                   
    n1 = iadj[indx];                                                          
    if (n1 == nl)
        goto L7;                                                            
    goto L3;                                                                  
                                                                              
L4: indx++;                                                                   
    n2 = iadj[indx];                                                          
    if (intp_left(ctx, x[n0],y[n0],x[n2],y[n2],xp,yp) == 0)                      
        goto L8;                                                              

    n1 = n2;                                                                  
    if (n1 != nl)
        goto L4;                                                  
    if (intp_left(ctx, x[n0],y[n0],x[nf],y[nf],xp,yp) == 0)                      
        goto L7;                                                              
    if (xp == x[n0] && yp == y[n0])
        goto L6;                           
                                                                              
L5: if (intp_left(ctx, x[n1],y[n1],x[n0],y[n0],xp,yp) == 0)                      
        goto L6;                                                              
    if (n1 == nf)
        goto L20;                                                 
    indx--;                                                                   
    n1 = iadj[indx];                                                          
    goto L5;                                                                  
                                                                              
L6: n0 = n1;                                                                  
    goto L1;                                                                  
                                                                              
L7: n2 = nf;                                                                  
                                                                              
L8: n3 = n0;                                                                  
L9: if (intp_left(ctx, x[n1],y[n1],x[n2],y[n2],xp,yp))
        goto L13;                      
                                                                              
    indx = iend[n2];                                                          
    if (iadj[indx] != n1)
        goto L10;                                         
                                                                              
    indx = 1;                                                                 
    if (n2 != 1)
        indx = iend[n2 - 1] + 1;                                     
    n4 = iadj[indx];                                                          
    goto L11;                                                                 
                                                                              
L10:
    indx--;                                                                   
    if (iadj[indx] != n1)
        goto L10;                                         
    n4 = iadj[indx + 1];                                                      
    if (n4 != 0)
        goto L11;                                                  
                                                                              
    nf = n2;                                                                  
    nl = n1;                                                                  
    goto L16;                                                                 
                                                                              
L11:
    if (intp_left(ctx, x[n0],y[n0],x[n4],y[n4],xp,yp))
        goto L12;                      
    n3 = n2;                                                                  
    n2 = n4;                                                                  
    goto L9;                                                                  

L12:
    n3 = n1;                                                                  
    n1 = n4;                                                                  
    goto L9;                                                                  
                                                                              
L13:
    indx = iend[n1];                                                          
    if (iadj[indx] != 0)
        goto L15;                                          
                                                                              
    if (n3 != iadj[indx - 1])
        goto L14;                                       
                                                                              
    if (intp_left(ctx, x[n1],y[n1],x[n3],y[n3],xp,yp) == 0)
        goto L14;                                                             
                                                                              
    *i1 = n1;                                                                 
    *i2 = n3;                                                                 
    *i3 = 0;                                                                  
    return;                                                                  
                                                                              
L14:
    indx = 1;                                                                 
    if (n1 != 1)
        indx = iend[n1 - 1] + 1;                                     
    if (n2 != iadj[indx])
        goto L15;                                         
                                                                              
    if (intp_left(ctx, x[n2],y[n2],x[n1],y[n1],xp,yp) == 0)                      
        goto L15;                                                             
                                                                              
    *i1 = n2;                                                                 
    *i2 = n1;                                                                 
    *i3 = 0;                                                                  
    return;                                                                   
                                                                              
L15:
    *i1 = n1;                                                                 
    *i2 = n2;                                                                 
    *i3 = n3;                                                                 
    return;                                                                   
                                                                              
L16:
    indx = 1;                                                                 
    if (nf != 1)
        indx = iend[nf - 1] + 1;                                     
    next = iadj[indx];                                                        

    if (intp_left(ctx, x[nf],y[nf],x[next],y[next],xp,yp))                            
        goto L17;                                                             
    nf = next;                                                                
    goto L16;                                                                 
                                                                              
L17:
    *i1 = nf;                                                                 
                                                                              
L18:
    indx = iend[nl] - 1;                                                      
    next = iadj[indx];                                                        

    if (intp_left(ctx, x[next],y[next],x[nl],y[nl],xp,yp))                            
        goto L19;                                                             
    nl = next;                                                                
    goto L18;                                                                 
                                                                              
L19:
    *i2 = nl;                                                                 
    *i3 = 0;                                                                  
    return;                                                                   
                                                                              
L20:
    *i1 = *i2 = *i3 = 0;                                                      
    return;                                                                   
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_swptst(in1,in2,io1,io2,x,y)                                        */
/*                                                                          */
/*  THIS FUNCTION DECIDES WHETHER OR NOT TO REPLACE A                       */
/*  DIAGONAL ARC IN A QUADRILATERAL WITH THE OTHER DIAGONAL.                */
/*  THE DETERMINATION IS BASED ON THE SIZES OF THE ANGLES                   */
/*  CONTAINED IN THE 2 TRIANGLES DEFINED BY THE DIAGONAL.                   */
/*  THE DIAGONAL IS CHOSEN TO MAXIMIZE THE SMALLEST OF THE                  */
/*  SIX ANGLES OVER THE TWO PAIRS OF TRIANGLES.                             */
/*                                                                          */
/*  INPUT PARAMETERS -  In1,In2,IO1,IO2 - NODE INDICES OF THE               */
/*                             FOUR POINTS DEFINING THE                     */
/*                             QUADRILATERAL.  IO1 AND IO2                  */
/*                             ARE CURRENTLY COnnECTED BY A                 */
/*                             DIAGONAL ARC.  THIS ARC                      */
/*                             SHOULD BE REPLACED BY AN ARC                 */
/*                             COnnECTING In1, In2 IF THE                   */
/*                             DECISION IS MADE TO SWAP.                    */
/*                             In1,IO1,IO2 MUST BE IN                       */
/*                             COUNTERCLOCKWISE ORDER.                      */
/*                                                                          */
/*                       X,Y - VECTORS OF NODAL COORDINATES.                */
/*                             (X(I),Y(I)) ARE THE COORD-                   */
/*                             INATES OF NODE I FOR I = In1,                */
/*                             In2, IO1, OR IO2.                            */
/*                                                                          */
/*  Return 1 if the arc connecting io1 and io2 is to be replaced,           */
/*  otherwise return 0.                                                     */
/*                                                                          */
/*  Functions called: none                                                  */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */
 
int intp_swptst(TDAContext *ctx, int in1,int in2,int io1,int io2,double *x,double *y)             
{
    (void)ctx;        /* unused: the signature is shared */
    double dx11,dx12,dx22,dx21,dy11,dy12,dy22,dy21,sin1,sin2,cos1,cos2,sin12; 
                                                                              
    dx11 = x[io1] - x[in1];                                                   
    dx12 = x[io2] - x[in1];                                                   
    dx22 = x[io2] - x[in2];                                                   
    dx21 = x[io1] - x[in2];                                                   
                                                                              
    dy11 = y[io1] - y[in1];                                                   
    dy12 = y[io2] - y[in1];                                                   
    dy22 = y[io2] - y[in2];                                                   
    dy21 = y[io1] - y[in2];                                                   
                                                                              
    cos1 = dx11 * dx12 + dy11 * dy12;                                         
    cos2 = dx22 * dx21 + dy22 * dy21;                                         
                                                                              
    if (cos1 >= 0.0 && cos2 >= 0.0)
        return(0);                 
                                   
    if (cos1 < 0.0 && cos2 < 0.0)
        return(1);               
                                 
    sin1 = dx11 * dy12 - dx12 * dy11;                                         
    sin2 = dx22 * dy21 - dx21 * dy22;                                         
    sin12 = sin1 * cos2 + cos1 * sin2;                                        
    if (sin12 >= 0.0)
        return(0);   
    return(1);       
}                                                                             

/* -##--------------------------------------------------------------------- */
/*  intp_index(nvertx,nabor,iadj,iend)                                      */
/*                                                                          */
/*  THIS FUNCTION RETURNS THE INDEX OF NABOR IN THE                         */
/*  ADJACENCY LIST FOR NVERTX.                                              */
/*                                                                          */
/*  INPUT PARAMETERS - NVERTX - NODE WHOSE ADJACENCY LIST IS                */
/*                              TO BE SEARCHED.                             */
/*                     NABOR - NODE WHOSE INDEX IS TO BE                    */
/*                             RETURNED.  NABOR MUST BE                     */
/*                             COnnECTED TO NVERTX.                         */
/*                     iadj - SET OF ADJACENCY LISTS.                       */
/*                     iend - POINTERS TO THE ENDS OF                       */
/*                            ADJACENCY LISTS IN iadj.                      */
/*                                                                          */
/*  Return: index - iadj[index] = nabor.                                    */
/*                                                                          */
/*  Functions called: none                                                  */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_index(TDAContext *ctx, int nvertx,int nabor,int *iadj,int *iend)             
{
    (void)ctx;        /* unused: the signature is shared */
    int indx;
                                                                              
    indx = iend[nvertx] + 1;                                                  
    while (1) {
        if (iadj[--indx] == nabor)
            break;
    }
    return(indx);
}

/* -##--------------------------------------------------------------------- */
/*  intp_left(x1,y1,x2,y2,x0,y0)                                            */
/*                                                                          */
/*  This function computes the sign of a cross product (z-component).       */
/*  Return 1 if (x0,y0) is on or to the left of the vector from (x1,y1) to  */
/*  (x2,y2), otherwise return 0.                                            */

int intp_left(TDAContext *ctx, double x1,double y1,double x2,double y2,double x0,double y0)
{
    (void)ctx;        /* unused: the signature is shared */
    if ((x2 - x1) * (y0 - y1) >= (x0 - x1) * (y2 - y1))                       
        return(1);
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  intp_trmtst(n,x,y,iadj,iend,tol)                                        */
/*                                                                          */
/*  THIS ROUTINE TESTS THE VALIDITY OF THE DATA STRUCTURE                   */
/*  REPRESENTING A THIESSEN TRIANGULATION CREATED BY SUBROU-                */
/*  TINE TRMESH.  THE FOLLOWING PROPERTIES ARE TESTED --                    */
/*  1)  IEND(1) .GE. 3 AND IEND(K) .GE. IEND(K-1)+3 FOR K =                 */
/*      2,...,N (EACH NODE HAS AT LEAST THREE NEIGHBORS).                   */
/*  2)  0 .LE. IADJ(K) .LE. N FOR K = 1,...,IEND(N) (IADJ                   */
/*      ENTRIES ARE NODAL INDICES OR ZEROS REPRESENTING THE                 */
/*      BOUNDARY).                                                          */
/*  3)  NB .GE. 3, NT = 2N-NB-2, AND NA = 3N-NB-3 WHERE NB,                 */
/*      NT, AND NA ARE THE NUMBERS OF BOUNDARY NODES, TRI-                  */
/*      ANGLES, AND ARCS, RESPECTIVELY.                                     */
/*  4)  EACH CIRCUMCIRCLE DEFINED BY THE VERTICES OF A TRI-                 */
/*      ANGLE CONTAINS NO NODES IN ITS INTERIOR.  THIS PROP-                */
/*      ERTY DISTINGUISHES A THIESSEN TRIANGULATION FROM AN                 */
/*      ARBITRARY TRIANGULATION OF THE NODES.                               */
/*                                                                          */
/*  NOTE THAT NO TEST IS MADE FOR THE PROPERTY THAT A TRIANGU-              */
/*  LATION COVERS THE CONVEX HULL OF THE NODES, AND THUS A                  */
/*  TEST ON A DATA STRUCTURE ALTERED BY SUBROUTINE DELETE                   */
/*  SHOULD NOT RESULT IN AN ERROR.                                          */
/*                                                                          */
/*  INPUT PARAMETERS -                                                      */
/*      N - NUMBER OF NODES.  N .GE. 3.                                     */
/*      X,Y - NODAL COORDINATES.                                            */
/*      IADJ,IEND - TRIANGULATION DATA STRUCTURE.  SEE SUB-                 */
/*                  ROUTINE TRMESH.                                         */
/*      TOL - NONNEGATIVE TOLERANCE TO ALLOW FOR FLOATING-                  */
/*            POINT ERRORS IN THE CIRCUMCIRCLE TEST.  AN                    */
/*            ERROR SITUATION IS DEFINED AS R**2 - D**2 .GT.                */
/*            TOL WHERE R IS THE RADIUS OF A CIRCUMCIRCLE                   */
/*            AND D IS THE DISTANCE FROM THE CIRCUMCENTER                   */
/*            TO THE NEAREST NODE.  A REASONABLE VALUE                      */
/*            FOR TOL IS 10*EPS WHERE EPS IS THE MACHINE                    */
/*            PRECISION.  THE TEST IS EFFECTIVELY BYPASSED                  */
/*            BY MAKING TOL LARGER THAN THE DIAMETER OF THE                 */
/*            CONVEX HULL OF THE NODES.                                     */
/*                                                                          */
/*  Return:   -1 IF ONE OR MORE NULL TRIANGLES (AREA = 0)                   */
/*               ARE PRESENT BUT NO (OTHER) ERRORS                          */
/*               WERE ENCOUNTERED.  A NULL TRIANGLE IS                      */
/*               AN ERROR ONLY IF IT OCCURS IN THE                          */
/*               THE INTERIOR.                                              */
/*             0 IF NO ERRORS OR NULL TRIANGLES WERE ENCOUNTERED            */
/*             1 IF N .LT. 3 OR TOL .LT. 0.                                 */
/*             2 IF AN IEND OR IADJ ENTRY IS OUT OF RANGE                   */
/*             3 IF THE TRIANGULATION PARAMETERS (NB, NT, and NA)           */
/*               ARE INCONSISTENT.                                          */
/*             4 IF A TRIANGLE CONTAINS A NODE INTERIOR                     */
/*               TO ITS CIRCUMCIRCLE.                                       */
/*                                                                          */
/*  Functions called: intp_circum().                                        */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_trmtst(TDAContext *ctx, int n,double *x,double *y,int *iadj,int *iend,double tol)  
{
    register int i = 0;
    int nn = 0,indx = 0,indf = 0,indl = 0,na = 0,nb = 0,nt = 0,nfail = 0,n1 = 0,n2 = 0,n3 = 0,ierr = 0;
    double rtol = 0.0,r = 0.0,tmp1 = 0.0,tmp2 = 0.0,cx = 0.0,cy = 0.0;

    nn = n;
    rtol = tol;
    if (nn < 3 || rtol < 0.0)
        return(1);           

    nb = nt = nfail = 0;
 
    /* LOOP ON TRIANGLES (N1,N2,N3) SUCH THAT N2 AND N3 INDEX
       ADJACENT NEIGHBORS OF N1 AND ARE BOTH LARGER THAN N1
      (TRIANGLES ARE ASSOCIATED WITH THEIR SMALLEST INDEX). */
 
    indf = 1;
    for (n1 = 1; n1 <= nn; ++n1) {

        indl = iend[n1];
        if (indl < indf + 2)
            return(2); 

        if (iadj[indl] == 0)
            nb++;                              

        /* LOOP ON NEIGHBORS OF N1 */
 
        for (indx = indf; indx <= indl; ++indx) {

            n2 = iadj[indx];
            if (n2 < 0 || n2 > nn || (indx < indl && n2 == 0))
                return(2);  

            if (indx < indl)
                n3 = iadj[indx + 1];
        
            else if (indx == indl)
                n3 = iadj[indf];

            if (n2 < n1 || n3 < n1)  
                continue;
                                     
            nt++;             
 
            /* COMPUTE THE COORDINATES OF THE CIRCUMCENTER OF (N1,N2,N3) */
  
            ierr = intp_circum(ctx, x[n1],x[n2],x[n3],y[n1],y[n2],y[n3],&cx,&cy);
            if (ierr)   
                continue;

            /* TEST FOR NODES WITHIN THE CIRCUMCIRCLE. */
 
            tmp1 = cx - x[n1];
            tmp2 = cy - y[n1];           
            r = tmp1 * tmp1 + tmp2 * tmp2 - rtol;
            if (r <= 0.0)
                continue;           

            for (i = 1; i <= nn; ++i) {
                if (i == n1 || i == n2 || i == n3)
                    continue;                                  
                tmp1 = cx - x[i];
                tmp2 = cy - y[i];                           
                if (tmp1 * tmp1 + tmp2 * tmp2 < r) {
                    nfail++;                          /* node i is interior */
                    break;
                }               
            }
        }
        indf = indl + 1;
    }

    /* CHECK PARAMETERS FOR CONSISTENCY AND TEST FOR NFAIL = 0 */
 
    na = (iend[nn] - nb) / 2;
    if (nb < 3 || nt != (2 * nn - nb - 2) || na != (3 * nn - nb - 3))
        return(3);                    
                               
    if (nfail != 0)
        return(4);
 
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  intp_trmtst_emsg(r)   Report error message.                             */

void intp_trmtst_emsg(TDAContext *ctx, int r)
{
    printf1(ctx, "Error %d in triangulation data structure.\n",r);
    if (r == 1)  
        printf1(ctx, "Unvalid number of points or tolerance.\n");
    else if (r == 4)
        printf1(ctx, "A triangle contains a node interior to its circumcircle.\n");
}

/* -##--------------------------------------------------------------------- */
/*  intp_circum(x1,x2,x3,y1,y2,y3,cx,cy)                                    */
/*                                                                          */
/*  THIS SUBROUTINE COMPUTES THE COORDINATES OF THE CENTER                  */
/*  OF A CIRCLE DEFINED BY THREE POINTS IN THE PLANE.                       */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*  X1,...,Y3 - X AND Y COORDINATES OF THREE POINTS IN THE PLANE.           */
/*                                                                          */
/*  OUTPUT PARAMETERS -                                                     */
/*      CX,CY - COORDINATES OF THE CENTER OF THE CIRCLE                     */
/*                                                                          */
/*  Return: 0 if successful, -1 if error.                                   */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_circum(TDAContext *ctx, double x1,double x2,double x3, double y1,double y2,double y3,double *cx,double *cy)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double a,fx,fy,u[4],v[4],ds[4];
 
    /* SET U(K) AND V(K) TO THE X AND Y COMPONENTS OF THE EDGE
       OPPOSITE VERTEX K, TREATING THE POINTS AS VERTICES OF A TRIANGLE. */
  
    u[1] = x3 - x2;
    u[2] = x1 - x3;
    u[3] = x2 - x1;
    v[1] = y3 - y2;
    v[2] = y1 - y3;
    v[3] = y2 - y1;
 
    /* SET A TO TWICE THE SIGNED AREA OF THE TRIANGLE.  A .GT. 0
       IFF (X3,Y3) IS STRICTLY TO THE LEFT OF THE EDGE FROM (X1,Y1) TO (X2,Y2) */
 
    a = u[1] * v[2] - u[2] * v[1];
    if (a == 0.0)
        return(-1);
 
    /* SET DS(K) TO THE SQUARED DISTANCE FROM THE ORIGIN TO VERTEX K. */
 
    ds[1] = x1 * x1 + y1 * y1;
    ds[2] = x2 * x2 + y2 * y2;
    ds[3] = x3 * x3 + y3 * y3;
 
    /* COMPUTE FACTORS OF CX AND CY. */
 
    fx = fy = 0.0;
    for (i = 1; i <= 3; ++i) {
        fx -= ds[i] * v[i];
        fy += ds[i] * u[i];
    }
    *cx = fx / 2.0 / a;
    *cy = fy / 2.0 / a;
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  intp_bnodes(n,iadj,iend,nb,na,nt,nodes)                                 */
/*                                                                          */
/*  GIVEN A TRIANGULATION OF N POINTS IN THE PLANE, THIS                    */
/*  ROUTINE RETURNS A VECTOR CONTAINING THE INDICES, IN                     */
/*  COUNTERCLOCKWISE ORDER, OF THE NODES ON THE BOUNDARY OF                 */
/*  THE CONVEX HULL OF THE SET OF POINTS.                                   */
/*                                                                          */
/*  INPUT PARAMETERS                                                        */
/*      N - NUMBER OF NODES IN THE MESH.                                    */
/*      IADJ - SET OF ADJACENCY LISTS OF NODES IN THE MESH.                 */
/*      IEND - POINTERS TO THE ENDS OF ADJACENCY LISTS IN IADJ FOR          */
/*             EACH NODE IN THE MESH.                                       */
/*                                                                          */
/*      NODES - VECTOR OF LENGTH .GE. NB.(NB .LE. N).                       */
/*                                                                          */
/*  IADJ AND IEND MAY BE CREATED BY TRMESH AND ARE NOT                      */
/*  ALTERED BY THIS ROUTINE.                                                */
/*                                                                          */
/*  OUTPUT PARAMETERS                                                       */
/*      NB - NUMBER OF BOUNDARY NODES (this value is returned)              */
/*      NA,NT - NUMBER OF ARCS AND TRIANGLES, RESPECTIVELY, IN THE MESH.    */
/*      NODES - VECTOR OF NB BOUNDARY NODE INDICES RANGING FROM 1 TO N.     */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_bnodes(TDAContext *ctx, int n,int *iadj,int *iend,int *na,int *nt,int *nodes)
{
    (void)ctx;        /* unused: the signature is shared */
    int k,nst,indl,n0,indf;
 
    nst = 1;
L1: indl = iend[nst];
    if (iadj[indl] == 0)
        goto L2;
    nst++;           
    goto L1;
 
L2: nodes[1] = nst;
    k = 1;
    n0 = nst;
 
    /* TRAVERSE THE BOUNDARY IN COUNTERCLOCKWISE ORDER */
 
L3: indf = 1;
    if (n0 > 1)
        indf = iend[n0 - 1] + 1;
    n0 = iadj[indf];
    if (n0 == nst)
        goto L4;
    k++;          
    nodes[k] = n0;
    goto L3;
 
L4: *nt = 2 * n - k - 2;
    *na = *nt + n - 1;
    return(k);
}

/* -##--------------------------------------------------------------------- */
/*  intp_coords(x,y,x1,x2,x3,y1,y2,y3,r)                                    */
/*                                                                          */
/*  THIS ROUTINE COMPUTES THE THREE BARYCENTRIC COORDINATES                 */
/*  OF A POINT IN THE PLANE FOR A GIVEN TRIANGLE.                           */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      X,Y   X AND Y COORDINATES OF THE POINT                              */
/*            WHOSE BARYCENTRIC COORDINATES ARE DESIRED.                    */
/*                                                                          */
/*      X1,X2,X3,Y1,Y2,Y3   COORDINATES OF THE VERTICES OF THE TRIANGLE.    */
/*                                                                          */
/*  OUTPUT PARAMETERS                                                       */
/*      R   3-VECTOR OF BARYCENTRIC COORDINATES UNLESS.                     */
/*          NOTE THAT R(I) .LT. 0. IFF (X,Y) IS TO THE                      */
/*          RIGHT OF THE VECTOR FROM VERTEX I+1 TO VERTEX I+2               */
/*          (CYCLICAL ARITHMETIC).                                          */
/*                                                                          */
/*  Return 0 if successful, -1 if error (points are collinear).             */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_coords(TDAContext *ctx, double x,double y,double x1,double x2,double x3, double y1,double y2,double y3,double *r)
{
    (void)ctx;        /* unused: the signature is shared */
    double xp,yp,area,u[4],v[4];

    u[1] = x3 - x2;
    u[2] = x1 - x3;
    u[3] = x2 - x1;
 
    v[1] = y3 - y2;
    v[2] = y1 - y3;
    v[3] = y2 - y1;
 
    area = u[1] * v[2] - u[2] * v[1];
    if (area == 0.0)            
        return(-1);
 
    r[1] = (u[1] * (y - y2) - v[1] * (x - x2)) / area;
    xp = x - x1;
    yp = y - y1;
    r[2] = (u[2] * yp - v[2] * xp) / area;
    r[3] = (u[3] * yp - v[3] * xp) / area;
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  intp_intrc0(n,px,py,x,y,z,iadj,iend,ist,pz)                             */
/*                                                                          */
/*  GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,                  */
/*  THIS ROUTINE COMPUTES THE VALUE AT (PX,PY) OF A PIECEWISE               */
/*  LINEAR SURFACE WHICH INTERPOLATES DATA VALUES AT THE                    */
/*  VERTICES OF THE TRIANGLES.  THE SURFACE IS EXTENDED IN A                */
/*  CONTINUOUS FASHION BEYOND THE BOUNDARY OF THE TRIANGULAR                */
/*  MESH, ALLOWING EXTRAPOLATION.                                           */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      N - NUMBER OF NODES IN THE MESH. N .GE. 3.                          */
/*      PX,PY - POINT AT WHICH THE INTERPOLATED VALUE IS DESIRED.           */
/*      X,Y - VECTORS OF COORDINATES OF THE NODES IN THE MESH.              */
/*      Z - VECTOR OF DATA VALUES AT THE NODES.                             */
/*      IADJ - SET OF ADJACENCY LISTS OF NODES IN THE MESH.                 */
/*      IEND - POINTERS TO THE ENDS OF ADJACENCY LISTS IN IADJ              */
/*      IST - INDEX OF THE STARTING NODE IN THE SEARCH FOR A TRIANGLE       */
/*            CONTAINING (PX,PY).  1 .LE. IST .LE. N.  THE OUTPUT VALUE     */
/*            OF IST FROM A PREVIOUS CALL MAY BE A GOOD CHOICE.             */
/*                                                                          */
/*      IADJ AND IEND MAY BE CREATED BY TRMESH.                             */
/*                                                                          */
/*      INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED.                    */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      IST - INDEX OF ONE OF THE VERTICES OF THE TRIANGLE CONTAINING       */
/*            (PX,PY) if successful return.                                 */
/*      PZ - VALUE OF THE INTERPOLATORY SURFACE AT (PX,PY) OR ZERO          */
/*           if no success.                                                 */
/*                                                                          */
/*  Return:  0  IF NO ERRORS WERE ENCOUNTERED.                              */
/*           1  IF NO ERRORS WERE ENCOUNTERED AND EXTRAPOLATION             */
/*              WAS PERFORMED.                                              */
/*          -1 IF N OR IST IS OUT OF RANGE.                                 */
/*          -2 IF THE NODES ARE COLLINEAR.                                  */
/*                                                                          */
/*  Functions called: intp_trfind(), intp_coords().                         */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_intrc0(TDAContext *ctx, int n,double px,double py,double *x,double *y,double *z, int *iadj,int *iend,int *ist,double *pz)
{
    int i1,i2,i3,n1,n2,indx,ier;
    double xp,yp,r[4],x1,y1,x2,y2,dp,tmp1,tmp2;
 
    *pz = 0.0;
    if (n < 3 || *ist < 1 || *ist > n)  
        return(-1);
                                      
    xp = px;
    yp = py;
 
    /* FIND A TRIANGLE CONTAINING P IF P IS WITHIN THE MESH BOUNDARY */
 
    intp_trfind(ctx, *ist,xp,yp,x,y,iadj,iend,&i1,&i2,&i3);
    if (i1 == 0)             
        return(-2);

    *ist = i1;
    if (i3 == 0)
        goto L1;
 
    /* COMPUTE BARYCENTRIC COORDINATES */
 
    ier = intp_coords(ctx, xp,yp,x[i1],x[i2],x[i3],y[i1],y[i2],y[i3],r);
    if (ier != 0)
        return(-2);
                   
    *pz = r[1] * z[i1] + r[2] * z[i2] + r[3] * z[i3];
    return(0);
 
    /* P IS OUTSIDE OF THE MESH BOUNDARY.  EXTRAPOLATE TO P BY
       EXTENDING THE INTERPOLATORY SURFACE AS A CONSTANT
       BEYOND THE BOUNDARY.  THUS PZ IS THE SURFACE FUNCTION
       VALUE AT Q WHERE Q IS THE CLOSEST BOUNDARY POINT TO P. */
 
    /* DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM THE
       RIGHTMOST VISIBLE NODE I1. */
 
L1: n2 = i1;
 
    /* SET N1 TO THE LAST NONZERO NEIGHBOR OF N2 AND COMPUTE DP */
 
L2:
    indx = iend[n2] - 1;
    n1 = iadj[indx];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];
    dp = (x1 - x2) * (xp - x2) + (y1 - y2) * (yp - y2);
    if (dp <= 0.0)
        goto L3;  
                              
    if ((xp - x1) * (x2 - x1) + (yp - y1) * (y2 - y1) > 0.0)
        goto L4;
    n2 = n1;
    goto L2;
 
    /* N2 IS THE CLOSEST BOUNDARY POINT TO P */
 
L3:
    *pz = z[n2];
    return(1);
 
    /* THE CLOSEST BOUNDARY POINT TO P LIES ON N2-N1.  COMPUTE
       ITS COORDINATES WITH RESPECT TO N2-N1. */
 
L4:
    tmp1 = x2 - x1;
    tmp2 = y2 - y1;
    r[1] = dp / (tmp1 * tmp1 + tmp2 * tmp2);
    r[2] = 1.0 - r[1];
    *pz = r[1] * z[n1] + r[2] * z[n2];
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  intp_intrc1(n,px,py,x,y,z,iadj,iend,ist,pz)                             */
/*                                                                          */
/*  GIVEN A TRIANGULATION OF A SET OF POINTS IN THE PLANE,                  */
/*  THIS ROUTINE DETERMINES A PIECEWISE CUBIC FUNCTION F(X,Y)               */
/*  WHICH INTERPOLATES A SET OF DATA VALUES AND PARTIAL                     */
/*  DERIVATIVES AT THE VERTICES.  F HAS CONTINUOUS FIRST                    */
/*  DERIVATIVES OVER THE MESH AND EXTENDS BEYOND THE MESH                   */
/*  BOUNDARY ALLOWING EXTRAPOLATION.  INTERPOLATION IS EXACT                */
/*  FOR QUADRATIC DATA.  THE VALUE OF F AT (PX,PY) IS                       */
/*  RETURNED.                                                               */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      N - NUMBER OF NODES IN THE MESH. N .GE. 3.                          */
/*      PX,PY - COORDINATES OF A POINT AT WHICH F IS TO BE EVALUATED.       */
/*      X,Y - VECTORS OF COORDINATES OF THE NODES IN THE MESH.              */
/*      Z - VECTOR OF DATA VALUES AT THE NODES.                             */
/*      IADJ - SET OF ADJACENCY LISTS OF NODES IN THE MESH.                 */
/*      IEND - POINTERS TO THE ENDS OF ADJACENCY LISTS IN IADJ              */
/*                                                                          */
/*      ZXZY - 2 BY N ARRAY WHOSE COLUMNS                                   */
/*             CONTAIN ESTIMATED PARTIAL DER-                               */
/*             IVATIVES AT THE NODES (X PARTIALS IN THE FIRST ROW) IF       */
/*                           IFLAG = 1, NOT USED IF IFLAG                   */
/*                           = 0.                                           */
/*                                                                          */
/*      IST - INDEX OF THE STARTING NODE IN THE SEARCH FOR A TRIANGLE       */
/*            CONTAINING (PX,PY).  1 .LE. IST .LE. N.  THE OUTPUT VALUE     */
/*            OF IST FROM A PREVIOUS CALL MAY BE A GOOD CHOICE.             */
/*                                                                          */
/*      IADJ AND IEND MAY BE CREATED BY TRMESH AND DERIVATIVE               */
/*      ESTIMATES MAY BE COMPUTED BY GRADL OR GRADG.                        */
/*                                                                          */
/*      INPUT PARAMETERS OTHER THAN IST ARE NOT ALTERED BY THIS             */
/*      ROUTINE.                                                            */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      IST - INDEX OF ONE OF THE VERTICES OF                               */
/*            THE TRIANGLE CONTAINING (PX,PY) UNLESS IER .LT. 0.            */
/*      PZ - VALUE OF F AT (PX,PY), OR 0 IF IER .LT. 0.                     */
/*                                                                          */
/*  Return: 0 IF NO ERRORS WERE ENCOUNTERED.                                */
/*          1 IF NO ERRORS WERE ENCOUNTERED AND EXTRAPOLATION               */
/*            WAS PERFORMED.                                                */
/*         -1 IF N, IFLAG, OR IST IS OUT OF RANGE.                          */
/*         -2 IF THE NODES ARE COLLINEAR.                                   */
/*                                                                          */
/*  Functions called: intp_trfind(), intp_gradl(), intp_tval()              */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_intrc1(TDAContext *ctx, int n,double px,double py,double *x,double *y,double *z, int *iadj,int *iend,int *ist,double *pz)
{
    int nn,i1,i2,i3,n1,n2,indx,ierr;
    double xp,yp,zx1,zy1,zx2,zy2,zx3,zy3,x1,y1,x2,y2,x3,y3,z1,z2,z3;
    double dum,dp,u,v,xq,yq,r1,r2,a1,a2,b1,b2,c1,c2,f1,f2;
 
    nn = n;
    *pz = 0.0;
    if (nn < 3 || *ist < 1 || *ist > nn)
        return(-1);                     

    xp = px;
    yp = py;
 
    /* FIND A TRIANGLE CONTAINING P IF P IS WITHIN THE MESH BOUNDARY */
 
    intp_trfind(ctx, *ist,xp,yp,x,y,iadj,iend,&i1,&i2,&i3);

    if (i1 == 0)             
        return(-2);

    *ist = i1;
    if (i3 == 0)
        goto L3;
 
    /* COMPUTE DERIVATIVE ESTIMATES AT THE VERTICES */ 

    intp_gradl(ctx, nn,i1,x,y,z,iadj,iend,&zx1,&zy1);
    intp_gradl(ctx, nn,i2,x,y,z,iadj,iend,&zx2,&zy2);
    intp_gradl(ctx, nn,i3,x,y,z,iadj,iend,&zx3,&zy3);

    /* SET LOCAL PARAMETERS FOR CALL TO TVAL */ 
 
/* L2 */ x1 = x[i1];
    y1 = y[i1];
    x2 = x[i2];
    y2 = y[i2];
    x3 = x[i3];
    y3 = y[i3];
    z1 = z[i1];
    z2 = z[i2];
    z3 = z[i3];

    ierr = intp_tval(ctx, xp,yp,x1,x2,x3,y1,y2,y3,z1,z2,z3,zx1,zx2,zx3,zy1,zy2,zy3,
            pz,&dum,&dum,0);
    if (ierr != 0)
        return(-2);               
    return(0);
 
L3: /* P IS OUTSIDE OF THE MESH BOUNDARY.  EXTRAPOLATE TO P BY
       PASSING A LINEAR FUNCTION OF ONE VARIABLE THROUGH THE
       VALUE AND DIRECTIONAL DERIVATIVE (IN THE DIRECTION
       P-Q) OF THE INTERPOLATORY SURFACE (TVAL) AT Q WHERE
       Q IS THE CLOSEST BOUNDARY POINT TO P. */
 
    /* DETERMINE Q BY TRAVERSING THE BOUNDARY STARTING FROM
       THE RIGHTMOST VISIBLE NODE I1. */
 
    n2 = i1;
 
    /* SET N1 TO THE LAST NONZERO NEIGHBOR OF N2 AND COMPUTE DP */
 
L4: indx = iend[n2] - 1;
    n1 = iadj[indx];
    x1 = x[n1];
    y1 = y[n1];
    x2 = x[n2];
    y2 = y[n2];
    dp = (x1 - x2) * (xp - x2) + (y1 - y2) * (yp - y2);
    if (dp <= 0.0)
        goto L5;  

    if ((xp - x1) * (x2 - x1) + (yp - y1) * (y2 - y1) > 0.0)
        goto L8;

    n2 = n1;
    goto L4;
 
    /* N2 IS THE CLOSEST BOUNDARY POINT TO P.  COMPUTE PARTIAL
       DERIVATIVES AT N2. */ 
L5:
    intp_gradl(ctx, nn,n2,x,y,z,iadj,iend,&zx2,&zy2);

    /* COMPUTE EXTRAPOLATED VALUE AT P */
 
    *pz = z[n2] + zx2 * (xp - x2) + zy2 * (yp - y2);
    return(1);
 
    /* THE CLOSEST BOUNDARY POINT Q LIES ON N2-N1.  COMPUTE
       PARTIALS AT N1 AND N2. */
 
L8: intp_gradl(ctx, nn,n1,x,y,z,iadj,iend,&zx1,&zy1);
    intp_gradl(ctx, nn,n2,x,y,z,iadj,iend,&zx2,&zy2);

    /* COMPUTE Q, ITS BARYCENTRIC COORDINATES, AND THE CARDINAL
       FUNCTIONS FOR EXTRAPOLATION */
 
/* L10 */
    u = x2 - x1;
    v = y2 - y1; 
    r1 = dp / (u * u + v * v);
    r2 = 1.0 - r1;
    xq = r1 * x1 + r2 * x2;
    yq = r1 * y1 + r2 * y2;
    f1 = r1 * r1 * r2;
    f2 = r1 * r2 * r2;
    a1 = r1 + (f1 - f2);
    a2 = r2 - (f1 - f2);
    b1 = u * f1;
    b2 = -u * f2;
    c1 = v * f1;
    c2 = -v * f2;
 
    /* COMPUTE THE VALUE OF THE INTERPOLATORY SURFACE (TVAL) AT Q */
 
    *pz = a1 * z[n1] + a2 * z[n2] + b1 * zx1 + b2 * zx2 + c1 * zy1 + c2 * zy2;
 
    /* COMPUTE THE EXTRAPOLATED VALUE AT P */
 
    *pz += (r1 * zx1 + r2 * zx2) * (xp - xq) + (r1 * zy1 + r2 * zy2)*(yp - yq);
    return(1);
}

/* -##--------------------------------------------------------------------- */
/*  intp_rotate(n,c,s,x,y)                                                  */
/*                                           ( C  S)                        */
/*  THIS ROUTINE APPLIES THE GIVENS ROTATION (     ) TO THE                 */
/*                                           (-S  C)                        */
/*                (X(1) ... X(N))                                           */
/*  2 BY N MATRIX (             ).  THIS ROUTINE WAS TAKEN FROM             */
/*                (Y(1) ... Y(N))                                           */
/*  LINPACK.                                                                */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      N - NUMBER OF COLUMNS TO BE ROTATED.                                */
/*      C,S - ELEMENTS OF THE GIVENS ROTATION.                              */
/*            THESE MAY BE DETERMINED BY SUBROUTINE GIVENS.                 */
/*      X,Y - VECTORS OF LENGTH .GE. N CONTAINING THE 2-VECTORS TO BE       */
/*            ROTATED.                                                      */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      X,Y - ROTATED VECTORS                                               */
/*                                                                          */
/*  Functions called: none                                                  */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_rotate(TDAContext *ctx, int n,double c,double s,double *x,double *y) 
{
    (void)ctx;        /* unused: the signature is shared */
    register int i; 
    double xi,yi;

    if (n <=  0 || (c == 1.0 && s == 0.0))
        return;                                           

    for (i = 1; i <= n; ++i) {
        xi = x[i];
        yi = y[i];
        x[i] = c * xi + s * yi;
        y[i] = -s * xi + c * yi;
    }
}

/* -##--------------------------------------------------------------------- */
/*  intp_givens(a,b,c,s)                                                    */
/*                                                                          */
/*  THIS ROUTINE CONSTRUCTS THE GIVENS PLANE ROTATION --                    */
/*      ( C  S)                                                             */
/*  G = (     ) WHERE C*C + S*S = 1 -- WHICH ZEROS THE SECOND               */
/*      (-S  C)                                                             */
/*  ENTRY OF THE 2-VECTOR (A B)-TRANSPOSE.  A CALL TO GIVENS                */
/*  IS NORMALLY FOLLOWED BY A CALL TO ROTATE WHICH APPLIES                  */
/*  THE TRANSFORMATION TO A 2 BY N MATRIX.  THIS ROUTINE WAS                */
/*  TAKEN FROM LINPACK.                                                     */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      A,B - COMPONENTS OF THE 2-VECTOR TO BE ROTATED.                     */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      A - OVERWRITTEN BY R = +/-SQRT(A*A + B*B)                           */
/*      B - OVERWRITTEN BY A VALUE Z WHICH ALLOWS C AND S TO BE RECOVERED   */
/*          AS FOLLOWS -                                                    */
/*                         C = SQRT(1-Z*Z), S=Z IF ABS(Z)                   */
/*                             .LE. 1.                                      */
/*                         C = 1/Z, S = SQRT(1-C*C) IF                      */
/*                             ABS(Z) .GT. 1.                               */
/*                                                                          */
/*                         C - +/-(A/R)                                     */
/*                                                                          */
/*                         S - +/-(B/R)                                     */
/*                                                                          */
/*  Functions called: none                                                  */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_givens(TDAContext *ctx, double *a,double *b,double *c,double *s) 
{
    (void)ctx;        /* unused: the signature is shared */
    double aa,bb,r,u,v;

    aa = *a;
    bb = *b;
    if (fabs(aa) <= fabs(bb))
        goto L1;             
 
    u = aa + aa;
    v = bb / u;
    r = sqrt(0.25 + v * v) * u;
    *c = aa / r;
    *s = v * (*c + *c);
 
    *b = *s;
    *a = r;
    return;  
 
L1:
    if (bb == 0.0)
        goto L2;
    u = bb + bb;
    v = aa / u;
 
    *a = sqrt(0.25 + v * v) * u;
    *s = bb / *a;
    *c = v * (*s + *s);
 
    *b = 1.0;
    if (*c != 0.0)
        *b = 1.0 / *c;
    return;
 
L2:
    *c = 1.0;
    *s = 0.0;
    return;  
}

/* -##--------------------------------------------------------------------- */
/*  intp_setup(xk,yk,zk,xi,yi,zi,s1,s2,r,row)                               */
/*                                                                          */
/*  THIS ROUTINE SETS UP THE I-TH ROW OF AN AUGMENTED RE-                   */
/*  GRESSION MATRIX FOR A WEIGHTED LEAST-SQUARES FIT OF A                   */
/*  QUADRATIC FUNCTION Q(X,Y) TO A SET OF DATA VALUES Z WHERE               */
/*  Q(XK,YK) = ZK.  THE FIRST 3 COLUMNS (QUADRATIC TERMS) ARE               */
/*  SCALED BY 1/S2 AND THE FOURTH AND FIFTH COLUMNS (LINEAR                 */
/*  TERMS) ARE SCALED BY 1/S1.  THE WEIGHT IS (R-D)/(R*D) IF                */
/*  R .GT. D AND 0 IF R .LE. D, WHERE D IS THE DISTANCE                     */
/*  BETWEEN NODES I AND K.                                                  */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      XK,YK,ZK - COORDINATES AND DATA VALUE                               */
/*                 AT NODE K -- INTERPOLATED BY Q.                          */
/*      XI,YI,ZI - COORDINATES AND DATA VALUE AT NODE I.                    */
/*      S1,S2 - INVERSE SCALE FACTORS.                                      */
/*      R - RADIUS OF INFLUENCE ABOUT NODE K DEFINING THE WEIGHT.           */
/*      ROW - VECTOR OF LENGTH 6.                                           */
/*                                                                          */
/*  OUTPUT PARAMETER:                                                       */
/*      ROW - VECTOR CONTAINING A ROW OF THE                                */
/*            AUGMENTED REGRESSION MATRIX.                                  */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

void intp_setup(TDAContext *ctx, double xk,double yk,double zk,double xi,double yi,double zi, double s1,double s2,double r,double *row)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double dx,dy,dxsq,dysq,d,w,w1,w2;

    dx = xi - xk;
    dy = yi - yk;
    dxsq = dx * dx;
    dysq = dy * dy;
    d = sqrt(dxsq + dysq);

    if (d <= 0.0 || d >= r) {
        for (i = 1; i <= 6; ++i)  
            row[i] = 0.0;
        return;
    }                        
    w = (r - d) / r / d;
    w1 = w / s1;
    w2 = w / s2;
    row[1] = dxsq * w2;
    row[2] = dx * dy * w2;
    row[3] = dysq * w2;
    row[4] = dx * w1;
    row[5] = dy * w1;
    row[6] = (zi - zk) * w;
}

/* -##--------------------------------------------------------------------- */
/* intp_tval(x,y,x1,x2,x3,y1,y2,y3,z1,z2,z3,zx1,zx2,zx3,zy1,zy2,zy3,w,wx,wy)*/
/*                                                                          */
/*  GIVEN FUNCTION VALUES AND FIRST PARTIAL DERIVATIVES AT                  */
/*  THE THREE VERTICES OF A TRIANGLE, THIS ROUTINE DETERMINES               */
/*  A FUNCTION W WHICH AGREES WITH THE GIVEN DATA, RETURNING                */
/*  THE VALUE AND (OPTIONALLY) FIRST PARTIAL DERIVATIVES OF W               */
/*  AT A POINT (X,Y) IN THE TRIANGLE.  THE INTERPOLATION                    */
/*  METHOD IS EXACT FOR QUADRATIC POLYNOMIAL DATA.  THE                     */
/*  TRIANGLE IS PARTITIONED INTO THREE SUBTRIANGLES WITH                    */
/*  EQUAL AREAS.  W IS CUBIC IN EACH SUBTRIANGLE AND ALONG                  */
/*  THE EDGES, BUT HAS ONLY ONE CONTINUOUS DERIVATIVE ACROSS                */
/*  EDGES.  THE NORMAL DERIVATIVE OF W VARIES LINEARLY ALONG                */
/*  EACH OUTER EDGE.  THE VALUES AND PARTIAL DERIVATIVES OF W               */
/*  ALONG A TRIANGLE EDGE DEPEND ONLY ON THE DATA VALUES AT                 */
/*  THE ENDPOINTS OF THE EDGE.  THUS THE METHOD YIELDS C-1                  */
/*  CONTINUITY WHEN USED TO INTERPOLATE OVER A TRIANGULAR                   */
/*  GRID.  THIS ALGORITHM IS DUE TO C. L. LAWSON.                           */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      X,Y - COORDINATES OF A POINT AT WHICH W IS TO BE EVALUATED.         */
/*      X1,X2,X3,Y1,Y2,Y3 - COORDINATES OF THE VERTICES OF                  */
/*                          A TRIANGLE CONTAINING (X,Y).                    */
/*      Z1,Z2,Z3 - FUNCTION VALUES AT THE VERTICES TO BE INTERPOLATED.      */
/*      ZX1,ZX2,ZX3 - X-DERIVATIVE VALUES AT THE VERTICES.                  */
/*      ZY1,ZY2,ZY3 - Y-DERIVATIVE VALUES AT THE VERTICES.                  */
/*      IFLAG - OPTION INDICATOR                                            */
/*              IFLAG = 0 IF ONLY W IS TO BE COMPUTED.                      */
/*              IFLAG = 1 IF W, WX, AND WY ARE TO BE RETURNED.              */
/*                                                                          */
/*  INPUT PARAMETERS ARE NOT ALTERED BY THIS ROUTINE.                       */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      W - ESTIMATED VALUE OF THE INTERPOLATORY FUNCTION AT (X,Y) IF       */
/*          IER = 0.  OTHERWISE W = 0.                                      */
/*      WX,WY - PARTIAL DERIVATIVES OF W AT (X,Y) IF IER = 0 AND            */
/*              IFLAG = 1, UNCHANGED IF IFLAG .NE. 1, ZERO                  */
/*                         IF IER .NE. 0 AND IFLAG = 1.                     */
/*                                                                          */
/*  Return:  0 IF NO ERRORS WERE ENCOUNTERED.                               */
/*           1 IF THE VERTICES OF THE TRIANGLE ARE COLLINEAR.               */
/*                                                                          */
/*  Functions called: none.                                                 */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_tval(TDAContext *ctx, double x,double y,double x1,double x2,double x3, double y1,double y2,double y3,double z1,double z2,double z3, double zx1,double zx2,double zx3,double zy1,double zy2,double zy3, double *w,double *wx,double *wy,int iflag)
{
    int i,ip1,ip2,ip3;
    double u[4],v[4],sl[4],area,xp,yp,r[4],rx[4],ry[4],phi[4],phix[4];
    double phiy[4],rmin,c1,c2,ro[4],rox[4],roy[4],f[4],g[4],gx[4],gy[4];
    double p[4],px[4],py[4],q[4],qx[4],qy[4],a[4],ax[4],ay[4],b[4],bx[4];
    double by[4],c[4],cx[4],cy[4];
 
    u[1] = x3 - x2;
    u[2] = x1 - x3;
    u[3] = x2 - x1;
 
    v[1] = y3 - y2;
    v[2] = y1 - y3;
    v[3] = y2 - y1;
 
    for (i = 1; i <= 3; ++i) {
        sl[i] = u[i] * u[i] + v[i] * v[i];
    }                
    area = u[1] * v[2] - u[2] * v[1];
    if (area == 0.0) {            
        *w = 0.0;
        if (iflag != 1)
            return(1);
        *wx = *wy = 0.0;
        return(1);
    }
    r[1] = (u[1] * (y - y2) - v[1] * (x - x2)) / area;
    xp = x - x1;
    yp = y - y1;
    r[2] = (u[2] * yp - v[2] * xp) / area;
    r[3] = (u[3] * yp - v[3] * xp) / area;
 
    phi[1] = r[2] * r[3];
    phi[2] = r[3] * r[1];
    phi[3] = r[1] * r[2];
 
    rmin = dmin(ctx, r[1],dmin(ctx, r[2],r[3]));
    if (rmin != r[1])
        goto L3;
    ip1 = 1;
    ip2 = 2;
    ip3 = 3;
    goto L5;

L3: if (rmin != r[2])
        goto L4;
    ip1 = 2;
    ip2 = 3;
    ip3 = 1;
    goto L5;

L4: ip1 = 3;
    ip2 = 1;
    ip3 = 2;
 
L5: c1 = rmin * rmin / 2.0;
    c2 = rmin / 3.0;
    ro[ip1] = (phi[ip1] + 5.0 * c1 / 3.0) * r[ip1] - c1;
    ro[ip2] = c1 * (r[ip3] - c2);
    ro[ip3] = c1 * (r[ip2] - c2);
 
    f[1] = 3.0 * (sl[2] - sl[3]) / sl[1];
    f[2] = 3.0 * (sl[3] - sl[1]) / sl[2];
    f[3] = 3.0 * (sl[1] - sl[2]) / sl[3];
 
    g[1] = (r[2] - r[3]) * phi[1] + f[1] * ro[1] - ro[2] + ro[3];
    g[2] = (r[3] - r[1]) * phi[2] + f[2] * ro[2] - ro[3] + ro[1];
    g[3] = (r[1] - r[2]) * phi[3] + f[3] * ro[3] - ro[1] + ro[2];
 
    for (i = 1; i <= 3; ++i) {
        p[i] = g[i] + phi[i];
        q[i] = g[i] - phi[i];
    }              
    a[1] = r[1] + g[3] - g[2];
    a[2] = r[2] + g[1] - g[3];
    a[3] = r[3] + g[2] - g[1];
 
    b[1] = u[3] * p[3] + u[2] * q[2];
    b[2] = u[1] * p[1] + u[3] * q[3];
    b[3] = u[2] * p[2] + u[1] * q[1];
 
    c[1] = v[3] * p[3] + v[2] * q[2];
    c[2] = v[1] * p[1] + v[3] * q[3];
    c[3] = v[2] * p[2] + v[1] * q[1];
 
    *w = a[1] * z1 + a[2] * z2 + a[3] * z3 + (b[1] * zx1 + b[2] * zx2 
        + b[3] * zx3 + c[1] * zy1 + c[2] * zy2 + c[3] * zy3) / 2.0;
    if (iflag != 1)
        return(0);               
 
    for (i = 1; i <= 3; ++i) {
        rx[i] = -v[i] / area;
        ry[i] = u[i] / area;
    }               
    phix[1] = r[2] * rx[3] + rx[2] * r[3];
    phiy[1] = r[2] * ry[3] + ry[2] * r[3];
    phix[2] = r[3] * rx[1] + rx[3] * r[1];
    phiy[2] = r[3] * ry[1] + ry[3] * r[1];
    phix[3] = r[1] * rx[2] + rx[1] * r[2];
    phiy[3] = r[1] * ry[2] + ry[1] * r[2];
 
    rox[ip1] = rx[ip1] * (phi[ip1] + 5.0 *c1) + r[ip1] * (phix[ip1] - rx[ip1]);
    roy[ip1] = ry[ip1] * (phi[ip1] + 5.0 *c1) + r[ip1] * (phiy[ip1] - ry[ip1]);
    rox[ip2] = rx[ip1] * (phi[ip2] - c1) + c1 * rx[ip3];
    roy[ip2] = ry[ip1] * (phi[ip2] - c1) + c1 * ry[ip3];
    rox[ip3] = rx[ip1] * (phi[ip3] - c1) + c1 * rx[ip2];
    roy[ip3] = ry[ip1] * (phi[ip3] - c1) + c1 * ry[ip2];
 
    gx[1] = (rx[2] - rx[3]) * phi[1] + (r[2] - r[3]) * phix[1] + f[1] * rox[1] - rox[2] + rox[3];
    gy[1] = (ry[2] - ry[3]) * phi[1] + (r[2] - r[3]) * phiy[1] + f[1] * roy[1] - roy[2] + roy[3];
    gx[2] = (rx[3] - rx[1]) * phi[2] + (r[3] - r[1]) * phix[2] + f[2] * rox[2] - rox[3] + rox[1];
    gy[2] = (ry[3] - ry[1]) * phi[2] + (r[3] - r[1]) * phiy[2] + f[2] * roy[2] - roy[3] + roy[1];
    gx[3] = (rx[1] - rx[2]) * phi[3] + (r[1] - r[2]) * phix[3] + f[3] * rox[3] - rox[1] + rox[2];
    gy[3] = (ry[1] - ry[2]) * phi[3] + (r[1] - r[2]) * phiy[3] + f[3] * roy[3] - roy[1] + roy[2];
 
    for (i = 1; i <= 3; ++i) {
        px[i] = gx[i] + phix[i];
        py[i] = gy[i] + phiy[i];
        qx[i] = gx[i] - phix[i];
        qy[i] = gy[i] - phiy[i];
    }               
    ax[1] = rx[1] + gx[3] - gx[2];
    ay[1] = ry[1] + gy[3] - gy[2];
    ax[2] = rx[2] + gx[1] - gx[3];
    ay[2] = ry[2] + gy[1] - gy[3];
    ax[3] = rx[3] + gx[2] - gx[1];
    ay[3] = ry[3] + gy[2] - gy[1];
 
    bx[1] = u[3] * px[3] + u[2] * qx[2];
    by[1] = u[3] * py[3] + u[2] * qy[2];
    bx[2] = u[1] * px[1] + u[3] * qx[3];
    by[2] = u[1] * py[1] + u[3] * qy[3];
    bx[3] = u[2] * px[2] + u[1] * qx[1];
    by[3] = u[2] * py[2] + u[1] * qy[1];
 
    cx[1] = v[3] * px[3] + v[2] * qx[2];
    cy[1] = v[3] * py[3] + v[2] * qy[2];
    cx[2] = v[1] * px[1] + v[3] * qx[3];
    cy[2] = v[1] * py[1] + v[3] * qy[3];
    cx[3] = v[2] * px[2] + v[1] * qx[1];
    cy[3] = v[2] * py[2] + v[1] * qy[1];
 
    *wx = ax[1] * z1 + ax[2] * z2 + ax[3] * z3 + (bx[1] * zx1 +
          bx[2] * zx2 + bx[3] * zx3 + cx[1] * zy1 + cx[2] * zy2 + cx[3] * zy3) / 2.0;
    *wy = ay[1] * z1 + ay[2] * z2 + ay[3] * z3 + (by[1] * zx1 +
          by[2] * zx2 + by[3] * zx3 + cy[1] * zy1 + cy[2] * zy2 + cy[3] * zy3) / 2.0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  intp_gradl(n,k,x,y,z,iadj,iend,dx,dy)                                   */
/*                                                                          */
/*  GIVEN A THIESSEN TRIANGULATION OF N POINTS IN THE PLANE                 */
/*  WITH ASSOCIATED DATA VALUES Z, THIS SUBROUTINE ESTIMATES                */
/*  X AND Y PARTIAL DERIVATIVES AT NODE K.  THE DERIVATIVES                 */
/*  ARE TAKEN TO BE THE PARTIALS AT K OF A QUADRATIC FUNCTION               */
/*  WHICH INTERPOLATES Z(K) AND FITS THE DATA VALUES AT A SET               */
/*  OF NEARBY NODES IN A WEIGHTED LEAST SQUARES SENSE. A MAR-               */
/*  QUARDT STABILIZATION FACTOR IS USED IF NECESSARY TO ENSURE              */
/*  A WELL-CONDITIONED SYSTEM AND A LINEAR FITTING FUNCTION IS              */
/*  USED IF N .LT. 6.  THUS, A UNIQUE SOLUTION EXISTS UNLESS                */
/*  THE NODES ARE COLLINEAR.                                                */
/*  AN ALTERNATIVE ROUTINE, GRADG, EMPLOYS A GLOBAL METHOD                  */
/*  TO COMPUTE THE PARTIAL DERIVATIVES AT ALL OF THE NODES AT               */
/*  ONCE.  THAT METHOD IS MORE EFFICIENT (WHEN ALL PARTIALS                 */
/*  ARE NEEDED) AND MAY BE MORE ACCURATE, DEPENDING ON THE DATA.            */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      N - NUMBER OF NODES IN THE TRIANGULATION.  N .GE. 3.                */
/*      K - NODE AT WHICH DERIVATIVES ARE SOUGHT.  1 .LE. K .LE. N.         */
/*      X,Y - N-VECTORS CONTAINING THE CARTESIAN COORDINATES OF THE NODES.  */
/*      Z - N-VECTOR CONTAINING THE DATA VALUES ASSOCIATED WITH THE NODES.  */
/*      IADJ - SET OF ADJACENCY LISTS.                                      */
/*      IEND - POINTERS TO THE ENDS OF ADJACENCY LISTS FOR EACH NODE.       */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      DX,DY - ESTIMATED PARTIAL DERIVATIVES AT NODE K UNLESS IER .LT. 0.  */
/*                                                                          */
/*  Return: r > 0 IF NO ERRORS WERE ENCOUNTERED. r contains the number      */
/*                of nodes (including k) used in the fit. r = 3,4 or 5      */
/*                implies a linear fit.                                     */
/*          -1  IF N OR K IS OUT OF RANGE.                                  */
/*          -2  IF ALL NODES ARE COLLINEAR.                                 */
/*                                                                          */
/*  Functions called: intp_getnp(), intp_setup(), intp_givens(),            */
/*  intp_rotate().                                                          */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */
/*  Note: we change the order of the indices in the doubly-indexed array    */
/*  a[][].                                                                  */

int intp_gradl(TDAContext *ctx, int n,int k,double *x,double *y,double *z,int *iadj,int *iend, double *dx,double *dy)
{
    int nn = 0,kk = 0,lmn = 0,lmx = 0,lmin = 0,lmax = 0,lm1 = 0,lnp = 0,npts[32],np = 0,i = 0,j = 0,im1 = 0,jp1 = 0,ip1 = 0,l = 0;
    double sum = 0.0,ds = 0.0,r = 0.0,rs = 0.0,rtol = 0.0,avsq = 0.0,av = 0.0,xk = 0.0,yk = 0.0,zk = 0.0,c = 0.0,s = 0.0,ddmin = 0.0,dtol = 0.0,sf = 0.0;
    double a[50],*apt = NULL,*bpt = NULL;

    dtol = 0.01;
    rtol = 1.e-5;
    sf = 1.0;
    lmn = 10;
    lmx = 30;

    *dx = *dy = 0.0;
    nn = n;
    kk = k;
 
    if (nn < 3 || kk < 1 || kk > nn)  
        return(-1);
    lmin = imin(ctx, lmn,nn);
    lmax = imin(ctx, lmx,nn);
 
    sum = 0.0;
    npts[1] = kk;
    lm1 = lmin - 1;

    for (lnp = 2; lnp <= lm1; ++lnp) {
        intp_getnp(ctx, x,y,iadj,iend,lnp,npts,&ds);   
        sum += ds;             
    }
    for (lnp = lmin; lnp <= lmax; ++lnp) {
        intp_getnp(ctx, x,y,iadj,iend,lnp,npts,&rs);   

        if ((rs - ds) / ds <= rtol)
            goto L2;
        if (lnp > 6)
            goto L3;
L2:     sum += rs;          
    }
    rs = 1.1 * rs;
    lnp = lmax + 1;
 
L3:
    avsq = sum / (double)(lnp - 2);
    av = sqrt(avsq);
    r = sqrt(rs);
    xk = x[kk];
    yk = y[kk];
    zk = z[kk];
    if (lnp < 7)
        goto L12;
 
    for (i = 1; i <= 5; ++i) {
        np = npts[i + 1];
        /** CALL SETUP (XK,YK,ZK,X(NP),Y(NP),Z(NP),AV,AVSQ,R, A(1,I)) **/

        apt = a + i * 7;
        intp_setup(ctx, xk,yk,zk,x[np],y[np],z[np],av,avsq,r,apt);

        if (i == 1)
            continue;

        im1 = i - 1;
        for (j = 1; j <= im1; ++j) {
            jp1 = j + 1;
            l = 6 - j;
            /** CALL GIVENS (A(J,J),A(J,I),C,S) **/
            intp_givens(ctx, a + j * 7 + j,a + i * 7 + j,&c,&s);

            /** CALL ROTATE(L,C,S,A(JP1,J),A(JP1,I)) **/
            apt = a + j * 7 + jp1 - 1;
            bpt = a + i * 7 + jp1 - 1;
            intp_rotate(ctx, l,c,s,apt,bpt);
        }
    }
    i = 7;
L6:
    if (i == lnp)
        goto L8;
    np = npts[i];
    /** CALL SETUP (XK,YK,ZK,X(NP),Y(NP),Z(NP),AV,AVSQ,R, A(1,6)) **/

    apt = a + 6 * 7;

    intp_setup(ctx, xk,yk,zk,x[np],y[np],z[np],av,avsq,r,apt);


    for (j = 1; j <= 5; ++j) {
        jp1 = j + 1;
        l = 6 - j;
        /** CALL GIVENS (A(J,J),A(J,6),C,S) **/
        intp_givens(ctx, a + j * 7 + j,a + 6 * 7 + j,&c,&s);

        /** CALL ROTATE(L,C,S,A(JP1,J),A(JP1,6)) **/

        apt = a + j * 7 + jp1 - 1;
        bpt = a + 6 * 7 + jp1 - 1;
        intp_rotate(ctx, l,c,s,apt,bpt);
    }
    i++;          
    goto L6;
 
L8: ddmin = dmin(ctx, fabs(a[7 + 1]),fabs(a[2 * 7 + 2]));
    ddmin = dmin(ctx, ddmin,fabs(a[3 * 7 + 3]));
    ddmin = dmin(ctx, ddmin,fabs(a[4 * 7 + 4]));
    ddmin = dmin(ctx, ddmin,fabs(a[5 * 7 + 5]));

    if (ddmin >= dtol)
        goto L15;
    if (lnp > lmax)
        goto L9;
 
    lnp++;          
    if (lnp <= lmax) {
        intp_getnp(ctx, x,y,iadj,iend,lnp,npts,&rs);   
    }
    r = sqrt(1.1 * rs);
    goto L6;
 
L9: for (i = 1; i <= 3; ++i) {
        a[6 * 7 + i] = sf;
        ip1 = i + 1;
        for (j = ip1; j <= 6; ++j) {
            a[6 * 7 + j] = 0.0;
        }
        for (j = i; j <= 5; ++j) {
            jp1 = j + 1;
            l = 6 - j;
            /** CALL GIVENS (A(J,J),A(J,6),C,S) **/
            intp_givens(ctx, a + j * 7 + j,a + 6 * 7 + j,&c,&s);

            /** CALL ROTATE(L,C,S,A(JP1,J),A(JP1,6)) **/
            apt = a + j * 7 + jp1 - 1;
            bpt = a + 6 * 7 + jp1 - 1;
            intp_rotate(ctx, l,c,s,apt,bpt);
        }
    }
    goto L14;
 
L12:
    np = npts[2];

    /** CALL SETUP (XK,YK,ZK,X(NP),Y(NP),Z(NP),AV,AVSQ,R, A(1,4)) **/

    apt = a + 4 * 7;
    intp_setup(ctx, xk,yk,zk,x[np],y[np],z[np],av,avsq,r,apt);

    np = npts[3];

    /** CALL SETUP (XK,YK,ZK,X(NP),Y(NP),Z(NP),AV,AVSQ,R,A(1,5)) **/
    apt = a + 5 * 7;
    intp_setup(ctx, xk,yk,zk,x[np],y[np],z[np],av,avsq,r,apt);


    /** CALL GIVENS (A(4,4),A(4,5),C,S) **/
    intp_givens(ctx, a + 4 * 7 + 4,a + 5 * 7 + 4,&c,&s);

    /** CALL ROTATE(2,C,S,A(5,4),A(5,5)) **/
    apt = a + 4 * 7 + 4;
    bpt = a + 5 * 7 + 4;
    intp_rotate(ctx, 2,c,s,apt,bpt);
    if (lnp == 4)
        goto L14;
 
    lm1 = lnp - 1;
    for (i = 4; i <= lm1; ++i) {
        np = npts[i];
        /** CALL SETUP (XK,YK,ZK,X(NP),Y(NP),Z(NP),AV,AVSQ,R, A(1,6)) **/

        apt = a + 6 * 7;
        intp_setup(ctx, xk,yk,zk,x[np],y[np],z[np],av,avsq,r,apt);

        /** CALL GIVENS (A(4,4),A(4,6),C,S) **/
        intp_givens(ctx, a + 4 * 7 + 4,a + 6 * 7 + 4,&c,&s);

        /** CALL ROTATE (2,C,S,A(5,4),A(5,6)) **/
        apt = a + 4 * 7 + 4;
        bpt = a + 6 * 7 + 4;
        intp_rotate(ctx, 2,c,s,apt,bpt);

        /** CALL GIVENS (A(5,5),A(5,6),C,S) **/
        intp_givens(ctx, a + 5 * 7 + 5,a + 6 * 7 + 5,&c,&s);

        /** CALL ROTATE (1,C,S,A(6,5),A(6,6)) **/
        apt = a + 5 * 7 + 5;
        bpt = a + 6 * 7 + 5;
        intp_rotate(ctx, 1,c,s,apt,bpt);
    }
L14:
    ddmin = dmin(ctx, fabs(a[4 * 7 + 4]),fabs(a[5 * 7 + 5]));   
    if (ddmin < dtol)    
        return(-2);
 
L15:
    *dy = a[5 * 7 + 6] / a[5 * 7 + 5];
    *dx = (a[4 * 7 + 6] - a[4 * 7 + 5] * *dy) / a[4 * 7 + 4] / av;
    *dy = *dy / av;
    return(lnp - 1);
}

/* -##--------------------------------------------------------------------- */
/*  intp_getnp(x,y,iadj,iend,l,npts,ds)                                     */
/*                                                                          */
/*  GIVEN A THIESSEN TRIANGULATION OF N NODES AND AN ARRAY                  */
/*  NPTS CONTAINING THE INDICES OF L-1 NODES ORDERED BY                     */
/*  EUCLIDEAN DISTANCE FROM NPTS(1), THIS SUBROUTINE SETS                   */
/*  NPTS(L) TO THE INDEX OF THE NEXT NODE IN THE SEQUENCE --                */
/*  THE NODE, OTHER THAN NPTS(1),...,NPTS(L-1), WHICH IS                    */
/*  CLOSEST TO NPTS(1).  THUS, THE ORDERED SEQUENCE OF K                    */
/*  CLOSEST NODES TO N1 (INCLUDING N1) MAY BE DETERMINED BY                 */
/*  K-1 CALLS TO GETNP WITH NPTS(1) = N1 AND L = 2,3,...,K                  */
/*  FOR K .GE. 2.                                                           */
/*  THE ALGORITHM USES THE FACT THAT, IN A THIESSEN TRIAN-                  */
/*  GULATION, THE K-TH CLOSEST NODE TO A GIVEN NODE N1 IS A                 */
/*  NEIGHBOR OF ONE OF THE K-1 CLOSEST NODES TO N1.                         */
/*                                                                          */
/*  INPUT PARAMETERS:                                                       */
/*      X,Y - VECTORS OF LENGTH N CONTAINING                                */
/*            THE CARTESIAN COORDINATES OF THE NODES.                       */
/*      IADJ - SET OF ADJACENCY LISTS OF NODES IN THE TRIANGULATION.        */
/*      IEND - POINTERS TO THE ENDS OF ADJACENCY LISTS FOR EACH NODE        */
/*      L - NUMBER OF NODES IN THE SEQUENCE ON OUTPUT. 2 .LE. L .LE. N.     */
/*      NPTS - ARRAY OF LENGTH .GE. L CONTAINING THE INDICES OF THE L-1     */
/*             CLOSEST NODES TO NPTS(1) IN THE FIRST L-1 LOCATIONS.         */
/*                                                                          */
/*  IADJ AND IEND MAY BE CREATED BY SUBROUTINE TRMESH.                      */
/*                                                                          */
/*  INPUT PARAMETERS OTHER THAN NPTS ARE NOT ALTERED BY THIS                */
/*  ROUTINE.                                                                */
/*                                                                          */
/*  OUTPUT PARAMETERS:                                                      */
/*      NPTS - UPDATED WITH THE INDEX OF THE L-TH CLOSEST NODE TO           */
/*             NPTS(1) IN POSITION L UNLESS IER = 1.                        */
/*      DS - SQUARED EUCLIDEAN DISTANCE BETWEEN NPTS(1) AND NPTS(L)         */
/*                                                                          */
/*  Return:  0 IF NO ERRORS WERE ENCOUNTERED.                               */
/*           1 IF L IS OUT OF RANGE.                                        */
/*                                                                          */
/*  Functions called: intp_getnp(), intp_setup(), intp_givens(),            */
/*  intp_rotate().                                                          */
/*                                                                          */
/*  The original FORTRAN code is CALGO 624 (Robert Renka).                  */
/*                                                                          */

int intp_getnp(TDAContext *ctx, double *x,double *y,int *iadj,int *iend,int l,int *npts, double *ds)
{
    int lm1,n1,i,ni,np,indf,indl,indx,nb;
    double x1,y1,dnp,dnb,tmp1,tmp2;
 
    lm1 = l - 1;
    if (lm1 < 1)
        return(1);
                
    n1 = npts[1];
    x1 = x[n1];
    y1 = y[n1];
 
    for (i = 1; i <= lm1; ++i) {
        ni = npts[i];
        iend[ni] = -iend[ni];
    }
    np = 0;
    dnp = 0.0;
 
    for (i = 1; i <= lm1; ++i) {
        ni = npts[i];
        indf = 1;
        if (ni > 1)
            indf = iabs(ctx, iend[ni - 1]) + 1;
        indl = -iend[ni];
 
        for (indx = indf; indx <= indl; ++indx) {
            nb = iadj[indx];
            if (nb == 0 || iend[nb] < 0)
                continue;        
 
            tmp1 = x[nb] - x1;                      
            tmp2 = y[nb] - y1;        

            dnb = tmp1 * tmp1 + tmp2 * tmp2;        

            if (np != 0 && dnb >= dnp) 
                continue;

            np = nb;
            dnp = dnb;
        }
    }
    npts[l] = np;
    *ds = dnp;
 
    for (i = 1; i <= lm1; ++i) {
        ni = npts[i];
        iend[ni] = -iend[ni];
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  ak_idsfft(ncp,ndp,xd,yd,zd,nxi,nyi,xi,yi,zi,iwk,wk)                     */
/*                                                                          */
/*  SUBROUTINE  IDSFFT(MD,MCP,NDP,XD,YD,ZD,NXI,NYI,XI,YI,ZI,IWK,WK)         */
/*                                                                          */
/*  THIS SUBROUTINE PERFORMS SMOOTH SURFACE FITTING WHEN THE PRO-           */
/*  JECTIONS OF THE DATA POINTS IN THE X-Y PLANE ARE IRREGULARLY            */
/*  DISTRIBUTED IN THE PLANE.                                               */
/*                                                                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    ncp = NUMBER OF ADDITIONAL DATA POINTS USED FOR ESTI-                 */
/*          MATING PARTIAL DERIVATIVES AT EACH DATA POINT                   */
/*          (MUST BE 2 OR GREATER, BUT SMALLER THAN ndp),                   */
/*    ndp = NUMBER OF DATA POINTS (MUST BE 4 OR GREATER),                   */
/*    xd  = ARRAY OF DIMENSION ndp CONTAINING THE X                         */
/*          COORDINATES OF THE DATA POINTS,                                 */
/*    yd  = ARRAY OF DIMENSION ndp CONTAINING THE Y                         */
/*          COORDINATES OF THE DATA POINTS,                                 */
/*    zd  = ARRAY OF DIMENSION ndp CONTAINING THE Z                         */
/*          COORDINATES OF THE DATA POINTS,                                 */
/*    nxi = NUMBER OF OUTPUT GRID POINTS IN THE X COORDINATE                */
/*          (MUST BE 1 OR GREATER),                                         */
/*    nyi = NUMBER OF OUTPUT GRID POINTS IN THE Y COORDINATE                */
/*          (MUST BE 1 OR GREATER),                                         */
/*    xi  = ARRAY OF DIMENSION NXI CONTAINING THE X                         */
/*          COORDINATES OF THE OUTPUT GRID POINTS,                          */
/*    yi  = ARRAY OF DIMENSION NYI CONTAINING THE Y                         */
/*          COORDINATES OF THE OUTPUT GRID POINTS.                          */
/*                                                                          */
/*  THE OUTPUT PARAMETER IS                                                 */
/*    zi  = DOUBLY-DIMENSIONED ARRAY OF DIMENSION (NXI,NYI),                */
/*          WHERE THE INTERPOLATED Z VALUES AT THE OUTPUT                   */
/*          GRID POINTS ARE TO BE STORED.                                   */
/*                                                                          */
/*  THE OTHER PARAMETERS ARE                                                */
/*    iwk = INTEGER ARRAY OF DIMENSION                                      */
/*             MAX0(31,27+ncp)*ndp+NXI*NYI                                  */
/*          USED INTERNALLY AS A WORK AREA,                                 */
/*    wk  = ARRAY OF DIMENSION 5*ndp USED INTERNALLY AS A                   */
/*          WORK AREA.                                                      */
/*                                                                          */
/*  USE OF A VALUE BETWEEN 3 AND 5 (INCLUSIVE) FOR ncp IS RECOM-            */
/*  MENDED UNLESS THERE ARE EVIDENCES THAT DICTATE OTHERWISE.               */
/*                                                                          */
/*  THIS SUBROUTINE CALLS THE IDCLDP, IDGRID, IDPDRV, IDPTIP, AND           */
/*  IDTANG SUBROUTINES.                                                     */
/*                                                                          */
/*  Note: MD is not used. We always assume MD = 1.                          */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */


int ak_idsfft(TDAContext *ctx, int ncp,int ndp,double * xd,double *yd,double *zd, int nxi,int nyi,double * xi,double *yi,double *zi,int *iwk,double *wk)
{

    int ncp0,ndp0,nxi0,nyi0,jwipt,jwipl,jwiwp,jwipc,jwngp0,jwigp0,nt,nl;
    int jig0mx,jig1mn,nngp,jngp,iti,il1,il2,jwngp,ngp0,jig0mn,jigp;
    int jwigp,ixi,iyi,izi,izi1,ngp1,jig1mx,jwiwl,err;

    /* SETTING OF SOME INPUT PARAMETERS TO LOCAL VARIABLES. */
                     
    ncp0 = ncp;
    ndp0 = ndp;
    nxi0 = nxi;
    nyi0 = nyi;

    if (ncp0 < 2 || ncp0 >= ndp0 || ndp0 < 4 || nxi0 < 1 || nyi0 < 1)
        return(-1);

    iwk[1] = ncp0;
    iwk[2] = ndp0;
    iwk[3] = nxi0;
    iwk[4] = nyi0;

    /* ALLOCATION OF STORAGE AREAS IN THE IWK ARRAY. */

    jwipt = 16;
    jwiwl = 6 * ndp0 + 1;
    jwngp0 = jwiwl - 1;
    jwipl = 24 * ndp0 + 1;
    jwiwp = 30 * ndp0 + 1;
    jwipc = 27 * ndp0 + 1;
    jwigp0 = imax(ctx, 31,27 + ncp0) * ndp0;

    /* TRIANGULATES THE X-Y PLANE. */                 

    err = ak_idtang(ctx, ndp0,xd,yd,&nt,iwk + jwipt - 1,&nl,iwk + jwipl - 1,
        iwk + jwiwl - 1,iwk + jwiwp - 1,wk);                     

    if (err)
        return(-1);

    iwk[5] = nt;
    iwk[6] = nl;
    if(nt == 0)
        return(-2);           

    /* DETERMINES ncp POINTS CLOSEST TO EACH DATA POINT. */

    err = ak_idcldp(ctx, ndp0,xd,yd,ncp0,iwk + jwipc - 1);
    if (err)
        return(-1);

    if (iwk[jwipc] == 0)
        return(-3);

    /* SORTS OUTPUT GRID POINTS IN ASCENDING ORDER OF THE TRIANGLE
       NUMBER AND THE BORDER LINE SEGMENT NUMBER. */

    ak_idgrid(ctx, xd,yd,nt,iwk + jwipt - 1,nl,iwk + jwipl - 1,nxi0,nyi0,
            xi,yi,iwk + jwngp0,iwk + jwigp0);
   
    /* ESTIMATES PARTIAL DERIVATIVES AT ALL DATA POINTS. */
         
    ak_idpdrv(ctx, ndp0,xd,yd,zd,ncp0,iwk + jwipc - 1,wk);
       
    /* INTERPOLATES THE ZI VALUES. */

    ctx->AK_itpv = 0;
    jig0mx = 0;
    jig1mn = nxi0 * nyi0 + 1;
    nngp = nt + 2 * nl;

    for (jngp = 1; jngp <= nngp; ++jngp) {
        iti = jngp; 
        if (jngp <= nt)
            goto L81;

        il1 = (jngp - nt + 1) / 2;
        il2 = (jngp - nt + 2) / 2;
        if (il2 > nl)
            il2 = 1;
        iti = il1 * (nt + nl) + il2;
L81:
        jwngp = jwngp0 + jngp;
        ngp0 = iwk[jwngp];
        if (ngp0 == 0)
            goto L86;

        jig0mn = jig0mx + 1;
        jig0mx = jig0mx + ngp0;

        for (jigp = jig0mn; jigp <= jig0mx; ++jigp) {
            jwigp = jwigp0 + jigp;
            izi = iwk[jwigp];
            iyi = (izi - 1) / nxi0 + 1;
            ixi = izi - nxi0 * (iyi-1);

            izi1 = (ixi - 1) * nyi0 + iyi;

            ak_idptip(ctx, xd,yd,zd,nt,iwk + jwipt - 1,nl,iwk + jwipl - 1,wk,
                iti,xi[ixi],yi[iyi],zi + izi1);                  
        }

L86:
        jwngp = jwngp0 + 2 * nngp + 1 - jngp;
        ngp1 = iwk[jwngp];
        if (ngp1 == 0)
            goto L89;

        jig1mx = jig1mn - 1;
        jig1mn = jig1mn - ngp1;

        for (jigp = jig1mn; jigp <= jig1mx; ++jigp) {
          jwigp = jwigp0 + jigp;
          izi = iwk[jwigp];
          iyi = (izi - 1) / nxi0 + 1;
          ixi = izi - nxi0 * (iyi - 1);

          izi1 = (ixi - 1) * nyi0 + iyi;

          ak_idptip(ctx, xd,yd,zd,nt,iwk + jwipt - 1,nl,iwk + jwipl - 1,wk,
                iti,xi[ixi],yi[iyi],zi + izi1);                  
        }
L89:    ;
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  ak_idtang(ndp,xd,yd,nt,ipt,nl,ipl,iwl,iwp,wk)                           */
/*                                                                          */
/*  SUBROUTINE  IDTANG(ndp,XD,YD,NT,IPT,NL,IPL,IWL,IWP,WK)                  */
/*                                                                          */
/*  THIS SUBROUTINE PERFORMS TRIANGULATION.  IT DIVIDES THE X-Y             */
/*  PLANE INTO A NUMBER OF TRIANGLES ACCORDING TO GIVEN DATA                */
/*  POINTS IN THE PLANE, DETERMINES LINE SEGMENTS THAT FORM THE             */
/*  BORDER OF DATA AREA, AND DETERMINES THE TRIANGLE NUMBERS                */
/*  CORRESPONDING TO THE BORDER LINE SEGMENTS.                              */
/*  AT COMPLETION, POINT NUMBERS OF THE VERTEXES OF EACH TRIANGLE           */
/*  ARE LISTED COUNTER-CLOCKWISE.  POINT NUMBERS OF THE END POINTS          */
/*  OF EACH BORDER LINE SEGMENT ARE LISTED COUNTER-CLOCKWISE,               */
/*  LISTING ORDER OF THE LINE SEGMENTS BEING COUNTER-CLOCKWISE.             */
/*                                                                          */
/*  THIS SUBROUTINE CALLS THE IDXCHG FUNCTION.                              */
/*                                                                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    ndp = NUMBER OF DATA POINTS,                                          */
/*    XD  = ARRAY OF DIMENSION ndp CONTAINING THE                           */
/*          X COORDINATES OF THE DATA POINTS,                               */
/*    YD  = ARRAY OF DIMENSION ndp CONTAINING THE                           */
/*          Y COORDINATES OF THE DATA POINTS.                               */
/*                                                                          */
/*  THE OUTPUT PARAMETERS ARE                                               */
/*    NT  = NUMBER OF TRIANGLES,                                            */
/*    IPT = INTEGER ARRAY OF DIMENSION 6*ndp-15, WHERE THE                  */
/*          POINT NUMBERS OF THE VERTEXES OF THE (IT)TH                     */
/*          TRIANGLE ARE TO BE STORED AS THE (3*IT-2)ND,                    */
/*          (3*IT-1)ST, AND (3*IT)TH ELEMENTS,                              */
/*          IT=1,2,...,NT,                                                  */
/*    NL  = NUMBER OF BORDER LINE SEGMENTS,                                 */
/*    IPL = INTEGER ARRAY OF DIMENSION 6*ndp, WHERE THE                     */
/*          POINT NUMBERS OF THE END POINTS OF THE (IL)TH                   */
/*          BORDER LINE SEGMENT AND ITS RESPECTIVE TRIANGLE                 */
/*          NUMBER ARE TO BE STORED AS THE (3*IL-2)ND,                      */
/*          (3*IL-1)ST, AND (3*IL)TH ELEMENTS,                              */
/*          IL=1,2,..., NL.                                                 */
/*                                                                          */
/*  THE OTHER PARAMETERS ARE                                                */
/*    IWL = INTEGER ARRAY OF DIMENSION 18*ndp USED                          */
/*          INTERNALLY AS A WORK AREA,                                      */
/*    IWP = INTEGER ARRAY OF DIMENSION ndp USED                             */
/*          INTERNALLY AS A WORK AREA,                                      */
/*    WK  = ARRAY OF DIMENSION ndp USED INTERNALLY AS A                     */
/*          WORK AREA.                                                      */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int ak_idtang(TDAContext *ctx, int ndp,double *xd,double *yd,int *nt,int *ipt,int *nl,int *ipl, int *iwl,int *iwp,double *wk)
{
    int ndp0,ndpm1,ipmn1,ipmn2,ip1,jp1,jpmn,its,ip,jp,jpmx,nt0,ntt3,nl0;
    int nlt3,nsh,nsht3,jp2t3,jp3t3,jwl,nln = 0,nlnt3 = 0,nlf,irep,ilf,ipl1,ipl2,nlft2;
    int ipti,ntt3p3,ilft2,ntf,itt3r,ipt1,ipt2,ipt3,it1t3,it2t3,iplj1,iplj2;
    int it,jwl1,jwl1mn,ip2,ip3,jp2,jpc,ipti1,ipti2,ip1p1,jlt3,nlfc,itt3,itf[4];
    double dsqmn,dsqi,dsq12,xdmp,ydmp,ar,armn,armx,dx21,dy21,x1,y1,dxmn,dymn;
    double dxmx,dymx,dsqmx,dx,dy;

    int nrep = 100;
    double ratio = 1.e-6;

    ndp0 = ndp;
    ndpm1 = ndp0 - 1;

    if (ndp0 < 4)
        return(-1);

    /* DETERMINES THE CLOSEST PAIR OF DATA POINTS AND THEIR MIDPOINT. */

/* L20 */
    dsqmn = ak_dsqf(ctx, xd[1],yd[1],xd[2],yd[2]);
    ipmn1 = 1;
    ipmn2 = 2;

    for (ip1 = 1; ip1 <= ndpm1; ++ip1) {
        x1 = xd[ip1];
        y1 = yd[ip1];
        ip1p1 = ip1 + 1;

        for (ip2 = ip1p1; ip2 <= ndp0; ++ip2) {
            dsqi = ak_dsqf(ctx, x1,y1,xd[ip2],yd[ip2]);
            if (dsqi == 0.0)
                return(-2);

            if (dsqi >= dsqmn)
                continue;     

            dsqmn = dsqi;
            ipmn1 = ip1;
            ipmn2 = ip2;
        }
    }
    dsq12 = dsqmn;
    xdmp = (xd[ipmn1] + xd[ipmn2]) / 2.0;
    ydmp = (yd[ipmn1] + yd[ipmn2]) / 2.0;

    /* SORTS THE OTHER (ndp-2) DATA POINTS IN ASCENDING ORDER OF
       DISTANCE FROM THE MIDPOINT AND STORES THE SORTED DATA POINT
       NUMBERS IN THE IWP ARRAY. */

/* L30 */
    jp1 = 2;
    for (ip1 = 1; ip1 <= ndp0; ++ip1) {
        if (ip1 == ipmn1 || ip1 == ipmn2)
            continue;
        jp1++;       
        iwp[jp1] = ip1;
        wk[jp1] = ak_dsqf(ctx, xdmp,ydmp,xd[ip1],yd[ip1]);
    }
    for (jp1 = 3; jp1 <= ndpm1; ++jp1) {
        dsqmn = wk[jp1];
        jpmn = jp1;
        for (jp2 = jp1; jp2 <= ndp0; ++jp2) {
            if (wk[jp2] >= dsqmn)
                continue;
            dsqmn = wk[jp2];
            jpmn = jp2;
        }
        its = iwp[jp1];
        iwp[jp1] = iwp[jpmn];
        iwp[jpmn] = its;
        wk[jpmn] = wk[jp1];
    }
 
    /* IF NECESSarY, MODIFIES THE ORDERING IN SUCH A WAY THAT THE
       FIRST THREE DATA POINTS arE NOT COLLINEar. */

/* L35 */
    ar = dsq12 * ratio;

    x1 = xd[ipmn1];
    y1 = yd[ipmn1];
    dx21 = xd[ipmn2] - x1;
    dy21 = yd[ipmn2] - y1;

    for (jp = 3; jp <= ndp0; ++jp) {
        ip = iwp[jp];
        if (fabs((yd[ip] - y1) * dx21 - (xd[ip] - x1) * dy21) > ar)
            goto L37;
    }
    return(-3);
                 
L37:
    if (jp == 3)
        goto L40;
    jpmx = jp;
    jp = jpmx + 1;

    for (jpc = 4; jpc <= jpmx; ++jpc) {
        jp--;       
        iwp[jp] = iwp[jp-1];
    }
    iwp[3] = ip;

    /* FORMS THE FIRST TRIANGLE.  STORES POINT NUMBERS OF THE VER-
       TEXES OF THE TRIANGLE IN THE ipT ARRAY, AND STORES POINT NUM-
       BERS OF THE BORDER LINE SEGMENTS AND THE TRIANGLE NUMBER IN
       THE ipl ARRAY. */

L40:
    ip1 = ipmn1;
    ip2 = ipmn2;
    ip3 = iwp[3];
    if (ak_side(ctx, xd[ip1],yd[ip1],xd[ip2],yd[ip2],xd[ip3],yd[ip3]) >= 0.0)
        goto L41;

    ip1 = ipmn2;
    ip2 = ipmn1;
L41:
    nt0 = 1;
    ntt3 = 3;
    ipt[1] = ip1;
    ipt[2] = ip2;
    ipt[3] = ip3;
    nl0 = 3;
    nlt3 = 9;
    ipl[1] = ip1;
    ipl[2] = ip2;
    ipl[3] = 1;
    ipl[4] = ip2;
    ipl[5] = ip3;
    ipl[6] = 1;
    ipl[7] = ip3;
    ipl[8] = ip1;
    ipl[9] = 1;

    /* ADDS THE REMAINING (ndp-3) DATA POINTS, ONE BY ONE. */

/* L50 */
    for (jp1 = 4; jp1 <= ndp0; ++jp1) {
        ip1 = iwp[jp1];
        x1 = xd[ip1];
        y1 = yd[ip1];

        /* DETERMINES THE VISIBLE BORDER LINE SEGMENTS. */

        ip2 = ipl[1]; 
        jpmn = 1;
        dxmn = xd[ip2] - x1;
        dymn = yd[ip2] - y1;
        dsqmn = pow(dxmn,2.0) + pow(dymn,2.0);
        armn= dsqmn * ratio;
        jpmx = 1;
        dxmx = dxmn;
        dymx = dymn;
        dsqmx = dsqmn;
        armx = armn;

        for (jp2 = 2; jp2 <= nl0; ++jp2) {
            ip2 = ipl[3 * jp2 - 2];
            dx = xd[ip2] - x1;
            dy = yd[ip2] - y1;
            ar = dy * dxmn - dx * dymn;
            if (ar > armn)
                goto L51;

            dsqi= dx* dx + dy * dy;
            if (ar >= (-armn) && dsqi >= dsqmn)
                goto L51;

            jpmn = jp2;
            dxmn = dx;
            dymn = dy;
            dsqmn = dsqi;
            armn = dsqmn * ratio;
L51:
            ar = dy * dxmx - dx * dymx;
            if (ar < (-armx))
                continue;

            dsqi = dx * dx + dy * dy;
            if (ar <= armx && dsqi >= dsqmx)
                continue;                   

            jpmx = jp2;
            dxmx = dx;
            dymx = dy;
            dsqmx = dsqi;
            armx = dsqmx * ratio;
        }
        if (jpmx < jpmn)
            jpmx = jpmx + nl0;

        nsh = jpmn - 1;
        if (nsh <= 0)
            goto L60;

        /* SHIFTS (ROTATES) THE ipl ARRAY TO HAVE THE INVISIBLE BORDER
           LINE SEGMENTS CONTAINED IN THE FIRST ParT OF THE ipl ARRAY. */

        nsht3 = nsh * 3;

        for (jp2t3 = 3; jp2t3 <= nsht3; jp2t3 += 3) {
            jp3t3 = jp2t3 + nlt3;
            ipl[jp3t3 - 2] = ipl[jp2t3 - 2];
            ipl[jp3t3 - 1] = ipl[jp2t3 - 1];
            ipl[jp3t3]     = ipl[jp2t3];
        }
        for (jp2t3 = 3; jp2t3 <= nlt3; jp2t3 += 3) {
            jp3t3 = jp2t3 + nsht3;
            ipl[jp2t3 - 2] = ipl[jp3t3 - 2];
            ipl[jp2t3 - 1] = ipl[jp3t3 - 1];
            ipl[jp2t3]     = ipl[jp3t3];
        }
        jpmx = jpmx - nsh;

        /* ADDS TRIANGLES TO THE ipt ARRAY, UPDATES BORDER LINE
           SEGMENTS IN THE ipl ARRAY, AND SETS FLAGS FOR THE BORDER
           LINE SEGMENTS TO BE REEXAMINED IN THE IWL ARRAY. */

L60:    jwl = 0;
        for (jp2 = jpmx; jp2 <= nl0; ++jp2) {
            jp2t3 = jp2 * 3;
            ipl1 = ipl[jp2t3 - 2];
            ipl2 = ipl[jp2t3 - 1];
            it   = ipl[jp2t3];

            /* ADDS A TRIANGLE TO THE ipt ARRAY. */

            nt0 = nt0 + 1;
            ntt3 = ntt3 + 3;
            ipt[ntt3 - 2] = ipl2;
            ipt[ntt3 - 1] = ipl1;
            ipt[ntt3]     = ip1;

            /* UPDATES BORDER LINE SEGMENTS IN THE ipl ARRAY. */

            if (jp2 != jpmx)
                goto L61;

            ipl[jp2t3 - 1] = ip1;
            ipl[jp2t3]     = nt0;
L61:
            if (jp2 != nl0)
                goto L62;

            nln = jpmx + 1;
            nlnt3 = nln * 3;
            ipl[nlnt3 - 2] = ip1;
            ipl[nlnt3 - 1] = ipl[1];
            ipl[nlnt3]     = nt0;

            /* DETERMINES THE VERTEX THAT DOES NOT LIE ON THE BORDER
               LINE SEGMENTS. */

L62:
            itt3 = it * 3;
            ipti = ipt[itt3 - 2];
            if (ipti != ipl1 && ipti != ipl2)
                goto L63;

            ipti = ipt[itt3 - 1];
            if (ipti != ipl1 && ipti != ipl2)
                goto L63;
            ipti = ipt[itt3];

            /* CHECKS IF THE EXCHANGE IS NECESSARY. */

L63:
            if (ak_idxchg(ctx, xd,yd,ip1,ipti,ipl1,ipl2) == 0)  
                continue;

            /* MODIFIES THE ipt ARRAY WHEN NECESSARY. */

            ipt[itt3 - 2] = ipti;
            ipt[itt3 - 1] = ipl1;
            ipt[itt3]     = ip1;
            ipt[ntt3 - 1] = ipti;
            if (jp2 == jpmx)
                ipl[jp2t3] = it; 

            if (jp2 == nl0 && ipl[3] == it)
                ipl[3] = nt0;

            /* SETS FLAGS IN THE IWL ARRAY. */

            jwl = jwl + 4;
            iwl[jwl - 3] = ipl1;
            iwl[jwl - 2] = ipti;
            iwl[jwl - 1] = ipti;
            iwl[jwl]     = ipl2;
        }
        nl0 = nln;
        nlt3 = nlnt3;
        nlf = jwl / 2;
        if (nlf == 0)
            goto L79;

        /* IMPROVES TRIANGULATION. */

/* L70 */
        ntt3p3 = ntt3 + 3;
        for (irep = 1; irep <= nrep; ++irep) {

            for (ilf = 1; ilf <= nlf; ++ilf) {
                ilft2 = ilf * 2;
                ipl1 = iwl[ilft2 - 1];
                ipl2 = iwl[ilft2];

                /* LOCATES IN THE IPT ARRAY TWO TRIANGLES ON BOTH SIDES OF
                   THE FLAGGED LINE SEGMENT. */

                ntf = 0;
                for (itt3r = 3; itt3r <= ntt3; itt3r += 3) {
                    itt3 = ntt3p3 - itt3r;
                    ipt1 = ipt[itt3 - 2];
                    ipt2 = ipt[itt3 - 1];
                    ipt3 = ipt[itt3];
                    if (ipl1 != ipt1 && ipl1 != ipt2 && ipl1 != ipt3)
                        continue;
                    if (ipl2 != ipt1 && ipl2 != ipt2 && ipl2 != ipt3)
                        continue;                                     
                    ntf = ntf + 1;
                    itf[ntf] = itt3 / 3;
                    if (ntf == 2)
                        goto L72;
                }
                if (ntf < 2)
                    goto L76;

                /* DETERMINES THE VERTEXES OF THE TRIANGLES THAT DO NOT LIE
                   ON THE LINE SEGMENT. */

L72:
                it1t3 = itf[1] * 3;
                ipti1 = ipt[it1t3 - 2];
                if (ipti1 != ipl1 && ipti1 != ipl2)
                    goto L73;
                ipti1 = ipt[it1t3 - 1];
                if (ipti1 != ipl1 && ipti1 != ipl2)
                    goto L73;
                ipti1 = ipt[it1t3];
L73:
                it2t3 = itf[2] * 3;
                ipti2 = ipt[it2t3 - 2];
                if (ipti2 != ipl1 && ipti2 != ipl2)
                    goto L74;
                ipti2 = ipt[it2t3 - 1];
                if (ipti2 != ipl1 && ipti2 != ipl2)
                    goto L74;
                ipti2 = ipt[it2t3];

                /* CHECKS IF THE EXCHANGE IS NECESSARY. */

L74:            if (ak_idxchg(ctx, xd,yd,ipti1,ipti2,ipl1,ipl2) == 0)
                    goto L76;

                /* MODIFIES THE ipt ARRAY WHEN NECESSARY. */

                ipt[it1t3 - 2] = ipti1;
                ipt[it1t3 - 1] = ipti2;
                ipt[it1t3]     = ipl1;
                ipt[it2t3 - 2] = ipti2;
                ipt[it2t3 - 1] = ipti1;
                ipt[it2t3]     = ipl2;

                /* SETS NEW FLAGS. */

                jwl = jwl + 8;
                iwl[jwl - 7] = ipl1;
                iwl[jwl - 6] = ipti1;
                iwl[jwl - 5] = ipti1;
                iwl[jwl - 4] = ipl2;
                iwl[jwl - 3] = ipl2;
                iwl[jwl - 2] = ipti2;
                iwl[jwl - 1] = ipti2;
                iwl[jwl]     = ipl1;

                for (jlt3 = 3; jlt3 <= nlt3; jlt3 += 3) {
                    iplj1 = ipl[jlt3 - 2];
                    iplj2 = ipl[jlt3 - 1];
                    if ((iplj1 == ipl1 && iplj2 == ipti2) || 
                        (iplj2 == ipl1 && iplj1 == ipti2))
                        ipl[jlt3] = itf[1];

                    if ((iplj1 == ipl2 && iplj2 == ipti1) || 
                        (iplj2 == ipl2 && iplj1 == ipti1))
                        ipl[jlt3] = itf[2];
                }
L76:            ;  
            }
            nlfc = nlf;
            nlf = jwl / 2;
            if (nlf == nlfc)
                goto L79;

            /* RESETS THE IWL ARRAY FOR THE NEXT ROUND. */

            jwl = 0;
            jwl1mn = (nlfc + 1) * 2;
            nlft2 = nlf * 2;

            for (jwl1 = jwl1mn; jwl1 <= nlft2; jwl1 += 2) {
                jwl = jwl + 2;
                iwl[jwl - 1] = iwl[jwl1 - 1];
                iwl[jwl]     = iwl[jwl1];
            }
            nlf = jwl / 2;
/* L78 */        ;
        }

L79:    ;          
    }

    /* REARRANGES THE ipt ARRAY SO THAT THE VERTEXES OF EACH TRIANGLE
       ARE LISTED COUNTER-CLOCKWISE. */

/* L80 */
    for (itt3 = 3; itt3 <= ntt3; itt3 += 3) {
        ip1 = ipt[itt3 - 2];
        ip2 = ipt[itt3 - 1];
        ip3 = ipt[itt3];
        if (ak_side(ctx, xd[ip1],yd[ip1],xd[ip2],yd[ip2],xd[ip3],yd[ip3]) >= 0.0)
            continue;                                 
        ipt[itt3 - 2] = ip2;
        ipt[itt3 - 1] = ip1;
    }
    *nt = nt0;
    *nl = nl0;
    return(0);
}

double ak_dsqf(TDAContext *ctx, double u1,double v1,double u2,double v2)
{
    (void)ctx;        /* unused: the signature is shared */
    double r;
    r = pow(u2 - u1,2.0) + pow(v2 - v1,2.0);
    return(r);
}

double ak_side(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3)
{
    (void)ctx;        /* unused: the signature is shared */
    double r;
    r = (v3 - v1) * (u2 - u1) - (u3 - u1) * (v2 - v1);
    return(r);
}

/* -###-------------------------------------------------------------------- */
/*  ak_idxchg(x,y,i1,i2,i3,i4)                                              */
/*                                                                          */
/*  FUNCTION  IDXCHG(X,Y,I1,I2,I3,I4)                                       */
/*                                                                          */
/*  THIS FUNCTION DETERMINES WHETHER OR NOT THE EXCHANGE OF TWO             */
/*  TRIANGLES IS NECESSARY ON THE BASIS OF MAX-MIN-ANGLE CRITERION          */
/*  BY C. L. LAWSON.                                                        */
/*                                                                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    X,Y = ARRAYS CONTAINING THE COORDINATES OF THE DATA                   */
/*          POINTS,                                                         */
/*    I1,I2,I3,I4 = POINT NUMBERS OF FOUR POINTS P1, P2,                    */
/*          P3, AND P4 THAT FORM A QUADRILATERAL WITH P3                    */
/*          AND P4 CONNECTED DIAGONALLY.                                    */
/*                                                                          */
/*  THIS FUNCTION RETURNS AN INTEGER VALUE 1 (ONE) WHEN AN EX-              */
/*  CHANGE IS NECESSARY, AND 0 (ZERO) OTHERWISE.                            */

int ak_idxchg(TDAContext *ctx, double *x,double *y,int i1,int i2,int i3,int i4)  
{
    int idx;
    double x1,x2,x3,x4,y1,y2,y3,y4,u1,u2,u3,u4;
    double c1sq,b2sq,a1sq,b1sq,a2sq,c3sq,s1sq,s2sq,s3sq,s4sq;

    x1 = x[i1];
    y1 = y[i1];
    x2 = x[i2];
    y2 = y[i2];
    x3 = x[i3];
    y3 = y[i3];
    x4 = x[i4];
    y4 = y[i4];
               
    idx = 0;
    u3 = (y2 - y3) * (x1 - x3) - (x2 - x3) * (y1 - y3);
    u4 = (y1 - y4) * (x2 - x4) - (x1 - x4) * (y2 - y4);
    if(u3 * u4 <= 0.0)
        goto L30;

    u1 = (y3 - y1) * (x4 - x1) - (x3 - x1) * (y4 - y1);
    u2 = (y4 - y2) * (x3 - x2) - (x4 - x2) * (y3 - y2);

    a1sq = pow((x1 - x3),2.0) + pow((y1 - y3),2.0);
    b1sq = pow((x4 - x1),2.0) + pow((y4 - y1),2.0);
    c1sq = pow((x3 - x4),2.0) + pow((y3 - y4),2.0);
    a2sq = pow((x2 - x4),2.0) + pow((y2 - y4),2.0);
    b2sq = pow((x3 - x2),2.0) + pow((y3 - y2),2.0);
    c3sq = pow((x2 - x1),2.0) + pow((y2 - y1),2.0);

    s1sq = u1 * u1 / (c1sq * dmax(ctx, a1sq,b1sq));
    s2sq = u2 * u2 / (c1sq * dmax(ctx, a2sq,b2sq));
    s3sq = u3 * u3 / (c3sq * dmax(ctx, b2sq,a1sq));
    s4sq = u4 * u4 / (c3sq * dmax(ctx, b1sq,a2sq));

    if(dmin(ctx, s1sq,s2sq) <  dmax(ctx, s3sq,s4sq))
        idx = 1;

L30:
    return(idx);
}

/* -###-------------------------------------------------------------------- */
/*  ak_idptip(xd,yd,zd,nt,ipt,nl,ipl,pdd,iti,xii,yii,zii)                   */
/*                                                                          */
/*  SUBROUTINE  IDPTIP(XD,YD,ZD,NT,IPT,NL,IPL,PDD,ITI,XII,YII,ZII)          */
/*                                                                          */
/*  THIS SUBROUTINE PERFORMS PUNCTUAL INTERPOLATION OR EXTRAPOLA-           */
/*  TION, I.E., DETERMINES THE Z VALUE AT A POINT.                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    XD,YD,ZD = ARRAYS OF DIMENSION ndp CONTAINING THE X,                  */
/*          Y, AND Z COORDINATES OF THE DATA POINTS, WHERE                  */
/*          ndp IS THE NUMBER OF THE DATA POINTS,                           */
/*    NT  = NUMBER OF TRIANGLES,                                            */
/*    IPT = INTEGER ARRAY OF DIMENSION 3*NT CONTAINING THE                  */
/*          POINT NUMBERS OF THE VERTEXES OF THE TRIANGLES,                 */
/*    NL  = NUMBER OF BORDER LINE SEGMENTS,                                 */
/*    IPL = INTEGER ARRAY OF DIMENSION 3*NL CONTAINING THE                  */
/*          POINT NUMBERS OF THE END POINTS OF THE BORDER                   */
/*          LINE SEGMENTS AND THEIR RESPECTIVE TRIANGLE                     */
/*          NUMBERS,                                                        */
/*    PDD = ARRAY OF DIMENSION 5*ndp CONTAINING THE PARTIAL                 */
/*          DERIVATIVES AT THE DATA POINTS,                                 */
/*    ITI = TRIANGLE NUMBER OF THE TRIANGLE IN WHICH LIES                   */
/*          THE POINT FOR WHICH INTERPOLATION IS TO BE                      */
/*          PERFORMED,                                                      */
/*    XII,YII = X AND Y COORDINATES OF THE POINT FOR WHICH                  */
/*          INTERPOLATION IS TO BE PERFORMED.                               */
/*                                                                          */
/*  THE OUTPUT PARAMETER IS                                                 */
/*    ZII = INTERPOLATED Z VALUE.                                           */

void ak_idptip(TDAContext *ctx, double *xd,double *yd,double *zd,int nt,int *ipt,int nl, int *ipl,double *pdd,int iti,double xii,double yii,double *zii)
{
    int i = 0,idp = 0,ntl = 0,it0 = 0,jpdd = 0,kpd = 0,jipl = 0,jpd = 0,il1 = 0,il2 = 0,jipt = 0;
    double lu = 0.0,lv = 0.0,adbc = 0.0,aa = 0.0,ab = 0.0,ad = 0.0,ac = 0.0,bb = 0.0,bc = 0.0,bdt2 = 0.0,cc = 0.0,cd = 0.0,dx = 0.0,dy = 0.0,dd = 0.0;
    double p0 = 0.0,p1 = 0.0,p2 = 0.0,p3 = 0.0,p4 = 0.0,g1 = 0.0,g2 = 0.0;
    double h1 = 0.0,h2 = 0.0,h3 = 0.0,a = 0.0,b = 0.0,c = 0.0,d = 0.0,u = 0.0,v = 0.0,dlt = 0.0,thxu = 0.0,thuv = 0.0,csuv = 0.0,thus = 0.0,thsv = 0.0,act2 = 0.0;
    double x[4],y[4],z[4] = {0.0},pd[16],zu[4],zv[4],zuu[4],zuv[4],zvv[4];

/*
 AK_itpv = 0;

*/

    it0 = iti;
    ntl = nt + nl;
    if (it0 <= ntl)
        goto L20;

    il1 = it0 / ntl;
    il2 = it0 - il1 * ntl;
    if (il1 == il2)
        goto L40;

    goto L60;

    /* CALCULATION OF ZII BY INTERPOLATION.
       CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED. */

L20:
    if (it0 == ctx->AK_itpv)
        goto L30;

    /* LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE VERTEXES. */

/* L21 */
    jipt = 3 * (it0 - 1);
    jpd = 0;

    for (i = 1; i <= 3; ++i) {
        jipt = jipt + 1;
        idp = ipt[jipt];
        x[i] = xd[idp];
        y[i] = yd[idp];
        z[i] = zd[idp];
        jpdd = 5 * (idp - 1);

        for (kpd = 1; kpd <= 5; ++kpd) {
            jpd = jpd + 1;
            jpdd = jpdd + 1;
            pd[jpd] = pdd[jpdd];
        }
    }

    /* DETERMINES THE COEFFICIENTS FOR THE COORDINATE SYSTEM
       TRANSFORMATION FROM THE X-Y SYSTEM TO THE U-V SYSTEM AND VICE VERSA. */

/* L24 */
    ctx->s_ak_idptip_x0 = x[1];
    ctx->s_ak_idptip_y0 = y[1];
    a = x[2] - ctx->s_ak_idptip_x0;
    b = x[3] - ctx->s_ak_idptip_x0;
    c = y[2] - ctx->s_ak_idptip_y0;
    d = y[3] - ctx->s_ak_idptip_y0;
    ad = a * d;
    bc = b * c;
    dlt = ad - bc;
    ctx->s_ak_idptip_ap =  d / dlt;
    ctx->s_ak_idptip_bp = -b / dlt;
    ctx->s_ak_idptip_cp = -c / dlt;
    ctx->s_ak_idptip_dp =  a / dlt;

    /* CONVERTS THE PARTIAL DERIVATIVES AT THE VERTEXES OF THE
       TRIANGLE FOR THE U-V COORDINATE SYSTEM. */

/* L25 */
    aa = a * a;
    act2 = 2.0 * a * c;
    cc = c * c;
    ab = a * b;
    adbc = ad + bc;
    cd = c * d;
    bb = b * b;
    bdt2 = 2.0 * b * d;
    dd = d * d;

    for (i = 1; i <= 3; ++i) {
        jpd = 5  *  i;
        zu[i] = a * pd[jpd - 4] + c * pd[jpd - 3];
        zv[i] = b * pd[jpd - 4] + d * pd[jpd - 3];
        zuu[i] = aa * pd[jpd - 2] + act2 * pd[jpd - 1] + cc * pd[jpd];
        zuv[i] = ab * pd[jpd - 2] + adbc * pd[jpd - 1] + cd * pd[jpd];
        zvv[i] = bb * pd[jpd - 2] + bdt2 * pd[jpd - 1] + dd * pd[jpd];
    }

    /*  CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL.  */

/* L27 */
    ctx->s_ak_idptip_p00 = z[1];
    ctx->s_ak_idptip_p10 = zu[1];
    ctx->s_ak_idptip_p01 = zv[1];
    ctx->s_ak_idptip_p20 = 0.5 * zuu[1];
    ctx->s_ak_idptip_p11 = zuv[1];
    ctx->s_ak_idptip_p02 = 0.5 * zvv[1];
    h1 = z[2] - ctx->s_ak_idptip_p00 - ctx->s_ak_idptip_p10 - ctx->s_ak_idptip_p20;
    h2 = zu[2] - ctx->s_ak_idptip_p10 - zuu[1];
    h3 = zuu[2] - zuu[1];
    ctx->s_ak_idptip_p30 =  10.0 * h1 - 4.0 * h2 + 0.5 * h3;
    ctx->s_ak_idptip_p40 = -15.0 * h1 + 7.0 * h2 - h3;
    ctx->s_ak_idptip_p5 =   6.0 * h1 - 3.0 * h2 + 0.5 * h3;
    h1 = z[3] - ctx->s_ak_idptip_p00 - ctx->s_ak_idptip_p01 - ctx->s_ak_idptip_p02;
    h2 = zv[3] - ctx->s_ak_idptip_p01 - zvv[1];
    h3 = zvv[3] - zvv[1];
    ctx->s_ak_idptip_p03 =  10.0 * h1 - 4.0 * h2 + 0.5 * h3;
    ctx->s_ak_idptip_p04 = -15.0 * h1 + 7.0 * h2 - h3;
    ctx->s_ak_idptip_p05 =   6.0 * h1 - 3.0 * h2 + 0.5 * h3;
    lu = sqrt(aa + cc);
    lv = sqrt(bb + dd);
    thxu = atan2(c,a);
    thuv = atan2(d,b) - thxu;
    csuv = cos(thuv);
    ctx->s_ak_idptip_p41 = 5.0 * lv * csuv / lu * ctx->s_ak_idptip_p5;
    ctx->s_ak_idptip_p14 = 5.0 * lu * csuv / lv * ctx->s_ak_idptip_p05;
    h1 = zv[2] - ctx->s_ak_idptip_p01 - ctx->s_ak_idptip_p11 - ctx->s_ak_idptip_p41;
    h2 = zuv[2] - ctx->s_ak_idptip_p11 - 4.0 * ctx->s_ak_idptip_p41;
    ctx->s_ak_idptip_p21 =  3.0 * h1 - h2;
    ctx->s_ak_idptip_p31 = -2.0 * h1 + h2;
    h1 = zu[3] - ctx->s_ak_idptip_p10 - ctx->s_ak_idptip_p11 - ctx->s_ak_idptip_p14;
    h2 = zuv[3] - ctx->s_ak_idptip_p11 - 4.0 * ctx->s_ak_idptip_p14;
    ctx->s_ak_idptip_p12 =  3.0 * h1 - h2;
    ctx->s_ak_idptip_p13 = -2.0 * h1 + h2;
    thus = atan2(d - c,b - a) - thxu;
    thsv = thuv - thus;
    aa =  sin(thsv) / lu;
    bb = -cos(thsv) / lu;
    cc =  sin(thus) / lv;
    dd =  cos(thus) / lv;
    ac = aa * cc;
    ad = aa * dd;
    bc = bb * cc;
    g1 = aa * ac * (3.0 * bc + 2.0 * ad);
    g2 = cc * ac * (3.0 * ad + 2.0 * bc);
    h1 = -aa * aa * aa * (5.0 * aa * bb * ctx->s_ak_idptip_p5  + (4.0 * bc + ad) * ctx->s_ak_idptip_p41) - 
          cc * cc * cc * (5.0 * cc * dd * ctx->s_ak_idptip_p05 + (4.0 * ad + bc) * ctx->s_ak_idptip_p14);
    h2 = 0.5 * zvv[2] - ctx->s_ak_idptip_p02 - ctx->s_ak_idptip_p12;
    h3 = 0.5 * zuu[3] - ctx->s_ak_idptip_p20 - ctx->s_ak_idptip_p21;
    ctx->s_ak_idptip_p22 = (g1 * h2 + g2 * h3 - h1)/(g1 + g2);
    ctx->s_ak_idptip_p32 = h2 - ctx->s_ak_idptip_p22;
    ctx->s_ak_idptip_p23 = h3 - ctx->s_ak_idptip_p22;
    ctx->AK_itpv = it0;

    /*  CONVERTS xii AND yii TO U-V SYSTEM.  */

L30:
    dx = xii - ctx->s_ak_idptip_x0;
    dy = yii - ctx->s_ak_idptip_y0;
    u = ctx->s_ak_idptip_ap * dx + ctx->s_ak_idptip_bp * dy;
    v = ctx->s_ak_idptip_cp * dx + ctx->s_ak_idptip_dp * dy;

    /*  EVALUATES THE POLYNOMIAL.  */

/* L31 */
    p0 = ctx->s_ak_idptip_p00 + v * (ctx->s_ak_idptip_p01 + v * (ctx->s_ak_idptip_p02 + v * (ctx->s_ak_idptip_p03 + v * (ctx->s_ak_idptip_p04 + v * ctx->s_ak_idptip_p05))));
    p1 = ctx->s_ak_idptip_p10 + v * (ctx->s_ak_idptip_p11 + v * (ctx->s_ak_idptip_p12 + v * (ctx->s_ak_idptip_p13 + v * ctx->s_ak_idptip_p14)));
    p2 = ctx->s_ak_idptip_p20 + v * (ctx->s_ak_idptip_p21 + v * (ctx->s_ak_idptip_p22 + v * ctx->s_ak_idptip_p23));
    p3 = ctx->s_ak_idptip_p30 + v * (ctx->s_ak_idptip_p31 + v * ctx->s_ak_idptip_p32);
    p4 = ctx->s_ak_idptip_p40 + v * ctx->s_ak_idptip_p41;
    *zii = p0 + u * (p1 + u * (p2 + u * (p3 + u * (p4 + u * ctx->s_ak_idptip_p5))));
    return;

    /* CALCULATION OF ZII BY EXTRapOLATION IN THE RECTANGLE.
       CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED. */

L40:
    if (it0 == ctx->AK_itpv)
        goto L50;

    /* LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE END
       POINTS OF THE BORDER LINE SEGMENT. */

/* L41 */
    jipl = 3 * (il1 - 1);
    jpd = 0;

    for (i = 1; i <= 2; ++i) {
        jipl = jipl + 1;
        idp = ipl[jipl];
        x[i] = xd[idp];
        y[i] = yd[idp];
        z[i] = zd[idp];
        jpdd = 5 * (idp - 1);

        for (kpd = 1; kpd <= 5; ++kpd) {
            jpd = jpd + 1;
            jpdd = jpdd + 1;
            pd[jpd] = pdd[jpdd];
        }
    }

    /* DETERMINES THE COEFFICIENTS FOR THE COORDINATE SYSTEM
       TRANSFORMATION FROM THE X-Y SYSTEM TO THE U-V SYSTEM AND VICE VERSA.  */

/* L44 */
    ctx->s_ak_idptip_x0 = x[1];
    ctx->s_ak_idptip_y0 = y[1];
    a = y[2] - y[1];
    b = x[2] - x[1];
    c = -b;
    d = a;
    ad = a * d;
    bc = b * c;
    dlt = ad - bc;
    ctx->s_ak_idptip_ap =  d / dlt;
    ctx->s_ak_idptip_bp = -b / dlt;
    ctx->s_ak_idptip_cp = -ctx->s_ak_idptip_bp;
    ctx->s_ak_idptip_dp =  ctx->s_ak_idptip_ap;

    /* CONVERTS THE PARTIAL DERIVATIVES AT THE END POINTS OF THE
       BORDER LINE SEGMENT FOR THE U-V COORDINATE SYSTEM.  */

/* L45 */
    aa = a * a;
    act2 = 2.0 * a * c;
    cc = c * c;
    ab = a * b;
    adbc = ad + bc;
    cd = c * d;
    bb = b * b;
    bdt2 = 2.0 * b * d;
    dd = d * d;

    for (i = 1; i <= 2; ++i) {
        jpd = 5 * i;
        zu[i] = a * pd[jpd - 4] + c * pd[jpd - 3];
        zv[i] = b * pd[jpd - 4] + d * pd[jpd - 3];
        zuu[i] = aa * pd[jpd - 2] + act2 * pd[jpd - 1] + cc * pd[jpd];
        zuv[i] = ab * pd[jpd - 2] + adbc * pd[jpd - 1] + cd * pd[jpd];
        zvv[i] = bb * pd[jpd - 2] + bdt2 * pd[jpd - 1] + dd * pd[jpd];
    }

    /* CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL. */

/* L47 */
    ctx->s_ak_idptip_p00 = z[1];
    ctx->s_ak_idptip_p10 = zu[1];
    ctx->s_ak_idptip_p01 = zv[1];
    ctx->s_ak_idptip_p20 = 0.5 * zuu[1];
    ctx->s_ak_idptip_p11 = zuv[1];
    ctx->s_ak_idptip_p02 = 0.5 * zvv[1];
    h1 = z[2] - ctx->s_ak_idptip_p00 - ctx->s_ak_idptip_p01 - ctx->s_ak_idptip_p02;
    h2 = zv[2] - ctx->s_ak_idptip_p01 - zvv[1];
    h3 = zvv[2] - zvv[1];
    ctx->s_ak_idptip_p03 =  10.0 * h1 - 4.0 * h2 + 0.5 * h3;
    ctx->s_ak_idptip_p04 = -15.0 * h1 + 7.0 * h2 - h3;
    ctx->s_ak_idptip_p05 =   6.0 * h1 - 3.0 * h2 + 0.5 * h3;
    h1 = zu[2] - ctx->s_ak_idptip_p10 - ctx->s_ak_idptip_p11;
    h2 = zuv[2] - ctx->s_ak_idptip_p11;
    ctx->s_ak_idptip_p12 =  3.0 * h1 - h2;
    ctx->s_ak_idptip_p13 = -2.0 * h1 + h2;
    ctx->s_ak_idptip_p21 = 0.0;
    ctx->s_ak_idptip_p23 = -zuu[2] + zuu[1];
    ctx->s_ak_idptip_p22 = -1.5 * ctx->s_ak_idptip_p23;
    ctx->AK_itpv = it0;

    /* CONVERTS xii AND yii TO U-V SYSTEM. */

L50:
    dx = xii - ctx->s_ak_idptip_x0;
    dy = yii - ctx->s_ak_idptip_y0;
    u = ctx->s_ak_idptip_ap * dx + ctx->s_ak_idptip_bp * dy;
    v = ctx->s_ak_idptip_cp * dx + ctx->s_ak_idptip_dp * dy;

    /* EVALUATES THE POLYNOMIAL. */

/* L51 */
    p0 = ctx->s_ak_idptip_p00 + v * (ctx->s_ak_idptip_p01 + v * (ctx->s_ak_idptip_p02 + v * (ctx->s_ak_idptip_p03 + v * (ctx->s_ak_idptip_p04 + v * ctx->s_ak_idptip_p05))));
    p1 = ctx->s_ak_idptip_p10 + v * (ctx->s_ak_idptip_p11 + v * (ctx->s_ak_idptip_p12 + v * ctx->s_ak_idptip_p13));
    p2 = ctx->s_ak_idptip_p20 + v * (ctx->s_ak_idptip_p21 + v * (ctx->s_ak_idptip_p22 + v * ctx->s_ak_idptip_p23));
    *zii = p0 + u * (p1 + u * p2);
    return;

    /* CALCULATION OF ZII BY EXTRAPOLATION IN THE TRIANGLE.
       CHECKS IF THE NECESSARY COEFFICIENTS HAVE BEEN CALCULATED. */

L60:
    if (it0 == ctx->AK_itpv)
        goto L70;

    /* LOADS COORDINATE AND PARTIAL DERIVATIVE VALUES AT THE VERTEX
       OF THE TRIANGLE. */

/* L61 */
    jipl = 3 * il2 - 2;
    idp = ipl[jipl];
    ctx->s_ak_idptip_x0 = xd[idp];
    ctx->s_ak_idptip_y0 = yd[idp];
    ctx->s_ak_idptip_z0 = zd[idp];
    jpdd = 5 * (idp - 1);

    for (kpd = 1; kpd <= 5; ++kpd) {
        jpdd = jpdd + 1;
        pd[kpd] = pdd[jpdd];
    }

    /* CALCULATES THE COEFFICIENTS OF THE POLYNOMIAL. */

/* L67 */
    ctx->s_ak_idptip_p00 = z[1];
    ctx->s_ak_idptip_p10 = pd[1];
    ctx->s_ak_idptip_p01 = pd[2];
    ctx->s_ak_idptip_p20 = 0.5 * pd[3];
    ctx->s_ak_idptip_p11 = pd[4];
    ctx->s_ak_idptip_p02 = 0.5 * pd[5];
    ctx->AK_itpv = it0;

    /* CONVERTS xii AND yii TO U-V SYSTEM. */

L70:
    u = xii - ctx->s_ak_idptip_x0;
    v = yii - ctx->s_ak_idptip_y0;

    /* EVALUATES THE POLYNOMIAL. */

/* L71 */
    p0 = ctx->s_ak_idptip_p00 + v * (ctx->s_ak_idptip_p01 + v * ctx->s_ak_idptip_p02);
    p1 = ctx->s_ak_idptip_p10 + v * ctx->s_ak_idptip_p11;
    *zii = p0 + u * (p1 + u * ctx->s_ak_idptip_p20);
    return;
}           

/* -#####------------------------------------------------------------------ */
/*  ak_idpdrv(ndp,xd,yd,zd,ncp,ipc,pd)                                      */
/*                                                                          */
/*  SUBROUTINE  IDPDRV(ndp,XD,YD,ZD,ncp,ipc,PD)                             */
/*                                                                          */
/*  THIS SUBROUTINE ESTIMATES PARTIAL DERIVATIVES OF THE FIRST AND          */
/*  SECOND ORDER AT THE DATA POINTS.                                        */
/*  THE INPUT PARAMETERS ARE                                                */
/*    ndp = NUMBER OF DATA POINTS,                                          */
/*    XD,YD,ZD = ARRAYS OF DIMENSION ndp CONTAINING THE X,                  */
/*          Y, AND Z COORDINATES OF THE DATA POINTS,                        */
/*    ncp = NUMBER OF ADDITIONAL DATA POINTS USED FOR ESTI-                 */
/*          MATING PARTIAL DERIVATIVES AT EACH DATA POINT,                  */
/*    ipc = INTEGER ARRAY OF DIMENSION ncp*ndp CONTAINING                   */
/*          THE POINT NUMBERS OF ncp DATA POINTS CLOSEST TO                 */
/*          EACH OF THE ndp DATA POINTS.                                    */
/*  THE OUTPUT PARAMETER IS                                                 */
/*    PD  = ARRAY OF DIMENSION 5*ndp, WHERE THE ESTIMATED                   */
/*          ZX, ZY, ZXX, ZXY, AND ZYY VALUES AT THE DATA                    */
/*          POINTS ARE TO BE STORED.                                        */
/*                                                                          */

void ak_idpdrv(TDAContext *ctx, int ndp,double *xd,double *yd,double *zd,int ncp, int *ipc,double *pd)
{
    (void)ctx;        /* unused: the signature is shared */
    int ndp0,ncp0,ncpm1,jipc,ipi,jpd,ip0,ic1,ic2,ic2mn,jpd0,jipc0;
    double nmx,nmy,nmz,nmxx,nmxy,nmyx,nmyy,dx1,dy1,dx2,dy2,zx0,zy0,x0,y0,z0;
    double dzx1,dzy1,dz1,dnmz,dz2,dnmx,dnmy,dzx2,dzy2;
    double dnmxx,dnmxy,dnmyx,dnmyy;

    ndp0 = ndp;
    ncp0 = ncp;
    ncpm1 = ncp0 - 1;

    /* ESTIMATION OF ZX AND ZY */ 

    for (ip0 = 1; ip0 <= ndp0; ++ip0) {

        x0 = xd[ip0];
        y0 = yd[ip0];
        z0 = zd[ip0];
        nmx = 0.0;
        nmy = 0.0;
        nmz = 0.0;
        jipc0 = ncp0 * (ip0 - 1);

        for (ic1 = 1; ic1 <= ncpm1; ++ic1) {
            jipc = jipc0 + ic1;
            ipi = ipc[jipc];
            dx1 = xd[ipi] - x0;
            dy1 = yd[ipi] - y0;
            dz1 = zd[ipi] - z0;
            ic2mn = ic1 + 1;

            for (ic2 = ic2mn; ic2 <= ncp0; ++ic2) {
                jipc = jipc0 + ic2;
                ipi = ipc[jipc];
                dx2 = xd[ipi] - x0;
                dy2 = yd[ipi] - y0;
                dnmz = dx1 * dy2 - dy1 * dx2;
                if (dnmz == 0.0)
                    goto L22;
                dz2 = zd[ipi] - z0;
                dnmx = dy1 * dz2 - dz1 * dy2;
                dnmy = dz1 * dx2 - dx1 * dz2;
                if (dnmz >= 0.0)
                    goto L21;
                dnmx =  -dnmx;
                dnmy =  -dnmy;
                dnmz =  -dnmz;
L21:
                nmx = nmx + dnmx;
                nmy = nmy + dnmy;
                nmz = nmz + dnmz;
L22:            ;
            }
        }
        jpd0 = 5 * ip0;
        pd[jpd0 - 4] =  -nmx / nmz;
        pd[jpd0 - 3] =  -nmy / nmz;
    }

    /* ESTIMATION OF ZXX, ZXY, AND ZYY */

/* L30 */
    for (ip0 = 1; ip0 <= ndp0; ++ip0) {
        jpd0 = jpd0 + 5;
        x0 = xd[ip0];
        jpd0 = 5 * ip0;
        y0 = yd[ip0];
        zx0 = pd[jpd0 - 4];
        zy0 = pd[jpd0 - 3];
        nmxx = 0.0;
        nmxy = 0.0;
        nmyx = 0.0;
        nmyy = 0.0;
        nmz  = 0.0;
        jipc0 = ncp0 * (ip0 - 1);

        for (ic1 = 1; ic1 <= ncpm1; ++ic1) {
            jipc = jipc0 + ic1;
            ipi = ipc[jipc];
            dx1 = xd[ipi] - x0;
            dy1 = yd[ipi] - y0;
            jpd = 5 * ipi;
            dzx1 = pd[jpd - 4] - zx0;
            dzy1 = pd[jpd - 3] - zy0;
            ic2mn = ic1 + 1;

            for (ic2 = ic2mn; ic2 <= ncp0; ++ic2) {
                jipc = jipc0 + ic2;
                ipi = ipc[jipc];
                dx2 = xd[ipi] - x0;
                dy2 = yd[ipi] - y0;
                dnmz  = dx1 * dy2  - dy1 * dx2;
                if (dnmz == 0.0)
                    goto L32;

                jpd = 5 * ipi;
                dzx2 = pd[jpd - 4] - zx0;
                dzy2 = pd[jpd - 3] - zy0;
                dnmxx = dy1 * dzx2 - dzx1 * dy2;
                dnmxy = dzx1 * dx2 - dx1 * dzx2;
                dnmyx = dy1 * dzy2 - dzy1 * dy2;
                dnmyy = dzy1 * dx2 - dx1 * dzy2;
                if (dnmz >= 0.0)
                    goto L31;  

                dnmxx =  -dnmxx;
                dnmxy =  -dnmxy;
                dnmyx =  -dnmyx;
                dnmyy =  -dnmyy;
                dnmz  =  -dnmz;
L31:
                nmxx = nmxx + dnmxx;
                nmxy = nmxy + dnmxy;
                nmyx = nmyx + dnmyx;                                            
                nmyy = nmyy + dnmyy;
                nmz  = nmz  + dnmz;
L32:            ;
            }
/* L33 */        ;  
        }
        pd[jpd0 - 2] =  -nmxx / nmz;
        pd[jpd0 - 1] =  -(nmxy + nmyx) / (2.0 * nmz);
        pd[jpd0]     =  -nmyy / nmz;
/* L34 */    ;          
    }
    return;
}           

/* -###-------------------------------------------------------------------- */
/*  ak_idgrid(xd,yd,nt,ipt,nl,ipl,nxi,nyi,xi,yi,ngp,igp)                    */
/*                                                                          */
/*  SUBROUTINE IDGRID(XD, YD, NT, IPT, NL, IPL, NXI, NYI, XI, YI,NGP, IGP)  */
/*                                                                          */
/*  THIS SUBROUTINE ORGANIZES GRID POINTS FOR SURFACE FITTING BY            */
/*  SORTING THEM IN ASCENDING ORDER OF TRIANGLE NUMBERS AND OF THE          */
/*  BORDER LINE SEGMENT NUMBER.                                             */
/*                                                                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    XD,YD = ARRAYS OF DIMENSION ndp CONTAINING THE X AND Y                */
/*          COORDINATES OF THE DATA POINTS, WHERE ndp IS THE                */
/*          NUMBER OF THE DATA POINTS,                                      */
/*    NT  = NUMBER OF TRIANGLES,                                            */
/*    IPT = INTEGER ARRAY OF DIMENSION 3*NT CONTAINING THE                  */
/*          POINT NUMBERS OF THE VERTEXES OF THE TRIANGLES,                 */
/*    NL  = NUMBER OF BORDER LINE SEGMENTS,                                 */
/*    IPL = INTEGER ARRAY OF DIMENSION 3*NL CONTAINING THE                  */
/*          POINT NUMBERS OF THE END POINTS OF THE BORDER                   */
/*          LINE SEGMENTS AND THEIR RESPECTIVE TRIANGLE                     */
/*          NUMBERS,                                                        */
/*    NXI = NUMBER OF GRID POINTS IN THE X COORDINATE,                      */
/*    NYI = NUMBER OF GRID POINTS IN THE Y COORDINATE,                      */
/*    XI,YI = ARRAYS OF DIMENSION NXI AND NYI CONTAINING                    */
/*          THE X AND Y COORDINATES OF THE GRID POINTS,                     */
/*          RESPECTIVELY.                                                   */
/*                                                                          */
/*  THE OUTPUT PARAMETERS ARE                                               */
/*    NGP = INTEGER ARRAY OF DIMENSION 2*(NT+2*NL) WHERE THE                */
/*          NUMBER OF GRID POINTS THAT BELONG TO EACH OF THE                */
/*          TRIANGLES OR OF THE BORDER LINE SEGMENTS ARE TO                 */
/*          BE STORED,                                                      */
/*    IGP = INTEGER ARRAY OF DIMENSION NXI*NYI WHERE THE                    */
/*          GRID POINT NUMBERS ARE TO BE STORED IN ASCENDING                */
/*          ORDER OF THE TRIANGLE NUMBER AND THE BORDER LINE                */
/*          SEGMENT NUMBER.                                                 */

void ak_idgrid(TDAContext *ctx, double *xd,double *yd,int nt,int *ipt,int nl,int *ipl, int nxi,int nyi,double *xi,double *yi,int *ngp,int *igp)
{
    int l,nt0,nl0,nxi0,nyi0,nxinyi,jngp0,jngp1,jigp0,jigp1,ngp0,ngp1,it0t3;
    int it0,il0,insd,ixi,iximx,ip1,ip2,ip3,iyi,izi,il0t3,ilp1t3,ilp1,iximn,jigp1i;
    double rr,xmn,ymn,xmx,ymx,ximn,yimn,ximx,yimx,x1,x2,x3,y1,y2,y3,xii,yii;

    nt0 = nt;
    nl0 = nl;
    nxi0 = nxi;
    nyi0 = nyi;
    nxinyi = nxi0 * nyi0;
    ximn = dmin(ctx, xi[1],xi[nxi0]);
    ximx = dmax(ctx, xi[1],xi[nxi0]);
    yimn = dmin(ctx, yi[1],yi[nyi0]);
    yimx = dmax(ctx, yi[1],yi[nyi0]);
 
    /* DETERMINES GRID POINTS INSIDE THE DATA AREA. */

    jngp0 = 0;
    jngp1 = 2 * (nt0 + 2 * nl0)  +  1;
    jigp0 = 0;
    jigp1 = nxinyi  +  1;

    for (it0 = 1; it0 <= nt0; ++it0) {
        ngp0 = 0;
        ngp1 = 0;
        it0t3 = it0 * 3;
        ip1 = ipt[it0t3 - 2];
        ip2 = ipt[it0t3 - 1];
        ip3 = ipt[it0t3];
        x1 = xd[ip1];
        y1 = yd[ip1];
        x2 = xd[ip2];
        y2 = yd[ip2];
        x3 = xd[ip3];
        y3 = yd[ip3]; 
        xmn = dmin(ctx, x1,dmin(ctx, x2,x3));
        xmx = dmax(ctx, x1,dmax(ctx, x2,x3));
        ymn = dmin(ctx, y1,dmin(ctx, y2,y3));
        ymx = dmax(ctx, y1,dmax(ctx, y2,y3));
        insd = 0;

        for (ixi = 1; ixi <= nxi0; ++ixi) {
            if (xi[ixi] >= xmn && xi[ixi] <= xmx)
                goto L10;

            if (insd == 0)
                continue;                
            iximx = ixi - 1;
            goto L30;
L10:
            if (insd == 1)
                continue;
            insd = 1;
            iximn = ixi;
        }

        if (insd == 0)
            goto L150;

        iximx = nxi0;

L30:
        for (iyi = 1; iyi <= nyi0; ++iyi) {
            yii = yi[iyi];
            if (yii < ymn || yii > ymx)
                goto L140;

            for (ixi = iximn; ixi <= iximx; ++ixi) {
                xii = xi[ixi];
                l = 0;
                if ((rr = ak_side1(ctx, x1,y1,x2,y2,xii,yii)) < 0.0)
                    goto L130;                                 
                else if (rr > 0.0)                             
                    goto L50;                                  
                l = 1;
L50:
                if ((rr = ak_side1(ctx, x2,y2,x3,y3,xii,yii)) < 0.0)
                    goto L130;                                 
                else if (rr > 0.0)
                    goto L70;
                l = 1;
L70:            if ((rr = ak_side1(ctx, x3,y3,x1,y1,xii,yii)) < 0.0)
                    goto L130;                          
                else if (rr > 0.0)
                    goto L90;
                l = 1;
L90:            izi = nxi0 * (iyi - 1)  +  ixi;

                if (l == 1)
                    goto L100;

                ngp0 = ngp0  +  1;
                jigp0 = jigp0 + 1;
                igp[jigp0] = izi;
                goto L130; 
L100:
                if (jigp1 > nxinyi)
                    goto L120;

                for (jigp1i = jigp1; jigp1i <= nxinyi; ++jigp1i) {
                    if (izi == igp[jigp1i])
                        goto L130;
                }
L120:
                ngp1 = ngp1 + 1;
                jigp1 = jigp1 - 1;
                igp[jigp1] = izi;
L130:           ;
            }
L140:       ;    
        }
L150:
        jngp0 = jngp0 + 1;
        ngp[jngp0] = ngp0;
        jngp1 = jngp1 - 1;
        ngp[jngp1] = ngp1;
/* L160 */   ;         
    }

    /* DETERMINES GRID POINTS OUTSIDE THE DATA AREA.
       IN SEMI - INFINITE RECTANGULAR AREA. */

    for (il0 = 1; il0 <= nl0; ++il0) {
        ngp0 = 0;
        ngp1 = 0;
        il0t3 = il0 * 3;
        ip1 = ipl[il0t3 - 2];
        ip2 = ipl[il0t3 - 1];
        x1 = xd[ip1];
        y1 = yd[ip1];
        x2 = xd[ip2];
        y2 = yd[ip2];
        xmn = ximn;
        xmx = ximx;
        ymn = yimn;
        ymx = yimx;

        if (y2 >= y1)
            xmn = dmin(ctx, x1,x2);
        if (y2 <= y1)
            xmx = dmax(ctx, x1,x2);
        if (x2 <= x1)
            ymn = dmin(ctx, y1,y2);
        if (x2 >= x1)
            ymx = dmax(ctx, y1,y2);
        insd = 0;

        for (ixi = 1; ixi <= nxi0; ++ixi) {
            if (xi[ixi] >= xmn && xi[ixi] <= xmx)
                goto L170;
            if (insd == 0)
                continue;  
            iximx = ixi - 1;
            goto L190;
L170:
            if (insd == 1)             
                continue;

            insd = 1;
            iximn = ixi;
        }
        if (insd == 0)
            goto L310;

        iximx = nxi0;
L190:
    
        for (iyi = 1; iyi <= nyi0; ++iyi) {
            yii = yi[iyi];
            if (yii < ymn || yii > ymx)
                continue;                                   

            for (ixi = iximn; ixi <= iximx; ++ixi) {
                xii = xi[ixi];
                l = 0;
                if ((rr = ak_side1(ctx, x1,y1,x2,y2,xii,yii)) < 0.0)
                    goto L210;                          
                else if (rr > 0.0)
                    goto L290;
                l = 1;
L210:   
                if ((rr = ak_spdt1(ctx, x2,y2,x1,y1,xii,yii)) < 0.0)
                    goto L290;                          
                else if (rr > 0.0)                      
                    goto L230;
                l = 1;
L230:
                if ((rr = ak_spdt1(ctx, x1,y1,x2,y2,xii,yii)) < 0.0)
                    goto L290;                          
                else if (rr > 0.0)                      
                    goto L250;                          
                l = 1;
L250:
                izi = nxi0 * (iyi - 1) + ixi;

                if (l == 1)
                    goto L260;
                ngp0 = ngp0 + 1;
                jigp0 = jigp0 + 1;
                igp[jigp0] = izi;
                goto L290;
L260:
                if (jigp1 > nxinyi)
                    goto L280;

                for (jigp1i = jigp1; jigp1i <= nxinyi; ++jigp1i) {
                    if (izi == igp[jigp1i])
                        goto L290;
                }
L280:
                ngp1 = ngp1 + 1;
                jigp1 = jigp1 - 1;
                igp[jigp1] = izi;
L290:           ;
            }
/* L300 */       ;  
        }

L310:
        jngp0 = jngp0 + 1;
        ngp[jngp0] = ngp0;
        jngp1 = jngp1 - 1;
        ngp[jngp1] = ngp1;

        /* IN SEMI - INFINITE TRIANGULAR AREA. */

        ngp0 = 0;
        ngp1 = 0;
        ilp1 = (il0 % nl0) + 1;
        ilp1t3 = ilp1 * 3;
        ip3 = ipl[ilp1t3 - 1];
        x3 = xd[ip3];
        y3 = yd[ip3];
        xmn = ximn;
        xmx = ximx;
        ymn = yimn;
        ymx = yimx;
        if (y3 >= y2 && y2 >= y1)
            xmn = x2;
        if (y3 <= y2 && y2 <= y1)
            xmx = x2;
        if (x3 <= x2 && x2 <= x1)
            ymn = y2;
        if (x3 >= x2 && x2 >= x1)
            ymx = y2;
        insd = 0;

        for (ixi = 1; ixi <= nxi0; ++ixi) {
            if (xi[ixi] >= xmn && xi[ixi] <= xmx)
                goto L320;

            if (insd == 0)
                continue;             
            iximx = ixi - 1;
            goto L340;
L320:
            if (insd == 1)              
                continue;

            insd = 1;
            iximn = ixi;
        }
        if (insd == 0)
            goto L440;

        iximx = nxi0;
L340:
        for (iyi = 1; iyi <= nyi0; ++iyi) {
            yii = yi[iyi];
            if (yii < ymn || yii > ymx)
                goto L430;

            for (ixi = iximn; ixi <= iximx; ++ixi) {
                xii = xi[ixi];
                l = 0;
                if ((rr = ak_spdt1(ctx, x1,y1,x2,y2,xii,yii)) < 0.0)
                    goto L360;                          
                else if (rr > 0.0)                      
                    goto L420;
                l = 1;
L360:
                if ((rr = ak_spdt1(ctx, x3,y3,x2,y2,xii,yii)) < 0.0)
                    goto L380;                          
                else if (rr > 0.0)                      
                    goto L420;
                l = 1;
L380:
                izi = nxi0 * (iyi - 1) + ixi;

                if (l == 1)
                    goto L390;
                ngp0 = ngp0 + 1;
                jigp0 = jigp0 + 1;
                igp[jigp0] = izi;
                goto L420;
L390:
                if (jigp1 > nxinyi)
                    goto L410;

                for (jigp1i = jigp1; jigp1i <= nxinyi; ++jigp1i) {
                    if (izi == igp[jigp1i])
                        goto L420;
                }
L410:
                ngp1 = ngp1 + 1;
                jigp1 = jigp1 - 1;
                igp[jigp1] = izi;
L420:           ; 
            }
L430:       ;   
        }

L440:
        jngp0 = jngp0 + 1;
        ngp[jngp0] = ngp0;
        jngp1 = jngp1 - 1;
        ngp[jngp1] = ngp1;
/* L450 */   ;          
    }
    return;
}           

double ak_side1(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3)
{
    (void)ctx;        /* unused: the signature is shared */
    double r;
    r = (u1 - u3) * (v2 - v3) - (v1 - v3) * (u2 - u3);
    return(r);
}

double ak_spdt1(TDAContext *ctx, double u1,double v1,double u2,double v2,double u3,double v3)
{
    (void)ctx;        /* unused: the signature is shared */
    double r;

    r = (u1 - u2) * (u3 - u2) + (v1 - v2) * (v3 - v2);
    return(r);
}

/* -###-------------------------------------------------------------------- */
/*  ak_idcldp(ndp,xd,yd,ncp,ipc)                                            */
/*                                                                          */
/*  SUBROUTINE  IDCLDP(ndp,XD,YD,ncp,ipc)                                   */
/*                                                                          */
/*  THIS SUBROUTINE SELECTS SEVERAL DATA POINTS THAT ARE CLOSEST            */
/*  TO EACH OF THE DATA POINT.                                              */
/*                                                                          */
/*  THE INPUT PARAMETERS ARE                                                */
/*    ndp = NUMBER OF DATA POINTS,                                          */
/*    xd,yd = ARRAYS OF DIMENSION ndp CONTAINING THE X AND Y                */
/*          COORDINATES OF THE DATA POINTS,                                 */
/*    ncp = NUMBER OF DATA POINTS CLOSEST TO EACH DATA                      */
/*          POINTS.                                                         */
/*                                                                          */
/*  THE OUTPUT PARAMETER IS                                                 */
/*    ipc = INTEGER ARRAY OF DIMENSION ncp*ndp, WHERE THE                   */
/*          POINT NUMBERS OF ncp DATA POINTS CLOSEST TO                     */
/*          EACH OF THE ndp DATA POINTS ARE TO BE STORED.                   */
/*                                                                          */
/*  THIS SUBROUTINE ARBITRARILY SETS A RESTRICTION THAT ncp MUST            */
/*  NOT EXCEED 25.                                                          */
/*                                                                          */
/*  Return 0 if successful, -2 if all points collinear, -1 if error in      */
/*  input parameters.                                                       */

int ak_idcldp(TDAContext *ctx, int ndp,double *xd,double *yd,int ncp,int *ipc)
{
    int ndp0,ncp0,ncpmx,ip1,ip2,ip3,j1,j2,j3,j4,jmx = 0,ip2mn,nclpt,ip3mn;
    double x1,y1,dsqmx,dsqi,dx12,dy12,dx13,dy13,dsqmn;
    double dsq0[26];
    int ipc0[26];

    ncpmx = 25;
    ndp0 = ndp;
    ncp0 = ncp;
    if (ndp0 < 2)
        return(-1);
                 
    if (ncp0 < 1 || ncp0 > ncpmx || ncp0 >= ndp0)
        return(-1);                                                 
                  
    for (ip1 = 1; ip1 <= ndp0; ++ip1) {

        /* SELECTS ncp POINTS. */

        x1 = xd[ip1];
        y1 = yd[ip1];
        j1 = 0;
        dsqmx = 0.0;

        for (ip2 = 1; ip2 <= ndp0; ++ip2) {
            if (ip2 == ip1)
                continue;                 
            dsqi = ak_dsqf(ctx, x1,y1,xd[ip2],yd[ip2]);
            j1++;         
            dsq0[j1] = dsqi;
            ipc0[j1] = ip2;
            if (dsqi <= dsqmx)
                goto L21;
            dsqmx = dsqi;
            jmx = j1;
L21:
            if (j1 >= ncp0)
                goto L23;
        }
L23:
        ip2mn = ip2 + 1;
        if (ip2mn > ndp0)
            goto L30;

        for (ip2 = ip2mn; ip2 <= ndp0; ++ip2) {
            if (ip2 == ip1)
                continue;

            dsqi = ak_dsqf(ctx, x1,y1,xd[ip2],yd[ip2]);
            if (dsqi >= dsqmx) 
                continue;

            dsq0[jmx] = dsqi;
            ipc0[jmx] = ip2;
            dsqmx = 0.0;

            for (j1 = 1; j1 <= ncp0; ++j1) {
                if (dsq0[j1] <= dsqmx)                  
                    continue;
                dsqmx = dsq0[j1];
                jmx = j1;
            }
        }

        /* CHECKS IF ALL THE ncp+1 POINTS ARE COLLINEAR. */

L30:
        ip2 = ipc0[1];
        dx12 = xd[ip2] - x1;
        dy12 = yd[ip2] - y1;

        for (j3 = 2; j3 <= ncp0; ++j3) {
            ip3 = ipc0[j3];
            dx13 = xd[ip3] - x1;
            dy13 = yd[ip3] - y1;
            if ((dy13 * dx12 - dx13 * dy12) != 0.0)
                goto L50;
        }

        /* SEARCHES FOR THE CLOSEST NONCOLLINEAR POINT. */

/* L40 */
        nclpt = 0;

        for (ip3 = 1; ip3 <= ndp0; ++ip3) {

            if (ip3 == ip1)
                continue;

            for (j4 = 1; j4 <= ncp0; ++j4) {
                if (ip3 == ipc0[j4])
                    goto L43;
            }
            dx13 = xd[ip3] - x1;
            dy13 = yd[ip3] - y1;
            if ((dy13 * dx12 - dx13 * dy12) == 0.0)
                goto L43;

            dsqi = ak_dsqf(ctx, x1,y1,xd[ip3],yd[ip3]);
            if (nclpt == 0)
                goto L42;
            if (dsqi >= dsqmn)
                goto L43;
L42:
            nclpt = 1;
            dsqmn = dsqi;
            ip3mn = ip3;
L43:        ; 
        }
        if (nclpt == 0)
            return(-1);   

        dsqmx = dsqmn;
        ipc0[jmx] = ip3mn;
   
        /* REPLACES THE LOCAL ARRAY FOR THE OUTPUT ARRAY. */
L50:
        j1 = (ip1 - 1) * ncp0;
        for (j2 = 1; j2 <= ncp0; ++j2) {
            j1++;        
            ipc[j1] = ipc0[j2];
        }
    }
    return(0);
}


