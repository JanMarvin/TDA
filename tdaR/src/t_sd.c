/****************************************************************************/
/*  t_sd                                                                    */
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
#include "t_gm.h"
#include "tda_context.h"

/*  functions in t_sd.c */

int sd_getdata(TDAContext *ctx, int i,int order,int cflag,int opt);

int sdplot(TDAContext *ctx); 
void sdplot_draw(TDAContext *ctx, int typ,int n,double *x,double *y,int nc);
int sdplot31(TDAContext *ctx);
int sdplot32(TDAContext *ctx);
int sdplot33(TDAContext *ctx);
int sdplot3_prism(TDAContext *ctx, double zval,int iz,int *npol,int *nsn,int *nsf);
int sdpdata(TDAContext *ctx);
int sdinf(TDAContext *ctx);
int sdinf_prn(TDAContext *ctx, int typ,int opt);
int sdencl(TDAContext *ctx);
int sdsel(TDAContext *ctx);
int sdsel_rec(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax, int *n1,int *n2,int *n3);
int sdsel_inside(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax,double x, double y);
int sdsel_prn(TDAContext *ctx, int typ,int n,double *x,double *y,int id,int id1,int sdid);        
int sdsel_prn_poly(TDAContext *ctx, int n,double *x,double *y,char *d,int *nf,double xmin, double ymin,double xmax,double ymax,int *id,int sdid,int *n3);

/* ------------------------------------------------------------------------ */
/*  sd_getdata(i,order,cflag,opt)                                           */
/*                                                                          */
/*  Get data for object i (data matrix record i) and put the data into      */
/*  standard arrays SDVarX, SDVarY. If order != 0 reverse order.            */
/*  If cflag != 0 create a closed polygon. If opt != 0 error message.       */
/*                                                                          */
/*  Return n = number of points in object, or -1 if error.                  */

int sd_getdata(TDAContext *ctx, int i,int order,int cflag,int opt)
{
    register int j,k;
    int fptr,n;
    char *p,buf[121];
    double tmp;

    if (ctx->SDVarFDef == 0) {
        printf1(ctx, "Fatal error in sd_getdata().\n");
        return(-1);
    }
    fptr = (int)get_data(ctx, ctx->SDVarSDPtr,i);
    n    = (int)get_data(ctx, ctx->SDVarSDN,i);

    if (fseek(ctx->SDVarFd,(long)fptr,SEEK_SET)) {
        n = -1;
        goto SDGDFin;
    }
    for (j = 0; j < n; ++j) {
        if (fgets(buf,120,ctx->SDVarFd) == NULL) {
            n = -1;
            goto SDGDFin;
        }                    
        p = skip_b(ctx, buf);
        if (sscanf(p,"%lg",&tmp) != 1) {
            n = -1;
            goto SDGDFin;
        }
        ctx->SDVarX[j] = tmp;

        p = skip_dbl(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%lg",&tmp) != 1) {
            n = -1;
            goto SDGDFin;
        }
        ctx->SDVarY[j] = tmp;
    }
    if (cflag) {
        if (ctx->SDVarX[n - 1] != ctx->SDVarX[0] || ctx->SDVarY[n - 1] != ctx->SDVarY[0]) {
            ctx->SDVarX[n] = ctx->SDVarX[0];
            ctx->SDVarY[n] = ctx->SDVarY[0];
            n++;
        }
    }
    if (order != 0) {
        j = 0; k = n - 1;
        while (j < k) {
            tmp = ctx->SDVarX[j];
            ctx->SDVarX[j] = ctx->SDVarX[k];
            ctx->SDVarX[k] = tmp;
            tmp = ctx->SDVarY[j];
            ctx->SDVarY[j] = ctx->SDVarY[k];
            ctx->SDVarY[k] = tmp;
            j++;
            k--;
        }
    }

SDGDFin:
    if (n < 1 && opt)
        printf1(ctx, "Error: cannot read data for spatial object %d\n",i + 1);
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  sdplot      Plot spatial data in two-dimensional coordinates.           */
/*                                                                          */  
/*              sdplot(                                                     */
/*                  opt=...,    option, def. 1                              */
/*                              1 = plot spatial objects                    */
/*                              2 = convex hull for all points              */
/*                              3 = convex hull for each polygon            */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def. 0.2 mm                     */
/*                  s=...,      number of marker symbol, def. not used      */
/*                  fs=...,     size of marker symbol, def. 2 mm            */
/*                  gs=...,     greyscale value, def. 1 (white)             */
/*                  nc=...,     if nc=1 then no clipping, def. 0            */
/*              );                                                          */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command and a valid PostScript 2d coordinate system.             */
/*                                                                          */
/*  Option 1.                                                               */
/*  The operation depends on the type of the spatial objects.               */
/*                                                                          */
/*  Type 1: Points.     This requires the specification of a marker symbol  */
/*                      Otherwise nothing is done.                          */
/*                                                                          */
/*  Type 2: Lines.      If the line width is positive, a line with the      */
/*                      specified line type and width is drawn. Independent */
/*                      of the line width, if a marker symbol is specified, */
/*                      it is plotted for each point in the line.           */
/*                                                                          */
/*  Type 3: Polygons.   Basically the same procedure as for type 2, but     */
/*                      an implicit closepath is added at the end of the    */  
/*                      data points. Also, if a valid grey scale value is   */
/*                      specified, the polygon is filled accordingly.       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdplot(TDAContext *ctx)
{
    register int i,j,k,ik;
    int err,n,typ,first,nn,np;
    double x,y,d;

    err = -1;
    printf1(ctx, "Plot spatial data (2d). Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (check_pcmd(ctx, 0,2))  /* check for PostScript 2d coordinate system */ 
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,1,0))       /* get parameters */
        goto PLSDFin;

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 1;

    if (ctx->PMS < 1 || ctx->PMS > 17)
        ctx->PMS = 0;           

    fprintf(ctx->PSFd,"\n%%#%d: sdplot\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);

    nn  = 0;    /* number of objects plotted */

    if (ctx->PMOPT == 1) {

        for (i = 0; i < ctx->NOC; ++i) {

            if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            if (typ == 1 && ctx->PMS == 0)  
                continue;
               
            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto PLSDFin;
             
            if (ctx->PMLW > 0.0 || ctx->PMGSFlg) {
                if (typ != 1)
                    sdplot_draw(ctx, typ,n,ctx->SDVarX,ctx->SDVarY,ctx->PMNC);
            }
            if (ctx->PMS) {  /* plot symbols */

                fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

                for (j = 0; j < n; ++j) {
                    x = ctx->SDVarX[j];
                    y = ctx->SDVarY[j];
                    x = ps_2dx(ctx, x);
                    y = ps_2dy(ctx, y);

                    d = ctx->PtMM * ctx->PMFS / 2.0;
                    ps_sym(ctx, ctx->PMS,x,y,d);

                    if (ctx->PMNC) {
                        d *= 2;
                        upd_bbox(ctx, 1,x + d,y + d);
                        upd_bbox(ctx, 1,x + d,y - d);
                        upd_bbox(ctx, 1,x - d,y + d);
                        upd_bbox(ctx, 1,x - d,y - d);
                    }
                }
            }
            nn++;
        }
    }
    else if (ctx->PMOPT == 2) {              /* convex hull for all points */
        np = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            np += (int)get_data(ctx, ctx->SDVarSDN,i);
        }
        if (np == 0)
            goto PLSDFin;
          
        if (alloc_acxf(ctx, np + 2))
            return(-1);   
        if (alloc_acyf(ctx, np + 2))
            return(-1);   
        if (alloc_acn(ctx, np + 2))
            return(-1);   
        if (alloc_acj(ctx, np + 2))
            return(-1);   
        if (alloc_aci(ctx, np + 2))
            return(-1);   

        k = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto PLSDFin;

            for (j = 0; j < n; ++j) {
                k++;
                ctx->AcXF[k] = (float)(ctx->SDVarX[j]);
                ctx->AcYF[k] = (float)(ctx->SDVarY[j]);
                ctx->AcN[k] = k;
            }
        }
        if (k != np) {
            printf1(ctx, "sdplot: k=%d np=%d\n",k,np);
            exit(0);
        }
        n = g_chull(ctx, ctx->AcXF,ctx->AcYF,np,ctx->AcN,ctx->AcJ,ctx->AcI);
        if (n < 0) {
            printf1(ctx, "Error: insufficient memory for convex hull calculation.\n");
            goto PLSDFin;
        }
        ik = ctx->AcI[1];
        first = 0;
        for (i = 1; i <= n; ++i) {
            j = ctx->AcJ[ik];
            x = (double)ctx->AcXF[j];
            y = (double)ctx->AcYF[j];
            ps_2dplot(ctx, x,y,first); 
            if (ctx->PMNC)
                upd_bbox(ctx, 0,x,y);
            first = 1;
            ik = ctx->AcI[ik];
        }
        j = ctx->AcJ[ctx->AcI[1]];
        x = (double)ctx->AcXF[j];
        y = (double)ctx->AcYF[j];
        ps_2dplot(ctx, x,y,1);     /* close path */
     
        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            fprintf(ctx->PSFd,"gsave\n%4.2f setgray\n",ctx->PMGS);    
            fprintf(ctx->PSFd,"fill\ngrestore\n");
        }
        if (ctx->PMLW > 0.0)
            fprintf(ctx->PSFd,"stroke\n");
        nn = ctx->NOC;
    }
    else if (ctx->PMOPT == 3) {              /* convex hulls for polygons */
                    
        for (i = 0; i < ctx->NOC; ++i) {
      
            if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 3)  
                continue;

            if ((np = sd_getdata(ctx, i,0,0,1)) < 1)
                goto PLSDFin;
          
            if (alloc_acxf(ctx, np + 2))
                return(-1);   
            if (alloc_acyf(ctx, np + 2))
                return(-1);   
            if (alloc_acn(ctx, np + 2))
                return(-1);   
            if (alloc_acj(ctx, np + 2))
                return(-1);   
            if (alloc_aci(ctx, np + 2))
                return(-1);   
      
            k = 0;
            for (j = 0; j < np; ++j) {
                k++;
                ctx->AcXF[k] = (float)(ctx->SDVarX[j]);
                ctx->AcYF[k] = (float)(ctx->SDVarY[j]);
                ctx->AcN[k] = k;
            }
            n = g_chull(ctx, ctx->AcXF,ctx->AcYF,np,ctx->AcN,ctx->AcJ,ctx->AcI);
            if (n < 0) {
                printf1(ctx, "Error: insufficient memory for convex hull calculation.\n");
                goto PLSDFin;
            }
            ik = ctx->AcI[1];
            first = 0;
            for (k = 1; k <= n; ++k) {
                j = ctx->AcJ[ik];
                x = (double)ctx->AcXF[j];
                y = (double)ctx->AcYF[j];
                ps_2dplot(ctx, x,y,first); 
                if (ctx->PMNC)
                    upd_bbox(ctx, 0,x,y);
                first = 1;
                ik = ctx->AcI[ik];
            }
            j = ctx->AcJ[ctx->AcI[1]];
            x = (double)ctx->AcXF[j];
            y = (double)ctx->AcYF[j];
            ps_2dplot(ctx, x,y,1);     /* close path */
     
            if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
                fprintf(ctx->PSFd,"gsave\n%4.2f setgray\n",ctx->PMGS);    
                fprintf(ctx->PSFd,"fill\ngrestore\n");
            }
            if (ctx->PMLW > 0.0)
                fprintf(ctx->PSFd,"stroke\n");
            nn++;
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
    printf1(ctx, "PostScript output for %d object(s) written to: %s\n",nn,ctx->PSFName);
    err = 0;

PLSDFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdplot_draw(typ,n,x,y,nc)                                               */
/*                                                                          */
/*  x[i], y[i] (i=0,n-1) contain the points of a line (typ 2) or            */
/*  polygon (typ 3). This function draws the line, or polygon, assuming     */
/*  2d cartesian coordinates. For type 3 there will be an implicit          */
/*  closepath. If nc != 0 update bounding box.                              */
    
void sdplot_draw(TDAContext *ctx, int typ,int n,double *x,double *y,int nc)
{
    register int i,first;

    first = 0;
    for (i = 0; i < n; ++i) {
        ps_2dplot(ctx, x[i],y[i],first);
        first = 1;
        if (nc)
            upd_bbox(ctx, 0,x[i],y[i]);
    }
    if (typ == 3) {
        fprintf(ctx->PSFd,"closepath\n");    
        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            fprintf(ctx->PSFd,"gsave\n%4.2f setgray\n",ctx->PMGS);    
            fprintf(ctx->PSFd,"fill\ngrestore\n");
        }
    }
    fprintf(ctx->PSFd,"stroke\n");
}

/* -##--------------------------------------------------------------------- */
/*  sdplot31    Plot type 1 objects in three-dimensional coordinates.       */
/*                                                                          */  
/*              sdplot31(                                                   */
/*                  zval = ..., z value of x-y plane, def. 0.0              */
/*                  zvar = ..., optional z coordinates                      */
/*                  s=...,      number of marker symbol, def. 0             */
/*                  fs=...,     size of marker symbol, def. 2 mm            */
/*                  nc=...,     if nc=1 then no clipping, def. 0            */
/*                                                                          */
/*              );                                                          */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command and a valid PostScript 3d coordinate system.             */
/*                                                                          */
/*  The command considers all type 1 objects (points) in the spatial data   */
/*  structure. The output depends on the zvar and s parameters. To explain  */
/*  the options, assume that a point has coordinates (x,y) and zval         */
/*  specifies a z coordinate z.                                             */
/*                                                                          */
/*  a)  If neither the zvar nor the s parameter is used, nothing is done.   */
/*                                                                          */
/*  b)  If the s parameter specifies a valid marker symbol, but the zvar    */
/*      parameter is not used, the symbol is plotted at (x,y,z).            */
/*                                                                          */
/*  c)  If the s parameter is not used but the zvar parameter specifies     */
/*      an attribute variable, the value of this variable is taken as a     */
/*      z coordinate for the currently selected point, say z1, and a        */
/*      straight line is drawn from (x,y,z) to (x,y,z1).                    */
/*                                                                          */
/*  d)  If both the s and zvar parameters are used, the command first draws */
/*      a straight line from (x,y,z) to (x,y,z1) and then plots the symbol  */
/*      at position (x,y,z1).                                               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdplot31(TDAContext *ctx)
{
    int err,i,n,nn,np,mode;
    double x,y,z,z1,d;

    err = -1;
    if (check_pcmd(ctx, 0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1(ctx, "Plot type 1 objects (3d). Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 8,1,0))       /* get parameters */
        goto PLS31Fin;

    if (ctx->PMS < 1 || ctx->PMS > 17)
        ctx->PMS = 0;           

    mode = 0;
    if (ctx->PMZVar >= 0)  
        mode = 1;

    if (mode == 0 && ctx->PMS == 0) {
        printf1(ctx, "Nothing done.\n");
        err = 0;
        goto PLS31Fin;
    }
    nn = 0;    /* number of objects plotted */

    fprintf(ctx->PSFd,"\n%%#%d: sdplot31\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMS != 0)
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

    np = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 1)
            continue;

        np++;
        if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
            goto PLS31Fin;
                      
        if (mode == 1) {
            z1 = get_data(ctx, ctx->PMZVar,i);
            ps_3dplot(ctx, ctx->SDVarX[0],ctx->SDVarY[0],ctx->PMZVal,0,0,ctx->PMNC);
            ps_3dplot(ctx, ctx->SDVarX[0],ctx->SDVarY[0],z1,1,0,ctx->PMNC);
            fprintf(ctx->PSFd,"stroke\n");
        }
        if (ctx->PMS) {
            if (mode == 1)
                z = z1;
            else
                z = ctx->PMZVal;

            ps_3dprj(ctx, ctx->SDVarX[0],ctx->SDVarY[0],z,&x,&y);
            x = ps_2dx(ctx, x);
            y = ps_2dy(ctx, y);

            d = ctx->PtMM * ctx->PMFS / 2.0;
            ps_sym(ctx, ctx->PMS,x,y,d);

            if (ctx->PMNC) {
                d *= 2;
                upd_bbox(ctx, 1,x + d,y + d);
                upd_bbox(ctx, 1,x + d,y - d);
                upd_bbox(ctx, 1,x - d,y + d);
                upd_bbox(ctx, 1,x - d,y - d);
            }
        }
        nn++;
    }
    fprintf(ctx->PSFd,"grestore\n");
    printf1(ctx, "Number of type 1 objects: %d\n",np);
    printf1(ctx, "PostScript output for %d object(s) written to: %s\n",nn,ctx->PSFName);
    err = 0;

PLS31Fin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdplot32    Plot type 2 objects in three-dimensional coordinates.       */
/*                                                                          */  
/*              sdplot32(                                                   */
/*                  zval = ..., z value of x-y plane, def. 0.0              */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def. 0.2 mm                     */
/*                  s=...,      number of marker symbol, def. not used      */
/*                  fs=...,     size of marker symbol, def. 2 mm            */
/*                  nc=...,     if nc=1 then no clipping, def. 0            */
/*                                                                          */
/*              );                                                          */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command and a valid PostScript 3d coordinate system.             */
/*                                                                          */
/*  The command recognizes all type 2 objects (lines) in the spatial data   */
/*  structure and, if line width > 0, plots the line in the x-y plane       */
/*  defined by the z coordinate specified with the zval parameter.          */
/*                                                                          */  
/*  In addition, if the s parameter specfies a valid marker symbol, this    */
/*  symbol is plotted at all end points of the segments of the line.        */  
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdplot32(TDAContext *ctx)
{
    register int i,j;
    int err,n,first,nn,np,mode;
    double x,y,z,d;

    err = -1;
    if (check_pcmd(ctx, 0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1(ctx, "Plot type 2 objects (3d). Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 8,1,0))       /* get parameters */
        goto PLS32Fin;

    if (ctx->PMS < 1 || ctx->PMS > 171)
        ctx->PMS = 0;           

    mode = 0;
    if (ctx->PMZVar >= 0)  
        mode = 1;

    nn = 0;    /* number of objects plotted */

    if (ctx->PMS == 0 && ctx->PMLW <= 0.0) {
        printf1(ctx, "Nothing done.\n");
        err = 0;
        goto PLS32Fin;
    }

    fprintf(ctx->PSFd,"\n%%#%d: sdplot32\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);

    if (ctx->PMS != 0)
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

    np = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 2)
            continue;

        np++;
        if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
            goto PLS32Fin;

        if (ctx->PMLW > 0.0) {
            first = 0;
            for (j = 0; j < n; ++j) {
                ps_3dplot(ctx, ctx->SDVarX[j],ctx->SDVarY[j],ctx->PMZVal,first,0,ctx->PMNC);
                first = 1;
            }
            fprintf(ctx->PSFd,"stroke\n");
        }
        if (ctx->PMS) {
            if (mode == 1)
                z = get_data(ctx, ctx->PMZVar,i);
            else
                z = ctx->PMZVal;

            for (j = 0; j < n; ++j) {
                ps_3dprj(ctx, ctx->SDVarX[j],ctx->SDVarY[j],z,&x,&y);
                x = ps_2dx(ctx, x);
                y = ps_2dy(ctx, y);

                d = ctx->PtMM * ctx->PMFS / 2.0;
                ps_sym(ctx, ctx->PMS,x,y,d);

                if (ctx->PMNC) {
                    d *= 2;
                    upd_bbox(ctx, 1,x + d,y + d);
                    upd_bbox(ctx, 1,x + d,y - d);
                    upd_bbox(ctx, 1,x - d,y + d);
                    upd_bbox(ctx, 1,x - d,y - d);
                }
            }
        }
        nn++;
    }
    fprintf(ctx->PSFd,"grestore\n");
    printf1(ctx, "Number of type 2 objects: %d\n",np);
    printf1(ctx, "PostScript output for %d object(s) written to: %s\n",nn,ctx->PSFName);
    err = 0;

PLS32Fin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdplot33    Plot type 3 objects in three-dimensional coordinates.       */
/*                                                                          */  
/*              sdplot33(                                                   */
/*                  zval = ..., z value of x-y plane, def. 0.0              */
/*                  zvar = ..., optional z coordinates                      */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def. 0.2 mm                     */
/*                  s=...,      number of marker symbol, def. not used      */
/*                  fs=...,     size of marker symbol, def. 2 mm            */
/*                  gs=...,     greyscale value, def. 1 (white)             */
/*                  nc=...,     if nc=1 then no clipping, def. 0            */
/*                                                                          */
/*              );                                                          */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command and a valid PostScript 3d coordinate system.             */
/*                                                                          */
/*  The command recognizes all type 3 objects (polygons) in the spatial     */
/*  data structure. The operation depends on the zvar parameter.            */
/*                                                                          */
/*  a) If the zvar parameter is not used, the polygons are plotted in the   */
/*     x-y plane defined by the z-coordinate defined by zval. Given a       */
/*     positive line width, the polygon is plotted with the specified line  */
/*     type and line width. In this case, the gs parameter can be used to   */
/*     specify a grey scale value for filling the area of the polygons.     */
/*                                                                          */
/*     In addition, if the s parameter specfies a valid marker symbol, this */
/*     symbol is plotted at all end points of the segments of the polygon.  */  
/*                                                                          */
/*  b) If the zvar parameter is used to specify an attribute variable, the  */
/*     values of this variable are taken to define the heights of prisms    */  
/*     having the polygons as their bases. The command then tries to plot   */
/*     the prisms in a way that non-visible parts are not shown.            */
/*                                                                          */
/*     Note 1: The algorithm is based on the possibility to overdraw filled */
/*     regions with the PostScript fill command. In order to work properly, */
/*     it is required that the view parameter in the psetup3 command        */  
/*     specifies a positive latitude.                                       */
/*                                                                          */
/*     Note 2: Also, in order to work properly, it is required that all     */
/*     prisms have a positive height. Otherwise, they will not be drawn.    */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdplot33(TDAContext *ctx)
{
    register int i,j;
    int err,n,first,nn,mode,ns,npol,nsn,nsf;
    double x,y,d;

    err = -1;
    if (check_pcmd(ctx, 0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1(ctx, "Plot type 3 objects (3d). Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 8,1,0))       /* get parameters */
        goto PLS33Fin;

    if (ctx->PMS < 1 || ctx->PMS > 17)
        ctx->PMS = 0;           

    if (ctx->PMZVar >= 0) {
        mode = 1;
        if (ctx->PMGSFlg == 0 || ctx->PMGS < 0.0 || ctx->PMGS > 1.0) {
            ctx->PMGSFlg = 1;
            ctx->PMGS = 1.0;
        }
        ctx->PMS = 0;
    }
    else {
        mode = 0;
        if (ctx->PMGS < 0.0 || ctx->PMGS > 1.0)
            ctx->PMGSFlg = 0;
    }
    nn = 0;    /* number of objects plotted */
    nsn = nsf = npol = 0;

    if (mode == 0 && ctx->PMS == 0 && ctx->PMLW <= 0.0 && ctx->PMGSFlg == 0) {
        printf1(ctx, "Nothing done.\n");
        err = 0;
        goto PLS33Fin;
    }

    fprintf(ctx->PSFd,"\n%%#%d: sdplot33\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);

    if (mode == 0) {
        if (ctx->PMS != 0)
            fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

        for (i = 0; i < ctx->NOC; ++i) {
            if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 3)
                continue;

            npol++;
            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto PLS33Fin;

            ns = 0;
            if (ctx->PMLW > 0.0 || ctx->PMGSFlg != 0) {
                first = 0;
                for (j = 0; j < n; ++j) {
                    ps_3dplot(ctx, ctx->SDVarX[j],ctx->SDVarY[j],ctx->PMZVal,first,0,ctx->PMNC);
                    first = 1;
                }
                fprintf(ctx->PSFd,"closepath\n");    
                if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
                    fprintf(ctx->PSFd,"gsave\n%4.2f setgray\n",ctx->PMGS);    
                    fprintf(ctx->PSFd,"fill\ngrestore\n");
                }
                fprintf(ctx->PSFd,"stroke\n");
                ns = 1;
            }
            if (ctx->PMS) {
                for (j = 0; j < n; ++j) {
                    ps_3dprj(ctx, ctx->SDVarX[j],ctx->SDVarY[j],ctx->PMZVal,&x,&y);
                    x = ps_2dx(ctx, x);
                    y = ps_2dy(ctx, y);

                    d = ctx->PtMM * ctx->PMFS / 2.0;
                    ps_sym(ctx, ctx->PMS,x,y,d);

                    if (ctx->PMNC) {
                        d *= 2;
                        upd_bbox(ctx, 1,x + d,y + d);
                        upd_bbox(ctx, 1,x + d,y - d);
                        upd_bbox(ctx, 1,x - d,y + d);
                        upd_bbox(ctx, 1,x - d,y - d);
                    }
                    ns = 1;
                }
            }
            if (ns)
                nn++;
        }
        fprintf(ctx->PSFd,"grestore\n");
    }
    else {
        nn = sdplot3_prism(ctx, ctx->PMZVal,ctx->PMZVar,&npol,&nsn,&nsf);
        if (nn < 0)
            goto PLS33Fin;

        fprintf(ctx->PSFd,"grestore\n");
    }
    printf1(ctx, "Number of type 3 objects: %d\n",npol);
    if (nn > 0)
        printf1(ctx, "PostScript output for %d object(s) written to: %s\n",nn,ctx->PSFName);
    if (nsn > 0)
        printf1(ctx, "Omitted %d prism(s) with non-positive height.\n",nsn);
    if (nsf > 0)
        printf1(ctx, "Omitted %d faces with unknown priority.\n",nsf);
    err = 0;

PLS33Fin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdplot3_prism(zval,iz,npol,nsn,nsf)                                     */
/*                                                                          */
/*  Plot the polygons in the current spatial data as prisms. Return         */
/*  number of objects plotted, or -1 if an error occurred. Return the       */
/*  number of polygons in npol, the number of prisms with non-positive      */ 
/*  height in nsn, and the number of omitted faces in nsf.                  */

int sdplot3_prism(TDAContext *ctx, double zval,int iz,int *npol,int *nsn,int *nsf)
{
    register int i,j,k,l;
    int npp,nn,np,nx,nx0,ng,n0,npl;
    double a,zval1,x,y,z,xia,xib,xja,xjb,yia,yib,yja,yjb,xa,ya,xb,yb,yi,yj;
    double tmpi,tmpj,tol;

    tol = 100.0 * ctx->EPSI1;

    *npol = *nsn = *nsf;

    /* always with fill */

    if (ctx->PMGSFlg == 0 || ctx->PMGS < 0.0 || ctx->PMGS > 1.0) {
        ctx->PMGSFlg = 1;
        ctx->PMGS = 1.0;
    }

    /* AcC[i] = 1 for positive orientation (counterclockwise), -1 for
       negative orientation, 0 if not a polygon. */

    if (alloc_acc(ctx, ctx->NOC + 1))
        return(-1);   

    npp = 0;                /* number of polygons */
    nn  = 0;                /* number of points */

    for (i = 0; i < ctx->NOC; ++i) {

        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 3)  
            continue;

        *npol += 1;
        zval1 = get_data(ctx, iz,i);          
        if (zval1 <= zval) {
            *nsn += 1;
            continue;
        }

        if ((np = sd_getdata(ctx, i,0,1,1)) < 1)
            return(-1);       

        a = g_pol_area(ctx, np - 1,ctx->SDVarX,ctx->SDVarY);
        npp++;
        nn += np;
        if (a >= 0.0)
            ctx->AcC[i] = 1;
        else
            ctx->AcC[i] = -1;
    }
    if (npp == 0)
        return(0);

    /* create a list containing all segments. AcD has the following meaning:
       AcD = 0 do not draw vertical lines
             1 only draw left vertical line
             2 draw both vertical lines
             3 only draw right vertical line
    */

    if (alloc_acxf(ctx, nn + 1))
        return(-1);   
    if (alloc_acuf(ctx, nn + 1))
        return(-1);   
    if (alloc_acyf(ctx, nn + 1))
        return(-1);   
    if (alloc_acvf(ctx, nn + 1))
        return(-1);   
    if (alloc_aczf(ctx, nn + 1))
        return(-1);   
    if (alloc_acd(ctx, nn + 1))
        return(-1);   

    nx = 0;                                
    for (i = 0; i < ctx->NOC; ++i) {
        if (ctx->AcC[i] == 0)
            continue;

        l = 0;
        if (ctx->AcC[i] == -1)
            l = 1;

        if ((np = sd_getdata(ctx, i,l,1,1)) < 1)
            return(-1);       
                 
        zval1 = get_data(ctx, iz,i);          
        if (zval1 <= zval)  
            continue;

        nx0 = nx;
        for (j = 1; j < np; ++j) {
            ps_3dprj1(ctx, ctx->SDVarX[j - 1],ctx->SDVarY[j - 1],zval,&xia,&y,&z);
            ctx->AcXF[nx] = (float)xia;
            ctx->AcYF[nx] = (float)y;
            ps_3dprj1(ctx, ctx->SDVarX[j],ctx->SDVarY[j],zval,&xa,&ya,&z);
            ctx->AcUF[nx] = (float)xa;
            ctx->AcVF[nx] = (float)ya;
            ps_3dprj1(ctx, ctx->SDVarX[j],ctx->SDVarY[j],zval1,&x,&yb,&z);
            ctx->AcZF[nx] = (float)(yb - ya);
      
            if (xia <= xa) {
                ctx->AcD[nx] = 1;
                if (j > 1) {
                    if (ctx->AcD[nx - 1] != 0) {
                        if (ctx->AcD[nx - 1] == 1)
                            ctx->AcD[nx - 1] = 2;
                    }
                    else {
                        xia = (double)ctx->AcXF[nx - 1];
                        yia = (double)ctx->AcYF[nx - 1];
                        xib = (double)ctx->AcUF[nx - 1];
                        yib = (double)ctx->AcVF[nx - 1];
       
                        if (g_left(ctx, xia,yia,xib,yib,xa,ya) == 0)  
                            ctx->AcD[nx] = 3;
                    }
                }
            }
            else {   
                ctx->AcYF[nx] = (float)((double)(ctx->AcYF[nx]) - (tol));
                ctx->AcVF[nx] = (float)((double)(ctx->AcVF[nx]) - (tol));

                if (j > 1 && ctx->AcD[nx - 1] == 1) {
                    xia = (double)ctx->AcXF[nx - 1];
                    yia = (double)ctx->AcYF[nx - 1];
                    xib = (double)ctx->AcUF[nx - 1];
                    yib = (double)ctx->AcVF[nx - 1];
                    if (g_left(ctx, xia,yia,xib,yib,xa,ya))   
                        ctx->AcD[nx - 1] = 2;
                }
            }             
            nx++;
        }
        if (ctx->AcD[nx - 1] == 1 && ctx->AcD[nx0] == 0)  
            ctx->AcD[nx - 1] = 2;

    }
        
    /* create an edge list for a graph that records the priority   
       relations: i -> j if i must be drawn before j. */
     
    ng = 0;         /* number of nodes */

    for (l = 0; l < 2; ++l) {

        k = 0;
        for (i = 0; i < nx; ++i) {
            xia = (double)ctx->AcXF[i];
            xib = (double)ctx->AcUF[i];
            yia = (double)ctx->AcYF[i];
            yib = (double)ctx->AcVF[i];

            for (j = i + 1; j < nx; ++j) {
                xja = (double)ctx->AcXF[j];
                xjb = (double)ctx->AcUF[j];
                yja = (double)ctx->AcYF[j];
                yjb = (double)ctx->AcVF[j];

                if (xia > xib) {
                    g_xchxy(ctx, &xia,&xib);
                    g_xchxy(ctx, &yia,&yib);
                }
                if (xja > xjb) {
                    g_xchxy(ctx, &xja,&xjb);
                    g_xchxy(ctx, &yja,&yjb);
                }
                if ((xa = xja) < xia)  
                    xa = xia;

                if ((xb = xjb) > xib)  
                    xb = xib;

                if (xa >= xb)    
                    continue;
 
                if ((tmpi = xib - xia) == 0.0)
                    tmpi = ctx->EPSI1;

                if ((tmpj = xjb - xja) == 0.0)
                    tmpj = ctx->EPSI1;

                x = (xa + xb) / 2.0;
                yi = (yib - yia) * (x - xia) / tmpi + yia;
                yj = (yjb - yja) * (x - xja) / tmpj + yja;
                /***
                if (fabs(yi - yj) <= EPSI1) {
                    continue;
                }
                else 
                **/                  
                if (yi <= yj) {
                    if (l) {
                        ctx->AcN[k] = j;
                        ctx->AcM[k] = i;
                        k++;
                    }
                }
                else {
                    if (l) {
                        ctx->AcN[k] = i;
                        ctx->AcM[k] = j;
                        k++;
                    }
                }
                if (l == 0)
                    ng++;
            }
        }
        if (l)
            break;

        if (alloc_acn(ctx, ng + 1))
            return(-1);   
        if (alloc_acm(ctx, ng + 1))
            return(-1);   
    }

    /* sort the edge list according to nodes in AcN and AcM */                
    
    if (ng > 0) {
        if (sorti2(ctx, ng,ctx->AcN,ctx->AcM))
            return(-1);
    }

    /* create pointer to begin of new nodes in AcK */

    if (alloc_ack(ctx, nx + 1))
        return(-1);   

    for (i = 0; i < nx; ++i)
        ctx->AcK[i] = -1;

    n0 = -1;
    for (i = 0; i < ng; ++i) {
        if (ctx->AcN[i] != n0) {
            n0 = ctx->AcN[i];
            ctx->AcK[n0] = i;
        }           
    }

    /* first fill all polygons at the base heigth zval */

    for (i = 0; i < ctx->NOC; ++i) {
        if (ctx->AcC[i] == 0)
            continue;

        l = 0;
        if (ctx->AcC[i] == -1)
            l = 1;

        if ((np = sd_getdata(ctx, i,l,1,1)) < 1)
            return(-1);       
                 
        fprintf(ctx->PSFd,"gsave\n");    

        l = 0;
        for (j = 0; j < np; ++j) {
            ps_3dprj1(ctx, ctx->SDVarX[j],ctx->SDVarY[j],zval,&x,&y,&z);
            ps_2dplot(ctx, x,y,l);
            l = 1;
        }
        fprintf(ctx->PSFd,"%4.2f setgray\n",ctx->PMGS);    
        fprintf(ctx->PSFd,"fill\n");
        fprintf(ctx->PSFd,"stroke\n");
        fprintf(ctx->PSFd,"grestore\n");
    }

    /* plot the faces according to the priorities given by the graph */

    if (alloc_acns(ctx, nx + 1))     /* counts indegrees of nodes */
        return(-1);   

    for (i = 0; i < ng; ++i)  
        ctx->AcNS[ctx->AcM[i]] += 1;

    npl = 0;
    while (1) {
        j = 0;
        for (i = 0; i < nx; ++i) {
            if (ctx->AcNS[i] == 0) {
                xia = (double)ctx->AcXF[i];
                yia = (double)ctx->AcYF[i];
                xib = (double)ctx->AcUF[i];
                yib = (double)ctx->AcVF[i];
                z   = (double)ctx->AcZF[i];

                fprintf(ctx->PSFd,"gsave\n");    
                ps_2dplot(ctx, xia,yia,0);
                ps_2dplot(ctx, xib,yib,1);
                ps_2dplot(ctx, xib,yib + z,1);
                ps_2dplot(ctx, xia,yia + z,1);
                ps_2dplot(ctx, xia,yia,1);
                fprintf(ctx->PSFd,"%4.2f setgray\n",ctx->PMGS);    
                fprintf(ctx->PSFd,"fill\n");
                fprintf(ctx->PSFd,"stroke\n");
                fprintf(ctx->PSFd,"grestore\n");

                ps_2dplot(ctx, xia,yia + z,0);
                ps_2dplot(ctx, xib,yib + z,1);
                fprintf(ctx->PSFd,"stroke\n");

                if (ctx->AcD[i] != 0) {
                    ps_2dplot(ctx, xia,yia,0);
                    ps_2dplot(ctx, xib,yib,1);
                    fprintf(ctx->PSFd,"stroke\n");
                      
                    if (ctx->AcD[i] <= 2) {
                        ps_2dplot(ctx, xia,yia,0);
                        ps_2dplot(ctx, xia,yia + z,1);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    if (ctx->AcD[i] >= 2) {
                        ps_2dplot(ctx, xib,yib,0);
                        ps_2dplot(ctx, xib,yib + z,1);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                }
                npl++;
        
                j++;
                ctx->AcNS[i] = -1;
                k = ctx->AcK[i];
                if (k < 0)
                    continue;

                for (l = k; l < ng; ++l) {      /* adjust indegrees */
                    if (ctx->AcN[l] != i)
                        break;
                    ctx->AcNS[ctx->AcM[l]] -= 1;
                }
            }
        }
        if (j == 0)   
            break;
    }   
    if (npl < nx) {
        *nsf = nx - npl;
        printf1(ctx, "Warning: omitted %d faces.\n",nx - npl);
    }

    /* free memory */
    alloc_acc(ctx, 0);
    alloc_acd(ctx, 0);
    alloc_acxf(ctx, 0);
    alloc_acyf(ctx, 0);
    alloc_aczf(ctx, 0);
    alloc_acuf(ctx, 0);
    alloc_acvf(ctx, 0);
    alloc_acn(ctx, 0);
    alloc_acm(ctx, 0);
    alloc_ack(ctx, 0);
    alloc_acns(ctx, 0);

    return(npp);
}

/* -##--------------------------------------------------------------------- */
/*  sdpdata     Write spatial data to output file.                          */
/*                                                                          */  
/*              sdpdata(                                                    */
/*                  opt=...,    option, def. 1                              */
/*                              1 = points                                  */
/*                              2 = line segments                           */
/*                              3 = line segments (closed polygons)         */
/*                  fmt=...,    print format for coordinates, def. 10.4     */
/*              ) = output_file;                                            */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command.                                                         */
/*                                                                          */  
/*  Option 1 writes one record for each point in all type 1, type 2 and     */
/*  type 3 objects. Each record has three entries:                          */
/*                                                                          */
/*  SDID  x-coordinate  y-coordinate                                        */
/*                                                                          */
/*  Option 2 writes one record for each line segment in all type 2 and      */
/*  type 3 objects. Each record has five entries:                           */
/*                                                                          */
/*  SDID  x1-coordinate  y1-coordinate  x1-coordinate  y2-coordinate        */
/*                                                                          */
/*  Option 3 is identical with option 2 but an additional segment from the  */
/*  last to the first point of each polygon is added (if not already in     */
/*  the data set).                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdpdata(TDAContext *ctx)
{
    register int i,j;
    int err,n,typ,sdid,nrec,cflag;

    err = -1;
    printf1(ctx, "Write spatial data. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 7,1,1))       /* get parameters */
        goto SDPDFin;

    if (ctx->PMOPT > 3)
        ctx->PMOPT = 1;

    cflag = 0;
    if (ctx->PMOPT == 3) {
        cflag = 1;
        ctx->PMOPT = 2;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        if (ctx->PMOPT == 2 && typ == 1)
            continue;

        sdid = (int)get_data(ctx, ctx->SDVarSDID,i);                     

        if (ctx->PMOPT == 1) {
            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto SDPDFin;

            for (j = 0; j < n; ++j) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[ctx->SDVarSDID],(double)sdid);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarX[j]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarY[j]);
                fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                /* SDVarX/Y are double; the file is written at fmt=,
                   default 10.4, so real geographic coordinates lost
                   everything past the fourth decimal on the way back */
                {
                    double erow[3];
                    erow[0] = (double)sdid;
                    erow[1] = ctx->SDVarX[j];
                    erow[2] = ctx->SDVarY[j];
                    tda_export_row(ctx, "sdpdata.table", erow, 3);
                }
#endif
                nrec++;
            }
        }    
        else if (ctx->PMOPT == 2) {
            if ((n = sd_getdata(ctx, i,0,cflag,1)) < 1)
                goto SDPDFin;

            for (j = 1; j < n; ++j) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[ctx->SDVarSDID],(double)sdid);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarX[j - 1]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarY[j - 1]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarX[j]);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->SDVarY[j]);
                fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
                {
                    double erow[5];
                    erow[0] = (double)sdid;
                    erow[1] = ctx->SDVarX[j - 1];
                    erow[2] = ctx->SDVarY[j - 1];
                    erow[3] = ctx->SDVarX[j];
                    erow[4] = ctx->SDVarY[j];
                    tda_export_row(ctx, "sdpdata.table", erow, 5);
                }
#endif
                nrec++;
            }
        }    
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDPDFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdinf   Elementary inforamtion about spatial data.                      */
/*                                                                          */  
/*          sdinf(                                                          */
/*              opt=...,    option, def. 1                                  */  
/*                          1 = basic information                           */
/*                          2 = additional information about orientation    */  
/*                              of polygons                                 */
/*              eps=...,    tolerance for area check, def. 1.e-6            */
/*          );                                                              */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command and provides some basic information about the spatial    */
/*  objects.                                                                */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdinf(TDAContext *ctx)
{
    int err;

    err = -1;
    printf1(ctx, "Information about spatial data. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,0))       /* get parameters */
        goto SDINFFin;

    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;

    if (sdinf_prn(ctx, 1,ctx->PMOPT))
        goto SDINFFin;

    if (sdinf_prn(ctx, 2,ctx->PMOPT))
        goto SDINFFin;

    if (sdinf_prn(ctx, 3,ctx->PMOPT))
        goto SDINFFin;

    err = 0;

SDINFFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdinf_prn(typ,opt)                                                      */

int sdinf_prn(TDAContext *ctx, int typ,int opt)
{
    (void)opt;        /* unused: the signature is shared */
    register int i,j;
    int np,nt,t,n,first,cflag,npp,npn,np0;
    double xmin,xmax,ymin,ymax,pdmin,pdmax,d;

    cflag = 0;
    if (typ == 3)
        cflag = 1;

    pdmin = ctx->DBLMAX;
    pdmax = 0.0;

    first = 1;
    np0 = npp = npn = np = nt = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        t = (int)get_data(ctx, ctx->SDVarSDTyp,i);                       
        if (t != typ)
            continue;
        if ((n = sd_getdata(ctx, i,0,cflag,1)) < 1)
            return(-1);             
        if (typ == 3)
            n--;
        nt++;
        np += n;

        if (first) {
            xmin = xmax = ctx->SDVarX[0];
            ymin = ymax = ctx->SDVarY[0];
            first = 0;
        }
        for (j = 0; j < n; ++j) {
            xmin = dmin(ctx, xmin,ctx->SDVarX[j]);
            xmax = dmax(ctx, xmax,ctx->SDVarX[j]);
            ymin = dmin(ctx, ymin,ctx->SDVarY[j]);
            ymax = dmax(ctx, ymax,ctx->SDVarY[j]);
        }
        if (typ == 2 || typ == 3) {
            for (j = 1; j < n; ++j) {
                d = g_len(ctx, ctx->SDVarX[j],ctx->SDVarY[j],ctx->SDVarX[j-1],ctx->SDVarY[j-1]);
                pdmin = dmin(ctx, pdmin,d);
                pdmax = dmax(ctx, pdmax,d);
            }
        }
        if (ctx->PMOPT == 2 && typ == 3) {
            d = g_pol_area(ctx, n,ctx->SDVarX,ctx->SDVarY);
            if (d >= ctx->PMEPS)
                npp++;
            else if (d <= -ctx->PMEPS)
                npn++;
            else
                np0++;
        }
    }
    printf1(ctx, "\nNumber of ");
    if (typ == 1)
        printf1(ctx, "points:");
    else if (typ == 2)
        printf1(ctx, "lines:");
    else            
        printf1(ctx, "polygons:");
    printf1(ctx, " %d",nt);
    if (nt > 0 && (typ == 2 || typ == 3))
        printf1(ctx, "  (number of points: %d)",np);
    printf1(ctx, "\n");
    if (nt == 0)
        return(0);

    printf1(ctx, "XMin: %18.12f   YMin: %18.12f\n",xmin,ymin);
    printf1(ctx, "XMax: %18.12f   YMax: %18.12f\n",xmax,ymax);

    if (typ == 2 || typ == 3) {
        printf1(ctx, "\nMinimal segment length: %18.12f\n",pdmin);
        printf1(ctx, "Maximal segment length: %18.12f\n",pdmax);
    }
    if (ctx->PMOPT == 2 && typ == 3) {
        printf1(ctx, "\nPolygons with positive orientation: %d\n",npp);
        printf1(ctx, "Polygons with negative orientation: %d\n",npn);
        printf1(ctx, "Polygons with almost zero (less than %g) area: %d\n",ctx->PMEPS,np0);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sdencl      Creating enclosing boundaries.                              */
/*                                                                          */  
/*              sdencl(                                                     */
/*                  opt=...,        type of boundary, def. 1                */
/*                                  1 = rectangle                           */
/*                                  2 = convex hull                         */
/*                  fmt=...,        print format for coordinates, def. 10.4 */
/*              ) = output_file;                                            */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command.                                                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdencl(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,nn,nrec,typ,first,np;
    double xmin = 0.0,xmax = 0.0,ymin = 0.0,ymax = 0.0;                

    err = -1;
    nrec = 0;
    printf1(ctx, "Enclosing boundaries. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto SDENCLFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    printf1(ctx, "Type of boundary: ");
    if (ctx->PMOPT == 2)
        printf1(ctx, "convex hull.\n\n");
    else {
        printf1(ctx, "rectangle.\n\n");
        ctx->PMOPT = 1;
    }
    nn = imax(ctx, ctx->SDVarNT,4);

    if (ctx->PMOPT == 1) {                   /* rectangle */
        first = 1;
        for (i = 0; i < ctx->NOC; ++i) {
            typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i);                       
            if (typ < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto SDENCLFin;         

            for (j = 0; j < n; ++j) {
                if (first) {
                    xmin = xmax = ctx->SDVarX[j];
                    ymin = ymax = ctx->SDVarY[j];
                    first = 0;
                }
                else {
                    xmin = dmin(ctx, xmin,ctx->SDVarX[j]);
                    xmax = dmax(ctx, xmax,ctx->SDVarX[j]);
                    ymin = dmin(ctx, ymin,ctx->SDVarY[j]);
                    ymax = dmax(ctx, ymax,ctx->SDVarY[j]);
                }
            }
        }
        fprintf(ctx->PMFd,"1 3 4\n");
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,xmin);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ymin);
        fprintf(ctx->PMFd,"\n");
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,xmax);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ymin);
        fprintf(ctx->PMFd,"\n");
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,xmax);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ymax);
        fprintf(ctx->PMFd,"\n");
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,xmin);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ymax);
        fprintf(ctx->PMFd,"\n");
        nrec = 5;
    }
    else if (ctx->PMOPT == 2) {              /* convex hull */

        if (alloc_acxf(ctx, nn + 2))
            goto SDENCLFin;
        if (alloc_acyf(ctx, nn + 2))
            goto SDENCLFin;
        if (alloc_acn(ctx, nn + 2))
            return(-1);   
        if (alloc_acj(ctx, nn + 2))
            return(-1);   
        if (alloc_aci(ctx, nn + 2))
            return(-1);   

        np = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i);                       
            if (typ < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto SDENCLFin;

            for (j = 0; j < n; ++j) {
                np++;
                ctx->AcXF[np] = (float)(ctx->SDVarX[j]);
                ctx->AcYF[np] = (float)(ctx->SDVarY[j]);
                ctx->AcN[np] = np;
            }
        }
        n = g_chull(ctx, ctx->AcXF,ctx->AcYF,np,ctx->AcN,ctx->AcJ,ctx->AcI);
        if (n < 0) {
            printf1(ctx, "Error: insufficient memory for convex hull calculation.\n");
            goto SDENCLFin;
        }
        fprintf(ctx->PMFd,"1 3 %d\n",n);
        nrec = 1;
        k = ctx->AcI[1];
        for (i = 1; i <= n; ++i) {
            j = ctx->AcJ[k];
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->AcXF[j]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->AcYF[j]);
            fprintf(ctx->PMFd,"\n");
            nrec++;
            k = ctx->AcI[k];
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDENCLFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdsel       Selection of spatial objects.                               */
/*                                                                          */  
/*              sdsel(                                                      */
/*                  rec=...,    rectangle: xmin,ymin,xmax,ymax              */
/*                  dxa=...,    offset for geographical coordinates         */
/*                  fmt=...,    print format for coordinates, def. 10.4     */
/*              ) = output_file;                                            */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. The output file is again a spatial data file containing */
/*  all points, lines and polygons that fall in the specified rectangle.    */  
/*  lines and polygons are properly clipped.                                */
/*                                                                          */
/*  The dxa parameter is an option for geographical coordinates. If         */
/*  dxa = x, then for all x coordinates will be added the value 360.        */
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int sdsel(TDAContext *ctx)
{
    int err,nrec,n1,n2,n3;

    err = -1;
    printf1(ctx, "Selection of spatial objects. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto SDLSELFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    nrec = n1  = n2 = n3 = 0;

    if (ctx->PMRECFlg) {
        printf1(ctx, "Region defined by rectangle: %g, %g, %g, %g\n",
            ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax);
        if (ctx->PMRECXMin >= ctx->PMRECXMax - ctx->EPSI1 || ctx->PMRECYMin >= ctx->PMRECYMax - ctx->EPSI1) {
            printf1(ctx, "This is not a valid rectangle.\n");
            goto SDLSELFin;
        }
        nrec = sdsel_rec(ctx, ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax,&n1,&n2,&n3);
        if (nrec < 0)
            goto SDLSELFin;
    }
    else {
        printf1(ctx, "Error: need rec parameter.\n");
        goto SDLSELFin;
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    printf1(ctx, "%d points, %d lines, %d polygons.\n",n1,n2,n3);
    err = 0;

SDLSELFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_rec(xmin,ymin,xmax,ymax,n1,n2,n3)                                 */  
/*                                                                          */  
/*  Write points and lines inside rectangle into spatial output file.       */
/*  Return number of records, or -1 if error.                               */
/*                                                                          */ 
/*  Return the number of points in n1, number of lines in n2, and number    */
/*  of polygons in n3, or -1 if error.                                      */

int sdsel_rec(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax, int *n1,int *n2,int *n3)
{
    register int i,j,k;
    int ii,n,nrec,typ,sdid,id,id1,r,r0,r1,f,i1,i2,i3,cflag,ni;
    double x0,y0,x1,y1,xxa,yya,xxb,yyb,xs,ys;
    double xa[4],ya[4],xb[4],yb[4];

    *n1 = *n2 = *n3 = 0;

    if (alloc_acx(ctx, 2 * ctx->SDVarMax + 10))
        return(-1);   
    if (alloc_acy(ctx, 2 * ctx->SDVarMax + 10))
        return(-1);   
    if (alloc_acn(ctx, 2 * ctx->SDVarMax + 10))
        return(-1);   
    if (alloc_acd(ctx, 2 * ctx->SDVarMax + 10))
        return(-1);   

    xa[0] = xmin;
    ya[0] = ymin;
    xb[0] = xmax;
    yb[0] = ymin;

    xa[1] = xmax;
    ya[1] = ymin;
    xb[1] = xmax;
    yb[1] = ymax;

    xa[2] = xmin;
    ya[2] = ymax;
    xb[2] = xmax;
    yb[2] = ymax;

    xa[3] = xmin;
    ya[3] = ymin;
    xb[3] = xmin;
    yb[3] = ymax;

    id = nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) {

        typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i);                       
        if (typ < 1 || typ > 3)
            continue;

        sdid = (int)get_data(ctx, ctx->SDVarSDID,i);                     

        if (typ == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(ctx, i,0,cflag,1)) < 1)
            return(-1);             

        if (ctx->PMDXAFlg) {                         /* adjust x coordinates */
            for (j = 0; j < n; ++j) {
                if (ctx->SDVarX[j] < ctx->PMDXA)
                    ctx->SDVarX[j] += 360.0;
            }
        }


        if (typ == 1) {                         /* point */
            x0 = ctx->SDVarX[0];
            y0 = ctx->SDVarY[0];

            if (sdsel_inside(ctx, xmin,ymin,xmax,ymax,x0,y0)) {         
                id++;
                nrec += sdsel_prn(ctx, 1,1,ctx->SDVarX,ctx->SDVarY,id,1,sdid);
                *n1 += 1;
            }
        }
        else if (typ == 2) {                    /* line */
            id1 = 0;
            k = 0;
            x0 = ctx->SDVarX[0];
            y0 = ctx->SDVarY[0];
            r0 = sdsel_inside(ctx, xmin,ymin,xmax,ymax,x0,y0);           
            if (r0) {
                ctx->AcX[k] = x0;
                ctx->AcY[k] = y0;
                k++;
            }
            for (j = 1; j < n; ++j) {
                x1 = ctx->SDVarX[j];
                y1 = ctx->SDVarY[j];
                r1 = sdsel_inside(ctx, xmin,ymin,xmax,ymax,x1,y1);           

                if (r0) {               /* (x0,y0) inside */
                    if (r1) {           /* (x1,y1) inside */
                        ctx->AcX[k] = x1;
                        ctx->AcY[k] = y1;
                        k++;
                    }
                    else {              /* (x1,y1) outside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(ctx, x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r == 2 || (r == 1 && (xxa != x0 || yya != y0)))
                                break;
                        }
                        if (r == 0) {
                            printf1(ctx, "Error (1) in sdsel algorithm.\n");
                            exit(0);
                        }
                        else if (r == 1) {      /* intersection in a point */
                            xs = xxa;
                            ys = yya;
                        }
                        else {                               /* overlap */
                            if (x0 == xxa && y0 == yya) {       
                                xs = xxb;
                                ys = yyb;
                            }
                            else {
                                xs = xxa;
                                ys = yya;
                            }
                        }
                        ctx->AcX[k] = xs;
                        ctx->AcY[k] = ys;
                        k++;
                        id++;
                        id1++;
                        nrec += sdsel_prn(ctx, 2,k,ctx->AcX,ctx->AcY,id,id1,sdid);
                        *n2 += 1;
                        k = 0;
                    }
                    x0 = x1;
                    y0 = y1;
                    r0 = r1;
                }
                else {                  /* (x0,y0) outside */
                    if (r1) {           /* (x1,y1) inside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(ctx, x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1(ctx, "Error (2) in sdsel algorithm.\n");
                            exit(0);
                        }
                        else if (r == 1) {      /* intersection in a point */
                            xs = xxa;
                            ys = yya;
                        }
                        else {                               /* overlap */
                            if (x1 == xxa && y1 == yya) {
                                xs = xxb;
                                ys = yyb;
                            }
                            else {
                                xs = xxa;
                                ys = yya;
                            }
                        }
                        ctx->AcX[k] = xs;
                        ctx->AcY[k] = ys;
                        k++;
                        ctx->AcX[k] = x1;
                        ctx->AcY[k] = y1;
                        k++;
                        x0 = x1;
                        y0 = y1;
                        r0 = r1;
                    }
                    else {                          /* (x1,y1) outside */
                        f = 0;
                        if (x0 < xmin) {
                            if (x1 >= xmin) {
                                i1 = 3; i2 = 0; i3 = 2;
                                f = 1;
                            }
                        }
                        else if (x0 > xmax) {
                            if (x1 <= xmax) {
                                i1 = 1; i2 = 0; i3 = 2;
                                f = 1;
                            }
                        }
                        else if (y0 < ymin) {
                            if (y1 >= ymin) {
                                i1 = 0; i2 = 3; i3 = 1;
                                f = 1;
                            }
                        }
                        else if (y0 > ymax) {
                            if (y1 <= ymax) {
                                i1 = 2; i2 = 3; i3 = 1;
                                f = 1;
                            }
                        }
                        if (f == 1) {
                            if (g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i1],ya[i1],xb[i1],yb[i1],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i2],ya[i2],xb[i2],yb[i2],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i3],ya[i3],xb[i3],yb[i3],&xxa,&yya,&xxb,&yyb)) {

                                ctx->AcX[k] = x0 = xxa;
                                ctx->AcY[k] = y0 = yya;
                                k++;
                                r0 = 1;
                                j--;
                                f = 2;
                            }
                        }
                        if (f != 2) {
                            x0 = x1;
                            y0 = y1;
                            r0 = r1;
                        }
                    }
                }
            }
            if (k >= 2) {
                id++;
                id1++;
                nrec += sdsel_prn(ctx, 2,k,ctx->AcX,ctx->AcY,id,id1,sdid);
                *n2 += 1;
                k = 0;
            }
        }
        else {                              /* ### polygon */

            /* make positive (ccw) orientation */

            x0 = g_pol_area(ctx, n,ctx->SDVarX,ctx->SDVarY);

            if (x0 < 0.0) {
                j = 0;
                k = n - 1;
                while (j < k) {
                    x0 = ctx->SDVarX[j];
                    ctx->SDVarX[j] = ctx->SDVarX[k];
                    ctx->SDVarX[k] = x0;
                    y0 = ctx->SDVarY[j];
                    ctx->SDVarY[j] = ctx->SDVarY[k];
                    ctx->SDVarY[k] = y0;
                    j++; 
                    k--;
                }
            }
            ni = 0;
            id1 = 0;
            k = 0;
            x0 = ctx->SDVarX[k];
            y0 = ctx->SDVarY[k];
            r0 = sdsel_inside(ctx, xmin,ymin,xmax,ymax,x0,y0);           
            if (r0) {
                ctx->AcX[k] = x0;
                ctx->AcY[k] = y0;
                ctx->AcD[k] = 0;
                ctx->AcN[k] = -1;
                k++;
                ni++;
            }
            for (j = 1; j < n; ++j) {
                x1 = ctx->SDVarX[j];
                y1 = ctx->SDVarY[j];
                r1 = sdsel_inside(ctx, xmin,ymin,xmax,ymax,x1,y1);           

                if (r0) {               /* (x0,y0) inside */
                    if (r1) {           /* (x1,y1) inside */
                        ctx->AcX[k] = x1;
                        ctx->AcY[k] = y1;
                        ctx->AcD[k] = 0;
                        ctx->AcN[k] = -1;
                        if (k > 0)  
                            ctx->AcN[k - 1] = k;
                        k++;
                        ni++;
                    }
                    else {              /* (x1,y1) outside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(ctx, x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);

/********
tda_out("r=%d xxa=%g %g x0=%g,%g x1=%g,%g\n",r,xxa,yya,x0,y0,x1,y1);

                            if (r == 2 || (r == 1 && (xxa != x0 || yya != y0)))
                                break;
********/
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1(ctx, "Error (3) in sdsel algorithm.\n");
                            printf1(ctx, "[x0=%20.12f,%20.12f x1=%20.12f,%20.12f]\n",x0,y0,x1,y1);
                            exit(0);
                        }
                        else if (r == 1) {      /* intersection in a point */
                            xs = xxa;
                            ys = yya;
                        }
                        else {                               /* overlap */
                            if (x0 == xxa && y0 == yya) {       
                                xs = xxb;
                                ys = yyb;
                            }
                            else {
                                xs = xxa;
                                ys = yya;
                            }
                        }
                        ctx->AcX[k] = xs;
                        ctx->AcY[k] = ys;
                        ctx->AcD[k] = 1;
                        ctx->AcN[k] = -1;
                        if (k > 0)  
                            ctx->AcN[k - 1] = k;
                        k++;
                    }
                    x0 = x1;
                    y0 = y1;
                    r0 = r1;
                }
                else {                  /* (x0,y0) outside */
                    if (r1) {           /* (x1,y1) inside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(ctx, x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1(ctx, "Error (4) in sdsel algorithm.\n");
                            printf1(ctx, "[x0=%20.16f,%20.16f x1=%20.16f,%20.16f]\n",x0,y0,x1,y1);

                            ii=0;
                            r = g_intersect(ctx, x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            tda_out("r=%d\n",r);

                            printf1(ctx, "[x0=%20.16f,%20.16f x1=%20.16f,%20.16f]\n",
                                   xa[ii],ya[ii],xb[ii],yb[ii]);



                            exit(0);
                        }
                        else if (r == 1) {      /* intersection in a point */
                            xs = xxa;
                            ys = yya;
                        }
                        else {                               /* overlap */
                            if (x1 == xxa && y1 == yya) {
                                xs = xxb;
                                ys = yyb;
                            }
                            else {
                                xs = xxa;
                                ys = yya;
                            }
                        }
                        ctx->AcX[k] = xs;
                        ctx->AcY[k] = ys;
                        ctx->AcD[k] = 1;
                        ctx->AcN[k] = -1;
                        k++;
                        ctx->AcX[k] = x1;
                        ctx->AcY[k] = y1;
                        ctx->AcD[k] = 0;
                        ctx->AcN[k] = -1;
                        ctx->AcN[k - 1] = k;
                        k++;
                        x0 = x1;
                        y0 = y1;
                        r0 = r1;
                    }
                    else {                          /* (x1,y1) outside */
                        f = 0;
                        if (x0 < xmin) {
                            if (x1 >= xmin) {
                                i1 = 3; i2 = 0; i3 = 2;
                                f = 1;
                            }
                        }
                        else if (x0 > xmax) {
                            if (x1 <= xmax) {
                                i1 = 1; i2 = 0; i3 = 2;
                                f = 1;
                            }
                        }
                        else if (y0 < ymin) {
                            if (y1 >= ymin) {
                                i1 = 0; i2 = 3; i3 = 1;
                                f = 1;
                            }
                        }
                        else if (y0 > ymax) {
                            if (y1 <= ymax) {
                                i1 = 2; i2 = 3; i3 = 1;
                                f = 1;
                            }
                        }
                        if (f == 1) {
                            if (g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i1],ya[i1],xb[i1],yb[i1],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i2],ya[i2],xb[i2],yb[i2],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(ctx, x0,y0,x1,y1,
                                  xa[i3],ya[i3],xb[i3],yb[i3],&xxa,&yya,&xxb,&yyb)) {

                                ctx->AcX[k] = x0 = xxa;
                                ctx->AcY[k] = y0 = yya;
                                ctx->AcD[k] = 1;
                                ctx->AcN[k] = -1;
                                k++;
                                r0 = 1;
                                j--;
                                f = 2;
                            }
                        }
                        if (f != 2) {
                            x0 = x1;
                            y0 = y1;
                            r0 = r1;
                        }
                    }
                }
            }
            if (ni == n) {          /* all points inside rectangle */
tda_out("sdid=%d inside\n",sdid);
                id++;
                nrec += sdsel_prn(ctx, 3,n,ctx->AcX,ctx->AcY,id,1,sdid);        
                *n3 += 1;
            }
            else if (ni > 0) {          
tda_out("sdid=%d ni=%d\n",sdid,ni);

                /* Here the result might consist of several distinct
                   polygons, this is handled by the function sdsel_prn_poly   
                   by completing the graph construction and processing ist
                   components. First close the graph. */

                if (fabs(ctx->AcX[k-1] - ctx->AcX[0]) < ctx->EPSI1 &&
                    fabs(ctx->AcY[k-1] - ctx->AcY[0]) < ctx->EPSI1)
                    ctx->AcN[k - 1] = 0;
                else
                    ctx->AcD[k - 1] = 1;
                                 
                f = sdsel_prn_poly(ctx, k,ctx->AcX,ctx->AcY,ctx->AcD,ctx->AcN,xmin,ymin,xmax,ymax,&id,sdid,n3);
                if (f < 0)
                    return(-1);
                nrec += f;
            }
/**
else
tda_out("sdid=%d outside\n",sdid);
**/
        }
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_inside(xmin,ymin,xmax,ymax,x,y)                                   */  
/*                                                                          */  
/*  Return 1 if (x,y) is inside the rectangle, otherwise return 0           */

int sdsel_inside(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax,double x, double y)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)    
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_prn(typ,n,x,y,id,id1,sdid)   print to output file.                */

int sdsel_prn(TDAContext *ctx, int typ,int n,double *x,double *y,int id,int id1,int sdid)         
{
    register int i,nrec;

    fprintf(ctx->PMFd,"%8d %d %6d %8d %6d\n",id,typ,n,sdid,id1);
    nrec = 1;
    for (i = 0; i < n; ++i) {
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,x[i]);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y[i]);
        fprintf(ctx->PMFd,"\n");
        nrec++;
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_prn_poly(n,x,y,d,nf,id,sdid,n3)  Print intersecting polygons.     */
/*                                                                          */

int sdsel_prn_poly(TDAContext *ctx, int n,double *x,double *y,char *d,int *nf,double xmin, double ymin,double xmax,double ymax,int *id,int sdid,int *n3)
{
    register int j,k,l;
    int m,nn,nd,nd1,nrec,ef,id1,idx[4];

    /* Add extremal points of bounding box (if not already present). Also
       check that only points lying on the bounding box (d = 1) might
       have an undefined follower. */

    idx[0] = idx[1] = idx[2] = idx[3] = 1;
    for (j = 0; j < n; ++j) {
        if (nf[j] < 0 && d[j] != 1) {
            printf1(ctx, "Error (5) in sdsel algorithm.\n");
            return(-1);
        }
        if (x[j] == xmin && y[j] == ymin)
            idx[0] = 0;
        else if (x[j] == xmax && y[j] == ymin)
            idx[1] = 0;
        else if (x[j] == xmax && y[j] == ymax)
            idx[2] = 0;
        else if (x[j] == xmin && y[j] == ymax)
            idx[3] = 0;
    }
    nn = n;
    if (idx[0] == 1) {
        x[nn] = xmin;
        y[nn] = ymin;
        nf[nn] = -1;
        d[nn] = 1;
        nn++;     
    }
    if (idx[1] == 1) {
        x[nn] = xmax;
        y[nn] = ymin;
        nf[nn] = -1;
        d[nn] = 1;
        nn++;     
    }
    if (idx[2] == 1) {
        x[nn] = xmax;
        y[nn] = ymax;
        nf[nn] = -1;
        d[nn] = 1;
        nn++;     
    }
    if (idx[3] == 1) {
        x[nn] = xmin;
        y[nn] = ymax;
        nf[nn] = -1;
        d[nn] = 1;
        nn++;     
    }
/**       
    tda_out("VORHER  n=%d nn=%d \n",n,nn);
    for (j = 0; j < nn; ++j) {
        tda_out("j=%4d d=%d nf=%3d x=%f y=%f\n",j,d[j],nf[j],x[j],y[j]);
    }
**/     
    if (alloc_acu(ctx, nn + 1))
        return(-1);   
    if (alloc_acv(ctx, nn + 1))
        return(-1);   
    if (alloc_ack(ctx, nn + 1))
        return(-1);   
    if (alloc_acl(ctx, nn + 1))
        return(-1);   

    nd = 0;
    for (j = 0; j < nn; ++j) {
        if (d[j] == 1) {
            ctx->AcU[nd] = x[j];
            ctx->AcV[nd] = y[j];
            ctx->AcL[nd] = j;
            nd++;
        }
    }
/**
tda_out("VOR SORT\n");
    for (j = 0; j < nd; ++j) {
        tda_out("j=%4d u=%f v=%f l=%d nf=%3d\n",j,AcU[j],AcV[j],AcL[j],nf[AcL[j]]);
    }
**/ 
        
    if (sortdp2c(ctx, nd,ctx->AcU,ctx->AcV,ctx->AcK,xmin,ymin,xmax,ymax))
        return(-1);

     
    newline(ctx);
         
            
    /* for all nodes with nf[j] < 0, find adjacent node */

    nd1 = nd - 1;
    for (j = 0; j < nd1; ++j) {
        k = ctx->AcL[ctx->AcK[j]];
        if (nf[k] < 0)  
            nf[k] = ctx->AcL[ctx->AcK[j + 1]];
    }
    k = ctx->AcL[ctx->AcK[nd1]];
    if (nf[k] < 0)
        nf[k] = ctx->AcL[ctx->AcK[0]];
      
    tda_out("\n\nNach einfuegen  nn=%d  n=%d\n",nn,n);
    for (j = 0; j < nn; ++j) {
        tda_out("j=%4d d=%d nf=%3d x=%f y=%f\n",j,d[j],nf[j],x[j],y[j]);
    }
              
    /* print to output file */

tda_out("Writing polygon %d n=%d \n",*id + 1,n );

    for (j = 0; j < n; ++j) 
        d[j] = 1;

    nrec = id1 = 0;
    while (1) {
        k = -1;
        for (j = 0; j < n; ++j) {
            if (d[j] == 1) {
                k = j;
                break;
            }
        }
        if (k < 0)
            break;
        l = k;

tda_out("Next part k=%d\n",k);
        ef = m = 0;
        while (1) {
            m++;
            d[l] = 0;
            l = nf[l];             
            if (l == k)
                break;
            if (d[l] == 0) {
                printf1(ctx, "Error (6) in sdsel algorithm.\n");
                ef = 1;        
                break;
            }
        }
tda_out("m=%d\n",m);
        if (ef || m < 3)
            continue;

        id1++;
        *id += 1;
        fprintf(ctx->PMFd,"%8d 3 %6d %8d %6d ",*id,m,sdid,id1);
        fprintf(ctx->PMFd,"\n");
        nrec++;
        *n3 += 1;

        l = k;
        while (1) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,x[l]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y[l]);
            fprintf(ctx->PMFd,"\n");
            nrec++;
            l = nf[l];
            if (l == k)
                break;
        }
    }
tda_out("finished\n");
/* return(err);     */
    return(nrec);
}









