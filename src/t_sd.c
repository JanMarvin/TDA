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

/*  functions in t_sd.c */

int sd_getdata(int i,int order,int cflag,int opt);

int sdplot(void); 
void sdplot_draw(int typ,int n,double *x,double *y,int nc);
int sdplot31(void);
int sdplot32(void);
int sdplot33(void);
int sdplot3_prism(double zval,int iz,int *npol,int *nsn,int *nsf);
int sdpdata(void);
int sdinf(void);
int sdinf_prn(int typ,int opt);
int sdencl(void);
int sdsel(void);
int sdsel_rec(double xmin,double ymin,double xmax,double ymax,
    int *n1,int *n2,int *n3);
int sdsel_inside(double xmin,double ymin,double xmax,double ymax,double x,
    double y);
int sdsel_prn(int typ,int n,double *x,double *y,int id,int id1,int sdid);        
int sdsel_prn_poly(int n,double *x,double *y,char *d,int *nf,double xmin,
    double ymin,double xmax,double ymax,int *id,int sdid,int *n3);

/* ------------------------------------------------------------------------ */
/*  sd_getdata(i,order,cflag,opt)                                           */
/*                                                                          */
/*  Get data for object i (data matrix record i) and put the data into      */
/*  standard arrays SDVarX, SDVarY. If order != 0 reverse order.            */
/*  If cflag != 0 create a closed polygon. If opt != 0 error message.       */
/*                                                                          */
/*  Return n = number of points in object, or -1 if error.                  */

int sd_getdata(int i,int order,int cflag,int opt)
{
    register int j,k;
    int fptr,n;
    char *p,buf[121];
    double tmp;

    if (SDVarFDef == 0) {
        printf1("Fatal error in sd_getdata().\n");
        return(-1);
    }
    fptr = (int)get_data(SDVarSDPtr,i);
    n    = (int)get_data(SDVarSDN,i);

    if (fseek(SDVarFd,(long)fptr,SEEK_SET)) {
        n = -1;
        goto SDGDFin;
    }
    for (j = 0; j < n; ++j) {
        if (fgets(buf,120,SDVarFd) == NULL) {
            n = -1;
            goto SDGDFin;
        }                    
        p = skip_b(buf);
        if (sscanf(p,"%lg",&tmp) != 1) {
            n = -1;
            goto SDGDFin;
        }
        SDVarX[j] = tmp;

        p = skip_dbl(p);
        p = skip_b(p);
        if (sscanf(p,"%lg",&tmp) != 1) {
            n = -1;
            goto SDGDFin;
        }
        SDVarY[j] = tmp;
    }
    if (cflag) {
        if (SDVarX[n - 1] != SDVarX[0] || SDVarY[n - 1] != SDVarY[0]) {
            SDVarX[n] = SDVarX[0];
            SDVarY[n] = SDVarY[0];
            n++;
        }
    }
    if (order != 0) {
        j = 0; k = n - 1;
        while (j < k) {
            tmp = SDVarX[j];
            SDVarX[j] = SDVarX[k];
            SDVarX[k] = tmp;
            tmp = SDVarY[j];
            SDVarY[j] = SDVarY[k];
            SDVarY[k] = tmp;
            j++;
            k--;
        }
    }

SDGDFin:
    if (n < 1 && opt)
        printf1("Error: cannot read data for spatial object %d\n",i + 1);
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

int sdplot(void)
{
    register int i,j,k,ik;
    int err,n,r,typ,first,nn,np;
    double x,y,d;

    err = -1;
    printf1("Plot spatial data (2d). Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (check_pcmd(0,2))  /* check for PostScript 2d coordinate system */ 
        return(-1);

    if (parm(CmdBuf + 6,1,0))       /* get parameters */
        goto PLSDFin;

    if (PMOPT > 3)
        PMOPT = 1;

    if (PMS < 1 || PMS > 17)
        PMS = 0;           

    fprintf(PSFd,"\n%%#%d: sdplot\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);

    nn  = 0;    /* number of objects plotted */

    if (PMOPT == 1) {

        for (i = 0; i < NOC; ++i) {

            if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            if (typ == 1 && PMS == 0)  
                continue;
               
            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto PLSDFin;
             
            if (PMLW > 0.0 || PMGSFlg) {
                if (typ != 1)
                    sdplot_draw(typ,n,SDVarX,SDVarY,PMNC);
            }
            if (PMS) {  /* plot symbols */

                fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

                for (j = 0; j < n; ++j) {
                    x = SDVarX[j];
                    y = SDVarY[j];
                    x = ps_2dx(x);
                    y = ps_2dy(y);

                    d = PtMM * PMFS / 2.0;
                    ps_sym(PMS,x,y,d);

                    if (PMNC) {
                        d *= 2;
                        upd_bbox(1,x + d,y + d);
                        upd_bbox(1,x + d,y - d);
                        upd_bbox(1,x - d,y + d);
                        upd_bbox(1,x - d,y - d);
                    }
                }
            }
            nn++;
        }
    }
    else if (PMOPT == 2) {              /* convex hull for all points */
        np = 0;
        for (i = 0; i < NOC; ++i) {
            if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            np += (int)get_data(SDVarSDN,i);
        }
        if (np == 0)
            goto PLSDFin;
          
        if (alloc_acxf(np + 2))
            return(-1);   
        if (alloc_acyf(np + 2))
            return(-1);   
        if (alloc_acn(np + 2))
            return(-1);   
        if (alloc_acj(np + 2))
            return(-1);   
        if (alloc_aci(np + 2))
            return(-1);   

        k = 0;
        for (i = 0; i < NOC; ++i) {
            if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto PLSDFin;

            for (j = 0; j < n; ++j) {
                k++;
                AcXF[k] = SDVarX[j];
                AcYF[k] = SDVarY[j];
                AcN[k] = k;
            }
        }
        if (k != np) {
            printf1("sdplot: k=%d np=%d\n",k,np);
            exit(0);
        }
        n = g_chull(AcXF,AcYF,np,AcN,AcJ,AcI);
        if (n < 0) {
            printf1("Error: insufficient memory for convex hull calculation.\n");
            goto PLSDFin;
        }
        ik = AcI[1];
        first = 0;
        for (i = 1; i <= n; ++i) {
            j = AcJ[ik];
            x = (double)AcXF[j];
            y = (double)AcYF[j];
            ps_2dplot(x,y,first); 
            if (PMNC)
                upd_bbox(0,x,y);
            first = 1;
            ik = AcI[ik];
        }
        j = AcJ[AcI[1]];
        x = (double)AcXF[j];
        y = (double)AcYF[j];
        ps_2dplot(x,y,1);     /* close path */
     
        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            fprintf(PSFd,"gsave\n%4.2f setgray\n",PMGS);    
            fprintf(PSFd,"fill\ngrestore\n");
        }
        if (PMLW > 0.0)
            fprintf(PSFd,"stroke\n");
        nn = NOC;
    }
    else if (PMOPT == 3) {              /* convex hulls for polygons */
                    
        for (i = 0; i < NOC; ++i) {
      
            if ((int)get_data(SDVarSDTyp,i) != 3)  
                continue;

            if ((np = sd_getdata(i,0,0,1)) < 1)
                goto PLSDFin;
          
            if (alloc_acxf(np + 2))
                return(-1);   
            if (alloc_acyf(np + 2))
                return(-1);   
            if (alloc_acn(np + 2))
                return(-1);   
            if (alloc_acj(np + 2))
                return(-1);   
            if (alloc_aci(np + 2))
                return(-1);   
      
            k = 0;
            for (j = 0; j < np; ++j) {
                k++;
                AcXF[k] = SDVarX[j];
                AcYF[k] = SDVarY[j];
                AcN[k] = k;
            }
            n = g_chull(AcXF,AcYF,np,AcN,AcJ,AcI);
            if (n < 0) {
                printf1("Error: insufficient memory for convex hull calculation.\n");
                goto PLSDFin;
            }
            ik = AcI[1];
            first = 0;
            for (k = 1; k <= n; ++k) {
                j = AcJ[ik];
                x = (double)AcXF[j];
                y = (double)AcYF[j];
                ps_2dplot(x,y,first); 
                if (PMNC)
                    upd_bbox(0,x,y);
                first = 1;
                ik = AcI[ik];
            }
            j = AcJ[AcI[1]];
            x = (double)AcXF[j];
            y = (double)AcYF[j];
            ps_2dplot(x,y,1);     /* close path */
     
            if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
                fprintf(PSFd,"gsave\n%4.2f setgray\n",PMGS);    
                fprintf(PSFd,"fill\ngrestore\n");
            }
            if (PMLW > 0.0)
                fprintf(PSFd,"stroke\n");
            nn++;
        }
    }
    fprintf(PSFd,"grestore\n");
    printf1("PostScript output for %d object(s) written to: %s\n",nn,PSFName);
    err = 0;

PLSDFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdplot_draw(typ,n,x,y,nc)                                               */
/*                                                                          */
/*  x[i], y[i] (i=0,n-1) contain the points of a line (typ 2) or            */
/*  polygon (typ 3). This function draws the line, or polygon, assuming     */
/*  2d cartesian coordinates. For type 3 there will be an implicit          */
/*  closepath. If nc != 0 update bounding box.                              */
    
void sdplot_draw(int typ,int n,double *x,double *y,int nc)
{
    register int i,first;

    first = 0;
    for (i = 0; i < n; ++i) {
        ps_2dplot(x[i],y[i],first);
        first = 1;
        if (nc)
            upd_bbox(0,x[i],y[i]);
    }
    if (typ == 3) {
        fprintf(PSFd,"closepath\n");    
        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            fprintf(PSFd,"gsave\n%4.2f setgray\n",PMGS);    
            fprintf(PSFd,"fill\ngrestore\n");
        }
    }
    fprintf(PSFd,"stroke\n");
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

int sdplot31(void)
{
    int err,i,n,nn,np,mode;
    double x,y,z,z1,d;

    err = -1;
    if (check_pcmd(0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1("Plot type 1 objects (3d). Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 8,1,0))       /* get parameters */
        goto PLS31Fin;

    if (PMS < 1 || PMS > 17)
        PMS = 0;           

    mode = 0;
    if (PMZVar >= 0)  
        mode = 1;

    if (mode == 0 && PMS == 0) {
        printf1("Nothing done.\n");
        err = 0;
        goto PLS31Fin;
    }
    nn = 0;    /* number of objects plotted */

    fprintf(PSFd,"\n%%#%d: sdplot31\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);
    if (PMS != 0)
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

    np = 0;
    for (i = 0; i < NOC; ++i) {
        if ((int)get_data(SDVarSDTyp,i) != 1)
            continue;

        np++;
        if ((n = sd_getdata(i,0,0,1)) < 1)
            goto PLS31Fin;
                      
        if (mode == 1) {
            z1 = get_data(PMZVar,i);
            ps_3dplot(SDVarX[0],SDVarY[0],PMZVal,0,0,PMNC);
            ps_3dplot(SDVarX[0],SDVarY[0],z1,1,0,PMNC);
            fprintf(PSFd,"stroke\n");
        }
        if (PMS) {
            if (mode == 1)
                z = z1;
            else
                z = PMZVal;

            ps_3dprj(SDVarX[0],SDVarY[0],z,&x,&y);
            x = ps_2dx(x);
            y = ps_2dy(y);

            d = PtMM * PMFS / 2.0;
            ps_sym(PMS,x,y,d);

            if (PMNC) {
                d *= 2;
                upd_bbox(1,x + d,y + d);
                upd_bbox(1,x + d,y - d);
                upd_bbox(1,x - d,y + d);
                upd_bbox(1,x - d,y - d);
            }
        }
        nn++;
    }
    fprintf(PSFd,"grestore\n");
    printf1("Number of type 1 objects: %d\n",np);
    printf1("PostScript output for %d object(s) written to: %s\n",nn,PSFName);
    err = 0;

PLS31Fin:
    p_clean();
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

int sdplot32(void)
{
    register int i,j;
    int err,n,first,nn,np,mode;
    double x,y,z,z1,d;

    err = -1;
    if (check_pcmd(0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1("Plot type 2 objects (3d). Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 8,1,0))       /* get parameters */
        goto PLS32Fin;

    if (PMS < 1 || PMS > 171)
        PMS = 0;           

    mode = 0;
    if (PMZVar >= 0)  
        mode = 1;

    nn = 0;    /* number of objects plotted */

    if (PMS == 0 && PMLW <= 0.0) {
        printf1("Nothing done.\n");
        err = 0;
        goto PLS32Fin;
    }

    fprintf(PSFd,"\n%%#%d: sdplot32\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);

    if (PMS != 0)
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

    np = 0;
    for (i = 0; i < NOC; ++i) {

        if ((int)get_data(SDVarSDTyp,i) != 2)
            continue;

        np++;
        if ((n = sd_getdata(i,0,0,1)) < 1)
            goto PLS32Fin;

        if (PMLW > 0.0) {
            first = 0;
            for (j = 0; j < n; ++j) {
                ps_3dplot(SDVarX[j],SDVarY[j],PMZVal,first,0,PMNC);
                first = 1;
            }
            fprintf(PSFd,"stroke\n");
        }
        if (PMS) {
            if (mode == 1)
                z = z1;
            else
                z = PMZVal;

            for (j = 0; j < n; ++j) {
                ps_3dprj(SDVarX[j],SDVarY[j],z,&x,&y);
                x = ps_2dx(x);
                y = ps_2dy(y);

                d = PtMM * PMFS / 2.0;
                ps_sym(PMS,x,y,d);

                if (PMNC) {
                    d *= 2;
                    upd_bbox(1,x + d,y + d);
                    upd_bbox(1,x + d,y - d);
                    upd_bbox(1,x - d,y + d);
                    upd_bbox(1,x - d,y - d);
                }
            }
        }
        nn++;
    }
    fprintf(PSFd,"grestore\n");
    printf1("Number of type 2 objects: %d\n",np);
    printf1("PostScript output for %d object(s) written to: %s\n",nn,PSFName);
    err = 0;

PLS32Fin:
    p_clean();
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

int sdplot33(void)
{
    register int i,j;
    int err,n,first,nn,mode,ns,npol,nsn,nsf;
    double x,y,z,z1,d;

    err = -1;
    if (check_pcmd(0,3))  /* check for PostScript 3d coordinate system */ 
        return(-1);

    printf1("Plot type 3 objects (3d). Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 8,1,0))       /* get parameters */
        goto PLS33Fin;

    if (PMS < 1 || PMS > 17)
        PMS = 0;           

    if (PMZVar >= 0) {
        mode = 1;
        if (PMGSFlg == 0 || PMGS < 0.0 || PMGS > 1.0) {
            PMGSFlg = 1;
            PMGS = 1.0;
        }
        PMS = 0;
    }
    else {
        mode = 0;
        if (PMGS < 0.0 || PMGS > 1.0)
            PMGSFlg = 0;
    }
    nn = 0;    /* number of objects plotted */
    nsn = nsf = npol = 0;

    if (mode == 0 && PMS == 0 && PMLW <= 0.0 && PMGSFlg == 0) {
        printf1("Nothing done.\n");
        err = 0;
        goto PLS33Fin;
    }

    fprintf(PSFd,"\n%%#%d: sdplot33\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);

    if (mode == 0) {
        if (PMS != 0)
            fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

        for (i = 0; i < NOC; ++i) {
            if ((int)get_data(SDVarSDTyp,i) != 3)
                continue;

            npol++;
            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto PLS33Fin;

            ns = 0;
            if (PMLW > 0.0 || PMGSFlg != 0) {
                first = 0;
                for (j = 0; j < n; ++j) {
                    ps_3dplot(SDVarX[j],SDVarY[j],PMZVal,first,0,PMNC);
                    first = 1;
                }
                fprintf(PSFd,"closepath\n");    
                if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
                    fprintf(PSFd,"gsave\n%4.2f setgray\n",PMGS);    
                    fprintf(PSFd,"fill\ngrestore\n");
                }
                fprintf(PSFd,"stroke\n");
                ns = 1;
            }
            if (PMS) {
                for (j = 0; j < n; ++j) {
                    ps_3dprj(SDVarX[j],SDVarY[j],PMZVal,&x,&y);
                    x = ps_2dx(x);
                    y = ps_2dy(y);

                    d = PtMM * PMFS / 2.0;
                    ps_sym(PMS,x,y,d);

                    if (PMNC) {
                        d *= 2;
                        upd_bbox(1,x + d,y + d);
                        upd_bbox(1,x + d,y - d);
                        upd_bbox(1,x - d,y + d);
                        upd_bbox(1,x - d,y - d);
                    }
                    ns = 1;
                }
            }
            if (ns)
                nn++;
        }
        fprintf(PSFd,"grestore\n");
    }
    else {
        nn = sdplot3_prism(PMZVal,PMZVar,&npol,&nsn,&nsf);
        if (nn < 0)
            goto PLS33Fin;

        fprintf(PSFd,"grestore\n");
    }
    printf1("Number of type 3 objects: %d\n",npol);
    if (nn > 0)
        printf1("PostScript output for %d object(s) written to: %s\n",nn,PSFName);
    if (nsn > 0)
        printf1("Omitted %d prism(s) with non-positive height.\n",nsn);
    if (nsf > 0)
        printf1("Omitted %d faces with unknown priority.\n",nsf);
    err = 0;

PLS33Fin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdplot3_prism(zval,iz,npol,nsn,nsf)                                     */
/*                                                                          */
/*  Plot the polygons in the current spatial data as prisms. Return         */
/*  number of objects plotted, or -1 if an error occurred. Return the       */
/*  number of polygons in npol, the number of prisms with non-positive      */ 
/*  height in nsn, and the number of omitted faces in nsf.                  */

int sdplot3_prism(double zval,int iz,int *npol,int *nsn,int *nsf)
{
    register int i,j,k,l;
    int npp,nn,np,nx,nx0,ng,m,n,n0,npl,np1;
    double a,zval1,x,y,z,xia,xib,xja,xjb,yia,yib,yja,yjb,xa,ya,xb,yb,yi,yj;
    double tmpi,tmpj,tol;

    tol = 100.0 * EPSI1;

    *npol = *nsn = *nsf;

    /* always with fill */

    if (PMGSFlg == 0 || PMGS < 0.0 || PMGS > 1.0) {
        PMGSFlg = 1;
        PMGS = 1.0;
    }

    /* AcC[i] = 1 for positive orientation (counterclockwise), -1 for
       negative orientation, 0 if not a polygon. */

    if (alloc_acc(NOC + 1))
        return(-1);   

    npp = 0;                /* number of polygons */
    nn  = 0;                /* number of points */

    for (i = 0; i < NOC; ++i) {

        if ((int)get_data(SDVarSDTyp,i) != 3)  
            continue;

        *npol += 1;
        zval1 = get_data(iz,i);          
        if (zval1 <= zval) {
            *nsn += 1;
            continue;
        }

        if ((np = sd_getdata(i,0,1,1)) < 1)
            return(-1);       

        a = g_pol_area(np - 1,SDVarX,SDVarY);
        npp++;
        nn += np;
        if (a >= 0.0)
            AcC[i] = 1;
        else
            AcC[i] = -1;
    }
    if (npp == 0)
        return(0);

    /* create a list containing all segments. AcD has the following meaning:
       AcD = 0 do not draw vertical lines
             1 only draw left vertical line
             2 draw both vertical lines
             3 only draw right vertical line
    */

    if (alloc_acxf(nn + 1))
        return(-1);   
    if (alloc_acuf(nn + 1))
        return(-1);   
    if (alloc_acyf(nn + 1))
        return(-1);   
    if (alloc_acvf(nn + 1))
        return(-1);   
    if (alloc_aczf(nn + 1))
        return(-1);   
    if (alloc_acd(nn + 1))
        return(-1);   

    nx = 0;                                
    for (i = 0; i < NOC; ++i) {
        if (AcC[i] == 0)
            continue;

        l = 0;
        if (AcC[i] == -1)
            l = 1;

        if ((np = sd_getdata(i,l,1,1)) < 1)
            return(-1);       
                 
        zval1 = get_data(iz,i);          
        if (zval1 <= zval)  
            continue;

        nx0 = nx;
        for (j = 1; j < np; ++j) {
            ps_3dprj1(SDVarX[j - 1],SDVarY[j - 1],zval,&xia,&y,&z);
            AcXF[nx] = (float)xia;
            AcYF[nx] = (float)y;
            ps_3dprj1(SDVarX[j],SDVarY[j],zval,&xa,&ya,&z);
            AcUF[nx] = (float)xa;
            AcVF[nx] = (float)ya;
            ps_3dprj1(SDVarX[j],SDVarY[j],zval1,&x,&yb,&z);
            AcZF[nx] = yb - ya;
      
            if (xia <= xa) {
                AcD[nx] = 1;
                if (j > 1) {
                    if (AcD[nx - 1] != 0) {
                        if (AcD[nx - 1] == 1)
                            AcD[nx - 1] = 2;
                    }
                    else {
                        xia = (double)AcXF[nx - 1];
                        yia = (double)AcYF[nx - 1];
                        xib = (double)AcUF[nx - 1];
                        yib = (double)AcVF[nx - 1];
       
                        if (g_left(xia,yia,xib,yib,xa,ya) == 0)  
                            AcD[nx] = 3;
                    }
                }
            }
            else {   
                AcYF[nx] -= tol;
                AcVF[nx] -= tol;

                if (j > 1 && AcD[nx - 1] == 1) {
                    xia = (double)AcXF[nx - 1];
                    yia = (double)AcYF[nx - 1];
                    xib = (double)AcUF[nx - 1];
                    yib = (double)AcVF[nx - 1];
                    if (g_left(xia,yia,xib,yib,xa,ya))   
                        AcD[nx - 1] = 2;
                }
            }             
            nx++;
        }
        if (AcD[nx - 1] == 1 && AcD[nx0] == 0)  
            AcD[nx - 1] = 2;

    }
        
    /* create an edge list for a graph that records the priority   
       relations: i -> j if i must be drawn before j. */
     
    ng = 0;         /* number of nodes */

    for (l = 0; l < 2; ++l) {

        k = 0;
        for (i = 0; i < nx; ++i) {
            xia = (double)AcXF[i];
            xib = (double)AcUF[i];
            yia = (double)AcYF[i];
            yib = (double)AcVF[i];

            for (j = i + 1; j < nx; ++j) {
                xja = (double)AcXF[j];
                xjb = (double)AcUF[j];
                yja = (double)AcYF[j];
                yjb = (double)AcVF[j];

                if (xia > xib) {
                    g_xchxy(&xia,&xib);
                    g_xchxy(&yia,&yib);
                }
                if (xja > xjb) {
                    g_xchxy(&xja,&xjb);
                    g_xchxy(&yja,&yjb);
                }
                if ((xa = xja) < xia)  
                    xa = xia;

                if ((xb = xjb) > xib)  
                    xb = xib;

                if (xa >= xb)    
                    continue;
 
                if ((tmpi = xib - xia) == 0.0)
                    tmpi = EPSI1;

                if ((tmpj = xjb - xja) == 0.0)
                    tmpj = EPSI1;

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
                        AcN[k] = j;
                        AcM[k] = i;
                        k++;
                    }
                }
                else {
                    if (l) {
                        AcN[k] = i;
                        AcM[k] = j;
                        k++;
                    }
                }
                if (l == 0)
                    ng++;
            }
        }
        if (l)
            break;

        if (alloc_acn(ng + 1))
            return(-1);   
        if (alloc_acm(ng + 1))
            return(-1);   
    }

    /* sort the edge list according to nodes in AcN and AcM */                
    
    if (ng > 0) {
        if (sorti2(ng,AcN,AcM))
            return(-1);
    }

    /* create pointer to begin of new nodes in AcK */

    if (alloc_ack(nx + 1))
        return(-1);   

    for (i = 0; i < nx; ++i)
        AcK[i] = -1;

    n0 = -1;
    for (i = 0; i < ng; ++i) {
        if (AcN[i] != n0) {
            n0 = AcN[i];
            AcK[n0] = i;
        }           
    }

    /* first fill all polygons at the base heigth zval */

    for (i = 0; i < NOC; ++i) {
        if (AcC[i] == 0)
            continue;

        l = 0;
        if (AcC[i] == -1)
            l = 1;

        if ((np = sd_getdata(i,l,1,1)) < 1)
            return(-1);       
                 
        fprintf(PSFd,"gsave\n");    

        l = 0;
        for (j = 0; j < np; ++j) {
            ps_3dprj1(SDVarX[j],SDVarY[j],zval,&x,&y,&z);
            ps_2dplot(x,y,l);
            l = 1;
        }
        fprintf(PSFd,"%4.2f setgray\n",PMGS);    
        fprintf(PSFd,"fill\n");
        fprintf(PSFd,"stroke\n");
        fprintf(PSFd,"grestore\n");
    }

    /* plot the faces according to the priorities given by the graph */

    if (alloc_acns(nx + 1))     /* counts indegrees of nodes */
        return(-1);   

    for (i = 0; i < ng; ++i)  
        AcNS[AcM[i]] += 1;

    npl = 0;
    while (1) {
        j = 0;
        for (i = 0; i < nx; ++i) {
            if (AcNS[i] == 0) {
                xia = (double)AcXF[i];
                yia = (double)AcYF[i];
                xib = (double)AcUF[i];
                yib = (double)AcVF[i];
                z   = (double)AcZF[i];

                fprintf(PSFd,"gsave\n");    
                ps_2dplot(xia,yia,0);
                ps_2dplot(xib,yib,1);
                ps_2dplot(xib,yib + z,1);
                ps_2dplot(xia,yia + z,1);
                ps_2dplot(xia,yia,1);
                fprintf(PSFd,"%4.2f setgray\n",PMGS);    
                fprintf(PSFd,"fill\n");
                fprintf(PSFd,"stroke\n");
                fprintf(PSFd,"grestore\n");

                ps_2dplot(xia,yia + z,0);
                ps_2dplot(xib,yib + z,1);
                fprintf(PSFd,"stroke\n");

                if (AcD[i] != 0) {
                    ps_2dplot(xia,yia,0);
                    ps_2dplot(xib,yib,1);
                    fprintf(PSFd,"stroke\n");
                      
                    if (AcD[i] <= 2) {
                        ps_2dplot(xia,yia,0);
                        ps_2dplot(xia,yia + z,1);
                        fprintf(PSFd,"stroke\n");
                    }
                    if (AcD[i] >= 2) {
                        ps_2dplot(xib,yib,0);
                        ps_2dplot(xib,yib + z,1);
                        fprintf(PSFd,"stroke\n");
                    }
                }
                npl++;
        
                j++;
                AcNS[i] = -1;
                k = AcK[i];
                if (k < 0)
                    continue;

                for (l = k; l < ng; ++l) {      /* adjust indegrees */
                    if (AcN[l] != i)
                        break;
                    AcNS[AcM[l]] -= 1;
                }
            }
        }
        if (j == 0)   
            break;
    }   
    if (npl < nx) {
        *nsf = nx - npl;
        printf1("Warning: omitted %d faces.\n",nx - npl);
    }

    /* free memory */
    alloc_acc(0);
    alloc_acd(0);
    alloc_acxf(0);
    alloc_acyf(0);
    alloc_aczf(0);
    alloc_acuf(0);
    alloc_acvf(0);
    alloc_acn(0);
    alloc_acm(0);
    alloc_ack(0);
    alloc_acns(0);

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

int sdpdata(void)
{
    register int i,j;
    int err,n,typ,sdid,nrec,cflag;

    err = -1;
    printf1("Write spatial data. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 7,1,1))       /* get parameters */
        goto SDPDFin;

    if (PMOPT > 3)
        PMOPT = 1;

    cflag = 0;
    if (PMOPT == 3) {
        cflag = 1;
        PMOPT = 2;
    }
    if (PMFmtF == 0)
        pmfmt(10,4);

    nrec = 0;
    for (i = 0; i < NOC; ++i) {

        if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        if (PMOPT == 2 && typ == 1)
            continue;

        sdid = (int)get_data(SDVarSDID,i);                     

        if (PMOPT == 1) {
            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto SDPDFin;

            for (j = 0; j < n; ++j) {
                fprintf(PMFd,VPFmtS[SDVarSDID],(double)sdid);
                fprintf(PMFd,PMFmtS,SDVarX[j]);
                fprintf(PMFd,PMFmtS,SDVarY[j]);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }    
        else if (PMOPT == 2) {
            if ((n = sd_getdata(i,0,cflag,1)) < 1)
                goto SDPDFin;

            for (j = 1; j < n; ++j) {
                fprintf(PMFd,VPFmtS[SDVarSDID],(double)sdid);
                fprintf(PMFd,PMFmtS,SDVarX[j - 1]);
                fprintf(PMFd,PMFmtS,SDVarY[j - 1]);
                fprintf(PMFd,PMFmtS,SDVarX[j]);
                fprintf(PMFd,PMFmtS,SDVarY[j]);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }    
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDPDFin:
    p_clean();
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

int sdinf(void)
{
    int err;

    err = -1;
    printf1("Information about spatial data. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 5,1,0))       /* get parameters */
        goto SDINFFin;

    if (PMOPT != 2)
        PMOPT = 1;

    if (sdinf_prn(1,PMOPT))
        goto SDINFFin;

    if (sdinf_prn(2,PMOPT))
        goto SDINFFin;

    if (sdinf_prn(3,PMOPT))
        goto SDINFFin;

    err = 0;

SDINFFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdinf_prn(typ,opt)                                                      */

int sdinf_prn(int typ,int opt)
{
    register int i,j;
    int np,nt,t,n,first,cflag,npp,npn,np0;
    double xmin,xmax,ymin,ymax,pdmin,pdmax,d;

    cflag = 0;
    if (typ == 3)
        cflag = 1;

    pdmin = DBLMAX;
    pdmax = 0.0;

    first = 1;
    np0 = npp = npn = np = nt = 0;
    for (i = 0; i < NOC; ++i) {
        t = (int)get_data(SDVarSDTyp,i);                       
        if (t != typ)
            continue;
        if ((n = sd_getdata(i,0,cflag,1)) < 1)
            return(-1);             
        if (typ == 3)
            n--;
        nt++;
        np += n;

        if (first) {
            xmin = xmax = SDVarX[0];
            ymin = ymax = SDVarY[0];
            first = 0;
        }
        for (j = 0; j < n; ++j) {
            xmin = dmin(xmin,SDVarX[j]);
            xmax = dmax(xmax,SDVarX[j]);
            ymin = dmin(ymin,SDVarY[j]);
            ymax = dmax(ymax,SDVarY[j]);
        }
        if (typ == 2 || typ == 3) {
            for (j = 1; j < n; ++j) {
                d = g_len(SDVarX[j],SDVarY[j],SDVarX[j-1],SDVarY[j-1]);
                pdmin = dmin(pdmin,d);
                pdmax = dmax(pdmax,d);
            }
        }
        if (PMOPT == 2 && typ == 3) {
            d = g_pol_area(n,SDVarX,SDVarY);
            if (d >= PMEPS)
                npp++;
            else if (d <= -PMEPS)
                npn++;
            else
                np0++;
        }
    }
    printf1("\nNumber of ");
    if (typ == 1)
        printf1("points:");
    else if (typ == 2)
        printf1("lines:");
    else            
        printf1("polygons:");
    printf1(" %d",nt);
    if (nt > 0 && (typ == 2 || typ == 3))
        printf1("  (number of points: %d)",np);
    printf1("\n");
    if (nt == 0)
        return(0);

    printf1("XMin: %18.12f   YMin: %18.12f\n",xmin,ymin);
    printf1("XMax: %18.12f   YMax: %18.12f\n",xmax,ymax);

    if (typ == 2 || typ == 3) {
        printf1("\nMinimal segment length: %18.12f\n",pdmin);
        printf1("Maximal segment length: %18.12f\n",pdmax);
    }
    if (PMOPT == 2 && typ == 3) {
        printf1("\nPolygons with positive orientation: %d\n",npp);
        printf1("Polygons with negative orientation: %d\n",npn);
        printf1("Polygons with almost zero (less than %g) area: %d\n",PMEPS,np0);
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

int sdencl(void)
{
    register int i,j,k;
    int err,n,nn,nrec,typ,first,np;
    double xmin,xmax,ymin,ymax;                

    err = -1;
    printf1("Enclosing boundaries. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto SDENCLFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    printf1("Type of boundary: ");
    if (PMOPT == 2)
        printf1("convex hull.\n\n");
    else {
        printf1("rectangle.\n\n");
        PMOPT = 1;
    }
    nn = imax(SDVarNT,4);

    if (PMOPT == 1) {                   /* rectangle */
        first = 1;
        for (i = 0; i < NOC; ++i) {
            typ  = (int)get_data(SDVarSDTyp,i);                       
            if (typ < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto SDENCLFin;         

            for (j = 0; j < n; ++j) {
                if (first) {
                    xmin = xmax = SDVarX[j];
                    ymin = ymax = SDVarY[j];
                    first = 0;
                }
                else {
                    xmin = dmin(xmin,SDVarX[j]);
                    xmax = dmax(xmax,SDVarX[j]);
                    ymin = dmin(ymin,SDVarY[j]);
                    ymax = dmax(ymax,SDVarY[j]);
                }
            }
        }
        fprintf(PMFd,"1 3 4\n");
        fprintf(PMFd,PMFmtS,xmin);
        fprintf(PMFd,PMFmtS,ymin);
        fprintf(PMFd,"\n");
        fprintf(PMFd,PMFmtS,xmax);
        fprintf(PMFd,PMFmtS,ymin);
        fprintf(PMFd,"\n");
        fprintf(PMFd,PMFmtS,xmax);
        fprintf(PMFd,PMFmtS,ymax);
        fprintf(PMFd,"\n");
        fprintf(PMFd,PMFmtS,xmin);
        fprintf(PMFd,PMFmtS,ymax);
        fprintf(PMFd,"\n");
        nrec = 5;
    }
    else if (PMOPT == 2) {              /* convex hull */

        if (alloc_acxf(nn + 2))
            goto SDENCLFin;
        if (alloc_acyf(nn + 2))
            goto SDENCLFin;
        if (alloc_acn(nn + 2))
            return(-1);   
        if (alloc_acj(nn + 2))
            return(-1);   
        if (alloc_aci(nn + 2))
            return(-1);   

        np = 0;
        for (i = 0; i < NOC; ++i) {
            typ  = (int)get_data(SDVarSDTyp,i);                       
            if (typ < 1 || typ > 3)
                continue;

            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto SDENCLFin;

            for (j = 0; j < n; ++j) {
                np++;
                AcXF[np] = SDVarX[j];
                AcYF[np] = SDVarY[j];
                AcN[np] = np;
            }
        }
        n = g_chull(AcXF,AcYF,np,AcN,AcJ,AcI);
        if (n < 0) {
            printf1("Error: insufficient memory for convex hull calculation.\n");
            goto SDENCLFin;
        }
        fprintf(PMFd,"1 3 %d\n",n);
        nrec = 1;
        k = AcI[1];
        for (i = 1; i <= n; ++i) {
            j = AcJ[k];
            fprintf(PMFd,PMFmtS,(double)AcXF[j]);
            fprintf(PMFd,PMFmtS,(double)AcYF[j]);
            fprintf(PMFd,"\n");
            nrec++;
            k = AcI[k];
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDENCLFin:
    p_clean();
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

int sdsel(void)
{
    register int i;
    int err,n,nrec,n1,n2,n3;

    err = -1;
    printf1("Selection of spatial objects. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto SDLSELFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    nrec = n1  = n2 = n3 = 0;

    if (PMRECFlg) {
        printf1("Region defined by rectangle: %g, %g, %g, %g\n",
            PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);
        if (PMRECXMin >= PMRECXMax - EPSI1 || PMRECYMin >= PMRECYMax - EPSI1) {
            printf1("This is not a valid rectangle.\n");
            goto SDLSELFin;
        }
        nrec = sdsel_rec(PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax,&n1,&n2,&n3);
        if (nrec < 0)
            goto SDLSELFin;
    }
    else {
        printf1("Error: need rec parameter.\n");
        goto SDLSELFin;
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    printf1("%d points, %d lines, %d polygons.\n",n1,n2,n3);
    err = 0;

SDLSELFin:
    p_clean();
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

int sdsel_rec(double xmin,double ymin,double xmax,double ymax,
    int *n1,int *n2,int *n3)
{
    register int i,j,k;
    int ii,n,nrec,typ,sdid,id,id1,r,r0,r1,f,i1,i2,i3,cflag,ni;
    double x0,y0,x1,y1,xxa,yya,xxb,yyb,xs,ys;
    double xa[4],ya[4],xb[4],yb[4];

    *n1 = *n2 = *n3 = 0;

    if (alloc_acx(2 * SDVarMax + 10))
        return(-1);   
    if (alloc_acy(2 * SDVarMax + 10))
        return(-1);   
    if (alloc_acn(2 * SDVarMax + 10))
        return(-1);   
    if (alloc_acd(2 * SDVarMax + 10))
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
    for (i = 0; i < NOC; ++i) {

        typ  = (int)get_data(SDVarSDTyp,i);                       
        if (typ < 1 || typ > 3)
            continue;

        sdid = (int)get_data(SDVarSDID,i);                     

        if (typ == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(i,0,cflag,1)) < 1)
            return(-1);             

        if (PMDXAFlg) {                         /* adjust x coordinates */
            for (j = 0; j < n; ++j) {
                if (SDVarX[j] < PMDXA)
                    SDVarX[j] += 360.0;
            }
        }


        if (typ == 1) {                         /* point */
            x0 = SDVarX[0];
            y0 = SDVarY[0];

            if (sdsel_inside(xmin,ymin,xmax,ymax,x0,y0)) {         
                id++;
                nrec += sdsel_prn(1,1,SDVarX,SDVarY,id,1,sdid);
                *n1 += 1;
            }
        }
        else if (typ == 2) {                    /* line */
            id1 = 0;
            k = 0;
            x0 = SDVarX[0];
            y0 = SDVarY[0];
            r0 = sdsel_inside(xmin,ymin,xmax,ymax,x0,y0);           
            if (r0) {
                AcX[k] = x0;
                AcY[k] = y0;
                k++;
            }
            for (j = 1; j < n; ++j) {
                x1 = SDVarX[j];
                y1 = SDVarY[j];
                r1 = sdsel_inside(xmin,ymin,xmax,ymax,x1,y1);           

                if (r0) {               /* (x0,y0) inside */
                    if (r1) {           /* (x1,y1) inside */
                        AcX[k] = x1;
                        AcY[k] = y1;
                        k++;
                    }
                    else {              /* (x1,y1) outside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r == 2 || (r == 1 && (xxa != x0 || yya != y0)))
                                break;
                        }
                        if (r == 0) {
                            printf1("Error (1) in sdsel algorithm.\n");
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
                        AcX[k] = xs;
                        AcY[k] = ys;
                        k++;
                        id++;
                        id1++;
                        nrec += sdsel_prn(2,k,AcX,AcY,id,id1,sdid);
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
                            r = g_intersect(x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1("Error (2) in sdsel algorithm.\n");
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
                        AcX[k] = xs;
                        AcY[k] = ys;
                        k++;
                        AcX[k] = x1;
                        AcY[k] = y1;
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
                            if (g_intersect(x0,y0,x1,y1,
                                  xa[i1],ya[i1],xb[i1],yb[i1],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(x0,y0,x1,y1,
                                  xa[i2],ya[i2],xb[i2],yb[i2],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(x0,y0,x1,y1,
                                  xa[i3],ya[i3],xb[i3],yb[i3],&xxa,&yya,&xxb,&yyb)) {

                                AcX[k] = x0 = xxa;
                                AcY[k] = y0 = yya;
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
                nrec += sdsel_prn(2,k,AcX,AcY,id,id1,sdid);
                *n2 += 1;
                k = 0;
            }
        }
        else {                              /* ### polygon */

            /* make positive (ccw) orientation */

            x0 = g_pol_area(n,SDVarX,SDVarY);

            if (x0 < 0.0) {
                j = 0;
                k = n - 1;
                while (j < k) {
                    x0 = SDVarX[j];
                    SDVarX[j] = SDVarX[k];
                    SDVarX[k] = x0;
                    y0 = SDVarY[j];
                    SDVarY[j] = SDVarY[k];
                    SDVarY[k] = y0;
                    j++; 
                    k--;
                }
            }
            ni = 0;
            id1 = 0;
            k = 0;
            x0 = SDVarX[k];
            y0 = SDVarY[k];
            r0 = sdsel_inside(xmin,ymin,xmax,ymax,x0,y0);           
            if (r0) {
                AcX[k] = x0;
                AcY[k] = y0;
                AcD[k] = 0;
                AcN[k] = -1;
                k++;
                ni++;
            }
            for (j = 1; j < n; ++j) {
                x1 = SDVarX[j];
                y1 = SDVarY[j];
                r1 = sdsel_inside(xmin,ymin,xmax,ymax,x1,y1);           

                if (r0) {               /* (x0,y0) inside */
                    if (r1) {           /* (x1,y1) inside */
                        AcX[k] = x1;
                        AcY[k] = y1;
                        AcD[k] = 0;
                        AcN[k] = -1;
                        if (k > 0)  
                            AcN[k - 1] = k;
                        k++;
                        ni++;
                    }
                    else {              /* (x1,y1) outside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);

/********
printf("r=%d xxa=%g %g x0=%g,%g x1=%g,%g\n",r,xxa,yya,x0,y0,x1,y1);

                            if (r == 2 || (r == 1 && (xxa != x0 || yya != y0)))
                                break;
********/
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1("Error (3) in sdsel algorithm.\n");
                            printf1("[x0=%20.12f,%20.12f x1=%20.12f,%20.12f]\n",x0,y0,x1,y1);
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
                        AcX[k] = xs;
                        AcY[k] = ys;
                        AcD[k] = 1;
                        AcN[k] = -1;
                        if (k > 0)  
                            AcN[k - 1] = k;
                        k++;
                    }
                    x0 = x1;
                    y0 = y1;
                    r0 = r1;
                }
                else {                  /* (x0,y0) outside */
                    if (r1) {           /* (x1,y1) inside */

                        for (ii = 0; ii < 4; ++ii) {
                            r = g_intersect(x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            if (r)
                                break;
                        }
                        if (r == 0) {
                            printf1("Error (4) in sdsel algorithm.\n");
                            printf1("[x0=%20.16f,%20.16f x1=%20.16f,%20.16f]\n",x0,y0,x1,y1);

                            ii=0;
                            r = g_intersect(x0,y0,x1,y1,
                                xa[ii],ya[ii],xb[ii],yb[ii],&xxa,&yya,&xxb,&yyb);
                            printf("r=%d\n",r);

                            printf1("[x0=%20.16f,%20.16f x1=%20.16f,%20.16f]\n",
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
                        AcX[k] = xs;
                        AcY[k] = ys;
                        AcD[k] = 1;
                        AcN[k] = -1;
                        k++;
                        AcX[k] = x1;
                        AcY[k] = y1;
                        AcD[k] = 0;
                        AcN[k] = -1;
                        AcN[k - 1] = k;
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
                            if (g_intersect(x0,y0,x1,y1,
                                  xa[i1],ya[i1],xb[i1],yb[i1],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(x0,y0,x1,y1,
                                  xa[i2],ya[i2],xb[i2],yb[i2],&xxa,&yya,&xxb,&yyb) ||
                                g_intersect(x0,y0,x1,y1,
                                  xa[i3],ya[i3],xb[i3],yb[i3],&xxa,&yya,&xxb,&yyb)) {

                                AcX[k] = x0 = xxa;
                                AcY[k] = y0 = yya;
                                AcD[k] = 1;
                                AcN[k] = -1;
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
printf("sdid=%d inside\n",sdid);
                id++;
                nrec += sdsel_prn(3,n,AcX,AcY,id,1,sdid);        
                *n3 += 1;
            }
            else if (ni > 0) {          
printf("sdid=%d ni=%d\n",sdid,ni);

                /* Here the result might consist of several distinct
                   polygons, this is handled by the function sdsel_prn_poly   
                   by completing the graph construction and processing ist
                   components. First close the graph. */

                if (fabs(AcX[k-1] - AcX[0]) < EPSI1 &&
                    fabs(AcY[k-1] - AcY[0]) < EPSI1)
                    AcN[k - 1] = 0;
                else
                    AcD[k - 1] = 1;
                                 
                f = sdsel_prn_poly(k,AcX,AcY,AcD,AcN,xmin,ymin,xmax,ymax,&id,sdid,n3);
                if (f < 0)
                    return(-1);
                nrec += f;
            }
/**
else
printf("sdid=%d outside\n",sdid);
**/
        }
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_inside(xmin,ymin,xmax,ymax,x,y)                                   */  
/*                                                                          */  
/*  Return 1 if (x,y) is inside the rectangle, otherwise return 0           */

int sdsel_inside(double xmin,double ymin,double xmax,double ymax,double x,
    double y)
{
    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)    
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_prn(typ,n,x,y,id,id1,sdid)   print to output file.                */

int sdsel_prn(int typ,int n,double *x,double *y,int id,int id1,int sdid)         
{
    register int i,nrec;

    fprintf(PMFd,"%8d %d %6d %8d %6d\n",id,typ,n,sdid,id1);
    nrec = 1;
    for (i = 0; i < n; ++i) {
        fprintf(PMFd,PMFmtS,x[i]);
        fprintf(PMFd,PMFmtS,y[i]);
        fprintf(PMFd,"\n");
        nrec++;
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  sdsel_prn_poly(n,x,y,d,nf,id,sdid,n3)  Print intersecting polygons.     */
/*                                                                          */

int sdsel_prn_poly(int n,double *x,double *y,char *d,int *nf,double xmin,
    double ymin,double xmax,double ymax,int *id,int sdid,int *n3)
{
    register int j,k,l;
    int m,nn,nd,nd1,nrec,ef,id1,idx[4];
    double xa,ya,xb,yb;
int err = 0;

    /* Add extremal points of bounding box (if not already present). Also
       check that only points lying on the bounding box (d = 1) might
       have an undefined follower. */

    idx[0] = idx[1] = idx[2] = idx[3] = 1;
    for (j = 0; j < n; ++j) {
        if (nf[j] < 0 && d[j] != 1) {
            printf1("Error (5) in sdsel algorithm.\n");
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
    printf("VORHER  n=%d nn=%d \n",n,nn);
    for (j = 0; j < nn; ++j) {
        printf("j=%4d d=%d nf=%3d x=%f y=%f\n",j,d[j],nf[j],x[j],y[j]);
    }
**/     
    if (alloc_acu(nn + 1))
        return(-1);   
    if (alloc_acv(nn + 1))
        return(-1);   
    if (alloc_ack(nn + 1))
        return(-1);   
    if (alloc_acl(nn + 1))
        return(-1);   

    nd = 0;
    for (j = 0; j < nn; ++j) {
        if (d[j] == 1) {
            AcU[nd] = x[j];
            AcV[nd] = y[j];
            AcL[nd] = j;
            nd++;
        }
    }
/**
printf("VOR SORT\n");
    for (j = 0; j < nd; ++j) {
        printf("j=%4d u=%f v=%f l=%d nf=%3d\n",j,AcU[j],AcV[j],AcL[j],nf[AcL[j]]);
    }
**/ 
        
    if (sortdp2c(nd,AcU,AcV,AcK,xmin,ymin,xmax,ymax))
        return(-1);

     
    newline();
printf("SORT\n");
    for (j = 0; j < nd; ++j) {
        k = AcK[j];
        printf("j=%4d u=%20.12f v=%20.12f l=%d nf=%3d k=%d \n",j,AcU[k],AcV[k],AcL[k],nf[AcL[k]],k );
    }
         
            
    /* for all nodes with nf[j] < 0, find adjacent node */

    nd1 = nd - 1;
    for (j = 0; j < nd1; ++j) {
        k = AcL[AcK[j]];
        if (nf[k] < 0)  
            nf[k] = AcL[AcK[j + 1]];
    }
    k = AcL[AcK[nd1]];
    if (nf[k] < 0)
        nf[k] = AcL[AcK[0]];
      
    printf("\n\nNach einfuegen  nn=%d  n=%d\n",nn,n);
    for (j = 0; j < nn; ++j) {
        printf("j=%4d d=%d nf=%3d x=%f y=%f\n",j,d[j],nf[j],x[j],y[j]);
    }
              
    /* print to output file */

printf("Writing polygon %d n=%d \n",*id + 1,n );

    for (j = 0; j < n; ++j) 
        d[j] = 1;

    nrec = id1 = 0;
    while (1) {
printf("new while\n");

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

printf("Next part k=%d\n",k);
        ef = m = 0;
        while (1) {
            m++;
            d[l] = 0;
            l = nf[l];             
            if (l == k)
                break;
            if (d[l] == 0) {
                printf1("Error (6) in sdsel algorithm.\n");
                ef = 1;        
err = -1;
                break;
            }
        }
printf("m=%d\n",m);
        if (ef || m < 3)
            continue;

        id1++;
        *id += 1;
        fprintf(PMFd,"%8d 3 %6d %8d %6d ",*id,m,sdid,id1);
        fprintf(PMFd,"\n");
        nrec++;
        *n3 += 1;

        l = k;
        while (1) {
            fprintf(PMFd,PMFmtS,x[l]);
            fprintf(PMFd,PMFmtS,y[l]);
            fprintf(PMFd,"\n");
            nrec++;
            l = nf[l];
            if (l == k)
                break;
        }
    }
printf("finished\n");
/* return(err);     */
    return(nrec);
}









