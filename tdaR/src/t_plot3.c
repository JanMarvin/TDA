/****************************************************************************/
/*  t_plot3                                                                 */
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
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_gf.h"
#include "t_int.h"
#include "t_var.h"
#include "t_sort.h"
#include "t_eval3.h"
#include "t_plot.h"
#include "t_intp.h"
#include "t_gmin.h"
#include "tda_context.h"

/*  functions in t_plot3.c */

double arc_to_degree(TDAContext *ctx, double x);
double degree_to_arc(TDAContext *ctx, double x);
void ps_3dprj(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y);
void ps_3dprj1(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dprj_inv(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dnorm(TDAContext *ctx, double *x1,double *x2,double *x3);
double ps_3dist(TDAContext *ctx, double x1,double x2,double x3,double y1,double y2,double y3);
void ps_3dgeo(TDAContext *ctx, double *x1,double *x2,double *x3);
void ps_3dgeo_inv(TDAContext *ctx, double *x1,double *x2,double *x3);
void ps_3dplot(TDAContext *ctx, double x1,double x2,double x3,int opt,int geo,int nc);
void ps_3dplot1(TDAContext *ctx, double x,double y,int opt,int nc);
int plotp3(TDAContext *ctx, int typ);
int plcurv3(TDAContext *ctx);
int pltext3(TDAContext *ctx);
int plcirc3(TDAContext *ctx);
void pl_circ3(TDAContext *ctx, double x,double y,double z,double r,double s,double a,double b, int hide,double xp,double yp,double zp,int nc);
int plglob3(TDAContext *ctx);
int pl_lon3(TDAContext *ctx, double x,double y,double z,double r,double lon,double lat1,double lat2,int nc);
int pl_lat3(TDAContext *ctx, double x,double y,double z,double r,double lat,double lon1,double lon2,int nc);
int pl_meridian(TDAContext *ctx, double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_meridian1(TDAContext *ctx, double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_parallel(TDAContext *ctx, double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);
int pl_parallel1(TDAContext *ctx, double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);

int plsurf3(TDAContext *ctx);
int plsurf3d(TDAContext *ctx);
int plsurf3p(TDAContext *ctx, int nu,int nv,int mu,int mv);
int plsurf3t(TDAContext *ctx, double x,double y,double z,double xtmax,double ytmax);
int plsurf3t1(TDAContext *ctx, double x,double xtmax);
int plsurf3t2(TDAContext *ctx, double x,double y,double z,int iu,int iv);
int plsurf3dt(TDAContext *ctx, double ua,double va,double du,double dv,float *x,float *y,float *z);
int plsurf3c(TDAContext *ctx, double ax,double ay,double bx,double by,double cx,double cy,double x,double y);
void plsurf3n(TDAContext *ctx, int opt,double x,double y,double z,double x1,double y1,double z1, double xtmax,double ytmax,double *x0,double *y0);
int plsurf3o(TDAContext *ctx, int iu,int iv,int mu,int mv);
void plsurf3ct(TDAContext *ctx, double xtmax,double ytmax);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */

#define SHORTMIN -32760
#define SHORTMAX  32760



/* ------------------------------------------------------------------------ */
/*  arc_to_degree(x)                                                        */

double arc_to_degree(TDAContext *ctx, double x)
{
    (void)ctx;        /* unused: the signature is shared */
    return(180.0 * x / Pi);
}

/* ------------------------------------------------------------------------ */
/*  degree_to_arc(x)                                                        */

double degree_to_arc(TDAContext *ctx, double x)
{
    (void)ctx;        /* unused: the signature is shared */
    return(x * Pi / 180.0);
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj(x1,x2,x3,x,y)  Return projection in x and y.                   */

void ps_3dprj(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y)
{
    *x = ctx->PSPR11 * x1 + ctx->PSPR12 * x2;
    *y = ctx->PSPR21 * x1 + ctx->PSPR22 * x2 + ctx->PSPR23 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj1(x1,x2,x3,x,y,z)  Return projection in x, y and z.             */

void ps_3dprj1(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z)
{
    *x = ctx->PSPR11 * x1 + ctx->PSPR12 * x2 + ctx->PSPR13 * x3;
    *y = ctx->PSPR21 * x1 + ctx->PSPR22 * x2 + ctx->PSPR23 * x3;
    *z = ctx->PSPR31 * x1 + ctx->PSPR32 * x2 + ctx->PSPR33 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj_inv(x1,x2,x3,x,y,z)  Return inverse projection in x, y and z.  */

void ps_3dprj_inv(TDAContext *ctx, double x1,double x2,double x3,double *x,double *y,double *z)
{
    *x = ctx->PSPI11 * x1 + ctx->PSPI12 * x2 + ctx->PSPI13 * x3;
    *y = ctx->PSPI21 * x1 + ctx->PSPI22 * x2 + ctx->PSPI23 * x3;
    *z = ctx->PSPI31 * x1 + ctx->PSPI32 * x2 + ctx->PSPI33 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dnorm(x1,x2,x3)     Return normalized vector.                       */

void ps_3dnorm(TDAContext *ctx, double *x1,double *x2,double *x3)
{
    (void)ctx;        /* unused: the signature is shared */
    double u;

    u = sqrt(*x1 * *x1 + *x2 * *x2 + *x3 * *x3);         
    if (fabs(u) > 0.0) {
        *x1 /= u;
        *x2 /= u;
        *x3 /= u;
    }
}

/* ------------------------------------------------------------------------ */
/*  ps_3dist(x1,x2,x3,y1,y2,y3)   Return distance.                          */

double ps_3dist(TDAContext *ctx, double x1,double x2,double x3,double y1,double y2,double y3)
{
    (void)ctx;        /* unused: the signature is shared */
    x1 -= y1;
    x2 -= y2;
    x3 -= y3;
    return(sqrt(x1 * x1 + x2 * x2 + x3 * x3));         
}

/* ------------------------------------------------------------------------ */
/*  ps_3dgeo(x1,x2,x3)  Change from geographic to cartesian coordinates.    */
/*                      Assume geographic coordinates in degrees.           */

void ps_3dgeo(TDAContext *ctx, double *x1,double *x2,double *x3)
{
    double p,lambda,beta,cosl,cosb,sinl,sinb;

    if (*x3 < ctx->EPSI1) {
        *x1 = *x2 = *x3 = 0.0;
        return;
    }
    if (*x1 > 180)
        *x1 -= 360.0;
    else if (*x1 < -180.0)
        *x1 += 360;

    if (*x2 < -90.0 || *x2 > 90.0) {
        if (*x2 < -90.0)
            *x2 = -180.0 - *x2;
        else if (*x2 > 90.0)
            *x2 = 180.0 - *x2;

        if (*x1 <= 0.0)
            *x1 += 180.0;
        else    
            *x1 -= 180.0;
    }


    if (*x1 < -180.0 || *x1 > 180.0 || *x2 < -90.0 || *x2 > 90.0) {
        printf1(ctx, "Error in geographical coordinates: %g, %g, %g.\n",*x1,*x2,*x3);
        printf1(ctx, "Cannot continue.\n");
        exit(0);        
    }
    p = Pi / 180.0;
    lambda = *x1 * p;
    beta   = *x2 * p;
    cosl = cos(lambda);
    cosb = cos(beta);
    sinl = sin(lambda);
    sinb = sin(beta);

    *x1 = cosl * cosb * *x3;
    *x2 = sinl * cosb * *x3;
    *x3 *= sinb;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dgeo_inv(x1,x2,x3)  Change from cartesian to geographic coordinates */
/*                          Return geographic coordinates in degrees.       */

void ps_3dgeo_inv(TDAContext *ctx, double *x1,double *x2,double *x3)
{
    double r;

    r = *x1 * *x1 + *x2 * *x2;
    if (r < ctx->EPSI1) {
        *x2 = 90.0;
        return;
    }
    *x1 = arc_to_degree(ctx, atan2(*x2,*x1));
    *x2 = arc_to_degree(ctx, Pi / 2.0 - atan2(sqrt(r),*x3));
    *x3 = sqrt(r + *x3 * *x3);
}

/* ------------------------------------------------------------------------ */
/*  ps_3dplot   PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values.                                       */
/*                                                                          */
/*              if geo != 0 translate from geographic into cartesian        */
/*              coordinates. If nc != 0 update bounding box.                */

void ps_3dplot(TDAContext *ctx, double x1,double x2,double x3,int opt,int geo,int nc)
{
    double x,y,px,py;

    if (geo)  
        ps_3dgeo(ctx, &x1,&x2,&x3);  
    ps_3dprj(ctx, x1,x2,x3,&x,&y);

    px = (x - ctx->PA1[0]) / ctx->UXLen;
    py = (y - ctx->PA1[1]) / ctx->UYLen;
    px *= ctx->PSXLen;
    py *= ctx->PSYLen;

    if (!opt)  
        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(ctx->PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(ctx->PSFd,"%5.2f %5.2f ",px,py);

    if (nc)
        upd_bbox(ctx, 1,px,py); 
}

/* ------------------------------------------------------------------------ */
/*  ps_3dplot1  PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values. No projection. The function           */
/*              assumes x,y coordinates in projection space.                */
/*              If nc != 0 update bounding box.                             */

void ps_3dplot1(TDAContext *ctx, double x,double y,int opt,int nc)
{
    double px,py;

    px = (x - ctx->PA1[0]) / ctx->UXLen;
    py = (y - ctx->PA1[1]) / ctx->UYLen;
    px *= ctx->PSXLen;
    py *= ctx->PSYLen;

    if (!opt)  
        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(ctx->PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(ctx->PSFd,"%5.2f %5.2f ",px,py);

    if (nc)
        upd_bbox(ctx, 1,px,py); 
}

/*--------------------------------------------------------------------------*/
/*  plotp3(typ) Plot polygon.                                               */
/*                                                                          */
/*  typ0:       plotp(                                                      */
/*                  geo=...,        1 use geographic coordinates, def. 0    */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  gs=...,         grey scale value                        */
/*                  s=...,          marker symbol                           */
/*                  fs=...,         font size, def. 2mm                     */
/*                  a=...,          end in arrow                            */
/*                  nc=...,         1 if no clipping                        */
/*              ) = x1,y1,z1,...;   points                                  */
/*                                                                          */
/*  typ1:       plot(...) = VX,VY,VZ;                                       */
/*                                                                          */  
/*              same parameters, in addition: sel=..., case selection.      */
/*                                                                          */
/*              geo = 0 : cartesian coordinates                             */  
/*              geo = 1 : geographic coordinates: lon,lat,r                 */
/*                        (-180 <= lon,lat <= 180, r >= 0)                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plotp3(TDAContext *ctx, int typ)
{
    int i,n,nn,err,first,ix1 = 0,ix2 = 0,ix3 = 0;
    double x1,x2,x3,x,y,xa,ya;

    err = -1;
    if (check_pcmd(ctx, 2,3))
        return(-1);

    if (typ == 0) {
        if (parm(ctx, ctx->CmdBuf + 6,7,1))    /* get parameters */
            goto PLP3Fin;

        n = ctx->PMRHSN;
        if (n < 6 || (n / 3) * 3 != n) {
            p_err(ctx, -45,1);
            goto PLP3Fin;
        }
    }
    else {
        if (parm(ctx, ctx->CmdBuf + 5,4,1))    /* get parameters */
            goto PLP3Fin;

        if (ctx->NOC < 2)
            goto PLP3Fin;

        n = ctx->NOC;
        if (ctx->PMNV != 3) {
            p_err(ctx, -1,1);
            goto PLP3Fin;
        }
        ix1 = ctx->PMVIdx[0];
        ix2 = ctx->PMVIdx[1];
        ix3 = ctx->PMVIdx[2];
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);
    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC == 0)
        set_clip(ctx);

    /* first process fill option */

    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        ps_lwidth(ctx, 0.0);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = ctx->PMRHSX[i];
                x2 = ctx->PMRHSX[i + 1];
                x3 = ctx->PMRHSX[i + 2];
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ctx, ix1,i);       
                x2 = get_data(ctx, ix2,i);       
                x3 = get_data(ctx, ix3,i);       
            }
            if (ctx->PMGEO) {
                ps_3dgeo(ctx, &x1,&x2,&x3);
            }
            ps_3dprj(ctx, x1,x2,x3,&x,&y);
            x = ps_2dx(ctx, x);
            y = ps_2dy(ctx, y);

            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);    

            nn++;
            if (first) {
            /*  fprintf(PSFd,"%5.2f %5.2f m\n",x,0.0);  */
                fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                if (ctx->PMDIR == 1)
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,ya);
                else if (ctx->PMDIR == 2)
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",xa,y);

                fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
            }
            xa = x;
            ya = y;

            i++;   
            if (typ == 0)
                i += 2;
        }
        if (nn) {
            /** fprintf(PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);   **/
            fprintf(ctx->PSFd,"gsave\n");
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\nnewpath\n");
        }
    }

    /* now the standard polygon */

    if (ctx->PMLW > 0.0)
        ps_lwidth(ctx, ctx->PMLW);

    if (ctx->PMLW > 0.0 && ctx->PMLT > 0) {
        ps_ltyp(ctx, ctx->PMLT);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = ctx->PMRHSX[i];
                x2 = ctx->PMRHSX[i + 1];
                x3 = ctx->PMRHSX[i + 2];
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ctx, ix1,i);       
                x2 = get_data(ctx, ix2,i);       
                x3 = get_data(ctx, ix3,i);       
            }
            if (ctx->PMGEO) {
                ps_3dgeo(ctx, &x1,&x2,&x3);
            }
            ps_3dprj(ctx, x1,x2,x3,&x,&y);
            x = ps_2dx(ctx, x);
            y = ps_2dy(ctx, y);

            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);    

            nn++;
            if (first) {
                fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
            }
            i++;   
            if (typ == 0)
                i += 2;
            if (i >= n)
                break;
            xa = x;
            ya = y;
        }
        if (nn >= 2 && ctx->PMAFlg && ctx->PMA1 > 0.0 && ctx->PMA2 > 0.0) {       /* plot arrow */
            x -= xa;
            y -= ya;
            if (fabs(y) > ctx->EPSI1) {
                fprintf(ctx->PSFd,"gsave\ncurrentpoint\nstroke m\n");
                fprintf(ctx->PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
                fprintf(ctx->PSFd,"%5.2f %5.2f scale\n",ctx->PMA1,ctx->PMA2);
                fprintf(ctx->PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
            }
        }
        if (nn)
            fprintf(ctx->PSFd,"stroke\n");
    }
    if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

        i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = ctx->PMRHSX[i];
                x2 = ctx->PMRHSX[i + 1];
                x3 = ctx->PMRHSX[i + 2];
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ctx, ix1,i);       
                x2 = get_data(ctx, ix2,i);       
                x3 = get_data(ctx, ix3,i);       
            }
            if (ctx->PMGEO) {
                ps_3dgeo(ctx, &x1,&x2,&x3);
            }
            ps_3dprj(ctx, x1,x2,x3,&x,&y);
            x = ps_2dx(ctx, x);
            y = ps_2dy(ctx, y);

            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);    

            ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            i++;   
            if (typ == 0)
                i += 2;
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLP3Fin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  plcurv3.        Plot curve.                                             */
/*                                                                          */
/*  plcurv3(                                                                */
/*      rx=...,     range of argument: a(d)b                                */
/*      f1=...,     x1 = f(x)                                               */
/*      f2=...,     x2 = f(x)                                               */
/*      f3=...,     x3 = f(x)                                               */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      nc=...,     1 if no clipping                                        */
/*  );                                                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plcurv3(TDAContext *ctx)
{
    register int i;
    int err,r,f,n,np;      
    double x,y,tmp;
           
    err = -1;
    if (check_pcmd(ctx, 0,3))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,1,0))    /* get parameters */
        goto PLCFin;

    if (ctx->PMF1A <= 0 || ctx->PMF2A <= 0 || ctx->PMF3A <= 0 || ctx->PMRXFlg == 0) {
        printf1(ctx, "Need arguments: rx, f1, f2, f3.\n");
        goto PLCFin;
    }
    np = 0;
    x = (double)(ctx->PMRXA);
    while ((double)(x) <= (double)((ctx->PMRXB)) + (double)(ctx->PMRXD)) {
        np++;
        x += (double)(ctx->PMRXD);
    }
    if (alloc_acx(ctx, np + 1))
        goto PLCFin;
    if (alloc_acy(ctx, np + 1))
        goto PLCFin;
    if (alloc_acz(ctx, np + 1))
        goto PLCFin;

    for (f = 1; f <= 3; ++f) {
        if (f == 1)
            r = get_func(ctx, ctx->PMF1,1,&n,0,"f1");
        else if (f == 2)
            r = get_func(ctx, ctx->PMF2,1,&n,0,"f2");
        else                          
            r = get_func(ctx, ctx->PMF3,1,&n,0,"f3");
        if (r)
            goto PLCFin;
        if (n) {
            get_func(ctx, NULL,0,&n,0,NULL);
            p_err(ctx, -40,1);
            goto PLCFin;
        }
        if (ctx->FNArgN > 1) {
            printf1(ctx, "Error: at most one function argument (x).\n");
            goto PLCFin;
        }
        x = (double)(ctx->PMRXA);
        for (i = 0; i < np; ++i) {
            if (ctx->FNArgN == 1)
                ctx->FNArgVal[0] = x;    
            r = get_flval(ctx, &y,ctx->FNArgN,ctx->FNArgVal,0,1,&tmp,&tmp,&tmp);
            if (r) {      /* r from v_eval1(ctx) */
                printf1(ctx, "Can't evaluate function.\n");
                prn_emsg2(ctx, r);
                goto PLCFin;
            }
            if (f == 1)
                ctx->AcX[i] = y;
            else if (f == 2)
                ctx->AcY[i] = y;
            else
                ctx->AcZ[i] = y;

            if ((double)((((x += (double)(ctx->PMRXD))))) > (double)(ctx->PMRXB))
                x = (double)(ctx->PMRXB);
        }
        get_func(ctx, NULL,0,&n,0,NULL);
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC == 0)
        set_clip(ctx);

    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    ps_3dplot(ctx, ctx->AcX[0],ctx->AcY[0],ctx->AcZ[0],0,0,0);

    for (i = 1; i < np; ++i)  
        ps_3dplot(ctx, ctx->AcX[i],ctx->AcY[i],ctx->AcZ[i],1,0,ctx->PMNC);
       
    fprintf(ctx->PSFd,"stroke\ngrestore\n");
    err = 0;

PLCFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pltext3.    Plot text string.                                           */
/*                                                                          */
/*              pltext3(                                                    */
/*                  xyz=            position: x,y,z                         */
/*                  fs =...,        font size, def. 2mm                     */
/*                  sc =...,        1 = centered, def. 0                    */
/*                  r  =...,        rotation, def. 0                        */
/*                  s  =...,        marker symbol                           */
/*                  wf = ...,       if wf=1 put string into white           */  
/*                                  bounding box, def. 0                    */
/*              ) = string;                                                 */ 
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pltext3(TDAContext *ctx)
{
    int cen,err = -1;
    double x,y,px,py,l;

    if (check_pcmd(ctx, 0,3))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,8,1))    /* get parameters */
        goto PLT3Fin;

    if (ctx->PMXYZFlg == 0) {
        p_err(ctx, -1,1);
        goto PLT3Fin;
    }
    if (ctx->PMWF != 0)
        ctx->PMWF = 1;

    if (ctx->PMR)  
        cen = 0;
    else  
        cen = (int)ctx->PMSC;
       
    ps_3dprj(ctx, ctx->PM3X,ctx->PM3Y,ctx->PM3Z,&x,&y);

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    px = ps_2dx(ctx, x);
    py = ps_2dy(ctx, y);

    if (ctx->PMS) {
        ctx->PMR = cen = 0;
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);
        ps_sym(ctx, ctx->PMS,px - 1.5 * ctx->PtMM * ctx->PMFS,py + 0.5 * ctx->PtMM * ctx->PMFS,ctx->PtMM * ctx->PMFS / 2.0);
    }
    plot_str(ctx, x,y,ctx->PMRHSTR,1.5 * ctx->PMFS * ctx->PtMM,0,cen,ctx->PMR,0,ctx->PMWF);

    /* update bounding box */

    l = 0.6 * (double)strlen(ctx->PMRHSTR) * ctx->PMFS * ctx->PtMM;
    upd_bbox(ctx, 1,px,py);     
    upd_bbox(ctx, 1,px + l,py + l);     

    upd_bbox(ctx, 1,px,py + 1.3 * ctx->PMFS * ctx->PtMM);
    upd_bbox(ctx, 1,px,py - 0.5 * ctx->PMFS * ctx->PtMM);

    if (ctx->PMR > 0)  
        upd_bbox(ctx, 1,px - 1.3 * ctx->PMFS * ctx->PtMM,y);

    if (ctx->PMR < 0)  
        upd_bbox(ctx, 1,px + 1.3 * ctx->PMFS * ctx->PtMM,y);

    err = 0;

PLT3Fin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  plcirc3()   Plot circle.                                                */
/*                                                                          */
/*              plcirc3(                                                    */
/*                  xyz=...,        center                                  */
/*                  dvec=...,       direction vector                        */
/*                  geo=...,        if 1 dir vector in geogr. coordinates   */
/*                  hide=...,        0 show all, def. 0                     */
/*                                  -1 don't show hidden lines              */
/*                                   i > 0 line type for hidden parts       */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  nc=...,         1 if no clipping                        */
/*              ) = r [a,b];        radius                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plcirc3(TDAContext *ctx)
{
    int err;
    double a,b,r,s;
           
    err = -1;
    if (check_pcmd(ctx, 0,3))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,7,1))    /* get parameters */
        goto PLCIFin;

    if (ctx->PMXYZFlg == 0) {
        p_err(ctx, -1,1);
        goto PLCIFin;
    }

    if (ctx->PMRHSN == 1) {
        a = 0.0;
        b = 360.0;
    }
    else if (ctx->PMRHSN == 3) {
        a = ctx->PMRHSX[1];
        b = ctx->PMRHSX[2];
    }      
    else {
        p_err(ctx, -1,1);
        goto PLCIFin;
    }
    r = ctx->PMRHSX[0];
    if (r <= 0.0 || a < -360.0 || a > b || b > 360.0) {
        p_err(ctx, -1,1);
        goto PLCIFin;
    }
    if (!ctx->PMDVECFlg) {
        p_err(ctx, -1,1);
        goto PLCIFin;
    }
    if (ctx->PMGEO) {
        ps_3dgeo(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);     
    }
    ps_3dnorm(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);        
    s = sqrt(ctx->PMDVECX * ctx->PMDVECX + ctx->PMDVECY * ctx->PMDVECY);

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);
    pl_circ3(ctx, ctx->PM3X,ctx->PM3Y,ctx->PM3Z,r,s,a,b,ctx->PMHIDE,ctx->PM3X,ctx->PM3Y,ctx->PM3Z,ctx->PMNC);

    err = 0;

PLCIFin:
    p_clean(ctx);
    return(err);
}
                                              

/*--------------------------------------------------------------------------*/
/*  pl_circ3(x,y,z,r,s,a,b,hide,xp,yp,zp,nc)                                */
/*                                                                          */
/*  Plot circle with radius r and center x,y,z from a to b. xp,yp,zp is the */   
/*  center of the projection to decide what is hidden. If nc != 0 update    */
/*  bounding box (no clipping).                                             */

void pl_circ3(TDAContext *ctx, double x,double y,double z,double r,double s,double a,double b, int hide,double xp,double yp,double zp,int nc)
{
    int i,cflag,hflag,first,lt;
    double u,x1,x2,x3,h,hh,cosu,sinu;

    hh = ctx->PSPR31 * xp + ctx->PSPR32 * yp + ctx->PSPR33 * zp;
    lt = ctx->PMLT;
    cflag = hflag = first = 0;
    for (i = -360; i < 360; ++i) {
        if ((double)i < a)   
            continue;
        else if ((double)i > b)
            break;
        u = (double)i * Pi / 180.0;
        cosu = cos(u);
        sinu = sin(u);
        if (s > ctx->EPSI1) {
            x1 = x - (ctx->PMDVECY * cosu + ctx->PMDVECX * ctx->PMDVECZ * sinu) * r / s;
            x2 = y + (ctx->PMDVECX * cosu - ctx->PMDVECY * ctx->PMDVECZ * sinu) * r / s;
            x3 = z + r * s * sinu;          
        }
        else {
            x1 = x + r * cosu;
            x2 = y + r * sinu;
            x3 = z;
        }
        if (hide < 0) {
            h = ctx->PSPR31 * x1 + ctx->PSPR32 * x2 + ctx->PSPR33 * x3;
            if (h < hh -ctx->EPSI1) {
                if (first) {
                    fprintf(ctx->PSFd,"stroke\ngrestore\n");
                    hflag = first = 0;
                }
                continue;
            }
        }
        else if (hide > 0) {
            h = ctx->PSPR31 * x1 + ctx->PSPR32 * x2 + ctx->PSPR33 * x3;
            if (h < hh -ctx->EPSI1) {
                if (cflag == 0) {
                    if (first) {
                        fprintf(ctx->PSFd,"stroke\ngrestore\n");
                        hflag = first = 0;
                    }
                    lt = hide;
                    cflag = 1;
                }
            }
            else if (cflag) {
                fprintf(ctx->PSFd,"stroke\ngrestore\n");
                lt = ctx->PMLT;
                first = cflag = hflag = 0;
            }
        }
        if (first == 0) {
            fprintf(ctx->PSFd,"gsave\n");
            if (nc == 0)
                set_clip(ctx);
            ps_ltyp(ctx, lt);
            ps_lwidth(ctx, ctx->PMLW);
            hflag = 1;
        }
        ps_3dplot(ctx, x1,x2,x3,first,0,nc);
        first = 1;
    }
    if (hflag)
        fprintf(ctx->PSFd,"stroke\ngrestore\n");
}

/*--------------------------------------------------------------------------*/
/*  pl_meridian(lon,lata,latb,r,n,nc,lt,lw)                                 */
/*                                                                          */
/*  plot meridian at longitude lon from latitude lata to latb. r is radius. */
/*  n is the number of steps, nc is the no clipping option. lt and lw are,  */
/*  respectively, the line typ and line with for drawing the lines.         */
/*                                                                          */
/*  -180 <= lata < latb <= 180. Max range from lata to latb is 180 degrees. */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_meridian(TDAContext *ctx, double lon,double lata,double latb,double r,int n,int nc,int lt,double lw)
{
    fprintf(ctx->PSFd,"gsave\n");
    if (nc == 0)
        set_clip(ctx);
    ps_ltyp(ctx, lt);
    ps_lwidth(ctx, lw);

      
    /** tda_out("meridian lon=%g lata=%g latb=%g\n",lon,lata,latb); **/


    if (latb > 90.0) {
        if (pl_meridian1(ctx, lon,lata,90.0,r,n,nc,lt,lw))
            return(-1);

        if ((lon += 180.0) > 180)
            lon -= 360.0;

        if (pl_meridian1(ctx, lon,180.0 - latb,90.0,r,n,nc,lt,lw))
            return(-1);
    }
    else if (lata < -90.0) {

        if (pl_meridian1(ctx, lon,latb,-90.0,r,n,nc,lt,lw))
            return(-1);

        if ((lon += 180.0) > 180)
            lon -= 360.0;

        if (pl_meridian1(ctx, lon,180.0 + lata,-90.0,r,n,nc,lt,lw))
            return(-1);
    }
    else {
        if (pl_meridian1(ctx, lon,lata,latb,r,n,nc,lt,lw))
            return(-1);
    }
    fprintf(ctx->PSFd,"grestore\n");
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_meridian1(lon,lata,latb,r,n,nc,lt,lw)                                */
/*                                                                          */
/*  plot meridian at longitude lon from latitude lata to latb. r is radius. */
/*  n is the number of steps, nc is the no clipping option. lt and lw are,  */
/*  respectively, the line typ and line with for drawing the lines.         */
/*                                                                          */
/*  -90 <= lata < latb <= 90.                                               */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_meridian1(TDAContext *ctx, double lon,double lata,double latb,double r,int n,int nc,int lt,double lw)
{
    (void)lt; (void)lw;        /* unused: the signature is shared */
    int i,first;
    double delta,x,y,z,h;               

    delta = (latb - lata) / (double)n;

    first = 0;
    for (i = 0; i <= n; ++i) {
        x = lon;
        y = lata + (double)i * delta;
        z = r;
        ps_3dgeo(ctx, &x,&y,&z);  

        h = ctx->PSPR31 * x + ctx->PSPR32 * y + ctx->PSPR33 * z;

        if (h < -ctx->EPSI1) {           /* plot only if visible */
            if (first) {
                fprintf(ctx->PSFd,"stroke\n");
                first = 0;
            }
        }
        else {
            ps_3dplot(ctx, x,y,z,first,0,nc);
            first = 1;
        }
    }
    if (first)
        fprintf(ctx->PSFd,"stroke\n");
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_parallel(lat,lona,lonb,r,n,nc,lt,lw)                                 */
/*                                                                          */
/*  plot parallel at latitude lat from longitude lona to lonb. r is radius. */
/*  n is the number of steps, nc is the no clipping option. lt and lw are,  */
/*  respectively, the line type and line width for drawing the lines.       */
/*                                                                          */
/*  -270 <- lona < lonb <= +270.                                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_parallel(TDAContext *ctx, double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw)
{


    fprintf(ctx->PSFd,"gsave\n");
    if (nc == 0)
        set_clip(ctx);
    ps_ltyp(ctx, lt);
    ps_lwidth(ctx, lw);

    if (lona <= lonb) {
        if (pl_parallel1(ctx, lat,lona,lonb,r,n,nc,lt,lw))
            return(-1);
    }
    else {
        if (pl_parallel1(ctx, lat,lona,180.0,r,n,nc,lt,lw))
            return(-1);

        if (pl_parallel1(ctx, lat,-180.0,lonb,r,n,nc,lt,lw))
            return(-1);
    }
    fprintf(ctx->PSFd,"grestore\n");
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_parallel1(lat,lona,lonb,r,n,nc,lt,lw)                                */
/*                                                                          */
/*  plot parallel at latitude lat from longitude lona to lonb. r is radius. */
/*  n is the number of steps, nc is the no clipping option. lt and lw are,  */
/*  respectively, the line type and line width for drawing the lines.       */

int pl_parallel1(TDAContext *ctx, double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw)
{
    (void)lt; (void)lw;        /* unused: the signature is shared */
    int i,first;
    double delta,x,y,z,h;               

    delta = (lonb - lona) / (double)n;

    first = 0;
    for (i = 0; i <= n; ++i) {

        x = lona + (double)i * delta;
        y = lat;
        z = r;
        ps_3dgeo(ctx, &x,&y,&z);  

        h = ctx->PSPR31 * x + ctx->PSPR32 * y + ctx->PSPR33 * z;

        if (h < -ctx->EPSI1) {           /* plot only if visible */
            if (first) {
                fprintf(ctx->PSFd,"stroke\n");
                first = 0;
            }
        }
        else {
            ps_3dplot(ctx, x,y,z,first,0,nc);
            first = 1;
        }
    }
    if (first)
        fprintf(ctx->PSFd,"stroke\n");
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  plglob3()   Plot globe.                                                 */
/*                                                                          */
/*              plglob3(                                                    */
/*                  xyz=...,        center, def. 0,0,0                      */
/*                  lon=...,        sequence of longitudes                  */
/*                  lat=...,        sequence of latitudes                   */
/*                  lt=...,         line type, def. 1                       */
/*                  lw=...,         line width, def. 0.2 mm                 */
/*                  cont=...,       1 add contour, def. 0                   */
/*              ) = r;              radius                                  */
/*                                                                          */
/*  lon and lat can be given as any sequence of floating point values. Also */
/*  possible is the ,, notation. However, values foor lon must be in the    */
/*  range -180,...,180, and values for lat must be in the range -90,...,90. */ 
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plglob3(TDAContext *ctx)
{
    int err,i;
    double r,s,x;  
           
    err = -1;
    if (check_pcmd(ctx, 0,3))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,7,1))    /* get parameters */
        goto PLGLFin;

    if ((r = ctx->PMRHSX[0]) <= 0.0) {
        p_err(ctx, -1,1);
        goto PLGLFin;
    }

    /* plot contour of sphere */
      
    ctx->PMDVECX = ctx->PSLon;
    ctx->PMDVECY = ctx->PSLat;
    ctx->PMDVECZ = r;       
    ps_3dgeo(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);          
    ps_3dnorm(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);        
    s = sqrt(ctx->PMDVECX * ctx->PMDVECX + ctx->PMDVECY * ctx->PMDVECY);
                      
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    if (ctx->PMCONT)
        pl_circ3(ctx, ctx->PM3X,ctx->PM3Y,ctx->PM3Z,r,s,0.0,360.0,0,ctx->PM3X,ctx->PM3Y,ctx->PM3Z,0);
     
    for (i = 0; i < ctx->PMNTP; ++i) {
        x = ctx->PMTP[i];
        if (x < -180.0 || x > 180.0)
            continue;

        if (pl_lon3(ctx, ctx->PM3X,ctx->PM3Y,ctx->PM3Z,r,x,0.0,360.0,0))
            goto PLGLFin;
    }
    for (i = 0; i < ctx->PMNTP1; ++i) {
        x = ctx->PMTP1[i];
        if (x < -90.0 || x > 90.0)
            continue;

        if (pl_lat3(ctx, ctx->PM3X,ctx->PM3Y,ctx->PM3Z,r,x,0.0,360.0,0))
            goto PLGLFin;
    }
    err = 0;

PLGLFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_lon3(x,y,z,r,lon,lat1,lat2,nc)                                       */
/*                                                                          */
/*  Plot circle with radius r and longitude lon from lat1 to lat2.          */
/*  x,y,z is center of globe. if nc != 0 then no clipping.                  */
/*  Return 0 if OK, -1 if error.                                            */

int pl_lon3(TDAContext *ctx, double x,double y,double z,double r,double lon,double lat1,double lat2,int nc)
{
    double s;

    if (lon > 90.0)
        lon -= 90.0;
    else
        lon += 90.0;

    ctx->PMDVECX = lon;
    ctx->PMDVECY = 0.0;
    ctx->PMDVECZ = r;
    ps_3dgeo(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);    
    ps_3dnorm(ctx, &ctx->PMDVECX,&ctx->PMDVECY,&ctx->PMDVECZ);        
    s = sqrt(ctx->PMDVECX * ctx->PMDVECX + ctx->PMDVECY * ctx->PMDVECY);
/*  pl_circ3(x,y,z,r,s,lat1,lat2,-1,PM3X,PM3Y,PM3Z,nc);     */
    pl_circ3(ctx, x,y,z,r,s,lat1,lat2,-1,x,y,z,nc);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_lat3(x,y,z,r,lat,lon1,lon2,nc)                                       */
/*                                                                          */
/*  Plot circle with radius r and latitude lat (from lon1 to lon2)          */
/*  x,y,z is center of globe. If nc != 0 then no clipping.                  */
/*  Return 0 if OK, -1 if error.                                            */

int pl_lat3(TDAContext *ctx, double x,double y,double z,double r,double lat,double lon1,double lon2,int nc)
{
    double s,u,z1;

    u = lat * Pi / 180.0;
    z1 = z + sin(u) * r;
    r *= cos(u);
    ctx->PMDVECX = 0.0;    
    ctx->PMDVECY = 0.0;   
    ctx->PMDVECZ = 1;   
    s = sqrt(ctx->PMDVECX * ctx->PMDVECX + ctx->PMDVECY * ctx->PMDVECY);
    pl_circ3(ctx, x,y,z1,r,s,lon1,lon2,-1,x,y,z,nc);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  plsurf3         Plot a parametrically defined surface.                  */
/*                                                                          */
/*  plsurf3(                                                                */
/*      ru=...,     range of first argument: ua,ub,nu,mu                    */
/*      rv=...,     range of second argument: va,vb,nv,mv                   */
/*      f1=...,     x1 = f(u,v)                                             */
/*      f2=...,     x2 = f(u,v)                                             */
/*      f3=...,     x3 = f(u,v)                                             */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      cont=...,   1 add internal contour, def. 0                          */
/*      gs=...,     grey scale value                                        */
/*      nc=...,     1 if no clipping                                        */
/*  );                                                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plsurf3(TDAContext *ctx)
{
    int err,nu,mu,nv,mv;
    double ua,ub,va,vb,du,dv;

    err = -1;
    if (check_pcmd(ctx, 0,3))
        return(-1);

    printf1(ctx, "Plot parametrically defined surface. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,1,0))    /* get parameters */
        goto PLSFin;

    if (ctx->PMF1A <= 0 || ctx->PMF2A <= 0 || ctx->PMF3A <= 0 || ctx->PMRUFlg == 0 || ctx->PMRVFlg == 0) {
        printf1(ctx, "Need arguments: ru, rv, f1, f2, f3.\n");
        goto PLSFin;
    }
    if (ctx->PMNC != 1)
        ctx->PMNC = 0;

    nu = ctx->PMRUN;
    nv = ctx->PMRVN;
    mu = ctx->PMRUM;
    mv = ctx->PMRVM;
    ua = ctx->PMRUA;
    ub = ctx->PMRUB;
    va = ctx->PMRVA;
    vb = ctx->PMRVB;
    if (nu < 2 || nv < 2) {
        printf1(ctx, "Error: nu and nv must be greater than 1.\n");
        goto PLSFin;
    }
    ctx->NNU = (nu - 1) * mu + 1;
    ctx->NNV = (nv - 1) * mv + 1;
    ctx->NNUV = ctx->NNU * ctx->NNV;

    printf1(ctx, "\nRange of u and v parameters.\n");
    printf1(ctx, "u: %16.6f %16.6f\n",ua,ub);
    printf1(ctx, "v: %16.6f %16.6f\n",va,vb);
    printf1(ctx, "\nSize of evalutation grid: %d x %d = %d points.\n",ctx->NNU,ctx->NNV,ctx->NNUV);
    printf1(ctx, "Size of visible grid: %d x %d = %d points.\n\n",nu,nv,nu * nv);

    /* get coordinates at evaluation grid into AcXF, AcYF, AcZF */

    if (alloc_acxf(ctx, ctx->NNUV + 1))
        goto PLSFin;
    if (alloc_acyf(ctx, ctx->NNUV + 1))
        goto PLSFin;
    if (alloc_aczf(ctx, ctx->NNUV + 1))
        goto PLSFin;

    du = (ub - ua) / (ctx->NNU - 1);
    dv = (vb - va) / (ctx->NNV - 1);

    if (plsurf3dt(ctx, ua,va,du,dv,ctx->AcXF,ctx->AcYF,ctx->AcZF))         /* get data */
        goto PLSFin;
           
    printf1(ctx, "\n             u              v        f1(u,v)        f2(u,v)        f3(u,v)\n");
    printf1(ctx, "%14.4f %14.4f %14.4f %14.4f %14.4f\n",ua,va,(double)(ctx->AcXF[0]),(double)(ctx->AcYF[0]),(double)(ctx->AcZF[0]));
    printf1(ctx, "%14.4f %14.4f %14.4f %14.4f %14.4f\n",ua,vb,(double)(ctx->AcXF[ctx->NNV-1]),(double)(ctx->AcYF[ctx->NNV-1]),(double)(ctx->AcZF[ctx->NNV-1]));
    printf1(ctx, "%14.4f %14.4f %14.4f %14.4f %14.4f\n",ub,va,(double)(ctx->AcXF[(ctx->NNU-1)*ctx->NNV]),(double)(ctx->AcYF[(ctx->NNU-1)*ctx->NNV]),(double)(ctx->AcZF[(ctx->NNU-1)*ctx->NNV]));
    printf1(ctx, "%14.4f %14.4f %14.4f %14.4f %14.4f\n\n",ub,vb,(double)(ctx->AcXF[ctx->NNUV-1]),(double)(ctx->AcYF[ctx->NNUV-1]),(double)(ctx->AcZF[ctx->NNUV-1]));

    if (plsurf3p(ctx, nu,nv,mu,mv))                         /* perform the plot */
        goto PLSFin;

    err = 0;

PLSFin:
    p_clean(ctx);
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  plsurf3d        Plot surface of a set of data points.                   */
/*                                                                          */
/*  plsurf3d(                                                               */
/*      alg=...,    algorithm, def. 1                                       */


/*      ru=...,     range of first argument: ua,ub,nu,mu                    */
/*      rv=...,     range of second argument: va,vb,nv,mv                   */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      cont=...,   1 add internal contour, def. 0                          */
/*      gs=...,     grey scale value                                        */
/*      nc=...,     1 if no clipping                                        */
/*  ) = X,Y,Z;      variables providing the data points.                    */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plsurf3d(TDAContext *ctx)
{
    register int i,k;
    int err,nu,mu,nv,mv,r,iu,iv,ncp;
    double ua,ub,va,vb,du,dv,ix,iy,iz,u,v;

    ncp = 3;
    err = -1;
    if (check_pcmd(ctx, 0,3))
        return(-1);

    printf1(ctx, "Surface of a point set. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 8,4,1))     /* get parameters */
        goto PLSDFin;

    if (ctx->PMNV != 3) {
        printf1(ctx, "Error: need exactly three variables on right-hand side.\n");
        goto PLSDFin;
    }
    if (ctx->NOC < 3) {
        printf1(ctx, "Error: need at least three data points.\n");
        goto PLSDFin;
    }
    if (ctx->PMRUFlg == 0 || ctx->PMRVFlg == 0) {
        printf1(ctx, "Error: need ru and rv parameters.\n");
        goto PLSDFin;
    }
    if (ctx->PMALG < 1 || ctx->PMALG > 3)
        ctx->PMALG = 1;

    if (ctx->PMNC != 1)
        ctx->PMNC = 0;

    nu = ctx->PMRUN;
    nv = ctx->PMRVN;
    mu = ctx->PMRUM;
    mv = ctx->PMRVM;
    ua = ctx->PMRUA;
    ub = ctx->PMRUB;
    va = ctx->PMRVA;
    vb = ctx->PMRVB;
    if (nu < 2 || nv < 2) {
        printf1(ctx, "Error: nu and nv must be greater than 1.\n");
        goto PLSDFin;
    }
    ctx->NNU = (nu - 1) * mu + 1;
    ctx->NNV = (nv - 1) * mv + 1;
    ctx->NNUV = ctx->NNU * ctx->NNV;

    printf1(ctx, "Algorithm: %d\n\n",ctx->PMALG);
    printf1(ctx, "Range of u and v parameters.\n");
    printf1(ctx, "u: %16.6f %16.6f\n",ua,ub);
    printf1(ctx, "v: %16.6f %16.6f\n",va,vb);
    printf1(ctx, "\nSize of evalutation grid: %d x %d = %d points.\n",ctx->NNU,ctx->NNV,ctx->NNUV);
    printf1(ctx, "Size of visible grid: %d x %d = %d points.\n\n",nu,nv,nu * nv);

    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];
    iz = ctx->PMVIdx[2];

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto PLSDFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto PLSDFin;
    if (alloc_acz(ctx, ctx->NOC + 1))
        goto PLSDFin;

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcX[i + 1] = get_data(ctx, (int)ix,i);
        ctx->AcY[i + 1] = get_data(ctx, (int)iy,i);
        ctx->AcZ[i + 1] = get_data(ctx, (int)iz,i);
    }
    intp_dbox3(ctx, ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcZ);     /* write bounding box of input data */
    newline(ctx);

    if (alloc_acu(ctx, ctx->NNU + 2))
        goto PLSDFin;
    if (alloc_acv(ctx, ctx->NNV + 2))
        goto PLSDFin;
    if (alloc_acw(ctx, ctx->NNUV + 2))
        goto PLSDFin;

    du = (ub - ua) / (ctx->NNU - 1);
    dv = (vb - va) / (ctx->NNV - 1);

    for (iu = 1; iu <= ctx->NNU; ++iu) {
        u = ua + (double)(iu - 1) * du;
        ctx->AcU[iu] = u;
    }
    for (iv = 1; iv <= ctx->NNV; ++iv) {
        v = va + (double)(iv - 1) * dv;
        ctx->AcV[iv] = v;
    }
    r = intp_grid(ctx, ctx->PMALG,ctx->NOC,ctx->AcX,ctx->AcY,ctx->AcZ,ctx->NNU,ctx->NNV,ctx->AcU,ctx->AcV,ctx->AcW,ncp);
    if (r)
        goto PLSDFin;

    /* get coordinates at evaluation grid into AcXF, AcYF, AcZF */
             
    if (alloc_acxf(ctx, ctx->NNUV + 1))
        goto PLSDFin;
    if (alloc_acyf(ctx, ctx->NNUV + 1))
        goto PLSDFin;
    if (alloc_aczf(ctx, ctx->NNUV + 1))
        goto PLSDFin;

    for (iu = 0; iu < ctx->NNU; ++iu) {
        u = ctx->AcU[iu + 1];
        for (iv = 0; iv < ctx->NNV; ++iv) {
            v = ctx->AcV[iv + 1];
            k = iu * ctx->NNV + iv;
            ctx->AcXF[k] = (float)(u);            
            ctx->AcYF[k] = (float)(v);                
            ctx->AcZF[k] = (float)(ctx->AcW[k + 1]);
        }
    }
    if (plsurf3p(ctx, nu,nv,mu,mv))                  /* perform the plot */
        goto PLSDFin;
                           
    err = 0;

PLSDFin:
    p_clean(ctx);
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  plsurf3p    Plot the surface. This function is called by plsurf3 and    */
/*              plsurf3d.                                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plsurf3p(TDAContext *ctx, int nu,int nv,int mu,int mv)
{
    register int i,j,k;
    int iu,iv,k1,k2,h,h1,h2,nflag; 
    double fxa = 0.0,fxb = 0.0,fya = 0.0,fyb = 0.0,fza = 0.0,fzb = 0.0;
    double xtmax,ytmax,x,y,z,x1,y1,z1,x2,y2,z2,x0,y0,dir;
    double ax,bx,cx,dx,ay,by,cy,dy;

    /* Calculate size of frame buffer. Also substitute AcXF, AcYF and
       AcZF by the coordinates of the projection. */
         
    for (k = 0; k < ctx->NNUV; ++k) {
        ps_3dprj1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],(double)ctx->AcZF[k],&x,&y,&z);
        ctx->AcXF[k] = (float)x;
        ctx->AcYF[k] = (float)y;
        ctx->AcZF[k] = (float)z;
      
        if (k == 0) {
            fxa = fxb = x;
            fya = fyb = y;
            fxa = fxb = z;
        }
        else {
            fxa = dmin(ctx, fxa,x);
            fxb = dmax(ctx, fxb,x);
            fya = dmin(ctx, fya,y);
            fyb = dmax(ctx, fyb,y);
            fza = dmin(ctx, fza,x);
            fzb = dmax(ctx, fzb,z);
        }
    }
    if (alloc_ack(ctx, ctx->NNUV + 1))                /* pointer for sorting */
        return(-1);   
       
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    if (ctx->PMGSFlg) {
        if (ctx->PMGS < 0.0)
            ctx->PMGS = 0.0;
        else if (ctx->PMGS > 1.0)
            ctx->PMGS = 1.0;

        if (ctx->PMGSFlg == 2) {
            if (ctx->PMGS1 < 0.0)
                ctx->PMGS1 = 0.0;
            else if (ctx->PMGS1 > 1.0)
                ctx->PMGS1 = 1.0;
        }
        else
            ctx->PMGS1 = ctx->PMGS;

        if (sortdp2f(ctx, ctx->NNUV,ctx->AcZF,ctx->AcXF,ctx->AcK))
            return(-1);   
        
        fprintf(ctx->PSFd,"gsave\n");

        for (i = 0; i < ctx->NNUV; ++i) {
            k = ctx->AcK[i];
            iu = k / ctx->NNV;
            iv = k % ctx->NNV;
            if (iu >= ctx->NNU - 1 || iv >= ctx->NNV - 1)
                continue;

            k1 = iu * ctx->NNV + iv;
            k2 = (iu + 1) * ctx->NNV + iv + 1;
       
            ax = (double)ctx->AcXF[k1];
            ay = (double)ctx->AcYF[k1];
            bx = (double)ctx->AcXF[k2];
            by = (double)ctx->AcYF[k2];
            cx = (double)ctx->AcXF[k2 - 1];
            cy = (double)ctx->AcYF[k2 - 1];
            dx = (double)ctx->AcXF[k1 + 1];
            dy = (double)ctx->AcYF[k1 + 1];

            ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
            ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
            ps_3dplot1(ctx, cx,cy,1,ctx->PMNC);
            fprintf(ctx->PSFd,"gsave\n");

            dir = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            if (dir >= 0.0)
                ps_fill(ctx, ctx->PMGS);
            else
                ps_fill(ctx, ctx->PMGS1);
            fprintf(ctx->PSFd,"grestore\nnewpath\n");

            ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
            ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
            ps_3dplot1(ctx, dx,dy,1,ctx->PMNC);
            fprintf(ctx->PSFd,"gsave\n");

            dir = (dx - ax) * (by - ay) - (dy - ay) * (bx - ax);
            if (dir >= 0.0)
                ps_fill(ctx, ctx->PMGS);
            else
                ps_fill(ctx, ctx->PMGS1);
            fprintf(ctx->PSFd,"grestore\nnewpath\n");
        }
        fprintf(ctx->PSFd,"grestore\n");
    }

    /* Calculate max lenght of the triangle sides. */
        
    xtmax = ytmax = 0.0;

    for (iu = 0; iu < ctx->NNU - 1; ++iu) {
        for (iv = 0; iv < ctx->NNV - 1; ++iv) {

            k1 = iu * ctx->NNV + iv;
            k2 = (iu + 1) * ctx->NNV + iv + 1;
       
            ax = (double)ctx->AcXF[k1];
            ay = (double)ctx->AcYF[k1];
            bx = (double)ctx->AcXF[k2];
            by = (double)ctx->AcYF[k2];
            cx = (double)ctx->AcXF[k2 - 1];
            cy = (double)ctx->AcYF[k2 - 1];
            dx = (double)ctx->AcXF[k1 + 1];
            dy = (double)ctx->AcYF[k1 + 1];

            xtmax = dmax(ctx, xtmax,fabs(ax - cx));
            xtmax = dmax(ctx, xtmax,fabs(bx - cx));
            xtmax = dmax(ctx, xtmax,fabs(ax - dx));
            xtmax = dmax(ctx, xtmax,fabs(bx - dx));
            xtmax = dmax(ctx, xtmax,fabs(ax - bx));
            xtmax = dmax(ctx, xtmax,fabs(dx - cx));

            ytmax = dmax(ctx, ytmax,fabs(ay - cy));
            ytmax = dmax(ctx, ytmax,fabs(by - cy));
            ytmax = dmax(ctx, ytmax,fabs(ay - dy));
            ytmax = dmax(ctx, ytmax,fabs(by - dy));
            ytmax = dmax(ctx, ytmax,fabs(ay - by));
            ytmax = dmax(ctx, ytmax,fabs(dy - cy));
        }
    }

    /* sort points first wrt AcXF, then AcYF */

    if (sortdp2f(ctx, ctx->NNUV,ctx->AcXF,ctx->AcYF,ctx->AcK))
        return(-1);   

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);

    for (i = 0; i < nu; ++i) {                      /* plot u-lines */
        nflag = 0;
        iu = i * mu;
        k  = iu * ctx->NNV;
        x = (double)ctx->AcXF[k];
        y = (double)ctx->AcYF[k];
        z = (double)ctx->AcZF[k];
        h = plsurf3t(ctx, x,y,z,xtmax,ytmax);

        if (h == 0) {
            ps_3dplot1(ctx, x,y,0,ctx->PMNC);
            nflag = 1;
        }
        for (j = 1; j < ctx->NNV; ++j) {
            x1 = (double)ctx->AcXF[k + j];
            y1 = (double)ctx->AcYF[k + j];
            z1 = (double)ctx->AcZF[k + j];
            h1 = plsurf3t(ctx, x1,y1,z1,xtmax,ytmax);

            if (h == 0 && h1 == 0) {        /* check midpoint */
                x2 = (x + x1) / 2.0;
                y2 = (y + y1) / 2.0;
                z2 = (z + z1) / 2.0;
                h2 = plsurf3t(ctx, x2,y2,z2,xtmax,ytmax);
                if (h2 == 0) {
                    ps_3dplot1(ctx, x1,y1,1,ctx->PMNC);
                }
                else {
                    plsurf3n(ctx, 1,x,y,z,x2,y2,z2,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(ctx, x0,y0,1,ctx->PMNC);
                    fprintf(ctx->PSFd,"stroke\n");
                    nflag = 0;

                    plsurf3n(ctx, 0,x2,y2,z2,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(ctx, x0,y0,0,ctx->PMNC);
                    nflag = 1;
                }
            }
            else if (h == 0 && h1 == 1) {
                plsurf3n(ctx, 1,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(ctx, x0,y0,1,ctx->PMNC);
                fprintf(ctx->PSFd,"stroke\n");
                nflag = 0;
            }
            else if (h != 0 && h1 == 0) {
                plsurf3n(ctx, 0,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(ctx, x0,y0,0,ctx->PMNC);
                nflag = 1;
            }
            x = x1;
            y = y1;
            z = z1;
            h = h1;
        }
        if (nflag) {
            fprintf(ctx->PSFd,"stroke\n");
            nflag = 0;
        }
    }
    for (i = 0; i < nv; ++i) {                      /* plot v-lines */
        nflag = 0;
        iv = i * mv;
        x = (double)ctx->AcXF[iv];
        y = (double)ctx->AcYF[iv];
        z = (double)ctx->AcZF[iv];
        h = plsurf3t(ctx, x,y,z,xtmax,ytmax);

        if (h == 0) {
            ps_3dplot1(ctx, x,y,0,ctx->PMNC);
            nflag = 1;
        }
        for (j = 1; j < ctx->NNU; ++j) {
            k = j * ctx->NNV + iv;
            x1 = (double)ctx->AcXF[k];
            y1 = (double)ctx->AcYF[k];
            z1 = (double)ctx->AcZF[k];
            h1 = plsurf3t(ctx, x1,y1,z1,xtmax,ytmax);
       
            if (h == 0 && h1 == 0) {        /* check midpoint */
                x2 = (x + x1) / 2.0;
                y2 = (y + y1) / 2.0;
                z2 = (z + z1) / 2.0;
                h2 = plsurf3t(ctx, x2,y2,z2,xtmax,ytmax);
                if (h2 == 0) {
                    ps_3dplot1(ctx, x1,y1,1,ctx->PMNC);
                }
                else {
                    plsurf3n(ctx, 1,x,y,z,x2,y2,z2,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(ctx, x0,y0,1,ctx->PMNC);
                    fprintf(ctx->PSFd,"stroke\n");
                    nflag = 0;

                    plsurf3n(ctx, 0,x2,y2,z2,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(ctx, x0,y0,0,ctx->PMNC);
                    nflag = 1;
                }
            }
            else if (h == 0 && h1 == 1) {
                plsurf3n(ctx, 1,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(ctx, x0,y0,1,ctx->PMNC);
                fprintf(ctx->PSFd,"stroke\n");
                nflag = 0;
            }
            else if (h != 0 && h1 == 0) {
                plsurf3n(ctx, 0,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(ctx, x0,y0,0,ctx->PMNC);
                nflag = 1;
            }
            x = x1;
            y = y1;
            z = z1;
            h = h1;
        }
        if (nflag) {
            fprintf(ctx->PSFd,"stroke\n");
            nflag = 0;
        }
    }
    fprintf(ctx->PSFd,"grestore\n");

    if (ctx->PMCONT == 1)                        /* add internal contour line */
        plsurf3ct(ctx, xtmax,ytmax);

    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  plsurf3t(x,y,z)     Return 1 if (x,y,z) is hidden by some triangle.     */
/*                      Otherwise return 0.                                 */

int plsurf3t(TDAContext *ctx, double x,double y,double z,double xtmax,double ytmax)
{
    register int k,kp;
    int iu,iv,k0;

    k0 = plsurf3t1(ctx, x,xtmax);
    k = k0 + 1;
    while (--k >= 0) {
        kp = ctx->AcK[k];
        if (fabs(x - (double)ctx->AcXF[kp]) > xtmax)
            break;                            
        if (fabs(y - (double)ctx->AcYF[kp]) > ytmax)
            continue;

        iu = kp / ctx->NNV;
        iv = kp % ctx->NNV;

        if (iu >= ctx->NNU - 1 || iv >= ctx->NNV - 1)
            continue;
        if (plsurf3t2(ctx, x,y,z,iu,iv))  
            return(1);
    }
    k = k0;
    while (++k < ctx->NNUV) {
        kp = ctx->AcK[k];
        if (fabs(x - (double)ctx->AcXF[kp]) > xtmax)
            break;                            
        if (fabs(y - (double)ctx->AcYF[kp]) > ytmax)
            continue;

        iu = kp / ctx->NNV;
        iv = kp % ctx->NNV;

        if (iu >= ctx->NNU - 1 || iv >= ctx->NNV - 1)
            continue;
        if (plsurf3t2(ctx, x,y,z,iu,iv))  
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3t1(x,y,z)    Return index of nearest neighbor.                   */

int plsurf3t1(TDAContext *ctx, double x,double xtmax)
{
    register int k0,k1,k2;
    double x0,x1,x2;

    k0 = 0;
    k1 = ctx->NNUV - 1;
    x0 = (double)ctx->AcXF[ctx->AcK[k0]];
    if (fabs(x - x0) <= xtmax)
        return(k0);

    x1 = (double)ctx->AcXF[ctx->AcK[k1]];
    if (fabs(x - x1) <= xtmax)
        return(k1);

    while (x0 < x && x < x1) {
        k2 = (k0 + k1) / 2;
        x2 = (double)ctx->AcXF[ctx->AcK[k2]];
        if (x <= x2) {
            k1 = k2;
            x1 = x2;
        }
        else {
            k0 = k2;
            x0 = x2;
        }
        if (k0 >= k1 - 1)
            break;
    }
    if (x <= x0)
        return(k0);
    else
        return(k1);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3t2(x,y,z,iu,iv)   Return 1 if (x,y,z) is hidden by the triangle  */
/*                           defined by iu and iv. Otherwise return 0.      */

int plsurf3t2(TDAContext *ctx, double x,double y,double z,int iu,int iv)
{
    int k1,k2;
    double ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,a,b,c,d,zz;
    double tol = 0.001;

    k1 = iu * ctx->NNV + iv;
    k2 = (iu + 1) * ctx->NNV + iv + 1;
       
    ax = (double)ctx->AcXF[k1];
    ay = (double)ctx->AcYF[k1];
    az = (double)ctx->AcZF[k1];
    bx = (double)ctx->AcXF[k2];
    by = (double)ctx->AcYF[k2];
    bz = (double)ctx->AcZF[k2];
    cx = (double)ctx->AcXF[k2 - 1];
    cy = (double)ctx->AcYF[k2 - 1];
    cz = (double)ctx->AcZF[k2 - 1];
    dx = (double)ctx->AcXF[k1 + 1];
    dy = (double)ctx->AcYF[k1 + 1];
    dz = (double)ctx->AcZF[k1 + 1];

    if (plsurf3c(ctx, ax,ay,bx,by,cx,cy,x,y) == 1) {

        /* get coefficients of plane equation */
        a = (by - ay) * (cz - az) - (bz - az) * (cy - ay); 
        b = (bz - az) * (cx - ax) - (bx - ax) * (cz - az); 
        c = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax); 
        d = -(a * ax + b * ay + c * az);
        if (c == 0.0)
            c = ctx->EPSI1;
        zz = -(a * x + b * y + d) / c;
        if (z < zz - tol)  
            return(1);
    }  
    if (plsurf3c(ctx, ax,ay,bx,by,dx,dy,x,y) == 1) {
        a = (by - ay) * (dz - az) - (bz - az) * (dy - ay); 
        b = (bz - az) * (dx - ax) - (bx - ax) * (dz - az); 
        c = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax); 
        d = -(a * ax + b * ay + c * az);
        if (c == 0.0)
            c = ctx->EPSI1;
        zz = -(a * x + b * y + d) / c;
        if (z < zz - tol)  
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3dt() Get data from the three functions.                          */
/*              Return 0 if OK, -1 if error.                                */

int plsurf3dt(TDAContext *ctx, double ua,double va,double du,double dv,float *x,float *y,float *z)
{
    register int iu,iv,k;
    int fi,n,r,atyp;
    double tmp,f,u,v;

    for (fi = 1; fi <= 3; ++fi) {
        if (fi == 1)  
            r = get_func(ctx, ctx->PMF1,1,&n,0,"f1");
        else if (fi == 2)
            r = get_func(ctx, ctx->PMF2,1,&n,0,"f2");
        else                          
            r = get_func(ctx, ctx->PMF3,1,&n,0,"f3");
        if (r)
            return(-1);  
        if (n) {
            get_func(ctx, NULL,0,&n,0,NULL);
            p_err(ctx, -40,1);
            return(-1);  
        }
       
        /* note: function arguments alphabetically ordered: u, then v */

        if (ctx->FNArgN == 0)                                            
            atyp = 0;
        else if (ctx->FNArgN == 1) {
            if (!strncmp(ctx->FNArgDef[0],"u",(size_t)(ctx->FNArgLen[0])))
                atyp = 1;
            else if (!strncmp(ctx->FNArgDef[0],"v",(size_t)(ctx->FNArgLen[0])))
                atyp = 2;
            else {
                printf1(ctx, "Error: arguments must be u or v.\n");
                return(-1);
            }
        }
        else if (ctx->FNArgN == 2)  
            atyp = 3;
        else {
            printf1(ctx, "Error: at most two function arguments (u,v).\n");
            return(-1);   
        }
        for (iu = 0; iu < ctx->NNU; ++iu) {
            u = ua + (double)iu * du;
            for (iv = 0; iv < ctx->NNV; ++iv) {
                v = va + (double)iv * dv;
                if (atyp == 1)  
                    ctx->FNArgVal[0] = u;
                else if (atyp == 2)  
                    ctx->FNArgVal[0] = v;
                else if (atyp == 3) {
                    ctx->FNArgVal[0] = u;
                    ctx->FNArgVal[1] = v;
                }
                r = get_flval(ctx, &f,ctx->FNArgN,ctx->FNArgVal,0,1,&tmp,&tmp,&tmp);
                if (r) {      /* r from v_eval1(ctx) */
                    printf1(ctx, "Can't evaluate function.\n");
                    prn_emsg2(ctx, r);
                    return(-1);   
                }
                k = iu * ctx->NNV + iv;

                if (fi == 1)  
                    x[k] = (float)f;
                else if (fi == 2)  
                    y[k] = (float)f;
                else                                     
                    z[k] = (float)f;
            }
        }
        get_func(ctx, NULL,0,&n,0,NULL);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3c(ax,ay,bx,by,cx,cy,x,y)                                         */
/*                                                                          */
/*  Return 1 if (x,y) inside triangle, otherwise 0.                         */
/*  Return -1 if all three points on a line.                                */

int plsurf3c(TDAContext *ctx, double ax,double ay,double bx,double by,double cx,double cy, double x,double y)
{
    (void)ctx;        /* unused: the signature is shared */
    double b,c,d;

    x -= ax;
    y -= ay;
    bx -= ax;
    by -= ay;
    cx -= ax;
    cy -= ay;

    d = bx * cy - by * cx;

    if (d == 0.0)
        return(-1);

    b = x * cy - y * cx;
    c = bx * y - by * x;
        
    return(d >  0.0 ? b > 0.0 && c > 0.0 && b + c <  d:
                      b < 0.0 && c < 0.0 && b + c >  d);
}
 
/* ------------------------------------------------------------------------ */
/*  plsurf3n(opt,x,y,z,x1,y1,z1,xtmax,ytmax,x0,y0)                          */
/*                                                                          */
/*  opt = 0: hidden to not hidden.                                          */
/*  opt = 1: not hidden to hidden.                                          */

void plsurf3n(TDAContext *ctx, int opt,double x,double y,double z,double x1,double y1,double z1, double xtmax,double ytmax,double *x0,double *y0)
{
    int i,h;
    double x2,y2,z2;

    for (i = 0; i < 5; ++i) {
        x2 = (x + x1) / 2.0;
        y2 = (y + y1) / 2.0;
        z2 = (z + z1) / 2.0;
        h  = plsurf3t(ctx, x2,y2,z2,xtmax,ytmax);
        if (h == opt) {
            x1 = x2;
            y1 = y2;
            z1 = z2;
        }
        else {
            x = x2;
            y = y2;
            z = z2;
        }
    }
    *x0 = x2;
    *y0 = y2;
}

/* ------------------------------------------------------------------------ */
/*  plsurf3o()  Check faces. Return 0 if all triangles have same positive   */
/*              direction, 1 if all triangles have same negative direction. */
/*              Return -1 if some triangles have different directions.      */

int plsurf3o(TDAContext *ctx, int iu,int iv,int mu,int mv)
{
    register int i,j,k1,k2;
    int first,idir = 0;
    double ax,ay,bx,by,cx,cy,dx,dy,dir;

    iu *= mu;
    iv *= mv;

    first = 1;
    for (i = 0; i < mu; ++i) {
        for (j = 0; j < mv; ++j) {

            k1 = (iu + i) * ctx->NNV + iv + j;
            k2 = (iu + i + 1) * ctx->NNV + iv + j + 1;
       
            ax = (double)ctx->AcXF[k1];
            ay = (double)ctx->AcYF[k1];
            bx = (double)ctx->AcXF[k2];
            by = (double)ctx->AcYF[k2];
            cx = (double)ctx->AcXF[k2 - 1];
            cy = (double)ctx->AcYF[k2 - 1];
            dx = (double)ctx->AcXF[k1 + 1];
            dy = (double)ctx->AcYF[k1 + 1];

            dir = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            if (first) {
                if (dir >= 0.0)
                    idir = 1;
                else
                    idir = 0;
            }
            if (dir >= 0.0 && idir == 0)
                return(-1);
            else if (dir < 0.0 && idir == 1)
                return(-1);

            dir = (dx - ax) * (by - ay) - (dy - ay) * (bx - ax);
            if (dir >= 0.0 && idir == 0)
                return(-1);
            else if (dir < 0.0 && idir == 1)
                return(-1);
                 
        }
    }
    return(idir);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3r(iu,iv)     Fill face beginning at iu,iv                        */

void plsurf3r(TDAContext *ctx, int iu,int iv,int mu,int mv,int *p,double gs)
{
    (void)p;        /* unused: the signature is shared */
    register int i,k;

    iu *= mu;
    iv *= mv;
    k = iu * ctx->NNV + iv;

    fprintf(ctx->PSFd,"gsave\n");
    ps_3dplot1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],0,ctx->PMNC);
    for (i = 1; i <= mv; ++i) {
        k++;
        ps_3dplot1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],1,ctx->PMNC);
    }
    for (i = 1; i <= mu; ++i) {
        k = (iu + i) * ctx->NNV + iv + mv;
        ps_3dplot1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],1,ctx->PMNC);
    }
    for (i = 0; i < mv; ++i) {
        k--;                            
        ps_3dplot1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],1,ctx->PMNC);
    }
    for (i = mu; i > 0; --i) {
        k = (iu + i) * ctx->NNV + iv;
        ps_3dplot1(ctx, (double)ctx->AcXF[k],(double)ctx->AcYF[k],1,ctx->PMNC);
    }
    ps_fill(ctx, gs);
    fprintf(ctx->PSFd,"grestore\n");
}

/* ------------------------------------------------------------------------ */
/*  plsurf3ct()     Add internal contour line.                              */

void plsurf3ct(TDAContext *ctx, double xtmax,double ytmax)
{
    register int k,iu,iv;
    int ra,rb;
    double ax,ay,az,bx,by,bz,cx,cy,dx,dy,dir,dir1;

    fprintf(ctx->PSFd,"\n%%#%d: internal contour\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);

    for (iu = 0; iu < ctx->NNU - 1; ++iu) {
        for (iv = 0; iv < ctx->NNV - 1; ++iv) {
            if (iu == 0 && iv == 0)
                continue;

            k  = iu * ctx->NNV + iv;
            ax = (double)ctx->AcXF[k];
            ay = (double)ctx->AcYF[k];
            az = (double)ctx->AcZF[k];

            if (iv > 0) {
                k  = (iu + 1) * ctx->NNV + iv;
                bx = (double)ctx->AcXF[k];
                by = (double)ctx->AcYF[k];
                bz = (double)ctx->AcZF[k];

                k++;                     
                cx = (double)ctx->AcXF[k];
                cy = (double)ctx->AcYF[k];
                k  = iu * ctx->NNV + iv - 1;
                dx = (double)ctx->AcXF[k];
                dy = (double)ctx->AcYF[k];

                dir  = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
                dir1 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);

                if ((dir >= 0.0 && dir1 >= 0.0) || (dir <= 0.0 && dir1 <= 0.0)) {
                    ra = plsurf3t(ctx, ax,ay,az,xtmax,ytmax);
                    rb = plsurf3t(ctx, bx,by,bz,xtmax,ytmax);             

                    if (ra == 0 && rb == 0) {
                        ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
                        ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    else if (ra == 0 && rb == 1) {
                        plsurf3n(ctx, 1,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
                        ps_3dplot1(ctx, cx,cy,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    else if (ra == 1 && rb == 0) {
                        plsurf3n(ctx, 0,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ctx, cx,cy,0,ctx->PMNC);
                        ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                }
            }
            if (iu > 0) {
                k  = iu * ctx->NNV + iv + 1;
                bx = (double)ctx->AcXF[k];
                by = (double)ctx->AcYF[k];
                bz = (double)ctx->AcZF[k];

                k  = (iu + 1) * ctx->NNV + iv + 1;
                cx = (double)ctx->AcXF[k];
                cy = (double)ctx->AcYF[k];
                k  = (iu - 1) * ctx->NNV + iv;
                dx = (double)ctx->AcXF[k];
                dy = (double)ctx->AcYF[k];

                dir  = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
                dir1 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);

                if ((dir >= 0.0 && dir1 >= 0.0) || (dir <= 0.0 && dir1 <= 0.0)) {

                    ra = plsurf3t(ctx, ax,ay,az,xtmax,ytmax);
                    rb = plsurf3t(ctx, bx,by,bz,xtmax,ytmax);             

                    if (ra == 0 && rb == 0) {
                        ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
                        ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    else if (ra == 0 && rb == 1) {
                        plsurf3n(ctx, 1,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ctx, ax,ay,0,ctx->PMNC);
                        ps_3dplot1(ctx, cx,cy,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                    else if (ra == 1 && rb == 0) {
                        plsurf3n(ctx, 0,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ctx, cx,cy,0,ctx->PMNC);
                        ps_3dplot1(ctx, bx,by,1,ctx->PMNC);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                }
            }
         
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
}


