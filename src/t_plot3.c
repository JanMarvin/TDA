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

/*  functions in t_plot3.c */

double arc_to_degree(double x);
double degree_to_arc(double x);
void ps_3dprj(double x1,double x2,double x3,double *x,double *y);
void ps_3dprj1(double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dprj_inv(double x1,double x2,double x3,double *x,double *y,double *z);
void ps_3dnorm(double *x1,double *x2,double *x3);
double ps_3dist(double x1,double x2,double x3,double y1,double y2,double y3);
void ps_3dgeo(double *x1,double *x2,double *x3);
void ps_3dgeo_inv(double *x1,double *x2,double *x3);
void ps_3dplot(double x1,double x2,double x3,int opt,int geo,int nc);
void ps_3dplot1(double x,double y,int opt,int nc);
int plotp3(int typ);
int plcurv3(void);
int pltext3(void);
int plcirc3(void);
void pl_circ3(double x,double y,double z,double r,double s,double a,double b,
    int hide,double xp,double yp,double zp,int nc);
int plglob3(void);
int pl_lon3(double x,double y,double z,double r,double lon,double lat1,double lat2,int nc);
int pl_lat3(double x,double y,double z,double r,double lat,double lon1,double lon2,int nc);
int pl_meridian(double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_meridian1(double lon,double lata,double latb,double r,int n,int nc,int lt,double lw);
int pl_parallel(double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);
int pl_parallel1(double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw);

int plsurf3(void);
int plsurf3d(void);
int plsurf3p(int nu,int nv,int mu,int mv);
int plsurf3t(double x,double y,double z,double xtmax,double ytmax);
int plsurf3t1(double x,double xtmax);
int plsurf3t2(double x,double y,double z,int iu,int iv);
int plsurf3dt(double ua,double va,double du,double dv,float *x,float *y,float *z);
int plsurf3c(double ax,double ay,double bx,double by,double cx,double cy,double x,double y);
void plsurf3n(int opt,double x,double y,double z,double x1,double y1,double z1,
    double xtmax,double ytmax,double *x0,double *y0);
int plsurf3o(int iu,int iv,int mu,int mv);
void plsurf3ct(double xtmax,double ytmax);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */

#define SHORTMIN -32760
#define SHORTMAX  32760

int NNU = 0;            /* size of evaluation grid: NNU, NNV, NNUV          */
int NNV = 0;
int NNUV = 0;


/* ------------------------------------------------------------------------ */
/*  arc_to_degree(x)                                                        */

double arc_to_degree(double x)
{
    return(180.0 * x / Pi);
}

/* ------------------------------------------------------------------------ */
/*  degree_to_arc(x)                                                        */

double degree_to_arc(double x)
{
    return(x * Pi / 180.0);
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj(x1,x2,x3,x,y)  Return projection in x and y.                   */

void ps_3dprj(double x1,double x2,double x3,double *x,double *y)
{
    *x = PSPR11 * x1 + PSPR12 * x2;
    *y = PSPR21 * x1 + PSPR22 * x2 + PSPR23 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj1(x1,x2,x3,x,y,z)  Return projection in x, y and z.             */

void ps_3dprj1(double x1,double x2,double x3,double *x,double *y,double *z)
{
    *x = PSPR11 * x1 + PSPR12 * x2 + PSPR13 * x3;
    *y = PSPR21 * x1 + PSPR22 * x2 + PSPR23 * x3;
    *z = PSPR31 * x1 + PSPR32 * x2 + PSPR33 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dprj_inv(x1,x2,x3,x,y,z)  Return inverse projection in x, y and z.  */

void ps_3dprj_inv(double x1,double x2,double x3,double *x,double *y,double *z)
{
    *x = PSPI11 * x1 + PSPI12 * x2 + PSPI13 * x3;
    *y = PSPI21 * x1 + PSPI22 * x2 + PSPI23 * x3;
    *z = PSPI31 * x1 + PSPI32 * x2 + PSPI33 * x3;
}

/* ------------------------------------------------------------------------ */
/*  ps_3dnorm(x1,x2,x3)     Return normalized vector.                       */

void ps_3dnorm(double *x1,double *x2,double *x3)
{
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

double ps_3dist(double x1,double x2,double x3,double y1,double y2,double y3)
{
    x1 -= y1;
    x2 -= y2;
    x3 -= y3;
    return(sqrt(x1 * x1 + x2 * x2 + x3 * x3));         
}

/* ------------------------------------------------------------------------ */
/*  ps_3dgeo(x1,x2,x3)  Change from geographic to cartesian coordinates.    */
/*                      Assume geographic coordinates in degrees.           */

void ps_3dgeo(double *x1,double *x2,double *x3)
{
    double p,lambda,beta,cosl,cosb,sinl,sinb;

    if (*x3 < EPSI1) {
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
        printf1("Error in geographical coordinates: %g, %g, %g.\n",*x1,*x2,*x3);
        printf1("Cannot continue.\n");
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

void ps_3dgeo_inv(double *x1,double *x2,double *x3)
{
    double r;

    r = *x1 * *x1 + *x2 * *x2;
    if (r < EPSI1) {
        *x2 = 90.0;
        return;
    }
    *x1 = arc_to_degree(atan2(*x2,*x1));
    *x2 = arc_to_degree(Pi / 2.0 - atan2(sqrt(r),*x3));
    *x3 = sqrt(r + *x3 * *x3);
}

/* ------------------------------------------------------------------------ */
/*  ps_3dplot   PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values.                                       */
/*                                                                          */
/*              if geo != 0 translate from geographic into cartesian        */
/*              coordinates. If nc != 0 update bounding box.                */

void ps_3dplot(double x1,double x2,double x3,int opt,int geo,int nc)
{
    double x,y,px,py;

    if (geo)  
        ps_3dgeo(&x1,&x2,&x3);  
    ps_3dprj(x1,x2,x3,&x,&y);

    px = (x - PA1[0]) / UXLen;
    py = (y - PA1[1]) / UYLen;
    px *= PSXLen;
    py *= PSYLen;

    if (!opt)  
        fprintf(PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(PSFd,"%5.2f %5.2f ",px,py);

    if (nc)
        upd_bbox(1,px,py); 
}

/* ------------------------------------------------------------------------ */
/*  ps_3dplot1  PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values. No projection. The function           */
/*              assumes x,y coordinates in projection space.                */
/*              If nc != 0 update bounding box.                             */

void ps_3dplot1(double x,double y,int opt,int nc)
{
    double px,py;

    px = (x - PA1[0]) / UXLen;
    py = (y - PA1[1]) / UYLen;
    px *= PSXLen;
    py *= PSYLen;

    if (!opt)  
        fprintf(PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(PSFd,"%5.2f %5.2f ",px,py);

    if (nc)
        upd_bbox(1,px,py); 
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

int plotp3(int typ)
{
    int i,n,nn,err,first,ix1,ix2,ix3;
    double x1,x2,x3,x,y,xa,ya;

    err = -1;
    if (check_pcmd(2,3))
        return(-1);

    if (typ == 0) {
        if (parm(CmdBuf + 6,7,1))    /* get parameters */
            goto PLP3Fin;

        n = PMRHSN;
        if (n < 6 || (n / 3) * 3 != n) {
            p_err(-45,1);
            goto PLP3Fin;
        }
    }
    else {
        if (parm(CmdBuf + 5,4,1))    /* get parameters */
            goto PLP3Fin;

        if (NOC < 2)
            goto PLP3Fin;

        n = NOC;
        if (PMNV != 3) {
            p_err(-1,1);
            goto PLP3Fin;
        }
        ix1 = PMVIdx[0];
        ix2 = PMVIdx[1];
        ix3 = PMVIdx[2];
    }

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);
    fprintf(PSFd,"gsave\n");
    if (PMNC == 0)
        set_clip();

    /* first process fill option */

    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        ps_lwidth(0.0);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = PMRHSX[i];
                x2 = PMRHSX[i + 1];
                x3 = PMRHSX[i + 2];
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ix1,i);       
                x2 = get_data(ix2,i);       
                x3 = get_data(ix3,i);       
            }
            if (PMGEO) {
                ps_3dgeo(&x1,&x2,&x3);
            }
            ps_3dprj(x1,x2,x3,&x,&y);
            x = ps_2dx(x);
            y = ps_2dy(y);

            if (PMNC)
                upd_bbox(1,x,y);    

            nn++;
            if (first) {
            /*  fprintf(PSFd,"%5.2f %5.2f m\n",x,0.0);  */
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                if (PMDIR == 1)
                    fprintf(PSFd,"%5.2f %5.2f l\n",x,ya);
                else if (PMDIR == 2)
                    fprintf(PSFd,"%5.2f %5.2f l\n",xa,y);

                fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
            }
            xa = x;
            ya = y;

            i++;   
            if (typ == 0)
                i += 2;
        }
        if (nn) {
            /** fprintf(PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);   **/
            fprintf(PSFd,"gsave\n");
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\nnewpath\n");
        }
    }

    /* now the standard polygon */

    if (PMLW > 0.0)
        ps_lwidth(PMLW);

    if (PMLW > 0.0 && PMLT > 0) {
        ps_ltyp(PMLT);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = PMRHSX[i];
                x2 = PMRHSX[i + 1];
                x3 = PMRHSX[i + 2];
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ix1,i);       
                x2 = get_data(ix2,i);       
                x3 = get_data(ix3,i);       
            }
            if (PMGEO) {
                ps_3dgeo(&x1,&x2,&x3);
            }
            ps_3dprj(x1,x2,x3,&x,&y);
            x = ps_2dx(x);
            y = ps_2dy(y);

            if (PMNC)
                upd_bbox(1,x,y);    

            nn++;
            if (first) {
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
            }
            i++;   
            if (typ == 0)
                i += 2;
            if (i >= n)
                break;
            xa = x;
            ya = y;
        }
        if (nn >= 2 && PMAFlg && PMA1 > 0.0 && PMA2 > 0.0) {       /* plot arrow */
            x -= xa;
            y -= ya;
            if (fabs(y) > EPSI1) {
                fprintf(PSFd,"gsave\ncurrentpoint\nstroke m\n");
                fprintf(PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
                fprintf(PSFd,"%5.2f %5.2f scale\n",PMA1,PMA2);
                fprintf(PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
            }
        }
        if (nn)
            fprintf(PSFd,"stroke\n");
    }
    if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

        i = 0;
        while (i < n) {

            if (typ == 0) {
                x1 = PMRHSX[i];
                x2 = PMRHSX[i + 1];
                x3 = PMRHSX[i + 2];
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x1 = get_data(ix1,i);       
                x2 = get_data(ix2,i);       
                x3 = get_data(ix3,i);       
            }
            if (PMGEO) {
                ps_3dgeo(&x1,&x2,&x3);
            }
            ps_3dprj(x1,x2,x3,&x,&y);
            x = ps_2dx(x);
            y = ps_2dy(y);

            if (PMNC)
                upd_bbox(1,x,y);    

            ps_sym(PMS,x,y,PtMM * PMFS / 2.0);
            i++;   
            if (typ == 0)
                i += 2;
        }
    }
    fprintf(PSFd,"grestore\n");
    err = 0;

PLP3Fin:
    p_clean();
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

int plcurv3(void)
{
    register int i;
    int err,r,f,n,np;      
    double x,y,tmp;
           
    err = -1;
    if (check_pcmd(0,3))
        return(-1);

    if (parm(CmdBuf + 7,1,0))    /* get parameters */
        goto PLCFin;

    if (PMF1A <= 0 || PMF2A <= 0 || PMF3A <= 0 || PMRXFlg == 0) {
        printf1("Need arguments: rx, f1, f2, f3.\n");
        goto PLCFin;
    }
    np = 0;
    x = PMRXA;
    while (x <= PMRXB + PMRXD) {
        np++;
        x += PMRXD;
    }
    if (alloc_acx(np + 1))
        goto PLCFin;
    if (alloc_acy(np + 1))
        goto PLCFin;
    if (alloc_acz(np + 1))
        goto PLCFin;

    for (f = 1; f <= 3; ++f) {
        if (f == 1)
            r = get_func(PMF1,1,&n,0,"f1");
        else if (f == 2)
            r = get_func(PMF2,1,&n,0,"f2");
        else                          
            r = get_func(PMF3,1,&n,0,"f3");
        if (r)
            goto PLCFin;
        if (n) {
            get_func(NULL,0,&n,0,NULL);
            p_err(-40,1);
            goto PLCFin;
        }
        if (FNArgN > 1) {
            printf1("Error: at most one function argument (x).\n");
            goto PLCFin;
        }
        x = PMRXA;
        for (i = 0; i < np; ++i) {
            if (FNArgN == 1)
                FNArgVal[0] = x;    
            r = get_flval(&y,FNArgN,FNArgVal,0,1,&tmp,&tmp,&tmp);
            if (r) {      /* r from v_eval1() */
                printf1("Can't evaluate function.\n");
                prn_emsg2(r);
                goto PLCFin;
            }
            if (f == 1)
                AcX[i] = y;
            else if (f == 2)
                AcY[i] = y;
            else
                AcZ[i] = y;

            if ((x += PMRXD) > PMRXB)
                x = PMRXB;
        }
        get_func(NULL,0,&n,0,NULL);
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC == 0)
        set_clip();

    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    ps_3dplot(AcX[0],AcY[0],AcZ[0],0,0,0);

    for (i = 1; i < np; ++i)  
        ps_3dplot(AcX[i],AcY[i],AcZ[i],1,0,PMNC);
       
    fprintf(PSFd,"stroke\ngrestore\n");
    err = 0;

PLCFin:
    p_clean();
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

int pltext3(void)
{
    int cen,err = -1;
    double x,y,px,py,l;

    if (check_pcmd(0,3))
        return(-1);

    if (parm(CmdBuf + 7,8,1))    /* get parameters */
        goto PLT3Fin;

    if (PMXYZFlg == 0) {
        p_err(-1,1);
        goto PLT3Fin;
    }
    if (PMWF != 0)
        PMWF = 1;

    if (PMR)  
        cen = 0;
    else  
        cen = (int)PMSC;
       
    ps_3dprj(PM3X,PM3Y,PM3Z,&x,&y);

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    px = ps_2dx(x);
    py = ps_2dy(y);

    if (PMS) {
        PMR = cen = 0;
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);
        ps_sym(PMS,px - 1.5 * PtMM * PMFS,py + 0.5 * PtMM * PMFS,PtMM * PMFS / 2.0);
    }
    plot_str(x,y,PMRHSTR,1.5 * PMFS * PtMM,0,cen,PMR,0,PMWF);

    /* update bounding box */

    l = 0.6 * (double)strlen(PMRHSTR) * PMFS * PtMM;
    upd_bbox(1,px,py);     
    upd_bbox(1,px + l,py + l);     

    upd_bbox(1,px,py + 1.3 * PMFS * PtMM);
    upd_bbox(1,px,py - 0.5 * PMFS * PtMM);

    if (PMR > 0)  
        upd_bbox(1,px - 1.3 * PMFS * PtMM,y);

    if (PMR < 0)  
        upd_bbox(1,px + 1.3 * PMFS * PtMM,y);

    err = 0;

PLT3Fin:
    p_clean();
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

int plcirc3(void)
{
    int err;
    double a,b,r,s;
           
    err = -1;
    if (check_pcmd(0,3))
        return(-1);

    if (parm(CmdBuf + 7,7,1))    /* get parameters */
        goto PLCIFin;

    if (PMXYZFlg == 0) {
        p_err(-1,1);
        goto PLCIFin;
    }

    if (PMRHSN == 1) {
        a = 0.0;
        b = 360.0;
    }
    else if (PMRHSN == 3) {
        a = PMRHSX[1];
        b = PMRHSX[2];
    }      
    else {
        p_err(-1,1);
        goto PLCIFin;
    }
    r = PMRHSX[0];
    if (r <= 0.0 || a < -360.0 || a > b || b > 360.0) {
        p_err(-1,1);
        goto PLCIFin;
    }
    if (!PMDVECFlg) {
        p_err(-1,1);
        goto PLCIFin;
    }
    if (PMGEO) {
        ps_3dgeo(&PMDVECX,&PMDVECY,&PMDVECZ);     
    }
    ps_3dnorm(&PMDVECX,&PMDVECY,&PMDVECZ);        
    s = sqrt(PMDVECX * PMDVECX + PMDVECY * PMDVECY);

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);
    pl_circ3(PM3X,PM3Y,PM3Z,r,s,a,b,PMHIDE,PM3X,PM3Y,PM3Z,PMNC);

    err = 0;

PLCIFin:
    p_clean();
    return(err);
}
                                              

/*--------------------------------------------------------------------------*/
/*  pl_circ3(x,y,z,r,s,a,b,hide,xp,yp,zp,nc)                                */
/*                                                                          */
/*  Plot circle with radius r and center x,y,z from a to b. xp,yp,zp is the */   
/*  center of the projection to decide what is hidden. If nc != 0 update    */
/*  bounding box (no clipping).                                             */

void pl_circ3(double x,double y,double z,double r,double s,double a,double b,
    int hide,double xp,double yp,double zp,int nc)
{
    int i,cflag,hflag,first,lt;
    double u,x1,x2,x3,h,hh,cosu,sinu;

    hh = PSPR31 * xp + PSPR32 * yp + PSPR33 * zp;
    lt = PMLT;
    cflag = hflag = first = 0;
    for (i = -360; i < 360; ++i) {
        if ((double)i < a)   
            continue;
        else if ((double)i > b)
            break;
        u = (double)i * Pi / 180.0;
        cosu = cos(u);
        sinu = sin(u);
        if (s > EPSI1) {
            x1 = x - (PMDVECY * cosu + PMDVECX * PMDVECZ * sinu) * r / s;
            x2 = y + (PMDVECX * cosu - PMDVECY * PMDVECZ * sinu) * r / s;
            x3 = z + r * s * sinu;          
        }
        else {
            x1 = x + r * cosu;
            x2 = y + r * sinu;
            x3 = z;
        }
        if (hide < 0) {
            h = PSPR31 * x1 + PSPR32 * x2 + PSPR33 * x3;
            if (h < hh -EPSI1) {
                if (first) {
                    fprintf(PSFd,"stroke\ngrestore\n");
                    hflag = first = 0;
                }
                continue;
            }
        }
        else if (hide > 0) {
            h = PSPR31 * x1 + PSPR32 * x2 + PSPR33 * x3;
            if (h < hh -EPSI1) {
                if (cflag == 0) {
                    if (first) {
                        fprintf(PSFd,"stroke\ngrestore\n");
                        hflag = first = 0;
                    }
                    lt = hide;
                    cflag = 1;
                }
            }
            else if (cflag) {
                fprintf(PSFd,"stroke\ngrestore\n");
                lt = PMLT;
                first = cflag = hflag = 0;
            }
        }
        if (first == 0) {
            fprintf(PSFd,"gsave\n");
            if (nc == 0)
                set_clip();
            ps_ltyp(lt);
            ps_lwidth(PMLW);
            hflag = 1;
        }
        ps_3dplot(x1,x2,x3,first,0,nc);
        first = 1;
    }
    if (hflag)
        fprintf(PSFd,"stroke\ngrestore\n");
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

int pl_meridian(double lon,double lata,double latb,double r,int n,int nc,int lt,double lw)
{
    fprintf(PSFd,"gsave\n");
    if (nc == 0)
        set_clip();
    ps_ltyp(lt);
    ps_lwidth(lw);

      
    /** printf("meridian lon=%g lata=%g latb=%g\n",lon,lata,latb); **/


    if (latb > 90.0) {
        if (pl_meridian1(lon,lata,90.0,r,n,nc,lt,lw))
            return(-1);

        if ((lon += 180.0) > 180)
            lon -= 360.0;

        if (pl_meridian1(lon,180.0 - latb,90.0,r,n,nc,lt,lw))
            return(-1);
    }
    else if (lata < -90.0) {

        if (pl_meridian1(lon,latb,-90.0,r,n,nc,lt,lw))
            return(-1);

        if ((lon += 180.0) > 180)
            lon -= 360.0;

        if (pl_meridian1(lon,180.0 + lata,-90.0,r,n,nc,lt,lw))
            return(-1);
    }
    else {
        if (pl_meridian1(lon,lata,latb,r,n,nc,lt,lw))
            return(-1);
    }
    fprintf(PSFd,"grestore\n");
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

int pl_meridian1(double lon,double lata,double latb,double r,int n,int nc,int lt,double lw)
{
    int i,first;
    double delta,x,y,z,h;               

    delta = (latb - lata) / (double)n;

    first = 0;
    for (i = 0; i <= n; ++i) {
        x = lon;
        y = lata + (double)i * delta;
        z = r;
        ps_3dgeo(&x,&y,&z);  

        h = PSPR31 * x + PSPR32 * y + PSPR33 * z;

        if (h < -EPSI1) {           /* plot only if visible */
            if (first) {
                fprintf(PSFd,"stroke\n");
                first = 0;
            }
        }
        else {
            ps_3dplot(x,y,z,first,0,nc);
            first = 1;
        }
    }
    if (first)
        fprintf(PSFd,"stroke\n");
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

int pl_parallel(double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw)
{
    double delta;                          

    delta = (lonb - lona) / (double)n;

    fprintf(PSFd,"gsave\n");
    if (nc == 0)
        set_clip();
    ps_ltyp(lt);
    ps_lwidth(lw);

    if (lona <= lonb) {
        if (pl_parallel1(lat,lona,lonb,r,n,nc,lt,lw))
            return(-1);
    }
    else {
        if (pl_parallel1(lat,lona,180.0,r,n,nc,lt,lw))
            return(-1);

        if (pl_parallel1(lat,-180.0,lonb,r,n,nc,lt,lw))
            return(-1);
    }
    fprintf(PSFd,"grestore\n");
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_parallel1(lat,lona,lonb,r,n,nc,lt,lw)                                */
/*                                                                          */
/*  plot parallel at latitude lat from longitude lona to lonb. r is radius. */
/*  n is the number of steps, nc is the no clipping option. lt and lw are,  */
/*  respectively, the line type and line width for drawing the lines.       */

int pl_parallel1(double lat,double lona,double lonb,double r,int n,int nc,int lt,double lw)
{
    int i,first;
    double delta,x,y,z,h;               

    delta = (lonb - lona) / (double)n;

    first = 0;
    for (i = 0; i <= n; ++i) {

        x = lona + (double)i * delta;
        y = lat;
        z = r;
        ps_3dgeo(&x,&y,&z);  

        h = PSPR31 * x + PSPR32 * y + PSPR33 * z;

        if (h < -EPSI1) {           /* plot only if visible */
            if (first) {
                fprintf(PSFd,"stroke\n");
                first = 0;
            }
        }
        else {
            ps_3dplot(x,y,z,first,0,nc);
            first = 1;
        }
    }
    if (first)
        fprintf(PSFd,"stroke\n");
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

int plglob3(void)
{
    int err,i;
    double r,s,x;  
           
    err = -1;
    if (check_pcmd(0,3))
        return(-1);

    if (parm(CmdBuf + 7,7,1))    /* get parameters */
        goto PLGLFin;

    if ((r = PMRHSX[0]) <= 0.0) {
        p_err(-1,1);
        goto PLGLFin;
    }

    /* plot contour of sphere */
      
    PMDVECX = PSLon;
    PMDVECY = PSLat;
    PMDVECZ = r;       
    ps_3dgeo(&PMDVECX,&PMDVECY,&PMDVECZ);          
    ps_3dnorm(&PMDVECX,&PMDVECY,&PMDVECZ);        
    s = sqrt(PMDVECX * PMDVECX + PMDVECY * PMDVECY);
                      
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    if (PMCONT)
        pl_circ3(PM3X,PM3Y,PM3Z,r,s,0.0,360.0,0,PM3X,PM3Y,PM3Z,0);
     
    for (i = 0; i < PMNTP; ++i) {
        x = PMTP[i];
        if (x < -180.0 || x > 180.0)
            continue;

        if (pl_lon3(PM3X,PM3Y,PM3Z,r,x,0.0,360.0,0))
            goto PLGLFin;
    }
    for (i = 0; i < PMNTP1; ++i) {
        x = PMTP1[i];
        if (x < -90.0 || x > 90.0)
            continue;

        if (pl_lat3(PM3X,PM3Y,PM3Z,r,x,0.0,360.0,0))
            goto PLGLFin;
    }
    err = 0;

PLGLFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_lon3(x,y,z,r,lon,lat1,lat2,nc)                                       */
/*                                                                          */
/*  Plot circle with radius r and longitude lon from lat1 to lat2.          */
/*  x,y,z is center of globe. if nc != 0 then no clipping.                  */
/*  Return 0 if OK, -1 if error.                                            */

int pl_lon3(double x,double y,double z,double r,double lon,double lat1,double lat2,int nc)
{
    double s;

    if (lon > 90.0)
        lon -= 90.0;
    else
        lon += 90.0;

    PMDVECX = lon;
    PMDVECY = 0.0;
    PMDVECZ = r;
    ps_3dgeo(&PMDVECX,&PMDVECY,&PMDVECZ);    
    ps_3dnorm(&PMDVECX,&PMDVECY,&PMDVECZ);        
    s = sqrt(PMDVECX * PMDVECX + PMDVECY * PMDVECY);
/*  pl_circ3(x,y,z,r,s,lat1,lat2,-1,PM3X,PM3Y,PM3Z,nc);     */
    pl_circ3(x,y,z,r,s,lat1,lat2,-1,x,y,z,nc);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  pl_lat3(x,y,z,r,lat,lon1,lon2,nc)                                       */
/*                                                                          */
/*  Plot circle with radius r and latitude lat (from lon1 to lon2)          */
/*  x,y,z is center of globe. If nc != 0 then no clipping.                  */
/*  Return 0 if OK, -1 if error.                                            */

int pl_lat3(double x,double y,double z,double r,double lat,double lon1,double lon2,int nc)
{
    double s,u,z1;

    u = lat * Pi / 180.0;
    z1 = z + sin(u) * r;
    r *= cos(u);
    PMDVECX = 0.0;    
    PMDVECY = 0.0;   
    PMDVECZ = 1;   
    s = sqrt(PMDVECX * PMDVECX + PMDVECY * PMDVECY);
    pl_circ3(x,y,z1,r,s,lon1,lon2,-1,x,y,z,nc);
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

int plsurf3(void)
{
    int err,nu,mu,nv,mv;
    double ua,ub,va,vb,du,dv;

    err = -1;
    if (check_pcmd(0,3))
        return(-1);

    printf1("Plot parametrically defined surface. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,1,0))    /* get parameters */
        goto PLSFin;

    if (PMF1A <= 0 || PMF2A <= 0 || PMF3A <= 0 || PMRUFlg == 0 || PMRVFlg == 0) {
        printf1("Need arguments: ru, rv, f1, f2, f3.\n");
        goto PLSFin;
    }
    if (PMNC != 1)
        PMNC = 0;

    nu = PMRUN;
    nv = PMRVN;
    mu = PMRUM;
    mv = PMRVM;
    ua = PMRUA;
    ub = PMRUB;
    va = PMRVA;
    vb = PMRVB;
    if (nu < 2 || nv < 2) {
        printf1("Error: nu and nv must be greater than 1.\n");
        goto PLSFin;
    }
    NNU = (nu - 1) * mu + 1;
    NNV = (nv - 1) * mv + 1;
    NNUV = NNU * NNV;

    printf1("\nRange of u and v parameters.\n");
    printf1("u: %16.6f %16.6f\n",ua,ub);
    printf1("v: %16.6f %16.6f\n",va,vb);
    printf1("\nSize of evalutation grid: %d x %d = %d points.\n",NNU,NNV,NNUV);
    printf1("Size of visible grid: %d x %d = %d points.\n\n",nu,nv,nu * nv);

    /* get coordinates at evaluation grid into AcXF, AcYF, AcZF */

    if (alloc_acxf(NNUV + 1))
        goto PLSFin;
    if (alloc_acyf(NNUV + 1))
        goto PLSFin;
    if (alloc_aczf(NNUV + 1))
        goto PLSFin;

    du = (ub - ua) / (NNU - 1);
    dv = (vb - va) / (NNV - 1);

    if (plsurf3dt(ua,va,du,dv,AcXF,AcYF,AcZF))         /* get data */
        goto PLSFin;
           
    printf1("\n             u              v        f1(u,v)        f2(u,v)        f3(u,v)\n");
    printf1("%14.4f %14.4f %14.4f %14.4f %14.4f\n",ua,va,AcXF[0],AcYF[0],AcZF[0]);
    printf1("%14.4f %14.4f %14.4f %14.4f %14.4f\n",ua,vb,AcXF[NNV-1],AcYF[NNV-1],AcZF[NNV-1]);
    printf1("%14.4f %14.4f %14.4f %14.4f %14.4f\n",ub,va,AcXF[(NNU-1)*NNV],AcYF[(NNU-1)*NNV],AcZF[(NNU-1)*NNV]);
    printf1("%14.4f %14.4f %14.4f %14.4f %14.4f\n\n",ub,vb,AcXF[NNUV-1],AcYF[NNUV-1],AcZF[NNUV-1]);

    if (plsurf3p(nu,nv,mu,mv))                         /* perform the plot */
        goto PLSFin;

    err = 0;

PLSFin:
    p_clean();
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

int plsurf3d(void)
{
    register int i,k;
    int err,nu,mu,nv,mv,r,iu,iv,ncp;
    double ua,ub,va,vb,du,dv,ix,iy,iz,u,v;

    ncp = 3;
    err = -1;
    if (check_pcmd(0,3))
        return(-1);

    printf1("Surface of a point set. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 8,4,1))     /* get parameters */
        goto PLSDFin;

    if (PMNV != 3) {
        printf1("Error: need exactly three variables on right-hand side.\n");
        goto PLSDFin;
    }
    if (NOC < 3) {
        printf1("Error: need at least three data points.\n");
        goto PLSDFin;
    }
    if (PMRUFlg == 0 || PMRVFlg == 0) {
        printf1("Error: need ru and rv parameters.\n");
        goto PLSDFin;
    }
    if (PMALG < 1 || PMALG > 3)
        PMALG = 1;

    if (PMNC != 1)
        PMNC = 0;

    nu = PMRUN;
    nv = PMRVN;
    mu = PMRUM;
    mv = PMRVM;
    ua = PMRUA;
    ub = PMRUB;
    va = PMRVA;
    vb = PMRVB;
    if (nu < 2 || nv < 2) {
        printf1("Error: nu and nv must be greater than 1.\n");
        goto PLSDFin;
    }
    NNU = (nu - 1) * mu + 1;
    NNV = (nv - 1) * mv + 1;
    NNUV = NNU * NNV;

    printf1("Algorithm: %d\n\n",PMALG);
    printf1("Range of u and v parameters.\n");
    printf1("u: %16.6f %16.6f\n",ua,ub);
    printf1("v: %16.6f %16.6f\n",va,vb);
    printf1("\nSize of evalutation grid: %d x %d = %d points.\n",NNU,NNV,NNUV);
    printf1("Size of visible grid: %d x %d = %d points.\n\n",nu,nv,nu * nv);

    ix = PMVIdx[0];
    iy = PMVIdx[1];
    iz = PMVIdx[2];

    if (alloc_acx(NOC + 1))
        goto PLSDFin;
    if (alloc_acy(NOC + 1))
        goto PLSDFin;
    if (alloc_acz(NOC + 1))
        goto PLSDFin;

    for (i = 0; i < NOC; ++i) {
        AcX[i + 1] = get_data(ix,i);
        AcY[i + 1] = get_data(iy,i);
        AcZ[i + 1] = get_data(iz,i);
    }
    intp_dbox3(NOC,AcX,AcY,AcZ);     /* write bounding box of input data */
    newline();

    if (alloc_acu(NNU + 2))
        goto PLSDFin;
    if (alloc_acv(NNV + 2))
        goto PLSDFin;
    if (alloc_acw(NNUV + 2))
        goto PLSDFin;

    du = (ub - ua) / (NNU - 1);
    dv = (vb - va) / (NNV - 1);

    for (iu = 1; iu <= NNU; ++iu) {
        u = ua + (double)(iu - 1) * du;
        AcU[iu] = u;
    }
    for (iv = 1; iv <= NNV; ++iv) {
        v = va + (double)(iv - 1) * dv;
        AcV[iv] = v;
    }
    r = intp_grid(PMALG,NOC,AcX,AcY,AcZ,NNU,NNV,AcU,AcV,AcW,ncp);
    if (r)
        goto PLSDFin;

    /* get coordinates at evaluation grid into AcXF, AcYF, AcZF */
             
    if (alloc_acxf(NNUV + 1))
        goto PLSDFin;
    if (alloc_acyf(NNUV + 1))
        goto PLSDFin;
    if (alloc_aczf(NNUV + 1))
        goto PLSDFin;

    for (iu = 0; iu < NNU; ++iu) {
        u = AcU[iu + 1];
        for (iv = 0; iv < NNV; ++iv) {
            v = AcV[iv + 1];
            k = iu * NNV + iv;
            AcXF[k] = u;            
            AcYF[k] = v;                
            AcZF[k] = AcW[k + 1];
        }
    }
    if (plsurf3p(nu,nv,mu,mv))                  /* perform the plot */
        goto PLSDFin;
                           
    err = 0;

PLSDFin:
    p_clean();
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  plsurf3p    Plot the surface. This function is called by plsurf3 and    */
/*              plsurf3d.                                                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int plsurf3p(int nu,int nv,int mu,int mv)
{
    register int i,j,k;
    int iu,iv,k1,k2,h,h1,h2,nflag; 
    double fxa,fxb,fya,fyb,fza,fzb;
    double xtmax,ytmax,x,y,z,x1,y1,z1,x2,y2,z2,x0,y0,dir;
    double ax,bx,cx,dx,ay,by,cy,dy;

    /* Calculate size of frame buffer. Also substitute AcXF, AcYF and
       AcZF by the coordinates of the projection. */
         
    for (k = 0; k < NNUV; ++k) {
        ps_3dprj1((double)AcXF[k],(double)AcYF[k],(double)AcZF[k],&x,&y,&z);
        AcXF[k] = (float)x;
        AcYF[k] = (float)y;
        AcZF[k] = (float)z;
      
        if (k == 0) {
            fxa = fxb = x;
            fya = fyb = y;
            fxa = fxb = z;
        }
        else {
            fxa = dmin(fxa,x);
            fxb = dmax(fxb,x);
            fya = dmin(fya,y);
            fyb = dmax(fyb,y);
            fza = dmin(fza,x);
            fzb = dmax(fzb,z);
        }
    }
    if (alloc_ack(NNUV + 1))                /* pointer for sorting */
        return(-1);   
       
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    if (PMGSFlg) {
        if (PMGS < 0.0)
            PMGS = 0.0;
        else if (PMGS > 1.0)
            PMGS = 1.0;

        if (PMGSFlg == 2) {
            if (PMGS1 < 0.0)
                PMGS1 = 0.0;
            else if (PMGS1 > 1.0)
                PMGS1 = 1.0;
        }
        else
            PMGS1 = PMGS;

        if (sortdp2f(NNUV,AcZF,AcXF,AcK))
            return(-1);   
        
        fprintf(PSFd,"gsave\n");

        for (i = 0; i < NNUV; ++i) {
            k = AcK[i];
            iu = k / NNV;
            iv = k % NNV;
            if (iu >= NNU - 1 || iv >= NNV - 1)
                continue;

            k1 = iu * NNV + iv;
            k2 = (iu + 1) * NNV + iv + 1;
       
            ax = (double)AcXF[k1];
            ay = (double)AcYF[k1];
            bx = (double)AcXF[k2];
            by = (double)AcYF[k2];
            cx = (double)AcXF[k2 - 1];
            cy = (double)AcYF[k2 - 1];
            dx = (double)AcXF[k1 + 1];
            dy = (double)AcYF[k1 + 1];

            ps_3dplot1(ax,ay,0,PMNC);
            ps_3dplot1(bx,by,1,PMNC);
            ps_3dplot1(cx,cy,1,PMNC);
            fprintf(PSFd,"gsave\n");

            dir = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
            if (dir >= 0.0)
                ps_fill(PMGS);
            else
                ps_fill(PMGS1);
            fprintf(PSFd,"grestore\nnewpath\n");

            ps_3dplot1(ax,ay,0,PMNC);
            ps_3dplot1(bx,by,1,PMNC);
            ps_3dplot1(dx,dy,1,PMNC);
            fprintf(PSFd,"gsave\n");

            dir = (dx - ax) * (by - ay) - (dy - ay) * (bx - ax);
            if (dir >= 0.0)
                ps_fill(PMGS);
            else
                ps_fill(PMGS1);
            fprintf(PSFd,"grestore\nnewpath\n");
        }
        fprintf(PSFd,"grestore\n");
    }

    /* Calculate max lenght of the triangle sides. */
        
    xtmax = ytmax = 0.0;

    for (iu = 0; iu < NNU - 1; ++iu) {
        for (iv = 0; iv < NNV - 1; ++iv) {

            k1 = iu * NNV + iv;
            k2 = (iu + 1) * NNV + iv + 1;
       
            ax = (double)AcXF[k1];
            ay = (double)AcYF[k1];
            bx = (double)AcXF[k2];
            by = (double)AcYF[k2];
            cx = (double)AcXF[k2 - 1];
            cy = (double)AcYF[k2 - 1];
            dx = (double)AcXF[k1 + 1];
            dy = (double)AcYF[k1 + 1];

            xtmax = dmax(xtmax,fabs(ax - cx));
            xtmax = dmax(xtmax,fabs(bx - cx));
            xtmax = dmax(xtmax,fabs(ax - dx));
            xtmax = dmax(xtmax,fabs(bx - dx));
            xtmax = dmax(xtmax,fabs(ax - bx));
            xtmax = dmax(xtmax,fabs(dx - cx));

            ytmax = dmax(ytmax,fabs(ay - cy));
            ytmax = dmax(ytmax,fabs(by - cy));
            ytmax = dmax(ytmax,fabs(ay - dy));
            ytmax = dmax(ytmax,fabs(by - dy));
            ytmax = dmax(ytmax,fabs(ay - by));
            ytmax = dmax(ytmax,fabs(dy - cy));
        }
    }

    /* sort points first wrt AcXF, then AcYF */

    if (sortdp2f(NNUV,AcXF,AcYF,AcK))
        return(-1);   

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);

    for (i = 0; i < nu; ++i) {                      /* plot u-lines */
        nflag = 0;
        iu = i * mu;
        k  = iu * NNV;
        x = (double)AcXF[k];
        y = (double)AcYF[k];
        z = (double)AcZF[k];
        h = plsurf3t(x,y,z,xtmax,ytmax);

        if (h == 0) {
            ps_3dplot1(x,y,0,PMNC);
            nflag = 1;
        }
        for (j = 1; j < NNV; ++j) {
            x1 = (double)AcXF[k + j];
            y1 = (double)AcYF[k + j];
            z1 = (double)AcZF[k + j];
            h1 = plsurf3t(x1,y1,z1,xtmax,ytmax);

            if (h == 0 && h1 == 0) {        /* check midpoint */
                x2 = (x + x1) / 2.0;
                y2 = (y + y1) / 2.0;
                z2 = (z + z1) / 2.0;
                h2 = plsurf3t(x2,y2,z2,xtmax,ytmax);
                if (h2 == 0) {
                    ps_3dplot1(x1,y1,1,PMNC);
                }
                else {
                    plsurf3n(1,x,y,z,x2,y2,z2,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(x0,y0,1,PMNC);
                    fprintf(PSFd,"stroke\n");
                    nflag = 0;

                    plsurf3n(0,x2,y2,z2,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(x0,y0,0,PMNC);
                    nflag = 1;
                }
            }
            else if (h == 0 && h1 == 1) {
                plsurf3n(1,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(x0,y0,1,PMNC);
                fprintf(PSFd,"stroke\n");
                nflag = 0;
            }
            else if (h != 0 && h1 == 0) {
                plsurf3n(0,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(x0,y0,0,PMNC);
                nflag = 1;
            }
            x = x1;
            y = y1;
            z = z1;
            h = h1;
        }
        if (nflag) {
            fprintf(PSFd,"stroke\n");
            nflag = 0;
        }
    }
    for (i = 0; i < nv; ++i) {                      /* plot v-lines */
        nflag = 0;
        iv = i * mv;
        x = (double)AcXF[iv];
        y = (double)AcYF[iv];
        z = (double)AcZF[iv];
        h = plsurf3t(x,y,z,xtmax,ytmax);

        if (h == 0) {
            ps_3dplot1(x,y,0,PMNC);
            nflag = 1;
        }
        for (j = 1; j < NNU; ++j) {
            k = j * NNV + iv;
            x1 = (double)AcXF[k];
            y1 = (double)AcYF[k];
            z1 = (double)AcZF[k];
            h1 = plsurf3t(x1,y1,z1,xtmax,ytmax);
       
            if (h == 0 && h1 == 0) {        /* check midpoint */
                x2 = (x + x1) / 2.0;
                y2 = (y + y1) / 2.0;
                z2 = (z + z1) / 2.0;
                h2 = plsurf3t(x2,y2,z2,xtmax,ytmax);
                if (h2 == 0) {
                    ps_3dplot1(x1,y1,1,PMNC);
                }
                else {
                    plsurf3n(1,x,y,z,x2,y2,z2,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(x0,y0,1,PMNC);
                    fprintf(PSFd,"stroke\n");
                    nflag = 0;

                    plsurf3n(0,x2,y2,z2,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                    ps_3dplot1(x0,y0,0,PMNC);
                    nflag = 1;
                }
            }
            else if (h == 0 && h1 == 1) {
                plsurf3n(1,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(x0,y0,1,PMNC);
                fprintf(PSFd,"stroke\n");
                nflag = 0;
            }
            else if (h != 0 && h1 == 0) {
                plsurf3n(0,x,y,z,x1,y1,z1,xtmax,ytmax,&x0,&y0);
                ps_3dplot1(x0,y0,0,PMNC);
                nflag = 1;
            }
            x = x1;
            y = y1;
            z = z1;
            h = h1;
        }
        if (nflag) {
            fprintf(PSFd,"stroke\n");
            nflag = 0;
        }
    }
    fprintf(PSFd,"grestore\n");

    if (PMCONT == 1)                        /* add internal contour line */
        plsurf3ct(xtmax,ytmax);

    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  plsurf3t(x,y,z)     Return 1 if (x,y,z) is hidden by some triangle.     */
/*                      Otherwise return 0.                                 */

int plsurf3t(double x,double y,double z,double xtmax,double ytmax)
{
    register int k,kp;
    int iu,iv,k0;

    k0 = plsurf3t1(x,xtmax);
    k = k0 + 1;
    while (--k >= 0) {
        kp = AcK[k];
        if (fabs(x - (double)AcXF[kp]) > xtmax)
            break;                            
        if (fabs(y - (double)AcYF[kp]) > ytmax)
            continue;

        iu = kp / NNV;
        iv = kp % NNV;

        if (iu >= NNU - 1 || iv >= NNV - 1)
            continue;
        if (plsurf3t2(x,y,z,iu,iv))  
            return(1);
    }
    k = k0;
    while (++k < NNUV) {
        kp = AcK[k];
        if (fabs(x - (double)AcXF[kp]) > xtmax)
            break;                            
        if (fabs(y - (double)AcYF[kp]) > ytmax)
            continue;

        iu = kp / NNV;
        iv = kp % NNV;

        if (iu >= NNU - 1 || iv >= NNV - 1)
            continue;
        if (plsurf3t2(x,y,z,iu,iv))  
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3t1(x,y,z)    Return index of nearest neighbor.                   */

int plsurf3t1(double x,double xtmax)
{
    register int k0,k1,k2;
    double x0,x1,x2;

    k0 = 0;
    k1 = NNUV - 1;
    x0 = (double)AcXF[AcK[k0]];
    if (fabs(x - x0) <= xtmax)
        return(k0);

    x1 = (double)AcXF[AcK[k1]];
    if (fabs(x - x1) <= xtmax)
        return(k1);

    while (x0 < x && x < x1) {
        k2 = (k0 + k1) / 2;
        x2 = (double)AcXF[AcK[k2]];
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

int plsurf3t2(double x,double y,double z,int iu,int iv)
{
    int k1,k2;
    double ax,ay,az,bx,by,bz,cx,cy,cz,dx,dy,dz,a,b,c,d,zz;
    double tol = 0.001;

    k1 = iu * NNV + iv;
    k2 = (iu + 1) * NNV + iv + 1;
       
    ax = (double)AcXF[k1];
    ay = (double)AcYF[k1];
    az = (double)AcZF[k1];
    bx = (double)AcXF[k2];
    by = (double)AcYF[k2];
    bz = (double)AcZF[k2];
    cx = (double)AcXF[k2 - 1];
    cy = (double)AcYF[k2 - 1];
    cz = (double)AcZF[k2 - 1];
    dx = (double)AcXF[k1 + 1];
    dy = (double)AcYF[k1 + 1];
    dz = (double)AcZF[k1 + 1];

    if (plsurf3c(ax,ay,bx,by,cx,cy,x,y) == 1) {

        /* get coefficients of plane equation */
        a = (by - ay) * (cz - az) - (bz - az) * (cy - ay); 
        b = (bz - az) * (cx - ax) - (bx - ax) * (cz - az); 
        c = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax); 
        d = -(a * ax + b * ay + c * az);
        if (c == 0.0)
            c = EPSI1;
        zz = -(a * x + b * y + d) / c;
        if (z < zz - tol)  
            return(1);
    }  
    if (plsurf3c(ax,ay,bx,by,dx,dy,x,y) == 1) {
        a = (by - ay) * (dz - az) - (bz - az) * (dy - ay); 
        b = (bz - az) * (dx - ax) - (bx - ax) * (dz - az); 
        c = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax); 
        d = -(a * ax + b * ay + c * az);
        if (c == 0.0)
            c = EPSI1;
        zz = -(a * x + b * y + d) / c;
        if (z < zz - tol)  
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3dt() Get data from the three functions.                          */
/*              Return 0 if OK, -1 if error.                                */

int plsurf3dt(double ua,double va,double du,double dv,float *x,float *y,float *z)
{
    register int iu,iv,k;
    int fi,n,r,atyp;
    double tmp,f,u,v;

    for (fi = 1; fi <= 3; ++fi) {
        if (fi == 1)  
            r = get_func(PMF1,1,&n,0,"f1");
        else if (fi == 2)
            r = get_func(PMF2,1,&n,0,"f2");
        else                          
            r = get_func(PMF3,1,&n,0,"f3");
        if (r)
            return(-1);  
        if (n) {
            get_func(NULL,0,&n,0,NULL);
            p_err(-40,1);
            return(-1);  
        }
       
        /* note: function arguments alphabetically ordered: u, then v */

        if (FNArgN == 0)                                            
            atyp = 0;
        else if (FNArgN == 1) {
            if (!strncmp(FNArgDef[0],"u",FNArgLen[0]))
                atyp = 1;
            else if (!strncmp(FNArgDef[0],"v",FNArgLen[0]))
                atyp = 2;
            else {
                printf1("Error: arguments must be u or v.\n");
                return(-1);
            }
        }
        else if (FNArgN == 2)  
            atyp = 3;
        else {
            printf1("Error: at most two function arguments (u,v).\n");
            return(-1);   
        }
        for (iu = 0; iu < NNU; ++iu) {
            u = ua + (double)iu * du;
            for (iv = 0; iv < NNV; ++iv) {
                v = va + (double)iv * dv;
                if (atyp == 1)  
                    FNArgVal[0] = u;
                else if (atyp == 2)  
                    FNArgVal[0] = v;
                else if (atyp == 3) {
                    FNArgVal[0] = u;
                    FNArgVal[1] = v;
                }
                r = get_flval(&f,FNArgN,FNArgVal,0,1,&tmp,&tmp,&tmp);
                if (r) {      /* r from v_eval1() */
                    printf1("Can't evaluate function.\n");
                    prn_emsg2(r);
                    return(-1);   
                }
                k = iu * NNV + iv;

                if (fi == 1)  
                    x[k] = (float)f;
                else if (fi == 2)  
                    y[k] = (float)f;
                else                                     
                    z[k] = (float)f;
            }
        }
        get_func(NULL,0,&n,0,NULL);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  plsurf3c(ax,ay,bx,by,cx,cy,x,y)                                         */
/*                                                                          */
/*  Return 1 if (x,y) inside triangle, otherwise 0.                         */
/*  Return -1 if all three points on a line.                                */

int plsurf3c(double ax,double ay,double bx,double by,double cx,double cy,
    double x,double y)
{
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

void plsurf3n(int opt,double x,double y,double z,double x1,double y1,double z1,
    double xtmax,double ytmax,double *x0,double *y0)
{
    int i,h;
    double x2,y2,z2;

    for (i = 0; i < 5; ++i) {
        x2 = (x + x1) / 2.0;
        y2 = (y + y1) / 2.0;
        z2 = (z + z1) / 2.0;
        h  = plsurf3t(x2,y2,z2,xtmax,ytmax);
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

int plsurf3o(int iu,int iv,int mu,int mv)
{
    register int i,j,k1,k2;
    int first,idir;
    double ax,ay,bx,by,cx,cy,dx,dy,dir;

    iu *= mu;
    iv *= mv;

    first = 1;
    for (i = 0; i < mu; ++i) {
        for (j = 0; j < mv; ++j) {

            k1 = (iu + i) * NNV + iv + j;
            k2 = (iu + i + 1) * NNV + iv + j + 1;
       
            ax = (double)AcXF[k1];
            ay = (double)AcYF[k1];
            bx = (double)AcXF[k2];
            by = (double)AcYF[k2];
            cx = (double)AcXF[k2 - 1];
            cy = (double)AcYF[k2 - 1];
            dx = (double)AcXF[k1 + 1];
            dy = (double)AcYF[k1 + 1];

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

void plsurf3r(int iu,int iv,int mu,int mv,int *p,double gs)
{
    register int i,k;

    iu *= mu;
    iv *= mv;
    k = iu * NNV + iv;

    fprintf(PSFd,"gsave\n");
    ps_3dplot1((double)AcXF[k],(double)AcYF[k],0,PMNC);
    for (i = 1; i <= mv; ++i) {
        k++;
        ps_3dplot1((double)AcXF[k],(double)AcYF[k],1,PMNC);
    }
    for (i = 1; i <= mu; ++i) {
        k = (iu + i) * NNV + iv + mv;
        ps_3dplot1((double)AcXF[k],(double)AcYF[k],1,PMNC);
    }
    for (i = 0; i < mv; ++i) {
        k--;                            
        ps_3dplot1((double)AcXF[k],(double)AcYF[k],1,PMNC);
    }
    for (i = mu; i > 0; --i) {
        k = (iu + i) * NNV + iv;
        ps_3dplot1((double)AcXF[k],(double)AcYF[k],1,PMNC);
    }
    ps_fill(gs);
    fprintf(PSFd,"grestore\n");
}

/* ------------------------------------------------------------------------ */
/*  plsurf3ct()     Add internal contour line.                              */

void plsurf3ct(double xtmax,double ytmax)
{
    register int k,iu,iv;
    int ra,rb;
    double ax,ay,az,bx,by,bz,cx,cy,dx,dy,dir,dir1;

    fprintf(PSFd,"\n%%#%d: internal contour\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);

    for (iu = 0; iu < NNU - 1; ++iu) {
        for (iv = 0; iv < NNV - 1; ++iv) {
            if (iu == 0 && iv == 0)
                continue;

            k  = iu * NNV + iv;
            ax = (double)AcXF[k];
            ay = (double)AcYF[k];
            az = (double)AcZF[k];

            if (iv > 0) {
                k  = (iu + 1) * NNV + iv;
                bx = (double)AcXF[k];
                by = (double)AcYF[k];
                bz = (double)AcZF[k];

                k++;                     
                cx = (double)AcXF[k];
                cy = (double)AcYF[k];
                k  = iu * NNV + iv - 1;
                dx = (double)AcXF[k];
                dy = (double)AcYF[k];

                dir  = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
                dir1 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);

                if ((dir >= 0.0 && dir1 >= 0.0) || (dir <= 0.0 && dir1 <= 0.0)) {
                    ra = plsurf3t(ax,ay,az,xtmax,ytmax);
                    rb = plsurf3t(bx,by,bz,xtmax,ytmax);             

                    if (ra == 0 && rb == 0) {
                        ps_3dplot1(ax,ay,0,PMNC);
                        ps_3dplot1(bx,by,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                    else if (ra == 0 && rb == 1) {
                        plsurf3n(1,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ax,ay,0,PMNC);
                        ps_3dplot1(cx,cy,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                    else if (ra == 1 && rb == 0) {
                        plsurf3n(0,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(cx,cy,0,PMNC);
                        ps_3dplot1(bx,by,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                }
            }
            if (iu > 0) {
                k  = iu * NNV + iv + 1;
                bx = (double)AcXF[k];
                by = (double)AcYF[k];
                bz = (double)AcZF[k];

                k  = (iu + 1) * NNV + iv + 1;
                cx = (double)AcXF[k];
                cy = (double)AcYF[k];
                k  = (iu - 1) * NNV + iv;
                dx = (double)AcXF[k];
                dy = (double)AcYF[k];

                dir  = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
                dir1 = (bx - ax) * (dy - ay) - (by - ay) * (dx - ax);

                if ((dir >= 0.0 && dir1 >= 0.0) || (dir <= 0.0 && dir1 <= 0.0)) {

                    ra = plsurf3t(ax,ay,az,xtmax,ytmax);
                    rb = plsurf3t(bx,by,bz,xtmax,ytmax);             

                    if (ra == 0 && rb == 0) {
                        ps_3dplot1(ax,ay,0,PMNC);
                        ps_3dplot1(bx,by,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                    else if (ra == 0 && rb == 1) {
                        plsurf3n(1,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(ax,ay,0,PMNC);
                        ps_3dplot1(cx,cy,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                    else if (ra == 1 && rb == 0) {
                        plsurf3n(0,ax,ay,az,bx,by,bz,xtmax,ytmax,&cx,&cy);
                        ps_3dplot1(cx,cy,0,PMNC);
                        ps_3dplot1(bx,by,1,PMNC);
                        fprintf(PSFd,"stroke\n");
                    }
                }
            }
         
        }
    }
    fprintf(PSFd,"grestore\n");
}


