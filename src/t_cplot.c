/****************************************************************************/
/*  t_cplot                                                                 */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_gf.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_mat.h"

/*  functions in t_cplot.c */

int pl_plotc(void);
void cont_plot(int nlev,double *flev,int opt);
void cont_draw(double x, double y, int n, int level,int nlev,double *flev,
    int opt);
int pl_plotcm(void);
int pl_plotr(void);
void pl_plotr_c(int im,int row,int col,double lev,double dx,double dy);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

#define PCMAX 1.e6; /* max function value in contour plots                  */

int PCNX = 0;       /* number of x axis intervals                           */
int PCNY = 0;       /* number of y axis intervals                           */
double PCDX = 0.0;  /* width of x axis intervals                            */
double PCDY = 0.0;  /* width of y axis intervals                            */

float *PCFV;        /* array with function values                           */
char *PCBitM;       /* array used in cont_plot()                            */
unsigned
char PBMsk[8] = {        /* Bitmasks for bit-wise stored values in PCBitM   */
    0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01
};

/*--------------------------------------------------------------------------*/
/*  pl_plotc()          Contour plot.                                       */
/*                      CmdBuf: plotc(n,lt,lw,nc,x) = function;             */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotc(void)
{
    register int i,j;
    int nn,err,r,pcbita,pcfva;
    double tmp;
           
    pcbita = pcfva = 0;
    err = -1;

    if (check_pcmd(2,2))
        return(-1);

    printf1("Contour plot.\n");

    if (parm(CmdBuf + 5,9,1))    /* get parameters */
        goto PLCFin;

    if (FNFlg == 0 || PMNTP < 1) {
        p_err(-1,1);
        goto PLCFin;
    }
    if (FNArgN != 2) {
        printf1("Error: there must be exactly two function arguments.\n");
        goto PLCFin;
    }
    if (PMFTYP5) {
        p_err(-40,1);
        goto PLCFin;
    }
    if (FVFlg)
        printf1("Sum over %d data matrix cases.\n",NOC);

    NINTMUsed = -1;

    if (PMNNFlg == 0 || PMNN1 < 1 || PMNN2 < 1)
        PMNN1 = PMNN2 = 10;
      
    PCDX = UXLen / (double) PMNN1;
    PCDY = UYLen / (double) PMNN2;
    PCNX = PMNN1 + 1;
    PCNY = PMNN2 + 1;

    if (!(PCFV = (float *) calloc(PCNX * PCNY + 1,sizeof(float)))) {
        p_err(-2,1);
        goto PLCFin;
    }
    memrq(PCNX * PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * PCNX * PCNY * PMNTP / 8 + 1;
    if (!(PCBitM = (char *)calloc(nn,sizeof(char)))) {
        p_err(-2,1);
        goto PLCFin;
    }
    memrq(nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= PCNX; ++i) {
        for (j = 1; j <= PCNY; ++j) {

            FNArgVal[0] = PA1[0] + (double)(i - 1) * PCDX;
            FNArgVal[1] = PA1[1] + (double)(j - 1) * PCDY;

            r = get_flval(&tmp,FNArgN,FNArgVal,0,0,&tmp,&tmp,&tmp);
            if (r) {      /* r from v_eval1() */
                printf1("Can't evaluate function.\n");
                prn_emsg2(r);
                goto PLCFin;
            }
            PCFV[(i - 1) * PCNY + j] = (float)tmp;
        }
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    cont_plot(PMNTP,PMTP,0);
    err = 0;

PLCFin:
    if (pcfva) {
        free((char *)PCFV);
        memrq(-PCNX * PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)PCBitM);
        memrq(-nn,sizeof(char));
    }
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cont_plot   Contour plot routine.                                       */
/*                                                                          */
/*    The code is adopted from                                              */
/*    Snyder, ACM algorithm 531, Trans. in Math. Software 4(3), 290 - 294   */

void cont_plot(int nlev,double *flev,int opt)
{
    register int k,l,ii,jj,q;
    int ix,icv,ni,ks,jump;
    int idir,nxidir,ibkey,iflag,icur,jcur,iedge;
    int ij[3],l1[5],l2[5],i1[3],i2[3],i3[7];
    double tmp,cval,dmax,z1,z2,zz,xint[5],xy[3];

    ks = 0;
    l1[1] = PCNX;
    l1[2] = PCNY;
    l1[3] = l1[4] = -1;
    i1[1] = i2[1] =  1;
    i1[2] = 0;
    i2[2] = -1;
    i3[1] = i3[4] = i3[5] = 1;
    i3[2] = i3[3] = i3[6] = 0;

    dmax = PCMAX;

    /* set current pen position, default corresponds to z(1,1). */
    
    xy[1] = xy[2] = 1.0;
    
    icur = (int)xy[1];
    if (icur > PCNX) icur = PCNX;
    if (icur < 1)    icur = 1;
    jcur = (int)xy[2];
    if (jcur > PCNY) jcur = PCNY;
    if (jcur < 1)    jcur = 1;
    ibkey = 0;

PCL10:
    ij[1] = icur;
    ij[2] = jcur;
PCL20:
    l2[1] = ij[1];
    l2[3] = -ij[1];
    l2[2] = ij[2];
    l2[4] = -ij[2];
    idir = 0;

PCL30:
    nxidir = idir + 1;
    k = nxidir;
    if (nxidir > 3)
        nxidir = 0;
PCL40:
    if (ij[1] < 0) ij[1] = -ij[1];
    if (ij[2] < 0) ij[2] = -ij[2];

    if (PCFV[(ij[1] - 1) * PCNY + ij[2]] > dmax)
        goto PCL140;

    l = 1;
    while (ij[l] >= l1[l]) {
        if (++l > 2)
            goto PCL140;
    }
    while (1) {
        ii = ij[1] + i1[l];
        jj = ij[2] + i1[3 - l];

        while (PCFV[(ii - 1) * PCNY + jj] > dmax) {
            if (++l > 2)
                goto PCL140;
            else {
                while (ij[l] >= l1[l]) {
                    if (++l > 2)
                        goto PCL140;
                }
                ii = ij[1] + i1[l];
                jj = ij[2] + i1[3 - l];
            }
        }
        jump = 1;

PCL60:
        ix = 1;
        if (ij[3 - l] != 1) {
            ii = ij[1] - i1[3 - l];
            jj = ij[2] - i1[l];
    
            if (PCFV[(ii - 1) * PCNY + jj] <= dmax) {
                ii = ij[1] + i2[l];
                jj = ij[2] + i2[3 - l];

                if (PCFV[(ii - 1) * PCNY + jj] < dmax)
                    ix = 0;
            }
            if (ij[3 - l] >= l1[3 - l])
                goto PCL90;
        }
        ii = ij[1] + i1[3 - l];
        jj = ij[2] + i1[l];

        if (PCFV[(ii - 1) * PCNY + jj] <= dmax &&
                                     PCFV[ij[1] * PCNY + ij[2] + 1] < dmax) {
            if (jump)
                goto PCL100;
            else {
                if (ix != 0)
                    iflag = 4;
                iedge = ks + 2;
                if (iedge > 4)
                    iedge -= 4;
                xint[iedge] = xint[ks];
                goto PCL200;
            }
        }
PCL90:
        ix += 2;
        if (!jump) {
            if (ix != 0)
                iflag = 4;
            iedge = ks + 2;
            if (iedge > 4)
                iedge -= 4;
            xint[iedge] = xint[ks];
            goto PCL200;
        }

PCL100:
        if (ix != 3 && (ix + ibkey) != 0) {

            ii = ij[1] + i1[l];
            jj = ij[2] + i1[3 - l];
            z1 = PCFV[(ij[1] - 1) * PCNY + ij[2]];
            z2 = PCFV[(ii - 1) * PCNY + jj];

            for (icv = 1; icv <= nlev; ++icv) {

                q = 2 * (PCNX * (PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
                if (!(*(PCBitM + q / 8) & PBMsk[q % 8])) {
                    tmp = z1;
                    if (tmp > z2)
                        tmp = z2;
                    if (flev[icv - 1] > tmp) {
                        tmp = z1;
                        if (tmp < z2)
                            tmp = z2;
                        if (flev[icv - 1] <= tmp) {
                            iedge = l;
                            cval = flev[icv - 1];
                            if (ix != 1)
                                iedge += 2;
                            iflag = 2 + ibkey;
                            xint[iedge] = (cval - z1) / (z2 - z1);
                            goto PCL200;
                        }
                    }
                    q = 2 * (PCNX * (PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
                    *(PCBitM + q / 8) |= PBMsk[q % 8];
                }
            }
        }
        if (++l > 2)  
           break;

        while (ij[l] >= l1[l]) {
            if (++l > 2)
                goto PCL140;
        }
    }

PCL140:
    l = (idir % 2) + 1;
    if (l1[k] < 0)
        ij[l] = -ij[l];

    while (1) {
        if (ij[l] < l1[k]) {
            ij[l] += 1;
            if (ij[l] <= l2[k])
                goto PCL40;

            l2[k] = ij[l];
            idir = nxidir;
            goto PCL30;
        }
        if (idir != nxidir) {
            nxidir++;
            ij[l] = l1[k];
            k = nxidir;
            l = 3 - l;
            ij[l] = l2[k];
            if (nxidir > 3)
                nxidir = 0;
        }
        else if (ibkey != 0)  
            return;
        else 
            break;
    }
    ibkey = 1;
    goto PCL10;

PCL200:
    while (1) {
    
        xy[l] = (double)ij[l] + xint[iedge];
        xy[3 - l] = (double) ij[3 - l];
        q = 2 * (PCNX * (PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
        *(PCBitM + q / 8) |= PBMsk[q % 8];

        cont_draw(xy[1],xy[2],iflag,icv,nlev,flev,opt);
    
        if (iflag >= 4) {
            icur = ij[1];
            jcur = ij[2];
            goto PCL20;
        }
        ni = 1;
        if (iedge >= 3) {
            ij[1] -= i3[iedge];
            ij[2] -= i3[iedge + 2];
        }
        for (k = 1; k <= 4; ++k) {
    
            if (k != iedge) {
                ii = ij[1] + i3[k];
                jj = ij[2] + i3[k + 1];
                z1 = PCFV[(ii - 1) * PCNY + jj];
                ii = ij[1] + i3[k + 1];
                jj = ij[2] + i3[k + 2];
                z2 = PCFV[(ii - 1) * PCNY + jj];
    
                tmp = z1;
                if (tmp > z2)
                    tmp = z2;
                if (cval > tmp) {
                    tmp = z1;
                    if (tmp < z2)
                        tmp = z2;
                    if (cval <= tmp) {
                        if (k == 1 || k == 4) {
                            zz = z1;
                            z1 = z2;
                            z2 = zz;
                        }
                        xint[k] = (cval - z1) / (z2 - z1);
                        ni++;
                        ks = k;
                    }
                }
            }
        }
        if (ni != 2) {
            ks = 5 - iedge;
            if (xint[3] >= xint[1]) {
                ks = 3 - iedge;
                if (ks <= 0)
                    ks += 4;
            }
        }
        l = ks;
        iflag = 1;
        jump = 0;
        if (ks >= 3) {
            ij[1] += i3[ks];
            ij[2] += i3[ks + 2];
            l = ks - 2;
        }
        q = 2 * (PCNX * (PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
        if (!(*(PCBitM + q / 8) & PBMsk[q % 8]))  
            goto PCL60;
    
        iflag = 5;
        iedge = ks + 2;
        if (iedge > 4)
            iedge -= 4;
    
        xint[iedge] = xint[ks];
    }
}   

/* ------------------------------------------------------------------------ */
/*  cont_draw   Draw function to be used with cont_plot.                    */

void cont_draw(double x, double y, int n, int level,int nlev,double *flev,
    int opt)
{
    register int i,j,fflag;
    double g,xx,yy,px,py;    
    static double xf,yf;
    static double eps = 1.e-3;

    i = (int)x;
    j = (int)y;
    xx = PA1[0] + (double)(i - 1) * PCDX;
    yy = PA1[1] + (double)(j - 1) * PCDY;
    xx += (x - (double)i) * PCDX;
    yy += (y - (double)j) * PCDY;

    if (opt) {
        xx += PCDX / 2.0;
        yy += PCDY / 2.0;
    }


    px = PtMM * (xx - PA1[0]) * PXLen / UXLen;
    py = PtMM * (yy - PA1[1]) * PYLen / UYLen;

    switch (n) {
        case  1:  fprintf(PSFd,"%4.2f %4.2f l\n",px,py);
                  break;
        case  2:
        case  3:  fprintf(PSFd,"gsave %% level(%d) = %g\n",level,flev[level - 1]);

                  ps_lwidth(PMLW);
                  ps_ltyp(PMLT);

                  fprintf(PSFd,"%4.2f %4.2f m\n",px,py);
                  xf = px;
                  yf = py; 
                  break;
        case  4:
        case  5:  fprintf(PSFd,"%4.2f %4.2f l\n",px,py);

                  /* fill only curves which can be closed */

                  if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {

                      fflag = 0;
                      if (fabs(px - xf) <= eps || fabs(py - yf) <= eps)  
                          fflag = 1;
                      else if ((fabs(px) <= eps && fabs(yf) <= eps) ||
                          (fabs(py) <= eps && fabs(xf) <= eps)) {
                          fprintf(PSFd,"0 0 l %% aaa\n");
                          fflag = 1;
                      }
                      else if ((fabs(px) <= eps && fabs(yf - PSYLen) <= eps) ||
                               (fabs(py - PSYLen) <= eps && fabs(xf) <= eps)) {
                          fprintf(PSFd,"0 %5.2f l %% bbb\n",PSYLen);
                          fflag = 1;
                      }
                      else if ((fabs(px - PSXLen) <= eps && fabs(yf) <= eps) ||
                               (fabs(py) <= eps && fabs(xf - PSXLen) <= eps)) {
                          fprintf(PSFd,"%5.2f 0 l %% ccc\n",PSXLen);
                          fflag = 1;
                      }
                      else if ((fabs(px - PSXLen) <= eps && fabs(yf - PSYLen) <= eps) ||
                             (fabs(py - PSYLen) <= eps && fabs(xf - PSYLen) <= eps)) {
                          fprintf(PSFd,"%5.2f %5.2f l %% ddd\n",PSXLen,PSYLen);
                          fflag = 1;
                      }
                      if (fflag) {
                          if (PMLW > 0.0)
                              fprintf(PSFd,"gsave\n");

                          if (PMGSFlg == 1)
                              g = PMGS * (1.0 + (double)(1 - level) / (double)nlev);
                          else if (nlev > 1)  
                              g = PMGS - (double)(level - 1) * (PMGS - PMGS1) / (double)(nlev - 1);

                          ps_fill(g);
                          if (PMLW > 0.0)
                              fprintf(PSFd,"grestore\n");
                      }
                  }
                  fprintf(PSFd,"stroke\n");
                  fprintf(PSFd,"grestore\n");
                  break;
        default:  fprintf(PSFd,"%% cont_draw n=%d level=%d\n",n,level);
                  break;
    }
}


/*--------------------------------------------------------------------------*/
/*  pl_plotcm   Contour plot.                                               */
/*                                                                          */
/*              plotcm(                                                     */
/*                  opt=...,    option, def. 1                              */
/*                              1 use all matrix elements                   */
/*                              2 use only non-negative matrix elements     */
/*                  dopt=...,   option for matrix rows, def. 1              */
/*                              1 use matrix rows in standard order         */
/*                              2 use matrix rows in reverse order          */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def.                            */
/*                  x=...,      level                                       */
/*                  gs=...,     grey-scale option                           */
/*              ) = matrix_name;                                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_plotcm(void)
{
    register int i,j,k;
    int nn,err,r,pcbita,pcfva,im;
    double tmp;
           
    pcbita = pcfva = 0;
    err = -1;

    if (check_pcmd(2,2))
        return(-1);

    printf1("Contour plot.\n");
    if (parm(CmdBuf + 6,14,1))    /* get parameters */
        goto PLCMFin;

    if (PMNVTyp != 1 || PMNV != 1) {
        printf1("Error: need exactly one matrix on right-hand side.\n");
        goto PLCMFin;
    }
    if (PMOPT != 2)
        PMOPT = 1;
    if (PMDOPT != 2)
        PMDOPT = 1;

    im = PMVIdx[0];
    PCNY = MatRow[im];
    PCNX = MatCol[im];

    PCDX = UXLen / (double)PCNX;
    PCDY = UYLen / (double)PCNY;

    if (!(PCFV = (float *) calloc(PCNX * PCNY + 1,sizeof(float)))) {
        p_err(-2,1);
        goto PLCMFin;
    }
    memrq(PCNX * PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * PCNX * PCNY * PMNTP / 8 + 1;
    if (!(PCBitM = (char *)calloc(nn,sizeof(char)))) {
        p_err(-2,1);
        goto PLCMFin;
    }
    memrq(nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= PCNX; ++i) {
        for (j = 0; j < PCNY; ++j) {
            tmp = (float)MatVal[im][j * PCNX + i];       
            if (PMOPT == 2 && tmp < 0.0)
                tmp = PCMAX - 1.0;  

            if (PMDOPT == 2)
                k = j + 1;
            else
                k = PCNY - j;

            PCFV[(i - 1) * PCNY + k] = (float)tmp;
        }
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    cont_plot(PMNTP,PMTP,1);
    err = 0;

PLCMFin:
    if (pcfva) {
        free((char *)PCFV);
        memrq(-PCNX * PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)PCBitM);
        memrq(-nn,sizeof(char));
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotr    Plot grey-scaled relief.                                    */
/*                                                                          */
/*              plotr(                                                      */
/*                  dopt=...,   option for matrix rows, def. 1              */
/*                              1 use matrix rows in standard order         */
/*                              2 use matrix rows in reverse order          */
/*                  gs = min,max, def. 0,1                                  */
/*              ) = matrix_name;                                            */

int pl_plotr(void)
{
    register int i,j,k;
    int err,im,row,col,nm;
    double a,b,d,dx,dy,x,y,tmp,vmin,vmax,g;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,14,1))    /* get parameters */
        goto PLRFin;

    if (PMNVTyp != 1 || PMNV != 1) {
        printf1("Error: need exactly one matrix on right-hand side.\n");
        goto PLRFin;
    }
    if (PMDOPT != 2)
        PMDOPT = 1;

    a = 0.0;
    b = 1.0;
    if (PMGSFlg == 2) {
        a = PMGS;
        b = PMGS1;
    }
    im = PMVIdx[0];
    row = MatRow[im];
    col = MatCol[im];

    dx = UXLen / (double)col;
    dy = UYLen / (double)row;

    vmin = DBLMAX;
    vmax = -1.0;      
    nm = 0;
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = MatVal[im][i * col + j];
            if (tmp < 0.0)
                nm++;
            else {
                vmin = dmin(vmin,tmp);
                vmax = dmax(vmax,tmp);
            }
        }
    }
    printf1("Number of cells: %d x %d = %d\n",row,col,row * col);
    printf1("Number of empty cells: %d\n",nm);
    if (nm == row * col) {
        err = 0;
        goto PLRFin;
    }
    printf("Minimum: %g\n",vmin);
    printf("Maximum: %g\n",vmax);
    d = vmax - vmin;
    if (d < EPSI1) {
        printf1("Error: values almost identical.\n");
        goto PLRFin;
    }

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_lwidth(PMLW);
    ps_ltyp(PMLT);
   
    for (i = row - 1; i >= 0; --i) {

        if (PMDOPT == 2)
            k = row - i - 1;
        else
            k = i;

        for (j = 1; j <= col; ++j) {

            tmp = MatVal[im][k * col + j];
            if (tmp < 0.0)
                continue;

            g = (vmax - tmp) / d;
            g = a + g * (b - a);
            x = PA1[0] + dx * (double)(j - 1);
            y = PA1[1] + dy * (double)(row - i - 1);

            fprintf(PSFd,"gsave %4.2f setgray\n",g);    
            ps_2dplot(x,y,0);    
            ps_2dplot(x,y + dy,1);    
            ps_2dplot(x + dx,y + dy,1);    
            ps_2dplot(x + dx,y,1);    
            ps_2dplot(x,y,1);    
            fprintf(PSFd,"fill\ngrestore\n");
        }
    }
    fprintf(PSFd,"grestore\n");
/**
    pl_plotr_c(im,row,col,30.0,dx,dy);
    pl_plotr_c(im,row,col,50.0,dx,dy);
    pl_plotr_c(im,row,col,80.0,dx,dy);
    pl_plotr_c(im,row,col,100.0,dx,dy);
    pl_plotr_c(im,row,col,120.0,dx,dy);
**/ 

    err = 0;

PLRFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotr_c(                                                             */

void pl_plotr_c(int im,int row,int col,double lev,double dx,double dy)
{
    register int i,j;
    double x,y,a,b;                                

    fprintf(PSFd,"gsave\n");    
    for (i = row - 1; i >= 0; --i) {
        for (j = 1; j <= col; ++j) {
            a = MatVal[im][i * col + j];
            x = PA1[0] + dx * (double)(j - 1);
            y = PA1[1] + dy * (double)(row - i - 1);

            if (j < col) {
                b = MatVal[im][i * col + j + 1];
                if ((a < lev && b >= lev) || (b < lev && a >= lev)) {
                    ps_2dplot(x + dx,y,0);    
                    ps_2dplot(x + dx,y + dy,1);    
                    fprintf(PSFd,"stroke\n");
                }
            }
            if (i > 0) {
                b = MatVal[im][(i - 1) * col + j];
                if ((a < lev && b >= lev) || (b < lev && a >= lev)) {
                    ps_2dplot(x,y + dy,0);    
                    ps_2dplot(x + dx,y + dy,1);    
                    fprintf(PSFd,"stroke\n");
                }
            }
        }
    }
    fprintf(PSFd,"grestore\n");
}
