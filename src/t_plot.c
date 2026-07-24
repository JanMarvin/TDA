/****************************************************************************/
/*  t_plot                                                                  */
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
#include "t_graph.h"     
#include "t_spl.h"
#include "t_int.h"
#include "t_gmin.h"
#include "t_var.h"
#include "t_areg.h"
#include "t_sort.h"
#include "t_ml.h"
#include "t_eval3.h"
#include "t_gm.h"
#include "t_mat.h"

/*  functions in t_plot.c */

int check_pcmd(int opt,int dim);              
int check_ps(int dim);              
void upd_bbox(int opt,double x,double y); 
void set_clip(void); 
void ps_ltyp(int typ);
void ps_lwidth(double lw);
void ps_setfs(double fs);
void ps_fill(double g); 
int ps_cclip(double x,double y);
void ps_2dplot(double x,double y, int opt);
double ps_2dx(double x);
double ps_2dy(double y);
void ps_sym(int typ, double px, double py,double siz);
void plot_str(double x,double y,char *s, double siz, int opt,int adj,
    int r,int xopt,int cflag);
void ps_nlab(int n,double x,double y,double fs,int opt);
int ps_grid(int typ);
int pl_label(int typ);
int pl_text(void);
int pl_frame(void);
int pl_rec(void);
int pl_plotp(int typ);
int pl_plotm(void);
int pl_plotf(void);
int pl_ploth(void);
int pl_plotd(void);
int pl_ploto(int typ);
int pl_plotk(void); 
int pl_plotch(void);
int pl_plots(int typ);

/*--------------------------------------------------------------------------*/
int PSFFlg = 0;             /* set if PostScript output file defined        */
int PS3DFlg = 0;            /* set by psetup3().                            */
int PSPROJ = 0;             /* set by psetupg()                             */

char PSFName[FNMaxLen+1];   /* name of PostScript output file               */
FILE *PSFd;                 /* file handle for PostScript output file       */
int BBLX,BBLY,BBUX,BBUY;    /* bounding box                                 */
int PSONUM = 0;             /* plot object number                           */ 

double PtMM = 2.835;        /* points per mm                                */

double LWCSMax = 0.0;       /* max line width used for axes                 */
double LWDef = 0.2;         /* default line width                           */
double LW1Def = 0.05;       /* default second line width (grids etc.)       */
double FSDef = 2.0;         /* default font size                            */
double PTLen = 1.8;         /* tick length in mm                            */

int PSXOrg = 0;             /* origin of PostScript figure                  */
int PSYOrg = 0;
int XOrg = 150;             /* origin of PostScript figure, redefined       */
int YOrg = 460;             /* by psorg commands.                           */
double SCALXFac = 1.0;      /* PostScript x scaling factor                  */
double SCALYFac = 1.0;      /* PostScript y scaling factor                  */
double ROTFac = 0.0;        /* PostScript rotation factor                   */

double PXLen = 120.0;       /* length of X axis in mm                       */
double PYLen =  80.0;       /* length of Y axis in mm                       */
double PSXLen = 0.0;        /* lenght of X axis in points                   */
double PSYLen = 0.0;        /* lenght of Y axis in points                   */

int PSLog[2];               /* flags for logarithmic axes                   */
double PA1[2];              /* logical start of axis                        */
double PA2[2];              /* logical end of axis                          */

double UXLen;               /* length of x axis in user units               */
double UYLen;               /* length of y axis in user units               */
double UZLen;               /* length of z axis in user units               */

double PX3D[2];             /* logical x axis                               */
double PY3D[2];             /* logical y axis                               */
double PZ3D[2];             /* logical z axis                               */
double PX3C = 0.0;          /* center of x axis                             */
double PY3C = 0.0;          /* center of y axis                             */
double PZ3C = 0.0;          /* center of z axis                             */

double PSLon = 30.0;        /* direction of projection, longitude           */
double PSLat = 30.0;        /* direction of projection, latitude            */

double PSPR11,PSPR12,PSPR13;    /* projection matrix                        */
double PSPR21,PSPR22,PSPR23;
double PSPR31,PSPR32,PSPR33;

double PSPI11,PSPI12,PSPI13;    /* inverse of projection matrix             */
double PSPI21,PSPI22,PSPI23;
double PSPI31,PSPI32,PSPI33;

/* ------------------------------------------------------------------------ */
/*  check_pcmd(opt,dim)   Check plot command in CmdBuf. First call          */  
/*                        check_cmd(opt), then check for valid PostScript   */
/*                        output file.                                      */  
/*                        Return 0 if OK, -1 if error.                      */

int check_pcmd(int opt,int dim)       
{
    int err;
    err = check_cmd(opt);
    if (err == 0)  
        err = check_ps(dim);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_ps(dim)       Check for valid PostScript output file.             */
/*                      Return 0 if OK, -1 if error.                        */

int check_ps(int dim)            
{
    if (PSFFlg != 2) {
        printf1("Error: need PostScript output file and coordinates.\n");
        return(-1);
    }
    if (dim == 2 && PS3DFlg) {
        printf1("Error: current coordinate system is 3-dimensional.\n");
        return(-1);
    }
    else if (dim == 3 && PS3DFlg == 0) {
        printf1("Error: current coordinate system is 2-dimensional.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  upd_bbox(opt,x,y)   Update bounding box for point (x,y). If opt == 0    */
/*                      user coordinates, otherwise PostScript coordinates. */

void upd_bbox(int opt,double x,double y)
{
    int ix,iy;

    if (opt == 0) {
        ix = (int)ps_2dx(x) + XOrg;
        iy = (int)ps_2dy(y) + YOrg;
    }
    else {
        ix = (int)x + XOrg;
        iy = (int)y + YOrg;
    }
    if (BBLY > iy)
        BBLY = iy;
    if (BBUY < iy)
        BBUY = iy;
    if (BBLX > ix)
        BBLX = ix;
    if (BBUX < ix)
        BBUX = ix;
}

/* ------------------------------------------------------------------------ */
/*  set_clip()      Set default clipping region. Observe LWCSMax            */

void set_clip(void)
{
    double x;  
    x = LWCSMax * PtMM;
    fprintf(PSFd,"%5.2f %5.2f m\n%5.2f %5.2f l\n%5.2f %5.2f l\n%5.2f %5.2f l\nclosepath\nclip\nnewpath\n",
                                  x,x,x,PSYLen,PSXLen,PSYLen,PSXLen,x);
}

/* ------------------------------------------------------------------------ */
/*  ps_ltyp(typ)  Select line type for PostScript output.                   */
/*      Note: only some standard predefined formats are used here.          */

void ps_ltyp(int typ)
{
    switch (typ) {
        case  2:  fprintf(PSFd,"[1 2]"); break;
        case  3:  fprintf(PSFd,"[1 3]"); break;
        case  4:  fprintf(PSFd,"[1 5]"); break;
        case  5:  fprintf(PSFd,"[2 2]"); break;
        case  6:  fprintf(PSFd,"[4 2]"); break;
        case  7:  fprintf(PSFd,"[6 4]"); break;
        case  8:  fprintf(PSFd,"[6 2]"); break;
        case  9:  fprintf(PSFd,"[9 4]"); break;
        default:  fprintf(PSFd,"[]"); break;
    }
    fprintf(PSFd," 0 setdash\n");
}

/* ------------------------------------------------------------------------ */
/*  ps_lwidth(lw)   Set line width to lw.                                   */

void ps_lwidth(double lw)
{
    fprintf(PSFd,"%6.4f setlinewidth\n",PtMM * lw);
}

/* ------------------------------------------------------------------------ */
/*  ps_setfs(fs)    Set font size.                                          */

void ps_setfs(double fs)
{
    fprintf(PSFd,"/fsiz %7.4f def FT (",fs); 
}

/* ------------------------------------------------------------------------ */
/*  ps_fill     PostScript fill with greyscale value command.               */

void ps_fill(double g)
{
    if (g >= 0.0 && g <= 1.0)  
        fprintf(PSFd," %6.4f setgray fill ",g);
}

/* ------------------------------------------------------------------------ */
/*  ps_cclip(x,y)   x,y in PostScript coordinates. Return 1 if (x,y) is     */
/*                  inside coordinate system, otherwise 0.                  */

int ps_cclip(double x,double y)
{
    if (x < (double)XOrg || x > (double)XOrg + PSXLen)
        return(0);
    if (y < (double)YOrg || y > (double)YOrg + PSYLen)
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dplot   PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values.                                       */

void ps_2dplot(double x,double y, int opt)
{
    double px,py;

    if (!PSLog[0])
        px = (x - PA1[0]) / UXLen;
    else if (x <= 0.0)
        px = 0.0;
    else
        px = rlog(x / PA1[0]) / rlog(PA2[0] / PA1[0]);

    if (!PSLog[1])
        py = (y - PA1[1]) / UYLen;
    else if (y <= 0.0)
        py = 0.0;
    else
        py = rlog(y / PA1[1]) / rlog(PA2[1] / PA1[1]);

    px *= PtMM * PXLen;
    py *= PtMM * PYLen;

    if (!opt)  
        fprintf(PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(PSFd,"%5.2f %5.2f ",px,py);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dx(x)   return x in PostScript coordinates.                         */

double ps_2dx(double x)
{
    double px;

    if (!PSLog[0])
        px = (x - PA1[0]) / UXLen;
    else if (x <= 0.0)
        px = 0.0;
    else
        px = rlog(x / PA1[0]) / rlog(PA2[0] / PA1[0]);

    px *= PtMM * PXLen;
    return(px);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dy(y)   return y in PostScript coordinates.                         */

double ps_2dy(double y)
{
    double py;

    if (!PSLog[1])
        py = (y - PA1[1]) / UYLen;
    else if (y <= 0.0)
        py = 0.0;
    else
        py = rlog(y / PA1[1]) / rlog(PA2[1] / PA1[1]);

    py *= PtMM * PYLen;
    return(py);
}

/* ------------------------------------------------------------------------ */
/*  ps_sym(typ,x,y,siz)     plot symbol (typ) at PostScript coordinates x,y */

void ps_sym(int typ, double px, double py,double siz)
{
            
    if (typ >= 1 && typ <= 17)
        fprintf(PSFd,"%%#symbol: %d %5.2f %5.2f %5.2f\n",typ,px,py,siz);
            
    switch (typ) {
        case  1:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m cross\n",px,py,px,py);
                  break;
        case  2:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m xsym\n",px,py,px,py);
                  break;
        case  3:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m cross xsym\n",px,py,px,py);
                  break;
        case  4:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle stroke\n",px,py,px,py);
                  break;
        case  5:  fprintf(PSFd,"%5.2f %5.2f circle fill\n",px,py);
                  break;
        case  6:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle ",px,py,px,py);
                  fprintf(PSFd,"stroke %5.2f %5.2f m cross\n",px,py);
                  break;
        case  7:  fprintf(PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle ",px,py,px,py);
                  fprintf(PSFd,"stroke %5.2f %5.2f m xsym\n",px,py);
                  break;
        case  8:  fprintf(PSFd,"%5.2f %5.2f m square stroke grestore\n",px,py);
                  break;
        case  9:  fprintf(PSFd,"%5.2f %5.2f m square fill grestore\n",px,py);
                  break;
        case 10:  fprintf(PSFd,"%5.2f %5.2f m trian1 stroke grestore\n",px,py);
                  break;
        case 11:  fprintf(PSFd,"%5.2f %5.2f m trian1 fill grestore\n",px,py);
                  break;
        case 12:  fprintf(PSFd,"%5.2f %5.2f m trian2 stroke grestore\n",px,py);
                  break;
        case 13:  fprintf(PSFd,"%5.2f %5.2f m trian2 fill grestore\n",px,py);
                  break;
        case 14:  fprintf(PSFd,"%5.2f %5.2f m square stroke grestore cross\n",px,py);
                  break;
        case 15:  fprintf(PSFd,"%5.2f %5.2f m rhomb stroke grestore\n",px,py);
                  break;
        case 16:  fprintf(PSFd,"%5.2f %5.2f m rhomb fill grestore\n",px,py);
                  break;
        case 17:  fprintf(PSFd,"%5.2f %5.2f m rhomb stroke grestore cross\n",px,py);
                  break;
        default:  break;
    }
}

/* ------------------------------------------------------------------------ */
/*  plot_str (x,y,s,siz,opt,adj,r,xopt,cflag)                               */
/*                                                                          */
/*  Plot string s in size at x,y.                                           */
/*  If opt == 1 don't end with show.                                        */
/*  If adj = 1 center, 2 right justified                                    */
/*  If r != 0 rotate.                                                       */
/*  If xopt != 0 PostScript coordinates                                     */
/*  If cflag != 0 put string in white bounding box                          */

void plot_str(double x,double y,char *s, double siz, int opt,int adj,int r,
    int xopt,int cflag)
{
    register int sflg,scnt;
    register char c, *p;

    if (xopt == 0) {
        x = ps_2dx(x);
        y = ps_2dy(y);
    }
    fprintf(PSFd,"%%#text: %5.2f %5.2f %5.2f %d %d %s\n",x,y,siz,adj,r,s);

    fprintf(PSFd,"gsave\n");
    fprintf(PSFd,"%5.2f %5.2f translate\n0 0 m\n",x,y);

    if (r != 0)
        fprintf(PSFd,"%d rotate\n",r);

    p = s;
    sflg = 0;
    if (*p != '@') {
        fprintf(PSFd,"/fsiz %5.2f def FT (",siz);  /* default is TFont */
        sflg = 1;
    }
    while (*p) {
        c = *p;
        if (c == '@' && isdigit((int) *(p+1))) {  /* switch to symbol font */
                          
            if (sflg == 1)
                fprintf(PSFd,")\nshow\n"); 
            fprintf(PSFd,"/fsiz %5.2f def FS (\\",siz); 
            sflg = 2;
            scnt = 0;
        }
        else {
            if (sflg == 2) {
                if (scnt >= 3 || isdigit((int)c) == 0) {
                    fprintf(PSFd,")\nshow\n/fsiz %5.2f def FT (",siz); 
                    sflg = 1;
                    scnt = 0;
                }
                else
                    scnt++;
            }
            if (c == '(' || c == ')')  
                fprintf(PSFd,"\\");
            fprintf(PSFd,"%c",c);
        }
        p++;
    }
    fprintf(PSFd,") ");
    if (opt == 0) {
        if (adj == 1)
            fprintf(PSFd,"center\n");
        else if (adj == 2)
            fprintf(PSFd,"aright\n");
        if (cflag)
            fprintf(PSFd,"sclear ");
        fprintf(PSFd,"show\n");
    }
    else
        fprintf(PSFd,"\n");
    fprintf(PSFd,"grestore\n");
}

/* ------------------------------------------------------------------------ */
/*  ps_nlab(n,x,y,fs,opt)                                                   */
/*                                                                          */
/*  Plot numerical label n at (x,y) with font size fs                       */
/*  opt = 0 : centered for x and y direction                                */  
/*  opt = 1 : centered and white rectangle around character                 */
/*  opt = 2 : right justified                                               */

void ps_nlab(int n,double x,double y,double fs,int opt)
{
    double gx,gy,ax,ay;
    int adj;
    char buf[20];

    ax = fs * UXLen / PXLen;
    gx = 0.6 * ax;                   
    gy = 0.9 * fs * UYLen / PYLen;
    if (n > 9)
        gx += 0.7 * ax;
    if (n > 99)
        gx += 0.7 * ax;
    gx /= 1.9;
    gy /= 1.9;
    ax = 1.6 * gx;
    ay = 1.6 * gy;

    if (opt == 1) {
        fprintf(PSFd,"gsave\n");       
        ps_lwidth(0.0);
        ps_2dplot(x - ax,y - ay,0);    
        ps_2dplot(x - ax,y + ay,1);    
        ps_2dplot(x + ax,y + ay,1);    
        ps_2dplot(x + ax,y - ay,1);    
        ps_fill(1.0);
        fprintf(PSFd,"grestore\n");   
    }
    sprintf(buf,"%d",n);
    if (opt <= 1)
        adj = 1;
    else
        adj = 2;
             
    plot_str(x,y - gy,buf,1.5 * fs * PtMM,0,adj,0,0,0);
}

/*--------------------------------------------------------------------------*/
/*  ps_grid(typ)        plot grid lines. Syntax:                            */
/*                      plxgrid(lt=,lw=) = x1,x2,...                        */
/*                      typ: 0 = x-axis, 1 = y-axis.                        */
/*                      Return 0 if OK, -1 if error.                        */

int ps_grid(int typ)
{
    int i,err;
    double x;

    err = -1;

    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 7,7,1))    /* get parameters */
        goto PSGFin;

    if (PMRHSN > 0) {
        if (PMLTFlg == 0)
            PMLT = 3;           /* default line type */

        if (PMLWFlg == 0)
            PMLW = 0.05;        /* default line width */

        fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

        ps_ltyp(PMLT);
        ps_lwidth(PMLW);
           
        for (i = 0; i < PMRHSN; ++i) {
            x = PMRHSX[i];
            if (typ == 0) {
                if (x < PA1[1] || x > PA2[1])
                    continue;

                ps_2dplot(PA1[0],x,0);
                ps_2dplot(PA2[0],x,1);
            }
            else {
                if (x < PA1[0] || x > PA2[0])
                    continue;

                ps_2dplot(x,PA1[1],0);
                ps_2dplot(x,PA2[1],1);
            }
            fprintf(PSFd,"stroke\n");
        }
    }
    err = 0;

PSGFin:
    p_clean();
    return(err);
}
   
/*--------------------------------------------------------------------------*/
/*  pl_label(typ)       plot label.                                         */
/*                                                                          */
/*              typ 0   plabel(sc=,fs=) = string;                           */
/*              typ 1   pxlabel(sc=,fs=) = string;                          */
/*              typ 2   pylabel(sc=,fs=) = string;                          */
/*                      Return 0 if OK, -1 if error.                        */

int pl_label(int typ)
{
    int r,err = -1;
    double x,y,d;

    if (check_pcmd(1,2))
        return(-1);

    r = 6;
    if (typ)
        r++;

    if (parm(CmdBuf + r,8,1))    /* get parameters */
         goto PLLFin;

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    if (typ == 0) {

        d = 3.0 + PMSC;
        if (PMFSFlg == 0)
            PMFS = 3.0;

        x = PSXLen / 2.0;
        y = PSYLen + d * PtMM;         
        plot_str(x,y,PMRHSTR,1.5 * PMFS * PtMM,0,1,0,1,0);
        y += 1.5 * PMFS * PtMM;
        upd_bbox(1,x,y);
    }
    else if (typ == 1) {    /* x axis */

        if (PMFSFlg == 0)
            PMFS = 2.4;

        d = 4.8 + PMSC + 1.6 * PMFS;

        x = PSXLen / 2.0;
        y = -d * PtMM;
        plot_str(x,y,PMRHSTR,1.5 * PMFS * PtMM,0,1,0,1,0);
        y -= 0.5 * PMFS * PtMM;
        upd_bbox(1,x,y);
    }
    else if (typ == 2) {        /* y axis */

        d = 7.0 + PMSC;
        if (PMFSFlg == 0)
            PMFS = 2.4;

        x = -d * PtMM;           
        y = PSYLen / 2.0;
        plot_str(x,y,PMRHSTR,1.5 * PMFS * PtMM,0,1,90,1,0);
        x -= 1.5 * PMFS * PtMM;
        upd_bbox(1,x,y);
    }
    err = 0;

PLLFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_text()   Plot a string.                                              */
/*                      plttext(pos=,fs=,sc=,r=,s=) = string;               */
/*              pltext(                                                     */
/*                  xy=...,     x,y position                                */
/*                  fs=...,     font size, def. 2 (mm)                      */
/*                  sc=...,     if sc=1 centered, def. 0                    */
/*                  s=...,      symbol                                      */
/*                  r=...,      rotation, def. 0                            */ 
/*                  wf=...,     if wf=1, put string in white bounding       */
/*                              box, def. 0                                 */
/*              ) = string;                                                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_text(void)
{
    int cen,len,err = -1;
    double x,y;

    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 6,8,1))    /* get parameters */
        goto PLTFin;

    if (PMXYFlg == 0) {
        p_err(-1,1);
        goto PLTFin;
    }
    if (PMWF != 0)
        PMWF = 1;

    if (PMR)  
        cen = 0;
    else  
        cen = (int)PMSC;
       
    x = ps_2dx(PMX);
    y = ps_2dy(PMY);

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    if (PMS) {
        PMR = cen = 0;
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);
        ps_sym(PMS,x - 1.5 * PtMM * PMFS,y + 0.5 * PtMM * PMFS,PtMM * PMFS / 2.0);
    }
    plot_str(PMX,PMY,PMRHSTR,1.5 * PMFS * PtMM,0,cen,PMR,0,PMWF);

    /* update bounding box */

    len = strlen(PMRHSTR) + 1;
    upd_bbox(1,x,y);     
    upd_bbox(1,x,y + 1.5 * PMFS * PtMM);
    upd_bbox(1,x,y - 0.8 * PMFS * PtMM);

    if (cen == 0) {
        upd_bbox(1,x - 1.5 * PMFS * PtMM,y);
        upd_bbox(1,x + len * PMFS * PtMM,y);
    }
    else {
        len /= 2.0;
        upd_bbox(1,x + len * PMFS * PtMM,y);
        upd_bbox(1,x - len * PMFS * PtMM,y);
    }
    err = 0;

PLTFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_frame()          plot frame.                                         */
/*                      plframe(lt,lw,gs)                                   */
/*                      Return 0 if OK, -1 if error.                        */

int pl_frame(void)
{
    int err = -1;

    if (check_pcmd(1,0))
        return(-1);

    if (parm(CmdBuf + 7,0,0))    /* get parameters */
        goto PLFRFin;

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);


    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave %4.2f setgray\n",PMGS);    
        fprintf(PSFd,"0 0 m\n");
        fprintf(PSFd,"%5.2f 0 l\n",PSXLen);
        fprintf(PSFd,"%5.2f %5.2f l\n",PSXLen,PSYLen);
        fprintf(PSFd,"0 %5.2f l\n",PSYLen);
        fprintf(PSFd,"0 0 l\n");
        fprintf(PSFd,"fill\ngrestore\n");
    }
    if (PMLW > 0.0) {
        fprintf(PSFd,"0 0 m\n");
        fprintf(PSFd,"%5.2f 0 l\n",PSXLen);
        fprintf(PSFd,"%5.2f %5.2f l\n",PSXLen,PSYLen);
        fprintf(PSFd,"0 %5.2f l\n",PSYLen);
        fprintf(PSFd,"0 0 l\n");
        fprintf(PSFd,"closepath\nstroke\n");
    }
    fprintf(PSFd,"grestore\n");
   
    err = 0;

PLFRFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_rec()            plot rectangle. Syntax:                             */
/*                      plrec(lt,lw,gs,r) = x,y,xlen,ylen;                  */
/*                      gs = grey value, r = rotate                         */
/*                      Return 0 if OK, -1 if error.                        */

int pl_rec(void)
{
    int err = -1;
    double x,y,x1,y1;

    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,7,1))    /* get parameters */
        goto PLRFin;

    if (PMRHSN != 4) {
        p_err(-1,1);
        goto PLRFin;
    }
    x = ps_2dx(PMRHSX[0]);
    y = ps_2dy(PMRHSX[1]);
    x1 = ps_2dx(PMRHSX[0] + PMRHSX[2]);
    y1 = ps_2dy(PMRHSX[1] + PMRHSX[3]);
  
    upd_bbox(1,x,y);        /* update bounding box */
    upd_bbox(1,x1,y1);

    if (fabs(PMRHSX[2]) < EPSI || fabs(PMRHSX[3]) < EPSI) {
        err = 0;
        goto PLRFin;
    }

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);

    fprintf(PSFd,"%5.2f %5.2f translate\n",x,y);
    x1 -= x;
    y1 -= y;
    x = y = 0;

    if (PMR != 0 && PMR < 360 && PMR > -360)
        fprintf(PSFd,"%d rotate\n",PMR);

    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n%4.2f setgray\n",PMGS);    
        fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
        fprintf(PSFd,"%5.2f %5.2f l\n",x1,y);
        fprintf(PSFd,"%5.2f %5.2f l\n",x1,y1);
        fprintf(PSFd,"%5.2f %5.2f l\n",x,y1);
        fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
        fprintf(PSFd,"fill\ngrestore\n");
    }
    if (PMLW > 0.0) {
        fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
        fprintf(PSFd,"%5.2f %5.2f l\n",x1,y);
        fprintf(PSFd,"%5.2f %5.2f l\n",x1,y1);
        fprintf(PSFd,"%5.2f %5.2f l\n",x,y1);
        fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
        fprintf(PSFd,"closepath\nstroke\n");
    }
    fprintf(PSFd,"grestore\n");
   
    err = 0;

PLRFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotp(typ)       Plot polygon.                                       */
/*                                                                          */
/*      typ 0:  plotp(lt,lw,gs,s,fs,r,dir,a,nc,ns) = x1,y1,x2,y2,...;       */
/*      typ 1:  plot (lt,lw,gs,s,fs,r,dir,a,nc,sel,ns) = VX,VY;             */
/*              gs = grey value, s = symbol, fs=symbol size                 */
/*              r = connect                                                 */
/*              a = arrow                                                   */  
/*              dir = 1 (plotr) 2 (plotl)                                   */
/*              ns=1 only horizontal lines                                  */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int pl_plotp(int typ)
{
    int i,ix,iy,n,nn,err,first;
    double x,y,xa,ya;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (typ == 0) {
        if (parm(CmdBuf + 5,7,1))    /* get parameters */
            goto PLPFin;
    
        n = PMRHSN;
        if (n < 4 || (n / 2) * 2 != n) {
            p_err(-25,1);
            goto PLPFin;
        }
    }
    else {
        if (parm(CmdBuf + 4,4,1))    /* get parameters */
            goto PLPFin;

        if (NOC < 2)
            goto PLPFin;

        n = NOC;
        if (PMNV != 2) {
            p_err(-1,1);
            goto PLPFin;
        }
        ix = PMVIdx[0];
        iy = PMVIdx[1];
    }


    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();

    /* first process fill option */

    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        ps_lwidth(0.0);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x = ps_2dx(PMRHSX[i]);
                y = ps_2dy(PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(get_data(ix,i));       
                y = ps_2dy(get_data(iy,i));       
            }
            if (PMNC)
                upd_bbox(1,x,y);              /* update bounding box */
            nn++;
            if (first) {
                fprintf(PSFd,"%5.2f %5.2f m\n",x,0.0);
                fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
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
                i++;
        }
        if (nn) {
            fprintf(PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);
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
                x = ps_2dx(PMRHSX[i]);
                y = ps_2dy(PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(get_data(ix,i));       
                y = ps_2dy(get_data(iy,i));       
            }
            if (PMNC)
                upd_bbox(1,x,y);              /* update bounding box */
            nn++;

            if (first) {
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                if (PMDIR == 1) {
                    fprintf(PSFd,"%5.2f %5.2f l\n",x,ya);
                    if (PMNS != 1)
                        fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
                    else
                        fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
                }
                else if (PMDIR == 2) {
                    if (PMNS != 1) 
                        fprintf(PSFd,"%5.2f %5.2f l\n",xa,y);
                    else
                        fprintf(PSFd,"%5.2f %5.2f m\n",xa,y);
                    fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
                }
                else
                    fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
            }
            if (PMR)  
                fprintf(PSFd,"0 %5.2f rl 0 %5.2f rl\n",-y,y);

            i++;   
            if (typ == 0)
                i++;
            if (i >= n)
                break;
            xa = x;
            ya = y;
        }
        if (nn >= 2 && PMAFlg && PMA1 > 0.0 && PMA2 > 0.0) {       /* plot arrow */
        
            fprintf(PSFd,"gsave\ncurrentpoint\nstroke m\n");
            x -= xa;
            y -= ya;
            fprintf(PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
            fprintf(PSFd,"%5.2f %5.2f scale\n",PMA1,PMA2);
            fprintf(PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
        }
        if (nn)
            fprintf(PSFd,"stroke\n");
    }
    if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

        i = 0;
        while (i < n) {

            if (typ == 0) {
                x = ps_2dx(PMRHSX[i]);
                y = ps_2dy(PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(get_data(ix,i));       
                y = ps_2dy(get_data(iy,i));       
            }
            if (PMNC)
                upd_bbox(1,x,y);                 /* update bounding box */
            ps_sym(PMS,x,y,PtMM * PMFS / 2.0);
            i++;   
            if (typ == 0)
                i++;
        }
    }
    fprintf(PSFd,"grestore\n");
    err = 0;

PLPFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotm()          Plot polygon, separately for each data matrix row.  */
/*                                                                          */
/*                      plotm(lt,lw,gs,s,fs,r,dir,a,nc,sel) = VX,VY;        */
/*                      gs = grey value, s = symbol, fs=symbol size         */
/*                      r = connect                                         */
/*                      a = arrow                                           */  
/*                      dir = 1 (plotr) 2 (plotl)                           */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotm(void)
{
    int j,icase,ix,iy,n,err,first;
    double x,y,xa,ya;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,4,1))    /* get parameters */
        goto PLMFin;

    n = PMNV / 2;
    if (PMNV < 2 || 2 * n != PMNV) {
        p_err(-1,1);
        goto PLMFin;
    }
       
    for (icase = 0; icase < NOC; ++icase) {

        fprintf(PSFd,"\n%%#%d: %s [case %d]\n",++PSONUM,CmdBuf,icase+1);

        fprintf(PSFd,"gsave\n");
        if (PMNC != 1)
            set_clip();

        /* first process fill option */

        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            ps_lwidth(0.0);
            xa = ya = 0.0;
            first = 1;
            for (j = 0; j < PMNV; j += 2) {
                ix = PMVIdx[j];
                iy = PMVIdx[j + 1];
                x = ps_2dx(get_data(ix,icase));       
                y = ps_2dy(get_data(iy,icase));       

                if (PMNC)
                    upd_bbox(1,x,y);              /* update bounding box */

                if (first) {
                    fprintf(PSFd,"%5.2f %5.2f m\n",x,0.0);
                    fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
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
            }
            fprintf(PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\nnewpath\n");
        }

        /* now the standard polygon */

        if (PMLW > 0.0)
            ps_lwidth(PMLW);

        if (PMLW > 0.0 && PMLT > 0) {
            ps_ltyp(PMLT);
            xa = ya = 0.0;

            first = 1;
            for (j = 0; j < PMNV; j += 2) {
                ix = PMVIdx[j];
                iy = PMVIdx[j + 1];
                x = ps_2dx(get_data(ix,icase));       
                y = ps_2dy(get_data(iy,icase));       

                if (PMNC)
                    upd_bbox(1,x,y);              /* update bounding box */

                if (first) {
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
                if (PMR)  
                    fprintf(PSFd,"0 %5.2f rl 0 %5.2f rl\n",-y,y);

            }
            if (n >= 2 && PMAFlg && PMA1 > 0.0 && PMA2 > 0.0) {       /* plot arrow */
        
                fprintf(PSFd,"gsave\ncurrentpoint\nstroke m\n");
                x -= xa;
                y -= ya;
                fprintf(PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
                fprintf(PSFd,"%5.2f %5.2f scale\n",PMA1,PMA2);
                fprintf(PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
            }
            fprintf(PSFd,"stroke\n");
        }
        if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

            fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

            for (j = 0; j < PMNV; j += 2) {
                ix = PMVIdx[j];
                iy = PMVIdx[j + 1];
                x = ps_2dx(get_data(ix,icase));       
                y = ps_2dy(get_data(iy,icase));       

                if (PMNC)
                    upd_bbox(1,x,y);                 /* update bounding box */
                ps_sym(PMS,x,y,PtMM * PMFS / 2.0);
            }
        }
        fprintf(PSFd,"grestore\n");
    }
    err = 0;

PLMFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotf()          Plot function.                                      */
/*                                                                          */
/*  plotf(                                                                  */
/*      rx=...,     range of function argument                              */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      gs=...,     grey scale value, def. 1 (white)                        */
/*      nc=...,     no clipping option, def. 0                              */
/*  ) = function;                                                           */
/*                                                                          */
/*  The right-hand side must provide the definition of a function, depending*/
/*  on a single argument. The function is evaluated in the range [a,b] on   */
/*  the x axis, with increments of size d; these values must be specified   */
/*  with the parameter:                                                     */
/*                                                                          */
/*  rx = a (d) b,                                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_plotf(void)
{
    int err,r,first,n,deriv;
    double x,y,y1;
    char *p;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    deriv = 0;
    n = 5;
    p = CmdBuf + 5;
    if (*p == '1') {
        deriv = 1;
        n++;
    }
    else if (*p == '2') {
        deriv = 2;
        n++;
    }
    if (parm(CmdBuf + n,9,1))    /* get parameters */
        goto PLFFin;

    if (FNFlg == 0) {
        p_err(-1,1);
        goto PLFFin;
    }
    if (FNArgN != 1) {
        printf1("Error: there must be exactly one function argument.\n");
        goto PLFFin;
    }
    if (PMFTYP5) {
        p_err(-40,1);
        goto PLFFin;
    }
    if (FVFlg)
        printf1("Sum over %d data matrix cases.\n",NOC);

    if (deriv) {        /* need memory for gradient and hessian */
          
        if (fnd_alloc(1,deriv,FNArgN,FNPN,0))
            goto PLFFin;

        if (alloc_acu(FNArgN + 1))
            goto PLFFin;

        if (alloc_acv(FNArgN * FNArgN + 1))
            goto PLFFin;
    }
    NINTMUsed = -1;

    if (PMRXFlg == 0) {
        PMRXA = PA1[0];
        PMRXB = PA2[0];
        PMRXD = (PMRXB - PMRXA) / 100.0;
    }

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    if (PMNC != 1)
        set_clip();

    x = PMRXA;
    first = 1;
    while (x <= PMRXB + PMRXD) {

        if (x > PMRXB)
            x = PMRXB;

        FNArgVal[0] = x;

        r = get_flval(&y,FNArgN,FNArgVal,deriv,1,AcU,AcV,AcV);
        if (r) {      /* r from v_eval1() */

            printf1("Can't evaluate function or derivatives.\n");
            prn_emsg2(r);
            goto PLFFin;
        }
        if (deriv == 1)
            y = AcU[1];
        else if (deriv > 1)
            y = AcV[1];
       
        if (first) {
            ps_2dplot(x,y,0);
            y1 = y;
            first = 0;
        }
        ps_2dplot(x,y,1);

        if (PMNC)           /* update bounding box */  
            upd_bbox(0,x,y);     
    
        if (fabs(x - PMRXB) <= EPSI)
            break;

        x += PMRXD;
    }
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0 && first == 0) {
        fprintf(PSFd,"gsave\n");
        ps_2dplot(PMRXB,PA1[1],1);    
        ps_2dplot(PMRXA,PA1[1],1);    
        ps_2dplot(PMRXA,y1,1);    
        ps_lwidth(0.0);
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    fprintf(PSFd,"stroke\ngrestore\n");

    if (NINTMUsed >= 0)  
        ni_info();

    err = 0;

PLFFin:
    if (deriv)
        fnd_alloc(0,0,0,0,0);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_ploth()      Plot histogram                                          */
/*                                                                          */
/*  ploth(                                                                  */
/*      x=...,      intervals, required                                     */
/*      w=...,      variable for weights                                    */
/*      dscal=...,  scaling factor, def. 1.0                                */
/*      s=...,      s=0 (default) classes are left closed                   */
/*                  s=1 classes are left open                               */
/*      ns=...,     ns=1, don't plot vertical lines, def. ns=0              */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      gs=...,     grey scale value, def. 1 (white)                        */
/*      nc=...,     no clipping option, def. 0                              */
/*      df=...,     print table to output file                              */
/*      fmt=...,    print format for output file, def. 10.4                 */
/*  ) = variable;                                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_ploth(void)
{
    register int i,j,k;
    int err,n,ix,gs;
    double tmp,xmin,xmax,wt,wsum,scal;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,4,1))    /* get parameters */
        goto PLHFin;

    if (PMNV != 1) {
        printf1("Error: there must be exactly one variable name.\n");
        goto PLHFin;
    }
    n = PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (n < 1) {
        printf1("Error: need classes defined with x parameter.\n");
        goto PLHFin;
    }
    if (DScalFlg && DScal > 0.0)
        scal = DScal;
    else
        scal = 1.0;

    for (j = 1; j < n; ++j) {
        if (PMTP[j] <= PMTP[j - 1]) {
            printf1("Error: need non-empty classes.\n");
            goto PLHFin;
        }
    }
    ix = PMVIdx[0];
    xmin = PMTP[0];
    xmax = PMTP[n - 1];

    if (alloc_acw(n + 1))
        goto PLHFin;
    if (alloc_acy(n + 1))
        goto PLHFin;

    if (PMF1Def) {
        fprintf(PMF1d,"# Frequency distribution: %s\n",VName[ix]);
        if (PMWVar >= 0)  
            fprintf(PMF1d,"# Using weights defined by: %s\n",VName[PMWVar]);
        if (scal != 1.0)  
            fprintf(PMF1d,"# Scaling factor: %g\n",scal);
    }
    if (PMFmtF == 0)            /* default print format */
        pmfmt(10,4);

    wsum = 0.0;
    wt = 1.0;
    for (i = 0; i < NOC; ++i) {

        if (PMWVar >= 0) {
            wt = get_data(PMWVar,i);
            if (wt < 0.0) {
                printf1("Error: found negative weight in case %d\n",i + 1);
                goto PLHFin;
            }
        }
        tmp = get_data(ix,i);

        if (tmp < xmin) {
            printf1("Error: value in case %d is %g (less than %g)\n",    
                                                    i + 1,tmp,xmin);
            goto PLHFin;
        }
        if (tmp > xmax) {
            printf1("Error: value in case %d is %g (greater than %g)\n",    
                                                    i + 1,tmp,xmax);
            goto PLHFin;
        }
        k = n - 1;
        for (j = 1; j < n; ++j) {
            if (PMS == 1) {
                if (tmp <= PMTP[j]) {
                    k = j;
                    break;
                }
            }
            else {
                if (tmp < PMTP[j]) {
                    k = j;
                    break;
                }
            }
        }
        AcW[k] += wt;
        wsum += wt;
    }
    if (wsum <= 0.0) {
        printf1("Error: sum of weights is zero.\n");
        goto PLHFin;
    }   
    for (j = 1; j < n; ++j)  
        AcY[j] = scal * AcW[j] / (wsum * (PMTP[j] - PMTP[j - 1]));

    if (PMF1Def) {
        for (j = 1; j < n; ++j) {
            fprintf(PMF1d,PMFmtS,PMTP[j - 1]);
            fprintf(PMF1d,PMFmtS,PMTP[j]);
            fprintf(PMF1d,PMFmtS,AcW[j]);
            fprintf(PMF1d,PMFmtS,AcW[j] / wsum);
            fprintf(PMF1d,PMFmtS,AcY[j]);
            fprintf(PMF1d,"\n");
        }
    }
    gs = 0;
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0)  
        gs = 1;

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    if (PMNC != 1)
        set_clip();

    ps_2dplot(PMTP[0],0.0,0);
    if (PMNC)   
        upd_bbox(0,PMTP[0],0.0);     

    for (j = 1; j < n; ++j) {
        ps_2dplot(PMTP[j - 1],AcY[j],1);
        if (PMNC)   
            upd_bbox(0,PMTP[j - 1],AcY[j]);     

        ps_2dplot(PMTP[j],AcY[j],1);
        if (PMNC)   
            upd_bbox(0,PMTP[j],AcY[j]);     
    }
    ps_2dplot(PMTP[n - 1],0.0,1);
      
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");
        ps_2dplot(PMTP[0],0.0,1);    
        ps_lwidth(0.0);
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    if (PMNS != 1) {
        for (j = 1; j < n - 1; ++j) {
            ps_2dplot(PMTP[j],0.0,0);
            ps_2dplot(PMTP[j],AcY[j],1);
        }   
    }
    fprintf(PSFd,"stroke\ngrestore\n");
    err = 0;

PLHFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotd()      Plot density                                            */
/*                                                                          */
/*  plotd(                                                                  */
/*      x=...,      sequence of evaluation points, required                 */
/*      d=...,      bandwidth, def. 1.0                                     */
/*      k=...,      kernel, def. 1                                          */
/*                  1 : uniform                                             */
/*                  2 : triangle                                            */
/*                  3 : quartic                                             */
/*                  4 : Epanechnikov                                        */
/*      dscal=...,  scaling factor, def. 1.0                                */
/*      lt=...,     line type, def. 1 (solid line)                          */
/*      lw=...,     line width, def. 0.2 (mm)                               */
/*      gs=...,     grey scale value, def. 1 (white)                        */
/*      nc=...,     no clipping option, def. 0                              */
/*      df=...,     print table to output file                              */
/*      fmt=...,    print format for output file, def. 10.4                 */
/*  ) = variable;                                                           */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_plotd(void)
{
    register int i,j;
    int err,nx,ix,m,gs;
    double w,a,b,ymax,scal;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,4,1))    /* get parameters */
        goto PLDFin;

    if (PMFmtF == 0)  
        pmfmt(10,4);

    if (PMK < 2 || PMK > 4)
        PMK = 1;

    if (PMNV != 1) {
        printf1("Error: there must be exactly one variable name.\n");
        goto PLDFin;
    }
    if (DScalFlg && DScal > 0.0)
        scal = DScal;
    else
        scal = 1.0;

    nx = PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (nx < 1) {
        printf1("Error: need x parameter.\n");
        goto PLDFin;
    }
    ix = PMVIdx[0];

    if (PMD < EPSI1)
        PMD = 1.0;

    if (alloc_acx(NOC + 1))
        goto PLDFin;
    if (alloc_acu(NOC + 1))
        goto PLDFin; 
    if (alloc_acy(nx + 1))
        goto PLDFin;

    for (i = 0; i < NOC; ++i)
        AcX[i] = get_data(ix,i);

    if (sortd(NOC,AcX,0))     
        goto PLDFin;      

    if (PMF1Def) {
        fprintf(PMF1d,"# Density plot: %s\n",VName[ix]);
        fprintf(PMF1d,"# Kernel: ");
        switch (PMK) {
            case 2:  fprintf(PMF1d,"triangle\n"); break;
            case 3:  fprintf(PMF1d,"quartic\n"); break;
            case 4:  fprintf(PMF1d,"Epanechnikov\n"); break;
            default: fprintf(PMF1d,"uniform\n"); break;
        }
    }
    w = PMD / 2.0;
    ymax = 0.0;
    for (j = 0; j < nx; ++j) {
  
        a = PMTP[j] - w;
        b = PMTP[j] + w;
        m = 0;
        for (i = 0; i < NOC; ++i) {
            if (AcX[i] < a)
                continue;
            if (AcX[i] > b)
                break;

            AcU[m++] = AcX[i];
        }
        if (m > 0) {
            AcY[j] = mkern1(PMK,m,PMD,PMTP[j],AcU) / (double)NOC;
            ymax = dmax(ymax,AcY[j]);
        }
        else
            AcY[j] = 0.0;

        AcY[j] *= scal;

        if (PMF1Def) {
            fprintf(PMF1d,PMFmtS,PMTP[j]);
            fprintf(PMF1d,PMFmtS,AcY[j]);
            fprintf(PMF1d,"\n");
        }
    }
    printf1("Max height of density: %g\n",ymax);

    gs = 0;
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0)  
        gs = 1;

    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    if (PMNC != 1)
        set_clip();

    if (gs) {
        ps_2dplot(PMTP[0],0.0,0);
        ps_2dplot(PMTP[0],AcY[0],1);
    }
    else
        ps_2dplot(PMTP[0],AcY[0],0);
    if (PMNC)   
        upd_bbox(0,PMTP[0],AcY[0]);     

    for (j = 1; j < nx; ++j) {
        ps_2dplot(PMTP[j],AcY[j],1);
        if (PMNC)   
            upd_bbox(0,PMTP[j],AcY[j]);     
    }
    if (gs) {
        ps_2dplot(PMTP[nx - 1],0.0,1);
        ps_2dplot(PMTP[0],0.0,1);
        fprintf(PSFd,"gsave\n");
        ps_lwidth(0.0);
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    fprintf(PSFd,"stroke\ngrestore\n");
    err = 0;

PLDFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_ploto(typ)       Plot cirlce or ellipse                              */
/*                                                                          */
/*          typ 0 :     ploto(pos,lt,lw,gs,nc) = r [,alpha,beta];           */
/*          typ 1 :     plote(pos,lt,lw,gs,nc) = a,b[,r];                   */
/*                      Return 0 if OK, -1 if error.                        */

int pl_ploto(int typ)
{
    int err;
    double x,y,a,b,r;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,7,1))    /* get parameters */
        goto PLOFin;

    if (PMXYFlg == 0) {
        p_err(-1,1);
        goto PLOFin;
    }
    if (typ == 0) {
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
            goto PLOFin;
        }
        r = PMRHSX[0];
        if (r <= 0.0 || a < 0.0 || b <= a || b > 360.0) {
            p_err(-1,1);
            goto PLOFin;
        }
    }
    else {
        if (PMRHSN != 2 && PMRHSN != 3) {
            p_err(-1,1);
            goto PLOFin;
        }
        a = PMRHSX[0];
        b = PMRHSX[1];
        if (PMRHSN == 3)
            r = PMRHSX[2];
        else
            r = 0.0;

        if (a <= 0.0 || b <= 0.0 || r < -360.0 || r >= 360.0) {
            p_err(-1,1);
            goto PLOFin;
        }
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    if (PMNC != 1)
        set_clip();

    x = ps_2dx(PMX);      /* translate to PostScript coordinates */
    y = ps_2dy(PMY);

    if (typ == 0) {
        r = ps_2dx(r) - ps_2dx(0.0);
            
        fprintf(PSFd,"%5.2f %5.2f m\nstroke\n",x,y);
        fprintf(PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f arc\n",x,y,r,a,b);

        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            if (a > 0.0 || b < 360.0)  
                fprintf(PSFd,"%5.2f %5.2f l\nclosepath\n",x,y);
            fprintf(PSFd,"gsave\n");
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\n");
        }
        a = b = r;
    }
    else {

        a = ps_2dx(a) - ps_2dx(0.0);
        b = ps_2dy(b) - ps_2dy(0.0);
        fprintf(PSFd,"%5.2f %5.2f translate\n",x,y);
        fprintf(PSFd,"%5.2f rotate\n",r);

        if (a != 0.0)
            fprintf(PSFd," 1.00 %5.2f scale\n",b / a);
        fprintf(PSFd,"0 0 %5.2f 0 360 arc\n",a);

        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            fprintf(PSFd,"gsave\n");
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\n");
        }
    }
    fprintf(PSFd,"stroke\ngrestore\n");
    if (PMNC) {
        upd_bbox(1,x + a,y);
        upd_bbox(1,x - a,y);
        upd_bbox(1,x,y + b);
        upd_bbox(1,x,y - b);
    }
    err = 0;

PLOFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotk()          Plot arc.                                           */
/*                                                                          */
/*                      plotk(lt,lw,a) = x1,y1,x2,y2,...;                   */
/*                      a = arrow                                           */  
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotk(void)
{
    int i,n,err,first,ptyp;
    double x,y,xl,yl,tmp,tmp1,tmp2,rot,phik,phil,r;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 5,7,1))  /* get parameters */
        goto PLKFin;

    if (PMSC >= 180.0 || PMSC <= -180.0) {
        p_err(-1,1);
        goto PLKFin;
    }
    n = PMRHSN;
    if (n < 4 || (n / 2) * 2 != n) {
        p_err(-25,1);
        goto PLKFin;
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
    if (PMNC != 1)
        set_clip();

    if (fabs(PMSC) > EPSI1)
        n -= 2;
    first = 1;
    for (i = 0; i < n; i += 2) {

        x = ps_2dx(PMRHSX[i]);
        y = ps_2dy(PMRHSX[i + 1]);
        if (PMNC) 
            upd_bbox(1,x,y);              /* update bounding box */

        if (first) {
            fprintf(PSFd,"%5.2f %5.2f m\n",x,y);
            first = 0;
        }
        else
            fprintf(PSFd,"%5.2f %5.2f l\n",x,y);
    }
    if (fabs(PMSC) <= EPSI1)
        goto PLKF;
     
    /* plot the final arc */

    x  = ps_2dx(PMRHSX[n]);                           
    xl = ps_2dx(PMRHSX[n - 2]);
    y  = ps_2dy(PMRHSX[n + 1]);                           
    yl = ps_2dy(PMRHSX[n - 1]);
    tmp1 = x - xl;
    tmp2 = y - yl;
    if (PMNC) 
        upd_bbox(1,x,y);              /* update bounding box */

    if (fabs(tmp1) > EPSI1 || fabs(tmp2) > EPSI1) {

        if (fabs(tmp1) > EPSI1)
            tmp = atan(tmp2 / tmp1);
        else
            tmp = Pi / 2.0;

        if (PMSC < 0.0) {
            PMSC = -PMSC;
            ptyp = 1;
        }
        else
            ptyp = 0;

        rot = tmp * 180 / Pi - PMSC;
        if (tmp1 > 0.0 || (tmp1 == 0.0 && tmp2 > 0.0)) 
            tmp += 0.5 * Pi;
        else {
            tmp += 1.5 * Pi;
            rot += 180;
        }
        phik = tmp - PMSC * Pi / 180.0;
        phil = tmp + PMSC * Pi / 180.0;
        r = 0.0;
        tmp = cos(phik) - cos(phil);
        if (fabs(tmp) > EPSI1)  
            r = tmp1 / tmp;

        if (fabs(r) <= EPSI1) {
            tmp = sin(phik) - sin(phil);
            if (fabs(tmp) > EPSI1)  
                r = tmp2 / tmp;
        }
        if (fabs(r) > EPSI1) {

            if (r < 0.0)  
                r = -r;
                           
            tmp1 = x - r * cos(phik);
            tmp2 = y - r * sin(phik);
            phik *= 180.0 / Pi;
            phil *= 180.0 / Pi;

            if (ptyp == 1) {
                tmp1 = x + xl - tmp1;
                tmp2 = y + yl - tmp2;
                phik += 180.0;
                phil += 180.0;
                rot += 2.0 * PMSC;
            }
            fprintf(PSFd,"stroke\n%5.2f %5.2f ",tmp1,tmp2);
            fprintf(PSFd,"%5.2f %5.2f %5.2f arc\n",r,phik,phil);

            if (PMAFlg && PMA1 > 0.0 && PMA2 > 0.0) {       /* plot arrow */
                fprintf(PSFd,"%5.2f %5.2f m\ncurrentpoint\nstroke\nm\n",x,y);
                fprintf(PSFd,"%5.2f rotate\n%5.2f %5.2f scale\n",rot,PMA1,PMA2);
                fprintf(PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\n");
                fprintf(PSFd,"closepath\nfill\n");
            }
        }
    }
PLKF:
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
    err = 0;

PLKFin:
    p_clean();
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  pl_plotch()         Plot convex hull.                                   */
/*                                                                          */
/*                      plotch(lt,lw,gs,s,fs,sel,nc,ns,ic) = X,Y;           */
/*                      if ns >= 2 use Akima for smoothing                  */
/*                      if ic > 0 expand convex hull                        */
/*                                                                          */
/*                      Correced July 2001: ns option only if n > 2.        */
/*                                                                          */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotch(void)
{
    register int i,j,ik,l;
    int ix,iy,m,n,nn,r,err,ns;
    double x,y,dx,dy,xm,ym,xa,xb,ya,yb;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (parm(CmdBuf + 6,4,1))    /* get parameters */
        goto PLHFin;

    if (NOC < 2)
        goto PLHFin;

    if (PMNV != 2) {
        p_err(-1,1);
        goto PLHFin;
    }
    ix = PMVIdx[0];
    iy = PMVIdx[1];

    if (alloc_acxf(NOC + 1))
        goto PLHFin;
    if (alloc_acyf(NOC + 1))
        goto PLHFin;
    if (alloc_acn(NOC + 1))
        goto PLHFin;
    if (alloc_aci(NOC + 1))
        goto PLHFin;
    if (alloc_acj(NOC + 1))
        goto PLHFin;

    nn = 0;
    for (i = 0; i < NOC; ++i) {
        if (eval_sve(i) == 0)   
            continue;
        nn++;
        AcXF[nn] = (float)get_data(ix,i);       
        AcYF[nn] = (float)get_data(iy,i);       
        AcN[nn] = nn;
    }
    n = g_chull(AcXF,AcYF,nn,AcN,AcJ,AcI);   
    if (n <= 0) {
        if (n < 0)
            p_err(-2,1);
        goto PLHFin;
    }
    if (PMNS >= 2 && n < 3)  
        PMNS = 1;

    if (PMNS >= 2 && n > 2) {        /* use Akima algorithm for smoothing */

        ns = n + 4;
        if (alloc_acx(ns + 1))
            goto PLHFin;
        if (alloc_acy(ns + 1))
            goto PLHFin;

        xm = ym = 0.0;

        ik = AcI[1];                        /* get points of convex hull */
        for (i = 1; i <= n; ++i) {
            j = AcJ[ik];
            AcX[i + 2] = (double)AcXF[j];
            AcY[i + 2] = (double)AcYF[j];
            ik = AcI[ik];
        }
        AcX[n + 3] = AcX[3];
        AcY[n + 3] = AcY[3];
        AcX[n + 4] = AcX[4];
        AcY[n + 4] = AcY[4];

        AcX[1] = AcX[n + 1];
        AcY[1] = AcY[n + 1];
        AcX[2] = AcX[n + 2];
        AcY[2] = AcY[n + 2];
             
        if (PMIC > 0.0) {
            for (i = 1; i <= ns; ++i) {
                xm += AcX[i];
                ym += AcY[i];
            }
            xm /= (double)ns;
            ym /= (double)ns;

            xa = ya = xb = yb = 0.0;
            for (i = 1; i <= ns; ++i) {
                xa = dmax(xa,xm - AcX[i]);
                xb = dmax(xb,AcX[i] - xm);
                ya = dmax(ya,ym - AcY[i]);
                yb = dmax(yb,AcY[i] - ym);
            }
            dx = PMIC * UXLen / PXLen;
            dy = PMIC * UYLen / PYLen;

            for (i = 1; i <= ns; ++i) {
                x = AcX[i];
                y = AcY[i];
                if (x < xm - EPSI1) {
                    x -= dx * (xm - x) / xa;
                    if (y > ym + EPSI1)  
                        y += dy * (y - ym) / yb;
                    else if (y < ym - EPSI1)
                        y -= dy * (ym - y) / ya;
                }
                else if (x > xm + EPSI1) {
                    x += dx * (x - xm) / xb;
                    if (y > ym + EPSI1)
                        y += dy * (y - ym) / yb;
                    else if (y < ym - EPSI1)
                        y -= dy * (ym - y) / ya;
                }
                AcX[i] = x;
                AcY[i] = y;
            }
        }
        m = (ns - 1) * PMNS + 1;

        if (alloc_acu(m + 1))
            goto PLHFin;
        if (alloc_acv(m + 1))
            goto PLHFin;

        r = splakima(2,ns,AcX,AcY,PMNS,m,AcU,AcV);

        if (r) {
            printf1("Error (%d) in smoothing algorithm.\n",r);
            PMNS = -1;      /* continue with standard method */
        }
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_lwidth(PMLW);
    ps_ltyp(PMLT);
   
    if (PMNS >= 2) {

        l = 1;
        for (i = 2 * PMNS + 1; i <= m - PMNS; ++i) {
            if (l) {
                ps_2dplot(AcU[i],AcV[i],0);
                l = 0;
            }
            else
                ps_2dplot(AcU[i],AcV[i],1);

            if (PMNC)
                upd_bbox(0,AcU[i],AcV[i]);
        }
    }
    else {

        ik = AcI[1];
        l = 0;
        for (i = 1; i <= n; ++i) {
            j = AcJ[ik];
            x = (double)AcXF[j];
            y = (double)AcYF[j];
            ps_2dplot(x,y,l); 
            if (PMNC)
                upd_bbox(0,x,y);
            l = 1;
            ik = AcI[ik];
        }
        ik = AcI[1];
        j = AcJ[ik];
        ps_2dplot((double)AcXF[j],(double)AcYF[j],1);     /* close path */
    }
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");       
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");   
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");

    if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

        for (i = 0; i < NOC; ++i) {
            if (eval_sve(i) == 0)   
                continue;

            x = ps_2dx(get_data(ix,i));       
            y = ps_2dy(get_data(iy,i));       
            ps_sym(PMS,x,y,PtMM * PMFS / 2.0);
        }
    }
    fprintf(PSFd,"grestore\n");
    err = 0;

PLHFin:
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plots(typ)       Smoothing with Akima algorithm.                     */
/*                                                                          */
/*              typ 0:  plotsp(ns,lt,lw,gs,nc,sc,s,fs) = x1,y1,x2,y2,...;   */
/*              typ 1:  plots (ns,lt,lw,gs,nc,sc,sel,s,fs) = VX,VY;         */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plots(int typ)
{
    register int i,j;
    int ix,iy,m,n,nn,r,err,first;
    double x,y,xa,ya;
           
    err = -1;
    if (check_pcmd(1,2))
        return(-1);

    if (typ == 0) {
        if (parm(CmdBuf + 6,7,1))    /* get parameters */
            goto PLSFin;

        n = PMRHSN / 2;
        if (n < 2 || 2 * n != PMRHSN) {
            p_err(-25,1);
            goto PLSFin;
        }
    }
    else {
        if (parm(CmdBuf + 5,4,1))    /* get parameters */
            goto PLSFin;

        if (NOC < 2)
            goto PLSFin;

        n = NOC;
        if (PMNV != 2) {
            p_err(-1,1);
            goto PLSFin;
        }
        ix = PMVIdx[0];
        iy = PMVIdx[1];
    }
    if (PMNS < 2)                   /* number of subintervals */
        PMNS = 2;

    if (alloc_acx(n + 2))
        goto PLSFin;
    if (alloc_acy(n + 2))
        goto PLSFin;

    nn = j = 0;                     /* get data points */
    xa = ya = 0.0;
    for (i = 0; i < n; ++i) {

        if (typ == 0) {
            x = PMRHSX[j++];
            y = PMRHSX[j++];
        }
        else {
            if (eval_sve(i) == 0)  
                continue;
               
            x = get_data(ix,i);       
            y = get_data(iy,i);       
        }
        if (nn == 0 || x != xa || y != ya) {
            nn++;
            AcX[nn] = x;
            AcY[nn] = y;
        }   
        xa = x;
        ya = y;
    }
    if (nn < 2) {
        p_err(-26,1);
        goto PLSFin;
    }
    if ((int)PMSC == 1) {           /* closed curve */

        if (AcX[nn] != AcX[1] || AcY[nn] != AcY[1]) {
            nn++;
            AcX[nn] = AcX[1];
            AcY[nn] = AcY[1];
        }
    }
    m = (nn - 1) * PMNS + 1;

    if (alloc_acu(m + 1))
        goto PLSFin;
    if (alloc_acv(m + 1))
        goto PLSFin;

    r = splakima(2,nn,AcX,AcY,PMNS,m,AcU,AcV);
    if (r) {
        printf1("Error in smoothing algorithm.\n");
        goto PLSFin;
    }
    if (PMF1Def) {      /* write into output file */
        for (i = 1; i <= m; ++i) {
            fprintf(PMF1d,"%6d ",i);
            fprintf(PMF1d,PMFmtS,AcU[i]);
            fprintf(PMF1d,PMFmtS,AcV[i]);
            fprintf(PMF1d,"\n");
        }
    }
    fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);

    fprintf(PSFd,"gsave\n");
    if (PMNC != 1)
        set_clip();
    ps_lwidth(PMLW);
    ps_ltyp(PMLT);
   
    first = 1;
    for (i = 1; i <= m; ++i) {
        if (first) {
            ps_2dplot(AcU[i],AcV[i],0);
            first = 0;
        }
        else
            ps_2dplot(AcU[i],AcV[i],1);

        if (PMNC)
            upd_bbox(0,AcU[i],AcV[i]);
    }

    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0 && first == 0) {
        if ((int)PMSC == 1) {
            fprintf(PSFd,"gsave\n");       
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\n");   
        }
        else if (first == 0) {
            fprintf(PSFd,"gsave\n");
            ps_2dplot(AcU[m],PA1[1],1);    
            ps_2dplot(AcU[1],PA1[1],1);    
            ps_2dplot(AcU[1],AcV[1],1);    
            ps_lwidth(0.0);
            ps_fill(PMGS);
            fprintf(PSFd,"grestore\n");
        }
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");

    if (PMS >= 1 && PMS <= 17 && PMFS > 0.0) {  /* plot symbols */

        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

        ps_sym(PMS,ps_2dx(AcU[1]),ps_2dy(AcV[1]),PtMM * PMFS / 2.0);
        ps_sym(PMS,ps_2dx(AcU[m]),ps_2dy(AcV[m]),PtMM * PMFS / 2.0);
    }
    fprintf(PSFd,"grestore\n");
    err = 0;

PLSFin:
    p_clean();
    return(err);
}

