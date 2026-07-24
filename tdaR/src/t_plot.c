#include "tda_rhooks.h"
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
#include "tda_context.h"

/*  functions in t_plot.c */

int check_pcmd(TDAContext *ctx, int opt,int dim);              
int check_ps(TDAContext *ctx, int dim);              
void upd_bbox(TDAContext *ctx, int opt,double x,double y); 
void set_clip(TDAContext *ctx); 
void ps_ltyp(TDAContext *ctx, int typ);
void ps_lwidth(TDAContext *ctx, double lw);
void ps_setfs(TDAContext *ctx, double fs);
void ps_fill(TDAContext *ctx, double g); 
int ps_cclip(TDAContext *ctx, double x,double y);
void ps_2dplot(TDAContext *ctx, double x,double y, int opt);
double ps_2dx(TDAContext *ctx, double x);
double ps_2dy(TDAContext *ctx, double y);
void ps_sym(TDAContext *ctx, int typ, double px, double py,double siz);
void plot_str(TDAContext *ctx, double x,double y,char *s, double siz, int opt,int adj, int r,int xopt,int cflag);
void ps_nlab(TDAContext *ctx, int n,double x,double y,double fs,int opt);
int ps_grid(TDAContext *ctx, int typ);
int pl_label(TDAContext *ctx, int typ);
int pl_text(TDAContext *ctx);
int pl_frame(TDAContext *ctx);
int pl_rec(TDAContext *ctx);
int pl_plotp(TDAContext *ctx, int typ);
int pl_plotm(TDAContext *ctx);
int pl_plotf(TDAContext *ctx);
int pl_ploth(TDAContext *ctx);
int pl_plotd(TDAContext *ctx);
int pl_ploto(TDAContext *ctx, int typ);
int pl_plotk(TDAContext *ctx); 
int pl_plotch(TDAContext *ctx);
int pl_plots(TDAContext *ctx, int typ);

/*--------------------------------------------------------------------------*/

int BBLX,BBLY,BBUX,BBUY;    /* bounding box                                 */









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

int check_pcmd(TDAContext *ctx, int opt,int dim)       
{
    int err;
    err = check_cmd(ctx, opt);
    if (err == 0)  
        err = check_ps(ctx, dim);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_ps(dim)       Check for valid PostScript output file.             */
/*                      Return 0 if OK, -1 if error.                        */

int check_ps(TDAContext *ctx, int dim)            
{
    if (ctx->PSFFlg != 2) {
        printf1(ctx, "Error: need PostScript output file and coordinates.\n");
        return(-1);
    }
    if (dim == 2 && ctx->PS3DFlg) {
        printf1(ctx, "Error: current coordinate system is 3-dimensional.\n");
        return(-1);
    }
    else if (dim == 3 && ctx->PS3DFlg == 0) {
        printf1(ctx, "Error: current coordinate system is 2-dimensional.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  upd_bbox(opt,x,y)   Update bounding box for point (x,y). If opt == 0    */
/*                      user coordinates, otherwise PostScript coordinates. */

void upd_bbox(TDAContext *ctx, int opt,double x,double y)
{
    int ix,iy;

    if (opt == 0) {
        ix = (int)ps_2dx(ctx, x) + ctx->XOrg;
        iy = (int)ps_2dy(ctx, y) + ctx->YOrg;
    }
    else {
        ix = (int)x + ctx->XOrg;
        iy = (int)y + ctx->YOrg;
    }
    if (ctx->BBLY > iy)
        ctx->BBLY = iy;
    if (ctx->BBUY < iy)
        ctx->BBUY = iy;
    if (ctx->BBLX > ix)
        ctx->BBLX = ix;
    if (ctx->BBUX < ix)
        ctx->BBUX = ix;
}

/* ------------------------------------------------------------------------ */
/*  set_clip()      Set default clipping region. Observe LWCSMax            */

void set_clip(TDAContext *ctx)
{
    double x;  
    x = ctx->LWCSMax * ctx->PtMM;
    fprintf(ctx->PSFd,"%5.2f %5.2f m\n%5.2f %5.2f l\n%5.2f %5.2f l\n%5.2f %5.2f l\nclosepath\nclip\nnewpath\n",
                                  x,x,x,ctx->PSYLen,ctx->PSXLen,ctx->PSYLen,ctx->PSXLen,x);
}

/* ------------------------------------------------------------------------ */
/*  ps_ltyp(typ)  Select line type for PostScript output.                   */
/*      Note: only some standard predefined formats are used here.          */

void ps_ltyp(TDAContext *ctx, int typ)
{
    switch (typ) {
        case  2:  fprintf(ctx->PSFd,"[1 2]"); break;
        case  3:  fprintf(ctx->PSFd,"[1 3]"); break;
        case  4:  fprintf(ctx->PSFd,"[1 5]"); break;
        case  5:  fprintf(ctx->PSFd,"[2 2]"); break;
        case  6:  fprintf(ctx->PSFd,"[4 2]"); break;
        case  7:  fprintf(ctx->PSFd,"[6 4]"); break;
        case  8:  fprintf(ctx->PSFd,"[6 2]"); break;
        case  9:  fprintf(ctx->PSFd,"[9 4]"); break;
        default:  fprintf(ctx->PSFd,"[]"); break;
    }
    fprintf(ctx->PSFd," 0 setdash\n");
}

/* ------------------------------------------------------------------------ */
/*  ps_lwidth(lw)   Set line width to lw.                                   */

void ps_lwidth(TDAContext *ctx, double lw)
{
    fprintf(ctx->PSFd,"%6.4f setlinewidth\n",ctx->PtMM * lw);
}

/* ------------------------------------------------------------------------ */
/*  ps_setfs(fs)    Set font size.                                          */

void ps_setfs(TDAContext *ctx, double fs)
{
    fprintf(ctx->PSFd,"/fsiz %7.4f def FT (",fs); 
}

/* ------------------------------------------------------------------------ */
/*  ps_fill     PostScript fill with greyscale value command.               */

void ps_fill(TDAContext *ctx, double g)
{
    if (g >= 0.0 && g <= 1.0)  
        fprintf(ctx->PSFd," %6.4f setgray fill ",g);
}

/* ------------------------------------------------------------------------ */
/*  ps_cclip(x,y)   x,y in PostScript coordinates. Return 1 if (x,y) is     */
/*                  inside coordinate system, otherwise 0.                  */

int ps_cclip(TDAContext *ctx, double x,double y)
{
    if (x < (double)ctx->XOrg || x > (double)ctx->XOrg + ctx->PSXLen)
        return(0);
    if (y < (double)ctx->YOrg || y > (double)ctx->YOrg + ctx->PSYLen)
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dplot   PostScript: 2D moveto or lineto, with translation           */
/*              opt == 0 moveto, opt = 1 lineto, opt = 2 translate          */
/*              opt == 3 only values.                                       */

void ps_2dplot(TDAContext *ctx, double x,double y, int opt)
{
    double px,py;

    if (!ctx->PSLog[0])
        px = (x - ctx->PA1[0]) / ctx->UXLen;
    else if (x <= 0.0)
        px = 0.0;
    else
        px = rlog(ctx, x / ctx->PA1[0]) / rlog(ctx, ctx->PA2[0] / ctx->PA1[0]);

    if (!ctx->PSLog[1])
        py = (y - ctx->PA1[1]) / ctx->UYLen;
    else if (y <= 0.0)
        py = 0.0;
    else
        py = rlog(ctx, y / ctx->PA1[1]) / rlog(ctx, ctx->PA2[1] / ctx->PA1[1]);

    px *= ctx->PtMM * ctx->PXLen;
    py *= ctx->PtMM * ctx->PYLen;

    if (!opt)  
        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",px,py);
    else if (opt == 1)  
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",px,py);
    else if (opt == 2)
        fprintf(ctx->PSFd,"%5.2f %5.2f translate\n",px,py);
    else 
        fprintf(ctx->PSFd,"%5.2f %5.2f ",px,py);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dx(x)   return x in PostScript coordinates.                         */

double ps_2dx(TDAContext *ctx, double x)
{
    double px;

    if (!ctx->PSLog[0])
        px = (x - ctx->PA1[0]) / ctx->UXLen;
    else if (x <= 0.0)
        px = 0.0;
    else
        px = rlog(ctx, x / ctx->PA1[0]) / rlog(ctx, ctx->PA2[0] / ctx->PA1[0]);

    px *= ctx->PtMM * ctx->PXLen;
    return(px);
}

/* ------------------------------------------------------------------------ */
/*  ps_2dy(y)   return y in PostScript coordinates.                         */

double ps_2dy(TDAContext *ctx, double y)
{
    double py;

    if (!ctx->PSLog[1])
        py = (y - ctx->PA1[1]) / ctx->UYLen;
    else if (y <= 0.0)
        py = 0.0;
    else
        py = rlog(ctx, y / ctx->PA1[1]) / rlog(ctx, ctx->PA2[1] / ctx->PA1[1]);

    py *= ctx->PtMM * ctx->PYLen;
    return(py);
}

/* ------------------------------------------------------------------------ */
/*  ps_sym(typ,x,y,siz)     plot symbol (typ) at PostScript coordinates x,y */

void ps_sym(TDAContext *ctx, int typ, double px, double py,double siz)
{
            
    if (typ >= 1 && typ <= 17)
        fprintf(ctx->PSFd,"%%#symbol: %d %5.2f %5.2f %5.2f\n",typ,px,py,siz);
            
    switch (typ) {
        case  1:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m cross\n",px,py,px,py);
                  break;
        case  2:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m xsym\n",px,py,px,py);
                  break;
        case  3:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f m cross xsym\n",px,py,px,py);
                  break;
        case  4:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle stroke\n",px,py,px,py);
                  break;
        case  5:  fprintf(ctx->PSFd,"%5.2f %5.2f circle fill\n",px,py);
                  break;
        case  6:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle ",px,py,px,py);
                  fprintf(ctx->PSFd,"stroke %5.2f %5.2f m cross\n",px,py);
                  break;
        case  7:  fprintf(ctx->PSFd,"%5.2f %5.2f circle invers %5.2f %5.2f circle ",px,py,px,py);
                  fprintf(ctx->PSFd,"stroke %5.2f %5.2f m xsym\n",px,py);
                  break;
        case  8:  fprintf(ctx->PSFd,"%5.2f %5.2f m square stroke grestore\n",px,py);
                  break;
        case  9:  fprintf(ctx->PSFd,"%5.2f %5.2f m square fill grestore\n",px,py);
                  break;
        case 10:  fprintf(ctx->PSFd,"%5.2f %5.2f m trian1 stroke grestore\n",px,py);
                  break;
        case 11:  fprintf(ctx->PSFd,"%5.2f %5.2f m trian1 fill grestore\n",px,py);
                  break;
        case 12:  fprintf(ctx->PSFd,"%5.2f %5.2f m trian2 stroke grestore\n",px,py);
                  break;
        case 13:  fprintf(ctx->PSFd,"%5.2f %5.2f m trian2 fill grestore\n",px,py);
                  break;
        case 14:  fprintf(ctx->PSFd,"%5.2f %5.2f m square stroke grestore cross\n",px,py);
                  break;
        case 15:  fprintf(ctx->PSFd,"%5.2f %5.2f m rhomb stroke grestore\n",px,py);
                  break;
        case 16:  fprintf(ctx->PSFd,"%5.2f %5.2f m rhomb fill grestore\n",px,py);
                  break;
        case 17:  fprintf(ctx->PSFd,"%5.2f %5.2f m rhomb stroke grestore cross\n",px,py);
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

void plot_str(TDAContext *ctx, double x,double y,char *s, double siz, int opt,int adj,int r, int xopt,int cflag)
{
    register int sflg,scnt = 0;
    register char c, *p;

    if (xopt == 0) {
        x = ps_2dx(ctx, x);
        y = ps_2dy(ctx, y);
    }
    fprintf(ctx->PSFd,"%%#text: %5.2f %5.2f %5.2f %d %d %s\n",x,y,siz,adj,r,s);

    fprintf(ctx->PSFd,"gsave\n");
    fprintf(ctx->PSFd,"%5.2f %5.2f translate\n0 0 m\n",x,y);

    if (r != 0)
        fprintf(ctx->PSFd,"%d rotate\n",r);

    p = s;
    sflg = 0;
    if (*p != '@') {
        fprintf(ctx->PSFd,"/fsiz %5.2f def FT (",siz);  /* default is TFont */
        sflg = 1;
    }
    while (*p) {
        c = *p;
        if (c == '@' && isdigit((int) *(p+1))) {  /* switch to symbol font */
                          
            if (sflg == 1)
                fprintf(ctx->PSFd,")\nshow\n"); 
            fprintf(ctx->PSFd,"/fsiz %5.2f def FS (\\",siz); 
            sflg = 2;
            scnt = 0;
        }
        else {
            if (sflg == 2) {
                if (scnt >= 3 || isdigit((int)c) == 0) {
                    fprintf(ctx->PSFd,")\nshow\n/fsiz %5.2f def FT (",siz); 
                    sflg = 1;
                    scnt = 0;
                }
                else
                    scnt++;
            }
            if (c == '(' || c == ')')  
                fprintf(ctx->PSFd,"\\");
            fprintf(ctx->PSFd,"%c",c);
        }
        p++;
    }
    fprintf(ctx->PSFd,") ");
    if (opt == 0) {
        if (adj == 1)
            fprintf(ctx->PSFd,"center\n");
        else if (adj == 2)
            fprintf(ctx->PSFd,"aright\n");
        if (cflag)
            fprintf(ctx->PSFd,"sclear ");
        fprintf(ctx->PSFd,"show\n");
    }
    else
        fprintf(ctx->PSFd,"\n");
    fprintf(ctx->PSFd,"grestore\n");
}

/* ------------------------------------------------------------------------ */
/*  ps_nlab(n,x,y,fs,opt)                                                   */
/*                                                                          */
/*  Plot numerical label n at (x,y) with font size fs                       */
/*  opt = 0 : centered for x and y direction                                */  
/*  opt = 1 : centered and white rectangle around character                 */
/*  opt = 2 : right justified                                               */

void ps_nlab(TDAContext *ctx, int n,double x,double y,double fs,int opt)
{
    double gx,gy,ax,ay;
    int adj;
    char buf[20];

    ax = fs * ctx->UXLen / ctx->PXLen;
    gx = 0.6 * ax;                   
    gy = 0.9 * fs * ctx->UYLen / ctx->PYLen;
    if (n > 9)
        gx += 0.7 * ax;
    if (n > 99)
        gx += 0.7 * ax;
    gx /= 1.9;
    gy /= 1.9;
    ax = 1.6 * gx;
    ay = 1.6 * gy;

    if (opt == 1) {
        fprintf(ctx->PSFd,"gsave\n");       
        ps_lwidth(ctx, 0.0);
        ps_2dplot(ctx, x - ax,y - ay,0);    
        ps_2dplot(ctx, x - ax,y + ay,1);    
        ps_2dplot(ctx, x + ax,y + ay,1);    
        ps_2dplot(ctx, x + ax,y - ay,1);    
        ps_fill(ctx, 1.0);
        fprintf(ctx->PSFd,"grestore\n");   
    }
    snprintf(buf,sizeof(buf),"%d",n);
    if (opt <= 1)
        adj = 1;
    else
        adj = 2;
             
    plot_str(ctx, x,y - gy,buf,1.5 * fs * ctx->PtMM,0,adj,0,0,0);
}

/*--------------------------------------------------------------------------*/
/*  ps_grid(typ)        plot grid lines. Syntax:                            */
/*                      plxgrid(lt=,lw=) = x1,x2,...                        */
/*                      typ: 0 = x-axis, 1 = y-axis.                        */
/*                      Return 0 if OK, -1 if error.                        */

int ps_grid(TDAContext *ctx, int typ)
{
    int i,err;
    double x;

    err = -1;

    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,7,1))    /* get parameters */
        goto PSGFin;

    if (ctx->PMRHSN > 0) {
        if (ctx->PMLTFlg == 0)
            ctx->PMLT = 3;           /* default line type */

        if (ctx->PMLWFlg == 0)
            ctx->PMLW = 0.05;        /* default line width */

        fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

        ps_ltyp(ctx, ctx->PMLT);
        ps_lwidth(ctx, ctx->PMLW);
           
        for (i = 0; i < ctx->PMRHSN; ++i) {
            x = ctx->PMRHSX[i];
            if (typ == 0) {
                if (x < ctx->PA1[1] || x > ctx->PA2[1])
                    continue;

                ps_2dplot(ctx, ctx->PA1[0],x,0);
                ps_2dplot(ctx, ctx->PA2[0],x,1);
            }
            else {
                if (x < ctx->PA1[0] || x > ctx->PA2[0])
                    continue;

                ps_2dplot(ctx, x,ctx->PA1[1],0);
                ps_2dplot(ctx, x,ctx->PA2[1],1);
            }
            fprintf(ctx->PSFd,"stroke\n");
        }
    }
    err = 0;

PSGFin:
    p_clean(ctx);
    return(err);
}
   
/*--------------------------------------------------------------------------*/
/*  pl_label(typ)       plot label.                                         */
/*                                                                          */
/*              typ 0   plabel(sc=,fs=) = string;                           */
/*              typ 1   pxlabel(sc=,fs=) = string;                          */
/*              typ 2   pylabel(sc=,fs=) = string;                          */
/*                      Return 0 if OK, -1 if error.                        */

int pl_label(TDAContext *ctx, int typ)
{
    int r,err = -1;
    double x,y,d;

    if (check_pcmd(ctx, 1,2))
        return(-1);

    r = 6;
    if (typ)
        r++;

    if (parm(ctx, ctx->CmdBuf + r,8,1))    /* get parameters */
         goto PLLFin;

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    if (typ == 0) {

        d = 3.0 + ctx->PMSC;
        if (ctx->PMFSFlg == 0)
            ctx->PMFS = 3.0;

        x = ctx->PSXLen / 2.0;
        y = ctx->PSYLen + d * ctx->PtMM;         
        plot_str(ctx, x,y,ctx->PMRHSTR,1.5 * ctx->PMFS * ctx->PtMM,0,1,0,1,0);
        y += 1.5 * ctx->PMFS * ctx->PtMM;
        upd_bbox(ctx, 1,x,y);
    }
    else if (typ == 1) {    /* x axis */

        if (ctx->PMFSFlg == 0)
            ctx->PMFS = 2.4;

        d = 4.8 + ctx->PMSC + 1.6 * ctx->PMFS;

        x = ctx->PSXLen / 2.0;
        y = -d * ctx->PtMM;
        plot_str(ctx, x,y,ctx->PMRHSTR,1.5 * ctx->PMFS * ctx->PtMM,0,1,0,1,0);
        y -= 0.5 * ctx->PMFS * ctx->PtMM;
        upd_bbox(ctx, 1,x,y);
    }
    else if (typ == 2) {        /* y axis */

        d = 7.0 + ctx->PMSC;
        if (ctx->PMFSFlg == 0)
            ctx->PMFS = 2.4;

        x = -d * ctx->PtMM;           
        y = ctx->PSYLen / 2.0;
        plot_str(ctx, x,y,ctx->PMRHSTR,1.5 * ctx->PMFS * ctx->PtMM,0,1,90,1,0);
        x -= 1.5 * ctx->PMFS * ctx->PtMM;
        upd_bbox(ctx, 1,x,y);
    }
    err = 0;

PLLFin:
    p_clean(ctx);
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

int pl_text(TDAContext *ctx)
{
    int cen,len,err = -1;
    double x,y;

    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,8,1))    /* get parameters */
        goto PLTFin;

    if (ctx->PMXYFlg == 0) {
        p_err(ctx, -1,1);
        goto PLTFin;
    }
    if (ctx->PMWF != 0)
        ctx->PMWF = 1;

    if (ctx->PMR)  
        cen = 0;
    else  
        cen = (int)ctx->PMSC;
       
    x = ps_2dx(ctx, ctx->PMX);
    y = ps_2dy(ctx, ctx->PMY);

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    if (ctx->PMS) {
        ctx->PMR = cen = 0;
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);
        ps_sym(ctx, ctx->PMS,x - 1.5 * ctx->PtMM * ctx->PMFS,y + 0.5 * ctx->PtMM * ctx->PMFS,ctx->PtMM * ctx->PMFS / 2.0);
    }
    plot_str(ctx, ctx->PMX,ctx->PMY,ctx->PMRHSTR,1.5 * ctx->PMFS * ctx->PtMM,0,cen,ctx->PMR,0,ctx->PMWF);

    /* update bounding box */

    len = (int)(strlen(ctx->PMRHSTR) + 1);
    upd_bbox(ctx, 1,x,y);     
    upd_bbox(ctx, 1,x,y + 1.5 * ctx->PMFS * ctx->PtMM);
    upd_bbox(ctx, 1,x,y - 0.8 * ctx->PMFS * ctx->PtMM);

    if (cen == 0) {
        upd_bbox(ctx, 1,x - 1.5 * ctx->PMFS * ctx->PtMM,y);
        upd_bbox(ctx, 1,x + len * ctx->PMFS * ctx->PtMM,y);
    }
    else {
        len = (int)(len / 2.0);
        upd_bbox(ctx, 1,x + len * ctx->PMFS * ctx->PtMM,y);
        upd_bbox(ctx, 1,x - len * ctx->PMFS * ctx->PtMM,y);
    }
    err = 0;

PLTFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_frame()          plot frame.                                         */
/*                      plframe(lt,lw,gs)                                   */
/*                      Return 0 if OK, -1 if error.                        */

int pl_frame(TDAContext *ctx)
{
    int err = -1;

    if (check_pcmd(ctx, 1,0))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 7,0,0))    /* get parameters */
        goto PLFRFin;

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);


    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave %4.2f setgray\n",ctx->PMGS);    
        fprintf(ctx->PSFd,"0 0 m\n");
        fprintf(ctx->PSFd,"%5.2f 0 l\n",ctx->PSXLen);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",ctx->PSXLen,ctx->PSYLen);
        fprintf(ctx->PSFd,"0 %5.2f l\n",ctx->PSYLen);
        fprintf(ctx->PSFd,"0 0 l\n");
        fprintf(ctx->PSFd,"fill\ngrestore\n");
    }
    if (ctx->PMLW > 0.0) {
        fprintf(ctx->PSFd,"0 0 m\n");
        fprintf(ctx->PSFd,"%5.2f 0 l\n",ctx->PSXLen);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",ctx->PSXLen,ctx->PSYLen);
        fprintf(ctx->PSFd,"0 %5.2f l\n",ctx->PSYLen);
        fprintf(ctx->PSFd,"0 0 l\n");
        fprintf(ctx->PSFd,"closepath\nstroke\n");
    }
    fprintf(ctx->PSFd,"grestore\n");
   
    err = 0;

PLFRFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_rec()            plot rectangle. Syntax:                             */
/*                      plrec(lt,lw,gs,r) = x,y,xlen,ylen;                  */
/*                      gs = grey value, r = rotate                         */
/*                      Return 0 if OK, -1 if error.                        */

int pl_rec(TDAContext *ctx)
{
    int err = -1;
    double x,y,x1,y1;

    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,7,1))    /* get parameters */
        goto PLRFin;

    if (ctx->PMRHSN != 4) {
        p_err(ctx, -1,1);
        goto PLRFin;
    }
    x = ps_2dx(ctx, ctx->PMRHSX[0]);
    y = ps_2dy(ctx, ctx->PMRHSX[1]);
    x1 = ps_2dx(ctx, ctx->PMRHSX[0] + ctx->PMRHSX[2]);
    y1 = ps_2dy(ctx, ctx->PMRHSX[1] + ctx->PMRHSX[3]);
  
    upd_bbox(ctx, 1,x,y);        /* update bounding box */
    upd_bbox(ctx, 1,x1,y1);

    if (fabs(ctx->PMRHSX[2]) < ctx->EPSI || fabs(ctx->PMRHSX[3]) < ctx->EPSI) {
        err = 0;
        goto PLRFin;
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);

    fprintf(ctx->PSFd,"%5.2f %5.2f translate\n",x,y);
    x1 -= x;
    y1 -= y;
    x = y = 0;

    if (ctx->PMR != 0 && ctx->PMR < 360 && ctx->PMR > -360)
        fprintf(ctx->PSFd,"%d rotate\n",ctx->PMR);

    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n%4.2f setgray\n",ctx->PMGS);    
        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x1,y);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x1,y1);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y1);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
        fprintf(ctx->PSFd,"fill\ngrestore\n");
    }
    if (ctx->PMLW > 0.0) {
        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x1,y);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x1,y1);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y1);
        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
        fprintf(ctx->PSFd,"closepath\nstroke\n");
    }
    fprintf(ctx->PSFd,"grestore\n");
   
    err = 0;

PLRFin:
    p_clean(ctx);
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

int pl_plotp(TDAContext *ctx, int typ)
{
    int i,ix = 0,iy = 0,n,nn,err,first;
    double x,y,xa,ya;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (typ == 0) {
        if (parm(ctx, ctx->CmdBuf + 5,7,1))    /* get parameters */
            goto PLPFin;
    
        n = ctx->PMRHSN;
        if (n < 4 || (n / 2) * 2 != n) {
            p_err(ctx, -25,1);
            goto PLPFin;
        }
    }
    else {
        if (parm(ctx, ctx->CmdBuf + 4,4,1))    /* get parameters */
            goto PLPFin;

        if (ctx->NOC < 2)
            goto PLPFin;

        n = ctx->NOC;
        if (ctx->PMNV != 2) {
            p_err(ctx, -1,1);
            goto PLPFin;
        }
        ix = ctx->PMVIdx[0];
        iy = ctx->PMVIdx[1];
    }


    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);

    /* first process fill option */

    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        ps_lwidth(ctx, 0.0);
        xa = ya = 0.0;

        first = 1;
        nn = i = 0;
        while (i < n) {

            if (typ == 0) {
                x = ps_2dx(ctx, ctx->PMRHSX[i]);
                y = ps_2dy(ctx, ctx->PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(ctx, get_data(ctx, ix,i));       
                y = ps_2dy(ctx, get_data(ctx, iy,i));       
            }
            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);              /* update bounding box */
            nn++;
            if (first) {
                fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,0.0);
                fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
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
                i++;
        }
        if (nn) {
            fprintf(ctx->PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);
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
                x = ps_2dx(ctx, ctx->PMRHSX[i]);
                y = ps_2dy(ctx, ctx->PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(ctx, get_data(ctx, ix,i));       
                y = ps_2dy(ctx, get_data(ctx, iy,i));       
            }
            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);              /* update bounding box */
            nn++;

            if (first) {
                fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
                first = 0;
            }
            else {         
                if (ctx->PMDIR == 1) {
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,ya);
                    if (ctx->PMNS != 1)
                        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
                    else
                        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
                }
                else if (ctx->PMDIR == 2) {
                    if (ctx->PMNS != 1) 
                        fprintf(ctx->PSFd,"%5.2f %5.2f l\n",xa,y);
                    else
                        fprintf(ctx->PSFd,"%5.2f %5.2f m\n",xa,y);
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
                }
                else
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
            }
            if (ctx->PMR)  
                fprintf(ctx->PSFd,"0 %5.2f rl 0 %5.2f rl\n",-y,y);

            i++;   
            if (typ == 0)
                i++;
            if (i >= n)
                break;
            xa = x;
            ya = y;
        }
        if (nn >= 2 && ctx->PMAFlg && ctx->PMA1 > 0.0 && ctx->PMA2 > 0.0) {       /* plot arrow */
        
            fprintf(ctx->PSFd,"gsave\ncurrentpoint\nstroke m\n");
            x -= xa;
            y -= ya;
            fprintf(ctx->PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
            fprintf(ctx->PSFd,"%5.2f %5.2f scale\n",ctx->PMA1,ctx->PMA2);
            fprintf(ctx->PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
        }
        if (nn)
            fprintf(ctx->PSFd,"stroke\n");
    }
    if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

        i = 0;
        while (i < n) {

            if (typ == 0) {
                x = ps_2dx(ctx, ctx->PMRHSX[i]);
                y = ps_2dy(ctx, ctx->PMRHSX[i + 1]);
            }
            else {
                if (eval_sve(ctx, i) == 0) {
                    i++;
                    continue;
                }
                x = ps_2dx(ctx, get_data(ctx, ix,i));       
                y = ps_2dy(ctx, get_data(ctx, iy,i));       
            }
            if (ctx->PMNC)
                upd_bbox(ctx, 1,x,y);                 /* update bounding box */
            ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            i++;   
            if (typ == 0)
                i++;
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLPFin:
    p_clean(ctx);
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

int pl_plotm(TDAContext *ctx)
{
    int j,icase,ix,iy,n,err,first;
    double x = 0.0,y = 0.0,xa,ya;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))    /* get parameters */
        goto PLMFin;

    n = ctx->PMNV / 2;
    if (ctx->PMNV < 2 || 2 * n != ctx->PMNV) {
        p_err(ctx, -1,1);
        goto PLMFin;
    }
       
    for (icase = 0; icase < ctx->NOC; ++icase) {

        fprintf(ctx->PSFd,"\n%%#%d: %s [case %d]\n",++ctx->PSONUM,ctx->CmdBuf,icase+1);

        fprintf(ctx->PSFd,"gsave\n");
        if (ctx->PMNC != 1)
            set_clip(ctx);

        /* first process fill option */

        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            ps_lwidth(ctx, 0.0);
            xa = ya = 0.0;
            first = 1;
            for (j = 0; j < ctx->PMNV; j += 2) {
                ix = ctx->PMVIdx[j];
                iy = ctx->PMVIdx[j + 1];
                x = ps_2dx(ctx, get_data(ctx, ix,icase));       
                y = ps_2dy(ctx, get_data(ctx, iy,icase));       

                if (ctx->PMNC)
                    upd_bbox(ctx, 1,x,y);              /* update bounding box */

                if (first) {
                    fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,0.0);
                    fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
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
            }
            fprintf(ctx->PSFd,"%5.2f %5.2f l\ngsave\n",x,0.0);
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\nnewpath\n");
        }

        /* now the standard polygon */

        if (ctx->PMLW > 0.0)
            ps_lwidth(ctx, ctx->PMLW);

        if (ctx->PMLW > 0.0 && ctx->PMLT > 0) {
            ps_ltyp(ctx, ctx->PMLT);
            xa = ya = 0.0;

            first = 1;
            for (j = 0; j < ctx->PMNV; j += 2) {
                ix = ctx->PMVIdx[j];
                iy = ctx->PMVIdx[j + 1];
                x = ps_2dx(ctx, get_data(ctx, ix,icase));       
                y = ps_2dy(ctx, get_data(ctx, iy,icase));       

                if (ctx->PMNC)
                    upd_bbox(ctx, 1,x,y);              /* update bounding box */

                if (first) {
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
                if (ctx->PMR)  
                    fprintf(ctx->PSFd,"0 %5.2f rl 0 %5.2f rl\n",-y,y);

            }
            if (n >= 2 && ctx->PMAFlg && ctx->PMA1 > 0.0 && ctx->PMA2 > 0.0) {       /* plot arrow */
        
                fprintf(ctx->PSFd,"gsave\ncurrentpoint\nstroke m\n");
                x -= xa;
                y -= ya;
                fprintf(ctx->PSFd,"%5.2f %5.2f\natan\nrotate\n",y,x);
                fprintf(ctx->PSFd,"%5.2f %5.2f scale\n",ctx->PMA1,ctx->PMA2);
                fprintf(ctx->PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\nclosepath\nfill\ngrestore\n");
            }
            fprintf(ctx->PSFd,"stroke\n");
        }
        if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

            fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

            for (j = 0; j < ctx->PMNV; j += 2) {
                ix = ctx->PMVIdx[j];
                iy = ctx->PMVIdx[j + 1];
                x = ps_2dx(ctx, get_data(ctx, ix,icase));       
                y = ps_2dy(ctx, get_data(ctx, iy,icase));       

                if (ctx->PMNC)
                    upd_bbox(ctx, 1,x,y);                 /* update bounding box */
                ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
            }
        }
        fprintf(ctx->PSFd,"grestore\n");
    }
    err = 0;

PLMFin:
    p_clean(ctx);
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

int pl_plotf(TDAContext *ctx)
{
    int err,r,first,n,deriv;
    double x,y,y1 = 0.0;
    char *p;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    deriv = 0;
    n = 5;
    p = ctx->CmdBuf + 5;
    if (*p == '1') {
        deriv = 1;
        n++;
    }
    else if (*p == '2') {
        deriv = 2;
        n++;
    }
    if (parm(ctx, ctx->CmdBuf + n,9,1))    /* get parameters */
        goto PLFFin;

    if (ctx->FNFlg == 0) {
        p_err(ctx, -1,1);
        goto PLFFin;
    }
    if (ctx->FNArgN != 1) {
        printf1(ctx, "Error: there must be exactly one function argument.\n");
        goto PLFFin;
    }
    if (ctx->PMFTYP5) {
        p_err(ctx, -40,1);
        goto PLFFin;
    }
    if (ctx->FVFlg)
        printf1(ctx, "Sum over %d data matrix cases.\n",ctx->NOC);

    if (deriv) {        /* need memory for gradient and hessian */
          
        if (fnd_alloc(ctx, 1,deriv,ctx->FNArgN,ctx->FNPN,0))
            goto PLFFin;

        if (alloc_acu(ctx, ctx->FNArgN + 1))
            goto PLFFin;

        if (alloc_acv(ctx, ctx->FNArgN * ctx->FNArgN + 1))
            goto PLFFin;
    }
    ctx->NINTMUsed = -1;

    if (ctx->PMRXFlg == 0) {
        ctx->PMRXA = (float)(ctx->PA1[0]);
        ctx->PMRXB = (float)(ctx->PA2[0]);
        ctx->PMRXD = (float)((double)((ctx->PMRXB - ctx->PMRXA)) / 100.0);
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMNC != 1)
        set_clip(ctx);

    x = (double)(ctx->PMRXA);
    first = 1;
    while ((double)(x) <= (double)((ctx->PMRXB)) + (double)(ctx->PMRXD)) {

        if ((double)(x) > (double)(ctx->PMRXB))
            x = (double)(ctx->PMRXB);

        ctx->FNArgVal[0] = x;

        r = get_flval(ctx, &y,ctx->FNArgN,ctx->FNArgVal,deriv,1,ctx->AcU,ctx->AcV,ctx->AcV);
        if (r) {      /* r from v_eval1(ctx) */

            printf1(ctx, "Can't evaluate function or derivatives.\n");
            prn_emsg2(ctx, r);
            goto PLFFin;
        }
        if (deriv == 1)
            y = ctx->AcU[1];
        else if (deriv > 1)
            y = ctx->AcV[1];
       
        if (first) {
            ps_2dplot(ctx, x,y,0);
            y1 = y;
            first = 0;
        }
        ps_2dplot(ctx, x,y,1);

        if (ctx->PMNC)           /* update bounding box */  
            upd_bbox(ctx, 0,x,y);     
    
        if (fabs((double)(x) - (double)(ctx->PMRXB)) <= ctx->EPSI)
            break;

        x += (double)(ctx->PMRXD);
    }
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0 && first == 0) {
        fprintf(ctx->PSFd,"gsave\n");
        ps_2dplot(ctx,(double)(ctx->PMRXB),ctx->PA1[1],1);    
        ps_2dplot(ctx,(double)(ctx->PMRXA),ctx->PA1[1],1);    
        ps_2dplot(ctx,(double)(ctx->PMRXA),y1,1);    
        ps_lwidth(ctx, 0.0);
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    fprintf(ctx->PSFd,"stroke\ngrestore\n");

    if (ctx->NINTMUsed >= 0)  
        ni_info(ctx);

    err = 0;

PLFFin:
    if (deriv)
        fnd_alloc(ctx, 0,0,0,0,0);
    p_clean(ctx);
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

int pl_ploth(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,ix;
    double tmp,xmin,xmax,wt,wsum,scal;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))    /* get parameters */
        goto PLHFin;

    if (ctx->PMNV != 1) {
        printf1(ctx, "Error: there must be exactly one variable name.\n");
        goto PLHFin;
    }
    n = ctx->PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (n < 1) {
        printf1(ctx, "Error: need classes defined with x parameter.\n");
        goto PLHFin;
    }
    if (ctx->DScalFlg && ctx->DScal > 0.0)
        scal = ctx->DScal;
    else
        scal = 1.0;

    for (j = 1; j < n; ++j) {
        if (ctx->PMTP[j] <= ctx->PMTP[j - 1]) {
            printf1(ctx, "Error: need non-empty classes.\n");
            goto PLHFin;
        }
    }
    ix = ctx->PMVIdx[0];
    xmin = ctx->PMTP[0];
    xmax = ctx->PMTP[n - 1];

    if (alloc_acw(ctx, n + 1))
        goto PLHFin;
    if (alloc_acy(ctx, n + 1))
        goto PLHFin;

    if (ctx->PMF1Def) {
        fprintf(ctx->PMF1d,"# Frequency distribution: %s\n",ctx->VName[ix]);
        if (ctx->PMWVar >= 0)  
            fprintf(ctx->PMF1d,"# Using weights defined by: %s\n",ctx->VName[ctx->PMWVar]);
        if (scal != 1.0)  
            fprintf(ctx->PMF1d,"# Scaling factor: %g\n",scal);
    }
    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    wsum = 0.0;
    wt = 1.0;
    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->PMWVar >= 0) {
            wt = get_data(ctx, ctx->PMWVar,i);
            if (wt < 0.0) {
                printf1(ctx, "Error: found negative weight in case %d\n",i + 1);
                goto PLHFin;
            }
        }
        tmp = get_data(ctx, ix,i);

        if (tmp < xmin) {
            printf1(ctx, "Error: value in case %d is %g (less than %g)\n",    
                                                    i + 1,tmp,xmin);
            goto PLHFin;
        }
        if (tmp > xmax) {
            printf1(ctx, "Error: value in case %d is %g (greater than %g)\n",    
                                                    i + 1,tmp,xmax);
            goto PLHFin;
        }
        k = n - 1;
        for (j = 1; j < n; ++j) {
            if (ctx->PMS == 1) {
                if (tmp <= ctx->PMTP[j]) {
                    k = j;
                    break;
                }
            }
            else {
                if (tmp < ctx->PMTP[j]) {
                    k = j;
                    break;
                }
            }
        }
        ctx->AcW[k] += wt;
        wsum += wt;
    }
    if (wsum <= 0.0) {
        printf1(ctx, "Error: sum of weights is zero.\n");
        goto PLHFin;
    }   
    for (j = 1; j < n; ++j)  
        ctx->AcY[j] = scal * ctx->AcW[j] / (wsum * (ctx->PMTP[j] - ctx->PMTP[j - 1]));

    if (ctx->PMF1Def) {
        for (j = 1; j < n; ++j) {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->PMTP[j - 1]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->PMTP[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcW[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcW[j] / wsum);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[j]);
            fprintf(ctx->PMF1d,"\n");
        }
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMNC != 1)
        set_clip(ctx);

    ps_2dplot(ctx, ctx->PMTP[0],0.0,0);
    if (ctx->PMNC)   
        upd_bbox(ctx, 0,ctx->PMTP[0],0.0);     

    for (j = 1; j < n; ++j) {
        ps_2dplot(ctx, ctx->PMTP[j - 1],ctx->AcY[j],1);
        if (ctx->PMNC)   
            upd_bbox(ctx, 0,ctx->PMTP[j - 1],ctx->AcY[j]);     

        ps_2dplot(ctx, ctx->PMTP[j],ctx->AcY[j],1);
        if (ctx->PMNC)   
            upd_bbox(ctx, 0,ctx->PMTP[j],ctx->AcY[j]);     
    }
    ps_2dplot(ctx, ctx->PMTP[n - 1],0.0,1);
      
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n");
        ps_2dplot(ctx, ctx->PMTP[0],0.0,1);    
        ps_lwidth(ctx, 0.0);
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    if (ctx->PMNS != 1) {
        for (j = 1; j < n - 1; ++j) {
            ps_2dplot(ctx, ctx->PMTP[j],0.0,0);
            ps_2dplot(ctx, ctx->PMTP[j],ctx->AcY[j],1);
        }   
    }
    fprintf(ctx->PSFd,"stroke\ngrestore\n");
    err = 0;

PLHFin:
    p_clean(ctx);
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

int pl_plotd(TDAContext *ctx)
{
    register int i,j;
    int err,nx,ix,m,gs;
    double w,a,b,ymax,scal;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,4,1))    /* get parameters */
        goto PLDFin;

    if (ctx->PMFmtF == 0)  
        pmfmt(ctx, 10,4);

    if (ctx->PMK < 2 || ctx->PMK > 4)
        ctx->PMK = 1;

    if (ctx->PMNV != 1) {
        printf1(ctx, "Error: there must be exactly one variable name.\n");
        goto PLDFin;
    }
    if (ctx->DScalFlg && ctx->DScal > 0.0)
        scal = ctx->DScal;
    else
        scal = 1.0;

    nx = ctx->PMNTP;      /* number of arguments in PMTP[] (t=...) */
    if (nx < 1) {
        printf1(ctx, "Error: need x parameter.\n");
        goto PLDFin;
    }
    ix = ctx->PMVIdx[0];

    if (ctx->PMD < ctx->EPSI1)
        ctx->PMD = 1.0;

    if (alloc_acx(ctx, ctx->NOC + 1))
        goto PLDFin;
    if (alloc_acu(ctx, ctx->NOC + 1))
        goto PLDFin; 
    if (alloc_acy(ctx, nx + 1))
        goto PLDFin;

    for (i = 0; i < ctx->NOC; ++i)
        ctx->AcX[i] = get_data(ctx, ix,i);

    if (sortd(ctx, ctx->NOC,ctx->AcX,0))     
        goto PLDFin;      

    if (ctx->PMF1Def) {
        fprintf(ctx->PMF1d,"# Density plot: %s\n",ctx->VName[ix]);
        fprintf(ctx->PMF1d,"# Kernel: ");
        switch (ctx->PMK) {
            case 2:  fprintf(ctx->PMF1d,"triangle\n"); break;
            case 3:  fprintf(ctx->PMF1d,"quartic\n"); break;
            case 4:  fprintf(ctx->PMF1d,"Epanechnikov\n"); break;
            default: fprintf(ctx->PMF1d,"uniform\n"); break;
        }
    }
    w = ctx->PMD / 2.0;
    ymax = 0.0;
    for (j = 0; j < nx; ++j) {
  
        a = ctx->PMTP[j] - w;
        b = ctx->PMTP[j] + w;
        m = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            if (ctx->AcX[i] < a)
                continue;
            if (ctx->AcX[i] > b)
                break;

            ctx->AcU[m++] = ctx->AcX[i];
        }
        if (m > 0) {
            ctx->AcY[j] = mkern1(ctx, ctx->PMK,m,ctx->PMD,ctx->PMTP[j],ctx->AcU) / (double)ctx->NOC;
            ymax = dmax(ctx, ymax,ctx->AcY[j]);
        }
        else
            ctx->AcY[j] = 0.0;

        ctx->AcY[j] *= scal;

        if (ctx->PMF1Def) {
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->PMTP[j]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[j]);
            fprintf(ctx->PMF1d,"\n");
        }
    }
    printf1(ctx, "Max height of density: %g\n",ymax);

    gs = 0;
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0)  
        gs = 1;

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMNC != 1)
        set_clip(ctx);

    if (gs) {
        ps_2dplot(ctx, ctx->PMTP[0],0.0,0);
        ps_2dplot(ctx, ctx->PMTP[0],ctx->AcY[0],1);
    }
    else
        ps_2dplot(ctx, ctx->PMTP[0],ctx->AcY[0],0);
    if (ctx->PMNC)   
        upd_bbox(ctx, 0,ctx->PMTP[0],ctx->AcY[0]);     

    for (j = 1; j < nx; ++j) {
        ps_2dplot(ctx, ctx->PMTP[j],ctx->AcY[j],1);
        if (ctx->PMNC)   
            upd_bbox(ctx, 0,ctx->PMTP[j],ctx->AcY[j]);     
    }
    if (gs) {
        ps_2dplot(ctx, ctx->PMTP[nx - 1],0.0,1);
        ps_2dplot(ctx, ctx->PMTP[0],0.0,1);
        fprintf(ctx->PSFd,"gsave\n");
        ps_lwidth(ctx, 0.0);
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    fprintf(ctx->PSFd,"stroke\ngrestore\n");
    err = 0;

PLDFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_ploto(typ)       Plot cirlce or ellipse                              */
/*                                                                          */
/*          typ 0 :     ploto(pos,lt,lw,gs,nc) = r [,alpha,beta];           */
/*          typ 1 :     plote(pos,lt,lw,gs,nc) = a,b[,r];                   */
/*                      Return 0 if OK, -1 if error.                        */

int pl_ploto(TDAContext *ctx, int typ)
{
    int err;
    double x,y,a,b,r;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,7,1))    /* get parameters */
        goto PLOFin;

    if (ctx->PMXYFlg == 0) {
        p_err(ctx, -1,1);
        goto PLOFin;
    }
    if (typ == 0) {
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
            goto PLOFin;
        }
        r = ctx->PMRHSX[0];
        if (r <= 0.0 || a < 0.0 || b <= a || b > 360.0) {
            p_err(ctx, -1,1);
            goto PLOFin;
        }
    }
    else {
        if (ctx->PMRHSN != 2 && ctx->PMRHSN != 3) {
            p_err(ctx, -1,1);
            goto PLOFin;
        }
        a = ctx->PMRHSX[0];
        b = ctx->PMRHSX[1];
        if (ctx->PMRHSN == 3)
            r = ctx->PMRHSX[2];
        else
            r = 0.0;

        if (a <= 0.0 || b <= 0.0 || r < -360.0 || r >= 360.0) {
            p_err(ctx, -1,1);
            goto PLOFin;
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMNC != 1)
        set_clip(ctx);

    x = ps_2dx(ctx, ctx->PMX);      /* translate to PostScript coordinates */
    y = ps_2dy(ctx, ctx->PMY);

    if (typ == 0) {
        r = ps_2dx(ctx, r) - ps_2dx(ctx, 0.0);
            
        fprintf(ctx->PSFd,"%5.2f %5.2f m\nstroke\n",x,y);
        fprintf(ctx->PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f arc\n",x,y,r,a,b);

        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            if (a > 0.0 || b < 360.0)  
                fprintf(ctx->PSFd,"%5.2f %5.2f l\nclosepath\n",x,y);
            fprintf(ctx->PSFd,"gsave\n");
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\n");
        }
        a = b = r;
    }
    else {

        a = ps_2dx(ctx, a) - ps_2dx(ctx, 0.0);
        b = ps_2dy(ctx, b) - ps_2dy(ctx, 0.0);
        fprintf(ctx->PSFd,"%5.2f %5.2f translate\n",x,y);
        fprintf(ctx->PSFd,"%5.2f rotate\n",r);

        if (a != 0.0)
            fprintf(ctx->PSFd," 1.00 %5.2f scale\n",b / a);
        fprintf(ctx->PSFd,"0 0 %5.2f 0 360 arc\n",a);

        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            fprintf(ctx->PSFd,"gsave\n");
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\n");
        }
    }
    fprintf(ctx->PSFd,"stroke\ngrestore\n");
    if (ctx->PMNC) {
        upd_bbox(ctx, 1,x + a,y);
        upd_bbox(ctx, 1,x - a,y);
        upd_bbox(ctx, 1,x,y + b);
        upd_bbox(ctx, 1,x,y - b);
    }
    err = 0;

PLOFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotk()          Plot arc.                                           */
/*                                                                          */
/*                      plotk(lt,lw,a) = x1,y1,x2,y2,...;                   */
/*                      a = arrow                                           */  
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotk(TDAContext *ctx)
{
    int i,n,err,first,ptyp;
    double x,y,xl,yl,tmp,tmp1,tmp2,rot,phik,phil,r;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,7,1))  /* get parameters */
        goto PLKFin;

    if (ctx->PMSC >= 180.0 || ctx->PMSC <= -180.0) {
        p_err(ctx, -1,1);
        goto PLKFin;
    }
    n = ctx->PMRHSN;
    if (n < 4 || (n / 2) * 2 != n) {
        p_err(ctx, -25,1);
        goto PLKFin;
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMNC != 1)
        set_clip(ctx);

    if (fabs(ctx->PMSC) > ctx->EPSI1)
        n -= 2;
    first = 1;
    for (i = 0; i < n; i += 2) {

        x = ps_2dx(ctx, ctx->PMRHSX[i]);
        y = ps_2dy(ctx, ctx->PMRHSX[i + 1]);
        if (ctx->PMNC) 
            upd_bbox(ctx, 1,x,y);              /* update bounding box */

        if (first) {
            fprintf(ctx->PSFd,"%5.2f %5.2f m\n",x,y);
            first = 0;
        }
        else
            fprintf(ctx->PSFd,"%5.2f %5.2f l\n",x,y);
    }
    if (fabs(ctx->PMSC) <= ctx->EPSI1)
        goto PLKF;
     
    /* plot the final arc */

    x  = ps_2dx(ctx, ctx->PMRHSX[n]);                           
    xl = ps_2dx(ctx, ctx->PMRHSX[n - 2]);
    y  = ps_2dy(ctx, ctx->PMRHSX[n + 1]);                           
    yl = ps_2dy(ctx, ctx->PMRHSX[n - 1]);
    tmp1 = x - xl;
    tmp2 = y - yl;
    if (ctx->PMNC) 
        upd_bbox(ctx, 1,x,y);              /* update bounding box */

    if (fabs(tmp1) > ctx->EPSI1 || fabs(tmp2) > ctx->EPSI1) {

        if (fabs(tmp1) > ctx->EPSI1)
            tmp = atan(tmp2 / tmp1);
        else
            tmp = Pi / 2.0;

        if (ctx->PMSC < 0.0) {
            ctx->PMSC = -ctx->PMSC;
            ptyp = 1;
        }
        else
            ptyp = 0;

        rot = tmp * 180 / Pi - ctx->PMSC;
        if (tmp1 > 0.0 || (tmp1 == 0.0 && tmp2 > 0.0)) 
            tmp += 0.5 * Pi;
        else {
            tmp += 1.5 * Pi;
            rot += 180;
        }
        phik = tmp - ctx->PMSC * Pi / 180.0;
        phil = tmp + ctx->PMSC * Pi / 180.0;
        r = 0.0;
        tmp = cos(phik) - cos(phil);
        if (fabs(tmp) > ctx->EPSI1)  
            r = tmp1 / tmp;

        if (fabs(r) <= ctx->EPSI1) {
            tmp = sin(phik) - sin(phil);
            if (fabs(tmp) > ctx->EPSI1)  
                r = tmp2 / tmp;
        }
        if (fabs(r) > ctx->EPSI1) {

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
                rot += 2.0 * ctx->PMSC;
            }
            fprintf(ctx->PSFd,"stroke\n%5.2f %5.2f ",tmp1,tmp2);
            fprintf(ctx->PSFd,"%5.2f %5.2f %5.2f arc\n",r,phik,phil);

            if (ctx->PMAFlg && ctx->PMA1 > 0.0 && ctx->PMA2 > 0.0) {       /* plot arrow */
                fprintf(ctx->PSFd,"%5.2f %5.2f m\ncurrentpoint\nstroke\nm\n",x,y);
                fprintf(ctx->PSFd,"%5.2f rotate\n%5.2f %5.2f scale\n",rot,ctx->PMA1,ctx->PMA2);
                fprintf(ctx->PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\n");
                fprintf(ctx->PSFd,"closepath\nfill\n");
            }
        }
    }
PLKF:
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLKFin:
    p_clean(ctx);
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

int pl_plotch(TDAContext *ctx)
{
    register int i,j,ik,l;
    int ix,iy,m = 0,n,nn,r,err,ns;
    double x,y,dx,dy,xm,ym,xa,xb,ya,yb;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,4,1))    /* get parameters */
        goto PLHFin;

    if (ctx->NOC < 2)
        goto PLHFin;

    if (ctx->PMNV != 2) {
        p_err(ctx, -1,1);
        goto PLHFin;
    }
    ix = ctx->PMVIdx[0];
    iy = ctx->PMVIdx[1];

    if (alloc_acxf(ctx, ctx->NOC + 1))
        goto PLHFin;
    if (alloc_acyf(ctx, ctx->NOC + 1))
        goto PLHFin;
    if (alloc_acn(ctx, ctx->NOC + 1))
        goto PLHFin;
    if (alloc_aci(ctx, ctx->NOC + 1))
        goto PLHFin;
    if (alloc_acj(ctx, ctx->NOC + 1))
        goto PLHFin;

    nn = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if (eval_sve(ctx, i) == 0)   
            continue;
        nn++;
        ctx->AcXF[nn] = (float)get_data(ctx, ix,i);       
        ctx->AcYF[nn] = (float)get_data(ctx, iy,i);       
        ctx->AcN[nn] = nn;
    }
    n = g_chull(ctx, ctx->AcXF,ctx->AcYF,nn,ctx->AcN,ctx->AcJ,ctx->AcI);   
    if (n <= 0) {
        if (n < 0)
            p_err(ctx, -2,1);
        goto PLHFin;
    }
    if (ctx->PMNS >= 2 && n < 3)  
        ctx->PMNS = 1;

    if (ctx->PMNS >= 2 && n > 2) {        /* use Akima algorithm for smoothing */

        ns = n + 4;
        if (alloc_acx(ctx, ns + 1))
            goto PLHFin;
        if (alloc_acy(ctx, ns + 1))
            goto PLHFin;

        xm = ym = 0.0;

        ik = ctx->AcI[1];                        /* get points of convex hull */
        for (i = 1; i <= n; ++i) {
            j = ctx->AcJ[ik];
            ctx->AcX[i + 2] = (double)ctx->AcXF[j];
            ctx->AcY[i + 2] = (double)ctx->AcYF[j];
            ik = ctx->AcI[ik];
        }
        ctx->AcX[n + 3] = ctx->AcX[3];
        ctx->AcY[n + 3] = ctx->AcY[3];
        ctx->AcX[n + 4] = ctx->AcX[4];
        ctx->AcY[n + 4] = ctx->AcY[4];

        ctx->AcX[1] = ctx->AcX[n + 1];
        ctx->AcY[1] = ctx->AcY[n + 1];
        ctx->AcX[2] = ctx->AcX[n + 2];
        ctx->AcY[2] = ctx->AcY[n + 2];
             
        if (ctx->PMIC > 0.0) {
            for (i = 1; i <= ns; ++i) {
                xm += ctx->AcX[i];
                ym += ctx->AcY[i];
            }
            xm /= (double)ns;
            ym /= (double)ns;

            xa = ya = xb = yb = 0.0;
            for (i = 1; i <= ns; ++i) {
                xa = dmax(ctx, xa,xm - ctx->AcX[i]);
                xb = dmax(ctx, xb,ctx->AcX[i] - xm);
                ya = dmax(ctx, ya,ym - ctx->AcY[i]);
                yb = dmax(ctx, yb,ctx->AcY[i] - ym);
            }
            dx = ctx->PMIC * ctx->UXLen / ctx->PXLen;
            dy = ctx->PMIC * ctx->UYLen / ctx->PYLen;

            for (i = 1; i <= ns; ++i) {
                x = ctx->AcX[i];
                y = ctx->AcY[i];
                if (x < xm - ctx->EPSI1) {
                    x -= dx * (xm - x) / xa;
                    if (y > ym + ctx->EPSI1)  
                        y += dy * (y - ym) / yb;
                    else if (y < ym - ctx->EPSI1)
                        y -= dy * (ym - y) / ya;
                }
                else if (x > xm + ctx->EPSI1) {
                    x += dx * (x - xm) / xb;
                    if (y > ym + ctx->EPSI1)
                        y += dy * (y - ym) / yb;
                    else if (y < ym - ctx->EPSI1)
                        y -= dy * (ym - y) / ya;
                }
                ctx->AcX[i] = x;
                ctx->AcY[i] = y;
            }
        }
        m = (ns - 1) * ctx->PMNS + 1;

        if (alloc_acu(ctx, m + 1))
            goto PLHFin;
        if (alloc_acv(ctx, m + 1))
            goto PLHFin;

        r = splakima(ctx, 2,ns,ctx->AcX,ctx->AcY,ctx->PMNS,m,ctx->AcU,ctx->AcV);

        if (r) {
            printf1(ctx, "Error (%d) in smoothing algorithm.\n",r);
            ctx->PMNS = -1;      /* continue with standard method */
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, ctx->PMLT);
   
    if (ctx->PMNS >= 2) {

        l = 1;
        for (i = 2 * ctx->PMNS + 1; i <= m - ctx->PMNS; ++i) {
            if (l) {
                ps_2dplot(ctx, ctx->AcU[i],ctx->AcV[i],0);
                l = 0;
            }
            else
                ps_2dplot(ctx, ctx->AcU[i],ctx->AcV[i],1);

            if (ctx->PMNC)
                upd_bbox(ctx, 0,ctx->AcU[i],ctx->AcV[i]);
        }
    }
    else {

        ik = ctx->AcI[1];
        l = 0;
        for (i = 1; i <= n; ++i) {
            j = ctx->AcJ[ik];
            x = (double)ctx->AcXF[j];
            y = (double)ctx->AcYF[j];
            ps_2dplot(ctx, x,y,l); 
            if (ctx->PMNC)
                upd_bbox(ctx, 0,x,y);
            l = 1;
            ik = ctx->AcI[ik];
        }
        ik = ctx->AcI[1];
        j = ctx->AcJ[ik];
        ps_2dplot(ctx, (double)ctx->AcXF[j],(double)ctx->AcYF[j],1);     /* close path */
    }
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n");       
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");   
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");

    if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

        for (i = 0; i < ctx->NOC; ++i) {
            if (eval_sve(ctx, i) == 0)   
                continue;

            x = ps_2dx(ctx, get_data(ctx, ix,i));       
            y = ps_2dy(ctx, get_data(ctx, iy,i));       
            ps_sym(ctx, ctx->PMS,x,y,ctx->PtMM * ctx->PMFS / 2.0);
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLHFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plots(typ)       Smoothing with Akima algorithm.                     */
/*                                                                          */
/*              typ 0:  plotsp(ns,lt,lw,gs,nc,sc,s,fs) = x1,y1,x2,y2,...;   */
/*              typ 1:  plots (ns,lt,lw,gs,nc,sc,sel,s,fs) = VX,VY;         */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plots(TDAContext *ctx, int typ)
{
    register int i,j;
    int ix = 0,iy = 0,m,n,nn,r,err,first;
    double x,y,xa,ya;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (typ == 0) {
        if (parm(ctx, ctx->CmdBuf + 6,7,1))    /* get parameters */
            goto PLSFin;

        n = ctx->PMRHSN / 2;
        if (n < 2 || 2 * n != ctx->PMRHSN) {
            p_err(ctx, -25,1);
            goto PLSFin;
        }
    }
    else {
        if (parm(ctx, ctx->CmdBuf + 5,4,1))    /* get parameters */
            goto PLSFin;

        if (ctx->NOC < 2)
            goto PLSFin;

        n = ctx->NOC;
        if (ctx->PMNV != 2) {
            p_err(ctx, -1,1);
            goto PLSFin;
        }
        ix = ctx->PMVIdx[0];
        iy = ctx->PMVIdx[1];
    }
    if (ctx->PMNS < 2)                   /* number of subintervals */
        ctx->PMNS = 2;

    if (alloc_acx(ctx, n + 2))
        goto PLSFin;
    if (alloc_acy(ctx, n + 2))
        goto PLSFin;

    nn = j = 0;                     /* get data points */
    xa = ya = 0.0;
    for (i = 0; i < n; ++i) {

        if (typ == 0) {
            x = ctx->PMRHSX[j++];
            y = ctx->PMRHSX[j++];
        }
        else {
            if (eval_sve(ctx, i) == 0)  
                continue;
               
            x = get_data(ctx, ix,i);       
            y = get_data(ctx, iy,i);       
        }
        if (nn == 0 || x != xa || y != ya) {
            nn++;
            ctx->AcX[nn] = x;
            ctx->AcY[nn] = y;
        }   
        xa = x;
        ya = y;
    }
    if (nn < 2) {
        p_err(ctx, -26,1);
        goto PLSFin;
    }
    if ((int)ctx->PMSC == 1) {           /* closed curve */

        if (ctx->AcX[nn] != ctx->AcX[1] || ctx->AcY[nn] != ctx->AcY[1]) {
            nn++;
            ctx->AcX[nn] = ctx->AcX[1];
            ctx->AcY[nn] = ctx->AcY[1];
        }
    }
    m = (nn - 1) * ctx->PMNS + 1;

    if (alloc_acu(ctx, m + 1))
        goto PLSFin;
    if (alloc_acv(ctx, m + 1))
        goto PLSFin;

    r = splakima(ctx, 2,nn,ctx->AcX,ctx->AcY,ctx->PMNS,m,ctx->AcU,ctx->AcV);
    if (r) {
        printf1(ctx, "Error in smoothing algorithm.\n");
        goto PLSFin;
    }
    if (ctx->PMF1Def) {      /* write into output file */
        for (i = 1; i <= m; ++i) {
            fprintf(ctx->PMF1d,"%6d ",i);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcU[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcV[i]);
            fprintf(ctx->PMF1d,"\n");
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, ctx->PMLT);
   
    first = 1;
    for (i = 1; i <= m; ++i) {
        if (first) {
            ps_2dplot(ctx, ctx->AcU[i],ctx->AcV[i],0);
            first = 0;
        }
        else
            ps_2dplot(ctx, ctx->AcU[i],ctx->AcV[i],1);

        if (ctx->PMNC)
            upd_bbox(ctx, 0,ctx->AcU[i],ctx->AcV[i]);
    }

    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0 && first == 0) {
        if ((int)ctx->PMSC == 1) {
            fprintf(ctx->PSFd,"gsave\n");       
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\n");   
        }
        else if (first == 0) {
            fprintf(ctx->PSFd,"gsave\n");
            ps_2dplot(ctx, ctx->AcU[m],ctx->PA1[1],1);    
            ps_2dplot(ctx, ctx->AcU[1],ctx->PA1[1],1);    
            ps_2dplot(ctx, ctx->AcU[1],ctx->AcV[1],1);    
            ps_lwidth(ctx, 0.0);
            ps_fill(ctx, ctx->PMGS);
            fprintf(ctx->PSFd,"grestore\n");
        }
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");

    if (ctx->PMS >= 1 && ctx->PMS <= 17 && ctx->PMFS > 0.0) {  /* plot symbols */

        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

        ps_sym(ctx, ctx->PMS,ps_2dx(ctx, ctx->AcU[1]),ps_2dy(ctx, ctx->AcV[1]),ctx->PtMM * ctx->PMFS / 2.0);
        ps_sym(ctx, ctx->PMS,ps_2dx(ctx, ctx->AcU[m]),ps_2dy(ctx, ctx->AcV[m]),ctx->PtMM * ctx->PMFS / 2.0);
    }
    fprintf(ctx->PSFd,"grestore\n");
    err = 0;

PLSFin:
    p_clean(ctx);
    return(err);
}



void tda_reset_t_plot(void)
{
    BBLX = BBLY = BBUX = BBUY = 0;
    PSPR11 = PSPR12 = PSPR13 = PSPR21 = PSPR22 = PSPR23 = 0.0;
    PSPR31 = PSPR32 = PSPR33 = 0.0;
    PSPI11 = PSPI12 = PSPI13 = PSPI21 = PSPI22 = PSPI23 = 0.0;
    PSPI31 = PSPI32 = PSPI33 = 0.0;
}
