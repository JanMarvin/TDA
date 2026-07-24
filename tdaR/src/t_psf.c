/****************************************************************************/
/*  t_psf                                                                   */
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
#include "t_plot.h"
#include "t_parm.h"
#include "t_dplot.h"
#include "t_psinit.h"
#include "t_cplot.h"
#include "t_xplot.h"
#include "t_plot3.h"
#include "t_gf.h"
#include "t_plot3.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_psf.c                                                    */

int psfile(TDAContext *ctx);
int psetup(TDAContext *ctx);
void setup_2dsys(TDAContext *ctx);
void prn_psfe(TDAContext *ctx, char *s);
int ps_close(TDAContext *ctx, int opt,int scx,int kflag);
void rot_adjust(TDAContext *ctx);
int pl_axis(TDAContext *ctx, int typ,int opt,int opt1,int lt,double lw);
int psetup3(TDAContext *ctx);
void setup_3dsys(TDAContext *ctx);
void setup_3proj(TDAContext *ctx, double lon,double lat);

/*--------------------------------------------------------------------------*/
/*  psfile()    Open PostScript file.                                       */
/*              CmdBuf: psf = fname.                                        */
/*              Return 0 if OK, -1 if error.                                */

int psfile(TDAContext *ctx)
{
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->PSFFlg) {
        printf1(ctx, "Current PostScript file will be closed.\n");
        ps_close(ctx, 0,1,0);
    }
    strcpy(ctx->PSFName,ctx->CmdBuf + 7);

    if (!(ctx->PSFd = fopen(ctx->PSFName,OPEN_WR))) {
        printf1(ctx, "Error: can't create: %s\n",ctx->PSFName);
        return(-1);
    }
    ps_init(ctx);          /* write header to output file */
    ctx->PSFFlg = 1;
    ctx->PSPROJ = ctx->PS3DFlg = ctx->PSONUM = 0;
    printf1(ctx, "Opened new PostScript file: %s\n",ctx->PSFName);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  psetup()    Setup plot environment.                                     */
/*              CmdBuf: psetup(...)                                         */
/*              Return 0 if OK, -1 if error.                                */

int psetup(TDAContext *ctx)
{
    int n,m;           
    double x,y;
    register char *p;

    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->PSFFlg == 0) {
        printf1(ctx, "Error: need a PostScript file.\n");
        return(-1);
    }
    ctx->PSLog[0] = ctx->PSLog[1] = -1;                           

    ctx->XOrg = 150;         /* origin of PostScript figure */
    ctx->YOrg = 460;  
    ctx->SCALXFac = 1.0;     /* PostScript x scaling factor */
    ctx->SCALYFac = 1.0;     /* PostScript y scaling factor */
    ctx->ROTFac = 0.0;       /* PostScript rotation factor */

    ctx->PXLen = 120.0;      /* phys length of X axis in mm */
    ctx->PYLen =  80.0;      /* phys length of Y axis in mm */

    p = ctx->CmdBuf + 6;
    while (*++p) {
        if (sscanf(p,"pxlen=%lg",&x) == 1 && x > 0) {
            ctx->PXLen = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"pylen=%lg",&y) == 1 && y > 0) {
            ctx->PYLen = y;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            ctx->XOrg = n;
            ctx->YOrg = m;
            p = skip_int(ctx, p + 6);
            p = skip_int(ctx, p + 1);
        }
        else if (sscanf(p,"psscal=%lg,%lg",&x,&y) == 2 && x > 0.0 && y > 0.0) {
            ctx->SCALXFac = x;
            ctx->SCALYFac = y;
            p = skip_dbl(ctx, p + 7);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"psrot=%lg",&x) == 1 && x >= 0.0 && x < 360.0) {
            ctx->ROTFac = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"pxa(log=%d)=%lg,%lg",&n,&x,&y) == 3 && x < y &&
            (n == 0 || n == 1)) {      
            ctx->PA1[0] = x;
            ctx->PA2[0] = y;
            ctx->PSLog[0] = n;
            p = skip_int(ctx, p + 8);
            p = skip_dbl(ctx, p + 2);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"pxa=%lg,%lg",&x,&y) == 2 && x < y) {   
            ctx->PA1[0] = x;
            ctx->PA2[0] = y;
            ctx->PSLog[0] = 0;
            p = skip_dbl(ctx, p + 4);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"pya(log=%d)=%lg,%lg",&n,&x,&y) == 3 && x < y &&
            (n == 0 || n == 1)) {      
            ctx->PA1[1] = x;
            ctx->PA2[1] = y;
            ctx->PSLog[1] = n;
            p = skip_int(ctx, p + 8);
            p = skip_dbl(ctx, p + 2);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"pya=%lg,%lg",&x,&y) == 2 && x < y) {   
            ctx->PA1[1] = x;
            ctx->PA2[1] = y;
            ctx->PSLog[1] = 0;
            p = skip_dbl(ctx, p + 4);
            p = skip_dbl(ctx, p + 1);
        }
        if (*p != ',' && *p != ')') {
            prn_psfe(ctx, p);
            return(-1);    
        }
        *p = '\0';
    }
    if (ctx->PSLog[0] < 0 || ctx->PSLog[1] < 0) {
        printf1(ctx, "Error: need logical coordinates (pxa and pya parameters).\n");
        return(-1);
    }
    setup_2dsys(ctx);
    return(0);
}

/*--##----------------------------------------------------------------------*/
/*  setup_2dsys(lon,lat)        Set up 2d coordinate system.                */

void setup_2dsys(TDAContext *ctx)
{

    ctx->UXLen = ctx->PA2[0] - ctx->PA1[0];
    ctx->UYLen = ctx->PA2[1] - ctx->PA1[1];

    ctx->PSXLen = ctx->PtMM * ctx->PXLen;   /* phys. x axis length points */
    ctx->PSYLen = ctx->PtMM * ctx->PYLen;   /* phys. y axis length points */

    printf1(ctx, "Creating a new 2d coordinate system.\n");
    printf1(ctx, "Size (width x height in mm): %g x %g\n",ctx->PXLen,ctx->PYLen);
    printf1(ctx, "PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    ctx->XOrg,ctx->XOrg + (int)(ctx->PSXLen + 0.5),ctx->YOrg,ctx->YOrg + (int)(ctx->PSYLen + 0.5));   
    printf1(ctx, "Origin (PostScript x and y coordinates): %d,%d\n",ctx->XOrg,ctx->YOrg);
    printf1(ctx, "Scaling (x and y axis): %g,%g\n",ctx->SCALXFac,ctx->SCALYFac);
    printf1(ctx, "Rotation (degrees): %g\n\n",ctx->ROTFac);

    if (ctx->PSFFlg == 1) {      /* init bounding box parameters */
        ctx->BBLX = ctx->XOrg - 5;
        ctx->BBLY = ctx->YOrg - 5;
        ctx->BBUX = (int)(ctx->BBLX + (int)ctx->PSXLen + 10.0);
        ctx->BBUY = (int)(ctx->BBLY + (int)ctx->PSYLen + 10.0);
    }
    ps_inis(ctx, 1);

    printf1(ctx, "User coordinates X axis: %16.4f, %16.4f ",ctx->PA1[0],ctx->PA2[0]);
    if (ctx->PSLog[0])
        printf1(ctx, "(logarithmic)");
    printf1(ctx, "\nUser coordinates Y axis: %16.4f, %16.4f ",ctx->PA1[1],ctx->PA2[1]);
    if (ctx->PSLog[1])
        printf1(ctx, "(logarithmic)");
    printf1(ctx, "\n");

    ctx->PSFFlg = 2;              /* set to 2 if setup successful */
    ctx->PSPROJ = ctx->PS3DFlg = 0;    /* 2-d coordinate system */
}

/* ------------------------------------------------------------------------ */
/*  prn_nve.    Print error message.                                        */

void prn_psfe(TDAContext *ctx, char *s)
{
    register char *p = s;

    printf1(ctx, "Syntax error: ");
    if (!*p)
        printf1(ctx, "check brackets and semicolon.\n");
    else {     
        while (*p && p < s + 20)
            printf1(ctx, "%c",*p++);
        if (*p)
            printf1(ctx, " ...");
        printf1(ctx, "\n");
    }
}

/*--------------------------------------------------------------------------*/
/*  ps_close(opt,scx,kflag)  Close PostScript file.                         */
/*                           opt != 0 make message.                         */
/*                           if scx also free scx.                          */
/*                           if kflag != 0 keep parameters.                 */
/*                           Return 0 if OK, -1 if error.                   */

int ps_close(TDAContext *ctx, int opt,int scx,int kflag)
{
    int err = 0;                
    long pos;
    char buf[132];
         
    err = 0;
    if (ctx->PSFFlg) {
        if (opt)  
            printf1(ctx, "> Closed current PostScript file: %s\n",ctx->PSFName);

        fprintf(ctx->PSFd,"showpage\n");
        fclose(ctx->PSFd);
    }
    if (ctx->PSFFlg == 2) {     /* reprint bounding box */

        if (!(ctx->PSFd = fopen(ctx->PSFName,OPEN_RU)))  
            err = -1;
        else {

            if (ctx->PS3DFlg) {
                ctx->BBLX -= 4;
                ctx->BBLY -= 4;
                ctx->BBUX += 4;
                ctx->BBUY += 4;
            }
            if (ctx->ROTFac != 0.0) {        /* adjust for rotation */
                rot_adjust(ctx);
            }

            while (fgets(buf,131,ctx->PSFd)) {
                pos = ftell(ctx->PSFd);
                if (!strncmp(buf + 2,"CreationDate",12)) {
                    fflush(ctx->PSFd);
                    fseek(ctx->PSFd,pos,0);
                    fprintf(ctx->PSFd,"%%%%BoundingBox: %5d %5d %5d %5d\n",ctx->BBLX,ctx->BBLY,ctx->BBUX,ctx->BBUY);
                    break;
                }
            }
            fclose(ctx->PSFd);
        }
        if (err)  
            printf1(ctx, "\nError: can't update bounding box information in %s.\n",ctx->PSFName);
    }
    if (kflag == 0) {
        ctx->PSONUM = ctx->PSYOrg = ctx->PSXOrg = ctx->BBLX = ctx->BBLY = ctx->BBUX = ctx->BBUY = ctx->PSPROJ = ctx->PS3DFlg = ctx->PSFFlg = 0;
    }
    if (scx)
        alloc_scx(ctx, 0);
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  rot_adjust.  Adjust bounding box for rotation.                          */

void rot_adjust(TDAContext *ctx)
{
    double phi,sphi,cphi,xa,ya,xb,yb,x,y,xmin,xmax,ymin,ymax;

    phi = degree_to_arc(ctx, ctx->ROTFac);
    sphi = sin(phi);
    cphi = cos(phi);

    xa = (double)(ctx->BBLX - ctx->XOrg);
    xb = (double)(ctx->BBUX - ctx->XOrg);
    ya = (double)(ctx->BBLY - ctx->YOrg);
    yb = (double)(ctx->BBUY - ctx->YOrg);

    x = xa * cphi - ya * sphi;
    y = xa * sphi + ya * cphi;
    xmin = xmax = x;
    ymin = ymax = y;

    x = xa * cphi - yb * sphi;
    y = xa * sphi + yb * cphi;
    xmin = dmin(ctx, xmin,x);
    xmax = dmax(ctx, xmax,x);
    ymin = dmin(ctx, ymin,y);
    ymax = dmax(ctx, ymax,y);

    x = xb * cphi - ya * sphi;
    y = xb * sphi + ya * cphi;
    xmin = dmin(ctx, xmin,x);
    xmax = dmax(ctx, xmax,x);
    ymin = dmin(ctx, ymin,y);
    ymax = dmax(ctx, ymax,y);

    x = xb * cphi - yb * sphi;
    y = xb * sphi + yb * cphi;
    xmin = dmin(ctx, xmin,x);
    xmax = dmax(ctx, xmax,x);
    ymin = dmin(ctx, ymin,y);
    ymax = dmax(ctx, ymax,y);

    ctx->BBLX = ctx->XOrg + (int)floor(xmin) - 1;
    ctx->BBUX = ctx->XOrg + (int)ceil(xmax)  + 1;
    ctx->BBLY = ctx->YOrg + (int)floor(ymin) - 1;
    ctx->BBUY = ctx->YOrg + (int)ceil(ymax)  + 1;
}

/*--------------------------------------------------------------------------*/
/*  pl_axis(typ,opt,opt1,lt,lw)                                             */
/*                                                                          */
/*                      plot axis: plxa and plya command.                   */
/*                      update LWCSMax.                                     */
/*                      If opt=1 direct call from xplot()                   */
/*                      if opt1 = 0 use PMLT and PMLW, otherwise lt, lw     */
/*                                                                          */
/*                      Return 0 if OK, -1 if error.                        */

int pl_axis(TDAContext *ctx, int typ,int opt,int opt1,int lt,double lw)
{
    int err,n,i,len;
    double xa,xb,ya,yb,x,xx,y,yy,d = 0.0,s,tlen,flen;
    char buf[100];

    err = -1;

    if (opt == 0) {
        if (check_pcmd(ctx, 1,2))
            return(-1);

        /* printf1("Plot %c axis in user's coordinate system.\n",'X' + typ); */

        if (parm(ctx, ctx->CmdBuf + 4,7,0))    /* get parameters */
            goto PLAFin;
    }
    if (ctx->PMRHSN > 0 && opt == 0) {
        if (ctx->PMRHSN != 4) {
            p_err(ctx, -1,1);
            goto PLAFin;
        }
        xa = ctx->PMRHSX[0];
        ya = ctx->PMRHSX[1];
        xb = ctx->PMRHSX[2];
        yb = ctx->PMRHSX[3];
    }
    else if (typ == 0) {
        xa = ctx->PA1[0];
        ya = ctx->PA1[1];
        xb = ctx->PA2[0];
        yb = ctx->PA1[1];
    }
    else  {
        xa = ctx->PA1[0];
        ya = ctx->PA1[1];
        xb = ctx->PA1[0];
        yb = ctx->PA2[1];
    }
    if (opt == 1) {
        if (typ == 0)
            fprintf(ctx->PSFd,"\n%%#%d: plxa\n",++ctx->PSONUM);
        else
            fprintf(ctx->PSFd,"\n%%#%d: plya\n",++ctx->PSONUM);
    }
    else
        fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);
    fprintf(ctx->PSFd,"%%#%c-axis\n",'x' + typ);
    if (ctx->PMFS > 0.0) {
        fprintf(ctx->PSFd,"/fsiz %5.2f def FT\n",1.5 * ctx->PMFS * ctx->PtMM);
        makefmt(ctx, &ctx->PMFmt1,&ctx->PMFmt2,ctx->PMFmtS,sizeof(ctx->PMFmtS),0,ctx->SEPC,1);
         
        if (typ) {
            fprintf(ctx->PSFd,"/lpt %5.2f def\n",1.5 * ctx->PMFS * ctx->PtMM);
            fprintf(ctx->PSFd,"/xmin 0 def\n");
        }
    }
    if (opt1 == 0) {
        ps_ltyp(ctx, ctx->PMLT);
        ps_lwidth(ctx, ctx->PMLW);
        if (ctx->LWCSMax < ctx->PMLW)
            ctx->LWCSMax = ctx->PMLW;
    }
    else {
        ps_ltyp(ctx, lt);
        ps_lwidth(ctx, lw);
        if (ctx->LWCSMax < lw)
            ctx->LWCSMax = lw;
    }

    ps_2dplot(ctx, xa,ya,0);
    ps_2dplot(ctx, xb,yb,1);
    fprintf(ctx->PSFd,"stroke\n");

    if ((ctx->PSLog[typ] == 0 && ctx->PMSC > 0.0) || (ctx->PSLog[typ] && ctx->PMSC > 1.0)) {
        len = 1;
        n = (int)ctx->PMIC;
        x = xa;
        y = ya;
        if (ctx->PMTL <= 0.0)    /* default tick length: 1.8 mm */
            ctx->PMTL = ctx->PTLen;

        if (n > 0) {
            if (!ctx->PSLog[typ])
                d = ctx->PMSC / (double)n;
            else
                d = (ctx->PMSC - 1.0) / (double)n;
        }

        if (typ == 0) {         /* x axis */

            if (ctx->PMDIR) {
                tlen = ctx->PMTL;
                flen = 0.5;
            }   
            else {
                tlen = -ctx->PMTL;
                flen = -1.5;
            }
            while (x <= xb + ctx->EPSI1) {

                if (fabs(x) <= ctx->EPSI1)
                    x = 0.0;

                ps_2dplot(ctx, x,y,0);
                fprintf(ctx->PSFd,"0 %5.2f rl\n",tlen * ctx->PtMM);
                fprintf(ctx->PSFd,"%%#stroke\n");

                if (ctx->PMFS > 0.0) {
                    xx = x + ctx->PMOFF;
                    fprintf(ctx->PSFd,"%%#text: %5.2f %5.2f %5.2f 1 0 %g\n",ps_2dx(ctx, x),ps_2dy(ctx, y) + (flen * ctx->PMFS + tlen) * ctx->PtMM,1.5 * ctx->PMFS * ctx->PtMM,xx);
    
                    fprintf(ctx->PSFd,"gsave\n0 %5.2f rm\n",flen * ctx->PMFS * ctx->PtMM);
                    rt_fprintf_d(ctx, ctx->PSFd,ctx->PMFmtS,xx);
                    fprintf(ctx->PSFd,"center\nshow\ngrestore\n");
                    rt_snprintf_d(buf,sizeof(buf),ctx->PMFmtS,xx);
                    if ((size_t)len < strlen(buf))
                        len = (int)(strlen(buf));
                }
                fprintf(ctx->PSFd,"stroke\n");

                if (n > 0) {
                    for (i = 1; i < n; ++i) {
                        if (!ctx->PSLog[typ])
                            s = x + d * (double)i;
                        else
                            s = x + x * d * (double)i;

                        if (s >= xb + ctx->EPSI1)
                            break;

                        ps_2dplot(ctx, s,y,0);
                        fprintf(ctx->PSFd,"0 %5.2f rl\n",0.7 * tlen * ctx->PtMM);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                }
                if (ctx->PSLog[typ] == 0)
                    x += ctx->PMSC;
                else
                    x *= ctx->PMSC;
            }

            /* update bounding box */

            x = ps_2dx(ctx, xa);                 
            y = ps_2dy(ctx, ya) + 2.0 * tlen * ctx->PtMM;
            if (ctx->PMFS > 0.0)
                y += 2.0 * flen * ctx->PMFS * ctx->PtMM;
            upd_bbox(ctx, 1,x,y);
            y = ps_2dy(ctx, yb) + 2.0 * tlen * ctx->PtMM;
            if (ctx->PMFS > 0.0)
                y += 2.0 * (flen + 1.0) * ctx->PMFS * ctx->PtMM;
            upd_bbox(ctx, 1,x,y);

            if (ctx->PMFS > 0.0) {
                len -= 3;
                if (len < 1)
                    len = 1;
                len = (len + 1) / 2;
                upd_bbox(ctx, 1,x - 2.0 * (double)len * ctx->PMFS * ctx->PtMM,y);
                x = ps_2dx(ctx, xb);                 
                upd_bbox(ctx, 1,x + 2.0 * (double)len * ctx->PMFS * ctx->PtMM,y);
            }
        }
        else {      /* y axis */

            if (ctx->PMDIR) {
                tlen = ctx->PMTL;
                flen = 0.6;
            }   
            else {
                tlen = -ctx->PMTL;
                flen = -0.5;
            }
            while (y <= yb + ctx->EPSI1) {

                if (fabs(y) <= ctx->EPSI1)
                    y = 0.0;

                ps_2dplot(ctx, x,y,0);
                fprintf(ctx->PSFd,"%5.2f 0 rl\n",tlen * ctx->PtMM);
                fprintf(ctx->PSFd,"%%#stroke\n");

                if (ctx->PMFS > 0.0) {
                    yy = y + ctx->PMOFF;
                    if (ctx->PMDIR == 0) {
                        fprintf(ctx->PSFd,"%%#text: %5.2f %5.2f %5.2f 2 0 %g\n",
                        (flen * ctx->PMFS + tlen) * ctx->PtMM,ps_2dy(ctx, y) - 0.3 * ctx->PMFS * ctx->PtMM,1.5 * ctx->PMFS * ctx->PtMM,yy);

                        fprintf(ctx->PSFd,"gsave\n%5.2f 0 rm\n",flen * ctx->PMFS * ctx->PtMM);
                    }
                    else {
                        fprintf(ctx->PSFd,"%%#text: %5.2f %5.2f %5.2f 0 0 %g\n",
                        ps_2dx(ctx, x) + (flen * ctx->PMFS + tlen) * ctx->PtMM,ps_2dy(ctx, y) - 0.3 * ctx->PMFS * ctx->PtMM,1.5 * ctx->PMFS * ctx->PtMM,yy);

                        fprintf(ctx->PSFd,"gsave\n%5.2f %5.2f rm\n",
                                     flen * ctx->PMFS * ctx->PtMM,-0.4 * ctx->PMFS * ctx->PtMM);
                    }
                    rt_fprintf_d(ctx, ctx->PSFd,ctx->PMFmtS,yy);
                    if (ctx->PMDIR == 0)
                        fprintf(ctx->PSFd,"adjust\n");
                    fprintf(ctx->PSFd,"show\ngrestore\n");
                    rt_snprintf_d(buf,sizeof(buf),ctx->PMFmtS,yy);
                    if ((size_t)len < strlen(buf))
                        len = (int)(strlen(buf));
                }
                fprintf(ctx->PSFd,"stroke\n");

                if (n > 0) {
                    for (i = 1; i < n; ++i) {
                        if (!ctx->PSLog[typ])
                            s = y + d * (double)i;
                        else
                            s = y + y * d * (double)i;

                        if (s >= yb + ctx->EPSI1)
                            break;

                        ps_2dplot(ctx, x,s,0);
                        fprintf(ctx->PSFd,"%5.2f 0 rl\n",0.7 * tlen * ctx->PtMM);
                        fprintf(ctx->PSFd,"stroke\n");
                    }
                }
                if (ctx->PSLog[typ] == 0)
                    y += ctx->PMSC;
                else
                    y *= ctx->PMSC;
            }

            /* update bounding box */

            x = ps_2dx(ctx, xa) + 2.0 * tlen * ctx->PtMM;   
            y = ps_2dy(ctx, ya);                 
            if (ctx->PMFS > 0.0) {
                x += 2.0 * flen * ctx->PMFS * ctx->PtMM;
                len -= 3;
                if (len < 1)
                    len = 1;
                if (ctx->PMDIR == 0)
                    x -= 2.0 * (double)len * ctx->PMFS * ctx->PtMM;
                else             
                    x += 2.0 * (double)len * ctx->PMFS * ctx->PtMM;
            }
            upd_bbox(ctx, 1,x,y);

            x = ps_2dx(ctx, xb) + 2.0 * tlen * ctx->PtMM;
            if (ctx->PMFS > 0.0) {
                x += 2.0 * flen * ctx->PMFS * ctx->PtMM;
                len -= 3;
                if (len < 1)
                    len = 1;
                if (ctx->PMDIR == 0)
                    x -= 2.0 * (double)len * ctx->PMFS * ctx->PtMM;
                else             
                    x += 2.0 * (double)len * ctx->PMFS * ctx->PtMM;
            }
            upd_bbox(ctx, 1,x,y);

            if (ctx->PMFS > 0.0) {
                y = ps_2dy(ctx, ya) - 2.0 * ctx->PMFS * ctx->PtMM;   
                upd_bbox(ctx, 1,x,y);
                y = ps_2dy(ctx, yb) + 2.0 * ctx->PMFS * ctx->PtMM;   
                upd_bbox(ctx, 1,x,y);
            }
        }
    }
    err = 0;

PLAFin:
    if (opt == 0)
        p_clean(ctx);
    return(err);
}

/*--##----------------------------------------------------------------------*/
/*  psetup3()   Setup 3-d plot environment.                                 */
/*                                                                          */
/*              psetup3(                                                    */
/*                  view=lon,lat,     direction of projection, def. 30,30   */
/*                  pxa=xa,xb,        logical x axis                        */
/*                  pya=ya,yb,        logical y axis                        */
/*                  pza=za,zb,        logical z axis                        */
/*                  pxlen=...,        horizontal size of plot (mm), def 100 */
/*                  psorg=...,        origin of PostScript CS, def. 100,100 */
/*                  psscal=...,       scaling, def. 1,1                     */
/*                  psrot=...,        rotation, def. 0                      */
/*              );                                                          */
/*                                                                          */
/*              -180 <= lon <= 180, -90 <= lat <= 90                        */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int psetup3(TDAContext *ctx)
{
    int n,m,xs,ys,zs;  
    double x,y;
    register char *p;

    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->PSFFlg == 0) {
        printf1(ctx, "Error: need a PostScript file.\n");
        return(-1);
    }
    ctx->PSLog[0] = ctx->PSLog[1] = 0;                            
    xs = ys = zs = 0;

    ctx->XOrg = 100;         /* origin of PostScript figure */
    ctx->YOrg = 100;  
    ctx->SCALXFac = 1.0;     /* PostScript x scaling factor */
    ctx->SCALYFac = 1.0;     /* PostScript y scaling factor */
    ctx->ROTFac = 0.0;       /* PostScript rotation factor */

    ctx->PXLen = 100.0;      /* phys length of X axis in mm */

    ctx->PSLon = ctx->PSLat = 30.0;   /* direction of projection */

    p = ctx->CmdBuf + 7;
    while (*++p) {
        if (sscanf(p,"pxlen=%lg",&x) == 1 && x > 0) {
            ctx->PXLen = x;
            p = skip_dbl(ctx, p + 6);
        }
        /***
        else if (sscanf(p,"pylen=%lg",&y) == 1 && y > 0) {
            PYLen = y;
            p = skip_dbl(ctx, p + 6);
        }
        ***/
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            ctx->XOrg = n;
            ctx->YOrg = m;
            p = skip_int(ctx, p + 6);
            p = skip_int(ctx, p + 1);
        }
        else if (sscanf(p,"psscal=%lg,%lg",&x,&y) == 2 && x > 0.0 && y > 0.0) {
            ctx->SCALXFac = x;
            ctx->SCALYFac = y;
            p = skip_dbl(ctx, p + 7);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"psrot=%lg",&x) == 1 && x >= 0.0 && x < 360.0) {
            ctx->ROTFac = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"pxa=%lg,%lg",&x,&y) == 2 && x < y) {   
            ctx->PX3D[0] = x;
            ctx->PX3D[1] = y;
            ctx->PX3C = (x + y) / 2.0;
            xs = 1;
            p = skip_dbl(ctx, p + 4);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"pya=%lg,%lg",&x,&y) == 2 && x < y) {   
            ctx->PY3D[0] = x;
            ctx->PY3D[1] = y;
            ctx->PY3C = (x + y) / 2.0;
            ys = 1;          
            p = skip_dbl(ctx, p + 4);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"pza=%lg,%lg",&x,&y) == 2 && x < y) {   
            ctx->PZ3D[0] = x;
            ctx->PZ3D[1] = y;
            ctx->PZ3C = (x + y) / 2.0;
            zs = 1;        
            p = skip_dbl(ctx, p + 4);
            p = skip_dbl(ctx, p + 1);
        }
        else if (sscanf(p,"view=%lg,%lg",&x,&y) == 2 &&   
                 -180.0 <= x && x <= 180.0 && -90.0 <= y && y <= 90.0) {  
            ctx->PSLon = x;
            ctx->PSLat = y;
            p = skip_dbl(ctx, p + 5);
            p = skip_dbl(ctx, p + 1);
        }
        if (*p != ',' && *p != ')') {
            prn_psfe(ctx, p);
            return(-1);    
        }
        *p = '\0';
    }
    if (xs == 0 || ys == 0 || zs == 0) {
        printf1(ctx, "Error: need pxa, pya and pza.\n");
        return(-1);
    }
    setup_3dsys(ctx);
    return(0);
}

/*--##----------------------------------------------------------------------*/
/*  setup_3dsys()               Set up 3d coordinate system.                */

void setup_3dsys(TDAContext *ctx)
{
    register int i,j,k;
    double x,y;

    setup_3proj(ctx, ctx->PSLon,ctx->PSLat);   /* set up projection matrix */

    printf1(ctx, "Creating a new 3d coordinate system.\n\n");
    printf1(ctx, "                 x axis           y-axis           z axis\n");
    prnchar(ctx, '-',57,1);
    printf1(ctx, "Begin  %16.4f %16.4f %16.4f\n",ctx->PX3D[0],ctx->PY3D[0],ctx->PZ3D[0]);
    printf1(ctx, "End    %16.4f %16.4f %16.4f\n",ctx->PX3D[1],ctx->PY3D[1],ctx->PZ3D[1]);
    printf1(ctx, "Center %16.4f %16.4f %16.4f\n",ctx->PX3C,ctx->PY3C,ctx->PZ3C);

    ps_3dprj(ctx, ctx->PX3D[0],ctx->PY3D[0],ctx->PZ3D[0],&x,&y);
    ctx->PA1[0] = ctx->PA2[0] = x;
    ctx->PA1[1] = ctx->PA2[1] = y;

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            for (k = 0; k < 2; ++k) {
                ps_3dprj(ctx, ctx->PX3D[i],ctx->PY3D[j],ctx->PZ3D[k],&x,&y);
                ctx->PA1[0] = dmin(ctx, ctx->PA1[0],x);
                ctx->PA2[0] = dmax(ctx, ctx->PA2[0],x);
                ctx->PA1[1] = dmin(ctx, ctx->PA1[1],y);
                ctx->PA2[1] = dmax(ctx, ctx->PA2[1],y);
            }
        }
    }
    ctx->UXLen = ctx->PA2[0] - ctx->PA1[0];
    ctx->UYLen = ctx->PA2[1] - ctx->PA1[1];
    ctx->PYLen = ctx->PXLen * ctx->UYLen / ctx->UXLen;
    ctx->PSXLen = ctx->PtMM * ctx->PXLen;
    ctx->PSYLen = ctx->PtMM * ctx->PYLen; 

    printf1(ctx, "\nDirection of projection. Longitude: %lg  Latitude: %lg (degrees)\n",ctx->PSLon,ctx->PSLat);
    printf1(ctx, "Coordinates of projection plane.\n");
    printf1(ctx, "X direction: %16.4f %16.4f  Size: %6.2f (mm)\n",ctx->PA1[0],ctx->PA2[0],ctx->PXLen);
    printf1(ctx, "Y direction: %16.4f %16.4f  Size: %6.2f (mm)\n",ctx->PA1[1],ctx->PA2[1],ctx->PYLen);
    printf1(ctx, "Origin of PostScript coordinates: %d,%d\n",ctx->XOrg,ctx->YOrg);
    printf1(ctx, "Scaling: %g,%g\n",ctx->SCALXFac,ctx->SCALYFac);
    printf1(ctx, "Rotation (degrees): %g\n\n",ctx->ROTFac);

    if (ctx->PSFFlg == 1) {      /* init bounding box parameters */
        ctx->BBLX = ctx->XOrg;
        ctx->BBLY = ctx->YOrg;
        ctx->BBUX = ctx->BBLX + (int)(ctx->PSXLen + 0.5);
        ctx->BBUY = ctx->BBLY + (int)(ctx->PSYLen + 0.5);
    }
    ps_inis(ctx, 1);

    ctx->PSFFlg = 2;     /* set to 2 if setup successful */
    ctx->PS3DFlg = 1;
    ctx->PSPROJ = 0;
}

/*--##----------------------------------------------------------------------*/
/*  setup_3dproj(lon,lat)       Set up 3d projection matrix.                */

void setup_3proj(TDAContext *ctx, double lon,double lat)
{
    lon *= Pi / 180.0;
    lat *= Pi / 180.0;
    ctx->PSPR11 = -sin(lon);
    ctx->PSPR12 =  cos(lon);
    ctx->PSPR13 =  0.0;          
    ctx->PSPR21 = -cos(lon) * sin(lat);
    ctx->PSPR22 = -sin(lon) * sin(lat);
    ctx->PSPR23 =  cos(lat);   
    ctx->PSPR31 =  cos(lon) * cos(lat);
    ctx->PSPR32 =  sin(lon) * cos(lat);
    ctx->PSPR33 =  sin(lat);   
/**
tda_out("PSPR %lg %lg %lg\n",PSPR11,PSPR12,PSPR13);
tda_out("PSPR %lg %lg %lg\n",PSPR21,PSPR22,PSPR23);
tda_out("PSPR %lg %lg %lg\n",PSPR31,PSPR32,PSPR33);
**/

    ctx->PSPI11 = -sin(lon);
    ctx->PSPI12 = -cos(lon) * sin(lat);
    ctx->PSPI13 =  cos(lon) * cos(lat);
    ctx->PSPI21 =  cos(lon);              
    ctx->PSPI22 = -sin(lon) * sin(lat);
    ctx->PSPI23 =  sin(lon) * cos(lat);
    ctx->PSPI31 =  0.0;                  
    ctx->PSPI32 =  cos(lat);               
    ctx->PSPI33 =  sin(lat);   
/**
tda_out("PSPI %lg %lg %lg\n",PSPI11,PSPI12,PSPI13);
tda_out("PSPI %lg %lg %lg\n",PSPI21,PSPI22,PSPI23);
tda_out("PSPI %lg %lg %lg\n",PSPI31,PSPI32,PSPI33);
**/

}
 
