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

/* ------------------------------------------------------------------------ */
/*  functions in t_psf.c                                                    */

int psfile(void);
int psetup(void);
void setup_2dsys(void);
void prn_psfe(char *s);
int ps_close(int opt,int scx,int kflag);
void rot_adjust(void);
int pl_axis(int typ,int opt,int opt1,int lt,double lw);
int psetup3(void);
void setup_3dsys();
void setup_3proj(double lon,double lat);

/*--------------------------------------------------------------------------*/
/*  psfile()    Open PostScript file.                                       */
/*              CmdBuf: psf = fname.                                        */
/*              Return 0 if OK, -1 if error.                                */

int psfile(void)
{
    if (check_cmd(1))
        return(-1);

    if (PSFFlg) {
        printf1("Current PostScript file will be closed.\n");
        ps_close(0,1,0);
    }
    strcpy(PSFName,CmdBuf + 7);

    if (!(PSFd = fopen(PSFName,OPEN_WR))) {
        printf1("Error: can't create: %s\n",PSFName);
        return(-1);
    }
    ps_init();          /* write header to output file */
    PSFFlg = 1;
    PSPROJ = PS3DFlg = PSONUM = 0;
    printf1("Opened new PostScript file: %s\n",PSFName);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  psetup()    Setup plot environment.                                     */
/*              CmdBuf: psetup(...)                                         */
/*              Return 0 if OK, -1 if error.                                */

int psetup(void)
{
    int n,m;           
    double x,y;
    register char *p;

    if (check_cmd(0))
        return(-1);

    if (PSFFlg == 0) {
        printf1("Error: need a PostScript file.\n");
        return(-1);
    }
    PSLog[0] = PSLog[1] = -1;                           

    XOrg = 150;         /* origin of PostScript figure */
    YOrg = 460;  
    SCALXFac = 1.0;     /* PostScript x scaling factor */
    SCALYFac = 1.0;     /* PostScript y scaling factor */
    ROTFac = 0.0;       /* PostScript rotation factor */

    PXLen = 120.0;      /* phys length of X axis in mm */
    PYLen =  80.0;      /* phys length of Y axis in mm */

    p = CmdBuf + 6;
    while (*++p) {
        if (sscanf(p,"pxlen=%lg",&x) == 1 && x > 0) {
            PXLen = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"pylen=%lg",&y) == 1 && y > 0) {
            PYLen = y;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            XOrg = n;
            YOrg = m;
            p = skip_int(p + 6);
            p = skip_int(p + 1);
        }
        else if (sscanf(p,"psscal=%lg,%lg",&x,&y) == 2 && x > 0.0 && y > 0.0) {
            SCALXFac = x;
            SCALYFac = y;
            p = skip_dbl(p + 7);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"psrot=%lg",&x) == 1 && x >= 0.0 && x < 360.0) {
            ROTFac = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"pxa(log=%d)=%lg,%lg",&n,&x,&y) == 3 && x < y &&
            (n == 0 || n == 1)) {      
            PA1[0] = x;
            PA2[0] = y;
            PSLog[0] = n;
            p = skip_int(p + 8);
            p = skip_dbl(p + 2);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"pxa=%lg,%lg",&x,&y) == 2 && x < y) {   
            PA1[0] = x;
            PA2[0] = y;
            PSLog[0] = 0;
            p = skip_dbl(p + 4);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"pya(log=%d)=%lg,%lg",&n,&x,&y) == 3 && x < y &&
            (n == 0 || n == 1)) {      
            PA1[1] = x;
            PA2[1] = y;
            PSLog[1] = n;
            p = skip_int(p + 8);
            p = skip_dbl(p + 2);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"pya=%lg,%lg",&x,&y) == 2 && x < y) {   
            PA1[1] = x;
            PA2[1] = y;
            PSLog[1] = 0;
            p = skip_dbl(p + 4);
            p = skip_dbl(p + 1);
        }
        if (*p != ',' && *p != ')') {
            prn_psfe(p);
            return(-1);    
        }
        *p = '\0';
    }
    if (PSLog[0] < 0 || PSLog[1] < 0) {
        printf1("Error: need logical coordinates (pxa and pya parameters).\n");
        return(-1);
    }
    setup_2dsys();
    return(0);
}

/*--##----------------------------------------------------------------------*/
/*  setup_2dsys(lon,lat)        Set up 2d coordinate system.                */

void setup_2dsys(void)
{

    UXLen = PA2[0] - PA1[0];
    UYLen = PA2[1] - PA1[1];

    PSXLen = PtMM * PXLen;   /* phys. x axis length points */
    PSYLen = PtMM * PYLen;   /* phys. y axis length points */

    printf1("Creating a new 2d coordinate system.\n");
    printf1("Size (width x height in mm): %g x %g\n",PXLen,PYLen);
    printf1("PostScript coordinates. X-axis: %d,%d  Y-axis: %d,%d\n",
                    XOrg,XOrg + (int)(PSXLen + 0.5),YOrg,YOrg + (int)(PSYLen + 0.5));   
    printf1("Origin (PostScript x and y coordinates): %d,%d\n",XOrg,YOrg);
    printf1("Scaling (x and y axis): %g,%g\n",SCALXFac,SCALYFac);
    printf1("Rotation (degrees): %g\n\n",ROTFac);

    if (PSFFlg == 1) {      /* init bounding box parameters */
        BBLX = XOrg - 5;
        BBLY = YOrg - 5;
        BBUX = BBLX + (int)PSXLen + 10.0;
        BBUY = BBLY + (int)PSYLen + 10.0;
    }
    ps_inis(1);

    printf1("User coordinates X axis: %16.4f, %16.4f ",PA1[0],PA2[0]);
    if (PSLog[0])
        printf1("(logarithmic)");
    printf1("\nUser coordinates Y axis: %16.4f, %16.4f ",PA1[1],PA2[1]);
    if (PSLog[1])
        printf1("(logarithmic)");
    printf1("\n");

    PSFFlg = 2;              /* set to 2 if setup successful */
    PSPROJ = PS3DFlg = 0;    /* 2-d coordinate system */
}

/* ------------------------------------------------------------------------ */
/*  prn_nve.    Print error message.                                        */

void prn_psfe(char *s)
{
    register char *p = s;

    printf1("Syntax error: ");
    if (!*p)
        printf1("check brackets and semicolon.\n");
    else {     
        while (*p && p < s + 20)
            printf1("%c",*p++);
        if (*p)
            printf1(" ...");
        printf1("\n");
    }
}

/*--------------------------------------------------------------------------*/
/*  ps_close(opt,scx,kflag)  Close PostScript file.                         */
/*                           opt != 0 make message.                         */
/*                           if scx also free scx.                          */
/*                           if kflag != 0 keep parameters.                 */
/*                           Return 0 if OK, -1 if error.                   */

int ps_close(int opt,int scx,int kflag)
{
    int err = 0;                
    long pos;
    char buf[132];
         
    err = 0;
    if (PSFFlg) {
        if (opt)  
            printf1("> Closed current PostScript file: %s\n",PSFName);

        fprintf(PSFd,"showpage\n");
        fclose(PSFd);
    }
    if (PSFFlg == 2) {     /* reprint bounding box */

        if (!(PSFd = fopen(PSFName,OPEN_RU)))  
            err = -1;
        else {

            if (PS3DFlg) {
                BBLX -= 4;
                BBLY -= 4;
                BBUX += 4;
                BBUY += 4;
            }
            if (ROTFac != 0.0) {        /* adjust for rotation */
                rot_adjust();
            }

            while (fgets(buf,131,PSFd)) {
                pos = ftell(PSFd);
                if (!strncmp(buf + 2,"CreationDate",12)) {
                    fflush(PSFd);
                    fseek(PSFd,pos,0);
                    fprintf(PSFd,"%%%%BoundingBox: %5d %5d %5d %5d\n",BBLX,BBLY,BBUX,BBUY);
                    break;
                }
            }
            fclose(PSFd);
        }
        if (err)  
            printf1("\nError: can't update bounding box information in %s.\n",PSFName);
    }
    if (kflag == 0) {
        PSONUM = PSYOrg = PSXOrg = BBLX = BBLY = BBUX = BBUY = PSPROJ = PS3DFlg = PSFFlg = 0;
    }
    if (scx)
        alloc_scx(0);
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  rot_adjust.  Adjust bounding box for rotation.                          */

void rot_adjust(void)
{
    double phi,sphi,cphi,xa,ya,xb,yb,x,y,xmin,xmax,ymin,ymax;

    phi = degree_to_arc(ROTFac);
    sphi = sin(phi);
    cphi = cos(phi);

    xa = (double)(BBLX - XOrg);
    xb = (double)(BBUX - XOrg);
    ya = (double)(BBLY - YOrg);
    yb = (double)(BBUY - YOrg);

    x = xa * cphi - ya * sphi;
    y = xa * sphi + ya * cphi;
    xmin = xmax = x;
    ymin = ymax = y;

    x = xa * cphi - yb * sphi;
    y = xa * sphi + yb * cphi;
    xmin = dmin(xmin,x);
    xmax = dmax(xmax,x);
    ymin = dmin(ymin,y);
    ymax = dmax(ymax,y);

    x = xb * cphi - ya * sphi;
    y = xb * sphi + ya * cphi;
    xmin = dmin(xmin,x);
    xmax = dmax(xmax,x);
    ymin = dmin(ymin,y);
    ymax = dmax(ymax,y);

    x = xb * cphi - yb * sphi;
    y = xb * sphi + yb * cphi;
    xmin = dmin(xmin,x);
    xmax = dmax(xmax,x);
    ymin = dmin(ymin,y);
    ymax = dmax(ymax,y);

    BBLX = XOrg + (int)floor(xmin) - 1;
    BBUX = XOrg + (int)ceil(xmax)  + 1;
    BBLY = YOrg + (int)floor(ymin) - 1;
    BBUY = YOrg + (int)ceil(ymax)  + 1;
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

int pl_axis(int typ,int opt,int opt1,int lt,double lw)
{
    int err,n,i,len;
    double xa,xb,ya,yb,x,xx,y,yy,d,s,tlen,flen;
    char buf[100];

    err = -1;

    if (opt == 0) {
        if (check_pcmd(1,2))
            return(-1);

        /* printf1("Plot %c axis in user's coordinate system.\n",'X' + typ); */

        if (parm(CmdBuf + 4,7,0))    /* get parameters */
            goto PLAFin;
    }
    if (PMRHSN > 0 && opt == 0) {
        if (PMRHSN != 4) {
            p_err(-1,1);
            goto PLAFin;
        }
        xa = PMRHSX[0];
        ya = PMRHSX[1];
        xb = PMRHSX[2];
        yb = PMRHSX[3];
    }
    else if (typ == 0) {
        xa = PA1[0];
        ya = PA1[1];
        xb = PA2[0];
        yb = PA1[1];
    }
    else  {
        xa = PA1[0];
        ya = PA1[1];
        xb = PA1[0];
        yb = PA2[1];
    }
    if (opt == 1) {
        if (typ == 0)
            fprintf(PSFd,"\n%%#%d: plxa\n",++PSONUM);
        else
            fprintf(PSFd,"\n%%#%d: plya\n",++PSONUM);
    }
    else
        fprintf(PSFd,"\n%%#%d: %s\n",++PSONUM,CmdBuf);
    fprintf(PSFd,"%%#%c-axis\n",'x' + typ);
    if (PMFS > 0.0) {
        fprintf(PSFd,"/fsiz %5.2f def FT\n",1.5 * PMFS * PtMM);
        makefmt(&PMFmt1,&PMFmt2,PMFmtS,0,SEPC,1);
         
        if (typ) {
            fprintf(PSFd,"/lpt %5.2f def\n",1.5 * PMFS * PtMM);
            fprintf(PSFd,"/xmin 0 def\n");
        }
    }
    if (opt1 == 0) {
        ps_ltyp(PMLT);
        ps_lwidth(PMLW);
        if (LWCSMax < PMLW)
            LWCSMax = PMLW;
    }
    else {
        ps_ltyp(lt);
        ps_lwidth(lw);
        if (LWCSMax < lw)
            LWCSMax = lw;
    }

    ps_2dplot(xa,ya,0);
    ps_2dplot(xb,yb,1);
    fprintf(PSFd,"stroke\n");

    if ((PSLog[typ] == 0 && PMSC > 0.0) || (PSLog[typ] && PMSC > 1.0)) {
        len = 1;
        n = (int)PMIC;
        x = xa;
        y = ya;
        if (PMTL <= 0.0)    /* default tick length: 1.8 mm */
            PMTL = PTLen;

        if (n > 0) {
            if (!PSLog[typ])
                d = PMSC / (double)n;
            else
                d = (PMSC - 1.0) / (double)n;
        }

        if (typ == 0) {         /* x axis */

            if (PMDIR) {
                tlen = PMTL;
                flen = 0.5;
            }   
            else {
                tlen = -PMTL;
                flen = -1.5;
            }
            while (x <= xb + EPSI1) {

                if (fabs(x) <= EPSI1)
                    x = 0.0;

                ps_2dplot(x,y,0);
                fprintf(PSFd,"0 %5.2f rl\n",tlen * PtMM);
                fprintf(PSFd,"%%#stroke\n");

                if (PMFS > 0.0) {
                    xx = x + PMOFF;
                    fprintf(PSFd,"%%#text: %5.2f %5.2f %5.2f 1 0 %g\n",ps_2dx(x),ps_2dy(y) + (flen * PMFS + tlen) * PtMM,1.5 * PMFS * PtMM,xx);
    
                    fprintf(PSFd,"gsave\n0 %5.2f rm\n",flen * PMFS * PtMM);
                    fprintf(PSFd,PMFmtS,xx);
                    fprintf(PSFd,"center\nshow\ngrestore\n");
                    sprintf(buf,PMFmtS,xx);
                    if (len < strlen(buf))
                        len = strlen(buf);
                }
                fprintf(PSFd,"stroke\n");

                if (n > 0) {
                    for (i = 1; i < n; ++i) {
                        if (!PSLog[typ])
                            s = x + d * (double)i;
                        else
                            s = x + x * d * (double)i;

                        if (s >= xb + EPSI1)
                            break;

                        ps_2dplot(s,y,0);
                        fprintf(PSFd,"0 %5.2f rl\n",0.7 * tlen * PtMM);
                        fprintf(PSFd,"stroke\n");
                    }
                }
                if (PSLog[typ] == 0)
                    x += PMSC;
                else
                    x *= PMSC;
            }

            /* update bounding box */

            x = ps_2dx(xa);                 
            y = ps_2dy(ya) + 2.0 * tlen * PtMM;
            if (PMFS > 0.0)
                y += 2.0 * flen * PMFS * PtMM;
            upd_bbox(1,x,y);
            y = ps_2dy(yb) + 2.0 * tlen * PtMM;
            if (PMFS > 0.0)
                y += 2.0 * (flen + 1.0) * PMFS * PtMM;
            upd_bbox(1,x,y);

            if (PMFS > 0.0) {
                len -= 3;
                if (len < 1)
                    len = 1;
                len = (len + 1) / 2;
                upd_bbox(1,x - 2.0 * (double)len * PMFS * PtMM,y);
                x = ps_2dx(xb);                 
                upd_bbox(1,x + 2.0 * (double)len * PMFS * PtMM,y);
            }
        }
        else {      /* y axis */

            if (PMDIR) {
                tlen = PMTL;
                flen = 0.6;
            }   
            else {
                tlen = -PMTL;
                flen = -0.5;
            }
            while (y <= yb + EPSI1) {

                if (fabs(y) <= EPSI1)
                    y = 0.0;

                ps_2dplot(x,y,0);
                fprintf(PSFd,"%5.2f 0 rl\n",tlen * PtMM);
                fprintf(PSFd,"%%#stroke\n");

                if (PMFS > 0.0) {
                    yy = y + PMOFF;
                    if (PMDIR == 0) {
                        fprintf(PSFd,"%%#text: %5.2f %5.2f %5.2f 2 0 %g\n",
                        (flen * PMFS + tlen) * PtMM,ps_2dy(y) - 0.3 * PMFS * PtMM,1.5 * PMFS * PtMM,yy);

                        fprintf(PSFd,"gsave\n%5.2f 0 rm\n",flen * PMFS * PtMM);
                    }
                    else {
                        fprintf(PSFd,"%%#text: %5.2f %5.2f %5.2f 0 0 %g\n",
                        ps_2dx(x) + (flen * PMFS + tlen) * PtMM,ps_2dy(y) - 0.3 * PMFS * PtMM,1.5 * PMFS * PtMM,yy);

                        fprintf(PSFd,"gsave\n%5.2f %5.2f rm\n",
                                     flen * PMFS * PtMM,-0.4 * PMFS * PtMM);
                    }
                    fprintf(PSFd,PMFmtS,yy);
                    if (PMDIR == 0)
                        fprintf(PSFd,"adjust\n");
                    fprintf(PSFd,"show\ngrestore\n");
                    sprintf(buf,PMFmtS,yy);
                    if (len < strlen(buf))
                        len = strlen(buf);
                }
                fprintf(PSFd,"stroke\n");

                if (n > 0) {
                    for (i = 1; i < n; ++i) {
                        if (!PSLog[typ])
                            s = y + d * (double)i;
                        else
                            s = y + y * d * (double)i;

                        if (s >= yb + EPSI1)
                            break;

                        ps_2dplot(x,s,0);
                        fprintf(PSFd,"%5.2f 0 rl\n",0.7 * tlen * PtMM);
                        fprintf(PSFd,"stroke\n");
                    }
                }
                if (PSLog[typ] == 0)
                    y += PMSC;
                else
                    y *= PMSC;
            }

            /* update bounding box */

            x = ps_2dx(xa) + 2.0 * tlen * PtMM;   
            y = ps_2dy(ya);                 
            if (PMFS > 0.0) {
                x += 2.0 * flen * PMFS * PtMM;
                len -= 3;
                if (len < 1)
                    len = 1;
                if (PMDIR == 0)
                    x -= 2.0 * (double)len * PMFS * PtMM;
                else             
                    x += 2.0 * (double)len * PMFS * PtMM;
            }
            upd_bbox(1,x,y);

            x = ps_2dx(xb) + 2.0 * tlen * PtMM;
            if (PMFS > 0.0) {
                x += 2.0 * flen * PMFS * PtMM;
                len -= 3;
                if (len < 1)
                    len = 1;
                if (PMDIR == 0)
                    x -= 2.0 * (double)len * PMFS * PtMM;
                else             
                    x += 2.0 * (double)len * PMFS * PtMM;
            }
            upd_bbox(1,x,y);

            if (PMFS > 0.0) {
                y = ps_2dy(ya) - 2.0 * PMFS * PtMM;   
                upd_bbox(1,x,y);
                y = ps_2dy(yb) + 2.0 * PMFS * PtMM;   
                upd_bbox(1,x,y);
            }
        }
    }
    err = 0;

PLAFin:
    if (opt == 0)
        p_clean();
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

int psetup3(void)
{
    int n,m,xs,ys,zs;  
    double x,y;
    register char *p;

    if (check_cmd(0))
        return(-1);

    if (PSFFlg == 0) {
        printf1("Error: need a PostScript file.\n");
        return(-1);
    }
    PSLog[0] = PSLog[1] = 0;                            
    xs = ys = zs = 0;

    XOrg = 100;         /* origin of PostScript figure */
    YOrg = 100;  
    SCALXFac = 1.0;     /* PostScript x scaling factor */
    SCALYFac = 1.0;     /* PostScript y scaling factor */
    ROTFac = 0.0;       /* PostScript rotation factor */

    PXLen = 100.0;      /* phys length of X axis in mm */

    PSLon = PSLat = 30.0;   /* direction of projection */

    p = CmdBuf + 7;
    while (*++p) {
        if (sscanf(p,"pxlen=%lg",&x) == 1 && x > 0) {
            PXLen = x;
            p = skip_dbl(p + 6);
        }
        /***
        else if (sscanf(p,"pylen=%lg",&y) == 1 && y > 0) {
            PYLen = y;
            p = skip_dbl(p + 6);
        }
        ***/
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            XOrg = n;
            YOrg = m;
            p = skip_int(p + 6);
            p = skip_int(p + 1);
        }
        else if (sscanf(p,"psscal=%lg,%lg",&x,&y) == 2 && x > 0.0 && y > 0.0) {
            SCALXFac = x;
            SCALYFac = y;
            p = skip_dbl(p + 7);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"psrot=%lg",&x) == 1 && x >= 0.0 && x < 360.0) {
            ROTFac = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"pxa=%lg,%lg",&x,&y) == 2 && x < y) {   
            PX3D[0] = x;
            PX3D[1] = y;
            PX3C = (x + y) / 2.0;
            xs = 1;
            p = skip_dbl(p + 4);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"pya=%lg,%lg",&x,&y) == 2 && x < y) {   
            PY3D[0] = x;
            PY3D[1] = y;
            PY3C = (x + y) / 2.0;
            ys = 1;          
            p = skip_dbl(p + 4);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"pza=%lg,%lg",&x,&y) == 2 && x < y) {   
            PZ3D[0] = x;
            PZ3D[1] = y;
            PZ3C = (x + y) / 2.0;
            zs = 1;        
            p = skip_dbl(p + 4);
            p = skip_dbl(p + 1);
        }
        else if (sscanf(p,"view=%lg,%lg",&x,&y) == 2 &&   
                 -180.0 <= x && x <= 180.0 && -90.0 <= y && y <= 90.0) {  
            PSLon = x;
            PSLat = y;
            p = skip_dbl(p + 5);
            p = skip_dbl(p + 1);
        }
        if (*p != ',' && *p != ')') {
            prn_psfe(p);
            return(-1);    
        }
        *p = '\0';
    }
    if (xs == 0 || ys == 0 || zs == 0) {
        printf1("Error: need pxa, pya and pza.\n");
        return(-1);
    }
    setup_3dsys();
    return(0);
}

/*--##----------------------------------------------------------------------*/
/*  setup_3dsys()               Set up 3d coordinate system.                */

void setup_3dsys(void)
{
    register int i,j,k;
    double x,y;

    setup_3proj(PSLon,PSLat);   /* set up projection matrix */

    printf1("Creating a new 3d coordinate system.\n\n");
    printf1("                 x axis           y-axis           z axis\n");
    prnchar('-',57,1);
    printf1("Begin  %16.4f %16.4f %16.4f\n",PX3D[0],PY3D[0],PZ3D[0]);
    printf1("End    %16.4f %16.4f %16.4f\n",PX3D[1],PY3D[1],PZ3D[1]);
    printf1("Center %16.4f %16.4f %16.4f\n",PX3C,PY3C,PZ3C);

    ps_3dprj(PX3D[0],PY3D[0],PZ3D[0],&x,&y);
    PA1[0] = PA2[0] = x;
    PA1[1] = PA2[1] = y;

    for (i = 0; i < 2; ++i) {
        for (j = 0; j < 2; ++j) {
            for (k = 0; k < 2; ++k) {
                ps_3dprj(PX3D[i],PY3D[j],PZ3D[k],&x,&y);
                PA1[0] = dmin(PA1[0],x);
                PA2[0] = dmax(PA2[0],x);
                PA1[1] = dmin(PA1[1],y);
                PA2[1] = dmax(PA2[1],y);
            }
        }
    }
    UXLen = PA2[0] - PA1[0];
    UYLen = PA2[1] - PA1[1];
    PYLen = PXLen * UYLen / UXLen;
    PSXLen = PtMM * PXLen;
    PSYLen = PtMM * PYLen; 

    printf1("\nDirection of projection. Longitude: %lg  Latitude: %lg (degrees)\n",PSLon,PSLat);
    printf1("Coordinates of projection plane.\n");
    printf1("X direction: %16.4f %16.4f  Size: %6.2f (mm)\n",PA1[0],PA2[0],PXLen);
    printf1("Y direction: %16.4f %16.4f  Size: %6.2f (mm)\n",PA1[1],PA2[1],PYLen);
    printf1("Origin of PostScript coordinates: %d,%d\n",XOrg,YOrg);
    printf1("Scaling: %g,%g\n",SCALXFac,SCALYFac);
    printf1("Rotation (degrees): %g\n\n",ROTFac);

    if (PSFFlg == 1) {      /* init bounding box parameters */
        BBLX = XOrg;
        BBLY = YOrg;
        BBUX = BBLX + (int)(PSXLen + 0.5);
        BBUY = BBLY + (int)(PSYLen + 0.5);
    }
    ps_inis(1);

    PSFFlg = 2;     /* set to 2 if setup successful */
    PS3DFlg = 1;
    PSPROJ = 0;
}

/*--##----------------------------------------------------------------------*/
/*  setup_3dproj(lon,lat)       Set up 3d projection matrix.                */

void setup_3proj(double lon,double lat)
{
    lon *= Pi / 180.0;
    lat *= Pi / 180.0;
    PSPR11 = -sin(lon);
    PSPR12 =  cos(lon);
    PSPR13 =  0.0;          
    PSPR21 = -cos(lon) * sin(lat);
    PSPR22 = -sin(lon) * sin(lat);
    PSPR23 =  cos(lat);   
    PSPR31 =  cos(lon) * cos(lat);
    PSPR32 =  sin(lon) * cos(lat);
    PSPR33 =  sin(lat);   
/**
printf("PSPR %lg %lg %lg\n",PSPR11,PSPR12,PSPR13);
printf("PSPR %lg %lg %lg\n",PSPR21,PSPR22,PSPR23);
printf("PSPR %lg %lg %lg\n",PSPR31,PSPR32,PSPR33);
**/

    PSPI11 = -sin(lon);
    PSPI12 = -cos(lon) * sin(lat);
    PSPI13 =  cos(lon) * cos(lat);
    PSPI21 =  cos(lon);              
    PSPI22 = -sin(lon) * sin(lat);
    PSPI23 =  sin(lon) * cos(lat);
    PSPI31 =  0.0;                  
    PSPI32 =  cos(lat);               
    PSPI33 =  sin(lat);   
/**
printf("PSPI %lg %lg %lg\n",PSPI11,PSPI12,PSPI13);
printf("PSPI %lg %lg %lg\n",PSPI21,PSPI22,PSPI23);
printf("PSPI %lg %lg %lg\n",PSPI31,PSPI32,PSPI33);
**/

}
 
