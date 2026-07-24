/****************************************************************************/
/*  t_dplot                                                                 */
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

/*  functions in t_dplot.c */

int dplot(void);

/* ------------------------------------------------------------------------ */
/*  dplot    arrange several PostScript plots into a single figure.         */
/*           Return 0 if OK, -1 if error                                    */

#define DPMAXFN  50

int dplot(void)
{
    FILE *fd,*fd1;
    register int i,j; 
    int err,n,m,fn,fflag;
    int x,x1,x2,y,y1,y2,xorg1,yorg1,gxorg,gyorg,xlen,ylen,dx,dy,bx,by;
    register char *p;
    char *ofn,*fnp[DPMAXFN][DPMAXFN],tbuf[150];
    int fnn[DPMAXFN],xorg[DPMAXFN][DPMAXFN],yorg[DPMAXFN][DPMAXFN];
    double xscal,yscal,tmp;

    err = -1;
    if (check_cmd(2))
        return(-1);
   
    printf1("Combining PostScript files. Current memory: %d bytes.\n",MemReq);

    fflag = 0;
    xlen = 120;
    ylen =  80;
    gxorg = gyorg = 10;

    p = CmdBuf + 5;
    if (*p != '(')
        goto DPFin;

    ofn = skip_blev(p);
    if (*ofn++ != '=' || !*ofn)
        goto DPFin;

    if (!(fd = fopen(ofn,OPEN_WR))) {
        printf1("Error: can't open %s\n",ofn);
        goto DPFin;
    }
    fflag = 1;

    fn = 0;
         
    while (*p) {
        if (*p == ',' || *p == '(')
            p++;
        if (*p == ')')
            break;

        if (!strncmp(p,"fn=",3)) {
            if (fn >= DPMAXFN) {
                err = -2;
                goto DPFin;
            }
            i = 0;
            p += 3;
            while (i < DPMAXFN) {
                fnp[fn][i++] = p;
                p = skip_nc(p);
                if (*p == ')') {
                    *p = '\0';
                    break;
                }
                if (*p != ',')
                    goto DPFin;
                *p++ = '\0';
                if (*p == ')' || !strncmp(p,"fn=",3) || !strncmp(p,"pxlen",5) || !strncmp(p,"pylen",5) || !strncmp(p,"psorg",5)) {
                    break;
                }                   
            }
            fnn[fn] = i;
            fn++;
        }
        else if (sscanf(p,"pxlen=%d",&n) == 1 && n > 0) {
            xlen = n;
            p = skip_int(p + 6);
        }
        else if (sscanf(p,"pylen=%d",&n) == 1 && n > 0) {
            ylen = n;
            p = skip_int(p + 6);
        }
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            gxorg = n;
            gyorg = m;
            p = skip_int(p + 6);
            p = skip_int(p + 1);
        }
        else
            goto DPFin;
    }
    err = -3;
    if (fn == 0) {
        printf1("Error: no input files.\n");
        goto DPFin;
    }
    bx = by = yorg1 = 0; 

    for (i = fn - 1; i >= 0; --i) {

        printf1("\nRow %d\n",fn - i);

        xorg1 = 0;    
        dx = dy = 0;

        for (j = 0; j < fnn[i]; ++j) {

            printf1("%s  ",fnp[i][j]);

            if (!(fd1 = fopen(fnp[i][j],OPEN_RD))) {
                printf1("--  can't open this file.\n");
                goto DPFin;
            }
            n = 0;
            while (n < 3 && fgets(tbuf,130,fd1)) {
                if (sscanf(tbuf,"%%%%BoundingBox: %d %d %d %d",
                                              &x1,&y1,&x2,&y2) == 4) {
                    n++;
                }
                else if (sscanf(tbuf,"/xorg %d",&x) == 1)
                    n++; 
                else if (sscanf(tbuf,"/yorg %d",&y) == 1)
                    n++; 
            }
            fclose (fd1);
            if (n != 3) {
                printf1("--  can't find bounding box, xorg or yorg coordinates.\n");
                goto DPFin;
            }
            printf1("bounding box: %3d %3d %3d %3d",x1,y1,x2,y2);
            printf1("  xorg%4d  yorg%4d\n",x,y);

            dx += x2 - x1;
            if (dy < y2 - y1)
                dy = y2 - y1;
         
            xorg[i][j] = xorg1 + x - x1;
            yorg[i][j] = yorg1 + y - y1;
            xorg1 += x2 - x1;
        }
        if (bx < dx)
            bx = dx;
        by += dy;
        yorg1 += dy;
    }
    if (bx <= 0.0 || by <= 0.0)
        goto DPFin;

    xscal = ((double)xlen * PtMM)/ (double)bx;
    yscal = ((double)ylen * PtMM)/ (double)by;

    bx    = gxorg + (int)(bx * xscal);
    by    = gyorg + (int)(by * yscal);
  
    printf1("\nNew bounding box: %d %d %d %d",gxorg,gyorg,bx,by);
    printf1("\nScaled with xscale = %g  yscale = %g\n",xscal,yscal);
             
    /****************************************************************
       now read the files again and write with modifications into
       the output file.    

       default header:

       %!PS-Adobe-2.0 EPSF-1.2
       %%BoundingBox: 124 442 556 751

       /xorg 150 def  % X and Y coordinates of origin
       /yorg 460 def
    ***************************************************************/
   
    printf1("\nCopying input file(s) into output file: %s\n",ofn);

    fprintf(fd,"%%!PS-Adobe-2.0 EPSF-1.2\n");
    fprintf(fd,"%%%%BoundingBox: %d %d %d %d\n\n",gxorg,gyorg,bx,by);
    fprintf(fd,"%% changed: global scaling\n");
    fprintf(fd,"/sfx %5.2f def\n",xscal);
    fprintf(fd,"/sfy %5.2f def\n",yscal);

    for (i = fn - 1; i >= 0; --i) {
        for (j = 0; j < fnn[i]; ++j) {

            if (!(fd1 = fopen(fnp[i][j],OPEN_RD))) {
                printf1("Error: can't reopen: %s\n",fnp[i][j]);
                goto DPFin;
            }
            fprintf(fd,"\n%% input file: %s\n",fnp[i][j]);
            fprintf(fd,"gsave %% added\n");

            while (fgets(tbuf,130,fd1)) {
                if (!strncmp(tbuf,"showpage",8)) ;
                else if (!strncmp(tbuf + 1,"!PS-Adobe",9)) ;
                else if (!strncmp(tbuf + 2,"BoundingBox",11)) ;
                else if (!strncmp(tbuf,"/sfx",3))  ;
                else if (!strncmp(tbuf,"/sfy",3))  ;
                else if (!strncmp(tbuf,"/xorg",5)) { 
                    tmp = (double)gxorg + (xscal * (double)xorg[i][j]);
                    fprintf(fd,"/xorg %d def %% changed\n",(int)tmp);
                }
                else if (!strncmp(tbuf,"/yorg",5)) { 
                    tmp = (double)gyorg + (yscal * (double)yorg[i][j]);
                    fprintf(fd,"/yorg %d def %% changed\n",(int)tmp);
                }
                else   
                    fprintf(fd,"%s",tbuf);      
            }
            fprintf(fd,"grestore %% added\n");
            fclose (fd1);
        }
    }
    fprintf(fd,"showpage\n");
    err = 0;

DPFin:  
    if (fflag)
        fclose(fd);
    if (err == -1)
        p_err(-1,1);
    else if (err == -2)
        printf1("Error: exceeded max number of input files.\n");
    return(err);
}

