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
#include "tda_context.h"

/*  functions in t_dplot.c */

int dplot(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  dplot    arrange several PostScript plots into a single figure.         */
/*           Return 0 if OK, -1 if error                                    */

#define DPMAXFN  50

int dplot(TDAContext *ctx)
{
    FILE *fd,*fd1;
    register int i,j; 
    int err,n,m,fn,fflag;
    int x,x1,x2,y,y1,y2,yorg1,gxorg,gyorg,xlen,ylen,bx,by;
    register char *p;
    char *ofn,*fnp[DPMAXFN][DPMAXFN],tbuf[150];
    int fnn[DPMAXFN],xorg[DPMAXFN][DPMAXFN],yorg[DPMAXFN][DPMAXFN];
    int bxo[DPMAXFN][DPMAXFN],byo[DPMAXFN][DPMAXFN];
    int bxl[DPMAXFN][DPMAXFN],byl[DPMAXFN][DPMAXFN];
    int cw,ch,nc;
    double xscal,yscal,tmp;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);
   
    printf1(ctx, "Combining PostScript files. Current memory: %d bytes.\n",ctx->MemReq);

    fflag = 0;
    xlen = 120;
    ylen =  80;
    gxorg = gyorg = 10;

    p = ctx->CmdBuf + 5;
    if (*p != '(')
        goto DPFin;

    ofn = skip_blev(ctx, p);
    if (*ofn++ != '=' || !*ofn)
        goto DPFin;

    if (!(fd = fopen(ofn,OPEN_WR))) {
        printf1(ctx, "Error: can't open %s\n",ofn);
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
                p = skip_nc(ctx, p);
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
            p = skip_int(ctx, p + 6);
        }
        else if (sscanf(p,"pylen=%d",&n) == 1 && n > 0) {
            ylen = n;
            p = skip_int(ctx, p + 6);
        }
        else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
            gxorg = n;
            gyorg = m;
            p = skip_int(ctx, p + 6);
            p = skip_int(ctx, p + 1);
        }
        else
            goto DPFin;
    }
    err = -3;
    if (fn == 0) {
        printf1(ctx, "Error: no input files.\n");
        goto DPFin;
    }
    /****************************************************************
       Read every input file's bounding box and origin first, then
       place them.  The cells are all the size of the largest input,
       so the plots line up in columns and rows whatever their own
       boxes happen to be.  Sizing each cell from its own file --
       what this did before -- makes the columns of one row disagree
       with the next as soon as two inputs differ in width, which
       they do whenever one carries a longer pltext label than
       another: TDA's bounding box reserves room for the text, and a
       neighbour without it is then narrower.  Inputs of equal size,
       as in the manual's own example, give the same layout either
       way.
    ****************************************************************/

    cw = ch = nc = 0;

    for (i = fn - 1; i >= 0; --i) {

        printf1(ctx, "\nRow %d\n",fn - i);

        if (nc < fnn[i])
            nc = fnn[i];

        for (j = 0; j < fnn[i]; ++j) {

            printf1(ctx, "%s  ",fnp[i][j]);

            if (!(fd1 = fopen(fnp[i][j],OPEN_RD))) {
                printf1(ctx, "--  can't open this file.\n");
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
                printf1(ctx, "--  can't find bounding box, xorg or yorg coordinates.\n");
                goto DPFin;
            }
            printf1(ctx, "bounding box: %3d %3d %3d %3d",x1,y1,x2,y2);
            printf1(ctx, "  xorg%4d  yorg%4d\n",x,y);

            bxo[i][j] = x - x1;
            byo[i][j] = y - y1;
            bxl[i][j] = x2 - x1;
            byl[i][j] = y2 - y1;

            if (cw < bxl[i][j])
                cw = bxl[i][j];
            if (ch < byl[i][j])
                ch = byl[i][j];
        }
    }
    if (cw <= 0 || ch <= 0 || nc <= 0)
        goto DPFin;

    yorg1 = 0;
    for (i = fn - 1; i >= 0; --i) {
        for (j = 0; j < fnn[i]; ++j) {
            xorg[i][j] = j * cw + bxo[i][j];
            yorg[i][j] = yorg1 + byo[i][j];
        }
        yorg1 += ch;
    }
    bx = nc * cw;
    by = fn * ch;

    printf1(ctx, "\nCell size: %d %d\n",cw,ch);

    xscal = ((double)xlen * ctx->PtMM)/ (double)bx;
    yscal = ((double)ylen * ctx->PtMM)/ (double)by;

    bx    = gxorg + (int)(bx * xscal);
    by    = gyorg + (int)(by * yscal);
  
    printf1(ctx, "\nNew bounding box: %d %d %d %d",gxorg,gyorg,bx,by);
    printf1(ctx, "\nScaled with xscale = %g  yscale = %g\n",xscal,yscal);
             
    /****************************************************************
       now read the files again and write with modifications into
       the output file.    

       default header:

       %!PS-Adobe-2.0 EPSF-1.2
       %%BoundingBox: 124 442 556 751

       /xorg 150 def  % X and Y coordinates of origin
       /yorg 460 def
    ***************************************************************/
   
    printf1(ctx, "\nCopying input file(s) into output file: %s\n",ofn);

    fprintf(fd,"%%!PS-Adobe-2.0 EPSF-1.2\n");
    fprintf(fd,"%%%%BoundingBox: %d %d %d %d\n\n",gxorg,gyorg,bx,by);
    fprintf(fd,"%% changed: global scaling\n");
    fprintf(fd,"/sfx %5.2f def\n",xscal);
    fprintf(fd,"/sfy %5.2f def\n",yscal);

    for (i = fn - 1; i >= 0; --i) {
        for (j = 0; j < fnn[i]; ++j) {

            if (!(fd1 = fopen(fnp[i][j],OPEN_RD))) {
                printf1(ctx, "Error: can't reopen: %s\n",fnp[i][j]);
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
        p_err(ctx, -1,1);
    else if (err == -2)
        printf1(ctx, "Error: exceeded max number of input files.\n");
    return(err);
}

