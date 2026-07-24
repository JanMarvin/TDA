/****************************************************************************/
/*  t_top                                                                   */
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
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_gm.h"
#include "t_sd.h"

/*  functions in t_top.c */

int sdrel(void);
int sdrel_edge(int nit,int ni,double *xi,double *yi,   
    int njt,int nj,double *xj,double *yj,int opt);



/* ------------------------------------------------------------------------ */
/*  sdrel   Elementary topological relations for spatial data.              */
/*                                                                          */  
/*          sdrel(                                                          */
/*              opt=...,    option, def. 1                                  */
/*                          1 = containment relations (inside)              */
/*                          2 = containment relations (strictly inside)     */
/*              nfmt=...,   integer print format, def. 4                    */
/*                                                                          */
/*          ) = file_name;                                                  */
/*                                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdrel(void)
{
    register int i,j;
    int err,nrec,nn,ni,nj,nit,njt,typ,r,idi,idj,cflag,iflag;

    err = -1;
    printf1("Elementary topological relations. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto SDRELFin;

    if (PMOPT < 1 || PMOPT > 3)
        PMOPT = 1;

    iflag = 0;      
    if (PMOPT == 2)
        iflag = 1;

    nrec = 0;

    /* create node list */

    if (alloc_acn(NOC + 1))
        goto SDRELFin;

    nn = 0; /* number of nodes */

    for (i = 0; i < NOC; ++i) {
        if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
            continue;
        AcN[nn++] = i;
    }
    if (nn == 0) {
        printf1("Error: node list is empty.\n");
        goto SDRELFin;
    }
printf("nn=%d\n",nn);

 
    /* allocate arrays for first object */

    if (alloc_acx(SDVarMax + 3))
        goto SDRELFin;
    if (alloc_acy(SDVarMax + 3))
        goto SDRELFin;

    /* create edge list */

    for (i = 0; i < nn; ++i) {
        if ((nit = (int)get_data(SDVarSDTyp,i)) < 1 || nit > 3)
            continue;
        idi = (int)get_data(SDVarSDID,i);                        

        if (nit == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((ni = sd_getdata(i,0,cflag,1)) < 1)
            goto SDRELFin;             

        for (j = 0; j < ni; ++j) {
            AcX[j] = SDVarX[j];
            AcY[j] = SDVarY[j];
        }

        for (j = 0; j < nn; ++j) {
            if ((njt = (int)get_data(SDVarSDTyp,j)) < 1 || njt > 3)
                continue;

            if (j == i)
                continue;

            if (njt == 3)
                cflag = 1;
            else
                cflag = 0;

            if ((nj = sd_getdata(j,0,cflag,1)) < 1)
                goto SDRELFin;             

            r = sdrel_edge(nit,ni,AcX,AcY,njt,nj,SDVarX,SDVarY,iflag);
            if (r < 0)
                goto SDRELFin;

            if (r == 1) {
                idj = (int)get_data(SDVarSDID,j);                        

                fprintf(PMFd,PMNFmtS,i + 1);
                fprintf(PMFd,PMNFmtS,j + 1);
                fprintf(PMFd,PMNFmtS,idi);
                fprintf(PMFd,PMNFmtS,idj);
                fprintf(PMFd,"\n");
                nrec++;
            }
        }

    }

/**

        if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

               
            if ((n = sd_getdata(i,0,0,1)) < 1)
                goto SDRELFin;
**/          
        

    printf1("\n%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDRELFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdrel_edge(nit,ni,xi,yi,njt,nj,xj,yj,opt)                               */
/*                                                                          */
/*  Return 1 if object i is contained in object j, otherwise return 0.      */
/*  Note: for polygons it is assumed that there is a last point that        */
/*  equals the first one.                                                   */
/*                                                                          */
/*  If opt == 0 include boundaries, otherwise strictly inside.              */
/*                                                                          */

int sdrel_edge(int nit,int ni,double *xi,double *yi,   
    int njt,int nj,double *xj,double *yj,int opt)
{
    register int i,j;
    int r,r1,iflag,tflag;
    double x0,y0,xa,ya,xb,yb,xxa,yya,xxb,yyb;

printf("\n nit=%d ni=%d njt=%d nj=%d\n",nit,ni,njt,nj);
printf("FIRST: ");
for (i = 0; i < ni; ++i)
printf("%g,%g ",xi[i],yi[i]);
printf("\nSECOND: ");
for (i = 0; i < nj; ++i)
printf("%g,%g ",xj[i],yj[i]);
newline();


    if (nit == 1) {         /* first object is a point */
        x0 = xi[0];
        y0 = yi[0];

        if (njt == 1) {     /* second object is a point */
            if (opt)
                return(0);
                       
            if (x0 == xj[0] && y0 == yj[0])
                return(1);
            return(0);
        }
        else if (njt == 2) {                /* second object is a polyline */
            for (j = 1; j < nj; ++j) {
                xa = xj[j - 1];
                ya = yj[j - 1];
                xb = xj[j];
                yb = yj[j];
                if (x0 < dmin(xa,xb) || x0 > dmax(xa,xb) ||
                    y0 < dmin(ya,yb) || y0 > dmax(ya,yb))
                    continue;  

                if (opt) {
                    if ((x0 == xa && y0 == ya) || (x0 == xb && y0 == yb))
                        return(0);
                }           
                r = g_intersect(x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb);
                if (r != 0)  
                    return(1);
            }
            return(0);
        }
        else {                      /* second object is a polygon */

            /* check whether (x0,y0) is on the boundary */

            for (j = 1; j < nj; ++j) {
                xa = xj[j - 1];
                ya = yj[j - 1];
                xb = xj[j];
                yb = yj[j];
                r = g_intersect(x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb);
                if (r != 0) {
                    if (opt)
                        return(0);
                    else
                        return(1);
                }
            }

            /* check whether (x0,y0) is inside */

            r = g_pol_checkp(nj - 1,xj,yj,x0,y0);
            return(r);
        }
    }
    else if (nit == 2) {        /* first object is a polyline */

        if (njt == 1)           /* second is point */
            return(0);

        else if (njt == 2) {    /* second is polyline */

            if (opt)
                return(0);

            iflag = 0;
            for (i = 1; i < ni; ++i) {
                xa = xi[i - 1];
                ya = yi[i - 1];
                xb = xi[i];
                yb = yi[i];

                if (xa == xb && ya == yb)
                    tflag = 1;
                else
                    tflag = 0;

                iflag = 0;
                for (j = 1; j < nj; ++j) {
                    r = g_intersect(xa,ya,xb,yb,xj[j-1],yj[j-1],xj[j],yj[j],&xxa,&yya,&xxb,&yyb);
                    if (r == 0)
                        continue;
                    if (tflag == 1) {
                        iflag = 1;      /* i inside j */
                        break;
                    }
                    if (r != 2)
                        continue;

                    if (xxa == xa && yya == ya && xxb == xb && yyb == yb) {
                        iflag = 1;
                        break;
                    }
                }
                if (iflag == 0)
                    return(0);
            }
            return(iflag);
        }
        else {                  /* second is polygon */

            /* check whether all points of the polyline are inside  
               of the polygon */

            for (i = 0; i < ni; ++i) {
                x0 = xi[i];
                y0 = yi[i];

                r = g_pol_checkp(nj - 1,xj,yj,x0,y0);
                if (r == 0 && opt)
                    return(0);

                /* check whether (x0,y0) is on the boundary */

                for (j = 1; j < nj; ++j) {
                    xa = xj[j - 1];
                    ya = yj[j - 1];
                    xb = xj[j];
                    yb = yj[j];
                    r1 = g_intersect(x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb);
                    if (r1 == 0) {
                        if (r == 0)
                            return(0);
                    }
                    else {              /* on the boundary */
                        if (opt)
                            return(0);
                    }
                }
            }

            /* all points of first object are inside polygon or on its
               boundary. check for crossing segments. */

            for (i = 1; i < ni; ++i) {
                xa = xi[i - 1];
                ya = yi[i - 1];
                xb = xi[i];
                yb = yi[i];

                for (j = 1; j < nj; ++j) {
                    r = g_intersect(xa,ya,xb,yb,xj[j-1],yj[j-1],xj[j],yj[j],&xxa,&yya,&xxb,&yyb);
                    if (r != 0)
                        return(0);
                }
            }
            return(1);
        }
    }
    return(0);
}



/*  g_intersect(xia,yia,xib,yib,xja,yja,xjb,yjb,xxa,yya,xxb,yyb)            */
/*                                                                          */
/*  Input is two line segments given by [(xia,yia),(xib,yib)] and           */
/*  [(xja,yja),(xjb,yjb)]. The functions determines whether and where the   */
/*  two segments intersect. Return code:                                    */
/*                                                                          */
/*  0  if the segments do not intersect                                     */
/*  1  if the segments intersect in a single points. In this case the       */
/*     coordinates of the intersection point are returned in (xxa,yya).     */
/*  2  if the segments overlap in a segment. In this case the coordinates   */
/*     of the overlap are returned in [(xxa,yya),(xxb,yyb)].                */
/*                                                                          */
/*  The code is based on a proposal by Dan Sunday (April-B 2001 Algorithm). */
/**
int g_intersect(double xia,double yia,double xib,double yib,
    double xja,double yja,double xjb,double yjb,double *xxa,double *yya,
    double *xxb,double *yyb)
**/ 
