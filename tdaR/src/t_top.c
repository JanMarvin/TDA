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
#include "tda_context.h"

/*  functions in t_top.c */

int sdrel(TDAContext *ctx);
int sdrel_edge(TDAContext *ctx, int nit,int ni,double *xi,double *yi,    int njt,int nj,double *xj,double *yj,int opt);



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

int sdrel(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,nn,ni,nj,nit,njt,typ,r,idi,idj,cflag,iflag;

    err = -1;
    printf1(ctx, "Elementary topological relations. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto SDRELFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
        ctx->PMOPT = 1;

    iflag = 0;      
    if (ctx->PMOPT == 2)
        iflag = 1;

    nrec = 0;

    /* create node list */

    if (alloc_acn(ctx, ctx->NOC + 1))
        goto SDRELFin;

    nn = 0; /* number of nodes */

    for (i = 0; i < ctx->NOC; ++i) {
        if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
            continue;
        ctx->AcN[nn++] = i;
    }
    if (nn == 0) {
        printf1(ctx, "Error: node list is empty.\n");
        goto SDRELFin;
    }
#if 0   /* leftover debug output, disabled -- see README */
tda_out("nn=%d\n",nn);
#endif

 
    /* allocate arrays for first object */

    if (alloc_acx(ctx, ctx->SDVarMax + 3))
        goto SDRELFin;
    if (alloc_acy(ctx, ctx->SDVarMax + 3))
        goto SDRELFin;

    /* create edge list */

    for (i = 0; i < nn; ++i) {
        if ((nit = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || nit > 3)
            continue;
        idi = (int)get_data(ctx, ctx->SDVarSDID,i);                        

        if (nit == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((ni = sd_getdata(ctx, i,0,cflag,1)) < 1)
            goto SDRELFin;             

        for (j = 0; j < ni; ++j) {
            ctx->AcX[j] = ctx->SDVarX[j];
            ctx->AcY[j] = ctx->SDVarY[j];
        }

        for (j = 0; j < nn; ++j) {
            if ((njt = (int)get_data(ctx, ctx->SDVarSDTyp,j)) < 1 || njt > 3)
                continue;

            if (j == i)
                continue;

            if (njt == 3)
                cflag = 1;
            else
                cflag = 0;

            if ((nj = sd_getdata(ctx, j,0,cflag,1)) < 1)
                goto SDRELFin;             

            r = sdrel_edge(ctx, nit,ni,ctx->AcX,ctx->AcY,njt,nj,ctx->SDVarX,ctx->SDVarY,iflag);
            if (r < 0)
                goto SDRELFin;

            if (r == 1) {
                idj = (int)get_data(ctx, ctx->SDVarSDID,j);                        

                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j + 1);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idi);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idj);
                fprintf(ctx->PMFd,"\n");
                nrec++;
            }
        }

    }

/**

        if ((typ  = (int)get_data(ctx, SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

               
            if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
                goto SDRELFin;
**/          
        

    printf1(ctx, "\n%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDRELFin:
    p_clean(ctx);
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

int sdrel_edge(TDAContext *ctx, int nit,int ni,double *xi,double *yi,    int njt,int nj,double *xj,double *yj,int opt)
{
    register int i,j;
    int r,r1,iflag,tflag;
    double x0,y0,xa,ya,xb,yb,xxa,yya,xxb,yyb;

#if 0   /* leftover debug output, disabled -- see README */
tda_out("\n nit=%d ni=%d njt=%d nj=%d\n",nit,ni,njt,nj);
tda_out("FIRST: ");
for (i = 0; i < ni; ++i)
tda_out("%g,%g ",xi[i],yi[i]);
tda_out("\nSECOND: ");
for (i = 0; i < nj; ++i)
tda_out("%g,%g ",xj[i],yj[i]);
newline(ctx);
#endif


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
                if (x0 < dmin(ctx, xa,xb) || x0 > dmax(ctx, xa,xb) ||
                    y0 < dmin(ctx, ya,yb) || y0 > dmax(ctx, ya,yb))
                    continue;  

                if (opt) {
                    if ((x0 == xa && y0 == ya) || (x0 == xb && y0 == yb))
                        return(0);
                }           
                r = g_intersect(ctx, x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb);
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
                r = g_intersect(ctx, x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb);
                if (r != 0) {
                    if (opt)
                        return(0);
                    else
                        return(1);
                }
            }

            /* check whether (x0,y0) is inside */

            r = g_pol_checkp(ctx, nj - 1,xj,yj,x0,y0);
            return(r);
        }
    }
    else if (nit == 2 || nit == 3) {    /* first object is a polyline,
                                           or a polygon: closed, so the
                                           same containment test applies
                                           (the missing nit == 3 branch
                                           fell through to return(0), so
                                           no polygon was ever reported
                                           inside anything) */

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
                    r = g_intersect(ctx, xa,ya,xb,yb,xj[j-1],yj[j-1],xj[j],yj[j],&xxa,&yya,&xxb,&yyb);
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

                r = g_pol_checkp(ctx, nj - 1,xj,yj,x0,y0);
                if (r == 0 && opt)
                    return(0);

                /* check whether (x0,y0) is on the boundary.  As
                   shipped, this loop returned "not contained" at the
                   first edge the point was NOT on whenever the point
                   was not strictly inside -- so a point lying on a
                   later edge was wrongly rejected.  Decide only after
                   all edges are seen. */

                r1 = 0;
                for (j = 1; j < nj; ++j) {
                    xa = xj[j - 1];
                    ya = yj[j - 1];
                    xb = xj[j];
                    yb = yj[j];
                    if (g_intersect(ctx, x0,y0,x0,y0,xa,ya,xb,yb,&xxa,&yya,&xxb,&yyb)) {
                        r1 = 1;         /* on the boundary */
                        break;
                    }
                }
                if (r1) {
                    if (opt)
                        return(0);
                }
                else if (r == 0)
                    return(0);
            }

            /* all points of first object are inside polygon or on its
               boundary. check for crossing segments. */

            for (i = 1; i < ni; ++i) {
                xa = xi[i - 1];
                ya = yi[i - 1];
                xb = xi[i];
                yb = yi[i];

                for (j = 1; j < nj; ++j) {
                    r = g_intersect(ctx, xa,ya,xb,yb,xj[j-1],yj[j-1],xj[j],yj[j],&xxa,&yya,&xxb,&yyb);
                    if (r == 0)
                        continue;
                    if (opt)            /* strictly inside: any contact */
                        return(0);

                    /* boundaries included: running along the boundary
                       (r == 2) or touching it at a vertex of either
                       segment is not a crossing; only a transversal
                       intersection at an interior point is.  As
                       shipped, any contact was rejected, so an object
                       with an edge on the boundary was never reported
                       contained even with opt = 1. */

                    if (r == 2)
                        continue;
                    if ((xxa == xa && yya == ya) || (xxa == xb && yya == yb) ||
                        (xxa == xj[j-1] && yya == yj[j-1]) ||
                        (xxa == xj[j] && yya == yj[j]))
                        continue;
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
int g_intersect(double xia,double yia,double xib,double yib, double xja,double yja,double xjb,double yjb,double *xxa,double *yya, double *xxb,double *yyb)
**/ 
