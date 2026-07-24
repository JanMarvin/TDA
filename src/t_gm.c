/****************************************************************************/
/*  t_gm                                                                    */
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

/*  functions in t_gm.c */

void g_spat(double ax,double ay,double az,double bx,double by,double bz,
    double *rx,double *ry,double *rz);
int g_intersect(double xia,double yia,double xib,double yib,
    double xja,double yja,double xjb,double yjb,double *xxa,double *yya,
    double *xxb,double *yyb);
int g_inseg(double px,double py,double xa,double ya,double xb,double yb);
double g_perp(double ux,double uy,double vx,double vy);
int g_left(double xa,double ya,double xb,double yb,double x,double y);
int g_ccw(double x0,double y0,double x1,double y1,double x2,double y2);
int g_isect(double x1,double y1,double x2,double y2,double x3,double y3,
    double x4,double y4);
double g_line_len(int n,double *x,double *y);
double g_pol_area(int n,double *x,double *y);
double g_pol_xmin(int n,double *x,double *y);
double g_pol_xmax(int n,double *x,double *y);
double g_pol_ymin(int n,double *x,double *y);
double g_pol_ymax(int n,double *x,double *y);
void g_pol_rect(int n,double *x,double *y,double *xmin,double *xmax,
    double *ymin,double *ymax);
int g_inpoly(int n, double *x,double *y,double xa,double ya);
int g_l_inpoly(int n, double *x,double *y,double xa,double ya,double xb,
    double yb);
int g_p_inpoly(int na, double *xa,double *ya,int nb,double *xb,double *yb);
int g_pol_checkp(int n,double *x,double *y,double px,double py);
int g_pol_isect(int n, double *x,double *y,double xa,double ya,
    double xb,double yb);
int g_chull(float *x,float *y,int m,int *in,int *ih,int *il);
void g_chull_split(float *x,float *y,int m,int *in,int ii,int jj,int s,
    int *iabv,int *na,int *maxa,int *ibel,int *nb,int *maxb);
void g_xchxy(double *x,double *y);
double g_len(double xa,double ya,double xb,double yb);
double g_vlen(double x,double y);
void g_poly_ccw(int n,double *x,double *y);
int g_line_intersect(double p1x,double p1y,double p2x,double p2y,
    double q1x,double q1y,double q2x,double q2y,
    double *s1x,double *s1y,double *s2x,double *s2y,double eps);


/* ------------------------------------------------------------------------ */
/*  g_spat(ax,ay,az,bx,by,bz,rx,ry,rz)                                      */
/*                                                                          */
/*  Return a x b in (rx,ry,rz)                                              */

void g_spat(double ax,double ay,double az,double bx,double by,double bz,
    double *rx,double *ry,double *rz)
{
    *rx = ay * bz - az * by;
    *ry = az * bx - ax * bz;
    *rz = ax * by - ay * bx;
}

/* ------------------------------------------------------------------------ */
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

int g_intersect(double xia,double yia,double xib,double yib,
    double xja,double yja,double xjb,double yjb,double *xxa,double *yya,
    double *xxb,double *yyb)
{
    float ux,uy,vx,vy,wx,wy,d,du,dv,w2x,w2y,t,t0,t1,s1;

    ux = xib - xia;
    uy = yib - yia;
    vx = xjb - xja;
    vy = yjb - yja;
    wx = xia - xja;
    wy = yia - yja;
    d = g_perp(ux,uy,vx,vy);

    /*  test if segments are parallel (includes either being a point) */

    if (fabs(d) < EPSI1) {       /* parallel */
        if (g_perp(ux,uy,wx,wy) != 0.0 || g_perp(vx,vy,wx,wy) != 0.0) /* not collinear */
            return(0);
     
        /* they are collinear or degenerate */

        du = ux * ux + uy * uy;
        dv = vx * vx + vy * vy;

        if (du == 0.0 && dv == 0.0) {       /* both segments are points */
            if (xia != xja || yia != yja)   /* distinct points */
                return(0);
            *xxa = xia;
            *yya = yia;
            return(1);
        }
        if (du == 0.0) {                /* first segment is single point */

            if (g_inseg(xia,yia,xja,yja,xjb,yjb) == 0)
                return(0);

            *xxa = xia;
            *yya = yia;
            return(1);
        }
        if (dv == 0.0) {                /* second segment is single point */

            if (g_inseg(xja,yja,xia,yia,xib,yib) == 0)
                return(0);

            *xxa = xja;
            *yya = yja;
            return(1);
        }

        /* segments are collinear. overlap? */

        w2x = xib - xja;
        w2y = yib - yja;

        if (vx != 0.0) {
            t0 = wx / vx;
            t1 = w2x / vx;
        }
        else {
            t0 = wy / vy;
            t1 = w2y / vy;
        }
        if (t0 > t1) {
            t = t0;
            t0 = t1;
            t1 = t;
        }
        if (t0 > 1.0 || t1 < 0.0)       /* no overlap */
            return(0);

        t0 = dmax(t0,0.0);
        t1 = dmin(t1,1.0);
        if (t0 == t1) {                 /* intersection is a point */
            *xxa = xja + t0 * vx;
            *yya = yja + t0 * vy;
            return(1);
        }

        /* they overlap in a valid segment */

        *xxa = xja + t0 * vx;
        *yya = yja + t0 * vy;
        *xxb = xja + t1 * vx;
        *yyb = yja + t1 * vy;
        return(2);
    }

    /* the segments are skew and may intersect in a point */

    /* get the intersect parameter for first segment */

    s1 = g_perp(vx,vy,wx,wy) / d;
    if (s1 < 0.0 || s1 > 1.0)       /* no intersect */
        return(0);

    /* get the intersect parameter for second segment */

    t1 = g_perp(ux,uy,wx,wy) / d;
    if (t1 < 0.0 || t1 > 1.0)       /* no intersect */
        return(0);

    *xxa = xia + s1 * ux;
    *yya = yia + s1 * uy;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  g_inseq(px,py,xa,ya,xb,yb)                                              */
/*                                                                          */
/*  Given a point (px,py) and a collinear line segment [(xa,ya),(xb,yb)],   */
/*  the function returns 1 if the point is inside the segment, otherwise    */
/*  the function returns 0.                                                 */
/*                                                                          */
/*  The code is based on a proposal by Dan Sunday (April-B 2001 Algorithm). */

int g_inseg(double px,double py,double xa,double ya,double xb,double yb)
{
    if (xa != xb) {     /* segment is not vertical */
        if (xa <= px && px <= xb)
            return(1);
        if (xa >= px && px >= xb)
            return(1);
    }
    else {
        if (ya <= py && py <= yb)
            return(1);
        if (ya >= py && py >= yb)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_perp(ux,uy,vx,vy) Returns the perp product of (ux,uy) and (vx,vy)     */

double g_perp(double ux,double uy,double vx,double vy)
{
    return(ux * vy - uy * vx);
}

/* ------------------------------------------------------------------------ */
/*  g_left(xa,ya,xb,yb,x,y)                                                 */
/*                                                                          */
/*  Return 1 if (x,y) is on or to the left of the segment from (xa,ya)      */
/*  to (xb,yb), otherwise return 0.                                         */  

int g_left(double xa,double ya,double xb,double yb,double x,double y)
{
    if ((xb - xa) * (y - ya) >= (x - xa) * (yb - ya))                       
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_ccw(x0,y0,x1,y1,x2,y2)                                                */
/*                                                                          */
/*  Return 1 if when travelling from the first to the second to the third   */
/*  point, the direction is counterclockwise, otherwise return -1.          */
    
int g_ccw(double x0,double y0,double x1,double y1,double x2,double y2)
{
    double dx1,dx2,dy1,dy2;

    dx1 = x1 - x0;
    dx2 = x2 = x0;
    dy1 = y1 - y0;
    dy2 = y2 - y0;

    if ((dx1 * dy2) > (dy1 * dx2))
        return(1);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  g_isect(x1,y1,x2,y2,x3,y3,x4,y4)                                        */
/*                                                                          */
/*  Given two line segments, (x1,y1,x2,y2) and (x3,y3,x4,y4), return 1 if   */
/*  they intersect, otherwise return 0.                                     */
    
int g_isect(double x1,double y1,double x2,double y2,double x3,double y3,
    double x4,double y4)
{
    if (g_ccw(x1,y1,x2,y2,x3,y3) * g_ccw(x1,y1,x2,y2,x4,y4) <= 0 &&
        g_ccw(x3,y3,x4,y4,x1,y1) * g_ccw(x3,y3,x4,y4,x2,y2) <= 0)       
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_line_len(n,x,y)                                                       */
/*                                                                          */
/*  Given a sequence of points x[i],y[i], i = 0,...,n-1, the function       */
/*  returns the sum of the euclidean length of the segments.                */
 
double g_line_len(int n,double *x,double *y)
{
    register int i;
    double len,xx,yy;

    len = 0.0;
    for (i = 1; i < n; ++i) {
        xx = x[i] - x[i - 1];
        yy = y[i] - y[i - 1];
        len += sqrt(xx * xx + yy * yy);
    }
    return(len);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_area(n,x,y)                                                       */
/*                                                                          */
/*  Given a polygon x[i],y[i] (i = 0,...,n-1), this function returns its    */
/*  area. The value is >= 0 if the polygon is positively oriented           */
/*  (counterclockwise), otherwise the value is <= 0.                        */
 
double g_pol_area(int n,double *x,double *y)
{
    register int i;
    double a,x0,y0,dx1,dy1,dx2,dy2;

    if (n < 3)
        return(0.0);

    a  = 0.0;
    x0 = x[0];
    y0 = y[0];
    dx1 = x[1] - x0;
    dy1 = y[1] - y0;

    for (i = 2; i < n; ++i) {
        dx2 = x[i] - x0; 
        dy2 = y[i] - y0;
        a += dx1 * dy2 - dx2 * dy1;
        dx1 = dx2;
        dy1 = dy2;
    }
    return(a / 2.0);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_xmin(n,x,y)                                                       */
/*                                                                          */
/*  Return xmin for the smallest rectangle that contains the polygon given  */
/*  by x[i],y[i] (i=0,n-1).                                                 */
    
double g_pol_xmin(int n,double *x,double *y)
{
    register int i;
    double xmin;

    xmin = x[0];
    for (i = 1; i < n; ++i) {
        if (x[i] < xmin)
            xmin = x[i];
    }
    return(xmin);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_xmax(n,x,y)                                                       */
/*                                                                          */
/*  Return xmax for the smallest rectangle that contains the polygon given  */
/*  by x[i],y[i] (i=0,n-1).                                                 */
    
double g_pol_xmax(int n,double *x,double *y)
{
    register int i;
    double xmax;

    xmax = x[0];
    for (i = 1; i < n; ++i) {
        if (x[i] > xmax)
            xmax = x[i];
    }
    return(xmax);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_ymin(n,x,y)                                                       */
/*                                                                          */
/*  Return ymin for the smallest rectangle that contains the polygon given  */
/*  by x[i],y[i] (i=0,n-1).                                                 */
    
double g_pol_ymin(int n,double *x,double *y)
{
    register int i;
    double ymin;

    ymin = y[0];
    for (i = 1; i < n; ++i) {
        if (y[i] < ymin)
            ymin = y[i];
    }
    return(ymin);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_ymax(n,x,y)                                                       */
/*                                                                          */
/*  Return ymax for the smallest rectangle that contains the polygon given  */
/*  by x[i],y[i] (i=0,n-1).                                                 */
    
double g_pol_ymax(int n,double *x,double *y)
{
    register int i;
    double ymax;

    ymax = y[0];
    for (i = 1; i < n; ++i) {
        if (y[i] > ymax)
            ymax = y[i];
    }
    return(ymax);
}

/* ------------------------------------------------------------------------ */
/*  g_pol_rect(n,x,y,xmin,xmax,ymin,ymax)                                   */
/*                                                                          */
/*  Calculate the smallest rectangle that contains the polygon given by     */
/*  x[i],y[i] (i=0,n-1).                                                    */
    
void g_pol_rect(int n,double *x,double *y,double *xmin,double *xmax,
    double *ymin,double *ymax)
{
    register int i;

    *xmin = *xmax = x[0];
    *ymin = *ymax = y[0];
    for (i = 1; i < n; ++i) {
        if (x[i] < *xmin)
            *xmin = x[i];
        if (x[i] > *xmax)
            *xmax = x[i];
        if (y[i] < *ymin)
            *ymin = y[i];
        if (y[i] > *ymax)
            *ymax = y[i];
    }
}

/* ------------------------------------------------------------------------ */
/*  g_inpoly(n,x,y,xa,ya)                                                   */
/*                                                                          */
/*  Return 1 if the point (xa,ya) is inside the polygon defined by          */
/*  by x[i],y[i] (i=0,n-1), including its boundary, otherwise return 0.     */  
/*                                                                          */
/*  It is assumed that (x[0],y[0] != x[n-1],y[n-1]. The edge from           */
/*  (x[n-1],y[n-1]) to (x[0],y[0]) is assumed implicitly.                   */
/*                                                                          */
/*  The code is a modified version of the inpoly function written by        */
/*  Bob Stein (Linux Journal, March 1997).                                  */

int g_inpoly(int n, double *x,double *y,double xa,double ya)
{
    register int i,inside;
    double xold,xnew,yold,ynew,x1,x2,y1,y2,tmp1,tmp2;

    if (n < 3)
        return(0);

    inside = 0;
    xold = x[n - 1];
    yold = y[n - 1];
    for (i = 0; i < n; ++i) {
        xnew = x[i];
        ynew = y[i];
        if (xnew > xold) {
            x1 = xold;
            x2 = xnew;
            y1 = yold;
            y2 = ynew;
        }
        else {
            x1 = xnew;
            x2 = xold;
            y1 = ynew;
            y2 = yold;
        }
        tmp1 = (ya - y1) * (x2 - x1);
        tmp2 = (y2 - y1) * (xa - x1);
        if (fabs(tmp1 - tmp2) <= EPSI1)
            return(1);

        if ((xnew < xa) == (xa <= xold) && tmp1 < tmp2)
            inside = !inside;

        xold = xnew;
        yold = ynew;
    }
    return(inside);
} 

/* ------------------------------------------------------------------------ */
/*  g_l_inpoly(n,x,y,xa,ya,xb,yb)                                           */
/*                                                                          */
/*  Return 1 if the line segment [(xa,ya),(xb,yb)] is inside the polygon    */  
/*  defined by x[i],y[i] (i=0,n-1), including its boundary,                 */
/*  otherwise return 0.                                                     */  
/*                                                                          */
/*  It is assumed that (x[0],y[0] != x[n-1],y[n-1]. The edge from           */
/*  (x[n-1],y[n-1]) to (x[0],y[0]) is assumed implicitly.                   */

int g_l_inpoly(int n, double *x,double *y,double xa,double ya,double xb,
    double yb)
{
    register int i;
    double xold,xnew,yold,ynew,xxa,yya,xxb,yyb;

    if (n < 3)
        return(0);

    if (g_inpoly(n,x,y,xa,ya) == 0)
        return(0);
    if (g_inpoly(n,x,y,xb,yb) == 0)
        return(0);

    xold = x[n - 1];
    yold = y[n - 1];
    for (i = 0; i < n; ++i) {
        xnew = x[i];
        ynew = y[i];
        if (g_intersect(xa,ya,xb,yb,xold,yold,xnew,ynew,&xxa,&yya,&xxb,&yyb) == 1) {
            if ((fabs(xa - xxa) > EPSI1 || fabs(ya - yya) > EPSI) &&
                (fabs(xb - xxa) > EPSI1 || fabs(yb - yya) > EPSI) &&
                (fabs(xold - xxa) > EPSI1 || fabs(yold - yya) > EPSI) &&
                (fabs(xnew - xxa) > EPSI1 || fabs(ynew - yya) > EPSI))  
                return(0);
        }
        xold = xnew;
        yold = ynew;
    }
    return(1);
} 

/* ------------------------------------------------------------------------ */
/*  g_p_inpoly(na,xa,ya,nb,xb,yb)                                           */
/*                                                                          */
/*  Return 1 if polygon a is included in polygon b, otherwise return 0.     */ 
/*                                                                          */
/*  It is assumed that (xa[0],ya[0]) = (xa[na-1],ya[na-1]) and              */
/*                     (xb[0],yb[0]) = (xb[nb-1],yb[nb-1]).                 */

int g_p_inpoly(int na, double *xa,double *ya,int nb,double *xb,double *yb)
{
    register int i;

    for (i = 1; i < na; ++i) {
        if (g_l_inpoly(nb - 1,xb,yb,xa[i-1],ya[i-1],xa[i],ya[i]) == 0)
            return(0);
    }
    return(1);
} 

/* ------------------------------------------------------------------------ */
/*  g_pol_checkp(n,x,y,px,py)                                               */
/*                                                                          */
/*  Return 1 if the point (px,py) is inside the polygon defined by          */  
/*  x[i],y[i] (i = 0,n-1), otherwise return 0.                              */
/*                                                                          */  
/*  Note: points that are on the line segments are taken to be outside      */
/*  the polygon.                                                            */
/*                                                                          */
/*  The code follows an algorithm proposed by Randolph Franklin as          */
/*  presented in a discussion by Paul Bourke.                               */
    
int g_pol_checkp1(int n,double *x,double *y,double px,double py)
{
    register int i,j;
    int c = 0;

    for (i = 0, j = n - 1; i < n; j = i++) {
        if ((((y[i] <= py) && (py < y[j])) ||
             ((y[j] <= py) && (py < y[i]))) &&
             (px < (x[j] - x[i]) * (py - y[i]) / (y[j] - y[i]) + x[i]))  
            c = !c;
    }
    return(c);
}


int g_pol_checkp(int n,double *x,double *y,double px,double py)
{
    register int i,c;
    double p1x,p1y,p2x,p2y,xt;

    c = 0;
    p1x = x[0];
    p1y = y[0];
    for (i = 1; i <= n; ++i) {
        p2x = x[i % n];
        p2y = y[i % n];
        if (py > dmin(p1y,p2y)) {
            if (py <= dmax(p1y,p2y)) {
                if (px <= dmax(p1x,p2x)) {
                    if (p1y != p2y) {
                        xt = (py - p1y) * (p2x - p1x) / (p2y - p1y) + p1x;
                        if (p1x == p2x || px <= xt)
                            c++;   
                    }   
                }
            }
        }
        p1x = p2x;
        p1y = p2y;
    }
    if (c % 2 == 0)
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  g_isectpol(n,x,y,xa,ya,xb,yb)                                           */
/*                                                                          */
/*  Return 1 if the line segment defined by [(xa,ya),(xb,yb)] has a point   */
/*  in common with the polygon defined by x[i],y[i] (i=0,n-1), otherwise    */
/*  return 0.                                                               */

int g_pol_isect(int n, double *x,double *y,double xa,double ya,
    double xb,double yb)
{
    register int i,r;
    double xxa,yya,xxb,yyb;

    for (i = 1; i < n; ++i) {
        r = g_intersect(xa,ya,xb,yb,x[i-1],y[i-1],x[i],y[i],&xxa,&yya,&xxb,&yyb);
        if (r != 0)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_chull(x,y,m,in,ih,il)                                                 */
/*                                                                          */
/*  Convex hull of a set of points. Adapted from Algorithm 523,             */
/*  Collected Alg. from ACM. See W.F. Eddy, A New Convex Hull Algorithm     */
/*  for Planar Sets, ACM Transactions on Mathematical Software 3 (1977),    */
/*  pp. 398 - 403.                                                          */
/*                                                                          */
/*  float *x    x coordinate of data points,x[i], i = 1,2,3,...             */
/*  float *y    y coordinate of data points                                 */
/*  int m       number of points in input subset                            */
/*  int *in     subscripts for array x[] of the points in input subset      */
/*  int *ih     subscripts of array x[] of vertices of convex hull          */
/*  int *il     linked list giving the elements of array ih[]               */
/*                                                                          */
/*  Return: -1 if insufficient memory, otherwise number of points in        */
/*                                                      ih[] and il[].      */
/*                                                                          */
/*  Note: all subscripts begin with 1.                                      */

int g_chull(float *x,float *y,int m,int *in,int *ih,int *il) 
{
    register int i,j;
    int kn,kx,mp1,min,mx,maxe,mine,nib,ma,inh,mm,mb,mbb,mxa,mxb,mxbb,ilinh;
    int nh,nia,*ia,*ib;

    if (m < 1)
        return(0);

    if (!(ia = (int *)calloc(m + 1,sizeof(int))))
        return(-1);       

    if (!(ib = (int *)calloc(m + 1,sizeof(int)))) {
        free((char *)ia);
        return(-1);       
    }
    if (m == 1)          
        goto CONV22;

    il[1] = 2;
    il[2] = 1;
    kn = in[1];
    kx = in[2];
    if (m == 2)          
        goto CONV21;
       
    mp1 = m + 1;
    min = 1;
    mx = 1;
    kx = in[1];
    maxe = mine = 0;

    /* FIND TWO VERTICES OF THE CONVEX HULL FOR THE INITIAL
       PARTITION */

    for (i = 2; i <= m; ++i) {
        j = in[i];
        if (x[j] > x[kx]) {
            maxe = 0;
            mx = i;
            kx = j;
        }
        else if (x[j] == x[kx]) 
            maxe = 1;

        if (x[j] < x[kn]) {
            mine = 0;
            min = i;
            kn = j;
        }
        else if (x[j] == x[kn])
            mine = 1;
    }

    /* IF THE MAX AND MIN ARE EQUAL, ALL M POINTS LIE ON A
       VERTICAL LINE */

    if (kx == kn)          
        goto CONV18;

    /* IF MAXE (OR MINE) HAS THE VALUE TRUE THERE ARE SEVERAL
       MAXIMA (OR MINIMA) WITH EQUAL FIRST COORDINATES */

    if (maxe || mine)           
        goto CONV23;

CONV7:
    ih[1] = kx;
    ih[2] = kn;
    nh = 3;
    inh = 1;
    nib = 1;
    ma = m;
    in[mx] = in[m];
    in[m] = kx;
    mm = m - 2;
    if (min == m)
        min = mx;
    in[min] = in[m - 1];
    in[m - 1] = kn;

    /* BEGIN BY PARTITIONING THE ROOT OF THE TREE */

    g_chull_split(x,y,mm,in,ih[1],ih[2],0,ia,&mb,&mxa,ib,&ia[ma],&mxbb);

    /* FIRST TRAVERSE THE LEFT HALF OF THE TREE
       START WITH THE LEFT SON */

    while (1) {
        nib += ia[ma];
        ma--;   
CONV9:  
        if (mxa == 0)                
            break;
        il[nh] = il[inh];
        il[inh] = nh;
        ih[nh] = ia[mxa]; 
        ia[mxa] = ia[mb];
        mb--;   
        nh += 1;
        if (mb == 0) {         
            inh = il[inh];
            break;  
        }
        ilinh = il[inh];

        g_chull_split(x,y,mb,ia,ih[inh],ih[ilinh],1,ia,&mbb,&mxa,ib + nib - 1,&ia[ma],&mxb);
        mb = mbb;
    }

    /* THEN THE RIGHT SON */

    while (1) {
        inh = il[inh];
        ma++;   
        nib -= ia[ma];
        if (ma >= m)           
            break;

        if (ia[ma] != 0) {           
            ilinh = il[inh];

            /* ON THE LEFT SIDE OF THE TREE, THE RIGHT SON OF A RIGHT SON
               MUST REPRESENT A SUBSET OF POINTS WHICH IS INSIDE A
               TRIANGLE WITH VERTICES WHICH ARE ALSO VERTICES OF THE
               CONVEX POLYGON AND HENCE THE SUBSET MAY BE NEGLECTED. */

            g_chull_split(x,y,ia[ma],ib + nib - 1,ih[inh],ih[ilinh],2,ia,&mb,&mxa,ib + nib - 1,&mbb,&mxb);
            ia[ma] = mbb;
            goto CONV9;
        }
    }             

    /* NOW TRAVERSE THE RIGHT HALF OF THE TREE */

    mxb = mxbb;
    ma = m;
    mb = ia[ma];
    nia = 1;
    ia[ma] = 0;

    /* START WITH THE RIGHT SON */

    while (1) {
        nia += ia[ma];
        ma--;   
CONV14: 
        if (mxb == 0)            
            break;

        il[nh] = il[inh];
        il[inh] = nh;
        ih[nh] = ib[mxb];
        ib[mxb] = ib[mb];
        mb--;   
        nh += 1;
        if (mb == 0) {           
            inh = il[inh];
            break;
        }

        ilinh = il[inh];

        g_chull_split(x,y,mb,ib + nib - 1,ih[inh],ih[ilinh],-1,ia + nia - 1,&ia[ma],&mxa,ib + nib - 1,&mbb,&mxb);
        mb = mbb;
    }

    /* THEN THE LEFT SON */

    while (1) {
        inh = il[inh];
        ma++;   
        nia -= ia[ma];
        if (ma == mp1)          
            goto CONVFin;

        if (ia[ma] != 0)            
            break;
    }
    ilinh = il[inh];

    /* ON THE RIGHT SIDE OF THE TREE, THE LEFT SON OF A LEFT SON
       MUST REPRESENT A SUBSET OF POINTS WHICH IS INSIDE A
       TRIANGLE WITH VERTICES WHICH ARE ALSO VERTICES OF THE
       CONVEX POLYGON AND HENCE THE SUBSET MAY BE NEGLECTED. */

    g_chull_split(x,y,ia[ma],ia + nia - 1,ih[inh],ih[ilinh],-2,ia + nia - 1,&mbb,&mxa,ib + nib - 1,&mb,&mxb);
    goto CONV14;

    /* ALL THE SPECIAL CASES ARE HANDLED DOWN HERE
       IF ALL THE POINTS LIE ON A VERTICAL LINE */

CONV18: 
    kx = in[1];
    kn = in[1];

    for (i = 1; i <= m; ++i) {
        j = in[i];
        if (y[j] > y[kx]) {        
            mx = i;
            kx = j;
        }
        if (y[j] < y[kn]) {         
            min = i;
            kn = j;
        }
    }

    if (kx == kn)          
        goto CONV22;

    /* IF THERE ARE ONLY TWO POINTS */

CONV21: 
    ih[1] = kx;
    ih[2] = kn;
    nh = 3;
    if (x[kn] == x[kx] && y[kn] == y[kx])
        nh = 2;
    goto CONVFin;

    /* IF THERE IS ONLY ONE POINT */

CONV22:
    nh = 2;
    ih[1] = in[1];
    il[1] = 1;
    goto CONVFin;

    /* MULTIPLE EXTREMES ARE HANDLED HERE
       IF THERE ARE SEVERAL POINTS WITH THE (SAME) LARGEST
       FIRST COORDINATE */

CONV23:
    if (maxe != 0) {         

        for (i = 1; i <= m; ++i) {
            j = in[i];
            if (x[j] == x[kx] && y[j] > y[kx]) {         
                mx = i;
                kx = j;
            }
        }
    }

    /* IF THERE ARE SEVERAL POINTS WITH THE (SAME) SMALLEST
       FIRST COORDINATE */
    
    if (mine != 0) {             
        for (i = 1; i <= m; ++i) {
            j = in[i];
            if (x[j] == x[kn] && y[j] < y[kn]) {       
                min = i;
                kn = j;
            }
        }
    }
    goto CONV7;

CONVFin:  
    free((char *)ia);
    free((char *)ib);
    nh -= 1;

    return(nh);
}           

/*--------------------------------------------------------------------------*/
/*  g_chull_split(x,y,m,in,ii,jj,s,iabv,na,maxa,ibel,nb,maxb)               */
/*                                                                          */
/*  Function used by convex().                                              */
/*                                                                          */
/*  THIS SUBROUTINE TAKES THE M POINTS OF ARRAY X,Y WHOSE                   */
/*  SUBSCRIPTS ARE IN ARRAY IN AND PARTITIONS THEM BY THE                   */
/*  LINE JOINING THE TWO POINTS IN ARRAY X,Y WHOSE SUBSCRIPTS               */
/*  ARE II AND JJ. THE SUBSCRIPTS OF THE POINTS ABOVE THE                   */
/*  LINE ARE PUT INTO ARRAY IABV, AND THE SUBSCRIPTS OF THE                 */
/*  POINTS BELOW ARE PUT INTO ARRAY IBEL. NA AND NB ARE,                    */
/*  RESPECTIVELY, THE NUMBER OF POINTS ABOVE THE LINE AND THE               */
/*  NUMBER BELOW. MAXA AND MAXB ARE THE SUBSCRIPTS FOR ARRAY                */
/*  X,Y OF THE POINT FURTHEST ABOVE THE LINE AND THE POINT                  */
/*  FURTHEST BELOW, RESPECTIVELY. IF EITHER SUBSET IS NULL                  */
/*  THE CORRESPONDING SUBSCRIPT (MAXA OR MAXB) IS SET TO ZERO               */
  
void g_chull_split(float *x,float *y,int m,int *in,int ii,int jj,int s,
    int *iabv,int *na,int *maxa,int *ibel,int *nb,int *maxb) 
{
    register int i,is;
    int t;
    float z,dir,a,b,xt,up,down;

    t = 0;
    *na = 0;
    *nb = 0;
    *maxa = 0;
    *maxb = 0;
    up = down = 0.0;

    /* CHECK TO SEE IF THE LINE IS VERTICAL */

    if (x[jj] != x[ii]) {           
        a = (y[jj] - y[ii]) / (x[jj] - x[ii]);
        b = y[ii] - a * x[ii];
    }
    else {
        xt = x[ii];
        if (y[jj] < y[ii])
            dir = -1.0;
        else
            dir = 1.0;

        if (s < 0)
            dir *= -1.0;

        t = 1;   
    }

    for (i = 1; i <= m; ++i) {
        is = in[i];
        if (t)           
            z = dir * (x[is] - xt);
        else  
            z = y[is] - a * x[is] - b;

        if (z <= 0.0) {

            if (s != 2 && z < 0.0) { /* THE POINT IS BELOW THE LINE */
                *nb += 1;
                ibel[*nb] = is;
                if (z <= down) {       
                    down = z;
                    *maxb = *nb;
                }
            }
        }
        else {  /* THE POINT IS ABOVE THE LINE */
            if (s != -2) {       
                *na += 1;
                iabv[*na] = is;
                if (z >= up) {       
                    up = z;
                    *maxa = *na;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  g_xchxy(x,y)    Exchange x and y                                        */

void g_xchxy(double *x,double *y)
{
    double tmp;
    tmp = *x;
    *x = *y;
    *y = tmp;
}

/* ------------------------------------------------------------------------ */
/*  g_len(xa,ya,xb,yb)  Return euclidean length of the line segment         */
/*                      [(xa,ya),(xb,yb)].                                  */

double g_len(double xa,double ya,double xb,double yb)
{
    xa -= xb;
    ya -= yb;
    return(sqrt(xa * xa + ya * ya));
}

/* ------------------------------------------------------------------------ */
/*  g_vlen(x,y)     Return euclidean length of the vector (x,y).            */

double g_vlen(double x,double y)
{
    return(sqrt(x * x + y * y));
}

/* ------------------------------------------------------------------------ */
/*  g_poly_ccw(n,x,y)                                                       */
/*                                                                          */
/*  Check whether the polygon defined by x[i], y[i] has a positive (ccw)    */
/*  orientation. If not, reverse order.                                     */

void g_poly_ccw(int n,double *x,double *y)
{
    register int j,k;
    double t;

    if (g_pol_area(n,x,y) >= 0.0)
        return;

    j = 0;
    k = n - 1;
    while (j < k) {
        t = x[j];
        x[j] = x[k];
        x[k] = t;
        t = y[j];
        y[j] = y[k];
        y[k] = t;
        j++; 
        k--;
    }
}

/* ------------------------------------------------------------------------ */
/*  g_line_intersect(p1x,p1y,p2x,p2y,q1x,q1y,q2x,q2y,s1x,s1y,s2x,s2y,eps)   */
/*                                                                          */
/*  Given two line segments [(p1x,p1y),(p2x,p2y)] and [(q1x,q1y),(q2x,q2y)] */
/*  this functions checks their intersection.                               */
/*                                                                          */
/*  Return 0    if no intersection                                          */
/*         1    single intersection point returned in (s1x,s1y).            */
/*         2    overlap in a line segment, in this case the overlap is      */
/*              returned in [(s1x,s1y),(s2x,s2y)]                           */
/*                                                                          */
/*  The code is adapted from the MultClip program written and copyrighted   */
/*  (according to GPL) by Michael Leonov and Alexey Nikitin.                */
/*                                                                          */
/*  Note: two points are treated as equal if their distance is less         */
/*  than eps.                                                               */

typedef double Vector[2];

#define Vsub2(r,a,b)    {(r)[0] = (a)[0] - (b)[0]; (r)[1] = (a)[1] - (b)[1];}
#define Vadd2(r,a,b)    {(r)[0] = (a)[0] + (b)[0]; (r)[1] = (a)[1] + (b)[1];}

#define Vcpy2(r,a)      {(r)[0] = (a)[0]; (r)[1] = (a)[1];}
#define Vequ2(a,b)      ((a)[0] == (b)[0] && (a)[1] == (b)[1])

#define Vadds(r,a,b,s)  {(r)[0] = ((a)[0] + (b)[0]) * (s); (r)[1] = ((a)[1] + (b)[1]) * (s);}

#define Vswp2(a,b) { float t; \
    t = (a)[0], (a)[0] = (b)[0], (b)[0] = t; \
    t = (a)[1], (a)[1] = (b)[1], (b)[1] = t; \
}

double vect_det2(Vector v1,Vector v2)
{
    return ((v1[0] * v2[1]) - (v2[0] * v1[1]));
}

double vect_m_dist(Vector v1,Vector v2)
{
    double dx = v1[0] - v2[0];
    double dy = v1[1] - v2[1];
    double dd = sqrt(dx * dx + dy * dy);
    
    if (dx > 0)  
        return(+dd);
    if (dx < 0)
        return(-dd);
    if (dy > 0)
        return(+dd);
    return(-dd);
}
 
void RecursiveIsect(Vector p1,Vector p2,Vector q1,Vector q2,Vector ret);

int g_line_intersect(double p1x,double p1y,double p2x,double p2y,
    double q1x,double q1y,double q2x,double q2y,
    double *s1x,double *s1y,double *s2x,double *s2y,double eps)
{                  
    double rpx,rpy,rqx,rqy,t,deel;
    Vector p1,p2,q1,q2,s1,s2,p,pq1,pq2,q,qp1,qp2,s,qs;
    int c1,c2,c1s,c2s;

    p1[0] = p1x; p1[1] = p1y;
    p2[0] = p2x; p2[1] = p2y;
    q1[0] = q1x; q1[1] = q1y;
    q2[0] = q2x; q2[1] = q2y;

    if (dmax(p1[0],p2[0]) < dmin(q1[0],q2[0]) ||
        dmax(q1[0],q2[0]) < dmin(p1[0],p2[0]))
        return(0);

    rpx = p2[0] - p1[0];
    rpy = p2[1] - p1[1];
    rqx = q2[0] - q1[0];
    rqy = q2[1] - q1[1];

    deel = rpx * rqy - rpy * rqx; /* vect_det(rp,rq); */

    /* every value below eps is considered as being 0. */

    if (fabs(deel) < eps) {           /* parallel */
         
        double inpr, dc1, dc2, d1, d2, h;   /* Check too see whether p1-p2 and q1-q2 are on the same line */
        Vector hp1, hq1, hp2, hq2, q1p1, q1q2;

        Vsub2(q1p1,q1,p1);
        Vsub2(q1q2,q1,q2);

        inpr = vect_det2(q1p1, q1q2);

        /* If this product is not zero then p1 is not on q1-q2! */

        if (fabs(inpr) > eps)
            return(0);

        dc1 = 0.0;  /* m_len(p1 - p1) */

        dc2 = vect_m_dist(p1,p2);
        d1 = vect_m_dist(p1,q1);
        d2 = vect_m_dist(p1,q2);

        /* Sorting the independent points from small to large: */

        Vcpy2(hp1,p1);
        Vcpy2(hp2,p2);
        Vcpy2(hq1,q1);
        Vcpy2(hq2,q2);

        if (dc1 > dc2) {        /* hv and h are used as help-variable. */
            Vswp2(hp1,hp2);
            h = dc1, dc1 = dc2, dc2 = h;
        }
        if (d1 > d2) {
            Vswp2(hq1, hq2);
            h = d1, d1 = d2, d2 = h;
        }

        /* Now the line-pieces are compared: */

        if (dc1 < d1) {
            if (dc2 < d1)
                return(0);
            if (dc2 < d2) {
                Vcpy2(s1,hp2);
                Vcpy2(s2,hq1);
            }
            else {
                Vcpy2(s1,hq1);
                Vcpy2(s2,hq2);
            };
        }
        else {
            if (dc1 > d2)
                return(0);

            if (dc2 < d2) {
                Vcpy2(s1,hp1);
                Vcpy2(s2,hp2);
            }
            else {
                Vcpy2(s1,hp1);
                Vcpy2(s2,hq2);
            };
        }
        *s1x = s1[0]; *s1y = s1[1];
        if (Vequ2(s1,s2))
            return(1);

        *s2x = s2[0]; *s2y = s2[1];
        return(2);
    }
    else {                             /* not parallel */
             
        /* We have the lines:
         * l1: p1 + s(p2 - p1)
         * l2: q1 + t(q2 - q1)
         * And we want to know the intersection point.
         * Calculate t:
         * p1 + s(p2-p1) = q1 + t(q2-q1)
         * which is similar to the two equations:
         * p1x + s * rpx = q1x + t * rqx
         * p1y + s * rpy = q1y + t * rqy
         * Multiplying these by rpy resp. rpx gives:
         * rpy * p1x + s * rpx * rpy = rpy * q1x + t * rpy * rqx
         * rpx * p1y + s * rpx * rpy = rpx * q1y + t * rpx * rqy
         * Subtracting these gives:
         * rpy * p1x - rpx * p1y = rpy * q1x - rpx * q1y + t * ( rpy * rqx - rpx * rqy )
         * So t can be isolated:
         * t = (rpy * ( p1x - q1x ) + rpx * ( - p1y + q1y )) / ( rpy * rqx - rpx * rqy )
         * and deel = rpx * rqy - rpy * rqx
         */

        double inp1, inp2, inp3, inp4, inp1s, inp2s, inp3s, inp4s;


        if (Vequ2(q1,p1) || Vequ2(q1,p2))
            Vcpy2(s1,q1)

        else if (Vequ2(q2,p1) || Vequ2(q2,p2))
            Vcpy2(s1,q2)

        else {
            t = -(rpy * (-q1[0] + p1[0]) + rpx * (q1[1] - p1[1])) / deel;
            s1[0] = (q1[0] + t * rqx);
            s1[1] = (q1[1] + t * rqy);
        }

        /*
         * The intersection point is valid if it is
         * 1) on q1-q2 --> t >= 0 && t <= 1
         * 2) on p1-p2 --> p1 must be on the other side of q1-q2 as p2
         *    This is so if the difference of the x coordinate of p1-s1 has the
         *    opposite sign as the x coordinate of p2-s2. So the multiplication of
         *    these two must be negative. This might fail if p1-p2 is a vertical line;
         *    this can be solved by adding the same product for the y coordinates
         */

        Vsub2(p,p2,p1);
        Vsub2(pq1,q1,p1);
        Vsub2(pq2,q2,p1);
        inp1 = vect_det2(p,pq1);
        inp2 = vect_det2(p,pq2);

        c1 = (inp1 * inp2 <= 0.0);

        Vsub2(q,q2,q1);
        Vsub2(qp1,p1,q1);
        Vsub2(qp2,p2,q1);

        inp3 = vect_det2(q, qp1);
        inp4 = vect_det2(q, qp2);
        c2 = (inp3 * inp4 <= 0.0);

        {
            double l_q, l_p;

            /* Say that *s1 equals one of the points if the relative
               distance is smaller than eps */

            l_q = g_len(q1[0],q1[1],q2[0],q2[1]);
            l_p = g_len(p1[0],p1[1],p2[0],p2[1]);

            if (eps > g_len(s1[0],s1[1],p2[0],p2[1]) / l_p) {
                *s1 = *p2;
                c2 = 2;
            }
            else if (eps > g_len(s1[0],s1[1],p1[0],p1[1]) / l_p) {
                *s1 = *p1;
                c2 = 2;
            }
            else {
                if (eps > g_len(s1[0],s1[1],q2[0],q2[1]) / l_q) {
                    *s1 = *q2;
                    c1 = 2;
                }
                else if (eps > g_len(s1[0],s1[1],q1[0],q1[1]) / l_q) {
                    *s1 = *q1;
                    c1 = 2;
                }
            }
        }
        Vsub2(s,s1,p1);
        inp1s = vect_det2(s,pq1);
        inp2s = vect_det2(s,pq2);
        c1s = Vequ2(s1, p1) ? -1 : inp1s * inp2s <= 0.0;

        Vsub2(qs,s1,q1);
        inp3s = vect_det2(qs,qp1);
        inp4s = vect_det2(qs,qp2);
        c2s = Vequ2(s1, q1) ? -1 : inp3s * inp4s <= 0.0;
        {
            int failed = 0; /* Roundig errors might make the statements below untrue */

            if (!((c1 == 0) || (c2 == 0) || (c1s == (c1 != 0)) || (c1s == -1)))
                failed = 1;
            else if (!((c1 == 0) || (c2 == 0) || (c2s == (c2 != 0)) || (c2s == -1)))
                failed = 2;
            else if (c1 && c2 && (g_len(s1[0],s1[1],q1[0],q1[1]) > g_len(q2[0],q2[1],q1[0],q1[1])))
                failed = 3;
            else if (c1 && c2 && (g_len(s1[0],s1[1],q2[0],q2[1]) > g_len(q2[0],q2[1],q1[0],q1[1])))
                failed = 4;
            else if (c1 && c2 && (g_len(s1[0],s1[1],p1[0],p1[1]) > g_len(p2[0],p2[1],p1[0],p1[1])))
                failed = 5;
            else if (c1 && c2 && (g_len(s1[0],s1[1],p2[0],p2[1]) > g_len(p2[0],p2[1],p1[0],p1[1])))
                failed = 6;
            if ((failed >= 3) && (c1 == 2 || c2 == 2)) {
                failed = -1;
                c1 = c2 = 0;
            }
            if (failed > 0)
                RecursiveIsect(p1,p2,q1,q2,s1);

            *s1x = s1[0]; *s1y = s1[1];
            *s2x = s2[0]; *s2y = s2[1];
            return(c1 && c2);
        }
    }
}
             
void RecursiveIsect(Vector p1,Vector p2,Vector q1,Vector q2,Vector ret)
{
    /* Find intersection point of p1-p2 and q1-q2 by iteratively taking
       the middle of p1-p2 */

    Vector q,s1,s2,m,tmp;
    double inp;

    Vsub2(q,q2,q1);

    if (g_vlen(q[0],q[1]) < g_len(p2[0],p2[1],p1[0],p1[1])) {
        RecursiveIsect(q1,q2,p1,p2,ret);
        return;
    }
    Vsub2(s1,p1,q1);
    Vsub2(s2,p2,q1);
    Vadds(m,s1,s2,0.5f);

    inp = vect_det2(q,s1);

    if (inp > 0.0) {
        Vcpy2(tmp,s1);
        Vcpy2(s1,s2);
        Vcpy2(s2,tmp);
    }

    while (!Vequ2(m,s1) && !Vequ2(m,s2)) {
/**
        assert(vect_det2(q, s1) <= 0);
        assert(vect_det2(q, s2) >= 0);
**/
        if (vect_det2(q,m) <= 0.0)
            Vcpy2(s1, m)
        else
            Vcpy2(s2, m)

        Vadds(m,s1,s2,0.5f);
    }
/*  assert(vect_len2(m) <= vect_dist2(q2, q1));     */

    Vadd2(ret,m,q1);
}

