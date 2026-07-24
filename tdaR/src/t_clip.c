/****************************************************************************/
/*  t_clip                                                                  */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_gm.h"
#include "tda_context.h"
#include "t_sd.h"

/*  functions in t_clip.c */

int sdclip(TDAContext *ctx);
void sdclip_err_msg(TDAContext *ctx, int n,int id);
void sdclip_write_header(TDAContext *ctx, FILE *fd,int id,int typ,int n,int sdid,int id1);
void sdclip_write_point(TDAContext *ctx, FILE *fd,double x,double y);
int sdclip_point_rec(TDAContext *ctx, double x,double y,int *id,int sdid,int *nrec);
int clip_inside_rec(TDAContext *ctx, double x,double y,double xmin,double ymin,double xmax, double ymax);
int sdclip_line_rec(TDAContext *ctx, int n,double *x,double *y,double *sx,double *sy, int *id,int sdid,int *nrec);
void sdclip_line_prn(TDAContext *ctx, int n,double *x,double *y,int id,int sdid,int id1);

int sdclip_poly_rec(TDAContext *ctx, int n,double *x,double *y,double *sx,double *sy,int *nf, int *d,double xmin,double ymin,double xmax,double ymax);


void sdclip_clip_prn(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax, int id,int sdid);
void sdclip_poly_prn(TDAContext *ctx, int n,double *x,double *y,int id,int sdid,int id1);
int sdclip_poly_list(TDAContext *ctx, int n,double *x,double *y,int *nf,int *d,int *id, int sdid,int *nr);

int clip_simp(TDAContext *ctx, int n,double *x,double *y,double xmin,double ymin, double xmax,double ymax);
int clip_simp_typ(TDAContext *ctx, double x,double y,double xmin,double ymin, double xmax,double ymax);
int clip_input(TDAContext *ctx, double *xmin,double *ymin,double *xmax,double *ymax);
int clip_line(TDAContext *ctx, double xa,double ya,double xb,double yb,double xmin,double ymin, double xmax,double ymax,double *sax,double *say,double *sbx,double *sby, int *ca,int *cb,int *ra,int *rb,int *dir);
int clip_line_code(TDAContext *ctx, double x,double y,double xmin,double ymin,double xmax, double ymax);


/* -##--------------------------------------------------------------------- */
/*  sdclip      Clipping of spatial objects.                                */
/*                                                                          */  
/*              sdclip(                                                     */
/*                  rec=...,    rectangle: xmin,ymin,xmax,ymax              */
/*                  fmt=...,    print format for coordinates, def. 10.4     */
/*              ) = output_file;                                            */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. The output file is again a spatial data file containing */
/*  all points, lines and polygons intersecting the specified rectangle.    */  
/*                                                                          */
/*  For polygon clipping, it is assumed that polygons are simple and also   */
/*  that their intersection with the clip region results in one or more     */ 
/*  simple polygons.                                                        */  
/*                                                                          */  
/*  Return: 0 if OK, -1 if error.                                           */

int sdclip(TDAContext *ctx)
{
    register int i;
    int err,r,rp,n,nn,nr,nrec,typ,id,sdid,cflag,nmax,nmax2,nmax3;
    int n1,n1p,n2,n2p,n3,n3p;

    err = -1;
    nrec = id = n1 = n2 = n3 = n1p = n2p = n3p = 0;

    printf1(ctx, "Clipping spatial objects. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto SDCLIPFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    nrec = 0;
    if (ctx->PMRECFlg) {

        printf1(ctx, "Region defined by rectangle: %g, %g, %g, %g\n",
            ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax);
        if (ctx->PMRECXMin >= ctx->PMRECXMax - ctx->EPSI1 || ctx->PMRECYMin >= ctx->PMRECYMax - ctx->EPSI1) {
            printf1(ctx, "This is not a valid rectangle.\n");
            goto SDCLIPFin;
        }
        nn = 4;
        if (alloc_acx(ctx, nn))
            goto SDCLIPFin;
        if (alloc_acy(ctx, nn))
            goto SDCLIPFin;

        ctx->AcX[0] = ctx->PMRECXMin;
        ctx->AcY[0] = ctx->PMRECYMin;
        ctx->AcX[1] = ctx->PMRECXMax;
        ctx->AcY[1] = ctx->PMRECYMin;
        ctx->AcX[2] = ctx->PMRECXMax;
        ctx->AcY[2] = ctx->PMRECYMax;
        ctx->AcX[3] = ctx->PMRECXMin;
        ctx->AcY[3] = ctx->PMRECYMax;
    }
    /************************
    else if (PMF2Def) {
        printf1(ctx, "option if is currently not supported.\n");
        goto SDCLIPFin;

        nn = clip_input(ctx, &PMRECXMin,&PMRECYMin,&PMRECXMax,&PMRECYMax);
        if (nn < 3)
            goto SDCLIPFin;
      
        printf1(ctx, "Number of points: %d\n",nn);
        printf1(ctx, "Bounding box: %g, %g, %g, %g\n",
            PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);
    }
    *******************/
    else {
        printf1(ctx, "Error: need rec parameter.\n");
        goto SDCLIPFin;
    }
/**
    tda_out("nn=%d\n",nn);
    for (i = 0; i < nn; ++i)
        tda_out("i=%3d x=%f y=%f\n",i,AcX[i],AcY[i]);
**/
    /* allocate temporary storage for line and polygon clipping */

    nmax2 = nmax3 = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i);                        
        if (typ == 2 || typ == 3) {
            n = (int)get_data(ctx, ctx->SDVarSDN,i);
            if (typ == 2)
                nmax2 = imax(ctx, nmax2,n);
            else
                nmax3 = imax(ctx, nmax3,n);
        }
    }
    nmax = imax(ctx, nmax2,nmax3);
    if (nmax > 0) {
        if (alloc_acu(ctx, 2 * nmax + 10))
            goto SDCLIPFin;          
        if (alloc_acv(ctx, 2 * nmax + 10))
            goto SDCLIPFin;          

        if (nmax3 > 0) {
            if (alloc_acn(ctx, 2 * nmax3 + 10))
                goto SDCLIPFin;          
            if (alloc_acr(ctx, 2 * nmax3 + 10))
                goto SDCLIPFin;          
        }
    }

    /* do for all spatial objects */

    for (i = 0; i < ctx->NOC; ++i) {

        if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        sdid = (int)get_data(ctx, ctx->SDVarSDID,i);

        if (typ == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(ctx, i,0,cflag,1)) < 1)
            goto SDCLIPFin;
             
        if (typ == 1) {             /* points */
            n1++;
            r = sdclip_point_rec(ctx, ctx->SDVarX[0],ctx->SDVarY[0],&id,sdid,&nr);
            if (r) {
                n1p++;
                nrec += nr;
            }
        }
        else if (typ == 2) {        /* lines */
            n2++;
            r = sdclip_line_rec(ctx, n,ctx->SDVarX,ctx->SDVarY,ctx->AcU,ctx->AcV,&id,sdid,&nr);
            if (r) {
                n2p++;
                nrec += nr;
            }
        }
        else if (typ == 3) {        /* polygons */
            n3++;
            rp = sdclip_poly_rec(ctx, n - 1,ctx->SDVarX,ctx->SDVarY,ctx->AcU,ctx->AcV,ctx->AcN,ctx->AcR,
                 ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax);

            if (rp < 0)
                goto SDCLIPFin;

            if (rp == 1) {          /* polygon completely inside clip region */

                id += 1;
                sdclip_poly_prn(ctx, n - 1,ctx->SDVarX,ctx->SDVarY,id,sdid,1);
                n3p += 1;
                nrec += n;
            }
            else if (rp == 2) {     /* clip region completely inside polygon */
                id += 1;
                sdclip_clip_prn(ctx, ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax,id,sdid);
                n3p += 1;
                nrec += 5;
            }
            else if (rp > 2) {      /* some non-trivial intersection */

                r = sdclip_poly_list(ctx, rp,ctx->AcU,ctx->AcV,ctx->AcN,ctx->AcR,&id,sdid,&nr);
                if (r < 0) {
                    sdclip_err_msg(ctx, 1,sdid);
                }
                else if (r > 0) {
                    n3p += r;
                    nrec += nr;
                }
            }

        }
    }
    printf1(ctx, "\nNumber of type 1 objects: %d. Written: %d\n",n1,n1p);
    printf1(ctx, "Number of type 2 objects: %d. Written: %d\n",n2,n2p);
    printf1(ctx, "Number of type 3 objects: %d. Written: %d\n",n3,n3p);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDCLIPFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_err_msg(n,id)    Error message.                                  */

void sdclip_err_msg(TDAContext *ctx, int n,int id)
{
    if (n == 1) {
        printf1(ctx, "Error in clipping polygon with ID %d. The polygon, or its intersection\n",id);
        printf1(ctx, "with the clip region, might not be simple. Will be skipped.\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  sdclip_write_header(fd,id,typ,n,sdid,id1)                               */
/*                                                                          */  
/*  Write object description record to output file fd.                      */

void sdclip_write_header(TDAContext *ctx, FILE *fd,int id,int typ,int n,int sdid,int id1)
{
    (void)ctx;
    fprintf(fd,"%8d %d %6d %6d %6d\n",id,typ,n,sdid,id1);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_write_point(fd,x,y)   Write (x,y) to output file.                */

void sdclip_write_point(TDAContext *ctx, FILE *fd,double x,double y)
{
    rt_fprintf_d(ctx, fd,ctx->PMFmtS,x);
    rt_fprintf_d(ctx, fd,ctx->PMFmtS,y);
    fprintf(fd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  sdclip_point_rec(x,y,id,sdid)                                           */
/*                                                                          */
/*  If point (x,y) fits into the clip region write to output file and       */
/*  return 1, otherwise return 0. Increase id, return number of records     */ 
/*  written to output file in nrec.                                         */
/*                                                                          */
/*  This function assumes a rectangular clip region defined by              */
/*  PMRECXMin, PMRECYMin, PMRECXMax, PMRECYMax.                             */

int sdclip_point_rec(TDAContext *ctx, double x,double y,int *id,int sdid,int *nrec)
{
    *nrec = 0;

    if (clip_inside_rec(ctx, x,y,ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax) == 0)      
        return(0);

    *id += 1;
    sdclip_write_header(ctx, ctx->PMFd,*id,1,1,sdid,1);
    /* fprintf(PMFd,"%8d 1 %6d %6d %6d\n",*id,1,sdid,1); */
    sdclip_write_point(ctx, ctx->PMFd,x,y);
    *nrec = 2;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  clip_inside_rec(x,y,xmin,ymin,xmax,ymax)                                */
/*                                                                          */
/*  Return 1 if (x,y) is inside the rectangle xmin,ymin,xmax,ymax,          */
/*  otherwise return 0.                                                     */

int clip_inside_rec(TDAContext *ctx, double x,double y,double xmin,double ymin,double xmax, double ymax)
{
    (void)ctx;        /* unused: the signature is shared */
    if (x < xmin || x > xmax || y < ymin || y > ymax)  
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_line_rec(n,x,y,sx,sy,id,sdid,nrec)                               */
/*                                                                          */
/*  Clip the line given by x[i], y[i] (i = 0,n-1) to the rectangle defined  */
/*  by PMRECXMin, PMRECYMin, PMRECXMax, PMRECYMax.                          */
/*                                                                          */
/*  Return number of lines written to the output file, number of records    */
/*  written in nrec.                                                        */

int sdclip_line_rec(TDAContext *ctx, int n,double *x,double *y,double *sx,double *sy, int *id,int sdid,int *nrec)
{
    register int i;
    int r,np,rt,id1,ca,cb,ra,rb,dir;
    double sax,say,sbx,sby;

    *nrec = 0;
    np = rt = id1 = 0;

    for (i = 1; i < n; ++i) {
        r = clip_line(ctx, x[i-1],y[i-1],x[i],y[i],ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,
            ctx->PMRECYMax,&sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir);

        if (r == 0 || g_len(ctx, sax,say,sbx,sby) < ctx->EPSI1) {
            if (np > 0) {
                *id += 1;
                id1 += 1;
                sdclip_line_prn(ctx, np,sx,sy,*id,sdid,id1);
                *nrec += np + 1;
                np = 0;
                rt = 1;
            }
        }
        else {
            if (ca != 0) {
                if (np > 0) {
                    *id += 1;
                    id1 += 1;
                    sdclip_line_prn(ctx, np,sx,sy,*id,sdid,id1);
                    *nrec += np + 1;
                    np = 0;
                    rt = 1;
                }
                sx[np] = sax;
                sy[np] = say;
                np++;
                sx[np] = sbx;
                sy[np] = sby;
                np++;
            }
            else {
                if (np == 0) {
                    sx[np] = sax;
                    sy[np] = say;
                    np++;
                }
                sx[np] = sbx;
                sy[np] = sby;
                np++;

                if (cb != 0) {
                    *id += 1;
                    id1 += 1;
                    sdclip_line_prn(ctx, np,sx,sy,*id,sdid,id1);
                    *nrec += np + 1;
                    np = 0;
                    rt = 1;
                }
            }          
        }
    }
    if (np > 0) {
        *id += 1;
        id1 += 1;
        sdclip_line_prn(ctx, np,sx,sy,*id,sdid,id1);
        *nrec += np + 1;
        np = 0;
        rt = 1;
    }
    return(rt);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_line_prn  Write part of the clipped line to output file.         */

void sdclip_line_prn(TDAContext *ctx, int n,double *x,double *y,int id,int sdid,int id1)
{
    register int i;

    if (n == 0)
        return;

    sdclip_write_header(ctx, ctx->PMFd,id,2,n,sdid,id1);
    /* fprintf(PMFd,"%8d 2 %6d %6d %6d\n",id,n,sdid,id1); */
    for (i = 0; i < n; ++i)  
        sdclip_write_point(ctx, ctx->PMFd,x[i],y[i]);
}

/* -##--------------------------------------------------------------------- */
/*  sdclip_poly_rec(n,x,y,sx,sy,nf,d,xmin,ymin,xmax,ymax)                   */
/*                                                                          */
/*  Clip the polygon given by x[i], y[i] (i = 0,n-1) to the rectangle       */
/*  defined by xmin, ymin, xmax, ymax.                                      */
/*                                                                          */
/*  Note: the function uses AcTmp, AcTmp1, AcK and AcI as temporary         */
/*  storage.                                                                */
/*                                                                          */
/*  Return:   -1 if insufficient memory                                     */
/*            -2 if error in algorithm                                      */
/*             0 if no valid intersection                                   */
/*             1 if polygon completely inside clip region                   */
/*             2 if clip region completely inside polygon                   */
/*           > 2 number of points in edge list (sx,sy,nf,d)                 */

int sdclip_poly_rec(TDAContext *ctx, int n,double *x,double *y,double *sx,double *sy,int *nf, int *d,double xmin,double ymin,double xmax,double ymax)
{
    register int i,k,ip;
    int ia,ib,r,np,nb,nnb,ni,ra,rb,k1,ip1,n1,n2,n3,n4,dir,bflag[6];
    int ca,cb,nfcnt,rflag;
    double sax,say,sbx,sby;

    g_poly_ccw(ctx, n,x,y);      /* make ccw orientation */

    /* preliminary check for inside/outside clip region */



    ni = n1 = n2 = n3 = n4 = 0;
    for (i = 0; i < n; ++i) {

/*
tda_out("ii=%3d x=%12.6f %12.6f\n",i,x[i],y[i]);
*/

        if (x[i] < xmin + ctx->EPSI1)
            n1++;
        else if (x[i] > xmax - ctx->EPSI1)
            n2++;

        if (y[i] < ymin + ctx->EPSI1)
            n3++;
        else if (y[i] > ymax - ctx->EPSI1)
            n4++;

        if (clip_inside_rec(ctx, x[i],y[i],xmin,ymin,xmax,ymax))      
            ni++;
    }
tda_out("ni=%d n=%d\n",ni,n);



    if (ni == n) {          /* completely inside */ 
        return(1);
    }
    if (n1 == n || n2 == n || n3 == n || n4 == n) {
        tda_out("OUTSIDE\n");
        return(0);
    }
    rflag = 1;
    nfcnt = nb = np = 0;
    bflag[1] = bflag[2] = bflag[3] = bflag[4] = bflag[5] = 0;

    for (i = 1; i <= n; ++i) {

        ia = i - 1;
        if (i == n)
            ib = 0;
        else
            ib = i;

        /***
        tda_out("\ni=%d xa=%g,%g xb=%g,%g\n", i, x[ia],y[ia],x[ib],y[ib]);
        ***/

        r = clip_line(ctx, x[ia],y[ia],x[ib],y[ib],xmin,ymin,xmax,ymax,
            &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir);

        /**
        tda_out("r=%d ra=%d rb=%d dir=%d  \n",r,ra,rb,dir);
        **/

        if (r == 0)  
            continue;
                       
        /**
        tda_out("sax=%g,%g sbx=%g,%g\n",sax,say,sbx,sby);
        **/

        if (ra == rb && dir == -1) {    
            tda_out("X-continue\n");
            continue;
        }

        if (dir > 0 && ra == rb) {
            tda_out("Z-continue\n");
            rflag = 1;
            continue;
        }

        /* flag for orientation */

        if (((ra == -1 || ra == 2) && (rb == -1 || rb == 1) && sax > sbx) ||
            ((ra == -3 || ra == 4) && (rb == -3 || rb == 3) && sax < sbx) ||
            ((ra == -4 || ra == 1) && (rb == -4 || rb == 4) && say < sby) ||
            ((ra == -2 || ra == 3) && (rb == -2 || rb == 2) && say > sby)) { 

            tda_out("RDIR-continue\n");
            continue;
        }
        /**
        tda_out("rflag = %d\n",rflag);
        **/


        if (rflag) {            /* keep first point */

            sx[np] = sax;
            sy[np] = say;
            nf[np] = -1;          
            d[np] = 0;
            if (ra != 0) {      /* on boundary */
                d[np] = 1;   
                nb++;
                if (ra > 0)
                    bflag[ra] = 1;
            }
            np++;
            rflag = 0;
        }
        if (i == n) {
            if (np > 0 && (d[0] == 0 || (fabs(sbx - x[0]) <= ctx->EPSI1 &&
                                         fabs(sby - y[0]) <= ctx->EPSI1))) {
                nf[np - 1] = 0;
                nfcnt++;
                break;
            }
        }
        if (np > 0) {
            nf[np - 1] = np;
            nfcnt++;
        }
        sx[np] = sbx;
        sy[np] = sby;
        nf[np] = -1;           
        d[np] = 0;
        if (rb != 0) {
            d[np] = 1;
            nb++;
            if (rb > 0)
                bflag[rb] = 1;
        }
        np++;

        if (dir == -1 || dir == 2)
            rflag = 1;

    }
/**        
    tda_out("\nVOR  np=%d ni=%d nb=%d nfcnt=%d     \n",np,ni,nb,nfcnt);
    for (i = 0; i < np; ++i) {
        tda_out("i=%4d x=%f y=%f nf=%3d d=%d\n",i,sx[i],sy[i],nf[i],d[i]  );
    }
**/          

    /* check whether clip region is inside the polygon */

    if (np > 0 && nfcnt == 0) {
        tda_out("Y-return\n");
        return(0);
    }
    if (np == 0) {  
        sax = (xmin + xmax) / 2.0;
        say = (ymin + ymax) / 2.0;

        if (g_inpoly(ctx, n,x,y,sax,say)) {      
            tda_out("CLIP COMPLETELY INSIDE \n");
            return(2);
        }
        tda_out("NOT INTERSECT\n");
        return(0);
    }
           
    tda_out("\nnp=%d nb=%d    bflag: ",np,nb );
    for (i = 1; i <= 4; ++i)
        tda_out("%d ",bflag[i]);
    newline(ctx);
           

    if (nb < 2) {           /* no valid intersection */
        tda_out("no valid intersection nb=%d\n",nb);
        return(0);
    }

    /* add corner points */

    if (bflag[1] == 0) {
        sx[np] = xmin;
        sy[np] = ymin;
        nf[np] = -1;
        d[np] = 2;
        nb++;
        np++;
    }
    if (bflag[2] == 0) {
        sx[np] = xmax;
        sy[np] = ymin;
        nf[np] = -1;
        d[np] = 2;
        nb++;
        np++;
    }
    if (bflag[3] == 0) {
        sx[np] = xmax;
        sy[np] = ymax;
        nf[np] = -1;
        d[np] = 2;
        nb++;
        np++;
    }
    if (bflag[4] == 0) {
        sx[np] = xmin;
        sy[np] = ymax;
        nf[np] = -1;
        d[np] = 2;
        nb++;
        np++;
    }

/**       
    tda_out("\nNach corner add   np=%d ni=%d nb=%d   \n",np,ni,nb );
    for (i = 0; i < np; ++i) {
        tda_out("i=%4d x=%f y=%f nf=%3d d=%d\n",i,sx[i],sy[i],nf[i],d[i]  );
    }
**/        

    /* sort points on boundary */

    if (alloc_actmp(ctx, nb + 1))
        return(-1);              
    if (alloc_actmp1(ctx, nb + 1))
        return(-1);              
    if (alloc_ack(ctx, nb + 1))
        return(-1);              
    if (alloc_aci(ctx, nb + 1))
        return(-1);              

    nnb = 0;
    for (i = 0; i < np; ++i) {
        if (d[i] != 0) {
            ctx->AcTmp[nnb] = sx[i];
            ctx->AcTmp1[nnb] = sy[i];
            ctx->AcI[nnb] = i;
            nnb++;
        }
    }

    if (sortdp2c(ctx, nnb,ctx->AcTmp,ctx->AcTmp1,ctx->AcK,xmin,ymin,xmax,ymax))
        return(-1);

/*         
    tda_out("\nafter sorting nbb=%d\n",nnb);
    for (i = 0; i < nnb; ++i) {
        k = AcK[i];
        ip = AcI[k];
        tda_out("aci=%4d x=%f y=%f nf=%3d d=%d\n",ip,AcTmp[k],AcTmp1[k],nf[ip],d[ip]);
    }
*/        

    /* complete the edge list */

    for (i = 0; i < nnb; ++i) {
        k = ctx->AcK[i];
        ip = ctx->AcI[k];
        if (nf[ip] < 0) {
            if (i == nnb - 1) 
                k1 = ctx->AcK[0];
            else
                k1 = ctx->AcK[i + 1];
          
            ip1 = ctx->AcI[k1];
            nf[ip] = ip1;
            if (nf[ip1] >= 0 && nf[ip1] == ip) { /* no valid intersection */
                tda_out("NOT VALID\n");
                return(0);
            }
        }
    }
/**       
    tda_out("\nnach completed    nbb=%d\n",nnb);
    for (i = 0; i < nnb; ++i) {
        k = AcK[i];
        ip = AcI[k];
        tda_out("aci=%4d x=%f y=%f nf=%3d d=%d\n",ip,AcTmp[k],AcTmp1[k],nf[ip],d[ip]);
    }
**/        

    return(np);
}


/* ------------------------------------------------------------------------ */
/*  sdclip_clip_prn(xmin,ymin,xmax,ymax,id,sdid)                            */
/*                                                                          */
/*  Write rectangular clip region to output file.                           */

void sdclip_clip_prn(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax, int id,int sdid)
{
    sdclip_write_header(ctx, ctx->PMFd,id,3,4,sdid,1);
    /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",id,4,sdid,1); */
    sdclip_write_point(ctx, ctx->PMFd,xmin,ymin);
    sdclip_write_point(ctx, ctx->PMFd,xmax,ymin);
    sdclip_write_point(ctx, ctx->PMFd,xmax,ymax);
    sdclip_write_point(ctx, ctx->PMFd,xmin,ymax);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_poly_prn(n,x,y,id,sdid,id1)                                      */
/*                                                                          */
/*  Write polygon x[i], y[i] (i=0,ln-1) to output file.                     */

void sdclip_poly_prn(TDAContext *ctx, int n,double *x,double *y,int id,int sdid,int id1)
{
    register int i;

    sdclip_write_header(ctx, ctx->PMFd,id,3,n,sdid,id1);
    /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",id,n,sdid,id1); */
    for (i = 0; i < n; ++i)  
        sdclip_write_point(ctx, ctx->PMFd,x[i],y[i]);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_poly_list(n,x,y,nf,d,id,sdid,nr)                                 */
/*                                                                          */
/*  This function gets an edge list created by sdclip_poly_rec and writes   */
/*  all polygons to the output file. Return in nr the number of records     */
/*  written to the output file.                                             */
/*                                                                          */
/*  Return number of polygons written.                                      */

int sdclip_poly_list(TDAContext *ctx, int n,double *x,double *y,int *nf,int *d,int *id, int sdid,int *nr)
{
    register int i,j,k;
    int m,np,id1;

    /*************
    tda_out("sdclip_poly_list\n");
    for (i = 0; i < n; ++i) {
        tda_out("i=%3d x=%f y=%f nf=%3d d=%d\n",i,x[i],y[i],nf[i],d[i]);
    }
    *********/

    *nr = np = id1 = 0;
    for (i = 0; i < n; ++i) {
        if (d[i] == 2)
            d[i] = -2;
        else
            d[i] = -1;
    }

    /* first check whether the edge list represents a set of valid
       simple polygons. */

    while (1) {
        j = -1;
        for (i = 0; i < n; ++i) {
            if (d[i] == -1) {
                j = i;
                break;
            }
        }
        if (j < 0)
            break;

        m = 1;
        k = j;
        d[k] = j;

        while (1) {
            k = nf[k];
            if (k == j)
                break;
            if (++m > n || d[k] == j) {
                m = 0;
                break;
            }
            d[k] = j;
        }
        if (m == 0)             /* no valid simple polygon */
            return(-1);
    }

    /* now write to output file */

    for (i = 0; i < n; ++i) {
        if (d[i] != -2)
            d[i] = -1;
    }

    while (1) {
        j = -1;
        for (i = 0; i < n; ++i) {
            if (d[i] == -1) {
                j = i;
                break;
            }
        }
        if (j < 0)
            break;

        m = 1;
        k = j;
        d[k] = j;

        while (1) {
            k = nf[k];
            if (k == j)
                break;
            if (++m > n || d[k] == j) {
                m = 0;
                break;
            }
            d[k] = j;
        }
        if (m == 0)         
            return(-1);

        if (m < 3)
            continue;

        *id += 1;
        id1 += 1;
        sdclip_write_header(ctx, ctx->PMFd,*id,3,m,sdid,id1);
        /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",*id,m,sdid,id1); */
        *nr += 1;
        k = j;
        while (1) {
            sdclip_write_point(ctx, ctx->PMFd,x[k],y[k]);
            *nr += 1;
            k = nf[k];
            if (k == j)
                break;
        }
        np++;
    }
    return(np);
}

/* ------------------------------------------------------------------------ */
/*  clip_simp(n,x,y,xmin,ymin,xmax,ymax)                                    */

int clip_simp(TDAContext *ctx, int n,double *x,double *y,double xmin,double ymin, double xmax,double ymax)
{
    int i,k,t0,t1,t2;

    t0 = clip_simp_typ(ctx, x[0],y[0],xmin,ymin,xmax,ymax);
    t1 = clip_simp_typ(ctx, x[1],y[1],xmin,ymin,xmax,ymax);

    k = 2;
    for (i = 2; i < n; ++i) {
        t2 = clip_simp_typ(ctx, x[i],y[i],xmin,ymin,xmax,ymax);
        if (t2 == 0 || t0 != t1 || t1 != t2) {
            x[k] = x[i];
            y[k] = y[i];
            t0 = t1;
            t1 = t2;
            k++;
        }
        else if (t0 == t1 && t1 == t2) {
            x[k - 1] = x[i];
            y[k - 1] = y[i];
            t1 = t2;
        }
    }
tda_out("INPUT n=%d output k=%d\n",n,k);

    return(k);
}

int clip_simp_typ(TDAContext *ctx, double x,double y,double xmin,double ymin, double xmax,double ymax)
{
    if (x < xmin - ctx->EPSI1) {
        if (y < ymin)
            return(8);
        else if (y > ymax)
            return(2);
        else
            return(1);
    }
    else if (x > xmax + ctx->EPSI1) {
        if (y < ymin)
            return(6);
        else if (y > ymax)
            return(4);
        else
            return(5);
    }
    else if (y < ymin - ctx->EPSI1)
        return(7);
    else if (y > ymax + ctx->EPSI1)
        return(3);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  clip_input      Read polygon from input file PMF2d.                     */
/*                  Return number of points, or -1 if error.                */

int clip_input(TDAContext *ctx, double *xmin,double *ymin,double *xmax,double *ymax) 
{
    int i,n,err;
    char buf[101],*p;
    double x,y;

    printf1(ctx, "Reading input file: %s\n",ctx->PMF2dName);

    err = 0;
    if (!fgets(buf,100,ctx->PMF2d)) { 
        printf1(ctx, "Error: cannot read object description record.\n");
        return(-1);
    }               
    if (err == 0) {
        p = skip_b(ctx, buf);
        if (sscanf(p,"%d",&n) != 1)
            err = -1;
    }
    if (err == 0) {
        p = skip_int(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%d",&n) != 1 || n != 3)
            err = -1;
    }
    if (err == 0) {
        p = skip_int(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%d",&n) != 1 || n < 3)
            err = -1;
    }
    if (err == -1) {
        printf1(ctx, "Error: invalid object description record.\n");
        return(-1);
    }
    if (alloc_acx(ctx, n))
        return(-1);             
    if (alloc_acy(ctx, n))
        return(-1);           

    for (i = 0; i < n; ++i) {
        if (!fgets(buf,100,ctx->PMF2d)) {
            printf1(ctx, "Error: cannot read coordinates record %d.\n",i + 1);
            return(-1);
        }
        p = skip_b(ctx, buf);
        if (sscanf(p,"%lf",&x) != 1) {
            printf1(ctx, "Error: invalid x coordinate in record %d\n",i + 1);
            return(-1);
        }
        p = skip_dbl(ctx, p);
        p = skip_b(ctx, p);
        if (sscanf(p,"%lf",&y) != 1) {
            printf1(ctx, "Error: invalid x coordinate in record %d\n",i + 1);
            return(-1);
        }
        ctx->AcX[i] = x;
        ctx->AcY[i] = y; 

        if (i == 0) {
            *xmin = *xmax = x;
            *ymin = *ymax = y; 
        }
        else {
            *xmin = dmin(ctx, *xmin,x);
            *xmax = dmax(ctx, *xmax,x);
            *ymin = dmin(ctx, *ymin,y);
            *ymax = dmax(ctx, *ymax,y);
        }
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  clip_line(xa,ya,xb,yb,xmin,ymin,xmax,ymax,sax,say,sbx,sby,ca,cb,        */
/*                                                            ra,rb,dir)    */
/*                                                                          */
/*  Given a line segment [(xa,ya),xb,yb)] find intersection with the        */
/*  rectangle defined by xmin, ymin, xmax,ymax. Uses a version of the       */
/*  Cohen-Sutherland algorithm.                                             */  
/*                                                                          */
/*  Return 0 if no intersection, 1 if intersection.                         */
/*  If intersection, return intersection points in (sax,say) and (sbx,sby)  */
/*  and codes in ca and cb.                                                 */
/*                                                                          */
/*  Return direction in dir, based on original ca and cb codes:             */ 
/*                                                                          */
/*  dir = 0   if ca = 0 and cb = 0                                          */
/*       -1   if ca = 0 and cb = 1                                          */
/*        1   if ca = 1 and cb = 0                                          */
/*        2   if ca = 1 and cb = 1                                          */
/*                                                                          */
/*  ra, rb =  0  inside rectangle (without boundary)                        */
/*            1  lower left corner                                          */  
/*            2  lower right corner                                         */
/*            3  upper right corner                                         */ 
/*            4  upper left corner                                          */
/*           -1  lower boundary                                             */
/*           -2  right boundary                                             */
/*           -3  upper boundary                                             */
/*           -4  left boundary                                              */
/*                                                                          */

int clip_line(TDAContext *ctx, double xa,double ya,double xb,double yb,double xmin,double ymin, double xmax,double ymax,double *sax,double *say,double *sbx,double *sby, int *ca,int *cb,int *ra,int *rb,int *dir)
{
    double m;

    *ca = clip_line_code(ctx, xa,ya,xmin,ymin,xmax,ymax);
    *cb = clip_line_code(ctx, xb,yb,xmin,ymin,xmax,ymax);

    if (*ca & *cb)
        return(0);

    if (*ca == 0) {
        if (*cb == 0)
            *dir = 0;
        else
            *dir = -1;
    }
    else if (*cb == 0)
        *dir = 1;
    else
        *dir = 2;

    *ra = *rb = 0;
    *sax = xa;
    *say = ya;
    *sbx = xb;
    *sby = yb;

    if (*ca == 0) {
        if (fabs(xa - xmin) < ctx->EPSI1) {
            if (fabs(ya - ymin) < ctx->EPSI1)
                *ra = 1;
            else if (fabs(ya - ymax) < ctx->EPSI1)
                *ra = 4;
            else
                *ra = -4;
        }
        else if (fabs(xa - xmax) < ctx->EPSI1) {
            if (fabs(ya - ymin) < ctx->EPSI1)
                *ra = 2;
            else if (fabs(ya - ymax) < ctx->EPSI1)
                *ra = 3;
            else
                *ra = -2;
        }
        else if (fabs(ya - ymin) < ctx->EPSI1)
            *ra = -1;
        else if (fabs(ya - ymax) < ctx->EPSI1)
            *ra = -3;
    }
    if (*cb == 0) {
        if (fabs(xb - xmin) < ctx->EPSI1) {
            if (fabs(yb - ymin) < ctx->EPSI1)
                *rb = 1;
            else if (fabs(yb - ymax) < ctx->EPSI1)
                *rb = 4;
            else
                *rb = -4;
        }
        else if (fabs(xb - xmax) < ctx->EPSI1) {
            if (fabs(yb - ymin) < ctx->EPSI1)
                *rb = 2;
            else if (fabs(yb - ymax) < ctx->EPSI1)
                *rb = 3;
            else
                *rb = -2;
        }
        else if (fabs(yb - ymin) < ctx->EPSI1)
            *rb = -1;
        else if (fabs(yb - ymax) < ctx->EPSI1)
            *rb = -3;
    }
    if (*ca == 0 && *cb == 0)
        return(1);
      
    if (fabs(xa - xb) < ctx->EPSI1) {
        if (xa < xmin || xb > xmax)
            return(0);

        if (*ca != 0) {
            if (ya < ymin) {
                *say = ymin;
                                                        
                if (fabs(xa - xmin) < ctx->EPSI1)
                    *ra = 1;
                else if (fabs(xa - xmax) < ctx->EPSI1)
                    *ra = 2;
                else
                    *ra = -1;
            }
            else if (ya > ymax) {
                *say = ymax;

                if (fabs(xa - xmin) < ctx->EPSI1)
                    *ra = 4;
                else if (fabs(xa - xmax) < ctx->EPSI1)
                    *ra = 3;
                else
                    *ra = -3;
            }
        }
        if (*cb != 0) {
            if (yb < ymin) {
                *sby = ymin;

                if (fabs(xb - xmin) < ctx->EPSI1)
                    *rb = 1;
                else if (fabs(xb - xmax) < ctx->EPSI1)
                    *rb = 2;
                else
                    *rb = -1;
            }
            else if (yb > ymax) {
                *sby = ymax;

                if (fabs(xb - xmin) < ctx->EPSI1)
                    *rb = 4;
                else if (fabs(xa - xmax) < ctx->EPSI1)
                    *rb = 3;
                else
                    *rb = -3;
            }
        }
        if (*ra == *rb && *ra > 0 && *dir == 2)   
            return(0);
        return(1);
    }

    m = (ya - yb) / (xa - xb);
         
    if (*ca) {
        if (xa < xmin) {
            *sax = xmin;
            *say = yb + m * (xmin - xb);

            if (*say < ymin) {
                *say = ymin;
                *sax = xb + (ymin - yb) / m;
            }
            else if (*say > ymax) {
                *say = ymax;
                *sax = xb + (ymax - yb) / m;
            }
        }
        else if (xa > xmax) {
            *sax = xmax;
            *say = yb + m * (xmax - xb);

            if (*say < ymin) {
                *say = ymin;
                *sax = xb + (ymin - yb) / m;
            }
            else if (*say > ymax) {
                *say = ymax;
                *sax = xb + (ymax - yb) / m;
            }
        }
        else if (ya < ymin) {
            *say = ymin;
            *sax = xb + (ymin - yb) / m;
        }
        else if (ya > ymax) {
            *say = ymax;
            *sax = xb + (ymax - yb) / m;
        }
        if (*sax < xmin || *sax > xmax || *say < ymin || *say > ymax)
            return(0);

        if (fabs(*sax - xmin) < ctx->EPSI1) {
            if (fabs(*say - ymin) < ctx->EPSI1)
                *ra = 1;
            else if (fabs(*say - ymax) < ctx->EPSI1)
                *ra = 4;
            else
                *ra = -4;
        }
        else if (fabs(*sax - xmax) < ctx->EPSI1) {
            if (fabs(*say - ymin) < ctx->EPSI1)
                *ra = 2;
            else if (fabs(*say - ymax) < ctx->EPSI1)
                *ra = 3;
            else
                *ra = -2;
        }
        else if (fabs(*say - ymin) < ctx->EPSI1)
            *ra = -1;
        else if (fabs(*say - ymax) < ctx->EPSI1)
            *ra = -3;
    }
    if (*cb) {
        if (xb < xmin) {
            *sbx = xmin;
            *sby = yb + m * (xmin - xb);

            if (*sby < ymin) {
                *sby = ymin;
                *sbx = xb + (ymin - yb) / m;
            }
            else if (*sby > ymax) {
                *sby = ymax;
                *sbx = xb + (ymax - yb) / m;
            }
        }
        else if (xb > xmax) {
            *sbx = xmax;
            *sby = yb + m * (xmax - xb);

            if (*sby < ymin) {
                *sby = ymin;
                *sbx = xb + (ymin - yb) / m;
            }
            else if (*sby > ymax) {
                *sby = ymax;
                *sbx = xb + (ymax - yb) / m;
            }
        }
        else if (yb < ymin) {
            *sby = ymin;
            *sbx = xb + (ymin - yb) / m;
        }
        else if (yb > ymax) {
            *sby = ymax;
            *sbx = xb + (ymax - yb) / m;
        }
        if (*sbx < xmin || *sbx > xmax || *sby < ymin || *sby > ymax)
            return(0);

        if (fabs(*sbx - xmin) < ctx->EPSI1) {
            if (fabs(*sby - ymin) < ctx->EPSI1)
                *rb = 1;
            else if (fabs(*sby - ymax) < ctx->EPSI1)
                *rb = 4;
            else
                *rb = -4;
        }
        else if (fabs(*sbx - xmax) < ctx->EPSI1) {
            if (fabs(*sby - ymin) < ctx->EPSI1)
                *rb = 2;
            else if (fabs(*sby - ymax) < ctx->EPSI1)
                *rb = 3;
            else
                *rb = -2;
        }
        else if (fabs(*sby - ymin) < ctx->EPSI1)
            *rb = -1;
        else if (fabs(*sby - ymax) < ctx->EPSI1)
            *rb = -3;
    }
    if (*ra == *rb && *ra > 0 && *dir == 2)   
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  clip_line_code. Return code for clip_line function.                     */

int clip_line_code(TDAContext *ctx, double x,double y,double xmin,double ymin,double xmax, double ymax)
{
    (void)ctx;        /* unused: the signature is shared */
    int c = 0;
    if (x < xmin)
        c |= 0x08;
    else if (x > xmax)
        c |= 0x04;
    if (y < ymin)
        c |= 0x02;
    else if (y > ymax)
        c |= 0x01;
    return(c);
}



