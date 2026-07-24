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

/*  functions in t_clip.c */

int sdclip(void);
void sdclip_err_msg(int n,int id);
void sdclip_write_header(FILE *fd,int id,int typ,int n,int sdid,int id1);
void sdclip_write_point(FILE *fd,double x,double y);
int sdclip_point_rec(double x,double y,int *id,int sdid,int *nrec);
int clip_inside_rec(double x,double y,double xmin,double ymin,double xmax,
    double ymax);
int sdclip_line_rec(int n,double *x,double *y,double *sx,double *sy,
    int *id,int sdid,int *nrec);
void sdclip_line_prn(int n,double *x,double *y,int id,int sdid,int id1);

int sdclip_poly_rec(int n,double *x,double *y,double *sx,double *sy,int *nf,
    int *d,double xmin,double ymin,double xmax,double ymax);


void sdclip_clip_prn(double xmin,double ymin,double xmax,double ymax,
    int id,int sdid);
void sdclip_poly_prn(int n,double *x,double *y,int id,int sdid,int id1);
int sdclip_poly_list(int n,double *x,double *y,int *nf,int *d,int *id,
    int sdid,int *nr);

int clip_simp(int n,double *x,double *y,double xmin,double ymin,
    double xmax,double ymax);
int clip_simp_typ(double x,double y,double xmin,double ymin,
    double xmax,double ymax);
int clip_input(double *xmin,double *ymin,double *xmax,double *ymax);
int clip_line(double xa,double ya,double xb,double yb,double xmin,double ymin,
    double xmax,double ymax,double *sax,double *say,double *sbx,double *sby, 
    int *ca,int *cb,int *ra,int *rb,int *dir);
int clip_line_code(double x,double y,double xmin,double ymin,double xmax,
    double ymax);


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

int sdclip(void)
{
    register int i,j;
    int err,r,rp,n,nn,nr,ns,nrec,typ,id,sdid,cflag,nmax,nmax2,nmax3;
    int n1,n1p,n2,n2p,n3,n3p;

    err = -1;
    nrec = id = n1 = n2 = n3 = n1p = n2p = n3p = 0;

    printf1("Clipping spatial objects. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto SDCLIPFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    nrec = 0;
    if (PMRECFlg) {

        printf1("Region defined by rectangle: %g, %g, %g, %g\n",
            PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);
        if (PMRECXMin >= PMRECXMax - EPSI1 || PMRECYMin >= PMRECYMax - EPSI1) {
            printf1("This is not a valid rectangle.\n");
            goto SDCLIPFin;
        }
        nn = 4;
        if (alloc_acx(nn))
            goto SDCLIPFin;
        if (alloc_acy(nn))
            goto SDCLIPFin;

        AcX[0] = PMRECXMin;
        AcY[0] = PMRECYMin;
        AcX[1] = PMRECXMax;
        AcY[1] = PMRECYMin;
        AcX[2] = PMRECXMax;
        AcY[2] = PMRECYMax;
        AcX[3] = PMRECXMin;
        AcY[3] = PMRECYMax;
    }
    /************************
    else if (PMF2Def) {
        printf1("option if is currently not supported.\n");
        goto SDCLIPFin;

        nn = clip_input(&PMRECXMin,&PMRECYMin,&PMRECXMax,&PMRECYMax);
        if (nn < 3)
            goto SDCLIPFin;
      
        printf1("Number of points: %d\n",nn);
        printf1("Bounding box: %g, %g, %g, %g\n",
            PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);
    }
    *******************/
    else {
        printf1("Error: need rec parameter.\n");
        goto SDCLIPFin;
    }
/**
    printf("nn=%d\n",nn);
    for (i = 0; i < nn; ++i)
        printf("i=%3d x=%f y=%f\n",i,AcX[i],AcY[i]);
**/
    /* allocate temporary storage for line and polygon clipping */

    nmax2 = nmax3 = 0;
    for (i = 0; i < NOC; ++i) {
        typ  = (int)get_data(SDVarSDTyp,i);                        
        if (typ == 2 || typ == 3) {
            n = (int)get_data(SDVarSDN,i);
            if (typ == 2)
                nmax2 = imax(nmax2,n);
            else
                nmax3 = imax(nmax3,n);
        }
    }
    nmax = imax(nmax2,nmax3);
    if (nmax > 0) {
        if (alloc_acu(2 * nmax + 10))
            goto SDCLIPFin;          
        if (alloc_acv(2 * nmax + 10))
            goto SDCLIPFin;          

        if (nmax3 > 0) {
            if (alloc_acn(2 * nmax3 + 10))
                goto SDCLIPFin;          
            if (alloc_acr(2 * nmax3 + 10))
                goto SDCLIPFin;          
        }
    }

    /* do for all spatial objects */

    for (i = 0; i < NOC; ++i) {

        if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        sdid = (int)get_data(SDVarSDID,i);

        if (typ == 3)
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(i,0,cflag,1)) < 1)
            goto SDCLIPFin;
             
        if (typ == 1) {             /* points */
            n1++;
            r = sdclip_point_rec(SDVarX[0],SDVarY[0],&id,sdid,&nr);
            if (r) {
                n1p++;
                nrec += nr;
            }
        }
        else if (typ == 2) {        /* lines */
            n2++;
            r = sdclip_line_rec(n,SDVarX,SDVarY,AcU,AcV,&id,sdid,&nr);
            if (r) {
                n2p++;
                nrec += nr;
            }
        }
        else if (typ == 3) {        /* polygons */
            n3++;
            rp = sdclip_poly_rec(n - 1,SDVarX,SDVarY,AcU,AcV,AcN,AcR,
                 PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);

            if (rp < 0)
                goto SDCLIPFin;

            if (rp == 1) {          /* polygon completely inside clip region */

                id += 1;
                sdclip_poly_prn(n - 1,SDVarX,SDVarY,id,sdid,1);
                n3p += 1;
                nrec += n;
            }
            else if (rp == 2) {     /* clip region completely inside polygon */
                id += 1;
                sdclip_clip_prn(PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax,id,sdid);
                n3p += 1;
                nrec += 5;
            }
            else if (rp > 2) {      /* some non-trivial intersection */

                r = sdclip_poly_list(rp,AcU,AcV,AcN,AcR,&id,sdid,&nr);
                if (r < 0) {
                    sdclip_err_msg(1,sdid);
                }
                else if (r > 0) {
                    n3p += r;
                    nrec += nr;
                }
            }

        }
    }
    printf1("\nNumber of type 1 objects: %d. Written: %d\n",n1,n1p);
    printf1("Number of type 2 objects: %d. Written: %d\n",n2,n2p);
    printf1("Number of type 3 objects: %d. Written: %d\n",n3,n3p);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDCLIPFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_err_msg(n,id)    Error message.                                  */

void sdclip_err_msg(int n,int id)
{
    if (n == 1) {
        printf1("Error in clipping polygon with ID %d. The polygon, or its intersection\n",id);
        printf1("with the clip region, might not be simple. Will be skipped.\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  sdclip_write_header(fd,id,typ,n,sdid,id1)                               */
/*                                                                          */  
/*  Write object description record to output file fd.                      */

void sdclip_write_header(FILE *fd,int id,int typ,int n,int sdid,int id1)
{
    fprintf(fd,"%8d %d %6d %6d %6d\n",id,typ,n,sdid,id1);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_write_point(fd,x,y)   Write (x,y) to output file.                */

void sdclip_write_point(FILE *fd,double x,double y)
{
    fprintf(fd,PMFmtS,x);
    fprintf(fd,PMFmtS,y);
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

int sdclip_point_rec(double x,double y,int *id,int sdid,int *nrec)
{
    *nrec = 0;

    if (clip_inside_rec(x,y,PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax) == 0)      
        return(0);

    *id += 1;
    sdclip_write_header(PMFd,*id,1,1,sdid,1);
    /* fprintf(PMFd,"%8d 1 %6d %6d %6d\n",*id,1,sdid,1); */
    sdclip_write_point(PMFd,x,y);
    *nrec = 2;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  clip_inside_rec(x,y,xmin,ymin,xmax,ymax)                                */
/*                                                                          */
/*  Return 1 if (x,y) is inside the rectangle xmin,ymin,xmax,ymax,          */
/*  otherwise return 0.                                                     */

int clip_inside_rec(double x,double y,double xmin,double ymin,double xmax,
    double ymax)
{
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

int sdclip_line_rec(int n,double *x,double *y,double *sx,double *sy,
    int *id,int sdid,int *nrec)
{
    register int i;
    int r,np,rt,id1,first,ca,cb,ra,rb,dir;
    double sax,say,sbx,sby;

    *nrec = 0;
    np = rt = id1 = 0;

    for (i = 1; i < n; ++i) {
        r = clip_line(x[i-1],y[i-1],x[i],y[i],PMRECXMin,PMRECYMin,PMRECXMax,
            PMRECYMax,&sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir);

        if (r == 0 || g_len(sax,say,sbx,sby) < EPSI1) {
            if (np > 0) {
                *id += 1;
                id1 += 1;
                sdclip_line_prn(np,sx,sy,*id,sdid,id1);
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
                    sdclip_line_prn(np,sx,sy,*id,sdid,id1);
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
                    sdclip_line_prn(np,sx,sy,*id,sdid,id1);
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
        sdclip_line_prn(np,sx,sy,*id,sdid,id1);
        *nrec += np + 1;
        np = 0;
        rt = 1;
    }
    return(rt);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_line_prn  Write part of the clipped line to output file.         */

void sdclip_line_prn(int n,double *x,double *y,int id,int sdid,int id1)
{
    register int i;

    if (n == 0)
        return;

    sdclip_write_header(PMFd,id,2,n,sdid,id1);
    /* fprintf(PMFd,"%8d 2 %6d %6d %6d\n",id,n,sdid,id1); */
    for (i = 0; i < n; ++i)  
        sdclip_write_point(PMFd,x[i],y[i]);
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

int sdclip_poly_rec(int n,double *x,double *y,double *sx,double *sy,int *nf,
    int *d,double xmin,double ymin,double xmax,double ymax)
{
    register int i,k,ip;
    int ia,ib,r,np,na,nb,nnb,ni,ra,rb,k1,ip1,n1,n2,n3,n4,dir,bflag[6];
    int ca,cb,nfcnt,nflag,rflag,rdir;
    double sax,say,sbx,sby;

    g_poly_ccw(n,x,y);      /* make ccw orientation */

    /* preliminary check for inside/outside clip region */

printf("xmin=%g,%g xmax=%g,%g\n",xmin,ymin,xmax,ymax);


    ni = n1 = n2 = n3 = n4 = 0;
    for (i = 0; i < n; ++i) {

/*
printf("ii=%3d x=%12.6f %12.6f\n",i,x[i],y[i]);
*/

        if (x[i] < xmin + EPSI1)
            n1++;
        else if (x[i] > xmax - EPSI1)
            n2++;

        if (y[i] < ymin + EPSI1)
            n3++;
        else if (y[i] > ymax - EPSI1)
            n4++;

        if (clip_inside_rec(x[i],y[i],xmin,ymin,xmax,ymax))      
            ni++;
    }
printf("ni=%d n=%d\n",ni,n);



    if (ni == n) {          /* completely inside */ 
        printf("INSIDE\n");
        return(1);
    }
    if (n1 == n || n2 == n || n3 == n || n4 == n) {
        printf("OUTSIDE\n");
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
        printf("\ni=%d xa=%g,%g xb=%g,%g\n", i, x[ia],y[ia],x[ib],y[ib]);
        ***/

        r = clip_line(x[ia],y[ia],x[ib],y[ib],xmin,ymin,xmax,ymax,
            &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir);

        /**
        printf("r=%d ra=%d rb=%d dir=%d  \n",r,ra,rb,dir);
        **/

        if (r == 0)  
            continue;
                       
        /**
        printf("sax=%g,%g sbx=%g,%g\n",sax,say,sbx,sby);
        **/

        if (ra == rb && dir == -1) {    
            printf("X-continue\n");
            continue;
        }

        if (dir > 0 && ra == rb) {
            printf("Z-continue\n");
            rflag = 1;
            continue;
        }

        /* flag for orientation */

        rdir = 0;
        if (((ra == -1 || ra == 2) && (rb == -1 || rb == 1) && sax > sbx) ||
            ((ra == -3 || ra == 4) && (rb == -3 || rb == 3) && sax < sbx) ||
            ((ra == -4 || ra == 1) && (rb == -4 || rb == 4) && say < sby) ||
            ((ra == -2 || ra == 3) && (rb == -2 || rb == 2) && say > sby)) { 
            rdir = 1;

            printf("RDIR-continue\n");
            continue;
        }
        /**
        printf("rflag = %d\n",rflag);
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
            if (np > 0 && (d[0] == 0 || (fabs(sbx - x[0]) <= EPSI1 &&
                                         fabs(sby - y[0]) <= EPSI1))) {
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
    printf("\nVOR  np=%d ni=%d nb=%d nfcnt=%d     \n",np,ni,nb,nfcnt);
    for (i = 0; i < np; ++i) {
        printf("i=%4d x=%f y=%f nf=%3d d=%d\n",i,sx[i],sy[i],nf[i],d[i]  );
    }
**/          

    /* check whether clip region is inside the polygon */

    if (np > 0 && nfcnt == 0) {
        printf("Y-return\n");
        return(0);
    }
    if (np == 0) {  
        sax = (xmin + xmax) / 2.0;
        say = (ymin + ymax) / 2.0;

        if (g_inpoly(n,x,y,sax,say)) {      
            printf("CLIP COMPLETELY INSIDE \n");
            return(2);
        }
        printf("NOT INTERSECT\n");
        return(0);
    }
           
    printf("\nnp=%d nb=%d    bflag: ",np,nb );
    for (i = 1; i <= 4; ++i)
        printf("%d ",bflag[i]);
    newline();
           

    if (nb < 2) {           /* no valid intersection */
        printf("no valid intersection nb=%d\n",nb);
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
    printf("\nNach corner add   np=%d ni=%d nb=%d   \n",np,ni,nb );
    for (i = 0; i < np; ++i) {
        printf("i=%4d x=%f y=%f nf=%3d d=%d\n",i,sx[i],sy[i],nf[i],d[i]  );
    }
**/        

    /* sort points on boundary */

    if (alloc_actmp(nb + 1))
        return(-1);              
    if (alloc_actmp1(nb + 1))
        return(-1);              
    if (alloc_ack(nb + 1))
        return(-1);              
    if (alloc_aci(nb + 1))
        return(-1);              

    nnb = 0;
    for (i = 0; i < np; ++i) {
        if (d[i] != 0) {
            AcTmp[nnb] = sx[i];
            AcTmp1[nnb] = sy[i];
            AcI[nnb] = i;
            nnb++;
        }
    }

    if (sortdp2c(nnb,AcTmp,AcTmp1,AcK,xmin,ymin,xmax,ymax))
        return(-1);

/*         
    printf("\nafter sorting nbb=%d\n",nnb);
    for (i = 0; i < nnb; ++i) {
        k = AcK[i];
        ip = AcI[k];
        printf("aci=%4d x=%f y=%f nf=%3d d=%d\n",ip,AcTmp[k],AcTmp1[k],nf[ip],d[ip]);
    }
*/        

    /* complete the edge list */

    for (i = 0; i < nnb; ++i) {
        k = AcK[i];
        ip = AcI[k];
        if (nf[ip] < 0) {
            if (i == nnb - 1) 
                k1 = AcK[0];
            else
                k1 = AcK[i + 1];
          
            ip1 = AcI[k1];
            nf[ip] = ip1;
            if (nf[ip1] >= 0 && nf[ip1] == ip) { /* no valid intersection */
                printf("NOT VALID\n");
                return(0);
            }
        }
    }
/**       
    printf("\nnach completed    nbb=%d\n",nnb);
    for (i = 0; i < nnb; ++i) {
        k = AcK[i];
        ip = AcI[k];
        printf("aci=%4d x=%f y=%f nf=%3d d=%d\n",ip,AcTmp[k],AcTmp1[k],nf[ip],d[ip]);
    }
**/        

    return(np);
}


/* ------------------------------------------------------------------------ */
/*  sdclip_clip_prn(xmin,ymin,xmax,ymax,id,sdid)                            */
/*                                                                          */
/*  Write rectangular clip region to output file.                           */

void sdclip_clip_prn(double xmin,double ymin,double xmax,double ymax,
    int id,int sdid)
{
    sdclip_write_header(PMFd,id,3,4,sdid,1);
    /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",id,4,sdid,1); */
    sdclip_write_point(PMFd,xmin,ymin);
    sdclip_write_point(PMFd,xmax,ymin);
    sdclip_write_point(PMFd,xmax,ymax);
    sdclip_write_point(PMFd,xmin,ymax);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_poly_prn(n,x,y,id,sdid,id1)                                      */
/*                                                                          */
/*  Write polygon x[i], y[i] (i=0,ln-1) to output file.                     */

void sdclip_poly_prn(int n,double *x,double *y,int id,int sdid,int id1)
{
    register int i;

    sdclip_write_header(PMFd,id,3,n,sdid,id1);
    /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",id,n,sdid,id1); */
    for (i = 0; i < n; ++i)  
        sdclip_write_point(PMFd,x[i],y[i]);
}

/* ------------------------------------------------------------------------ */
/*  sdclip_poly_list(n,x,y,nf,d,id,sdid,nr)                                 */
/*                                                                          */
/*  This function gets an edge list created by sdclip_poly_rec and writes   */
/*  all polygons to the output file. Return in nr the number of records     */
/*  written to the output file.                                             */
/*                                                                          */
/*  Return number of polygons written.                                      */

int sdclip_poly_list(int n,double *x,double *y,int *nf,int *d,int *id,
    int sdid,int *nr)
{
    register int i,j,k;
    int m,np,id1;

    /*************
    printf("sdclip_poly_list\n");
    for (i = 0; i < n; ++i) {
        printf("i=%3d x=%f y=%f nf=%3d d=%d\n",i,x[i],y[i],nf[i],d[i]);
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
        sdclip_write_header(PMFd,*id,3,m,sdid,id1);
        /* fprintf(PMFd,"%8d 3 %6d %6d %6d\n",*id,m,sdid,id1); */
        *nr += 1;
        k = j;
        while (1) {
            sdclip_write_point(PMFd,x[k],y[k]);
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

int clip_simp(int n,double *x,double *y,double xmin,double ymin,
    double xmax,double ymax)
{
    int i,k,t0,t1,t2;

    t0 = clip_simp_typ(x[0],y[0],xmin,ymin,xmax,ymax);
    t1 = clip_simp_typ(x[1],y[1],xmin,ymin,xmax,ymax);

    k = 2;
    for (i = 2; i < n; ++i) {
        t2 = clip_simp_typ(x[i],y[i],xmin,ymin,xmax,ymax);
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
printf("INPUT n=%d output k=%d\n",n,k);

    return(k);
}

int clip_simp_typ(double x,double y,double xmin,double ymin,
    double xmax,double ymax)
{
    if (x < xmin - EPSI1) {
        if (y < ymin)
            return(8);
        else if (y > ymax)
            return(2);
        else
            return(1);
    }
    else if (x > xmax + EPSI1) {
        if (y < ymin)
            return(6);
        else if (y > ymax)
            return(4);
        else
            return(5);
    }
    else if (y < ymin - EPSI1)
        return(7);
    else if (y > ymax + EPSI1)
        return(3);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  clip_input      Read polygon from input file PMF2d.                     */
/*                  Return number of points, or -1 if error.                */

int clip_input(double *xmin,double *ymin,double *xmax,double *ymax) 
{
    int i,n,err;
    char buf[101],*p;
    double x,y;

    printf1("Reading input file: %s\n",PMF2dName);

    err = 0;
    if (!fgets(buf,100,PMF2d)) { 
        printf1("Error: cannot read object description record.\n");
        return(-1);
    }               
    if (err == 0) {
        p = skip_b(buf);
        if (sscanf(p,"%d",&n) != 1)
            err = -1;
    }
    if (err == 0) {
        p = skip_int(p);
        p = skip_b(p);
        if (sscanf(p,"%d",&n) != 1 || n != 3)
            err = -1;
    }
    if (err == 0) {
        p = skip_int(p);
        p = skip_b(p);
        if (sscanf(p,"%d",&n) != 1 || n < 3)
            err = -1;
    }
    if (err == -1) {
        printf1("Error: invalid object description record.\n");
        return(-1);
    }
    if (alloc_acx(n))
        return(-1);             
    if (alloc_acy(n))
        return(-1);           

    for (i = 0; i < n; ++i) {
        if (!fgets(buf,100,PMF2d)) {
            printf1("Error: cannot read coordinates record %d.\n",i + 1);
            return(-1);
        }
        p = skip_b(buf);
        if (sscanf(p,"%lf",&x) != 1) {
            printf1("Error: invalid x coordinate in record %d\n",i + 1);
            return(-1);
        }
        p = skip_dbl(p);
        p = skip_b(p);
        if (sscanf(p,"%lf",&y) != 1) {
            printf1("Error: invalid x coordinate in record %d\n",i + 1);
            return(-1);
        }
        AcX[i] = x;
        AcY[i] = y; 

        if (i == 0) {
            *xmin = *xmax = x;
            *ymin = *ymax = y; 
        }
        else {
            *xmin = dmin(*xmin,x);
            *xmax = dmax(*xmax,x);
            *ymin = dmin(*ymin,y);
            *ymax = dmax(*ymax,y);
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

int clip_line(double xa,double ya,double xb,double yb,double xmin,double ymin,
    double xmax,double ymax,double *sax,double *say,double *sbx,double *sby, 
    int *ca,int *cb,int *ra,int *rb,int *dir)
{
    double m;

    *ca = clip_line_code(xa,ya,xmin,ymin,xmax,ymax);
    *cb = clip_line_code(xb,yb,xmin,ymin,xmax,ymax);

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
        if (fabs(xa - xmin) < EPSI1) {
            if (fabs(ya - ymin) < EPSI1)
                *ra = 1;
            else if (fabs(ya - ymax) < EPSI1)
                *ra = 4;
            else
                *ra = -4;
        }
        else if (fabs(xa - xmax) < EPSI1) {
            if (fabs(ya - ymin) < EPSI1)
                *ra = 2;
            else if (fabs(ya - ymax) < EPSI1)
                *ra = 3;
            else
                *ra = -2;
        }
        else if (fabs(ya - ymin) < EPSI1)
            *ra = -1;
        else if (fabs(ya - ymax) < EPSI1)
            *ra = -3;
    }
    if (*cb == 0) {
        if (fabs(xb - xmin) < EPSI1) {
            if (fabs(yb - ymin) < EPSI1)
                *rb = 1;
            else if (fabs(yb - ymax) < EPSI1)
                *rb = 4;
            else
                *rb = -4;
        }
        else if (fabs(xb - xmax) < EPSI1) {
            if (fabs(yb - ymin) < EPSI1)
                *rb = 2;
            else if (fabs(yb - ymax) < EPSI1)
                *rb = 3;
            else
                *rb = -2;
        }
        else if (fabs(yb - ymin) < EPSI1)
            *rb = -1;
        else if (fabs(yb - ymax) < EPSI1)
            *rb = -3;
    }
    if (*ca == 0 && *cb == 0)
        return(1);
      
    if (fabs(xa - xb) < EPSI1) {
        if (xa < xmin || xb > xmax)
            return(0);

        if (*ca != 0) {
            if (ya < ymin) {
                *say = ymin;
                                                        
                if (fabs(xa - xmin) < EPSI1)
                    *ra = 1;
                else if (fabs(xa - xmax) < EPSI1)
                    *ra = 2;
                else
                    *ra = -1;
            }
            else if (ya > ymax) {
                *say = ymax;

                if (fabs(xa - xmin) < EPSI1)
                    *ra = 4;
                else if (fabs(xa - xmax) < EPSI1)
                    *ra = 3;
                else
                    *ra = -3;
            }
        }
        if (*cb != 0) {
            if (yb < ymin) {
                *sby = ymin;

                if (fabs(xb - xmin) < EPSI1)
                    *rb = 1;
                else if (fabs(xb - xmax) < EPSI1)
                    *rb = 2;
                else
                    *rb = -1;
            }
            else if (yb > ymax) {
                *sby = ymax;

                if (fabs(xb - xmin) < EPSI1)
                    *rb = 4;
                else if (fabs(xa - xmax) < EPSI1)
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

        if (fabs(*sax - xmin) < EPSI1) {
            if (fabs(*say - ymin) < EPSI1)
                *ra = 1;
            else if (fabs(*say - ymax) < EPSI1)
                *ra = 4;
            else
                *ra = -4;
        }
        else if (fabs(*sax - xmax) < EPSI1) {
            if (fabs(*say - ymin) < EPSI1)
                *ra = 2;
            else if (fabs(*say - ymax) < EPSI1)
                *ra = 3;
            else
                *ra = -2;
        }
        else if (fabs(*say - ymin) < EPSI1)
            *ra = -1;
        else if (fabs(*say - ymax) < EPSI1)
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

        if (fabs(*sbx - xmin) < EPSI1) {
            if (fabs(*sby - ymin) < EPSI1)
                *rb = 1;
            else if (fabs(*sby - ymax) < EPSI1)
                *rb = 4;
            else
                *rb = -4;
        }
        else if (fabs(*sbx - xmax) < EPSI1) {
            if (fabs(*sby - ymin) < EPSI1)
                *rb = 2;
            else if (fabs(*sby - ymax) < EPSI1)
                *rb = 3;
            else
                *rb = -2;
        }
        else if (fabs(*sby - ymin) < EPSI1)
            *rb = -1;
        else if (fabs(*sby - ymax) < EPSI1)
            *rb = -3;
    }
    if (*ra == *rb && *ra > 0 && *dir == 2)   
        return(0);
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  clip_line_code. Return code for clip_line function.                     */

int clip_line_code(double x,double y,double xmin,double ymin,double xmax,
    double ymax)
{
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



