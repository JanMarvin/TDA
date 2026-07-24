/****************************************************************************/
/*  t_sdr                                                                   */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gf.h"
#include "t_psf.h"
#include "t_plot.h"
#include "t_plot3.h"
#include "t_gm.h"
#include "t_sd.h"
#include "tda_context.h"

/*  functions in t_sdr.c */

int sdnl(TDAContext *ctx);
int sdppol(TDAContext *ctx);
int sdlpol(TDAContext *ctx);
int sdipol(TDAContext *ctx);

int sdcpol(TDAContext *ctx);
int sdp_id(TDAContext *ctx, double xa,double ya,double xb,double yb,double tol);
void sdp_sort(TDAContext *ctx, int i,int k,int n,int *nptr,float *x,float *y,double *a);
int sdp_fpol(TDAContext *ctx, int nn,int nc,int *cn,int rn,float *xc,float *yc,char *idx, int *bnodes,int nbmax,int *pnodes,int npmax,int *sdid,int *nrec);

void sdp_comp(TDAContext *ctx, int k,int nc,int *cn);
int sdp_deg(TDAContext *ctx, int i,int j);       
int sdcpol_edge(TDAContext *ctx, int i,int j,int n,int *bnodes,char *idx);
int sdp_clear(TDAContext *ctx, int nn);
int sdp_cclear(TDAContext *ctx, int nn,int nc,int *cn);
void sdp_remove(TDAContext *ctx, int i);
void sdp_cut(TDAContext *ctx, int nn,int *vn,char *idx);
int sdp_visit(TDAContext *ctx, int i,int ir,int *vn,char *idx);
int sdp_find_idx(TDAContext *ctx, int i,int j);
int sdp_find_nxt(TDAContext *ctx, int i,int l,int opt);
int sdp_write(TDAContext *ctx, int id,int np,int *pnodes,float *xc,float *yc);
int sdp_list(TDAContext *ctx, int n,float *xa,float *ya,float *xb,float *yb,int *ptr,char *idx, int nmax,float tol);
int sdp_list_find(TDAContext *ctx, float xa0,float ya0,int n,float *xa,float *ya,float tol);
void sdp_list_free(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  sdnl    Creating a network from lines.                                  */
/*                                                                          */  
/*          sdnl(                                                           */
/*              opt=...,    option for output file, def. 1                  */
/*                          1 = only node numbers and length                */
/*                          2 = add coordinates                             */
/*              tol=...,    tolerance for identical points, def. 1.e-4      */
/*              max=...,    maximum of nodes, def. 2 * number of segments   */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 10.4          */
/*          ) = output_file;                                                */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. Only type 2 objects (lines) will be used. The command   */
/*  extracts all line segments from these lines and creates an undirected   */
/*  graph. Nodes are created from the end points of the line segments,      */
/*  from crossings, and from overlaps. Two nodes are connected by an edge   */
/*  if they are connected by (at least) one line segment. The edges are     */
/*  valued with their euclidean length.                                     */
/*                                                                          */
/*  The output file is written in the form of an edge list:                 */
/*                                                                          */
/*      i  j  v  xi yi xj yj                                                */
/*                                                                          */
/*  where i and j refer to nodes, v is the edge value and, if opt = 2,      */
/*  (xi,yi) and (xj,yj) are the coordinates of node i and node j,           */
/*  respectively.                                                           */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdnl(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,r,n,nl,ns,nrec;                 
    double len;

    err = -1;
    printf1(ctx, "Creating a network from lines. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto SDNLFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 2)
        ctx->PMOPT = 1;
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    nl = ns = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 2)
            continue;

        n = (int)get_data(ctx, ctx->SDVarSDN,i);             
        nl += 1;
        ns += n - 1;
    }
    newline(ctx);
    printf1(ctx, "Number of lines: %d\n",nl);
    printf1(ctx, "Number of segments: %d\n",ns);
    ctx->PMTol = dmax(ctx, ctx->PMTol,ctx->EPSI1);
    ctx->PMMax = imax(ctx, 2 * ns,ctx->PMMax);
    printf1(ctx, "Tolerance: %g\n",ctx->PMTol);
    printf1(ctx, "Maximal number of segments: %d\n",ctx->PMMax);

    if (alloc_acxf(ctx, ctx->PMMax + 2))
        goto SDNLFin;     
    if (alloc_acyf(ctx, ctx->PMMax + 2))
        goto SDNLFin;     
    if (alloc_acuf(ctx, ctx->PMMax + 2))
        goto SDNLFin;     
    if (alloc_acvf(ctx, ctx->PMMax + 2))
        goto SDNLFin;     
    if (alloc_ack(ctx, ctx->PMMax + 2))
        goto SDNLFin;     
    if (alloc_acd(ctx, ctx->PMMax + 2))
        goto SDNLFin;     

    k = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 2)
            continue;

        if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
            goto SDNLFin;

        for (j = 1; j < n; ++j) {
            ctx->AcXF[k] = (float)ctx->SDVarX[j - 1];
            ctx->AcYF[k] = (float)ctx->SDVarY[j - 1];
            ctx->AcUF[k] = (float)ctx->SDVarX[j];
            ctx->AcVF[k] = (float)ctx->SDVarY[j];
            k++;
        }
    }
    r = sdp_list(ctx, ns,ctx->AcXF,ctx->AcYF,ctx->AcUF,ctx->AcVF,ctx->AcK,ctx->AcD,ctx->PMMax,(float)ctx->PMTol);         
    if (r != 0)  
        goto SDNLFin;

    printf1(ctx, "\nNumber of nodes: %d\n",ctx->SD_NN);
    printf1(ctx, "Number of edges: %d\n",ctx->SD_NE / 2);

    nrec = 0;
    for (i = 1; i <= ctx->SD_NN; ++i) {
        k = ctx->SD_EPtr[i];
        n = ctx->SD_NA[i];
        for (l = 0; l < n; ++l) {
            j = ctx->SD_EList[k + l];
            if (j > i) {
                len = g_len(ctx, (double)ctx->SD_NLX[i],(double)ctx->SD_NLY[i],(double)ctx->SD_NLX[j],(double)ctx->SD_NLY[j]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j);        
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,len);        
                if (ctx->PMOPT == 2) {
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->SD_NLX[i]);        
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->SD_NLY[i]);        
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->SD_NLX[j]);        
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,(double)ctx->SD_NLY[j]);        
                }
                fprintf(ctx->PMFd,"\n");        
                nrec++;
            }
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDNLFin:
    sdp_list_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdppol  Find points in polygons                                         */
/*                                                                          */  
/*          sdppol(                                                         */
/*              opt=...,    option for output file, def. 1                  */
/*                          1 = standard data file                          */
/*                          2 = spatial data file                           */
/*                          3 = edge list of bi-modal graph                 */
/*              n=...,      number of polygon IDs, def. 1                   */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 10.4          */
/*          ) = output_file;                                                */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. There must be least one point and one polygon.          */
/*  For each point i and polygon j, the command checks whether the point    */
/*  is inside the polygon (including its boundaries). There are three       */
/*  output options.                                                         */
/*                                                                          */
/*  Option 1. The output file contains one record for each point withe      */
/*  the following record structure:                                         */
/*                                                                          */
/*      SDID  x  y  np  P1  P2 ...                                          */
/*                                                                          */
/*  SDID is the ID of the point in the spatial data structure, x and y are  */
/*  its coordinates. np is the number of polygons that contain the points.  */
/*  Then follow up to n columns containing the SDIDs of the polygons, or -1 */
/*  if the point is not inside. n can be specified with the n parameter.    */  
/*                                                                          */
/*  Option 2. The same information is provided in the form of a spatial     */
/*  data file:                                                              */
/*                                                                          */
/*      SDID  1  1  np  P1  P2  ...                                         */
/*          x y                                                             */
/*                                                                          */
/*  Option 3. The output file is an edge list for a bi-modal graph. There   */
/*  is one record for each point i and polygon j, given that i is inside j. */
/*                                                                          */
/*      i  j  x  y                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdppol(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,n,np,npol,nt,nrec,typ;
    double xmin,xmax,ymin,ymax;
    double x,y;

    err = -1;
    nrec = 0;
    printf1(ctx, "Find points in polygons. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto SDPPFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
        ctx->PMOPT = 1;
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);
    if (ctx->PMN < 1)
        ctx->PMN = 1;

    np = npol = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        typ = (int)get_data(ctx, ctx->SDVarSDTyp,i);           
        if (typ == 1)
            np++;
        else if (typ == 3)
            npol++;
    }
    newline(ctx);
    printf1(ctx, "Number of points: %d\n",np);
    printf1(ctx, "Number of polygons: %d\n",npol);
    if (np < 1 || npol < 1) {
        err = 0;
        goto SDPPFin;
    }   
    if (ctx->PMN > npol)
        ctx->PMN = npol;

    /* save points in AcX, AcY, id in AcI; boundaries for polygons in 
       AcXF,AxYF,AcUF,AcVF, index in AcJ. AcK is used to save polygon ids 
       for the points. */

    if (alloc_acx(ctx, np))
        goto SDPPFin;     
    if (alloc_acy(ctx, np))
        goto SDPPFin;     
    if (alloc_aci(ctx, np))
        goto SDPPFin;     

    if (alloc_acxf(ctx, npol))
        goto SDPPFin;     
    if (alloc_acyf(ctx, npol))
        goto SDPPFin;     
    if (alloc_acuf(ctx, npol))
        goto SDPPFin;     
    if (alloc_acvf(ctx, npol))
        goto SDPPFin;     
    if (alloc_acj(ctx, npol))
        goto SDPPFin;     

    if (alloc_ack(ctx, npol))
        goto SDPPFin;     

    i = j = 0;
    for (k = 0; k < ctx->NOC; ++k) {
        typ = (int)get_data(ctx, ctx->SDVarSDTyp,k);           
        if (typ != 1 && typ != 3)  
            continue;

        if ((n = sd_getdata(ctx, k,0,0,1)) < 1)
            goto SDPPFin;

        if (typ == 1) {
            ctx->AcI[i] = (int)get_data(ctx, ctx->SDVarSDID,k);           
            ctx->AcX[i] = ctx->SDVarX[0];
            ctx->AcY[i] = ctx->SDVarY[0];
            i++;
        }
        else {
            ctx->AcJ[j] = k;                                    
            xmax = ctx->SDVarX[0];
                                    
            xmin = (double)xmax;
            ymax = ctx->SDVarY[0];

            ymin = (double)ymax;
            for (l = 1; l < n; ++l) {
                xmin = (double)((float)dmin(ctx, (double)xmin,ctx->SDVarX[l]));
                xmax = (double)((float)dmax(ctx, (double)xmax,ctx->SDVarX[l]));
                ymin = (double)((float)dmin(ctx, (double)ymin,ctx->SDVarY[l]));
                ymax = (double)((float)dmax(ctx, (double)ymax,ctx->SDVarY[l]));
            }   
            ctx->AcXF[j] = (float)(xmin - ctx->EPSI1);
            ctx->AcUF[j] = (float)(xmax + ctx->EPSI1);
            ctx->AcYF[j] = (float)(ymin - ctx->EPSI1);
            ctx->AcVF[j] = (float)(ymax + ctx->EPSI1);
            j++;
        }
    }
    for (i = 0; i < np; ++i) {

        nt = 0;                     /* number of insides */
        x = ctx->AcX[i];
        y = ctx->AcY[i];

        for (j = 0; j < npol; ++j) {
            xmin = (double)(ctx->AcXF[j]);
            xmax = (double)(ctx->AcUF[j]);
            ymin = (double)(ctx->AcYF[j]);
            ymax = (double)(ctx->AcVF[j]);

            if ((double)((float)(x)) >= (double)(xmin) && (double)((float)(x)) <= (double)(xmax) &&
                (double)((float)(y)) >= (double)(ymin) && (double)((float)(y)) <= (double)(ymax)) {  

                if ((n = sd_getdata(ctx, ctx->AcJ[j],0,1,1)) < 1)
                    goto SDPPFin;
                              
                if (g_inpoly(ctx, n - 1,ctx->SDVarX,ctx->SDVarY,x,y)) {

                    k = (int)get_data(ctx, ctx->SDVarSDID,ctx->AcJ[j]);           

                    if (ctx->PMOPT == 3) {
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcI[i]);        
                        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,k);        
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);        
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[i]);        
                        fprintf(ctx->PMFd,"\n");        
                        nrec++;
                    }
                    else 
                        ctx->AcK[nt++] = k;                                         
                }
            }
        }
        if (ctx->PMOPT == 3)
            continue;

        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcI[i]);        
        if (ctx->PMOPT == 1) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);        
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[i]);        
        }
        else {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,1);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,1);        
        }                                
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nt);        
        n = imin(ctx, nt,ctx->PMN);
        for (j = 0; j < n; ++j)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcK[j]);        
        for (; j < ctx->PMN; ++j)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,-1);        
        fprintf(ctx->PMFd,"\n");        
        nrec++;
        if (ctx->PMOPT == 2) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[i]);        
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[i]);        
            fprintf(ctx->PMFd,"\n");        
            nrec++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDPPFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdlpol  Find lines in polygons.                                         */
/*                                                                          */  
/*          sdlpol(                                                         */
/*              opt=...,    option for output file, def. 1                  */
/*                          1 = standard data file                          */
/*                          2 = spatial data file                           */
/*                          3 = edge list of bi-modal graph                 */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 10.4          */
/*          ) = output_file;                                                */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. There must be least one point and one polygon.          */
/*  For each point i and polygon j, the command checks whether the point    */
/*  is inside the polygon (including its boundaries). There are three       */
/*  output options.                                                         */
/*                                                                          */
/*  Option 1. The output file contains one record for each point withe      */
/*  the following record structure:                                         */
/*                                                                          */
/*      SDID  x  y  np  P1  P2 ...                                          */
/*                                                                          */
/*  SDID is the ID of the point in the spatial data structure, x and y are  */
/*  its coordinates. np is the number of polygons that contain the points.  */
/*  Then follow up to n columns containing the SDIDs of the polygons, or -1 */
/*  if the point is not inside. n can be specified with the n parameter.    */  
/*                                                                          */
/*  Option 2. The same information is provided in the form of a spatial     */
/*  data file:                                                              */
/*                                                                          */
/*      SDID  1  1  np  P1  P2  ...                                         */
/*          x y                                                             */
/*                                                                          */
/*  Option 3. The output file is an edge list for a bi-modal graph. There   */
/*  is one record for each point i and polygon j, given that i is inside j. */
/*                                                                          */
/*      i  j  x  y                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdlpol(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,r,n,nl,npol,npl,nt,nrec,typ,nlmax,idl,idp;
    double xmin,xmax,ymin,ymax;

    err = -1;
    nrec = 0;
    printf1(ctx, "Find lines in polygons. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto SDLPFin;

    if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
        ctx->PMOPT = 1;
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    nlmax = nl = npol = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        typ = (int)get_data(ctx, ctx->SDVarSDTyp,i);           
        if (typ == 2) {
            nl++;
            nlmax = imax(ctx, nlmax,(int)get_data(ctx, ctx->SDVarSDN,i));
        }
        else if (typ == 3)
            npol++;
    }
    newline(ctx);
    printf1(ctx, "Number of lines: %d\n",nl);
    printf1(ctx, "Number of polygons: %d\n",npol);
    if (nl < 1 || npol < 1) {
        err = 0;
        goto SDLPFin;
    }   

    /* save boundaries for polygons in AcXF,AxYF,AcUF,AcVF, index in AcJ.
       AcX, AcY is used to save the points of a line. AcI is used for the 
       indices of the lines. */

    if (alloc_acx(ctx, nlmax + 1))
        goto SDLPFin;     
    if (alloc_acy(ctx, nlmax + 1))
        goto SDLPFin;     
    if (alloc_aci(ctx, nl))
        goto SDLPFin;     

    if (alloc_acxf(ctx, npol))
        goto SDLPFin;     
    if (alloc_acyf(ctx, npol))
        goto SDLPFin;     
    if (alloc_acuf(ctx, npol))
        goto SDLPFin;     
    if (alloc_acvf(ctx, npol))
        goto SDLPFin;     
    if (alloc_acj(ctx, npol))
        goto SDLPFin;     

    i = j = 0;
    for (k = 0; k < ctx->NOC; ++k) {
        typ = (int)get_data(ctx, ctx->SDVarSDTyp,k);           
        if (typ == 2) 
            ctx->AcI[i++] = k;
        if (typ != 3)  
            continue;

        if ((n = sd_getdata(ctx, k,0,0,1)) < 1)
            goto SDLPFin;

        ctx->AcJ[j] = k;                                    
        xmax = ctx->SDVarX[0];
                                    
        xmin = (double)xmax;
        ymax = ctx->SDVarY[0];

        ymin = (double)ymax;
        for (l = 1; l < n; ++l) {
            xmin = (double)((float)dmin(ctx, (double)xmin,ctx->SDVarX[l]));
            xmax = (double)((float)dmax(ctx, (double)xmax,ctx->SDVarX[l]));
            ymin = (double)((float)dmin(ctx, (double)ymin,ctx->SDVarY[l]));
            ymax = (double)((float)dmax(ctx, (double)ymax,ctx->SDVarY[l]));
        }   
        ctx->AcXF[j] = (float)(xmin - ctx->EPSI1);
        ctx->AcUF[j] = (float)(xmax + ctx->EPSI1);
        ctx->AcYF[j] = (float)(ymin - ctx->EPSI1);
        ctx->AcVF[j] = (float)(ymax + ctx->EPSI1);
        j++;
    }
    for (i = 0; i < nl; ++i) {

        idl = (int)get_data(ctx, ctx->SDVarSDID,ctx->AcI[i]);           

        if ((npl = sd_getdata(ctx, ctx->AcI[i],0,0,1)) < 1)
            goto SDLPFin;

        for (j = 0; j < npl; ++j) {
            ctx->AcX[j] = ctx->SDVarX[j];
            ctx->AcY[j] = ctx->SDVarY[j];
        }
        nt = 0;                     /* number of insides */
        r = 0;
        for (k = 0; k < npol; ++k) {
            xmin = (double)(ctx->AcXF[k]);
            xmax = (double)(ctx->AcUF[k]);
            ymin = (double)(ctx->AcYF[k]);
            ymax = (double)(ctx->AcVF[k]);

            r = 1;
            for (j = 0; j < npl; ++j) {
                if ((double)((float)(ctx->AcX[j])) < (double)(xmin) || (double)((float)(ctx->AcX[j])) > (double)(xmax) ||
                    (double)((float)(ctx->AcY[j])) < (double)(ymin) || (double)((float)(ctx->AcY[j])) > (double)(ymax)) {   
                    r = 0;    
                    break;
                }
            }
            if (r == 0)
                continue;

            if ((n = sd_getdata(ctx, ctx->AcJ[k],0,1,1)) < 1)
                goto SDLPFin;
                              
            for (j = 1; j < npl; ++j) {
                if (g_l_inpoly(ctx, n - 1,ctx->SDVarX,ctx->SDVarY,ctx->AcX[j-1],ctx->AcY[j-1],ctx->AcX[j],ctx->AcY[j]) == 0) {
                    r = 0;
                    break;
                }
            }
            if (r == 0)
                continue;

            idp = (int)get_data(ctx, ctx->SDVarSDID,ctx->AcJ[k]);           

            if (ctx->PMOPT == 3) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idl);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idp);        
                fprintf(ctx->PMFd,"\n");        
                nrec++;
            }
            break;
        }
        if (ctx->PMOPT == 3)
            continue;

        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idl);        
        if (ctx->PMOPT == 2) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,2);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,npl);        
        }                                
        if (r == 0)
            idp = -1;


        n = imin(ctx, nt,ctx->PMN);
        for (j = 0; j < n; ++j)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcK[j]);        
        for (; j < ctx->PMN; ++j)
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,-1);        
        fprintf(ctx->PMFd,"\n");        
        nrec++;
        if (ctx->PMOPT == 2) {
            for (j = 0; j < npl; ++j) {
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j]);        
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[j]);        
                fprintf(ctx->PMFd,"\n");        
                nrec++;
            }
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDLPFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  sdipol  Hierarchical inclusion of polygons.                             */
/*                                                                          */  
/*          sdipol(                                                         */
/*              opt=...,    option for output file, def. 1                  */
/*                          1 = standard data file                          */
/*                          2 = spatial data file                           */
/*                          3 = edge list of bi-modal graph                 */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 10.4          */
/*          ) = output_file;                                                */
/*                                                                          */
/*  This command requires a valid spatial data structure defined with the   */
/*  sdnvar command. There must be least one point and one polygon.          */
/*  For each point i and polygon j, the command checks whether the point    */
/*  is inside the polygon (including its boundaries). There are three       */
/*  output options.                                                         */
/*                                                                          */
/*  Option 1. The output file contains one record for each point withe      */
/*  the following record structure:                                         */
/*                                                                          */
/*      SDID  x  y  np  P1  P2 ...                                          */
/*                                                                          */
/*  SDID is the ID of the point in the spatial data structure, x and y are  */
/*  its coordinates. np is the number of polygons that contain the points.  */
/*  Then follow up to n columns containing the SDIDs of the polygons, or -1 */
/*  if the point is not inside. n can be specified with the n parameter.    */  
/*                                                                          */
/*  Option 2. The same information is provided in the form of a spatial     */
/*  data file:                                                              */
/*                                                                          */
/*      SDID  1  1  np  P1  P2  ...                                         */
/*          x y                                                             */
/*                                                                          */
/*  Option 3. The output file is an edge list for a bi-modal graph. There   */
/*  is one record for each point i and polygon j, given that i is inside j. */
/*                                                                          */
/*      i  j  x  y                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdipol(TDAContext *ctx)
{
    register int i,j,l;
    int err,r,n,npol,nrec,npmax,idi,idj,npi,npj,ii,jj;
    double xmin,xmax,ymin,ymax;

    err = -1;
    nrec = 0;
    printf1(ctx, "Hierarchical inclusion of polygons. Current memory: %d bytes.\n",ctx->MemReq);
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto SDIPFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    npmax = npol = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 3)      
            continue;
        npol++;
        npmax = imax(ctx, npmax,(int)get_data(ctx, ctx->SDVarSDN,i));
    }
    printf1(ctx, "\nNumber of polygons: %d\n",npol);
    if (npol < 1) {
        err = 0;
        goto SDIPFin;
    }   

    /* save boundaries for polygons in AcXF,AxYF,AcUF,AcVF, index in AcN. */

    if (alloc_acxf(ctx, npol))
        goto SDIPFin;     
    if (alloc_acyf(ctx, npol))
        goto SDIPFin;     
    if (alloc_acuf(ctx, npol))
        goto SDIPFin;     
    if (alloc_acvf(ctx, npol))
        goto SDIPFin;     
    if (alloc_acn(ctx, npol))
        goto SDIPFin;     
    if (alloc_acd(ctx, ctx->NOC + 1))
        goto SDIPFin;     

    j = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 3)      
            continue;

        if ((n = sd_getdata(ctx, i,0,0,1)) < 1)
            goto SDIPFin;

        ctx->AcN[j] = i;                                    
        xmax = ctx->SDVarX[0];
                                    
        xmin = (double)xmax;
        ymax = ctx->SDVarY[0];

        ymin = (double)ymax;
        for (l = 1; l < n; ++l) {
            xmin = (double)((float)dmin(ctx, (double)xmin,ctx->SDVarX[l]));
            xmax = (double)((float)dmax(ctx, (double)xmax,ctx->SDVarX[l]));
            ymin = (double)((float)dmin(ctx, (double)ymin,ctx->SDVarY[l]));
            ymax = (double)((float)dmax(ctx, (double)ymax,ctx->SDVarY[l]));
        }   
        ctx->AcXF[j] = (float)(xmin);
        ctx->AcUF[j] = (float)(xmax);
        ctx->AcYF[j] = (float)(ymin);
        ctx->AcVF[j] = (float)(ymax);
        j++;
    }

    if (alloc_acx(ctx, npmax + 2))
        goto SDIPFin;     
    if (alloc_acy(ctx, npmax + 2))
        goto SDIPFin;     

    for (i = 0; i < npol; ++i) {
        ii = ctx->AcN[i];
        idi = (int)get_data(ctx, ctx->SDVarSDID,ii);           

        if ((npi = sd_getdata(ctx, ii,0,1,1)) < 1)
            goto SDIPFin;

        for (l = 0; l < npi; ++l) {
            ctx->AcX[l] = ctx->SDVarX[l];
            ctx->AcY[l] = ctx->SDVarY[l];
        }
        xmin = (double)(ctx->AcXF[i]);
        xmax = (double)(ctx->AcUF[i]);
        ymin = (double)(ctx->AcYF[i]);
        ymax = (double)(ctx->AcVF[i]);

        for (j = i + 1; j < npol; ++j) {
            if ((double)(xmax) < (double)(ctx->AcXF[j]) || (double)(xmin) > (double)(ctx->AcUF[j]) ||
                (double)(ymax) < (double)(ctx->AcYF[j]) || (double)(ymin) > (double)(ctx->AcVF[j]))    
                continue;

            jj = ctx->AcN[j];
            if ((npj = sd_getdata(ctx, jj,0,1,1)) < 1)
                goto SDIPFin;

            r = g_p_inpoly(ctx, npi,ctx->AcX,ctx->AcY,npj,ctx->SDVarX,ctx->SDVarY); 
            if (r == 0) {
                r = g_p_inpoly(ctx, npj,ctx->SDVarX,ctx->SDVarY,npi,ctx->AcX,ctx->AcY);
                if (r == 1)
                    r = 2;
            }
            if (r == 0)
                continue;

            idj = (int)get_data(ctx, ctx->SDVarSDID,jj);           

            if (r == 1) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idi);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idj);        
            }
            else {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idj);        
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idi);        
            }
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,r);        
            fprintf(ctx->PMFd,"\n");        
            nrec++;

            ctx->AcD[ii] = 1;
            ctx->AcD[jj] = 1;
        }
    }

    /* add isolated nodes */

    for (i = 0; i < npol; ++i) {
        ii = ctx->AcN[i];
        if (ctx->AcD[ii] == 0) {
            idi = (int)get_data(ctx, ctx->SDVarSDID,ii);           
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idi);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,idi);        
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,-1);        
            fprintf(ctx->PMFd,"\n");        
            nrec++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

SDIPFin:
    p_clean(ctx);
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  sdcpol      Construct polygons from line segments.                      */
/*                                                                          */  
/*              sdcpol(                                                     */
/*                  opt=...,    option, def. 1                              */
/*                              1 print node list of graph                  */
/*                              2 print edge list of graph                  */
/*                              3 construction of polygons                  */
/*                  df=...,     output file (required)                      */
/*                  tol=...,    tolerance for identical points, def. 1.e-4  */
/*                  max=...,    maximum of nodes, def. 2 * NOC              */
/*                  nfmt=...,   integer print format, def. 4                */
/*                  fmt=...,    print format for coordinates, def. 10.4     */
/*              ) = X1,Y1,X2,Y2;                                            */
/*                                                                          */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdcpol(TDAContext *ctx)
{
    register int i,j,k,l;
    int err,r,ixa,iya,ixb,iyb,n,m,nn,nrec,nc,ne,sdid,npmax,bnmax;
    double x,y;

    err = -1;
    nrec = 0;
    printf1(ctx, "Construction of polygons from line segments. Current memory: %d bytes.\n",ctx->MemReq);
         
    if (parm(ctx, ctx->CmdBuf + 6,4,1))       /* get parameters */
        goto SDCPFin;

    if (ctx->PMNV != 4) {
        printf1(ctx, "Error: need four variables on right-hand side.\n");
        goto SDCPFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto SDCPFin;
    }
    if (ctx->NOC < 3) {
        printf1(ctx, "Error: need at least three line segments.\n");
        goto SDCPFin;
    }
    if (ctx->PMOPT < 1 || ctx->PMOPT > 3)
        ctx->PMOPT = 1;
    if (ctx->PMMax < 1)
        ctx->PMMax = 2 * ctx->NOC;

    ctx->PMTol = dmax(ctx, ctx->PMTol,ctx->EPSI1);
    ctx->PMMax = imax(ctx, 2 * ctx->NOC,ctx->PMMax);
    printf1(ctx, "Tolerance: %g\n",ctx->PMTol);
    printf1(ctx, "Maximal number of segments: %d\n",ctx->PMMax);

    /* Create the graph. First read segments, then call sdp_list(). */

    if (alloc_acxf(ctx, ctx->PMMax + 2))
        goto SDCPFin;     
    if (alloc_acyf(ctx, ctx->PMMax + 2))
        goto SDCPFin;     
    if (alloc_acuf(ctx, ctx->PMMax + 2))
        goto SDCPFin;     
    if (alloc_acvf(ctx, ctx->PMMax + 2))
        goto SDCPFin;     
    if (alloc_ack(ctx, ctx->PMMax + 2))
        goto SDCPFin;     
    if (alloc_acd(ctx, ctx->PMMax + 2))
        goto SDCPFin;     

    ixa = ctx->PMVIdx[0];
    iya = ctx->PMVIdx[1];
    ixb = ctx->PMVIdx[2];
    iyb = ctx->PMVIdx[3];

    for (i = 0; i < ctx->NOC; ++i) {
        ctx->AcXF[i] = (float)get_data(ctx, ixa,i);
        ctx->AcYF[i] = (float)get_data(ctx, iya,i);
        ctx->AcUF[i] = (float)get_data(ctx, ixb,i);
        ctx->AcVF[i] = (float)get_data(ctx, iyb,i);
    }
    r = sdp_list(ctx, ctx->NOC,ctx->AcXF,ctx->AcYF,ctx->AcUF,ctx->AcVF,ctx->AcK,ctx->AcD,ctx->PMMax,(float)ctx->PMTol);         
    if (r != 0)  
        goto SDCPFin;

    alloc_acxf(ctx, 0);   /* free memory */
    alloc_acyf(ctx, 0);
    alloc_acuf(ctx, 0);
    alloc_acvf(ctx, 0);
    alloc_ack(ctx, 0);
    alloc_acd(ctx, 0);
      
    nn = ctx->SD_NN;     /* number of nodes */
    ne = ctx->SD_NE;     /* number of edges */

    printf1(ctx, "Number of nodes: %d\n",nn);
    printf1(ctx, "Number of nodes: %d\n",ne);
    if (nn < 3 || ne < 3) {
        printf1(ctx, "Will not continue.\n");
        err = 0;
        goto SDCPFin;
    }
    if (ctx->PMOPT == 1) {                       /* print node list */
        for (i = 1; i <= nn; ++i) {
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLX[i]);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLY[i]);
            fprintf(ctx->PMF1d,"\n");
            nrec++;
        }
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
        err = 0;
        goto SDCPFin;
    }

    /* for each node i = 1,...,nn, sort the adjacent nodes in
       counterclockwise direction */

    printf1(ctx, "Sorting nodes in counterclockwise orientation.\n");

    if (alloc_actmp(ctx, nn + 1))
        goto SDCPFin;     

    for (i = 1; i <= nn; ++i) {
        if (ctx->SD_NA[i] <= 2)                /* no sorting necessary */
            continue;
        sdp_sort(ctx, i,ctx->SD_EPtr[i],ctx->SD_NA[i],ctx->SD_EList,ctx->SD_NLX,ctx->SD_NLY,ctx->AcTmp);    
    }
    alloc_actmp(ctx, 0);     /* free */

    /* check for multiple edges (segments) */

    n = 0;
    for (i = 1; i <= nn; ++i) {
        k = ctx->SD_EPtr[i];
        j = ctx->SD_EList[k];
        for (l = 1; l < ctx->SD_NA[i]; ++l) {
            if (ctx->SD_EList[k + l] == j) {
                ctx->SD_EList[k + l] = -1;
                n++;
            }
            else
                j = ctx->SD_EList[k + l];
        }
    }
    if (n > 0) {
        printf1(ctx, "Dropped %d segment(s) occurring at least twice.\n",n);
        printf1(ctx, "Remaining number of segments: %d\n",ne - n);
    }


    if (ctx->PMOPT == 2) {                   /* write edge list */
        for (i = 1; i <= nn; ++i) {
            k = ctx->SD_EPtr[i];
            for (l = 0; l < ctx->SD_NA[i]; ++l) {
                j = ctx->SD_EList[k + l];
                if (j < 1)
                    continue;

                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,i);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,j);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLX[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLY[i]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLX[j]);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)ctx->SD_NLY[j]);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
        }
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
        err = 0;
        goto SDCPFin;
    }
                   
    /* Prepare search for polygons. First remove nodes with degree less
       than two. SD_Deg will be used to store the degrees. */

    printf1(ctx, "Searching for polygons.\n");
    printf1(ctx, "Removing nodes with degree less than 2.\n");

    if (alloc_aci(ctx, nn + 1))
        goto SDCPFin;     

    /* create the initial list of degrees in SD_Deg */

    for (i = 1; i <= nn; ++i) {
        k = ctx->SD_EPtr[i];
        j = ctx->SD_EList[k];
        n = 0;
        for (l = 0; l < ctx->SD_NA[i]; ++l) {
            if (ctx->SD_EList[k + l] >= 1)   
                n++;
        }
        ctx->SD_Deg[i] = n;
    }
    if (sdp_clear(ctx, nn) < 0) {
        printf1(ctx, "Fatal error in degree calculation.\n");
        goto SDCPFin;
    }
    n = 0;
    for (i = 1; i <= nn; ++i) {
        if (ctx->SD_Deg[i] >= 2)
            n++;             
    }
    printf1(ctx, "Remaining number of nodes: %d\n",n);
    if (n < 3)
        goto SDCPFin;

    /* ### Remove all cut nodes with degree less than 3. */

    printf1(ctx, "Removing cut nodes with degree less than 3.\n");

    if (alloc_acl(ctx, nn + 1))      /* work space */
        goto SDCPFin;     
    if (alloc_acd(ctx, nn + 1))      
        goto SDCPFin;     
                   
    sdp_cut(ctx, nn,ctx->AcL,ctx->AcD);    

    /* Remove all nodes with degree less than 2. */

    if (sdp_clear(ctx, nn) < 0) {
        printf1(ctx, "Fatal error in degree calculation.\n");
        goto SDCPFin;
    }
    /****************************
    tda_out("X NEUE DEGREES\n");
    for (i = 1; i <= nn; ++i) {
        tda_out("i=%3d aci=%4d\n",i,SD_Deg[i]);
    }   
    tda_out("NEW SD_EList\n");
    for (i = 1; i <= nn; ++i) {
        k = SD_EPtr[i];
        for (l = 0; l < SD_NA[i]; ++l) {
            j = SD_EList[k + l];
            tda_out("i=%3d k=%3d j=%4d\n",i,k,j);
        }
    }
    *************/


    /* Loop over all components, nc is counter for components. The index
       number of the component is stored in AcL. */

    npmax = 2 * nn;             /* max number of nodes in a polygon */
    bnmax = 2 * nn;             /* max number of nodes on boundary */

    if (alloc_acm(ctx, npmax + 2))   /* used to save nodes of polygons */
        goto SDCPFin;     
    if (alloc_acs(ctx, bnmax + 2))   /* used to save nodes on boundary */
        goto SDCPFin;     

    n = 0;
    for (i = 1; i <= nn; ++i) {
        if (ctx->SD_Deg[i] < 2)
            ctx->AcL[i] = -1;
        else {
            ctx->AcL[i] = 0;
            n++;
        }
    }
    printf1(ctx, "Remaining number of nodes: %d.\n\n",n);                  
    if (n < 3) {
        err = 0;
        goto SDCPFin;
    }           
                  
    /* Create the components and do further work for each 
       component separately. */

    sdid = 0;               /* counts the polygons */
    nc = 0;                 /* counter for components */
    while (1) {
        k = -1;                             /* find first unprocessed node */
        for (i = 1; i <= nn; ++i) {
            if (ctx->AcL[i] == 0) {
                k = i;
                break;
            }
        }
        if (k < 0)              /* no more components */
            break;

        ctx->SDC_N = 0;              /* number of elements */
        nc++;                   /* counts the components */
        sdp_comp(ctx, k,nc,ctx->AcL);     /* find next component */

        printf1(ctx, "Component: %d. Nodes: %d. Polygons: ",nc,ctx->SDC_N);
        if (ctx->SDC_N < 3) {
            printf1(ctx, "0\n");
            continue;
        }

        /* find left-most node in current component to begin the search. */

        k = -1;
        x = y = 0.0;
        for (i = 1; i <= nn; ++i) {
            if (ctx->AcL[i] != nc)
                continue;

            if (k < 0 || ctx->SD_NLX[i] < (float)x    
                      || (ctx->SD_NLX[i] == (float)x && ctx->SD_NLY[i] > (float)y)) {
                x = (double)ctx->SD_NLX[i];
                y = (double)ctx->SD_NLY[i];
                k = i;
            }
        }
        tda_out("k====%d\n",k);

        r = sdp_fpol(ctx, nn,nc,ctx->AcL,k,ctx->SD_NLX,ctx->SD_NLY,ctx->AcD,ctx->AcS,bnmax,ctx->AcM,npmax,&sdid,&m);
        if (r < 0) {
            printf1(ctx, "Error: %d\n",r);
            goto SDCPFin;
        }
        else {
            printf1(ctx, "Polygons: %d\n",r);
            nrec += m;
        }

    }
    printf1(ctx, "Found %d polygon(s).\n",sdid);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    err = 0;

SDCPFin:
    sdp_list_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdp_id(xa,ya,xb,yb,tol)                                                 */
/*                                                                          */
/*  Return 1 if (xa,ya) and (xb,yb) are identical (up to tol), otherwise    */  
/*  return 0.                                                               */

int sdp_id(TDAContext *ctx, double xa,double ya,double xb,double yb,double tol)
{
    (void)ctx;        /* unused: the signature is shared */
    xb -= xa;
    yb -= ya;
    if (xb * xb + yb * yb <= tol)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  sdp_sort(i,k,n,nptr,x,y,a)                                              */
/*                                                                          */
/*  Sort the nodes that are adjacent to node i in counterclockwise          */
/*  direction. Note: it is assumed that no edge has zero length.            */

void sdp_sort(TDAContext *ctx, int i,int k,int n,int *nptr,float *x,float *y,double *a)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j,l;
    int fin;
    double x0,y0,xx,yy,tmp;                 

    x0 = (double)x[i];
    y0 = (double)y[i];

    for (j = 0; j < n; ++j) {
        l = nptr[k + j];
        xx = (double)x[l] - x0;
        yy = (double)y[l] - y0;
        a[j] = atan2(yy,xx);
        if (a[j] < 0.0)
            a[j] += 2.0 * Pi;
    }
    fin = 0;
    while (fin == 0) {
        fin = 1;
        for (j = 1; j < n; ++j) {
            if (a[j] < a[j - 1]) {
                tmp = a[j];
                a[j] = a[j - 1];
                a[j - 1] = tmp;
                l = nptr[k + j - 1];
                nptr[k + j - 1] = nptr[k + j];
                nptr[k + j] = l;
                fin = 0;
            }
        }
    }
}    

/* -###-------------------------------------------------------------------- */
/*  sdp_fpol(nn,nc,cn,rn,xc,yc,idx,bnodes,nbmax,pnodes,npmax,sdid,nrec)     */
/*                                                                          */
/*  nn is number of nodes in the graph.                                     */
/*  nc is level of current component. the component consists of all nodes   */
/*  i for which cn[i] = nc.                                                 */
/*  rn is the root node for the boundary.                                   */
/*  x[] and y[] contain the coordinates of the nodes                        */
/*  idx[] is work array of length nn.                                       */  
/*  bnodes[] is work array of length nbmax.                                 */  
/*  pnodes[] is work array of length npmax.                                 */  
/*                                                                          */  
/*  sdid counts the number of polygons and is used as the SDID for the      */
/*  spatial data file written to: PMF1d. Return number of records in nrec.  */
/*                                                                          */
/*  It is assumed that each node in the component has degree >= 2 and that  */
/*  the edge list is sorted in counterclockwise orientation.                */
/*                                                                          */
/*  Return >= 0 number of polygons, or if error:                            */
/*  -1 if a node has degree less than 2                                     */
/*  -2 error in algorithm.                                                  */
/*  -6 overflow in pnodes                                                   */
/*  -8 overflow in bnodes                                                   */

int sdp_fpol(TDAContext *ctx, int nn,int nc,int *cn,int rn,float *xc,float *yc,char *idx, int *bnodes,int nbmax,int *pnodes,int npmax,int *sdid,int *nrec)
{
    register int i,j,k,l,ll;
    int r,ii,n,nnc = 0,i0,j0,j1,j2,nb,np,rnl;
    double x,y,xm,ym,a0,a;

    r = 0;
    *nrec = 0;


    while (1) {             /* do until all polygons have been found */

        /* find the boundary beginning at root node rn in downward
           direction. save the boundary nodes in bnodes[] and flag in idx[]. */

        nb = 1;
        bnodes[nb] = rn;
        idx[rn] = 1;

        /* first find the left-most adjacent node to rn. we need to come
           back from this node. */

        k = ctx->SD_EPtr[rn];
        n = ctx->SD_NA[rn];

        rnl = -1;
        a0 = 0.0;
        for (ll = 0; ll < n; ++ll) {
            if ((j = ctx->SD_EList[k + ll]) >= 1) { 
                a = atan2((double)(yc[j] - yc[rn]),(double)(xc[j] - xc[rn]));
                if (rnl < 0 || a0 < a) {
                    rnl = j;
                    a0 = a;
                }
            }
        }
tda_out("rn=%d     most left rnl=%d a0=%f\n",rn, rnl,a0);

        if (rnl < 0) {
            tda_out("ERROR??? rn=%d  rnl=%d\n",rn, rnl);
            break;
        }


        /* find second node j0 */

        j0 = -1;
        xm = (double)(xc[rn]);
        ym = (double)(yc[rn]);
        a0 = 0.0;
        k = ctx->SD_EPtr[rn];            /* consider all nodes adjacent to rn */
        n = ctx->SD_NA[rn];
        for (ll = 0; ll < n; ++ll) {
            j = ctx->SD_EList[k + ll];        
            if (j < 1)
                continue;

            a = atan2((double)(yc[j]) - ym,(double)(xc[j]) - xm);
            if (j0 < 0 || a < a0) {
                j0 = j;
                a0 = a;
            }
        }
        nb++;
        bnodes[nb] = j0;
        idx[j0] = 1;

               
        tda_out("found second j0=%d\n",j0);

        if (j0 < 1) {
            tda_out("ERROR???\n");
            break;
        }


        /* find remaining nodes on the boundary */

        j1 = j0;
        j0 = rn;
        while (1) {
            if ((l = sdp_find_idx(ctx, j1,j0)) < 0)      /* index to edge list */
                return(-2);
            if ((j2 = sdp_find_nxt(ctx, j1,l,0)) < 0)    /* next node */
                return(-3);

            tda_out("found j2=%d\n",j2);

            /* must not come back to the same node */

            if (j2 == j0)  
                return(-10);

            if (j2 == rn) {
                if (j1 == rnl)          /* no need to go further */
                    break;
            }
            if (nb >= nbmax) 
                return(-8);

            nb++;
            bnodes[nb] = j2;
            idx[j2] = 1;
            j0 = j1;
            j1 = j2;
        }
        tda_out("\nBNodes: ");
        for (i = 1; i <= nb; ++i)
            tda_out("%d ",bnodes[i]);
        newline(ctx);


        /* find all polygons that have at least one segment on the boundary */
        /* store nodes of polygons in pnodes[], counter is np. */
       
        for (ii = 1; ii < nb; ++ii) {

            j0 = i0 = bnodes[ii];
            if (ctx->SD_Deg[i0] < 2) {
                tda_out("xcontinue\n");
                continue;
            }
            j1 = bnodes[ii + 1];
            k = ctx->SD_EPtr[j0];
            n = ctx->SD_NA[j0];
            l = -1;
            for (ll = 0; ll < n; ++ll) {
                if (ctx->SD_EList[k + ll] == j1) {
                    l = ll;
                    break;
                }
            }
            if (l < 0) {
                tda_out("continue\n");
                continue;
            }
            np = 0;
            pnodes[++np] = j0;
            pnodes[++np] = j1;

            tda_out("found j0=%d j1=%d\n",j0,j1);

            /* find additional points until arriving back at i0 */

            while (1) {

                if ((l = sdp_find_idx(ctx, j1,j0)) < 0) {   /* index to edge list */

tda_out("ERROR -4 j0=%d j1=%d\n",j0,j1 );

k = ctx->SD_EPtr[j0];
n = ctx->SD_NA[j0];
for (ll = 0; ll < n; ++ll)  
    tda_out("j0=%d k=%d n=%d elist=%d\n",j0,k,n,ctx->SD_EList[k + ll]);

k = ctx->SD_EPtr[j1];
n = ctx->SD_NA[j1];
for (ll = 0; ll < n; ++ll)  
    tda_out("j1=%d k=%d n=%d elist=%d\n",j1,k,n,ctx->SD_EList[k + ll]);



                    return(-4);
             }
                if ((j2 = sdp_find_nxt(ctx, j1,l,1)) < 0)   /* next node */
                    return(-5);
      
                tda_out("Next point: %d\n",j2);
                       
                if (j2 == i0) {
                    tda_out("found new polygon. np=%d\n",np);

                    /* write polygon to PMF1d. spatial data file format.   
                       count number of records in nrec. count polygon IDs
                       in sdid. */

                    *sdid += 1;
                    *nrec += sdp_write(ctx, *sdid,np,pnodes,xc,yc);
                    r++;

                    /* adjust degrees according to boundary edges. */

                    for (i = 1; i <= np; ++i) {
                        if (i < np) {
                            j0 = pnodes[i];
                            j1 = pnodes[i + 1];
                        }
                        else {
                            j0 = pnodes[np];
                            j1 = pnodes[1];
                        }
                        if (sdcpol_edge(ctx, j0,j1,nb,bnodes,idx) == 1) {
                            sdp_deg(ctx, j0,j1);       
                            sdp_deg(ctx, j1,j0);       
                        }
                    }
                    break;
                }

                /* check whether j2 already in the node list. */

                j = -1;
                for (i = 1; i <= np; ++i) {
                    if (pnodes[i] == j2) {
                        j = i;
                        break;
                    }   
                }
                if (j >= 1) {   /* not a valid polygon */
tda_out("NOT VALID j=%d j2=%d    \n",j,j2   );

                    /* adjust degrees for the first two nodes. */

                    if (sdcpol_edge(ctx, pnodes[1],pnodes[2],nb,bnodes,idx) == 1) {
                        sdp_deg(ctx, pnodes[1],pnodes[2]);       
                        sdp_deg(ctx, pnodes[2],pnodes[1]);       
                    }
                    break;
                }
                /* save new node in pnodes. */

                if (np >= npmax)    /* this might occur if the polygon */
                    return(-6);     /* is not simple */
                pnodes[++np] = j2;
                j0 = j1;
                j1 = j2;
            }
                      
            /* clear the current component from all nodes with degree < 2. */

            nnc = sdp_cclear(ctx, nn,nc,cn);
            if (nnc < 0)
                return(-9);

                   
            tda_out("FINAL degrees nc=%d nnc=%d     \n",nc,nnc  );
            for (i = 1; i <= nn; ++i) 
                tda_out("i=%3d deg=%3d cn=%4d  \n",i,ctx->SD_Deg[i],cn[i]  );
            
            if (nnc < 3)        /* Number of nodes is less than 3 */
                break;
        }
        if (nnc < 3)
            break;

        /* find root node for new reduced component. */

        rn = -1;
        x = y = 0.0;
        for (i = 1; i <= nn; ++i) {
            if (cn[i] != nc || ctx->SD_Deg[i] < 2)
                continue;

            if (rn < 0 || (double)(xc[i]) < x || ((double)(xc[i]) == x && (double)(yc[i]) > y)) {
                x = (double)(xc[i]);
                y = (double)(yc[i]);
                rn = i;
            }
        }
        tda_out("NEW rn = %d\n",rn);

        if (rn < 0) {
            tda_out("WARUM nnc=%d\n",nnc);

            break;
        }
    }
    return(r);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_comp(k,nc,cn)  Create next component with level nc. Use all nodes   */
/*                     with cn[i] = 0. Count number of elements in SDC_N    */

void sdp_comp(TDAContext *ctx, int k,int nc,int *cn)
{  
    register int l,j,n;

    ctx->SDC_N++;
    cn[k] = nc;
    n = ctx->SD_NA[k];
    if (n <= 0)
        return;

    for (l = 0; l < n; ++l) {
        j = ctx->SD_EList[ctx->SD_EPtr[k] + l];            /* edge from k to j */
        if (j >= 1 && cn[j] == 0)  
            sdp_comp(ctx, j,nc,cn);
    }
}

/* -###-------------------------------------------------------------------- */
/*  sdp_deg(i,j)   Remove node j from the followers of i in edge list.      */
/*                 Return 1 if removed, otherwise 0                         */

int sdp_deg(TDAContext *ctx, int i,int j)
{
    register int k,l,n;

    k = ctx->SD_EPtr[i];
    n = ctx->SD_NA[i];
    for (l = 0; l < n; ++l) {
        if (ctx->SD_EList[k + l] == j) {
            ctx->SD_EList[k + l] = -1;
            ctx->SD_Deg[i] -= 1;
            return(1);
        }
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  sdcpol_edge(i,j,n,bnodes,idx)                                           */
/*                                                                          */
/*  Return 1 if i and j are connected by an edge on the boundary,           */  
/*  otherwise return 0.                                                     */

int sdcpol_edge(TDAContext *ctx, int i,int j,int n,int *bnodes,char *idx)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k;

    if (idx[i] == 0 || idx[j] == 0)        /* not on boundary */
        return(0);

    for (k = 1; k <= n; ++k) {
        if (i == bnodes[k]) {
            if (k == 1) {
                if (j == bnodes[k + 1] || j == bnodes[n]) 
                    return(1);
            }
            else if (k == n) {
                if (j == bnodes[k - 1] || j == bnodes[1]) 
                    return(1);
            }
            else if (j == bnodes[k - 1] || j == bnodes[k + 1])
                return(1);
        }
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_clear(nn)                                                           */
/*                                                                          */
/*  Remove all nodes with degree less than 2. Return number of remaining    */  
/*  nodes. Or -1 if error.                                                  */

int sdp_clear(TDAContext *ctx, int nn)
{
    register int i,j,k,ll;
    int n,nd,fin;
                
    nd = fin = 0;
    while (fin == 0) {
        fin = 1;
        nd = 0;
        for (i = 1; i <= nn; ++i) {
            if (ctx->SD_Deg[i] < 1)
                continue;
            else if (ctx->SD_Deg[i] > 1)
                nd++;
            else {
                k = ctx->SD_EPtr[i];
                n = ctx->SD_NA[i];
                j = -1;
                for (ll = 0; ll < n; ++ll) {
                    if ((j = ctx->SD_EList[k + ll]) >= 1)  
                        break;
                }
                if (j < 0)
                    return(-1);

                sdp_deg(ctx, j,i);       
                ctx->SD_EList[k + ll] = -1;
                ctx->SD_Deg[i] -= 1;
                fin = 0;
            }
        }
    }
    return(nd);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_cclear(nn,nc,cn)                                                    */
/*                                                                          */
/*  Remove all nodes with degree less than 2 in the current component       */
/*  containing all nodes i with cn[i] = nc.                                 */
/*  Return number of remaining nodes. Or -1 if error.                       */

int sdp_cclear(TDAContext *ctx, int nn,int nc,int *cn)
{
    register int i,j,k,ll;
    int n,nd,fin;
                
    nd = fin = 0;
    while (fin == 0) {
        fin = 1;
        nd = 0;
        for (i = 1; i <= nn; ++i) {
            if (cn[i] != nc || ctx->SD_Deg[i] < 1)
                continue;
            else if (ctx->SD_Deg[i] > 1)
                nd++;
            else {
                k = ctx->SD_EPtr[i];
                n = ctx->SD_NA[i];
                j = -1;
                for (ll = 0; ll < n; ++ll) {
                    if ((j = ctx->SD_EList[k + ll]) >= 1)  
                        break;
                }
                if (j < 0)
                    return(-1);

                sdp_deg(ctx, j,i);       
                ctx->SD_EList[k + ll] = -1;
                ctx->SD_Deg[i] = 0;
                fin = 0;
            }
        }
    }
    return(nd);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_remove(i)    Remove node i from graph, adjust SD_Deg[].             */

void sdp_remove(TDAContext *ctx, int i)
{
    register int k,l,j,n;

    k = ctx->SD_EPtr[i];
    n = ctx->SD_NA[i];
    for (l = 0; l < n; ++l) {
        if ((j = ctx->SD_EList[k + l]) >= 1) {
            ctx->SD_EList[k + l] = -1;
            ctx->SD_Deg[i] -= 1;
            sdp_deg(ctx, j,i);
        }
    }
}


/* -###-------------------------------------------------------------------- */
/*  sdp_cut(nn,vn,idx)                                                      */
/*                                                                          */  
/*  nn is number of nodes. vn[] and idx[] are work arrays of length nn.     */

void sdp_cut(TDAContext *ctx, int nn,int *vn,char *idx)
{
    register int i,j;
    int k,nid; 

    ctx->SDC_ID = 0;

    for (i = 1; i <= nn; ++i)
        vn[i] = 0;

    while (1) {

        /* find an unprocessed node */

        k = -1;
        for (i = 1; i <= nn; ++i) {
            if (ctx->SD_Deg[i] >= 2 && vn[i] == 0) {
                k = i;                    
                break;
            }
        }
        if (k < 0)      /* no more components */
            break;

        /* Find the next component and, at the same time, find the
           cut nodes. */

        ctx->SDC_N = 0;                  /* number of elements in component */
        ctx->SDC_CNT = 0;                /* number of path from root */
        nid = ctx->SDC_ID;               /* nodes with vn[i] > nid belong to   
                                       current component. */
        for (i = 1; i <= nn; ++i)
            idx[i] = 0;

        sdp_visit(ctx, k,k,vn,idx);

        if (ctx->SDC_CNT < 2)
            idx[k] = 0;

        tda_out("\nNext component nnn_nn=%d nid=%d\n",ctx->SDC_N,nid );
/**
        for (j = 1; j <= nn; ++j)
            tda_out("j=%3d nnn=%4d nnns=%4d  SD_Deg=%3d\n",j,vn[j],idx[j],SD_Deg[j] );
**/ 

        for (j = 1; j <= nn; ++j) {
            if (ctx->SD_Deg[j] != 2)   
                idx[j] = 0;
        }     
        for (j = 1; j <= nn; ++j) {
            if (idx[j]) {

tda_out("REMOVE %d\n",j);
                sdp_remove(ctx, j);
            }
        }
    }
    tda_out("NO MORE\n");
}

/* ------------------------------------------------------------------------ */
/*  sdp_visit(i,ir,vn,idx)    Called by sdp_cut().                          */

int sdp_visit(TDAContext *ctx, int i,int ir,int *vn,char *idx)
{
    register int j,k,ll;
    int n,m,min;

    vn[i] = ++ctx->SDC_ID;
    min = ctx->SDC_ID;
    ctx->SDC_N++; 

    k = ctx->SD_EPtr[i];
    n = ctx->SD_NA[i];
    for (ll = 0; ll < n; ++ll) {
        if ((j = ctx->SD_EList[k + ll]) >= 1) {  /* edge from i to j */
            if (vn[j] == 0) {
                if (i == ir)
                    ctx->SDC_CNT++;

                m = sdp_visit(ctx, j,ir,vn,idx);
                if (m < min)
                    min = m;
                if (m >= vn[i])
                    idx[i] = 1;
            }
            else if (vn[j] < min)
                min = vn[j];
        }
    }
    return(min);
}

/* ------------------------------------------------------------------------ */
/*  sdp_find_idx(i,j)   If node j is not adjacent to node i return -1,      */
/*                      otherwise the index to the edge list SD_EList.      */

int sdp_find_idx(TDAContext *ctx, int i,int j)
{
    register int l,k,n;

    k = ctx->SD_EPtr[i];
    n = ctx->SD_NA[i];

    for (l = 0; l < n; ++l) {
        if (ctx->SD_EList[k + l] == j)    
            return(l);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  sdp_find_nxt(i,l,opt)   i is node number, l is index to its adjacent    */
/*                          nodes in the edge list. If opt = 0 return next  */
/*                          higher node number, otherwise next lower node   */
/*                          number. Or -1 if there is no adjacent node.     */

int sdp_find_nxt(TDAContext *ctx, int i,int l,int opt)
{
    register int ll,k,n,j;

    k = ctx->SD_EPtr[i];
    n = ctx->SD_NA[i];

    if (opt == 0) {
        for (ll = l + 1; ll < n; ++ll) {
            if ((j = ctx->SD_EList[k + ll]) >= 1)  
                return(j);           
        }
        for (ll = 0; ll < l; ++ll) {
            if ((j = ctx->SD_EList[k + ll]) >= 1)  
                return(j);               
        }
    }                 
    else {
        for (ll = l - 1; ll >= 0; --ll) {
            if ((j = ctx->SD_EList[k + ll]) >= 1)  
                return(j);           
        }
        for (ll = n - 1; ll > l; --ll) {
            if ((j = ctx->SD_EList[k + ll]) >= 1)  
                return(j);               
        }
    }
    return(-1);
}                 

/* ------------------------------------------------------------------------ */
/*  sdp_write(id,np,pnodes,xc,yc)                                           */
/*                                                                          */  
/*  Write polygon to PMF1d. Return number of records written.               */

int sdp_write(TDAContext *ctx, int id,int np,int *pnodes,float *xc,float *yc)
{
    register int i,j;
    int nrec = 0;

    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,id);
    fprintf(ctx->PMF1d,"3 ");
    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,np);
    fprintf(ctx->PMF1d,"\n");
    nrec++;         
    for (i = 1; i <= np; ++i) {
        j = pnodes[i];
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)xc[j]);
        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,(double)yc[j]);
        fprintf(ctx->PMF1d,"\n");
        nrec++;         
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  sdp_list(n,xa,ya,xb,yb,ptr,idx,nmax,tol)                                */
/*                                                                          */  
/*  Create the node and edge list for the graph.                            */
/*  [(xa[i],ya[i]),(xb[i],yb[i])] for i = 0,...,n - 1 is the given list     */
/*  of coordinates. All arrays are dimensioned for i = 0,...,nmax - 1.      */
/*                                                                          */
/*  Set: SD_NN = number of nodes                                            */
/*  Set: SD_NE = number of edges                                            */
/*                                                                          */
/*  Return 0 if successful, or negative value for error.                    */
/*  -1 if insufficient memory.                                              */
/*  -2 if nmax is too small.                                                */
/*  -3 if error in the algorithm.                                           */

int sdp_list(TDAContext *ctx, int n,float *xa,float *ya,float *xb,float *yb,int *ptr,char *idx, int nmax,float tol)
{
    register int i,j,k;
    int err,m,nne,nnd,nnl,r,ri,rj,fin,nfree,xfa,yfa,ii,jj;
    float t,xa0,ya0,*xf,*yf;
    double d,xai,yai,xbi,ybi,xaj,yaj,xbj,ybj,xxa,yya,xxb,yyb;

    err = -2;
    xfa = yfa = 0;

    for (i = 0; i < n; ++i) {
        tda_out("i=%3d %f %f %f %f\n",i,(double)(xa[i]),(double)(ya[i]),(double)(xb[i]),(double)(yb[i]));
    }
    newline(ctx);

    /* order all segments in the same x direction */

    for (i = 0; i < n; ++i) {
        if (xa[i] > xb[i]) {
            t = xa[i];
            xa[i] = xb[i];
            xb[i] = t;
            t = ya[i];
            ya[i] = yb[i];
            yb[i] = t;
        }
    }

    /* drop segment with length <= tol */

    nne = 0;
    for (i = 0; i < n; ++i) {
        if (sdp_id(ctx, (double)xa[i],(double)ya[i],(double)xb[i],(double)yb[i],(double)tol))
            continue;

        xa[nne] = xa[i];
        ya[nne] = ya[i];
        xb[nne] = xb[i];
        yb[nne] = yb[i];
        nne++;
    }

    if (nne == 0)
        return(0);
             
    /* crossings. ptr[] will be used to keep track of free places. */

    for (i = 0; i < nmax; ++i)  
        idx[i] = 0;
       
    nfree = fin = 0;
    while (fin == 0) {
        fin = 1;
        nnd = nne;

        for (i = 0; i < nne; ++i) {
            if (idx[i])
                continue;
            xai = (double)xa[i];
            yai = (double)ya[i];
            xbi = (double)xb[i];
            ybi = (double)yb[i];

            if (sdp_id(ctx, (double)xai,(double)yai,(double)xbi,(double)ybi,(double)tol)) {       /* zero length */

            /**
            tda_out("found zero lengh i=%d  will free i=%d nfree=%d   \n",i,i,nfree  );
            **/

                ptr[nfree] = i;
                idx[i] = 1;
                nfree++;
                continue;
            }
            for (j = i + 1; j < nne; ++j) {
                if (idx[j])
                    continue;

                xaj = (double)xa[j];
                yaj = (double)ya[j];
                xbj = (double)xb[j];
                ybj = (double)yb[j];

                if (sdp_id(ctx, (double)xaj,(double)yaj,(double)xbj,(double)ybj,(double)tol)) {       /* zero length */
                    /**
                    tda_out("found zero lengh j=%d  will free j=%d nfree=%d   \n",j,j,nfree  );
                    **/
                    ptr[nfree] = j;
                    idx[j] = 1;
                    nfree++;
                    continue;
                }
                /**
                tda_out("i=%d j=%d xi=%g %g %g %g xj=%g %g %g %g\n",i,j,xai,yai,xbi,ybi,
                        xaj,yaj,xbj,ybj);
                **/

                if ((sdp_id(ctx, (double)xai,(double)yai,(double)xaj,(double)yaj,(double)tol)  && 
                     sdp_id(ctx, xbi,ybi,xbj,ybj,(double)(tol))) ||    
                    (sdp_id(ctx, xai,yai,xbj,ybj,(double)(tol))  && 
                     sdp_id(ctx, xaj,yaj,xbi,ybi,(double)(tol)))) {     
                    /**
                    tda_out("found identical. will free j=%d nfree=%d   \n",j,nfree  );
                    **/
                    ptr[nfree] = j;
                    idx[j] = 1;
                    nfree++;
                    continue;
                }
                r = g_intersect(ctx, xai,yai,xbi,ybi,xaj,yaj,xbj,ybj,&xxa,&yya,&xxb,&yyb);
                /**
                tda_out("rrrr=%d\n",r);
                **/
                if (r == 0)
                    continue;

                if (r == 1) {           /* crossing in single point */

                    if (sdp_id(ctx, (double)xai,(double)yai,(double)xxa,(double)yya,(double)tol) || 
                        sdp_id(ctx, xbi,ybi,xxa,yya,(double)(tol)))       
                        ri = 0;
                    else
                        ri = 1;

                    if (sdp_id(ctx, (double)xaj,(double)yaj,(double)xxa,(double)yya,(double)tol) || 
                        sdp_id(ctx, xbj,ybj,xxa,yya,(double)(tol)))       
                        rj = 0;
                    else
                        rj = 1;
                    /**
                    tda_out("ri=%d rj=%d\n",ri,rj);
                    **/
                    if (ri) {           /* i needs splitting */

                        /**
                        tda_out("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                        **/
                        xb[i] = (float)xxa;
                        yb[i] = (float)yya;

                        if (nfree > 0) {
                            k = ptr[nfree - 1];
                            idx[k] = 0;
                        }
                        else {
                            if (nnd >= nmax - 1)
                                goto SDPListFin;
                            k = nnd;
                        }
                        xa[k] = (float)xxa;
                        ya[k] = (float)yya;
                        xb[k] = (float)xbi;
                        yb[k] = (float)ybi;
                        /**
                        tda_out("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xbi,ybi);
                        **/
                        xbi = xxa;
                        ybi = yya;

                        if (nfree > 0)
                            nfree--;
                        else
                            nnd++;
                        fin = 0;
                    }
                    if (rj) {           /* j needs splitting */
                        /**
                        tda_out("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                        **/
                        xb[j] = (float)xxa;
                        yb[j] = (float)yya;

                        if (nfree > 0) {
                            k = ptr[nfree - 1];
                            idx[k] = 0;
                        }
                        else {
                            if (nnd >= nmax - 1)
                                goto SDPListFin;
                            k = nnd;
                        }
                        xa[k] = (float)xxa;
                        ya[k] = (float)yya;
                        xb[k] = (float)xbj;
                        yb[k] = (float)ybj;
                        /**
                        tda_out("setzt bei k=%d ein: %g %g %g %g\n",k, xxa,yya,xbj,ybj);
                        **/
                        if (nfree > 0)
                            nfree--;
                        else
                            nnd++;
                        fin = 0;
                    }
                }      
                else if (r == 2) {      /* overlap */
            
                    if (xxa > xxb) {
                        d = xxa;
                        xxa = xxb;
                        xxb = d;
                        d = yya;  
                        yya = yyb;
                        yyb = d;
                    }
                    if (sdp_id(ctx, (double)xxa,(double)yya,(double)xxb,(double)yyb,(double)tol) == 0) {
                        /**
                        tda_out("xxa=%g %g xxb=%g %g nicht identisch\n",xxa,yya,xxb,yyb);
                        **/
                        if ((sdp_id(ctx, (double)xxa,(double)yya,(double)xai,(double)yai,(double)tol)  &&
                             sdp_id(ctx, xxb,yyb,xbi,ybi,(double)(tol))) ||
                            (sdp_id(ctx, xxb,yyb,xai,yai,(double)(tol))  &&
                             sdp_id(ctx, xxa,yya,xbi,ybi,(double)(tol)))) {
                            /**
                            tda_out("xxa usw identisch mit i.\n");
                            tda_out("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                            **/
                            xb[j] = (float)xxa;
                            yb[j] = (float)yya;

                            if (nfree > 0) {
                                k = ptr[nfree - 1];
                                idx[k] = 0;
                            }
                            else {
                                if (nnd >= nmax - 1)
                                    goto SDPListFin;
                                k = nnd;
                            }
                            /**
                            tda_out("setzt bei k=%d ein: %g %g %g %g\n",k,xxb,yyb,xbj,ybj);
                            **/
                            xa[k] = (float)xxb;
                            ya[k] = (float)yyb;
                            xb[k] = (float)xbj;
                            yb[k] = (float)ybj;

                            if (nfree > 0)
                                nfree--;
                            else
                                nnd++;
                            fin = 0;
                        }
                        else if ((sdp_id(ctx, xxa,yya,xaj,yaj,(double)(tol))  &&
                                  sdp_id(ctx, xxb,yyb,xbj,ybj,(double)(tol))) ||
                                 (sdp_id(ctx, xxb,yyb,xaj,yaj,(double)(tol))  &&
                                  sdp_id(ctx, xxa,yya,xbj,ybj,(double)(tol)))) {

                            /**
                            tda_out("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                            **/

                            xb[i] = (float)xxa;
                            yb[i] = (float)yya;

                            if (nfree > 0) {
                                k = ptr[nfree - 1];
                                idx[k] = 0;
                            }
                            else {
                                if (nnd >= nmax - 1)
                                    goto SDPListFin;
                                k = nnd;
                            }
                            /**
                            tda_out("setzt bei k=%d ein: %g %g %g %g\n",k,xxb,yyb,xbi,ybi);
                            **/
                            xa[k] = (float)xxb;
                            ya[k] = (float)yyb;
                            xb[k] = (float)xbi;
                            yb[k] = (float)ybi;

                            if (nfree > 0)
                                nfree--;
                            else
                                nnd++;
                            fin = 0;
                        }
                        else if (sdp_id(ctx, xxa,0.0,xxb,0.0,(double)(tol))) {
                            if (yai < yya) {
                                /**
                                tda_out("yai kleiner als yya\n");
                                tda_out("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                                tda_out("setzt bei xa[j] ein: %g %g\n",xxb,yyb);
                                **/

                                xb[i] = (float)xxa;
                                yb[i] = (float)yya;
                                xa[j] = (float)xxb;
                                ya[j] = (float)yyb;
                                xbi = xxa;
                                ybi = yya;
                            }
                            else {
                                /*
                                tda_out("yai nicht  kleiner als yya\n");
                                tda_out("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                                tda_out("setzt bei xa[i] ein: %g %g\n",xxb,yyb);
                                **/
                                xb[j] = (float)xxa;
                                yb[j] = (float)yya;
                                xa[i] = (float)xxb;
                                ya[i] = (float)yyb;
                                xai = xxb;
                                yai = yyb;
                            }
                            if (nfree > 0) {
                                k = ptr[nfree - 1];
                                idx[k] = 0;
                            }
                            else {
                                if (nnd >= nmax - 1)
                                    goto SDPListFin;
                                k = nnd;
                            }
                            /*
                            tda_out("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xxb,yyb);
                            */
                            xa[k] = (float)xxa;
                            ya[k] = (float)yya;
                            xb[k] = (float)xxb;
                            yb[k] = (float)yyb;

                            if (nfree > 0)
                                nfree--;
                            else
                                nnd++;
                            fin = 0;
                        }
                        else {
                            if (xai < xxa) {
                                /*
                                tda_out("xai kleiner als xxa\n");
                                tda_out("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                                tda_out("setzt bei xa[j] ein: %g %g\n",xxb,yyb);
                                **/

                                xb[i] = (float)xxa;
                                yb[i] = (float)yya;
                                xa[j] = (float)xxb;
                                ya[j] = (float)yyb;
                                xbi = xxa;
                                ybi = yya;
                            }
                            else {
                                /**
                                tda_out("xai nicht  kleiner als xxa\n");
                                tda_out("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                                tda_out("setzt bei xa[i] ein: %g %g\n",xxb,yyb);
                                **/
                                xb[j] = (float)xxa;
                                yb[j] = (float)yya;
                                xa[i] = (float)xxb;
                                ya[i] = (float)yyb;
                                xai = xxb;
                                yai = yyb;
                            }
                            if (nfree > 0) {
                                k = ptr[nfree - 1];
                                idx[k] = 0;
                            }
                            else {
                                if (nnd >= nmax - 1)
                                    goto SDPListFin;
                                k = nnd;
                            }
                            /*
                            tda_out("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xxb,yyb);
                            */
                            xa[k] = (float)xxa;
                            ya[k] = (float)yya;
                            xb[k] = (float)xxb;
                            yb[k] = (float)yyb;

                            if (nfree > 0)
                                nfree--;
                            else
                                nnd++;
                            fin = 0;
                        }
                    }
                }
            }
        }
        nne = nnd;
        /**
        tda_out("After nne=%d nfree=%d: ",nne,nfree );
        for (k = 0; k < nfree; ++k)
            tda_out("%d ",ptr[k]);
        newline(ctx);

        for (i = 0; i < nne; ++i) {
            tda_out("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]);
        }
        newline(ctx);
        **/
        if (nfree > 0) {        /* reorganize */
            k = 0;
            for (i = 0; i < nne; ++i) {
                if (idx[i] == 0) {
                    xa[k] = xa[i];
                    ya[k] = ya[i];
                    xb[k] = xb[i];
                    yb[k] = yb[i];
                    k++;
                }
                else  
                    idx[i] = 0;
            }
            nfree = 0;
            nne = k;
            /**
            tda_out("After reorg nne=%d nfree=%d: ",nne,nfree );
            for (i = 0; i < nne; ++i) {
                tda_out("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]);
            }
            newline(ctx);
            **/
        }
    }

    tda_out("After splitting nne=%d nnd=%d nfree=%d: ",nne,nnd,nfree );
    for (k = 0; k < nfree; ++k)
        tda_out("%d ",ptr[k]);
    newline(ctx);

    for (i = 0; i < nne; ++i) {
        tda_out("i=%3d %f %f %f %f idx=%d\n",i,(double)(xa[i]),(double)(ya[i]),(double)(xb[i]),(double)(yb[i]),idx[i]  );
    }
    newline(ctx);

    if (nfree > 0) {        /* reorganize */
        k = 0;
        for (i = 0; i < nne; ++i) {
            if (idx[i] == 0) {
                xa[k] = xa[i];
                ya[k] = ya[i];
                xb[k] = xb[i];
                yb[k] = yb[i];
                k++;
            }
            else  
                idx[i] = 0;
        }
        nfree = 0;
        nne = k;

        tda_out("After reorg nne=%d nfree=%d: ",nne,nfree );
        for (i = 0; i < nne; ++i) {
            tda_out("i=%3d %f %f %f %f idx=%d\n",i,(double)(xa[i]),(double)(ya[i]),(double)(xb[i]),(double)(yb[i]),idx[i]);
        }
        newline(ctx);
    }

    err = -1;
/***
    if (sortdp4f(ctx, nne,xa,ya,xb,yb,ptr))
        goto SDPListFin;

    for (i = 0; i < nne; ++i) {
        k = ptr[i];
        tda_out("i=%3d %f %f %f %f\n",i,xa[k],ya[k],xb[k],yb[k]);
    }
    newline(ctx);
**/        
    /* create list of unique coordinates */

    m = 2 * nne;
    if (!(ctx->SD_EPtr = (int *)calloc((size_t)(m),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, m,sizeof(int));
    ctx->SD_EPtrA = m;     

    if (!(xf = (float *)calloc((size_t)(m),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, m,sizeof(float));
    xfa = m;             

    if (!(yf = (float *)calloc((size_t)(m),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, m,sizeof(float));
    yfa = m;             

    k = 0;
    for (i = 0; i < nne; ++i) {
        xf[k] = xa[i];
        yf[k] = ya[i];
        k++;
        xf[k] = xb[i];
        yf[k] = yb[i];
        k++;
    }
    if (sortdp2f(ctx, m,xf,yf,ctx->SD_EPtr))
        goto SDPListFin;

    nnd = 0;                /* number of unique coordinates (nodes) */
    k = ctx->SD_EPtr[0];
    xa0 = (float)((double)(xf[k]) - 1000.0);
    ya0 = (float)((double)(yf[k]) - 1000.0);
    for (i = 0; i < m; ++i) {
        k = ctx->SD_EPtr[i];
        if (sdp_id(ctx, (double)xf[k],(double)yf[k],(double)xa0,(double)ya0,(double)tol) == 0) {
            nnd++;
            xa0 = xf[k];
            ya0 = yf[k];
        }
    }

    tda_out("nnd=%d\n",nnd);

    /* save unique list of coordinates in SD_NLX and SD_NLY */

    if (!(ctx->SD_NLX = (float *)calloc((size_t)(nnd + 2),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnd + 2,sizeof(float));
    ctx->SD_NLXA = nnd + 2;       

    if (!(ctx->SD_NLY = (float *)calloc((size_t)(nnd + 2),sizeof(float)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnd + 2,sizeof(float));
    ctx->SD_NLYA = nnd + 2;       

    j = 0;
    k = ctx->SD_EPtr[0];
    xa0 = (float)((double)(xf[k]) - 1000.0);
    ya0 = (float)((double)(yf[k]) - 1000.0);
    for (i = 0; i < m; ++i) {
        k = ctx->SD_EPtr[i];
        if (sdp_id(ctx, (double)xf[k],(double)yf[k],(double)xa0,(double)ya0,(double)tol) == 0) {
            j++;
            ctx->SD_NLX[j] = xa0 = xf[k];
            ctx->SD_NLY[j] = ya0 = yf[k];
        }
    }
    if (j != nnd) {
        err = -3;
        goto SDPListFin;
    }

    tda_out("LIST of node coordinates nnd=%d\n",nnd);

    for (i = 1; i <= nnd; ++i) {
        tda_out("i=%3d x=%f y=%f\n",i,(double)(ctx->SD_NLX[i]),(double)(ctx->SD_NLY[i]));
    }

    free((char *)xf);
    memrq(ctx, -xfa,sizeof(float));
    xfa = 0;
    free((char *)yf);
    memrq(ctx, -yfa,sizeof(float));
    yfa = 0;
    free((char *)ctx->SD_EPtr);
    memrq(ctx, -ctx->SD_EPtrA,sizeof(int));
    ctx->SD_EPtrA = 0;

    /* create forward-linked edge list in SD_EList, SD_EPtr, SD_NA */

    nnl = 2 * nne;
    if (!(ctx->SD_EList = (int *)calloc((size_t)(nnl + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnl + 1,sizeof(int));
    ctx->SD_EListA = nnl + 1;     

    if (!(ctx->SD_NA = (int *)calloc((size_t)(nnl + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnl + 1,sizeof(int));
    ctx->SD_NAA = nnl + 1;     

    j = 0;
    for (i = 0; i < nne; ++i) {
        tda_out("i=%3d %f %f %f %f\n",i,(double)(xa[i]),(double)(ya[i]),(double)(xb[i]),(double)(yb[i]));
        
        ii = sdp_list_find(ctx, xa[i],ya[i],nnd,ctx->SD_NLX,ctx->SD_NLY,tol);
        jj = sdp_list_find(ctx, xb[i],yb[i],nnd,ctx->SD_NLX,ctx->SD_NLY,tol);
        if (ii < 1 || jj < 1) {
            err = -3;
            goto SDPListFin;
        }
        j++;         
        ctx->SD_NA[j] = ii;
        ctx->SD_EList[j] = jj;
        j++;         
        ctx->SD_NA[j] = jj;
        ctx->SD_EList[j] = ii;
    }
    if (sorti2(ctx, nnl,ctx->SD_NA + 1,ctx->SD_EList + 1))
        goto SDPListFin;

    if (!(ctx->SD_EPtr = (int *)calloc((size_t)(nnd + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnd + 1,sizeof(int));
    ctx->SD_EPtrA = nnd + 1;     

    k = -1;
    j = 0;
    for (i = 1; i <= nnl; ++i) {
        if (ctx->SD_NA[i] != k) {
            ctx->SD_EPtr[++j] = i;
            k = ctx->SD_NA[i];
        }
    }
    if (j != nnd) {
        err = -3;
        goto SDPListFin;
    }
    free((char *)ctx->SD_NA);
    memrq(ctx, -ctx->SD_NAA,sizeof(int));
    ctx->SD_NAA = 0;

    if (!(ctx->SD_NA = (int *)calloc((size_t)(nnd + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnd + 1,sizeof(int));
    ctx->SD_NAA = nnd + 1;     

    for (i = 1; i < nnd; ++i)
        ctx->SD_NA[i] = ctx->SD_EPtr[i + 1] - ctx->SD_EPtr[i];
    ctx->SD_NA[nnd] = nnl - ctx->SD_EPtr[nnd] + 1;

    /* create degrees in SD_Deg */

    if (!(ctx->SD_Deg = (int *)calloc((size_t)(nnd + 1),sizeof(int)))) { 
        p_err(ctx, -2,1);
        goto SDPListFin;
    }   
    memrq(ctx, nnd + 1,sizeof(int));
    ctx->SD_DegA = nnd + 1;     

    for (i = 1; i <= nnd; ++i)
        ctx->SD_Deg[i] = ctx->SD_NA[i];


    ctx->SD_NN = nnd;
    ctx->SD_NE = nnl;
    err = 0;

SDPListFin:
    if (err == -2)  
        printf1(ctx, "Error: exceeded the maximum number of segments (%d).\n",nmax);
    else if (err == -3) {
        printf1(ctx, "Error while creating the edge list.\n");
        printf1(ctx, "Might be caused by rounding errors.\n");
        printf1(ctx, "Cannot continue.\n");
    }
    if (xfa > 0) {
        free((char *)xf);
        memrq(ctx, -xfa,sizeof(float));
    }
    if (yfa > 0) {
        free((char *)yf);
        memrq(ctx, -yfa,sizeof(float));
    }
    if (err)
        sdp_list_free(ctx);

    return(err);
}


/* ------------------------------------------------------------------------ */
/*  sdp_list_find(xa0,ya0,n,xa,ya,tol)                                      */
/*                                                                          */
/*  Return index of (xa0,ya0) in list xa[i],ya[i] (i = 1,n) or -1 if        */
/*  not found.                                                              */
 
int sdp_list_find(TDAContext *ctx, float xa0,float ya0,int n,float *xa,float *ya,float tol)
{
    register int i,l,r,k;

    l = 1;
    r = n;
    while (l <= r) {
        k = (l + r) / 2;
        if ((float)fabs((double)(xa0 - xa[k])) <= tol) {
            for (i = k; i <= n; ++i) {
                if ((float)fabs((double)(xa[i] - xa[k])) > tol)
                    break;
                if (sdp_id(ctx, (double)xa0,(double)ya0,(double)xa[i],(double)ya[i],(double)tol))        
                    return(i);
            }           
            for (i = k - 1; i >= 1; --i) {
                if ((float)fabs((double)(xa[i] - xa[k])) > tol)
                    break;
                if (sdp_id(ctx, (double)xa0,(double)ya0,(double)xa[i],(double)ya[i],(double)tol))        
                    return(i);
            }           
        }
        if (xa0 < xa[k])
            r = k - 1;
        else
            l = k + 1;
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  sdp_list_free. Free previously allocated memory for SD data structure.  */
 
void sdp_list_free(TDAContext *ctx)
{
    if (ctx->SD_NLXA > 0) {
        free((char *)ctx->SD_NLX);
        memrq(ctx, -ctx->SD_NLXA,sizeof(float));
        ctx->SD_NLXA = 0;
    }
    if (ctx->SD_NLYA > 0) {
        free((char *)ctx->SD_NLY);
        memrq(ctx, -ctx->SD_NLYA,sizeof(float));
        ctx->SD_NLYA = 0;
    }
    if (ctx->SD_EPtrA > 0) {
        free((char *)ctx->SD_EPtr);
        memrq(ctx, -ctx->SD_EPtrA,sizeof(int));
        ctx->SD_EPtrA = 0;
    }
    if (ctx->SD_NAA > 0) {
        free((char *)ctx->SD_NA);
        memrq(ctx, -ctx->SD_NAA,sizeof(int));
        ctx->SD_NAA = 0;
    }
    if (ctx->SD_DegA > 0) {
        free((char *)ctx->SD_Deg);
        memrq(ctx, -ctx->SD_DegA,sizeof(int));
        ctx->SD_DegA = 0;
    }
    if (ctx->SD_EListA > 0) {
        free((char *)ctx->SD_EList);
        memrq(ctx, -ctx->SD_EListA,sizeof(int));
        ctx->SD_EListA = 0;
    }
}



