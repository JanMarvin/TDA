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

/*  functions in t_sdr.c */

int sdnl(void);
int sdppol(void);
int sdlpol(void);
int sdipol(void);

int sdcpol(void);
int sdp_id(double xa,double ya,double xb,double yb,double tol);
void sdp_sort(int i,int k,int n,int *nptr,float *x,float *y,double *a);
int sdp_fpol(int nn,int nc,int *cn,int rn,float *xc,float *yc,char *idx,
    int *bnodes,int nbmax,int *pnodes,int npmax,int *sdid,int *nrec);

void sdp_comp(int k,int nc,int *cn);
int sdp_deg(int i,int j);       
int sdcpol_edge(int i,int j,int n,int *bnodes,char *idx);
int sdp_clear(int nn);
int sdp_cclear(int nn,int nc,int *cn);
void sdp_remove(int i);
void sdp_cut(int nn,int *vn,char *idx);
int sdp_visit(int i,int ir,int *vn,char *idx);
int sdp_find_idx(int i,int j);
int sdp_find_nxt(int i,int l,int opt);
int sdp_write(int id,int np,int *pnodes,float *xc,float *yc);
int sdp_list(int n,float *xa,float *ya,float *xb,float *yb,int *ptr,char *idx,
    int nmax,float tol);
int sdp_list_find(float xa0,float ya0,int n,float *xa,float *ya,float tol);
void sdp_list_free(void);

/* ------------------------------------------------------------------------ */
/*  Global variables                                                        */

int SD_NN = 0;              /* number of nodes                              */
int SD_NE = 0;              /* number of edges                              */
float *SD_NLX;              /* node list x coordinate                       */
float *SD_NLY;              /* node list y coordinate                       */
int *SD_EList;              /* forward-linked edge list                     */
int *SD_EPtr;               /* pointer to SD_EList                          */
int *SD_NA;                 /* number of adjacent nodes                     */
int *SD_Deg;                /* degrees of the nodes                         */
int SD_NLXA = 0;            /* used for allocation                          */
int SD_NLYA = 0;            /* used for allocation                          */
int SD_EListA = 0;          /* used for allocation                          */
int SD_EPtrA = 0;           /* used for allocation                          */
int SD_NAA = 0;             /* used for allocation                          */
int SD_DegA = 0;            /* used for allocation                          */
int SDC_ID = 0;
int SDC_CNT = 0;
int SDC_N = 0;

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

int sdnl(void)
{
    register int i,j,k,l;
    int err,r,n,nl,ns,nrec;                 
    double len;

    err = -1;
    printf1("Creating a network from lines. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto SDNLFin;

    if (PMOPT < 1 || PMOPT > 2)
        PMOPT = 1;
    if (PMFmtF == 0)
        pmfmt(10,4);

    nl = ns = 0;
    for (i = 0; i < NOC; ++i) {
        if ((int)get_data(SDVarSDTyp,i) != 2)
            continue;

        n = (int)get_data(SDVarSDN,i);             
        nl += 1;
        ns += n - 1;
    }
    newline();
    printf1("Number of lines: %d\n",nl);
    printf1("Number of segments: %d\n",ns);
    PMTol = dmax(PMTol,EPSI1);
    PMMax = imax(2 * ns,PMMax);
    printf1("Tolerance: %g\n",PMTol);
    printf1("Maximal number of segments: %d\n",PMMax);

    if (alloc_acxf(PMMax + 2))
        goto SDNLFin;     
    if (alloc_acyf(PMMax + 2))
        goto SDNLFin;     
    if (alloc_acuf(PMMax + 2))
        goto SDNLFin;     
    if (alloc_acvf(PMMax + 2))
        goto SDNLFin;     
    if (alloc_ack(PMMax + 2))
        goto SDNLFin;     
    if (alloc_acd(PMMax + 2))
        goto SDNLFin;     

    k = 0;
    for (i = 0; i < NOC; ++i) {
        if ((int)get_data(SDVarSDTyp,i) != 2)
            continue;

        if ((n = sd_getdata(i,0,0,1)) < 1)
            goto SDNLFin;

        for (j = 1; j < n; ++j) {
            AcXF[k] = (float)SDVarX[j - 1];
            AcYF[k] = (float)SDVarY[j - 1];
            AcUF[k] = (float)SDVarX[j];
            AcVF[k] = (float)SDVarY[j];
            k++;
        }
    }
    r = sdp_list(ns,AcXF,AcYF,AcUF,AcVF,AcK,AcD,PMMax,(float)PMTol);         
    if (r != 0)  
        goto SDNLFin;

    printf1("\nNumber of nodes: %d\n",SD_NN);
    printf1("Number of edges: %d\n",SD_NE / 2);

    nrec = 0;
    for (i = 1; i <= SD_NN; ++i) {
        k = SD_EPtr[i];
        n = SD_NA[i];
        for (l = 0; l < n; ++l) {
            j = SD_EList[k + l];
            if (j > i) {
                len = g_len((double)SD_NLX[i],(double)SD_NLY[i],(double)SD_NLX[j],(double)SD_NLY[j]);
                fprintf(PMFd,PMNFmtS,i);        
                fprintf(PMFd,PMNFmtS,j);        
                fprintf(PMFd,PMFmtS,len);        
                if (PMOPT == 2) {
                    fprintf(PMFd,PMFmtS,(double)SD_NLX[i]);        
                    fprintf(PMFd,PMFmtS,(double)SD_NLY[i]);        
                    fprintf(PMFd,PMFmtS,(double)SD_NLX[j]);        
                    fprintf(PMFd,PMFmtS,(double)SD_NLY[j]);        
                }
                fprintf(PMFd,"\n");        
                nrec++;
            }
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDNLFin:
    sdp_list_free();
    p_clean();
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

int sdppol(void)
{
    register int i,j,k,l;
    int err,r,n,np,npol,nt,nrec,typ;             
    float xmin,xmax,ymin,ymax;
    double x,y;

    err = -1;
    nrec = 0;
    printf1("Find points in polygons. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto SDPPFin;

    if (PMOPT < 1 || PMOPT > 3)
        PMOPT = 1;
    if (PMFmtF == 0)
        pmfmt(10,4);
    if (PMN < 1)
        PMN = 1;

    np = npol = 0;
    for (i = 0; i < NOC; ++i) {
        typ = (int)get_data(SDVarSDTyp,i);           
        if (typ == 1)
            np++;
        else if (typ == 3)
            npol++;
    }
    newline();
    printf1("Number of points: %d\n",np);
    printf1("Number of polygons: %d\n",npol);
    if (np < 1 || npol < 1) {
        err = 0;
        goto SDPPFin;
    }   
    if (PMN > npol)
        PMN = npol;

    /* save points in AcX, AcY, id in AcI; boundaries for polygons in 
       AcXF,AxYF,AcUF,AcVF, index in AcJ. AcK is used to save polygon ids 
       for the points. */

    if (alloc_acx(np))
        goto SDPPFin;     
    if (alloc_acy(np))
        goto SDPPFin;     
    if (alloc_aci(np))
        goto SDPPFin;     

    if (alloc_acxf(npol))
        goto SDPPFin;     
    if (alloc_acyf(npol))
        goto SDPPFin;     
    if (alloc_acuf(npol))
        goto SDPPFin;     
    if (alloc_acvf(npol))
        goto SDPPFin;     
    if (alloc_acj(npol))
        goto SDPPFin;     

    if (alloc_ack(npol))
        goto SDPPFin;     

    i = j = 0;
    for (k = 0; k < NOC; ++k) {
        typ = (int)get_data(SDVarSDTyp,k);           
        if (typ != 1 && typ != 3)  
            continue;

        if ((n = sd_getdata(k,0,0,1)) < 1)
            goto SDPPFin;

        if (typ == 1) {
            AcI[i] = (int)get_data(SDVarSDID,k);           
            AcX[i] = SDVarX[0];
            AcY[i] = SDVarY[0];
            i++;
        }
        else {
            AcJ[j] = k;                                    
            xmin = xmax = (float)SDVarX[0];
            ymin = ymax = (float)SDVarY[0];
            for (l = 1; l < n; ++l) {
                xmin = (float)dmin((double)xmin,SDVarX[l]);
                xmax = (float)dmax((double)xmax,SDVarX[l]);
                ymin = (float)dmin((double)ymin,SDVarY[l]);
                ymax = (float)dmax((double)ymax,SDVarY[l]);
            }   
            AcXF[j] = xmin - EPSI1;
            AcUF[j] = xmax + EPSI1;
            AcYF[j] = ymin - EPSI1;
            AcVF[j] = ymax + EPSI1;
            j++;
        }
    }
    for (i = 0; i < np; ++i) {

        nt = 0;                     /* number of insides */
        x = AcX[i];
        y = AcY[i];

        for (j = 0; j < npol; ++j) {
            xmin = AcXF[j];
            xmax = AcUF[j];
            ymin = AcYF[j];
            ymax = AcVF[j];

            if ((float)x >= xmin && (float)x <= xmax &&
                (float)y >= ymin && (float)y <= ymax) {  

                if ((n = sd_getdata(AcJ[j],0,1,1)) < 1)
                    goto SDPPFin;
                              
                if (g_inpoly(n - 1,SDVarX,SDVarY,x,y)) {

                    k = (int)get_data(SDVarSDID,AcJ[j]);           

                    if (PMOPT == 3) {
                        fprintf(PMFd,PMNFmtS,AcI[i]);        
                        fprintf(PMFd,PMNFmtS,k);        
                        fprintf(PMFd,PMFmtS,AcX[i]);        
                        fprintf(PMFd,PMFmtS,AcY[i]);        
                        fprintf(PMFd,"\n");        
                        nrec++;
                    }
                    else 
                        AcK[nt++] = k;                                         
                }
            }
        }
        if (PMOPT == 3)
            continue;

        fprintf(PMFd,PMNFmtS,AcI[i]);        
        if (PMOPT == 1) {
            fprintf(PMFd,PMFmtS,AcX[i]);        
            fprintf(PMFd,PMFmtS,AcY[i]);        
        }
        else {
            fprintf(PMFd,PMNFmtS,1);        
            fprintf(PMFd,PMNFmtS,1);        
        }                                
        fprintf(PMFd,PMNFmtS,nt);        
        n = imin(nt,PMN);
        for (j = 0; j < n; ++j)
            fprintf(PMFd,PMNFmtS,AcK[j]);        
        for (; j < PMN; ++j)
            fprintf(PMFd,PMNFmtS,-1);        
        fprintf(PMFd,"\n");        
        nrec++;
        if (PMOPT == 2) {
            fprintf(PMFd,PMFmtS,AcX[i]);        
            fprintf(PMFd,PMFmtS,AcY[i]);        
            fprintf(PMFd,"\n");        
            nrec++;
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDPPFin:
    p_clean();
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

int sdlpol(void)
{
    register int i,j,k,l;
    int err,r,n,nl,npol,npl,nt,nrec,typ,nlmax,idl,idp;
    float xmin,xmax,ymin,ymax;
    double x,y;

    err = -1;
    nrec = 0;
    printf1("Find lines in polygons. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto SDLPFin;

    if (PMOPT < 1 || PMOPT > 3)
        PMOPT = 1;
    if (PMFmtF == 0)
        pmfmt(10,4);

    nlmax = nl = npol = 0;
    for (i = 0; i < NOC; ++i) {
        typ = (int)get_data(SDVarSDTyp,i);           
        if (typ == 2) {
            nl++;
            nlmax = imax(nlmax,(int)get_data(SDVarSDN,i));
        }
        else if (typ == 3)
            npol++;
    }
    newline();
    printf1("Number of lines: %d\n",nl);
    printf1("Number of polygons: %d\n",npol);
    if (nl < 1 || npol < 1) {
        err = 0;
        goto SDLPFin;
    }   

    /* save boundaries for polygons in AcXF,AxYF,AcUF,AcVF, index in AcJ.
       AcX, AcY is used to save the points of a line. AcI is used for the 
       indices of the lines. */

    if (alloc_acx(nlmax + 1))
        goto SDLPFin;     
    if (alloc_acy(nlmax + 1))
        goto SDLPFin;     
    if (alloc_aci(nl))
        goto SDLPFin;     

    if (alloc_acxf(npol))
        goto SDLPFin;     
    if (alloc_acyf(npol))
        goto SDLPFin;     
    if (alloc_acuf(npol))
        goto SDLPFin;     
    if (alloc_acvf(npol))
        goto SDLPFin;     
    if (alloc_acj(npol))
        goto SDLPFin;     

    i = j = 0;
    for (k = 0; k < NOC; ++k) {
        typ = (int)get_data(SDVarSDTyp,k);           
        if (typ == 2) 
            AcI[i++] = k;
        if (typ != 3)  
            continue;

        if ((n = sd_getdata(k,0,0,1)) < 1)
            goto SDLPFin;

        AcJ[j] = k;                                    
        xmin = xmax = (float)SDVarX[0];
        ymin = ymax = (float)SDVarY[0];
        for (l = 1; l < n; ++l) {
            xmin = (float)dmin((double)xmin,SDVarX[l]);
            xmax = (float)dmax((double)xmax,SDVarX[l]);
            ymin = (float)dmin((double)ymin,SDVarY[l]);
            ymax = (float)dmax((double)ymax,SDVarY[l]);
        }   
        AcXF[j] = xmin - EPSI1;
        AcUF[j] = xmax + EPSI1;
        AcYF[j] = ymin - EPSI1;
        AcVF[j] = ymax + EPSI1;
        j++;
    }
    for (i = 0; i < nl; ++i) {

        idl = (int)get_data(SDVarSDID,AcI[i]);           

        if ((npl = sd_getdata(AcI[i],0,0,1)) < 1)
            goto SDLPFin;

        for (j = 0; j < npl; ++j) {
            AcX[j] = SDVarX[j];
            AcY[j] = SDVarY[j];
        }
        nt = 0;                     /* number of insides */
        r = 0;
        for (k = 0; k < npol; ++k) {
            xmin = AcXF[k];
            xmax = AcUF[k];
            ymin = AcYF[k];
            ymax = AcVF[k];

            r = 1;
            for (j = 0; j < npl; ++j) {
                if ((float)AcX[j] < xmin || (float)AcX[j] > xmax ||
                    (float)AcY[j] < ymin || (float)AcY[j] > ymax) {   
                    r = 0;    
                    break;
                }
            }
            if (r == 0)
                continue;

            if ((n = sd_getdata(AcJ[k],0,1,1)) < 1)
                goto SDLPFin;
                              
            for (j = 1; j < npl; ++j) {
                if (g_l_inpoly(n - 1,SDVarX,SDVarY,AcX[j-1],AcY[j-1],AcX[j],AcY[j]) == 0) {
                    r = 0;
                    break;
                }
            }
            if (r == 0)
                continue;

            idp = (int)get_data(SDVarSDID,AcJ[k]);           

            if (PMOPT == 3) {
                fprintf(PMFd,PMNFmtS,idl);        
                fprintf(PMFd,PMNFmtS,idp);        
                fprintf(PMFd,"\n");        
                nrec++;
            }
            break;
        }
        if (PMOPT == 3)
            continue;

        fprintf(PMFd,PMNFmtS,idl);        
        if (PMOPT == 2) {
            fprintf(PMFd,PMNFmtS,2);        
            fprintf(PMFd,PMNFmtS,npl);        
        }                                
        if (r == 0)
            idp = -1;


        n = imin(nt,PMN);
        for (j = 0; j < n; ++j)
            fprintf(PMFd,PMNFmtS,AcK[j]);        
        for (; j < PMN; ++j)
            fprintf(PMFd,PMNFmtS,-1);        
        fprintf(PMFd,"\n");        
        nrec++;
        if (PMOPT == 2) {
            for (j = 0; j < npl; ++j) {
                fprintf(PMFd,PMFmtS,AcX[j]);        
                fprintf(PMFd,PMFmtS,AcY[j]);        
                fprintf(PMFd,"\n");        
                nrec++;
            }
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDLPFin:
    p_clean();
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

int sdipol(void)
{
    register int i,j,k,l;
    int err,r,n,npol,nrec,npmax,idi,idj,npi,npj,ii,jj;
    float xmin,xmax,ymin,ymax;
    double x,y;

    err = -1;
    nrec = 0;
    printf1("Hierarchical inclusion of polygons. Current memory: %d bytes.\n",MemReq);
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (parm(CmdBuf + 6,1,1))       /* get parameters */
        goto SDIPFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    npmax = npol = 0;
    for (i = 0; i < NOC; ++i) {
        if ((int)get_data(SDVarSDTyp,i) != 3)      
            continue;
        npol++;
        npmax = imax(npmax,(int)get_data(SDVarSDN,i));
    }
    printf1("\nNumber of polygons: %d\n",npol);
    if (npol < 1) {
        err = 0;
        goto SDIPFin;
    }   

    /* save boundaries for polygons in AcXF,AxYF,AcUF,AcVF, index in AcN. */

    if (alloc_acxf(npol))
        goto SDIPFin;     
    if (alloc_acyf(npol))
        goto SDIPFin;     
    if (alloc_acuf(npol))
        goto SDIPFin;     
    if (alloc_acvf(npol))
        goto SDIPFin;     
    if (alloc_acn(npol))
        goto SDIPFin;     
    if (alloc_acd(NOC + 1))
        goto SDIPFin;     

    j = 0;
    for (i = 0; i < NOC; ++i) {
        if ((int)get_data(SDVarSDTyp,i) != 3)      
            continue;

        if ((n = sd_getdata(i,0,0,1)) < 1)
            goto SDIPFin;

        AcN[j] = i;                                    
        xmin = xmax = (float)SDVarX[0];
        ymin = ymax = (float)SDVarY[0];
        for (l = 1; l < n; ++l) {
            xmin = (float)dmin((double)xmin,SDVarX[l]);
            xmax = (float)dmax((double)xmax,SDVarX[l]);
            ymin = (float)dmin((double)ymin,SDVarY[l]);
            ymax = (float)dmax((double)ymax,SDVarY[l]);
        }   
        AcXF[j] = xmin;
        AcUF[j] = xmax;
        AcYF[j] = ymin;
        AcVF[j] = ymax;
        j++;
    }

    if (alloc_acx(npmax + 2))
        goto SDIPFin;     
    if (alloc_acy(npmax + 2))
        goto SDIPFin;     

    for (i = 0; i < npol; ++i) {
        ii = AcN[i];
        idi = (int)get_data(SDVarSDID,ii);           

        if ((npi = sd_getdata(ii,0,1,1)) < 1)
            goto SDIPFin;

        for (l = 0; l < npi; ++l) {
            AcX[l] = SDVarX[l];
            AcY[l] = SDVarY[l];
        }
        xmin = AcXF[i];
        xmax = AcUF[i];
        ymin = AcYF[i];
        ymax = AcVF[i];

        for (j = i + 1; j < npol; ++j) {
            if (xmax < AcXF[j] || xmin > AcUF[j] ||
                ymax < AcYF[j] || ymin > AcVF[j])    
                continue;

            jj = AcN[j];
            if ((npj = sd_getdata(jj,0,1,1)) < 1)
                goto SDIPFin;

            r = g_p_inpoly(npi,AcX,AcY,npj,SDVarX,SDVarY); 
            if (r == 0) {
                r = g_p_inpoly(npj,SDVarX,SDVarY,npi,AcX,AcY);
                if (r == 1)
                    r = 2;
            }
            if (r == 0)
                continue;

            idj = (int)get_data(SDVarSDID,jj);           

            if (r == 1) {
                fprintf(PMFd,PMNFmtS,idi);        
                fprintf(PMFd,PMNFmtS,idj);        
            }
            else {
                fprintf(PMFd,PMNFmtS,idj);        
                fprintf(PMFd,PMNFmtS,idi);        
            }
            fprintf(PMFd,PMNFmtS,r);        
            fprintf(PMFd,"\n");        
            nrec++;

            AcD[ii] = 1;
            AcD[jj] = 1;
        }
    }

    /* add isolated nodes */

    for (i = 0; i < npol; ++i) {
        ii = AcN[i];
        if (AcD[ii] == 0) {
            idi = (int)get_data(SDVarSDID,ii);           
            fprintf(PMFd,PMNFmtS,idi);        
            fprintf(PMFd,PMNFmtS,idi);        
            fprintf(PMFd,PMNFmtS,-1);        
            fprintf(PMFd,"\n");        
            nrec++;
        }
    }
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

SDIPFin:
    p_clean();
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

int sdcpol(void)
{
    register int i,j,k,l;
    int err,r,ixa,iya,ixb,iyb,n,m,nn,nrec,nc,ncn,ne,sdid,npmax,bnmax,nid;
    double x,y,x1,y1;

    err = -1;
    nrec = 0;
    printf1("Construction of polygons from line segments. Current memory: %d bytes.\n",MemReq);
         
    if (parm(CmdBuf + 6,4,1))       /* get parameters */
        goto SDCPFin;

    if (PMNV != 4) {
        printf1("Error: need four variables on right-hand side.\n");
        goto SDCPFin;
    }
    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMF1Def == 0) {
        printf1("Error: need an output file.\n");
        goto SDCPFin;
    }
    if (NOC < 3) {
        printf1("Error: need at least three line segments.\n");
        goto SDCPFin;
    }
    if (PMOPT < 1 || PMOPT > 3)
        PMOPT = 1;
    if (PMMax < 1)
        PMMax = 2 * NOC;

    PMTol = dmax(PMTol,EPSI1);
    PMMax = imax(2 * NOC,PMMax);
    printf1("Tolerance: %g\n",PMTol);
    printf1("Maximal number of segments: %d\n",PMMax);

    /* Create the graph. First read segments, then call sdp_list(). */

    if (alloc_acxf(PMMax + 2))
        goto SDCPFin;     
    if (alloc_acyf(PMMax + 2))
        goto SDCPFin;     
    if (alloc_acuf(PMMax + 2))
        goto SDCPFin;     
    if (alloc_acvf(PMMax + 2))
        goto SDCPFin;     
    if (alloc_ack(PMMax + 2))
        goto SDCPFin;     
    if (alloc_acd(PMMax + 2))
        goto SDCPFin;     

    ixa = PMVIdx[0];
    iya = PMVIdx[1];
    ixb = PMVIdx[2];
    iyb = PMVIdx[3];

    for (i = 0; i < NOC; ++i) {
        AcXF[i] = (float)get_data(ixa,i);
        AcYF[i] = (float)get_data(iya,i);
        AcUF[i] = (float)get_data(ixb,i);
        AcVF[i] = (float)get_data(iyb,i);
    }
    r = sdp_list(NOC,AcXF,AcYF,AcUF,AcVF,AcK,AcD,PMMax,(float)PMTol);         
    if (r != 0)  
        goto SDCPFin;

    alloc_acxf(0);   /* free memory */
    alloc_acyf(0);
    alloc_acuf(0);
    alloc_acvf(0);
    alloc_ack(0);
    alloc_acd(0);
      
    nn = SD_NN;     /* number of nodes */
    ne = SD_NE;     /* number of edges */

    printf1("Number of nodes: %d\n",nn);
    printf1("Number of nodes: %d\n",ne);
    if (nn < 3 || ne < 3) {
        printf1("Will not continue.\n");
        err = 0;
        goto SDCPFin;
    }
    if (PMOPT == 1) {                       /* print node list */
        for (i = 1; i <= nn; ++i) {
            fprintf(PMF1d,PMNFmtS,i);
            fprintf(PMF1d,PMFmtS,(double)SD_NLX[i]);
            fprintf(PMF1d,PMFmtS,(double)SD_NLY[i]);
            fprintf(PMF1d,"\n");
            nrec++;
        }
        printf1("%d records written to: %s\n",nrec,PMF1dName);
        err = 0;
        goto SDCPFin;
    }

    /* for each node i = 1,...,nn, sort the adjacent nodes in
       counterclockwise direction */

    printf1("Sorting nodes in counterclockwise orientation.\n");

    if (alloc_actmp(nn + 1))
        goto SDCPFin;     

    for (i = 1; i <= nn; ++i) {
        if (SD_NA[i] <= 2)                /* no sorting necessary */
            continue;
        sdp_sort(i,SD_EPtr[i],SD_NA[i],SD_EList,SD_NLX,SD_NLY,AcTmp);    
    }
    alloc_actmp(0);     /* free */

    /* check for multiple edges (segments) */

    n = 0;
    for (i = 1; i <= nn; ++i) {
        k = SD_EPtr[i];
        j = SD_EList[k];
        for (l = 1; l < SD_NA[i]; ++l) {
            if (SD_EList[k + l] == j) {
                SD_EList[k + l] = -1;
                n++;
            }
            else
                j = SD_EList[k + l];
        }
    }
    if (n > 0) {
        printf1("Dropped %d segment(s) occurring at least twice.\n",n);
        printf1("Remaining number of segments: %d\n",ne - n);
    }


    if (PMOPT == 2) {                   /* write edge list */
        for (i = 1; i <= nn; ++i) {
            k = SD_EPtr[i];
            for (l = 0; l < SD_NA[i]; ++l) {
                j = SD_EList[k + l];
                if (j < 1)
                    continue;

                fprintf(PMF1d,PMNFmtS,i);
                fprintf(PMF1d,PMNFmtS,j);
                fprintf(PMF1d,PMFmtS,(double)SD_NLX[i]);
                fprintf(PMF1d,PMFmtS,(double)SD_NLY[i]);
                fprintf(PMF1d,PMFmtS,(double)SD_NLX[j]);
                fprintf(PMF1d,PMFmtS,(double)SD_NLY[j]);
                fprintf(PMF1d,"\n");
                nrec++;
            }
        }
        printf1("%d records written to: %s\n",nrec,PMF1dName);
        err = 0;
        goto SDCPFin;
    }
                   
    /* Prepare search for polygons. First remove nodes with degree less
       than two. SD_Deg will be used to store the degrees. */

    printf1("Searching for polygons.\n");
    printf1("Removing nodes with degree less than 2.\n");

    if (alloc_aci(nn + 1))
        goto SDCPFin;     

    /* create the initial list of degrees in SD_Deg */

    for (i = 1; i <= nn; ++i) {
        k = SD_EPtr[i];
        j = SD_EList[k];
        n = 0;
        for (l = 0; l < SD_NA[i]; ++l) {
            if (SD_EList[k + l] >= 1)   
                n++;
        }
        SD_Deg[i] = n;
    }
    if (sdp_clear(nn) < 0) {
        printf1("Fatal error in degree calculation.\n");
        goto SDCPFin;
    }
    /***********/
    printf("NEUE DEGREES\n");
    for (i = 1; i <= nn; ++i) {
        printf("i=%3d aci=%4d\n",i,SD_Deg[i]);
    }   

    printf("NEW SD_EList\n");
    for (i = 1; i <= nn; ++i) {
        k = SD_EPtr[i];
        for (l = 0; l < SD_NA[i]; ++l) {
            j = SD_EList[k + l];
            printf("i=%3d k=%3d j=%4d\n",i,k,j);
        }
    }
    /*******/

    n = 0;
    for (i = 1; i <= nn; ++i) {
        if (SD_Deg[i] >= 2)
            n++;             
    }
    printf1("Remaining number of nodes: %d\n",n);
    if (n < 3)
        goto SDCPFin;

    /* ### Remove all cut nodes with degree less than 3. */

    printf1("Removing cut nodes with degree less than 3.\n");

    if (alloc_acl(nn + 1))      /* work space */
        goto SDCPFin;     
    if (alloc_acd(nn + 1))      
        goto SDCPFin;     
                   
    sdp_cut(nn,AcL,AcD);    

    /* Remove all nodes with degree less than 2. */

    if (sdp_clear(nn) < 0) {
        printf1("Fatal error in degree calculation.\n");
        goto SDCPFin;
    }
    /****************************
    printf("X NEUE DEGREES\n");
    for (i = 1; i <= nn; ++i) {
        printf("i=%3d aci=%4d\n",i,SD_Deg[i]);
    }   
    printf("NEW SD_EList\n");
    for (i = 1; i <= nn; ++i) {
        k = SD_EPtr[i];
        for (l = 0; l < SD_NA[i]; ++l) {
            j = SD_EList[k + l];
            printf("i=%3d k=%3d j=%4d\n",i,k,j);
        }
    }
    *************/


    /* Loop over all components, nc is counter for components. The index
       number of the component is stored in AcL. */

    npmax = 2 * nn;             /* max number of nodes in a polygon */
    bnmax = 2 * nn;             /* max number of nodes on boundary */

    if (alloc_acm(npmax + 2))   /* used to save nodes of polygons */
        goto SDCPFin;     
    if (alloc_acs(bnmax + 2))   /* used to save nodes on boundary */
        goto SDCPFin;     

    n = 0;
    for (i = 1; i <= nn; ++i) {
        if (SD_Deg[i] < 2)
            AcL[i] = -1;
        else {
            AcL[i] = 0;
            n++;
        }
    }
    printf1("Remaining number of nodes: %d.\n\n",n);                  
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
            if (AcL[i] == 0) {
                k = i;
                break;
            }
        }
        if (k < 0)              /* no more components */
            break;

        SDC_N = 0;              /* number of elements */
        nc++;                   /* counts the components */
        sdp_comp(k,nc,AcL);     /* find next component */

        printf1("Component: %d. Nodes: %d. Polygons: ",nc,SDC_N);
        if (SDC_N < 3) {
            printf1("0\n");
            continue;
        }

        /* find left-most node in current component to begin the search. */

        k = -1;
        x = y = 0.0;
        for (i = 1; i <= nn; ++i) {
            if (AcL[i] != nc)
                continue;

            if (k < 0 || SD_NLX[i] < (float)x    
                      || (SD_NLX[i] == (float)x && SD_NLY[i] > (float)y)) {
                x = (double)SD_NLX[i];
                y = (double)SD_NLY[i];
                k = i;
            }
        }
        printf("k====%d\n",k);

        r = sdp_fpol(nn,nc,AcL,k,SD_NLX,SD_NLY,AcD,AcS,bnmax,AcM,npmax,&sdid,&m);
        if (r < 0) {
            printf1("Error: %d\n",r);
            goto SDCPFin;
        }
        else {
            printf1("Polygons: %d\n",r);
            nrec += m;
        }

    }
    printf1("Found %d polygon(s).\n",sdid);
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    err = 0;

SDCPFin:
    sdp_list_free();
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdp_id(xa,ya,xb,yb,tol)                                                 */
/*                                                                          */
/*  Return 1 if (xa,ya) and (xb,yb) are identical (up to tol), otherwise    */  
/*  return 0.                                                               */

int sdp_id(double xa,double ya,double xb,double yb,double tol)
{
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

void sdp_sort(int i,int k,int n,int *nptr,float *x,float *y,double *a)
{
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

int sdp_fpol(int nn,int nc,int *cn,int rn,float *xc,float *yc,char *idx,
    int *bnodes,int nbmax,int *pnodes,int npmax,int *sdid,int *nrec)
{
    register int i,j,k,l,ll;
    int r,ii,n,m,nnc,i0,j0,j1,j2,nb,np,fin,rnl;
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

        k = SD_EPtr[rn];
        n = SD_NA[rn];

        rnl = -1;
        a0 = 0.0;
        for (ll = 0; ll < n; ++ll) {
            if ((j = SD_EList[k + ll]) >= 1) { 
                a = atan2((double)(yc[j] - yc[rn]),(double)(xc[j] - xc[rn]));
                if (rnl < 0 || a0 < a) {
                    rnl = j;
                    a0 = a;
                }
            }
        }
printf("rn=%d     most left rnl=%d a0=%f\n",rn, rnl,a0);

        if (rnl < 0) {
            printf("ERROR??? rn=%d  rnl=%d\n",rn, rnl);
            break;
        }


        /* find second node j0 */

        j0 = -1;
        xm = xc[rn];
        ym = yc[rn];
        a0 = 0.0;
        k = SD_EPtr[rn];            /* consider all nodes adjacent to rn */
        n = SD_NA[rn];
        for (ll = 0; ll < n; ++ll) {
            j = SD_EList[k + ll];        
            if (j < 1)
                continue;

            a = atan2(yc[j] - ym,xc[j] - xm);
            if (j0 < 0 || a < a0) {
                j0 = j;
                a0 = a;
            }
        }
        nb++;
        bnodes[nb] = j0;
        idx[j0] = 1;

               
        printf("found second j0=%d\n",j0);

        if (j0 < 1) {
            printf("ERROR???\n");
            break;
        }


        /* find remaining nodes on the boundary */

        j1 = j0;
        j0 = rn;
        while (1) {
            if ((l = sdp_find_idx(j1,j0)) < 0)      /* index to edge list */
                return(-2);
            if ((j2 = sdp_find_nxt(j1,l,0)) < 0)    /* next node */
                return(-3);

            printf("found j2=%d\n",j2);

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
        printf("\nBNodes: ");
        for (i = 1; i <= nb; ++i)
            printf("%d ",bnodes[i]);
        newline();


        /* find all polygons that have at least one segment on the boundary */
        /* store nodes of polygons in pnodes[], counter is np. */
       
        for (ii = 1; ii < nb; ++ii) {

            j0 = i0 = bnodes[ii];
            if (SD_Deg[i0] < 2) {
                printf("xcontinue\n");
                continue;
            }
            j1 = bnodes[ii + 1];
            k = SD_EPtr[j0];
            n = SD_NA[j0];
            l = -1;
            for (ll = 0; ll < n; ++ll) {
                if (SD_EList[k + ll] == j1) {
                    l = ll;
                    break;
                }
            }
            if (l < 0) {
                printf("continue\n");
                continue;
            }
            np = 0;
            pnodes[++np] = j0;
            pnodes[++np] = j1;

            printf("found j0=%d j1=%d\n",j0,j1);

            /* find additional points until arriving back at i0 */

            while (1) {

                if ((l = sdp_find_idx(j1,j0)) < 0) {   /* index to edge list */

printf("ERROR -4 j0=%d j1=%d\n",j0,j1 );

k = SD_EPtr[j0];
n = SD_NA[j0];
for (ll = 0; ll < n; ++ll)  
    printf("j0=%d k=%d n=%d elist=%d\n",j0,k,n,SD_EList[k + ll]);

k = SD_EPtr[j1];
n = SD_NA[j1];
for (ll = 0; ll < n; ++ll)  
    printf("j1=%d k=%d n=%d elist=%d\n",j1,k,n,SD_EList[k + ll]);



                    return(-4);
             }
                if ((j2 = sdp_find_nxt(j1,l,1)) < 0)   /* next node */
                    return(-5);
      
                printf("Next point: %d\n",j2);
                       
                if (j2 == i0) {
                    printf("found new polygon. np=%d\n",np);

                    /* write polygon to PMF1d. spatial data file format.   
                       count number of records in nrec. count polygon IDs
                       in sdid. */

                    *sdid += 1;
                    *nrec += sdp_write(*sdid,np,pnodes,xc,yc);
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
                        if (sdcpol_edge(j0,j1,nb,bnodes,idx) == 1) {
                            sdp_deg(j0,j1);       
                            sdp_deg(j1,j0);       
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
printf("NOT VALID j=%d j2=%d    \n",j,j2   );

                    /* adjust degrees for the first two nodes. */

                    if (sdcpol_edge(pnodes[1],pnodes[2],nb,bnodes,idx) == 1) {
                        sdp_deg(pnodes[1],pnodes[2]);       
                        sdp_deg(pnodes[2],pnodes[1]);       
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
                      
            printf("ende des polygons np=%d\n",np );
            for (i = 1; i <= np; ++i)
                printf("%d ",pnodes[i]);
            newline();

                   
            printf("NEW degrees\n");
            for (i = 1; i <= nn; ++i) 
                printf("i=%3d deg=%3d\n",i,SD_Deg[i]);
                   
           

            /* clear the current component from all nodes with degree < 2. */

            nnc = sdp_cclear(nn,nc,cn);
            if (nnc < 0)
                return(-9);

                   
            printf("FINAL degrees nc=%d nnc=%d     \n",nc,nnc  );
            for (i = 1; i <= nn; ++i) 
                printf("i=%3d deg=%3d cn=%4d  \n",i,SD_Deg[i],cn[i]  );
            
            if (nnc < 3)        /* Number of nodes is less than 3 */
                break;
        }
        if (nnc < 3)
            break;

        /* find root node for new reduced component. */

        rn = -1;
        x = y = 0.0;
        for (i = 1; i <= nn; ++i) {
            if (cn[i] != nc || SD_Deg[i] < 2)
                continue;

            if (rn < 0 || xc[i] < x || (xc[i] == x && yc[i] > y)) {
                x = xc[i];
                y = yc[i];
                rn = i;
            }
        }
        printf("NEW rn = %d\n",rn);

        if (rn < 0) {
            printf("WARUM nnc=%d\n",nnc);

            break;
        }
    }
    return(r);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_comp(k,nc,cn)  Create next component with level nc. Use all nodes   */
/*                     with cn[i] = 0. Count number of elements in SDC_N    */

void sdp_comp(int k,int nc,int *cn)
{  
    register int l,j,n;

    SDC_N++;
    cn[k] = nc;
    n = SD_NA[k];
    if (n <= 0)
        return;

    for (l = 0; l < n; ++l) {
        j = SD_EList[SD_EPtr[k] + l];            /* edge from k to j */
        if (j >= 1 && cn[j] == 0)  
            sdp_comp(j,nc,cn);
    }
}

/* -###-------------------------------------------------------------------- */
/*  sdp_deg(i,j)   Remove node j from the followers of i in edge list.      */
/*                 Return 1 if removed, otherwise 0                         */

int sdp_deg(int i,int j)
{
    register int k,l,n;

    k = SD_EPtr[i];
    n = SD_NA[i];
    for (l = 0; l < n; ++l) {
        if (SD_EList[k + l] == j) {
            SD_EList[k + l] = -1;
            SD_Deg[i] -= 1;
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

int sdcpol_edge(int i,int j,int n,int *bnodes,char *idx)
{
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

int sdp_clear(int nn)
{
    register int i,j,k,ll;
    int n,nd,fin;
                
    nd = fin = 0;
    while (fin == 0) {
        fin = 1;
        nd = 0;
        for (i = 1; i <= nn; ++i) {
            if (SD_Deg[i] < 1)
                continue;
            else if (SD_Deg[i] > 1)
                nd++;
            else {
                k = SD_EPtr[i];
                n = SD_NA[i];
                j = -1;
                for (ll = 0; ll < n; ++ll) {
                    if ((j = SD_EList[k + ll]) >= 1)  
                        break;
                }
                if (j < 0)
                    return(-1);

                sdp_deg(j,i);       
                SD_EList[k + ll] = -1;
                SD_Deg[i] -= 1;
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

int sdp_cclear(int nn,int nc,int *cn)
{
    register int i,j,k,ll;
    int n,nd,fin;
                
    nd = fin = 0;
    while (fin == 0) {
        fin = 1;
        nd = 0;
        for (i = 1; i <= nn; ++i) {
            if (cn[i] != nc || SD_Deg[i] < 1)
                continue;
            else if (SD_Deg[i] > 1)
                nd++;
            else {
                k = SD_EPtr[i];
                n = SD_NA[i];
                j = -1;
                for (ll = 0; ll < n; ++ll) {
                    if ((j = SD_EList[k + ll]) >= 1)  
                        break;
                }
                if (j < 0)
                    return(-1);

                sdp_deg(j,i);       
                SD_EList[k + ll] = -1;
                SD_Deg[i] = 0;
                fin = 0;
            }
        }
    }
    return(nd);
}

/* -###-------------------------------------------------------------------- */
/*  sdp_remove(i)    Remove node i from graph, adjust SD_Deg[].             */

void sdp_remove(int i)
{
    register int k,l,j,n;

    k = SD_EPtr[i];
    n = SD_NA[i];
    for (l = 0; l < n; ++l) {
        if ((j = SD_EList[k + l]) >= 1) {
            SD_EList[k + l] = -1;
            SD_Deg[i] -= 1;
            sdp_deg(j,i);
        }
    }
}


/* -###-------------------------------------------------------------------- */
/*  sdp_cut(nn,vn,idx)                                                      */
/*                                                                          */  
/*  nn is number of nodes. vn[] and idx[] are work arrays of length nn.     */

void sdp_cut(int nn,int *vn,char *idx)
{
    register int i,j;
    int k,nid; 

    SDC_ID = 0;

    for (i = 1; i <= nn; ++i)
        vn[i] = 0;

    while (1) {

        /* find an unprocessed node */

        k = -1;
        for (i = 1; i <= nn; ++i) {
            if (SD_Deg[i] >= 2 && vn[i] == 0) {
                k = i;                    
                break;
            }
        }
        if (k < 0)      /* no more components */
            break;

        /* Find the next component and, at the same time, find the
           cut nodes. */

        SDC_N = 0;                  /* number of elements in component */
        SDC_CNT = 0;                /* number of path from root */
        nid = SDC_ID;               /* nodes with vn[i] > nid belong to   
                                       current component. */
        for (i = 1; i <= nn; ++i)
            idx[i] = 0;

        sdp_visit(k,k,vn,idx);

        if (SDC_CNT < 2)
            idx[k] = 0;

        printf("\nNext component nnn_nn=%d nid=%d\n",SDC_N,nid );
/**
        for (j = 1; j <= nn; ++j)
            printf("j=%3d nnn=%4d nnns=%4d  SD_Deg=%3d\n",j,vn[j],idx[j],SD_Deg[j] );
**/ 

        for (j = 1; j <= nn; ++j) {
            if (SD_Deg[j] != 2)   
                idx[j] = 0;
        }     
        for (j = 1; j <= nn; ++j) {
            if (idx[j]) {

printf("REMOVE %d\n",j);
                sdp_remove(j);
            }
        }
    }
    printf("NO MORE\n");
}

/* ------------------------------------------------------------------------ */
/*  sdp_visit(i,ir,vn,idx)    Called by sdp_cut().                          */

int sdp_visit(int i,int ir,int *vn,char *idx)
{
    register int j,k,ll;
    int n,m,min;

    vn[i] = ++SDC_ID;
    min = SDC_ID;
    SDC_N++; 

    k = SD_EPtr[i];
    n = SD_NA[i];
    for (ll = 0; ll < n; ++ll) {
        if ((j = SD_EList[k + ll]) >= 1) {  /* edge from i to j */
            if (vn[j] == 0) {
                if (i == ir)
                    SDC_CNT++;

                m = sdp_visit(j,ir,vn,idx);
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

int sdp_find_idx(int i,int j)
{
    register int l,k,n;

    k = SD_EPtr[i];
    n = SD_NA[i];

    for (l = 0; l < n; ++l) {
        if (SD_EList[k + l] == j)    
            return(l);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  sdp_find_nxt(i,l,opt)   i is node number, l is index to its adjacent    */
/*                          nodes in the edge list. If opt = 0 return next  */
/*                          higher node number, otherwise next lower node   */
/*                          number. Or -1 if there is no adjacent node.     */

int sdp_find_nxt(int i,int l,int opt)
{
    register int ll,k,n,j;

    k = SD_EPtr[i];
    n = SD_NA[i];

    if (opt == 0) {
        for (ll = l + 1; ll < n; ++ll) {
            if ((j = SD_EList[k + ll]) >= 1)  
                return(j);           
        }
        for (ll = 0; ll < l; ++ll) {
            if ((j = SD_EList[k + ll]) >= 1)  
                return(j);               
        }
    }                 
    else {
        for (ll = l - 1; ll >= 0; --ll) {
            if ((j = SD_EList[k + ll]) >= 1)  
                return(j);           
        }
        for (ll = n - 1; ll > l; --ll) {
            if ((j = SD_EList[k + ll]) >= 1)  
                return(j);               
        }
    }
    return(-1);
}                 

/* ------------------------------------------------------------------------ */
/*  sdp_write(id,np,pnodes,xc,yc)                                           */
/*                                                                          */  
/*  Write polygon to PMF1d. Return number of records written.               */

int sdp_write(int id,int np,int *pnodes,float *xc,float *yc)
{
    register int i,j;
    int nrec = 0;

    fprintf(PMF1d,PMNFmtS,id);
    fprintf(PMF1d,"3 ");
    fprintf(PMF1d,PMNFmtS,np);
    fprintf(PMF1d,"\n");
    nrec++;         
    for (i = 1; i <= np; ++i) {
        j = pnodes[i];
        fprintf(PMF1d,PMFmtS,(double)xc[j]);
        fprintf(PMF1d,PMFmtS,(double)yc[j]);
        fprintf(PMF1d,"\n");
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

int sdp_list(int n,float *xa,float *ya,float *xb,float *yb,int *ptr,char *idx,
    int nmax,float tol) 
{
    register int i,j,k;
    int err,m,nne,nnd,nnl,r,ri,rj,fin,nfree,xfa,yfa,ii,jj;
    float t,xa0,ya0,*xf,*yf;
    double d,xai,yai,xbi,ybi,xaj,yaj,xbj,ybj,xxa,yya,xxb,yyb;

    err = -2;
    xfa = yfa = 0;

    for (i = 0; i < n; ++i) {
        printf("i=%3d %f %f %f %f\n",i,xa[i],ya[i],xb[i],yb[i]);
    }
    newline();

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
        if (sdp_id(xa[i],ya[i],xb[i],yb[i],tol))   
            continue;

        xa[nne] = xa[i];
        ya[nne] = ya[i];
        xb[nne] = xb[i];
        yb[nne] = yb[i];
        nne++;
    }

    printf("drop length less than\n");
    for (i = 0; i < nne; ++i) {
        printf("i=%3d %f %f %f %f\n",i,xa[i],ya[i],xb[i],yb[i]);
    }
    newline();

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

            if (sdp_id(xai,yai,xbi,ybi,tol)) {       /* zero length */

            /**
            printf("found zero lengh i=%d  will free i=%d nfree=%d   \n",i,i,nfree  );
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

                if (sdp_id(xaj,yaj,xbj,ybj,tol)) {       /* zero length */
                    /**
                    printf("found zero lengh j=%d  will free j=%d nfree=%d   \n",j,j,nfree  );
                    **/
                    ptr[nfree] = j;
                    idx[j] = 1;
                    nfree++;
                    continue;
                }
                /**
                printf("i=%d j=%d xi=%g %g %g %g xj=%g %g %g %g\n",i,j,xai,yai,xbi,ybi,
                        xaj,yaj,xbj,ybj);
                **/

                if ((sdp_id(xai,yai,xaj,yaj,tol)  && 
                     sdp_id(xbi,ybi,xbj,ybj,tol)) ||    
                    (sdp_id(xai,yai,xbj,ybj,tol)  && 
                     sdp_id(xaj,yaj,xbi,ybi,tol))) {     
                    /**
                    printf("found identical. will free j=%d nfree=%d   \n",j,nfree  );
                    **/
                    ptr[nfree] = j;
                    idx[j] = 1;
                    nfree++;
                    continue;
                }
                r = g_intersect(xai,yai,xbi,ybi,xaj,yaj,xbj,ybj,&xxa,&yya,&xxb,&yyb);
                /**
                printf("rrrr=%d\n",r);
                **/
                if (r == 0)
                    continue;

                if (r == 1) {           /* crossing in single point */

                    if (sdp_id(xai,yai,xxa,yya,tol) || 
                        sdp_id(xbi,ybi,xxa,yya,tol))       
                        ri = 0;
                    else
                        ri = 1;

                    if (sdp_id(xaj,yaj,xxa,yya,tol) || 
                        sdp_id(xbj,ybj,xxa,yya,tol))       
                        rj = 0;
                    else
                        rj = 1;
                    /**
                    printf("ri=%d rj=%d\n",ri,rj);
                    **/
                    if (ri) {           /* i needs splitting */

                        /**
                        printf("setzt bei xb[i] ein: %g %g\n",xxa,yya);
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
                        printf("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xbi,ybi);
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
                        printf("setzt bei xb[j] ein: %g %g\n",xxa,yya);
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
                        printf("setzt bei k=%d ein: %g %g %g %g\n",k, xxa,yya,xbj,ybj);
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
                    if (sdp_id(xxa,yya,xxb,yyb,tol) == 0) {
                        /**
                        printf("xxa=%g %g xxb=%g %g nicht identisch\n",xxa,yya,xxb,yyb);
                        **/
                        if ((sdp_id(xxa,yya,xai,yai,tol)  &&
                             sdp_id(xxb,yyb,xbi,ybi,tol)) ||
                            (sdp_id(xxb,yyb,xai,yai,tol)  &&
                             sdp_id(xxa,yya,xbi,ybi,tol))) {
                            /**
                            printf("xxa usw identisch mit i.\n");
                            printf("setzt bei xb[j] ein: %g %g\n",xxa,yya);
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
                            printf("setzt bei k=%d ein: %g %g %g %g\n",k,xxb,yyb,xbj,ybj);
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
                        else if ((sdp_id(xxa,yya,xaj,yaj,tol)  &&
                                  sdp_id(xxb,yyb,xbj,ybj,tol)) ||
                                 (sdp_id(xxb,yyb,xaj,yaj,tol)  &&
                                  sdp_id(xxa,yya,xbj,ybj,tol))) {

                            /**
                            printf("setzt bei xb[i] ein: %g %g\n",xxa,yya);
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
                            printf("setzt bei k=%d ein: %g %g %g %g\n",k,xxb,yyb,xbi,ybi);
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
                        else if (sdp_id(xxa,0.0,xxb,0.0,tol)) {
                            if (yai < yya) {
                                /**
                                printf("yai kleiner als yya\n");
                                printf("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                                printf("setzt bei xa[j] ein: %g %g\n",xxb,yyb);
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
                                printf("yai nicht  kleiner als yya\n");
                                printf("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                                printf("setzt bei xa[i] ein: %g %g\n",xxb,yyb);
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
                            printf("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xxb,yyb);
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
                                printf("xai kleiner als xxa\n");
                                printf("setzt bei xb[i] ein: %g %g\n",xxa,yya);
                                printf("setzt bei xa[j] ein: %g %g\n",xxb,yyb);
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
                                printf("xai nicht  kleiner als xxa\n");
                                printf("setzt bei xb[j] ein: %g %g\n",xxa,yya);
                                printf("setzt bei xa[i] ein: %g %g\n",xxb,yyb);
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
                            printf("setzt bei k=%d ein: %g %g %g %g\n",k,xxa,yya,xxb,yyb);
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
        printf("After nne=%d nfree=%d: ",nne,nfree );
        for (k = 0; k < nfree; ++k)
            printf("%d ",ptr[k]);
        newline();

        for (i = 0; i < nne; ++i) {
            printf("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]);
        }
        newline();
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
            printf("After reorg nne=%d nfree=%d: ",nne,nfree );
            for (i = 0; i < nne; ++i) {
                printf("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]);
            }
            newline();
            **/
        }
    }

    printf("After splitting nne=%d nnd=%d nfree=%d: ",nne,nnd,nfree );
    for (k = 0; k < nfree; ++k)
        printf("%d ",ptr[k]);
    newline();

    for (i = 0; i < nne; ++i) {
        printf("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]  );
    }
    newline();

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

        printf("After reorg nne=%d nfree=%d: ",nne,nfree );
        for (i = 0; i < nne; ++i) {
            printf("i=%3d %f %f %f %f idx=%d\n",i,xa[i],ya[i],xb[i],yb[i],idx[i]);
        }
        newline();
    }

    err = -1;
/***
    if (sortdp4f(nne,xa,ya,xb,yb,ptr))
        goto SDPListFin;

    for (i = 0; i < nne; ++i) {
        k = ptr[i];
        printf("i=%3d %f %f %f %f\n",i,xa[k],ya[k],xb[k],yb[k]);
    }
    newline();
**/        
    /* create list of unique coordinates */

    m = 2 * nne;
    if (!(SD_EPtr = (int *)calloc(m,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(m,sizeof(int));
    SD_EPtrA = m;     

    if (!(xf = (float *)calloc(m,sizeof(float)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(m,sizeof(float));
    xfa = m;             

    if (!(yf = (float *)calloc(m,sizeof(float)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(m,sizeof(float));
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
    if (sortdp2f(m,xf,yf,SD_EPtr))
        goto SDPListFin;

    printf("xf yf array\n");

    for (i = 0; i < m; ++i) {
        k = SD_EPtr[i];
        printf("i=%3d %f %f\n",i,xf[k],yf[k]);
    }

    nnd = 0;                /* number of unique coordinates (nodes) */
    k = SD_EPtr[0];
    xa0 = xf[k] - 1000.0;
    ya0 = yf[k] - 1000.0;
    for (i = 0; i < m; ++i) {
        k = SD_EPtr[i];
        if (sdp_id((double)xf[k],(double)yf[k],(double)xa0,(double)ya0,tol) == 0) {
            nnd++;
            xa0 = xf[k];
            ya0 = yf[k];
        }
    }

    printf("nnd=%d\n",nnd);

    /* save unique list of coordinates in SD_NLX and SD_NLY */

    if (!(SD_NLX = (float *)calloc(nnd + 2,sizeof(float)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnd + 2,sizeof(float));
    SD_NLXA = nnd + 2;       

    if (!(SD_NLY = (float *)calloc(nnd + 2,sizeof(float)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnd + 2,sizeof(float));
    SD_NLYA = nnd + 2;       

    j = 0;
    k = SD_EPtr[0];
    xa0 = xf[k] - 1000.0;
    ya0 = yf[k] - 1000.0;
    for (i = 0; i < m; ++i) {
        k = SD_EPtr[i];
        if (sdp_id((double)xf[k],(double)yf[k],(double)xa0,(double)ya0,tol) == 0) {
            j++;
            SD_NLX[j] = xa0 = xf[k];
            SD_NLY[j] = ya0 = yf[k];
        }
    }
    if (j != nnd) {
        err = -3;
        goto SDPListFin;
    }

    printf("LIST of node coordinates nnd=%d\n",nnd);

    for (i = 1; i <= nnd; ++i) {
        printf("i=%3d x=%f y=%f\n",i,SD_NLX[i],SD_NLY[i]);
    }

    free((char *)xf);
    memrq(-xfa,sizeof(float));
    xfa = 0;
    free((char *)yf);
    memrq(-yfa,sizeof(float));
    yfa = 0;
    free((char *)SD_EPtr);
    memrq(-SD_EPtrA,sizeof(int));
    SD_EPtrA = 0;

    /* create forward-linked edge list in SD_EList, SD_EPtr, SD_NA */

    nnl = 2 * nne;
    if (!(SD_EList = (int *)calloc(nnl + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnl + 1,sizeof(int));
    SD_EListA = nnl + 1;     

    if (!(SD_NA = (int *)calloc(nnl + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnl + 1,sizeof(int));
    SD_NAA = nnl + 1;     

    j = 0;
    for (i = 0; i < nne; ++i) {
        printf("i=%3d %f %f %f %f\n",i,xa[i],ya[i],xb[i],yb[i]);
        
        ii = sdp_list_find(xa[i],ya[i],nnd,SD_NLX,SD_NLY,tol);
        jj = sdp_list_find(xb[i],yb[i],nnd,SD_NLX,SD_NLY,tol);
        if (ii < 1 || jj < 1) {
            err = -3;
            goto SDPListFin;
        }
        j++;         
        SD_NA[j] = ii;
        SD_EList[j] = jj;
        j++;         
        SD_NA[j] = jj;
        SD_EList[j] = ii;
    }
    if (sorti2(nnl,SD_NA + 1,SD_EList + 1))
        goto SDPListFin;

    printf("Kantenliste\n");
    for (i = 1; i <= nnl; ++i)
        printf("i=%3d ii=%3d jj=%3d\n",i,SD_NA[i],SD_EList[i]);


    if (!(SD_EPtr = (int *)calloc(nnd + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnd + 1,sizeof(int));
    SD_EPtrA = nnd + 1;     

    k = -1;
    j = 0;
    for (i = 1; i <= nnl; ++i) {
        if (SD_NA[i] != k) {
            SD_EPtr[++j] = i;
            k = SD_NA[i];
        }
    }
    if (j != nnd) {
        err = -3;
        goto SDPListFin;
    }
    free((char *)SD_NA);
    memrq(-SD_NAA,sizeof(int));
    SD_NAA = 0;

    if (!(SD_NA = (int *)calloc(nnd + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnd + 1,sizeof(int));
    SD_NAA = nnd + 1;     

    for (i = 1; i < nnd; ++i)
        SD_NA[i] = SD_EPtr[i + 1] - SD_EPtr[i];
    SD_NA[nnd] = nnl - SD_EPtr[nnd] + 1;

    /* create degrees in SD_Deg */

    if (!(SD_Deg = (int *)calloc(nnd + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto SDPListFin;
    }   
    memrq(nnd + 1,sizeof(int));
    SD_DegA = nnd + 1;     

    for (i = 1; i <= nnd; ++i)
        SD_Deg[i] = SD_NA[i];


    printf("Pointer und NE\n");
    for (i = 1; i <= nnd; ++i)
        printf("i=%3d ptr=%5d ne=%4d\n",i,SD_EPtr[i],SD_NA[i]);

    SD_NN = nnd;
    SD_NE = nnl;
    err = 0;

SDPListFin:
    if (err == -2)  
        printf1("Error: exceeded the maximum number of segments (%d).\n",nmax);
    else if (err == -3) {
        printf1("Error while creating the edge list.\n");
        printf1("Might be caused by rounding errors.\n");
        printf1("Cannot continue.\n");
    }
    if (xfa > 0) {
        free((char *)xf);
        memrq(-xfa,sizeof(float));
    }
    if (yfa > 0) {
        free((char *)yf);
        memrq(-yfa,sizeof(float));
    }
    if (err)
        sdp_list_free();

    return(err);
}


/* ------------------------------------------------------------------------ */
/*  sdp_list_find(xa0,ya0,n,xa,ya,tol)                                      */
/*                                                                          */
/*  Return index of (xa0,ya0) in list xa[i],ya[i] (i = 1,n) or -1 if        */
/*  not found.                                                              */
 
int sdp_list_find(float xa0,float ya0,int n,float *xa,float *ya,float tol)
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
                if (sdp_id((double)xa0,(double)ya0,(double)xa[i],(double)ya[i],(double)tol))        
                    return(i);
            }           
            for (i = k - 1; i >= 1; --i) {
                if ((float)fabs((double)(xa[i] - xa[k])) > tol)
                    break;
                if (sdp_id((double)xa0,(double)ya0,(double)xa[i],(double)ya[i],(double)tol))        
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
 
void sdp_list_free(void)
{
    if (SD_NLXA > 0) {
        free((char *)SD_NLX);
        memrq(-SD_NLXA,sizeof(float));
        SD_NLXA = 0;
    }
    if (SD_NLYA > 0) {
        free((char *)SD_NLY);
        memrq(-SD_NLYA,sizeof(float));
        SD_NLYA = 0;
    }
    if (SD_EPtrA > 0) {
        free((char *)SD_EPtr);
        memrq(-SD_EPtrA,sizeof(int));
        SD_EPtrA = 0;
    }
    if (SD_NAA > 0) {
        free((char *)SD_NA);
        memrq(-SD_NAA,sizeof(int));
        SD_NAA = 0;
    }
    if (SD_DegA > 0) {
        free((char *)SD_Deg);
        memrq(-SD_DegA,sizeof(int));
        SD_DegA = 0;
    }
    if (SD_EListA > 0) {
        free((char *)SD_EList);
        memrq(-SD_EListA,sizeof(int));
        SD_EListA = 0;
    }
}



