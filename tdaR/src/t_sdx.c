/****************************************************************************/
/*  t_sdx                                                                   */
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
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gf.h"
#include "t_stata.h"
#include "tda_context.h"

/*  functions in t_sdx.c */

int sdgen(TDAContext *ctx); 
void sdgenv(TDAContext *ctx, int irec); 
int sdshp(TDAContext *ctx); 
int shp_typ(TDAContext *ctx, int n,int opt);
void shp_dbfw(TDAContext *ctx, char *buf,int rlen);
int shp_dbfh(TDAContext *ctx, int opt,FILE *fd,char *buf,int *nrec,int *hlen,int *nvar,int *rlen,int ns);
int shp_dbff(TDAContext *ctx, FILE *fd,char *buf,int nvar);
int shp_dbf_a(TDAContext *ctx, int nvar,int opt);

int rdbf(TDAContext *ctx); 
int prn_dbf_var(TDAContext *ctx, int nvar);
int rcsv(TDAContext *ctx); 
int rplz(TDAContext *ctx); 

int gtopo(TDAContext *ctx); 
int gtopo_rpixel(TDAContext *ctx, int i,int j,int nc,int *err);
int sddcwp(TDAContext *ctx); 
int sdgshhs(TDAContext *ctx); 

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */
 
#define HLEN 100            /* index/shape file header length               */
#define RHLEN  8            /* index file record length                     */
#define DBFRLEN 32          /* dbf file record length                       */
#define DBFFLEN 11          /* dbf file field length                        */


/* ------------------------------------------------------------------------ */
/* Parameters for gshhs.                                                    */

#define GSHHSNP 1435084     /* max number of points in a polygon */     

#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) ((int)(((unsigned int)(i4) >> 24) + \
                    ((unsigned int)((i4) >> 8) & 65280u) + \
                    (((unsigned int)(i4) & 65280u) << 8) + \
                    (((unsigned int)(i4) & 255u) << 24)))

struct GSHHS {	
	int id;      		                /* Unique polygon id number, starting at 0 */
	int n;	                       	/* Number of points in this polygon */
	int level;	  	                 /* 1 land, 2 lake, 3 island_in_lake, 4 pond_in_island_in_lake */
	int west, east, south, north;	 /* min/max extent in micro-degrees */
	int area;		                    /* Area of polygon in 1/10 km^2 */
	short int greenwich;	          /* Greenwich is 1 if Greenwich is crossed */
	short int source;		            /* 0 = CIA WDBII, 1 = WVS */
};
struct	POINT {
	int	x;
	int	y;
};

/* ------------------------------------------------------------------------ */
/*  sdgen   Generating a spatial data file.                                 */
/*                                                                          */
/*          sdgen(                                                          */
/*              id=...,     name of SDID variable (required)                */
/*              xyv=...,    two variables with coordinates (required)       */
/*              av=...,     additional attribute variables (optional)       */
/*              typ=...,    type of spatial objects, def. 1                 */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 10.6          */
/*                                                                          */
/*          ) = output_file;                                                */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdgen(TDAContext *ctx)  
{
    register int i,j;
    int err,nrec,n,idx,idy,id,ida,ia;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Generating a spatial data file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto SDGENFin;
        
    if (ctx->PMNV != 2) {
        printf1(ctx, "Error: need two variables for coordinates.\n");
        goto SDGENFin;
    }
    if (ctx->PMID < 0) {
        printf1(ctx, "Error: need an ID variable.\n");
        goto SDGENFin;
    }
    if (ctx->PMTyp < 1)
        ctx->PMTyp = 1;

    if (ctx->PMTyp > 3) {
        printf1(ctx, "Error in type specification.\n");
        goto SDGENFin;
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,6);
     
    if (alloc_acx(ctx, ctx->NOC + 1))
        goto SDGENFin;
    if (alloc_acy(ctx, ctx->NOC + 1))
        goto SDGENFin;

    idx = ctx->PMVIdx[0];
    idy = ctx->PMVIdx[1];

    i = nrec = 0;
    while (i < ctx->NOC) {
        ia = i;
        ida = (int)get_data(ctx, ctx->PMID,ia);
        n = 0;
        while (i < ctx->NOC) {
            ctx->AcX[n] = get_data(ctx, idx,i);
            ctx->AcY[n] = get_data(ctx, idy,i);
            n++;
            i++;
            /* i has already been advanced, so the lookahead needs its
               own bound: on the last case it read element NOC */
            if (ctx->PMTyp == 1 || i >= ctx->NOC ||
                (id = (int)get_data(ctx, ctx->PMID,i)) != ida)
                break;
        }
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ida);
        if (n == 1)        
            j = 1;
        else if (n == 2)
            j = 2;
        else
            j = ctx->PMTyp;
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
        sdgenv(ctx, ia);
        fprintf(ctx->PMFd,"\n");
        nrec++; 

        for (j = 0; j < n; ++j) {
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcX[j]);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ctx->AcY[j]);
            fprintf(ctx->PMFd,"\n");
            nrec++; 
        }
        ida = id;
    }
      
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);

    err = 0;     

SDGENFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdgenv()    Print additional attribute variables. Called by sdgen().    */

void sdgenv(TDAContext *ctx, int irec)
{
    register int i,j;
    double tmp;

    for (i = 0; i < ctx->PM1NV; ++i) {
        j = ctx->PM1VIdx[i];
        if (ctx->VTyp[j] != 1) {
            tmp = get_data(ctx, j,irec);                  
            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[j],tmp);
        }
        else {          /* string variables */
            get_str(ctx, ctx->SVBuf,j,irec);
            fprintf(ctx->PMFd,"%s ",ctx->SVBuf);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  sdshp       Convert an ESRI shape file.                                 */
/*                                                                          */
/*              sdshp(                                                      */
/*                  df=...,     creates a spatial data file                 */
/*                  dtda=...,   TDA description file                        */
/*                  fmt=...,    print format data records, def. 12.4        */
/*              ) = path to shape file;                                     */
/*                                                                          */
/*  Reads a shape file according to the ESRI Shapefile Technical            */
/*  Description, July 1998.                                                 */
/*                                                                          */
/*  The right hand side should provide a path to the shape file without     */
/*  a file extension. It will be assumed that three files can be opened:    */
/*                                                                          */  
/*  *.shp   the shape file                                                  */
/*  *.shx   the index file                                                  */
/*  *.dbf   a dBASE file containing descriptions                            */
/*                                                                          */
/*  The command reads the files and prints some information to stdout.      */
/*  If the df parameter is used to specify an ouptut file, the command      */
/*  will write the data to this output file conforming to a standard        */
/*  TDA spatial data file.                                                  */
/*                                                                          */
/*  The following shape types are recognized:                               */
/*                                                                          */
/*  Type 0 : null shapes without geometrical data (ignored)                 */  
/*  Type 1 : point                                                          */
/*  Type 3 : polyline                                                       */
/*  Type 5 : polygon                                                        */
/*                                                                          */
/*  If the input file contains any other type of shape, the command will    */
/*  stop with an error message.                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdshp(TDAContext *ctx)  
{
    FILE *dbffd,*shpfd,*shxfd;
    register int i,j,k,l;
    int err,n,m,m1,mn,ns,np,typ,len,flen,nrec,hlen,nvar,rlen,n0,n1,drec,drec1;
    int dbfbufa,shpbufa,shpbufl,dbffd_open,shpfd_open,shxfd_open,npart,npnt;
    int id,id1,id2;
    double x,y,xmin,xmax,ymin,ymax;
    char *dbfbuf,*shpbuf,*p,*q;

    err = -1;
    nvar = shpbufa = dbfbufa = dbffd_open = shpfd_open = shxfd_open = 0;

    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Reading a shape file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,8,1)) {     /* get parameters */
        goto SDSHPFin;
    }   
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 12,4);

    printf1(ctx, "File: %s\n",ctx->PMRHSTR);

    shpbufl = 10000;
    if (!(shpbuf = (char *)calloc((size_t)(shpbufl + 2),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto SDSHPFin;
    }
    shpbufa = shpbufl + 2;
    memrq(ctx, shpbufa,sizeof(char));

    strcpy(shpbuf,ctx->PMRHSTR);                         /* open index file */
    strcat(shpbuf,".shx");
    if (!(shxfd = fopen(shpbuf,OPEN_RB))) {   
        printf1(ctx, "Error: can't open the index file.\n");         
        goto SDSHPFin;
    }
    shxfd_open = 1;

    strcpy(shpbuf,ctx->PMRHSTR);                         /* open shape file */
    strcat(shpbuf,".shp");
    if (!(shpfd = fopen(shpbuf,OPEN_RB))) {   
        printf1(ctx, "Error: can't open the shp file.\n");         
        goto SDSHPFin;
    }
    shpfd_open = 1; 

    strcpy(shpbuf,ctx->PMRHSTR);                         /* open dbf file */
    strcat(shpbuf,".dbf");
    if (!(dbffd = fopen(shpbuf,OPEN_RB))) {   
        printf1(ctx, "Error: can't open the dbf file.\n");         
        goto SDSHPFin;
    }
    dbffd_open = 1;
                        
    printf1(ctx, "\nReading the index file.\n");

    if (fread(shpbuf,sizeof(char),HLEN,shxfd) != HLEN) {
        p_err(ctx, -7,1);
        goto SDSHPFin;
    }
    ctx->HILO = 1;
    n = st_geti(ctx, shpbuf);
    if (n != 9994) {
        printf1(ctx, "Probably not a proper index/shape file.\n");
        goto SDSHPFin;
    }
    flen = 2 * st_geti(ctx, shpbuf + 24);
    flen -= HLEN;
    ns = flen / 8;                  /* number of shapes */ 
 
    ctx->HILO = 2;
    n = st_geti(ctx, shpbuf + 28);
    printf1(ctx, "Version: %d\n",n);
    typ = st_geti(ctx, shpbuf + 32);
    if (shp_typ(ctx, typ,1))             /* check shape type */
        goto SDSHPFin;

    tda_out("Number of shapes: %d\n",ns);

    xmin = st_getd(ctx, shpbuf + 36);
    ymin = st_getd(ctx, shpbuf + 44);
    xmax = st_getd(ctx, shpbuf + 52);
    ymax = st_getd(ctx, shpbuf + 60);

    printf1(ctx, "\nBounding box X: %20.8lf  %20.8lf\n",xmin,xmax);
    printf1(ctx, "Bounding box Y: %20.8lf  %20.8lf\n",ymin,ymax);

    printf1(ctx, "\nReading the dbf file.\n");

    if (shp_dbfh(ctx, 1,dbffd,shpbuf,&nrec,&hlen,&nvar,&rlen,ns))  /* read dbf file header */
        goto SDSHPFin;

    n = imax(ctx, rlen,DBFRLEN);
    if (!(dbfbuf = (char *)calloc((size_t)(n + 2),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto SDSHPFin;
    }
    dbfbufa = n + 2;
    memrq(ctx, dbfbufa,sizeof(char));
                  
    if (shp_dbf_a(ctx, nvar,1))              /* allocate memory for fields */
        goto SDSHPFin;
         
    if (shp_dbff(ctx, dbffd,dbfbuf,nvar))    /* read fields */
        goto SDSHPFin;

    if (prn_dbf_var(ctx, nvar))              /* info about dbf variables */
        goto SDSHPFin;

    if (ctx->PMF1Def == 0) {                 /* without output file end here */
        err = 0;
        goto SDSHPFin;
    }
    printf1(ctx, "Creating the output data file: %s\n",ctx->PMF1dName);
    drec1 = drec = n1 = n0 = 0;

    /* skip rest of dbf header */

    if ((hlen -= (nvar + 1) * DBFRLEN) < 1) {
        printf1(ctx, "Error in dbf header length.\n");
        goto SDSHPFin;
    }
    if (fseek(dbffd,(long)hlen,SEEK_CUR)) {
        printf1(ctx, "Error: cannot seek to begin of fields in dbf file.\n");
        goto SDSHPFin;
    }
    id = id1 = 0;
    for (i = 0; i < ns; ++i) {

        /* read current record from the dbf file */

        if (fread(dbfbuf,sizeof(char),(size_t)(rlen),dbffd) != (size_t)(rlen)) {
            p_err(ctx, -7,1);
            goto SDSHPFin;
        }

        /* read the index file */

        if (fread(shpbuf,sizeof(char),8,shxfd) != 8) {
            p_err(ctx, -7,1);
            goto SDSHPFin;
        }
        ctx->HILO = 1;
        n = st_geti(ctx, shpbuf);        /* n = pointer to shape file */
        len = st_geti(ctx, shpbuf + 4);  /* len = size of entries */ 

        /* read the shape file */

        if (fseek(shpfd,(long)(2 * n),0)) {
            printf1(ctx, "Can't seek to offset %d in shape file.\n",2 * n);
            goto SDSHPFin;
        }
        if (fread(shpbuf,sizeof(char),8,shpfd) != 8) {
            p_err(ctx, -7,1);
            goto SDSHPFin;
        }
        m = st_geti(ctx, shpbuf + 4);
        if (m != len) {
            printf1(ctx, "Error: inconsistent information in index and shape file.\n");
            goto SDSHPFin;
        }
        len *= 2;
        if (len > shpbufl) {            /* adjust buffer length */
            free(shpbuf);
            memrq(ctx, -shpbufa,sizeof(char));

            shpbufl = len;
            if (!(shpbuf = (char *)calloc((size_t)(shpbufl + 2),sizeof(char)))) {
                p_err(ctx, -2,1);
                goto SDSHPFin;
            }
            shpbufa = shpbufl + 2;
            memrq(ctx, shpbufa,sizeof(char));
        }
        if (fread(shpbuf,sizeof(char),(size_t)(len),shpfd) != (size_t)(len)) {
            p_err(ctx, -7,1);
            goto SDSHPFin;
        }
        ctx->HILO = 2;
        typ = st_geti(ctx, shpbuf);            /* shape type */
        id1++;
        switch (typ) {
            case 0:                                 /* null shape */
                n0++;
                break;   

            case 1:                                 /* point */
                id++;
                id2 = 1;
                np = mn = 1;
                fprintf(ctx->PMF1d,"%10d %1d %7d %10d %4d ",id,mn,np,id1,id2);
                shp_dbfw(ctx, dbfbuf,rlen);
                fprintf(ctx->PMF1d,"\n");
                drec++;

                x = st_getd(ctx, shpbuf + 4);
                y = st_getd(ctx, shpbuf + 12);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,x);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,y);
                fprintf(ctx->PMF1d,"\n");
                drec1++;
                break;

            case 3:                                 /* polyline */
                mn = 2;
                TDA_FALLTHROUGH;

            case 5:                                 /* polygon */
                if (typ == 5)
                    mn = 3;

                npart = st_geti(ctx, shpbuf + 36);
                npnt  = st_geti(ctx, shpbuf + 40);

                if (npart > 1)
                    n1++;

                p = shpbuf + 44;
                m = st_geti(ctx, p);
                p += 4;

                id2 = 0;
                for (j = 1; j <= npart; ++j) {
                    id++;
                    id2++;

                    if (j == npart)
                        m1 = npnt;
                    else {
                        m1 = st_geti(ctx, p);
                        p += 4;
                    }
                    np = m1 - m; 

                    fprintf(ctx->PMF1d,"%10d %1d %7d %10d %4d ",id,mn,np,id1,id2);
                    shp_dbfw(ctx, dbfbuf,rlen);
                    fprintf(ctx->PMF1d,"\n");
                    drec++;

                    q = shpbuf + 44 + 4 * npart + m * 16;

                    for (k = 0; k < np; ++k) {
                        x = st_getd(ctx, q);
                        q += 8;
                        y = st_getd(ctx, q);
                        q += 8;
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,x);
                        rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,y);
                        fprintf(ctx->PMF1d,"\n");
                        drec1++;
                    }
                    m = m1;
                }
                break;

            default:
                printf1(ctx, "Error: found unknown shape type %d.\n",typ);
                goto SDSHPFin;
        }
        prn_message(ctx, i,0,0);
    }
    prn_message(ctx, i,1,0);
         
    printf1(ctx, "Object-description records: %d. Data records: %d\n",drec,drec1);
    printf1(ctx, "Total number of records written to output file: %d\n",drec + drec1);
    if (n0 > 0)
        printf1(ctx, "Note: ignored %d records containing null type shapes.\n",n0);
    if (n1 > 0)
        printf1(ctx, "Note: %d shapes contain more than one part.\n",n1);

    if (ctx->PMTDAFDef) {                /* create TDA decription file */

        fprintf(ctx->PMTDAFd,"sdnvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMF1dName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",drec);
        fprintf(ctx->PMTDAFd,"  ffmt=c1(1-10),c2(12),c3(14-20),c4(22-31),c5(33-36),");
            
        n = 5; j = 6; l = 38;
        for (i = 0; i < nvar; ++i) {
            k = l + ctx->DBFFL[i] - 1;
            if (ctx->DBFFT[i] == 0) {
                if (++n > 5) {
                    fprintf(ctx->PMTDAFd,"\n       ");
                    n = 1;
                }
                fprintf(ctx->PMTDAFd,"c%d(%d-%d),",j,l,k);
                j++;
            }
            l = k + 2;
        }
        fprintf(ctx->PMTDAFd,"\n\n  SDID<5>[10.0] = c1,\n");
        fprintf(ctx->PMTDAFd,"  SDTyp<1>[1.0] = c2,\n");
        fprintf(ctx->PMTDAFd,"  SDN<5>[7.0] = c3,\n");
        fprintf(ctx->PMTDAFd,"  SDPtr<5>[10.0] = rd,\n");
        fprintf(ctx->PMTDAFd,"  SDID1<5>[10.0] = c4,\n");
        fprintf(ctx->PMTDAFd,"  SDID2<5>[ 4.0] = c5,\n");

        j = 6; l = 38;
        for (i = 0; i < nvar; ++i) {
            k = l + ctx->DBFFL[i] - 1;
            if (ctx->DBFFT[i] == 0) {
                if (ctx->DBFFD[i] > 0)
                    n = 8;
                else
                    n = 5;

                fprintf(ctx->PMTDAFd,"  %s<%d>[%d.%d] = c%d,\n",ctx->DBFFN + i * DBFFLEN,
                                            n,ctx->DBFFL[i],ctx->DBFFD[i],j);
                j++;
            }
            else if (ctx->DBFFT[i] == 1) {
                fprintf(ctx->PMTDAFd,"  %s = str(%d,%d),\n",ctx->DBFFN + i * DBFFLEN,l,k);
            }
            l = k + 2;
        }
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;
    newline(ctx);

SDSHPFin:
    if (dbfbufa > 0) {
        free(dbfbuf);
        memrq(ctx, -dbfbufa,sizeof(char));
    }
    if (shpbufa > 0) {
        free(shpbuf);
        memrq(ctx, -shpbufa,sizeof(char));
    }

    shp_dbf_a(ctx, nvar,0);

    if (dbffd_open)  
        fclose(dbffd);
    if (shpfd_open)  
        fclose(shpfd);
    if (shxfd_open)  
        fclose(shxfd);

    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  shp_typ(n,opt)  Check shape type. Return 0 if OK, -1 if not defined.    */
/*                  If opt != 0 print.                                      */

int shp_typ(TDAContext *ctx, int n,int opt)
{
    if (opt)
        printf1(ctx, "Shape type: ");

    switch (n) {
        case 1:     if (opt) printf1(ctx, "Point");
                    break;
        case 3:     if (opt) printf1(ctx, "PolyLine");
                    break;
        case 5:     if (opt) printf1(ctx, "Polygon");
                    break;
        /*****************************************
        case 8:     if (opt) printf1(ctx, "MultiPoint");
                    break;
        case 11:    if (opt) printf1(ctx, "PointZ");
                    break;
        case 13:    if (opt) printf1(ctx, "PolyLineZ");
                    break;
        case 15:    if (opt) printf1(ctx, "PolygonZ");
                    break;
        case 18:    if (opt) printf1(ctx, "MultiPointZ");
                    break;
        case 21:    if (opt) printf1(ctx, "PointM");
                    break;
        case 23:    if (opt) printf1(ctx, "PolyLineM");
                    break;
        case 25:    if (opt) printf1(ctx, "PolygonM");
                    break;
        case 28:    if (opt) printf1(ctx, "MultiPointM");
                    break;
        case 31:    if (opt) printf1(ctx, "MultiPatch");
                    break;
        ******************************************/
        default:    if (opt) printf1(ctx, "%d (unknown)\n",n);
                    return(-1);
    }
    if (opt)
        printf1(ctx, " (%d)\n",n);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbfw(buf,rlen)  Write current dbf record to output file.            */

void shp_dbfw(TDAContext *ctx, char *buf,int rlen)
{
    register int j,k,l;
    char *p;

    p = buf + 1;
    k = l = 0;
    for (j = 1; j < rlen; ++j) {
        if (*p == '\0')
            *p = ' ';
        fprintf(ctx->PMF1d,"%c",*p++);
        if (++l >= ctx->DBFFL[k]) {
            fprintf(ctx->PMF1d," ");
            k++;
            l = 0;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  shp_dbfh(opt,fd,buf,nrec,hlen,nvar,rlen,ns)                             */
/*                                                                          */
/*  Read dbf file header. If opt != 0, check for consistency with shp file. */
/*  ns is the number of shapes found in the shape file.                     */
/*  Return 0 if OK, -1 if error.                                            */

int shp_dbfh(TDAContext *ctx, int opt,FILE *fd,char *buf,int *nrec,int *hlen,int *nvar,int *rlen,int ns)
{
    if (fread(buf,sizeof(char),DBFRLEN,fd) != DBFRLEN) {
        printf1(ctx, "Error: can't read header of dbf file.\n");
        return(-1);
    }
    ctx->HILO = 2;
    *nrec = st_geti(ctx, buf + 4);
    *hlen = st_gets(ctx, buf + 8);
    *rlen = st_gets(ctx, buf + 10);
   
    if (opt && *nrec != ns) {
        printf1(ctx, "Error: number of records in dbf file is %d.\n",*nrec);
        printf1(ctx, "Inconsistent with information in index/shape file.\n");
        return(-1);
    }
    tda_out("Number of records: %d\n",*nrec);
    tda_out("Record length: %d\n",*rlen);
    if (*rlen < 1) {
        printf1(ctx, "This is probably not a dbf file.\n");
        return(-1);
    }
    *nvar = *hlen - 32;
    *nvar /= 32;

    tda_out("Number of fields: %d\n",*nvar);
    if (*nvar < 1) {
        printf1(ctx, "This is probably not a dbf file.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbff(nvar)        Read dbf fields and save information in gloabal   */
/*                        arrays DBFFN, DBFFT, DBFFL, DBFFD, DBFFP.         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int shp_dbff(TDAContext *ctx, FILE *fd,char *buf,int nvar)
{
    register int i;
    char *vtyp; 
    unsigned char u;

    for (i = 0; i < nvar; ++i) {
        if (fread(buf,sizeof(char),DBFRLEN,fd) != DBFRLEN) {
            printf1(ctx, "Error: can't read fields of dbf file.\n");
            return(-1);
        }
        strncpy(ctx->DBFFN + i * DBFFLEN,buf,DBFFLEN);
        vtyp = buf + DBFFLEN;
        if (*vtyp == 'N' || *vtyp == 'F') {
            ctx->DBFFT[i] = 0;
            ctx->DBFFD[i] = (int)buf[17];
        }
        else if (*vtyp == 'C')
            ctx->DBFFT[i] = 1;
        else
            ctx->DBFFT[i] = 2;

        u = (unsigned char)buf[16];
        ctx->DBFFL[i] = (int)u;        

        if (i == 0)
            ctx->DBFFP[i] = 1;
        else
            ctx->DBFFP[i] = ctx->DBFFL[i - 1] + ctx->DBFFP[i - 1];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbf_a(nvar,opt)   If opt != 0 allocate arrays for n dbf fields,     */
/*                        otherwise free previously allocated memory.       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int shp_dbf_a(TDAContext *ctx, int nvar,int opt)
{

    if (ctx->DBFFN_A > 0) {
        free((char *)ctx->DBFFN);
        memrq(ctx, -ctx->DBFFN_A,sizeof(char));
        ctx->DBFFN_A = 0;
    }
    if (ctx->DBFFT_A > 0) {
        free((char *)ctx->DBFFT);
        memrq(ctx, -ctx->DBFFT_A,sizeof(int));
        ctx->DBFFT_A = 0;
    }
    if (ctx->DBFFL_A > 0) {
        free((char *)ctx->DBFFL);
        memrq(ctx, -ctx->DBFFL_A,sizeof(int));
        ctx->DBFFL_A = 0;
    }
    if (ctx->DBFFD_A > 0) {
        free((char *)ctx->DBFFD);
        memrq(ctx, -ctx->DBFFD_A,sizeof(int));
        ctx->DBFFD_A = 0;
    }
    if (ctx->DBFFP_A > 0) {
        free((char *)ctx->DBFFP);
        memrq(ctx, -ctx->DBFFP_A,sizeof(int));
        ctx->DBFFP_A = 0;
    }
    if (opt == 0)
        return(0);

    if (!(ctx->DBFFN = (char *)calloc((size_t)(nvar) * (size_t)((DBFFLEN + 1)),sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    ctx->DBFFN_A = nvar * (DBFFLEN + 1);
    memrq(ctx, ctx->DBFFN_A,sizeof(char));

    if (!(ctx->DBFFT = (int *)calloc((size_t)(nvar),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    ctx->DBFFT_A = nvar;
    memrq(ctx, ctx->DBFFT_A,sizeof(int));

    if (!(ctx->DBFFL = (int *)calloc((size_t)(nvar),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    ctx->DBFFL_A = nvar;
    memrq(ctx, ctx->DBFFL_A,sizeof(int));

    if (!(ctx->DBFFD = (int *)calloc((size_t)(nvar),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    ctx->DBFFD_A = nvar;
    memrq(ctx, ctx->DBFFD_A,sizeof(int));

    if (!(ctx->DBFFP = (int *)calloc((size_t)(nvar),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    ctx->DBFFP_A = nvar;
    memrq(ctx, ctx->DBFFP_A,sizeof(int));
    return(0);
}


/* ------------------------------------------------------------------------ */
/*  rdbf            Read a dBase file.                                      */
/*                                                                          */
/*                  rdbf(                                                   */
/*                      df=...,     output file                             */
/*                                                                          */
/*                  ) = fname;      dbf input file.                         */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rdbf(TDAContext *ctx)  
{
    int err,nrec,nvar,rlen,hlen,i,j,k,l,wrec;
    char *p,buf[DBFRLEN + 2];

    err = -1;
    nvar = 0;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Reading a dBase file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {     /* get parameters */
        goto RDBFin;
    }   
    printf1(ctx, "File name: %s\n\n",ctx->PMFdName);

    if (shp_dbfh(ctx, 0,ctx->PMFd,buf,&nrec,&hlen,&nvar,&rlen,0))    /* read dbf file header */
        goto RDBFin;

    hlen -= (nvar + 1) * DBFRLEN;
    if (hlen < 1) {
        printf1(ctx, "Error: in dbf header length.\n");
        goto RDBFin;
    }
    if (shp_dbf_a(ctx, nvar,1))          /* allocate memory for fields */
        goto RDBFin;

    if (shp_dbff(ctx, ctx->PMFd,buf,nvar))    /* read fields */
        goto RDBFin;

    if (prn_dbf_var(ctx, nvar))          /* info about dbf variables */
        goto RDBFin;

    if (nrec > 0 && ctx->PMF1Def) {      /* read data and print to output file */

        wrec = 0;
        if (alloc_acc(ctx, rlen + 1))
            goto RDBFin;

        /* skip rest of dbf header */

        if (fseek(ctx->PMFd,(long)hlen,SEEK_CUR)) {
            printf1(ctx, "Error: cannot seek to begin of fields in dbf file.\n");
            goto RDBFin;
        }
        for (i = 0; i < nrec; ++i) {

            if (fread(ctx->AcC,sizeof(char),(size_t)(rlen),ctx->PMFd) != (size_t)(rlen)) {
                p_err(ctx, -7,1);
                goto RDBFin;
            }
            p = ctx->AcC + 1;
            k = l = 0;
#ifdef TDA_R_PACKAGE
            {
                /*  Each field of the record as its own string, so the R
                    side gets the attribute table without re-splitting
                    the fixed-width file by column position.  A dbf
                    holds character fields ('C', type 1 above) as well
                    as numbers, so the values travel as text and are
                    converted per field on the R side.  */
                char fbuf[256];
                int fi, fj, fp = 1;
                for (fi = 0; fi < nvar; ++fi) {
                    int fl = ctx->DBFFL[fi];
                    if (fl > (int)sizeof(fbuf) - 1)
                        fl = (int)sizeof(fbuf) - 1;
                    for (fj = 0; fj < fl && fp + fj < rlen; ++fj) {
                        char c = ctx->AcC[fp + fj];
                        fbuf[fj] = (c == '\0') ? ' ' : c;
                    }
                    fbuf[fj] = '\0';
                    tda_export_str_row(ctx, "rdbf.values", fbuf);
                    fp += ctx->DBFFL[fi];
                }
            }
#endif
            for (j = 1; j < rlen; ++j) {
                if (*p == '\0')
                    *p = ' ';
                fprintf(ctx->PMF1d,"%c",*p++);
                if (++l >= ctx->DBFFL[k]) {
                    fprintf(ctx->PMF1d," ");
                    k++;
                    l = 0;
                }
            }
            fprintf(ctx->PMF1d,"\n");
            wrec++;
        }
#ifdef TDA_R_PACKAGE
        /* field names and types, so the R side need not read the
           printed variable table either */
        {
            int fi;
            char nbuf[80];
            for (fi = 0; fi < nvar; ++fi) {
                snprintf(nbuf, sizeof(nbuf), "%.*s|%d",
                         DBFFLEN, ctx->DBFFN + fi * DBFFLEN,
                         ctx->DBFFT[fi]);
                tda_export_str_row(ctx, "rdbf.fields", nbuf);
            }
        }
#endif
        printf1(ctx, "%d records written to: %s\n",wrec,ctx->PMF1dName);
    }
    err = 0;     

RDBFin:
    shp_dbf_a(ctx, nvar,0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_dbf_var(nvar)       Print table with dbf variables.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int prn_dbf_var(TDAContext *ctx, int nvar)
{
    int i,n;

    printf1(ctx, "\nIdx  Type  Length  D.Count  Name\n");
    prnchar(ctx, '-',32,1);

    n = 0;                                  /* save var definitions */
    for (i = 0; i < nvar; ++i) {
        if (ctx->DBFFL[i] < 1) {
            printf1(ctx, "Error: in dbf field definitions.\n");
            return(-1);   
        }
        printf1(ctx, "%3d %4d   %5d  %7d   %s\n",i + 1,ctx->DBFFT[i],ctx->DBFFL[i],
                            ctx->DBFFD[i],ctx->DBFFN + i * DBFFLEN);  

        if (ctx->DBFFT[i] != 0 && ctx->DBFFT[i] != 1)
            n++;
    }
    newline(ctx);
    if (n > 0) {
        printf1(ctx, "Warning: the dbf file contains %d fields which are\n",n);
        printf1(ctx, "neither numerical nor character. Will be ignored.\n\n");
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rcsv            Read a CSV file.                                        */
/*                                                                          */
/*                  rcsv(                                                   */
/*                      df=...,     output file                             */
/*                                                                          */
/*                  ) = fname;      csv input file.                         */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rcsv(TDAContext *ctx)  
{
    int err,nrec,c,ns,nsmax,j,n;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Reading a csv file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {     /* get parameters */
        goto RCSVFin;
    }   
    printf1(ctx, "File name: %s\n\n",ctx->PMFdName);

    nsmax = ns = nrec = 0; 
    {
    int rlen = 0, rlenmax = 0;
    while ((c = fgetc(ctx->PMFd)) != EOF) {
        if (c == '\r')
            continue;
        if (c == ';')
            ns++;
        else if (c == '\n') {
            nrec++;
            nsmax = imax(ctx, nsmax,ns);
            ns = 0;
            rlenmax = imax(ctx, rlenmax,rlen);
            rlen = 0;
        }
        else
            rlen++;
        rlen += (c == ';');
    }
#ifdef TDA_R_PACKAGE
    ctx->RcsvLineMax = rlenmax + ns + rlen + 2;
#endif
    }
    tda_out("Number of records: %d\n",nrec);
    tda_out("Maximal number of entries per line: %d\n",nsmax);
    if (nsmax < 1)  
        goto RCSVFin;
                    
    if (alloc_acn(ctx, nsmax + 1))
        goto RCSVFin;

    if (fseek(ctx->PMFd,(long)0,0)) {             
        printf1(ctx, "Error: cannot seek to begin of file.\n");
        goto RCSVFin;
    }
    n = j = 0; 
    while ((c = fgetc(ctx->PMFd)) != EOF) {
        if (c == '\r')
            continue;
        if (c == ';') {
            if (j < nsmax)  
                ctx->AcN[j] = imax(ctx, ctx->AcN[j],n);
            j++;
            n = 0;
        }             
        else if (c == '\n')  
            j = 0;
        else
            n++;
    }
    if (fseek(ctx->PMFd,(long)0,0)) {             
        printf1(ctx, "Error: cannot seek to begin of file.\n");
        goto RCSVFin;
    }
    n = j = 0; 
    {
#ifdef TDA_R_PACKAGE
    /* the padded rendering below carries no separator guarantee, so
       the raw semicolon-bounded lines are exported alongside it and
       the R side splits them unambiguously */
    char *rline = (char *)calloc(ctx->RcsvLineMax > 0 ?
                                 (size_t)ctx->RcsvLineMax : 2, 1);
    int rpos = 0;
#endif
    while ((c = fgetc(ctx->PMFd)) != EOF) {
        if (c == '\r')
            continue;
#ifdef TDA_R_PACKAGE
        if (rline != NULL) {
            if (c == '\n') {
                rline[rpos] = '\0';
                tda_export_str_row(ctx, "rcsv.lines", rline);
                rpos = 0;
            }
            else if (rpos + 1 < ctx->RcsvLineMax)
                rline[rpos++] = (char)c;
        }
#endif
        if (c == ';') {
            if (j < nsmax) {
                while (n++ < ctx->AcN[j])
                    tda_out("%c",' ');
            }
            j++;
            n = 0;
        }             
        else if (c == '\n') {
            tda_out("%c",c);
            j = 0;
        }
        else {
            tda_out("%c",c);
            n++;
        }
    }
#ifdef TDA_R_PACKAGE
    free(rline);
#endif
    }
    /*  err was set to -1 at the top and never cleared, so rcsv reported
        failure on every run, including the ones that converted the file
        correctly.  Reaching here means the whole input was read.  */
    err = 0;

RCSVFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  rplz            Read special version of a CSV file (PLZ data)           */
/*                                                                          */
/*                  rplz(                                                   */
/*                      df=...,     output file                             */
/*                                                                          */
/*                  ) = fname;      csv input file.                         */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rplz(TDAContext *ctx)  
{
    int err,nrec,c,n,ns,rl,rlmax,id,nt,nl;
    char *p,*q;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Reading a plz file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,10,1)) {     /* get parameters */
        goto RPLZFin;
    }   
    printf1(ctx, "File name: %s\n\n",ctx->PMFdName);

    rlmax = rl = nrec = 0; 
    while ((c = fgetc(ctx->PMFd)) != EOF) {
        rl++;
        if (c == '\n') {
            nrec++;
            rlmax = imax(ctx, rlmax,rl);
            rl = 0;
        }
    }
    printf1(ctx, "Number of records: %d\n",nrec);
    printf1(ctx, "Maximal record length: %d\n",rlmax);
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Need an output file.\n");
        err = 0;
        goto RPLZFin;
    }
    if (fseek(ctx->PMFd,(long)0,0)) {             
        printf1(ctx, "Error: cannot seek to begin of file.\n");
        goto RPLZFin;
    }
    if (alloc_acc(ctx, rlmax + 10))
        goto RPLZFin;

    id = nrec = 0;
    while (fgets(ctx->AcC,rlmax,ctx->PMFd)) {
        nrec++;
        if (nrec == 1)
            continue;

        ns = 0; 
        p = ctx->AcC;
        while (*p) {
            if (*p == ';') {
                ns++;
                *p = ' ';
            }
            else if (*p == ',')
                *p = ' ';
            if (ns == 3) {
                *p = '\0';
                q = p + 1;
                p = q;
                n = 0;
                while (*p) {
                    if (*p == ',')
                        n++; 
                    p++;
                }
                nt = n;
                if (nt > 3)
                    nt = 3;

                id++;
                fprintf(ctx->PMF1d,"%12d %d %8d %s\n",id,nt,n,ctx->AcC);

                p = q;
                nl = 0;
                while (*p) {
                    if (*p == ',') {
                        fprintf(ctx->PMF1d,"  ");
                        nl = 1;
                        p++;
                    }
                    else if (*p == ' ') {
                        fprintf(ctx->PMF1d,"\n");
                        nl = 0;
                        p = skip_b(ctx, p);
                    }
                    else if (*p == ':') {
                        fprintf(ctx->PMF1d,"\n");
                        nl = 0;
                        p++;
                    }
                    else if (*p == '\n') {
                        fprintf(ctx->PMF1d,"\n");
                        nl = 0;
                        break;            
                    }
                    else {
                        fprintf(ctx->PMF1d,"%c",*p);
                        nl = 1;
                        p++;
                    }
                }
                if (nl)
                    fprintf(ctx->PMF1d,"\n");
                break;
            }
            p++; 
        }
        /************ 
        if (id > 10)
        break;
        ************/
    }
    printf1(ctx, "%d records written to: %s\n",id,ctx->PMF1dName);
    err = 0;

RPLZFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  gtopo   Read GTOPO30 elevation files.                                   */
/*                                                                          */
/*          gtopo(                                                          */
/*              rows=...,   number of rows                                  */
/*              cols=...,   number of columns                               */  
/*              ulx=...,    x coordinate of upper left corner               */
/*              uly=...,    y coordinate of upper left corner               */
/*              dx=...,     x dimension of a pixel                          */
/*              dy=...,     y dimension of a pixel                          */
/*                                                                          */
/*              df=...,     output file                                     */
/*              fmt=...,    print format, def. 10.6                         */
/*              opt=...,    option, def. 1 (see below)                      */
/*              lon=...,    select range of longitudes                      */ 
/*              lat=...,    select range of latitudes                       */
/*                                                                          */
/*          ) = fname;      path to binary GTOPO30 file.                    */
/*                                                                          */
/*                                                                          */
/*  Note: ulx, uly, dx and dy must be given in decimal degrees.             */
/*                                                                          */
/*  Options                                                                 */
/*                                                                          */
/*  Option 1. In this case the user must provide a specific longitude       */
/*      and a range of latitudes with the parameters                        */
/*                                                                          */
/*          lon = lon_value,                                                */
/*          lat = lata_value,latb_value,                                    */
/*                                                                          */
/*      or alternatively a specific latitude and a range of longitudes      */
/*      with the parameters                                                 */
/*                                                                          */
/*          lat = lat_value,                                                */
/*          lon = lona_value,lonb_value,                                    */
/*                                                                          */
/*      For the given range, the output file will contain records each      */
/*      having three columns: (1) longitude, (2) latitude, (3) elevation    */
/*      value.                                                              */


/*  Return 0 if OK, -1 if error.                                            */

int gtopo(TDAContext *ctx)  
{
    register int i,j; 
    int err,err1,nrec,v,r,i0,i1,j0,j1;                 
    double dxa,dxb,dya,dyb,lona,lonb,lata,latb;

    err = -1;
    err1 = nrec = 0;

    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Reading a GTOPO30 elevation file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,10,1)) {     /* get parameters */
        goto GTOPOFin;
    }   
    if (ctx->PMOPT != 1)
        ctx->PMOPT = 1;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,6);

    printf1(ctx, "File name: %s\n\n",ctx->PMFdName);
    printf1(ctx, "Number of rows: %d\n",ctx->PMRows);
    printf1(ctx, "Number of columns: %d\n",ctx->PMCols);
    if (ctx->PMRows < 1 || ctx->PMCols < 1) {
        printf1(ctx, "Error: need positive values.\n");
        goto GTOPOFin;
    }
    printf1(ctx, "Coordinates: ulx = %18.14lf uly=%18.14lf\n",ctx->PMULX,ctx->PMULY);
    printf1(ctx, "Grid spacing: dx = %18.14lf  dy=%18.14lf\n",ctx->PMDX,ctx->PMDY);
    if (ctx->PMDX < ctx->EPSI1 || ctx->PMDY < ctx->EPSI1) {
        printf1(ctx, "Error: dx and dy must be strictly positive.\n");
        goto GTOPOFin;
    }
    dxa = ctx->PMULX;
    dyb = ctx->PMULY;
    dxb = dxa + (ctx->PMCols - 1) * ctx->PMDX;
    dya = dyb - (ctx->PMRows - 1) * ctx->PMDY;

    printf1(ctx, "\nRange of the grid in decimal degrees.\n");
    printf1(ctx, "Longitude: %18.14lf %18.14lf\n",dxa,dxb);
    printf1(ctx, "Latitude:  %18.14lf %18.14lf\n\n",dya,dyb);

    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Without an output file: simply reading the file.\n");
        for (i = 0; i < ctx->PMRows; i++) {
            for (j = 0; j < ctx->PMCols; ++j) {
                v = gtopo_rpixel(ctx, i,j,ctx->PMCols,&r);
                if (r) {
                    printf1(ctx, "Error: cannot read entry in row %d, column %d.\n",i + 1,j + 1);
                    goto GTOPOFin;
                }
            }
            prn_message(ctx, i + 1,0,0);
        }
        prn_message(ctx, ctx->PMRows,1,0);
        printf1(ctx, "No errors occurred.\n");
        err = 0;
        goto GTOPOFin;
    }
    printf1(ctx, "Output file: %s\n",ctx->PMF1dName);
    printf1(ctx, "Option: %d\n\n",ctx->PMOPT);

    err1 = 1;
    if (ctx->PMNTP < 1 || ctx->PMNTP > 2 || ctx->PMNTP1 < 1 || ctx->PMNTP1 > 2)  
        goto GTOPOFin;
      
    if ((ctx->PMNTP == 2 && ctx->PMTP[0] > ctx->PMTP[1]) || (ctx->PMNTP1 == 2 && ctx->PMTP1[0] > ctx->PMTP1[1]))  
        goto GTOPOFin;
        
    for (i = 0; i < ctx->PMNTP; ++i) {
        if (ctx->PMTP[i] < -180.0 || ctx->PMTP[i] > 180.0)  
            goto GTOPOFin;
    }
    for (j = 0; j < ctx->PMNTP1; ++j) {
        if (ctx->PMTP1[j] < -90.0 || ctx->PMTP1[j] > 90.0)  
            goto GTOPOFin;
    }
    err1 = 0;

    if (ctx->PMOPT == 1) {
        if (ctx->PMNTP == 1 && ctx->PMNTP1 == 2) {
            lona = lonb = ctx->PMTP[0];
            lata = ctx->PMTP1[0];
            latb = ctx->PMTP1[1];
            if (lona < dxa || lonb > dxb || lata > dyb || latb < dya) {
                err1 = 2;
                goto GTOPOFin;
            }
            j0 = j1 = (int)((lona - dxa) / ctx->PMDX + 0.5);
            lona = lonb = dxa + j0 * ctx->PMDX;

            if (lata < dya)  
                i0 = ctx->PMRows - 1;
            else  
                i0 = (int)((dyb - lata) / ctx->PMDY + 0.5);
            lata = dyb - i0 * ctx->PMDY;

            if (latb > dyb)  
                i1 = 0;             
            else  
                i1 = (int)((dyb - latb) / ctx->PMDY + 0.5);
            latb = dyb - i1 * ctx->PMDY;
        }
        else if (ctx->PMNTP == 2 && ctx->PMNTP1 == 1) {
            lona = ctx->PMTP[0];
            lonb = ctx->PMTP[1];
            lata = latb = ctx->PMTP1[0];

            if (lata < dya || latb > dyb || lona > dxb || lonb < dxa) {
                err1 = 2;
                goto GTOPOFin;
            }
            i0 = i1 = (int)((lata - dya) / ctx->PMDY + 0.5);
            lata = latb = dya + i0 * ctx->PMDY;

            if (lona < dxa)  
                j0 = 0;
            else  
                j0 = (int)((lona - dxa) / ctx->PMDX + 0.5);
            lona = dxa + j0 * ctx->PMDX;

            if (lonb > dxb)  
                j1 = ctx->PMCols - 1;
            else  
                j1 = (int)((lonb - dxa) / ctx->PMDX + 0.5);
            lonb = dxa + j1 * ctx->PMDX;

        }
        else {   
            printf1(ctx, "Syntax error in lon and/or lat parameters.\n");
            goto GTOPOFin;
        }
        if (i0 < i1 || j0 > j1) {
            err1 = 2;
            goto GTOPOFin;
        }
        printf1(ctx, "Will use the following range of coordinates.\n");
        printf1(ctx, "Longitude         Pixel  Latitude          Pixel\n");
        printf1(ctx, "%16.12lf %6d  %16.12lf %6d\n",lona,j0,lata,i0);
        printf1(ctx, "%16.12lf %6d  %16.12lf %6d\n\n",lonb,j1,latb,i1);
                
        for (i = i0; i >= i1; i--) {
            for (j = j0; j <= j1; ++j) {
                v = gtopo_rpixel(ctx, i,j,ctx->PMCols,&r);
                if (r) {
                    printf1(ctx, "Error: cannot read entry in row %d, column %d.\n",i + 1,j + 1);
                    goto GTOPOFin;
                }
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dxa + j * ctx->PMDX);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,dyb - i * ctx->PMDY);
                fprintf(ctx->PMF1d,"%8d\n",v);
                nrec++;

            }
        }
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    }
    err = 0;     

GTOPOFin:
    if (err1 == 1)  
        printf1(ctx, "Invalid values for lon or lat parameter.\n");
    else if (err1 == 2)  
        printf1(ctx, "Selected lon/lat range is empty or outside of the grid region.\n");

    p_clean(ctx);
    return(err);
}


/* -##--------------------------------------------------------------------- */
/*  gtopo_rpixel(i,j,nc,*err)                                               */
/*                                                                          */
/*  Return value of Pixel in row i, column j. nc is the number of columns   */
/*  in the file.  Return err = 0 if OK, -1 if an error occured.             */  
/*                                                                          */
/*  Note: row and column counting begins with 0.                            */

int gtopo_rpixel(TDAContext *ctx, int i,int j,int nc,int *err)
{
    int n;
    char buf[3];

    *err = -1;
    n = 2 * (i * nc + j);

    if (n < 0 || fseek(ctx->PMFd,(long)n,0))  
        return(0.0);
                                        
    if (fread(buf,sizeof(char),2,ctx->PMFd) != 2)  
        return(0.0);
                                               
    ctx->HILO = 1;
    n = st_gets(ctx, buf);
    *err = 0;
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  sddcwp  Conversion of DCW polygon point data file.                      */
/*                                                                          */
/*          sddcwp(                                                         */
/*              df=...,     output file (required)                          */
/*              nc=...,     if nc=1 do not include centroids, def. nc=0     */
/*              nfmt=...,   integer print format, def. 4                    */
/*              fmt=...,    floating point print format, def. 12.6          */
/*                                                                          */
/*          ) = input_file;                                                 */
/*                                                                          */
/*  This command is intended to transform a Polygon Point ASCII file from   */
/*  the DCW (Digital Chart of the Worl) project into a TDA spatial data     */
/*  file. The input file must conform to the file format provided by        */
/*  www.maproom.psu.edu/dcw/  (option: DCW Polygon Point Generation).       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sddcwp(TDAContext *ctx)  
{
    register int i; 
    int err,rec,nrec,blen,m,ii,np,npmax,npt,nobj,nobje,id,nt1,nt3;
    char *p,name[200];
    double x,y;

    blen = 1000;                /* read buffer length */
    err = -1;

    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Conversion of DCW polygon point data file. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,10,1)) {     /* get parameters */
        goto SDDCWPFin;
    }   
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 12,6);

    printf1(ctx, "Input file: %s\n",ctx->PMFdName);
    if (ctx->PMF1Def == 0) {
        printf1(ctx, "Error: need an output file.\n");
        goto SDDCWPFin;
    }
    if (alloc_acc(ctx, blen + 1))
        goto SDDCWPFin;
 
    nt1 = nt3 = id = npmax = nrec = 0;
    for (ii = 0; ii < 2; ++ii) {
        nobje = nobj = npt = rec = 0;
        while (fgets(ctx->AcC,blen,ctx->PMFd)) {
            rec++;
            prn_message(ctx, rec,0,0);
            if (rec == 1) {
                if (ii == 0) {
                    printf1(ctx, "First record: %s\n",ctx->AcC);
                    strcpy(name,ctx->AcC);
                    p = name + strlen(name);
                    while (--p >= name) {
                        if (*p != ' ' && *p != LF && *p != CR)
                            break;
                    }
                    *++p = '\0';
                }
                continue;
            }
            p = skip_b(ctx, ctx->AcC);
            if (!strncmp(p,"END",3))
                break;

            if (sscanf(p,"%d",&m) != 1) {
                printf1(ctx, "Error: cannot read polygon ID in record %d.\n",rec);
                goto SDDCWPFin;
            }
            nobj++;

            np = 0;
            while (fgets(ctx->AcC,blen,ctx->PMFd)) {
                rec++;
                p = skip_b(ctx, ctx->AcC);
                if (!strncmp(p,"END",3))
                    break;

                if (sscanf(p,"%lg",&x) != 1) {
                    printf1(ctx, "Error: cannot read x coordinate in record %d.\n",rec);
                    goto SDDCWPFin;
                }
                p = skip_dbl(ctx, p);
                p = skip_b(ctx, p);
                if (sscanf(p,"%lg",&y) != 1) {
                    printf1(ctx, "Error: cannot read y coordinate in record %d.\n",rec);
                    goto SDDCWPFin;
                }
                if (ii && np < npmax) {
                    ctx->AcX[np] = x;
                    ctx->AcY[np] = y;
                }
                np++;
                npt++;
            }
            if (np == 0) {
                nobj++;
                continue;
            }
            if (ii == 0)
                npmax = imax(ctx, np,npmax);
            else {
                if (ctx->PMNC == 0) {
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,++id);
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,1);
                    rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,1);
                    fprintf(ctx->PMF1d,"%s\n",name);
                    nrec++;
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[0]);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[0]);
                    fprintf(ctx->PMF1d,"\n");
                    nrec++;
                    nt1++;
                }
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,++id);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,3);
                rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,np - 1);
                fprintf(ctx->PMF1d,"%s\n",name);
                nrec++;
                for (i = 1; i < np; ++i) {
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcX[i]);
                    rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,ctx->AcY[i]);
                    fprintf(ctx->PMF1d,"\n");
                    nrec++;
                }
                nt3++;
            }
        }
        prn_message(ctx, rec,1,0);
        if (ii == 0) {
            printf1(ctx, "Number of records: %d\n",rec);
            printf1(ctx, "Number of data blocks: %d\n",nobj);
            printf1(ctx, "Number of empty blocks: %d\n",nobje);
            printf1(ctx, "Number of data points: %d\n",npt);
            printf1(ctx, "Maximal number of points per object: %d\n",npmax);

            if (alloc_acx(ctx, npmax + 1))
                goto SDDCWPFin;
            if (alloc_acy(ctx, npmax + 1))
                goto SDDCWPFin;

            if (fseek(ctx->PMFd,(long)0,0)) {             
                printf1(ctx, "Error: cannot seek to begin of file.\n");
                goto SDDCWPFin;
            }
        }
    }
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);
    printf1(ctx, "Number of type 1 objects: %d\n",nt1);
    printf1(ctx, "Number of type 3 objects: %d\n",nt3);

    err = 0;     

SDDCWPFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdgshhs     Conversion of GSHHS data file.                              */
/*                                                                          */
/*  sdgshhs(                                                                */
/*      df=...,     output file (required)                                  */
/*      opt=...,    option, def. 1                                          */
/*                  1 = standard coding of longitudes                       */
/*                  2 = longitudes unchanged                                */
/*      level=...,  select polygons with specified level, def. 0            */
/*      fmt=...,    floating point print format, def. 10.5                  */
/*      dtda=...,   create TDA description file                             */
/*                                                                          */
/*  ) = gshhs_input_file;                                                   */
/*                                                                          */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int sdgshhs(TDAContext *ctx)  
{
    register int i;
   	int	max_east,err,nmax,rec,nrec,nst,np,nps,ii,nseg,nobj;
    double lon,lat,xmin = 0.0,xmax = 0.0,ymin = 0.0,ymax = 0.0;
   	struct	POINT p;
   	struct GSHHS h;
    int hb[11];                      /* raw header words, see below */
    short int hs[2];
    int gver;                        /* GSHHS data format version */

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Conversion of GSHHS files. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,10,1)) {     /* get parameters */
        goto SDGSHHSFin;
    }   
    if (ctx->PMOPT < 1 || ctx->PMOPT > 2)
        ctx->PMOPT = 1;

    /*  The upper bound was 4, which silently reset level=5 and level=6
        to 0 -- "all levels" -- so asking for the Antarctic ice front
        drew the whole file instead.  GSHHG added those two levels in
        2.3.0 (5 ice front, 6 grounding line) and this reader predates
        them.

        Levels are also selectable as a SET, encoded as a bitmask above
        64: bit 0 set means level 1, bit 1 level 2, and so on.  A plain
        1..6 still means exactly that one level, so nothing that worked
        before changes.  */
    if (ctx->PMLEVEL > 64) {
        ctx->PMLEVELMask = ctx->PMLEVEL - 64;
        ctx->PMLEVEL = 0;
    }
    else {
        ctx->PMLEVELMask = 0;
        if (ctx->PMLEVEL < 1 || ctx->PMLEVEL > 6)
            ctx->PMLEVEL = 0;
    }

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,5);

    if (ctx->PMF1Def == 0) { 
        printf1(ctx, "Error: need an output file.\n");
        goto SDGSHHSFin;
    }
    printf1(ctx, "Input file: %s\n",ctx->PMFdName);
    if (ctx->PMLEVEL)
        printf1(ctx, "Level selection: %d\n",ctx->PMLEVEL);
    else if (ctx->PMLEVELMask)
        printf1(ctx, "Level selection (set): %d\n",ctx->PMLEVELMask);
    newline(ctx);

    if (alloc_acxf(ctx, GSHHSNP + 2))
        goto SDGSHHSFin;
    if (alloc_acyf(ctx, GSHHSNP + 2))
        goto SDGSHHSFin;

    rec = nrec = nmax = np = nps = 0;
    /*  nobj numbers the OBJECTS written, nps counts the POLYGONS
        selected.  They differ once a polygon is split at the
        dateline: using one counter for both reported "2 polygons"
        and "4 polygons with selected level" for the same file.  */
    nobj = 0;
   	max_east = 270000000;
    nst = sizeof(struct POINT);

    /*  The header is read WORD BY WORD rather than as a struct, because
        GSHHS has two incompatible layouts and the file says which one it
        is only implicitly:

          v1.x  8 ints (id, n, level, west, east, south, north, area)
                followed by 2 shorts (greenwich, source)   -- 36 bytes
          v2.x  11 ints, where the third is a PACKED FLAG and three more
                (area_full, container, ancestor) follow    -- 44 bytes

        Both start with the same 8 ints, so those are read first and the
        third one decides: in v1 it is `level`, which is 1..4; in v2 it
        is the flag, whose version byte puts it at 256 or above.  That
        is the discriminator -- no guessing from the file size, which
        fails on a subset like gshhs_l.b.

        Reading v2 files with the v1 struct is what made TDA stop with
        "while reading the file": the first polygon's points were read
        from 8 bytes into its own header.  */
   	while (fread((void *)hb,sizeof(int),(size_t)8,ctx->PMFd) == 8) {

        prn_message(ctx, ++rec,0,0);
        np++;
        if (ARCHTyp == 2) {
            for (i = 0; i < 8; ++i)
                hb[i] = swabi4 ((unsigned int)hb[i]);
        }
        h.id = hb[0];
        h.n = hb[1];
        h.west = hb[3];
        h.east = hb[4];
        h.south = hb[5];
        h.north = hb[6];
        h.area = hb[7];

        if (hb[2] < 0 || hb[2] > 4) {         /* v2.x: packed flag */
            if (fread((void *)(hb + 8),sizeof(int),(size_t)3,ctx->PMFd) != 3) {
                printf1(ctx, "Error: while reading the file.\n");
                goto SDGSHHSFin;
            }
            if (ARCHTyp == 2) {
                for (i = 8; i < 11; ++i)
                    hb[i] = swabi4 ((unsigned int)hb[i]);
            }
            gver = (hb[2] >> 8) & 255;
            h.level = hb[2] & 255;
            h.greenwich = (short int)((hb[2] >> 16) & 3);
            h.source = (short int)((hb[2] >> 24) & 1);
        }
        else {                                 /* v1.x: level + 2 shorts */
            if (fread((void *)hs,sizeof(short int),(size_t)2,ctx->PMFd) != 2) {
                printf1(ctx, "Error: while reading the file.\n");
                goto SDGSHHSFin;
            }
            if (ARCHTyp == 2) {
                hs[0] = (short int)swabi2 ((unsigned int)hs[0]);
                hs[1] = (short int)swabi2 ((unsigned int)hs[1]);
            }
            gver = 1;
            h.level = hb[2];
            h.greenwich = hs[0];
            h.source = hs[1];
        }
        if (np == 1)
            printf1(ctx, "GSHHS data format version: %d\n",gver);
        nmax = imax(ctx, nmax,h.n);

        /*  Grow the point buffers if this polygon needs more than the
            fixed GSHHSNP the initial allocation used.  That constant is
            whatever the largest polygon happened to be in the release
            Rohwer built against; GSHHG has grown since, and the loop
            below writes AcXF[i] for i < h.n with no check at all.  A
            polygon larger than the buffer would therefore write past the
            end of the allocation -- the same failure that made rxls
            abort and take the R session with it, which is why this is
            checked rather than assumed to be big enough.  */
        if (h.n + 2 > ctx->AcXFN) {
            if (alloc_acxf(ctx, h.n + 2))
                goto SDGSHHSFin;
            if (alloc_acyf(ctx, h.n + 2))
                goto SDGSHHSFin;
            printf1(ctx, "Point buffer grown to %d points.\n",h.n + 2);
        }
                 
        if ((ctx->PMLEVEL && h.level != ctx->PMLEVEL) ||
            (ctx->PMLEVELMask && (h.level < 1 || h.level > 30 ||
             !(ctx->PMLEVELMask & (1 << (h.level - 1)))))) {
            if (fseek(ctx->PMFd,(long)(h.n * nst),SEEK_CUR)) {
                printf1(ctx, "gshhs: seek error.\n");
                goto SDGSHHSFin;
            }
            continue;
        }
        nps++;

      		/*****************************************************
        source = (h.source == 1) ? 'W' : 'C'; 
      		w = h.west  * 1.0e-6;
      		e = h.east  * 1.0e-6;
      		s = h.south * 1.0e-6;
      		n = h.north * 1.0e-6;
        fprintf(PMF1d,"P%d %6d %8d %2d %1d %13.3lf %10.5lf %10.5lf %10.5lf %10.5lf\n",
                h.greenwich,h.id, h.n, h.level, h.source, area, w, e, s, n);
        ********************************************************************/

        for (i = 0; i < h.n; ++i) {
         			if (fread ((void *)&p,(size_t)sizeof(struct POINT),(size_t)1,ctx->PMFd) != 1) {
                printf1(ctx, "Error: while reading the file.\n");
                goto SDGSHHSFin; 
            }
            prn_message(ctx, ++rec,0,0);

            if (ARCHTyp == 2) {
                p.x = swabi4((unsigned int)p.x);
                p.y = swabi4((unsigned int)p.y);
            }
            /*  TWO wraps, and opt only ever gated the second.  The
                first fires whenever the polygon has the greenwich flag
                and the point is east of max_east -- which drops from
                270 to 180 degrees after the first polygon, so from the
                second polygon on, ANY point past 180 was moved to the
                far side however opt was set.  A polygon spanning the
                dateline then has consecutive points 350 degrees apart
                and is drawn straight across the map.

                opt=2 is "longitudes unchanged" in TDA's own help, so
                leaving one of the two wraps running contradicts the
                option's documented meaning.  Both are gated now.  */
            if (ctx->PMOPT == 1) {
         			    lon = (h.greenwich && p.x > max_east)
                        ? p.x * 1.0e-6 - 360.0 : p.x * 1.0e-6;
                if (lon > 180.0)
                    lon -= 360.0;
            }
            else {
                /*  opt=2 is "longitudes unchanged" in TDA's help, and
                    taken literally that cannot draw a world map: a real
                    GSHHS file does not use ONE convention.  In
                    gshhs_l.b 5201 polygons run past 180 on a 0..360
                    basis while the two Antarctic ones are negative, so
                    left raw they span -180..360 -- 540 degrees, with the
                    world squeezed into part of it.  Normalising the
                    negatives up gives one consistent 0..360 map, which
                    is what the option is for.  */
                lon = p.x * 1.0e-6;
                if (lon < 0.0)
                    lon += 360.0;
            }
         			lat = p.y * 1.0e-6;

            ctx->AcXF[i] = (float)lon;   
            ctx->AcYF[i] = (float)lat;

            if (i == 0) {
                xmin = xmax = lon;
                ymin = ymax = lat;
            }
            else {
                xmin = dmin(ctx, xmin,lon);
                xmax = dmax(ctx, xmax,lon);
                ymin = dmin(ctx, ymin,lat);
                ymax = dmax(ctx, ymax,lat);
            }
        }
      		max_east = 180000000;	      /* Only Eurasiafrica needs 270 */
        
        /*  SPLIT AT THE DATELINE.  With opt=1 a longitude past 180 is
            moved to the far side point by point, so a polygon that
            spans the dateline ends up with two consecutive points
            ~360 degrees apart.  Drawn as one polyline that is a
            straight streak across the whole map -- an artefact of the
            wrapping, carrying no data.

            The polygon is therefore emitted as one object per run of
            points that does NOT jump, so each piece is drawn where it
            belongs and nothing is drawn between them.  A polygon that
            does not cross produces exactly one object, byte for byte
            what it produced before.

            Only opt=1 needs this: opt=2 leaves longitudes on 0..360,
            where a spanning polygon is already continuous.  */
        /*  A polygon that is split is no longer a closed area: each
            piece is an open run of coastline, so the pieces go out as
            type 2 (line) rather than type 3 (polygon).  An unsplit
            polygon stays type 3, exactly as before.  TDA also refuses a
            type 3 object with fewer than three points, which a split
            can easily produce.  */
        nseg = 1;
        {
            /*  A jump of more than 180 degrees between consecutive
                points is a seam crossing in EITHER convention: the
                dateline at +-180 under opt=1, the 0/360 meridian under
                opt=2.  Both need the same treatment, so this is no
                longer gated on opt.  */
            for (i = 1; i < h.n; ++i) {
                double d = (double)ctx->AcXF[i] - (double)ctx->AcXF[i - 1];
                if (d > 180.0 || d < -180.0)
                    nseg++;
            }
            /*  A type 3 object is CLOSED by the renderer: it draws a
                final segment from the last point back to the first.
                Antarctica runs from +180 to -180, so that closing
                segment spans the whole map and is drawn as a streak
                even though no two consecutive points jump.  Such a
                polygon goes out as a line instead, leaving the two ends
                unjoined -- which is what they are, once the coastline
                has been cut at the dateline.  */
            {
                double dc = (double)ctx->AcXF[0] - (double)ctx->AcXF[h.n - 1];
                if (h.n > 1 && (dc > 180.0 || dc < -180.0))
                    nseg++;
            }
        }
        for (ii = 0; ii < h.n; ) {
            int seg_end = h.n;
            int styp = (nseg > 1) ? 2 : 3;
            double sxmin, sxmax, symin, symax;

            for (i = ii + 1; i < h.n; ++i) {
                double d = (double)ctx->AcXF[i] - (double)ctx->AcXF[i - 1];
                if (d > 180.0 || d < -180.0) {
                    seg_end = i;
                    break;
                }
            }
            sxmin = sxmax = (double)ctx->AcXF[ii];
            symin = symax = (double)ctx->AcYF[ii];
            for (i = ii; i < seg_end; ++i) {
                double lx = (double)ctx->AcXF[i];
                double ly = (double)ctx->AcYF[i];
                sxmin = dmin(ctx, sxmin,lx);
                sxmax = dmax(ctx, sxmax,lx);
                symin = dmin(ctx, symin,ly);
                symax = dmax(ctx, symax,ly);
            }
            /*  A run of one point is still a vertex of the original
                coastline, so it goes out as a type 1 (point) object
                rather than being dropped -- discarding it would lose
                data to make the drawing tidy.  */
            if (seg_end - ii < 2)
                styp = 1;
            nobj++;
            fprintf(ctx->PMF1d,"%10d %d %7d %d ",nobj,styp,
                    seg_end - ii,h.level);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,sxmin);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,sxmax);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,symin);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,symax);
            fprintf(ctx->PMF1d,"%10d %d %d\n",h.area,h.greenwich,h.source);
            nrec++;

            for (i = ii; i < seg_end; ++i) {
                lon = (double)ctx->AcXF[i];
                lat = (double)ctx->AcYF[i];
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,lon);
                rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,lat);
                fprintf(ctx->PMF1d,"\n");
                nrec++;
            }
            ii = seg_end;
        }
    }
    prn_message(ctx, rec,1,0);

    printf1(ctx, "Number of polygons: %d (max number of points: %d)\n",np,nmax);
    printf1(ctx, "Number of polygons with selected level: %d\n",nps);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMF1dName);

    if (ctx->PMTDAFDef) {                /* create TDA decription file */

        fprintf(ctx->PMTDAFd,"sdnvar(\n");
        fprintf(ctx->PMTDAFd,"  dfile = %s,\n",ctx->PMF1dName);
        fprintf(ctx->PMTDAFd,"  noc = %d,\n",nps);
        fprintf(ctx->PMTDAFd,"  SDID     [10.0] = c1,\n");
        fprintf(ctx->PMTDAFd,"  SDTyp <1>[ 1.0] = c2,\n");
        fprintf(ctx->PMTDAFd,"  SDN   <5>[ 7.0] = c3,\n");
        fprintf(ctx->PMTDAFd,"  SDPtr <5>[10.0] = rd,\n");
        fprintf(ctx->PMTDAFd,"  Level <1>[ 1.0] = c4,\n");
        fprintf(ctx->PMTDAFd,"  XMin  <4>[10.5] = c5,\n");
        fprintf(ctx->PMTDAFd,"  XMax  <4>[10.5] = c6,\n");
        fprintf(ctx->PMTDAFd,"  YMin  <4>[10.5] = c7,\n");
        fprintf(ctx->PMTDAFd,"  YMax  <4>[10.5] = c8,\n");
        fprintf(ctx->PMTDAFd,"  Area  <4>[10.0] = c9,\n");
        fprintf(ctx->PMTDAFd,"  GFlag <1>[ 1.0] = c10,\n");
        fprintf(ctx->PMTDAFd,"  Source<1>[ 1.0] = c11,\n");
        fprintf(ctx->PMTDAFd,");\n");
        printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
    }
    err = 0;

SDGSHHSFin:
    p_clean(ctx);
    return(err);
}


