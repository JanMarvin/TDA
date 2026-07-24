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

/*  functions in t_sdx.c */

int sdgen(void); 
void sdgenv(int irec); 
int sdshp(void); 
int shp_typ(int n,int opt);
void shp_dbfw(char *buf,int rlen);
int shp_dbfh(int opt,FILE *fd,char *buf,int *nrec,int *hlen,int *nvar,int *rlen,int ns);
int shp_dbff(FILE *fd,char *buf,int nvar);
int shp_dbf_a(int nvar,int opt);

int rdbf(void); 
int prn_dbf_var(int nvar);
int rcsv(void); 
int rplz(void); 

int gtopo(void); 
int gtopo_rpixel(int i,int j,int nc,int *err);
int sddcwp(void); 
int sdgshhs(void); 

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */
 
#define HLEN 100            /* index/shape file header length               */
#define RHLEN  8            /* index file record length                     */
#define DBFRLEN 32          /* dbf file record length                       */
#define DBFFLEN 11          /* dbf file field length                        */

char *DBFFN;                /* dbf: field name                              */
int *DBFFT;                 /* dbf: field type                              */
int *DBFFL;                 /* dbf: field length                            */
int *DBFFD;                 /* dbf: field decimal count                     */
int *DBFFP;                 /* dbf: field pointer                           */
int DBFFN_A = 0;            /* flag for memory allocation                   */
int DBFFT_A = 0;
int DBFFL_A = 0;
int DBFFD_A = 0;
int DBFFP_A = 0;

/* ------------------------------------------------------------------------ */
/* Parameters for gshhs.                                                    */

#define GSHHSNP 1435084     /* max number of points in a polygon */     

#define swabi2(i2) (((i2) >> 8) + (((i2) & 255) << 8))
#define swabi4(i4) (((i4) >> 24) + (((i4) >> 8) & 65280) + (((i4) & 65280) << 8) + (((i4) & 255) << 24))

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

int sdgen(void)  
{
    register int i,j;
    int err,nrec,n,idx,idy,id,ida,ia;
    double x,y,a,b;

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Generating a spatial data file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,1,1))       /* get parameters */
        goto SDGENFin;
        
    if (PMNV != 2) {
        printf1("Error: need two variables for coordinates.\n");
        goto SDGENFin;
    }
    if (PMID < 0) {
        printf1("Error: need an ID variable.\n");
        goto SDGENFin;
    }
    if (PMTyp < 1)
        PMTyp = 1;

    if (PMTyp > 3) {
        printf1("Error in type specification.\n");
        goto SDGENFin;
    }
    if (PMFmtF == 0)
        pmfmt(10,6);
     
    if (alloc_acx(NOC + 1))
        goto SDGENFin;
    if (alloc_acy(NOC + 1))
        goto SDGENFin;

    idx = PMVIdx[0];
    idy = PMVIdx[1];

    i = nrec = 0;
    while (i < NOC) {
        ia = i;
        ida = (int)get_data(PMID,ia);
        n = 0;
        while (i < NOC) {
            AcX[n] = get_data(idx,i);
            AcY[n] = get_data(idy,i);
            n++;
            i++;
            if (PMTyp == 1 || (id = (int)get_data(PMID,i)) != ida)
                break;
        }
        fprintf(PMFd,PMNFmtS,ida);
        if (n == 1)        
            j = 1;
        else if (n == 2)
            j = 2;
        else
            j = PMTyp;
        fprintf(PMFd,PMNFmtS,j);
        fprintf(PMFd,PMNFmtS,n);
        sdgenv(ia);
        fprintf(PMFd,"\n");
        nrec++; 

        for (j = 0; j < n; ++j) {
            fprintf(PMFd,PMFmtS,AcX[j]);
            fprintf(PMFd,PMFmtS,AcY[j]);
            fprintf(PMFd,"\n");
            nrec++; 
        }
        ida = id;
    }
      
    printf1("%d records written to: %s\n",nrec,PMFdName);

    err = 0;     

SDGENFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  sdgenv()    Print additional attribute variables. Called by sdgen().    */

void sdgenv(int irec)
{
    register int i,j;
    double tmp;

    for (i = 0; i < PM1NV; ++i) {
        j = PM1VIdx[i];
        if (VTyp[j] != 1) {
            tmp = get_data(j,irec);                  
            fprintf(PMFd,VPFmtS[j],tmp);
        }
        else {          /* string variables */
            get_str(SVBuf,j,irec);
            fprintf(PMFd,"%s ",SVBuf);
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

int sdshp(void)  
{
    FILE *dbffd,*shpfd,*shxfd;
    register int i,j,k,l;
    int err,n,m,m1,mn,ns,np,typ,len,flen,nrec,hlen,nvar,rlen,n0,n1,drec,drec1;
    int dbfbufa,shpbufa,shpbufl,dbffd_open,shpfd_open,shxfd_open,npart,npnt;
    int id,id1,id2;
    double x,y,xmin,xmax,ymin,ymax;
    char *dbfbuf,*shpbuf,*p,*q;

    err = -1;
    shpbufa = dbfbufa = dbffd_open = shpfd_open = shxfd_open = 0;

    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading a shape file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,8,1)) {     /* get parameters */
        goto SDSHPFin;
    }   
    if (PMFmtF == 0)
        pmfmt(12,4);

    printf1("File: %s\n",PMRHSTR);

    shpbufl = 10000;
    if (!(shpbuf = (char *)calloc(shpbufl + 2,sizeof(char)))) {
        p_err(-2,1);
        goto SDSHPFin;
    }
    shpbufa = shpbufl + 2;
    memrq(shpbufa,sizeof(char));

    strcpy(shpbuf,PMRHSTR);                         /* open index file */
    strcat(shpbuf,".shx");
    if (!(shxfd = fopen(shpbuf,OPEN_RD))) {   
        printf1("Error: can't open the index file.\n");         
        goto SDSHPFin;
    }
    shxfd_open = 1;

    strcpy(shpbuf,PMRHSTR);                         /* open shape file */
    strcat(shpbuf,".shp");
    if (!(shpfd = fopen(shpbuf,OPEN_RD))) {   
        printf1("Error: can't open the shp file.\n");         
        goto SDSHPFin;
    }
    shpfd_open = 1; 

    strcpy(shpbuf,PMRHSTR);                         /* open dbf file */
    strcat(shpbuf,".dbf");
    if (!(dbffd = fopen(shpbuf,OPEN_RD))) {   
        printf1("Error: can't open the dbf file.\n");         
        goto SDSHPFin;
    }
    dbffd_open = 1;
                        
    printf1("\nReading the index file.\n");

    if (fread(shpbuf,sizeof(char),HLEN,shxfd) != HLEN) {
        p_err(-7,1);
        goto SDSHPFin;
    }
    HILO = 1;
    n = st_geti(shpbuf);
    if (n != 9994) {
        printf1("Probably not a proper index/shape file.\n");
        goto SDSHPFin;
    }
    flen = 2 * st_geti(shpbuf + 24);
    flen -= HLEN;
    ns = flen / 8;                  /* number of shapes */ 
 
    HILO = 2;
    n = st_geti(shpbuf + 28);
    printf1("Version: %d\n",n);
    typ = st_geti(shpbuf + 32);
    if (shp_typ(typ,1))             /* check shape type */
        goto SDSHPFin;

    printf("Number of shapes: %d\n",ns);

    xmin = st_getd(shpbuf + 36);
    ymin = st_getd(shpbuf + 44);
    xmax = st_getd(shpbuf + 52);
    ymax = st_getd(shpbuf + 60);

    printf1("\nBounding box X: %20.8lf  %20.8lf\n",xmin,xmax);
    printf1("Bounding box Y: %20.8lf  %20.8lf\n",ymin,ymax);

    printf1("\nReading the dbf file.\n");

    if (shp_dbfh(1,dbffd,shpbuf,&nrec,&hlen,&nvar,&rlen,ns))  /* read dbf file header */
        goto SDSHPFin;

    n = imax(rlen,DBFRLEN);
    if (!(dbfbuf = (char *)calloc(n + 2,sizeof(char)))) {
        p_err(-2,1);
        goto SDSHPFin;
    }
    dbfbufa = n + 2;
    memrq(dbfbufa,sizeof(char));
                  
    if (shp_dbf_a(nvar,1))              /* allocate memory for fields */
        goto SDSHPFin;
         
    if (shp_dbff(dbffd,dbfbuf,nvar))    /* read fields */
        goto SDSHPFin;

    if (prn_dbf_var(nvar))              /* info about dbf variables */
        goto SDSHPFin;

    if (PMF1Def == 0) {                 /* without output file end here */
        err = 0;
        goto SDSHPFin;
    }
    printf1("Creating the output data file: %s\n",PMF1dName);
    drec1 = drec = n1 = n0 = 0;

    /* skip rest of dbf header */

    if ((hlen -= (nvar + 1) * DBFRLEN) < 1) {
        printf1("Error in dbf header length.\n");
        goto SDSHPFin;
    }
    if (fseek(dbffd,(long)hlen,SEEK_CUR)) {
        printf1("Error: cannot seek to begin of fields in dbf file.\n");
        goto SDSHPFin;
    }
    id = id1 = 0;
    for (i = 0; i < ns; ++i) {

        /* read current record from the dbf file */

        if (fread(dbfbuf,sizeof(char),rlen,dbffd) != rlen) {
            p_err(-7,1);
            goto SDSHPFin;
        }

        /* read the index file */

        if (fread(shpbuf,sizeof(char),8,shxfd) != 8) {
            p_err(-7,1);
            goto SDSHPFin;
        }
        HILO = 1;
        n = st_geti(shpbuf);        /* n = pointer to shape file */
        len = st_geti(shpbuf + 4);  /* len = size of entries */ 

        /* read the shape file */

        if (fseek(shpfd,(long)(2 * n),0)) {
            printf1("Can't seek to offset %d in shape file.\n",2 * n);
            goto SDSHPFin;
        }
        if (fread(shpbuf,sizeof(char),8,shpfd) != 8) {
            p_err(-7,1);
            goto SDSHPFin;
        }
        m = st_geti(shpbuf + 4);
        if (m != len) {
            printf1("Error: inconsistent information in index and shape file.\n");
            goto SDSHPFin;
        }
        len *= 2;
        if (len > shpbufl) {            /* adjust buffer length */
            free(shpbuf);
            memrq(-shpbufa,sizeof(char));

            shpbufl = len;
            if (!(shpbuf = (char *)calloc(shpbufl + 2,sizeof(char)))) {
                p_err(-2,1);
                goto SDSHPFin;
            }
            shpbufa = shpbufl + 2;
            memrq(shpbufa,sizeof(char));
        }
        if (fread(shpbuf,sizeof(char),len,shpfd) != len) {
            p_err(-7,1);
            goto SDSHPFin;
        }
        HILO = 2;
        typ = st_geti(shpbuf);            /* shape type */
        id1++;
        switch (typ) {
            case 0:                                 /* null shape */
                n0++;
                break;   

            case 1:                                 /* point */
                id++;
                id2 = 1;
                np = mn = 1;
                fprintf(PMF1d,"%10d %1d %7d %10d %4d ",id,mn,np,id1,id2);
                shp_dbfw(dbfbuf,rlen);
                fprintf(PMF1d,"\n");
                drec++;

                x = st_getd(shpbuf + 4);
                y = st_getd(shpbuf + 12);
                fprintf(PMF1d,PMFmtS,x);
                fprintf(PMF1d,PMFmtS,y);
                fprintf(PMF1d,"\n");
                drec1++;
                break;

            case 3:                                 /* polyline */
                mn = 2;

            case 5:                                 /* polygon */
                if (typ == 5)
                    mn = 3;

                npart = st_geti(shpbuf + 36);
                npnt  = st_geti(shpbuf + 40);

                if (npart > 1)
                    n1++;

                p = shpbuf + 44;
                m = st_geti(p);
                p += 4;

                id2 = 0;
                for (j = 1; j <= npart; ++j) {
                    id++;
                    id2++;

                    if (j == npart)
                        m1 = npnt;
                    else {
                        m1 = st_geti(p);
                        p += 4;
                    }
                    np = m1 - m; 

                    fprintf(PMF1d,"%10d %1d %7d %10d %4d ",id,mn,np,id1,id2);
                    shp_dbfw(dbfbuf,rlen);
                    fprintf(PMF1d,"\n");
                    drec++;

                    q = shpbuf + 44 + 4 * npart + m * 16;

                    for (k = 0; k < np; ++k) {
                        x = st_getd(q);
                        q += 8;
                        y = st_getd(q);
                        q += 8;
                        fprintf(PMF1d,PMFmtS,x);
                        fprintf(PMF1d,PMFmtS,y);
                        fprintf(PMF1d,"\n");
                        drec1++;
                    }
                    m = m1;
                }
                break;

            default:
                printf1("Error: found unknown shape type %d.\n",typ);
                goto SDSHPFin;
        }
        prn_message(i,0,0);
    }
    prn_message(i,1,0);
         
    printf1("Object-description records: %d. Data records: %d\n",drec,drec1);
    printf1("Total number of records written to output file: %d\n",drec + drec1);
    if (n0 > 0)
        printf1("Note: ignored %d records containing null type shapes.\n",n0);
    if (n1 > 0)
        printf1("Note: %d shapes contain more than one part.\n",n1);

    if (PMTDAFDef) {                /* create TDA decription file */

        fprintf(PMTDAFd,"sdnvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
        fprintf(PMTDAFd,"  noc = %d,\n",drec);
        fprintf(PMTDAFd,"  ffmt=c1(1-10),c2(12),c3(14-20),c4(22-31),c5(33-36),");
            
        n = 5; j = 6; l = 38;
        for (i = 0; i < nvar; ++i) {
            k = l + DBFFL[i] - 1;
            if (DBFFT[i] == 0) {
                if (++n > 5) {
                    fprintf(PMTDAFd,"\n       ");
                    n = 1;
                }
                fprintf(PMTDAFd,"c%d(%d-%d),",j,l,k);
                j++;
            }
            l = k + 2;
        }
        fprintf(PMTDAFd,"\n\n  SDID<5>[10.0] = c1,\n");
        fprintf(PMTDAFd,"  SDTyp<1>[1.0] = c2,\n");
        fprintf(PMTDAFd,"  SDN<5>[7.0] = c3,\n");
        fprintf(PMTDAFd,"  SDPtr<5>[10.0] = rd,\n");
        fprintf(PMTDAFd,"  SDID1<5>[10.0] = c4,\n");
        fprintf(PMTDAFd,"  SDID2<5>[ 4.0] = c5,\n");

        j = 6; l = 38;
        for (i = 0; i < nvar; ++i) {
            k = l + DBFFL[i] - 1;
            if (DBFFT[i] == 0) {
                if (DBFFD[i] > 0)
                    n = 8;
                else
                    n = 5;

                fprintf(PMTDAFd,"  %s<%d>[%d.%d] = c%d,\n",DBFFN + i * DBFFLEN,
                                            n,DBFFL[i],DBFFD[i],j);
                j++;
            }
            else if (DBFFT[i] == 1) {
                fprintf(PMTDAFd,"  %s = str(%d,%d),\n",DBFFN + i * DBFFLEN,l,k);
            }
            l = k + 2;
        }
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;
    newline();

SDSHPFin:
    if (dbfbufa > 0) {
        free(dbfbuf);
        memrq(-dbfbufa,sizeof(char));
    }
    if (shpbufa > 0) {
        free(shpbuf);
        memrq(-shpbufa,sizeof(char));
    }

    shp_dbf_a(nvar,0);

    if (dbffd_open)  
        fclose(dbffd);
    if (shpfd_open)  
        fclose(shpfd);
    if (shxfd_open)  
        fclose(shxfd);

    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  shp_typ(n,opt)  Check shape type. Return 0 if OK, -1 if not defined.    */
/*                  If opt != 0 print.                                      */

int shp_typ(int n,int opt)
{
    if (opt)
        printf1("Shape type: ");

    switch (n) {
        case 1:     if (opt) printf1("Point");
                    break;
        case 3:     if (opt) printf1("PolyLine");
                    break;
        case 5:     if (opt) printf1("Polygon");
                    break;
        /*****************************************
        case 8:     if (opt) printf1("MultiPoint");
                    break;
        case 11:    if (opt) printf1("PointZ");
                    break;
        case 13:    if (opt) printf1("PolyLineZ");
                    break;
        case 15:    if (opt) printf1("PolygonZ");
                    break;
        case 18:    if (opt) printf1("MultiPointZ");
                    break;
        case 21:    if (opt) printf1("PointM");
                    break;
        case 23:    if (opt) printf1("PolyLineM");
                    break;
        case 25:    if (opt) printf1("PolygonM");
                    break;
        case 28:    if (opt) printf1("MultiPointM");
                    break;
        case 31:    if (opt) printf1("MultiPatch");
                    break;
        ******************************************/
        default:    if (opt) printf1("%d (unknown)\n",n);
                    return(-1);
    }
    if (opt)
        printf1(" (%d)\n",n);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbfw(buf,rlen)  Write current dbf record to output file.            */

void shp_dbfw(char *buf,int rlen)
{
    register int j,k,l;
    char *p;

    p = buf + 1;
    k = l = 0;
    for (j = 1; j < rlen; ++j) {
        if (*p == '\0')
            *p = ' ';
        fprintf(PMF1d,"%c",*p++);
        if (++l >= DBFFL[k]) {
            fprintf(PMF1d," ");
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

int shp_dbfh(int opt,FILE *fd,char *buf,int *nrec,int *hlen,int *nvar,int *rlen,int ns)
{
    if (fread(buf,sizeof(char),DBFRLEN,fd) != DBFRLEN) {
        printf1("Error: can't read header of dbf file.\n");
        return(-1);
    }
    HILO = 2;
    *nrec = st_geti(buf + 4);
    *hlen = st_gets(buf + 8);
    *rlen = st_gets(buf + 10);
   
    if (opt && *nrec != ns) {
        printf1("Error: number of records in dbf file is %d.\n",*nrec);
        printf1("Inconsistent with information in index/shape file.\n");
        return(-1);
    }
    printf("Number of records: %d\n",*nrec);
    printf("Record length: %d\n",*rlen);
    if (*rlen < 1) {
        printf1("This is probably not a dbf file.\n");
        return(-1);
    }
    *nvar = *hlen - 32;
    *nvar /= 32;

    printf("Number of fields: %d\n",*nvar);
    if (*nvar < 1) {
        printf1("This is probably not a dbf file.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbff(nvar)        Read dbf fields and save information in gloabal   */
/*                        arrays DBFFN, DBFFT, DBFFL, DBFFD, DBFFP.         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int shp_dbff(FILE *fd,char *buf,int nvar)
{
    register int i;
    char *vtyp; 
    unsigned char u;

    for (i = 0; i < nvar; ++i) {
        if (fread(buf,sizeof(char),DBFRLEN,fd) != DBFRLEN) {
            printf1("Error: can't read fields of dbf file.\n");
            return(-1);
        }
        strncpy(DBFFN + i * DBFFLEN,buf,DBFFLEN);
        vtyp = buf + DBFFLEN;
        if (*vtyp == 'N' || *vtyp == 'F') {
            DBFFT[i] = 0;
            DBFFD[i] = (int)buf[17];
        }
        else if (*vtyp == 'C')
            DBFFT[i] = 1;
        else
            DBFFT[i] = 2;

        u = (unsigned char)buf[16];
        DBFFL[i] = (int)u;        

        if (i == 0)
            DBFFP[i] = 1;
        else
            DBFFP[i] = DBFFL[i - 1] + DBFFP[i - 1];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  shp_dbf_a(nvar,opt)   If opt != 0 allocate arrays for n dbf fields,     */
/*                        otherwise free previously allocated memory.       */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int shp_dbf_a(int nvar,int opt)
{

    if (DBFFN_A > 0) {
        free((char *)DBFFN);
        memrq(-DBFFN_A,sizeof(char));
        DBFFN_A = 0;
    }
    if (DBFFT_A > 0) {
        free((char *)DBFFT);
        memrq(-DBFFT_A,sizeof(int));
        DBFFT_A = 0;
    }
    if (DBFFL_A > 0) {
        free((char *)DBFFL);
        memrq(-DBFFL_A,sizeof(int));
        DBFFL_A = 0;
    }
    if (DBFFD_A > 0) {
        free((char *)DBFFD);
        memrq(-DBFFD_A,sizeof(int));
        DBFFD_A = 0;
    }
    if (DBFFP_A > 0) {
        free((char *)DBFFP);
        memrq(-DBFFP_A,sizeof(int));
        DBFFP_A = 0;
    }
    if (opt == 0)
        return(0);

    if (!(DBFFN = (char *)calloc(nvar * (DBFFLEN + 1),sizeof(char)))) {
        p_err(-2,1);
        return(-1);    
    }
    DBFFN_A = nvar * (DBFFLEN + 1);
    memrq(DBFFN_A,sizeof(char));

    if (!(DBFFT = (int *)calloc(nvar,sizeof(int)))) {
        p_err(-2,1);
        return(-1);    
    }
    DBFFT_A = nvar;
    memrq(DBFFT_A,sizeof(int));

    if (!(DBFFL = (int *)calloc(nvar,sizeof(int)))) {
        p_err(-2,1);
        return(-1);    
    }
    DBFFL_A = nvar;
    memrq(DBFFL_A,sizeof(int));

    if (!(DBFFD = (int *)calloc(nvar,sizeof(int)))) {
        p_err(-2,1);
        return(-1);    
    }
    DBFFD_A = nvar;
    memrq(DBFFD_A,sizeof(int));

    if (!(DBFFP = (int *)calloc(nvar,sizeof(int)))) {
        p_err(-2,1);
        return(-1);    
    }
    DBFFP_A = nvar;
    memrq(DBFFP_A,sizeof(int));
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

int rdbf(void)  
{
    int err,nrec,nvar,rlen,hlen,i,j,k,l,n,wrec;
    char *p,buf[DBFRLEN + 2];

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading a dBase file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,10,1)) {     /* get parameters */
        goto RDBFin;
    }   
    printf1("File name: %s\n\n",PMFdName);

    if (shp_dbfh(0,PMFd,buf,&nrec,&hlen,&nvar,&rlen,0))    /* read dbf file header */
        goto RDBFin;

    hlen -= (nvar + 1) * DBFRLEN;
    if (hlen < 1) {
        printf1("Error: in dbf header length.\n");
        goto RDBFin;
    }
    if (shp_dbf_a(nvar,1))          /* allocate memory for fields */
        goto RDBFin;

    if (shp_dbff(PMFd,buf,nvar))    /* read fields */
        goto RDBFin;

    if (prn_dbf_var(nvar))          /* info about dbf variables */
        goto RDBFin;

    if (nrec > 0 && PMF1Def) {      /* read data and print to output file */

        wrec = 0;
        if (alloc_acc(rlen + 1))
            goto RDBFin;

        /* skip rest of dbf header */

        if (fseek(PMFd,(long)hlen,SEEK_CUR)) {
            printf1("Error: cannot seek to begin of fields in dbf file.\n");
            goto RDBFin;
        }
        for (i = 0; i < nrec; ++i) {

            if (fread(AcC,sizeof(char),rlen,PMFd) != rlen) {
                p_err(-7,1);
                goto RDBFin;
            }
            p = AcC + 1;
            k = l = 0;
            for (j = 1; j < rlen; ++j) {
                if (*p == '\0')
                    *p = ' ';
                fprintf(PMF1d,"%c",*p++);
                if (++l >= DBFFL[k]) {
                    fprintf(PMF1d," ");
                    k++;
                    l = 0;
                }
            }
            fprintf(PMF1d,"\n");
            wrec++;
        }
        printf1("%d records written to: %s\n",wrec,PMF1dName);
    }
    err = 0;     

RDBFin:
    shp_dbf_a(nvar,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_dbf_var(nvar)       Print table with dbf variables.                 */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int prn_dbf_var(int nvar)
{
    int i,n;

    printf1("\nIdx  Type  Length  D.Count  Name\n");
    prnchar('-',32,1);

    n = 0;                                  /* save var definitions */
    for (i = 0; i < nvar; ++i) {
        if (DBFFL[i] < 1) {
            printf1("Error: in dbf field definitions.\n");
            return(-1);   
        }
        printf1("%3d %4d   %5d  %7d   %s\n",i + 1,DBFFT[i],DBFFL[i],
                            DBFFD[i],DBFFN + i * DBFFLEN);  

        if (DBFFT[i] != 0 && DBFFT[i] != 1)
            n++;
    }
    newline();
    if (n > 0) {
        printf1("Warning: the dbf file contains %d fields which are\n",n);
        printf1("neither numerical nor character. Will be ignored.\n\n");
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

int rcsv(void)  
{
    int err,nrec,c,ns,nsmax,j,n;

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading a csv file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,10,1)) {     /* get parameters */
        goto RCSVFin;
    }   
    printf1("File name: %s\n\n",PMFdName);

    nsmax = ns = nrec = 0; 
    while ((c = fgetc(PMFd)) != EOF) {
        if (c == ';')
            ns++;
        else if (c == '\n') {
            nrec++;
            nsmax = imax(nsmax,ns);
            ns = 0;
        }
    }
    printf("Number of records: %d\n",nrec);
    printf("Maximal number of entries per line: %d\n",nsmax);
    if (nsmax < 1)  
        goto RCSVFin;
                    
    if (alloc_acn(nsmax + 1))
        goto RCSVFin;

    if (fseek(PMFd,(long)0,0)) {             
        printf1("Error: cannot seek to begin of file.\n");
        goto RCSVFin;
    }
    n = j = 0; 
    while ((c = fgetc(PMFd)) != EOF) {
        if (c == ';') {
            if (j < nsmax)  
                AcN[j] = imax(AcN[j],n);
            j++;
            n = 0;
        }             
        else if (c == '\n')  
            j = 0;
        else
            n++;
    }
    if (fseek(PMFd,(long)0,0)) {             
        printf1("Error: cannot seek to begin of file.\n");
        goto RCSVFin;
    }
    n = j = 0; 
    while ((c = fgetc(PMFd)) != EOF) {
        if (c == ';') {
            if (j < nsmax) {
                while (n++ < AcN[j])
                    fputc(' ',stdout);
            }
            j++;
            n = 0;
        }             
        else if (c == '\n') {
            fputc(c,stdout);
            j = 0;
        }
        else {
            fputc(c,stdout);
            n++;
        }
    }

RCSVFin:
    p_clean();
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

int rplz(void)  
{
    int err,nrec,c,n,ns,rl,rlmax,id,nt,nl;
    double x,y;
    char *p,*q;

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading a plz file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,10,1)) {     /* get parameters */
        goto RPLZFin;
    }   
    printf1("File name: %s\n\n",PMFdName);

    rlmax = rl = nrec = 0; 
    while ((c = fgetc(PMFd)) != EOF) {
        rl++;
        if (c == '\n') {
            nrec++;
            rlmax = imax(rlmax,rl);
            rl = 0;
        }
    }
    printf1("Number of records: %d\n",nrec);
    printf1("Maximal record length: %d\n",rlmax);
    if (PMF1Def == 0) {
        printf1("Need an output file.\n");
        err = 0;
        goto RPLZFin;
    }
    if (fseek(PMFd,(long)0,0)) {             
        printf1("Error: cannot seek to begin of file.\n");
        goto RPLZFin;
    }
    if (alloc_acc(rlmax + 10))
        goto RPLZFin;

    id = nrec = 0;
    while (fgets(AcC,rlmax,PMFd)) {
        nrec++;
        if (nrec == 1)
            continue;

        ns = 0; 
        p = AcC;
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
                fprintf(PMF1d,"%12d %d %8d %s\n",id,nt,n,AcC);

                p = q;
                nl = 0;
                while (*p) {
                    if (*p == ',') {
                        fprintf(PMF1d,"  ");
                        nl = 1;
                        p++;
                    }
                    else if (*p == ' ') {
                        fprintf(PMF1d,"\n");
                        nl = 0;
                        p = skip_b(p);
                    }
                    else if (*p == ':') {
                        fprintf(PMF1d,"\n");
                        nl = 0;
                        p++;
                    }
                    else if (*p == '\n') {
                        fprintf(PMF1d,"\n");
                        nl = 0;
                        break;            
                    }
                    else {
                        fprintf(PMF1d,"%c",*p);
                        nl = 1;
                        p++;
                    }
                }
                if (nl)
                    fprintf(PMF1d,"\n");
                break;
            }
            p++; 
        }
        /************ 
        if (id > 10)
        break;
        ************/
    }
    printf1("%d records written to: %s\n",id,PMF1dName);
    err = 0;

RPLZFin:
    p_clean();
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

int gtopo(void)  
{
    register int i,j; 
    int err,err1,nrec,v,r,i0,i1,j0,j1;                 
    double dxa,dxb,dya,dyb,lona,lonb,lata,latb;

    err = -1;
    err1 = nrec = 0;

    if (check_cmd(0))
        return(-1);
                    
    printf1("Reading a GTOPO30 elevation file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,10,1)) {     /* get parameters */
        goto GTOPOFin;
    }   
    if (PMOPT != 1)
        PMOPT = 1;

    if (PMFmtF == 0)
        pmfmt(10,6);

    printf1("File name: %s\n\n",PMFdName);
    printf1("Number of rows: %d\n",PMRows);
    printf1("Number of columns: %d\n",PMCols);
    if (PMRows < 1 || PMCols < 1) {
        printf1("Error: need positive values.\n");
        goto GTOPOFin;
    }
    printf1("Coordinates: ulx = %18.14lf uly=%18.14lf\n",PMULX,PMULY);
    printf1("Grid spacing: dx = %18.14lf  dy=%18.14lf\n",PMDX,PMDY);
    if (PMDX < EPSI1 || PMDY < EPSI1) {
        printf1("Error: dx and dy must be strictly positive.\n");
        goto GTOPOFin;
    }
    dxa = PMULX;
    dyb = PMULY;
    dxb = dxa + (PMCols - 1) * PMDX;
    dya = dyb - (PMRows - 1) * PMDY;

    printf1("\nRange of the grid in decimal degrees.\n");
    printf1("Longitude: %18.14lf %18.14lf\n",dxa,dxb);
    printf1("Latitude:  %18.14lf %18.14lf\n\n",dya,dyb);

    if (PMF1Def == 0) {
        printf1("Without an output file: simply reading the file.\n");
        for (i = 0; i < PMRows; i++) {
            for (j = 0; j < PMCols; ++j) {
                v = gtopo_rpixel(i,j,PMCols,&r);
                if (r) {
                    printf1("Error: cannot read entry in row %d, column %d.\n",i + 1,j + 1);
                    goto GTOPOFin;
                }
            }
            prn_message(i + 1,0,0);
        }
        prn_message(PMRows,1,0);
        printf1("No errors occurred.\n");
        err = 0;
        goto GTOPOFin;
    }
    printf1("Output file: %s\n",PMF1dName);
    printf1("Option: %d\n\n",PMOPT);

    err1 = 1;
    if (PMNTP < 1 || PMNTP > 2 || PMNTP1 < 1 || PMNTP1 > 2)  
        goto GTOPOFin;
      
    if ((PMNTP == 2 && PMTP[0] > PMTP[1]) || (PMNTP1 == 2 && PMTP1[0] > PMTP1[1]))  
        goto GTOPOFin;
        
    for (i = 0; i < PMNTP; ++i) {
        if (PMTP[i] < -180.0 || PMTP[i] > 180.0)  
            goto GTOPOFin;
    }
    for (j = 0; j < PMNTP1; ++j) {
        if (PMTP1[j] < -90.0 || PMTP1[j] > 90.0)  
            goto GTOPOFin;
    }
    err1 = 0;

    if (PMOPT == 1) {
        if (PMNTP == 1 && PMNTP1 == 2) {
            lona = lonb = PMTP[0];
            lata = PMTP1[0];
            latb = PMTP1[1];
            if (lona < dxa || lonb > dxb || lata > dyb || latb < dya) {
                err1 = 2;
                goto GTOPOFin;
            }
            j0 = j1 = (int)((lona - dxa) / PMDX + 0.5);
            lona = lonb = dxa + j0 * PMDX;

            if (lata < dya)  
                i0 = PMRows - 1;
            else  
                i0 = (int)((dyb - lata) / PMDY + 0.5);
            lata = dyb - i0 * PMDY;

            if (latb > dyb)  
                i1 = 0;             
            else  
                i1 = (int)((dyb - latb) / PMDY + 0.5);
            latb = dyb - i1 * PMDY;
        }
        else if (PMNTP == 2 && PMNTP1 == 1) {
            lona = PMTP[0];
            lonb = PMTP[1];
            lata = latb = PMTP1[0];

            if (lata < dya || latb > dyb || lona > dxb || lonb < dxa) {
                err1 = 2;
                goto GTOPOFin;
            }
            i0 = i1 = (int)((lata - dya) / PMDY + 0.5);
            lata = latb = dya + i0 * PMDY;

            if (lona < dxa)  
                j0 = 0;
            else  
                j0 = (int)((lona - dxa) / PMDX + 0.5);
            lona = dxa + j0 * PMDX;

            if (lonb > dxb)  
                j1 = PMCols - 1;
            else  
                j1 = (int)((lonb - dxa) / PMDX + 0.5);
            lonb = dxa + j1 * PMDX;

        }
        else {   
            printf1("Syntax error in lon and/or lat parameters.\n");
            goto GTOPOFin;
        }
        if (i0 < i1 || j0 > j1) {
            err1 = 2;
            goto GTOPOFin;
        }
        printf1("Will use the following range of coordinates.\n");
        printf1("Longitude         Pixel  Latitude          Pixel\n");
        printf1("%16.12lf %6d  %16.12lf %6d\n",lona,j0,lata,i0);
        printf1("%16.12lf %6d  %16.12lf %6d\n\n",lonb,j1,latb,i1);
                
        for (i = i0; i >= i1; i--) {
            for (j = j0; j <= j1; ++j) {
                v = gtopo_rpixel(i,j,PMCols,&r);
                if (r) {
                    printf1("Error: cannot read entry in row %d, column %d.\n",i + 1,j + 1);
                    goto GTOPOFin;
                }
                fprintf(PMF1d,PMFmtS,dxa + j * PMDX);
                fprintf(PMF1d,PMFmtS,dyb - i * PMDY);
                fprintf(PMF1d,"%8d\n",v);
                nrec++;

            }
        }
        printf1("%d records written to: %s\n",nrec,PMF1dName);
    }
    err = 0;     

GTOPOFin:
    if (err1 == 1)  
        printf1("Invalid values for lon or lat parameter.\n");
    else if (err1 == 2)  
        printf1("Selected lon/lat range is empty or outside of the grid region.\n");

    p_clean();
    return(err);
}


/* -##--------------------------------------------------------------------- */
/*  gtopo_rpixel(i,j,nc,*err)                                               */
/*                                                                          */
/*  Return value of Pixel in row i, column j. nc is the number of columns   */
/*  in the file.  Return err = 0 if OK, -1 if an error occured.             */  
/*                                                                          */
/*  Note: row and column counting begins with 0.                            */

int gtopo_rpixel(int i,int j,int nc,int *err)
{
    int n;
    char buf[3];

    *err = -1;
    n = 2 * (i * nc + j);

    if (n < 0 || fseek(PMFd,(long)n,0))  
        return(0.0);
                                        
    if (fread(buf,sizeof(char),2,PMFd) != 2)  
        return(0.0);
                                               
    HILO = 1;
    n = st_gets(buf);
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

int sddcwp(void)  
{
    register int i; 
    int err,rec,nrec,blen,m,ii,np,npmax,npt,nobj,nobje,id,nt1,nt3;
    char *p,name[200];
    double x,y;

    blen = 1000;                /* read buffer length */
    err = -1;

    if (check_cmd(0))
        return(-1);
                    
    printf1("Conversion of DCW polygon point data file. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,10,1)) {     /* get parameters */
        goto SDDCWPFin;
    }   
    if (PMFmtF == 0)
        pmfmt(12,6);

    printf1("Input file: %s\n",PMFdName);
    if (PMF1Def == 0) {
        printf1("Error: need an output file.\n");
        goto SDDCWPFin;
    }
    if (alloc_acc(blen + 1))
        goto SDDCWPFin;
 
    nt1 = nt3 = id = npmax = nrec = 0;
    for (ii = 0; ii < 2; ++ii) {
        nobje = nobj = npt = rec = 0;
        while (fgets(AcC,blen,PMFd)) {
            rec++;
            prn_message(rec,0,0);
            if (rec == 1) {
                if (ii == 0) {
                    printf1("First record: %s\n",AcC);
                    strcpy(name,AcC);
                    p = name + strlen(name);
                    while (--p >= name) {
                        if (*p != ' ' && *p != LF && *p != CR)
                            break;
                    }
                    *++p = '\0';
                }
                continue;
            }
            p = skip_b(AcC);
            if (!strncmp(p,"END",3))
                break;

            if (sscanf(p,"%d",&m) != 1) {
                printf1("Error: cannot read polygon ID in record %d.\n",rec);
                goto SDDCWPFin;
            }
            nobj++;

            np = 0;
            while (fgets(AcC,blen,PMFd)) {
                rec++;
                p = skip_b(AcC);
                if (!strncmp(p,"END",3))
                    break;

                if (sscanf(p,"%lg",&x) != 1) {
                    printf1("Error: cannot read x coordinate in record %d.\n",rec);
                    goto SDDCWPFin;
                }
                p = skip_dbl(p);
                p = skip_b(p);
                if (sscanf(p,"%lg",&y) != 1) {
                    printf1("Error: cannot read y coordinate in record %d.\n",rec);
                    goto SDDCWPFin;
                }
                if (ii && np < npmax) {
                    AcX[np] = x;
                    AcY[np] = y;
                }
                np++;
                npt++;
            }
            if (np == 0) {
                nobj++;
                continue;
            }
            if (ii == 0)
                npmax = imax(np,npmax);
            else {
                if (PMNC == 0) {
                    fprintf(PMF1d,PMNFmtS,++id);
                    fprintf(PMF1d,PMNFmtS,1);
                    fprintf(PMF1d,PMNFmtS,1);
                    fprintf(PMF1d,"%s\n",name);
                    nrec++;
                    fprintf(PMF1d,PMFmtS,AcX[0]);
                    fprintf(PMF1d,PMFmtS,AcY[0]);
                    fprintf(PMF1d,"\n");
                    nrec++;
                    nt1++;
                }
                fprintf(PMF1d,PMNFmtS,++id);
                fprintf(PMF1d,PMNFmtS,3);
                fprintf(PMF1d,PMNFmtS,np - 1);
                fprintf(PMF1d,"%s\n",name);
                nrec++;
                for (i = 1; i < np; ++i) {
                    fprintf(PMF1d,PMFmtS,AcX[i]);
                    fprintf(PMF1d,PMFmtS,AcY[i]);
                    fprintf(PMF1d,"\n");
                    nrec++;
                }
                nt3++;
            }
        }
        prn_message(rec,1,0);
        if (ii == 0) {
            printf1("Number of records: %d\n",rec);
            printf1("Number of data blocks: %d\n",nobj);
            printf1("Number of empty blocks: %d\n",nobje);
            printf1("Number of data points: %d\n",npt);
            printf1("Maximal number of points per object: %d\n",npmax);

            if (alloc_acx(npmax + 1))
                goto SDDCWPFin;
            if (alloc_acy(npmax + 1))
                goto SDDCWPFin;

            if (fseek(PMFd,(long)0,0)) {             
                printf1("Error: cannot seek to begin of file.\n");
                goto SDDCWPFin;
            }
        }
    }
    printf1("%d records written to: %s\n",nrec,PMF1dName);
    printf1("Number of type 1 objects: %d\n",nt1);
    printf1("Number of type 3 objects: %d\n",nt3);

    err = 0;     

SDDCWPFin:
    p_clean();
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

int sdgshhs(void)  
{
    register int i;
   	int	max_east,err,nmax,rec,nrec,nst,np,nps;
    double lon,lat,xmin,xmax,ymin,ymax;
   	struct	POINT p;
   	struct GSHHS h;

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Conversion of GSHHS files. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,10,1)) {     /* get parameters */
        goto SDGSHHSFin;
    }   
    if (PMOPT < 1 || PMOPT > 2)
        PMOPT = 1;

    if (PMLEVEL < 1 || PMLEVEL > 4)
        PMLEVEL = 0;

    if (PMFmtF == 0)
        pmfmt(10,5);

    if (PMF1Def == 0) { 
        printf1("Error: need an output file.\n");
        goto SDGSHHSFin;
    }
    printf1("Input file: %s\n",PMFdName);
    if (PMLEVEL)
        printf1("Level selection: %d\n",PMLEVEL);
    newline();

    if (alloc_acxf(GSHHSNP + 2))
        goto SDGSHHSFin;
    if (alloc_acyf(GSHHSNP + 2))
        goto SDGSHHSFin;

    rec = nrec = nmax = np = nps = 0;
   	max_east = 270000000;
    nst = sizeof(struct POINT);

   	while (fread((void *)&h,(size_t)sizeof(struct GSHHS),(size_t)1,PMFd) == 1) {

        prn_message(++rec,0,0);
        np++;
        if (ARCHTyp == 2) {
          		h.id = swabi4 ((unsigned int)h.id);
          		h.n = swabi4 ((unsigned int)h.n);
          		h.level = swabi4 ((unsigned int)h.level);
          		h.west = swabi4 ((unsigned int)h.west);
          		h.east = swabi4 ((unsigned int)h.east);
          		h.south = swabi4 ((unsigned int)h.south);
          		h.north = swabi4 ((unsigned int)h.north);
          		h.area = swabi4 ((unsigned int)h.area);
          		h.greenwich = swabi2 ((unsigned int)h.greenwich);
          		h.source = swabi2 ((unsigned int)h.source);
        }
        nmax = imax(nmax,h.n);
                 
        if (PMLEVEL && h.level != PMLEVEL) {
            if (fseek(PMFd,(long)(h.n * nst),SEEK_CUR)) {
                printf1("gshhs: seek error.\n");
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
         			if (fread ((void *)&p,(size_t)sizeof(struct POINT),(size_t)1,PMFd) != 1) {
                printf1("Error: while reading the file.\n");
                goto SDGSHHSFin; 
            }
            prn_message(++rec,0,0);

            if (ARCHTyp == 2) {
                p.x = swabi4((unsigned int)p.x);
                p.y = swabi4((unsigned int)p.y);
            }
         			lon = (h.greenwich && p.x > max_east) ? p.x * 1.0e-6 - 360.0 : p.x * 1.0e-6;
         			lat = p.y * 1.0e-6;

            if (PMOPT == 1 && lon > 180.0)        
                lon -= 360.0;       

            AcXF[i] = (float)lon;   
            AcYF[i] = (float)lat;

            if (i == 0) {
                xmin = xmax = lon;
                ymin = ymax = lat;
            }
            else {
                xmin = dmin(xmin,lon);
                xmax = dmax(xmax,lon);
                ymin = dmin(ymin,lat);
                ymax = dmax(ymax,lat);
            }
        }
      		max_east = 180000000;	      /* Only Eurasiafrica needs 270 */
        
        fprintf(PMF1d,"%10d 3 %7d %d ",nps,h.n,h.level);
        fprintf(PMF1d,PMFmtS,xmin);
        fprintf(PMF1d,PMFmtS,xmax);
        fprintf(PMF1d,PMFmtS,ymin);
        fprintf(PMF1d,PMFmtS,ymax);
        fprintf(PMF1d,"%10d %d %d\n",h.area,h.greenwich,h.source);
        nrec++;

      		for (i = 0; i < h.n; ++i) {
            lon = (double)AcXF[i];
            lat = (double)AcYF[i];
            fprintf(PMF1d,PMFmtS,lon);
            fprintf(PMF1d,PMFmtS,lat);
            fprintf(PMF1d,"\n");
            nrec++;
        }
    }
    prn_message(rec,1,0);

    printf1("Number of polygons: %d (max number of points: %d)\n",np,nmax);
    printf1("Number of polygons with selected level: %d\n",nps);
    printf1("%d records written to: %s\n",nrec,PMF1dName);

    if (PMTDAFDef) {                /* create TDA decription file */

        fprintf(PMTDAFd,"sdnvar(\n");
        fprintf(PMTDAFd,"  dfile = %s,\n",PMF1dName);
        fprintf(PMTDAFd,"  noc = %d,\n",nps);
        fprintf(PMTDAFd,"  SDID     [10.0] = c1,\n");
        fprintf(PMTDAFd,"  SDTyp <1>[ 1.0] = c2,\n");
        fprintf(PMTDAFd,"  SDN   <5>[ 7.0] = c3,\n");
        fprintf(PMTDAFd,"  SDPtr <5>[10.0] = rd,\n");
        fprintf(PMTDAFd,"  Level <1>[ 1.0] = c4,\n");
        fprintf(PMTDAFd,"  XMin  <4>[10.5] = c5,\n");
        fprintf(PMTDAFd,"  XMax  <4>[10.5] = c6,\n");
        fprintf(PMTDAFd,"  YMin  <4>[10.5] = c7,\n");
        fprintf(PMTDAFd,"  YMax  <4>[10.5] = c8,\n");
        fprintf(PMTDAFd,"  Area  <4>[10.0] = c9,\n");
        fprintf(PMTDAFd,"  GFlag <1>[ 1.0] = c10,\n");
        fprintf(PMTDAFd,"  Source<1>[ 1.0] = c11,\n");
        fprintf(PMTDAFd,");\n");
        printf1("TDA description written to: %s\n",PMTDAFName);
    }
    err = 0;

SDGSHHSFin:
    p_clean();
    return(err);
}


