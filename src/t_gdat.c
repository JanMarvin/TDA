/****************************************************************************/
/*  t_gdat                                                                  */
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
#include "t_eval.h"  
#include "t_eval1.h"  
#include "t_eval2.h"  
#include "t_eval3.h"  
#include "t_var.h"  
#include "t_zoo.h"  
#include "t_rzoo.h"  
#include "t_gf.h"  
#include "t_sort.h"  
#include "t_mdat.h"  
#include "t_pdat.h"  
#include "t_alloc.h"  
#include "t_rand.h"  

/*  functions in t_gdat.c */

int new_var(void);
int check_vsel(void);               
void prn_dfd(int nv,short *vidx);
int check_matchv(int iv,int idx);
void prn_mval(void);
int check_nv(void);
void prn_nve(char *s);
void prn_nve1(char *s);
void make_vfmt(int nidx);
void make_afmt(void); 
int save_isel(int opt,int idx);
int check_iref(int cnt,int *ptyp,int idx);
int save_vsel(int opt,int idx);
int check_vref(int cnt,int *ptyp,int idx);
int save_bsel(int opt);
int save_break(int opt,int idx);
int check_vcj(int idx,int opt);
int check_ffmt(int opt);
double get_data(int j,int i);
void put_data(double x,int j,int i);
void put_str(char *buf,int blen,int j,int i,int bflag);
void get_str(char *buf,int j,int i);
void clear_str(char *buf,int blen,int j);
double dscan(char *p,int len,int *mval);
int sep_char(char c);
char *skip_sep(char *p);
char *skip_dval(char *p);
int gen_dm(int idx,int nvar,int sflag,int sdflag);
void prn_gdf(int n);
int fnd_match(double *xm,int *n);
int put_t4var(int idx,int r1,int r2,int bn);
void cpy_dmrow(int ir,int ij,int idx);
int tst_bsel(int r1,int r2,int idx);
int tst_vsel(int i,int rec);
int tst_break(int r1,int r2);
int read_dfi(int opt,int sflag);
int read_df(int sflag,int idx);
int skip_df(int n,int opt);         
int get_nsdxy(char *buf,double *x,double *y);
int check_mval(int mval);

int sdnvar(void);
int check_nvsd(void);
int check_sd(int opt);
void sdnvar_close(void);
int sdnvar_alloc(int opt,int n);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int SDVarDef = 0;           /* set to 1 if spatial data defined             */
int SDVarFDef = 0;          /* set to 1 if SDVarFd open                     */
FILE *SDVarFd;              /* file handle                                  */
int SDVarSDID = 0;          /* index of SDID variable                       */
int SDVarSDTyp = 1;         /* index of SDTyp variable                      */
int SDVarSDN = 2;           /* index of SDN variable                        */
int SDVarSDPtr = 3;         /* index of SDPtr variable                      */
int SDVarNP = 0;            /* number of points                             */
int SDVarNL = 0;            /* number of lines                              */
int SDVarNPol = 0;          /* number of polygons                           */
int SDVarNU = 0;            /* number of unknown objects                    */
int SDVarMax = 0;           /* max number of points in object               */
int SDVarNT = 0;            /* total number of points                       */
double SDVarXMin = 0.0;     /* bounding box                                 */
double SDVarXMax = 0.0;
double SDVarYMin = 0.0;
double SDVarYMax = 0.0;

double *SDVarX;             /* standard array for spatial objects           */
double *SDVarY;
int SDVarN = 0;             /* number of allocated entries                  */

/*--------------------------------------------------------------------------*/
#define NVPFmtSLen 10   /* max length of print format strings               */
int DMDef = 0;          /* set to 1 if data matrix is available.            */
int NOCMaxA = 0;        /* max number of cases.                             */
int NOCDM = 0;          /* actual number of cases in data matrix.           */
int NOC = 0;            /* number of data matrix rows to be used, this      */
                        /* depends on the actual tsel command.              */
/*--------------------------------------------------------------------------*/
int NVIdx = 0;          /* start of new variables                           */
int NVNum = 0;          /* number of new variables                          */
int VFmtFlg = 0;        /* set if print format defined with fmt             */
int VFmt1 = 0;     
int VFmt2 = 0;
int VTNOC = 0;          /* set by noc parameter                             */
int BSIZE = 0;          /* set by bsize parameter                           */
/*--------------------------------------------------------------------------*/
char ISEPC = '\0';      /* separation character in data files               */
char *ISelPtr = NULL;   /* pointer to isel expression                       */
short ISESCnt = 0;      /* parser stack for isel expression                 */
int *ISESTyp;
double *ISESVal;

char *VSelPtr = NULL;   /* pointer to vsel expression                       */
short VSESCnt = 0;      /* parser stack for vsel expression                 */
int *VSESTyp;
double *VSESVal;

char *BSelPtr = NULL;   /* pointer to bsel expression                       */
short BSESCnt = 0;      /* parser stack for bsel expression                 */
int *BSESTyp;
double *BSESVal;

char *BRKPtr = NULL;    /* pointer to break expression                      */
short BRKESCnt = 0;     /* parser stack for break expression                */
int *BRKESTyp;
double *BRKESVal;

#define DFILMax 100     /* max number of data files                         */
int DFILN = 0;          /* number of data files                             */
char *DFILFN[DFILMax];  /* pointer to file names                            */
FILE *RFd;              /* file pointer                                     */
int REOF = 1;           /* if file opened                                   */
int DFILNI = 0;         /* data file counter                                */
int DRecLen = 0;        /* set for fixed record length                      */
int DRecLen1 = 0;           
int NMRec = 1;          /* number of multiple records                       */
int NRec = 0;           /* count number of records                          */
int NRec1 = 0;          /* count number of records for each file            */

char *RBuf;             /* read buffer                                      */
int RBufA = 0;          /* if allocated                                     */
char *RSBuf;            /* second read buffer (for strings)                 */
int RSBufA = 0;         /* if allocated                                     */
int RBufLen = 0;        /* length of read buffer                            */
int RBufALen = 0;       /* actual buffer length                             */

int VCJMax = 0;         /* Highest Cj index                                 */
char *VCJFlg;           /* array with flags for Cj required                 */
double *VCJVal;         /* array with values of required Cj                 */

char *FFMTPtr = NULL;   /* pointer to ffmt parameter                        */
int FFMTCnt = 0;        /* number of ffmt terms                             */
int *FFMTC1;            /* first column                                     */
int *FFMTC2;            /* second column                                    */
int FFMTA = 0;          /* if allocated                                     */

int DBLKVar = -1;       /* variable defined by dblock option                */

int MatchNV = 0;        /* number of matching variables                     */
short MatchVA[MaxMatchV];   /* first variables in match command             */
short MatchVB[MaxMatchV];   /* second variables in match command            */

unsigned
char BMsk[8] = {            /* Bitmasks for bit-wise stored variables       */
    0x80, 0x40, 0x20, 0x10, 0x08, 0x04, 0x02, 0x01
};
/*--------------------------------------------------------------------------*/
double MStarVal  = -1.0;    /* default missing value code: star             */
double MPntVal   = -1.0;    /* default missing value code: point            */
double MBlnkVal  = -1.0;    /* default missing value code: blank            */
double MGenVal   = -1.0;    /* general missing value code                   */
double MMatchVal = -3.0;    /* missing value code for no matching           */

int MGenFlg = 0;
int MStarN  = 0;            /* number of mstar missing values               */
int MPntN   = 0;            /* number of mpnt missing values                */
int MBlnkN  = 0;            /* number of mblnk missing values               */
int MGenN   = 0;            /* number of mgen missing values                */
/*--------------------------------------------------------------------------*/
double *AVVAL;              /* temporary storage of archive variables       */
int ARCDICFlg = 0;          /* set by arcdic option                         */
/*--------------------------------------------------------------------------*/
int WIVar = -1;             /* internal variable number of weights          */
double WSum = 0.0;          /* sum of weights                               */
double WSumS = 0.0;         /* factor defined by wnorm option               */
double WNorm = 1.0;         /* sum W(i) * WNOrm = (required) WSum           */
int WNormFlag = 0;          /* 1 if wnorm used in cwt command               */
/*--------------------------------------------------------------------------*/
int TSelFlg = 0;            /* set if tsel command active                   */
int *TSelect;               /* array with select indices                    */
int TSelectA = 0;           /* if allocated                                 */
/*--------------------------------------------------------------------------*/
FILE *GDFd;                 /* output file with df option                   */
int GDFlg = 0;              /* set by df parameter.                         */
int GDKFlg = 0;             /* set by df1 parameter                         */
char *GDFNPtr = NULL;       /* file name                                    */
int GDNRec = 0;             /* number of records written to GDFd            */
char *KeepPtr = NULL;       /* pointer to keep expression                   */
char *DropPtr = NULL;       /* pointer to drop expression                   */
short *NVarIdx;             /* list of new variables for df option          */
int NVarNV = 0;             /* number of variables allocated                */
char **NVPFmtS;             /* format strings                               */
int NVPFmtSA = 0;           /* allocated                                    */
int NVPFmtSA1 = 0;          /* allocated                                    */
char NVSEPC = ' ';          /* sepc=...                                     */
int NVL0 = 0;               /* l0=...                                       */

/* ------------------------------------------------------------------------ */
/*  new_var.    Command: nvar(...)                                          */
/*              Create new variables. If a data matrix already exists,      */
/*              add the new variables. Otherwise create a new data matrix   */
/*              implying that tsel = off and no case weights.               */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int new_var(void)
{
    register int i,j,k;
    int m,n,err,dmdef,nb,w,d,sflag;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Creating new variables. ");
    prn_mem();
    tsel_off(1);        /* turn off temporary case selection */

    dmdef = DMDef;
    DFILN = MatchNV = VFmtFlg = NVNum = 0;
    DRecLen = 0;            /* variable record length */
    NMRec = 1;              /* logical = physical record */
    DBLKVar = -1;           /* no block mode */
    AVDFN = -1;             /* no archive file */
    ARCDICFlg = 0;          /* set by arcdic option */
    RDMNInit = 0;           /* init random number generator */
    ISEPC = '\0';           /* separation character */
    VTNOC = 0;
    BSIZE = 0;

    EDVALSn = EDVALOrg = EDVALDes = 0;
    EDVALTs = EDVALTf = 0.0;
         
    NVIdx = get_nidx();     /* index to first new variable */

    /* set default missing value codes */

    MStarVal  = -1.0; 
    MPntVal   = -1.0; 
    MBlnkVal  = -1.0; 
    MGenVal   = -1.0; 
    MMatchVal = -3.0; 
    GDNRec = GDKFlg = GDFlg = MGenFlg = MStarN = MPntN = MBlnkN = MGenN = 0; 

    ISelPtr  = NULL;        /* set by isel option */
    VSelPtr  = NULL;        /* set by vsel option */
    BSelPtr  = NULL;        /* set by bsel option */
    GDFNPtr  = NULL;        /* file name by df option */
    KeepPtr  = NULL;            
    DropPtr  = NULL;            

    NVArc = NVArc1 = 0;     /* number of archive variables, counted by save_var() */

    if (check_nv())         /* check command etc */
        goto NVFin;

    if (NVNum == 0) {
        printf1("No variables defined.\n");
        goto NVFin;
    }
    newline();      

    if (VLabelLen > 0 && VLabelLen < 8)      
        VLabelLen = 8;

    sflag = 0;
    i = NVIdx;                  /* count string variables */
    while (i >= 0) {
        if (VTyp[i] == 1 && VTypA[i] == 0)
            sflag++;
        i = VNxt[i];
    }

    if (NVArc1 > 0) {
        printf1("Number of archive variables: %d\n",NVArc1);
        if (ARCDef == 0) {
            p_err(-9,1);
            goto NVFin;
        }
        NVArc = NVArc1;
        printf1("Searching for these variables in archive: %s\n",ZOONam);
        if (alloc_avar(NVIdx,NVArc,1)) {
            p_err(-2,1);
            goto NVFin;
        }

        /* check whether all variables belong to same data file */

        AVDFN = check_avar(NVArc,ARCDICFlg);
        if (AVDFN < 0)
            goto NVFin;

        make_afmt();        /* make new print formats for archive variables */
        newline();
    }
    if (VFmtFlg)
        make_vfmt(NVIdx);   /* make new print formats if requested with fmt */

    prn_var(NVIdx);         /* print list of new variables */

    if (DMDef == 0) {
        if (GDFlg == 0)
            NOCMaxA = VTNOC;
        else {
            if (BSIZE > 0)
                NOCMaxA = BSIZE;
            else if (VTNOC > 0)
                NOCMaxA = VTNOC;
            else
                NOCMaxA = VTNOC = NOCDef;
        }
        if (NOCMaxA <= 0) {
            if (AVDFN >= 0 && GDFlg == 0)
                NOCMaxA = ZAFNRec[AVDFN];
            else
                NOCMaxA = NOCDef;
        }
         
        if (GDFlg == 0) {
            printf1("\nCreating a new data matrix.\n");
            printf1("Maximum number of cases: %d\n",NOCMaxA);
        }
        else {
            printf1("\nData will be directly written to: %s\n",GDFNPtr);
            printf1("Will not add variables to data matrix.\n");

            if (!(GDFd = fopen(GDFNPtr,OPEN_WR))) {
                printf1("Error: can't open: %s\n",GDFNPtr);
                goto NVFin;
            }
            GDFlg = 2;
   
            /* create list of variables */

            if (KeepPtr != NULL) {
                printf1("Checking keep=%s\n",KeepPtr); 
                get_nvia(KeepPtr,&n,1,&nb);
                if (n < 1)
                    goto NVFin;
                if (nb) {
                    p_err(-42,1);
                    goto NVFin;
                }

                if (!(NVarIdx = (short *)calloc(n,sizeof(short)))) {
                    p_err(-2,1);
                    goto NVFin;
                }
                memrq(n,sizeof(short));
                NVarNV = n;       
                for (i = 0; i < n; ++i)
                    NVarIdx[i] = VLVIdx[i];
            }
            else {
                if (!(NVarIdx = (short *)calloc(NVNum,sizeof(short)))) {
                    p_err(-2,1);
                    goto NVFin;
                }
                memrq(NVNum,sizeof(short));
                NVarNV = NVNum;

                j = 0;
                i = NVIdx;          
                while (i >= 0) {
                    NVarIdx[j++] = i;
                    i = VNxt[i];
                }
                if (DropPtr != NULL) {
                    printf1("Checking drop=%s\n",DropPtr); 
                    get_nvia(DropPtr,&n,1,&nb);
                    if (n < 1)
                        goto NVFin;
                    if (nb) {
                        p_err(-42,1);
                        goto NVFin;
                    }

                    m = 0;
                    for (i = 0; i < NVarNV; ++i) {
                        k = NVarIdx[i];
                        for (j = 0; j < n; ++j) {
                            if (k == VLVIdx[j]) {
                                NVarIdx[i] = k = -1;
                                break;
                            }
                        }
                        if (k >= 0)
                            m++;
                    }
                    if (m == 0) {
                        printf1("Error: all variables dropped.\n");
                        goto NVFin;
                    }
                }
            }
                                    /* make print format */
            if (NVarNV > 0) {
                if (!(NVPFmtS = (char **)calloc(NVarNV,sizeof(char *)))) {
                    p_err(-2,1);             
                    goto NVFin;
                }         
                memrq(NVarNV,sizeof(char *));
                NVPFmtSA = NVarNV;

                if (NVL0 != 0)
                    NVL0 = -1;

                for (i = 0; i < NVarNV; ++i) {  

                    if (!(NVPFmtS[i] = (char *)calloc(NVPFmtSLen,sizeof(char)))) {
                        p_err(-2,1);             
                        goto NVFin; 
                    }         
                    memrq(NVPFmtSLen,sizeof(char));
                    NVPFmtSA1++;  

                    k = NVarIdx[i];
                    w = (int)VPFmt1[k];             
                    d = (int)VPFmt2[k];             
                    if (VTyp[k] != 1)
                        makefmt(&w,&d,NVPFmtS[i],0,NVSEPC,NVL0);
                }
            }
        }
    }
    else if (GDFlg) {
        printf1("\nError: df parameter cannot be used if a data matrix exists.\n");
        goto NVFin;
    }

    if (ISelPtr != NULL) {              /* isel option */
        if (check_vsel())
            goto NVFin;
        if (save_isel(1,NVIdx))
            goto NVFin;
    }
    if (VSelPtr != NULL) {              /* vsel option */
        if (check_vsel())
            goto NVFin;

        if (DBLKVar >= 0)  
            printf1("Block mode: vsel will be ignored.\n");
        else {
            if (save_vsel(1,NVIdx))
                goto NVFin;
        }
    }
    if (BRKPtr != NULL) {               /* break option */
        if (check_vsel())
            goto NVFin;
        if (DBLKVar >= 0)  
            printf1("Block mode: break will be ignored.\n");
        else {
            if (save_break(1,NVIdx))
                goto NVFin;
        }
    }
    if (check_vcj(NVIdx,1))             /* check cj references */
        goto NVFin;

    if (VCJMax > 0 || sflag) {          /* need a data file */

        if (NVArc > 0) {
            printf1("Error: referring to a data file with c terms or string\n");
            printf1("variables is not compatible with using archive variables.\n");
            goto NVFin;
        }

        if (DFILN > 0) {
            printf1("\nUsing data file(s): %s",DFILFN[0]);
            for (i = 1; i < DFILN; ++i)
                printf1(",%s",DFILFN[i]);
            printf1("\n");
            if (DRecLen > 0)
                printf1("Fixed record length: %d\n",DRecLen);
            if (NMRec > 1)
                printf1("Logical records consist of %d physical records.\n",NMRec);
        }
        else {
            printf1("\nError: need a data file.\n");
            goto NVFin;
        }
        if (FFMTCnt > 0) {
            if (check_ffmt(1))
                goto NVFin;
        }
        else {
            printf1("Free format. Separation character(s): ");
            if (ISEPC != '\0')
                printf1("%04x [hex]\n",ISEPC);
            else 
                printf1("default.\n");
        }
        if (VTNOC > 0)  
            printf1("Reading maximal %d records.\n",VTNOC);

    }
    else if (DFILN > 0) {
        printf1("\nWarning: dfile parameter(s) will be ignored.\n");
        DFILN = 0;
    }
    if (DMDef) {            /* if a data matrix already exists */
   
        printf1("\nNew variables will be added to existing data matrix.\n");

        if (MatchNV > 0) {
            if (NVArc == 0 && DFILN == 0) {
                printf1("Always trivial matching without data files.\n");
                MatchNV = 0;
            }
            else {
                printf1("\nMatching:\nNew variable    Existing variable\n");
                prnchar('-',33,1);
                     
                for (i = 0; i < MatchNV; ++i) {
                    j = MatchVA[i];
                    printf1("%s",VName[j]);
                    prnchar(' ',16 - strlen(VName[j]),0);
                    printf1("%s\n",VName[MatchVB[i]]);
                    if (VTyp[j] != 2 && VTyp[j] != 3) {
                        printf1("\nError: %s is not of type 2 or 3.\n",VName[j]);
                        goto NVFin;
                    }
                    if (check_matchv(MatchVB[i],NVIdx)) {
                        printf1("\nError: %s is not an existing variable.\n",VName[MatchVB[i]]);
                        goto NVFin;
                    }
                    j = MatchVB[i];
                    if (VTyp[j] == 1) {
                        printf1("\nError: %s must not be of type 1 (string variable).\n",VName[j]);
                        goto NVFin;
                    }
                }
                printf1("\n");
                if (DBLKVar >= 0) {
                    printf1("Error: dblock and match not compatible.\n");
                    goto NVFin;
                }
                if (BRKPtr != NULL) {
                    printf1("Break will be ignored.\n");
                    save_break(0,0);
                }
                if (vsort(MatchNV,MatchVB,1,1,1))       /* sort */
                    goto NVFin;
            }
        }
        else
            printf1("Trivial matching.\n");
    }
    else if (MatchNV > 0) {
        printf1("Match command will be ignored.\n");
        MatchNV = 0;
    }
            
    if (DBLKVar >= 0) {
        printf1("\nBlock mode defined by variable: %s\n",VName[DBLKVar]);
        i = VTyp[DBLKVar];
        if (DBLKVar >= NVIdx && i != 2 && i != 3) {
            printf1("Error: type of block mode variable must be 2 or 3.\n");
            goto NVFin;
        }
        if (GDFlg)
            printf1("Maximum block size: %d\n",NOCMaxA);

        if (BSelPtr != NULL) {              /* bsel option */
            if (check_vsel())
                goto NVFin;

            if (save_bsel(1))
                goto NVFin;
        }
    }
    else if (BSelPtr != NULL) {
        printf1("Block select ignored.\n");
        BSelPtr = NULL;
    }

    if (GDFlg && DBLKVar < 0) {     /* check type of new variables */

        i = NVIdx;
        while (i >= 0) {
            if (VTyp[i] > 3) {
                printf1("Error: if not in block mode, df option allows only variable types 1,2,3.\n");
                goto NVFin;
            }
            i = VNxt[i];
        }
    }

    if (NVArc > 0 || DFILN > 0)
        printf1("\n");         
  
    /* allocate memory for new variables */

    if (alloc_vdat(NVIdx,1)) {
        printf1("Insufficient memory for new variables.\n");
        goto NVFin;
    }
    err = gen_dm(NVIdx,NVNum,sflag,0);      /* create data */
    if (err == 0) {
        prn_mval();                 /* info about missing values */
    
        if (GDFlg && NVNum > 0 && (PMTDAFDef || PMSPSSFDef))
            prn_dfd(NVarNV,NVarIdx);
    }

NVFin:
    if (GDFlg > 1)                 
        fclose(GDFd);

    if (NVarNV > 0) {
        free((char *)NVarIdx);
        memrq(-NVarNV,sizeof(short));
        NVarNV = 0;
    }
    if (NVPFmtSA1 > 0) {
        for (i = 0; i < NVPFmtSA1; ++i) {  
            free((char *)NVPFmtS[i]); 
            memrq(-NVPFmtSLen,sizeof(char));
        }
        NVPFmtSA1 = 0;
    }
    if (NVPFmtSA > 0) {
        free((char *)NVPFmtS); 
        memrq(-NVPFmtSA,sizeof(char *));
        NVPFmtSA = 0;
    }
    alloc_vl(0);

    save_isel(0,0);
    save_vsel(0,0);
    save_bsel(0);
    save_break(0,0);
    check_vcj(NVIdx,0);
    check_ffmt(0);
    if (NVArc)
        alloc_avar(NVIdx,NVArc,0);      
    vsort(MatchNV,MatchVB,0,0,1);

    if (err || GDFlg) {
        if (NVNum > 0)
            clear_avar(NVIdx);
        printf1("No new variables created. ");
    }
    else {
        if (dmdef == 0) {
            WIVar = -1;
            WSum = (double)NOC;
            WSumS = 0.0;
            WNorm = 1.0;
            WNormFlag = 0;
        }
        printf1("\nEnd of creating new variables. ");
    }
    GDNRec = GDFlg = NVArc = 0;
    prn_mem();
    /* prnchar('-',LLEN,1); */
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_vsel()    Check vsel, bsel and break                              */ 
/*                  Return 0 if OK, -1 if error.                            */

int check_vsel(void)                
{
    if (NVArc <= 0 && DFILN <= 0) {
        printf1("\nError: isel, vsel, bsel and break can only be used when reading\n"); 
        printf1("from a data file or archive.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_dfd(nv,vidx)    Print TDA and SPSS description files.               */

void prn_dfd(int nv,short *vidx)
{
    register int i,j,k;

    j = 0;
    for (i = 0; i < nv; ++i) {
        k = vidx[i];
        if (k >= 0)
            vidx[j++] = k;
    }
    newline();
    if (PMTDAFDef)  
        dtda(GDFNPtr,GDNRec,j,vidx,1);
    if (PMSPSSFDef)  
        dspss(GDFNPtr,j,vidx,1);
    newline();
}

/* ------------------------------------------------------------------------ */
/*  check_matchv(iv,idx)    Check whether variable iv is a new variable.    */
/*                          Return 1 if yes, else 0.                        */

int check_matchv(int iv,int idx)
{
    while (idx >= 0) {
        if (iv == idx)
            return(1);
        idx = VNxt[idx];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_mval    Info about missing values.                                  */

void prn_mval(void)
{
    int i,n,k,l,ll;

    if (NVArc > 0) {

        printf1("\nMissing values in numerical variables");
        n = 0;
        for (i = 0; i < NVArc; ++i) {
            if (AVMBlnk[i] || AVMStar[i] || AVMPnt[i] || AVMGen[i]) {
                if (n == 0) {
                    ll = 8;
                    for (k = 0; k < NVArc; ++k) {
                        l = strlen(VDef[AVIdx[k]]) - 2;
                        if (ll < l)
                            ll = l;
                    }
                    printf1("\nVariable  ");
                    prnchar(' ',ll - 8,0);
                    printf1("Blanks  Stars  Points  General\n");
                    prnchar('-',32 + ll,1);
                }
                printf1("%s  ",VDef[AVIdx[i]] + 2);
                prnchar(' ',ll - strlen(VDef[AVIdx[i]]) + 2,0);
                printf1("%6d ",AVMBlnk[i]); 
                printf1("%6d ",AVMStar[i]); 
                printf1("%7d ",AVMPnt[i]); 
                printf1("%8d\n",AVMGen[i]);
                n++;
            }
        }
        if (n == 0)  
            printf1(" (blank,star,point,general): none\n");
        else  
            printf1("\nSubstitution: MBlnk=%g MPnt=%g MStar=%g MGen=%g\n",
                                MBlnkVal,MPntVal,MStarVal,MGenVal);
    }
    else if (DFILN > 0) {
   
        if (MStarN + MPntN + MBlnkN + MGenN == 0) {
            printf1("Missing values in data file(s): none.\n");
            return;
        }
        printf1("\nMissing values  numerical code\n");
        prnchar('-',30,1);        
        printf1("Blank %8d      %lg\n",MBlnkN,MBlnkVal);
        printf1("Star %9d      %lg\n",MStarN,MStarVal);
        printf1("Point %8d      %lg\n",MPntN,MPntVal);
        if (MGenFlg)
            printf1("General %6d      %lg\n",MGenN,MGenVal);
    }
}

/* ------------------------------------------------------------------------ */
/*  check_nv    Check nvar command in CmdBuf.                               */
/*              Return 0 if OK, -1 if error.                                */

int check_nv(void)
{
    int err,r,i,n,l,m1,m2,nb;
    register char c,*p,*q;
    char sc,vname[VNLMax + 1];
    double x;

    NVSEPC = ' ';             
    NVL0 = VTNOC = BSIZE = 0;
    err = -1;
    p = CmdBuf + 4;
    while (*++p) {
        if (sscanf(p,"noc=%d",&n) == 1 && n > 0) {
            VTNOC = n;
            p = skip_int(p + 4);
        }
        else if (sscanf(p,"bsize=%d",&n) == 1 && n > 0) {
            BSIZE = n;
            p = skip_int(p + 6);
        }
        else if (sscanf(p,"dreclen=%d",&n) == 1 && n > 0) {
            DRecLen = n;
            p = skip_int(p + 8);
        }
        else if (sscanf(p,"nmrec=%d",&n) == 1 && n > 0) {
            NMRec = n;
            p = skip_int(p + 6);
        }
        else if (sscanf(p,"fmt=%d.%d",&VFmt1,&VFmt2) == 2) {
            VFmtFlg = 1;
            p = skip_int(p + 4);
            p = skip_int(p + 1);
        }
        else if (sscanf(p,"mstar=%lf",&x) == 1) {
            MStarVal = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"mpnt=%lf",&x) == 1) {
            MPntVal = x;
            p = skip_dbl(p + 5);
        }
        else if (sscanf(p,"mblnk=%lf",&x) == 1) {
            MBlnkVal = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"mgen=%lf",&x) == 1) {
            MGenVal = x;
            MGenFlg = 1;
            p = skip_dbl(p + 5);
        }
        else if (sscanf(p,"mmatch=%lf",&x) == 1) {
            MMatchVal = x;
            p = skip_dbl(p + 7);
        }
        else if (!strncmp(p,"arcdic",6)) {
            ARCDICFlg = 1;
            p += 6;
        }
        else if (!strncmp(p,"dfile=",6)) {
            if (DFILN >= DFILMax) {
                printf1("Exceeded max number of data files.\n");
                goto CNVFin;
            }
            DFILFN[DFILN++] = p + 6;
            p = skip_nc(p);
        }
        else if (!strncmp(p,"isel=",5)) {
            ISelPtr = p + 5;
            p = skip_expr(p + 5);
        }
        else if (!strncmp(p,"vsel=",5)) {
            VSelPtr = p + 5;
            p = skip_expr(p + 5);
        }
        else if (!strncmp(p,"bsel=",5)) {
            BSelPtr = p + 5;
            p = skip_expr(p + 5);
        }
        else if (!strncmp(p,"break=",6)) {
            BRKPtr = p + 6;
            p = skip_expr(p + 6);
        }
        else if (!strncmp(p,"df=",3)) {
            GDFNPtr = p + 3;
            GDFlg = 1;
            p = skip_nc(p + 3);
        }
        else if (!strncmp(p,"df1=",4)) {
            GDFNPtr = p + 4;
            GDKFlg = GDFlg = 1;
            p = skip_nc(p + 4);
        }
        else if (!strncmp(p,"keep=",5)) {
            KeepPtr = p + 5;
            p = get_nvi(p + 5,&n,0,VLVIdx,&nb);
            if (n < 1 || nb) {
                printf1("Error in keep expression.\n");
                goto CNVFin;
            }
        }
        else if (!strncmp(p,"drop=",5)) {
            DropPtr = p + 5;
            p = get_nvi(p + 5,&n,0,VLVIdx,&nb);
            if (n < 1 || nb) {
                printf1("Error in drop expression.\n");
                goto CNVFin;
            }
        }
        else if (!strncmp(p,"dtda",4)) {    /* PMTDA file */
            p = get_file(10,p + 4,&r);
            if (r)  
                 goto CNVFin;
        }
        else if (!strncmp(p,"dspss",5)) {   /* PMSPSS file */
            p = get_file(11,p + 5,&r);
            if (r)  
                goto CNVFin;
        }
        else if (!strncmp(p,"match=",6)) {
            q = p + 6;
            n = 0;
            for (i = 0; i < MaxMatchV; ++i) {
                if ((MatchVA[i] = get_vidx(q,vname)) < 0) {
                    prn_nve1(p);
                    goto CNVFin;
                }
                q += strlen(vname);

                if (*q++ != ',' || (MatchVB[i] = get_vidx(q,vname)) < 0) {
                    prn_nve1(p);
                    goto CNVFin;
                }
                q += strlen(vname);
                n++;       

                if (*q != ',' || check_vname(q + 1) == 0)
                    break;
                q++;
            }
            MatchNV = n;
            p = q;
        }
        else if (sscanf(p,"isc=%c",&sc) == 1) {
            if (sc == 't')
                sc = '\t';
            ISEPC = sc;
            p += 5;                
        }
        else if (!strncmp(p,"dblock=",7)) {
            DBLKVar = get_vidx(p + 7,vname);
            if (DBLKVar < 0) {
                prn_nve1(p);
                goto CNVFin;
            }
            p += 7 + strlen(vname);
        }
        else if (!strncmp(p,"ffmt=c",6)) {
            FFMTCnt = 0;       
            q = p + 5;
            while (1) {
                if ((sscanf(q,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d,%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)) {
                    q = skip_int(q + 1);
                    q = skip_blev(q);
                }
                else {
                    prn_nve(p);
                    goto CNVFin;
                }
                FFMTCnt++;
                if (*q != ',' || *(q + 1) != 'c')
                    break;
                q++;       
            }
            FFMTPtr = p + 5;
            p = q;
        }
        else if ((l = get_vnlen(p)) > 0) {       /* variable */
            q = skip_nc(p + l);
            c = *q;
            *q = '\0';
            if (save_var(p,0))
                goto CNVFin;
            NVNum++;
            *q = c;
            p = q; 
        }
        else if (!strncmp(p,"sepc=none",9)) {
            NVSEPC = '\0';
            p += 9;
        }
        else if (sscanf(p,"sepc=%c",&NVSEPC) == 1) {
            if (NVSEPC == 't')
                NVSEPC = '\t';
            p += 6;                
        }
        else if (sscanf(p,"l0=%d",&n) == 1) {
            NVL0 = n;
            p = skip_int(p + 3);
        }
        if (*p != ',' && *p != ')') {
            prn_nve(p);
            goto CNVFin;
        }
        *p = '\0';
    }
    err = 0;

CNVFin:
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_nve.    Print error message.                                        */

void prn_nve(char *s)
{
    register char *p = s;

    printf1("Syntax error: ");
    if (!*p)
        printf1("check brackets and semicolon.\n");
    else {     
        while (*p && p < s + 20)
            printf1("%c",*p++);
        if (*p)
            printf1(" ...");
        printf1("\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_nve1    Print error message.                                        */

void prn_nve1(char *s)
{
    register char *p = s;

    printf1("Syntax error or undefined variable: ");
    while (*p && p < s + 30)
        printf1("%c",*p++);
    if (*p)
        printf1(" ...");
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  make_vfmt(nidx)   Make new print format, beginning with nidx.           */
/*                    Only for variables with 0.0, not for string var.      */
   
void make_vfmt(int nidx)
{
    int i,w1,w2;

    i = nidx;   
    while (i >= 0) {              
           
        w1 = (int)VPFmt1[i];
        w2 = (int)VPFmt2[i];
                 
        if (w1 == 0 && w2 == 0 && VTyp[i] != 1) {
            w1 = VFmt1;
            w2 = VFmt2;
            makefmt(&w1,&w2,VPFmtS[i],0,' ',0);
            VPFmt1[i] = (short)w1;
            VPFmt2[i] = (short)w2;
        }
        i = VNxt[i];
    }
}

/* ------------------------------------------------------------------------ */
/*  make_afmt   Make new print format for archive variables.                */
/*              Only if current format is 0.0.                              */

void make_afmt(void)  
{
    int i,j,s,w1,w2;
   
    for (i = 0; i < NVArc; ++i) {
        j = AVIdx[i];
        if (VTyp[j] == 1)
            continue;

        if (VPFmt1[j] == 0 && VPFmt2[j] == 0) {
            w1 = AVFmt1[i];
            w2 = AVFmt2[i];
            s = get_afmt(&w1,&w2);

            makefmt(&w1,&w2,VPFmtS[j],0,' ',0);
            VPFmt1[j] = (short)w1;
            VPFmt2[j] = (short)w2;
            /** VSLen[j] = (short)s; **/
        }
    }
}
     
/* ------------------------------------------------------------------------ */
/*  save_isel(opt,idx)                                                      */
/*                                                                          */
/*                  If opt != 0 check whether the isel expression in        */
/*                  ISelPtr could be correctly parsed. If this is the case, */
/*                  allocate a stack in                                     */
/*                  ISESTyp[], ISESVal[], ISESCnt and copy the parser stack */
/*                  to this stack. Return 0.                                */
/*                  If expression cannot be correctly parsed, don't alloc   */
/*                  the stack and return -1. Also if expression contains    */
/*                  type 2 or 3 operators.                                  */
/*                                                                          */
/*                  Restriction: the isel expression may not refer          */
/*                  to any new variables beginning at index idx.            */
/*                                                                          */
/*                  If opt == 0 free the previously allocated stack and     */
/*                  set ISESCnt = 0.                                        */

int save_isel(int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (ISESCnt > 0) {
            free((char *)ISESTyp);
            free((char *)ISESVal);
            memrq(-ISESCnt,sizeof(double) + sizeof(int));
        }
        ISESCnt = 0;
        ISelPtr = NULL;
        return(0);
    }
    printf1("Input select (isel): %s\n",ISelPtr);

    if ((n = v_parse(ISelPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d) in isel expression.\n",n);
        if (n < 0)
            prn_emsg1(n);
        return(-1);
    }
       
    /* check for type 2 and 3 operators */

    check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,0);
    if (n2 > 0 || n3 > 0) {
        printf1("Error: isel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_iref(ESCnt,ESTyp,idx)) {
        printf1("Reference error in isel expression.\n");
        return(-1);
    }
    if (!(ISESTyp = (int *)calloc(ESCnt,sizeof(int)))) {
        p_err(-2,1);
        return(-1);
    }
    if (!(ISESVal = (double *)calloc(ESCnt,sizeof(double)))) {
        free((char *)ISESTyp);
        p_err(-2,1);
        return(-1);      
    }
    memrq(ESCnt,sizeof(double) + sizeof(int));

    ISESCnt = ESCnt;
    for (i = 0; i < ESCnt; ++i) {
        ISESTyp[i] = ESTyp[i];
        ISESVal[i] = ESVal[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_iref(cnt,ptyp,idx)                                                */
/*                                                                          */
/*      Check whether variable with parser stack cnt,ptyp refers to         */
/*      another variable in the variable list beginning at index idx;       */
/*      only recognized if not VTypA = 2 (archive variables)                */
/*                                                                          */
/*      Return 0 if OK                                                      */
/*      -1 if reference error                                               */
         
int check_iref(int cnt,int *ptyp,int idx)
{
    register int i,k,t;
                 
    for (i = 0; i < cnt; ++i) {
        t = iabs(ptyp[i]);
        if (t >= VOFFS && t < COFFS) {
            k = idx;
            while (k >= 0) {
                if (t == VOFFS + k && VTypA[k] != 2)
                    return(-1);
                k = VNxt[k];
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  save_vsel(opt,idx)                                                      */
/*                                                                          */
/*                  If opt != 0 check whether the vsel expression in        */
/*                  VSelPtr could be correctly parsed. If this is           */
/*                  the case allocate a stack in                            */
/*                  VSESTyp[], VSESVal[], VSESCnt and copy the parser stack */
/*                  to this stack. Return 0.                                */
/*                  If expression cannot be correctly parsed, don't alloc   */
/*                  the stack and return -1. Also if the expression         */
/*                  contains type 2 or 3 operators.                         */
/*                                                                          */
/*                  Restriction: the vsel expression may not refer          */
/*                  to any new variables beginning at index idx that        */
/*                  contain type 2 operators.                               */
/*                                                                          */
/*                  If opt == 0 free the previously allocated stack and     */
/*                  set VSESCnt = 0.                                        */

int save_vsel(int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (VSESCnt > 0) {
            free((char *)VSESTyp);
            free((char *)VSESVal);
            memrq(-VSESCnt,sizeof(double) + sizeof(int));
        }
        VSESCnt = 0;
        VSelPtr = NULL;
        return(0);
    }
    printf1("Input select (vsel): %s\n",VSelPtr);

    if ((n = v_parse(VSelPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d) in vsel expression.\n",n);
        if (n < 0)
            prn_emsg1(n);
        return(-1);
    }

    /* check for type 2 and 3 operators */

    check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,0);
    if (n2 > 0 || n3 > 0) {
        printf1("Error: vsel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_vref(ESCnt,ESTyp,idx)) {
        printf1("Error: vsel expression may only refer to variables of type 1 - 3.\n");
        return(-1);
    }
    if (!(VSESTyp = (int *)calloc(ESCnt,sizeof(int)))) {
        p_err(-2,1);
        return(-1);
    }
    if (!(VSESVal = (double *)calloc(ESCnt,sizeof(double)))) {
        free((char *)VSESTyp);
        p_err(-2,1);
        return(-1);
    }
    memrq(ESCnt,sizeof(double) + sizeof(int));

    VSESCnt = ESCnt;
    for (i = 0; i < ESCnt; ++i) {
        VSESTyp[i] = ESTyp[i];
        VSESVal[i] = ESVal[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_vref(cnt,ptyp,idx)                                                */
/*                                                                          */
/*      Check whether variable with parser stack cnt,ptyp refers to         */
/*      another variable in the variable list, beginning at index idx,      */
/*      containing type 2 operators.                                        */
/*                                                                          */
/*      Return 0 if OK                                                      */
/*      -1 if there are variables of type 4 or higher                       */

int check_vref(int cnt,int *ptyp,int idx)
{
    register int i,k,t;
                 
    for (i = 0; i < cnt; ++i) {

        t = iabs(ptyp[i]);

        if (t >= VOFFS && t < COFFS) {
            k = idx;
            while (k >= 0) {
                if (t == VOFFS + k && VTyp[k] >= 4)
                    return(-1);
                k = VNxt[k];
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  save_bsel(opt)                                                          */
/*                                                                          */
/*                  If opt != 0 check whether the bsel expression in        */
/*                  BSelPtr could be correctly parsed. If this is           */
/*                  the case allocate a stack in                            */
/*                  BSESTyp[], BSESVal[], BSESCnt and copy the parser stack */
/*                  to this stack. Return 0.                                */
/*                  If expression cannot be correctly parsed, don't alloc   */
/*                  the stack and return -1.                                */
/*                                                                          */
/*                  Note: expression must not contain type 2 or 3 operators */
/*                  but may refer to type 4 variables.                      */
/*                                                                          */
/*                  If opt == 0 free the previously allocated stack and     */
/*                  set BSESCnt = 0.                                        */

int save_bsel(int opt)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (BSESCnt > 0) {
            free((char *)BSESTyp);
            free((char *)BSESVal);
            memrq(-BSESCnt,sizeof(double) + sizeof(int));
        }
        BSESCnt = 0;
        BSelPtr = NULL;
        return(0);
    }
    printf1("Block select: %s\n",BSelPtr);

    if ((n = v_parse(BSelPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d) in vsel expression.\n",n);
        if (n < 0)
            prn_emsg1(n);
        return(-1);
    }

    /* check for type 2 and 3 operators and c terms */

    check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,1);
    if (n2 > 0 || n3 > 0) {
        printf1("Error: bsel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (nc > 0) {
        printf1("Error: bsel expressions must not refer to c terms.\n");
        return(-1);
    }
    if (!(BSESTyp = (int *)calloc(ESCnt,sizeof(int)))) {
        p_err(-2,1);
        return(-1);
    }
    if (!(BSESVal = (double *)calloc(ESCnt,sizeof(double)))) {
        free((char *)BSESTyp);
        p_err(-2,1);
        return(-1);
    }
    memrq(ESCnt,sizeof(double) + sizeof(int));

    BSESCnt = ESCnt;
    for (i = 0; i < ESCnt; ++i) {
        BSESTyp[i] = ESTyp[i];
        BSESVal[i] = ESVal[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  save_break(opt,idx)                                                     */
/*                                                                          */
/*                  If opt != 0 check whether the break expression in       */
/*                  BRKPtr could be correctly parsed. If this is            */
/*                  the case allocate a stack in                            */
/*                  BRKESTyp[], BRKESVal[], BRKESCnt and copy the parser    */
/*                  stack to this stack. Return 0.                          */
/*                  If expression cannot be correctly parsed, don't alloc   */
/*                  the stack and return -1. Also if the expression         */
/*                  contains type 2 or 3 operators.                         */
/*                                                                          */
/*                  Restriction: the break expression may not refer         */
/*                  to any new variables beginning at index idx that        */
/*                  contain type 2 operators.                               */
/*                                                                          */
/*                  If opt == 0 free the previously allocated stack and     */
/*                  set BRKESCnt = 0.                                       */

int save_break(int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (BRKESCnt > 0) {
            free((char *)BRKESTyp);
            free((char *)BRKESVal);
            memrq(-BRKESCnt,sizeof(double) + sizeof(int));
        }
        BRKESCnt = 0;
        BRKPtr = NULL;
        return(0);
    }
    printf1("Break expression: %s\n",BRKPtr);

    if ((n = v_parse(BRKPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d) in break expression.\n",n);
        if (n < 0)
            prn_emsg1(n);
        return(-1);
    }

    /* check for type 2 and 3 operators */

    check_expr(ESCnt,ESTyp,&nv,&nc,&n2,&n3,1);
    if (n2 > 0 || n3 > 0) {
        printf1("Error: break expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_vref(ESCnt,ESTyp,idx)) {
        printf1("Error: break expression may only refer to variables of type 1 - 3.\n");
        return(-1);
    }
    if (!(BRKESTyp = (int *)calloc(ESCnt,sizeof(int)))) {
        p_err(-2,1);
        return(-1);
    }
    if (!(BRKESVal = (double *)calloc(ESCnt,sizeof(double)))) {
        free((char *)BRKESTyp);
        p_err(-2,1);
        return(-1);
    }
    memrq(ESCnt,sizeof(double) + sizeof(int));

    BRKESCnt = ESCnt;
    for (i = 0; i < ESCnt; ++i) {
        BRKESTyp[i] = ESTyp[i];
        BRKESVal[i] = ESVal[i];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_vcj(idx,opt)      Check new variables, isel, vsel, break for      */
/*                          cj references. If opt != 0 create array         */
/*                          VCJFlg[j] = 1 for used cj references,           */
/*                          i = 1,...,VCJMax. Set VCJMax = highest cj index */
/*                          If opt == 0 free previously allocated           */
/*                          memory. Return -1 if error.                     */

int check_vcj(int idx,int opt)
{
    register int i,j,k;   
    int ii,n,cmax;

    if (opt == 0) {
        if (VCJMax > 0) {
            free((char *)VCJFlg);
            free((char *)VCJVal);
            memrq(-VCJMax - 1,sizeof(char) + sizeof(double));
            VCJMax = 0;
        }
        return(0);
    }

    cmax = ii = 0;
VCJREP:
    i = idx;
    while (i >= 0) {
        n = VESCnt[i];
        for (j = 0; j < n; ++j) {
            k = VESTyp[i][j];
            if (k >= COFFS && k < COFFMAX) {
                k -= COFFS;
                if (ii)
                    VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
        i = VNxt[i];
    }
    if (ISelPtr != NULL) {
        n = ISESCnt;
        for (j = 0; j < n; ++j) {
            k = ISESTyp[j];
            if (k >= COFFS && k < COFFMAX) {
                k -= COFFS;
                if (ii)
                    VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (VSelPtr != NULL) {
        n = VSESCnt;
        for (j = 0; j < n; ++j) {
            k = VSESTyp[j];
            if (k >= COFFS && k < COFFMAX) {
                k -= COFFS;
                if (ii)
                    VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (BRKPtr != NULL) {
        n = BRKESCnt;
        for (j = 0; j < n; ++j) {
            k = BRKESTyp[j];
            if (k >= COFFS && k < COFFMAX) {
                k -= COFFS;
                if (ii)
                    VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (ii || cmax == 0)
        return(0);

    VCJMax = cmax;
    if (!(VCJFlg = (char *)calloc(cmax + 1,sizeof(char)))) {
        p_err(-2,1);
        return(-1);      
    }
    if (!(VCJVal = (double *)calloc(cmax + 1,sizeof(double)))) {
        free((char *)VCJFlg);
        p_err(-2,1);
        return(-1);      
    }
    memrq(cmax + 1,sizeof(char) + sizeof(double));
    ii++;
    goto VCJREP;
}

/* ------------------------------------------------------------------------ */
/*  check_ffmt(opt)    Create arrays FFMTC1[i] and FFMTC2[i] for ffmt terms */
/*                     i = 1,...,VCJMax. Also, if opt = 1, print info       */
/*                     to standard output.                                  */
/*                     If opt == 0 free previously allocated memory.        */
/*                     Return 0 if OK, -1 if error.                         */

int check_ffmt(int opt)
{
    register int i;   
    int n,m1,m2,err; 
    register char *p;

    if (opt == 0) {
        if (FFMTA > 0) {
            free((char *)FFMTC1);
            free((char *)FFMTC2);
            memrq(-FFMTA,sizeof(int));
            FFMTCnt = FFMTA = 0;
            FFMTPtr = NULL;
        }
        return(0);
    }
    if (!(FFMTC1 = (int *)calloc(VCJMax + 1,sizeof(int)))) {
        p_err(-2,1);
        return(-1);      
    }
    if (!(FFMTC2 = (int *)calloc(VCJMax + 1,sizeof(int)))) {
        free((char *)FFMTC1);
        p_err(-2,1);
        return(-1);      
    }
    FFMTA = 2 * VCJMax + 2;
    memrq(FFMTA,sizeof(int));

    p = FFMTPtr;
             
    while (1) {
        if ((sscanf(p,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1))   
            ;
        else if (sscanf(p,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)    
            m2 = m1;
        else 
            goto CFFMTErr;

        if (n < 1)
            goto CFFMTErr;
        else if (n <= VCJMax) {
            FFMTC1[n] = m1;
            FFMTC2[n] = m2;
        }
        p = skip_int(p + 1);
        p = skip_blev(p);

        if (*p != ',' || *(p + 1) != 'c')
            break;
        p++;       
    }
    err = 0;
    printf1("\nUsing fixed format data\n");
    prnchar('-',23,1);
    for (i = 1; i <= VCJMax; ++i) {
        if (VCJFlg[i] == 0)
            continue;
        printf1("c%-4d ",i);
        if (FFMTC1[i] > 0) {
            printf1("column%4d",FFMTC1[i]);
            if (FFMTC2[i] > FFMTC1[i])
                printf1(" -%4d",FFMTC2[i]);
            printf1("\n");
        }
        else {
            printf1("error, need column information\n");
            err = -1;
        }
    }
    newline();
    return(err);

CFFMTErr:
    printf1("Error in ffmt parameter.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  get_data(j,i)       Return the value of the j.th variable for case i.   */
/*                      i = 0,...,NOC - 1.                                  */
/*                      Note that the data matrix contains NOCDM cases,     */
/*                      NOC is the number of cases after tsel selection.    */
/*                                                                          */
/*                      If TSelFlg == 1, apply tsel case selection.         */
          
double get_data(int j,int i)
{
    char *p0,*p1;
    short *p2;
    float *p4;
    int *p5;
    double *p8;

    if (VTyp[j] == 1)                   /* string variables */
        return(0.0);
    else if (VTyp[j] == 2) {            /* numerical constants */
        p8 = (double *)VDPtr[j];
        return((double)p8[0]);
    }
    if (TSelFlg)                        /* get index from TSelect */
        i = TSelect[i];     
    else if (REPSelFlg)                 /* get index from REPSelect */
        i = REPSelect[i];

    switch (VSLen[j]) {
        
        case 0: p0 = VDPtr[j] + i / 8;
  
                if (*p0 & BMsk[i % 8])
                    return(1.0);
                else
                    return(0.0);
           
        case 1:  p1 = (char *)VDPtr[j];
                 return((double)p1[i]);

        case 2:  p2 = (short *)VDPtr[j];
                 return((double)p2[i]);

        case 4:  p4 = (float *)VDPtr[j];
                 return((double)p4[i]);

        case 5:  p5 = (int *)VDPtr[j];
                 return((double)p5[i]);

        case 8:  p8 = (double *)VDPtr[j];
                 return((double)p8[i]);

        default: break;
    }
    return(0.0);
}

/*--------------------------------------------------------------------------*/
/*  put_data(x,j,i)     Put value x in data matrix for variable j in case i */
/*                                                                          */
/*                      If TSelFlg == 1, apply tsel case selection.         */
          
void put_data(double x,int j,int i)
{
    char *p0,*p1;
    short *p2;
    float *p4;
    int *p5;
    double *p8;

    if (VTyp[j] == 1)               /* string variables */
        return;
    else if (VTyp[j] == 2) {        /* numerical constants */
        p8 = (double *)VDPtr[j];
        p8[0] = (double)x;
        return;
    }
    if (TSelFlg)                        /* get index to next selected case */
        i = TSelect[i];     
    else if (REPSelFlg)                 /* get index from REPSelect */
        i = REPSelect[i];

    switch (VSLen[j]) {

        case 0:  p0 = VDPtr[j] + i / 8;
   
                 if ((int)x)
                    *p0 |= BMsk[i % 8];
                 else
                    *p0 &= ~BMsk[i % 8];
                 break;

        case 1:  p1 = (char *)VDPtr[j];
                 p1[i] = (char)x;
                 break;
        case 2:  p2 = (short *)VDPtr[j];
                 p2[i] = (short)x;
                 break;
        case 4:  p4 = (float *)VDPtr[j];
                 p4[i] = (float)x;
                 break;
        case 5:  p5 = (int *)VDPtr[j];
                 p5[i] = (int)x;
                 break;
        case 8:  p8 = (double *)VDPtr[j];
                 p8[i] = (double)x;
                 break;
        default: break;
    }
}

/*--------------------------------------------------------------------------*/
/*  put_str(buf,blen,j,i,bflag)                                             */
/*                                                                          */
/*  Put string for variable j into data matrix for case i.                  */
/*  If bflag != 0 use blanks.                                               */
/*  If TSelFlg == 1, apply tsel case selection.                             */
          
void put_str(char *buf,int blen,int j,int i,int bflag)
{
    register int k;
    register char *p,*q,*m;
    int len;

    if (TSelFlg)                        /* get index to next selected case */
        i = TSelect[i];     
    else if (REPSelFlg)                 /* get index from REPSelect */
        i = REPSelect[i];

    len = -VSLen[j];
    p = buf;
    q = VDPtr[j] + i * len;
    if (bflag)
        m = p; 
    else
        m = buf + blen;

    for (k = 0; k < len; ++k) {
        if (p < m)
            *q++ = *p++; 
        else
            *q++ = ' ';
    }   
}

/*--------------------------------------------------------------------------*/
/*  clear_str(buf,blen,j)                                                   */
/*  Clear string for variable j in buf.                                     */
          
void clear_str(char *buf,int blen,int j)
{
    register int k;
    register char *p,*m;
    int n,len;

    n = VStrN[j];               /* first column of string */
    len = -VSLen[j];

    p = buf + n;
    m = buf + blen;

    for (k = 0; k < len; ++k) {
        if (p < m)
            *p++ = ' '; 
        else
            break;           
    }   
}

/*--------------------------------------------------------------------------*/
/*  get_str(buf,j,i)                                                        */
/*                                                                          */
/*  Get string from variable j, case i, into buf.                           */
/*  If TSelFlg == 1, apply tsel case selection.                             */
          
void get_str(char *buf,int j,int i)
{
    register int k;
    register char *p,*q;
    int len;

    if (TSelFlg)                        /* get index to next selected case */
        i = TSelect[i];     
    else if (REPSelFlg)                 /* get index from REPSelect */
        i = REPSelect[i];

    len = -VSLen[j];
    p = VDPtr[j] + i * len;
    q = buf;

    for (k = 0; k < len; ++k)  
        *q++ = *p++; 
    *q = '\0';
}

/* ------------------------------------------------------------------------ */
/*  dscan.  Scan ASCII string and convert to double.                        */
/*          On return, mval has one of the following values:                */
/*          0 no missing value                                              */
/*          1 missing value: blank                                          */
/*          2 missing value: star                                           */
/*          3 missing value: point                                          */
/*          4 missing value: general                                        */
/*                                                                          */
/*          If len > 0 only this number of characters is used.              */

#define MaxMDig     9       /* Max number of digits in the mantissa,        */
                            /* otherwise we use sscanf().                   */

double dscan(char *p,int len,int *mval)
{
    long int itmp = 0;
    register int i = 0;
    register int k = 0;
    register int n = 0;
    register int neg = 0;
    register int pt = 0;
    register char c,*q;
    double tmp = 0.0;

    *mval = 0;
    q = p;

    if (*q == '-') {
        i = neg = 1;
        q++;
    }
    else if (*q == '+') {
        i = 1;
        q++;
    } 

    while (*q) {

        if (sep_char(*q))
            break;

        if (*q == '*') {
            *mval = 2;
            return(MStarVal);
        }
        else if (*q >= '0' && *q <= '9') {
            itmp *= 10L;
            itmp += (long)(*q - '0');
            k++;
        }
        else if (*q == '.') {
            if (pt) {
                *mval = 4;
                return(MGenVal);
            }
            pt = 1;
            n = k;
        }
        else if (*q != '\n') {

            if (*q == 'E' || *q == 'e') {
                if (len) {
                    c = *(p + len);
                    *(p + len) = '\0';
                }
                i = sscanf(p,"%lg",&tmp);           
                if (len)
                    *(p + len) = c;
                if (i == 1)
                    return(tmp);
            }
            *mval = 4;
            return(MGenVal);
        }
        q++;
        if (len && ++i >= len)
            break;
    }
    if (!k) {
        if (pt) {
            *mval = 3;
            return(MPntVal);
        }
        *mval = 1;
        return(MBlnkVal);
    }
    else if (k > MaxMDig) {

        if (len) {
            c = *(p + len);
            *(p + len) = '\0';
        }
        if (sscanf(p,"%lg",&tmp) != 1) {
            *mval = 4;
            tmp = MGenVal;
        }
        if (len)
            *(p + len) = c;
        return(tmp);
    }
    if (neg)
        tmp = (double)(-itmp);
    else
        tmp = (double)itmp;

    if (pt) {
        k -= n;
        while (k--)
            tmp /= 10.0;
    }
    return(tmp);
}

/*--------------------------------------------------------------------------*/
/*  sep_char(c)     return 1 if c is a separation character, otherwise 0    */

int sep_char(char c)
{
    if (ISEPC != '\0') {
        if (c == ISEPC)     
            return(1);
        else
            return(0);
    }
    else if (c == ' ' || c == ',' || c == ';' || c == '\t' || c == LF || c == CR)
        return(1);
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  skip_sep(p)     skip separation characters in a data file record.       */
/*                  return pointer to the next character.                   */

char *skip_sep(char *p)
{
    if (ISEPC != '\0') {
        if (sep_char(*p))
            p++;
    }
    else {
        while (sep_char(*p))
            p++;
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  skip_dval(p)    skip an entry in a data file's record.                  */
/*                  return pointer to the next character.                   */

char *skip_dval(char *p)
{
    while (*p) {
        if (*p == '"')      
            while (*++p && *p != '"') ;

        if (sep_char(*p))
            break;
        p++;
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  gen_dm(idx,nvar,sflag,sdflag)                                           */  
/*                                                                          */
/*  Create nvar new variables, begining at index idx.                       */
/*                                                                          */
/*  A)  If DMDef = 0, that is, if a data matrix does not already exist.     */
/*      Create a maximum of NOCMaxA cases.                                  */
/*                                                                          */
/*      a)  If DFilN > 0 assume data files DFILFN[i] and create cases       */
/*          according to data file records.                                 */  
/*                                                                          */
/*          Assume that each logical record consists of NMRec physical      */
/*          records.                                                        */
/*                                                                          */
/*          If DRecLen == 0 assume variable record length, if DRecLen > 0   */
/*          assume fixed record length.                                     */
/*                                                                          */
/*          If FFMTCnt > 0 assume fixed format in FFMTC1[],FFMTC2[].        */
/*                                                                          */
/*      b)  Alternatively, if NVArc > 0, assume NVArc archive variables.    */  
/*                                                                          */
/*      If DBLKVar >= 0 assume block mode.                                  */
/*                                                                          */
/*  B)  If DMDef != 0 a data matrix already exists.                         */
/*                                                                          */
/*      a) If MatchNV = 0 assume trivial matching.                          */
/*                                                                          */
/*      b) If MatchNV > 0 match according to MatchVA[],MatchVB[]. Assume    */
/*         that data are already sorted with respect to MatchVB. Pointer    */
/*         is available in VSORTPtr.                                        */
/*                                                                          */
/*  C)  If GDFlg != 0 write directly to output file GDFd without saving     */
/*      variables in internal data matrix. In this case there are only      */
/*      variables of type 1 -- 3 and a data matrix does not already exist.  */
/*      Or also type 4 variables if in block mode.                          */
/*                                                                          */
/*  If sdflag != 0 assume spatial data generation, i.e., this function      */
/*  is called by sdnvar().                                                  */
/*                                                                          */
/*  Return 0 if success, other values indicate error.                       */

int gen_dm(int idx,int nvar,int sflag,int sdflag)
{
    register int j,k,l,ll;
    int rm,rn,r,err,jm,pcnt,rsba,drec,dflag,nb,bcnt,bi,nm,bmin,bmax;
    int nr,wflag,first,ns,nss,nst,fptr;
    double x,id,lastid;
           
    if (MatchNV > 0) {              
        if (alloc_acx(MatchNV))     /* use AcX for match values */
            return(-1);
        if (alloc_acn(NOCDM))       /* use AcN to indicate filled rows */
            return(-1);
    }
    first = 1;
    rsba = wflag = bmax = err = dflag = 0;
    bmin = INTMAX;
    if (DFILN > 0 || NVArc > 0)
        dflag = 1;

    drec = 0;   /* count data file records */
    MBlnkN = MStarN = MGenN = MPntN = 0;    /* counter for miss values */

    if (NVArc > 0) {            /* init reading from archive */
        dbf_init(1);
        if (SILENTFlg < 2)
            printfe("Reading archive file: %s\n",ZAFNam[AVDFN]);

        j = VIFirst;            /* clear AVVAL */
        while (j >= 0) {
            AVVAL[j] = 0.0;
            j = VNxt[j];
        }
    }
    else if (DFILN > 0) {       /* use data files */

        if (read_dfi(1,sflag)) {    /* init read_df function */
            err = -1;
            goto GDMFin;
        }
    }      

    lastid = DBLMAX;
    nb = bi = bcnt = 0;

    pcnt = ns = k = 0;  /* count data matrix rows */

    while (1) {
                    
        if (GDFlg && k >= NOCMaxA) {
            printf1("Error: exceeded maximum block size.\n");
            err = -1;
            break;
        }

        if (NVArc > 0) {    /* ## get archive variables into AVVAL[] */
                            /* only numerical variables */
            r = get_avar();
            if (r)  
                prn_message(++pcnt,0,0);
        }
        else if (DFILN > 0) {       /* get next record from data file(s) */
            r = read_df(sflag,idx); /* values are in VCJVal[]. */
            if (r)  
                prn_message(++pcnt,0,0);
        }
        else
            r = 1;

        if (r != 1) {               /* r = 0 if EOF */
            err = r;
            break;
        }       
        if (dflag) {
            drec++;
            if (VTNOC > 0 && drec > VTNOC)
                break;
        }
        else
            drec = k + 1;

        /* Check isel. Note that v_eval1() is called with vflg = 1,
           that is: values are assumed to be in VCJVal[] and AVVAL[]. */

        if (ISESCnt) {

            r = v_eval1(k,ISESCnt,ISESTyp,ISESVal,ESIdx,&x,1,0,0,0,0);
             
            if (r) {
                printf1("Can't evaluate isel expression for record %d in case %d.\n",
                                                          drec,k + 1);
                prn_emsg2(r);
                err = -1;
                break;
            }
            if (fabs(x) < EPSI2)  
                goto DMCONT;
        }
        if (MatchNV > 0) {

            for (j = 0; j < MatchNV; ++j) {
                jm = MatchVA[j];
                r = v_eval1(k,VESCnt[jm],VESTyp[jm],VESVal[jm],ESIdx,&x,1,0,0,0,0);
                if (r) {
                    printf1("Can't evaluate variable %s for sorted case %d (record %d).\n",
                                                    VName[jm],k + 1,drec);
                    prn_emsg2(r);
                    err = -1;
                    break;
                }
                AcX[j] = x;
            }
            if (err)
                break;

            rm = fnd_match(AcX,&rn);  /* check match */
            if (rm < 0)
                goto DMCONT;

            /* matching is for rn rows beginning at row rm */

        }

        /* put string variables into data matrix */

        j = idx;
        while (j >= 0) {
            if (VTyp[j] == 1 && VTypA[j] == 0) {    
                if (MatchNV == 0)
                    put_str(RSBuf + VStrN[j],RBufALen - VStrN[j],j,k,0); 
                else {
                    for (l = 0; l < rn; ++l) {
                        ll = VSORTPtr[rm + l];
                        if (AcN[ll] == 0)  
                            put_str(RSBuf + VStrN[j],RBufALen - VStrN[j],j,ll,0); 
                        else
                            wflag = 1;
                    }
                }
            }
            j = VNxt[j];
        }

        /* ## put archive variables into data matrix */

        if (NVArc > 0) {
            for (j = 0; j < NVArc; ++j) {
                if (MatchNV == 0) {
                    if (VTyp[AVIdx[j]] == 1) {  /* string var */
                        get_astr(SVBuf,j);
                        put_str(SVBuf,SVBufLen,AVIdx[j],k,0);
                    }
                    else
                        put_data(AVVAL[AVIdx[j]],AVIdx[j],k);
                }
                else {
                    for (l = 0; l < rn; ++l) {
                        ll = VSORTPtr[rm + l];
                        if (AcN[ll] == 0) {
                            if (VTyp[AVIdx[j]] == 1) {  /* string var */
                                get_astr(SVBuf,j);
                                put_str(SVBuf,SVBufLen,AVIdx[j],ll,0);
                            }
                            else
                                put_data(AVVAL[AVIdx[j]],AVIdx[j],ll);
                        }
                        else
                            wflag = 1;
                    }
                }
            }
        }

        /* evaluate remaining variables and put values into data matrix */
        /* skip variables with type >= 4 */

        j = idx;
        while (j >= 0) {

            if (VTyp[j] == 3 || (VTyp[j] == 2 && first == 1)) {

                /* archive variables already in data matrix */
                       
                if (MatchNV == 0) {

                    r = v_eval1(k,VESCnt[j],VESTyp[j],VESVal[j],ESIdx,&x,0,0,0,0,0);
                    if (r) {
                        printf1("Can't evaluate variable %s for case %d.\n",
                                                         VName[j],k + 1);
                        if (VCJMax > 0 || NVArc > 0)  
                            printf1("Evaluation based on data file record %d.\n",drec);
                        prn_emsg2(r);
                        err = -1;
                        break;          
                    }
                    put_data(x,j,k);      
                }
                else {

                    for (l = 0; l < rn; ++l) {
                        ll = VSORTPtr[rm + l];
                        if (AcN[ll] == 0) {
                            r = v_eval1(ll,VESCnt[j],VESTyp[j],VESVal[j],ESIdx,&x,0,0,0,0,0);
                            if (r) {
                                printf1("Can't evaluate variable %s for sorted case %d.\n",
                                                            VName[j],ll + 1);
                                if (VCJMax > 0 || NVArc > 0)  
                                    printf1("Evaluation based on data file record %d.\n",drec);
                                prn_emsg2(r);
                                err = -1;
                                break;          
                            }
                            put_data(x,j,ll);      
                        }
                        else
                            wflag = 1;
                    }
                    if (err)
                        break;
                }
            }
            j = VNxt[j];
        }
        first = 0;
        if (err)
            break;

        /* if block mode check ID variable and evaluate type 4 variables */
        /* cannot be used if matching data */

        if (DBLKVar >= 0) {

            id = get_data(DBLKVar,k);
            if (fabs(id - lastid) > EPSI1) {

                if (k > bi) {                   /* new block */
                    err = put_t4var(idx,bi,k,nb);  /* evaluate type 4 variables */
                    if (err)
                        break;

                    if (BSESCnt) {              /* check bsel */
                        nr = tst_bsel(bi,k,idx);
                        if (nr <= 0) {
                            if (nr < 0) {
                                err = nr;
                                break;
                            }   
                            cpy_dmrow(k,bi,idx);
                            k = bi;             /* skip whole block */
                            lastid = id;
                            bcnt = 1;
                            goto DMCONTA;
                        }
                    }
                    else
                        nr = bcnt;

                    nb++;
                    if (bmin > nr)
                        bmin = nr;
                    if (bmax < nr)
                        bmax = nr;


                    if (GDFlg) {            /* write to output file */
                        prn_gdf(nr);
                        cpy_dmrow(k,bi,idx);
                        if (VTNOC > 0 && GDNRec >= VTNOC)
                            break;
                    }
                    else {
                        cpy_dmrow(k,bi + nr,idx);
                        bi += nr;
                    }
                }
                lastid = id;
                bcnt = 1;
                k = bi;
            }
            else 
                bcnt++;
        }

        if (VSESCnt && DBLKVar < 0) {   /* check vsel if not in block mode */

            if (MatchNV == 0) {

                r = tst_vsel(k,drec);
                if (r < 0) {
                    err = r;
                    break;
                }
                else if (r == 0)
                    goto DMCONT;
            }
            else {
                for (l = 0; l < rn; ++l) {
                    ll = VSORTPtr[rm + l];
                    if (AcN[ll] == 0) {
                        r = tst_vsel(ll,drec);
                        if (r < 0) {
                            err = -1;
                            break;
                        }
                        else if (r > 0) {
                            AcN[ll] = 1;
                            ns++;
                        }
                    }
                }
                if (err)
                    break;
            }
        }
        else if (MatchNV > 0) {     /* indicate values */

            for (l = 0; l < rn; ++l) {
                ll = VSORTPtr[rm + l];
                if (AcN[ll] == 0) {
                    AcN[ll] = 1;
                    ns++;
                }
            }
        }  

        if (DBLKVar < 0) {      /* if not in block mode */

            /* Check break. Note that v_eval1() is called with vflag=0, that is,
               values are assumed to be already in data matrix. Cannot be used
               if matching data. */

            if (BRKESCnt) {
  
                r = v_eval1(k,BRKESCnt,BRKESTyp,BRKESVal,ESIdx,&x,0,0,0,0,0);
                if (r) {
                    printf1("Can't evaluate break expression for record %d in case %d.\n",drec,k + 1);
                    prn_emsg2(r);
                    err = -1;
                    break;
                }
                if (fabs(x) >= EPSI2)  
                    break;         
            }
        }

DMCONTA:
        k++;                    /* count data matrix rows */

        if (GDFlg && DBLKVar < 0) {
            prn_gdf(k);
            k = 0;
            if (VTNOC > 0 && GDNRec >= VTNOC)
                break;
            goto DMCONT;
        }
        if (sdflag) {        /* special for spatial data */

            nst = get_data(SDVarSDTyp,k - 1);    /* type of object */
            if (nst == 1)
                SDVarNP++;
            else if (nst == 2)
                SDVarNL++;
            else if (nst == 3)
                SDVarNPol++;
            else 
                SDVarNU++;

            nss = get_data(SDVarSDN,k - 1);    /* number of following records */

            if (nss < 1) {
                printf1("Error: SDN contains a value less than 1 (object %d).\n",k);
                err = -1;
                break;
            }
            if (nst == 1) {
                if (nss != 1) {
                    printf1("Error: type 1 object (%d) has more than one point.\n",k);
                    err = -1;
                    break;
                }
                SDVarNT += 1;
            }
            else if (nst == 2) {
                if (nss < 2) {
                    printf1("Error: type 2 object (%d) has less than two points.\n",k);
                    err = -1;
                    break;
                }
                SDVarNT += nss;
            }
            else if (nst == 3) {
                if (nss < 3) {
                    printf1("Error: type 3 object (%d) has less than three points.\n",k);
                    err = -1;
                    break;
                }
                SDVarNT += nss;
            }
            if ((fptr = skip_df(nss,1)) < 0) {       
                printf1("Error: cannot read %d data records for object %d.\n",nss,k);
                err = -1;
                break;
            }
            put_data((double)fptr,SDVarSDPtr,k - 1);      
            SDVarMax = imax(SDVarMax,nss);
        }
        if (MatchNV == 0) {
            if (k >= NOCMaxA || (NOCDM > 0 && k >= NOCDM))
                break;
        }
        else if (ns >= NOCDM)
            break;

DMCONT: 
        if (VTNOC > 0 && k >= VTNOC)
            break;
    }
    if (err)
        goto GDMFin;

    if (dflag)
        prn_message(pcnt,1,0);

    if (k > bi) {
        err = put_t4var(idx,bi,k,nb);  /* put type 4 variables into data matrix */
        if (err)
            goto GDMFin;             

        if (DBLKVar >= 0) {                 /* block mode */
            if (BSESCnt) {                  /* check bsel */
                nr = tst_bsel(bi,k,idx);
                if (nr <= 0) {
                    if (nr < 0) {
                        err = nr;
                        goto GDMFin;
                    }   
                    k = bi;       
                }
            }
            else
                nr = bcnt;

            if (nr > 0) {
                nb++;
                if (bmin > nr)
                    bmin = nr;
                if (bmax < nr)
                    bmax = nr;

                if (GDFlg) {            /* write to output file */
                    prn_gdf(nr);
                    k = 0;
                }   
                else
                    k = bi + nr;    
            }
        }
    }
    nm = 0;                     /* number of cases with no match */
    if (NOCDM == 0)
        NOC = NOCDM = k;

    else if (MatchNV > 0) {     /* insert missing values for non-match */
        nm = NOCDM - ns;
        if (nm > 0) {
            for (l = 0; l < NOCDM; ++l) {
                if (AcN[l] == 0) {
                    j = idx;
                    while (j >= 0) {
                        if (first || VTyp[j] != 2) {
                            if (VTyp[j] != 1)
                                put_data(MMatchVal,j,l);
                            else if (VTyp[j] == 1)
                                put_str(RSBuf + VStrN[j],RBufALen - VStrN[j],j,l,1); 
                        }
                        j = VNxt[j];
                    }
                }
            }
        }
    }
    else if (k < NOCDM) {     /* fill with miss values */

        nm = NOCDM - k;
        while (k < NOCDM) {
            j = idx;
            while (j >= 0) {
                if (VTyp[j] != 1)
                    put_data(MMatchVal,j,k);
                else
                    put_str(RSBuf + VStrN[j],RBufALen - VStrN[j],j,k,1); 
                j = VNxt[j];
            }
            k++;
        }
    }
    if (dflag) {
        printf1("Read records: %d ",drec);
        if (VCJMax > 0 && NMRec > 1)  
            printf1("(%d physical records).",NRec);
        printf1("\n");
    }

GDMFin:  
    if (MatchNV > 0) {
        alloc_acx(0);
        alloc_acn(0);
    }
    if (DFILN > 0)     /* free read buffer and reset values */
        read_dfi(0,0);    

    if (err == 0) {
        if (GDFlg == 0) {
            if (ISESCnt || VSESCnt || BRKESCnt || BSESCnt)  
                printf1("Selected for data matrix: %d records.\n",NOCDM);
        }
    }                       
    if (err || NOCDM <= 0) {
        if (err)
            err = -1;
          
        if (GDFlg == 0) {
            if (DMDef == 0)
                printf1("No data matrix created.\n");
            else
                printf1("Nothing added to current data matrix.\n");
        }
        else if (err == 0) {
            printf1("%d records written to output file: %s\n",GDNRec,GDFNPtr);
            if (DBLKVar >= 0) {
                printf1("Number of blocks: %d\n",nb);
                printf1("Min number of cases per block: %d\n",bmin);
                printf1("Max number of cases per block: %d\n",bmax);
            }
        }
        return(err);
    }
    if (DMDef == 0) {
        printf1("\nCreated a new data matrix.\n");
        printf1("Number of cases: %d\n",NOCDM);
        printf1("Number of variables: %d\n",nvar);
        DMDef = 1;
    }
    else {
        printf1("\nAdded %d variable(s) to existing data matrix.\n",nvar);
        printf1("Number of cases with no match: %d\n",nm);
        if (nm > 0)  
            printf1("Substituted missing value code: %g\n",MMatchVal);
    }
    if (wflag) {
        printf1("Warning: values in variable(s) %s",VName[MatchVA[0]]);
        for (j = 1; j < MatchNV; ++j)  
            printf1(",%s",VName[MatchVA[j]]);
        printf1(" are not unique.\n");
    }
    if (DBLKVar >= 0) {
        printf1("\nNumber of blocks: %d\n",nb);
        printf1("Min number of cases per block: %d\n",bmin);
        printf1("Max number of cases per block: %d\n",bmax);
    }

    /* set VTypA = 1 for new variables */

    j = idx;    
    while (j >= 0) {
        VTypA[j] = 1;
        j = VNxt[j];
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  prn_gdf(n)          Print n records from data matrix to output file     */
/*                      GDFd. Use variables in NVarIdx[].                   */
/*                      Use print formats in NVPFmtS[]                      */
/*                                                                          */
/*                      If GDKFlg != 0 keep only first record in block.     */  

void prn_gdf(int n)
{     
    register int i,j,k;   
    double tmp;

    for (i = 0; i < n; ++i) {
        if (i > 0 && GDKFlg && DBLKVar >= 0)
            break;
        for (j = 0; j < NVarNV; ++j) {
            k = NVarIdx[j];
            if (k >= 0) {
                if (VTyp[k] != 1) {
                    tmp = get_data(k,i);
                    fprintf(GDFd,NVPFmtS[j],tmp);
                }
                else {  
                    get_str(SVBuf,k,i);
                    fprintf(GDFd,"%s",SVBuf);
                    if (NVSEPC)
                        fprintf(GDFd,"%c",NVSEPC);
                }
            }
        }
        fprintf(GDFd,"\n");
        GDNRec++;
    }
}

/* ------------------------------------------------------------------------ */
/*  fnd_match(xm,n)                                                         */
/*                  Search for a record in the data matrix where xm[]       */
/*                  is equal to values of the MatchVB[] variables.          */
/*                  If no match return -1, otherwise return index j such    */
/*                  that VSORTPtr[j] points to first data matrix row with   */
/*                  a matching value. Also, return in n the number of       */
/*                  rows with a matching value.                             */

int fnd_match(double *xm,int *n)
{     
    register int i,j,l,r;   
    int fnd,jj,nn;
    double tmp;                  

    *n = 0;
    l = 0;
    r = NOCDM - 1;

    while (r >= l) {

        j = (l + r) / 2;

        fnd = 1;
        tmp = get_data(MatchVB[0],VSORTPtr[j]);
        if (fabs(tmp - xm[0]) > EPSI2)  
            fnd = 0;

        if (fnd) {

            while (j > 0) {
                if (fabs(tmp - get_data(MatchVB[0],VSORTPtr[j - 1])) > EPSI2) 
                    break;
                j--;
            }
            jj = -1;
            nn = 0;
            fnd = 0;
            while (j < NOCDM) {
                if (fabs(tmp - get_data(MatchVB[0],VSORTPtr[j])) > EPSI2)  
                    break;
  
                fnd = 1;
                for (i = 1; i < MatchNV; ++i) {
                    if (fabs(xm[i] - get_data(MatchVB[i],VSORTPtr[j])) > EPSI2) {
                        fnd = 0;
                        break;
                    }   
                }
                if (fnd) {
                    if (jj < 0)
                        jj = j;
                    nn++;
                }
                j++;
            }
            *n = nn;
            return(jj);
        }
        else if (xm[0] < tmp)  
            r = j - 1;
        else      
            l = j + 1;
    }
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  put_t4var(idx,r1,r2,bn)   Put type 4 variables into data matrix. Begin  */  
/*                            with index idx. Use data matrix rows r1 - r2. */
/*                            Return 0 if OK, -1 if error.                  */  
          
int put_t4var(int idx,int r1,int r2,int bn)
{
    int j,r;
    double tmp;

    j = idx;
    while (j >= 0) {

        if (VTyp[j] == 4) {
            r = v_eval2(j,VESCnt[j],VESTyp[j],VESVal[j],ESIdx,&tmp,r1,r2,0,bn);
            if (r) {
                printf1("\nCan't evaluate variable %s in block containing cases %d - %d.\n",VName[j],r1 + 1,r2);
                prn_emsg2(r);
                return(-1);
            }
        }
        j = VNxt[j];
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  cpy_dmrow(ir,ij,idx)    Copy data matrix row from ir to ij, use         */
/*                          variables beginning at idx.                     */

void cpy_dmrow(int ir,int ij,int idx)
{
    register int j;
    double tmp;

    if (ir == ij)
        return;

    j = idx;
    while (j >= 0) {
        if (VTyp[j] != 1) {
            tmp = get_data(j,ir);
            put_data(tmp,j,ij);      
        }
        else {
            get_str(SVBuf,j,ir);
            put_str(SVBuf,SVBufLen,j,ij,0);
        }
        j = VNxt[j];
    }
}

/*--------------------------------------------------------------------------*/
/*  tst_bsel(r1,r2,idx)     Test bsel expression for cases r1 ... r2.       */
/*                          Copy selected rows into contiguous block.       */
/*                          Return -1 if error, or number in block.         */

int tst_bsel(int r1,int r2,int idx)
{
    register int i,k;
    int r,n;
    double tmp;

    n = 0;
    i = r1;
    for (k = r1; k < r2; ++k) {
        r = v_eval1(k,BSESCnt,BSESTyp,BSESVal,ESIdx,&tmp,0,0,0,0,0);
        if (r) {
            printf1("Can't evaluate vsel expression in case %d.\n",k + 1);
            prn_emsg2(r);
            return(-1);
        }
        if (fabs(tmp) >= EPSI2) {   /* copy row k into row i */
            cpy_dmrow(k,i,idx);
            i++;
            n++;
        }
    }
    return(n);
}

/* -------------------------------------------------------------------------*/
/*  tst_vsel(i,rec)     Test vsel for case i, record rec.                   */
/*                      Return  -1  if error                                */  
/*                               1  if expression is true,                  */  
/*                               0  if not true.                            */

int tst_vsel(int i,int rec)
{
    int r;
    double tmp;

    r = v_eval1(i,VSESCnt,VSESTyp,VSESVal,ESIdx,&tmp,0,0,0,0,0);
    if (r) {
        printf1("Can't evaluate vsel expression for record %d in case %d.\n",rec,i + 1);
        prn_emsg2(r);
        return(-1);
    }
    if (fabs(tmp) < EPSI2)  
        return(0);                 
    return(1);
}

/* -------------------------------------------------------------------------*/
/*  tst_break(r1,r2)    Test break expression for cases r1 ... r2.          */
/*                      Return  -1  if error                                */  
/*                               1  if true for at least one record.        */  
/*                               0  if to be used                           */

int tst_break(int r1,int r2)
{
    int k,r;
    double tmp;

    for (k = r1; k < r2; ++k) {
        r = v_eval1(k,BRKESCnt,BRKESTyp,BRKESVal,ESIdx,&tmp,0,0,0,0,0);
        if (r) {
            printf1("Can't evaluate break expression in case %d.\n",k + 1);
            prn_emsg2(r);
            return(-1);
        }
        if (fabs(tmp) >= EPSI2)  
            return(1);   
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  read_dfi(opt,sflag)                                                     */
/*                                                                          */
/*  Init for reading data files.                                            */
/*  sflag != 0 if string variables.                                         */
/*  if opt != 0 prepare calls of function read_df().                        */
/*  if opt == 0 free memory and reset values.                               */
/*                                                                          */
/*  return 1 if OK, -1 if error.                                            */
          
int read_dfi(int opt,int sflag)
{
    if (opt) {
        if (DRecLen == 0)  
            RBufLen = RLMaxDef + 1;
        else {
            DRecLen1 = DRecLen * NMRec;
            RBufLen = DRecLen1 + 1;
        }
        if (!(RBuf = (char *)calloc(RBufLen + 1,sizeof(char)))) {
            p_err(-2,1);
            RBufLen = 0;
            return(-1);
        }
        if (sflag) {
            if (!(RSBuf = (char *)calloc(RBufLen + 1,sizeof(char)))) {
                p_err(-2,1);
                free(RBuf);
                RBufLen = 0;
                return(-1);
            }
            RSBufA = RBufLen + 1;
            memrq(RSBufA,sizeof(char));          
        }
        memrq(RBufLen + 1,sizeof(char));          
        RBufA = RBufLen + 1;
        NRec1 = NRec = 0;
        DFILNI = 0;     /* first data file */
        REOF = 1;       /* not open */
    }
    else {
        if (REOF == 0)
            fclose(RFd);

        if (RBufA > 0) {
            free(RBuf);
            memrq(-RBufA,sizeof(char));          
            RBufA = 0;
        }
        if (RSBufA > 0) {
            free(RSBuf);
            memrq(-RSBufA,sizeof(char));          
            RSBufA = 0;
        }
        RBufLen = 0;
    }
    REOF = 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  read_df()       read the next record from input data file(s) and save   */
/*                  results in VCJVal[]. It is assumed that the read buffer */
/*                  RBuf is already allocated.                              */
/*                                                                          */
/*  Note: files are only read up to lines that begin with ^Z (hex 1a).      */
/*  Return 1 if OK, 0 if EOF, -1 if error.                                  */
          
int read_df(int sflag,int idx)
{
    int i,l,rc,mval,err,cnt;
    register char *p,*q;

    RBufALen = err = 0;

    while (1) {
        if (REOF) {
            if (DFILNI >= DFILN)    /* no more data files */
                break;       

            if (!(RFd = fopen(DFILFN[DFILNI],OPEN_RD))) {
                printf1("Can't open data file: %s\n",DFILFN[DFILNI]);
                return(-1);
            }
            NRec1 = REOF = 0;
            printf1("Reading data file: %s\n",DFILFN[DFILNI]);
            DFILNI++;
        }            
        while (1) {
            if (DRecLen == 0) {

                if (fgets(RBuf,RBufLen,RFd) == NULL || *RBuf == 0x1a)
                    break;
                rc = strlen(RBuf);

                if (NMRec > 1) {

                    p = RBuf + rc;

                    for (i = 1; i < NMRec; ++i) {
                        cnt = RBufLen - rc;
                        if (cnt < 10) {
                            printf1("Error: exceeded max read buffer length (%d bytes).\n",RLMaxDef);
                            err = -1;
                            break;
                        }
                        if (fgets(p,cnt,RFd) == NULL) {
                            printf1("Warning: input file does not contain a multiple of %d records.\n",NMRec);
                            printf1("Stopped with reading %d physical records from %s.\n",NRec1,DFILFN[DFILNI-1]);
                            break;
                        }
                        cnt = strlen(p);
                        rc += cnt;
                        p += cnt;
                    }
                    if (i < NMRec)
                        break;
                }
                NRec += NMRec;
                NRec1 += NMRec;
            }
            else {
                rc = fread(RBuf,sizeof(char),DRecLen1,RFd);
                if (rc != DRecLen1) {
                    if (rc > 0) {
                        printf1("Warning: stopped with reading a logical record of less than %d bytes.\n",DRecLen1);
                        printf1("Last physical record(s) (%d bytes) will be ignored.\n",rc);
                    }
                    break;
                }
                NRec += NMRec;
                NRec1 += NMRec;
            }
            if (check_drec(RBuf)) {      /* check for data records */

                p = RBuf + strlen(RBuf);
                while (--p >= RBuf && (*p == LF || *p == CR))
              		    rc--;
                *(RBuf + rc) = '\0';

                p = RBuf;
                RBufALen = rc;                  /* actual buffer length */
                for (i = 0; i < rc; ++i) {
                    if (*p == LF || *p == CR) {
                        if (ISEPC != '\0')
                            *p = ISEPC;
                        else
                            *p = ' ';
                    }
                    p++;
                }

                /* if string variables, copy RBuf into RSBuf and clear
                   strings in RBuf */

                if (sflag) {
                    p = RBuf;
                    q = RSBuf;
                    for (i = 0; i < rc; ++i)
                        *q++ = *p++;

                    if (FFMTCnt == 0) {
                        i = idx;
                        while (i >= 0) {
                            if (VTyp[i] == 1 && VTypA[i] == 0)
                                clear_str(RBuf,RBufALen,i); 
                            i = VNxt[i];
                        }
                    }
                }

                /* read required values into VCJVal[] */
    
                if (FFMTCnt > 0) {      /* fixed format */

                    for (i = 1; i <= VCJMax; ++i) {
                        if (VCJFlg[i]) {
                            if (FFMTC1[i] > rc) {
                                VCJVal[i] = MBlnkVal;
                                MBlnkN++;  
                                continue;
                            }

                            if (FFMTC2[i] > rc) {
                                err = -1;
                                break;
                            }
                            p = RBuf + FFMTC1[i] - 1;
                            l = FFMTC2[i] - FFMTC1[i] + 1;
                            while (l > 0) {
                                if (*p != ' ')
                                    break;
                                p++;
                                l--;
                            }
                            if (l <= 0) {
                                VCJVal[i] = MBlnkVal;
                                MBlnkN++;  
                            }
                            else {
                                VCJVal[i] = dscan(p,l,&mval);
                                if (check_mval(mval)) {
                                    err = -1;
                                    break;
                                }
                            }
                        }
                    }
                }
                else {          /* free format */

                    p = RBuf;
                    for (i = 1; i <= VCJMax; ++i) { 

                        mval = 0;
                        q = skip_sep(p);    /* skip separator */

                        if (VCJFlg[i]) {
                            if (!*q || *q == LF || *q == CR) {

                                if (ISEPC != '\0' && *q == ISEPC) {
                                    VCJVal[i] = MBlnkVal;
                                    MBlnkN++;
                                }
                                else {
                                    err = -1;
                                    break;
                                }
                            }
                            else {
                                VCJVal[i] = dscan(q,0,&mval);
                                if (check_mval(mval)) {

                                    err = -1;
                                    break;
                                }
                            }
                        }
                        if (i == VCJMax)
                            break;

                        p = skip_dval(q);
                    }
                }
                if (err) {
                    printf1("Can't read c%d in data file %s, record %d.\n",i,DFILFN[DFILNI-1],NRec1);
                    break;
                }
                return(1);      /* return this record, values in VCJVal */
            }
        }
        fclose(RFd);
        REOF = 1;
        if (err)
            break;       
    }
    return(err);
}

/*--###---------------------------------------------------------------------*/
/*  skip_df(n,opt)  Skip n record in current data file. Return file         */
/*                  pointer to first of these records, or -1 if error.      */
/*                  If opt != 0 read the data records and update the        */
/*                  bounding box parameters.                                */

int skip_df(int n,int opt)  
{
    int i,fptr;
    double x,y; 

    fptr = (int)ftell(RFd);
    if (fptr < 0)
        return(-1);

    for (i = 0; i < n; ++i) {
        if (fgets(RBuf,RBufLen,RFd) == NULL || *RBuf == 0x1a)
            return(-1);

        if (opt) {
            if (get_nsdxy(RBuf,&x,&y))  
                return(-1);

            SDVarXMin = dmin(x,SDVarXMin);
            SDVarXMax = dmax(x,SDVarXMax);
            SDVarYMin = dmin(y,SDVarYMin);
            SDVarYMax = dmax(y,SDVarYMax);
        }
    }       
    return(fptr);
}

/* ------------------------------------------------------------------------ */
/*  get_nsdxy(buf,x,y)                                                      */
/*                                                                          */ 
/*  It is assumed that buf contains at least two floating point entries     */
/*  that are separated by at least one blank character. This function       */
/*  returns the two values in x and y, respectively.                        */
/*                                                                          */
/*  Return 0 if successful, -1 if error.                                    */

int get_nsdxy(char *buf,double *x,double *y)
{
    char *p;

    p = skip_b(buf);
    if (sscanf(p,"%lg",x) != 1)
        return(-1);

    p = skip_dbl(p);
    p = skip_b(p);
    if (sscanf(p,"%lg",y) != 1)
        return(-1);
    return(0);
}

/*--###---------------------------------------------------------------------*/
/*  check_mval(mval)        check missing values.                           */
          
int check_mval(int mval)
{
    if (mval == 0)
        return(0);
    if (mval == 1)
        MBlnkN++;
    else if (mval == 2)
        MStarN++;
    else if (mval == 3)
        MPntN++;
    else if (MGenFlg)
        MGenN++;
    else
        return(-1);
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  sdnvar()    Read a data file to create a spatial data structure.        */
/*                                                                          */
/*  This command is similar to the nvar command. It is required that a      */
/*  data matrix does not already exist. Also the first four variables       */
/*  must be defined as follows and in this order:                           */
/*                                                                          */
/*  SDID     = ...,     ID of spatial object                                */
/*  SDTyp    = ...,     Type of spatial object                              */
/*                      1 = point                                           */  
/*                      2 = line                                            */ 
/*                      3 = polygon                                         */
/*  SDN = ...,          number of coordinates of the object (must be 1 if   */
/*                      SDTyp = 1).                                         */
/*  SDPtr<5> = rd,      used as a file pointer for the objects              */
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdnvar(void)
{
    register int i;
    int err,sflag;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Reading a spatial data file. Current memory: %d bytes.\n",MemReq);
    if (DMDef) {
        printf1("Error: a data matrix already exists.\n");
        return(-1);
    }
    sdnvar_close();         /* close previous file (if any) */

    SDVarNT = SDVarNP = SDVarNL = SDVarNPol = SDVarNU = 0;
    SDVarXMin = SDVarYMin =  (DBLMAX / 10.0);
    SDVarXMax = SDVarYMax = -(DBLMAX / 10.0);

    DFILN = MatchNV = VFmtFlg = NVNum = 0;

    DRecLen = 0;            /* variable record length */
    NMRec = 1;              /* logical = physical record */
    DBLKVar = -1;           /* no block mode */
    AVDFN = -1;             /* no archive file */
    ARCDICFlg = 0;          /* set by arcdic option */
    RDMNInit = 0;           /* init random number generator */
    ISEPC = '\0';           /* separation character */
    VTNOC = 0;
    BSIZE = 0;

    EDVALSn = EDVALOrg = EDVALDes = 0;
    EDVALTs = EDVALTf = 0.0;
         
    NVIdx = get_nidx();     /* index to first new variable */

    /* set default missing value codes */

    MStarVal  = -1.0; 
    MPntVal   = -1.0; 
    MBlnkVal  = -1.0; 
    MGenVal   = -1.0; 
    MMatchVal = -3.0; 
    GDNRec = GDKFlg = GDFlg = MGenFlg = MStarN = MPntN = MBlnkN = MGenN = 0; 

    ISelPtr  = NULL;        /* set by isel option */
    VSelPtr  = NULL;        /* set by vsel option */
    BSelPtr  = NULL;        /* set by bsel option */
    GDFNPtr  = NULL;        /* file name by df option */
    KeepPtr  = NULL;            
    DropPtr  = NULL;            

    NVArc = NVArc1 = 0;     /* number of archive variables, counted by save_var() */

    if (check_nvsd())       /* check command etc */
        goto NVSDFin;

    if (NVNum == 0) {
        printf1("No variables defined.\n");
        goto NVSDFin;
    }
    if (DFILN == 0) {
        printf1("No data file defined.\n");
        goto NVSDFin;
    }
    else if (DFILN > 1) {
        printf1("Error: can only use a single data file.\n");
        goto NVSDFin;
    }
    newline();      

    if (VLabelLen > 0 && VLabelLen < 8)      
        VLabelLen = 8;

    sflag = 0;
    i = NVIdx;                  /* count string variables */
    while (i >= 0) {
        if (VTyp[i] == 1 && VTypA[i] == 0)
            sflag++;
        i = VNxt[i];
    }
    if (VFmtFlg)
        make_vfmt(NVIdx);   /* make new print formats if requested with fmt */

    prn_var(NVIdx);         /* print list of new variables */

    /* check for standard variables */

    if (check_sd(1))
        goto NVSDFin;

    NOCMaxA = VTNOC;
    if (NOCMaxA <= 0)  
        NOCMaxA = NOCDef;

    printf1("\nCreating a new data matrix.\n");
    printf1("Maximum number of cases: %d\n",NOCMaxA);

    if (check_vcj(NVIdx,1))             /* check cj references */
        goto NVSDFin;

    printf1("\nUsing data file(s): %s",DFILFN[0]);
    for (i = 1; i < DFILN; ++i)
        printf1(",%s",DFILFN[i]);
    printf1("\n");
    if (DRecLen > 0)
        printf1("Fixed record length: %d\n",DRecLen);

    if (FFMTCnt > 0) {
        if (check_ffmt(1))
            goto NVSDFin;
    }
    else {
        printf1("Free format. Separation character(s): ");
        if (ISEPC != '\0')
            printf1("%04x [hex]\n",ISEPC);
        else 
            printf1("default.\n");
    }
    if (VTNOC > 0)  
        printf1("Reading maximal %d records.\n",VTNOC);

    /* allocate memory for new variables */

    if (alloc_vdat(NVIdx,1)) {
        printf1("Insufficient memory for new variables.\n");
        goto NVSDFin;
    }
    err = gen_dm(NVIdx,NVNum,sflag,1);      /* create data */
    if (err)  
        goto NVSDFin;

    prn_mval();                           /* info about missing values */

    printf1("\nNumber of points: %d\n",SDVarNP);
    printf1("Number of lines: %d\n",SDVarNL);
    printf1("Number of polygons: %d\n",SDVarNPol);
    printf1("Number of unknown objects: %d\n",SDVarNU);

    printf1("\nMaximal number of points in spatial objects: %d\n",SDVarMax);
    printf1("Total number of points in spatial objects: %d\n",SDVarNT);
    if (SDVarMax > 0) {
        if (!(SDVarFd = fopen(DFILFN[0],OPEN_RD))) {
            printf1("Error: cannot re-open data file: %s\n",DFILFN[0]);
            err = -1;
            goto NVSDFin;
        }
        SDVarFDef = SDVarDef = 1;  

        printf1("Note: %s remains opened for further access.\n\n",DFILFN[0]);

        printf1("X values. Minimum: %14.4lf  Maximum: %14.4lf\n",SDVarXMin,SDVarXMax);
        printf1("Y values. Minimum: %14.4lf  Maximum: %14.4lf\n",SDVarYMin,SDVarYMax);

        if (sdnvar_alloc(1,SDVarMax + 20)) {
            printf1("\nError: Insufficient memory for spatial objects.\n");
            err = -1;
        }
    }
    else
        err = -1;

NVSDFin:
    alloc_vl(0);

    check_vcj(NVIdx,0);
    check_ffmt(0);

    if (err) {
        if (NVNum > 0)
            clear_avar(NVIdx);
        printf1("No new variables created. ");
        sdnvar_close();
    }
    else {
        WIVar = -1;
        WSum = (double)NOC;
        WSumS = 0.0;
        WNorm = 1.0;
        WNormFlag = 0;
        printf1("\nEnd of creating new variables. ");
    }
    prn_mem();
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  check_nvsd  Check sdnvar command in CmdBuf.                             */
/*              Return 0 if OK, -1 if error.                                */

int check_nvsd(void)
{
    int err,n,l,m1,m2;
    register char c,*p,*q;
    char sc;
    double x;

    NVSEPC = ' ';             
    NVL0 = VTNOC = BSIZE = 0;
    err = -1;
    p = CmdBuf + 6;
    while (*++p) {
        if (sscanf(p,"noc=%d",&n) == 1 && n > 0) {
            VTNOC = n;
            p = skip_int(p + 4);
        }
        else if (sscanf(p,"dreclen=%d",&n) == 1 && n > 0) {
            DRecLen = n;
            p = skip_int(p + 8);
        }
        else if (sscanf(p,"fmt=%d.%d",&VFmt1,&VFmt2) == 2) {
            VFmtFlg = 1;
            p = skip_int(p + 4);
            p = skip_int(p + 1);
        }
        else if (sscanf(p,"mstar=%lf",&x) == 1) {
            MStarVal = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"mpnt=%lf",&x) == 1) {
            MPntVal = x;
            p = skip_dbl(p + 5);
        }
        else if (sscanf(p,"mblnk=%lf",&x) == 1) {
            MBlnkVal = x;
            p = skip_dbl(p + 6);
        }
        else if (sscanf(p,"mgen=%lf",&x) == 1) {
            MGenVal = x;
            MGenFlg = 1;
            p = skip_dbl(p + 5);
        }
        else if (!strncmp(p,"dfile=",6)) {
            if (DFILN >= DFILMax) {
                printf1("Exceeded max number of data files.\n");
                goto CNVSDFin;
            }
            DFILFN[DFILN++] = p + 6;
            p = skip_nc(p);
        }
        else if (!strncmp(p,"break=",6)) {
            BRKPtr = p + 6;
            p = skip_expr(p + 6);
        }
        else if (sscanf(p,"isc=%c",&sc) == 1) {
            if (sc == 't')
                sc = '\t';
            ISEPC = sc;
            p += 5;                
        }
        else if (!strncmp(p,"ffmt=c",6)) {
            FFMTCnt = 0;       
            q = p + 5;
            while (1) {
                if ((sscanf(q,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d,%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)) {
                    q = skip_int(q + 1);
                    q = skip_blev(q);
                }
                else {
                    prn_nve(p);
                    goto CNVSDFin;
                }
                FFMTCnt++;
                if (*q != ',' || *(q + 1) != 'c')
                    break;
                q++;       
            }
            FFMTPtr = p + 5;
            p = q;
        }
        else if ((l = get_vnlen(p)) > 0) {       /* variable */
            q = skip_nc(p + l);
            c = *q;
            *q = '\0';
            if (save_var(p,0))
                goto CNVSDFin;
            NVNum++;
            *q = c;
            p = q; 
        }
        if (*p != ',' && *p != ')') {
            prn_nve(p);
            goto CNVSDFin;
        }
        *p = '\0';
    }
    err = 0;

CNVSDFin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  check_sd(opt)   Check whether the standard variables are correctly      */
/*                  defined. Return 0 if OK, -1 if error.                   */

int check_sd(int opt)
{
    if (VTyp[SDVarSDID]   < 2 || VTyp[SDVarSDID]  > 3 ||
        VTyp[SDVarSDPtr] != 3 ||
        VTyp[SDVarSDTyp]  < 2 || VTyp[SDVarSDTyp] > 3 ||
        VTyp[SDVarSDN]    < 2 || VTyp[SDVarSDN]   > 3 ||  
        strcmp(VName[SDVarSDID],"SDID")  || 
        strcmp(VName[SDVarSDTyp],"SDTyp") ||
        strcmp(VName[SDVarSDN],"SDN") ||
        strcmp(VName[SDVarSDPtr],"SDPtr")) {
        if (opt)
            printf1("\nError: need valid definition of standard variables.\n");
        return(-1);
    }
    if (VSLen[SDVarSDPtr] != 5) {
        if (opt)
            printf1("\nError: SDPtr must have storage type 5.\n");
        return(-1);
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  sdnvar_close()  If SDVarDef != 0 close the file.                        */

void sdnvar_close(void)
{
    sdnvar_alloc(0,0);
    if (SDVarFDef)
        fclose(SDVarFd);
    SDVarNT = SDVarMax = SDVarFDef = SDVarDef = 0;
}

/* -##--------------------------------------------------------------------- */
/*  sdnvar_alloc(opt,n)                                                     */
/*                                                                          */
/*  If opt != 0 allocate standard array for spatial objects with n          */
/*  elements, otherwise free previously allocated memory.                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdnvar_alloc(int opt,int n)
{

    if (SDVarN > 0) {
        free((char *)SDVarX);
        free((char *)SDVarY);
        memrq(-2 * SDVarN,sizeof(double));
        SDVarN = 0;
    }
    if (opt && n > 0) {
        if (!(SDVarX = (double *)calloc(n,sizeof(double)))) {
            return(-1);
        }
        if (!(SDVarY = (double *)calloc(n,sizeof(double)))) {
            free((char *)SDVarX);
            return(-1);
        }
        memrq(2 * n,sizeof(double));
        SDVarN = n;
    }
    return(0);
}


