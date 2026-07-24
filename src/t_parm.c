/****************************************************************************/
/*  t_parm                                                                  */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-94 Goetz Rohwer. All rights reserved.           */
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
#include "t_alloc.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_plot.h"
#include "t_ml.h"
#include "t_gmin.h"
#include "t_mat.h"
#include "t_rand.h"
#include "t_imat.h"
#include "t_gf.h"

/*  functions in t_parm.c */

void p_clean(void);
void p_clean1(void);
void p_fclose(void);
int parm(char *p,int opt,int rs);
void prn_perr(char *s);
void pmfn_free(void);
int sve_alloc(int n,char *s);
int rhstr_alloc(int n);
int pmstr_alloc(int n);
int pmf1_alloc(int n);
int pmf2_alloc(int n);
int pmf3_alloc(int n);
int rhs_alloc(int n);
int prc_alloc(int n);
int prcn_alloc(int n);
int pm_valloc(int n);
int pm_zalloc(int n);
int pm_v1alloc(int n);
int pm_v2alloc(int n);
int pm_v3alloc(int n);
char *get_tp(char *tp,int *err,int opt);
void free_tp(void);        
char *get_tp1(char *tp,int *err,int opt);
void free_tp1(void);        
char *get_cn(char *tp,int *err);
void free_cn(void);        
char *get_flags(char *tp,int *err,int opt);
int pmps_alloc(int n,int m);
char *get_pattern(char *p,int *err);
char *get_fn(char *p,char *fname,int *err);
char *get_var(char *p,int *err);
char *get_var1(char *p,int *err);
char *get_var2(char *p,int *err);
char *get_var3(char *p,int *err);
char *get_varx(char *p,int *err);
char *get_xp(char *s,int *err);
void free_xp(void);        
char *get_file(int typ,char *p,int *err);
char *get_box(char *s,int *err);
void free_box(void);        
char *get_mf(char *p,int *err);
char *get_mp(char *p,int l,int *mdef,char *name,int *err);
char *get_mod(char *p,int *err);
void free_mod(void);             
char *get_fmt(int k,char *p,int *err);
int get_fmt_alloc(int k,int n);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

/*  input parameter of commands                                             */

#define PMFNMax 100             /* max strings in fn parameter              */
char *PMFNam[PMFNMax];
int PMFNA[PMFNMax];
int PMFN = 0;

FILE *PMFd;                     /* right-hand side file                     */
char PMFdName[FNMaxLen + 1];    /* name of right-hand side file             */
int PMFDef = 0;                 /* set if file opened                       */

FILE *PMF1d;                    /* df= file (always output file)            */
char PMF1dName[FNMaxLen + 1];   
int PMF1Def = 0;                

FILE *PMPCFd;                   /* pcf= file (always output file)           */
char PMPCFName[FNMaxLen + 1];   
int PMPCFDef = 0;                

FILE *PMF2d;                    /* if= file (always input file)             */
char PMF2dName[FNMaxLen + 1];   
int PMF2Def = 0;                
int PMIFTyp = 0;                /* if(string)=...                           */

FILE *PMIF1d;                   /* if1= file (always input file)            */
char PMIF1Name[FNMaxLen + 1];   
int PMIF1Def = 0;                

FILE *PMIF2d;                   /* if2= file (always input file)            */
char PMIF2Name[FNMaxLen + 1];   
int PMIF2Def = 0;                

FILE *PMCovFd;                  /* pcov= or pcova=                          */
char PMCovFName[FNMaxLen + 1];  /* name of file                             */
int PMCovFDef = 0;              /* set if file opened                       */
int PMCovWFlg = 0;              /* set if file written                      */

FILE *PMResFd;                  /* pres= or presa=                          */
char PMResFName[FNMaxLen + 1];  /* name of file                             */
int PMResFDef = 0;              /* set if file opened                       */

FILE *PMPPFd;                   /* ppar= or ppara=                          */
char PMPPFName[FNMaxLen + 1];   /* name of file                             */
int PMPPFDef = 0;               /* set if file opened                       */
int PMPPWFlg = 0;               /* set if file written                      */

FILE *PMTabFd;                  /* ptab= or ptaba=                          */
char PMTabFName[FNMaxLen + 1];  /* name of file                             */
int PMTabFDef = 0;              /* set if file opened                       */

FILE *PMTab1Fd;                 /* ptab1= or ptab1a=                       */
char PMTab1FName[FNMaxLen + 1]; /* name of file                            */
int PMTab1FDef = 0;             /* set if file opened                      */

FILE *PMProtFd;                 /* prot=                                    */
char PMProtFName[FNMaxLen + 1]; /* name of file                             */
int PMProtFDef = 0;             /* set if file opened                       */

FILE *PMPFNFd;                  /* pfn()= or pfna()=                        */
char PMPFNFName[FNMaxLen + 1];  /* name of file                             */
int PMPFNFDef = 0;              /* set if file opened                       */

FILE *PMTDAFd;                  /* dtda=                                    */
char PMTDAFName[FNMaxLen + 1];  /* name of file                             */
int PMTDAFDef = 0;              /* set if file opened                       */

FILE *PMSPSSFd;                 /* dspss=                                   */
char PMSPSSFName[FNMaxLen + 1]; /* name of file                             */
int PMSPSSFDef = 0;             /* set if file opened                       */

FILE *PMDVARFd;                 /* dvar=                                    */
char PMDVARFName[FNMaxLen + 1]; /* name of file                             */
int PMDVARFDef = 0;             /* set if file opened                       */
int PMDVARFN = 1;               /* logical file number                      */
int PMDVARP = 0;                /* p option                                 */
int PMDVARPN = 0;               /* pn option                                */
int PMDVARVN = 0;               /* number of variables                      */
short *PMDVARVNP;               /* var option: 0 (p), 1 (pn)                */
short *PMDVARVNL;               /* size of substring                        */
char **PMDVARVName;             /* name of variable                         */
int PMDVARVNPA = 0;
int PMDVARVNLA = 0;
int PMDVARVNameA = 0;

FILE *PMARCFd;                  /* arcd=                                    */
char PMARCFName[FNMaxLen + 1];  /* name of file                             */
int PMARCFDef = 0;              /* set if file opened                       */
int PMARCFZOO = 0;             /* zoo = ...                                */
char PMARCFZOOF[FNMaxLen + 1];   
int PMARCFVDF = 0;              /* vdf = ...                                */
char PMARCFVDFF[FNMaxLen + 1];   

int MFN = 0;                    /* number of merge files, mf=...            */
int MFI1[MFMAX];
int MFI2[MFMAX];
int MFJ1[MFMAX];
int MFJ2[MFMAX];
FILE *MFFD[MFMAX];
char MFFNAME[MFMAX][FNMaxLen+1];

int PMMPCovDef = 0;             /* set for mpcov=...                        */
char PMMPCovName[VNLMax+1];     /* corresponding matrix name                */
int PMMPParDef = 0;             /* set for mppar=...                        */
char PMMPParName[VNLMax+1];     /* corresponding matrix name                */
int PMMPLogDef = 0;             /* set for mplog=...                        */
char PMMPLogName[VNLMax+1];     /* corresponding matrix name                */
int PMMPGradDef = 0;            /* set for mpgrad=...                       */
char PMMPGradName[VNLMax+1];    /* corresponding matrix name                */
int PMMPResDef = 0;             /* set for mpres=...                        */
char PMMPResName[VNLMax+1];     /* corresponding matrix name                */

int PMRHSFlg = 0;               /* if right-hand side given                 */
int PMRHSI = 0;                 /* right-hand side integer argument         */
double PMRHSA = 0.0;            /* right-hand side: first double            */  
double PMRHSB = 0.0;            /* right-hand side: second double           */  
double PMRHSD = 0.0;            /* right-hand side: increment               */  

int PMRHSN = 0;                 /* number of right-hand side values         */
double *PMRHSX;                 /* array of right-hand side values          */
int PMRHSXA = 0;                /* allocated elements in PMRHSX[]           */

char *PMRHSTR;                  /* right-hand side string                   */
int PMRHSTRA = 0;               /* allocated                                */

int PMOPT = 1;                  /* opt=                                     */
int PMDOPT = -1;                /* dopt=                                    */
int PMMETH = -1;                /* meth=                                    */
int PMPRNO = 0;                 /* prn=  print option                       */
int PMPLOT = 0;                 /* plot= plot option                        */
int PMCT = 0;                   /* ct=...                                   */

char SEPC = ' ';                /* Seperation character for print formats   */
char XSEPC = ' ';             
int PMFmt1 = 0;                 /* print format: fmt=                       */
int PMFmt2 = 0;
int PMFmtF = 0;
char PMFmtS[20];        

int PMTFmt1 = 10;               /* print format: tfmt=                      */
int PMTFmt2 = 4;
int PMTFmtF = 0;
char PMTFmtS[20];        

int PMMFmt1 = 12;               /* print format: mfmt=                      */
int PMMFmt2 =  4;
int PMMFmtF = 0;
char PMMFmtS[20];        

int PMNFmt  =  0;               /* print format: nfmt=                      */
int PMNFmtF = 0;
char PMNFmtS[20];        

int PMPFmt1 = -19;              /* print format: pfmt=                      */
int PMPFmt2 =  11;
int PMPFmtF = 0;
char PMPFmtS[20];        

int PMSDFmt1 = 12;              /* print format: sdfmt=                     */
int PMSDFmt2 = 4;  
int PMSDFmtF = 0;
char PMSDFmtS[20];        

int PMXFmtN[3];                 /* fmt0=..., fmt1=..., fmt2=...              */
int PMXFmtLen = 12;
int *PMXFmt1[3];
int *PMXFmt2[3];
char *PMXFmtS[3];
int PMXFmt1A[3];   
int PMXFmt2A[3];   
int PMXFmtSA[3];   

int PMAttr = 0;                 /* attr=...                                 */
int PMMaxCat = 0;               /* maxcat=  def NOC                         */
int PMMaxCatFlg = 0;            /* set if maxcat used                       */
int PMNOC = 1000;               /* noc=                                     */
int PMNOCFlg = 0;
int PMGLEN = 0;                 /* glen=                                    */
int PMLEN = -1;                 /* len= (sequence length)                   */
int PMSN = -1;                  /* sn=  (sequence number)                   */
int PMSN1 = -1;        
int PMMSG = -1;                 /* msg=...                                  */
int PMNS = -1;                  /* ns=                                      */
int PMNDIM = -1;                /* ndim=                                    */
int PMNC =  0;                  /* nc=                                      */
int PMMin = -1;                 /* min=                                     */
int PMMax = -1;                 /* max=                                     */
int PMN = 0;                    /* n=                                       */
int PMWF = 0;                   /* wf=                                      */  
int PMPCheck = 1;               /* pcheck=                                  */  

int PMNV = 0;                   /* number of variables in PMVIdx[]          */
int PMNVTyp = 0;                /* 0 if variable names, 1 if matrix names   */
short *PMVIdx;                  /* indices of variables: v=varlist          */
int PMNZ = 0;                   /* number of variables in PMZIdx[]          */
short *PMZIdx;                  /* indices of variables: v= (...),(...)     */
int PMNVLEN = 0;                /* max length of var names in varlist       */

int PM1NV = 0;                  /* number of variables in PM1VIdx[]         */
short *PM1VIdx;                 /* indices of variables: v=varlist          */
int PM1NVLEN = 0;               /* max length of var names in varlist       */

int PM2NV = 0;                  /* number of variables in PM2VIdx[]         */
short *PM2VIdx;                 /* indices of variables: v=varlist          */
int PM2NVLEN = 0;               /* max length of var names in varlist       */

int PM3NV = 0;                  /* number of variables in PM3VIdx[]         */
short *PM3VIdx;                 /* indices of variables: v=varlist          */
int PM3NVLEN = 0;               /* max length of var names in varlist       */

int PMDBlockV = -1;             /* variable defined with dblock=            */
int PMID = -1;                  /* variable defined with id=                */
int PMORG = -1;                 /* variable defined with org=               */
int PMDES = -1;                 /* variable defined with des=               */
int PMTS = -1;                  /* variable defined with ts=                */
int PMTF = -1;                  /* variable defined with tf=                */

int PMYL = -1;                  /* variable defined with yl=                */
int PMYH = -1;                  /* variable defined with yh=                */
int PMCEN = -1;                 /* variable defined with cen=               */
int PMTRUNC = -1;               /* variable defined with trunc=             */
int PMSCAL =  -1;               /* variable defined with scale=             */

                                /* t=, tp=, qo=, qt=, wt=                   */
int PMNTP = 0;                  /* number of time points in PMTP[]          */
double *PMTP;                   /* time points                              */
int PMNTP1 = 0;                 /* number of time points in PMTP1[]         */
double *PMTP1;                  /* time points                              */
int PMQOFlg = 0;                /* set for qo=                              */
int PMQTFlg = 0;                /* set for qt=                              */

int PMNX = 0;                   /* number of xp= parameter                  */
double *PMXX;                   
double *PMXA;                   
double *PMXB;                   

int PMBoxN = 0;                 /* number of box= parameter                 */
double *PMBoxL;                   
double *PMBoxU;                   

int PMNCN = 0;                  /* number of elements in PMCN[]             */
int *PMCN;                      /* cn=                                      */
int PMDM[11];                   /* dm=    flags                             */
int PMSM[11];                   /* sm=    flags                             */
int PMTST[11];                  /* tst=   flags                             */
int PMREL[11];                  /* rel=   flags                             */
int PMSK[11];                   /* sk=    flags                             */
int PMR = 0;                    /* r=                                       */
int PMS = 0;                    /* s=                                       */
int PMSD = 0;                   /* sd=                                      */
double PMSC = 0.0;              /* sc=                                      */
double PMIC = 0.0;              /* ic=                                      */
int PMSCFlg = 0;                /* set if sc=    used                       */
int PMICFlg = 0;                /* set if ic=    used                       */
int PMLEVEL = 0;                /* level=                                   */
int PMNLEV = 0;                 /* nlev=                                    */
int PMPROJ = 0;                 /* proj=                                    */
int PMViewFlg = 0;              /* set if view is specified                 */
double PMViewLon = 0.0;         /* view=                                    */
double PMViewLat = 0.0;
double PMRHem = 90.0;           /* rhem=                                    */  

int PMXOrg = 100;               /* psorg=                                   */
int PMYOrg = 100;
double PMPSRot = 0.0;           /* psrot=                                   */

double PMAlpha = 0.0;           /* alpha =                                  */
double PMBeta  = 0.0;           /* beta  =                                  */
double PMGamma = 0.0;           /* gamma =                                  */

double PMIDFA = -1.0;           /* idf=alpha,beta                           */
double PMIDFB = -1.0;

int PMNPS = 0;                  /* number of test patterns                  */
short *PMPS[PMPSMax];           /* test patterns                            */
int PMPSN[PMPSMax];             /* length of test patterns                  */

int PMGIdx = -1;                /* index of variable defined with g=        */
double PMSIG = -1.0;            /* sig=                                     */
double PMOFF =  0.0;            /* off=                                     */

int PMRXFlg = 0;                /* if rx = a (d) b defined                  */
float PMRXA = 0.0;           
float PMRXB = 0.0;
float PMRXD = 0.0;
int PMRYFlg = 0;                /* if ry = a (d) b defined                  */
float PMRYA = 0.0;           
float PMRYB = 0.0;
float PMRYD = 0.0;

int PMRRN = 0;                  /* rr=n or rr=n,m                           */
int PMRRM = 0;

double PMRDA = 0.0;             /* rd=a[,b]                                 */  
double PMRDB = 0.0;             

int PMPLFlg = 0;                /* pl=    (plot flag)                       */
int PMLT = 1;                   /* lt=    (line type)                       */
int PMLTFlg = 0; 
int PMLT1 = 1;                  /* lt1=   (second line type)                */
double PMLW = 0.2;              /* lw=    (line width)                      */
int PMLWFlg = 0;
double PMLW1 = 0.05;            /* lw1=   (second line width)               */
int PMLW1Flg = 0;

double PMFSX = 0.0;             /* fsx=   (font size)                       */
double PMFSY = 0.0;             /* fsx=   (font size)                       */

double PMFS = 2.0;              /* fs=    (font size)                       */
int PMFSFlg = 0;
double PMFSS = 2.0;             /* fss=   (font size)                       */
int PMFSSFlg = 0;
double PMTL = 0.0;              /* tl=    (tick length)                     */

int PMLOG = 0;                  /* log=                                     */
int PMDIR = 0;                  /* dir=                                     */

int PMSizeFlg = 0;              /* set for size=...                         */
double PMSize = 0.0;
int PMOrder = 0;                /* order=...                                */
int PMRegionFlg = 0;            /* set for region=...                       */
double PMRegion1 = 0.0;
double PMRegion2 = 0.0;
int PMXYFlg = 0;                /* set for xy= or ab= or pos=               */
double PMX = 0.0;
double PMY = 0.0;
int PMXYZFlg = 0;               /* set for xyz=                             */
double PM3X = 0.0;
double PM3Y = 0.0;
double PM3Z = 0.0;
int PMDVECFlg = 0;              /* set for dvec=                            */
double PMDVECX = 0.0;
double PMDVECY = 0.0;
double PMDVECZ = 0.0;
int PMGEO = 0;                  /* geo=...                                  */
int PMHIDE = 0;                 /* hide=...                                 */
int PMRUFlg = 0;                /* set for ru=...                           */
double PMRUA = 0.0;
double PMRUB = 0.0;
int PMRUN = 0;
int PMRUM = 0;
int PMRVFlg = 0;                /* set for rv=...                           */
double PMRVA = 0.0;
double PMRVB = 0.0;
int PMRVN = 0;
int PMRVM = 0;
int PMNP = 0;                   /* np = ...                                 */
double PMULX = 0.0;             /* ulx=...                                  */
double PMULY = 0.0;             /* uly=...                                  */
double PMDX = 0.0;              /* dx=...                                   */
double PMDY = 0.0;              /* dy=...                                   */
int PMDXFlg = 0;
double PMDXA = 0.0;             /* dxa=...                                  */
int PMDXAFlg = 0;
int PMRows = 0;                 /* rows=...                                 */
int PMCols = 0;                 /* cols=...                                 */
double PMZMin = 0.0;            /* zmin=...                                 */
double PMZVal = 0.0;            /* zval=...                                 */
int PMZVar = -1;                /* variable defined with zvar=              */
double PMTol = 1.e-4;           /* tol=...                                  */

double PMGS = 0.0;              /* gs=...[,...]                             */
double PMGS1 = 0.0;     
int PMGSFlg = 0;                /* 1 if one, 2 if two arguments             */
int PMCONT = 0;                 /* cont=...                                 */

double PMA1 = 0.0;              /* a=a1,a2                                  */
double PMA2 = 0.0;      
int PMAFlg = 0;    

int PMNNFlg = 0;    
int PMNN1 = 0;                  /* nn=n1,n2                                 */
int PMNN2 = 0;      
int PMAGE1 = -1;                /* age=a1,a2                                */
int PMAGE2 = -1;
int PMYEAR1 = -1;               /* year=y1,y2                               */
int PMYEAR2 = -1;

double PMD = 0.0;               /* d=                                       */

char *PRC;                      /* rc=      (recode information)            */
int PRCAlloc = 0;               /* if PRC allocated                         */

char *PRCN;                     /* rcn=     (recode information)            */
int PRCNAlloc = 0;              /* if PRCN allocated                        */

int SVEFlg = 0;                 /* set if sel expression available          */
short SVECnt = 0;               /* parser stack for sel expression          */
int *SVETyp;
double *SVEVal;
int SVETypA = 0;
int SVEValA = 0;
char *SVESTR;                   /* string with select expression            */
int SVESTRA = 0;                /* if allocated                             */
int PMatNameFlg = 0;            /* set if mdef=...                          */
char PMatName[VNLMax + 1];      /* name defined with mdef                   */
int PMALG = 0;                  /* select algorithm                         */
int PMPERM = 0;                 /* perm=...                                 */

int PMRECFlg = 0;               /* set for rec=...                          */
double PMRECXMin = 0.0;
double PMRECYMin = 0.0;
double PMRECXMax = 0.0;
double PMRECYMax = 0.0;

/* ------------------------------------------------------------------------ */
int PMNW = 1;                   /* nw=    (number of waves)                 */
int PMBlock = 0;                /* block=                                   */  
int PMNI = 0;                   /* ni=                                      */  
int PMNQ = 0;                   /* nq=                                      */  
int PMNConS = 0;                /* con=  number of con= expressions         */
int NCONSTR = 0;                /* number of constraint expressions         */
/* ------------------------------------------------------------------------ */
double PMCFrac = 0.5;           /* cfrac=                                   */
int PMCSF = 0;                  /* csf                                      */
int PMNXA = 0;                  /* number of xa(... strings                 */
int PMDSVFlg = 0;               /* set if dsv=                              */
char PMDSVName[FNMaxLen + 1];   /* name of dsv file                         */
/* ------------------------------------------------------------------------ */
int PMPMN = 0;                  /* number of panel miss value codes         */
double PMPMVal[10];             /* pmval= (panel miss value)                */
int PMPMin = 1;                 /* pmin= (min number of participation)      */
int PMRes[10];                  /* res = r1,r2,...                          */
int PMResN = 0;                 /* number of r1,r2,... in PMRes             */
/* ------------------------------------------------------------------------ */
double PMNDIGIT = 15.0;         /* ndigit                                   */
int PMNDIGFlg = 0;
/* ------------------------------------------------------------------------ */
int PMRRFlg = 0;                /* set by rrisk                             */
int PMPRN = 0;                  /* number of prate strings                  */
/* ------------------------------------------------------------------------ */
double PMRERR = 1.e-4;          /* rerr=                                    */
int PMRERRFlg = 0;
double PMAERR = 1.e-4;          /* aerr=                                    */
int PMAERRFlg = 0;
int PMM = 1;                    /* m=                                       */
int PMMFlg = 0;                 /* set if m=...                             */
int PMM1 = 1;                   /* m1=                                      */
int PMM1Flg = 0;                /* set if m1=...                            */
int PMM2 = 1;                   /* m2=                                      */
int PMM2Flg = 0;                /* set if m2=...                            */
int PMDEG = 0;                  /* deg=                                     */
double PMKGam = 1.0;            /* kgam=                                    */
int PMKGamFlg = 0;
/* ------------------------------------------------------------------------ */
int PMKeep = 0;                 /* set by keep=varlist                      */
int PMDrop = 0;                 /* set by drop=varlist                      */
int PMTransp = 0;               /* set by transp                            */
int PMArcDic = 0;               /* set by arcdic                            */
double PMMSYS = -5.0;           /* msys=...                                 */
int PMSORTFlg = 0;              /* set by sort,                             */
int PMAP = 0;                   /* ap=...                                   */
int DGRPFlg = 0;                /* set by dgrp option                       */
/* ------------------------------------------------------------------------ */
int PMGT = 0;                   /* gt=                                      */
int PMGTT = 1;                  /* gt(n)=                                   */
int PMRT = 0;                   /* rt=                                      */
double PMMR = 2.0;              /* mr=                                      */
int PMFTYP5 = 0;                /* set if function contains type 5 var      */
/* ------------------------------------------------------------------------ */
int PMYWVar = -1;               /* variable specified with yw=...           */
int PMWVar = -1;                /* variable specified with w=...            */
int PMNBOX = -1;                /* max number of boxes                      */
int PMGN = 0;                   /* number of points for g_min3()            */
int PMGN1 = 0;      
int PMGNK = 0;                  /* number of nearest neighbors              */
double PMGD = 1.0;              /* d for random search, II                  */
int PMEVSN = -1;                /* ev = [sn,j,k]                            */
int PMEV1 = -1; 
int PMEV2 = -1; 
int XEFlg = 0;                  /* xe=...                                   */
/* ------------------------------------------------------------------------ */
double PMICOSTA = -1.0;         /* icost = ...                              */
double PMICOSTB = -1.0;     
int PMICOSTMAT = -1;  
int PMSCOSTM = -1;              /* scost = ...                              */
int PMSCOSTMAT = -1;  
/* ------------------------------------------------------------------------ */
int PMNHP = 6;                  /* nhp = ...                                */
int PMNMPnt = 1;                /* nmp=...                                  */
int PMPTyp = 1;                 /* ptyp=                                    */
int PMTyp = 0;                  /* typ=                                     */
double PMEPS = 1.e-6;           /* eps=                                     */
int PMEPSFlg = 0;
/* ------------------------------------------------------------------------ */
char *PMSTR;                    /* str=...                                  */
int PMSTRA = 0;                 /* allocated                                */
/* ------------------------------------------------------------------------ */
int PMLINK = 0;                 /* link=                                    */
int PMMXCYC = -1;               /* mxcyc=...                                */
int PMIV = 0;                   /* iv=                                      */
int PMMIX = 0;                  /* mix=...                                  */
int PMK = 0;                    /* k=...                                    */
int PML0 = 0;                   /* l0=...                                   */
/* ------------------------------------------------------------------------ */
int PMXLenFlg = 0;
int PMYLenFlg = 0;
double PMXLen = 120.0;          /* pxlen=...                                */
double PMYLen =  80.0;          /* pylen=...                                */
int PMSEED = 0;                 /* seed=...                                 */  
int PMDIM = 0;                  /* dim=...                                  */  
int PMCG = 0;                   /* cg=...                                   */  
/* ------------------------------------------------------------------------ */
int PMModN = 0;                 /* number of model strings (mod=...)        */
char *PMModS[MaxLL];            /* strings                                  */
int PMModSA[MaxLL];             /* length if allocated                      */
int SCRNFlg = 0;                /* set if screen parameter                  */
/* ------------------------------------------------------------------------ */
char *PMF1;                     /* f1=...  function expression              */
char *PMF2;                     /* f2=...  function expression              */
char *PMF3;                     /* f3=...  function expression              */
int PMF1A = 0;                  /* allocated                                */
int PMF2A = 0;                  /* allocated                                */
int PMF3A = 0;                  /* allocated                                */

/* ------------------------------------------------------------------------ */
/*  p_clean()   Free all previously allocated memory and close files.       */

void p_clean(void)
{
    p_clean1();
    a_clean();          /* free all memory in t_alloc */
}

/* ------------------------------------------------------------------------ */
/*  p_clean1()   Free all previously allocated memory and close files.      */

void p_clean1(void)
{
    int n;

    p_fclose();
    pmfn_free();
    sve_alloc(0,NULL);
    rhstr_alloc(0);
    pmstr_alloc(0);
    pmf1_alloc(0);
    pmf2_alloc(0);
    pmf3_alloc(0);
    rhs_alloc(0);
    prc_alloc(0);
    prcn_alloc(0);
    pm_valloc(0);
    pm_zalloc(0);
    pm_v1alloc(0);
    pm_v2alloc(0);
    pm_v3alloc(0);
    free_tp();
    free_tp1();
    free_xp();
    free_box();
    free_cn();
    pmps_alloc(-1,0);
    get_func(NULL,0,&n,0,NULL);
    alloc_vl(0);
    rdmn_free();
    free_mod();             
    get_fmt_alloc(0,0);
    get_fmt_alloc(1,0);
    get_fmt_alloc(2,0);
}

/* ------------------------------------------------------------------------ */
/*  p_fclose()      Close files (if open)                                   */

void p_fclose(void)
{
    register int i,l;

    if (PMFDef) {
        fclose(PMFd);
        PMFDef = 0;
    }
    if (PMF1Def) {
        fclose(PMF1d);
        PMF1Def = 0;
    }
    if (PMF2Def) {
        fclose(PMF2d);
        PMIFTyp = PMF2Def = 0;
    }
    if (PMCovFDef) {
        fclose(PMCovFd);
        PMCovWFlg = PMCovFDef = 0;
    }
    if (PMResFDef) {
        fclose(PMResFd);
        PMResFDef = 0;
    }
    if (PMPPFDef) {
        fclose(PMPPFd);
        PMPPWFlg = PMPPFDef = 0;
    }
    if (PMProtFDef) {
        fclose(PMProtFd);
        PMProtFDef = 0;
    }
    if (PMPFNFDef) {
        fclose(PMPFNFd);
        PMPFNFDef = 0;
    }
    if (PMTabFDef) {
        fclose(PMTabFd);
        PMTabFDef = 0;
    }
    if (PMTab1FDef) {
        fclose(PMTab1Fd);
        PMTab1FDef = 0;
    }
    if (PMTDAFDef) {
        fclose(PMTDAFd);
        PMTDAFDef = 0;
    }
    if (PMSPSSFDef) {
        fclose(PMSPSSFd);
        PMSPSSFDef = 0;
    }
    if (PMPCFDef) {
        fclose(PMPCFd);
        PMPCFDef = 0;
    }
    if (PMDVARFDef) {
        fclose(PMDVARFd);
        PMDVARFDef = PMDVARP = PMDVARPN = 0;
        if (PMDVARVN > 0) {
            if (PMDVARVNameA > 0) {
                for (i = 0; i < PMDVARVN; ++i) {
                    if (PMDVARVNP[i]) {
                        l = strlen(PMDVARVName[i]) + 1;
                        free(PMDVARVName[i]);
                        memrq(-l,sizeof(char));
                    }
                }
                free((char *)PMDVARVName);
                memrq(-PMDVARVNameA,sizeof(char *));
                PMDVARVNameA = 0;
            }
            if (PMDVARVNPA > 0) {
                free((char *)PMDVARVNP);
                memrq(-PMDVARVNPA,sizeof(int));
                PMDVARVNPA = 0;
            }
            if (PMDVARVNLA > 0) {
                free((char *)PMDVARVNL);
                memrq(-PMDVARVNLA,sizeof(int));
                PMDVARVNLA = 0;
            }
        }
        PMDVARVN = 0;
    }
    if (PMARCFDef) {
        fclose(PMARCFd);
        PMARCFDef = PMARCFZOO = PMARCFVDF = 0;
    }
    for (i = 0; i < MFN; ++i)  
        fclose(MFFD[i]);
    MFN = 0;
}

/* ------------------------------------------------------------------------ */
/*  parm(p,opt,rs)      Get parameters for command p having the form        */
/*                      (.....)=...                                         */
/*                                                                          */
/*                      opt = 0     right-hand side not allowed             */
/*                      opt = 1     right-hand side output file             */
/*                      opt = 2     right-hand side input file              */
/*                      opt = 3     right-hand side integer                 */
/*                      opt = 4     right-hand side varlist                 */
/*                      opt = 5     right-hand = xa,xb (double)             */
/*                      opt = 6     right-hand = xa(d)xb (double)           */
/*                      opt = 7     right-hand = x1,x2,...,xn (double)      */
/*                      opt = 8     right-hand side: string                 */
/*                      opt = 9     right-hand side: function               */
/*                      opt = 10    right-hand side: read binary            */
/*                      opt = 11    right-hand side: expression             */
/*                      opt = 12    right-hand side: inclusion function     */
/*                      opt = 13    right-hand side: write binary           */
/*                      opt = 14    right-hand side: list of variable or    */
/*                                  matrix names (but not mixed).           */
/*                                                                          */
/*                      If rs = 1 there must be a right-hand side argument. */
/*                      Return: 0 if OK, otherwise error.                   */

int parm(char *p,int opt,int rs)
{
    register int i;
    int n,m,m1,err,r,fnd;
    register char c,*q;
    char vname[VNLMax + 1];
    char *p0;
    double a,b,d;

    p0 = p;
    err = -1;
    p_clean1();                 /* clean previously allocated arrays */

    PMIFTyp = 0;                /* if(string)=... */
    SEPC = ' ';                 /* separation character */
    XSEPC = ' ';    
    NCONSTR = 0;                /* number of constraint expressions */
    PMAttr = 0;                 /* attr=... */
    PMMaxCat = 1000;            /* maxcat= */
    PMMaxCatFlg = 0;
    PMRHSFlg = 0;               /* set if right-hand side given */
    PMRHSI = 0;                 /* right-hand side integer */
    PMRHSA = 0.0;               /* right-hand side: double */
    PMRHSB = 0.0;               /* right-hand side: double */
    PMRHSD = 0.0;               /* right-hand side: increment */
    PMOPT = 1;                  /* opt =  */
    PMDOPT = -1;                /* dopt =  */
    PMMETH = -1;                /* meth= */
    PMPRNO  = 0;                /* prn =  */
    PMPLOT = 0;                 /* plot= */
    PMCT = 0;                   /* ct=... */

    PMFmt1  = 0;                /* fmt */
    PMFmt2  = 0;           
    PMFmtF  = 0;    

    PMTFmt1 = 10;               /* tfmt */
    PMTFmt2 = 4;           
    PMTFmtF = 0;        

    PMMFmt1 = 12;               /* mfmt */
    PMMFmt2 = 4;           
    PMMFmtF = 0;        

    PMNFmt  = 4;                /* nfmt */
    PMNFmtF = 0;        

    PMPFmt1 = -19;              /* pfmt */
    PMPFmt2 =  11;         
    PMPFmtF = 0;        

    PMSDFmt1 = 12;              /* sdfmt */
    PMSDFmt2 = 4;           
    PMSDFmtF = 0;        

    PMNOC = 1000;               /* noc= */
    PMNOCFlg = 0;
    PMGLEN = 0;                 /* glen= */
    PMLEN = -1;                 /* len= sequence length */
    PMSN = -1;                  /* sn= sequence number */
    PMSN1 = -1; 
    PMMSG = -1;                 /* msg= */
    PMNS = -1;                  /* ns=  */
    PMNDIM = -1;                /* ndim=  */
    PMNC =  0;                  /* nc= */
    PMMin = -1;                 /* min= */
    PMMax = -1;                 /* max= */
    PMWF = 0;                   /* wf= */  
    PMPCheck = 1;               /* pcheck= */  

    PMDBlockV = -1;             /* dblock= */
    PMID  = -1;                 /* id= */
    PMORG = -1;                 /* org= */
    PMDES = -1;                 /* des= */
    PMTS  = -1;                 /* ts= */
    PMTF  = -1;                 /* tf= */

    PMYL = -1;                  /* yl= */
    PMYH = -1;                  /* yh= */
    PMCEN = -1;                 /* cen= */
    PMTRUNC = -1;               /* trunc= */
    PMSCAL  = -1;               /* scale= */

    PMQOFlg = 0;                /* set for qo= */
    PMQTFlg = 0;                /* set for qt= */
   
    PMN = 0;                    /* n= */
    PMM = 1;                    /* m= */
    PMMFlg = 0;                 /* set if m=... */
    PMM1 = 1;                   /* m1= */
    PMM1Flg = 0;                /* set if m1=... */
    PMM2 = 1;                   /* m2= */
    PMM2Flg = 0;                /* set if m2=... */
    PMR = 0;                    /* r= */
    PMS = 0;                    /* s= */
    PMSD = 0;                   /* sd= */
    PMSC = 0.0;                 /* sc= */
    PMIC = 0.0;                 /* ic= */
    PMSCFlg = PMICFlg = 0;
    PMNLEV = 0;                 /* nlev=                                    */
    PMLEVEL = 0;                /* level=                                   */
    PMPROJ = 0;                 /* proj= */
    PMViewFlg = 0;              /* set if view is specified */
    PMViewLon = 0.0;            /* view= */
    PMViewLat = 0.0;
    PMRHem = 90.0;              /* rhem=  */  
    PMXOrg = 100;               /* psorg= */
    PMYOrg = 100;
    PMPSRot = 0.0;              /* psrot= */

    PMAlpha = 0.0;              /* alpha = */
    PMBeta  = 0.0;              /* beta  = */
    PMGamma = 0.0;              /* gamma = */
    PMIDFA = PMIDFB = -1.0;     /* idf=alpha,beta */
    PMGIdx = -1;                /* index of variable defined with g= */
    PMSIG = -1.0;               /* sig= */
    PMOFF =  0.0;               /* off= */
    PMRXFlg = 0;                /* if rx = a (d) b defined */
    PMRXA = 0.0;           
    PMRXB = 0.0;
    PMRXD = 0.0;
    PMRYFlg = 0;                /* if rx = a (d) b defined */
    PMRYA = 0.0;           
    PMRYB = 0.0;
    PMRYD = 0.0;
    PMRRN = 0;                  /* rr=n or rr=n,m */
    PMRRM = 0;
    PMLOG = 0;                  /* log= */
    PMDIR = 0;                  /* dir= */
    PMSizeFlg = 0;              /* set for size=... */
    PMSize = 0.0;
    PMOrder = 0;                /* oder=... */
    PMRegionFlg = 0;            /* set for region=... */
    PMRegion1 = 0.0;
    PMRegion2 = 0.0;
    PMXYFlg = 0;                /* set for xy= */
    PMX = 0.0;
    PMY = 0.0;
    PMD = 0.0;                  /* d= */
    PMXYZFlg = 0;               /* set for xyz= */
    PM3X = 0.0;
    PM3Y = 0.0;
    PM3Z = 0.0;                                
    PMDVECFlg = 0;              /* set for dvec= */
    PMDVECX = 0.0;
    PMDVECY = 0.0;
    PMDVECZ = 0.0;                          
    PMGEO = 0;                  /* geo=... */
    PMHIDE = 0;                 /* hide=... */
    PMRUFlg = 0;                /* set for ru=... */
    PMRUA = 0.0;
    PMRUB = 0.0;
    PMRUN = 0;
    PMRUM = 0;
    PMRVFlg = 0;                /* set for rv=...*/
    PMRVA = 0.0;
    PMRVB = 0.0;
    PMRVN = 0;
    PMRVM = 0;
    PMNP = 0;                   /* np = ... */
    PMULX = 0.0;                /* ulx=... */
    PMULY = 0.0;                /* uly=... */
    PMDX = 0.0;                 /* dx=...  */
    PMDY = 0.0;                 /* dy=...  */
    PMDXFlg = 0;
    PMDXA = 0.0;                /* dxa=... */
    PMDXAFlg = 0;
    PMRows = 0;                 /* rows=... */
    PMCols = 0;                 /* cols=... */
    PMZMin = 0.0;               /* zmin=... */
    PMZVal = 0.0;               /* zval=... */
    PMZVar = -1;                /* zvar= */
    PMTol = 1.e-4;              /* tol=... */

    PMRDA = 0.0;                /* rd= */
    PMRDB = 0.0;                              
    PMPLFlg = 0;                /* pl= (plot flag) */
    PMLT = 1;                   /* lt= (line type) */
    PMLT1 = 1;                  /* lt1= (line type) */
    PMLW = LWDef;               /* lw= (line width) */
    PMLW1 = LW1Def;             /* lw1= (line width) */

    PMFSX = 0.0;                /* fsx= (font size) */
    PMFSY = 0.0;                /* fsx= (font size) */

    PMFSS = PMFS = FSDef;       /* fss= fs=... (font size) */
    PMTL = 0.0;                 /* tl= (tick length) */
    PMLTFlg = PMLWFlg = PMLW1Flg = PMFSFlg = PMFSSFlg = 0;

    PMGS = 0.0;                 /* gs=  (grey scale value) */
    PMGS1 = 0.0; 
    PMGSFlg = 0;
    PMCONT = 0;                 /* cont=... */

    PMA1 = PMA2 = 0.0;          /* a= a1,a2 */
    PMAFlg = 0;
    SVEFlg = 0;
    FNFlg = 0;
    PMNNFlg = PMNN1 = PMNN2 = 0;    /* nn=n1,n2 */
    PMAGE1 = PMAGE2 = -1;           /* age=a1,a2 */
    PMYEAR1 = PMYEAR2 = -1;         /* year=y1,y2 */

    PMNConS = 0;                /* number of con= expressions */
    PMNW = 1;                   /* nw= (number of waves) */
    PMBlock = 0;                /* block= */  
    PMNI = 0;                   /* ni= */
    PMNQ = 0;                   /* nq= */

    PMCSF = 0;                  /* csf */
    PMNXA = 0;                  /* number of xa(... strings */
    PMDSVFlg = 0;               /* set if dsv= */
    PMCFrac = 0.5;              /* cfrac= */
    PMPMN = 0;                  /* pmval= number of codes */
    PMPMin = 1;                 /* pmin= (min number of participation) */
    PMResN = 0;                 /* res=... */

    PMNDIGIT = 15.0;            /* ndigit */
    PMNDIGFlg = 0;
    PMRRFlg = 0;                /* set by rrisk */
    PMPRN = 0;                  /* number of prate strings */
    PMRERR = 1.e-4;             /* rerr= */
    PMRERRFlg = 0;
    PMAERR = 1.e-4;             /* aerr= */
    PMAERRFlg = 0;
    PMDEG = 0;                  /* deg= */
    PMKGam = 1.0;               /* kgam= */
    PMKGamFlg = 0;
    PMKeep = 0;                 /* set by keep=varlist */
    PMDrop = 0;                 /* set by drop=varlist */
    PMTransp = 0;               /* set by transp */
    PMArcDic = 0;               /* set by arcdic */
    PMMSYS = -5.0;              /* msys= */
    PMSORTFlg = 0;              /* set by sort, */
    PMAP = 0;                   /* ap= */
    DGRPFlg = 0;                /* set by dgrp option */
    PMGT = 0;                   /* gt= */
    PMGTT = 1;                  /* gtt */
    PMRT = 0;                   /* rt= */
    PMMR = 2.0;                 /* mr= */
    PMFTYP5 = 0;                /* set if function contains type 5 var */
    PMYWVar = -1;               /* variable specified with yw=... */
    PMWVar = -1;                /* variable specified with w=... */
    PMNBOX = -1;                /* max number of boxes */
    PMGN = 0;                   /* number of points for g_min3() */
    PMGN1 = 0;
    PMGNK = 0;                  /* number of nearest neighbors */
    PMGD = 1.0;                 /* gd=  d for random search, II */
    PMEVSN = -1;                /* ev = [sn,j,k] */
    PMEV1 = -1;     
    PMEV2 = -1; 
    XEFlg = 0;                  /* xe=... */
    PMICOSTA = -1.0;            /* icost = ...  */
    PMICOSTB = -1.0;     
    PMICOSTMAT = -1;  
    PMSCOSTM = -1;              /* scost = ... */
    PMSCOSTMAT = -1;  
    PMNHP = 6;                  /* nhp = ... */
    PMNMPnt = 1;                /* nmp=... */
    PMPTyp = 1;                 /* ptyp=... */
    PMTyp = 0;                  /* typ= ... */
    PMEPS = 1.e-6;              /* eps= */
    PMEPSFlg = 0;
    PMLINK = 0;                 /* link= */
    PMMXCYC = -1;               /* mxcyc=... */
    PMIV = 0;                   /* iv= */
    PMMIX = 0;                  /* mix=... */  
    PMK = 0;                    /* k=... */
    PML0 = 0;                   /* l0=... */
    PMXLenFlg = 0;
    PMYLenFlg = 0;
    PMXLen = 120.0;             /* pxlen=... */
    PMYLen =  80.0;             /* pylen=... */
    PMMPCovDef = 0;             /* set for mpcov=... */
    PMMPParDef = 0;             /* set for mppar=... */
    PMMPLogDef = 0;             /* set for mplog=... */
    PMMPGradDef = 0;            /* set for mpgrad=... */
    PMMPResDef = 0;             /* set for mpres=... */
    PMatNameFlg = 0;            /* set if mdef=...  */
    PMALG = 0;                  /* select algorithm */
    PMSEED = 0;                 /* seed=... */
    PMDIM = 0;                  /* dim=... */
    PMCG = 0;                   /* cg=... */
    PMPERM = 0;                 /* perm=... */
    SCRNFlg = 0;                /* screen */
    PMRECFlg = 0;               /* set for rec=... */


    for (i = 1; i <= 10; ++i)
        PMSK[i] = PMSM[i] = PMDM[i] = PMTST[i] = PMREL[i] = 0;

    if (*p == '(') {

        while (1) {
            if (*++p == ')')
                break;

            fnd = 1;
            err = -1;
            p0 = p;

            switch (*p) {

                case 'a':

                    if (sscanf(p,"a=%lf,%lf",&a,&b) == 2) {
                        PMA1 = a;
                        PMA2 = b;
                        PMAFlg = 1;
                        p = skip_dbl(p + 2);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"ab=%lf,%lf",&a,&b) == 2) {
                        PMX = a;
                        PMY = b;
                        PMXYFlg = 1;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"amue=%lg",&a) == 1 && a > 0.0 && a < 1.0) {
                        AMue = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"aerr=%lg",&a) == 1 && a >= 0.0) {
                        PMAERR = a;
                        PMAERRFlg = 1;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"ap=%d",&n) == 1) {
                        PMAP = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"alg=%d",&n) == 1 && n > 0) {
                        PMALG = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"arcdic",6)) {
                        PMArcDic = 1;
                        p += 6;
                    }
                    else if (!strncmp(p,"arcd",4)) {    /* PMARC file */
                        p = get_file(13,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"alpha=%lf",&a) == 1) {
                        PMAlpha = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (!strncmp(p,"av=",3)) {
                        p = get_var1(p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"attr=%d",&n) == 1 && n >= 1) {
                        PMAttr = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"age=%d,%d",&n,&m) == 2 && n >= 0 && m >= n) {
                        PMAGE1 = n;
                        PMAGE2 = m;
                        p = skip_int(p + 4);
                        p = skip_int(p + 1);
                    }   
                    else
                        fnd = 0;
                    break;

                case 'b':

                    if (!strncmp(p,"box=",4)) {
                        p = get_box(p + 4,&r);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"beta=%lf",&a) == 1) {
                        PMBeta = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"block=%d",&n) == 1 && n >= 0) {
                        PMBlock = n;
                        p = skip_int(p + 6);
                    }
                    else
                        fnd = 0;
                    break;

                case 'c':
                    if (!strncmp(p,"cn=",3)) {
                        p = get_cn(p + 3,&r);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"crit=%d",&n) == 1 && n >= 1 && n <= 3) {
                        Crite = n;
                        CritFlg = 1;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"ccov=%d",&n) == 1 && n >= 1 && n <= 3) {
                        CCTyp = n;
                        p = skip_int(p + 5);
                    }
                    else if (!strncmp(p,"csf,",4)) {     
                        PMCSF = 1;
                        p = skip_int(p + 3);
                    }   
                    else if (sscanf(p,"cfrac=%lg",&a) == 1 && a >= 0.0 && a <= 1.0) {
                        PMCFrac = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (!strncmp(p,"con=",4)) {     
                        PMNConS++; 
                        p = skip_com(p);
                    }
                    else if (!strncmp(p,"cen=",4)) {
                        if ((PMCEN = get_vidx1(p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (sscanf(p,"cg=%d",&n) == 1) {
                        PMCG = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"cont=%d",&n) == 1) {
                        PMCONT = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"cols=%d",&n) == 1 && n > 0) {
                        PMCols = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"ct=%d",&n) == 1 && n >= 0) {
                        PMCT = n;
                        p = skip_int(p + 3);
                    }
                    else
                        fnd = 0;
                    break;

                case 'd':
                    if (!strncmp(p,"des=",4)) {
                        if ((PMDES = get_vidx1(p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (!strncmp(p,"drop=",5)) {
                        p = get_var(p + 5,&r);
                        if (r)  
                            goto PARMFin;
                        PMDrop = 1;
                    }
                    else if (!strncmp(p,"dm=",3)) {
                        p = get_flags(p + 3,&r,1);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"deg=%d",&n) == 1 && n >= 0) {
                        PMDEG = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"df",2)) {     /* PMF1d output file */
                        p = get_file(1,p + 2,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"dir=%d",&n) == 1) {
                        PMDIR = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"d=%lg",&a) == 1) {
                        PMD = a;
                        p = skip_dbl(p + 2);
                    }
                    else if (sscanf(p,"dscal=%lg",&a) == 1 && a != 0.0) {
                        DScal = a;
                        DScalFlg = 1;
                        p = skip_dbl(p + 6);
                    }
                    else if (!strncmp(p,"dsv",3)) {    /* dsv file */
                        p = get_fn(p + 3,PMDSVName,&r);
                        if (r)  
                            goto PARMFin;
                        PMDSVFlg = 1;
                    }
                    else if (sscanf(p,"dopt=%d",&n) == 1 && n >= 0) {
                        PMDOPT = n;
                        p = skip_int(p + 5);
                    }
                    else if (!strncmp(p,"dtda",4)) {    /* PMTDA file */
                        p = get_file(10,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"dspss",5)) {   /* PMSPSS file */
                        p = get_file(11,p + 5,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"dvar",4)) {   /* PMDVAR file */
                        p = get_file(12,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"dim=%d",&n) == 1) {
                        PMDIM = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"dvec=%lf,%lf,%lf",&a,&b,&d) == 3) {
                        PMDVECX = a;
                        PMDVECY = b;
                        PMDVECZ = d;
                        PMDVECFlg = 1;
                        p = skip_dbl(p + 5);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                    }
                    else if (!strncmp(p,"dgrp=",5)) {
                        n = 0;
                        q = p + 4;
                        while (*(q + 1) == '[') {
                            q = skip_nc(q + 1);
                            if (*(q - 1) != ']')
                                goto PARMFin;
                            n++;
                        }
                        if (q <= p + 7)
                            goto PARMFin;

                        c = *q;
                        *q = '\0';
                        n = strlen(p + 5);
                        if (n == 0)
                            goto PARMFin;

                        if (rhstr_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(PMRHSTR,p + 5);
                        *q = c;
                        p = q;
                        DGRPFlg = n;
                    }
                    else if (sscanf(p,"dx=%lg",&a) == 1) {
                        PMDX = a;
                        PMDXFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"dy=%lg",&a) == 1) {
                        PMDY = a;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"dxa=%lg",&a) == 1) {
                        PMDXA = a;
                        PMDXAFlg = 1;
                        p = skip_dbl(p + 4);
                    }
                    else if (!strncmp(p,"dblock=",7)) {
                        if ((PMDBlockV = get_vidx1(p + 7,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 7;
                    }
                    else
                        fnd = 0;
                    break;

                case 'e':
                    if (sscanf(p,"eps=%lg",&a) == 1 && a > 0.0) {
                        PMEPS = a;
                        PMEPSFlg = 1;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"ev=[%d,%d,%d]",&n,&m,&m1) == 3 && n >= 1 && m >= 0 && m1 >= 0) {
                        PMEVSN = n;
                        PMEV1 = m;
                        PMEV2 = m1;
                        p = skip_int(p + 4);
                        p = skip_int(p + 1);
                        p = skip_int(p + 1) + 1;
                    }   
                    else
                        fnd = 0;
                    break;
                
                case 'f':
                    if (sscanf(p,"fmt=%d.%d",&n,&m) == 2) {
                        PMFmt1 = n;
                        PMFmt2 = m;
                        p = skip_dbl(p + 4);
                        PMFmtF = 1;
                    }   
                    else if (sscanf(p,"fs=%lg",&a) == 1 && a >= 0.0) {
                        PMFS = a;
                        PMFSFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"fss=%lg",&a) == 1 && a >= 0.0) {
                        PMFSS = a;
                        PMFSSFlg = 1;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"fsx=%lg",&a) == 1 && a >= 0.0) {
                        PMFSX = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"fsy=%lg",&a) == 1 && a >= 0.0) {
                        PMFSY = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (!strncmp(p,"fn=",3)) {
                        p += 3;
                        q = skip_nc(p);
                        c = *q;
                        *q = '\0';
                        n = strlen(p);
   
                        if (PMFN < PMFNMax) {
                            if (!(PMFNam[PMFN] = (char *)calloc(n + 1,sizeof(char)))) {
                                err = -2;    
                                goto PARMFin;
                            }
                            memrq(n + 1,sizeof(char));
                            PMFNA[PMFN] = n + 1;
                            strcpy(PMFNam[PMFN],p);   
                            PMFN++; 
                        }
                        *q = c;
                        p = q;
                    }
                    else if (!strncmp(p,"fmt0=",5)) {
                        p = get_fmt(0,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"fmt1=",5)) {
                        p = get_fmt(1,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"fmt2=",5)) {
                        p = get_fmt(2,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"f1=,",3)) {
                        p += 3;
                        q = skip_expr(p);        
                        n = (int)(q - p);
                        if (n == 0)   
                            goto PARMFin;

                        if (pmf1_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(PMF1,p,n);
                        p = q;
                    }
                    else if (!strncmp(p,"f2=,",3)) {
                        p += 3;
                        q = skip_expr(p);        
                        n = (int)(q - p);
                        if (n == 0)   
                            goto PARMFin;

                        if (pmf2_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(PMF2,p,n);
                        p = q;
                    }
                    else if (!strncmp(p,"f3=,",3)) {
                        p += 3;
                        q = skip_expr(p);        
                        n = (int)(q - p);
                        if (n == 0)   
                            goto PARMFin;

                        if (pmf3_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(PMF3,p,n);
                        p = q;
                    }
                    else
                        fnd = 0;
                    break;

                case 'g':
                    if (!strncmp(p,"g=",2)) {
                        p += 2;
                        PMGIdx = get_vidx1(p,vname);
                        if (PMGIdx < 0) {
                            err = -4;
                            goto PARMFin;
                        }
                        p = skip_com(p);
                    }
                    else if (sscanf(p,"gd=%lg",&a) == 1 && a > 0.0) {
                        PMGD = a;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"gs=%lg,%lg",&a,&b) == 2) {
                        PMGS = a;
                        PMGS1 = b;
                        PMGSFlg = 2;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"gs=%lg",&a) == 1) {
                        PMGS = a;
                        PMGSFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"gt(%d)=%d",&m,&n) == 2 &&
                                    m >= 1 && m <= 3 && n >= 1 && n <= 4) {
                        PMGT = n;
                        PMGTT = m;
                        p = skip_int(p + 3);
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"gt=%d",&n) == 1 && n >= 1 && n <= 4) {
                        PMGT = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"gn=%d,%d",&n,&m) == 2 && n >= 0 && m >= 0) {
                        PMGN = n;
                        PMGN1 = m;
                        p = skip_int(p + 3);
                        p = skip_int(p + 1);
                    }
                    else if (sscanf(p,"gn=%d",&n) == 1 && n >= 0) {
                        PMGN = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"gnk=%d",&n) == 1 && n >= 1) {
                        PMGNK = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"gss=",4)) {
                        p = get_tp(p + 4,&r,2);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"grp=",4)) {
                        p = get_var1(p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"geo=%d",&n) == 1) {
                        PMGEO = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"gamma=%lf",&a) == 1) {
                        PMGamma = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"glen=%d",&n) == 1 && n >= 0) {
                        PMGLEN = n;
                        p = skip_int(p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'h':
                    if (sscanf(p,"hide=%d",&n) == 1) {
                        PMHIDE = n;
                        p = skip_int(p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'i':
                    if (!strncmp(p,"id=",3)) {
                        if ((PMID = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"ic=%lf",&a) == 1) {
                        PMIC = a;
                        PMICFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"idf=%lf,%lf",&a,&b) == 2) {
                        PMIDFA = a;
                        PMIDFB = b;
                        p = skip_dbl(p + 4);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"iv=%d",&n) == 1) {
                        PMIV = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"idf=%lf",&a) == 1) {
                        PMIDFA = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (!strncmp(p,"if=",3) || !strncmp(p,"if(",3)) {     /* PMF2d (input file) */

                        p += 2;
                        if (*p == '(') {
                            p++;
                            if (!strncmp(p,"keep",4)) {
                                PMIFTyp = 1;
                                p += 4;
                            }
                            else if (!strncmp(p,"drop",4)) {
                                PMIFTyp = 2;
                                p += 4;
                            }
                            else  
                                goto PARMFin;
                            if (*p++ != ')')
                                goto PARMFin;
                        }
                        p = get_fn(p,PMF2dName,&r);
                        if (r)  
                            goto PARMFin;

                        if (!(PMF2d = fopen(PMF2dName,OPEN_RD))) {   
                            printf1("Error: can't open: %s\n",PMF2dName);
                            err = 1;
                            goto PARMFin;
                        }
                        PMF2Def = 1;
                    }
                    else if (!strncmp(p,"icost=",6)) {
                        if (sscanf(p,"icost=%lf,%lf",&a,&b) == 2 && a >= 0.0 && b >= 0.0) {
                            PMICOSTA = a;
                            PMICOSTB = b;
                            p = skip_dbl(p + 6);
                            p = skip_dbl(p + 1);
                        }
                        else if (sscanf(p,"icost=%lf",&a) == 1 && a >= 0.0) {
                            PMICOSTA = a;
                            p = skip_dbl(p + 6);
                        }
                        else {
                            p = get_mname(p + 6,vname,0);
                            if (p == NULL)
                                goto PARMFin;
                            n = mat_getidx(vname,0);
                            if (n < 0) {
                                printf1("Error: undefined matrix name.\n");
                                goto PARMFin;
                            }
                            PMICOSTMAT = n;  
                        }
                    }
                    else if (!strncmp(p,"if1",3)) {    /* PMIF1d input file */
                        p = get_file(14,p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"if2",3)) {    /* PMIF2d input file */
                        p = get_file(15,p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;

                case 'k':

                    if (sscanf(p,"k=%d",&n) == 1) {
                        PMK = n;
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"kgam=%lg",&a) == 1 && a > 0.0) {
                        PMKGam = a;
                        PMKGamFlg = 1;
                        p = skip_dbl(p + 5);
                    }
                    else if (!strncmp(p,"keep=",5)) {
                        p = get_var(p + 5,&r);
                        if (r)  
                            goto PARMFin;
                        PMKeep = 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 'l':
                    if (sscanf(p,"len=%d",&n) == 1) {
                        PMLEN = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"link=%d",&n) == 1 && n >= 0) {
                        PMLINK = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"lt=%d",&n) == 1 && n >= 0) {
                        PMLT = n;
                        PMLTFlg = 1;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"lt1=%d",&n) == 1 && n >= 0) {
                        PMLT1 = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"lw=%lg",&a) == 1 && a >= 0.0) {
                        PMLW = a;
                        PMLWFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"lw1=%lg",&a) == 1 && a >= 0.0) {
                        PMLW1 = a;
                        PMLW1Flg = 1;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"log=%d",&n) == 1) {
                        PMLOG = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"lsecon",6) || !strncmp(p,"lsicon",6)) {     
                        NCONSTR++; 
                        p = skip_com(p);
                    }
                    else if (sscanf(p,"l0=%d",&n) == 1) {
                        PML0 = n;
                        p = skip_int(p + 3);
                    }
                    else if (!strncmp(p,"lon=",4)) {
                        p = get_tp(p + 4,&r,0);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"lat=",4)) {
                        p = get_tp1(p + 4,&r,0);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"level=%d",&n) == 1) {
                        PMLEVEL = n;
                        p = skip_int(p + 6);
                    }
                    else
                        fnd = 0;
                    break;

                case 'm':
                    if (sscanf(p,"m=%d",&n) == 1) {
                        PMM = n;
                        PMMFlg = 1;
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"m1=%d",&n) == 1) {
                        PMM1 = n;
                        PMM1Flg = 1;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"m2=%d",&n) == 1) {
                        PMM2 = n;
                        PMM2Flg = 1;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"mix=%d",&n) == 1) {
                        PMMIX = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"mr=%lf",&a) == 1 && a >= 1.0) {
                        PMMR = a;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"maxcat=%d",&n) == 1 && n > 0) {
                        PMMaxCat = n;
                        PMMaxCatFlg = 1;
                        p = skip_int(p + 7);
                    }
                    else if (sscanf(p,"min=%d",&n) == 1 && n >= 0) {
                        PMMin = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"max=%d",&n) == 1 && n >= 0) {
                        PMMax = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"mfmt=%d.%d",&n,&m) == 2) {
                        PMMFmt1 = n;
                        PMMFmt2 = m;
                        p = skip_dbl(p + 5);
                        PMMFmtF = 1;
                    }   
                    else if (sscanf(p,"mina=%d",&n) == 1 && n >= 1 && n <= 8) {
                        MINA = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"mxit=%d",&n) == 1 && n >= 0) {
                        MxIter = n;
                        MxItFlg = 1;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"mxcyc=%d",&n) == 1 && n >= 0) {
                        PMMXCYC = n;
                        p = skip_int(p + 6);
                    }
                    else if (sscanf(p,"mxitl=%d",&n) == 1 && n >= 0) {
                        MxIt1 = n;
                        p = skip_int(p + 6);
                    }
                    else if (sscanf(p,"msys=%lg",&a) == 1) {
                        PMMSYS = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (!strncmp(p,"mf=",3)) {
                        p = get_mf(p + 3,&r);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpcov",5)) {  
                        p = get_mp(p,5,&PMMPCovDef,PMMPCovName,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mppar",5)) {  
                        p = get_mp(p,5,&PMMPParDef,PMMPParName,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mplog",5)) {  
                        p = get_mp(p,5,&PMMPLogDef,PMMPLogName,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpgrad",6)) {  
                        p = get_mp(p,6,&PMMPGradDef,PMMPGradName,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpres",5)) {  
                        p = get_mp(p,5,&PMMPResDef,PMMPResName,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mdef=",5)) {  
                        p = get_mname(p + 5,PMatName,0);
                        if (p == NULL || (*p != ',' && *p != ')'))
                            goto PARMFin;
                        PMatNameFlg = 1;
                    }
/* ## */            else if (!strncmp(p,"mod=",4)) {  
                        p = get_mod(p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"msg=%d",&n) == 1 && n >= 0) {
                        PMMSG = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"meth=%d",&n) == 1 && n >= 1) {
                        PMMETH = n;
                        p = skip_int(p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'n':
                    if (sscanf(p,"noc=%d",&n) == 1 && n >= 1) {
                        PMNOC = n;
                        PMNOCFlg = 1;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"nfmt=%d",&n) == 1 && n > 0) {
                        PMNFmt = n;
                        p = skip_int(p + 5);
                        PMNFmtF = 1;
                    }   
                    else if (sscanf(p,"n=%d",&n) == 1 && n >= 1) {
                        PMN = n;
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"nhp=%d",&n) == 1 && n >= 1) {
                        PMNHP = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"nmp=%d",&n) == 1 && n >= 1) {
                        PMNMPnt = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"ns=%d",&n) == 1 && n >= 0) {
                        PMNS = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"ndim=%d",&n) == 1 && n >= 1) {
                        PMNDIM = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"nc=%d",&n) == 1) {
                        PMNC = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"nw=%d",&n) == 1 && n >= 1) {
                        PMNW = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"nq=%d",&n) == 1 && n >= 1) {
                        PMNQ = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"ni=%d",&n) == 1 && n >= 0) {
                        PMNI = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"nbox=%d",&n) == 1 && n >= 1) {
                        PMNBOX = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"nn=%d,%d",&n,&m) == 2) {
                        PMNNFlg = 1;
                        PMNN1 = n;
                        PMNN2 = m;
                        p = skip_int(p + 3);
                        p = skip_int(p + 1);
                    }   
                    else if (sscanf(p,"ndigit=%lg",&a) == 1 && a >= 1.0) {
                        PMNDIGIT = a;
                        PMNDIGFlg = 1;
                        p = skip_dbl(p + 7);
                    }
                    else if (sscanf(p,"np=%d",&n) == 1) {
                        PMNP = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"nlev=%d",&n) == 1 && n >= 1) {
                        PMNLEV = n;
                        p = skip_int(p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'o':
                    if (sscanf(p,"opt=%d",&n) == 1 && n >= 1) {
                        PMOPT = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"org=",4)) {
                        if ((PMORG = get_vidx1(p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (sscanf(p,"order=%d",&n) == 1) {
                        PMOrder = n;
                        p = skip_int(p + 6);
                    }
                    else if (sscanf(p,"off=%lf",&a) == 1) {
                        PMOFF = a;
                        p = skip_dbl(p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'p':
                    if (!strncmp(p,"ps=",3)) {
                        p = get_pattern(p + 3,&r);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"prn=%d",&n) == 1 && n >= 0) {
                        PMPRNO = n;
                        p = skip_int(p + 4);
                    }
                    else if (sscanf(p,"plot=%d",&n) == 1 && n >= 0) {
                        PMPLOT = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"ptyp=%d",&n) == 1 && n >= 1) {
                        PMPTyp = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"pos=%lf,%lf",&a,&b) == 2) {
                        PMX = a;
                        PMY = b;
                        PMXYFlg = 1;
                        p = skip_dbl(p + 4);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"pl=%d",&n) == 1 && n >= 0) {
                        PMPLFlg = n;
                        p = skip_int(p + 3);
                    }
                    else if (!strncmp(p,"prate",5)) {     
                        p += 5;
                        if (*p == '(') 
                            p = skip_blev(p);
                             
                        if (*p++ != '=')
                            goto PARMFin;
                        p = skip_com(p);
                        PMPRN++; 
                    }
                    else if (!strncmp(p,"prot",4)) {  /* PMProt file */
                        p = get_file(5,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ppar",4)) {    /* PMPPFd file */
                        p = get_file(4,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pcov",4)) {    /* PMCov file */
                        p = get_file(2,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pres",4)) {    /* PMRes file */
                        p = get_file(3,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pfn",3)) {    /* PMPFN file */
                        p = get_file(7,p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"pfmt=%d.%d",&n,&m) == 2) {
                        PMPFmt1 = n;
                        PMPFmt2 = m;
                        p = skip_dbl(p + 5);
                        PMPFmtF = 1;
                    }   
                    else if (sscanf(p,"pmin=%d",&n) == 1 && n >= 1) {
                        PMPMin = n;
                        p = skip_int(p + 5);
                    }
                    else if (!strncmp(p,"pmval=",6)) {
                        p += 5;
                        while (sscanf(p + 1,"%lf",&a) == 1) {
                            PMPMVal[PMPMN] = a;
                            p = skip_dbl(p + 1);
                            if (++PMPMN >= 10)
                                break;
                        }
                    }
                    else if (!strncmp(p,"ptab1",5)) {    /* PMTab1 file */
                        p = get_file(9,p + 5,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ptab",4)) {    /* PMTab file */
                        p = get_file(8,p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"pxlen=%lg",&a) == 1 && a > 0.0) {
                        PMXLen = a;
                        PMXLenFlg = 1;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"pylen=%lg",&a) == 1 && a > 0.0) {
                        PMYLen = a;
                        PMYLenFlg = 1;
                        p = skip_dbl(p + 6);
                    }
                    else if (!strncmp(p,"pcf",3)) {     /* PMPCFd output file */
                        p = get_file(6,p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (sscanf(p,"perm=%d",&n) == 1) {
                        PMPERM = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"proj=%d",&n) == 1 && n >= 1) {
                        PMPROJ = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
                        PMXOrg = n;
                        PMYOrg = m;
                        p = skip_int(p + 6);
                        p = skip_int(p + 1);
                    }
                    else if (sscanf(p,"psrot=%lg",&a) == 1) { /*      && a >= 0.0 && a < 360.0) {   */
                        PMPSRot = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"pcheck=%d",&n) == 1 && n >= 0) {
                        PMPCheck = n;
                        p = skip_int(p + 7);
                    }
                    else
                        fnd = 0;
                    break;

                case 'q':
                    if (!strncmp(p,"qo=",3)) {
                        p = get_tp(p + 3,&r,3);      
                        if (r)  
                            goto PARMFin;
                        PMQOFlg = 1;
                    }
                    else if (!strncmp(p,"qt=",3)) {
                        p = get_tp(p + 3,&r,1);      
                        if (r)
                            goto PARMFin;
                        PMQTFlg = 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 'r':
                    if (sscanf(p,"r=%d",&n) == 1) {
                        PMR = n;
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"rd=%lf,%lf",&a,&b) == 2) {
                        PMRDA = a;
                        PMRDB = b;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"rd=%lf",&a) == 1) {
                        PMRDA = a;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"rt=%d",&n) == 1 && n >= 0) {
                        PMRT = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"rr=%d,%d",&n,&m) == 2 && n >= 1 && m >= n) {
                        PMRRN = n;
                        PMRRM = m;
                        p = skip_int(p + 3);
                        p = skip_int(p + 1);
                    }
                    else if (sscanf(p,"rr=%d",&n) == 1 && n >= 1) {
                        PMRRM = PMRRN = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"rx=%lf(%lf)%lf",&a,&d,&b) == 3 && a < b && d > 0.0) {
                        PMRXFlg = 1;
                        PMRXA = (float)a;
                        PMRXB = (float)b;
                        PMRXD = (float)d;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"ry=%lf(%lf)%lf",&a,&d,&b) == 3 && a < b && d > 0.0) {
                        PMRYFlg = 1;
                        PMRYA = (float)a;
                        PMRYB = (float)b;
                        PMRYD = (float)d;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"ru=%lf,%lf,%d,%d",&a,&b,&n,&m) == 4 && a < b && n > 0 && m > 0) {
                        PMRUFlg = 1;
                        PMRUA = a;
                        PMRUB = b;
                        PMRUN = n;          
                        PMRUM = m;          
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                        p = skip_int(p + 1);
                        p = skip_int(p + 1);
                    }
                    else if (sscanf(p,"rv=%lf,%lf,%d,%d",&a,&b,&n,&m) == 4 && a < b && n > 0 && m > 0) {
                        PMRVFlg = 1;
                        PMRVA = a;
                        PMRVB = b;
                        PMRVN = n;          
                        PMRVM = m;          
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                        p = skip_int(p + 1);
                        p = skip_int(p + 1);
                    }
                    else if (!strncmp(p,"rcn=",4)) {
                        p += 4;
                        q = p;
                        while (*q) {
                            q = skip_int(q);
                            if (*q++ != '[')  
                                goto PARMFin;

                            while (*q) {
                                q = skip_int(q);
                                if (*q == ']')  
                                    break;
                                if (*q++ != ',')  
                                    goto PARMFin;
                                if (*q == ',')
                                    q++;
                            }
                            if (*q++ != ']')  
                                goto PARMFin;

                            if (*q != ',' || sscanf(q + 1,"%d",&n) != 1)
                                break;
                            q++;
                        }
                        if (q <= p)
                            goto PARMFin;

                        n = (int)(q - p);
                        if (prcn_alloc(n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(PRCN,p,n);
                        *(PRCN + n) = '\0';
                        p = q;
                    }
                    else if (!strncmp(p,"rc=",3)) {
                        p += 3;
                        q = p;
                        while (*q) {
                            if (sscanf(q,"%d",&n) != 1)
                                goto PARMFin;
                            q = skip_int(q);
                            if (*q++ != '[')  
                                goto PARMFin;

                            while (*q) {
                                q = skip_int(q);
                                if (*q == ']')  
                                    break;
                                if (*q++ != ',')  
                                    goto PARMFin;
                            }
                            if (*q++ != ']')  
                                goto PARMFin;

                            if (*q != ',' || sscanf(q + 1,"%d",&n) != 1)
                                break;
                            q++;
                        }
                        if (q <= p)
                            goto PARMFin;

                        n = (int)(q - p);
                        if (prc_alloc(n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(PRC,p,n);
                        *(PRC + n) = '\0';
                        p = q;
                    }
                    else if (!strncmp(p,"res=",4)) {
                        p += 3;
                        while (sscanf(p + 1,"%d",&n) == 1) {
                            if (PMResN < 10)
                                PMRes[PMResN++] = n;
                            p = skip_int(p + 1);
                        }
                    }
                    else if (!strncmp(p,"rrisk,",6)) {
                        PMRRFlg = 1; 
                        p += 5;
                    }
                    else if (sscanf(p,"rerr=%lg",&a) == 1 && a >= 0.0) {
                        PMRERR = a;
                        PMRERRFlg = 1;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"region=%lf,%lf",&a,&b) == 2) {
                        PMRegion1 = a;
                        PMRegion2 = b;
                        PMRegionFlg = 1;
                        p = skip_dbl(p + 7);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"rows=%d",&n) == 1 && n > 0) {
                        PMRows = n;
                        p = skip_int(p + 5);
                    }
                    else if (sscanf(p,"rec=%lf,%lf,%lf,%lf",&PMRECXMin,&PMRECYMin,&PMRECXMax,&PMRECYMax) == 4) {
                        PMRECFlg = 1;
                        p = skip_dbl(p + 4);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                    }
                    else if (!strncmp(p,"rel=",4)) {
                        p = get_flags(p + 4,&r,5);      
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"rhem=%lf",&a) == 1 && a > 0.0) {
                        PMRHem = a;
                        p = skip_dbl(p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 's':
                    if (!strncmp(p,"sm=[",4)) {
                        p += 4;
                        q = p;
                        n = 0;
                        while (*q && *q != ']') {
                            q++;
                            n++;
                        }
                        if (n == 0 || !*q == ']')
                            goto PARMFin;

                        if (rhstr_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        *q = '\0';
                        strcpy(PMRHSTR,p);
                        *q = ']';
                        p += n + 1;
                    }
                    else if (!strncmp(p,"str=,",4)) {
                        q = p + 4;
                        m = n = 0;
                        if (*q == '"') {
                            q++;
                            m = 1;
                        }
                        p0 = q;
                        while (*q) {
                            if (m == 1) {
                                if (*q == '"')
                                    break;    
                            }                 
                            else if (*q == ',' || *q == ')')
                                break;
                            n++;
                            q++;
                        }
                        c = *q;
                        *q = '\0';
                        if (pmstr_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(PMSTR,p0);
                        *q = c;
                        p = q;
                        if (m)
                            p++;
                    }
                    else if (sscanf(p,"sc=%lf",&a) == 1) {
                        PMSC = a;
                        PMSCFlg = 1;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"s=%d",&n) == 1) {
                        PMS = n;
                        p = skip_int(p + 2);
                    }
                    else if (sscanf(p,"sd=%d",&n) == 1) {
                        PMSD = n;
                        p = skip_int(p + 3);
                    }
                    else if (sscanf(p,"sn=%d,%d",&n,&m) == 2 && n >= 0 && m >= 0) {
                        PMSN = n;
                        PMSN1 = m;
                        p = skip_int(p + 3);
                        p = skip_int(p + 1);
                    }
                    else if (sscanf(p,"sn=%d",&n) == 1 && n >= 0) {
                        PMSN = n;
                        p = skip_int(p + 3);
                    }
                    else if (!strncmp(p,"sm=",3)) {
                        p = get_flags(p + 3,&r,2);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sort=",5)) {
                        p = get_var1(p + 5,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sort,",5)) {
                        PMSORTFlg = 1;
                        p += 4;
                    }
                    else if (sscanf(p,"sig=%lf",&a) == 1) {
                        PMSIG = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (!strncmp(p,"scost=",6)) {
                        if (sscanf(p,"scost=%d",&n) == 1 && n >= 0) {
                            PMSCOSTM = n;
                            p = skip_int(p + 6);
                        }
                        else {
                            p = get_mname(p + 6,vname,0);
                            if (p == NULL)
                                goto PARMFin;
                            n = mat_getidx(vname,0);
                            if (n < 0) {
                                printf1("Error: undefined matrix name.\n");
                                goto PARMFin;
                            }
                            PMSCOSTMAT = n;  
                        }
                    }
                    else if (sscanf(p,"seed=%d",&n) == 1) {
                        PMSEED = n;
                        p = skip_int(p + 5);
                    }
                    else if (!strncmp(p,"sel=",4)) {

                        q = skip_expr(p + 4);
                        c = *q;
                        *q = '\0';

                        if ((n = v_parse(p + 4,0)) < 0 || ESCnt <= 0) {
                            printf1("Syntax error (%d) in sel expression.\n",n);
                            if (n < 0)
                                prn_emsg1(n);
                            err = 1;
                            goto PARMFin;
                        }
                        if (sve_alloc(ESCnt,p)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        SVECnt = ESCnt;
                        for (i = 0; i < ESCnt; ++i) {
                            SVETyp[i] = ESTyp[i];
                            SVEVal[i] = ESVal[i];
                        }
                        SVEFlg = 1;
                        *q = c;
                        p = q;
                    }
                    else if (sscanf(p,"smin=%lg",&a) == 1 && a >= 0.0) {
                        SMin = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"slen=%lg",&a) == 1 && a > 0.0) {
                        SLen = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"sred=%lg",&a) == 1 && a > 0.0 && a < 1.0) {
                        SRed = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"sst=%d",&n) == 1 && n >= 0 && n <= 1) {
                        STFlg = n;
                        p = skip_int(p + 4);
                    }
                    else if (!strncmp(p,"sk=",3)) {
                        p = get_flags(p + 3,&r,4);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sepc=none",9)) {
                        XSEPC = '\0';
                        p += 9;
                    }
                    else if (sscanf(p,"sepc=%c",&XSEPC) == 1) {
                        if (XSEPC == 't')
                            XSEPC = '\t';
                        p += 6;                
                    }
                    else if (!strncmp(p,"screen",6)) {     
                        SCRNFlg = 1;
                        p += 6;          
                    }
                    else if (!strncmp(p,"scale=",6)) {
                        if ((PMSCAL = get_vidx1(p + 6,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 6;
                    }
                    else if (sscanf(p,"size=%lg",&a) == 1) {
                        PMSize = a;
                        PMSizeFlg = 1;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"sdfmt=%d.%d",&n,&m) == 2) {
                        PMSDFmt1 = n;
                        PMSDFmt2 = m;
                        p = skip_dbl(p + 6);
                        PMSDFmtF = 1;
                    }   
                    else
                        fnd = 0;
                    break;

                case 't':
                    if (!strncmp(p,"t=",2)) {
                        p = get_tp(p + 2,&r,1);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"tp=",3)) {
                        p = get_tp(p + 3,&r,1);      
                        if (r) 
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"tst=",4)) {
                        p = get_flags(p + 4,&r,3);      
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ts=",3)) {
                        if ((PMTS = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"tf=",3)) {
                        if ((PMTF = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"tl=%lf",&a) == 1 && a >= 0.0) {
                        PMTL = a;
                        p = skip_dbl(p + 3);
                    }
                    else if (sscanf(p,"tfmt=%d.%d",&n,&m) == 2) {
                        PMTFmt1 = n;
                        PMTFmt2 = m;
                        p = skip_dbl(p + 5);
                        PMTFmtF = 1;
                    }   
                    else if (sscanf(p,"tolg=%lg",&a) == 1 && a > 0.0) {
                        TOLG = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"tolf=%lg",&a) == 1 && a > 0.0) {
                        TOLF = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"tolp=%lg",&a) == 1 && a > 0.0) {
                        TOLP = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"tolv=%lg",&a) == 1 && a > 0.0) {
                        TOLV = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"tols=%lg",&a) == 1 && a > 0.0) {
                        TOLS = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"tolsg=%lg",&a) == 1 && a > 0.0) {
                        TOLSG = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"tolsp=%lg",&a) == 1 && a > 0.0) {
                        TOLSP = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"tolbw=%lg",&a) == 1 && a > 0.0) {
                        TOLBW = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"tolfd=%lg",&a) == 1 && a > 0.0) {
                        TOLFD = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"tolfe=%lg",&a) == 1 && a > 0.0) {
                        TOLFE = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (sscanf(p,"tolbc=%lg",&a) == 1 && a > 0.0) {
                        TOLBC = a;
                        p = skip_dbl(p + 6);
                    }
                    else if (!strncmp(p,"transp",6)) {
                        PMTransp = 1;
                        p += 6;
                    }
                    else if (!strncmp(p,"trunc=",6)) {
                        if ((PMTRUNC = get_vidx1(p + 6,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 6;
                    }
                    else if (sscanf(p,"tol=%lg",&a) == 1 && a >= 0.0) {
                        PMTol = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"typ=%d",&n) == 1 && n >= 1) {
                        PMTyp = n;
                        p = skip_int(p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'u':
                    if (sscanf(p,"ulx=%lg",&a) == 1) {
                        PMULX = a;
                        p = skip_dbl(p + 4);
                    }
                    else if (sscanf(p,"uly=%lg",&a) == 1) {
                        PMULY = a;
                        p = skip_dbl(p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'v':
                    if (!strncmp(p,"v=",2)) {
                        p = get_var(p + 2,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"vsel=",5)) {

                        q = skip_expr(p + 5);
                        c = *q;
                        *q = '\0';
                        n = strlen(p);
                        if (rhstr_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(PMRHSTR,p);
                        *q = c;
                        p = q;
                    }
                    else if (sscanf(p,"view=%lg,%lg",&a,&b) == 2 &&
                          -180.0 <= a && a <= 180.0 && -90.0 <= b && b <= 90.0) {
                        PMViewLon = a;
                        PMViewLat = b;
                        PMViewFlg = 1; 
                        p = skip_dbl(p + 5);
                        p = skip_dbl(p + 1);
                    }
                    else
                        fnd = 0;
                    break;

                case 'w':
                    if (!strncmp(p,"w=",2)) {
                        if ((PMWVar = get_vidx1(p + 2,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 2;
                    }
                    else if (sscanf(p,"wf=%d",&n) == 1 && n >= 0) {
                        PMWF = n;
                        p = skip_int(p + 3);
                    }
                    else if (!strncmp(p,"wt=",3)) {
                        p = get_tp(p + 3,&r,1);      
                        if (r) 
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;

                case 'x':
                    if (sscanf(p,"xy=%lf,%lf",&a,&b) == 2) {
                        PMX = a;
                        PMY = b;
                        PMXYFlg = 1;
                        p = skip_dbl(p + 3);
                        p = skip_dbl(p + 1);
                    }
                    else if (sscanf(p,"xyz=%lf,%lf,%lf",&a,&b,&d) == 3) {
                        PM3X = a;
                        PM3Y = b;
                        PM3Z = d;
                        PMXYZFlg = 1;
                        p = skip_dbl(p + 4);
                        p = skip_dbl(p + 1);
                        p = skip_dbl(p + 1);
                    }
                    else if (!strncmp(p,"x=",2)) {
                        p = get_tp(p + 2,&r,0);      
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) || 
                             !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3)) { 
                        p = skip_xa(p);
                        PMNXA++;
                    }
                    else if (!strncmp(p,"xp=",3)) {
                        p = get_xp(p + 3,&r);      
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xe=",3)) {
                        m = 0;
                        q = p + 2;
                        while (*(q + 1) == '[') {
                            q = skip_nc(q + 1);
                            if (*(q - 1) != ']')
                                goto PARMFin;
                            m++;
                        }
                        /*****************
                        if (q <= p + 12)
                            goto PARMFin;
                        *****************/
                        c = *q;
                        *q = '\0';
                        n = strlen(p + 3);
                        if (n == 0)
                            goto PARMFin;

                        if (rhstr_alloc(n + 1)) {     
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(PMRHSTR,p + 3);
                        *q = c;
                        p = q;
                        XEFlg = m;
                    }
                    else if (!strncmp(p,"xv=",3)) {
                        p = get_var1(p + 3,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xyv=",4)) {
                        p = get_var(p + 4,&r);
                        if (r)  
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;
        
                case 'y':
                    if (!strncmp(p,"y=",2)) {
                        p = get_tp(p + 2,&r,1);      
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"yw=",3)) {
                        if ((PMYWVar = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"yl=",3)) {
                        if ((PMYL = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"yh=",3)) {
                        if ((PMYH = get_vidx1(p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"year=%d,%d",&n,&m) == 2 && n >= 0 && m >= n) {
                        PMYEAR1 = n;
                        PMYEAR2 = m;
                        p = skip_int(p + 5);
                        p = skip_int(p + 1);
                    }   
                    else
                        fnd = 0;
                    break;

                case 'z':
                    if (sscanf(p,"zmin=%lg",&a) == 1) {
                        PMZMin = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (sscanf(p,"zval=%lg",&a) == 1) {
                        PMZVal = a;
                        p = skip_dbl(p + 5);
                    }
                    else if (!strncmp(p,"zvar=",5)) {
                        if ((PMZVar = get_vidx1(p + 5,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 5;
                    }
                    else
                        fnd = 0;
                    break;
                
                
                default:
                    fnd = 0;
                    break;
            }
            if (fnd == 0)
                goto PARMFin;
                   
            if (*p != ',')
                break;
        }
        if (*p++ != ')')
            goto PARMFin;
    }
    err = -1;
    if (*p == '=') {
        if (!*++p || opt == 0)
            goto PARMFin;
  
        if (opt <= 2 || opt == 10 || opt == 13) {
            strcpy(PMFdName,p);
            i = 0;
            if (opt == 2) {
                if (!(PMFd = fopen(p,OPEN_RD)))     
                    i = 1;
            }
            else if (opt == 10) {
                if (!(PMFd = fopen(p,OPEN_RB)))     
                    i = 1;
            }
            else if (opt == 13) {
                if (!(PMFd = fopen(p,OPEN_WB)))     
                    i = 1;
            }
            else if (PMAP == 0) {
                if (!(PMFd = fopen(p,OPEN_WR)))     
                    i = 1;
            }
            else if (!(PMFd = fopen(p,OPEN_AP)))     
                i = 1;

            if (i) {
                printf1("Error: can't open: %s\n",p);
                err = 1;
                goto PARMFin;
            }
            PMFDef = 1;
        }
        else if (opt == 3) {
            if (sscanf(p,"%d",&n) != 1)  
                goto PARMFin;
            PMRHSI = n;
        }   
        else if (opt == 4) {        /* must be a varlist */
            p = get_var(p,&err);
            if (*p)
                err = -1;
            if (err)
                goto PARMFin;
        }
        else if (opt == 5) {
            if (sscanf(p,"%lf,%lf",&a,&b) != 2)  
                goto PARMFin;
            PMRHSA = a;
            PMRHSB = b;
        }   
        else if (opt == 6) {
            if (sscanf(p,"%lf(%lf)%lf",&a,&d,&b) != 3)  
                goto PARMFin;
            PMRHSA = a;
            PMRHSB = b;
            PMRHSD = d;
        }   
        else if (opt == 7) {        /* list of double values */
            q = p;
            n = 1;
            while (*q) {
                if (*q++ == ',')
                    n++;
            }
            if (rhs_alloc(n)) {     
                err = 2;
                goto PARMFin;
            }
            for (i = 0; i < n; ++i) {
                if (sscanf(p,"%lf",&a) != 1)  
                    goto PARMFin;
                PMRHSX[i] = a;
                p = skip_dbl(p);
                if (i == n - 1)
                    break;
                if (*p++ != ',')
                    goto PARMFin;
            }
            PMRHSN = n;
        }
        else if (opt == 8) {        /* string */
            n = strlen(p);
            if (rhstr_alloc(n + 1)) {     
                err = -2;
                goto PARMFin;
            }
            strcpy(PMRHSTR,p);
        }
        else if (opt == 9) {        /* function */

            err = get_func(p,1,&n,0,NULL);
            if (err) {
                err = 1;
                goto PARMFin;
            }
            PMFTYP5 = n;
        }
        else if (opt == 11) {        /* expression */

            if ((n = v_parse(p,0)) < 0 || ESCnt <= 0) {
                printf1("Syntax error (%d) in right-hand side expression.\n",n);
                if (n < 0)
                    prn_emsg1(n);
                err = 1;
                goto PARMFin;
            }
            p = skip_expr(p);
        }
        else if (opt == 12) {        /* inclusion function */

            err = get_func(p,1,&n,1,NULL);
            if (err) {
                err = 1;
                goto PARMFin;
            }
            PMFTYP5 = n;
        }
        else if (opt == 14) { /* must be a list of variable or matrix names */
            p = get_varx(p,&err);
            if (*p)
                err = -1;
            if (err)
                goto PARMFin;
        }
        PMRHSFlg = 1; 
    }
    else if (*p)
        goto PARMFin;
    else if (rs) {
        printf1("Error: need right-hand side of command.\n");
        err = 1;
        goto PARMFin;
    }

    makefmt(&PMFmt1,&PMFmt2,PMFmtS,0,SEPC,0);
    makefmt(&PMTFmt1,&PMTFmt2,PMTFmtS,0,SEPC,0);
    makefmt(&PMMFmt1,&PMMFmt2,PMMFmtS,0,SEPC,0);
    makefmt(&PMPFmt1,&PMPFmt2,PMPFmtS,0,SEPC,0);
    makefmt(&PMSDFmt1,&PMSDFmt2,PMSDFmtS,0,SEPC,0);
    makenfmt(&PMNFmt,PMNFmtS,1,SEPC);

    if (PMNVTyp == 0) {
        PMNVLEN = imax(get_mxvlen(PMNV,PMVIdx),get_mxvlen(PMNZ,PMZIdx));
        PM1NVLEN = get_mxvlen(PM1NV,PM1VIdx);
        PM2NVLEN = get_mxvlen(PM2NV,PM2VIdx);
        PM3NVLEN = get_mxvlen(PM3NV,PM3VIdx);
    }
    err = 0;

PARMFin:                        /* ## */
    if (err < 0) {
        /***************
        p = skip_nc(p0);
        *p = '\0';
        *************/
        prn_perr(p0);
        p_err(err,1);
    }
    if (err)
        err = -1;
    return(err);

PARMFin1:               /* ## */
    p = skip_nc(p0);
    *p = '\0';
    prn_perr(p0);
    printf1("Undefined variable.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_perr    Print error message.                                        */

void prn_perr(char *s)
{
    register char *p = s;

    printf1("Error: ");
    while (*p && p < s + 20)
        printf1("%c",*p++);
    if (*p)
        printf1(" ...");
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  pmfn_free()     Free previously allocated memory PMFNam[].              */

void pmfn_free(void)
{
    register int i;

    for (i = 0; i < PMFN; ++i) {
        if (PMFNA[i] > 0) {
            free(PMFNam[i]);
            memrq(-PMFNA[i],sizeof(char));
            PMFNA[i] = 0;
        }
    }
    PMFN = 0;
}

/* ------------------------------------------------------------------------ */
/*  sve_alloc(n,s)  If n > 0 allocate parser stack for sel expression,      */  
/*                  otherwise free. Save string s in SVESTR.                */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int sve_alloc(int n,char *s)
{
    int l;

    if (SVESTRA > 0) {
        free((char *)SVESTR);
        memrq(-SVESTRA,sizeof(char));
        SVESTRA = 0;
    }
    if (SVETypA > 0) {
        free((char *)SVETyp);
        memrq(-SVETypA,sizeof(int));
        SVETypA = 0;
    }
    if (SVEValA > 0) {
        free((char *)SVEVal);
        memrq(-SVEValA,sizeof(double));
        SVEValA = 0;
    }
    SVECnt = SVEFlg = 0;
    if (n <= 0)
        return(0);

    if (!(SVETyp = (int *)calloc(n,sizeof(int))))  
        return(-2);
    memrq(n,sizeof(int));
    SVETypA = n;

    if (!(SVEVal = (double *)calloc(n,sizeof(double))))  
        return(-2);
    memrq(n,sizeof(double));
    SVEValA = n;

    if (s != NULL) {
        l = strlen(s);
        if (l > 0) {
            if (!(SVESTR = (char *)calloc(l + 1,sizeof(char))))  
                return(-2);
            memrq(l + 1,sizeof(char));
            SVESTRA = l + 1;
            strcpy(SVESTR,s);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rhstr_alloc(n)  If n > 0 allocate PMRHSTR else free previously          */  
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int rhstr_alloc(int n)
{
    if (PMRHSTRA > 0) {
        free(PMRHSTR);
        memrq(-PMRHSTRA,sizeof(char));
        PMRHSTRA = 0;
        DGRPFlg = XEFlg = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMRHSTR = (char *)calloc(n,sizeof(char))))   
        return(-2);
    memrq(n,sizeof(char));
    PMRHSTRA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmstr_alloc(n)  If n > 0 allocate PMSTR else free previously            */  
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmstr_alloc(int n)
{
    if (PMSTRA > 0) {
        free(PMSTR);
        memrq(-PMSTRA,sizeof(char));
        PMSTRA = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMSTR = (char *)calloc(n,sizeof(char))))   
        return(-2);
    memrq(n,sizeof(char));
    PMSTRA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf1_alloc(n)   If n > 0 allocate PMF1 else free previously             */  
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf1_alloc(int n)
{
    if (PMF1A > 0) {
        free(PMF1);
        memrq(-PMF1A,sizeof(char));
        PMF1A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMF1 = (char *)calloc(n,sizeof(char))))   
        return(-2);
    memrq(n,sizeof(char));
    PMF1A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf2_alloc(n)   If n > 0 allocate PMF2 else free previously             */  
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf2_alloc(int n)
{
    if (PMF2A > 0) {
        free(PMF2);
        memrq(-PMF2A,sizeof(char));
        PMF2A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMF2 = (char *)calloc(n,sizeof(char))))   
        return(-2);
    memrq(n,sizeof(char));
    PMF2A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf3_alloc(n)   If n > 0 allocate PMF3 else free previously             */  
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf3_alloc(int n)
{
    if (PMF3A > 0) {
        free(PMF3);
        memrq(-PMF3A,sizeof(char));
        PMF3A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMF3 = (char *)calloc(n,sizeof(char))))   
        return(-2);
    memrq(n,sizeof(char));
    PMF3A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rhs_alloc(n)    If n > 0 allocate PMRHSX else free previously allocated */  
/*                  PMRHSX.                                                 */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int rhs_alloc(int n)
{
    if (PMRHSXA > 0) {
        free((char *)PMRHSX);
        memrq(-PMRHSXA,sizeof(double));
        PMRHSN = PMRHSXA = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PMRHSX = (double *)calloc(n,sizeof(double))))   
        return(-2);
    memrq(n,sizeof(double));
    PMRHSXA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prc_alloc(n)    If n > 0 allocate PRC else free previously allocated    */  
/*                  PRC.                                                    */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int prc_alloc(int n)
{
    if (PRCAlloc > 0) {
        free(PRC);
        memrq(-PRCAlloc,sizeof(char));
        PRCAlloc = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PRC = (char *)calloc(n,sizeof(char))))   
        return(-2);
       
    memrq(n,sizeof(char));
    PRCAlloc = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prcn_alloc(n)   If n > 0 allocate PRCN else free previously allocated   */  
/*                  PRCN.                                                   */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int prcn_alloc(int n)
{
    if (PRCNAlloc > 0) {
        free(PRCN);
        memrq(-PRCNAlloc,sizeof(char));
        PRCNAlloc = 0;
    }
    if (n <= 0)
        return(0);

    if (!(PRCN = (char *)calloc(n,sizeof(char))))   
        return(-2);
       
    memrq(n,sizeof(char));
    PRCNAlloc = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_valloc()     If n > 0 allocate n elements in PMVIdx[], otherwise     */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_valloc(int n)
{
    if (PMNV > 0) {
        free((char *)PMVIdx);
        memrq(-PMNV,sizeof(short));
        PMNVTyp = PMNV = 0;
    }
    if (n > 0) {
        if (!(PMVIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }                 
        memrq(n,sizeof(short));
        PMNV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_zalloc()     If n > 0 allocate n elements in PMZIdx[], otherwise     */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_zalloc(int n)
{
    if (PMNZ > 0) {
        free((char *)PMZIdx);
        memrq(-PMNZ,sizeof(short));
        PMNZ = 0;
    }
    if (n > 0) {
        if (!(PMZIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }                 
        memrq(n,sizeof(short));
        PMNZ = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v1alloc()    If n > 0 allocate n elements in PM1VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v1alloc(int n)
{
    if (PM1NV > 0) {
        free((char *)PM1VIdx);
        memrq(-PM1NV,sizeof(short));
        PM1NV = 0;
    }
    if (n > 0) {
        if (!(PM1VIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }                 
        memrq(n,sizeof(short));
        PM1NV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v2alloc()    If n > 0 allocate n elements in PM2VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v2alloc(int n)
{
    if (PM2NV > 0) {
        free((char *)PM2VIdx);
        memrq(-PM2NV,sizeof(short));
        PM2NV = 0;
    }
    if (n > 0) {
        if (!(PM2VIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }                 
        memrq(n,sizeof(short));
        PM2NV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v3alloc()    If n > 0 allocate n elements in PM3VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v3alloc(int n)
{
    if (PM3NV > 0) {
        free((char *)PM3VIdx);
        memrq(-PM3NV,sizeof(short));
        PM3NV = 0;
    }
    if (n > 0) {
        if (!(PM3VIdx = (short *)calloc(n,sizeof(short)))) { 
            p_err(-2,1);
            return(-1);
        }                 
        memrq(n,sizeof(short));
        PM3NV = n;
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_tp(tp,err,opt)                                                      */
/*                                                                          */
/*                  Get time points: t=... or x=... or x=...                */
/*                  Sets the number of time points PMNTP and allocates the  */
/*                  array PMTP[].                                           */
/*                                                                          */
/*                  Syntax check depends on opt                             */
/*                  opt = 0 : any order                                     */
/*                  opt = 1 : nonnegative and strictly ascending            */
/*                  opt = 2 : 0 <= t <= 1, any order                        */
/*                  opt = 3 : 1 > t > 0 and descending                      */
/*                  Return: pointer to next char in tp string.              */
/*                                                                          */
/*                  err =  0 if successful                                  */
/*                        -1 syntax error                                   */
/*                        -2 insuff memory                                  */
/*                         1 exceeded limits (local message)                */

char *get_tp(char *tp,int *err,int opt)
{
    register int i,j,n;
    register char *p;
    int m;
    double a,a1,b,d,tmp;

    free_tp();
    *err = -1;
    p = tp;
    if (!*p)
        goto GETTPErr;

    for (i = 0; i < 2; ++i) {

        a1 = 0.0;
        p = tp;
        n = 0;

        while (*p) {
            if (sscanf(p,"%lg",&a) != 1) {
                if (n)
                    p--;
                break;
            }
            if (n == 0) {
                if (opt == 1 && a < 0.0)                 
                    goto GETTPErr;
                if (opt == 2 && a < 0.0)                 
                    goto GETTPErr;
                if (opt == 3 && a >= 1.0)                 
                    goto GETTPErr;
            }      
            else {
                if (opt == 1 && a <= a1)                 
                    goto GETTPErr;
                if (opt == 3 && (a1 <= a || a <= 0.0))     
                    goto GETTPErr;
            }
            a1 = a;
            if (i)
                PMTP[n] = a;
            n++;
            p = skip_dbl(p);

            if (*p == '(') {
                if (sscanf(++p,"%lg",&d) != 1 || d <= EPSI1)   
                    goto GETTPErr;
  
                p = skip_dbl(p);
                if (*p++ != ')')
                    goto GETTPErr;
                if (sscanf(p,"%lg",&b) != 1 || fabs(b - a) <= EPSI1)
                    goto GETTPErr;
                p = skip_dbl(p);

                if (opt == 1 && b <= a + EPSI1)
                    goto GETTPErr;
                if (opt == 2 && b > 1.0)
                    goto GETTPErr;
                if (opt == 3 && (b >= a - EPSI1 || b < EPSI1))
                    goto GETTPErr;
                if (b < a)
                    d = -d;
                   
                tmp = fabs((b - a) / d);
                if (fabs(tmp - 1.0) <= EPSI1)
                    m = 1;
                else {
                    m = (int)floor(tmp);                   
                    if (m <= 0 || m > 100000) {
                        printf1("Error: exceeded limits of t/tp/x/gss/qo/qt/lon parameter.\n");
                        printf1("or there might be an error in the increment.\n");
                        *err = 1;
                        goto GETTPErr;
                    }
                }
                if (i == 0)  
                    n += m;
                else {
                    for (j = 0; j < m; ++j) {
                        PMTP[n] = PMTP[n - 1] + d;
                        n++;
                    }
                }
                a += (double)m * d;
                if (a < b - EPSI1) {
                    if (i)
                        PMTP[n] = b;
                    n++;
                }
                a1 = b;
            }
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETTPErr;

        if (i == 0) {
            if (!(PMTP = (double *)calloc(n,sizeof(double)))) { 
                *err = -2;
                goto GETTPErr;
            }
            memrq(n,sizeof(double));
            PMNTP = n;
        }
    }
    if (n != PMNTP)
        gerr_exit(61);

    *err = 0;

    /***
    printf1("PMNTP=%d : ",PMNTP);
    for (i = 0; i < PMNTP; ++i)
    printf1("%lg ",PMTP[i]);
    printf1("\n");
    ***/
    return(p);

GETTPErr:
    free_tp();
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_tp()    Free memory which was used for time points.                */

void free_tp(void)         
{
    if (PMNTP > 0) {
        free((char *)PMTP);
        memrq(-PMNTP,sizeof(double));
        PMNTP = 0;
    }
    PMQOFlg = PMQTFlg = 0;
}

/*--##----------------------------------------------------------------------*/
/*  get_tp1(tp,err,opt)                                                     */
/*                                                                          */
/*                  Get time points: t=... or x=...                         */
/*                  Sets the number of time points PMNTP1 and allocates the */
/*                  array PMTP1[].                                          */
/*                                                                          */
/*                  Syntax check depends on opt                             */
/*                  opt = 0 : any order                                     */
/*                  opt = 1 : nonnegative and strictly ascending            */
/*                  opt = 2 : 0 <= t <= 1, any order                        */
/*                  opt = 3 : 1 > t > 0 and descending                      */
/*                  Return: pointer to next char in tp string.              */
/*                                                                          */
/*                  err =  0 if successful                                  */
/*                        -1 syntax error                                   */
/*                        -2 insuff memory                                  */
/*                         1 exceeded limits (local message)                */

char *get_tp1(char *tp,int *err,int opt)
{
    register int i,j,n;
    register char *p;
    int m;
    double a,a1,b,d,tmp;

    free_tp1();
    *err = -1;
    p = tp;
    if (!*p)
        goto GETTPErr;

    for (i = 0; i < 2; ++i) {

        a1 = 0.0;
        p = tp;
        n = 0;

        while (*p) {
            if (sscanf(p,"%lg",&a) != 1) {
                if (n)
                    p--;
                break;
            }
            if (n == 0) {
                if (opt == 1 && a < 0.0)                 
                    goto GETTPErr;
                if (opt == 2 && a < 0.0)                 
                    goto GETTPErr;
                if (opt == 3 && a >= 1.0)                 
                    goto GETTPErr;
            }      
            else {
                if (opt == 1 && a <= a1)                 
                    goto GETTPErr;
                if (opt == 3 && (a1 <= a || a <= 0.0))     
                    goto GETTPErr;
            }
            a1 = a;
            if (i)
                PMTP1[n] = a;
            n++;
            p = skip_dbl(p);

            if (*p == '(') {
                if (sscanf(++p,"%lg",&d) != 1 || d <= EPSI1)   
                    goto GETTPErr;
  
                p = skip_dbl(p);
                if (*p++ != ')')
                    goto GETTPErr;
                if (sscanf(p,"%lg",&b) != 1 || fabs(b - a) <= EPSI1)
                    goto GETTPErr;
                p = skip_dbl(p);

                if (opt == 1 && b <= a + EPSI1)
                    goto GETTPErr;
                if (opt == 2 && b > 1.0)
                    goto GETTPErr;
                if (opt == 3 && (b >= a - EPSI1 || b < EPSI1))
                    goto GETTPErr;
                if (b < a)
                    d = -d;
                   
                tmp = fabs((b - a) / d);
                if (fabs(tmp - 1.0) <= EPSI1)
                    m = 1;
                else {
                    m = (int)floor(tmp);                   
                    if (m <= 0 || m > 100000) {
                        printf1("Error: exceeded limits of t/tp/x/gss/qo/qt/lat parameter.\n");
                        printf1("or there might be an error in the increment.\n");
                        *err = 1;
                        goto GETTPErr;
                    }
                }
                if (i == 0)  
                    n += m;
                else {
                    for (j = 0; j < m; ++j) {
                        PMTP1[n] = PMTP1[n - 1] + d;
                        n++;
                    }
                }
                a += (double)m * d;
                if (a < b - EPSI1) {
                    if (i)
                        PMTP1[n] = b;
                    n++;
                }
                a1 = b;
            }
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETTPErr;

        if (i == 0) {
            if (!(PMTP1 = (double *)calloc(n,sizeof(double)))) { 
                *err = -2;
                goto GETTPErr;
            }
            memrq(n,sizeof(double));
            PMNTP1 = n;
        }
    }
    if (n != PMNTP1)
        gerr_exit(61);

    *err = 0;

    /***
    printf1("PMNTP1=%d : ",PMNTP1);
    for (i = 0; i < PMNTP1; ++i)
    printf1("%lg ",PMTP1[i]);
    printf1("\n");
    ***/
    return(p);

GETTPErr:
    free_tp1();
    return(p);
}

/* -###-------------------------------------------------------------------- */
/*  free_tp1()   Free memory which was used for time points.                */

void free_tp1(void)         
{
    if (PMNTP1 > 0) {
        free((char *)PMTP1);
        memrq(-PMNTP1,sizeof(double));
        PMNTP1 = 0;
    }
    /***  PMQOFlg = PMQTFlg = 0; **/
}

/*--------------------------------------------------------------------------*/
/*  get_cn(tp,err)  Get cn=.... sequence of integer numbers.                */
/*                  These numbers are put into PMCN[] which is allocated    */
/*                  here. The number of elements is recorded in PMNCN.      */
/*                  Return pointer to next char in string.                  */
/*                                                                          */
/*  err =  0 if successful.                                                 */
/*        -1 if syntax error                                                */
/*        -2 insuff memory                                                  */

char *get_cn(char *tp,int *err)
{
    register int i,n;
    register char *p;
    int m;

    free_cn();
    *err = -1;
    p = tp;
    if (!*p)
        goto GETCNErr;

    for (i = 0; i < 2; ++i) {

        p = tp;
        n = 0;

        while (*p) {
            if (sscanf(p,"%d",&m) != 1) {
                if (n)
                    p--;
                break;
            }
            if (i)
                PMCN[n] = m;
            n++;

            p = skip_int(p);
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETCNErr;

        if (i == 0) {
            if (!(PMCN = (int *)calloc(n,sizeof(int)))) { 
                *err = -2;
                return(p);
            }
            memrq(n,sizeof(int));
            PMNCN = n;
        }
    }
    if (n != PMNCN)
        gerr_exit(61);

    *err = 0;
    return(p);

GETCNErr:
    free_cn();
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_cn()    Free memory which was used for PMCN[]                      */

void free_cn(void)         
{
    if (PMNCN > 0) {
        free((char *)PMCN);
        memrq(-PMNCN,sizeof(int));
        PMNCN = 0;
    }
}

/*--------------------------------------------------------------------------*/
/*  get_flags(tp,err,opt)                                                   */
/*                                                                          */
/*  If opt = 1: dm=i1,i2,...    put into PMDM[]                             */
/*  If opt = 2: sm=i1,i2,...    put into PMSM[]                             */
/*  If opt = 3: tst=i1,i2,...   put into PMTST[]                            */
/*  If opt = 4: sk =i1,i2,...   put into PMSK[]                             */
/*  If opt = 5: rel=i1,i2,...   put into PMREL[]                            */
/*                                                                          */
/*  Range is 1,...,10                                                       */
/*  Return pointer to next char in string. Set err = 0 if successful,       */
/*  otherwise -1 with local err message.                                    */
/*                                                                          */
/*  err =  0  if successful                                                 */
/*        -1  if syntax error                                               */

char *get_flags(char *tp,int *err,int opt)
{
    register char n;
    register char *p;
    int m,nc;

    *err = -1;
    p = tp;
    if (!*p)
        goto GETFLGErr;

    p = tp;
    nc = n = 0;

    while (*p) {
        if (sscanf(p,"%d",&m) != 1) {
            if (n)
                p--;
            break;
        }
        if (opt != 4 && (m < 1 || m > 10))
            goto GETFLGErr;
        if (opt == 1)
            PMDM[m] = 1;
        else if (opt == 2)
            PMSM[m] = 1;
        else if (opt == 3)
            PMTST[m] = 1;
        else if (opt == 5)
            PMREL[m] = 1;
        else if (opt == 4 && nc < 10)
            PMSK[++nc] = m;

        n++;

        p = skip_int(p);
        if (*p != ',')
            break;
        p++;
    }
    if (n == 0)
        goto GETFLGErr;

    *err = 0;

GETFLGErr:
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  pmps_alloc(n,m)     If n >= 0 allocate PMPS[n] for m elements;          */
/*                      otherwise free whole structure.                     */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int pmps_alloc(int n,int m)
{
    register int i;

    if (n >= 0) {
        if (!(PMPS[n] = (short *)calloc(m,sizeof(short))))   
            return(-2);
        memrq(m,sizeof(short));
        PMPSN[n] = m;        
        PMNPS++;
        return(0);
    }
    for (i = 0; i < PMPSMax; ++i) {
        m = PMPSN[i];
        if (m > 0) {
            free((char *)PMPS[i]);
            memrq(-m,sizeof(short));
            PMPSN[i] = 0;
        }
    }
    PMNPS = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  get_pattern(tp,&err)                                                    */
/*                                                                          */
/*  Get test pattern: tp=[x1,x2,...],[y1,y2,...],...                        */
/*  Max number is PMPSMax. Put into PMPS[][], length into PMPSN[].          */
/*  xi... may be >= 0 (valid state) or                                      */
/*  ? = -1                                                                  */
/*  + = -2                                                                  */
/*  - = -3                                                                  */
/*  * = -4                                                                  */
/*                                                                          */
/*  err = 0 if OK, -1 if syntax error, -2 if insufficient memory.           */
/*  Return pointer to next char after end of parameter string.              */

char *get_pattern(char *tp,int *err)
{
    register int i;
    int n,m,k;
    register char *p,*q;

    *err = -1;
    n = -1;
    p = tp;
    while (*p == '[') {
        if (++n >= PMPSMax)  
            return(p);
             
        q = ++p;
        m = 1;
        while (*q && *q != ']') {
            if (*q++ == ',')
                m++;
        }
        if (pmps_alloc(n,m)) {          /* insufficient memory */
            *err = -2;
            return(p);
        }
        PMPSN[n] = m;

        i = 0;
        while (*p && *p != ']') {
            if (*p == '?')
                k = -1;
            else if (*p == '+')
                k = -2;
            else if (*p == '-')
                k = -3;
            else if (*p == '*')
                k = -4;
            else if (sscanf(p,"%d",&k) == 1 && k >= 0)
                p = skip_int(p) - 1;
            else
                return(p);

            PMPS[n][i++] = (short)k;

            if (*++p == ']')
                break;
            if (*p++ != ',')
                return(p);
        }   
        if (i != m)
            return(p);

        if (*p++ != ']')
            return(p);

        if (*p != ',' || *(p + 1) != '[')
            break;
        p++;
    }
    *err = 0;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_fn(p,fname,&err)    get filename at p, put into fname.              */
/*  Return pointer to next character.                                       */
/*                                                                          */
/*  err =  0 if OK, -1 if syntax error.                                     */  

char *get_fn(char *p,char *fname,int *err)
{
    register char c,*q;

    *err = -1;
    if (*p++ != '=')
        return(p);

    q = skip_com(p);
    if (q == p || (int)(q - p) > FNMaxLen) 
        return(p);
    c = *q;
    *q = '\0';  
    strcpy(fname,p);
                
    *err = 0;
    *q = c;
    return(q);
}

/* ------------------------------------------------------------------------ */
/*  get_var(p,err)      p points to a variable list. get variables and      */
/*                      put into PMVIdx. Count number in PMNV.              */
/*                      return pointer to next character.                   */
/*                                                                          */
/*                      if varlist = ...,(...),...,(...),...                */
/*                      put variables in bracketes into separate array      */
/*                      PMZIdx[], number is PMNZ.                           */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var(char *p,int *err)
{
    register int i,j;     
    int n,w,na,nb,ia,ib;

    *err = 1;
    p = get_nvia(p,&n,1,&nb);
    if (n <= 0)  
        return(p);

    na = n - nb;

    if (na > 0) {
        if (pm_valloc(na)) 
            return(p);
    }
    if (nb > 0) {
        if (pm_zalloc(nb)) 
            return(p);
    }
    ia = ib = w = 0;
    for (i = 0; i < n; ++i) {
        j = VLVIdx[i];
        if (j >= 0)
            PMVIdx[ia++] = j;
        else
            PMZIdx[ib++] = -j - 1;
        if (VTyp[j] == 5)
            w = 1;
    }
    if (w) 
        p_warn(-3,1);
    
    *err = 0;
    alloc_vl(0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var1(p,err)     p points to a variable list. get variables and      */
/*                      put into PM1VIdx. Count number in PM1NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var1(char *p,int *err)
{
    register int i;       
    int n,w,nb;

    *err = 1;
    p = get_nvia(p,&n,1,&nb);

    if (n <= 0)  
        return(p);
    if (nb) {
        p_err(-42,1);
        return(p);
    }
    if (pm_v1alloc(n)) 
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        PM1VIdx[i] = VLVIdx[i];
        if (VTyp[VLVIdx[i]] == 5)
            w = 1;
    }
    if (w) 
        p_warn(-3,1);

    *err = 0;
    alloc_vl(0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var2(p,err)     p points to a variable list. get variables and      */
/*                      put into PM2VIdx. Count number in PM2NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var2(char *p,int *err)
{
    register int i;       
    int n,w,nb;

    *err = 1;
    p = get_nvia(p,&n,1,&nb);

    if (n <= 0)  
        return(p);
    if (nb) {
        p_err(-42,1);
        return(p);
    }
    if (pm_v2alloc(n)) 
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        PM2VIdx[i] = VLVIdx[i];
        if (VTyp[VLVIdx[i]] == 5)
            w = 1;
    }
    if (w) 
        p_warn(-3,1);

    *err = 0;
    alloc_vl(0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var3(p,err)     p points to a variable list. get variables and      */
/*                      put into PM3VIdx. Count number in PM3NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var3(char *p,int *err)
{
    register int i;       
    int n,w,nb;

    *err = 1;
    p = get_nvia(p,&n,1,&nb);

    if (n <= 0)  
        return(p);
    if (nb) {
        p_err(-42,1);
        return(p);
    }
    if (pm_v3alloc(n)) 
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        PM3VIdx[i] = VLVIdx[i];
        if (VTyp[VLVIdx[i]] == 5)
            w = 1;
    }
    if (w) 
        p_warn(-3,1);

    *err = 0;
    alloc_vl(0);
    return(p);
}

/* -###=------------------------------------------------------------------- */
/*  get_varx(p,err)     p points to a list of variable or matrix names.     */
/*                      Put indices into PMVIdx. Count numer in PMNV.       */
/*                      If matrix names, set PMNVTyp = 1, otherwise 0.      */
/*                      Return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_varx(char *p,int *err)
{
    register int i,j;     
    int n,w,na,nb,ia,ib,nm,r;
    char *q,mname[VNLMax + 1];
    short idx[MaxMatDef];

    q = p;
    *err = 1;
    p = get_nvia(p,&n,0,&nb);
                             
    if (n > 0) {                /* list of variable names */
        na = n - nb;

        if (na > 0) {
            if (pm_valloc(na)) 
                return(p);
        }
        if (nb > 0) {
            if (pm_zalloc(nb)) 
                return(p);
        }
        ia = ib = w = 0;
        for (i = 0; i < n; ++i) {
            j = VLVIdx[i];
            if (j >= 0)
                PMVIdx[ia++] = j;
            else
                PMZIdx[ib++] = -j - 1;
            if (VTyp[j] == 5)
                w = 1;
        }
        if (w) 
            p_warn(-3,1);
    }  
    else {                      /* list of matrix names */

        r = nm = 0;
        p = q;
        while (1) {
            if ((q = get_mname(p,mname,0)) == NULL) {
                if (nm > 0)
                    p--;
                break;
            }
            if ((i = mat_getidx(mname,0)) < 0) {
                r = 1;
                break;
            }
            if (nm >= MaxMatDef) {
                r = 2;
                break;
            }
            idx[nm++] = i;
            p = q;
            if (*p != ',')
                break;
            p++;
        }
        if (r > 0 || nm == 0) {
            printf1("Error: invalid list of variable or matrix names.\n");
            if (r == 2)  
                printf1("Error: exceeded max number of matrix names.\n");
            return(p);
        }
        if (pm_valloc(nm)) 
            return(p);

        for (i = 0; i < nm; ++i)
            PMVIdx[i] = (short)idx[i];
        PMNVTyp = 1;
    }
    *err = 0;
    alloc_vl(0);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  get_xp(s,err)                                                           */
/*                  Get values xp=x1,x2,... or x=x1[a1,b1],...              */
/*                  Set the number of value in PMNX and allocates the       */
/*                  arrays PMXX[], PMXA[] and PMXB.                         */
/*                                                                          */  
/*                  Set err = 0 if OK, -1 if syntax error, -2 if            */  
/*                  insufficient memory.                                    */
/*                                                                          */
/*                  Return: pointer to next char in string.                 */

char *get_xp(char *s,int *err)
{
    register int i,n;
    register char *p;
    double a,b,x;

    free_xp();
    *err = -1;
    p = s;
    if (!*p)
        goto GETXErr;

    for (i = 0; i < 2; ++i) {

        p = s;
        n = 0;
        while (*p) {
            if (sscanf(p,"%lg[%lg,%lg]",&x,&a,&b) == 3) {
                p = skip_dbl(p);
                p = skip_dbl(p + 1);
                p = skip_dbl(p + 1) + 1;
                if (i) {
                    PMXX[n] = x;
                    PMXA[n] = a;
                    PMXB[n] = b;
                }
            }     
            else if (sscanf(p,"%lg",&x) == 1) {
                p = skip_dbl(p);
                if (i) {
                    PMXX[n] = x;
                    PMXA[n] = x - 1.0;
                    PMXB[n] = x + 1.0;
                }
            }
            else if (n) {
                p--;
                break;
            }
            else
                goto GETXErr;
            n++;
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETXErr;

        if (i == 0) {
            if (!(PMXX = (double *)calloc(n,sizeof(double)))) { 
                *err = -2;
                return(p);
            }
            if (!(PMXA = (double *)calloc(n,sizeof(double)))) { 
                free((char *)PMXX);
                *err = -2;
                return(p);
            }
            if (!(PMXB = (double *)calloc(n,sizeof(double)))) { 
                free((char *)PMXX);
                free((char *)PMXA);
                *err = -2;
                return(p);
            }
            memrq(3 * n,sizeof(double));
            PMNX = n;
        }
    }
    if (n != PMNX)
        gerr_exit(61);

    *err = 0;

    /***
    printf1("PMNX=%d : ",PMNX);
    for (i = 0; i < PMNX; ++i)  
        printf1("%lg %lg %lg\n",PMXX[i],PMXA[i],PMXB[i]);
    ***/
    return(p);

GETXErr:
    free_xp();
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_xp()    Free memory which was used for xp parameter.               */

void free_xp(void)         
{
    if (PMNX > 0) {
        free((char *)PMXX);
        free((char *)PMXA);
        free((char *)PMXB);
        memrq(-3 * PMNX,sizeof(double));
        PMNX = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  get_file(typ,p,err)                                                     */
/*                                                                          */
/*  err =  0 if successful                                                  */
/*        -1 if syntax error                                                */
/*        -2 if insuff memory                                               */
/*         1 if file cannot be opened.                                      */

char *get_file(int typ,char *p,int *err)
{
    register int i;
    int n,aflag,nn,np,np1,l,ii;
    char fname[FNMaxLen + 1];
    register char *q;

    *err = -1;

    nn = 1;
    aflag = 0;
    if (*p == 'a') {
        p++;
        aflag = 1;
    }
    if (typ == 1 || typ == 5) {         /* df or df1, prot or prot1 */
        if (*p == '1') {
            nn = 2;
            p++;
        }
        else if (*p == '2') {
            nn = 3;
            p++;
        }
    }
    if (*p == 'a') {
        p++;
        aflag = 1;
    }
    if (typ == 7) {         /* pfn file */
        if (*p != '(')
            return(p);
        q = skip_blev(p);
        if (*--q != ')')
            return(p);
        *q = '\0';
        n = strlen(++p);
        if (n < 1)
            return(p);

        if (rhstr_alloc(n + 1)) {     
            *err = -2;
            return(p);
        }
        strcpy(PMRHSTR,p);
        *q = ')';
        p = q + 1;
    }
    else if (typ == 12 && *p == '(') {         /* dvar file */

        q = p + 1;
        np = 0;
        while (*q) {    
            if (sscanf(q,"fn=%d",&n) == 1 && n > 0) {
                PMDVARFN = n;
                q = skip_int(q + 3);
            }
            else if (sscanf(q,"pn=%d",&n) == 1 && n > 0) {
                PMDVARPN = n;
                q = skip_int(q + 3);
            }
            else if (sscanf(q,"p=%d",&n) == 1 && n > 0) {
                PMDVARP = n;
                q = skip_int(q + 2);
            }
            else if (!strncmp(q,"p(",2) || !strncmp(q,"pn(",3)) {
                np++;            
                if (!strncmp(q,"p(",2))  
                    q += 2;
                else                      
                    q += 3;
                l = get_vnlen(q);
                if (l < 1)
                    return(p);
                q += l;
                if (*q++ != ')' || *q++ != '=')
                    return(p);
                if (sscanf(q,"%d",&n) != 1 || n < 1)   
                    return(p);
                q = skip_int(q);
            }
            else  
                return(p);
 
            if (*q != ',')
                break;
            q++;
        }
        if (*q++ != ')')
            return(p);

        if (PMDVARP > 0 && PMDVARPN > 0)
            return(p);

        if (np > 0 && (PMDVARP > 0 || PMDVARPN > 0))
            return(p);

        if (np > 0 && PMDVARP == 0 && PMDVARPN == 0) {
/* ### */
            PMDVARVN = np;

            if (!(PMDVARVNP = (short *)calloc(np,sizeof(int)))) {
                *err = -2;    
                return(p);
            }
            memrq(np,sizeof(int));
            PMDVARVNPA = np;

            if (!(PMDVARVNL = (short *)calloc(np,sizeof(int)))) {
                *err = -2;    
                return(p);
            }
            memrq(np,sizeof(int));
            PMDVARVNLA = np;

            if (!(PMDVARVName = (char **)calloc(np,sizeof(char *)))) {
                *err = -2;    
                return(p);
            }
            memrq(np,sizeof(char *));
            PMDVARVNameA = np;

            ii = 0; 
            q = p + 1;

            while (*q) {    
                if (sscanf(q,"fn=%d",&n) == 1 && n > 0) {
                    q = skip_int(q + 3);
                }
                else if (sscanf(q,"pn=%d",&n) == 1 && n > 0) {
                    q = skip_int(q + 3);
                }
                else if (sscanf(q,"p=%d",&n) == 1 && n > 0) {
                    q = skip_int(q + 2);
                }
                else if (!strncmp(q,"p(",2) || !strncmp(q,"pn(",3)) {

                    if (!strncmp(q,"p(",2)) {
                        q += 2;
                        np1 = 1;
                    }
                    else {
                        q += 3;
                        np1 = 2;
                    }
                    l = get_vnlen(q);
                    if (l < 1)
                        return(p);

                    if (!(PMDVARVName[ii] = (char *)calloc(l + 1,sizeof(char)))) {
                        *err = -2;    
                        return(p);
                    }
                    memrq(l + 1,sizeof(char));
                    PMDVARVNP[ii] = np1;
                    strncpy(PMDVARVName[ii],q,l);
                    q += l;
                    if (*q++ != ')' || *q++ != '=')
                        return(p);
                    if (sscanf(q,"%d",&n) != 1 || n < 1)   
                        return(p);
                    q = skip_int(q);
                    PMDVARVNL[ii] = n;
                    ii++;
                }
                else  
                    return(p);
 
                if (*q != ',')
                    break;
                q++;
            }
            if (*q++ != ')' || PMDVARVN != ii)
                return(p);
        }
        p = q;
    }
    else if (typ == 13 && *p == '(') {         /* arcd file */
        p++;
        while (*p) {
            if (!strncmp(p,"zoo=",4)) {
                p += 4;
                q = PMARCFZOOF;
                for (i = 0; i < FNMaxLen; ++i) {
                    if (!*p || *p == ',' || *p == ')')
                        break;
                    *q++ = *p++;
                }
                *q = '\0';
                PMARCFZOO = 1;
            }
            else if (!strncmp(p,"vdf=",4)) {
                p += 4;
                q = PMARCFVDFF;
                for (i = 0; i < FNMaxLen; ++i) {
                    if (!*p || *p == ',' || *p == ')')
                        break;
                    *q++ = *p++;
                }
                *q = '\0';
                PMARCFVDF = 1;
            }
            if (*p != ',')
                break;
            p++;
        }
        if (*p++ != ')')
            return(p);
    }
    p = get_fn(p,fname,err);
    if (*err)  
        return(p);
         
    *err = 0;
    switch (typ) {
        case 1: if (PMF1Def) {
                    fclose(PMF1d);
                    PMF1Def = 0;
                }
                if ((aflag == 0 && !(PMF1d = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMF1d = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMF1Def = nn;
                    strcpy(PMF1dName,fname);
                }
                break;

        case 2: if (PMCovFDef) {
                    fclose(PMCovFd);
                    PMCovFDef = 0;
                }
                if ((aflag == 0 && !(PMCovFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMCovFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMCovFDef = 1;
                    strcpy(PMCovFName,fname);
                }
                break;

        case 3: if (PMResFDef) {
                    fclose(PMResFd);
                    PMResFDef = 0;
                }
                if ((aflag == 0 && !(PMResFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMResFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMResFDef = 1;
                    strcpy(PMResFName,fname);
                }
                break;

        case 4: if (PMPPFDef) {
                    fclose(PMPPFd);
                    PMPPFDef = 0;
                }
                if ((aflag == 0 && !(PMPPFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMPPFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMPPFDef = 1;
                    strcpy(PMPPFName,fname);
                }
                break;

        case 5: if (PMProtFDef) {
                    fclose(PMProtFd);
                    PMProtFDef = 0;
                }
                if ((aflag == 0 && !(PMProtFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMProtFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMProtFDef = nn;
                    strcpy(PMProtFName,fname);
                }
                break;

        case 6: if (PMPCFDef) {
                    fclose(PMPCFd);
                    PMPCFDef = 0;
                }
                if ((aflag == 0 && !(PMPCFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMPCFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMPCFDef = nn;
                    strcpy(PMPCFName,fname);
                }
                break;

        case 7: if (PMPFNFDef) {
                    fclose(PMPFNFd);
                    PMPFNFDef = 0;
                }
                if ((aflag == 0 && !(PMPFNFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMPFNFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMPFNFDef = 2;
                    strcpy(PMPFNFName,fname);
                }
                break;

        case 8: if (PMTabFDef) {
                    fclose(PMTabFd);
                    PMTabFDef = 0;
                }
                if ((aflag == 0 && !(PMTabFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMTabFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMTabFDef = 1;
                    strcpy(PMTabFName,fname);
                }
                break;

        case 9: if (PMTab1FDef) {
                    fclose(PMTab1Fd);
                    PMTab1FDef = 0;
                }
                if ((aflag == 0 && !(PMTab1Fd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMTab1Fd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMTab1FDef = 1;
                    strcpy(PMTab1FName,fname);
                }
                break;

       case 10: if (PMTDAFDef) {
                    fclose(PMTDAFd);
                    PMTDAFDef = 0;
                }
                if ((aflag == 0 && !(PMTDAFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMTDAFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMTDAFDef = 1;
                    strcpy(PMTDAFName,fname);
                }
                break;

       case 11: if (PMSPSSFDef) {
                    fclose(PMSPSSFd);
                    PMSPSSFDef = 0;
                }
                if ((aflag == 0 && !(PMSPSSFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMSPSSFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMSPSSFDef = 1;
                    strcpy(PMSPSSFName,fname);
                }
                break;

       case 12: if (PMDVARFDef) {
                    fclose(PMDVARFd);
                    PMDVARFDef = 0;
                }
                if ((aflag == 0 && !(PMDVARFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMDVARFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMDVARFDef = 1;
                    strcpy(PMDVARFName,fname);
                }
                break;

       case 13: if (PMARCFDef) {
                    fclose(PMARCFd);
                    PMARCFDef = 0;
                }
                if ((aflag == 0 && !(PMARCFd = fopen(fname,OPEN_WR))) ||   
                    (aflag == 1 && !(PMARCFd = fopen(fname,OPEN_AP))))    
                    *err = 1;
                else {
                    PMARCFDef = 1;
                    strcpy(PMARCFName,fname);
                }
                break;

       case 14: if (PMIF1Def) {
                    fclose(PMIF1d);
                    PMIF1Def = 0;
                }
                if (!(PMIF1d = fopen(fname,OPEN_RD)))      
                    *err = 1;
                else {
                    PMIF1Def = 1;
                    strcpy(PMIF1Name,fname);
                }
                break;

       case 15: if (PMIF2Def) {
                    fclose(PMIF2d);
                    PMIF2Def = 0;
                }
                if (!(PMIF2d = fopen(fname,OPEN_RD)))      
                    *err = 1;
                else {
                    PMIF2Def = 1;
                    strcpy(PMIF2Name,fname);
                }
                break;

        default: break;
    }
    if (*err)  
        printf1("Error: can't open: %s\n",fname);
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  get_box(s,err)                                                          */
/*                  Get values box=l1,u1,...,lm,um                          */
/*                  PMBoxL[], PMBoxN                                        */
/*                  PMBoxU[],                                               */
/*                                                                          */  
/*                  Set err = 0 if OK, -1 if syntax error, -2 if            */  
/*                  insufficient memory.                                    */
/*                                                                          */
/*                  Return: pointer to next char in string.                 */

char *get_box(char *s,int *err)
{
    register int i,n;
    register char *p;
    double l,u;   

    free_box();
    *err = -1;
    p = s;
    if (!*p)
        goto GETBOXErr;

    for (i = 0; i < 2; ++i) {

        p = s;
        n = 0;
        while (*p) {
            if (sscanf(p,"%lg,%lg",&l,&u) == 2 && l < u) {
                p = skip_dbl(p);
                p = skip_dbl(p + 1);
                if (i) {
                    PMBoxL[n] = l;
                    PMBoxU[n] = u;
                }
            }
            else if (n) {
                p--;
                break;
            }
            else
                goto GETBOXErr;
            n++;
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETBOXErr;

        if (i == 0) {
            if (!(PMBoxL = (double *)calloc(n,sizeof(double)))) { 
                *err = -2;
                return(p);
            }
            if (!(PMBoxU = (double *)calloc(n,sizeof(double)))) { 
                free((char *)PMBoxL);
                *err = -2;
                return(p);
            }
            memrq(2 * n,sizeof(double));
            PMBoxN = n;         
        }
    }
    if (n != PMBoxN)
        gerr_exit(61);

    *err = 0;
    return(p);

GETBOXErr:
    free_box();
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_box()    Free memory which was used for box parameter.             */

void free_box(void)         
{
    if (PMBoxN > 0) {
        free((char *)PMBoxL);
        free((char *)PMBoxU);
        memrq(-2 * PMBoxN,sizeof(double));
        PMBoxN = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  get_mf(p,err)   Get list of merge files beginning at p, return pointer  */
/*                  to next character. err = 0 if OK, -1 if error.          */

char *get_mf(char *p,int *err)
{
    int n; 
    register char *q;

    *err = -1;
    n = 0;
    while (1) {
        q = p;
        while (*q && *q != '[')
            q++;
        if (*q != '[') {
            if (n)
                *err = 0;
            return(p);
        }
        if (n)
            p++;
        *q++ = '\0';
        if (MFN >= MFMAX) {
            printf1("Error: exceeded max number of merge files.\n");
            return(p);
        }
        if (get_fname(p,MFFNAME[MFN]) == 0)
            return(p);
        if (sscanf(q,"%d,%d,%d,%d]",&MFI1[MFN],&MFI2[MFN],&MFJ1[MFN],&MFJ2[MFN]) != 4)
            return(p);
        while (*q && *q != ']')
            q++;
        if (*q++ != ']')
            return(p);

        if (!(MFFD[MFN] = fopen(MFFNAME[MFN],OPEN_RD))) {   
            printf1("Error: can't open: %s\n",MFFNAME[MFN]);
            return(p);
        }
        MFN++;
        p = q;
        if (*p != ',') {
            if (*p != ')')
                return(p);
            break;
        }
        n++;
    }
    *err = 0;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_mp(p,l,def,name,err)    get matrix name for                         */
/*                              mpcov                                       */
/*                              mppar                                       */
/*                              mplog                                       */
/*                              mpgrad(v=...)   PM2NV                       */
/*                              mpres(v=...)    PM3NV                       */
/*                                                                          */
/*  Return pointer to next character. err = 0 if OK, +/-1 if error.         */

char *get_mp(char *p,int l,int *mdef,char *name,int *err)
{
    register int i;
    register char *q,*r;

    *err = -1;
    r = p + l;

    if (*r == '(') {
        if (!strncmp(p,"mpgrad(v=",9))
            r = get_var2(r + 3,err);
        else if (!strncmp(p,"mpres(v=",8))  
            r = get_var3(r + 3,err);
        else
            return(p);
  
        if (*err)               
            return(p);                     
        *err = -1;
        if (*r++ != ')')   
            return(p);
    }
    if (*r++ != '=')
        return(p);

    q = name;
    for (i = 0; i < VNLMax; ++i) {
        if (!*r || *r == ',' || *r == ')')
            break;
        *q++ = *r++;
    }
    *q = '\0';
    if (mat_ncheck(name,0) || (*r != ',' && *r != ')')) { /* no valid matrix name */
        printf1("Error: %s\nNo valid matrix name.\n",p);
        *err = 1;
        return(p);
    }
    *mdef = 1;     
    *err = 0;
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  get_mod(p,err)      get string mod=..., into PMModS. Max number         */
/*                      of strings is MaxLL. PMModN counts strings.         */
/*                                                                          */
/*  Return pointer to next character. err = 0 if OK, -2 if error.           */

char *get_mod(char *p,int *err)
{
    int n;               
    register char *q;

    *err = -1;

    if (PMModN >= MaxLL) {
        printf1("Error: exceeded max number of mod strings.\n");
        return(p);
    }
    q = p;
    n = 0;
    while (*q && *q != ',') {
        n++;
        q++;
    }
    if (!(PMModS[PMModN] = (char *)calloc(n + 1,sizeof(char)))) {
        *err = -2;    
        return(p);
    }
    memrq(n + 1,sizeof(char));
    PMModSA[PMModN] = n;    

    q = PMModS[PMModN];
    while (*p && *p != ',')  
        *q++ = *p++;
    *q = '\0';
    PMModN++;
    *err = 0;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_mod()  free PMModS                                                 */

void free_mod(void)              
{
    register int i;   

    for (i = 0; i < PMModN; ++i) {
        if (PMModSA[i] > 0) {
            free(PMModS[i]);
            memrq(-PMModSA[i] - 1,sizeof(char));
            PMModSA[i] = 0;
        }
    }
    PMModN = 0;
}

/* ------------------------------------------------------------------------ */
/*  get_fmt(k,p,&err)   get print formats                                   */

char *get_fmt(int k,char *p,int *err)
{
    register int i;
    register char *q;
    int n,m;

    *err = -1;
    n = 0;
    q = p;
    while (*q) {
        if (sscanf(q,"%d",&m) == 1) {
            q = skip_int(q);
            if (*q == '.') {
                if (sscanf(++q,"%d",&m) == 1)  
                    q = skip_int(q);
            }
            n++;
        }
        if (!*q || *q != ',')
            break;
        q++;
    }
    if (n == 0)
        return(p);

    if (get_fmt_alloc(k,n))
        return(p);

    n = 0;
    while (*p) {
        if (sscanf(p,"%d",&m) == 1) {
            p = skip_int(p);
            PMXFmt1[k][n] = m;
            PMXFmt2[k][n] = 0;
            if (*p == '.') {
                if (sscanf(++p,"%d",&m) == 1) {
                    PMXFmt2[k][n] = m;
                    p = skip_int(p);
                }
            }
            n++;
        }
        else {
            if (n > 0)
                p--;
            break;
        }
        if (!*p || *p != ',')
            break;
        p++;
    }
    for (i = 0; i < n; ++i)  
        makefmt(&PMXFmt1[k][i],&PMXFmt2[k][i],PMXFmtS[k] + i * PMXFmtLen,0,SEPC,0);

    PMXFmtN[k] = n;
    *err = 0;
    return(p);
}


/* ------------------------------------------------------------------------ */
/*  get_fmt_alloc(k,n)  allocate mem for n format strings.                  */
/*                      if n = 0 free previously allocated memory.          */
/*                      return 0 if OK, -1 if error.                        */

int get_fmt_alloc(int k,int n) 
{
    if (PMXFmt1A[k] > 0) {
        free((char *)PMXFmt1[k]);
        memrq(-PMXFmt1A[k],sizeof(int));
        PMXFmt1A[k] = 0;
    }
    if (PMXFmt2A[k] > 0) {
        free((char *)PMXFmt2[k]);
        memrq(-PMXFmt2A[k],sizeof(int));
        PMXFmt2A[k] = 0;
    }
    if (PMXFmtSA[k] > 0) {
        free((char *)PMXFmtS[k]);
        memrq(-PMXFmtSA[k],sizeof(char));
        PMXFmtSA[k] = 0;
    }
    PMXFmtN[k] = 0;

    if (n > 0) {
        if (!(PMXFmt1[k] = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        PMXFmt1A[k] = n;

        if (!(PMXFmt2[k] = (int *)calloc(n,sizeof(int)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n,sizeof(int));
        PMXFmt2A[k] = n;

        if (!(PMXFmtS[k] = (char *)calloc(n * PMXFmtLen,sizeof(char)))) {
            p_err(-2,1);
            return(-1);
        }
        memrq(n * PMXFmtLen,sizeof(char));
        PMXFmtSA[k] = n * PMXFmtLen;
        PMXFmtN[k] = n;
    }
    return(0);
}




