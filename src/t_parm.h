/* t_parm.h */

#ifndef _TPARM_H
#define _TPARM_H

/*  functions in t_parm.c */

void p_clean(void);
int parm(char *p,int opt,int rs);
int pm_valloc(int n);
void free_tp(void);        
char *get_file(int typ,char *p,int *err);
char *get_var(char *p,int *err);

extern char *PMFNam[];          /* fn=fn1,fn2,...                           */
extern int PMFN;

extern FILE *PMFd;              /* right-hand side file                     */
extern char PMFdName[];         /* name of right-hand side file             */
extern int PMFDef;              /* set if file opened                       */

extern FILE *PMF1d;             /* df= file (always output file)            */
extern char PMF1dName[];   
extern int PMF1Def;                

extern FILE *PMPCFd;            /* pcf= file (always output file)           */
extern char PMPCFName[];   
extern int PMPCFDef;                

extern FILE *PMF2d;             /* if= file (always input file)             */
extern char PMF2dName[];   
extern int PMF2Def;                
extern int PMIFTyp;             /* if(string)=...                           */

extern FILE *PMIF1d;            /* if1= file (always input file)            */
extern char PMIF1Name[];   
extern int PMIF1Def;                

extern FILE *PMIF2d;            /* if2= file (always input file)            */
extern char PMIF2Name[];   
extern int PMIF2Def;                

extern FILE *PMCovFd;           /* pcov= or pcova=                          */
extern char PMCovFName[];       /* name of file                             */
extern int PMCovFDef;           /* set if file opened                       */
extern int PMCovWFlg;           /* set if file written                      */

extern FILE *PMResFd;           /* pres= or presa=                          */
extern char PMResFName[];       /* name of file                             */
extern int PMResFDef;           /* set if file opened                       */

extern FILE *PMPPFd;            /* pp= or ppa=                              */
extern char PMPPFName[];        /* name of file                             */
extern int PMPPFDef;            /* set if file opened                       */
extern int PMPPWFlg;            /* set if file written                      */

extern FILE *PMTabFd;           /* ptab= or ptaba=                          */
extern char PMTabFName[];       /* name of file                             */
extern int PMTabFDef;           /* set if file opened                       */

extern FILE *PMTab1Fd;          /* ptab1= or ptab1a=                        */
extern char PMTab1FName[];      /* name of file                             */
extern int PMTab1FDef;          /* set if file opened                       */

extern FILE *PMProtFd;          /* prot=                                    */
extern char PMProtFName[];      /* name of file                             */
extern int PMProtFDef;          /* set if file opened                       */

extern FILE *PMPFNFd;           /* pfn()= or pfna()=                        */
extern char PMPFNFName[];       /* name of file                             */
extern int PMPFNFDef;           /* set if file opened                       */

extern FILE *PMTDAFd;           /* dtda=                                    */
extern char PMTDAFName[];       /* name of file                             */
extern int PMTDAFDef;           /* set if file opened                       */

extern FILE *PMSPSSFd;          /* dspss=                                   */
extern char PMSPSSFName[];      /* name of file                             */
extern int PMSPSSFDef;          /* set if file opened                       */

extern FILE *PMDVARFd;          /* dvar=                                    */
extern char PMDVARFName[];      /* name of file                             */
extern int PMDVARFDef;          /* set if file opened                       */
extern int PMDVARFN;            /* logical file number                      */
extern int PMDVARP;             /* p option                                 */
extern int PMDVARPN;            /* pn option                                */
extern int PMDVARVN;            /* number of variables                      */
extern short *PMDVARVNP;        /* var option: 0 (p), 1 (pn)                */
extern short *PMDVARVNL;        /* size of substring                        */
extern char **PMDVARVName;      /* name of variable                         */

extern FILE *PMARCFd;           /* arcd=                                    */
extern char PMARCFName[];       /* name of file                             */
extern int PMARCFDef;           /* set if file opened                       */
extern int PMARCFZOO;          /* zoo = ...                                */
extern char PMARCFZOOF[];   
extern int PMARCFVDF;           /* vdf = ...                                */
extern char PMARCFVDFF[];   

extern int MFN;                 /* number of merge files, mf=...            */
extern int MFI1[];
extern int MFI2[];
extern int MFJ1[];
extern int MFJ2[];
extern FILE *MFFD[];
extern char MFFNAME[][FNMaxLen+1];

extern int PMMPCovDef;          /* set for mpcov=...                        */
extern char PMMPCovName[];      /* corresponding matrix name                */
extern int PMMPParDef;          /* set for mppar=...                        */
extern char PMMPParName[];      /* corresponding matrix name                */
extern int PMMPLogDef;          /* set for mplogl=...                       */
extern char PMMPLogName[];      /* corresponding matrix name                */
extern int PMMPGradDef;         /* set for mpgrad=...                       */
extern char PMMPGradName[];     /* corresponding matrix name                */
extern int PMMPResDef;          /* set for mpres=...                        */
extern char PMMPResName[];      /* corresponding matrix name                */

extern int PMRHSFlg;            /* if right-hand side given                 */
extern int PMRHSI;              /* right-hand side integer argument         */
extern double PMRHSA;           /* right-hand side: first double            */  
extern double PMRHSB;           /* right-hand side: second double           */  
extern double PMRHSD;           /* right-hand side: increment               */  

extern int PMRHSN;              /* number of right-hand side values         */
extern double *PMRHSX;          /* array of right-hand side values          */
extern int PMRHSXA;             /* allocated elements in PMRHSX[]           */

extern char *PMRHSTR;           /* right-hand side string                   */
extern int PMRHSTRA;            /* allocated                                */

extern int PMOPT;               /* opt=...                                  */
extern int PMDOPT;              /* dopt=                                    */
extern int PMMETH;              /* meth=                                    */
extern int PMPRNO;              /* prn=..., print option                    */
extern int PMPLOT;              /* plot= plot option                        */
extern int PMCT;                /* ct=...                                   */

extern char SEPC;               /* Seperation charactr for print formats    */
extern char XSEPC;             
extern int PMFmt1;              /* print format: fmt=...                    */
extern int PMFmt2;
extern int PMFmtF;
extern char PMFmtS[];        

extern int PMTFmt1;             /* print format: tfmt=...                   */
extern int PMTFmt2;
extern int PMTFmtF;
extern char PMTFmtS[];        

extern int PMMFmt1;             /* print format: mfmt=...                   */
extern int PMMFmt2;
extern int PMMFmtF;
extern char PMMFmtS[];        

extern int PMNFmt;              /* print format: nfmt=...                   */
extern int PMNFmtF;
extern char PMNFmtS[];        

extern int PMPFmt1;             /* print format: pfmt=...                   */
extern int PMPFmt2;
extern int PMPFmtF;
extern char PMPFmtS[];        

extern int PMSDFmt1;            /* print format: sdfmt=                     */
extern int PMSDFmt2;
extern int PMSDFmtF;
extern char PMSDFmtS[];        

extern int PMXFmtN[];           /* fmt0=..., fmt1=..., fmt2=...             */
extern int PMXFmtLen;
extern int *PMXFmt1[];
extern int *PMXFmt2[];
extern char *PMXFmtS[];

extern int PMAttr;              /* attr=...                                 */
extern int PMMaxCat;            /* maxcat=..., def NOC                      */
extern int PMMaxCatFlg;         /* set if maxcat used                       */
extern int PMNOC;               /* noc=...                                  */
extern int PMNOCFlg;

extern int PMN;                 /* n=                                       */
extern int PMGLEN;              /* glen=                                    */
extern int PMLEN;               /* len=... (sequence length)                */
extern int PMSN;                /* sn=... (sequence number)                 */
extern int PMSN1; 
extern int PMMSG;               /* msg=...                                  */
extern int PMNS;                /* ns=...                                   */
extern int PMNDIM;              /* ndim=...                                 */
extern int PMNC;                /* nc=...                                   */
extern int PMMin;               /* min=...                                  */
extern int PMMax;               /* max=...                                  */
extern int PMWF;                /* wf=                                      */  
extern int PMPCheck;            /* pcheck=                                  */  

extern int PMNV;                /* number of variables in PMVIdx[]          */
extern int PMNVTyp;             /* 0 if variable names, 1 if matrix names   */
extern short *PMVIdx;           /* indices of variables                     */
extern int PMNZ;                /* number of variables in PMZIdx[]          */
extern short *PMZIdx;           /* indices of variables: v= (...),(...)     */
extern int PMNVLEN;             /* max length of var names in varlist       */

extern int PM1NV;               /* number of variables in PM1VIdx[]          */
extern short *PM1VIdx;          /* indices of variables                     */
extern int PM1NVLEN;            /* max length of var names in varlist       */

extern int PM2NV;               /* number of variables in PM2VIdx[]          */
extern short *PM2VIdx;          /* indices of variables                     */
extern int PM2NVLEN;            /* max length of var names in varlist       */

extern int PM3NV;               /* number of variables in PM3VIdx[]          */
extern short *PM3VIdx;          /* indices of variables                     */
extern int PM3NVLEN;            /* max length of var names in varlist       */

extern int PMDBlockV;           /* variable defined with dblock=            */
extern int PMID;                /* variable defined with id=...             */
extern int PMORG;               /* variable defined with org=...            */
extern int PMDES;               /* variable defined with des=...            */
extern int PMTS;                /* variable defined with ts=...             */
extern int PMTF;                /* variable defined with tf=...             */

extern int PMYL;                /* variable defined with yl=                */
extern int PMYH;                /* variable defined with yh=                */
extern int PMCEN;               /* variable defined with cen=               */
extern int PMTRUNC;             /* variable defined with trunc=             */
extern int PMSCAL;              /* variable defined with scale=             */

extern int PMNTP;               /* number of time points in PMTP[]          */
extern double *PMTP;            /* time points                              */
extern int PMNTP1;              /* number of time points in PMTP1[]         */
extern double *PMTP1;           /* time points                              */
extern int PMQOFlg;             /* set for qo=                              */
extern int PMQTFlg;             /* set for qt=                              */

extern int PMNX;                /* number of xp= parameter                  */
extern double *PMXX;                   
extern double *PMXA;                   
extern double *PMXB;                   

extern int PMBoxN;              /* number of box= parameter                 */
extern double *PMBoxL;                   
extern double *PMBoxU;                   

extern int PMNCN;               /* number of elements in PMCN[]             */
extern int *PMCN;               /* cn=...                                   */
extern int PMDM[];              /* dm=... flags                             */
extern int PMSM[];              /* sm=... flags                             */
extern int PMTST[];             /* tst=... flags                            */
extern int PMREL[];             /* rel=... flags                            */
extern int PMSK[];              /* sk=... flags                             */
extern int PMR;                 /* r=...                                    */
extern int PMS;                 /* s=...                                    */
extern int PMSD;                /* sd=                                      */
extern double PMSC;             /* sc=...                                   */
extern double PMIC;             /* ic=...                                   */
extern int PMSCFlg;             /* set if sc=... used                       */
extern int PMICFlg;             /* set if ic=... used                       */
extern int PMNLEV;              /* nlev=                                    */
extern int PMLEVEL;             /* level=                                   */
extern int PMPROJ;              /* proj=                                    */
extern int PMViewFlg;           /* set if view is specified                 */
extern double PMViewLon;        /* view=                                    */
extern double PMViewLat;
extern double PMRHem;           /* rhem=                                    */  
extern int PMXOrg;              /* psorg=                                   */
extern int PMYOrg;
extern double PMPSRot;          /* psrot=                                   */

extern double PMAlpha;          /* alpha =                                  */
extern double PMBeta;           /* beta  =                                  */
extern double PMGamma;          /* gamma =                                  */
extern double PMIDFA;           /* idf=alpha,beta                           */
extern double PMIDFB;

extern int PMNPS;               /* number of test patterns                  */
extern short *PMPS[];           /* test patterns                            */
extern int PMPSN[];             /* length of test patterns                  */
extern int PMGIdx;              /* index of variable defined with g=...     */
extern double PMSIG;            /* sig=...                                  */
extern double PMOFF;            /* off=...                                  */
extern int PMRXFlg;             /* if rx = a (d) b defined                  */
extern float PMRXA;           
extern float PMRXB;
extern float PMRXD;
extern int PMRYFlg;             /* if ry = a (d) b defined                  */
extern float PMRYA;           
extern float PMRYB;
extern float PMRYD;

extern int PMRRN;               /* rr=n or rr=n,m                           */
extern int PMRRM;

extern double PMD;              /* d=...                                    */

extern double PMRDA;            /* rd=                                      */  
extern double PMRDB;
extern int PMPLFlg;             /* pl=    (plot flag)                       */
extern int PMLT;                /* lt=    (line type)                       */
extern int PMLTFlg; 
extern int PMLT1;               /* lt1=   (line type)                       */
extern double PMLW;             /* lw=    (line width)                      */
extern int PMLWFlg; 
extern double PMLW1;            /* lw1=   (line width)                      */
extern int PMLW1Flg; 
extern double PMFSX;            /* fsx=   (font size)                       */
extern double PMFSY;            /* fsx=   (font size)                       */
extern double PMFS;             /* fs=    (font size)                       */
extern int PMFSFlg; 
extern double PMFSS;            /* fss=   (font size)                       */
extern int PMFSSFlg; 
extern double PMTL;             /* tl=    (tick length)                     */
extern int PMLOG;               /* log=                                     */
extern int PMDIR;               /* dir=                                     */
extern int PMSizeFlg;           /* set for size=...                         */
extern double PMSize;
extern int PMOrder;             /* order=...                                */
extern int PMRegionFlg;         /* set for region=...                       */
extern double PMRegion1;
extern double PMRegion2;
extern int PMXYFlg;             /* set for xy=...                           */
extern double PMX;
extern double PMY;
extern int PMXYZFlg;            /* set for xyz=                             */
extern double PM3X;
extern double PM3Y;
extern double PM3Z;
extern int PMDVECFlg;           /* set for dvec=                            */
extern double PMDVECX;
extern double PMDVECY;
extern double PMDVECZ;
extern int PMGEO;               /* geo=...                                  */
extern int PMHIDE;              /* hide=...                                 */
extern int PMRUFlg;             /* set for ru=...                           */
extern double PMRUA;
extern double PMRUB;
extern int PMRUN;
extern int PMRUM;
extern int PMRVFlg;             /* set for rv=...                           */
extern double PMRVA;
extern double PMRVB;
extern int PMRVN;
extern int PMRVM;
extern int PMNP;                /* np = ...                                 */
extern double PMULX;            /* ulx=...                                  */
extern double PMULY;            /* uly=...                                  */
extern double PMDX;             /* dx=...                                   */
extern double PMDY;             /* dy=...                                   */
extern int PMDXFlg;
extern double PMDXA;            /* dxa=...                                  */
extern int PMDXAFlg;
extern int PMRows;              /* rows=...                                 */
extern int PMCols;              /* cols=...                                 */
extern double PMZMin;           /* zmin=...                                 */
extern double PMZVal;           /* zval=...                                 */
extern int PMZVar;              /* variable defined with zvar=              */
extern double PMTol;            /* tol=...                                  */

extern double PMGS;             /* gs=...                                   */
extern double PMGS1; 
extern int PMGSFlg;    
extern int PMCONT;              /* cont=...                                 */

extern double PMA1;             /* a=a1,a2                                  */
extern double PMA2;      
extern int PMAFlg;    

extern int PMNNFlg;    
extern int PMNN1;               /* nn=n1,n2                                 */
extern int PMNN2;      
extern int PMAGE1;              /* age=a1,a2                                */
extern int PMAGE2;
extern int PMYEAR1;             /* year=y1,y2                               */
extern int PMYEAR2;

extern char *PRC;               /* rc=..... (recode information)            */
extern int PRCAlloc;            /* if PRC allocated                         */

extern char *PRCN;              /* rcn=     (recode information)            */
extern int PRCNAlloc;           /* if PRCN allocated                        */

extern int SVEFlg;              /* set if sel expression available          */
extern short SVECnt;            /* parser stack for sel expression          */
extern int *SVETyp;
extern double *SVEVal;
extern char *SVESTR;            /* string with select expression            */
extern int SVESTRA;             /* if allocated                             */

/* ------------------------------------------------------------------------ */
extern int PMNW;                /* nw=... (number of waves)                 */
extern int PMBlock;             /* block=                                   */  
extern int PMNI;                /* ni=...                                   */  
extern int PMNQ;                /* nq=                                      */  
extern int PMNConS;             /* con=  number of con= expressions         */
extern int NCONSTR;             /* number of constraint expressions         */
/* ------------------------------------------------------------------------ */
extern double PMCFrac;          /* cfrac=...                                */
extern int PMCSF;               /* csf=...                                  */
extern int PMNXA;               /* number of xa(... strings                 */
extern int PMDSVFlg;            /* set if dsv=...                           */
extern char PMDSVName[];        /* name of dsv file                         */
/* ------------------------------------------------------------------------ */
extern int PMPMN;               /* number of panel miss value codes         */
extern double PMPMVal[];        /* pmval= (panel miss value)                */
extern int PMPMin;              /* pmin= (min number of participation)      */
extern int PMRes[];             /* res = r1,r2,...                          */
extern int PMResN;              /* number of r1,r2,... in PMRes             */
/* ------------------------------------------------------------------------ */
extern double PMNDIGIT;         /* ndigit                                   */
extern int PMNDIGFlg;
/* ------------------------------------------------------------------------ */
extern int PMRRFlg;             /* set by rrisk                             */
extern int PMPRN;               /* number of prate strings                  */
/* ------------------------------------------------------------------------ */
extern double PMRERR;           /* rerr=                                    */
extern int PMRERRFlg;
extern double PMAERR;           /* aerr=                                    */
extern int PMAERRFlg;
extern int PMM;                 /* m=                                       */
extern int PMMFlg;              /* set if m=...                             */
extern int PMM1;                /* m1=                                      */
extern int PMM1Flg;             /* set if m1=...                            */
extern int PMM2;                /* m2=                                      */
extern int PMM2Flg;             /* set if m2=...                            */
extern int PMDEG;               /* deg=                                     */
extern double PMKGam;           /* kgam=                                    */
extern int PMKGamFlg;
/* ------------------------------------------------------------------------ */
extern int PMKeep;              /* set by keep=varlist                      */
extern int PMDrop;              /* set by drop=varlist                      */
extern int PMTransp;            /* set by transp                            */
extern int PMArcDic;            /* set by arcdic                            */
extern double PMMSYS;           /* msys=...                                 */
extern int PMSORTFlg;           /* set by sort,                             */
extern int PMAP;                /* ap=...                                   */
extern int DGRPFlg;             /* set by dgrp option                       */
/* ------------------------------------------------------------------------ */
extern int PMGT;                /* gt=                                      */
extern int PMGTT;               /* gt(n)=                                   */
extern int PMRT;                /* rt=                                      */
extern double PMMR;             /* mr=                                      */
extern int PMFTYP5;             /* set if function contains type 5 var      */
/* ------------------------------------------------------------------------ */
extern int PMYWVar;             /* variable specified with yw=...           */
extern int PMWVar;              /* variable specified with w=...            */
extern int PMNBOX;              /* max number of boxes                      */
extern int PMGN;                /* number of points for g_min3()            */
extern int PMGN1;
extern int PMGNK;               /* number of nearest neighbors              */
extern double PMGD;             /* d for random search, II                  */

extern int PMEVSN;              /* ev = [sn,j,k]                            */
extern int PMEV1; 
extern int PMEV2; 
extern int XEFlg;               /* xe=...                                   */
/* ------------------------------------------------------------------------ */
extern double PMICOSTA;         /* icost = ...                              */
extern double PMICOSTB;     
extern int PMICOSTMAT;  
extern int PMSCOSTM;            /* scost = ...                              */
extern int PMSCOSTMAT;  
/* ------------------------------------------------------------------------ */
extern int PMNHP;               /* nhp = ...                                */
extern int PMNMPnt;             /* nmp=...                                  */
extern int PMPTyp;              /* ptyp=                                    */
extern int PMTyp;               /* typ=                                     */
extern double PMEPS;            /* eps=                                     */
extern int PMEPSFlg;
extern char *PMSTR;             /* str=...                                  */
extern int PMSTRA;              /* allocated                                */
extern int PMLINK;              /* link=                                    */
extern int PMMXCYC;             /* mxcyc=...                                */
extern int PMIV;                /* iv=                                      */
extern int PMMIX;               /* mix=...                                  */
extern int PMK;                 /* k= ...                                   */
extern int PML0;                /* l0=...                                   */
/* ------------------------------------------------------------------------ */
extern int PMXLenFlg;
extern int PMYLenFlg;
extern double PMXLen;           /* pxlen=...                                */
extern double PMYLen;           /* pylen=...                                */
extern int PMatNameFlg;         /* set if mdef=...                          */
extern char PMatName[];         /* mame defined with mdef                   */
extern int PMALG;               /* select algorithm                         */
extern int PMSEED;              /* seed=...                                 */  
extern int PMDIM;               /* dim=...                                  */  
extern int PMCG;                /* cg=...                                   */  
extern int PMPERM;              /* perm=...                                 */
/* ------------------------------------------------------------------------ */
extern int PMModN;              /* number of model strings (mod=...)        */
extern char *PMModS[];          /* strings                                  */
extern int PMModSA[];           /* length if allocated                      */
extern int SCRNFlg;             /* set if screen parameter                  */
/* ------------------------------------------------------------------------ */
extern char *PMF1;              /* f1=...  function expression              */
extern char *PMF2;              /* f2=...  function expression              */
extern char *PMF3;              /* f3=...  function expression              */
extern int PMF1A;               /* allocated                                */
extern int PMF2A;               /* allocated                                */
extern int PMF3A;               /* allocated                                */
extern int PMRECFlg;            /* set for rec=...                          */
extern double PMRECXMin;
extern double PMRECYMin;
extern double PMRECXMax;
extern double PMRECYMax;

#endif /* _TPARM_H */

/* end of t_parm.h */


