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
#include "tda_context.h"

/*  functions in t_parm.c */

void p_clean(TDAContext *ctx);
void p_clean1(TDAContext *ctx);
void p_fclose(TDAContext *ctx);
int parm(TDAContext *ctx, char *p,int opt,int rs);
void prn_perr(TDAContext *ctx, char *s);
void pmfn_free(TDAContext *ctx);
int sve_alloc(TDAContext *ctx, int n,char *s);
int rhstr_alloc(TDAContext *ctx, int n);
int pmstr_alloc(TDAContext *ctx, int n);
int pmf1_alloc(TDAContext *ctx, int n);
int pmf2_alloc(TDAContext *ctx, int n);
int pmf3_alloc(TDAContext *ctx, int n);
int rhs_alloc(TDAContext *ctx, int n);
int prc_alloc(TDAContext *ctx, int n);
int prcn_alloc(TDAContext *ctx, int n);
int pm_valloc(TDAContext *ctx, int n);
int pm_zalloc(TDAContext *ctx, int n);
int pm_v1alloc(TDAContext *ctx, int n);
int pm_v2alloc(TDAContext *ctx, int n);
int pm_v3alloc(TDAContext *ctx, int n);
char *get_tp(TDAContext *ctx, char *tp,int *err,int opt);
void free_tp(TDAContext *ctx);
char *get_tp1(TDAContext *ctx, char *tp,int *err,int opt);
void free_tp1(TDAContext *ctx);
char *get_cn(TDAContext *ctx, char *tp,int *err);
void free_cn(TDAContext *ctx);
char *get_flags(TDAContext *ctx, char *tp,int *err,int opt);
int pmps_alloc(TDAContext *ctx, int n,int m);
char *get_pattern(TDAContext *ctx, char *p,int *err);
char *get_fn(TDAContext *ctx, char *p,char *fname,int *err);
char *get_var(TDAContext *ctx, char *p,int *err);
char *get_var1(TDAContext *ctx, char *p,int *err);
char *get_var2(TDAContext *ctx, char *p,int *err);
char *get_var3(TDAContext *ctx, char *p,int *err);
char *get_varx(TDAContext *ctx, char *p,int *err);
char *get_xp(TDAContext *ctx, char *s,int *err);
void free_xp(TDAContext *ctx);
char *get_file(TDAContext *ctx, int typ,char *p,int *err);
char *get_box(TDAContext *ctx, char *s,int *err);
void free_box(TDAContext *ctx);
char *get_mf(TDAContext *ctx, char *p,int *err);
char *get_mp(TDAContext *ctx, char *p,int l,int *mdef,char *name,int *err);
char *get_mod(TDAContext *ctx, char *p,int *err);
void free_mod(TDAContext *ctx);
char *get_fmt(TDAContext *ctx, int k,char *p,int *err);
int get_fmt_alloc(TDAContext *ctx, int k,int n);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

/*  input parameter of commands                                             */

#define PMFNMax 100             /* max strings in fn parameter              */






































                                /* t=, tp=, qo=, qt=, wt=                   */

























/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */
/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------------------------ */
/*  p_clean()   Free all previously allocated memory and close files.       */

void p_clean(TDAContext *ctx)
{
    p_clean1(ctx);
    a_clean(ctx);          /* free all memory in t_alloc */
}

/* ------------------------------------------------------------------------ */
/*  p_clean1()   Free all previously allocated memory and close files.      */

void p_clean1(TDAContext *ctx)
{
    int n;

    p_fclose(ctx);
    pmfn_free(ctx);
    sve_alloc(ctx, 0,NULL);
    rhstr_alloc(ctx, 0);
    pmstr_alloc(ctx, 0);
    pmf1_alloc(ctx, 0);
    pmf2_alloc(ctx, 0);
    pmf3_alloc(ctx, 0);
    rhs_alloc(ctx, 0);
    prc_alloc(ctx, 0);
    prcn_alloc(ctx, 0);
    pm_valloc(ctx, 0);
    pm_zalloc(ctx, 0);
    pm_v1alloc(ctx, 0);
    pm_v2alloc(ctx, 0);
    pm_v3alloc(ctx, 0);
    free_tp(ctx);
    free_tp1(ctx);
    free_xp(ctx);
    free_box(ctx);
    free_cn(ctx);
    pmps_alloc(ctx, -1,0);
    get_func(ctx, NULL,0,&n,0,NULL);
    alloc_vl(ctx, 0);
    rdmn_free(ctx);
    free_mod(ctx);
    get_fmt_alloc(ctx, 0,0);
    get_fmt_alloc(ctx, 1,0);
    get_fmt_alloc(ctx, 2,0);
}

/* ------------------------------------------------------------------------ */
/*  p_fclose()      Close files (if open)                                   */

void p_fclose(TDAContext *ctx)
{
    register int i,l;

    if (ctx->PMFDef) {
        fclose(ctx->PMFd);
        ctx->PMFDef = 0;
    }
    if (ctx->PMF1Def) {
        fclose(ctx->PMF1d);
        ctx->PMF1Def = 0;
    }
    if (ctx->PMF2Def) {
        fclose(ctx->PMF2d);
        ctx->PMIFTyp = ctx->PMF2Def = 0;
    }
    if (ctx->PMCovFDef) {
        fclose(ctx->PMCovFd);
        ctx->PMCovWFlg = ctx->PMCovFDef = 0;
    }
    if (ctx->PMResFDef) {
        fclose(ctx->PMResFd);
        ctx->PMResFDef = 0;
    }
    if (ctx->PMPPFDef) {
        fclose(ctx->PMPPFd);
        ctx->PMPPWFlg = ctx->PMPPFDef = 0;
    }
    if (ctx->PMProtFDef) {
        fclose(ctx->PMProtFd);
        ctx->PMProtFDef = 0;
    }
    if (ctx->PMPFNFDef) {
        fclose(ctx->PMPFNFd);
        ctx->PMPFNFDef = 0;
    }
    if (ctx->PMTabFDef) {
        fclose(ctx->PMTabFd);
        ctx->PMTabFDef = 0;
    }
    if (ctx->PMTab1FDef) {
        fclose(ctx->PMTab1Fd);
        ctx->PMTab1FDef = 0;
    }
    if (ctx->PMTDAFDef) {
        fclose(ctx->PMTDAFd);
        ctx->PMTDAFDef = 0;
    }
    if (ctx->PMSPSSFDef) {
        fclose(ctx->PMSPSSFd);
        ctx->PMSPSSFDef = 0;
    }
    if (ctx->PMPCFDef) {
        fclose(ctx->PMPCFd);
        ctx->PMPCFDef = 0;
    }
    if (ctx->PMDVARFDef) {
        fclose(ctx->PMDVARFd);
        ctx->PMDVARFDef = ctx->PMDVARP = ctx->PMDVARPN = 0;
        if (ctx->PMDVARVN > 0) {
            if (ctx->PMDVARVNameA > 0) {
                for (i = 0; i < ctx->PMDVARVN; ++i) {
                    if (ctx->PMDVARVNP[i]) {
                        l = (int)(strlen(ctx->PMDVARVName[i]) + 1);
                        free(ctx->PMDVARVName[i]);
                        memrq(ctx, -l,sizeof(char));
                    }
                }
                free((char *)ctx->PMDVARVName);
                memrq(ctx, -ctx->PMDVARVNameA,sizeof(char *));
                ctx->PMDVARVNameA = 0;
            }
            if (ctx->PMDVARVNPA > 0) {
                free((char *)ctx->PMDVARVNP);
                memrq(ctx, -ctx->PMDVARVNPA,sizeof(int));
                ctx->PMDVARVNPA = 0;
            }
            if (ctx->PMDVARVNLA > 0) {
                free((char *)ctx->PMDVARVNL);
                memrq(ctx, -ctx->PMDVARVNLA,sizeof(int));
                ctx->PMDVARVNLA = 0;
            }
        }
        ctx->PMDVARVN = 0;
    }
    if (ctx->PMARCFDef) {
        fclose(ctx->PMARCFd);
        ctx->PMARCFDef = ctx->PMARCFZOO = ctx->PMARCFVDF = 0;
    }
    for (i = 0; i < ctx->MFN; ++i)
        fclose(ctx->MFFD[i]);
    ctx->MFN = 0;
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

int parm(TDAContext *ctx, char *p,int opt,int rs)
{
    register int i;
    int n,m,m1,err,r,fnd;
    register char c,*q;
    char vname[VNLMax + 1];
    char *p0;
    double a,b,d;

    p0 = p;
    err = -1;
    p_clean1(ctx);                 /* clean previously allocated arrays */

    ctx->PMIFTyp = 0;                /* if(string)=... */
    ctx->SEPC = ' ';                 /* separation character */
    ctx->XSEPC = ' ';
    ctx->NCONSTR = 0;                /* number of constraint expressions */
    ctx->PMAttr = 0;                 /* attr=... */
    ctx->PMMaxCat = 1000;            /* maxcat= */
    ctx->PMMaxCatFlg = 0;
    ctx->PMRHSFlg = 0;               /* set if right-hand side given */
    ctx->PMRHSI = 0;                 /* right-hand side integer */
    ctx->PMRHSA = 0.0;               /* right-hand side: double */
    ctx->PMRHSB = 0.0;               /* right-hand side: double */
    ctx->PMRHSD = 0.0;               /* right-hand side: increment */
    ctx->PMOPT = 1;                  /* opt =  */
    ctx->PMDOPT = -1;                /* dopt =  */
    ctx->PMMETH = -1;                /* meth= */
    ctx->PMPRNO  = 0;                /* prn =  */
    ctx->PMPLOT = 0;                 /* plot= */
    ctx->PMCT = 0;                   /* ct=... */

    ctx->PMFmt1  = 0;                /* fmt */
    ctx->PMFmt2  = 0;
    ctx->PMFmtF  = 0;

    ctx->PMTFmt1 = 10;               /* tfmt */
    ctx->PMTFmt2 = 4;
    ctx->PMTFmtF = 0;

    ctx->PMMFmt1 = 12;               /* mfmt */
    ctx->PMMFmt2 = 4;
    ctx->PMMFmtF = 0;

    ctx->PMNFmt  = 4;                /* nfmt */
    ctx->PMNFmtF = 0;

    ctx->PMPFmt1 = -19;              /* pfmt */
    ctx->PMPFmt2 =  11;
    ctx->PMPFmtF = 0;

    ctx->PMSDFmt1 = 12;              /* sdfmt */
    ctx->PMSDFmt2 = 4;
    ctx->PMSDFmtF = 0;

    ctx->PMNOC = 1000;               /* noc= */
    ctx->PMNOCFlg = 0;
    ctx->PMGLEN = 0;                 /* glen= */
    ctx->PMLEN = -1;                 /* len= sequence length */
    ctx->PMSN = -1;                  /* sn= sequence number */
    ctx->PMSN1 = -1;
    ctx->PMMSG = -1;                 /* msg= */
    ctx->PMNS = -1;                  /* ns=  */
    ctx->PMNDIM = -1;                /* ndim=  */
    ctx->PMNC =  0;                  /* nc= */
    ctx->PMMin = -1;                 /* min= */
    ctx->PMMax = -1;                 /* max= */
    ctx->PMWF = 0;                   /* wf= */
    ctx->PMPCheck = 1;               /* pcheck= */

    ctx->PMDBlockV = -1;             /* dblock= */
    ctx->PMID  = -1;                 /* id= */
    ctx->PMORG = -1;                 /* org= */
    ctx->PMDES = -1;                 /* des= */
    ctx->PMTS  = -1;                 /* ts= */
    ctx->PMTF  = -1;                 /* tf= */

    ctx->PMYL = -1;                  /* yl= */
    ctx->PMYH = -1;                  /* yh= */
    ctx->PMCEN = -1;                 /* cen= */
    ctx->PMTRUNC = -1;               /* trunc= */
    ctx->PMSCAL  = -1;               /* scale= */

    ctx->PMQOFlg = 0;                /* set for qo= */
    ctx->PMQTFlg = 0;                /* set for qt= */

    ctx->PMN = 0;                    /* n= */
    ctx->PMM = 1;                    /* m= */
    ctx->PMMFlg = 0;                 /* set if m=... */
    ctx->PMM1 = 1;                   /* m1= */
    ctx->PMM1Flg = 0;                /* set if m1=... */
    ctx->PMM2 = 1;                   /* m2= */
    ctx->PMM2Flg = 0;                /* set if m2=... */
    ctx->PMR = 0;                    /* r= */
    ctx->PMS = 0;                    /* s= */
    ctx->PMSD = 0;                   /* sd= */
    ctx->PMSC = 0.0;                 /* sc= */
    ctx->PMIC = 0.0;                 /* ic= */
    ctx->PMSCFlg = ctx->PMICFlg = 0;
    ctx->PMNLEV = 0;                 /* nlev=                                    */
    ctx->PMLEVEL = 0;                /* level=                                   */
    ctx->PMPROJ = 0;                 /* proj= */
    ctx->PMViewFlg = 0;              /* set if view is specified */
    ctx->PMViewLon = 0.0;            /* view= */
    ctx->PMViewLat = 0.0;
    ctx->PMRHem = 90.0;              /* rhem=  */
    ctx->PMXOrg = 100;               /* psorg= */
    ctx->PMYOrg = 100;
    ctx->PMPSRot = 0.0;              /* psrot= */

    ctx->PMAlpha = 0.0;              /* alpha = */
    ctx->PMBeta  = 0.0;              /* beta  = */
    ctx->PMGamma = 0.0;              /* gamma = */
    ctx->PMIDFA = ctx->PMIDFB = -1.0;     /* idf=alpha,beta */
    ctx->PMGIdx = -1;                /* index of variable defined with g= */
    ctx->PMSIG = -1.0;               /* sig= */
    ctx->PMOFF =  0.0;               /* off= */
    ctx->PMRXFlg = 0;                /* if rx = a (d) b defined */
    ctx->PMRXA = 0.0;
    ctx->PMRXB = 0.0;
    ctx->PMRXD = 0.0;
    ctx->PMRYFlg = 0;                /* if rx = a (d) b defined */
    ctx->PMRYA = 0.0;
    ctx->PMRYB = 0.0;
    ctx->PMRYD = 0.0;
    ctx->PMRRN = 0;                  /* rr=n or rr=n,m */
    ctx->PMRRM = 0;
    ctx->PMLOG = 0;                  /* log= */
    ctx->PMDIR = 0;                  /* dir= */
    ctx->PMSizeFlg = 0;              /* set for size=... */
    ctx->PMSize = 0.0;
    ctx->PMOrder = 0;                /* oder=... */
    ctx->PMRegionFlg = 0;            /* set for region=... */
    ctx->PMRegion1 = 0.0;
    ctx->PMRegion2 = 0.0;
    ctx->PMXYFlg = 0;                /* set for xy= */
    ctx->PMX = 0.0;
    ctx->PMY = 0.0;
    ctx->PMD = 0.0;                  /* d= */
    ctx->PMXYZFlg = 0;               /* set for xyz= */
    ctx->PM3X = 0.0;
    ctx->PM3Y = 0.0;
    ctx->PM3Z = 0.0;
    ctx->PMDVECFlg = 0;              /* set for dvec= */
    ctx->PMDVECX = 0.0;
    ctx->PMDVECY = 0.0;
    ctx->PMDVECZ = 0.0;
    ctx->PMGEO = 0;                  /* geo=... */
    ctx->PMHIDE = 0;                 /* hide=... */
    ctx->PMRUFlg = 0;                /* set for ru=... */
    ctx->PMRUA = 0.0;
    ctx->PMRUB = 0.0;
    ctx->PMRUN = 0;
    ctx->PMRUM = 0;
    ctx->PMRVFlg = 0;                /* set for rv=...*/
    ctx->PMRVA = 0.0;
    ctx->PMRVB = 0.0;
    ctx->PMRVN = 0;
    ctx->PMRVM = 0;
    ctx->PMNP = 0;                   /* np = ... */
    ctx->PMULX = 0.0;                /* ulx=... */
    ctx->PMULY = 0.0;                /* uly=... */
    ctx->PMDX = 0.0;                 /* dx=...  */
    ctx->PMDY = 0.0;                 /* dy=...  */
    ctx->PMDXFlg = 0;
    ctx->PMDXA = 0.0;                /* dxa=... */
    ctx->PMDXAFlg = 0;
    ctx->PMRows = 0;                 /* rows=... */
    ctx->PMCols = 0;                 /* cols=... */
    ctx->PMZMin = 0.0;               /* zmin=... */
    ctx->PMZVal = 0.0;               /* zval=... */
    ctx->PMZVar = -1;                /* zvar= */
    ctx->PMTol = 1.e-4;              /* tol=... */

    ctx->PMRDA = 0.0;                /* rd= */
    ctx->PMRDB = 0.0;
    ctx->PMPLFlg = 0;                /* pl= (plot flag) */
    ctx->PMLT = 1;                   /* lt= (line type) */
    ctx->PMLT1 = 1;                  /* lt1= (line type) */
    ctx->PMLW = ctx->LWDef;               /* lw= (line width) */
    ctx->PMLW1 = ctx->LW1Def;             /* lw1= (line width) */

    ctx->PMFSX = 0.0;                /* fsx= (font size) */
    ctx->PMFSY = 0.0;                /* fsx= (font size) */

    ctx->PMFSS = ctx->PMFS = ctx->FSDef;       /* fss= fs=... (font size) */
    ctx->PMTL = 0.0;                 /* tl= (tick length) */
    ctx->PMLTFlg = ctx->PMLWFlg = ctx->PMLW1Flg = ctx->PMFSFlg = ctx->PMFSSFlg = 0;

    ctx->PMGS = 0.0;                 /* gs=  (grey scale value) */
    ctx->PMGS1 = 0.0;
    ctx->PMGSFlg = 0;
    ctx->PMCONT = 0;                 /* cont=... */

    ctx->PMA1 = ctx->PMA2 = 0.0;          /* a= a1,a2 */
    ctx->PMAFlg = 0;
    ctx->SVEFlg = 0;
    ctx->FNFlg = 0;
    ctx->PMNNFlg = ctx->PMNN1 = ctx->PMNN2 = 0;    /* nn=n1,n2 */
    ctx->PMAGE1 = ctx->PMAGE2 = -1;           /* age=a1,a2 */
    ctx->PMYEAR1 = ctx->PMYEAR2 = -1;         /* year=y1,y2 */

    ctx->PMNConS = 0;                /* number of con= expressions */
    ctx->PMNW = 1;                   /* nw= (number of waves) */
    ctx->PMBlock = 0;                /* block= */
    ctx->PMNI = 0;                   /* ni= */
    ctx->PMNQ = 0;                   /* nq= */

    ctx->PMCSF = 0;                  /* csf */
    ctx->PMNXA = 0;                  /* number of xa(... strings */
    ctx->PMDSVFlg = 0;               /* set if dsv= */
    ctx->PMCFrac = 0.5;              /* cfrac= */
    ctx->PMPMN = 0;                  /* pmval= number of codes */
    ctx->PMPMin = 1;                 /* pmin= (min number of participation) */
    ctx->PMResN = 0;                 /* res=... */

    ctx->PMNDIGIT = 15.0;            /* ndigit */
    ctx->PMNDIGFlg = 0;
    ctx->PMRRFlg = 0;                /* set by rrisk */
    ctx->PMPRN = 0;                  /* number of prate strings */
    ctx->PMRERR = 1.e-4;             /* rerr= */
    ctx->PMRERRFlg = 0;
    ctx->PMAERR = 1.e-4;             /* aerr= */
    ctx->PMAERRFlg = 0;
    ctx->PMDEG = 0;                  /* deg= */
    ctx->PMKGam = 1.0;               /* kgam= */
    ctx->PMKGamFlg = 0;
    ctx->PMKeep = 0;                 /* set by keep=varlist */
    ctx->PMDrop = 0;                 /* set by drop=varlist */
    ctx->PMTransp = 0;               /* set by transp */
    ctx->PMArcDic = 0;               /* set by arcdic */
    ctx->PMMSYS = -5.0;              /* msys= */
    ctx->PMSORTFlg = 0;              /* set by sort, */
    ctx->PMAP = 0;                   /* ap= */
    ctx->DGRPFlg = 0;                /* set by dgrp option */
    ctx->PMGT = 0;                   /* gt= */
    ctx->PMGTT = 1;                  /* gtt */
    ctx->PMRT = 0;                   /* rt= */
    ctx->PMMR = 2.0;                 /* mr= */
    ctx->PMFTYP5 = 0;                /* set if function contains type 5 var */
    ctx->PMYWVar = -1;               /* variable specified with yw=... */
    ctx->PMWVar = -1;                /* variable specified with w=... */
    ctx->PMNBOX = -1;                /* max number of boxes */
    ctx->PMGN = 0;                   /* number of points for g_min3() */
    ctx->PMGN1 = 0;
    ctx->PMGNK = 0;                  /* number of nearest neighbors */
    ctx->PMGD = 1.0;                 /* gd=  d for random search, II */
    ctx->PMEVSN = -1;                /* ev = [sn,j,k] */
    ctx->PMEV1 = -1;
    ctx->PMEV2 = -1;
    ctx->XEFlg = 0;                  /* xe=... */
    ctx->PMICOSTA = -1.0;            /* icost = ...  */
    ctx->PMICOSTB = -1.0;
    ctx->PMICOSTMAT = -1;
    ctx->PMSCOSTM = -1;              /* scost = ... */
    ctx->PMSCOSTMAT = -1;
    ctx->PMNHP = 6;                  /* nhp = ... */
    ctx->PMNMPnt = 1;                /* nmp=... */
    ctx->PMPTyp = 1;                 /* ptyp=... */
    ctx->PMTyp = 0;                  /* typ= ... */
    ctx->PMEPS = 1.e-6;              /* eps= */
    ctx->PMEPSFlg = 0;
    ctx->PMLINK = 0;                 /* link= */
    ctx->PMMXCYC = -1;               /* mxcyc=... */
    ctx->PMIV = 0;                   /* iv= */
    ctx->PMMIX = 0;                  /* mix=... */
    ctx->PMK = 0;                    /* k=... */
    ctx->PML0 = 0;                   /* l0=... */
    ctx->PMXLenFlg = 0;
    ctx->PMYLenFlg = 0;
    ctx->PMXLen = 120.0;             /* pxlen=... */
    ctx->PMYLen =  80.0;             /* pylen=... */
    ctx->PMMPCovDef = 0;             /* set for mpcov=... */
    ctx->PMMPParDef = 0;             /* set for mppar=... */
    ctx->PMMPLogDef = 0;             /* set for mplog=... */
    ctx->PMMPGradDef = 0;            /* set for mpgrad=... */
    ctx->PMMPResDef = 0;             /* set for mpres=... */
    ctx->PMatNameFlg = 0;            /* set if mdef=...  */
    ctx->PMALG = 0;                  /* select algorithm */
    ctx->PMSEED = 0;                 /* seed=... */
    ctx->PMDIM = 0;                  /* dim=... */
    ctx->PMCG = 0;                   /* cg=... */
    ctx->PMPERM = 0;                 /* perm=... */
    ctx->SCRNFlg = 0;                /* screen */
    ctx->PMRECFlg = 0;               /* set for rec=... */


    for (i = 1; i <= 10; ++i)
        ctx->PMSK[i] = ctx->PMSM[i] = ctx->PMDM[i] = ctx->PMTST[i] = ctx->PMREL[i] = 0;

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
                        ctx->PMA1 = a;
                        ctx->PMA2 = b;
                        ctx->PMAFlg = 1;
                        p = skip_dbl(ctx, p + 2);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"ab=%lf,%lf",&a,&b) == 2) {
                        ctx->PMX = a;
                        ctx->PMY = b;
                        ctx->PMXYFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"amue=%lg",&a) == 1 && a > 0.0 && a < 1.0) {
                        ctx->AMue = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"aerr=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMAERR = a;
                        ctx->PMAERRFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"ap=%d",&n) == 1) {
                        ctx->PMAP = n;
                        p = skip_int(ctx, p + 3);
                    }
                    /* n >= 0: bfa documents alg=0 as its default; 0 is
                       also PMALG's reset value, so accepting it is
                       identical to omitting the option for every
                       command.  With n > 0 here, the documented
                       bfa(alg=0) fell through to the unknown-parameter
                       path and its raw text was echoed as the error. */
                    else if (sscanf(p,"alg=%d",&n) == 1 && n >= 0) {
                        ctx->PMALG = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"arcdic",6)) {
                        ctx->PMArcDic = 1;
                        p += 6;
                    }
                    else if (!strncmp(p,"arcd",4)) {    /* PMARC file */
                        p = get_file(ctx, 13,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"alpha=%lf",&a) == 1) {
                        ctx->PMAlpha = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (!strncmp(p,"av=",3)) {
                        p = get_var1(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"attr=%d",&n) == 1 && n >= 1) {
                        ctx->PMAttr = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"age=%d,%d",&n,&m) == 2 && n >= 0 && m >= n) {
                        ctx->PMAGE1 = n;
                        ctx->PMAGE2 = m;
                        p = skip_int(ctx, p + 4);
                        p = skip_int(ctx, p + 1);
                    }
                    else
                        fnd = 0;
                    break;

                case 'b':

                    if (!strncmp(p,"box=",4)) {
                        p = get_box(ctx, p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"beta=%lf",&a) == 1) {
                        ctx->PMBeta = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"block=%d",&n) == 1 && n >= 0) {
                        ctx->PMBlock = n;
                        p = skip_int(ctx, p + 6);
                    }
                    else
                        fnd = 0;
                    break;

                case 'c':
                    if (!strncmp(p,"cn=",3)) {
                        p = get_cn(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"crit=%d",&n) == 1 && n >= 1 && n <= 3) {
                        ctx->Crite = n;
                        ctx->CritFlg = 1;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"ccov=%d",&n) == 1 && n >= 1 && n <= 3) {
                        ctx->CCTyp = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (!strncmp(p,"csf,",4)) {
                        ctx->PMCSF = 1;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"cfrac=%lg",&a) == 1 && a >= 0.0 && a <= 1.0) {
                        ctx->PMCFrac = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (!strncmp(p,"con=",4)) {
                        ctx->PMNConS++;
                        p = skip_com(ctx, p);
                    }
                    else if (!strncmp(p,"cen=",4)) {
                        if ((ctx->PMCEN = get_vidx1(ctx, p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (sscanf(p,"cg=%d",&n) == 1) {
                        ctx->PMCG = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"cont=%d",&n) == 1) {
                        ctx->PMCONT = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"cols=%d",&n) == 1 && n > 0) {
                        ctx->PMCols = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"ct=%d",&n) == 1 && n >= 0) {
                        ctx->PMCT = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else
                        fnd = 0;
                    break;

                case 'd':
                    if (!strncmp(p,"des=",4)) {
                        if ((ctx->PMDES = get_vidx1(ctx, p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (!strncmp(p,"drop=",5)) {
                        p = get_var(ctx, p + 5,&r);
                        if (r)
                            goto PARMFin;
                        ctx->PMDrop = 1;
                    }
                    else if (!strncmp(p,"dm=",3)) {
                        p = get_flags(ctx, p + 3,&r,1);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"deg=%d",&n) == 1 && n >= 0) {
                        ctx->PMDEG = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"df",2)) {     /* PMF1d output file */
                        p = get_file(ctx, 1,p + 2,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"dir=%d",&n) == 1) {
                        ctx->PMDIR = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"d=%lg",&a) == 1) {
                        ctx->PMD = a;
                        p = skip_dbl(ctx, p + 2);
                    }
                    else if (sscanf(p,"dscal=%lg",&a) == 1 && a != 0.0) {
                        ctx->DScal = a;
                        ctx->DScalFlg = 1;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (!strncmp(p,"dsv",3)) {    /* dsv file */
                        p = get_fn(ctx, p + 3,ctx->PMDSVName,&r);
                        if (r)
                            goto PARMFin;
                        ctx->PMDSVFlg = 1;
                    }
                    else if (sscanf(p,"dopt=%d",&n) == 1 && n >= 0) {
                        ctx->PMDOPT = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (!strncmp(p,"dtda",4)) {    /* PMTDA file */
                        p = get_file(ctx, 10,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"dspss",5)) {   /* PMSPSS file */
                        p = get_file(ctx, 11,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"dvar",4)) {   /* PMDVAR file */
                        p = get_file(ctx, 12,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"dim=%d",&n) == 1) {
                        ctx->PMDIM = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"dvec=%lf,%lf,%lf",&a,&b,&d) == 3) {
                        ctx->PMDVECX = a;
                        ctx->PMDVECY = b;
                        ctx->PMDVECZ = d;
                        ctx->PMDVECFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (!strncmp(p,"dgrp=",5)) {
                        n = 0;
                        q = p + 4;
                        while (*(q + 1) == '[') {
                            q = skip_nc(ctx, q + 1);
                            if (*(q - 1) != ']')
                                goto PARMFin;
                            n++;
                        }
                        if (q <= p + 7)
                            goto PARMFin;

                        c = *q;
                        *q = '\0';
                        n = (int)(strlen(p + 5));
                        if (n == 0)
                            goto PARMFin;

                        if (rhstr_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(ctx->PMRHSTR,p + 5);
                        *q = c;
                        p = q;
                        ctx->DGRPFlg = n;
                    }
                    else if (sscanf(p,"dx=%lg",&a) == 1) {
                        ctx->PMDX = a;
                        ctx->PMDXFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"dy=%lg",&a) == 1) {
                        ctx->PMDY = a;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"dxa=%lg",&a) == 1) {
                        ctx->PMDXA = a;
                        ctx->PMDXAFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (!strncmp(p,"dblock=",7)) {
                        if ((ctx->PMDBlockV = get_vidx1(ctx, p + 7,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 7;
                    }
                    else
                        fnd = 0;
                    break;

                case 'e':
                    if (sscanf(p,"eps=%lg",&a) == 1 && a > 0.0) {
                        ctx->PMEPS = a;
                        ctx->PMEPSFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"ev=[%d,%d,%d]",&n,&m,&m1) == 3 && n >= 1 && m >= 0 && m1 >= 0) {
                        ctx->PMEVSN = n;
                        ctx->PMEV1 = m;
                        ctx->PMEV2 = m1;
                        p = skip_int(ctx, p + 4);
                        p = skip_int(ctx, p + 1);
                        p = skip_int(ctx, p + 1) + 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 'f':
                    if (sscanf(p,"fmt=%d.%d",&n,&m) == 2) {
                        ctx->PMFmt1 = n;
                        ctx->PMFmt2 = m;
                        p = skip_dbl(ctx, p + 4);
                        ctx->PMFmtF = 1;
                    }
                    else if (sscanf(p,"fs=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMFS = a;
                        ctx->PMFSFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"fss=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMFSS = a;
                        ctx->PMFSSFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"fsx=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMFSX = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"fsy=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMFSY = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (!strncmp(p,"fn=",3)) {
                        p += 3;
                        q = skip_nc(ctx, p);
                        c = *q;
                        *q = '\0';
                        n = (int)(strlen(p));

                        if (ctx->PMFN < PMFNMax) {
                            if (!(ctx->PMFNam[ctx->PMFN] = (char *)calloc((size_t)(n + 1),sizeof(char)))) {
                                err = -2;
                                goto PARMFin;
                            }
                            memrq(ctx, n + 1,sizeof(char));
                            ctx->PMFNA[ctx->PMFN] = n + 1;
                            strcpy(ctx->PMFNam[ctx->PMFN],p);
                            ctx->PMFN++;
                        }
                        *q = c;
                        p = q;
                    }
                    else if (!strncmp(p,"fmt0=",5)) {
                        p = get_fmt(ctx, 0,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"fmt1=",5)) {
                        p = get_fmt(ctx, 1,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"fmt2=",5)) {
                        p = get_fmt(ctx, 2,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"f1=,",3)) {
                        p += 3;
                        q = skip_expr(ctx, p);
                        n = (int)(q - p);
                        if (n == 0)
                            goto PARMFin;

                        if (pmf1_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(ctx->PMF1,p,(size_t)(n));
                        p = q;
                    }
                    else if (!strncmp(p,"f2=,",3)) {
                        p += 3;
                        q = skip_expr(ctx, p);
                        n = (int)(q - p);
                        if (n == 0)
                            goto PARMFin;

                        if (pmf2_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(ctx->PMF2,p,(size_t)(n));
                        p = q;
                    }
                    else if (!strncmp(p,"f3=,",3)) {
                        p += 3;
                        q = skip_expr(ctx, p);
                        n = (int)(q - p);
                        if (n == 0)
                            goto PARMFin;

                        if (pmf3_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(ctx->PMF3,p,(size_t)(n));
                        p = q;
                    }
                    else
                        fnd = 0;
                    break;

                case 'g':
                    if (!strncmp(p,"g=",2)) {
                        p += 2;
                        ctx->PMGIdx = get_vidx1(ctx, p,vname);
                        if (ctx->PMGIdx < 0) {
                            err = -4;
                            goto PARMFin;
                        }
                        p = skip_com(ctx, p);
                    }
                    else if (sscanf(p,"gd=%lg",&a) == 1 && a > 0.0) {
                        ctx->PMGD = a;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"gs=%lg,%lg",&a,&b) == 2) {
                        ctx->PMGS = a;
                        ctx->PMGS1 = b;
                        ctx->PMGSFlg = 2;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"gs=%lg",&a) == 1) {
                        ctx->PMGS = a;
                        ctx->PMGSFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"gt(%d)=%d",&m,&n) == 2 &&
                                    m >= 1 && m <= 3 && n >= 1 && n <= 4) {
                        ctx->PMGT = n;
                        ctx->PMGTT = m;
                        p = skip_int(ctx, p + 3);
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"gt=%d",&n) == 1 && n >= 1 && n <= 4) {
                        ctx->PMGT = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"gn=%d,%d",&n,&m) == 2 && n >= 0 && m >= 0) {
                        ctx->PMGN = n;
                        ctx->PMGN1 = m;
                        p = skip_int(ctx, p + 3);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"gn=%d",&n) == 1 && n >= 0) {
                        ctx->PMGN = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"gnk=%d",&n) == 1 && n >= 1) {
                        ctx->PMGNK = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"gss=",4)) {
                        p = get_tp(ctx, p + 4,&r,2);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"grp=",4)) {
                        p = get_var1(ctx, p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"geo=%d",&n) == 1) {
                        ctx->PMGEO = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"gamma=%lf",&a) == 1) {
                        ctx->PMGamma = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"glen=%d",&n) == 1 && n >= 0) {
                        ctx->PMGLEN = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'h':
                    if (sscanf(p,"hide=%d",&n) == 1) {
                        ctx->PMHIDE = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'i':
                    if (!strncmp(p,"id=",3)) {
                        if ((ctx->PMID = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"ic=%lf",&a) == 1) {
                        ctx->PMIC = a;
                        ctx->PMICFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"idf=%lf,%lf",&a,&b) == 2) {
                        ctx->PMIDFA = a;
                        ctx->PMIDFB = b;
                        p = skip_dbl(ctx, p + 4);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"iv=%d",&n) == 1) {
                        ctx->PMIV = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"idf=%lf",&a) == 1) {
                        ctx->PMIDFA = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (!strncmp(p,"if=",3) || !strncmp(p,"if(",3)) {     /* PMF2d (input file) */

                        p += 2;
                        if (*p == '(') {
                            p++;
                            if (!strncmp(p,"keep",4)) {
                                ctx->PMIFTyp = 1;
                                p += 4;
                            }
                            else if (!strncmp(p,"drop",4)) {
                                ctx->PMIFTyp = 2;
                                p += 4;
                            }
                            else
                                goto PARMFin;
                            if (*p++ != ')')
                                goto PARMFin;
                        }
                        p = get_fn(ctx, p,ctx->PMF2dName,&r);
                        if (r)
                            goto PARMFin;

                        if (!(ctx->PMF2d = fopen(ctx->PMF2dName,OPEN_RD))) {
                            printf1(ctx, "Error: can't open: %s\n",ctx->PMF2dName);
                            err = 1;
                            goto PARMFin;
                        }
                        ctx->PMF2Def = 1;
                    }
                    else if (!strncmp(p,"icost=",6)) {
                        if (sscanf(p,"icost=%lf,%lf",&a,&b) == 2 && a >= 0.0 && b >= 0.0) {
                            ctx->PMICOSTA = a;
                            ctx->PMICOSTB = b;
                            p = skip_dbl(ctx, p + 6);
                            p = skip_dbl(ctx, p + 1);
                        }
                        else if (sscanf(p,"icost=%lf",&a) == 1 && a >= 0.0) {
                            ctx->PMICOSTA = a;
                            p = skip_dbl(ctx, p + 6);
                        }
                        else {
                            p = get_mname(ctx, p + 6,vname,0);
                            if (p == NULL)
                                goto PARMFin;
                            n = mat_getidx(ctx, vname,0);
                            if (n < 0) {
                                printf1(ctx, "Error: undefined matrix name.\n");
                                goto PARMFin;
                            }
                            ctx->PMICOSTMAT = n;
                        }
                    }
                    else if (!strncmp(p,"if1",3)) {    /* PMIF1d input file */
                        p = get_file(ctx, 14,p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"if2",3)) {    /* PMIF2d input file */
                        p = get_file(ctx, 15,p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;

                case 'k':

                    if (sscanf(p,"k=%d",&n) == 1) {
                        ctx->PMK = n;
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"kgam=%lg",&a) == 1 && a > 0.0) {
                        ctx->PMKGam = a;
                        ctx->PMKGamFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (!strncmp(p,"keep=",5)) {
                        p = get_var(ctx, p + 5,&r);
                        if (r)
                            goto PARMFin;
                        ctx->PMKeep = 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 'l':
                    if (sscanf(p,"len=%d",&n) == 1) {
                        ctx->PMLEN = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"link=%d",&n) == 1 && n >= 0) {
                        ctx->PMLINK = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"lt=%d",&n) == 1 && n >= 0) {
                        ctx->PMLT = n;
                        ctx->PMLTFlg = 1;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"lt1=%d",&n) == 1 && n >= 0) {
                        ctx->PMLT1 = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"lw=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMLW = a;
                        ctx->PMLWFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"lw1=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMLW1 = a;
                        ctx->PMLW1Flg = 1;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"log=%d",&n) == 1) {
                        ctx->PMLOG = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"lsecon",6) || !strncmp(p,"lsicon",6)) {
                        ctx->NCONSTR++;
                        p = skip_com(ctx, p);
                    }
                    else if (sscanf(p,"l0=%d",&n) == 1) {
                        ctx->PML0 = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (!strncmp(p,"lon=",4)) {
                        p = get_tp(ctx, p + 4,&r,0);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"lat=",4)) {
                        p = get_tp1(ctx, p + 4,&r,0);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"level=%d",&n) == 1) {
                        ctx->PMLEVEL = n;
                        p = skip_int(ctx, p + 6);
                    }
                    else
                        fnd = 0;
                    break;

                case 'm':
                    if (sscanf(p,"m=%d",&n) == 1) {
                        ctx->PMM = n;
                        ctx->PMMFlg = 1;
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"m1=%d",&n) == 1) {
                        ctx->PMM1 = n;
                        ctx->PMM1Flg = 1;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"m2=%d",&n) == 1) {
                        ctx->PMM2 = n;
                        ctx->PMM2Flg = 1;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"mix=%d",&n) == 1) {
                        ctx->PMMIX = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"mr=%lf",&a) == 1 && a >= 1.0) {
                        ctx->PMMR = a;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"maxcat=%d",&n) == 1 && n > 0) {
                        ctx->PMMaxCat = n;
                        ctx->PMMaxCatFlg = 1;
                        p = skip_int(ctx, p + 7);
                    }
                    else if (sscanf(p,"min=%d",&n) == 1 && n >= 0) {
                        ctx->PMMin = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"max=%d",&n) == 1 && n >= 0) {
                        ctx->PMMax = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"mfmt=%d.%d",&n,&m) == 2) {
                        ctx->PMMFmt1 = n;
                        ctx->PMMFmt2 = m;
                        p = skip_dbl(ctx, p + 5);
                        ctx->PMMFmtF = 1;
                    }
                    else if (sscanf(p,"mina=%d",&n) == 1 && n >= 1 && n <= 8) {
                        ctx->MINA = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"mxit=%d",&n) == 1 && n >= 0) {
                        ctx->MxIter = n;
                        ctx->MxItFlg = 1;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"mxcyc=%d",&n) == 1 && n >= 0) {
                        ctx->PMMXCYC = n;
                        p = skip_int(ctx, p + 6);
                    }
                    else if (sscanf(p,"mxitl=%d",&n) == 1 && n >= 0) {
                        ctx->MxIt1 = n;
                        p = skip_int(ctx, p + 6);
                    }
                    else if (sscanf(p,"msys=%lg",&a) == 1) {
                        ctx->PMMSYS = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (!strncmp(p,"mf=",3)) {
                        p = get_mf(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpcov",5)) {
                        p = get_mp(ctx, p,5,&ctx->PMMPCovDef,ctx->PMMPCovName,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mppar",5)) {
                        p = get_mp(ctx, p,5,&ctx->PMMPParDef,ctx->PMMPParName,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mplog",5)) {
                        p = get_mp(ctx, p,5,&ctx->PMMPLogDef,ctx->PMMPLogName,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpgrad",6)) {
                        p = get_mp(ctx, p,6,&ctx->PMMPGradDef,ctx->PMMPGradName,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mpres",5)) {
                        p = get_mp(ctx, p,5,&ctx->PMMPResDef,ctx->PMMPResName,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"mdef=",5)) {
                        p = get_mname(ctx, p + 5,ctx->PMatName,0);
                        if (p == NULL || (*p != ',' && *p != ')'))
                            goto PARMFin;
                        ctx->PMatNameFlg = 1;
                    }
/* ## */            else if (!strncmp(p,"mod=",4)) {
                        p = get_mod(ctx, p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"msg=%d",&n) == 1 && n >= 0) {
                        ctx->PMMSG = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"meth=%d",&n) == 1 && n >= 1) {
                        ctx->PMMETH = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'n':
                    if (sscanf(p,"noc=%d",&n) == 1 && n >= 1) {
                        ctx->PMNOC = n;
                        ctx->PMNOCFlg = 1;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"nfmt=%d",&n) == 1 && n > 0) {
                        ctx->PMNFmt = n;
                        p = skip_int(ctx, p + 5);
                        ctx->PMNFmtF = 1;
                    }
                    else if (sscanf(p,"n=%d",&n) == 1 && n >= 1) {
                        ctx->PMN = n;
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"nhp=%d",&n) == 1 && n >= 1) {
                        ctx->PMNHP = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"nmp=%d",&n) == 1 && n >= 1) {
                        ctx->PMNMPnt = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"ns=%d",&n) == 1 && n >= 0) {
                        ctx->PMNS = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"ndim=%d",&n) == 1 && n >= 1) {
                        ctx->PMNDIM = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"nc=%d",&n) == 1) {
                        ctx->PMNC = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"nw=%d",&n) == 1 && n >= 1) {
                        ctx->PMNW = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"nq=%d",&n) == 1 && n >= 1) {
                        ctx->PMNQ = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"ni=%d",&n) == 1 && n >= 0) {
                        ctx->PMNI = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"nbox=%d",&n) == 1 && n >= 1) {
                        ctx->PMNBOX = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"nn=%d,%d",&n,&m) == 2) {
                        ctx->PMNNFlg = 1;
                        ctx->PMNN1 = n;
                        ctx->PMNN2 = m;
                        p = skip_int(ctx, p + 3);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"ndigit=%lg",&a) == 1 && a >= 1.0) {
                        ctx->PMNDIGIT = a;
                        ctx->PMNDIGFlg = 1;
                        p = skip_dbl(ctx, p + 7);
                    }
                    else if (sscanf(p,"np=%d",&n) == 1) {
                        ctx->PMNP = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"nlev=%d",&n) == 1 && n >= 1) {
                        ctx->PMNLEV = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else
                        fnd = 0;
                    break;

                case 'o':
                    if (sscanf(p,"opt=%d",&n) == 1 && n >= 1) {
                        ctx->PMOPT = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"org=",4)) {
                        if ((ctx->PMORG = get_vidx1(ctx, p + 4,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 4;
                    }
                    else if (sscanf(p,"order=%d",&n) == 1) {
                        ctx->PMOrder = n;
                        p = skip_int(ctx, p + 6);
                    }
                    else if (sscanf(p,"off=%lf",&a) == 1) {
                        ctx->PMOFF = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'p':
                    if (!strncmp(p,"ps=",3)) {
                        p = get_pattern(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"prn=%d",&n) == 1 && n >= 0) {
                        ctx->PMPRNO = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (sscanf(p,"plot=%d",&n) == 1 && n >= 0) {
                        ctx->PMPLOT = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"ptyp=%d",&n) == 1 && n >= 1) {
                        ctx->PMPTyp = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"pos=%lf,%lf",&a,&b) == 2) {
                        ctx->PMX = a;
                        ctx->PMY = b;
                        ctx->PMXYFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"pl=%d",&n) == 1 && n >= 0) {
                        ctx->PMPLFlg = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (!strncmp(p,"prate",5)) {
                        p += 5;
                        if (*p == '(')
                            p = skip_blev(ctx, p);

                        if (*p++ != '=')
                            goto PARMFin;
                        p = skip_com(ctx, p);
                        ctx->PMPRN++;
                    }
                    else if (!strncmp(p,"prot",4)) {  /* PMProt file */
                        p = get_file(ctx, 5,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ppar",4)) {    /* PMPPFd file */
                        p = get_file(ctx, 4,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pcov",4)) {    /* PMCov file */
                        p = get_file(ctx, 2,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pres",4)) {    /* PMRes file */
                        p = get_file(ctx, 3,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"pfn",3)) {    /* PMPFN file */
                        p = get_file(ctx, 7,p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"pfmt=%d.%d",&n,&m) == 2) {
                        ctx->PMPFmt1 = n;
                        ctx->PMPFmt2 = m;
                        p = skip_dbl(ctx, p + 5);
                        ctx->PMPFmtF = 1;
                    }
                    else if (sscanf(p,"pmin=%d",&n) == 1 && n >= 1) {
                        ctx->PMPMin = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (!strncmp(p,"pmval=",6)) {
                        p += 5;
                        while (sscanf(p + 1,"%lf",&a) == 1) {
                            ctx->PMPMVal[ctx->PMPMN] = a;
                            p = skip_dbl(ctx, p + 1);
                            if (++ctx->PMPMN >= 10)
                                break;
                        }
                    }
                    else if (!strncmp(p,"ptab1",5)) {    /* PMTab1 file */
                        p = get_file(ctx, 9,p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ptab",4)) {    /* PMTab file */
                        p = get_file(ctx, 8,p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"pxlen=%lg",&a) == 1 && a > 0.0) {
                        ctx->PMXLen = a;
                        ctx->PMXLenFlg = 1;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"pylen=%lg",&a) == 1 && a > 0.0) {
                        ctx->PMYLen = a;
                        ctx->PMYLenFlg = 1;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (!strncmp(p,"pcf",3)) {     /* PMPCFd output file */
                        p = get_file(ctx, 6,p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"perm=%d",&n) == 1) {
                        ctx->PMPERM = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"proj=%d",&n) == 1 && n >= 1) {
                        ctx->PMPROJ = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"psorg=%d,%d",&n,&m) == 2) {
                        ctx->PMXOrg = n;
                        ctx->PMYOrg = m;
                        p = skip_int(ctx, p + 6);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"psrot=%lg",&a) == 1) { /*      && a >= 0.0 && a < 360.0) {   */
                        ctx->PMPSRot = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"pcheck=%d",&n) == 1 && n >= 0) {
                        ctx->PMPCheck = n;
                        p = skip_int(ctx, p + 7);
                    }
                    else
                        fnd = 0;
                    break;

                case 'q':
                    if (!strncmp(p,"qo=",3)) {
                        p = get_tp(ctx, p + 3,&r,3);
                        if (r)
                            goto PARMFin;
                        ctx->PMQOFlg = 1;
                    }
                    else if (!strncmp(p,"qt=",3)) {
                        p = get_tp(ctx, p + 3,&r,1);
                        if (r)
                            goto PARMFin;
                        ctx->PMQTFlg = 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 'r':
                    if (sscanf(p,"r=%d",&n) == 1) {
                        ctx->PMR = n;
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"rd=%lf,%lf",&a,&b) == 2) {
                        ctx->PMRDA = a;
                        ctx->PMRDB = b;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"rd=%lf",&a) == 1) {
                        ctx->PMRDA = a;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"rt=%d",&n) == 1 && n >= 0) {
                        ctx->PMRT = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"rr=%d,%d",&n,&m) == 2 && n >= 1 && m >= n) {
                        ctx->PMRRN = n;
                        ctx->PMRRM = m;
                        p = skip_int(ctx, p + 3);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"rr=%d",&n) == 1 && n >= 1) {
                        ctx->PMRRM = ctx->PMRRN = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"rx=%lf(%lf)%lf",&a,&d,&b) == 3 && a < b && d > 0.0) {
                        ctx->PMRXFlg = 1;
                        ctx->PMRXA = (float)a;
                        ctx->PMRXB = (float)b;
                        ctx->PMRXD = (float)d;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"ry=%lf(%lf)%lf",&a,&d,&b) == 3 && a < b && d > 0.0) {
                        ctx->PMRYFlg = 1;
                        ctx->PMRYA = (float)a;
                        ctx->PMRYB = (float)b;
                        ctx->PMRYD = (float)d;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"ru=%lf,%lf,%d,%d",&a,&b,&n,&m) == 4 && a < b && n > 0 && m > 0) {
                        ctx->PMRUFlg = 1;
                        ctx->PMRUA = a;
                        ctx->PMRUB = b;
                        ctx->PMRUN = n;
                        ctx->PMRUM = m;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_int(ctx, p + 1);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"rv=%lf,%lf,%d,%d",&a,&b,&n,&m) == 4 && a < b && n > 0 && m > 0) {
                        ctx->PMRVFlg = 1;
                        ctx->PMRVA = a;
                        ctx->PMRVB = b;
                        ctx->PMRVN = n;
                        ctx->PMRVM = m;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_int(ctx, p + 1);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (!strncmp(p,"rcn=",4)) {
                        p += 4;
                        q = p;
                        while (*q) {
                            q = skip_int(ctx, q);
                            if (*q++ != '[')
                                goto PARMFin;

                            while (*q) {
                                q = skip_int(ctx, q);
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
                        if (prcn_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(ctx->PRCN,p,(size_t)(n));
                        *(ctx->PRCN + n) = '\0';
                        p = q;
                    }
                    else if (!strncmp(p,"rc=",3)) {
                        p += 3;
                        q = p;
                        while (*q) {
                            if (sscanf(q,"%d",&n) != 1)
                                goto PARMFin;
                            q = skip_int(ctx, q);
                            if (*q++ != '[')
                                goto PARMFin;

                            while (*q) {
                                q = skip_int(ctx, q);
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
                        if (prc_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strncpy(ctx->PRC,p,(size_t)(n));
                        *(ctx->PRC + n) = '\0';
                        p = q;
                    }
                    else if (!strncmp(p,"res=",4)) {
                        p += 3;
                        while (sscanf(p + 1,"%d",&n) == 1) {
                            if (ctx->PMResN < 10)
                                ctx->PMRes[ctx->PMResN++] = n;
                            p = skip_int(ctx, p + 1);
                        }
                    }
                    else if (!strncmp(p,"rrisk,",6)) {
                        ctx->PMRRFlg = 1;
                        p += 5;
                    }
                    else if (sscanf(p,"rerr=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMRERR = a;
                        ctx->PMRERRFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"region=%lf,%lf",&a,&b) == 2) {
                        ctx->PMRegion1 = a;
                        ctx->PMRegion2 = b;
                        ctx->PMRegionFlg = 1;
                        p = skip_dbl(ctx, p + 7);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"rows=%d",&n) == 1 && n > 0) {
                        ctx->PMRows = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (sscanf(p,"rec=%lf,%lf,%lf,%lf",&ctx->PMRECXMin,&ctx->PMRECYMin,&ctx->PMRECXMax,&ctx->PMRECYMax) == 4) {
                        ctx->PMRECFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (!strncmp(p,"rel=",4)) {
                        p = get_flags(ctx, p + 4,&r,5);
                        if (r)
                            goto PARMFin;
                    }
                    else if (sscanf(p,"rhem=%lf",&a) == 1 && a > 0.0) {
                        ctx->PMRHem = a;
                        p = skip_dbl(ctx, p + 5);
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
                        if (n == 0 || *q != ']')
                            goto PARMFin;

                        if (rhstr_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        *q = '\0';
                        strcpy(ctx->PMRHSTR,p);
                        *q = ']';
                        p += n + 1;
                    }
#ifdef TDA_R_PACKAGE
                    /* ivreg1's search box.  The standalone command uses a
                       box that is hard-coded in ivreg1() -- values from the
                       dataset the command was being developed against --
                       so any data whose solution lies outside it cannot be
                       found.  The R package passes the box explicitly:
                       four lower bounds then four upper bounds for
                       (alpha center, alpha radius, beta center, beta
                       radius).  Guarded so the standalone build is
                       unchanged. */
                    else if (sscanf(p,"sbox=%lg,%lg,%lg,%lg,%lg,%lg,%lg,%lg",
                                    ctx->IVSBox,ctx->IVSBox + 1,
                                    ctx->IVSBox + 2,ctx->IVSBox + 3,
                                    ctx->IVSBox + 4,ctx->IVSBox + 5,
                                    ctx->IVSBox + 6,ctx->IVSBox + 7) == 8) {
                        ctx->IVSBoxSet = 1;
                        p = skip_dbl(ctx, p + 5);
                        for (m = 0; m < 7; ++m)
                            p = skip_dbl(ctx, p + 1);
                    }
#endif
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
                        if (pmstr_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(ctx->PMSTR,p0);
                        *q = c;
                        p = q;
                        if (m)
                            p++;
                    }
                    else if (sscanf(p,"sc=%lf",&a) == 1) {
                        ctx->PMSC = a;
                        ctx->PMSCFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"s=%d",&n) == 1) {
                        ctx->PMS = n;
                        p = skip_int(ctx, p + 2);
                    }
                    else if (sscanf(p,"sd=%d",&n) == 1) {
                        ctx->PMSD = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (sscanf(p,"sn=%d,%d",&n,&m) == 2 && n >= 0 && m >= 0) {
                        ctx->PMSN = n;
                        ctx->PMSN1 = m;
                        p = skip_int(ctx, p + 3);
                        p = skip_int(ctx, p + 1);
                    }
                    else if (sscanf(p,"sn=%d",&n) == 1 && n >= 0) {
                        ctx->PMSN = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (!strncmp(p,"sm=",3)) {
                        p = get_flags(ctx, p + 3,&r,2);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sort=",5)) {
                        p = get_var1(ctx, p + 5,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sort,",5)) {
                        ctx->PMSORTFlg = 1;
                        p += 4;
                    }
                    else if (sscanf(p,"sig=%lf",&a) == 1) {
                        ctx->PMSIG = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (!strncmp(p,"scost=",6)) {
                        if (sscanf(p,"scost=%d",&n) == 1 && n >= 0) {
                            ctx->PMSCOSTM = n;
                            p = skip_int(ctx, p + 6);
                        }
                        else {
                            p = get_mname(ctx, p + 6,vname,0);
                            if (p == NULL)
                                goto PARMFin;
                            n = mat_getidx(ctx, vname,0);
                            if (n < 0) {
                                printf1(ctx, "Error: undefined matrix name.\n");
                                goto PARMFin;
                            }
                            ctx->PMSCOSTMAT = n;
                        }
                    }
                    else if (sscanf(p,"seed=%d",&n) == 1) {
                        ctx->PMSEED = n;
                        p = skip_int(ctx, p + 5);
                    }
                    else if (!strncmp(p,"sel=",4)) {

                        q = skip_expr(ctx, p + 4);
                        c = *q;
                        *q = '\0';

                        if ((n = v_parse(ctx, p + 4,0)) < 0 || ctx->ESCnt <= 0) {
                            printf1(ctx, "Syntax error (%d) in sel expression.\n",n);
                            if (n < 0)
                                prn_emsg1(ctx, n);
                            err = 1;
                            goto PARMFin;
                        }
                        if (sve_alloc(ctx, ctx->ESCnt,p)) {
                            err = -2;
                            goto PARMFin;
                        }
                        ctx->SVECnt = ctx->ESCnt;
                        for (i = 0; i < ctx->ESCnt; ++i) {
                            ctx->SVETyp[i] = ctx->ESTyp[i];
                            ctx->SVEVal[i] = ctx->ESVal[i];
                        }
                        ctx->SVEFlg = 1;
                        *q = c;
                        p = q;
                    }
                    else if (sscanf(p,"smin=%lg",&a) == 1 && a >= 0.0) {
                        ctx->SMin = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"slen=%lg",&a) == 1 && a > 0.0) {
                        ctx->SLen = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"sred=%lg",&a) == 1 && a > 0.0 && a < 1.0) {
                        ctx->SRed = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"sst=%d",&n) == 1 && n >= 0 && n <= 1) {
                        ctx->STFlg = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else if (!strncmp(p,"sk=",3)) {
                        p = get_flags(ctx, p + 3,&r,4);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"sepc=none",9)) {
                        ctx->XSEPC = '\0';
                        p += 9;
                    }
                    else if (sscanf(p,"sepc=%c",&ctx->XSEPC) == 1) {
                        if (ctx->XSEPC == 't')
                            ctx->XSEPC = '\t';
                        p += 6;
                    }
                    else if (!strncmp(p,"screen",6)) {
                        ctx->SCRNFlg = 1;
                        p += 6;
                    }
                    else if (!strncmp(p,"scale=",6)) {
                        if ((ctx->PMSCAL = get_vidx1(ctx, p + 6,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 6;
                    }
                    else if (sscanf(p,"size=%lg",&a) == 1) {
                        ctx->PMSize = a;
                        ctx->PMSizeFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"sdfmt=%d.%d",&n,&m) == 2) {
                        ctx->PMSDFmt1 = n;
                        ctx->PMSDFmt2 = m;
                        p = skip_dbl(ctx, p + 6);
                        ctx->PMSDFmtF = 1;
                    }
                    else
                        fnd = 0;
                    break;

                case 't':
                    if (!strncmp(p,"t=",2)) {
                        p = get_tp(ctx, p + 2,&r,1);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"tp=",3)) {
                        p = get_tp(ctx, p + 3,&r,1);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"tst=",4)) {
                        p = get_flags(ctx, p + 4,&r,3);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"ts=",3)) {
                        if ((ctx->PMTS = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"tf=",3)) {
                        if ((ctx->PMTF = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"tl=%lf",&a) == 1 && a >= 0.0) {
                        ctx->PMTL = a;
                        p = skip_dbl(ctx, p + 3);
                    }
                    else if (sscanf(p,"tfmt=%d.%d",&n,&m) == 2) {
                        ctx->PMTFmt1 = n;
                        ctx->PMTFmt2 = m;
                        p = skip_dbl(ctx, p + 5);
                        ctx->PMTFmtF = 1;
                    }
                    else if (sscanf(p,"tolg=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLG = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"tolf=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLF = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"tolp=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLP = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"tolv=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLV = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"tols=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLS = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"tolsg=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLSG = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"tolsp=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLSP = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"tolbw=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLBW = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"tolfd=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLFD = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"tolfe=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLFE = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (sscanf(p,"tolbc=%lg",&a) == 1 && a > 0.0) {
                        ctx->TOLBC = a;
                        p = skip_dbl(ctx, p + 6);
                    }
                    else if (!strncmp(p,"transp",6)) {
                        ctx->PMTransp = 1;
                        p += 6;
                    }
                    else if (!strncmp(p,"trunc=",6)) {
                        if ((ctx->PMTRUNC = get_vidx1(ctx, p + 6,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 6;
                    }
                    else if (sscanf(p,"tol=%lg",&a) == 1 && a >= 0.0) {
                        ctx->PMTol = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"typ=%d",&n) == 1 && n >= 1) {
                        ctx->PMTyp = n;
                        p = skip_int(ctx, p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'u':
                    if (sscanf(p,"ulx=%lg",&a) == 1) {
                        ctx->PMULX = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else if (sscanf(p,"uly=%lg",&a) == 1) {
                        ctx->PMULY = a;
                        p = skip_dbl(ctx, p + 4);
                    }
                    else
                        fnd = 0;
                    break;

                case 'v':
                    if (!strncmp(p,"v=",2)) {
                        p = get_var(ctx, p + 2,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"vsel=",5)) {

                        q = skip_expr(ctx, p + 5);
                        c = *q;
                        *q = '\0';
                        n = (int)(strlen(p));
                        if (rhstr_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(ctx->PMRHSTR,p);
                        *q = c;
                        p = q;
                    }
                    else if (sscanf(p,"view=%lg,%lg",&a,&b) == 2 &&
                          -180.0 <= a && a <= 180.0 && -90.0 <= b && b <= 90.0) {
                        ctx->PMViewLon = a;
                        ctx->PMViewLat = b;
                        ctx->PMViewFlg = 1;
                        p = skip_dbl(ctx, p + 5);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else
                        fnd = 0;
                    break;

                case 'w':
                    if (!strncmp(p,"w=",2)) {
                        if ((ctx->PMWVar = get_vidx1(ctx, p + 2,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 2;
                    }
                    else if (sscanf(p,"wf=%d",&n) == 1 && n >= 0) {
                        ctx->PMWF = n;
                        p = skip_int(ctx, p + 3);
                    }
                    else if (!strncmp(p,"wt=",3)) {
                        p = get_tp(ctx, p + 3,&r,1);
                        if (r)
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;

                case 'x':
                    if (sscanf(p,"xy=%lf,%lf",&a,&b) == 2) {
                        ctx->PMX = a;
                        ctx->PMY = b;
                        ctx->PMXYFlg = 1;
                        p = skip_dbl(ctx, p + 3);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (sscanf(p,"xyz=%lf,%lf,%lf",&a,&b,&d) == 3) {
                        ctx->PM3X = a;
                        ctx->PM3Y = b;
                        ctx->PM3Z = d;
                        ctx->PMXYZFlg = 1;
                        p = skip_dbl(ctx, p + 4);
                        p = skip_dbl(ctx, p + 1);
                        p = skip_dbl(ctx, p + 1);
                    }
                    else if (!strncmp(p,"x=",2)) {
                        p = get_tp(ctx, p + 2,&r,0);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) ||
                             !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3)) {
                        p = skip_xa(ctx, p);
                        ctx->PMNXA++;
                    }
                    else if (!strncmp(p,"xp=",3)) {
                        p = get_xp(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xe=",3)) {
                        m = 0;
                        q = p + 2;
                        while (*(q + 1) == '[') {
                            q = skip_nc(ctx, q + 1);
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
                        n = (int)(strlen(p + 3));
                        if (n == 0)
                            goto PARMFin;

                        if (rhstr_alloc(ctx, n + 1)) {
                            err = -2;
                            goto PARMFin;
                        }
                        strcpy(ctx->PMRHSTR,p + 3);
                        *q = c;
                        p = q;
                        ctx->XEFlg = m;
                    }
                    else if (!strncmp(p,"xv=",3)) {
                        p = get_var1(ctx, p + 3,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"xyv=",4)) {
                        p = get_var(ctx, p + 4,&r);
                        if (r)
                            goto PARMFin;
                    }
                    else
                        fnd = 0;
                    break;

                case 'y':
                    if (!strncmp(p,"y=",2)) {
                        p = get_tp(ctx, p + 2,&r,1);
                        if (r)
                            goto PARMFin;
                    }
                    else if (!strncmp(p,"yw=",3)) {
                        if ((ctx->PMYWVar = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"yl=",3)) {
                        if ((ctx->PMYL = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (!strncmp(p,"yh=",3)) {
                        if ((ctx->PMYH = get_vidx1(ctx, p + 3,vname)) < 0)
                            goto PARMFin1;
                        p += strlen(vname) + 3;
                    }
                    else if (sscanf(p,"year=%d,%d",&n,&m) == 2 && n >= 0 && m >= n) {
                        ctx->PMYEAR1 = n;
                        ctx->PMYEAR2 = m;
                        p = skip_int(ctx, p + 5);
                        p = skip_int(ctx, p + 1);
                    }
                    else
                        fnd = 0;
                    break;

                case 'z':
                    if (sscanf(p,"zmin=%lg",&a) == 1) {
                        ctx->PMZMin = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (sscanf(p,"zval=%lg",&a) == 1) {
                        ctx->PMZVal = a;
                        p = skip_dbl(ctx, p + 5);
                    }
                    else if (!strncmp(p,"zvar=",5)) {
                        if ((ctx->PMZVar = get_vidx1(ctx, p + 5,vname)) < 0)
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
            strcpy(ctx->PMFdName,p);
            i = 0;
            if (opt == 2) {
                if (!(ctx->PMFd = fopen(p,OPEN_RD)))
                    i = 1;
            }
            else if (opt == 10) {
                if (!(ctx->PMFd = fopen(p,OPEN_RB)))
                    i = 1;
            }
            else if (opt == 13) {
                if (!(ctx->PMFd = fopen(p,OPEN_WB)))
                    i = 1;
            }
            else if (ctx->PMAP == 0) {
                if (!(ctx->PMFd = fopen(p,OPEN_WR)))
                    i = 1;
            }
            else if (!(ctx->PMFd = fopen(p,OPEN_AP)))
                i = 1;

            if (i) {
                printf1(ctx, "Error: can't open: %s\n",p);
                err = 1;
                goto PARMFin;
            }
            ctx->PMFDef = 1;
        }
        else if (opt == 3) {
            if (sscanf(p,"%d",&n) != 1)
                goto PARMFin;
            ctx->PMRHSI = n;
        }
        else if (opt == 4) {        /* must be a varlist */
            p = get_var(ctx, p,&err);
            if (*p)
                err = -1;
            if (err)
                goto PARMFin;
        }
        else if (opt == 5) {
            if (sscanf(p,"%lf,%lf",&a,&b) != 2)
                goto PARMFin;
            ctx->PMRHSA = a;
            ctx->PMRHSB = b;
        }
        else if (opt == 6) {
            if (sscanf(p,"%lf(%lf)%lf",&a,&d,&b) != 3)
                goto PARMFin;
            ctx->PMRHSA = a;
            ctx->PMRHSB = b;
            ctx->PMRHSD = d;
        }
        else if (opt == 7) {        /* list of double values */
            q = p;
            n = 1;
            while (*q) {
                if (*q++ == ',')
                    n++;
            }
            if (rhs_alloc(ctx, n)) {
                err = 2;
                goto PARMFin;
            }
            for (i = 0; i < n; ++i) {
                if (sscanf(p,"%lf",&a) != 1)
                    goto PARMFin;
                ctx->PMRHSX[i] = a;
                p = skip_dbl(ctx, p);
                if (i == n - 1)
                    break;
                if (*p++ != ',')
                    goto PARMFin;
            }
            ctx->PMRHSN = n;
        }
        else if (opt == 8) {        /* string */
            n = (int)(strlen(p));
            if (rhstr_alloc(ctx, n + 1)) {
                err = -2;
                goto PARMFin;
            }
            strcpy(ctx->PMRHSTR,p);
        }
        else if (opt == 9) {        /* function */

            err = get_func(ctx, p,1,&n,0,NULL);
            if (err) {
                err = 1;
                goto PARMFin;
            }
            ctx->PMFTYP5 = n;
        }
        else if (opt == 11) {        /* expression */

            if ((n = v_parse(ctx, p,0)) < 0 || ctx->ESCnt <= 0) {
                printf1(ctx, "Syntax error (%d) in right-hand side expression.\n",n);
                if (n < 0)
                    prn_emsg1(ctx, n);
                err = 1;
                goto PARMFin;
            }
            p = skip_expr(ctx, p);
        }
        else if (opt == 12) {        /* inclusion function */

            err = get_func(ctx, p,1,&n,1,NULL);
            if (err) {
                err = 1;
                goto PARMFin;
            }
            ctx->PMFTYP5 = n;
        }
        else if (opt == 14) { /* must be a list of variable or matrix names */
            p = get_varx(ctx, p,&err);
            if (*p)
                err = -1;
            if (err)
                goto PARMFin;
        }
        ctx->PMRHSFlg = 1;
    }
    else if (*p)
        goto PARMFin;
    else if (rs) {
        printf1(ctx, "Error: need right-hand side of command.\n");
        err = 1;
        goto PARMFin;
    }

    makefmt(ctx, &ctx->PMFmt1,&ctx->PMFmt2,ctx->PMFmtS,sizeof(ctx->PMFmtS),0,ctx->SEPC,0);
    makefmt(ctx, &ctx->PMTFmt1,&ctx->PMTFmt2,ctx->PMTFmtS,sizeof(ctx->PMTFmtS),0,ctx->SEPC,0);
    makefmt(ctx, &ctx->PMMFmt1,&ctx->PMMFmt2,ctx->PMMFmtS,sizeof(ctx->PMMFmtS),0,ctx->SEPC,0);
    makefmt(ctx, &ctx->PMPFmt1,&ctx->PMPFmt2,ctx->PMPFmtS,sizeof(ctx->PMPFmtS),0,ctx->SEPC,0);
    makefmt(ctx, &ctx->PMSDFmt1,&ctx->PMSDFmt2,ctx->PMSDFmtS,sizeof(ctx->PMSDFmtS),0,ctx->SEPC,0);
    makenfmt(ctx, &ctx->PMNFmt,ctx->PMNFmtS,sizeof(ctx->PMNFmtS),1,ctx->SEPC);

    if (ctx->PMNVTyp == 0) {
        ctx->PMNVLEN = imax(ctx, get_mxvlen(ctx, ctx->PMNV,ctx->PMVIdx),get_mxvlen(ctx, ctx->PMNZ,ctx->PMZIdx));
        ctx->PM1NVLEN = get_mxvlen(ctx, ctx->PM1NV,ctx->PM1VIdx);
        ctx->PM2NVLEN = get_mxvlen(ctx, ctx->PM2NV,ctx->PM2VIdx);
        ctx->PM3NVLEN = get_mxvlen(ctx, ctx->PM3NV,ctx->PM3VIdx);
    }
    err = 0;

PARMFin:                        /* ## */
    if (err < 0) {
        /***************
        p = skip_nc(ctx, p0);
        *p = '\0';
        *************/
        prn_perr(ctx, p0);
        p_err(ctx, err,1);
    }
    if (err)
        err = -1;
    return(err);

PARMFin1:               /* ## */
    p = skip_nc(ctx, p0);
    *p = '\0';
    prn_perr(ctx, p0);
    printf1(ctx, "Undefined variable.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_perr    Print error message.                                        */

void prn_perr(TDAContext *ctx, char *s)
{
    register char *p = s;

    printf1(ctx, "Error: ");
    while (*p && p < s + 20)
        printf1(ctx, "%c",*p++);
    if (*p)
        printf1(ctx, " ...");
    printf1(ctx, "\n");
}

/* ------------------------------------------------------------------------ */
/*  pmfn_free()     Free previously allocated memory PMFNam[].              */

void pmfn_free(TDAContext *ctx)
{
    register int i;

    for (i = 0; i < ctx->PMFN; ++i) {
        if (ctx->PMFNA[i] > 0) {
            free(ctx->PMFNam[i]);
            memrq(ctx, -ctx->PMFNA[i],sizeof(char));
            ctx->PMFNA[i] = 0;
        }
    }
    ctx->PMFN = 0;
}

/* ------------------------------------------------------------------------ */
/*  sve_alloc(n,s)  If n > 0 allocate parser stack for sel expression,      */
/*                  otherwise free. Save string s in SVESTR.                */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int sve_alloc(TDAContext *ctx, int n,char *s)
{
    int l;

    if (ctx->SVESTRA > 0) {
        free((char *)ctx->SVESTR);
        memrq(ctx, -ctx->SVESTRA,sizeof(char));
        ctx->SVESTRA = 0;
    }
    if (ctx->SVETypA > 0) {
        free((char *)ctx->SVETyp);
        memrq(ctx, -ctx->SVETypA,sizeof(int));
        ctx->SVETypA = 0;
    }
    if (ctx->SVEValA > 0) {
        free((char *)ctx->SVEVal);
        memrq(ctx, -ctx->SVEValA,sizeof(double));
        ctx->SVEValA = 0;
    }
    ctx->SVECnt = (short)(ctx->SVEFlg = 0);
    if (n <= 0)
        return(0);

    if (!(ctx->SVETyp = (int *)calloc((size_t)(n),sizeof(int))))
        return(-2);
    memrq(ctx, n,sizeof(int));
    ctx->SVETypA = n;

    if (!(ctx->SVEVal = (double *)calloc((size_t)(n),sizeof(double))))
        return(-2);
    memrq(ctx, n,sizeof(double));
    ctx->SVEValA = n;

    if (s != NULL) {
        l = (int)(strlen(s));
        if (l > 0) {
            if (!(ctx->SVESTR = (char *)calloc((size_t)(l + 1),sizeof(char))))
                return(-2);
            memrq(ctx, l + 1,sizeof(char));
            ctx->SVESTRA = l + 1;
            strcpy(ctx->SVESTR,s);
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rhstr_alloc(n)  If n > 0 allocate PMRHSTR else free previously          */
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int rhstr_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMRHSTRA > 0) {
        free(ctx->PMRHSTR);
        memrq(ctx, -ctx->PMRHSTRA,sizeof(char));
        ctx->PMRHSTRA = 0;
        ctx->DGRPFlg = ctx->XEFlg = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMRHSTR = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);
    memrq(ctx, n,sizeof(char));
    ctx->PMRHSTRA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmstr_alloc(n)  If n > 0 allocate PMSTR else free previously            */
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmstr_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMSTRA > 0) {
        free(ctx->PMSTR);
        memrq(ctx, -ctx->PMSTRA,sizeof(char));
        ctx->PMSTRA = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMSTR = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);
    memrq(ctx, n,sizeof(char));
    ctx->PMSTRA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf1_alloc(n)   If n > 0 allocate PMF1 else free previously             */
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf1_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMF1A > 0) {
        free(ctx->PMF1);
        memrq(ctx, -ctx->PMF1A,sizeof(char));
        ctx->PMF1A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMF1 = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);
    memrq(ctx, n,sizeof(char));
    ctx->PMF1A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf2_alloc(n)   If n > 0 allocate PMF2 else free previously             */
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf2_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMF2A > 0) {
        free(ctx->PMF2);
        memrq(ctx, -ctx->PMF2A,sizeof(char));
        ctx->PMF2A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMF2 = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);
    memrq(ctx, n,sizeof(char));
    ctx->PMF2A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pmf3_alloc(n)   If n > 0 allocate PMF3 else free previously             */
/*                  allocated array.                                        */
/*                                                                          */
/*  Return: 0 if OK, -2 if insufficient memory.                             */

int pmf3_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMF3A > 0) {
        free(ctx->PMF3);
        memrq(ctx, -ctx->PMF3A,sizeof(char));
        ctx->PMF3A = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMF3 = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);
    memrq(ctx, n,sizeof(char));
    ctx->PMF3A = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rhs_alloc(n)    If n > 0 allocate PMRHSX else free previously allocated */
/*                  PMRHSX.                                                 */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int rhs_alloc(TDAContext *ctx, int n)
{
    if (ctx->PMRHSXA > 0) {
        free((char *)ctx->PMRHSX);
        memrq(ctx, -ctx->PMRHSXA,sizeof(double));
        ctx->PMRHSN = ctx->PMRHSXA = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PMRHSX = (double *)calloc((size_t)(n),sizeof(double))))
        return(-2);
    memrq(ctx, n,sizeof(double));
    ctx->PMRHSXA = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prc_alloc(n)    If n > 0 allocate PRC else free previously allocated    */
/*                  PRC.                                                    */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int prc_alloc(TDAContext *ctx, int n)
{
    if (ctx->PRCAlloc > 0) {
        free(ctx->PRC);
        memrq(ctx, -ctx->PRCAlloc,sizeof(char));
        ctx->PRCAlloc = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PRC = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);

    memrq(ctx, n,sizeof(char));
    ctx->PRCAlloc = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prcn_alloc(n)   If n > 0 allocate PRCN else free previously allocated   */
/*                  PRCN.                                                   */
/*                                                                          */
/*  Return: 0 if OK, or -2 if insufficient memory.                          */

int prcn_alloc(TDAContext *ctx, int n)
{
    if (ctx->PRCNAlloc > 0) {
        free(ctx->PRCN);
        memrq(ctx, -ctx->PRCNAlloc,sizeof(char));
        ctx->PRCNAlloc = 0;
    }
    if (n <= 0)
        return(0);

    if (!(ctx->PRCN = (char *)calloc((size_t)(n),sizeof(char))))
        return(-2);

    memrq(ctx, n,sizeof(char));
    ctx->PRCNAlloc = n;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_valloc()     If n > 0 allocate n elements in PMVIdx[], otherwise     */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_valloc(TDAContext *ctx, int n)
{
    if (ctx->PMNV > 0) {
        free((char *)ctx->PMVIdx);
        memrq(ctx, -ctx->PMNV,sizeof(short));
        ctx->PMNVTyp = ctx->PMNV = 0;
    }
    if (n > 0) {
        if (!(ctx->PMVIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->PMNV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_zalloc()     If n > 0 allocate n elements in PMZIdx[], otherwise     */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_zalloc(TDAContext *ctx, int n)
{
    if (ctx->PMNZ > 0) {
        free((char *)ctx->PMZIdx);
        memrq(ctx, -ctx->PMNZ,sizeof(short));
        ctx->PMNZ = 0;
    }
    if (n > 0) {
        if (!(ctx->PMZIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->PMNZ = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v1alloc()    If n > 0 allocate n elements in PM1VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v1alloc(TDAContext *ctx, int n)
{
    if (ctx->PM1NV > 0) {
        free((char *)ctx->PM1VIdx);
        memrq(ctx, -ctx->PM1NV,sizeof(short));
        ctx->PM1NV = 0;
    }
    if (n > 0) {
        if (!(ctx->PM1VIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->PM1NV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v2alloc()    If n > 0 allocate n elements in PM2VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v2alloc(TDAContext *ctx, int n)
{
    if (ctx->PM2NV > 0) {
        free((char *)ctx->PM2VIdx);
        memrq(ctx, -ctx->PM2NV,sizeof(short));
        ctx->PM2NV = 0;
    }
    if (n > 0) {
        if (!(ctx->PM2VIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->PM2NV = n;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  pm_v3alloc()    If n > 0 allocate n elements in PM3VIdx[], otherwise    */
/*                  free previously allocated memory.                       */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int pm_v3alloc(TDAContext *ctx, int n)
{
    if (ctx->PM3NV > 0) {
        free((char *)ctx->PM3VIdx);
        memrq(ctx, -ctx->PM3NV,sizeof(short));
        ctx->PM3NV = 0;
    }
    if (n > 0) {
        if (!(ctx->PM3VIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(short));
        ctx->PM3NV = n;
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

char *get_tp(TDAContext *ctx, char *tp,int *err,int opt)
{
    register int i,j,n;
    register char *p;
    int m;
    double a,a1,b,d,tmp;

    free_tp(ctx);
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
                ctx->PMTP[n] = a;
            n++;
            p = skip_dbl(ctx, p);

            if (*p == '(') {
                if (sscanf(++p,"%lg",&d) != 1 || d <= ctx->EPSI1)
                    goto GETTPErr;

                p = skip_dbl(ctx, p);
                if (*p++ != ')')
                    goto GETTPErr;
                if (sscanf(p,"%lg",&b) != 1 || fabs(b - a) <= ctx->EPSI1)
                    goto GETTPErr;
                p = skip_dbl(ctx, p);

                if (opt == 1 && b <= a + ctx->EPSI1)
                    goto GETTPErr;
                if (opt == 2 && b > 1.0)
                    goto GETTPErr;
                if (opt == 3 && (b >= a - ctx->EPSI1 || b < ctx->EPSI1))
                    goto GETTPErr;
                if (b < a)
                    d = -d;

                tmp = fabs((b - a) / d);
                if (fabs(tmp - 1.0) <= ctx->EPSI1)
                    m = 1;
                else {
                    m = (int)floor(tmp);
                    if (m <= 0 || m > 100000) {
                        printf1(ctx, "Error: exceeded limits of t/tp/x/gss/qo/qt/lon parameter.\n");
                        printf1(ctx, "or there might be an error in the increment.\n");
                        *err = 1;
                        goto GETTPErr;
                    }
                }
                if (i == 0)
                    n += m;
                else {
                    for (j = 0; j < m; ++j) {
                        ctx->PMTP[n] = ctx->PMTP[n - 1] + d;
                        n++;
                    }
                }
                a += (double)m * d;
                if (a < b - ctx->EPSI1) {
                    if (i)
                        ctx->PMTP[n] = b;
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
            if (!(ctx->PMTP = (double *)calloc((size_t)(n),sizeof(double)))) {
                *err = -2;
                goto GETTPErr;
            }
            memrq(ctx, n,sizeof(double));
            ctx->PMNTP = n;
        }
    }
    if (n != ctx->PMNTP)
        gerr_exit(ctx, 61);

    *err = 0;

    /***
    printf1(ctx, "PMNTP=%d : ",PMNTP);
    for (i = 0; i < PMNTP; ++i)
    printf1(ctx, "%lg ",PMTP[i]);
    printf1(ctx, "\n");
    ***/
    return(p);

GETTPErr:
    free_tp(ctx);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_tp()    Free memory which was used for time points.                */

void free_tp(TDAContext *ctx)
{
    if (ctx->PMNTP > 0) {
        free((char *)ctx->PMTP);
        memrq(ctx, -ctx->PMNTP,sizeof(double));
        ctx->PMNTP = 0;
    }
    ctx->PMQOFlg = ctx->PMQTFlg = 0;
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

char *get_tp1(TDAContext *ctx, char *tp,int *err,int opt)
{
    register int i,j,n;
    register char *p;
    int m;
    double a,a1,b,d,tmp;

    free_tp1(ctx);
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
                ctx->PMTP1[n] = a;
            n++;
            p = skip_dbl(ctx, p);

            if (*p == '(') {
                if (sscanf(++p,"%lg",&d) != 1 || d <= ctx->EPSI1)
                    goto GETTPErr;

                p = skip_dbl(ctx, p);
                if (*p++ != ')')
                    goto GETTPErr;
                if (sscanf(p,"%lg",&b) != 1 || fabs(b - a) <= ctx->EPSI1)
                    goto GETTPErr;
                p = skip_dbl(ctx, p);

                if (opt == 1 && b <= a + ctx->EPSI1)
                    goto GETTPErr;
                if (opt == 2 && b > 1.0)
                    goto GETTPErr;
                if (opt == 3 && (b >= a - ctx->EPSI1 || b < ctx->EPSI1))
                    goto GETTPErr;
                if (b < a)
                    d = -d;

                tmp = fabs((b - a) / d);
                if (fabs(tmp - 1.0) <= ctx->EPSI1)
                    m = 1;
                else {
                    m = (int)floor(tmp);
                    if (m <= 0 || m > 100000) {
                        printf1(ctx, "Error: exceeded limits of t/tp/x/gss/qo/qt/lat parameter.\n");
                        printf1(ctx, "or there might be an error in the increment.\n");
                        *err = 1;
                        goto GETTPErr;
                    }
                }
                if (i == 0)
                    n += m;
                else {
                    for (j = 0; j < m; ++j) {
                        ctx->PMTP1[n] = ctx->PMTP1[n - 1] + d;
                        n++;
                    }
                }
                a += (double)m * d;
                if (a < b - ctx->EPSI1) {
                    if (i)
                        ctx->PMTP1[n] = b;
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
            if (!(ctx->PMTP1 = (double *)calloc((size_t)(n),sizeof(double)))) {
                *err = -2;
                goto GETTPErr;
            }
            memrq(ctx, n,sizeof(double));
            ctx->PMNTP1 = n;
        }
    }
    if (n != ctx->PMNTP1)
        gerr_exit(ctx, 61);

    *err = 0;

    /***
    printf1(ctx, "PMNTP1=%d : ",PMNTP1);
    for (i = 0; i < PMNTP1; ++i)
    printf1(ctx, "%lg ",PMTP1[i]);
    printf1(ctx, "\n");
    ***/
    return(p);

GETTPErr:
    free_tp1(ctx);
    return(p);
}

/* -###-------------------------------------------------------------------- */
/*  free_tp1()   Free memory which was used for time points.                */

void free_tp1(TDAContext *ctx)
{
    if (ctx->PMNTP1 > 0) {
        free((char *)ctx->PMTP1);
        memrq(ctx, -ctx->PMNTP1,sizeof(double));
        ctx->PMNTP1 = 0;
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

char *get_cn(TDAContext *ctx, char *tp,int *err)
{
    register int i,n;
    register char *p;
    int m;

    free_cn(ctx);
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
                ctx->PMCN[n] = m;
            n++;

            p = skip_int(ctx, p);
            if (*p != ',')
                break;
            p++;
        }
        if (n == 0)
            goto GETCNErr;

        if (i == 0) {
            if (!(ctx->PMCN = (int *)calloc((size_t)(n),sizeof(int)))) {
                *err = -2;
                return(p);
            }
            memrq(ctx, n,sizeof(int));
            ctx->PMNCN = n;
        }
    }
    if (n != ctx->PMNCN)
        gerr_exit(ctx, 61);

    *err = 0;
    return(p);

GETCNErr:
    free_cn(ctx);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_cn()    Free memory which was used for PMCN[]                      */

void free_cn(TDAContext *ctx)
{
    if (ctx->PMNCN > 0) {
        free((char *)ctx->PMCN);
        memrq(ctx, -ctx->PMNCN,sizeof(int));
        ctx->PMNCN = 0;
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

char *get_flags(TDAContext *ctx, char *tp,int *err,int opt)
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
            ctx->PMDM[m] = 1;
        else if (opt == 2)
            ctx->PMSM[m] = 1;
        else if (opt == 3)
            ctx->PMTST[m] = 1;
        else if (opt == 5)
            ctx->PMREL[m] = 1;
        else if (opt == 4 && nc < 10)
            ctx->PMSK[++nc] = m;

        n++;

        p = skip_int(ctx, p);
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

int pmps_alloc(TDAContext *ctx, int n,int m)
{
    register int i;

    if (n >= 0) {
        if (!(ctx->PMPS[n] = (short *)calloc((size_t)(m),sizeof(short))))
            return(-2);
        memrq(ctx, m,sizeof(short));
        ctx->PMPSN[n] = m;
        ctx->PMNPS++;
        return(0);
    }
    for (i = 0; i < PMPSMax; ++i) {
        m = ctx->PMPSN[i];
        if (m > 0) {
            free((char *)ctx->PMPS[i]);
            memrq(ctx, -m,sizeof(short));
            ctx->PMPSN[i] = 0;
        }
    }
    ctx->PMNPS = 0;
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

char *get_pattern(TDAContext *ctx, char *tp,int *err)
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
        if (pmps_alloc(ctx, n,m)) {          /* insufficient memory */
            *err = -2;
            return(p);
        }
        ctx->PMPSN[n] = m;

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
                p = skip_int(ctx, p) - 1;
            else
                return(p);

            ctx->PMPS[n][i++] = (short)k;

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

char *get_fn(TDAContext *ctx, char *p,char *fname,int *err)
{
    register char c,*q;

    *err = -1;
    if (*p++ != '=')
        return(p);

    q = skip_com(ctx, p);
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

char *get_var(TDAContext *ctx, char *p,int *err)
{
    register int i,j;
    int n,w,na,nb,ia,ib;

    *err = 1;
    p = get_nvia(ctx, p,&n,1,&nb);
    if (n <= 0)
        return(p);

    na = n - nb;

    if (na > 0) {
        if (pm_valloc(ctx, na))
            return(p);
    }
    if (nb > 0) {
        if (pm_zalloc(ctx, nb))
            return(p);
    }
    ia = ib = w = 0;
    for (i = 0; i < n; ++i) {
        j = ctx->VLVIdx[i];
        if (j >= 0)
            ctx->PMVIdx[ia++] = (short)(j);
        else
            ctx->PMZIdx[ib++] = (short)(-j - 1);
        /*  j is NEGATIVE for a z-variable -- that is what the branch
            above encodes -- so VTyp[j] indexed before the start of the
            array.
            The test asks whether the variable is type 5, which only a
            real variable can be, so it belongs inside the j >= 0 case.  */
        if (j >= 0 && ctx->VTyp[j] == 5)
            w = 1;
    }
    if (w)
        p_warn(ctx, -3,1);

    *err = 0;
    alloc_vl(ctx, 0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var1(p,err)     p points to a variable list. get variables and      */
/*                      put into PM1VIdx. Count number in PM1NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var1(TDAContext *ctx, char *p,int *err)
{
    register int i;
    int n,w,nb;

    *err = 1;
    p = get_nvia(ctx, p,&n,1,&nb);

    if (n <= 0)
        return(p);
    if (nb) {
        p_err(ctx, -42,1);
        return(p);
    }
    if (pm_v1alloc(ctx, n))
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        ctx->PM1VIdx[i] = ctx->VLVIdx[i];
        if (ctx->VTyp[ctx->VLVIdx[i]] == 5)
            w = 1;
    }
    if (w)
        p_warn(ctx, -3,1);

    *err = 0;
    alloc_vl(ctx, 0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var2(p,err)     p points to a variable list. get variables and      */
/*                      put into PM2VIdx. Count number in PM2NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var2(TDAContext *ctx, char *p,int *err)
{
    register int i;
    int n,w,nb;

    *err = 1;
    p = get_nvia(ctx, p,&n,1,&nb);

    if (n <= 0)
        return(p);
    if (nb) {
        p_err(ctx, -42,1);
        return(p);
    }
    if (pm_v2alloc(ctx, n))
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        ctx->PM2VIdx[i] = ctx->VLVIdx[i];
        if (ctx->VTyp[ctx->VLVIdx[i]] == 5)
            w = 1;
    }
    if (w)
        p_warn(ctx, -3,1);

    *err = 0;
    alloc_vl(ctx, 0);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  get_var3(p,err)     p points to a variable list. get variables and      */
/*                      put into PM3VIdx. Count number in PM3NV.            */
/*                      return pointer to next character.                   */
/*                                                                          */
/*  err =  0 if successful, otherwise err = 1. Error messages already       */
/*  done here.                                                              */

char *get_var3(TDAContext *ctx, char *p,int *err)
{
    register int i;
    int n,w,nb;

    *err = 1;
    p = get_nvia(ctx, p,&n,1,&nb);

    if (n <= 0)
        return(p);
    if (nb) {
        p_err(ctx, -42,1);
        return(p);
    }
    if (pm_v3alloc(ctx, n))
        return(p);

    w = 0;
    for (i = 0; i < n; ++i) {
        ctx->PM3VIdx[i] = ctx->VLVIdx[i];
        if (ctx->VTyp[ctx->VLVIdx[i]] == 5)
            w = 1;
    }
    if (w)
        p_warn(ctx, -3,1);

    *err = 0;
    alloc_vl(ctx, 0);
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

char *get_varx(TDAContext *ctx, char *p,int *err)
{
    register int i,j;
    int n,w,na,nb,ia,ib,nm,r;
    char *q,mname[VNLMax + 1];
    short idx[MaxMatDef];

    q = p;
    *err = 1;
    p = get_nvia(ctx, p,&n,0,&nb);

    if (n > 0) {                /* list of variable names */
        na = n - nb;

        if (na > 0) {
            if (pm_valloc(ctx, na))
                return(p);
        }
        if (nb > 0) {
            if (pm_zalloc(ctx, nb))
                return(p);
        }
        ia = ib = w = 0;
        for (i = 0; i < n; ++i) {
            j = ctx->VLVIdx[i];
            if (j >= 0)
                ctx->PMVIdx[ia++] = (short)(j);
            else
                ctx->PMZIdx[ib++] = (short)(-j - 1);
            if (ctx->VTyp[j] == 5)
                w = 1;
        }
        if (w)
            p_warn(ctx, -3,1);
    }
    else {                      /* list of matrix names */

        r = nm = 0;
        p = q;
        while (1) {
            if ((q = get_mname(ctx, p,mname,0)) == NULL) {
                if (nm > 0)
                    p--;
                break;
            }
            if ((i = mat_getidx(ctx, mname,0)) < 0) {
                r = 1;
                break;
            }
            if (nm >= MaxMatDef) {
                r = 2;
                break;
            }
            idx[nm++] = (short)(i);
            p = q;
            if (*p != ',')
                break;
            p++;
        }
        if (r > 0 || nm == 0) {
            printf1(ctx, "Error: invalid list of variable or matrix names.\n");
            if (r == 2)
                printf1(ctx, "Error: exceeded max number of matrix names.\n");
            return(p);
        }
        if (pm_valloc(ctx, nm))
            return(p);

        for (i = 0; i < nm; ++i)
            ctx->PMVIdx[i] = (short)idx[i];
        ctx->PMNVTyp = 1;
    }
    *err = 0;
    alloc_vl(ctx, 0);
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

char *get_xp(TDAContext *ctx, char *s,int *err)
{
    register int i,n;
    register char *p;
    double a,b,x;

    free_xp(ctx);
    *err = -1;
    p = s;
    if (!*p)
        goto GETXErr;

    for (i = 0; i < 2; ++i) {

        p = s;
        n = 0;
        while (*p) {
            if (sscanf(p,"%lg[%lg,%lg]",&x,&a,&b) == 3) {
                p = skip_dbl(ctx, p);
                p = skip_dbl(ctx, p + 1);
                p = skip_dbl(ctx, p + 1) + 1;
                if (i) {
                    ctx->PMXX[n] = x;
                    ctx->PMXA[n] = a;
                    ctx->PMXB[n] = b;
                }
            }
            else if (sscanf(p,"%lg",&x) == 1) {
                p = skip_dbl(ctx, p);
                if (i) {
                    ctx->PMXX[n] = x;
                    ctx->PMXA[n] = x - 1.0;
                    ctx->PMXB[n] = x + 1.0;
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
            if (!(ctx->PMXX = (double *)calloc((size_t)(n),sizeof(double)))) {
                *err = -2;
                return(p);
            }
            if (!(ctx->PMXA = (double *)calloc((size_t)(n),sizeof(double)))) {
                free((char *)ctx->PMXX);
                *err = -2;
                return(p);
            }
            if (!(ctx->PMXB = (double *)calloc((size_t)(n),sizeof(double)))) {
                free((char *)ctx->PMXX);
                free((char *)ctx->PMXA);
                *err = -2;
                return(p);
            }
            memrq(ctx, 3 * n,sizeof(double));
            ctx->PMNX = n;
        }
    }
    if (n != ctx->PMNX)
        gerr_exit(ctx, 61);

    *err = 0;

    /***
    printf1(ctx, "PMNX=%d : ",PMNX);
    for (i = 0; i < PMNX; ++i)
        printf1(ctx, "%lg %lg %lg\n",PMXX[i],PMXA[i],PMXB[i]);
    ***/
    return(p);

GETXErr:
    free_xp(ctx);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_xp()    Free memory which was used for xp parameter.               */

void free_xp(TDAContext *ctx)
{
    if (ctx->PMNX > 0) {
        free((char *)ctx->PMXX);
        free((char *)ctx->PMXA);
        free((char *)ctx->PMXB);
        memrq(ctx, -3 * ctx->PMNX,sizeof(double));
        ctx->PMNX = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  get_file(typ,p,err)                                                     */
/*                                                                          */
/*  err =  0 if successful                                                  */
/*        -1 if syntax error                                                */
/*        -2 if insuff memory                                               */
/*         1 if file cannot be opened.                                      */

char *get_file(TDAContext *ctx, int typ,char *p,int *err)
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
        q = skip_blev(ctx, p);
        if (*--q != ')')
            return(p);
        *q = '\0';
        n = (int)(strlen(++p));
        if (n < 1)
            return(p);

        if (rhstr_alloc(ctx, n + 1)) {
            *err = -2;
            return(p);
        }
        strcpy(ctx->PMRHSTR,p);
        *q = ')';
        p = q + 1;
    }
    else if (typ == 12 && *p == '(') {         /* dvar file */

        q = p + 1;
        np = 0;
        while (*q) {
            if (sscanf(q,"fn=%d",&n) == 1 && n > 0) {
                ctx->PMDVARFN = n;
                q = skip_int(ctx, q + 3);
            }
            else if (sscanf(q,"pn=%d",&n) == 1 && n > 0) {
                ctx->PMDVARPN = n;
                q = skip_int(ctx, q + 3);
            }
            else if (sscanf(q,"p=%d",&n) == 1 && n > 0) {
                ctx->PMDVARP = n;
                q = skip_int(ctx, q + 2);
            }
            else if (!strncmp(q,"p(",2) || !strncmp(q,"pn(",3)) {
                np++;
                if (!strncmp(q,"p(",2))
                    q += 2;
                else
                    q += 3;
                l = get_vnlen(ctx, q);
                if (l < 1)
                    return(p);
                q += l;
                if (*q++ != ')' || *q++ != '=')
                    return(p);
                if (sscanf(q,"%d",&n) != 1 || n < 1)
                    return(p);
                q = skip_int(ctx, q);
            }
            else
                return(p);

            if (*q != ',')
                break;
            q++;
        }
        if (*q++ != ')')
            return(p);

        if (ctx->PMDVARP > 0 && ctx->PMDVARPN > 0)
            return(p);

        if (np > 0 && (ctx->PMDVARP > 0 || ctx->PMDVARPN > 0))
            return(p);

        if (np > 0 && ctx->PMDVARP == 0 && ctx->PMDVARPN == 0) {
/* ### */
            ctx->PMDVARVN = np;

            if (!(ctx->PMDVARVNP = (short *)calloc((size_t)(np),sizeof(int)))) {
                *err = -2;
                return(p);
            }
            memrq(ctx, np,sizeof(int));
            ctx->PMDVARVNPA = np;

            if (!(ctx->PMDVARVNL = (short *)calloc((size_t)(np),sizeof(int)))) {
                *err = -2;
                return(p);
            }
            memrq(ctx, np,sizeof(int));
            ctx->PMDVARVNLA = np;

            if (!(ctx->PMDVARVName = (char **)calloc((size_t)(np),sizeof(char *)))) {
                *err = -2;
                return(p);
            }
            memrq(ctx, np,sizeof(char *));
            ctx->PMDVARVNameA = np;

            ii = 0;
            q = p + 1;

            while (*q) {
                if (sscanf(q,"fn=%d",&n) == 1 && n > 0) {
                    q = skip_int(ctx, q + 3);
                }
                else if (sscanf(q,"pn=%d",&n) == 1 && n > 0) {
                    q = skip_int(ctx, q + 3);
                }
                else if (sscanf(q,"p=%d",&n) == 1 && n > 0) {
                    q = skip_int(ctx, q + 2);
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
                    l = get_vnlen(ctx, q);
                    if (l < 1)
                        return(p);

                    if (!(ctx->PMDVARVName[ii] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                        *err = -2;
                        return(p);
                    }
                    memrq(ctx, l + 1,sizeof(char));
                    ctx->PMDVARVNP[ii] = (short)(np1);
                    strncpy(ctx->PMDVARVName[ii],q,(size_t)(l));
                    q += l;
                    if (*q++ != ')' || *q++ != '=')
                        return(p);
                    if (sscanf(q,"%d",&n) != 1 || n < 1)
                        return(p);
                    q = skip_int(ctx, q);
                    ctx->PMDVARVNL[ii] = (short)(n);
                    ii++;
                }
                else
                    return(p);

                if (*q != ',')
                    break;
                q++;
            }
            if (*q++ != ')' || ctx->PMDVARVN != ii)
                return(p);
        }
        p = q;
    }
    else if (typ == 13 && *p == '(') {         /* arcd file */
        p++;
        while (*p) {
            if (!strncmp(p,"zoo=",4)) {
                p += 4;
                q = ctx->PMARCFZOOF;
                for (i = 0; i < FNMaxLen; ++i) {
                    if (!*p || *p == ',' || *p == ')')
                        break;
                    *q++ = *p++;
                }
                *q = '\0';
                ctx->PMARCFZOO = 1;
            }
            else if (!strncmp(p,"vdf=",4)) {
                p += 4;
                q = ctx->PMARCFVDFF;
                for (i = 0; i < FNMaxLen; ++i) {
                    if (!*p || *p == ',' || *p == ')')
                        break;
                    *q++ = *p++;
                }
                *q = '\0';
                ctx->PMARCFVDF = 1;
            }
            if (*p != ',')
                break;
            p++;
        }
        if (*p++ != ')')
            return(p);
    }
    p = get_fn(ctx, p,fname,err);
    if (*err)
        return(p);

    *err = 0;
    switch (typ) {
        case 1: if (ctx->PMF1Def) {
                    fclose(ctx->PMF1d);
                    ctx->PMF1Def = 0;
                }
                if ((aflag == 0 && !(ctx->PMF1d = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMF1d = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMF1Def = nn;
                    strcpy(ctx->PMF1dName,fname);
                }
                break;

        case 2: if (ctx->PMCovFDef) {
                    fclose(ctx->PMCovFd);
                    ctx->PMCovFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMCovFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMCovFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMCovFDef = 1;
                    strcpy(ctx->PMCovFName,fname);
                }
                break;

        case 3: if (ctx->PMResFDef) {
                    fclose(ctx->PMResFd);
                    ctx->PMResFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMResFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMResFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMResFDef = 1;
                    strcpy(ctx->PMResFName,fname);
                }
                break;

        case 4: if (ctx->PMPPFDef) {
                    fclose(ctx->PMPPFd);
                    ctx->PMPPFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMPPFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMPPFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMPPFDef = 1;
                    strcpy(ctx->PMPPFName,fname);
                }
                break;

        case 5: if (ctx->PMProtFDef) {
                    fclose(ctx->PMProtFd);
                    ctx->PMProtFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMProtFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMProtFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMProtFDef = nn;
                    strcpy(ctx->PMProtFName,fname);
                }
                break;

        case 6: if (ctx->PMPCFDef) {
                    fclose(ctx->PMPCFd);
                    ctx->PMPCFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMPCFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMPCFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMPCFDef = nn;
                    strcpy(ctx->PMPCFName,fname);
                }
                break;

        case 7: if (ctx->PMPFNFDef) {
                    fclose(ctx->PMPFNFd);
                    ctx->PMPFNFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMPFNFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMPFNFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMPFNFDef = 2;
                    strcpy(ctx->PMPFNFName,fname);
                }
                break;

        case 8: if (ctx->PMTabFDef) {
                    fclose(ctx->PMTabFd);
                    ctx->PMTabFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMTabFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMTabFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMTabFDef = 1;
                    strcpy(ctx->PMTabFName,fname);
                }
                break;

        case 9: if (ctx->PMTab1FDef) {
                    fclose(ctx->PMTab1Fd);
                    ctx->PMTab1FDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMTab1Fd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMTab1Fd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMTab1FDef = 1;
                    strcpy(ctx->PMTab1FName,fname);
                }
                break;

       case 10: if (ctx->PMTDAFDef) {
                    fclose(ctx->PMTDAFd);
                    ctx->PMTDAFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMTDAFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMTDAFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMTDAFDef = 1;
                    strcpy(ctx->PMTDAFName,fname);
                }
                break;

       case 11: if (ctx->PMSPSSFDef) {
                    fclose(ctx->PMSPSSFd);
                    ctx->PMSPSSFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMSPSSFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMSPSSFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMSPSSFDef = 1;
                    strcpy(ctx->PMSPSSFName,fname);
                }
                break;

       case 12: if (ctx->PMDVARFDef) {
                    fclose(ctx->PMDVARFd);
                    ctx->PMDVARFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMDVARFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMDVARFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMDVARFDef = 1;
                    strcpy(ctx->PMDVARFName,fname);
                }
                break;

       case 13: if (ctx->PMARCFDef) {
                    fclose(ctx->PMARCFd);
                    ctx->PMARCFDef = 0;
                }
                if ((aflag == 0 && !(ctx->PMARCFd = fopen(fname,OPEN_WR))) ||
                    (aflag == 1 && !(ctx->PMARCFd = fopen(fname,OPEN_AP))))
                    *err = 1;
                else {
                    ctx->PMARCFDef = 1;
                    strcpy(ctx->PMARCFName,fname);
                }
                break;

       case 14: if (ctx->PMIF1Def) {
                    fclose(ctx->PMIF1d);
                    ctx->PMIF1Def = 0;
                }
                if (!(ctx->PMIF1d = fopen(fname,OPEN_RD)))
                    *err = 1;
                else {
                    ctx->PMIF1Def = 1;
                    strcpy(ctx->PMIF1Name,fname);
                }
                break;

       case 15: if (ctx->PMIF2Def) {
                    fclose(ctx->PMIF2d);
                    ctx->PMIF2Def = 0;
                }
                if (!(ctx->PMIF2d = fopen(fname,OPEN_RD)))
                    *err = 1;
                else {
                    ctx->PMIF2Def = 1;
                    strcpy(ctx->PMIF2Name,fname);
                }
                break;

        default: break;
    }
    if (*err)
        printf1(ctx, "Error: can't open: %s\n",fname);
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

char *get_box(TDAContext *ctx, char *s,int *err)
{
    register int i,n;
    register char *p;
    double l,u;

    free_box(ctx);
    *err = -1;
    p = s;
    if (!*p)
        goto GETBOXErr;

    for (i = 0; i < 2; ++i) {

        p = s;
        n = 0;
        while (*p) {
            if (sscanf(p,"%lg,%lg",&l,&u) == 2 && l < u) {
                p = skip_dbl(ctx, p);
                p = skip_dbl(ctx, p + 1);
                if (i) {
                    ctx->PMBoxL[n] = l;
                    ctx->PMBoxU[n] = u;
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
            if (!(ctx->PMBoxL = (double *)calloc((size_t)(n),sizeof(double)))) {
                *err = -2;
                return(p);
            }
            if (!(ctx->PMBoxU = (double *)calloc((size_t)(n),sizeof(double)))) {
                free((char *)ctx->PMBoxL);
                *err = -2;
                return(p);
            }
            memrq(ctx, 2 * n,sizeof(double));
            ctx->PMBoxN = n;
        }
    }
    if (n != ctx->PMBoxN)
        gerr_exit(ctx, 61);

    *err = 0;
    return(p);

GETBOXErr:
    free_box(ctx);
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_box()    Free memory which was used for box parameter.             */

void free_box(TDAContext *ctx)
{
    if (ctx->PMBoxN > 0) {
        free((char *)ctx->PMBoxL);
        free((char *)ctx->PMBoxU);
        memrq(ctx, -2 * ctx->PMBoxN,sizeof(double));
        ctx->PMBoxN = 0;
    }
}

/* ------------------------------------------------------------------------ */
/*  get_mf(p,err)   Get list of merge files beginning at p, return pointer  */
/*                  to next character. err = 0 if OK, -1 if error.          */

char *get_mf(TDAContext *ctx, char *p,int *err)
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
        if (ctx->MFN >= MFMAX) {
            printf1(ctx, "Error: exceeded max number of merge files.\n");
            return(p);
        }
        if (get_fname(ctx, p,ctx->MFFNAME[ctx->MFN]) == 0)
            return(p);
        if (sscanf(q,"%d,%d,%d,%d]",&ctx->MFI1[ctx->MFN],&ctx->MFI2[ctx->MFN],&ctx->MFJ1[ctx->MFN],&ctx->MFJ2[ctx->MFN]) != 4)
            return(p);
        while (*q && *q != ']')
            q++;
        if (*q++ != ']')
            return(p);

        if (!(ctx->MFFD[ctx->MFN] = fopen(ctx->MFFNAME[ctx->MFN],OPEN_RD))) {
            printf1(ctx, "Error: can't open: %s\n",ctx->MFFNAME[ctx->MFN]);
            return(p);
        }
        ctx->MFN++;
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

char *get_mp(TDAContext *ctx, char *p,int l,int *mdef,char *name,int *err)
{
    register int i;
    register char *q,*r;

    *err = -1;
    r = p + l;

    if (*r == '(') {
        if (!strncmp(p,"mpgrad(v=",9))
            r = get_var2(ctx, r + 3,err);
        else if (!strncmp(p,"mpres(v=",8))
            r = get_var3(ctx, r + 3,err);
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
    if (mat_ncheck(ctx, name,0) || (*r != ',' && *r != ')')) { /* no valid matrix name */
        printf1(ctx, "Error: %s\nNo valid matrix name.\n",p);
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

char *get_mod(TDAContext *ctx, char *p,int *err)
{
    int n;
    register char *q;

    *err = -1;

    if (ctx->PMModN >= MaxLL) {
        printf1(ctx, "Error: exceeded max number of mod strings.\n");
        return(p);
    }
    q = p;
    n = 0;
    while (*q && *q != ',') {
        n++;
        q++;
    }
    if (!(ctx->PMModS[ctx->PMModN] = (char *)calloc((size_t)(n + 1),sizeof(char)))) {
        *err = -2;
        return(p);
    }
    memrq(ctx, n + 1,sizeof(char));
    ctx->PMModSA[ctx->PMModN] = n;

    q = ctx->PMModS[ctx->PMModN];
    while (*p && *p != ',')
        *q++ = *p++;
    *q = '\0';
    ctx->PMModN++;
    *err = 0;
    return(p);
}

/* ------------------------------------------------------------------------ */
/*  free_mod()  free PMModS                                                 */

void free_mod(TDAContext *ctx)
{
    register int i;

    for (i = 0; i < ctx->PMModN; ++i) {
        if (ctx->PMModSA[i] > 0) {
            free(ctx->PMModS[i]);
            memrq(ctx, -ctx->PMModSA[i] - 1,sizeof(char));
            ctx->PMModSA[i] = 0;
        }
    }
    ctx->PMModN = 0;
}

/* ------------------------------------------------------------------------ */
/*  get_fmt(k,p,&err)   get print formats                                   */

char *get_fmt(TDAContext *ctx, int k,char *p,int *err)
{
    register int i;
    register char *q;
    int n,m;

    *err = -1;
    n = 0;
    q = p;
    while (*q) {
        if (sscanf(q,"%d",&m) == 1) {
            q = skip_int(ctx, q);
            if (*q == '.') {
                if (sscanf(++q,"%d",&m) == 1)
                    q = skip_int(ctx, q);
            }
            n++;
        }
        if (!*q || *q != ',')
            break;
        q++;
    }
    if (n == 0)
        return(p);

    if (get_fmt_alloc(ctx, k,n))
        return(p);

    n = 0;
    while (*p) {
        if (sscanf(p,"%d",&m) == 1) {
            p = skip_int(ctx, p);
            ctx->PMXFmt1[k][n] = m;
            ctx->PMXFmt2[k][n] = 0;
            if (*p == '.') {
                if (sscanf(++p,"%d",&m) == 1) {
                    ctx->PMXFmt2[k][n] = m;
                    p = skip_int(ctx, p);
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
        makefmt(ctx, &ctx->PMXFmt1[k][i],&ctx->PMXFmt2[k][i],ctx->PMXFmtS[k] + i * ctx->PMXFmtLen,(size_t)ctx->PMXFmtLen,0,ctx->SEPC,0);

    ctx->PMXFmtN[k] = n;
    *err = 0;
    return(p);
}


/* ------------------------------------------------------------------------ */
/*  get_fmt_alloc(k,n)  allocate mem for n format strings.                  */
/*                      if n = 0 free previously allocated memory.          */
/*                      return 0 if OK, -1 if error.                        */

int get_fmt_alloc(TDAContext *ctx, int k,int n)
{
    if (ctx->PMXFmt1A[k] > 0) {
        free((char *)ctx->PMXFmt1[k]);
        memrq(ctx, -ctx->PMXFmt1A[k],sizeof(int));
        ctx->PMXFmt1A[k] = 0;
    }
    if (ctx->PMXFmt2A[k] > 0) {
        free((char *)ctx->PMXFmt2[k]);
        memrq(ctx, -ctx->PMXFmt2A[k],sizeof(int));
        ctx->PMXFmt2A[k] = 0;
    }
    if (ctx->PMXFmtSA[k] > 0) {
        free((char *)ctx->PMXFmtS[k]);
        memrq(ctx, -ctx->PMXFmtSA[k],sizeof(char));
        ctx->PMXFmtSA[k] = 0;
    }
    ctx->PMXFmtN[k] = 0;

    if (n > 0) {
        if (!(ctx->PMXFmt1[k] = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->PMXFmt1A[k] = n;

        if (!(ctx->PMXFmt2[k] = (int *)calloc((size_t)(n),sizeof(int)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n,sizeof(int));
        ctx->PMXFmt2A[k] = n;

        if (!(ctx->PMXFmtS[k] = (char *)calloc((size_t)(n) * (size_t)(ctx->PMXFmtLen),sizeof(char)))) {
            p_err(ctx, -2,1);
            return(-1);
        }
        memrq(ctx, n * ctx->PMXFmtLen,sizeof(char));
        ctx->PMXFmtSA[k] = n * ctx->PMXFmtLen;
        ctx->PMXFmtN[k] = n;
    }
    return(0);
}
