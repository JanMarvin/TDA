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
#include "tda_context.h"

/*  functions in t_gdat.c */

int new_var(TDAContext *ctx);
int check_vsel(TDAContext *ctx);               
void prn_dfd(TDAContext *ctx, int nv,short *vidx);
int check_matchv(TDAContext *ctx, int iv,int idx);
void prn_mval(TDAContext *ctx);
int check_nv(TDAContext *ctx);
void prn_nve(TDAContext *ctx, char *s);
void prn_nve1(TDAContext *ctx, char *s);
void make_vfmt(TDAContext *ctx, int nidx);
void make_afmt(TDAContext *ctx); 
int save_isel(TDAContext *ctx, int opt,int idx);
int check_iref(TDAContext *ctx, int cnt,int *ptyp,int idx);
int save_vsel(TDAContext *ctx, int opt,int idx);
int check_vref(TDAContext *ctx, int cnt,int *ptyp,int idx);
int save_bsel(TDAContext *ctx, int opt);
int save_break(TDAContext *ctx, int opt,int idx);
int check_vcj(TDAContext *ctx, int idx,int opt);
int check_ffmt(TDAContext *ctx, int opt);
double get_data(TDAContext *ctx, int j,int i);
void put_data(TDAContext *ctx, double x,int j,int i);
void put_str(TDAContext *ctx, char *buf,int blen,int j,int i,int bflag);
void get_str(TDAContext *ctx, char *buf,int j,int i);
void clear_str(TDAContext *ctx, char *buf,int blen,int j);
double dscan(TDAContext *ctx, char *p,int len,int *mval);
int sep_char(TDAContext *ctx, char c);
char *skip_sep(TDAContext *ctx, char *p);
char *skip_dval(TDAContext *ctx, char *p);
int gen_dm(TDAContext *ctx, int idx,int nvar,int sflag,int sdflag);
void prn_gdf(TDAContext *ctx, int n);
int fnd_match(TDAContext *ctx, double *xm,int *n);
int put_t4var(TDAContext *ctx, int idx,int r1,int r2,int bn);
void cpy_dmrow(TDAContext *ctx, int ir,int ij,int idx);
int tst_bsel(TDAContext *ctx, int r1,int r2,int idx);
int tst_vsel(TDAContext *ctx, int i,int rec);
int tst_break(TDAContext *ctx, int r1,int r2);
int read_dfi(TDAContext *ctx, int opt,int sflag);
int read_df(TDAContext *ctx, int sflag,int idx);
int skip_df(TDAContext *ctx, int n,int opt);         
int get_nsdxy(TDAContext *ctx, char *buf,double *x,double *y);
int check_mval(TDAContext *ctx, int mval);

int sdnvar(TDAContext *ctx);
int check_nvsd(TDAContext *ctx);
int check_sd(TDAContext *ctx, int opt);
void sdnvar_close(TDAContext *ctx);
int sdnvar_alloc(TDAContext *ctx, int opt,int n);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */



/*--------------------------------------------------------------------------*/
                        /* depends on the actual tsel command.              */
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/




#define DFILMax 100     /* max number of data files                         */






/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/* ------------------------------------------------------------------------ */
/*  new_var.    Command: nvar(...)                                          */
/*              Create new variables. If a data matrix already exists,      */
/*              add the new variables. Otherwise create a new data matrix   */
/*              implying that tsel = off and no case weights.               */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int new_var(TDAContext *ctx)
{
    register int i,j,k;
    int m,n,err,dmdef,nb,w,d,sflag;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Creating new variables. ");
    prn_mem(ctx);
    tsel_off(ctx, 1);        /* turn off temporary case selection */

    dmdef = ctx->DMDef;
    ctx->DFILN = ctx->MatchNV = ctx->VFmtFlg = ctx->NVNum = 0;
    ctx->DRecLen = 0;            /* variable record length */
    ctx->NMRec = 1;              /* logical = physical record */
    ctx->DBLKVar = -1;           /* no block mode */
    ctx->AVDFN = -1;             /* no archive file */
    ctx->ARCDICFlg = 0;          /* set by arcdic option */
    ctx->RDMNInit = 0;           /* init random number generator */
    ctx->ISEPC = '\0';           /* separation character */
    ctx->VTNOC = 0;
    ctx->BSIZE = 0;

    ctx->EDVALSn = ctx->EDVALOrg = ctx->EDVALDes = 0;
    ctx->EDVALTs = ctx->EDVALTf = 0.0;
         
    ctx->NVIdx = get_nidx(ctx);     /* index to first new variable */

    /* set default missing value codes */

    ctx->MStarVal  = -1.0; 
    ctx->MPntVal   = -1.0; 
    ctx->MBlnkVal  = -1.0; 
    ctx->MGenVal   = -1.0; 
    ctx->MMatchVal = -3.0; 
    ctx->GDNRec = ctx->GDKFlg = ctx->GDFlg = ctx->MGenFlg = ctx->MStarN = ctx->MPntN = ctx->MBlnkN = ctx->MGenN = 0; 

    ctx->ISelPtr  = NULL;        /* set by isel option */
    ctx->VSelPtr  = NULL;        /* set by vsel option */
    ctx->BSelPtr  = NULL;        /* set by bsel option */
    ctx->GDFNPtr  = NULL;        /* file name by df option */
    ctx->KeepPtr  = NULL;            
    ctx->DropPtr  = NULL;            

    ctx->NVArc = ctx->NVArc1 = 0;     /* number of archive variables, counted by save_var(ctx) */

    if (check_nv(ctx))         /* check command etc */
        goto NVFin;

    if (ctx->NVNum == 0) {
        printf1(ctx, "No variables defined.\n");
        goto NVFin;
    }
    newline(ctx);      

    if (ctx->VLabelLen > 0 && ctx->VLabelLen < 8)      
        ctx->VLabelLen = 8;

    sflag = 0;
    i = ctx->NVIdx;                  /* count string variables */
    while (i >= 0) {
        if (ctx->VTyp[i] == 1 && ctx->VTypA[i] == 0)
            sflag++;
        i = ctx->VNxt[i];
    }

    if (ctx->NVArc1 > 0) {
        printf1(ctx, "Number of archive variables: %d\n",ctx->NVArc1);
        if (ctx->ARCDef == 0) {
            p_err(ctx, -9,1);
            goto NVFin;
        }
        ctx->NVArc = ctx->NVArc1;
        printf1(ctx, "Searching for these variables in archive: %s\n",ctx->ZOONam);
        if (alloc_avar(ctx, ctx->NVIdx,ctx->NVArc,1)) {
            p_err(ctx, -2,1);
            goto NVFin;
        }

        /* check whether all variables belong to same data file */

        ctx->AVDFN = check_avar(ctx, ctx->NVArc,ctx->ARCDICFlg);
        if (ctx->AVDFN < 0)
            goto NVFin;

        make_afmt(ctx);        /* make new print formats for archive variables */
        newline(ctx);
    }
    if (ctx->VFmtFlg)
        make_vfmt(ctx, ctx->NVIdx);   /* make new print formats if requested with fmt */

    prn_var(ctx, ctx->NVIdx);         /* print list of new variables */

    if (ctx->DMDef == 0) {
        if (ctx->GDFlg == 0)
            ctx->NOCMaxA = ctx->VTNOC;
        else {
            if (ctx->BSIZE > 0)
                ctx->NOCMaxA = ctx->BSIZE;
            else if (ctx->VTNOC > 0)
                ctx->NOCMaxA = ctx->VTNOC;
            else
                ctx->NOCMaxA = ctx->VTNOC = NOCDef;
        }
        if (ctx->NOCMaxA <= 0) {
            if (ctx->AVDFN >= 0 && ctx->GDFlg == 0)
                ctx->NOCMaxA = ctx->ZAFNRec[ctx->AVDFN];
            else
                ctx->NOCMaxA = NOCDef;
        }
         
        if (ctx->GDFlg == 0) {
            printf1(ctx, "\nCreating a new data matrix.\n");
            printf1(ctx, "Maximum number of cases: %d\n",ctx->NOCMaxA);
        }
        else {
            printf1(ctx, "\nData will be directly written to: %s\n",ctx->GDFNPtr);
            printf1(ctx, "Will not add variables to data matrix.\n");

            if (!(ctx->GDFd = fopen(ctx->GDFNPtr,OPEN_WR))) {
                printf1(ctx, "Error: can't open: %s\n",ctx->GDFNPtr);
                goto NVFin;
            }
            ctx->GDFlg = 2;
   
            /* create list of variables */

            if (ctx->KeepPtr != NULL) {
                printf1(ctx, "Checking keep=%s\n",ctx->KeepPtr); 
                get_nvia(ctx, ctx->KeepPtr,&n,1,&nb);
                if (n < 1)
                    goto NVFin;
                if (nb) {
                    p_err(ctx, -42,1);
                    goto NVFin;
                }

                if (!(ctx->NVarIdx = (short *)calloc((size_t)(n),sizeof(short)))) {
                    p_err(ctx, -2,1);
                    goto NVFin;
                }
                memrq(ctx, n,sizeof(short));
                ctx->NVarNV = n;       
                for (i = 0; i < n; ++i)
                    ctx->NVarIdx[i] = ctx->VLVIdx[i];
            }
            else {
                if (!(ctx->NVarIdx = (short *)calloc((size_t)(ctx->NVNum),sizeof(short)))) {
                    p_err(ctx, -2,1);
                    goto NVFin;
                }
                memrq(ctx, ctx->NVNum,sizeof(short));
                ctx->NVarNV = ctx->NVNum;

                j = 0;
                i = ctx->NVIdx;          
                while (i >= 0) {
                    ctx->NVarIdx[j++] = (short)(i);
                    i = ctx->VNxt[i];
                }
                if (ctx->DropPtr != NULL) {
                    printf1(ctx, "Checking drop=%s\n",ctx->DropPtr); 
                    get_nvia(ctx, ctx->DropPtr,&n,1,&nb);
                    if (n < 1)
                        goto NVFin;
                    if (nb) {
                        p_err(ctx, -42,1);
                        goto NVFin;
                    }

                    m = 0;
                    for (i = 0; i < ctx->NVarNV; ++i) {
                        k = ctx->NVarIdx[i];
                        for (j = 0; j < n; ++j) {
                            if (k == ctx->VLVIdx[j]) {
                                ctx->NVarIdx[i] = (short)(k = -1);
                                break;
                            }
                        }
                        if (k >= 0)
                            m++;
                    }
                    if (m == 0) {
                        printf1(ctx, "Error: all variables dropped.\n");
                        goto NVFin;
                    }
                }
            }
                                    /* make print format */
            if (ctx->NVarNV > 0) {
                if (!(ctx->NVPFmtS = (char **)calloc((size_t)(ctx->NVarNV),sizeof(char *)))) {
                    p_err(ctx, -2,1);             
                    goto NVFin;
                }         
                memrq(ctx, ctx->NVarNV,sizeof(char *));
                ctx->NVPFmtSA = ctx->NVarNV;

                if (ctx->NVL0 != 0)
                    ctx->NVL0 = -1;

                for (i = 0; i < ctx->NVarNV; ++i) {  

                    if (!(ctx->NVPFmtS[i] = (char *)calloc(NVPFmtSLen,sizeof(char)))) {
                        p_err(ctx, -2,1);             
                        goto NVFin; 
                    }         
                    memrq(ctx, NVPFmtSLen,sizeof(char));
                    ctx->NVPFmtSA1++;  

                    k = ctx->NVarIdx[i];
                    w = (int)ctx->VPFmt1[k];             
                    d = (int)ctx->VPFmt2[k];             
                    if (ctx->VTyp[k] != 1)
                        makefmt(ctx, &w,&d,ctx->NVPFmtS[i],NVPFmtSLen,0,ctx->NVSEPC,ctx->NVL0);
                }
            }
        }
    }
    else if (ctx->GDFlg) {
        printf1(ctx, "\nError: df parameter cannot be used if a data matrix exists.\n");
        goto NVFin;
    }

    if (ctx->ISelPtr != NULL) {              /* isel option */
        if (check_vsel(ctx))
            goto NVFin;
        if (save_isel(ctx, 1,ctx->NVIdx))
            goto NVFin;
    }
    if (ctx->VSelPtr != NULL) {              /* vsel option */
        if (check_vsel(ctx))
            goto NVFin;

        if (ctx->DBLKVar >= 0)  
            printf1(ctx, "Block mode: vsel will be ignored.\n");
        else {
            if (save_vsel(ctx, 1,ctx->NVIdx))
                goto NVFin;
        }
    }
    if (ctx->BRKPtr != NULL) {               /* break option */
        if (check_vsel(ctx))
            goto NVFin;
        if (ctx->DBLKVar >= 0)  
            printf1(ctx, "Block mode: break will be ignored.\n");
        else {
            if (save_break(ctx, 1,ctx->NVIdx))
                goto NVFin;
        }
    }
    if (check_vcj(ctx, ctx->NVIdx,1))             /* check cj references */
        goto NVFin;

    if (ctx->VCJMax > 0 || sflag) {          /* need a data file */

        if (ctx->NVArc > 0) {
            printf1(ctx, "Error: referring to a data file with c terms or string\n");
            printf1(ctx, "variables is not compatible with using archive variables.\n");
            goto NVFin;
        }

        if (ctx->DFILN > 0) {
            printf1(ctx, "\nUsing data file(s): %s",ctx->DFILFN[0]);
            for (i = 1; i < ctx->DFILN; ++i)
                printf1(ctx, ",%s",ctx->DFILFN[i]);
            printf1(ctx, "\n");
            if (ctx->DRecLen > 0)
                printf1(ctx, "Fixed record length: %d\n",ctx->DRecLen);
            if (ctx->NMRec > 1)
                printf1(ctx, "Logical records consist of %d physical records.\n",ctx->NMRec);
        }
        else {
            printf1(ctx, "\nError: need a data file.\n");
            goto NVFin;
        }
        if (ctx->FFMTCnt > 0) {
            if (check_ffmt(ctx, 1))
                goto NVFin;
        }
        else {
            printf1(ctx, "Free format. Separation character(s): ");
            if (ctx->ISEPC != '\0')
                printf1(ctx, "%04x [hex]\n",ctx->ISEPC);
            else 
                printf1(ctx, "default.\n");
        }
        if (ctx->VTNOC > 0)  
            printf1(ctx, "Reading maximal %d records.\n",ctx->VTNOC);

    }
    else if (ctx->DFILN > 0) {
        printf1(ctx, "\nWarning: dfile parameter(s) will be ignored.\n");
        ctx->DFILN = 0;
    }
    if (ctx->DMDef) {            /* if a data matrix already exists */
   
        printf1(ctx, "\nNew variables will be added to existing data matrix.\n");

        if (ctx->MatchNV > 0) {
            if (ctx->NVArc == 0 && ctx->DFILN == 0) {
                printf1(ctx, "Always trivial matching without data files.\n");
                ctx->MatchNV = 0;
            }
            else {
                printf1(ctx, "\nMatching:\nNew variable    Existing variable\n");
                prnchar(ctx, '-',33,1);
                     
                for (i = 0; i < ctx->MatchNV; ++i) {
                    j = ctx->MatchVA[i];
                    printf1(ctx, "%s",ctx->VName[j]);
                    prnchar(ctx, ' ',(int)(16 - strlen(ctx->VName[j])),0);
                    printf1(ctx, "%s\n",ctx->VName[ctx->MatchVB[i]]);
                    if (ctx->VTyp[j] != 2 && ctx->VTyp[j] != 3) {
                        printf1(ctx, "\nError: %s is not of type 2 or 3.\n",ctx->VName[j]);
                        goto NVFin;
                    }
                    if (check_matchv(ctx, ctx->MatchVB[i],ctx->NVIdx)) {
                        printf1(ctx, "\nError: %s is not an existing variable.\n",ctx->VName[ctx->MatchVB[i]]);
                        goto NVFin;
                    }
                    j = ctx->MatchVB[i];
                    if (ctx->VTyp[j] == 1) {
                        printf1(ctx, "\nError: %s must not be of type 1 (string variable).\n",ctx->VName[j]);
                        goto NVFin;
                    }
                }
                printf1(ctx, "\n");
                if (ctx->DBLKVar >= 0) {
                    printf1(ctx, "Error: dblock and match not compatible.\n");
                    goto NVFin;
                }
                if (ctx->BRKPtr != NULL) {
                    printf1(ctx, "Break will be ignored.\n");
                    save_break(ctx, 0,0);
                }
                if (vsort(ctx, ctx->MatchNV,ctx->MatchVB,1,1,1))       /* sort */
                    goto NVFin;
            }
        }
        else
            printf1(ctx, "Trivial matching.\n");
    }
    else if (ctx->MatchNV > 0) {
        printf1(ctx, "Match command will be ignored.\n");
        ctx->MatchNV = 0;
    }
            
    if (ctx->DBLKVar >= 0) {
        printf1(ctx, "\nBlock mode defined by variable: %s\n",ctx->VName[ctx->DBLKVar]);
        i = ctx->VTyp[ctx->DBLKVar];
        if (ctx->DBLKVar >= ctx->NVIdx && i != 2 && i != 3) {
            printf1(ctx, "Error: type of block mode variable must be 2 or 3.\n");
            goto NVFin;
        }
        if (ctx->GDFlg)
            printf1(ctx, "Maximum block size: %d\n",ctx->NOCMaxA);

        if (ctx->BSelPtr != NULL) {              /* bsel option */
            if (check_vsel(ctx))
                goto NVFin;

            if (save_bsel(ctx, 1))
                goto NVFin;
        }
    }
    else if (ctx->BSelPtr != NULL) {
        printf1(ctx, "Block select ignored.\n");
        ctx->BSelPtr = NULL;
    }

    if (ctx->GDFlg && ctx->DBLKVar < 0) {     /* check type of new variables */

        i = ctx->NVIdx;
        while (i >= 0) {
            if (ctx->VTyp[i] > 3) {
                printf1(ctx, "Error: if not in block mode, df option allows only variable types 1,2,3.\n");
                goto NVFin;
            }
            i = ctx->VNxt[i];
        }
    }

    if (ctx->NVArc > 0 || ctx->DFILN > 0)
        printf1(ctx, "\n");         
  
    /* allocate memory for new variables */

    if (alloc_vdat(ctx, ctx->NVIdx,1)) {
        printf1(ctx, "Insufficient memory for new variables.\n");
        goto NVFin;
    }
    err = gen_dm(ctx, ctx->NVIdx,ctx->NVNum,sflag,0);      /* create data */
    if (err == 0) {
        prn_mval(ctx);                 /* info about missing values */
    
        if (ctx->GDFlg && ctx->NVNum > 0 && (ctx->PMTDAFDef || ctx->PMSPSSFDef))
            prn_dfd(ctx, ctx->NVarNV,ctx->NVarIdx);
    }

NVFin:
    if (ctx->GDFlg > 1)                 
        fclose(ctx->GDFd);

    if (ctx->NVarNV > 0) {
        free((char *)ctx->NVarIdx);
        memrq(ctx, -ctx->NVarNV,sizeof(short));
        ctx->NVarNV = 0;
    }
    if (ctx->NVPFmtSA1 > 0) {
        for (i = 0; i < ctx->NVPFmtSA1; ++i) {  
            free((char *)ctx->NVPFmtS[i]); 
            memrq(ctx, -NVPFmtSLen,sizeof(char));
        }
        ctx->NVPFmtSA1 = 0;
    }
    if (ctx->NVPFmtSA > 0) {
        free((char *)ctx->NVPFmtS); 
        memrq(ctx, -ctx->NVPFmtSA,sizeof(char *));
        ctx->NVPFmtSA = 0;
    }
    alloc_vl(ctx, 0);

    save_isel(ctx, 0,0);
    save_vsel(ctx, 0,0);
    save_bsel(ctx, 0);
    save_break(ctx, 0,0);
    check_vcj(ctx, ctx->NVIdx,0);
    check_ffmt(ctx, 0);
    if (ctx->NVArc)
        alloc_avar(ctx, ctx->NVIdx,ctx->NVArc,0);      
    vsort(ctx, ctx->MatchNV,ctx->MatchVB,0,0,1);

    if (err || ctx->GDFlg) {
        if (ctx->NVNum > 0)
            clear_avar(ctx, ctx->NVIdx);
        printf1(ctx, "No new variables created. ");
    }
    else {
        if (dmdef == 0) {
            ctx->WIVar = -1;
            ctx->WSum = (double)ctx->NOC;
            ctx->WSumS = 0.0;
            ctx->WNorm = 1.0;
            ctx->WNormFlag = 0;
        }
        printf1(ctx, "\nEnd of creating new variables. ");
    }
    ctx->GDNRec = ctx->GDFlg = ctx->NVArc = 0;
    prn_mem(ctx);
    /* prnchar('-',LLEN,1); */
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_vsel()    Check vsel, bsel and break                              */ 
/*                  Return 0 if OK, -1 if error.                            */

int check_vsel(TDAContext *ctx)                
{
    if (ctx->NVArc <= 0 && ctx->DFILN <= 0) {
        printf1(ctx, "\nError: isel, vsel, bsel and break can only be used when reading\n"); 
        printf1(ctx, "from a data file or archive.\n");
        return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_dfd(nv,vidx)    Print TDA and SPSS description files.               */

void prn_dfd(TDAContext *ctx, int nv,short *vidx)
{
    register int i,j,k;

    j = 0;
    for (i = 0; i < nv; ++i) {
        k = vidx[i];
        if (k >= 0)
            vidx[j++] = (short)(k);
    }
    newline(ctx);
    if (ctx->PMTDAFDef)  
        dtda(ctx, ctx->GDFNPtr,ctx->GDNRec,j,vidx,1);
    if (ctx->PMSPSSFDef)  
        dspss(ctx, ctx->GDFNPtr,j,vidx,1);
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  check_matchv(iv,idx)    Check whether variable iv is a new variable.    */
/*                          Return 1 if yes, else 0.                        */

int check_matchv(TDAContext *ctx, int iv,int idx)
{
    while (idx >= 0) {
        if (iv == idx)
            return(1);
        idx = ctx->VNxt[idx];
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_mval    Info about missing values.                                  */

void prn_mval(TDAContext *ctx)
{
    int i,n,k,l,ll;

    if (ctx->NVArc > 0) {

        printf1(ctx, "\nMissing values in numerical variables");
        n = 0;
        for (i = 0; i < ctx->NVArc; ++i) {
            if (ctx->AVMBlnk[i] || ctx->AVMStar[i] || ctx->AVMPnt[i] || ctx->AVMGen[i]) {
                if (n == 0) {
                    ll = 8;
                    for (k = 0; k < ctx->NVArc; ++k) {
                        l = (int)(strlen(ctx->VDef[ctx->AVIdx[k]]) - 2);
                        if (ll < l)
                            ll = l;
                    }
                    printf1(ctx, "\nVariable  ");
                    prnchar(ctx, ' ',ll - 8,0);
                    printf1(ctx, "Blanks  Stars  Points  General\n");
                    prnchar(ctx, '-',32 + ll,1);
                }
                printf1(ctx, "%s  ",ctx->VDef[ctx->AVIdx[i]] + 2);
                prnchar(ctx, ' ',ll - (int)strlen(ctx->VDef[ctx->AVIdx[i]]) + 2,0);
                printf1(ctx, "%6d ",ctx->AVMBlnk[i]); 
                printf1(ctx, "%6d ",ctx->AVMStar[i]); 
                printf1(ctx, "%7d ",ctx->AVMPnt[i]); 
                printf1(ctx, "%8d\n",ctx->AVMGen[i]);
                n++;
            }
        }
        if (n == 0)  
            printf1(ctx, " (blank,star,point,general): none\n");
        else  
            printf1(ctx, "\nSubstitution: MBlnk=%g MPnt=%g MStar=%g MGen=%g\n",
                                ctx->MBlnkVal,ctx->MPntVal,ctx->MStarVal,ctx->MGenVal);
    }
    else if (ctx->DFILN > 0) {
   
        if (ctx->MStarN + ctx->MPntN + ctx->MBlnkN + ctx->MGenN == 0) {
            printf1(ctx, "Missing values in data file(s): none.\n");
            return;
        }
        printf1(ctx, "\nMissing values  numerical code\n");
        prnchar(ctx, '-',30,1);        
        printf1(ctx, "Blank %8d      %lg\n",ctx->MBlnkN,ctx->MBlnkVal);
        printf1(ctx, "Star %9d      %lg\n",ctx->MStarN,ctx->MStarVal);
        printf1(ctx, "Point %8d      %lg\n",ctx->MPntN,ctx->MPntVal);
        if (ctx->MGenFlg)
            printf1(ctx, "General %6d      %lg\n",ctx->MGenN,ctx->MGenVal);
    }
}

/* ------------------------------------------------------------------------ */
/*  check_nv    Check nvar command in CmdBuf.                               */
/*              Return 0 if OK, -1 if error.                                */

int check_nv(TDAContext *ctx)
{
    int err,r,i,n,l,m1,m2,nb;
    register char c,*p,*q;
    char sc,vname[VNLMax + 1];
    double x;

    ctx->NVSEPC = ' ';             
    ctx->NVL0 = ctx->VTNOC = ctx->BSIZE = 0;
    err = -1;
    p = ctx->CmdBuf + 4;
    while (*++p) {
        if (sscanf(p,"noc=%d",&n) == 1 && n > 0) {
            ctx->VTNOC = n;
            p = skip_int(ctx, p + 4);
        }
        else if (sscanf(p,"bsize=%d",&n) == 1 && n > 0) {
            ctx->BSIZE = n;
            p = skip_int(ctx, p + 6);
        }
        else if (sscanf(p,"dreclen=%d",&n) == 1 && n > 0) {
            ctx->DRecLen = n;
            p = skip_int(ctx, p + 8);
        }
        else if (sscanf(p,"nmrec=%d",&n) == 1 && n > 0) {
            ctx->NMRec = n;
            p = skip_int(ctx, p + 6);
        }
        else if (sscanf(p,"fmt=%d.%d",&ctx->VFmt1,&ctx->VFmt2) == 2) {
            ctx->VFmtFlg = 1;
            p = skip_int(ctx, p + 4);
            p = skip_int(ctx, p + 1);
        }
        else if (sscanf(p,"mstar=%lf",&x) == 1) {
            ctx->MStarVal = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"mpnt=%lf",&x) == 1) {
            ctx->MPntVal = x;
            p = skip_dbl(ctx, p + 5);
        }
        else if (sscanf(p,"mblnk=%lf",&x) == 1) {
            ctx->MBlnkVal = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"mgen=%lf",&x) == 1) {
            ctx->MGenVal = x;
            ctx->MGenFlg = 1;
            p = skip_dbl(ctx, p + 5);
        }
        else if (sscanf(p,"mmatch=%lf",&x) == 1) {
            ctx->MMatchVal = x;
            p = skip_dbl(ctx, p + 7);
        }
        else if (!strncmp(p,"arcdic",6)) {
            ctx->ARCDICFlg = 1;
            p += 6;
        }
        else if (!strncmp(p,"dfile=",6)) {
            if (ctx->DFILN >= DFILMax) {
                printf1(ctx, "Exceeded max number of data files.\n");
                goto CNVFin;
            }
            ctx->DFILFN[ctx->DFILN++] = p + 6;
            p = skip_nc(ctx, p);
        }
        else if (!strncmp(p,"isel=",5)) {
            ctx->ISelPtr = p + 5;
            p = skip_expr(ctx, p + 5);
        }
        else if (!strncmp(p,"vsel=",5)) {
            ctx->VSelPtr = p + 5;
            p = skip_expr(ctx, p + 5);
        }
        else if (!strncmp(p,"bsel=",5)) {
            ctx->BSelPtr = p + 5;
            p = skip_expr(ctx, p + 5);
        }
        else if (!strncmp(p,"break=",6)) {
            ctx->BRKPtr = p + 6;
            p = skip_expr(ctx, p + 6);
        }
        else if (!strncmp(p,"df=",3)) {
            ctx->GDFNPtr = p + 3;
            ctx->GDFlg = 1;
            p = skip_nc(ctx, p + 3);
        }
        else if (!strncmp(p,"df1=",4)) {
            ctx->GDFNPtr = p + 4;
            ctx->GDKFlg = ctx->GDFlg = 1;
            p = skip_nc(ctx, p + 4);
        }
        else if (!strncmp(p,"keep=",5)) {
            ctx->KeepPtr = p + 5;
            p = get_nvi(ctx, p + 5,&n,0,ctx->VLVIdx,&nb);
            if (n < 1 || nb) {
                printf1(ctx, "Error in keep expression.\n");
                goto CNVFin;
            }
        }
        else if (!strncmp(p,"drop=",5)) {
            ctx->DropPtr = p + 5;
            p = get_nvi(ctx, p + 5,&n,0,ctx->VLVIdx,&nb);
            if (n < 1 || nb) {
                printf1(ctx, "Error in drop expression.\n");
                goto CNVFin;
            }
        }
        else if (!strncmp(p,"dtda",4)) {    /* PMTDA file */
            p = get_file(ctx, 10,p + 4,&r);
            if (r)  
                 goto CNVFin;
        }
        else if (!strncmp(p,"dspss",5)) {   /* PMSPSS file */
            p = get_file(ctx, 11,p + 5,&r);
            if (r)  
                goto CNVFin;
        }
        else if (!strncmp(p,"match=",6)) {
            q = p + 6;
            n = 0;
            for (i = 0; i < MaxMatchV; ++i) {
                if ((ctx->MatchVA[i] = (short)get_vidx(ctx, q,vname)) < 0) {
                    prn_nve1(ctx, p);
                    goto CNVFin;
                }
                q += strlen(vname);

                if (*q++ != ',' || (ctx->MatchVB[i] = (short)get_vidx(ctx, q,vname)) < 0) {
                    prn_nve1(ctx, p);
                    goto CNVFin;
                }
                q += strlen(vname);
                n++;       

                if (*q != ',' || check_vname(ctx, q + 1) == 0)
                    break;
                q++;
            }
            ctx->MatchNV = n;
            p = q;
        }
        else if (sscanf(p,"isc=%c",&sc) == 1) {
            if (sc == 't')
                sc = '\t';
            ctx->ISEPC = sc;
            p += 5;                
        }
        else if (!strncmp(p,"dblock=",7)) {
            ctx->DBLKVar = get_vidx(ctx, p + 7,vname);
            if (ctx->DBLKVar < 0) {
                prn_nve1(ctx, p);
                goto CNVFin;
            }
            p += 7 + strlen(vname);
        }
        else if (!strncmp(p,"ffmt=c",6)) {
            ctx->FFMTCnt = 0;       
            q = p + 5;
            while (1) {
                if ((sscanf(q,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d,%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)) {
                    q = skip_int(ctx, q + 1);
                    q = skip_blev(ctx, q);
                }
                else {
                    prn_nve(ctx, p);
                    goto CNVFin;
                }
                ctx->FFMTCnt++;
                if (*q != ',' || *(q + 1) != 'c')
                    break;
                q++;       
            }
            ctx->FFMTPtr = p + 5;
            p = q;
        }
        else if ((l = get_vnlen(ctx, p)) > 0) {       /* variable */
            q = skip_nc(ctx, p + l);
            c = *q;
            *q = '\0';
            if (save_var(ctx, p,0))
                goto CNVFin;
            ctx->NVNum++;
            *q = c;
            p = q; 
        }
        else if (!strncmp(p,"sepc=none",9)) {
            ctx->NVSEPC = '\0';
            p += 9;
        }
        else if (sscanf(p,"sepc=%c",&ctx->NVSEPC) == 1) {
            if (ctx->NVSEPC == 't')
                ctx->NVSEPC = '\t';
            p += 6;                
        }
        else if (sscanf(p,"l0=%d",&n) == 1) {
            ctx->NVL0 = n;
            p = skip_int(ctx, p + 3);
        }
        if (*p != ',' && *p != ')') {
            prn_nve(ctx, p);
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

void prn_nve(TDAContext *ctx, char *s)
{
    register char *p = s;

    printf1(ctx, "Syntax error: ");
    if (!*p)
        printf1(ctx, "check brackets and semicolon.\n");
    else {     
        while (*p && p < s + 20)
            printf1(ctx, "%c",*p++);
        if (*p)
            printf1(ctx, " ...");
        printf1(ctx, "\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn_nve1    Print error message.                                        */

void prn_nve1(TDAContext *ctx, char *s)
{
    register char *p = s;

    printf1(ctx, "Syntax error or undefined variable: ");
    while (*p && p < s + 30)
        printf1(ctx, "%c",*p++);
    if (*p)
        printf1(ctx, " ...");
    printf1(ctx, "\n");
}

/* ------------------------------------------------------------------------ */
/*  make_vfmt(nidx)   Make new print format, beginning with nidx.           */
/*                    Only for variables with 0.0, not for string var.      */
   
void make_vfmt(TDAContext *ctx, int nidx)
{
    int i,w1,w2;

    i = nidx;   
    while (i >= 0) {              
           
        w1 = (int)ctx->VPFmt1[i];
        w2 = (int)ctx->VPFmt2[i];
                 
        if (w1 == 0 && w2 == 0 && ctx->VTyp[i] != 1) {
            w1 = ctx->VFmt1;
            w2 = ctx->VFmt2;
            makefmt(ctx, &w1,&w2,ctx->VPFmtS[i],VPFmtSLen,0,' ',0);
            ctx->VPFmt1[i] = (short)w1;
            ctx->VPFmt2[i] = (short)w2;
        }
        i = ctx->VNxt[i];
    }
}

/* ------------------------------------------------------------------------ */
/*  make_afmt   Make new print format for archive variables.                */
/*              Only if current format is 0.0.                              */

void make_afmt(TDAContext *ctx)  
{
    int i,j,w1,w2;
   
    for (i = 0; i < ctx->NVArc; ++i) {
        j = ctx->AVIdx[i];
        if (ctx->VTyp[j] == 1)
            continue;

        if (ctx->VPFmt1[j] == 0 && ctx->VPFmt2[j] == 0) {
            w1 = ctx->AVFmt1[i];
            w2 = ctx->AVFmt2[i];
            get_afmt(ctx, &w1,&w2);

            makefmt(ctx, &w1,&w2,ctx->VPFmtS[j],VPFmtSLen,0,' ',0);
            ctx->VPFmt1[j] = (short)w1;
            ctx->VPFmt2[j] = (short)w2;
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

int save_isel(TDAContext *ctx, int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (ctx->ISESCnt > 0) {
            free((char *)ctx->ISESTyp);
            free((char *)ctx->ISESVal);
            memrq(ctx, -ctx->ISESCnt,sizeof(double) + sizeof(int));
        }
        ctx->ISESCnt = 0;
        ctx->ISelPtr = NULL;
        return(0);
    }
    printf1(ctx, "Input select (isel): %s\n",ctx->ISelPtr);

    if ((n = v_parse(ctx, ctx->ISelPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d) in isel expression.\n",n);
        if (n < 0)
            prn_emsg1(ctx, n);
        return(-1);
    }
       
    /* check for type 2 and 3 operators */

    check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,0);
    if (n2 > 0 || n3 > 0) {
        printf1(ctx, "Error: isel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_iref(ctx, ctx->ESCnt,ctx->ESTyp,idx)) {
        printf1(ctx, "Reference error in isel expression.\n");
        return(-1);
    }
    if (!(ctx->ISESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(ctx->ISESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        free((char *)ctx->ISESTyp);
        p_err(ctx, -2,1);
        return(-1);      
    }
    memrq(ctx, ctx->ESCnt,sizeof(double) + sizeof(int));

    ctx->ISESCnt = ctx->ESCnt;
    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->ISESTyp[i] = ctx->ESTyp[i];
        ctx->ISESVal[i] = ctx->ESVal[i];
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
         
int check_iref(TDAContext *ctx, int cnt,int *ptyp,int idx)
{
    register int i,k,t;
                 
    for (i = 0; i < cnt; ++i) {
        t = iabs(ctx, ptyp[i]);
        if (t >= ctx->VOFFS && t < ctx->COFFS) {
            k = idx;
            while (k >= 0) {
                if (t == ctx->VOFFS + k && ctx->VTypA[k] != 2)
                    return(-1);
                k = ctx->VNxt[k];
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

int save_vsel(TDAContext *ctx, int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (ctx->VSESCnt > 0) {
            free((char *)ctx->VSESTyp);
            free((char *)ctx->VSESVal);
            memrq(ctx, -ctx->VSESCnt,sizeof(double) + sizeof(int));
        }
        ctx->VSESCnt = 0;
        ctx->VSelPtr = NULL;
        return(0);
    }
    printf1(ctx, "Input select (vsel): %s\n",ctx->VSelPtr);

    if ((n = v_parse(ctx, ctx->VSelPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d) in vsel expression.\n",n);
        if (n < 0)
            prn_emsg1(ctx, n);
        return(-1);
    }

    /* check for type 2 and 3 operators */

    check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,0);
    if (n2 > 0 || n3 > 0) {
        printf1(ctx, "Error: vsel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_vref(ctx, ctx->ESCnt,ctx->ESTyp,idx)) {
        printf1(ctx, "Error: vsel expression may only refer to variables of type 1 - 3.\n");
        return(-1);
    }
    if (!(ctx->VSESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(ctx->VSESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        free((char *)ctx->VSESTyp);
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, ctx->ESCnt,sizeof(double) + sizeof(int));

    ctx->VSESCnt = ctx->ESCnt;
    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->VSESTyp[i] = ctx->ESTyp[i];
        ctx->VSESVal[i] = ctx->ESVal[i];
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

int check_vref(TDAContext *ctx, int cnt,int *ptyp,int idx)
{
    register int i,k,t;
                 
    for (i = 0; i < cnt; ++i) {

        t = iabs(ctx, ptyp[i]);

        if (t >= ctx->VOFFS && t < ctx->COFFS) {
            k = idx;
            while (k >= 0) {
                if (t == ctx->VOFFS + k && ctx->VTyp[k] >= 4)
                    return(-1);
                k = ctx->VNxt[k];
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

int save_bsel(TDAContext *ctx, int opt)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (ctx->BSESCnt > 0) {
            free((char *)ctx->BSESTyp);
            free((char *)ctx->BSESVal);
            memrq(ctx, -ctx->BSESCnt,sizeof(double) + sizeof(int));
        }
        ctx->BSESCnt = 0;
        ctx->BSelPtr = NULL;
        return(0);
    }
    printf1(ctx, "Block select: %s\n",ctx->BSelPtr);

    if ((n = v_parse(ctx, ctx->BSelPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d) in vsel expression.\n",n);
        if (n < 0)
            prn_emsg1(ctx, n);
        return(-1);
    }

    /* check for type 2 and 3 operators and c terms */

    check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,1);
    if (n2 > 0 || n3 > 0) {
        printf1(ctx, "Error: bsel expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (nc > 0) {
        printf1(ctx, "Error: bsel expressions must not refer to c terms.\n");
        return(-1);
    }
    if (!(ctx->BSESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(ctx->BSESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        free((char *)ctx->BSESTyp);
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, ctx->ESCnt,sizeof(double) + sizeof(int));

    ctx->BSESCnt = ctx->ESCnt;
    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->BSESTyp[i] = ctx->ESTyp[i];
        ctx->BSESVal[i] = ctx->ESVal[i];
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

int save_break(TDAContext *ctx, int opt,int idx)
{
    int i,n,nv,nc,n2,n3;

    if (opt == 0) {
        if (ctx->BRKESCnt > 0) {
            free((char *)ctx->BRKESTyp);
            free((char *)ctx->BRKESVal);
            memrq(ctx, -ctx->BRKESCnt,sizeof(double) + sizeof(int));
        }
        ctx->BRKESCnt = 0;
        ctx->BRKPtr = NULL;
        return(0);
    }
    printf1(ctx, "Break expression: %s\n",ctx->BRKPtr);

    if ((n = v_parse(ctx, ctx->BRKPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d) in break expression.\n",n);
        if (n < 0)
            prn_emsg1(ctx, n);
        return(-1);
    }

    /* check for type 2 and 3 operators */

    check_expr(ctx, ctx->ESCnt,ctx->ESTyp,&nv,&nc,&n2,&n3,1);
    if (n2 > 0 || n3 > 0) {
        printf1(ctx, "Error: break expressions may not contain type 2 and 3 operators.\n");
        return(-1);
    }
    if (check_vref(ctx, ctx->ESCnt,ctx->ESTyp,idx)) {
        printf1(ctx, "Error: break expression may only refer to variables of type 1 - 3.\n");
        return(-1);
    }
    if (!(ctx->BRKESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    if (!(ctx->BRKESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        free((char *)ctx->BRKESTyp);
        p_err(ctx, -2,1);
        return(-1);
    }
    memrq(ctx, ctx->ESCnt,sizeof(double) + sizeof(int));

    ctx->BRKESCnt = ctx->ESCnt;
    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->BRKESTyp[i] = ctx->ESTyp[i];
        ctx->BRKESVal[i] = ctx->ESVal[i];
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

int check_vcj(TDAContext *ctx, int idx,int opt)
{
    register int i,j,k;   
    int ii,n,cmax;

    if (opt == 0) {
        if (ctx->VCJMax > 0) {
            free((char *)ctx->VCJFlg);
            free((char *)ctx->VCJVal);
            memrq(ctx, -ctx->VCJMax - 1,sizeof(char) + sizeof(double));
            ctx->VCJMax = 0;
        }
        return(0);
    }

    cmax = ii = 0;
VCJREP:
    i = idx;
    while (i >= 0) {
        n = ctx->VESCnt[i];
        for (j = 0; j < n; ++j) {
            k = ctx->VESTyp[i][j];
            if (k >= ctx->COFFS && k < ctx->COFFMAX) {
                k -= ctx->COFFS;
                if (ii)
                    ctx->VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
        i = ctx->VNxt[i];
    }
    if (ctx->ISelPtr != NULL) {
        n = ctx->ISESCnt;
        for (j = 0; j < n; ++j) {
            k = ctx->ISESTyp[j];
            if (k >= ctx->COFFS && k < ctx->COFFMAX) {
                k -= ctx->COFFS;
                if (ii)
                    ctx->VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (ctx->VSelPtr != NULL) {
        n = ctx->VSESCnt;
        for (j = 0; j < n; ++j) {
            k = ctx->VSESTyp[j];
            if (k >= ctx->COFFS && k < ctx->COFFMAX) {
                k -= ctx->COFFS;
                if (ii)
                    ctx->VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (ctx->BRKPtr != NULL) {
        n = ctx->BRKESCnt;
        for (j = 0; j < n; ++j) {
            k = ctx->BRKESTyp[j];
            if (k >= ctx->COFFS && k < ctx->COFFMAX) {
                k -= ctx->COFFS;
                if (ii)
                    ctx->VCJFlg[k] = 1;
                else if (cmax < k)
                    cmax = k;
            }
        }
    }
    if (ii || cmax == 0)
        return(0);

    ctx->VCJMax = cmax;
    if (!(ctx->VCJFlg = (char *)calloc((size_t)(cmax + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        return(-1);      
    }
    if (!(ctx->VCJVal = (double *)calloc((size_t)(cmax + 1),sizeof(double)))) {
        free((char *)ctx->VCJFlg);
        p_err(ctx, -2,1);
        return(-1);      
    }
    memrq(ctx, cmax + 1,sizeof(char) + sizeof(double));
    ii++;
    goto VCJREP;
}

/* ------------------------------------------------------------------------ */
/*  check_ffmt(opt)    Create arrays FFMTC1[i] and FFMTC2[i] for ffmt terms */
/*                     i = 1,...,VCJMax. Also, if opt = 1, print info       */
/*                     to standard output.                                  */
/*                     If opt == 0 free previously allocated memory.        */
/*                     Return 0 if OK, -1 if error.                         */

int check_ffmt(TDAContext *ctx, int opt)
{
    register int i;   
    int n,m1,m2,err; 
    register char *p;

    if (opt == 0) {
        if (ctx->FFMTA > 0) {
            free((char *)ctx->FFMTC1);
            free((char *)ctx->FFMTC2);
            memrq(ctx, -ctx->FFMTA,sizeof(int));
            ctx->FFMTCnt = ctx->FFMTA = 0;
            ctx->FFMTPtr = NULL;
        }
        return(0);
    }
    if (!(ctx->FFMTC1 = (int *)calloc((size_t)(ctx->VCJMax + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);      
    }
    if (!(ctx->FFMTC2 = (int *)calloc((size_t)(ctx->VCJMax + 1),sizeof(int)))) {
        free((char *)ctx->FFMTC1);
        p_err(ctx, -2,1);
        return(-1);      
    }
    ctx->FFMTA = 2 * ctx->VCJMax + 2;
    memrq(ctx, ctx->FFMTA,sizeof(int));

    p = ctx->FFMTPtr;
             
    while (1) {
        if ((sscanf(p,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1))   
            ;
        else if (sscanf(p,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)    
            m2 = m1;
        else 
            goto CFFMTErr;

        if (n < 1)
            goto CFFMTErr;
        else if (n <= ctx->VCJMax) {
            ctx->FFMTC1[n] = m1;
            ctx->FFMTC2[n] = m2;
        }
        p = skip_int(ctx, p + 1);
        p = skip_blev(ctx, p);

        if (*p != ',' || *(p + 1) != 'c')
            break;
        p++;       
    }
    err = 0;
    printf1(ctx, "\nUsing fixed format data\n");
    prnchar(ctx, '-',23,1);
    for (i = 1; i <= ctx->VCJMax; ++i) {
        if (ctx->VCJFlg[i] == 0)
            continue;
        printf1(ctx, "c%-4d ",i);
        if (ctx->FFMTC1[i] > 0) {
            printf1(ctx, "column%4d",ctx->FFMTC1[i]);
            if (ctx->FFMTC2[i] > ctx->FFMTC1[i])
                printf1(ctx, " -%4d",ctx->FFMTC2[i]);
            printf1(ctx, "\n");
        }
        else {
            printf1(ctx, "error, need column information\n");
            err = -1;
        }
    }
    newline(ctx);
    return(err);

CFFMTErr:
    printf1(ctx, "Error in ffmt parameter.\n");
    return(-1);
}

/*--------------------------------------------------------------------------*/
/*  get_data(j,i)       Return the value of the j.th variable for case i.   */
/*                      i = 0,...,NOC - 1.                                  */
/*                      Note that the data matrix contains NOCDM cases,     */
/*                      NOC is the number of cases after tsel selection.    */
/*                                                                          */
/*                      If TSelFlg == 1, apply tsel case selection.         */
          
double get_data(TDAContext *ctx, int j,int i)
{
    char *p0,*p1;
    short *p2;
    float *p4;
    int *p5;
    double *p8;

    if (ctx->VTyp[j] == 1)                   /* string variables */
        return(0.0);
    else if (ctx->VTyp[j] == 2) {            /* numerical constants */
        p8 = (double *)(void *)ctx->VDPtr[j];
        return((double)p8[0]);
    }
    if (ctx->TSelFlg)                        /* get index from TSelect */
        i = ctx->TSelect[i];     
    else if (ctx->REPSelFlg)                 /* get index from REPSelect */
        i = ctx->REPSelect[i];

    switch (ctx->VSLen[j]) {
        
        case 0: p0 = ctx->VDPtr[j] + i / 8;
  
                if (*p0 & ctx->BMsk[i % 8])
                    return(1.0);
                else
                    return(0.0);
           
        case 1:  p1 = (char *)(void *)ctx->VDPtr[j];
                 return((double)p1[i]);

        case 2:  p2 = (short *)(void *)ctx->VDPtr[j];
                 return((double)p2[i]);

        case 4:  p4 = (float *)(void *)ctx->VDPtr[j];
                 return((double)p4[i]);

        case 5:  p5 = (int *)(void *)ctx->VDPtr[j];
                 return((double)p5[i]);

        case 8:  p8 = (double *)(void *)ctx->VDPtr[j];
                 return((double)p8[i]);

        default: break;
    }
    return(0.0);
}

/*--------------------------------------------------------------------------*/
/*  put_data(x,j,i)     Put value x in data matrix for variable j in case i */
/*                                                                          */
/*                      If TSelFlg == 1, apply tsel case selection.         */
          
void put_data(TDAContext *ctx, double x,int j,int i)
{
    char *p0,*p1;
    short *p2;
    float *p4;
    int *p5;
    double *p8;

    if (ctx->VTyp[j] == 1)               /* string variables */
        return;
    else if (ctx->VTyp[j] == 2) {        /* numerical constants */
        p8 = (double *)(void *)ctx->VDPtr[j];
        p8[0] = (double)x;
        return;
    }
    if (ctx->TSelFlg)                        /* get index to next selected case */
        i = ctx->TSelect[i];     
    else if (ctx->REPSelFlg)                 /* get index from REPSelect */
        i = ctx->REPSelect[i];

    switch (ctx->VSLen[j]) {

        case 0:  p0 = ctx->VDPtr[j] + i / 8;
   
                 if ((int)x)
                    *p0 |= ctx->BMsk[i % 8];
                 else
                    *p0 &= ~ctx->BMsk[i % 8];
                 break;

        case 1:  p1 = (char *)(void *)ctx->VDPtr[j];
                 p1[i] = (char)x;
                 break;
        case 2:  p2 = (short *)(void *)ctx->VDPtr[j];
                 p2[i] = (short)x;
                 break;
        case 4:  p4 = (float *)(void *)ctx->VDPtr[j];
                 p4[i] = (float)x;
                 break;
        case 5:  p5 = (int *)(void *)ctx->VDPtr[j];
                 p5[i] = (int)x;
                 break;
        case 8:  p8 = (double *)(void *)ctx->VDPtr[j];
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
          
void put_str(TDAContext *ctx, char *buf,int blen,int j,int i,int bflag)
{
    register int k;
    register char *p,*q,*m;
    int len;

    if (ctx->TSelFlg)                        /* get index to next selected case */
        i = ctx->TSelect[i];     
    else if (ctx->REPSelFlg)                 /* get index from REPSelect */
        i = ctx->REPSelect[i];

    len = -ctx->VSLen[j];
    p = buf;
    q = ctx->VDPtr[j] + i * len;
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
          
void clear_str(TDAContext *ctx, char *buf,int blen,int j)
{
    register int k;
    register char *p,*m;
    int n,len;

    n = ctx->VStrN[j];               /* first column of string */
    len = -ctx->VSLen[j];

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
          
void get_str(TDAContext *ctx, char *buf,int j,int i)
{
    register int k;
    register char *p,*q;
    int len;

    if (ctx->TSelFlg)                        /* get index to next selected case */
        i = ctx->TSelect[i];     
    else if (ctx->REPSelFlg)                 /* get index from REPSelect */
        i = ctx->REPSelect[i];

    len = -ctx->VSLen[j];
    p = ctx->VDPtr[j] + i * len;
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

double dscan(TDAContext *ctx, char *p,int len,int *mval)
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

        if (sep_char(ctx, *q))
            break;

        if (*q == '*') {
            *mval = 2;
            return(ctx->MStarVal);
        }
        else if (*q >= '0' && *q <= '9') {
            itmp *= 10L;
            itmp += (long)(*q - '0');
            k++;
        }
        else if (*q == '.') {
            if (pt) {
                *mval = 4;
                return(ctx->MGenVal);
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
            return(ctx->MGenVal);
        }
        q++;
        if (len && ++i >= len)
            break;
    }
    if (!k) {
        if (pt) {
            *mval = 3;
            return(ctx->MPntVal);
        }
        *mval = 1;
        return(ctx->MBlnkVal);
    }
    else if (k > MaxMDig) {

        if (len) {
            c = *(p + len);
            *(p + len) = '\0';
        }
        if (sscanf(p,"%lg",&tmp) != 1) {
            *mval = 4;
            tmp = ctx->MGenVal;
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

int sep_char(TDAContext *ctx, char c)
{
    if (ctx->ISEPC != '\0') {
        if (c == ctx->ISEPC)     
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

char *skip_sep(TDAContext *ctx, char *p)
{
    if (ctx->ISEPC != '\0') {
        if (sep_char(ctx, *p))
            p++;
    }
    else {
        while (sep_char(ctx, *p))
            p++;
    }
    return(p);
}

/*--------------------------------------------------------------------------*/
/*  skip_dval(p)    skip an entry in a data file's record.                  */
/*                  return pointer to the next character.                   */

char *skip_dval(TDAContext *ctx, char *p)
{
    while (*p) {
        if (*p == '"')      
            while (*++p && *p != '"') ;

        if (sep_char(ctx, *p))
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

int gen_dm(TDAContext *ctx, int idx,int nvar,int sflag,int sdflag)
{
    register int j = 0,k = 0,l = 0,ll = 0;
    int rm = 0,rn = 0,r = 0,err = 0,jm = 0,pcnt = 0,drec = 0,dflag = 0,nb = 0,bcnt = 0,bi = 0,nm = 0,bmin = 0,bmax = 0;
    int nr = 0,wflag = 0,first = 0,ns = 0,nss = 0,nst = 0,fptr = 0;
    double x = 0.0,id = 0.0,lastid = 0.0;
           
    if (ctx->MatchNV > 0) {              
        if (alloc_acx(ctx, ctx->MatchNV))     /* use AcX for match values */
            return(-1);
        if (alloc_acn(ctx, ctx->NOCDM))       /* use AcN to indicate filled rows */
            return(-1);
    }
    first = 1;
    wflag = bmax = err = dflag = 0;
    bmin = ctx->INTMAX;
    if (ctx->DFILN > 0 || ctx->NVArc > 0)
        dflag = 1;

    drec = 0;   /* count data file records */
    ctx->MBlnkN = ctx->MStarN = ctx->MGenN = ctx->MPntN = 0;    /* counter for miss values */

    if (ctx->NVArc > 0) {            /* init reading from archive */
        dbf_init(ctx, 1);
        if (ctx->SILENTFlg < 2)
            printfe(ctx, "Reading archive file: %s\n",ctx->ZAFNam[ctx->AVDFN]);

        j = ctx->VIFirst;            /* clear AVVAL */
        while (j >= 0) {
            ctx->AVVAL[j] = 0.0;
            j = ctx->VNxt[j];
        }
    }
    else if (ctx->DFILN > 0) {       /* use data files */

        if (read_dfi(ctx, 1,sflag)) {    /* init read_df function */
            err = -1;
            goto GDMFin;
        }
    }      

    lastid = ctx->DBLMAX;
    nb = bi = bcnt = 0;

    pcnt = ns = k = 0;  /* count data matrix rows */

    while (1) {
                    
        if (ctx->GDFlg && k >= ctx->NOCMaxA) {
            printf1(ctx, "Error: exceeded maximum block size.\n");
            err = -1;
            break;
        }

        if (ctx->NVArc > 0) {    /* ## get archive variables into AVVAL[] */
                            /* only numerical variables */
            r = get_avar(ctx);
            if (r)  
                prn_message(ctx, ++pcnt,0,0);
        }
        else if (ctx->DFILN > 0) {       /* get next record from data file(s) */
            r = read_df(ctx, sflag,idx); /* values are in VCJVal[]. */
            if (r)  
                prn_message(ctx, ++pcnt,0,0);
        }
        else
            r = 1;

        if (r != 1) {               /* r = 0 if EOF */
            err = r;
            break;
        }       
        if (dflag) {
            drec++;
            if (ctx->VTNOC > 0 && drec > ctx->VTNOC)
                break;
        }
        else
            drec = k + 1;

        /* Check isel. Note that v_eval1() is called with vflg = 1,
           that is: values are assumed to be in VCJVal[] and AVVAL[]. */

        if (ctx->ISESCnt) {

            r = v_eval1(ctx, k,ctx->ISESCnt,ctx->ISESTyp,ctx->ISESVal,ctx->ESIdx,&x,1,0,0,0,0);
             
            if (r) {
                printf1(ctx, "Can't evaluate isel expression for record %d in case %d.\n",
                                                          drec,k + 1);
                prn_emsg2(ctx, r);
                err = -1;
                break;
            }
            if (fabs(x) < ctx->EPSI2)  
                goto DMCONT;
        }
        if (ctx->MatchNV > 0) {

            for (j = 0; j < ctx->MatchNV; ++j) {
                jm = ctx->MatchVA[j];
                r = v_eval1(ctx, k,ctx->VESCnt[jm],ctx->VESTyp[jm],ctx->VESVal[jm],ctx->ESIdx,&x,1,0,0,0,0);
                if (r) {
                    printf1(ctx, "Can't evaluate variable %s for sorted case %d (record %d).\n",
                                                    ctx->VName[jm],k + 1,drec);
                    prn_emsg2(ctx, r);
                    err = -1;
                    break;
                }
                ctx->AcX[j] = x;
            }
            if (err)
                break;

            rm = fnd_match(ctx, ctx->AcX,&rn);  /* check match */
            if (rm < 0)
                goto DMCONT;

            /* matching is for rn rows beginning at row rm */

        }

        /* put string variables into data matrix */

        j = idx;
        while (j >= 0) {
            if (ctx->VTyp[j] == 1 && ctx->VTypA[j] == 0) {    
                if (ctx->MatchNV == 0)
                    put_str(ctx, ctx->RSBuf + ctx->VStrN[j],ctx->RBufALen - ctx->VStrN[j],j,k,0); 
                else {
                    for (l = 0; l < rn; ++l) {
                        ll = ctx->VSORTPtr[rm + l];
                        if (ctx->AcN[ll] == 0)  
                            put_str(ctx, ctx->RSBuf + ctx->VStrN[j],ctx->RBufALen - ctx->VStrN[j],j,ll,0); 
                        else
                            wflag = 1;
                    }
                }
            }
            j = ctx->VNxt[j];
        }

        /* ## put archive variables into data matrix */

        if (ctx->NVArc > 0) {
            for (j = 0; j < ctx->NVArc; ++j) {
                if (ctx->MatchNV == 0) {
                    if (ctx->VTyp[ctx->AVIdx[j]] == 1) {  /* string var */
                        get_astr(ctx, ctx->SVBuf,j);
                        put_str(ctx, ctx->SVBuf,ctx->SVBufLen,ctx->AVIdx[j],k,0);
                    }
                    else
                        put_data(ctx, ctx->AVVAL[ctx->AVIdx[j]],ctx->AVIdx[j],k);
                }
                else {
                    for (l = 0; l < rn; ++l) {
                        ll = ctx->VSORTPtr[rm + l];
                        if (ctx->AcN[ll] == 0) {
                            if (ctx->VTyp[ctx->AVIdx[j]] == 1) {  /* string var */
                                get_astr(ctx, ctx->SVBuf,j);
                                put_str(ctx, ctx->SVBuf,ctx->SVBufLen,ctx->AVIdx[j],ll,0);
                            }
                            else
                                put_data(ctx, ctx->AVVAL[ctx->AVIdx[j]],ctx->AVIdx[j],ll);
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

            if (ctx->VTyp[j] == 3 || (ctx->VTyp[j] == 2 && first == 1)) {

                /* archive variables already in data matrix */
                       
                if (ctx->MatchNV == 0) {

                    r = v_eval1(ctx, k,ctx->VESCnt[j],ctx->VESTyp[j],ctx->VESVal[j],ctx->ESIdx,&x,0,0,0,0,0);
                    if (r) {
                        printf1(ctx, "Can't evaluate variable %s for case %d.\n",
                                                         ctx->VName[j],k + 1);
                        if (ctx->VCJMax > 0 || ctx->NVArc > 0)  
                            printf1(ctx, "Evaluation based on data file record %d.\n",drec);
                        prn_emsg2(ctx, r);
                        err = -1;
                        break;          
                    }
                    put_data(ctx, x,j,k);      
                }
                else {

                    for (l = 0; l < rn; ++l) {
                        ll = ctx->VSORTPtr[rm + l];
                        if (ctx->AcN[ll] == 0) {
                            r = v_eval1(ctx, ll,ctx->VESCnt[j],ctx->VESTyp[j],ctx->VESVal[j],ctx->ESIdx,&x,0,0,0,0,0);
                            if (r) {
                                printf1(ctx, "Can't evaluate variable %s for sorted case %d.\n",
                                                            ctx->VName[j],ll + 1);
                                if (ctx->VCJMax > 0 || ctx->NVArc > 0)  
                                    printf1(ctx, "Evaluation based on data file record %d.\n",drec);
                                prn_emsg2(ctx, r);
                                err = -1;
                                break;          
                            }
                            put_data(ctx, x,j,ll);      
                        }
                        else
                            wflag = 1;
                    }
                    if (err)
                        break;
                }
            }
            j = ctx->VNxt[j];
        }
        first = 0;
        if (err)
            break;

        /* if block mode check ID variable and evaluate type 4 variables */
        /* cannot be used if matching data */

        if (ctx->DBLKVar >= 0) {

            id = get_data(ctx, ctx->DBLKVar,k);
            if (fabs(id - lastid) > ctx->EPSI1) {

                if (k > bi) {                   /* new block */
                    err = put_t4var(ctx, idx,bi,k,nb);  /* evaluate type 4 variables */
                    if (err)
                        break;

                    if (ctx->BSESCnt) {              /* check bsel */
                        nr = tst_bsel(ctx, bi,k,idx);
                        if (nr <= 0) {
                            if (nr < 0) {
                                err = nr;
                                break;
                            }   
                            cpy_dmrow(ctx, k,bi,idx);
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


                    if (ctx->GDFlg) {            /* write to output file */
                        prn_gdf(ctx, nr);
                        cpy_dmrow(ctx, k,bi,idx);
                        if (ctx->VTNOC > 0 && ctx->GDNRec >= ctx->VTNOC)
                            break;
                    }
                    else {
                        cpy_dmrow(ctx, k,bi + nr,idx);
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

        if (ctx->VSESCnt && ctx->DBLKVar < 0) {   /* check vsel if not in block mode */

            if (ctx->MatchNV == 0) {

                r = tst_vsel(ctx, k,drec);
                if (r < 0) {
                    err = r;
                    break;
                }
                else if (r == 0)
                    goto DMCONT;
            }
            else {
                for (l = 0; l < rn; ++l) {
                    ll = ctx->VSORTPtr[rm + l];
                    if (ctx->AcN[ll] == 0) {
                        r = tst_vsel(ctx, ll,drec);
                        if (r < 0) {
                            err = -1;
                            break;
                        }
                        else if (r > 0) {
                            ctx->AcN[ll] = 1;
                            ns++;
                        }
                    }
                }
                if (err)
                    break;
            }
        }
        else if (ctx->MatchNV > 0) {     /* indicate values */

            for (l = 0; l < rn; ++l) {
                ll = ctx->VSORTPtr[rm + l];
                if (ctx->AcN[ll] == 0) {
                    ctx->AcN[ll] = 1;
                    ns++;
                }
            }
        }  

        if (ctx->DBLKVar < 0) {      /* if not in block mode */

            /* Check break. Note that v_eval1() is called with vflag=0, that is,
               values are assumed to be already in data matrix. Cannot be used
               if matching data. */

            if (ctx->BRKESCnt) {
  
                r = v_eval1(ctx, k,ctx->BRKESCnt,ctx->BRKESTyp,ctx->BRKESVal,ctx->ESIdx,&x,0,0,0,0,0);
                if (r) {
                    printf1(ctx, "Can't evaluate break expression for record %d in case %d.\n",drec,k + 1);
                    prn_emsg2(ctx, r);
                    err = -1;
                    break;
                }
                if (fabs(x) >= ctx->EPSI2)  
                    break;         
            }
        }

DMCONTA:
        k++;                    /* count data matrix rows */

        if (ctx->GDFlg && ctx->DBLKVar < 0) {
            prn_gdf(ctx, k);
            k = 0;
            if (ctx->VTNOC > 0 && ctx->GDNRec >= ctx->VTNOC)
                break;
            goto DMCONT;
        }
        if (sdflag) {        /* special for spatial data */

            nst = (int)(get_data(ctx, ctx->SDVarSDTyp,k - 1));    /* type of object */
            if (nst == 1)
                ctx->SDVarNP++;
            else if (nst == 2)
                ctx->SDVarNL++;
            else if (nst == 3)
                ctx->SDVarNPol++;
            else 
                ctx->SDVarNU++;

            nss = (int)(get_data(ctx, ctx->SDVarSDN,k - 1));    /* number of following records */

            if (nss < 1) {
                printf1(ctx, "Error: SDN contains a value less than 1 (object %d).\n",k);
                err = -1;
                break;
            }
            if (nst == 1) {
                if (nss != 1) {
                    printf1(ctx, "Error: type 1 object (%d) has more than one point.\n",k);
                    err = -1;
                    break;
                }
                ctx->SDVarNT += 1;
            }
            else if (nst == 2) {
                if (nss < 2) {
                    printf1(ctx, "Error: type 2 object (%d) has less than two points.\n",k);
                    err = -1;
                    break;
                }
                ctx->SDVarNT += nss;
            }
            else if (nst == 3) {
                if (nss < 3) {
                    printf1(ctx, "Error: type 3 object (%d) has less than three points.\n",k);
                    err = -1;
                    break;
                }
                ctx->SDVarNT += nss;
            }
            if ((fptr = skip_df(ctx, nss,1)) < 0) {       
                printf1(ctx, "Error: cannot read %d data records for object %d.\n",nss,k);
                err = -1;
                break;
            }
            put_data(ctx, (double)fptr,ctx->SDVarSDPtr,k - 1);      
            ctx->SDVarMax = imax(ctx, ctx->SDVarMax,nss);
        }
        if (ctx->MatchNV == 0) {
            if (k >= ctx->NOCMaxA || (ctx->NOCDM > 0 && k >= ctx->NOCDM))
                break;
        }
        else if (ns >= ctx->NOCDM)
            break;

DMCONT: 
        if (ctx->VTNOC > 0 && k >= ctx->VTNOC)
            break;
    }
    if (err)
        goto GDMFin;

    if (dflag)
        prn_message(ctx, pcnt,1,0);

    if (k > bi) {
        err = put_t4var(ctx, idx,bi,k,nb);  /* put type 4 variables into data matrix */
        if (err)
            goto GDMFin;             

        if (ctx->DBLKVar >= 0) {                 /* block mode */
            if (ctx->BSESCnt) {                  /* check bsel */
                nr = tst_bsel(ctx, bi,k,idx);
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

                if (ctx->GDFlg) {            /* write to output file */
                    prn_gdf(ctx, nr);
                    k = 0;
                }   
                else
                    k = bi + nr;    
            }
        }
    }
    nm = 0;                     /* number of cases with no match */
    if (ctx->NOCDM == 0)
        ctx->NOC = ctx->NOCDM = k;

    else if (ctx->MatchNV > 0) {     /* insert missing values for non-match */
        nm = ctx->NOCDM - ns;
        if (nm > 0) {
            for (l = 0; l < ctx->NOCDM; ++l) {
                if (ctx->AcN[l] == 0) {
                    j = idx;
                    while (j >= 0) {
                        if (first || ctx->VTyp[j] != 2) {
                            if (ctx->VTyp[j] != 1)
                                put_data(ctx, ctx->MMatchVal,j,l);
                            else if (ctx->VTyp[j] == 1)
                                put_str(ctx, ctx->RSBuf + ctx->VStrN[j],ctx->RBufALen - ctx->VStrN[j],j,l,1); 
                        }
                        j = ctx->VNxt[j];
                    }
                }
            }
        }
    }
    else if (k < ctx->NOCDM) {     /* fill with miss values */

        nm = ctx->NOCDM - k;
        while (k < ctx->NOCDM) {
            j = idx;
            while (j >= 0) {
                if (ctx->VTyp[j] != 1)
                    put_data(ctx, ctx->MMatchVal,j,k);
                else
                    put_str(ctx, ctx->RSBuf + ctx->VStrN[j],ctx->RBufALen - ctx->VStrN[j],j,k,1); 
                j = ctx->VNxt[j];
            }
            k++;
        }
    }
    if (dflag) {
        printf1(ctx, "Read records: %d ",drec);
        if (ctx->VCJMax > 0 && ctx->NMRec > 1)  
            printf1(ctx, "(%d physical records).",ctx->NRec);
        printf1(ctx, "\n");

        /* a data file that yields zero records used to complete
           "successfully", leaving a matrix whose variables point at
           nothing -- the NEXT command touching it crashed (any
           command: dstat, gdd, ...).  Refuse cleanly instead. */
        if (drec < 1) {
            printf1(ctx, "Error: no records read from data file.\n");
            err = -1;
            goto GDMFin;
        }
    }

GDMFin:  
    if (ctx->MatchNV > 0) {
        alloc_acx(ctx, 0);
        alloc_acn(ctx, 0);
    }
    if (ctx->DFILN > 0)     /* free read buffer and reset values */
        read_dfi(ctx, 0,0);    

    if (err == 0) {
        if (ctx->GDFlg == 0) {
            if (ctx->ISESCnt || ctx->VSESCnt || ctx->BRKESCnt || ctx->BSESCnt)  
                printf1(ctx, "Selected for data matrix: %d records.\n",ctx->NOCDM);
        }
    }                       
    if (err || ctx->NOCDM <= 0) {
        if (err)
            err = -1;
          
        if (ctx->GDFlg == 0) {
            if (ctx->DMDef == 0)
                printf1(ctx, "No data matrix created.\n");
            else
                printf1(ctx, "Nothing added to current data matrix.\n");
        }
        else if (err == 0) {
            printf1(ctx, "%d records written to output file: %s\n",ctx->GDNRec,ctx->GDFNPtr);
            if (ctx->DBLKVar >= 0) {
                printf1(ctx, "Number of blocks: %d\n",nb);
                printf1(ctx, "Min number of cases per block: %d\n",bmin);
                printf1(ctx, "Max number of cases per block: %d\n",bmax);
            }
        }
        return(err);
    }
    if (ctx->DMDef == 0) {
        printf1(ctx, "\nCreated a new data matrix.\n");
        printf1(ctx, "Number of cases: %d\n",ctx->NOCDM);
        printf1(ctx, "Number of variables: %d\n",nvar);
        ctx->DMDef = 1;
    }
    else {
        printf1(ctx, "\nAdded %d variable(s) to existing data matrix.\n",nvar);
        printf1(ctx, "Number of cases with no match: %d\n",nm);
        if (nm > 0)  
            printf1(ctx, "Substituted missing value code: %g\n",ctx->MMatchVal);
    }
    if (wflag) {
        printf1(ctx, "Warning: values in variable(s) %s",ctx->VName[ctx->MatchVA[0]]);
        for (j = 1; j < ctx->MatchNV; ++j)  
            printf1(ctx, ",%s",ctx->VName[ctx->MatchVA[j]]);
        printf1(ctx, " are not unique.\n");
    }
    if (ctx->DBLKVar >= 0) {
        printf1(ctx, "\nNumber of blocks: %d\n",nb);
        printf1(ctx, "Min number of cases per block: %d\n",bmin);
        printf1(ctx, "Max number of cases per block: %d\n",bmax);
    }

    /* set VTypA = 1 for new variables */

    j = idx;    
    while (j >= 0) {
        ctx->VTypA[j] = 1;
        j = ctx->VNxt[j];
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  prn_gdf(n)          Print n records from data matrix to output file     */
/*                      GDFd. Use variables in NVarIdx[].                   */
/*                      Use print formats in NVPFmtS[]                      */
/*                                                                          */
/*                      If GDKFlg != 0 keep only first record in block.     */  

void prn_gdf(TDAContext *ctx, int n)
{     
    register int i,j,k;   
    double tmp;

    for (i = 0; i < n; ++i) {
        if (i > 0 && ctx->GDKFlg && ctx->DBLKVar >= 0)
            break;
        for (j = 0; j < ctx->NVarNV; ++j) {
            k = ctx->NVarIdx[j];
            if (k >= 0) {
                if (ctx->VTyp[k] != 1) {
                    tmp = get_data(ctx, k,i);
                    rt_fprintf_d(ctx, ctx->GDFd,ctx->NVPFmtS[j],tmp);
                }
                else {  
                    get_str(ctx, ctx->SVBuf,k,i);
                    fprintf(ctx->GDFd,"%s",ctx->SVBuf);
                    if (ctx->NVSEPC)
                        fprintf(ctx->GDFd,"%c",ctx->NVSEPC);
                }
            }
        }
        fprintf(ctx->GDFd,"\n");
        ctx->GDNRec++;
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

int fnd_match(TDAContext *ctx, double *xm,int *n)
{     
    register int i,j,l,r;   
    int fnd,jj,nn;
    double tmp;                  

    *n = 0;
    l = 0;
    r = ctx->NOCDM - 1;

    while (r >= l) {

        j = (l + r) / 2;

        fnd = 1;
        tmp = get_data(ctx, ctx->MatchVB[0],ctx->VSORTPtr[j]);
        if (fabs(tmp - xm[0]) > ctx->EPSI2)  
            fnd = 0;

        if (fnd) {

            while (j > 0) {
                if (fabs(tmp - get_data(ctx, ctx->MatchVB[0],ctx->VSORTPtr[j - 1])) > ctx->EPSI2) 
                    break;
                j--;
            }
            jj = -1;
            nn = 0;
            fnd = 0;
            while (j < ctx->NOCDM) {
                if (fabs(tmp - get_data(ctx, ctx->MatchVB[0],ctx->VSORTPtr[j])) > ctx->EPSI2)  
                    break;
  
                fnd = 1;
                for (i = 1; i < ctx->MatchNV; ++i) {
                    if (fabs(xm[i] - get_data(ctx, ctx->MatchVB[i],ctx->VSORTPtr[j])) > ctx->EPSI2) {
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
          
int put_t4var(TDAContext *ctx, int idx,int r1,int r2,int bn)
{
    int j,r;
    double tmp;

    j = idx;
    while (j >= 0) {

        if (ctx->VTyp[j] == 4) {
            r = v_eval2(ctx, j,ctx->VESCnt[j],ctx->VESTyp[j],ctx->VESVal[j],ctx->ESIdx,&tmp,r1,r2,0,bn);
            if (r) {
                printf1(ctx, "\nCan't evaluate variable %s in block containing cases %d - %d.\n",ctx->VName[j],r1 + 1,r2);
                prn_emsg2(ctx, r);
                return(-1);
            }
        }
        j = ctx->VNxt[j];
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  cpy_dmrow(ir,ij,idx)    Copy data matrix row from ir to ij, use         */
/*                          variables beginning at idx.                     */

void cpy_dmrow(TDAContext *ctx, int ir,int ij,int idx)
{
    register int j;
    double tmp;

    if (ir == ij)
        return;

    j = idx;
    while (j >= 0) {
        if (ctx->VTyp[j] != 1) {
            tmp = get_data(ctx, j,ir);
            put_data(ctx, tmp,j,ij);      
        }
        else {
            get_str(ctx, ctx->SVBuf,j,ir);
            put_str(ctx, ctx->SVBuf,ctx->SVBufLen,j,ij,0);
        }
        j = ctx->VNxt[j];
    }
}

/*--------------------------------------------------------------------------*/
/*  tst_bsel(r1,r2,idx)     Test bsel expression for cases r1 ... r2.       */
/*                          Copy selected rows into contiguous block.       */
/*                          Return -1 if error, or number in block.         */

int tst_bsel(TDAContext *ctx, int r1,int r2,int idx)
{
    register int i,k;
    int r,n;
    double tmp;

    n = 0;
    i = r1;
    for (k = r1; k < r2; ++k) {
        r = v_eval1(ctx, k,ctx->BSESCnt,ctx->BSESTyp,ctx->BSESVal,ctx->ESIdx,&tmp,0,0,0,0,0);
        if (r) {
            printf1(ctx, "Can't evaluate vsel expression in case %d.\n",k + 1);
            prn_emsg2(ctx, r);
            return(-1);
        }
        if (fabs(tmp) >= ctx->EPSI2) {   /* copy row k into row i */
            cpy_dmrow(ctx, k,i,idx);
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

int tst_vsel(TDAContext *ctx, int i,int rec)
{
    int r;
    double tmp;

    r = v_eval1(ctx, i,ctx->VSESCnt,ctx->VSESTyp,ctx->VSESVal,ctx->ESIdx,&tmp,0,0,0,0,0);
    if (r) {
        printf1(ctx, "Can't evaluate vsel expression for record %d in case %d.\n",rec,i + 1);
        prn_emsg2(ctx, r);
        return(-1);
    }
    if (fabs(tmp) < ctx->EPSI2)  
        return(0);                 
    return(1);
}

/* -------------------------------------------------------------------------*/
/*  tst_break(r1,r2)    Test break expression for cases r1 ... r2.          */
/*                      Return  -1  if error                                */  
/*                               1  if true for at least one record.        */  
/*                               0  if to be used                           */

int tst_break(TDAContext *ctx, int r1,int r2)
{
    int k,r;
    double tmp;

    for (k = r1; k < r2; ++k) {
        r = v_eval1(ctx, k,ctx->BRKESCnt,ctx->BRKESTyp,ctx->BRKESVal,ctx->ESIdx,&tmp,0,0,0,0,0);
        if (r) {
            printf1(ctx, "Can't evaluate break expression in case %d.\n",k + 1);
            prn_emsg2(ctx, r);
            return(-1);
        }
        if (fabs(tmp) >= ctx->EPSI2)  
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
          
int read_dfi(TDAContext *ctx, int opt,int sflag)
{
    if (opt) {
        if (ctx->DRecLen == 0)  
            ctx->RBufLen = RLMaxDef + 1;
        else {
            ctx->DRecLen1 = ctx->DRecLen * ctx->NMRec;
            ctx->RBufLen = ctx->DRecLen1 + 1;
        }
        if (!(ctx->RBuf = (char *)calloc((size_t)(ctx->RBufLen + 1),sizeof(char)))) {
            p_err(ctx, -2,1);
            ctx->RBufLen = 0;
            return(-1);
        }
        if (sflag) {
            if (!(ctx->RSBuf = (char *)calloc((size_t)(ctx->RBufLen + 1),sizeof(char)))) {
                p_err(ctx, -2,1);
                free(ctx->RBuf);
                ctx->RBufLen = 0;
                return(-1);
            }
            ctx->RSBufA = ctx->RBufLen + 1;
            memrq(ctx, ctx->RSBufA,sizeof(char));          
        }
        memrq(ctx, ctx->RBufLen + 1,sizeof(char));          
        ctx->RBufA = ctx->RBufLen + 1;
        ctx->NRec1 = ctx->NRec = 0;
        ctx->DFILNI = 0;     /* first data file */
        ctx->REOF = 1;       /* not open */
    }
    else {
        if (ctx->REOF == 0)
            fclose(ctx->RFd);

        if (ctx->RBufA > 0) {
            free(ctx->RBuf);
            memrq(ctx, -ctx->RBufA,sizeof(char));          
            ctx->RBufA = 0;
        }
        if (ctx->RSBufA > 0) {
            free(ctx->RSBuf);
            memrq(ctx, -ctx->RSBufA,sizeof(char));          
            ctx->RSBufA = 0;
        }
        ctx->RBufLen = 0;
    }
    ctx->REOF = 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  read_df()       read the next record from input data file(s) and save   */
/*                  results in VCJVal[]. It is assumed that the read buffer */
/*                  RBuf is already allocated.                              */
/*                                                                          */
/*  Note: files are only read up to lines that begin with ^Z (hex 1a).      */
/*  Return 1 if OK, 0 if EOF, -1 if error.                                  */
          
int read_df(TDAContext *ctx, int sflag,int idx)
{
    int i,l,rc,mval,err,cnt;
    register char *p,*q;

    ctx->RBufALen = err = 0;

    while (1) {
        if (ctx->REOF) {
            if (ctx->DFILNI >= ctx->DFILN)    /* no more data files */
                break;       

            if (!(ctx->RFd = fopen(ctx->DFILFN[ctx->DFILNI],OPEN_RD))) {
                printf1(ctx, "Can't open data file: %s\n",ctx->DFILFN[ctx->DFILNI]);
                return(-1);
            }
            ctx->NRec1 = ctx->REOF = 0;
            printf1(ctx, "Reading data file: %s\n",ctx->DFILFN[ctx->DFILNI]);
            ctx->DFILNI++;
        }            
        while (1) {
            if (ctx->DRecLen == 0) {

                if (fgets(ctx->RBuf,ctx->RBufLen,ctx->RFd) == NULL || *ctx->RBuf == 0x1a)
                    break;
                rc = (int)(strlen(ctx->RBuf));

                if (ctx->NMRec > 1) {

                    p = ctx->RBuf + rc;

                    for (i = 1; i < ctx->NMRec; ++i) {
                        cnt = ctx->RBufLen - rc;
                        if (cnt < 10) {
                            printf1(ctx, "Error: exceeded max read buffer length (%d bytes).\n",RLMaxDef);
                            err = -1;
                            break;
                        }
                        if (fgets(p,cnt,ctx->RFd) == NULL) {
                            printf1(ctx, "Warning: input file does not contain a multiple of %d records.\n",ctx->NMRec);
                            printf1(ctx, "Stopped with reading %d physical records from %s.\n",ctx->NRec1,ctx->DFILFN[ctx->DFILNI-1]);
                            break;
                        }
                        cnt = (int)(strlen(p));
                        rc += cnt;
                        p += cnt;
                    }
                    if (i < ctx->NMRec)
                        break;
                }
                ctx->NRec += ctx->NMRec;
                ctx->NRec1 += ctx->NMRec;
            }
            else {
                rc = (int)(fread(ctx->RBuf,sizeof(char),(size_t)(ctx->DRecLen1),ctx->RFd));
                if (rc != ctx->DRecLen1) {
                    if (rc > 0) {
                        printf1(ctx, "Warning: stopped with reading a logical record of less than %d bytes.\n",ctx->DRecLen1);
                        printf1(ctx, "Last physical record(s) (%d bytes) will be ignored.\n",rc);
                    }
                    break;
                }
                ctx->NRec += ctx->NMRec;
                ctx->NRec1 += ctx->NMRec;
            }
            if (check_drec(ctx, ctx->RBuf)) {      /* check for data records */

                p = ctx->RBuf + strlen(ctx->RBuf);
                while (--p >= ctx->RBuf && (*p == LF || *p == CR))
              		    rc--;
                *(ctx->RBuf + rc) = '\0';

                p = ctx->RBuf;
                ctx->RBufALen = rc;                  /* actual buffer length */
                for (i = 0; i < rc; ++i) {
                    if (*p == LF || *p == CR) {
                        if (ctx->ISEPC != '\0')
                            *p = ctx->ISEPC;
                        else
                            *p = ' ';
                    }
                    p++;
                }

                /* if string variables, copy RBuf into RSBuf and clear
                   strings in RBuf */

                if (sflag) {
                    p = ctx->RBuf;
                    q = ctx->RSBuf;
                    for (i = 0; i < rc; ++i)
                        *q++ = *p++;

                    if (ctx->FFMTCnt == 0) {
                        i = idx;
                        while (i >= 0) {
                            if (ctx->VTyp[i] == 1 && ctx->VTypA[i] == 0)
                                clear_str(ctx, ctx->RBuf,ctx->RBufALen,i); 
                            i = ctx->VNxt[i];
                        }
                    }
                }

                /* read required values into VCJVal[] */
    
                if (ctx->FFMTCnt > 0) {      /* fixed format */

                    for (i = 1; i <= ctx->VCJMax; ++i) {
                        if (ctx->VCJFlg[i]) {
                            if (ctx->FFMTC1[i] > rc) {
                                ctx->VCJVal[i] = ctx->MBlnkVal;
                                ctx->MBlnkN++;  
                                continue;
                            }

                            if (ctx->FFMTC2[i] > rc) {
                                err = -1;
                                break;
                            }
                            p = ctx->RBuf + ctx->FFMTC1[i] - 1;
                            l = ctx->FFMTC2[i] - ctx->FFMTC1[i] + 1;
                            while (l > 0) {
                                if (*p != ' ')
                                    break;
                                p++;
                                l--;
                            }
                            if (l <= 0) {
                                ctx->VCJVal[i] = ctx->MBlnkVal;
                                ctx->MBlnkN++;  
                            }
                            else {
                                ctx->VCJVal[i] = dscan(ctx, p,l,&mval);
                                if (check_mval(ctx, mval)) {
                                    err = -1;
                                    break;
                                }
                            }
                        }
                    }
                }
                else {          /* free format */

                    p = ctx->RBuf;
                    for (i = 1; i <= ctx->VCJMax; ++i) { 

                        mval = 0;
                        q = skip_sep(ctx, p);    /* skip separator */

                        if (ctx->VCJFlg[i]) {
                            if (!*q || *q == LF || *q == CR) {

                                if (ctx->ISEPC != '\0' && *q == ctx->ISEPC) {
                                    ctx->VCJVal[i] = ctx->MBlnkVal;
                                    ctx->MBlnkN++;
                                }
                                else {
                                    err = -1;
                                    break;
                                }
                            }
                            else {
                                ctx->VCJVal[i] = dscan(ctx, q,0,&mval);
                                if (check_mval(ctx, mval)) {

                                    err = -1;
                                    break;
                                }
                            }
                        }
                        if (i == ctx->VCJMax)
                            break;

                        p = skip_dval(ctx, q);
                    }
                }
                if (err) {
                    printf1(ctx, "Can't read c%d in data file %s, record %d.\n",i,ctx->DFILFN[ctx->DFILNI-1],ctx->NRec1);
                    break;
                }
                return(1);      /* return this record, values in VCJVal */
            }
        }
        fclose(ctx->RFd);
        ctx->REOF = 1;
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

int skip_df(TDAContext *ctx, int n,int opt)  
{
    int i,fptr;
    double x,y; 

    fptr = (int)ftell(ctx->RFd);
    if (fptr < 0)
        return(-1);

    for (i = 0; i < n; ++i) {
        if (fgets(ctx->RBuf,ctx->RBufLen,ctx->RFd) == NULL || *ctx->RBuf == 0x1a)
            return(-1);

        if (opt) {
            if (get_nsdxy(ctx, ctx->RBuf,&x,&y))  
                return(-1);

            ctx->SDVarXMin = dmin(ctx, x,ctx->SDVarXMin);
            ctx->SDVarXMax = dmax(ctx, x,ctx->SDVarXMax);
            ctx->SDVarYMin = dmin(ctx, y,ctx->SDVarYMin);
            ctx->SDVarYMax = dmax(ctx, y,ctx->SDVarYMax);
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

int get_nsdxy(TDAContext *ctx, char *buf,double *x,double *y)
{
    char *p;

    p = skip_b(ctx, buf);
    if (sscanf(p,"%lg",x) != 1)
        return(-1);

    p = skip_dbl(ctx, p);
    p = skip_b(ctx, p);
    if (sscanf(p,"%lg",y) != 1)
        return(-1);
    return(0);
}

/*--###---------------------------------------------------------------------*/
/*  check_mval(mval)        check missing values.                           */
          
int check_mval(TDAContext *ctx, int mval)
{
    if (mval == 0)
        return(0);
    if (mval == 1)
        ctx->MBlnkN++;
    else if (mval == 2)
        ctx->MStarN++;
    else if (mval == 3)
        ctx->MPntN++;
    else if (ctx->MGenFlg)
        ctx->MGenN++;
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

int sdnvar(TDAContext *ctx)
{
    register int i;
    int err,sflag;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Reading a spatial data file. Current memory: %d bytes.\n",ctx->MemReq);
    if (ctx->DMDef) {
        printf1(ctx, "Error: a data matrix already exists.\n");
        return(-1);
    }
    sdnvar_close(ctx);         /* close previous file (if any) */

    ctx->SDVarNT = ctx->SDVarNP = ctx->SDVarNL = ctx->SDVarNPol = ctx->SDVarNU = 0;
    ctx->SDVarXMin = ctx->SDVarYMin =  (ctx->DBLMAX / 10.0);
    ctx->SDVarXMax = ctx->SDVarYMax = -(ctx->DBLMAX / 10.0);

    ctx->DFILN = ctx->MatchNV = ctx->VFmtFlg = ctx->NVNum = 0;

    ctx->DRecLen = 0;            /* variable record length */
    ctx->NMRec = 1;              /* logical = physical record */
    ctx->DBLKVar = -1;           /* no block mode */
    ctx->AVDFN = -1;             /* no archive file */
    ctx->ARCDICFlg = 0;          /* set by arcdic option */
    ctx->RDMNInit = 0;           /* init random number generator */
    ctx->ISEPC = '\0';           /* separation character */
    ctx->VTNOC = 0;
    ctx->BSIZE = 0;

    ctx->EDVALSn = ctx->EDVALOrg = ctx->EDVALDes = 0;
    ctx->EDVALTs = ctx->EDVALTf = 0.0;
         
    ctx->NVIdx = get_nidx(ctx);     /* index to first new variable */

    /* set default missing value codes */

    ctx->MStarVal  = -1.0; 
    ctx->MPntVal   = -1.0; 
    ctx->MBlnkVal  = -1.0; 
    ctx->MGenVal   = -1.0; 
    ctx->MMatchVal = -3.0; 
    ctx->GDNRec = ctx->GDKFlg = ctx->GDFlg = ctx->MGenFlg = ctx->MStarN = ctx->MPntN = ctx->MBlnkN = ctx->MGenN = 0; 

    ctx->ISelPtr  = NULL;        /* set by isel option */
    ctx->VSelPtr  = NULL;        /* set by vsel option */
    ctx->BSelPtr  = NULL;        /* set by bsel option */
    ctx->GDFNPtr  = NULL;        /* file name by df option */
    ctx->KeepPtr  = NULL;            
    ctx->DropPtr  = NULL;            

    ctx->NVArc = ctx->NVArc1 = 0;     /* number of archive variables, counted by save_var(ctx) */

    if (check_nvsd(ctx))       /* check command etc */
        goto NVSDFin;

    if (ctx->NVNum == 0) {
        printf1(ctx, "No variables defined.\n");
        goto NVSDFin;
    }
    if (ctx->DFILN == 0) {
        printf1(ctx, "No data file defined.\n");
        goto NVSDFin;
    }
    else if (ctx->DFILN > 1) {
        printf1(ctx, "Error: can only use a single data file.\n");
        goto NVSDFin;
    }
    newline(ctx);      

    if (ctx->VLabelLen > 0 && ctx->VLabelLen < 8)      
        ctx->VLabelLen = 8;

    sflag = 0;
    i = ctx->NVIdx;                  /* count string variables */
    while (i >= 0) {
        if (ctx->VTyp[i] == 1 && ctx->VTypA[i] == 0)
            sflag++;
        i = ctx->VNxt[i];
    }
    if (ctx->VFmtFlg)
        make_vfmt(ctx, ctx->NVIdx);   /* make new print formats if requested with fmt */

    prn_var(ctx, ctx->NVIdx);         /* print list of new variables */

    /* check for standard variables */

    if (check_sd(ctx, 1))
        goto NVSDFin;

    ctx->NOCMaxA = ctx->VTNOC;
    if (ctx->NOCMaxA <= 0)  
        ctx->NOCMaxA = NOCDef;

    printf1(ctx, "\nCreating a new data matrix.\n");
    printf1(ctx, "Maximum number of cases: %d\n",ctx->NOCMaxA);

    if (check_vcj(ctx, ctx->NVIdx,1))             /* check cj references */
        goto NVSDFin;

    printf1(ctx, "\nUsing data file(s): %s",ctx->DFILFN[0]);
    for (i = 1; i < ctx->DFILN; ++i)
        printf1(ctx, ",%s",ctx->DFILFN[i]);
    printf1(ctx, "\n");
    if (ctx->DRecLen > 0)
        printf1(ctx, "Fixed record length: %d\n",ctx->DRecLen);

    if (ctx->FFMTCnt > 0) {
        if (check_ffmt(ctx, 1))
            goto NVSDFin;
    }
    else {
        printf1(ctx, "Free format. Separation character(s): ");
        if (ctx->ISEPC != '\0')
            printf1(ctx, "%04x [hex]\n",ctx->ISEPC);
        else 
            printf1(ctx, "default.\n");
    }
    if (ctx->VTNOC > 0)  
        printf1(ctx, "Reading maximal %d records.\n",ctx->VTNOC);

    /* allocate memory for new variables */

    if (alloc_vdat(ctx, ctx->NVIdx,1)) {
        printf1(ctx, "Insufficient memory for new variables.\n");
        goto NVSDFin;
    }
    err = gen_dm(ctx, ctx->NVIdx,ctx->NVNum,sflag,1);      /* create data */
    if (err)  
        goto NVSDFin;

    prn_mval(ctx);                           /* info about missing values */

    printf1(ctx, "\nNumber of points: %d\n",ctx->SDVarNP);
    printf1(ctx, "Number of lines: %d\n",ctx->SDVarNL);
    printf1(ctx, "Number of polygons: %d\n",ctx->SDVarNPol);
    printf1(ctx, "Number of unknown objects: %d\n",ctx->SDVarNU);

    printf1(ctx, "\nMaximal number of points in spatial objects: %d\n",ctx->SDVarMax);
    printf1(ctx, "Total number of points in spatial objects: %d\n",ctx->SDVarNT);
    if (ctx->SDVarMax > 0) {
        if (!(ctx->SDVarFd = fopen(ctx->DFILFN[0],OPEN_RD))) {
            printf1(ctx, "Error: cannot re-open data file: %s\n",ctx->DFILFN[0]);
            err = -1;
            goto NVSDFin;
        }
        ctx->SDVarFDef = ctx->SDVarDef = 1;  

        printf1(ctx, "Note: %s remains opened for further access.\n\n",ctx->DFILFN[0]);

        printf1(ctx, "X values. Minimum: %14.4lf  Maximum: %14.4lf\n",ctx->SDVarXMin,ctx->SDVarXMax);
        printf1(ctx, "Y values. Minimum: %14.4lf  Maximum: %14.4lf\n",ctx->SDVarYMin,ctx->SDVarYMax);

        if (sdnvar_alloc(ctx, 1,ctx->SDVarMax + 20)) {
            printf1(ctx, "\nError: Insufficient memory for spatial objects.\n");
            err = -1;
        }
    }
    else
        err = -1;

NVSDFin:
    alloc_vl(ctx, 0);

    check_vcj(ctx, ctx->NVIdx,0);
    check_ffmt(ctx, 0);

    if (err) {
        if (ctx->NVNum > 0)
            clear_avar(ctx, ctx->NVIdx);
        printf1(ctx, "No new variables created. ");
        sdnvar_close(ctx);
    }
    else {
        ctx->WIVar = -1;
        ctx->WSum = (double)ctx->NOC;
        ctx->WSumS = 0.0;
        ctx->WNorm = 1.0;
        ctx->WNormFlag = 0;
        printf1(ctx, "\nEnd of creating new variables. ");
    }
    prn_mem(ctx);
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  check_nvsd  Check sdnvar command in CmdBuf.                             */
/*              Return 0 if OK, -1 if error.                                */

int check_nvsd(TDAContext *ctx)
{
    int err,n,l,m1,m2;
    register char c,*p,*q;
    char sc;
    double x;

    ctx->NVSEPC = ' ';             
    ctx->NVL0 = ctx->VTNOC = ctx->BSIZE = 0;
    err = -1;
    p = ctx->CmdBuf + 6;
    while (*++p) {
        if (sscanf(p,"noc=%d",&n) == 1 && n > 0) {
            ctx->VTNOC = n;
            p = skip_int(ctx, p + 4);
        }
        else if (sscanf(p,"dreclen=%d",&n) == 1 && n > 0) {
            ctx->DRecLen = n;
            p = skip_int(ctx, p + 8);
        }
        else if (sscanf(p,"fmt=%d.%d",&ctx->VFmt1,&ctx->VFmt2) == 2) {
            ctx->VFmtFlg = 1;
            p = skip_int(ctx, p + 4);
            p = skip_int(ctx, p + 1);
        }
        else if (sscanf(p,"mstar=%lf",&x) == 1) {
            ctx->MStarVal = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"mpnt=%lf",&x) == 1) {
            ctx->MPntVal = x;
            p = skip_dbl(ctx, p + 5);
        }
        else if (sscanf(p,"mblnk=%lf",&x) == 1) {
            ctx->MBlnkVal = x;
            p = skip_dbl(ctx, p + 6);
        }
        else if (sscanf(p,"mgen=%lf",&x) == 1) {
            ctx->MGenVal = x;
            ctx->MGenFlg = 1;
            p = skip_dbl(ctx, p + 5);
        }
        else if (!strncmp(p,"dfile=",6)) {
            if (ctx->DFILN >= DFILMax) {
                printf1(ctx, "Exceeded max number of data files.\n");
                goto CNVSDFin;
            }
            ctx->DFILFN[ctx->DFILN++] = p + 6;
            p = skip_nc(ctx, p);
        }
        else if (!strncmp(p,"break=",6)) {
            ctx->BRKPtr = p + 6;
            p = skip_expr(ctx, p + 6);
        }
        else if (sscanf(p,"isc=%c",&sc) == 1) {
            if (sc == 't')
                sc = '\t';
            ctx->ISEPC = sc;
            p += 5;                
        }
        else if (!strncmp(p,"ffmt=c",6)) {
            ctx->FFMTCnt = 0;       
            q = p + 5;
            while (1) {
                if ((sscanf(q,"c%d(%d-%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d,%d)",&n,&m1,&m2) == 3 && n >= 1 && m1 >= 1 && m2 >= m1) ||
                    (sscanf(q,"c%d(%d)",&n,&m1) == 2 && n >= 1 && m1 >= 1)) {
                    q = skip_int(ctx, q + 1);
                    q = skip_blev(ctx, q);
                }
                else {
                    prn_nve(ctx, p);
                    goto CNVSDFin;
                }
                ctx->FFMTCnt++;
                if (*q != ',' || *(q + 1) != 'c')
                    break;
                q++;       
            }
            ctx->FFMTPtr = p + 5;
            p = q;
        }
        else if ((l = get_vnlen(ctx, p)) > 0) {       /* variable */
            q = skip_nc(ctx, p + l);
            c = *q;
            *q = '\0';
            if (save_var(ctx, p,0))
                goto CNVSDFin;
            ctx->NVNum++;
            *q = c;
            p = q; 
        }
        if (*p != ',' && *p != ')') {
            prn_nve(ctx, p);
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

int check_sd(TDAContext *ctx, int opt)
{
    if (ctx->VTyp[ctx->SDVarSDID]   < 2 || ctx->VTyp[ctx->SDVarSDID]  > 3 ||
        ctx->VTyp[ctx->SDVarSDPtr] != 3 ||
        ctx->VTyp[ctx->SDVarSDTyp]  < 2 || ctx->VTyp[ctx->SDVarSDTyp] > 3 ||
        ctx->VTyp[ctx->SDVarSDN]    < 2 || ctx->VTyp[ctx->SDVarSDN]   > 3 ||  
        strcmp(ctx->VName[ctx->SDVarSDID],"SDID")  || 
        strcmp(ctx->VName[ctx->SDVarSDTyp],"SDTyp") ||
        strcmp(ctx->VName[ctx->SDVarSDN],"SDN") ||
        strcmp(ctx->VName[ctx->SDVarSDPtr],"SDPtr")) {
        if (opt)
            printf1(ctx, "\nError: need valid definition of standard variables.\n");
        return(-1);
    }
    if (ctx->VSLen[ctx->SDVarSDPtr] != 5) {
        if (opt)
            printf1(ctx, "\nError: SDPtr must have storage type 5.\n");
        return(-1);
    }
    return(0);
}

/* -##--------------------------------------------------------------------- */
/*  sdnvar_close()  If SDVarDef != 0 close the file.                        */

void sdnvar_close(TDAContext *ctx)
{
    sdnvar_alloc(ctx, 0,0);
    if (ctx->SDVarFDef)
        fclose(ctx->SDVarFd);
    ctx->SDVarNT = ctx->SDVarMax = ctx->SDVarFDef = ctx->SDVarDef = 0;
}

/* -##--------------------------------------------------------------------- */
/*  sdnvar_alloc(opt,n)                                                     */
/*                                                                          */
/*  If opt != 0 allocate standard array for spatial objects with n          */
/*  elements, otherwise free previously allocated memory.                   */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdnvar_alloc(TDAContext *ctx, int opt,int n)
{

    if (ctx->SDVarN > 0) {
        free((char *)ctx->SDVarX);
        free((char *)ctx->SDVarY);
        memrq(ctx, -2 * ctx->SDVarN,sizeof(double));
        ctx->SDVarN = 0;
    }
    if (opt && n > 0) {
        if (!(ctx->SDVarX = (double *)calloc((size_t)(n),sizeof(double)))) {
            return(-1);
        }
        if (!(ctx->SDVarY = (double *)calloc((size_t)(n),sizeof(double)))) {
            free((char *)ctx->SDVarX);
            return(-1);
        }
        memrq(ctx, 2 * n,sizeof(double));
        ctx->SDVarN = n;
    }
    return(0);
}


