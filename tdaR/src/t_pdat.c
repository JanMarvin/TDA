/****************************************************************************/
/*  t_pdat                                                                  */
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
#include "t_gdat.h"  
#include "t_var.h"  
#include "t_gf.h"  
#include "t_alloc.h"  
#include "t_sort.h"  
#include "t_sd.h"  
#include "tda_context.h"

/*  functions in t_pdat.c */

int pdata(TDAContext *ctx);
void dtda(TDAContext *ctx, char *fname,int noc,int nv,short *vidx,int nq);
void dspss(TDAContext *ctx, char *fname,int nv,short *vidx,int nq);
int pdatr(TDAContext *ctx);
int pdatd(TDAContext *ctx);
void pdatdr(TDAContext *ctx, int n,double *x);
double pdatd1(TDAContext *ctx, int opt,int n,double *x,double *y,int nw,double *w);


/* ------------------------------------------------------------------------ */
/*  pdata()     pdata command.                                              */
/*                                                                          */
/*          pdata (                                                         */
/*              keep=varlist,                                               */
/*              drop=varlist,                                               */
/*              sort=varlist,                                               */
/*              noc=...,        only first noc cases                        */
/*              nn = n1,n2      range of cases                              */
/*              transp          transpose data                              */
/*              nc = ...        line feed after nc variables                */
/*              nq = ...        line feed after nq records, def. 1          */
/*              ap=1,           append                                      */
/*              prn=...,        output options, def. 0                      */
/*                              1 triangle                                  */
/*                              2 square matrix                             */
/*                              3 write data as single column               */
/*              dtda=...,       tda description file                        */
/*              dspss=...,      SPSS description file                       */
/*              sepc=' ',       separation character, def ' '               */
/*                              also possible: sepc=none                    */
/*              l0=...,         1 if leading zeros, def. 0                  */
/*                                                                          */
/*              sd=...,         if 1 write s-data records, def. 0           */
/*                              (only with prn=0)                           */
/*              sdfmt=...,      print format for s-data records, def. 12.4  */
/*          ) = fname;          output file                                 */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */
 
int pdata(TDAContext *ctx)
{
    register int i,j,k,ii;
    int err,nrec,n,nn,nv,na,nb,sflag,lcnt,w,d,r,sdrec;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1)) {
        p_clean(ctx);
        return(-1);
    }
    if (ctx->PMSD) {
        if (check_sd(ctx, 0)) {
            printf1(ctx, "Error: no spatial data defined.\n");
            goto PDATFin;
        }
        if (ctx->PMPRNO != 0) {
            printf1(ctx, "Error: sd option requires prn=0.\n");
            goto PDATFin;
        }
    }
    if (ctx->PMKeep && ctx->PMDrop) {
        p_err(ctx, -16,1);
        goto PDATFin;
    }
    if (ctx->PMNQ < 1)
        ctx->PMNQ = 1;

    /* get list of variables in AcNS[] */

    if (alloc_acns(ctx, imax(ctx, ctx->NVAR,ctx->PMNV)))
        goto PDATFin;

    nv = 0;
    if (ctx->PMKeep) {
        for (k = 0; k < ctx->PMNV; ++k)           
            ctx->AcNS[nv++] = ctx->PMVIdx[k];              
    }
    else {
        j = ctx->VIFirst;
        while (j >= 0) {
            i = 1;
            if (ctx->PMDrop) {
                for (k = 0; k < ctx->PMNV; ++k) {
                    if (ctx->PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                ctx->AcNS[nv++] = (short)(j);
            j = ctx->VNxt[j];
        }
    }
       
    if (nv == 0) {
        printf1(ctx, "No variables selected.\n");
        goto PDATFin;
    }
    if (ctx->PMPRNO == 1 || ctx->PMPRNO == 2) {
        if (nv != 1) {
            printf1(ctx, "Error: prn options can only be used with a single variable.\n");
            goto PDATFin;
        }
    }

    sflag = 0;
    if (ctx->PM1NV > 0) {        /* sort */
        err = vsort(ctx, ctx->PM1NV,ctx->PM1VIdx,1,0,1);
        if (err)
            goto PDATFin;
        sflag = 1;
    }
    na = 0; 
    nb = ctx->NOC;
    if (ctx->PMNNFlg) {
        na = ctx->PMNN1 - 1;
        nb = ctx->PMNN2;
    }
    else if (ctx->PMNOCFlg)
        nb = ctx->PMNOC;
    if (na < 0)
        na = 0;
    if (nb > ctx->NOC)
        nb = ctx->NOC;

    nn = 0;
    if (ctx->PMPRNO == 2) {
        nn = (int)sqrt((double)(nb - na));
        if (nn * nn != nb - na) {
            printf1(ctx, "Error: number of cases not compatible with prn=2 option.\n");
            goto PDATFin;
        }
    }
    if (ctx->PML0 != 0)
        ctx->PML0 = -1;

    if (ctx->PML0 != 0 || ctx->XSEPC != ctx->SEPC) {        /* change print format */
        for (k = 0; k < nv; ++k) {
            j = ctx->AcNS[k];
            if (ctx->VTyp[j] == 1)
                continue;
            w = (int)ctx->VPFmt1[j];             
            d = (int)ctx->VPFmt2[j];             
            makefmt(ctx, &w,&d,ctx->VPFmtS[j],VPFmtSLen,0,ctx->XSEPC,ctx->PML0);
        }
    }

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing: %s\n",ctx->PMFdName);

    sdrec = nrec = 0;

#ifdef TDA_R_PACKAGE
    /*  The variables pdata is about to print, in print order: name and
        whether it is a string.  Without this the R side had to read the
        order and the names back out of the printed variable table, so a
        frame built from the exports put its string columns in the wrong
        place and named them from the console.  One row per variable,
        "<name>|<0 numeric, 1 string>".  */
    {
        int kk;
        char vrow[96];
        for (kk = 0; kk < nv; ++kk) {
            int jj = ctx->AcNS[kk];
            snprintf(vrow, sizeof(vrow), "%s|%d", ctx->VName[jj],
                     ctx->VTyp[jj] == 1 ? 1 : 0);
            tda_export_str_row(ctx, "pdata.vars", vrow);
        }
    }
#endif

    n = 0;
    if (ctx->PMTransp == 0) {
        lcnt = 0;
        for (i = na; i < nb; ++i) {
            ii = i;
            if (sflag)
                ii = ctx->VSORTPtr[i];

            if (ctx->PMPRNO == 1 || ctx->PMPRNO == 2) {
                j = ctx->AcNS[0];       
                if (ctx->VTyp[j] == 1)
                    continue;
                tmp = get_data(ctx, j,ii);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[j],tmp);
                if (ctx->PMPRNO == 1) {
                    if (++n > nn) {
                        fprintf(ctx->PMFd,"\n");
                        nrec++;
                        prn_message(ctx, nrec,0,1);
                        n = 0;
                        nn++;
                    }
                }
                else if (++n >= nn) {
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                    prn_message(ctx, nrec,0,1);
                    n = 0;
                }
            }
            else {              /* standard printing */
                n = 0;
                for (k = 0; k < nv; ++k) {
                    j = ctx->AcNS[k];
                    if (ctx->VTyp[j] != 1) {
                        tmp = get_data(ctx, j,ii);                  
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[j],tmp);
#ifdef TDA_R_PACKAGE
                        /*  The numeric counterpart of pdata.str.<Var>
                            below.  Without it a reader built on pdata
                            still had to parse the printed file back --
                            tda_read_spss did, and could not produce its
                            frame when the file was unavailable.  One
                            cell per call, closed into rows by the
                            endrow after the variable loop.  */
                        tda_export_cell(ctx, "pdata.values", tmp);
#endif
                    }
                    else {          /* string variables */
                        get_str(ctx, ctx->SVBuf,j,ii);
                        fprintf(ctx->PMFd,"%s",ctx->SVBuf);
#ifdef TDA_R_PACKAGE
                        /*  A string cannot travel on the numeric export
                            channel, and reading it back out of the file
                            means re-parsing fixed-width text -- exactly
                            what the exports exist to avoid.  Each value
                            is handed over as itself, in the order the
                            row is printed, under one key per variable so
                            the R side does not have to work out which
                            column it came from.  */
                        {
                            char ekey[80];
                            snprintf(ekey, sizeof(ekey), "pdata.str.%s",
                                     ctx->VName[j]);
                            tda_export_str_row(ctx, ekey, ctx->SVBuf);
                        }
#endif
                        if (ctx->XSEPC)
                            fprintf(ctx->PMFd,"%c",ctx->XSEPC);
                    }
                    if (k == nv - 1)
                        break;

                    if (ctx->PMPRNO == 3 || (ctx->PMNC > 0 && ++n >= ctx->PMNC)) {
                        fprintf(ctx->PMFd,"\n");
                        nrec++;
                        prn_message(ctx, nrec,0,1);
                        n = 0;
                    }
                }
#ifdef TDA_R_PACKAGE
                /* one exported row per printed case */
                tda_export_endrow(ctx, "pdata.values");
#endif
                if (++lcnt >= ctx->PMNQ) {
                    fprintf(ctx->PMFd,"\n");
                    nrec++;
                    prn_message(ctx, nrec,0,1);
                    lcnt = 0;
                }
                if (ctx->PMSD) {
                    if ((r = sd_getdata(ctx, ii,0,0,0)) < 1) {
                        printf1(ctx, "Error: cannot read s-data records.\n");
                        goto PDATFin;
                    }
                    for (j = 0; j < r; ++j) {
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMSDFmtS,ctx->SDVarX[j]);
                        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMSDFmtS,ctx->SDVarY[j]);
                        fprintf(ctx->PMFd,"\n");
                        sdrec++;
                    }
                }
            }
        }
        if (ctx->PMNQ > 1 && lcnt > 0 && lcnt < ctx->PMNQ) {
            fprintf(ctx->PMFd,"\n");
            nrec++;
            prn_message(ctx, nrec,0,1);
        }
    }
    else {
        for (k = 0; k < nv; ++k) {
            j = ctx->AcNS[k];
            if (ctx->VTyp[j] == 1)               /* string variables */
                continue;

            for (i = na; i < nb; ++i) {
                ii = i;
                if (sflag)
                    ii = ctx->VSORTPtr[i];
                tmp = get_data(ctx, j,ii);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[j],tmp);
            }
            fprintf(ctx->PMFd,"\n");
            nrec++;
            prn_message(ctx, nrec,0,1);
        }
        nv = nb - na;
    }
    prn_message(ctx, nrec,1,1);

    if (ctx->PMPRNO == 3)
        printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    else {
        printf1(ctx, "%d records with %d variables written to: %s\n",nrec,nv,ctx->PMFdName);
        if (ctx->PMSD)
            printf1(ctx, "Added %d s-data records.\n",sdrec);
    }
    if (ctx->PMTransp == 0 && ctx->PMPRNO == 0) {

        if (ctx->PMTDAFDef)
            dtda(ctx, ctx->PMFdName,nrec,nv,ctx->AcNS,ctx->PMNQ);

        if (ctx->PMSPSSFDef)
            dspss(ctx, ctx->PMFdName,nv,ctx->AcNS,ctx->PMNQ);
    }
    if (ctx->PML0 != 0 || ctx->XSEPC != ctx->SEPC) {        /* reset format */
        for (k = 0; k < nv; ++k) {
            j = ctx->AcNS[k];
            w = (int)ctx->VPFmt1[j];             
            d = (int)ctx->VPFmt2[j];             
            makefmt(ctx, &w,&d,ctx->VPFmtS[j],VPFmtSLen,0,ctx->SEPC,0);
        }
    }
    err = 0;

PDATFin:
    if (ctx->PM1NV > 0)                     
        vsort(ctx, 0,ctx->PM1VIdx,0,0,1);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dtda(fname,noc,nv,vidx,nq)                                              */
 
void dtda(TDAContext *ctx, char *fname,int noc,int nv,short *vidx,int nq)
{
    register int i,j,n,ni;

    fprintf(ctx->PMTDAFd,"nvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",noc);
    ni = 1;
    for (n = 1; n <= nq; ++n) {
        for (i = 0; i < nv; ++i) {
            j = vidx[i];
            fprintf(ctx->PMTDAFd,"  %s",ctx->VName[j]);
            if (nq > 1)
                fprintf(ctx->PMTDAFd,"%d",n);                                                 

            fprnchar(ctx, ctx->PMTDAFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]),0);
            fprintf(ctx->PMTDAFd,"<%d>[%d.%d]",ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j]);

            if (ctx->VLabel[j] != NULL)                                                      
                fprintf(ctx->PMTDAFd,"(%s)",ctx->VLabel[j]);                                           
            fprintf(ctx->PMTDAFd," = c%d,\n",ni++);
        }
    }
    fprintf(ctx->PMTDAFd,");\n");
    
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
}

/* ------------------------------------------------------------------------ */
/*  dspss(fname,noc,nv,vidx,nq)                                             */
 
void dspss(TDAContext *ctx, char *fname,int nv,short *vidx,int nq)
{
    register int i,j,n,nn,ii;

    fprintf(ctx->PMSPSSFd,"DATA LIST FILE='%s' FREE/\n",fname);                                              
    ii = 0;
    for (nn = 1; nn <= nq; ++nn) {
        for (i = 0; i < nv; ++i) {                                                      
            j = vidx[i];                                                                
            fprintf(ctx->PMSPSSFd," %s",ctx->VName[j]);                                                 
            if (nq > 1)
                fprintf(ctx->PMSPSSFd,"%d",nn);                                                 
            ii++;
            n = ii / 8;
            if (8 * n == ii) {
                fprintf(ctx->PMSPSSFd,"\n");
                ii = 0;
            }     
        }
    }
    fprintf(ctx->PMSPSSFd,".\nVARIABLE LABELS\n");
    for (nn = 1; nn <= nq; ++nn) {
        for (i = 0; i < nv; ++i) {
            j = vidx[i];
            fprintf(ctx->PMSPSSFd," %s",ctx->VName[j]);
            if (nq > 1)
                fprintf(ctx->PMSPSSFd,"%d",nn);                                                 
            fprnchar(ctx, ctx->PMSPSSFd,' ',ctx->VNameLen - (int)strlen(ctx->VName[j]),0);
            if (ctx->VLabel[j] != NULL)                                                      
                fprintf(ctx->PMSPSSFd," '%s'",ctx->VLabel[j]);                                           
            else
                fprintf(ctx->PMSPSSFd," '%s'",ctx->VName[j]);                                           
            if (i == nv - 1 && nn == nq)
                fprintf(ctx->PMSPSSFd,".\n");
            else
                fprintf(ctx->PMSPSSFd,"/\n");
        }
    }                                                                               
    printf1(ctx, "SPSS description written to: %s\n",ctx->PMSPSSFName);
} 

/* -##--------------------------------------------------------------------- */
/*  pdatr() print data for all pairs (i,j) of cases.                        */
/*                                                                          */
/*          pdatr (                                                         */
/*              keep = varlist,                                             */
/*              drop = varlist,                                             */
/*              nfmt = ...,     print format for i,j, def. 4                */  
/*          ) = fname;          output file                                 */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */
 
int pdatr(TDAContext *ctx)
{
    register int i,j,k,jv;
    int err,nrec,nv;
    double tmp;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1)) {
        p_clean(ctx);
        return(-1);
    }
    if (ctx->PMKeep && ctx->PMDrop) {
        p_err(ctx, -16,1);
        goto PDATRFin;
    }

    /* get list of variables in AcNS[] */

    if (alloc_acns(ctx, imax(ctx, ctx->NVAR,ctx->PMNV)))
        goto PDATRFin;

    nv = 0;
    if (ctx->PMKeep) {
        for (k = 0; k < ctx->PMNV; ++k)           
            ctx->AcNS[nv++] = ctx->PMVIdx[k];              
    }
    else {
        j = ctx->VIFirst;
        while (j >= 0) {
            i = 1;
            if (ctx->PMDrop) {
                for (k = 0; k < ctx->PMNV; ++k) {
                    if (ctx->PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                ctx->AcNS[nv++] = (short)(j);
            j = ctx->VNxt[j];
        }
    }
       
    if (nv == 0) {
        printf1(ctx, "No variables selected.\n");
        goto PDATRFin;
    }
    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing: %s\n",ctx->PMFdName);

    nrec = 0;
    for (i = 0; i < ctx->NOC; ++i) {
        for (j = 0; j < ctx->NOC; ++j) {
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i+1);  
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,j+1);  

            for (k = 0; k < nv; ++k) {
                jv = ctx->AcNS[k];
                if (ctx->VTyp[jv] != 1) {
                    tmp = get_data(ctx, jv,i);                  
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[jv],tmp);
                }
                else {          /* string variables */
                    get_str(ctx, ctx->SVBuf,jv,i);
                    fprintf(ctx->PMFd,"%s",ctx->SVBuf);
                    if (ctx->XSEPC)
                        fprintf(ctx->PMFd,"%c",ctx->XSEPC);
                }
            }
            for (k = 0; k < nv; ++k) {
                jv = ctx->AcNS[k];
                if (ctx->VTyp[jv] != 1) {
                    tmp = get_data(ctx, jv,j);                  
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[jv],tmp);
                }
                else {          /* string variables */
                    get_str(ctx, ctx->SVBuf,jv,j);
                    fprintf(ctx->PMFd,"%s",ctx->SVBuf);
                    if (ctx->XSEPC)
                        fprintf(ctx->PMFd,"%c",ctx->XSEPC);
                }
            }
            fprintf(ctx->PMFd,"\n");
            nrec++;
        }
    }
    prn_message(ctx, nrec,1,1);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    err = 0;

PDATRFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  pdatd() Create a distance matrix.                                       */
/*                                                                          */
/*          pdatd (                                                         */
/*              v=...,      optional varlist, def. all variables            */
/*              opt=...,    option, def. 1                                  */
/*                           1 euclidean distances                          */
/*                           2 city-block distance                          */
/*                           3 number of different values                   */
/*                           4 dissimilarity index                          */

/*              wt=...,      optional weights, def. none                    */
/*              fmt = ...,   print format, def. 10.4                        */  
/*          ) = fname;          output file                                 */
/*                                                                          */
/*          Return 0 if OK, -1 if error.                                    */
 
int pdatd(TDAContext *ctx)
{
    register int i,j,k;
    int err,nv;
    double d;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Create a distance matrix. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1)) {
        p_clean(ctx);
        return(-1);
    }
    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    if (ctx->PMOPT == 2)  
        printf1(ctx, "City-block distance.\n");
    else if (ctx->PMOPT == 3)  
        printf1(ctx, "Number of different values.\n");
    else if (ctx->PMOPT == 4)  
        printf1(ctx, "Dissimilarity index.\n");
    else {
        ctx->PMOPT = 1;    
        printf1(ctx, "Euclidean distance.\n");
    }   
    if (ctx->PMNV > 0)
        nv = ctx->PMNV;
    else
        nv = ctx->NVAR;

    /*  nv is 0 when the command runs with no data matrix -- pdatd() on
        a bare session, which is exactly what examples/tests/pdatd.cf
        does.  alloc_aci(0) leaves AcI NULL and the "Variables: ..."
        line below dereferenced it: a NULL read and SIGSEGV, which in
        this package takes the R session with it.  The 491-case suite
        did not notice, because the reference file records the output
        printed BEFORE the crash and check.py compares only that.  */
    if (nv < 1) {
        printf1(ctx, "Error: no data matrix.\n");
        goto PDATDFin;
    }
    if (alloc_aci(ctx, nv))
        goto PDATDFin;

    if (ctx->PMNV > 0) {
        for (j = 0; j < ctx->PMNV; ++j)
            ctx->AcI[j] = ctx->PMVIdx[j];
    }
    else { 
        j = 0;
        i = ctx->VIFirst;
        while (i >= 0) {
            ctx->AcI[j++] = i;
            i = ctx->VNxt[i];
        }
    }
    printf1(ctx, "Variables: %s",ctx->VName[ctx->AcI[0]]);
    for (k = 1; k < nv; ++k)  
        printf1(ctx, ", %s",ctx->VName[ctx->AcI[k]]);
    newline(ctx);
    newline(ctx);

    if (ctx->PMNTP > 0) {
        if (ctx->PMNTP != nv) {
            printf1(ctx, "Error: would need %d weights.\n",nv);
            goto PDATDFin;
        }
    }
    if (alloc_acx(ctx, ctx->NOC * ctx->NOC + 1))
        goto PDATDFin;

    if (alloc_acu(ctx, nv + 1))
        goto PDATDFin;
    if (alloc_acv(ctx, nv + 1))
        goto PDATDFin;

    for (i = 1; i < ctx->NOC; ++i) {
        for (k = 0; k < nv; ++k)  
            ctx->AcU[k] = get_data(ctx, ctx->AcI[k],i);

        if (ctx->PMOPT == 4)         /* change into relative values */
            pdatdr(ctx, nv,ctx->AcU);

        for (j = 0; j < i; ++j) {
            for (k = 0; k < nv; ++k)  
                ctx->AcV[k] = get_data(ctx, ctx->AcI[k],j);

            if (ctx->PMOPT == 4)                                               
                pdatdr(ctx, nv,ctx->AcV);

            ctx->AcX[i * ctx->NOC + j] = pdatd1(ctx, ctx->PMOPT,nv,ctx->AcU,ctx->AcV,ctx->PMNTP,ctx->PMTP);
        }
    }

    for (i = 0; i < ctx->NOC; ++i) {
        for (j = 0; j < ctx->NOC; ++j) {
            if (j <= i)
                d = ctx->AcX[i * ctx->NOC + j];
            else           
                d = ctx->AcX[j * ctx->NOC + i];
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,d);  
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "pdatd.table", d);
#endif
        }
        fprintf(ctx->PMFd,"\n");  
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "pdatd.table");
#endif
    }
    printf1(ctx, "%d records written to: %s\n",ctx->NOC,ctx->PMFdName);
    err = 0;

PDATDFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  pdatdr(n,x)     Change x into relative values                           */

void pdatdr(TDAContext *ctx, int n,double *x)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    double s = 0.0;

    for (i = 0; i < n; ++i)
        s += x[i];
    if (s != 0.0) {
        for (i = 0; i < n; ++i)
            x[i] /= s;
    }
}

/* -##--------------------------------------------------------------------- */
/*  pdatd1(opt,n,x,y,nw,w)                                                  */
/*                                                                          */
/*  Return distance of n-vectors x and y, depending on opt.                 */
/*  opt = 1 : euclidean distance                                            */
/*        2 : city-block distance                                           */
/*        3 : number of different values                                    */
/*        4 : dissimilarity index (x and y already contain rel. values)     */
/*                                                                          */
/*  If nw > 0 use weights w[k], k = 0,1,...,n-1                             */
  
double pdatd1(TDAContext *ctx, int opt,int n,double *x,double *y,int nw,double *w)
{
    register int i;
    double tmp,d,wt;

    wt = 1.0;
    d = 0.0;
    for (i = 0; i < n; ++i) {
        if (nw > 0)
            wt = w[i];

        tmp = x[i] - y[i];

        if (opt == 1)  
            d += wt * tmp * tmp;
        else if (opt == 2 || opt == 4)
            d += wt * fabs(tmp);
        else if (opt == 3) {
            if (fabs(tmp) > ctx->EPSI1)
                d += wt;
        }
    }
    if (opt == 1)
        d = sqrt(d);
    else if (opt == 4)
        d /= 2.0;       

    return(d);
}







