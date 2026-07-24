/****************************************************************************/
/*  t_mdat                                                                  */
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
#include "t_gdat.h"  
#include "t_var.h"  
#include "t_gf.h"  
#include "t_edat.h"  
#include "t_seq.h"  
#include "t_alloc.h"  
#include "t_gdd.h"      
#include "t_mat.h"      
#include "t_matf.h"      
#include "t_matc.h"      
#include "t_rzoo.h"
#include "t_sd.h"
#include "tda_context.h"

/*  functions in t_mdat.c */

int clear(TDAContext *ctx);
void clear_a(TDAContext *ctx);
int clearnl(TDAContext *ctx);
int tsel(TDAContext *ctx);       
void tsel_off(TDAContext *ctx, int opt);       
int cwt(TDAContext *ctx);       
int set_cwt(TDAContext *ctx);        
int wr_sys(TDAContext *ctx);
int rd_sys(TDAContext *ctx);          
int dblock(TDAContext *ctx);          
int dblock_alloc(TDAContext *ctx, int n);   
int repsel(TDAContext *ctx);          
int repsel_off(TDAContext *ctx);   
int repsel_alloc(TDAContext *ctx, int n);   

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */


/* ------------------------------------------------------------------------ */
/*  clear()         clear variables: clear or clear=varlist.                */
/*                  return 0 if OK, -1 if error.                            */
 
int clear(TDAContext *ctx)
{
    register int i,j,k,l; 
    int err,n,cwtflg,edflg,gdflg,seqflg,arcflg,sdflg,nb;
    register char *p;

    err = -1;
    sdflg = seqflg = edflg = gdflg = cwtflg = arcflg = 0;

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        printf1(ctx, "No data matrix (command ignored).\n");
        return(0);
    }
    if (!strcmp(ctx->CmdBuf,"clear")) {
        repsel_off(ctx);
        clear_dm(ctx);         /* clear data matrix, variables, namelists */
        arcflg = cwtflg = edflg = gdflg = seqflg = 1;
        err = 0;
    }
    else {
        p = ctx->CmdBuf + 5;
        if (*p++ != '=') {
            p_err(ctx, -1,1);
            return(-1);
        }
        p = get_nvia(ctx, p,&n,1,&nb);
        if (nb) {
            p_err(ctx, -42,1);
            goto CLEARFin;
        }
        if (n < 1 || *p) {
            if (*p)
                p_err(ctx, -4,1);
            goto CLEARFin;
        }
        for (i = 0; i < ctx->VLNV; ++i) {
            j = ctx->VLVIdx[i];
            if (j == ctx->WIVar)    
                cwtflg = 1;
            if (edflg == 0 && check_edj(ctx, j))
                edflg = 1;
            if (seqflg == 0 && check_sdj(ctx, j))
                seqflg = 1;
            if (gdflg == 0 && check_gdj(ctx, j))
                gdflg = 1;

            if (sdflg == 7 && ctx->SDVarDef && j < 4)
                sdflg = 1;  

            if (ctx->NNL > 0) {
                for (k = 0; k < MaxNL; ++k) {
                    n = ctx->NLNV[k];
                    for (l = 0; l < n; ++l) {
                        if (j == ctx->NLVIdx[k][l]) {
                            printf1(ctx, "Removed namelist: %s\n",ctx->NLName[k]);
                            nl_free(ctx, k);
                            break;
                        }
                    }
                }
            }
            clear_var(ctx, j);
        }
        err = 0;
    }
    if (ctx->VIFirst == -1) { 
        printf1(ctx, "Deleted whole data matrix. ");
        tsel_off(ctx, 0);       
        sdnvar_close(ctx);
        ctx->NOC = ctx->NOCDM = ctx->NOCMaxA = ctx->DMDef = 0;
        arcflg = 1;
    }
    if (cwtflg && ctx->WIVar >= 0) {       /* turn off case weights */
        ctx->WIVar = -1;
        ctx->WSum = ctx->WSumS = 0.0;
        ctx->WNormFlag = 0;
        printf1(ctx, "Previously defined case weights turned off.\n");
    }
    if (gdflg && ctx->GD_PERM == 0)    /* free gdd data */
        gdd_free(ctx, 1);   

    if (edflg)          /* free episode data */
        edat_off(ctx, 1);    
         
    if (seqflg)  
        seq_afree(ctx, 1);   /* free sequence data */

    if (arcflg)
        arcd_off(ctx);           /* turn off archive */

    if (sdflg)
        sdnvar_close(ctx);

    prn_mem(ctx);

CLEARFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  clear_a()       clear variables: clear or clear=varlist.                */
/*                  return 0 if OK, -1 if error.                            */
 
void clear_a(TDAContext *ctx)
{
    repsel_off(ctx);
    clear_dm(ctx);
 
    tsel_off(ctx, 0);       
    ctx->NOC = ctx->NOCDM = ctx->NOCMaxA = ctx->DMDef = 0;

    if (ctx->WIVar >= 0) {       /* turn off case weights */
        ctx->WIVar = -1;
        ctx->WSum = ctx->WSumS = 0.0;
        ctx->WNormFlag = 0;
    }
    if (ctx->GD_PERM == 0)     /* free gdd data */
        gdd_free(ctx, 0);   

    edat_off(ctx, 0);          /* free episode data */
    seq_afree(ctx, 0);         /* free sequence data */
    arcd_off(ctx);           /* turn off archive */
}

/* ------------------------------------------------------------------------ */
/*  clearnl()       clear namelists                                         */
/*                  return 0 if OK, -1 if error.                            */
 
int clearnl(TDAContext *ctx)
{
    register int i;
    int fin,fnd;
    register char *p,*q;

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->NNL == 0) {
        printf1(ctx, "No namelists defined.\n");
        return(0);
    }
    p = ctx->CmdBuf + 7;
    if (*p++ != '=') {
        p_err(ctx, -1,1);
        return(-1);
    }
    fin = 0;
    while (*p) {
        q = p;
        while (*q && *q != ',')  
            q++;
        if (!*q)
            fin = 1;
        else if (*q != ',') {
            p_err(ctx, -1,1);
            return(-1);
        }
        *q++ = '\0';
        fnd = 0;
        for (i = 0; i < MaxNL; ++i) {
            if (!strcmp(p,ctx->NLName[i])) {
                printf1(ctx, "Removed namelist: %s\n",ctx->NLName[i]);
                nl_free(ctx, i);
                fnd = 1;
                break;
            }
        }
        if (fnd == 0)  
            printf1(ctx, "Cannot find: %s\n",p);

        if (fin)
            break;
        p = q;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  tsel()      tsel = expression.                                          */
/*              Create TSelect and set NOC accordingly.                     */
/*              Update case weights if defined.                             */
/*                                                                          */
/*              Standard tsel expression must not contain type 2 or 3       */
/*              operators. However, the command can now also be used        */
/*              with matrix expressions.                                    */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int tsel(TDAContext *ctx)        
{
    register int i;
    int err,m,n,opt,row,col,ivflg;
    double ws,tmp;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        printf1(ctx, "No data matrix (command ignored).\n");
        return(0);
    }
    if (!strcmp(ctx->CmdBuf,"tsel=off")) {
        if (ctx->TSelFlg == 0)  
            printf1(ctx, "No active case selection (command ignored).\n");
        else {
            tsel_off(ctx, 2);
            set_cwt(ctx);      /* adjust weights */
        }
        return(0);
    }
    if (ctx->TSelFlg)
        tsel_off(ctx, 1);        /* free previously active tsel command */
    else {
        if (ctx->EDAvail)  
            edat_off(ctx, 1);    /* turn off episode data */
        if (ctx->SeqDN)  
            seq_afree(ctx, 1);   /* turn off sequence data */
        if (ctx->GD_PERM == 0)
            gdd_free(ctx, 1);    /* free relational data */
        repsel_off(ctx);       /* turn off repsel */
    }

    opt = 0;        /* check for standard expression */
    ivflg = 0;

    if ((n = v_parse(ctx, ctx->CmdBuf + 5,0)) < 0 || ctx->ESCnt <= 0) {

        /* if not a standard expression check for matrix expression */
           
        if (get_mexpr(ctx, ctx->CmdBuf + 5,&row,&col,&ivflg)) {
            printf1(ctx, "Error: can't evaluate the tsel expression.\n");
            goto TSELFin;
        }
        if (row != ctx->NOCDM || col != 1) {
            printf1(ctx, "Error: matrix expression must result in NOC x 1 matrix.\n");
            printf1(ctx, "Found a %d x %d matrix.\n",row,col);
            goto TSELFin;
        }
        if (alloc_actmp(ctx, ctx->NOCDM + 1))
            goto TSELFin;   

        if (eval_mexpr(ctx, 0,row,col,ctx->AcTmp,0,ctx->AcTmp,NULL,NULL,-1,0.0)) {
            printf1(ctx, "Error: can't evaluate the tsel expression.\n");
            goto TSELFin;
        }   
        opt = 1;
    }

    if (opt == 0) {

        /* check for type 2 operators */

        for (i = 0; i < ctx->ESCnt; ++i) {
            n = iabs(ctx, ctx->ESTyp[i]);
            if (n >= ctx->OPT2A && n < ctx->OPT2B) {
                printf1(ctx, "Error: tsel expressions may not contain type 2 operators.\n");
                goto TSELFin; 
            }
        }
    }
    printf1(ctx, "New temporary case selection: %s\n",ctx->CmdBuf + 5);

    if (!(ctx->TSelect = (int *)calloc((size_t)(ctx->NOCDM),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto TSELFin;
    }
    ctx->TSelectA = ctx->NOCDM;
    memrq(ctx, ctx->TSelectA,sizeof(int));
           
    /* create the select indicator variables */

    m = 0;
    ws = 0.0;
    for (i = 0; i < ctx->NOCDM; ++i) {

        if (opt == 0) {
            n = v_eval1(ctx, i,ctx->ESCnt,ctx->ESTyp,ctx->ESVal,ctx->ESIdx,&tmp,0,0,0,0,0);
            if (n) {
                printf1(ctx, "Error: can't evaluate tsel expression in case %d.\n",i + 1);
                prn_emsg2(ctx, n);
                goto TSELFin;
            }
        }
        else  
            tmp = ctx->AcTmp[i + 1];

        if (fabs(tmp) > ctx->EPSI2) {
            ctx->TSelect[m] = i;
            if (ctx->WIVar >= 0)  
                ws += get_data(ctx, ctx->WIVar,i);
            m++;
        }
    }
    if (m == 0) {
        printf1(ctx, "Error: number of selected cases is zero.\n");
        goto TSELFin;
    }
    if (ctx->WIVar >= 0 && ws <= ctx->EPSI2) {
        printf1(ctx, "Error: sum of weights for selected cases is almost zero.\n");
        goto TSELFin;
    }
    printf1(ctx, "Number of selected cases: %d\n",m);
    ctx->NOC = m;
    ctx->TSelFlg = 1;

    if (ctx->WIVar >= 0) {
        printf1(ctx, "Corresponding sum of weights: %lg",ws);
        if (ctx->WNormFlag)
            ctx->WSumS = (double)ctx->NOC;

        if (ctx->WSumS > 0.0) {
            ctx->WSum = ctx->WSumS;
            ctx->WNorm = ctx->WSum / ws;
            printf1(ctx, " (adjusted to: %lg).\n",ctx->WSum);
        }
        else {
            ctx->WSum = ws;
            ctx->WNorm = 1.0;
            printf1(ctx, "\n");
        }
    }
    else {
        ctx->WSum = (double)ctx->NOC;
        ctx->WNorm = 1.0;
        ctx->WSumS = 0.0;
    }
    err = 0;  

TSELFin:
    alloc_actmp(ctx, 0);
    if (err) {
        if (ctx->TSelectA > 0) {
            free((char *)ctx->TSelect);
            memrq(ctx, -ctx->TSelectA,sizeof(int));
            ctx->TSelectA = 0;
        }
        set_cwt(ctx);
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  tsel_off(opt)   Turn off tsel, free memory. Reset NOC to NOCDM.         */
/*                  Also free the following data structures:                */
/*                  - episode data                                          */
/*                  - sequence data                                         */
/*                  - relational data                                       */

void tsel_off(TDAContext *ctx, int opt)        
{
    if (ctx->TSelFlg == 0)
        return;

    if (ctx->EDAvail)  
        edat_off(ctx, opt);      /* turn off episode data */
         
    if (ctx->SeqDN)  
        seq_afree(ctx, opt);     /* turn off sequence data */
         
    if (ctx->GD_PERM == 0)
        gdd_free(ctx, opt);      /* free relational data */

    if (ctx->TSelectA > 0) {
        free((char *)ctx->TSelect);
        memrq(ctx, -ctx->TSelectA,sizeof(int));
        ctx->TSelectA = 0;
    }
    ctx->TSelFlg = 0;
    ctx->NOC = ctx->NOCDM;
    if (opt) {
        printf1(ctx, "Temporary case selection turned off.\n");
        if (opt == 2)  
            printf1(ctx, "Number of cases: %d\n",ctx->NOC);
    }
}

/* ------------------------------------------------------------------------ */
/*  cwt()       cwt(wnorm=) = VName                                         */
/*              create case weights.                                        */
/*              Return 0 if OK, -1 if error.                                */

int cwt(TDAContext *ctx)        
{
    int iv;
    double wn;
    char *p,vname[VNLMax+1];

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        printf1(ctx, "No data matrix (command ignored).\n");
        return(0);
    }
    ctx->WSumS = 0.0;
    if (ctx->WIVar >= 0) {
        printf1(ctx, "Previously defined case weights turned off.\n");
        ctx->WIVar = -1;
        ctx->WSum = (double)ctx->NOC;
        ctx->WSumS = 0.0;
        ctx->WNorm = 1.0;
        ctx->WNormFlag = 0;
        iv = 0;
    }
    else
        iv = 1;

    if (!strcmp(ctx->CmdBuf,"cwt=off")) {
        if (iv)
            printf1(ctx, "No case weights defined (command ignored).\n");
        return(0);
    }  
    ctx->WNormFlag = 0;

    p = ctx->CmdBuf + 3;
    if (*p == '(') {
        if (!strncmp(p,"(wnorm)",7)) {
            p += 7;
            ctx->WNormFlag = 1;
        }
        else {
            if (sscanf(p + 1,"wnorm=%lf",&wn) != 1 || wn <= 0.0) {
                p_err(ctx, -1,1);
                return(-1);
            }
            p = skip_blev(ctx, p);
            ctx->WSumS = wn;
        }
    }
    if (*p++ != '=') {
        p_err(ctx, -1,1);
        return(-1);
    }
    if ((iv = get_vidx(ctx, p,vname)) < 0) {
        printf1(ctx, "Error: undefined variable.\n");
        return(-1);
    }
    ctx->WIVar = iv;
    set_cwt(ctx);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  set_cwt()       Set variables for case weights.                         */
/*                  a)  WIVar < 0 (no case weights), then                   */
/*                      WSum = NOC, WNorm = 1.                              */
/*                  b)  WIVar >= 0 and WSumS == 0                           */
/*                      WSum = sum of W(i) and WNorm = 1.                   */
/*                  c)  WIVar >= 0 and WSumS > 0.0 or WNormFlag == 1        */
/*                      WSum = WSumS = sum of W(i) * WNorm,                 */
/*                      WNorm = WSumS / sum of W(i).                        */
/*                                                                          */
/*  Return 0 if OK, -1 if error (command is then ignored).                  */

int set_cwt(TDAContext *ctx)         
{
    register int i;
    double tmp,ws;

    ctx->WSum = (double)ctx->NOC;
    ctx->WNorm = 1.0;

    if (ctx->WIVar < 0)  
        return(0);
       
    ws = 0.0;
    for (i = 0; i < ctx->NOC; ++i) {
        tmp = get_data(ctx, ctx->WIVar,i);
        if (tmp < 0.0) {
            printf1(ctx, "Error: found negative weight in case %d (command ignored).\n",i + 1);
            ctx->WIVar = -1;
            return(-1);
        }
        ws += tmp;
    }
    printf1(ctx, "New case weights. Sum of weights: %lg",ws);
    if (ws <= ctx->EPSI) {
        printf1(ctx, " (weights will not be used).\n");
        ctx->WIVar = -1;
        return(-1);
    }
    if (ctx->WNormFlag)
        ctx->WSumS = (double)ctx->NOC;

    if (ctx->WSumS > 0.0) {
        ctx->WSum = ctx->WSumS;
        ctx->WNorm = ctx->WSum / ws;
        printf1(ctx, ", adjusted to: %lg\n",ctx->WSum);
    }
    else {
        ctx->WSum = ws;
        ctx->WNorm = 1.0;            
        printf1(ctx, "\n");
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  wr_sys()        Execute wsys command in CmdBuf.                         */
/*                  Syntax: wsys or wsys = fname                            */
/*                  Return 0 if OK, -1 if err.                              */

int wr_sys(TDAContext *ctx)
{
    FILE *fd;
    register int j;
    int err,r,len,nv;
    char *p,*fname;
    char *sysname = "tda.sys";

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        printf1(ctx, "No data matrix (command ignored).\n");
        return(0);
    }
    p = ctx->CmdBuf + 4;

    if (!*p)
        fname = sysname;
    else if (*p++ == '=' && *p)
        fname = p;
    else {
        printf1(ctx, "Syntax error.\n");
        return(-1);    
    }
    nv = 0;
    j = ctx->VIFirst;
    while (j >= 0) {
        if (ctx->VTyp[j] != 5)      /* skip type 5 variables */
            nv++;
        j = ctx->VNxt[j];
    }
    if (nv == 0) {
        printf1(ctx, "No variables (command ignored).\n");
        return(0);
    }
    if (!(fd = fopen(fname,OPEN_WB))) {
        printf1(ctx, "Can't create: %s\n",fname);
        return(-1);   
    }
    printf1(ctx, "Writing system file: %s\n",fname);

    fprintf(fd,"TDA System File (%3.1f). ",TDA_Version);
#if TIME_ON
    prn_time(ctx, fd);
#else
    fprintf(fd,"\n");
#endif

    fprintf(fd,"%d %d\n",ctx->NOCDM,nv);

    j = ctx->VIFirst;
    while (j >= 0) {
        if (ctx->VTyp[j] != 5) {  
            fprintf(fd,"%s %s ",ctx->VName[j],ctx->VDef[j]);
            if (ctx->VLabel[j] != NULL)
                fprintf(fd,"%s ",ctx->VLabel[j]);
            fprintf(fd,"\n");
            j = ctx->VNxt[j];
        }
    }
    j = ctx->VIFirst;
    while (j >= 0) {
        if (ctx->VTyp[j] != 5)    
            fprintf(fd,"%d %d %d %d\n",ctx->VTyp[j],ctx->VSLen[j],ctx->VPFmt1[j],ctx->VPFmt2[j]);
        j = ctx->VNxt[j];
    }

    j = ctx->VIFirst;
    while (j >= 0) {
        if (ctx->VTyp[j] != 5) {  
            len = get_slen(ctx, j,ctx->NOCDM);
            r = (int)(fwrite(ctx->VDPtr[j],(size_t)(len),1,fd));
            if (r != 1) {
                printf1(ctx, "Error in writing system file %s\n",fname);
                goto WSYSFin;
            }
            j = ctx->VNxt[j];
        }
    }
    printf1(ctx, "Data matrix (%d cases, %d variables) written to: %s\n",ctx->NOCDM,nv,fname);
    err = 0;   

WSYSFin:
    fclose(fd);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rd_sys()    Read TDA system file. CmdBuf: rsys or rsys=fname.           */
/*              Return 0 if OK, -1 if error.                                */

int rd_sys(TDAContext *ctx)           
{
    FILE *fd;
    register int j;
    int nn,l,l1,l2,l3,nv,m,n,w1,w2,err;
    char *p,*q,*fname,*sysname = "tda.sys";
    char vname[VNLMax + 1];
    char vdef[RLMaxDef];
    char vlabel[VLLMax + 1];

    nn = 20000;                     /* length of read buffer */

    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef) {
        printf1(ctx, "Error: a data matrix already exists (command ignored).\n");
        return(0);
    }
    p = ctx->CmdBuf + 4;

    if (!*p)
        fname = sysname;
    else if (*p++ == '=' && *p)
        fname = p;
    else {
        printf1(ctx, "Syntax error.\n");
        return(-1);    
    }
    if (!(fd = fopen(fname,OPEN_RB))) {
        printf1(ctx, "Can't open: %s\n",fname);
        return(-1);   
    }
    printf1(ctx, "Reading system file: %s\n",fname);

    ctx->VLabelLen = 0;
    if (alloc_acc(ctx, nn + 1)) {
        err = -1;
        goto RSYSFin; 
    }
    err = 0;
    if (!fgets(ctx->AcC,nn,fd)) {
        err = 1;
        goto RSYSFin;
    }  
    if (strncmp(ctx->AcC,"TDA System File",15)) {
        err = -2;
        goto RSYSFin;
    }
    printf1(ctx, "Identification: %s",ctx->AcC);

    if (!fgets(ctx->AcC,nn,fd)) {
        err = 2;
        goto RSYSFin;
    }  
    if (sscanf(ctx->AcC,"%d %d",&n,&nv) != 2 || n <= 0 || nv <= 0) {
        err = 3;
        goto RSYSFin;
    }       
    if (nv > ctx->MaxNV) {
        err = -3;
        goto RSYSFin;
    }
    ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = n;
    ctx->NVAR = 0;

    /* get variable names and definitions */
            
    for (j = 0; j < nv; ++j) {

        if (!fgets(ctx->AcC,nn,fd)) {
            err = 4;
            goto RSYSFin;
        }
        p = ctx->AcC;
        q = vname;
        l1 = 0;
        for (l = 0; l < VNLMax; ++l) {
            if (!*p || *p == ' ')
                break;
            *q++ = *p++;
            l1++;
        }
        *q = '\0';

        if (*p++ != ' ' || !*p) {
            err = 5;
            goto RSYSFin;
        }
        q = vdef;
        l2 = 0;
        for (l = 0; l < nn; ++l) {
            if (!*p || *p == ' ')
                break;
            *q++ = *p++;
            l2++;
        }
        *q = '\0';

        if (*p++ != ' ') {
            err = 6;
            goto RSYSFin;
        }
        l3 = 0;
        if (*p && *p != LF && *p != CR) {          /* label */

            q = vlabel;
            for (l = 0; l < VLLMax; ++l) {
                if (!*p || *p == ' ')
                    break;
                *q++ = *p++;
                l3++;
            }
            *q = '\0';
        }
        l1++;
        if (!(ctx->VName[j] = (char *)calloc((size_t)(l1),sizeof(char)))) { 
            p_err(ctx, -2,1);
            err = -1;
            goto RSYSFin;
        }
        l2++;
        if (!(ctx->VDef[j] = (char *)calloc((size_t)(l2),sizeof(char)))) { 
            p_err(ctx, -2,1);
            free(ctx->VName[j]);
            err = -1;
            goto RSYSFin;
        }
        memrq(ctx, l1,1);
        strcpy(ctx->VName[j],vname);
        memrq(ctx, l2,1);
        strcpy(ctx->VDef[j],vdef);
     
        if (l3 > 0) {
            if (ctx->VLabelLen < l3)
                ctx->VLabelLen = l3;
            l3++;
            if (!(ctx->VLabel[j] = (char *)calloc((size_t)(l3),sizeof(char)))) { 
                p_err(ctx, -2,1);
                free(ctx->VName[j]);
                memrq(ctx, -l1,1);
                free(ctx->VDef[j]);
                memrq(ctx, -l2,1);
                ctx->VLabel[j] = NULL;
                err = -1;
                goto RSYSFin;
            }
            memrq(ctx, l3,1);
            strcpy(ctx->VLabel[j],vlabel);
        }
        ctx->VAlloc[j] = 1;
        if (ctx->VIFirst < 0)
            ctx->VILast = ctx->VIFirst = (short)j;
        else {
            ctx->VNxt[ctx->VILast] = (short)(j); 
            ctx->VILast = (short)(j);
        }
        ctx->NVAR++;
    }
    for (j = 0; j < ctx->NVAR; ++j) {

        if (!fgets(ctx->AcC,RLMaxDef,fd)) {
            err = 7;
            goto RSYSFin;
        }
        if (sscanf(ctx->AcC,"%d %d %d %d",&n,&m,&w1,&w2) != 4) {
            err = 8;
            goto RSYSFin;
        }
        if (m >= 0 && m != 0 && m != 1 && m != 2 && m != 4 && m != 5 && m != 8) {
            err = 9;
            goto RSYSFin;
        }
        ctx->VTyp[j] = (char)n;
        ctx->VTypA[j] = 1;
        ctx->VSLen[j] = (short)m;
        makefmt(ctx, &w1,&w2,ctx->VPFmtS[j],VPFmtSLen,0,' ',0);
        ctx->VPFmt1[j] = (short)w1;
        ctx->VPFmt2[j] = (short)w2;
    }
    newline(ctx);           
    prn_var(ctx, ctx->VIFirst);
    newline(ctx);

    /* allocate memory for data */

    if (alloc_vdat(ctx, 0,1)) {
        p_err(ctx, -2,1);
        err = -1;
        goto RSYSFin;
    }
        
    /* read data */

    for (j = 0; j < ctx->NVAR; ++j) {

        l = get_slen(ctx, j,ctx->NOCDM);
        if (((n = (int)fread(ctx->VDPtr[j],(size_t)(l),1,fd))) != 1) {
            err = 10;
            goto RSYSFin;
        }
    }
    ctx->DMDef = 1;
    printf1(ctx, "Created a data matrix with %d variables and %d cases.\n",ctx->NVAR,ctx->NOCDM);

RSYSFin:
    fclose(fd);
    if (err > 0) {
        printf1(ctx, "Error (%d) in reading system file.\n",err);
        if (err == 10)  
            printf1(ctx, "R=%d NOC=%d L=%d J=%d\n",n,ctx->NOCDM,l,j);
    }
    else if (err == -2)
        printf1(ctx, "Probably not a TDA system file.\n");
    else if (err == -3) {
        printf1(ctx, "Error: system file contains %d variables.\n",nv);    
        printf1(ctx, "More than the current maximum: %d.\n",ctx->MaxNV);
    }
    if (err) {
        clear_dm(ctx);            /* clear all variables */
        ctx->NOCMaxA = ctx->NOCDM = ctx->NOC = 0;
        ctx->NVAR = 0;
        err = -1;
    }
    alloc_acc(ctx, 0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dblock      dblock(mdef=MName) = varlist                                */
/*  ##          Creates a new block data structure according to varlist     */
/*              Command is in CmdBuf                                        */
/*              Creates DBlckPtr[i] for i = 0,...,NOC - 1.                  */
/*              where DBlckPtr[i] is the block number of case i.            */
/*              If dblock; without arguments, then each case is a           */
/*              separate block.                                             */  
/*              If mdef parameter, create vector with block info.           */
/*              Note: tsel = off is currently active.                       */
/*              Return 0 if OK, -1 if error.                                */

int dblock(TDAContext *ctx)           
{
    register int i,j,k,l;
    int err,idx;
    double tmp;

    idx = err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 6,4,0))  
        goto DBLOCKFin;
       
    tsel_off(ctx, 0);                            /* turn off tsel */
    repsel_off(ctx);                           /* turn off repsel */

    if (ctx->NOC != ctx->DBlckPtrA) {
        if (dblock_alloc(ctx, ctx->NOC))  
            goto DBLOCKFin;
    }
    if (ctx->PMNV > 0) {
        if (alloc_actmp(ctx, ctx->PMNV))
            goto DBLOCKFin;
        for (j = 0; j < ctx->PMNV; ++j)
            ctx->AcTmp[j] = get_data(ctx, ctx->PMVIdx[j],0);
        k = 0;
    }
    else
        k = -1;

    if (ctx->PMatNameFlg) {
        if ((idx = mat_newmat(ctx, ctx->PMatName,ctx->NOC,1)) >= 0)  
            mdefcpy(ctx, ctx->MatDef[idx],ctx->CmdBuf);
    }
    for (i = 0; i < ctx->NOC; ++i) {
        if (ctx->PMNV == 0) {                    /* block = case */
            k++;
        }
        else {                              /* create blocks */
            l = 0;
            for (j = 0; j < ctx->PMNV; ++j) {
                tmp = get_data(ctx, ctx->PMVIdx[j],i);
                if (tmp != ctx->AcTmp[j]) {
                    l = 1;
                    ctx->AcTmp[j] = tmp;
                }
            }
            if (l)
                k++;
        }
        ctx->DBlckPtr[i] = k;
        if (idx >= 0)
            ctx->MatVal[idx][i + 1] = (double)(k + 1);            
    }
    ctx->BNOC = k + 1;
    printf2(ctx, "Found %d block(s).\n",ctx->BNOC);
    err = 0;

DBLOCKFin:
    p_clean(ctx);
    if (err)
        ctx->BNOC = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dblock_alloc(n)     if n > 0 allocate DBlckPtr, else free.              */
/*                      Return 0 if OK, -1 if error.                        */

int dblock_alloc(TDAContext *ctx, int n)    
{
    if (n == 0 || n != ctx->DBlckPtrA) {
        if (ctx->DBlckPtrA > 0) {
            free((char *)ctx->DBlckPtr);
            memrq(ctx, -ctx->DBlckPtrA,sizeof(int));
            ctx->BNOC = ctx->DBlckPtrA = 0;
        }
    }
    if (n <= 0)
        return(0);

    if (!(ctx->DBlckPtr = (int *)calloc((size_t)(n),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);           
    }
    ctx->DBlckPtrA = n;
    memrq(ctx, n,sizeof(int));
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  repsel      repsel (                                                    */
/*  ##              id = V,                                                 */
/*              ) = m-expression.                                           */
/*                                                                          */
/*              m-expression must have dimension BNOC x 1. Call this        */
/*              expression S. There are two different forms of the command. */
/*                                                                          */
/*              a)  without id, S(i) is the number of times that block i    */
/*                  is repeated in data generation by get_data().           */
/*                                                                          */
/*              b)  if id=V is given, then those cases are selected from    */  
/*                  from block i where S(i) equals the value of V.          */
/*                                                                          */
/*              Note: if tsel is currently active, then tsel = off.         */
/*              Return 0 if OK, -1 if error.                                */

int repsel(TDAContext *ctx)           
{
    register int i,j,k,l,kk,k1;
    int err,row,col,n,m,mm = 0,nn,ivflg;
    char *p;

    err = -1;
    ivflg = 0;

    if (check_cmd(ctx, 1))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 6,8,1))  
        goto REPSELFin;

    if (ctx->PMRHSTRA < 1)
        goto REPSELFin;

    repsel_off(ctx);

    if (!strcmp(ctx->PMRHSTR,"off")) {
        err = 0;
        goto REPSELFin;
    }
    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        goto REPSELFin;
    }
    if (ctx->BNOC == 0) {
        printf1(ctx, "Error: need a previous dblock command.\n");
        goto REPSELFin;
    }
    p = ctx->PMRHSTR;                                    /* get expression */
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto REPSELFin;
    if (col != 1) {
        mat_err(ctx, 7);
        goto REPSELFin;
    }
    if (*p) {
        mat_err(ctx, 0);
        goto REPSELFin;
    }
    if (ctx->BNOC != row) {
        printf1(ctx, "Error: number of blocks is %d, found %d rows.\n",ctx->BNOC,row);
        goto REPSELFin;
    }
    tsel_off(ctx, 0);                /* turn off tsel */

    nn = 0;                     /* calc new number of cases */
    m = ctx->DBlckPtr[0];   
    kk = k = 0;
    for (i = 1; i <= ctx->BNOC; ++i) {
        n = (int)ctx->MX[0][i];
        while (k <= ctx->NOC) {
            if (k < ctx->NOC)
                mm = ctx->DBlckPtr[k];
            else 
                mm++;

            if (mm == m) {
                if (ctx->PMID >= 0) {        /* id variable defined */
                    if (fabs(get_data(ctx, ctx->PMID,k) - ctx->MX[0][i]) <= ctx->EPSI1)
                        nn++;
                }
            }
            else {
                if (ctx->PMID < 0) {
                    for (j = 1; j <= n; ++j) {
                        k1 = kk;
                        while (k1++ < k)
                            nn++;
                    }
                    kk = k;
                }
                m = mm;
                break;
            }
            k++;
        }
    }
    if (nn == 0) {
        printf1(ctx, "Error: selected number of cases is zero.\n");
        goto REPSELFin;
    }
    if (repsel_alloc(ctx, nn))  
        goto REPSELFin;

    m = ctx->DBlckPtr[0];
    kk = l = k = 0;
    for (i = 1; i <= ctx->BNOC; ++i) {
        n = (int)ctx->MX[0][i];
        while (k <= ctx->NOC) {
            if (k < ctx->NOC)
                mm = ctx->DBlckPtr[k];
            else
                mm++;

            if (mm == m) {
                if (ctx->PMID >= 0) {
                    if (fabs(get_data(ctx, ctx->PMID,k) - ctx->MX[0][i]) <= ctx->EPSI1)
                        ctx->REPSelect[l++] = k;
                }
            }
            else {
                if (ctx->PMID < 0) {
                    for (j = 1; j <= n; ++j) {
                        k1 = kk;
                        while (k1 < k) {
                            ctx->REPSelect[l++] = k1;
                            k1++;
                        }
                    }
                    kk = k;
                }
                m = mm;
                break;
            }
            k++;
        }
    }
    /*** 
    printf1(ctx, "REPSelect:\n");
    for (l = 0; l < nn; ++l)
        printf1(ctx, "l=%2d sel=%d\n",l,REPSelect[l]);
    if (PMID >= 0)
        printf1(ctx, "ID variable: %s\n",VName[PMID]);
    ****/

    ctx->REPNOCOld = ctx->NOC;            /* save old NOC */
    ctx->NOC = nn;                   /* set new NOC */
    ctx->REPSelFlg = 1;              /* make active */
    printf2(ctx, "Selected %d cases.\n",ctx->NOC);

    set_cwt(ctx);                  /* adjust weights if defined */
    err = 0;

REPSELFin:
    mx_free(ctx);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  repsel_off()        Turns off the repsel command.                       */
/*                      Return 0 if OK, -1 if error.                        */

int repsel_off(TDAContext *ctx)    
{
    if (ctx->REPSelFlg) {
        if (ctx->DMDef) {
            if (ctx->REPNOCOld != ctx->NOCDM)  
                gerr_exit(ctx, 100);
            ctx->NOC = ctx->REPNOCOld;
        }
        printf2(ctx, "Repsel turned off.\n");
    }
    ctx->REPNOCOld = 0;
    repsel_alloc(ctx, 0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  repsel_alloc(n)     if n > 0 allocate REPSelect, else free.             */
/*                      Return 0 if OK, -1 if error.                        */

int repsel_alloc(TDAContext *ctx, int n)    
{
    ctx->REPSelFlg = 0;
    if (n == 0 || n != ctx->REPSelectA) {
        if (ctx->REPSelectA > 0) {
            free((char *)ctx->REPSelect);
            memrq(ctx, -ctx->REPSelectA,sizeof(int));
            ctx->REPSelectA = 0;
        }
        if (n == 0)
            return(0);
    }
    else
        return(0);

    if (!(ctx->REPSelect = (int *)calloc((size_t)(n),sizeof(int)))) {
        p_err(ctx, -2,1);
        return(-1);           
    }
    ctx->REPSelectA = n;
    memrq(ctx, n,sizeof(int));
    return(0);
}


