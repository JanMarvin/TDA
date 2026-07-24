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

/*  functions in t_mdat.c */

int clear(void);
void clear_a(void);
int clearnl(void);
int tsel(void);       
void tsel_off(int opt);       
int cwt(void);       
int set_cwt(void);        
int wr_sys(void);
int rd_sys(void);          
int dblock(void);          
int dblock_alloc(int n);   
int repsel(void);          
int repsel_off(void);   
int repsel_alloc(int n);   

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

int *DBlckPtr;              /* maps cases to blocks                         */
int DBlckPtrA = 0;          /* number of allocated ints                     */
int BNOC = 0;               /* number of blocks                             */ 
int *REPSelect;             /* indices for repsel                           */
int REPSelectA = 0;         /* number of allocated ints                     */
int REPSelFlg = 0;          /* set if repsel active                         */
int REPNOCOld = 0;          /* saves old NOC                                */

/* ------------------------------------------------------------------------ */
/*  clear()         clear variables: clear or clear=varlist.                */
/*                  return 0 if OK, -1 if error.                            */
 
int clear(void)
{
    register int i,j,k,l; 
    int err,n,cwtflg,edflg,gdflg,seqflg,arcflg,sdflg,nb;
    register char *p;

    err = -1;
    sdflg = seqflg = edflg = gdflg = cwtflg = arcflg = 0;

    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        printf1("No data matrix (command ignored).\n");
        return(0);
    }
    if (!strcmp(CmdBuf,"clear")) {
        repsel_off();
        clear_dm();         /* clear data matrix, variables, namelists */
        arcflg = cwtflg = edflg = gdflg = seqflg = 1;
        err = 0;
    }
    else {
        p = CmdBuf + 5;
        if (*p++ != '=') {
            p_err(-1,1);
            return(-1);
        }
        p = get_nvia(p,&n,1,&nb);
        if (nb) {
            p_err(-42,1);
            goto CLEARFin;
        }
        if (n < 1 || *p) {
            if (*p)
                p_err(-4,1);
            goto CLEARFin;
        }
        for (i = 0; i < VLNV; ++i) {
            j = VLVIdx[i];
            if (j == WIVar)    
                cwtflg = 1;
            if (edflg == 0 && check_edj(j))
                edflg = 1;
            if (seqflg == 0 && check_sdj(j))
                seqflg = 1;
            if (gdflg == 0 && check_gdj(j))
                gdflg = 1;

            if (sdflg == 7 && SDVarDef && j < 4)
                sdflg = 1;  

            if (NNL > 0) {
                for (k = 0; k < MaxNL; ++k) {
                    n = NLNV[k];
                    for (l = 0; l < n; ++l) {
                        if (j == NLVIdx[k][l]) {
                            printf1("Removed namelist: %s\n",NLName[k]);
                            nl_free(k);
                            break;
                        }
                    }
                }
            }
            clear_var(j);
        }
        err = 0;
    }
    if (VIFirst == -1) { 
        printf1("Deleted whole data matrix. ");
        tsel_off(0);       
        sdnvar_close();
        NOC = NOCDM = NOCMaxA = DMDef = 0;
        arcflg = 1;
    }
    if (cwtflg && WIVar >= 0) {       /* turn off case weights */
        WIVar = -1;
        WSum = WSumS = 0.0;
        WNormFlag = 0;
        printf1("Previously defined case weights turned off.\n");
    }
    if (gdflg && GD_PERM == 0)    /* free gdd data */
        gdd_free(1);   

    if (edflg)          /* free episode data */
        edat_off(1);    
         
    if (seqflg)  
        seq_afree(1);   /* free sequence data */

    if (arcflg)
        arcd_off();           /* turn off archive */

    if (sdflg)
        sdnvar_close();

    prn_mem();

CLEARFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  clear_a()       clear variables: clear or clear=varlist.                */
/*                  return 0 if OK, -1 if error.                            */
 
void clear_a(void)
{
    repsel_off();
    clear_dm();
 
    tsel_off(0);       
    NOC = NOCDM = NOCMaxA = DMDef = 0;

    if (WIVar >= 0) {       /* turn off case weights */
        WIVar = -1;
        WSum = WSumS = 0.0;
        WNormFlag = 0;
    }
    if (GD_PERM == 0)     /* free gdd data */
        gdd_free(0);   

    edat_off(0);          /* free episode data */
    seq_afree(0);         /* free sequence data */
    arcd_off();           /* turn off archive */
}

/* ------------------------------------------------------------------------ */
/*  clearnl()       clear namelists                                         */
/*                  return 0 if OK, -1 if error.                            */
 
int clearnl(void)
{
    register int i;
    int fin,fnd;
    register char *p,*q;

    if (check_cmd(1))
        return(-1);

    if (NNL == 0) {
        printf1("No namelists defined.\n");
        return(0);
    }
    p = CmdBuf + 7;
    if (*p++ != '=') {
        p_err(-1,1);
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
            p_err(-1,1);
            return(-1);
        }
        *q++ = '\0';
        fnd = 0;
        for (i = 0; i < MaxNL; ++i) {
            if (!strcmp(p,NLName[i])) {
                printf1("Removed namelist: %s\n",NLName[i]);
                nl_free(i);
                fnd = 1;
                break;
            }
        }
        if (fnd == 0)  
            printf1("Cannot find: %s\n",p);

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

int tsel(void)        
{
    register int i;
    int err,m,n,opt,row,col,ivflg;
    double ws,tmp;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        printf1("No data matrix (command ignored).\n");
        return(0);
    }
    if (!strcmp(CmdBuf,"tsel=off")) {
        if (TSelFlg == 0)  
            printf1("No active case selection (command ignored).\n");
        else {
            tsel_off(2);
            set_cwt();      /* adjust weights */
        }
        return(0);
    }
    if (TSelFlg)
        tsel_off(1);        /* free previously active tsel command */
    else {
        if (EDAvail)  
            edat_off(1);    /* turn off episode data */
        if (SeqDN)  
            seq_afree(1);   /* turn off sequence data */
        if (GD_PERM == 0)
            gdd_free(1);    /* free relational data */
        repsel_off();       /* turn off repsel */
    }

    opt = 0;        /* check for standard expression */
    ivflg = 0;

    if ((n = v_parse(CmdBuf + 5,0)) < 0 || ESCnt <= 0) {

        /* if not a standard expression check for matrix expression */
           
        if (get_mexpr(CmdBuf + 5,&row,&col,&ivflg)) {
            printf1("Error: can't evaluate the tsel expression.\n");
            goto TSELFin;
        }
        if (row != NOCDM || col != 1) {
            printf1("Error: matrix expression must result in NOC x 1 matrix.\n");
            printf1("Found a %d x %d matrix.\n",row,col);
            goto TSELFin;
        }
        if (alloc_actmp(NOCDM + 1))
            goto TSELFin;   

        if (eval_mexpr(0,row,col,AcTmp,0,AcTmp,NULL,NULL,-1,0.0)) {
            printf1("Error: can't evaluate the tsel expression.\n");
            goto TSELFin;
        }   
        opt = 1;
    }

    if (opt == 0) {

        /* check for type 2 operators */

        for (i = 0; i < ESCnt; ++i) {
            n = iabs(ESTyp[i]);
            if (n >= OPT2A && n < OPT2B) {
                printf1("Error: tsel expressions may not contain type 2 operators.\n");
                goto TSELFin; 
            }
        }
    }
    printf1("New temporary case selection: %s\n",CmdBuf + 5);

    if (!(TSelect = (int *)calloc(NOCDM,sizeof(int)))) {
        p_err(-2,1);
        goto TSELFin;
    }
    TSelectA = NOCDM;
    memrq(TSelectA,sizeof(int));
           
    /* create the select indicator variables */

    m = 0;
    ws = 0.0;
    for (i = 0; i < NOCDM; ++i) {

        if (opt == 0) {
            n = v_eval1(i,ESCnt,ESTyp,ESVal,ESIdx,&tmp,0,0,0,0,0);
            if (n) {
                printf1("Error: can't evaluate tsel expression in case %d.\n",i + 1);
                prn_emsg2(n);
                goto TSELFin;
            }
        }
        else  
            tmp = AcTmp[i + 1];

        if (fabs(tmp) > EPSI2) {
            TSelect[m] = i;
            if (WIVar >= 0)  
                ws += get_data(WIVar,i);
            m++;
        }
    }
    if (m == 0) {
        printf1("Error: number of selected cases is zero.\n");
        goto TSELFin;
    }
    if (WIVar >= 0 && ws <= EPSI2) {
        printf1("Error: sum of weights for selected cases is almost zero.\n");
        goto TSELFin;
    }
    printf1("Number of selected cases: %d\n",m);
    NOC = m;
    TSelFlg = 1;

    if (WIVar >= 0) {
        printf1("Corresponding sum of weights: %lg",ws);
        if (WNormFlag)
            WSumS = (double)NOC;

        if (WSumS > 0.0) {
            WSum = WSumS;
            WNorm = WSum / ws;
            printf1(" (adjusted to: %lg).\n",WSum);
        }
        else {
            WSum = ws;
            WNorm = 1.0;
            printf1("\n");
        }
    }
    else {
        WSum = (double)NOC;
        WNorm = 1.0;
        WSumS = 0.0;
    }
    err = 0;  

TSELFin:
    alloc_actmp(0);
    if (err) {
        if (TSelectA > 0) {
            free((char *)TSelect);
            memrq(-TSelectA,sizeof(int));
            TSelectA = 0;
        }
        set_cwt();
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  tsel_off(opt)   Turn off tsel, free memory. Reset NOC to NOCDM.         */
/*                  Also free the following data structures:                */
/*                  - episode data                                          */
/*                  - sequence data                                         */
/*                  - relational data                                       */

void tsel_off(int opt)        
{
    if (TSelFlg == 0)
        return;

    if (EDAvail)  
        edat_off(opt);      /* turn off episode data */
         
    if (SeqDN)  
        seq_afree(opt);     /* turn off sequence data */
         
    if (GD_PERM == 0)
        gdd_free(opt);      /* free relational data */

    if (TSelectA > 0) {
        free((char *)TSelect);
        memrq(-TSelectA,sizeof(int));
        TSelectA = 0;
    }
    TSelFlg = 0;
    NOC = NOCDM;
    if (opt) {
        printf1("Temporary case selection turned off.\n");
        if (opt == 2)  
            printf1("Number of cases: %d\n",NOC);
    }
}

/* ------------------------------------------------------------------------ */
/*  cwt()       cwt(wnorm=) = VName                                         */
/*              create case weights.                                        */
/*              Return 0 if OK, -1 if error.                                */

int cwt(void)        
{
    int err,iv;
    double wn;
    char *p,vname[VNLMax+1];

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        printf1("No data matrix (command ignored).\n");
        return(0);
    }
    WSumS = 0.0;
    if (WIVar >= 0) {
        printf1("Previously defined case weights turned off.\n");
        WIVar = -1;
        WSum = (double)NOC;
        WSumS = 0.0;
        WNorm = 1.0;
        WNormFlag = 0;
        iv = 0;
    }
    else
        iv = 1;

    if (!strcmp(CmdBuf,"cwt=off")) {
        if (iv)
            printf1("No case weights defined (command ignored).\n");
        return(0);
    }  
    WNormFlag = 0;

    p = CmdBuf + 3;
    if (*p == '(') {
        if (!strncmp(p,"(wnorm)",7)) {
            p += 7;
            WNormFlag = 1;
        }
        else {
            if (sscanf(p + 1,"wnorm=%lf",&wn) != 1 || wn <= 0.0) {
                p_err(-1,1);
                return(-1);
            }
            p = skip_blev(p);
            WSumS = wn;
        }
    }
    if (*p++ != '=') {
        p_err(-1,1);
        return(-1);
    }
    if ((iv = get_vidx(p,vname)) < 0) {
        printf1("Error: undefined variable.\n");
        return(-1);
    }
    WIVar = iv;
    set_cwt();
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

int set_cwt(void)         
{
    register int i;
    double tmp,ws;

    WSum = (double)NOC;
    WNorm = 1.0;

    if (WIVar < 0)  
        return(0);
       
    ws = 0.0;
    for (i = 0; i < NOC; ++i) {
        tmp = get_data(WIVar,i);
        if (tmp < 0.0) {
            printf1("Error: found negative weight in case %d (command ignored).\n",i + 1);
            WIVar = -1;
            return(-1);
        }
        ws += tmp;
    }
    printf1("New case weights. Sum of weights: %lg",ws);
    if (ws <= EPSI) {
        printf1(" (weights will not be used).\n");
        WIVar = -1;
        return(-1);
    }
    if (WNormFlag)
        WSumS = (double)NOC;

    if (WSumS > 0.0) {
        WSum = WSumS;
        WNorm = WSum / ws;
        printf1(", adjusted to: %lg\n",WSum);
    }
    else {
        WSum = ws;
        WNorm = 1.0;            
        printf1("\n");
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  wr_sys()        Execute wsys command in CmdBuf.                         */
/*                  Syntax: wsys or wsys = fname                            */
/*                  Return 0 if OK, -1 if err.                              */

int wr_sys(void)
{
    FILE *fd;
    register int j;
    int err,r,len,nv;
    char *p,*fname;
    char *sysname = "tda.sys";

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        printf1("No data matrix (command ignored).\n");
        return(0);
    }
    p = CmdBuf + 4;

    if (!*p)
        fname = sysname;
    else if (*p++ == '=' && *p)
        fname = p;
    else {
        printf1("Syntax error.\n");
        return(-1);    
    }
    nv = 0;
    j = VIFirst;
    while (j >= 0) {
        if (VTyp[j] != 5)      /* skip type 5 variables */
            nv++;
        j = VNxt[j];
    }
    if (nv == 0) {
        printf1("No variables (command ignored).\n");
        return(0);
    }
    if (!(fd = fopen(fname,OPEN_WB))) {
        printf1("Can't create: %s\n",fname);
        return(-1);   
    }
    printf1("Writing system file: %s\n",fname);

    fprintf(fd,"TDA System File (%3.1f). ",TDA_Version);
#if TIME_ON
    prn_time(fd);
#else
    fprintf(fd,"\n");
#endif

    fprintf(fd,"%d %d\n",NOCDM,nv);

    j = VIFirst;
    while (j >= 0) {
        if (VTyp[j] != 5) {  
            fprintf(fd,"%s %s ",VName[j],VDef[j]);
            if (VLabel[j] != NULL)
                fprintf(fd,"%s ",VLabel[j]);
            fprintf(fd,"\n");
            j = VNxt[j];
        }
    }
    j = VIFirst;
    while (j >= 0) {
        if (VTyp[j] != 5)    
            fprintf(fd,"%d %d %d %d\n",VTyp[j],VSLen[j],VPFmt1[j],VPFmt2[j]);
        j = VNxt[j];
    }

    j = VIFirst;
    while (j >= 0) {
        if (VTyp[j] != 5) {  
            len = get_slen(j,NOCDM);
            r = fwrite(VDPtr[j],len,1,fd);
            if (r != 1) {
                printf1("Error in writing system file %s\n",fname);
                goto WSYSFin;
            }
            j = VNxt[j];
        }
    }
    printf1("Data matrix (%d cases, %d variables) written to: %s\n",NOCDM,nv,fname);
    err = 0;   

WSYSFin:
    fclose(fd);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rd_sys()    Read TDA system file. CmdBuf: rsys or rsys=fname.           */
/*              Return 0 if OK, -1 if error.                                */

int rd_sys(void)           
{
    FILE *fd;
    register int j;
    int nn,l,l1,l2,l3,nv,m,n,w1,w2,err;
    char *p,*q,*fname,*sysname = "tda.sys";
    char vname[VNLMax + 1];
    char vdef[RLMaxDef];
    char vlabel[VLLMax + 1];

    nn = 20000;                     /* length of read buffer */

    if (check_cmd(1))
        return(-1);

    if (DMDef) {
        printf1("Error: a data matrix already exists (command ignored).\n");
        return(0);
    }
    p = CmdBuf + 4;

    if (!*p)
        fname = sysname;
    else if (*p++ == '=' && *p)
        fname = p;
    else {
        printf1("Syntax error.\n");
        return(-1);    
    }
    if (!(fd = fopen(fname,OPEN_RB))) {
        printf1("Can't open: %s\n",fname);
        return(-1);   
    }
    printf1("Reading system file: %s\n",fname);

    VLabelLen = 0;
    if (alloc_acc(nn + 1)) {
        err = -1;
        goto RSYSFin; 
    }
    err = 0;
    if (!fgets(AcC,nn,fd)) {
        err = 1;
        goto RSYSFin;
    }  
    if (strncmp(AcC,"TDA System File",15)) {
        err = -2;
        goto RSYSFin;
    }
    printf1("Identification: %s",AcC);

    if (!fgets(AcC,nn,fd)) {
        err = 2;
        goto RSYSFin;
    }  
    if (sscanf(AcC,"%d %d",&n,&nv) != 2 || n <= 0 || nv <= 0) {
        err = 3;
        goto RSYSFin;
    }       
    if (nv > MaxNV) {
        err = -3;
        goto RSYSFin;
    }
    NOCMaxA = NOCDM = NOC = n;
    NVAR = 0;

    /* get variable names and definitions */
            
    for (j = 0; j < nv; ++j) {

        if (!fgets(AcC,nn,fd)) {
            err = 4;
            goto RSYSFin;
        }
        p = AcC;
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
        if (!(VName[j] = (char *)calloc(l1,sizeof(char)))) { 
            p_err(-2,1);
            err = -1;
            goto RSYSFin;
        }
        l2++;
        if (!(VDef[j] = (char *)calloc(l2,sizeof(char)))) { 
            p_err(-2,1);
            free(VName[j]);
            err = -1;
            goto RSYSFin;
        }
        memrq(l1,1);
        strcpy(VName[j],vname);
        memrq(l2,1);
        strcpy(VDef[j],vdef);
     
        if (l3 > 0) {
            if (VLabelLen < l3)
                VLabelLen = l3;
            l3++;
            if (!(VLabel[j] = (char *)calloc(l3,sizeof(char)))) { 
                p_err(-2,1);
                free(VName[j]);
                memrq(-l1,1);
                free(VDef[j]);
                memrq(-l2,1);
                VLabel[j] = NULL;
                err = -1;
                goto RSYSFin;
            }
            memrq(l3,1);
            strcpy(VLabel[j],vlabel);
        }
        VAlloc[j] = 1;
        if (VIFirst < 0)
            VILast = VIFirst = j;
        else {
            VNxt[VILast] = j; 
            VILast = j;
        }
        NVAR++;
    }
    for (j = 0; j < NVAR; ++j) {

        if (!fgets(AcC,RLMaxDef,fd)) {
            err = 7;
            goto RSYSFin;
        }
        if (sscanf(AcC,"%d %d %d %d",&n,&m,&w1,&w2) != 4) {
            err = 8;
            goto RSYSFin;
        }
        if (m >= 0 && m != 0 && m != 1 && m != 2 && m != 4 && m != 5 && m != 8) {
            err = 9;
            goto RSYSFin;
        }
        VTyp[j] = (char)n;
        VTypA[j] = 1;
        VSLen[j] = (short)m;
        makefmt(&w1,&w2,VPFmtS[j],0,' ',0);
        VPFmt1[j] = (short)w1;
        VPFmt2[j] = (short)w2;
    }
    newline();           
    prn_var(VIFirst);
    newline();

    /* allocate memory for data */

    if (alloc_vdat(0,1)) {
        p_err(-2,1);
        err = -1;
        goto RSYSFin;
    }
        
    /* read data */

    for (j = 0; j < NVAR; ++j) {

        l = get_slen(j,NOCDM);
        if ((n = fread(VDPtr[j],l,1,fd)) != 1) {
            err = 10;
            goto RSYSFin;
        }
    }
    DMDef = 1;
    printf1("Created a data matrix with %d variables and %d cases.\n",NVAR,NOCDM);

RSYSFin:
    fclose(fd);
    if (err > 0) {
        printf1("Error (%d) in reading system file.\n",err);
        if (err == 10)  
            printf1("R=%d NOC=%d L=%d J=%d\n",n,NOCDM,l,j);
    }
    else if (err == -2)
        printf1("Probably not a TDA system file.\n");
    else if (err == -3) {
        printf1("Error: system file contains %d variables.\n",nv);    
        printf1("More than the current maximum: %d.\n",MaxNV);
    }
    if (err) {
        clear_dm();            /* clear all variables */
        NOCMaxA = NOCDM = NOC = 0;
        NVAR = 0;
        err = -1;
    }
    alloc_acc(0);
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

int dblock(void)           
{
    register int i,j,k,l;
    int err,idx;
    double tmp;

    idx = err = -1;
    if (check_cmd(1))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    if (parm(CmdBuf + 6,4,0))  
        goto DBLOCKFin;
       
    tsel_off(0);                            /* turn off tsel */
    repsel_off();                           /* turn off repsel */

    if (NOC != DBlckPtrA) {
        if (dblock_alloc(NOC))  
            goto DBLOCKFin;
    }
    if (PMNV > 0) {
        if (alloc_actmp(PMNV))
            goto DBLOCKFin;
        for (j = 0; j < PMNV; ++j)
            AcTmp[j] = get_data(PMVIdx[j],0);
        k = 0;
    }
    else
        k = -1;

    if (PMatNameFlg) {
        if ((idx = mat_newmat(PMatName,NOC,1)) >= 0)  
            mdefcpy(MatDef[idx],CmdBuf);
    }
    for (i = 0; i < NOC; ++i) {
        if (PMNV == 0) {                    /* block = case */
            k++;
        }
        else {                              /* create blocks */
            l = 0;
            for (j = 0; j < PMNV; ++j) {
                tmp = get_data(PMVIdx[j],i);
                if (tmp != AcTmp[j]) {
                    l = 1;
                    AcTmp[j] = tmp;
                }
            }
            if (l)
                k++;
        }
        DBlckPtr[i] = k;
        if (idx >= 0)
            MatVal[idx][i + 1] = (double)(k + 1);            
    }
    BNOC = k + 1;
    printf2("Found %d block(s).\n",BNOC);
    err = 0;

DBLOCKFin:
    p_clean();
    if (err)
        BNOC = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dblock_alloc(n)     if n > 0 allocate DBlckPtr, else free.              */
/*                      Return 0 if OK, -1 if error.                        */

int dblock_alloc(int n)    
{
    if (n == 0 || n != DBlckPtrA) {
        if (DBlckPtrA > 0) {
            free((char *)DBlckPtr);
            memrq(-DBlckPtrA,sizeof(int));
            BNOC = DBlckPtrA = 0;
        }
    }
    if (n <= 0)
        return(0);

    if (!(DBlckPtr = (int *)calloc(n,sizeof(int)))) {
        p_err(-2,1);
        return(-1);           
    }
    DBlckPtrA = n;
    memrq(n,sizeof(int));
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

int repsel(void)           
{
    register int i,j,k,l,kk,k1;
    int err,row,col,n,m,mm,nn,ivflg;
    char *p;

    err = -1;
    ivflg = 0;

    if (check_cmd(1))
        return(-1);

    if (parm(CmdBuf + 6,8,1))  
        goto REPSELFin;

    if (PMRHSTRA < 1)
        goto REPSELFin;

    repsel_off();

    if (!strcmp(PMRHSTR,"off")) {
        err = 0;
        goto REPSELFin;
    }
    if (DMDef == 0) {
        p_err(-18,1);
        goto REPSELFin;
    }
    if (BNOC == 0) {
        printf1("Error: need a previous dblock command.\n");
        goto REPSELFin;
    }
    p = PMRHSTR;                                    /* get expression */
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto REPSELFin;
    if (col != 1) {
        mat_err(7);
        goto REPSELFin;
    }
    if (*p) {
        mat_err(0);
        goto REPSELFin;
    }
    if (BNOC != row) {
        printf1("Error: number of blocks is %d, found %d rows.\n",BNOC,row);
        goto REPSELFin;
    }
    tsel_off(0);                /* turn off tsel */

    nn = 0;                     /* calc new number of cases */
    m = DBlckPtr[0];   
    kk = k = 0;
    for (i = 1; i <= BNOC; ++i) {
        n = (int)MX[0][i];
        while (k <= NOC) {
            if (k < NOC)
                mm = DBlckPtr[k];
            else 
                mm++;

            if (mm == m) {
                if (PMID >= 0) {        /* id variable defined */
                    if (fabs(get_data(PMID,k) - MX[0][i]) <= EPSI1)
                        nn++;
                }
            }
            else {
                if (PMID < 0) {
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
        printf1("Error: selected number of cases is zero.\n");
        goto REPSELFin;
    }
    if (repsel_alloc(nn))  
        goto REPSELFin;

    m = DBlckPtr[0];
    kk = l = k = 0;
    for (i = 1; i <= BNOC; ++i) {
        n = (int)MX[0][i];
        while (k <= NOC) {
            if (k < NOC)
                mm = DBlckPtr[k];
            else
                mm++;

            if (mm == m) {
                if (PMID >= 0) {
                    if (fabs(get_data(PMID,k) - MX[0][i]) <= EPSI1)
                        REPSelect[l++] = k;
                }
            }
            else {
                if (PMID < 0) {
                    for (j = 1; j <= n; ++j) {
                        k1 = kk;
                        while (k1 < k) {
                            REPSelect[l++] = k1;
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
    printf1("REPSelect:\n");
    for (l = 0; l < nn; ++l)
        printf1("l=%2d sel=%d\n",l,REPSelect[l]);
    if (PMID >= 0)
        printf1("ID variable: %s\n",VName[PMID]);
    ****/

    REPNOCOld = NOC;            /* save old NOC */
    NOC = nn;                   /* set new NOC */
    REPSelFlg = 1;              /* make active */
    printf2("Selected %d cases.\n",NOC);

    set_cwt();                  /* adjust weights if defined */
    err = 0;

REPSELFin:
    mx_free();
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  repsel_off()        Turns off the repsel command.                       */
/*                      Return 0 if OK, -1 if error.                        */

int repsel_off(void)    
{
    if (REPSelFlg) {
        if (DMDef) {
            if (REPNOCOld != NOCDM)  
                gerr_exit(100);
            NOC = REPNOCOld;
        }
        printf2("Repsel turned off.\n");
    }
    REPNOCOld = 0;
    repsel_alloc(0);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  repsel_alloc(n)     if n > 0 allocate REPSelect, else free.             */
/*                      Return 0 if OK, -1 if error.                        */

int repsel_alloc(int n)    
{
    REPSelFlg = 0;
    if (n == 0 || n != REPSelectA) {
        if (REPSelectA > 0) {
            free((char *)REPSelect);
            memrq(-REPSelectA,sizeof(int));
            REPSelectA = 0;
        }
        if (n == 0)
            return(0);
    }
    else
        return(0);

    if (!(REPSelect = (int *)calloc(n,sizeof(int)))) {
        p_err(-2,1);
        return(-1);           
    }
    REPSelectA = n;
    memrq(n,sizeof(int));
    return(0);
}


