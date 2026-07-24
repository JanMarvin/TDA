/****************************************************************************/
/*  t_edat                                                                  */
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
#include "t_eval.h"
#include "t_eval1.h"
#include "t_gdat.h"
#include "t_freq.h"
#include "t_var.h" 
#include "t_parm.h" 
#include "t_alloc.h" 
#include "t_gf.h" 

/*  functions in t_edat.c */

int check_edj(int j);
int edef(void);
void prn_eed(char *s);
void edat_off(int opt);
int check_nevar(int nidx);
int def_edat(int mode);
void prn_edef(void);
void edef_info(void);
int check_edat(void);
int alloc_edat(int mode);
void prn_edat(void);
int get_edat(int i,int *sn,int *org,int *des,double *ts,double *tf);
int get_id(int icase,int *err);
int sort_ed(int opt);    
int tscomp(const void *arg1, const void *arg2);
int tfcomp(const void *arg1, const void *arg2);
int get_spell(int init,int *icase,int *sn,int *org,int *des,
    double *ts,double *tf,int *spl,int *nspl);
int check_nedat(void);
int epdat(void);
void edtda(char *fname,int noc);
int epsdat(void);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int EDAvail = 0;        /* 1 if episode data allocated and available        */
int MEFlg = 0;          /* set for multiepisode data                        */

char *PORGPtr = NULL;   /* pointer to org parameter                         */
char *PDESPtr = NULL;   /* pointer to des parameter                         */
char *PTSPtr  = NULL;   /* pointer to  ts parameter                         */
char *PTFPtr  = NULL;   /* pointer to  tf parameter                         */
char *PIDPtr  = NULL;   /* pointer to  id parameter                         */
char *PSNPtr  = NULL;   /* pointer to  sn parameter                         */

char PORGStr[120];
char PDESStr[120];
char PTSStr[120];
char PTFStr[120];
char PIDStr[120];
char PSNStr[120];

/*  The following variables have the following values:                      */
/*  >= 0, then it is the internal variable number                           */
/*    -1, then it is a numerical constant                                   */
/*    -2, it is a general expression                                        */

int IDVar = -2;         /* type of id expression                            */
int SNVar = -2;         /* type of id expression                            */
int TSVar = -2;         /* type of id expression                            */
int TFVar = -2;         /* type of id expression                            */
int ORGVar = -2;        /* type of id expression                            */
int DESVar = -2;        /* type of id expression                            */

short IDESCnt = 0;      /* parser stack for id command                      */
int *IDESTyp;
double *IDESVal;

short SNESCnt = 0;      /* parser stack for sn command                      */
int *SNESTyp;
double *SNESVal;

short TSESCnt = 0;      /* parser stack for ts command                      */
int *TSESTyp;
double *TSESVal;

short TFESCnt = 0;      /* parser stack for tf command                      */
int *TFESTyp;
double *TFESVal;

short ORGESCnt = 0;     /* parser stack for org command                     */
int *ORGESTyp;
double *ORGESVal;

short DESESCnt = 0;     /* parser stack for des command                     */
int *DESESTyp;
double *DESESVal;

/*  We assume spell numbers must be not less than 1, origin and destination */
/*  state number must be not less than 0.                                   */

int MaxTran = 0;        /* max number of (sn,org,des) combinations          */

int MaxSnn = 0;         /* highest spell number                             */
int MaxOrg = 0;         /* highest origin state number                      */
int MaxOrg1 = 0;        /* MaxOrg1 = MaxOrg + 1                             */
int MaxDes = 0;         /* highest destination state number                 */
int MaxDes1 = 0;        /* MaxDes1 = MaxDes + 1                             */
int NTran = 0;          /* number of (sn,org,des) combinations              */
int *SnTran;            /* this is an array of lengt NTran for all          */
int SnTranA = 0;
int *OrgTran;           /* combinations of (sn,org,des) in the input data,  */
int OrgTranA = 0;
int *DesTran;           /* sorted according to sn,org,des.                  */
int DesTranA = 0;
int NTran1 = 0;         /* number of (sn,org,des) transitions               */
int *SnTran1;           /* this is an array of lengt NTran1 for all         */
int SnTran1A = 0;
int *OrgTran1;          /* transitions (sn,org,des) in the input data,      */
int OrgTran1A = 0;
int *DesTran1;          /* sorted according to sn,org,des.                  */
int DesTran1A = 0;

int *TranPtr;           /* pointer to entries in SnTran1,OrgTran1,DesTran1  */
int TranPtrA = 0;

int *TranNE;            /* number of episodes                               */
int TranNEA = 0;
float *TranWE;          /* weighted number of episodes                      */
int TranWEA = 0;
float *TranMD;          /* weighted mean duration                           */
int TranMDA = 0;
float *TranTSM;         /* minimum of starting times                        */
int TranTSMA = 0;
float *TranTFM;         /* maximum of ending times                          */
int TranTFMA = 0;

float *TWFreq;          /* weighted number of uncensored episodes           */
int TWFreqA = 0;
float *TWDur;           /* summed durations                                 */
int TWDurA = 0;

int ESortFlg = 0;       /* 1 if sorted                                      */
int *TSIdx;             /* pointer for sorting according to ts              */
int *TFIdx;             /* pointer for sorting according to tf              */
int TSSortA = 0;        /* set to allocated memory for TSIdx                */
int TFSortA = 0;        /* set to allocated memory for TFIdx                */
/* ------------------------------------------------------------------------ */
int NSP = 0;            /* number of split variables                        */
short SPLITVar[MaxSPL]; /* indices of split variables                       */
int NSplits = 0;        /* number of splits, set in get_spell()             */

double SPTTV[MaxSPL+3]; /* array of times used in get_spell()               */
short SPIDX[MaxSPL+3];  /* array of indices used in get_spell()             */

int NVAR5 = 0;          /* number of type 5 variables                       */
int NSPLEvalErr = 0;    /* number of errors when evaluating type 5 vars     */

/* ------------------------------------------------------------------------ */
/*  check_edj(j)    Check whether variable j is needed for episode data.    */
/*                  Return 1 if needed, otherwise 0.                        */

int check_edj(int j)
{
    register int i;

    if (EDAvail == 0)
        return(0);
    if (j == IDVar || j == SNVar || j == TSVar || j == TFVar || j == ORGVar || j == DESVar)
        return(1);
    for (i = 0; i < NSP; ++i) {
        if (j == SPLITVar[i])
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  edef()      Setup episode data. Command in CmdBuf.                      */
/*              edef(org,des,ts,tf,id,sn,maxtran,vardef...).                */
/*              Return 0 if OK, -1 if error.                                */
 
int edef(void)
{
    register int i;
    int err,n,j,l,nidx,nb;
    char c,*p,*q;                     

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (!strcmp(CmdBuf,"edef")) {
        if (EDAvail == 0)
            printf1("No episode data structure defined.\n");
        else
            edef_info();
        return(0);
    }
    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    if (EDAvail)  
        edat_off(1);

    MaxTran = MaxTranDef;       /* default max number of transitions */
    NSP = 0;                    /* number of split variables */
    nidx = get_nidx();          /* index to first new variable */
    NVAR5 = 0;                  /* number of type 5 variables */

    PORGPtr = NULL; 
    PDESPtr = NULL; 
    PTSPtr  = NULL; 
    PTFPtr  = NULL; 
    PIDPtr  = NULL; 
    PSNPtr  = NULL; 

    p = CmdBuf + 4;
    while (*++p) {
        if (sscanf(p,"maxtran=%d",&n) == 1 && n > 0) {
            MaxTran = n;
            p = skip_int(p + 8);
        }
        else if (!strncmp(p,"org=",4) && *(p + 4)) {     
            PORGPtr = p + 4;
            q = skip_expr(PORGPtr);
            strncpy(PORGStr,PORGPtr,(int)(q - PORGPtr));
            p = skip_nc(p + 4);
        }
        else if (!strncmp(p,"des=",4) && *(p + 4)) {     
            PDESPtr = p + 4;
            q = skip_expr(PDESPtr);
            strncpy(PDESStr,PDESPtr,(int)(q - PDESPtr));
            p = skip_nc(p + 4);
        }
        else if (!strncmp(p,"ts=",3) && *(p + 3)) {     
            PTSPtr = p + 3;
            q = skip_expr(PTSPtr);
            strncpy(PTSStr,PTSPtr,(int)(q - PTSPtr));
            p = skip_nc(p + 3);
        }
        else if (!strncmp(p,"tf=",3) && *(p + 3)) {     
            PTFPtr = p + 3;
            q = skip_expr(PTFPtr);
            strncpy(PTFStr,PTFPtr,(int)(q - PTFPtr));
            p = skip_nc(p + 3);
        }
        else if (!strncmp(p,"id=",3) && *(p + 3)) {     
            PIDPtr = p + 3;
            q = skip_expr(PIDPtr);
            strncpy(PIDStr,PIDPtr,(int)(q - PIDPtr));
            p = skip_nc(p + 3);
        }
        else if (!strncmp(p,"sn=",3) && *(p + 3)) {     
            PSNPtr = p + 3;
            q = skip_expr(PSNPtr);
            strncpy(PSNStr,PSNPtr,(int)(q - PSNPtr));
            p = skip_nc(p + 3);
        }
        else if (!strncmp(p,"split=",6)) {     
            q = p;
            p = get_nvia(p + 6,&n,1,&nb);
            if (n <= 0 || nb) {
                printf1("Check: %s\n",q);
                goto EDEFFin;
            }
            if (n > MaxSPL) {
                printf1("Error: exceeded max number of split variables.\n");
                goto EDEFFin;
            }
            j = 0;
            for (i = 0; i < n; ++i) {
                SPLITVar[i] = VLVIdx[i];
                if (VTyp[VLVIdx[i]] == 5)
                    j = 1;
            }
            if (j) {
                printf1("Error: cannot use type 5 variables for episode splitting.\n");
                goto EDEFFin;
            }
            NSP = n;
        }
        else if ((l = get_vnlen(p)) > 0) {       /* variable */
            q = skip_nc(p + l);
            c = *q;
            *q = '\0';
            if (save_var(p,1))
                goto EDEFFin;
            NVAR5++;
            *q = c;
            p = q; 
        }

        if (*p != ',' && *p != ')') {
            prn_eed(p);
            goto EDEFFin;
        }
        *p = '\0';
    }

    if (PORGPtr == NULL || PDESPtr == NULL || PTSPtr == NULL || PTFPtr == 0) {
        printf1("Error: need org, des, ts, and tf parameters.\n");
        goto EDEFFin;
    }
    if (PIDPtr != NULL && PSNPtr != 0)
        MEFlg = 1;

    else if (PIDPtr != NULL || PSNPtr != 0) {
        printf1("Error: need id and sn parameters for multi-episode data.\n");
        goto EDEFFin;
    }
    else
        MEFlg = 0;
   
    err = def_edat(1);              /* check definition */ 
    if (err)
        goto EDEFFin;
    prn_edef();                     /* print definition */

    if (NSP > 0) {
        printf1("Episode splitting with variable(s): %s",VName[SPLITVar[0]]);
        for (j = 1; j < NSP; ++j)
            printf1(",%s",VName[SPLITVar[j]]);
        printf1("\n");
    }
    newline();
    err = check_edat();             /* check data */ 
    if (err)
        goto EDEFFin;
    prn_edat();                     /* print info */

    if (NVAR5 > 0) {                /* check new variables */
        err = check_nevar(nidx);
        if (err)
            goto EDEFFin;
    }
    err = check_nedat();            /* check correct data generation */
    if (err)
        goto EDEFFin;

    EDAvail = 1;

EDEFFin:
    p_clean();
    if (err) {
        edat_off(0);  
        printf1("Could not create new episode data.\n");
    }
    else
        printf1("Successfully created new episode data.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_eed.    Print error message.                                        */

void prn_eed(char *s)
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
/*  edat_off(opt)  free all data structures for episode data.               */
/*                 if opt != 0 print message.                               */

void edat_off(int opt)
{
    register int j;

    if (ESortFlg)
        sort_ed(0);
    alloc_edat(0);
    def_edat(0);

    j = VIFirst;        /* remove type 5 variables */
    while (j >= 0) {
        if (VTyp[j] == 5) {
            clear_var(j);
            j = VIFirst;
        }
        else
            j = VNxt[j];
    }
    if (opt && EDAvail)
        printf1("Deallocated episode data structures.\n");
    EDAvail = 0;
}

/* ------------------------------------------------------------------------ */
/*  check_nevar()   check new type 5 variables, begin with index vidx.      */
/*                  return 0 if OK, -1 if error.                            */

int check_nevar(int nidx)
{
    int err;

    err = -1;
    printf1("\nType 5 variables (dependent on current episode data).\n\n");

    prn_var(nidx);          /* print variables */

    if (alloc_vdat(nidx,1)) {
        printf1("\nInsufficient memory for new variables.\n");
        goto NEVARFin;
    }
    err = 0;

NEVARFin:
    if (err) {
        clear_avar(nidx);
        NVAR5 = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  def_edat(mode)  If mode != 0 check and allocate new definition of       */
/*                  episode data, otherwise remove current definition.      */
/*                                                                          */
/*                  This function also sets the type of the expressions.    */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */
   
int def_edat(int mode)
{
    register int i;
    int n,err;

    err = 0;
    if (mode == 0) {
        NSP = MEFlg = 0;
        IDVar = SNVar = TSVar = TFVar = ORGVar = DESVar = -2;
        goto EDATFin;
    }
    err = -1;

    /* save org */

    if ((n = v_parse(PORGPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d): %s\n",n,PORGPtr);
        if (n < 0)
            prn_emsg1(n);
        goto EDATFin;
    }
    if (!(ORGESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
        p_err(-2,1);
        goto EDATFin;
    }
    if (!(ORGESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
        p_err(-2,1);
        free((char *)ORGESTyp);
        goto EDATFin;
    }
    memrq(ESCnt,sizeof(int) + sizeof(double));
    ORGESCnt = ESCnt;

    for (i = 0; i < ESCnt; ++i) {
        ORGESTyp[i] = ESTyp[i];
        ORGESVal[i] = ESVal[i];
    }
    if (ESCnt == 1) {
        if (ESTyp[0] == 0)
            ORGVar = -1;
        else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
            ORGVar = ESTyp[0] - VOFFS;
    }
         
    /* save des */

    if ((n = v_parse(PDESPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d): %s\n",n,PDESPtr);
        if (n < 0)
            prn_emsg1(n);
        goto EDATFin;
    }
    if (!(DESESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
        p_err(-2,1);
        goto EDATFin;
    }
    if (!(DESESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
        p_err(-2,1);
        free((char *)DESESTyp);
        goto EDATFin;
    }
    memrq(ESCnt,sizeof(int) + sizeof(double));
    DESESCnt = ESCnt;

    for (i = 0; i < ESCnt; ++i) {
        DESESTyp[i] = ESTyp[i];
        DESESVal[i] = ESVal[i];
    }
    if (ESCnt == 1) {
        if (ESTyp[0] == 0)
            DESVar = -1;
        else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
            DESVar = ESTyp[0] - VOFFS;
    }
        
    /* save ts */

    if ((n = v_parse(PTSPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d): %s\n",n,PTSPtr);
        if (n < 0)
            prn_emsg1(n);
        goto EDATFin;
    }
    if (!(TSESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
        p_err(-2,1);
        goto EDATFin;
    }
    if (!(TSESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
        p_err(-2,1);
        free((char *)TSESTyp);
        goto EDATFin;
    }
    memrq(ESCnt,sizeof(int) + sizeof(double));
    TSESCnt = ESCnt;

    for (i = 0; i < ESCnt; ++i) {
        TSESTyp[i] = ESTyp[i];
        TSESVal[i] = ESVal[i];
    }
    if (ESCnt == 1) {
        if (ESTyp[0] == 0)
            TSVar = -1;
        else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
            TSVar = ESTyp[0] - VOFFS;
    }
        
    /* save tf */

    if ((n = v_parse(PTFPtr,0)) < 0 || ESCnt <= 0) {
        printf1("Syntax or reference error (%d): %s\n",n,PTFPtr);
        if (n < 0)
            prn_emsg1(n);
        goto EDATFin;
    }
    if (!(TFESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
        p_err(-2,1);
        goto EDATFin;
    }
    if (!(TFESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
        p_err(-2,1);
        free((char *)TFESTyp);
        goto EDATFin;
    }
    memrq(ESCnt,sizeof(int) + sizeof(double));
    TFESCnt = ESCnt;

    for (i = 0; i < ESCnt; ++i) {
        TFESTyp[i] = ESTyp[i];
        TFESVal[i] = ESVal[i];
    }
    if (ESCnt == 1) {
        if (ESTyp[0] == 0)
            TFVar = -1;
        else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
            TFVar = ESTyp[0] - VOFFS;
    }
        
    /* save id */

    if (PIDPtr != NULL) {
        if ((n = v_parse(PIDPtr,0)) < 0 || ESCnt <= 0) {
            printf1("Syntax or reference error (%d): %s\n",n,PIDPtr);
            if (n < 0)
                prn_emsg1(n);
            goto EDATFin;
        }
        if (!(IDESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
            p_err(-2,1);
            goto EDATFin;
        }
        if (!(IDESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
            p_err(-2,1);
            free((char *)IDESTyp);
            goto EDATFin;
        }
        memrq(ESCnt,sizeof(int) + sizeof(double));
        IDESCnt = ESCnt;

        for (i = 0; i < ESCnt; ++i) {
            IDESTyp[i] = ESTyp[i];
            IDESVal[i] = ESVal[i];
        }
        if (ESCnt == 1) {
            if (ESTyp[0] == 0)
                IDVar = -1;
            else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
                IDVar = ESTyp[0] - VOFFS;
        }
    }
         
    /* save sn */

    if (PSNPtr != NULL) {
        if ((n = v_parse(PSNPtr,0)) < 0 || ESCnt <= 0) {
            printf1("Syntax or reference error (%d): %s\n",n,PSNPtr);
            if (n < 0)
                prn_emsg1(n);
            goto EDATFin;
        }
        if (!(SNESTyp = (int *)calloc(ESCnt,sizeof(int)))) { 
            p_err(-2,1);
            goto EDATFin;
        }
        if (!(SNESVal = (double *)calloc(ESCnt,sizeof(double)))) { 
            p_err(-2,1);
            free((char *)SNESTyp);
            goto EDATFin;
        }
        memrq(ESCnt,sizeof(int) + sizeof(double));
        SNESCnt = ESCnt;

        for (i = 0; i < ESCnt; ++i) {
            SNESTyp[i] = ESTyp[i];
            SNESVal[i] = ESVal[i];
        }
        if (ESCnt == 1) {
            if (ESTyp[0] == 0)
                SNVar = -1;
            else if (ESTyp[0] >= VOFFS && ESTyp[0] < COFFS)
                SNVar = ESTyp[0] - VOFFS;
        }
    }
    return(0);

EDATFin:

    if (IDESCnt > 0) {
        free((char *)IDESTyp);
        free((char *)IDESVal);
        memrq(-IDESCnt,sizeof(int) + sizeof(double));
        IDESCnt = 0;
    }
    if (SNESCnt > 0) {
        free((char *)SNESTyp);
        free((char *)SNESVal);
        memrq(-SNESCnt,sizeof(int) + sizeof(double));
        SNESCnt = 0;
    }
    if (TSESCnt > 0) {
        free((char *)TSESTyp);
        free((char *)TSESVal);
        memrq(-TSESCnt,sizeof(int) + sizeof(double));
        TSESCnt = 0;
    }
    if (TFESCnt > 0) {
        free((char *)TFESTyp);
        free((char *)TFESVal);
        memrq(-TFESCnt,sizeof(int) + sizeof(double));
        TFESCnt = 0;
    }
    if (ORGESCnt > 0) {
        free((char *)ORGESTyp);
        free((char *)ORGESVal);
        memrq(-ORGESCnt,sizeof(int) + sizeof(double));
        ORGESCnt = 0;
    }
    if (DESESCnt > 0) {
        free((char *)DESESTyp);
        free((char *)DESESVal);
        memrq(-DESESCnt,sizeof(int) + sizeof(double));
        DESESCnt = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_edef()  print definition of episode data.                           */
   
void prn_edef(void)
{
    printf1("Creating new ");
    if (MEFlg)
        printf1("multi-");
    else
        printf1("single ");
    printf1("episode data. Max number of transitions: %d.\n",MaxTran);
    printf1("Definition: ");
    if (MEFlg)  
        printf1("id=%s, sn=%s, ",PIDPtr,PSNPtr);
    printf1("org=%s, des=%s, ts=%s, tf=%s\n",PORGPtr,PDESPtr,PTSPtr,PTFPtr);
}

/* ------------------------------------------------------------------------ */
/*  edef_info()  info about definition of episode data.                     */
   
void edef_info(void)
{
    if (MEFlg)
        printf1("Multi-");
    else
        printf1("Single ");
    printf1("episode data. Max number of transitions: %d.\n",MaxTran);
    printf1("Definition: ");
    if (MEFlg)  
        printf1("id=%s, sn=%s, ",PIDStr,PSNStr);
    printf1("org=%s, des=%s, ts=%s, tf=%s\n\n",PORGStr,PDESStr,PTSStr,PTFStr);
    prn_edat();
}

/* ------------------------------------------------------------------------ */
/*  check_edat()    Check for 0 <= TS < TF, and create basic data           */
/*                  structures for episodes.                                */
/*                  Return 0 if successful, otherwise                       */
/*                  -1 if syntax error in expressions,                      */
/*                  -2 if insufficient memory                               */
/*                  -3 if negative starting time                            */
/*                  -4 if zero or negative duration                         */
/*                  -5 if negative origin or destination state              */
/*                  -6 if spell number less than 1                          */
/*                  -8 if number of transitions exceeds MaxTran             */
/*                  -9 if error in making frequency distribution            */
/*                 -10 if number of transitions is zero                     */
/*                                                                          */

int check_edat(void)
{
    register int i,j,k;
    int err,r,wflag,sn,sn1,org,org1,des,*idat,*tran,*freq,*bf;
    int idata,freqa,trana,bfa,wta,wfreqa;
    double x,ts,tf,*wt,*wfreq;
    float w,d;

    EDAvail = 0;
    wflag = idata = freqa = trana = bfa = wta = wfreqa = 0;
    err = -2;

    /* allocate memory for frequency distribution of transitions */

    if (!(idat = (int *)calloc(3 * NOC + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto CDFin;           
    }
    idata = 3 * NOC + 1;
    memrq(idata,sizeof(int));
              
    if (!(freq = (int *)calloc(MaxTran + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto CDFin;           
    }
    freqa = MaxTran + 1;
    memrq(freqa,sizeof(int));

    if (!(tran = (int *)calloc(3 * MaxTran + 1,sizeof(int)))) { 
        p_err(-2,1);
        goto CDFin;           
    }
    trana = 3 * MaxTran + 1;
    memrq(trana,sizeof(int));
              
    if (!(bf = (int *)calloc(MaxTran + 1,sizeof(int)))) { 
        p_err(-2,1);  
        goto CDFin;           
    }
    bfa = MaxTran + 1;
    memrq(bfa,sizeof(int));

    if (WIVar >= 0) {       /* if case weights */

        wflag = 1;

        if (!(wt = (double *)calloc(NOC + 1,sizeof(double)))) { 
            p_err(-2,1);
            goto CDFin;           
        }
        wta = NOC + 1;
        memrq(wta,sizeof(double));

        if (!(wfreq = (double *)calloc(MaxTran + 1,sizeof(double)))) { 
            p_err(-2,1);
            goto CDFin;           
        }
        wfreqa = MaxTran + 1;
        memrq(wfreqa,sizeof(double));
    }
    err = 0;

    for (i = 0; i < NOC; ++i) {

        if (TSVar == -1)           /* get starting time */
            ts = TSESVal[0];
        else if (TSVar >= 0)
            ts = get_data(TSVar,i);
        else {
            r = v_eval1(i,TSESCnt,TSESTyp,TSESVal,ESIdx,&ts,0,0,0,0,0);
            if (r) {
                printf1("Error (%d): can't evaluate starting time in case %d.\n",r,i+1);
                prn_emsg2(r);
                err = -1;
                break;
            }
        }
        if (TFVar == -1)           /* get ending time */
            tf = TFESVal[0];
        else if (TFVar >= 0)
            tf = get_data(TFVar,i);
        else {
            r = v_eval1(i,TFESCnt,TFESTyp,TFESVal,ESIdx,&tf,0,0,0,0,0);
            if (r) {
                printf1("Error (%d): can't evaluate ending time in case %d.\n",r,i+1);
                prn_emsg2(r);
                err = -1;
                break;
            }
        }

        /* check starting and ending time */

        if (ts < 0.0) {
            printf1("Error: found negative starting time (%lg) in case %d.\n",ts,i+1);
            err = -3;
            break;
        }
        if (tf <= ts) {
            printf1("Error: found zero or negative duration (%lg) in case %d.\n",tf - ts,i+1);
            err = -4;
            break;
        }

        if (ORGVar == -1)           /* get origin state */
            org = (int)ORGESVal[0];
        else if (ORGVar >= 0)
            org = (int)get_data(ORGVar,i);
        else {
            r = v_eval1(i,ORGESCnt,ORGESTyp,ORGESVal,ESIdx,&x,0,0,0,0,0);
            if (r) {
                printf1("Error (%d): can't evaluate origin state in case %d.\n",r,i+1);
                prn_emsg2(r);
                err = -1;
                break;
            }
            org = (int)x;
        }

        if (DESVar == -1)           /* get destination state */
            des = (int)DESESVal[0];
        else if (DESVar >= 0)
            des = (int)get_data(DESVar,i);
        else {
            r = v_eval1(i,DESESCnt,DESESTyp,DESESVal,ESIdx,&x,0,0,0,0,0);
            if (r) {
                printf1("Error (%d): can't evaluate destination state in case %d.\n",r,i+1);
                prn_emsg2(r);
                err = -1;
                break;
            }
            des = (int)x;
        }
        if (org < 0 || des < 0) {
            printf1("Error: found negative origin (%d) or destination (%d) state in case %d.\n",org,des,i+1);
            err = -5;
            break;
        }
        idat[i * 3 + 1] = 1;
        idat[i * 3 + 2] = org;
        idat[i * 3 + 3] = des;

        if (MEFlg) {      /* multiepisode data */
   
            if (SNVar == -1)           /* get spell number */
                sn = (int)SNESVal[0];
            else if (SNVar >= 0)
                sn = (int)get_data(SNVar,i);
            else {
                r = v_eval1(i,SNESCnt,SNESTyp,SNESVal,ESIdx,&x,0,0,0,0,0);
                if (r) {
                    printf1("Error (%d): can't evaluate destination state in case %d.\n",r,i+1);
                    prn_emsg2(r);
                    err = -1;
                    break;
                }
                sn = (int)x;
            }
            if (sn < 1) {
                printf1("Error: found spell number (%d) less than 1 in case %d.\n",sn,i+1);
                err = -6;
                break;
            }
            idat[i * 3 + 1] = sn;
        }
        if (WIVar >= 0)                             /* get weights */
            wt[i + 1] = get_data(WIVar,i) * WNorm;
    }
    if (err)
        goto CDFin;

    /* make joint frequency table for spell number, origin and destination states */

    r = cfreq(3,NOC,idat,MaxTran,tran,freq,bf,wflag,wt,wfreq);
    if (r <= 0) {
        printf1("Error in checking transitions.\n");
        err = -9;
        if (r == -1) {
            err = -8;
            printf1("Number of transitions exceeds the defined maximum (%d).\n",MaxTran);
        }
        else if (r == -2)  
            printf1("Stack size is too small.\n");
        else if (r == -3)
            printf1("Insufficient memory for frequency distribution.\n");
        else 
            printf1("Error %d in frequency distribution.\n",r);
        goto CDFin;
    }

    /* -------------------------------------------------------------------- */
    /* NTran    = number of (sn,org,des) combinations                       */
    /* NTran1   = number of (sn,org,des) transitions                        */
    /* MaxSnn   = highest spell number                                      */
    /* MaxOrg   = highest origin state number                               */
    /* MaxOrg1  = MaxOrg + 1                                                */
    /* MaxDes   = highest destination state number                          */
    /* MaxDes1  = MaxDes + 1                                                */

    NTran = r;
    MaxSnn = MaxOrg = MaxDes = NTran1 = 0;
  
    for (i = 1; i <= NTran; ++i) {

        j = bf[i];

        sn  = tran[(j - 1) * 3 + 1];
        org = tran[(j - 1) * 3 + 2];
        des = tran[(j - 1) * 3 + 3];

        if (des != org)
            NTran1++;
        if (MaxSnn < sn)
            MaxSnn = sn;
        if (MaxOrg < org)
            MaxOrg = org;
        if (MaxDes < des)
            MaxDes = des;
    }
    if (NTran1 == 0) {
        printf1("Number of transitions is zero.\n");
        err = -10;
        goto CDFin;
    }
    MaxOrg1 = MaxOrg + 1;
    MaxDes1 = MaxDes + 1;

    if (alloc_edat(1)) {    /* allocate additional data structures */
        err = -2;
        goto CDFin;
    }

    /* -------------------------------------------------------------------- */
    /*  copy frequency distribution into                                    */
    /*  SnTran[i]   = spell number                                          */
    /*  OrgTran[i]  = origin state                                          */
    /*  DesTran[i]  = destination state                                     */
    /*  TranNE[i]   = number of episodes                                    */
    /*  TranWE[i]   = weighted number of episodes                           */
    /*  i = 0,...,NTran - 1                                                 */

    for (i = 0; i < NTran; ++i) {

        j = bf[i + 1];
        SnTran[i]  = tran[(j - 1) * 3 + 1];
        OrgTran[i] = tran[(j - 1) * 3 + 2];
        DesTran[i] = tran[(j - 1) * 3 + 3];
        TranNE[i]  = freq[j];

        if (WIVar >= 0)
            TranWE[i] = (float)wfreq[j];
        else
            TranWE[i] = (float)freq[j];
    }
   
    /* -------------------------------------------------------------------- */
    /*  create parallel arrays, but only for transitions, in                */
    /*  SnTran1[i]   = spell number                                         */
    /*  OrgTran1[i]  = origin state                                         */
    /*  DesTran1[i]  = destination state                                    */

    k = 0;
    for (i = 0; i < NTran; ++i) {

        sn  = SnTran[i];
        org = OrgTran[i];
        des = DesTran[i];

        if (des != org) {
            SnTran1[k] = sn;
            OrgTran1[k] = org;
            DesTran1[k] = des;
            k++;
        }
    }
    SnTran1[k] = OrgTran1[k] = DesTran1[k] = -1;


    /********************
    printf1("\nTran1...\n");
    for (i = 0; i < NTran1; ++i) 
        printf1("%2d %2d %2d \n",SnTran1[i],OrgTran1[i],DesTran1[i]);
    ***************/
    /* -------------------------------------------------------------------- */
    /*  create pointer to SnTran1,OrgTran1,DesTran1 in TranPtr.             */
    /*  TranPtr[sn,org] = TranPtr[(sn - 1) * MaxOrg1 + org]                 */
    /*  sn = 1,...,MaxSnn; org = 0,...,MaxOrg. (MaxOrg1 = MaxOrg + 1)       */

    j = MaxSnn * MaxOrg1 + 1;
    for (i = 0; i < j; ++i)
        TranPtr[i] = -1;

    sn = SnTran1[0] - 1;
    org = OrgTran1[0] - 1;

    for (i = 0; i < NTran1; ++i) {

        sn1 = SnTran1[i];
        org1 = OrgTran1[i];

        if (sn1 != sn || org1 != org) {
            TranPtr[(sn1 - 1) * MaxOrg1 + org1] = i;
            sn = sn1;
            org = org1;
        }
    }
    /***********************
    printf1("\nTranPtr\n");
    for (i = 1; i <= MaxSnn; ++i) {
        for (j = 0; j <= MaxOrg; ++j)
            printf1("%2d ",TranPtr[(i - 1) * MaxOrg1 + j]);
        printf1("\n");
    }
    *********/
    /* -------------------------------------------------------------------- */
    /*  TranMD[i]  = weighted mean duration for transition i                */
    /*  TranTSM[i] = minimum of starting times for transition i             */
    /*  TranTFM[i] = maximum of ending times for transition i               */
    /*  i = 0,...,NTran-1.                                                  */

    x = DBLMAX;
    for (j = 0; j < NTran; ++j) {
        TranTSM[j] = (float)x;
        TranTFM[j] = -1.0;
    }

    for (i = 0; i < NOC; ++i) {

        if (get_edat(i,&sn,&org,&des,&ts,&tf)) {
            p_err(-39,1);
            err = -7;
            goto CDFin;
        }   
        if (WIVar >= 0)                             /* get weights */
            x = get_data(WIVar,i) * WNorm;
        else
            x = 1.0;

        for (j = 0; j < NTran; ++j) {
            if (sn == SnTran[j] && org == OrgTran[j] && des == DesTran[j]) {
                TranMD[j] += (float)(x * (tf - ts));
                if (TranTSM[j] > (float)ts)
                    TranTSM[j] = (float)ts;
                if (TranTFM[j] < (float)tf)
                    TranTFM[j] = (float)tf;
                break;
            }
        }
    }

    for (j = 0; j < NTran; ++j)  
        TranMD[j] /= TranWE[j];

    /* -------------------------------------------------------------------- */
    /*  TWFreq[i]  = weighted number of uncensored episodes                 */
    /*               for transition i = 0,...,NTran1 - 1.                   */
    /*  TWDur[i]   = sum of weighted duration for uncensored episodes       */
    /*               in transition i.                                       */

    for (i = 0; i < NTran1; ++i) {

        sn = SnTran1[i];
        org = OrgTran1[i];
        des = DesTran1[i];
        w = d = 0.0;
        k = 0;
        for (j = 0; j < NTran; ++j) {
            if (sn == SnTran[j] && org == OrgTran[j]) {
                k++;
                d += (double)TranMD[j] * TranWE[j];
                if (des == DesTran[j])
                    w += TranWE[j];
            }
            else if (k)
                break;
        }
        TWDur[i] = d;
        TWFreq[i] = w;
    }
    /***********
    printf1("TWFreq TWDur\n");
    for (i = 0; i < NTran1; ++i)
        printf1("%lf %lf\n",TWFreq[i],TWDur[i]);
    printf1("\n");
    *********/

CDFin: 
    if (wfreqa) {
        free((char *)wfreq);
        memrq(-wfreqa,sizeof(double));
    }
    if (wta) {
        free((char *)wt);
        memrq(-wta,sizeof(double));
    }
    if (bfa) {
        free((char *)bf);
        memrq(-bfa,sizeof(int));
    }   
    if (trana) {
        free((char *)tran);
        memrq(-trana,sizeof(int));
    }   
    if (freqa) {
        free((char *)freq);
        memrq(-freqa,sizeof(int));
    }
    if (idata) {
        free((char *)idat);
        memrq(-idata,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_edat(mode)    If mode != 0 allocate memory for episode data       */
/*                      structures, else free previously allocated memory.  */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */
   
int alloc_edat(int mode)
{
    int err = 0;

    if (mode == 0)  
        goto AADFin;

    err = -1;
    if (!(SnTran = (int *)calloc(NTran,sizeof(int))))   
        goto AADFin;         
    memrq(NTran,sizeof(int));
    SnTranA = NTran;
              
    if (!(OrgTran = (int *)calloc(NTran,sizeof(int))))   
        goto AADFin;         
    memrq(NTran,sizeof(int));
    OrgTranA = NTran;
              
    if (!(DesTran = (int *)calloc(NTran,sizeof(int))))   
        goto AADFin;         
    memrq(NTran,sizeof(int));
    DesTranA = NTran;
              
    if (!(SnTran1 = (int *)calloc(NTran1 + 1,sizeof(int))))   
        goto AADFin;         
    SnTran1A = NTran1 + 1;
    memrq(SnTran1A,sizeof(int));
              
    if (!(OrgTran1 = (int *)calloc(NTran1 + 1,sizeof(int))))   
        goto AADFin;         
    OrgTran1A = NTran1 + 1;
    memrq(OrgTran1A,sizeof(int));
              
    if (!(DesTran1 = (int *)calloc(NTran1 + 1,sizeof(int))))   
        goto AADFin;         
    DesTran1A = NTran1 + 1;
    memrq(DesTran1A,sizeof(int));
              
    if (!(TranPtr = (int *)calloc(MaxSnn * MaxOrg1 + 1,sizeof(int))))   
        goto AADFin;         
    TranPtrA = MaxSnn * MaxOrg1 + 1;
    memrq(TranPtrA,sizeof(int));

    if (!(TranNE = (int *)calloc(NTran,sizeof(int))))   
        goto AADFin;         
    TranNEA = NTran;
    memrq(TranNEA,sizeof(int));

    if (!(TranWE = (float *)calloc(NTran,sizeof(float))))   
        goto AADFin;         
    TranWEA = NTran;
    memrq(TranWEA,sizeof(float));

    if (!(TranMD = (float *)calloc(NTran,sizeof(float))))   
        goto AADFin;          
    TranMDA = NTran;
    memrq(TranMDA,sizeof(float));

    if (!(TranTSM = (float *)calloc(NTran,sizeof(float))))   
        goto AADFin;          
    TranTSMA = NTran;
    memrq(TranTSMA,sizeof(float));

    if (!(TranTFM = (float *)calloc(NTran,sizeof(float))))   
        goto AADFin;          
    TranTFMA = NTran;
    memrq(TranTFMA,sizeof(float));

    if (!(TWFreq = (float *)calloc(NTran1,sizeof(float))))   
        goto AADFin;          
    TWFreqA = NTran1;
    memrq(TWFreqA,sizeof(float));

    if (!(TWDur = (float *)calloc(NTran1,sizeof(float))))   
        goto AADFin;          
    TWDurA = NTran1;
    memrq(TWDurA,sizeof(float));

    return(0);

AADFin:
    if (err)
        p_err(-2,1);

    if (TWDurA) {
        free((char *)TWDur);
        memrq(-TWDurA,sizeof(float));
        TWDurA = 0;
    } 
    if (TWFreqA) {
        free((char *)TWFreq);
        memrq(-TWFreqA,sizeof(float));
        TWFreqA = 0;
    }  
    if (TranTFMA) {
        free((char *)TranTFM);
        memrq(-TranTFMA,sizeof(float));
        TranTFMA = 0;
    }  
    if (TranTSMA) {
        free((char *)TranTSM);
        memrq(-TranTSMA,sizeof(float));
        TranTSMA = 0;
    }  
    if (TranMDA) {
        free((char *)TranMD);
        memrq(-TranMDA,sizeof(float));
        TranMDA = 0;
    }
    if (TranWEA) {
        free((char *)TranWE);
        memrq(-TranWEA,sizeof(float));
        TranWEA = 0;
    }
    if (TranNEA) {
        free((char *)TranNE);
        memrq(-TranNEA,sizeof(int));
        TranNEA = 0;
    }
    if (TranPtrA) {
        free((char *)TranPtr);
        memrq(-TranPtrA,sizeof(int));
        TranPtrA = 0;
    }
    if (DesTran1A) {
        free((char *)DesTran1);
        memrq(-DesTran1A,sizeof(int));
        DesTran1A = 0;
    }
    if (OrgTran1A) {
        free((char *)OrgTran1);
        memrq(-OrgTran1A,sizeof(int));
        OrgTran1A = 0;
    }
    if (SnTran1A) {
        free((char *)SnTran1);
        memrq(-SnTran1A,sizeof(int));
        SnTran1A = 0;
    }
    if (DesTranA) {
        free((char *)DesTran);
        memrq(-DesTranA,sizeof(int));
        DesTranA = 0;
    }
    if (OrgTranA) {
        free((char *)OrgTran);
        memrq(-OrgTranA,sizeof(int));
        OrgTranA = 0;
    }
    if (SnTranA) {
        free((char *)SnTran);
        memrq(-SnTranA,sizeof(int));
        SnTranA = 0;
    }
    NTran = NTran1 = MaxSnn = MaxOrg = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_edat()      print table with information about episodes.            */
   
void prn_edat(void)
{
    register int i,k;
    int sn,sn1,org,org1,nn,nna;
    double sum,suma;

    prnchar(' ',41,0);
    printf1("Mean\n");
    printf1("SN  Org Des   Episodes    Weighted     Duration     TS Min      TF Max  Excl\n");

    sn1 = SnTran[1] - 1;
    nna = nn = k = 0;
    org1 = -999;
    suma = sum = 0.0;

    for (i = 0; i < NTran; ++i) {

        sn = SnTran[i];
        org = OrgTran[i];
        if (sn != sn1 || org != org1) {
            if (nn > 0) {
                printf1("Sum         %9d %12.2lf\n",nn,sum);                   
                nn = 0;
                sum = 0.0;
            }
            prnchar('-',76,1);
            sn1 = sn;
            org1 = org;
        }
        printf1("%2d %4d %3d %9d %12.2lf %11.2lf %11.2lf %11.2lf    ",
            sn,org,DesTran[i],TranNE[i],(double)TranWE[i],
            (double)TranMD[i],(double)TranTSM[i],(double)TranTFM[i]);

        if (TranPtr[(sn - 1) * MaxOrg1 + org] >= 0)
            printf1("-\n");
        else {
            printf1("*\n");
            k++;
        }
        nn += TranNE[i];
        nna += TranNE[i];
        sum += (double)TranWE[i];
        suma += (double)TranWE[i];
    }
    printf1("Sum         %9d %12.2lf\n",nn,sum);                   
    if (nna > nn) {
        prnchar('-',76,1);
        printf1("Sum         %9d %12.2lf\n",nna,suma);                   
    }
    if (k) {
        printf1("\nWarning: episodes with * in column Excl do not belong to any transition\n");
        printf1("and will be excluded from all estimation procedures.\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  get_edat(i,sn,org,des,ts,tf)                                            */
/*                                                                          */
/*  get episode data for case i. Return 0 if OK, -1 if error.               */  
   
int get_edat(int i,int *sn,int *org,int *des,double *ts,double *tf)
{
    int r;  
    double x;

    if (TSVar == -1)           /* get starting time */
        *ts = TSESVal[0];
    else if (TSVar >= 0)
        *ts = get_data(TSVar,i);
    else {
        r = v_eval1(i,TSESCnt,TSESTyp,TSESVal,ESIdx,ts,0,0,0,0,0);
        if (r)  
            return(-1);
    }
    if (TFVar == -1)           /* get ending time */
        *tf = TFESVal[0];
    else if (TFVar >= 0)
        *tf = get_data(TFVar,i);
    else {
        r = v_eval1(i,TFESCnt,TFESTyp,TFESVal,ESIdx,tf,0,0,0,0,0);
        if (r)  
            return(-1);
    }
    if (ORGVar == -1)           /* get origin state */
        *org = (int)ORGESVal[0];
    else if (ORGVar >= 0)
        *org = (int)get_data(ORGVar,i);
    else {
        r = v_eval1(i,ORGESCnt,ORGESTyp,ORGESVal,ESIdx,&x,0,0,0,0,0);
        if (r)  
            return(-1);
        *org = (int)x;
    }
    if (DESVar == -1)           /* get destination state */
        *des = (int)DESESVal[0];
    else if (DESVar >= 0)
        *des = (int)get_data(DESVar,i);
    else {
        r = v_eval1(i,DESESCnt,DESESTyp,DESESVal,ESIdx,&x,0,0,0,0,0);
        if (r)  
            return(-1);
        *des = (int)x;
    }
    if (MEFlg) {      /* multiepisode data */

        if (SNVar == -1)           /* get spell number */
            *sn = (int)SNESVal[0];
        else if (SNVar >= 0)
            *sn = (int)get_data(SNVar,i);
        else {
            r = v_eval1(i,SNESCnt,SNESTyp,SNESVal,ESIdx,&x,0,0,0,0,0);
            if (r)  
                return(-1);
            *sn = (int)x;
        }
    }
    else
        *sn = 1;
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  get_id(icase)   get ID for icase.                                       */
/*                  only with multi-episode data. Otherwise id = icase + 1  */

int get_id(int icase,int *err)
{
    int r,id;
    double x;

    *err = 0;
    id = icase + 1;

    if (MEFlg) {
        if (IDVar == -1)           /* get spell number */
            id = (int)IDESVal[0];
        else if (IDVar >= 0)
            id = (int)get_data(IDVar,icase);
        else {
            r = v_eval1(icase,IDESCnt,IDESTyp,IDESVal,ESIdx,&x,0,0,0,0,0);
            if (r)  
                *err = r;
            else             
                id = (int)x;
        }
    }
    return(id);
}

/*--------------------------------------------------------------------------*/
/*  sort_ed(opt)    if opt == 1 then                                        */
/*                  sort episode data wrt starting and ending times         */
/*                  if successful, then                                     */
/*                  TSIdx[i] is pointer to data sorted according to ts, and */
/*                  TFIdx[i] is pointer to data sorted according to tf,     */
/*                  if opt == 0 then free memory.                           */
/*                  Return 0 if OK, -1 if error.                            */

int sort_ed(int opt)  
{
    register int i;
    int err,org,des,sn;
    double ts,tf;

    err = 0;
    if (opt == 0)
        goto SEDFin;

    if (!(TSIdx = (int *)calloc(NOC,sizeof(int))))   
        goto SEDFin;
    TSSortA = NOC;
    memrq(NOC,sizeof(int));

    if (!(TFIdx = (int *)calloc(NOC,sizeof(int))))   
        goto SEDFin;
    TFSortA = NOC;
    memrq(NOC,sizeof(int));
      
    for (i = 0; i < NOC; ++i)
        TSIdx[i] = i;

    if (alloc_acu(NOC))
        goto SEDFin;
                                /* do not sort if ts = constant */

    if (TSVar == -1 || (TSVar >= 0 && VTyp[TSVar] == 2))  
        ;
    else {
        printf1("Sorting episodes according to starting times.\n");

        for (i = 0; i < NOC; ++i) {
            get_edat(i,&sn,&org,&des,&ts,&tf);
            AcU[i] = ts;
        }
        qsort((char *)TSIdx,NOC,sizeof(int),tscomp);
    }
    for (i = 0; i < NOC; ++i)
        TFIdx[i] = i;
                                /* do not sort if tf = constant */

    if (TFVar == -1 || (TFVar >= 0 && VTyp[TFVar] == 2))  
        ;
    else {
        printf1("Sorting episodes according to ending times.\n");

        if (alloc_aci(NOC))
            goto SEDFin;
        if (alloc_acj(NOC))
            goto SEDFin;

        for (i = 0; i < NOC; ++i) {
            get_edat(i,&sn,&org,&des,&ts,&tf);
            AcU[i] = tf;
            AcI[i] = org;
            AcJ[i] = des;
        }
        qsort((char *)TFIdx,NOC,sizeof(int),tfcomp);
    }
    ESortFlg = 1;
    return(0);

SEDFin:
    alloc_acu(0);
    alloc_acv(0);
    if (err)  
        p_err(-2,1);

    if (TSSortA) {
        free((char *)TSIdx);
        memrq(-TSSortA,sizeof(int));
        TSSortA = 0;
    }
    if (TFSortA) {
        free((char *)TFIdx);
        memrq(-TFSortA,sizeof(int));
        TFSortA = 0;
    }
    ESortFlg = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  tscomp  Compare function for quicksort. The compare is according to the */
/*          starting times TS.                                              */

int tscomp(const void *arg1,const void *arg2) 
{     
    double ts1,ts2;

    ts1 = AcU[*(int *)arg1];
    ts2 = AcU[*(int *)arg2];
    if (ts1 > ts2)
        return(1);
    else if (ts1 < ts2)
        return(-1);
    else  
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  tfcomp  Compare function for quicksort. Compare according to ending     */
/*          times TF. If more than two episodes have the same ending time,  */
/*          then censored episodes comes after uncensored episodes.         */

int tfcomp(const void *arg1,const void *arg2)
{     
    int org1,org2,des1,des2;
    double tf1,tf2;
   
    tf1 = AcU[*(int *)arg1];
    tf2 = AcU[*(int *)arg2];

    if (tf1 > tf2)
        return(1);
    else if (tf1 < tf2)
        return(-1);
    else {

        org1 = AcI[*(int *)arg1];
        org2 = AcI[*(int *)arg2];
        des1 = AcJ[*(int *)arg1];
        des2 = AcJ[*(int *)arg2];

        if (org1 != org2)       /* if O states are different ordering   */
            return(0);          /* is not important.                    */
        else {
            if (org1 != des1) { 
                if (org2 != des2)
                    return(0);
                else
                    return(-1);
            }
            else { 
                if (org2 != des2)
                    return(1);
                else
                    return(0);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_spell   Get the next spell, or split.                               */
/*                                                                          */

int get_spell(int init,int *icase,int *sn,int *org,int *des,
               double *ts, double *tf,int *spl,int *nspl)
{
    static int isplit;
    static int nsplit;
    static int ic;
    static int stat;
    static int org1;
    static int des1;
    static int sn1;
    double t,ttf,tl,ts1,tf1;
    register int j,k,m;

    if (init) {
        NSPLEvalErr = NSplits = 0;
        isplit = nsplit = 0;
        ic = -1;
        return(0);
    }
    if (++isplit > nsplit) {

        if (++ic >= NOC)         /* no more cases */
            return(0);
           
        j = 0;
        get_edat(ic,&sn1,&org1,&des1,&ts1,&tf1);

        SPTTV[++j] = ts1;
        tl = ttf = tf1;
        m = 0;
  
        while (m >= 0) {
            m = -1;
            for (k = 0; k < NSP; ++k) {
                t = get_data(SPLITVar[k],ic);
                if (t > SPTTV[j] && t < tl) {
                    m = SPLITVar[k];
                    tl = t;
                }
            }
            if (m >= 0) {
                SPTTV[++j] = tl;
                SPIDX[j] = m;
                tl = ttf;
            }
        }
        nsplit = j;
        isplit = 1;

        SPTTV[++j] = ttf;
    }
    *icase = ic;

    *sn = sn1;
    *ts = SPTTV[isplit];
    *tf = SPTTV[isplit + 1];
    if (isplit == 1)
        *org = stat = org1;
    else
        *org = stat;

    if (isplit == nsplit)
        *des = des1;
    else
        *des = stat;

    *spl  = isplit;
    *nspl = nsplit;

    /* evaluate type 5 variables */

    EDVALSn  = *sn;
    EDVALOrg = *org;
    EDVALDes = *des;
    EDVALTs  = *ts;
    EDVALTf  = *tf;
    EDVALTime = *tf;

    if (NVAR5 > 0) {   
        j = VIFirst;
        while (j >= 0) {

            if (VTyp[j] == 5) {

                k = v_eval1(ic,VESCnt[j],VESTyp[j],VESVal[j],ESIdx,&t,0,0,0,0,0);
                if (k) {
                    NSPLEvalErr++;
                    NPFlgs[5] += 1;
                    t = 0.0;
                }
                put_data(t,j,ic);
                /*** printf1("put t=%lg into j=%d ic=%d\n",t,j,ic); **/
            }
            j = VNxt[j];
        }
    }
    NSplits++;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_nedat()   check episode data.                                     */
/*                  return 0 if OK, -1 if error.                            */
 
int check_nedat(void)
{
    int nrec,icase,sn,org,des,spl,nspl,id,id1,idn,r;
    double ts,tf;

    get_spell(1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);

    id1 = INTMAX;
    idn = nrec = 0;
    while (get_spell(0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (MEFlg) {
            id = get_id(icase,&r);
            if (r) {
                printf1("\nError in evaluating ID variable in case %d\n",icase + 1);
                return(-1);
            }
            if (id != id1) {
                idn++;
                id1 = id;
            }
        }
        nrec++;
    }
    newline();
    if (MEFlg)
        printf1("Number of individuals: %d\n",idn);
    printf1("Number of episodes: %d\n",NOC);
    if (NSP > 0)
        printf1("Number of splits: %d\n",nrec);

    if (NVAR5 > 0) {
        printf1("Errors in evaluating type 5 variables: %d\n",NSPLEvalErr);
        if (NSPLEvalErr > 0)  
            return(-1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  epdat()     epdat command in CmdBuf: epdat(...) = fname                 */
/*              options:                                                    */
/*              fmt= print format for ts, tf; def. 6.2                      */
/*              v = varlist                                                 */
/*              noc = ... only first noc cases                              */
/*              dtda=fname                                                  */
/*              return 0 if OK, -1 if error.                                */
 
int epdat(void)
{
    register int i,j;
    int err,nrec,icase,sn,org,des,spl,nspl,id,r,nr;
    double ts,tf;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (EDAvail == 0) {
        p_err(-15,1);
        return(0);
    }
    if (parm(CmdBuf + 5,1,1))  
        goto EPDATFin;
       
    if (SILENTFlg < 2)
        printfe("Writing episode data to: %s\n",PMFdName);
    nr = nrec = 0;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(6,2);

    get_spell(1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    
    while (get_spell(0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        id = get_id(icase,&r);
        if (r)
            nr++;

        fprintf(PMFd,"%6d %3d %3d %3d %3d %3d ",id,sn,nspl,spl,org,des);
        fprintf(PMFd,PMFmtS,ts);
        fprintf(PMFd,PMFmtS,tf);

        for (i = 0; i < PMNV; ++i) {
            j = PMVIdx[i];
            fprintf(PMFd,VPFmtS[j],get_data(j,icase));
        }
        fprintf(PMFd,"\n");
        nrec++;
        if (PMNOCFlg && nrec >= PMNOC)
            break;
    }
    printf1("Episode data: %d records written to: %s\n",nrec,PMFdName);

    if (PMTDAFDef)
        edtda(PMFdName,nrec);

    if (nr > 0)  
        printf1("Warning: %d errors in evaluating ID variable (ignored).\n",nr);

    err = 0;

EPDATFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  edtda(fname,noc)                                                        */
 
void edtda(char *fname,int noc)
{
    register int i,j;

    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(PMTDAFd,"  noc = %d,\n",noc);
    fprintf(PMTDAFd,"  ID   [6.0] = c1, # id number\n");
    fprintf(PMTDAFd,"  SN   [3.0] = c2, # spell number\n");
    fprintf(PMTDAFd,"  NSPL [3.0] = c3, # number of splits\n");
    fprintf(PMTDAFd,"  SPL  [3.0] = c4, # split number\n");
    fprintf(PMTDAFd,"  ORG  [3.0] = c5, # origin state\n");
    fprintf(PMTDAFd,"  DES  [3.0] = c6, # destination state\n");
    fprintf(PMTDAFd,"  TS   [%d.%d] = c7, # starting time\n",PMFmt1,PMFmt2);
    fprintf(PMTDAFd,"  TF   [%d.%d] = c8, # ending time\n",PMFmt1,PMFmt2);

    for (i = 0; i < PMNV; ++i) {
        j = PMVIdx[i];
        fprintf(PMTDAFd,"  %s [%d.%d] = c%d,\n",VName[j],VPFmt1[j],VPFmt2[j],i + 9);
    }
    fprintf(PMTDAFd,");\n");
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/* ------------------------------------------------------------------------ */
/*  epsdat()    epsdat command in CmdBuf: epsdat(...) = fname               */
/*              repuired parameter: t= time points.                         */
/*              write distribution of states to output file.                */
/*              return 0 if OK, -1 if error.                                */
 
int epsdat(void)
{
    register int i,j,k;
    int err,m,n,sn,org,des,r,id,id1,s;
    double ts,tf,t;

    err = -1;
    if (check_cmd(1))
        return(-1);

    if (EDAvail == 0) {
        p_err(-15,1);
        return(0);
    }
    if (parm(CmdBuf + 6,1,1))  
        goto EPSDATFin;
       
    if (PMNTP == 0) {
        p_err(-17,1);
        goto EPSDATFin;
    }
    if (SILENTFlg < 2)
        printfe("Writing state distributions to: %s\n",PMFdName);
    printf1("Writing state distributions to: %s\nUsing ",PMFdName);
    if (MEFlg)
        printf1("multi-");
    else
        printf1("single ");
    printf1("episode data.\n");

    if (alloc_acn(MaxOrg + 1))
        goto EPSDATFin;

    if (alloc_aci(MaxOrg + 1))
        goto EPSDATFin;

    for (i = 0; i < NTran; ++i)  
        AcI[OrgTran[i]] = 1;
             
    fprintf(PMFd,"# State distributions.\n");
    fprintf(PMFd,"# Time    ");
    for (i = 0; i <= MaxOrg; ++i) {
        if (AcI[i])
            fprintf(PMFd," %6d",i);
    }
    fprintf(PMFd,"    Total  Missing\n");


    for (k = 0; k < PMNTP; ++k) {        /* for all time points */

        t = PMTP[k];
        for (j = 0; j <= MaxOrg; ++j)
            AcN[j] = 0;

        m = n = 0;
        if (MEFlg == 0) {
            for (i = 0; i < NOC; ++i) {

                if (get_edat(i,&sn,&org,&des,&ts,&tf)) {
                    p_err(-39,1);
                    goto EPSDATFin;
                }   
                if (ts <= t && t < tf) {
                    AcN[org] += 1;
                    n++;
                }
                m++;
            }
        }
        else {
             
            s = -1;
            id1 = INTMAX;
            for (i = 0; i < NOC; ++i) {

                id = get_id(i,&r);
                if (id != id1) {
                    if (s >= 0) {
                        AcN[s] += 1;
                        n++;
                        s = -1;
                    }  
                    m++;
                    id1 = id;
                }
                if (r || get_edat(i,&sn,&org,&des,&ts,&tf)) {
                    p_err(-39,1);
                    goto EPSDATFin;
                }   
                if (s < 0 && ts <= t && t < tf)  
                    s = org;
            }
            if (s >= 0) {
                AcN[s] += 1;
                n++;
            }
        }
        fprintf(PMFd,"%10.4lf",t);
        for (j = 0; j <= MaxOrg; ++j) {
            if (AcI[j])
                fprintf(PMFd," %6d",AcN[j]);
        }
        fprintf(PMFd," %8d %8d\n",n,m - n);
    }
    if (MEFlg)
        printf1("Number of individuals: %d\n",m);
    else
        printf1("Number of single episodes: %d\n",m);

    err = 0;

EPSDATFin:
    p_clean();
    return(err);
}


