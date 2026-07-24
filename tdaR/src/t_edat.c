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
#include "tda_context.h"
#include "tda_compat.h"

/*  functions in t_edat.c */

int check_edj(TDAContext *ctx, int j);
int edef(TDAContext *ctx);
void prn_eed(TDAContext *ctx, char *s);
void edat_off(TDAContext *ctx, int opt);
int check_nevar(TDAContext *ctx, int nidx);
int def_edat(TDAContext *ctx, int mode);
void prn_edef(TDAContext *ctx);
void edef_info(TDAContext *ctx);
int check_edat(TDAContext *ctx);
int alloc_edat(TDAContext *ctx, int mode);
void prn_edat(TDAContext *ctx);
int get_edat(TDAContext *ctx, int i,int *sn,int *org,int *des,double *ts,double *tf);
int get_id(TDAContext *ctx, int icase,int *err);
int sort_ed(TDAContext *ctx, int opt);
int tscomp(const void *, const void *, void *);
int tfcomp(const void *, const void *, void *);
int get_spell(TDAContext *ctx, int init,int *icase,int *sn,int *org,int *des, double *ts,double *tf,int *spl,int *nspl);
int check_nedat(TDAContext *ctx);
int epdat(TDAContext *ctx);
void edtda(TDAContext *ctx, char *fname,int noc);
int epsdat(TDAContext *ctx);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */




/*  The following variables have the following values:                      */
/*  >= 0, then it is the internal variable number                           */
/*    -1, then it is a numerical constant                                   */
/*    -2, it is a general expression                                        */








/*  We assume spell numbers must be not less than 1, origin and destination */
/*  state number must be not less than 0.                                   */






/* ------------------------------------------------------------------------ */



/* ------------------------------------------------------------------------ */
/*  check_edj(j)    Check whether variable j is needed for episode data.    */
/*                  Return 1 if needed, otherwise 0.                        */

int check_edj(TDAContext *ctx, int j)
{
    register int i;

    if (ctx->EDAvail == 0)
        return(0);
    if (j == ctx->IDVar || j == ctx->SNVar || j == ctx->TSVar || j == ctx->TFVar || j == ctx->ORGVar || j == ctx->DESVar)
        return(1);
    for (i = 0; i < ctx->NSP; ++i) {
        if (j == ctx->SPLITVar[i])
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  edef()      Setup episode data. Command in CmdBuf.                      */
/*              edef(org,des,ts,tf,id,sn,maxtran,vardef...).                */
/*              Return 0 if OK, -1 if error.                                */

int edef(TDAContext *ctx)
{
    register int i;
    int err,n,j,l,nidx,nb;
    char c,*p,*q;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (!strcmp(ctx->CmdBuf,"edef")) {
        if (ctx->EDAvail == 0)
            printf1(ctx, "No episode data structure defined.\n");
        else
            edef_info(ctx);
        return(0);
    }
    if (ctx->DMDef == 0) {
        p_err(ctx, -18,1);
        return(0);
    }
    if (ctx->EDAvail)
        edat_off(ctx, 1);

    ctx->MaxTran = MaxTranDef;       /* default max number of transitions */
    ctx->NSP = 0;                    /* number of split variables */
    nidx = get_nidx(ctx);          /* index to first new variable */
    ctx->NVAR5 = 0;                  /* number of type 5 variables */

    ctx->PORGPtr = NULL;
    ctx->PDESPtr = NULL;
    ctx->PTSPtr  = NULL;
    ctx->PTFPtr  = NULL;
    ctx->PIDPtr  = NULL;
    ctx->PSNPtr  = NULL;

    p = ctx->CmdBuf + 4;
    while (*++p) {
        if (sscanf(p,"maxtran=%d",&n) == 1 && n > 0) {
            ctx->MaxTran = n;
            p = skip_int(ctx, p + 8);
        }
        else if (!strncmp(p,"org=",4) && *(p + 4)) {
            ctx->PORGPtr = p + 4;
            q = skip_expr(ctx, ctx->PORGPtr);
            strncpy(ctx->PORGStr,ctx->PORGPtr,(size_t)((int)(q - ctx->PORGPtr)));
            p = skip_nc(ctx, p + 4);
        }
        else if (!strncmp(p,"des=",4) && *(p + 4)) {
            ctx->PDESPtr = p + 4;
            q = skip_expr(ctx, ctx->PDESPtr);
            strncpy(ctx->PDESStr,ctx->PDESPtr,(size_t)((int)(q - ctx->PDESPtr)));
            p = skip_nc(ctx, p + 4);
        }
        else if (!strncmp(p,"ts=",3) && *(p + 3)) {
            ctx->PTSPtr = p + 3;
            q = skip_expr(ctx, ctx->PTSPtr);
            strncpy(ctx->PTSStr,ctx->PTSPtr,(size_t)((int)(q - ctx->PTSPtr)));
            p = skip_nc(ctx, p + 3);
        }
        else if (!strncmp(p,"tf=",3) && *(p + 3)) {
            ctx->PTFPtr = p + 3;
            q = skip_expr(ctx, ctx->PTFPtr);
            strncpy(ctx->PTFStr,ctx->PTFPtr,(size_t)((int)(q - ctx->PTFPtr)));
            p = skip_nc(ctx, p + 3);
        }
        else if (!strncmp(p,"id=",3) && *(p + 3)) {
            ctx->PIDPtr = p + 3;
            q = skip_expr(ctx, ctx->PIDPtr);
            strncpy(ctx->PIDStr,ctx->PIDPtr,(size_t)((int)(q - ctx->PIDPtr)));
            p = skip_nc(ctx, p + 3);
        }
        else if (!strncmp(p,"sn=",3) && *(p + 3)) {
            ctx->PSNPtr = p + 3;
            q = skip_expr(ctx, ctx->PSNPtr);
            strncpy(ctx->PSNStr,ctx->PSNPtr,(size_t)((int)(q - ctx->PSNPtr)));
            p = skip_nc(ctx, p + 3);
        }
        else if (!strncmp(p,"split=",6)) {
            q = p;
            p = get_nvia(ctx, p + 6,&n,1,&nb);
            if (n <= 0 || nb) {
                printf1(ctx, "Check: %s\n",q);
                goto EDEFFin;
            }
            if (n > MaxSPL) {
                printf1(ctx, "Error: exceeded max number of split variables.\n");
                goto EDEFFin;
            }
            j = 0;
            for (i = 0; i < n; ++i) {
                ctx->SPLITVar[i] = ctx->VLVIdx[i];
                if (ctx->VTyp[ctx->VLVIdx[i]] == 5)
                    j = 1;
            }
            if (j) {
                printf1(ctx, "Error: cannot use type 5 variables for episode splitting.\n");
                goto EDEFFin;
            }
            ctx->NSP = n;
        }
        else if ((l = get_vnlen(ctx, p)) > 0) {       /* variable */
            q = skip_nc(ctx, p + l);
            c = *q;
            *q = '\0';
            if (save_var(ctx, p,1))
                goto EDEFFin;
            ctx->NVAR5++;
            *q = c;
            p = q;
        }

        if (*p != ',' && *p != ')') {
            prn_eed(ctx, p);
            goto EDEFFin;
        }
        *p = '\0';
    }

    if (ctx->PORGPtr == NULL || ctx->PDESPtr == NULL || ctx->PTSPtr == NULL || ctx->PTFPtr == 0) {
        printf1(ctx, "Error: need org, des, ts, and tf parameters.\n");
        goto EDEFFin;
    }
    if (ctx->PIDPtr != NULL && ctx->PSNPtr != 0)
        ctx->MEFlg = 1;

    else if (ctx->PIDPtr != NULL || ctx->PSNPtr != 0) {
        printf1(ctx, "Error: need id and sn parameters for multi-episode data.\n");
        goto EDEFFin;
    }
    else
        ctx->MEFlg = 0;

    err = def_edat(ctx, 1);              /* check definition */
    if (err)
        goto EDEFFin;
    prn_edef(ctx);                     /* print definition */

    if (ctx->NSP > 0) {
        printf1(ctx, "Episode splitting with variable(s): %s",ctx->VName[ctx->SPLITVar[0]]);
        for (j = 1; j < ctx->NSP; ++j)
            printf1(ctx, ",%s",ctx->VName[ctx->SPLITVar[j]]);
        printf1(ctx, "\n");
    }
    newline(ctx);
    err = check_edat(ctx);             /* check data */
    if (err)
        goto EDEFFin;
    prn_edat(ctx);                     /* print info */

    if (ctx->NVAR5 > 0) {                /* check new variables */
        err = check_nevar(ctx, nidx);
        if (err)
            goto EDEFFin;
    }
    err = check_nedat(ctx);            /* check correct data generation */
    if (err)
        goto EDEFFin;

    ctx->EDAvail = 1;

EDEFFin:
    p_clean(ctx);
    if (err) {
        edat_off(ctx, 0);
        printf1(ctx, "Could not create new episode data.\n");
    }
    else
        printf1(ctx, "Successfully created new episode data.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_eed.    Print error message.                                        */

void prn_eed(TDAContext *ctx, char *s)
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
/*  edat_off(opt)  free all data structures for episode data.               */
/*                 if opt != 0 print message.                               */

void edat_off(TDAContext *ctx, int opt)
{
    register int j;

    if (ctx->ESortFlg)
        sort_ed(ctx, 0);
    alloc_edat(ctx, 0);
    def_edat(ctx, 0);

    j = ctx->VIFirst;        /* remove type 5 variables */
    while (j >= 0) {
        if (ctx->VTyp[j] == 5) {
            clear_var(ctx, j);
            j = ctx->VIFirst;
        }
        else
            j = ctx->VNxt[j];
    }
    if (opt && ctx->EDAvail)
        printf1(ctx, "Deallocated episode data structures.\n");
    ctx->EDAvail = 0;
}

/* ------------------------------------------------------------------------ */
/*  check_nevar()   check new type 5 variables, begin with index vidx.      */
/*                  return 0 if OK, -1 if error.                            */

int check_nevar(TDAContext *ctx, int nidx)
{
    int err;

    err = -1;
    printf1(ctx, "\nType 5 variables (dependent on current episode data).\n\n");

    prn_var(ctx, nidx);          /* print variables */

    if (alloc_vdat(ctx, nidx,1)) {
        printf1(ctx, "\nInsufficient memory for new variables.\n");
        goto NEVARFin;
    }
    err = 0;

NEVARFin:
    if (err) {
        clear_avar(ctx, nidx);
        ctx->NVAR5 = 0;
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

int def_edat(TDAContext *ctx, int mode)
{
    register int i;
    int n,err;

    err = 0;
    if (mode == 0) {
        ctx->NSP = ctx->MEFlg = 0;
        ctx->IDVar = ctx->SNVar = ctx->TSVar = ctx->TFVar = ctx->ORGVar = ctx->DESVar = -2;
        goto EDATFin;
    }
    err = -1;

    /* save org */

    if ((n = v_parse(ctx, ctx->PORGPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PORGPtr);
        if (n < 0)
            prn_emsg1(ctx, n);
        goto EDATFin;
    }
    if (!(ctx->ORGESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto EDATFin;
    }
    if (!(ctx->ORGESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        p_err(ctx, -2,1);
        free((char *)ctx->ORGESTyp);
        goto EDATFin;
    }
    memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
    ctx->ORGESCnt = ctx->ESCnt;

    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->ORGESTyp[i] = ctx->ESTyp[i];
        ctx->ORGESVal[i] = ctx->ESVal[i];
    }
    if (ctx->ESCnt == 1) {
        if (ctx->ESTyp[0] == 0)
            ctx->ORGVar = -1;
        else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
            ctx->ORGVar = ctx->ESTyp[0] - ctx->VOFFS;
    }

    /* save des */

    if ((n = v_parse(ctx, ctx->PDESPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PDESPtr);
        if (n < 0)
            prn_emsg1(ctx, n);
        goto EDATFin;
    }
    if (!(ctx->DESESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto EDATFin;
    }
    if (!(ctx->DESESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        p_err(ctx, -2,1);
        free((char *)ctx->DESESTyp);
        goto EDATFin;
    }
    memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
    ctx->DESESCnt = ctx->ESCnt;

    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->DESESTyp[i] = ctx->ESTyp[i];
        ctx->DESESVal[i] = ctx->ESVal[i];
    }
    if (ctx->ESCnt == 1) {
        if (ctx->ESTyp[0] == 0)
            ctx->DESVar = -1;
        else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
            ctx->DESVar = ctx->ESTyp[0] - ctx->VOFFS;
    }

    /* save ts */

    if ((n = v_parse(ctx, ctx->PTSPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PTSPtr);
        if (n < 0)
            prn_emsg1(ctx, n);
        goto EDATFin;
    }
    if (!(ctx->TSESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto EDATFin;
    }
    if (!(ctx->TSESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        p_err(ctx, -2,1);
        free((char *)ctx->TSESTyp);
        goto EDATFin;
    }
    memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
    ctx->TSESCnt = ctx->ESCnt;

    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->TSESTyp[i] = ctx->ESTyp[i];
        ctx->TSESVal[i] = ctx->ESVal[i];
    }
    if (ctx->ESCnt == 1) {
        if (ctx->ESTyp[0] == 0)
            ctx->TSVar = -1;
        else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
            ctx->TSVar = ctx->ESTyp[0] - ctx->VOFFS;
    }

    /* save tf */

    if ((n = v_parse(ctx, ctx->PTFPtr,0)) < 0 || ctx->ESCnt <= 0) {
        printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PTFPtr);
        if (n < 0)
            prn_emsg1(ctx, n);
        goto EDATFin;
    }
    if (!(ctx->TFESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto EDATFin;
    }
    if (!(ctx->TFESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
        p_err(ctx, -2,1);
        free((char *)ctx->TFESTyp);
        goto EDATFin;
    }
    memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
    ctx->TFESCnt = ctx->ESCnt;

    for (i = 0; i < ctx->ESCnt; ++i) {
        ctx->TFESTyp[i] = ctx->ESTyp[i];
        ctx->TFESVal[i] = ctx->ESVal[i];
    }
    if (ctx->ESCnt == 1) {
        if (ctx->ESTyp[0] == 0)
            ctx->TFVar = -1;
        else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
            ctx->TFVar = ctx->ESTyp[0] - ctx->VOFFS;
    }

    /* save id */

    if (ctx->PIDPtr != NULL) {
        if ((n = v_parse(ctx, ctx->PIDPtr,0)) < 0 || ctx->ESCnt <= 0) {
            printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PIDPtr);
            if (n < 0)
                prn_emsg1(ctx, n);
            goto EDATFin;
        }
        if (!(ctx->IDESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto EDATFin;
        }
        if (!(ctx->IDESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
            p_err(ctx, -2,1);
            free((char *)ctx->IDESTyp);
            goto EDATFin;
        }
        memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
        ctx->IDESCnt = ctx->ESCnt;

        for (i = 0; i < ctx->ESCnt; ++i) {
            ctx->IDESTyp[i] = ctx->ESTyp[i];
            ctx->IDESVal[i] = ctx->ESVal[i];
        }
        if (ctx->ESCnt == 1) {
            if (ctx->ESTyp[0] == 0)
                ctx->IDVar = -1;
            else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
                ctx->IDVar = ctx->ESTyp[0] - ctx->VOFFS;
        }
    }

    /* save sn */

    if (ctx->PSNPtr != NULL) {
        if ((n = v_parse(ctx, ctx->PSNPtr,0)) < 0 || ctx->ESCnt <= 0) {
            printf1(ctx, "Syntax or reference error (%d): %s\n",n,ctx->PSNPtr);
            if (n < 0)
                prn_emsg1(ctx, n);
            goto EDATFin;
        }
        if (!(ctx->SNESTyp = (int *)calloc((size_t)(ctx->ESCnt),sizeof(int)))) {
            p_err(ctx, -2,1);
            goto EDATFin;
        }
        if (!(ctx->SNESVal = (double *)calloc((size_t)(ctx->ESCnt),sizeof(double)))) {
            p_err(ctx, -2,1);
            free((char *)ctx->SNESTyp);
            goto EDATFin;
        }
        memrq(ctx, ctx->ESCnt,sizeof(int) + sizeof(double));
        ctx->SNESCnt = ctx->ESCnt;

        for (i = 0; i < ctx->ESCnt; ++i) {
            ctx->SNESTyp[i] = ctx->ESTyp[i];
            ctx->SNESVal[i] = ctx->ESVal[i];
        }
        if (ctx->ESCnt == 1) {
            if (ctx->ESTyp[0] == 0)
                ctx->SNVar = -1;
            else if (ctx->ESTyp[0] >= ctx->VOFFS && ctx->ESTyp[0] < ctx->COFFS)
                ctx->SNVar = ctx->ESTyp[0] - ctx->VOFFS;
        }
    }
    return(0);

EDATFin:

    if (ctx->IDESCnt > 0) {
        free((char *)ctx->IDESTyp);
        free((char *)ctx->IDESVal);
        memrq(ctx, -ctx->IDESCnt,sizeof(int) + sizeof(double));
        ctx->IDESCnt = 0;
    }
    if (ctx->SNESCnt > 0) {
        free((char *)ctx->SNESTyp);
        free((char *)ctx->SNESVal);
        memrq(ctx, -ctx->SNESCnt,sizeof(int) + sizeof(double));
        ctx->SNESCnt = 0;
    }
    if (ctx->TSESCnt > 0) {
        free((char *)ctx->TSESTyp);
        free((char *)ctx->TSESVal);
        memrq(ctx, -ctx->TSESCnt,sizeof(int) + sizeof(double));
        ctx->TSESCnt = 0;
    }
    if (ctx->TFESCnt > 0) {
        free((char *)ctx->TFESTyp);
        free((char *)ctx->TFESVal);
        memrq(ctx, -ctx->TFESCnt,sizeof(int) + sizeof(double));
        ctx->TFESCnt = 0;
    }
    if (ctx->ORGESCnt > 0) {
        free((char *)ctx->ORGESTyp);
        free((char *)ctx->ORGESVal);
        memrq(ctx, -ctx->ORGESCnt,sizeof(int) + sizeof(double));
        ctx->ORGESCnt = 0;
    }
    if (ctx->DESESCnt > 0) {
        free((char *)ctx->DESESTyp);
        free((char *)ctx->DESESVal);
        memrq(ctx, -ctx->DESESCnt,sizeof(int) + sizeof(double));
        ctx->DESESCnt = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_edef()  print definition of episode data.                           */

void prn_edef(TDAContext *ctx)
{
    printf1(ctx, "Creating new ");
    if (ctx->MEFlg)
        printf1(ctx, "multi-");
    else
        printf1(ctx, "single ");
    printf1(ctx, "episode data. Max number of transitions: %d.\n",ctx->MaxTran);
    printf1(ctx, "Definition: ");
    if (ctx->MEFlg)
        printf1(ctx, "id=%s, sn=%s, ",ctx->PIDPtr,ctx->PSNPtr);
    printf1(ctx, "org=%s, des=%s, ts=%s, tf=%s\n",ctx->PORGPtr,ctx->PDESPtr,ctx->PTSPtr,ctx->PTFPtr);
}

/* ------------------------------------------------------------------------ */
/*  edef_info()  info about definition of episode data.                     */

void edef_info(TDAContext *ctx)
{
    if (ctx->MEFlg)
        printf1(ctx, "Multi-");
    else
        printf1(ctx, "Single ");
    printf1(ctx, "episode data. Max number of transitions: %d.\n",ctx->MaxTran);
    printf1(ctx, "Definition: ");
    if (ctx->MEFlg)
        printf1(ctx, "id=%s, sn=%s, ",ctx->PIDStr,ctx->PSNStr);
    printf1(ctx, "org=%s, des=%s, ts=%s, tf=%s\n\n",ctx->PORGStr,ctx->PDESStr,ctx->PTSStr,ctx->PTFStr);
    prn_edat(ctx);
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

int check_edat(TDAContext *ctx)
{
    register int i,j,k;
    int err,r,wflag,sn,sn1,org,org1,des,*idat,*tran,*freq,*bf;
    int idata,freqa,trana,bfa,wta,wfreqa;
    double x,ts,tf,*wt = NULL,*wfreq = NULL;
    double w,d;

    ctx->EDAvail = 0;
    wflag = idata = freqa = trana = bfa = wta = wfreqa = 0;
    err = -2;

    /* allocate memory for frequency distribution of transitions */

    if (!(idat = (int *)calloc((size_t)(3) * (size_t)(ctx->NOC) + 1,sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CDFin;
    }
    idata = 3 * ctx->NOC + 1;
    memrq(ctx, idata,sizeof(int));

    if (!(freq = (int *)calloc((size_t)(ctx->MaxTran + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CDFin;
    }
    freqa = ctx->MaxTran + 1;
    memrq(ctx, freqa,sizeof(int));

    if (!(tran = (int *)calloc((size_t)(3) * (size_t)(ctx->MaxTran) + 1,sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CDFin;
    }
    trana = 3 * ctx->MaxTran + 1;
    memrq(ctx, trana,sizeof(int));

    if (!(bf = (int *)calloc((size_t)(ctx->MaxTran + 1),sizeof(int)))) {
        p_err(ctx, -2,1);
        goto CDFin;
    }
    bfa = ctx->MaxTran + 1;
    memrq(ctx, bfa,sizeof(int));

    if (ctx->WIVar >= 0) {       /* if case weights */

        wflag = 1;

        if (!(wt = (double *)calloc((size_t)(ctx->NOC + 1),sizeof(double)))) {
            p_err(ctx, -2,1);
            goto CDFin;
        }
        wta = ctx->NOC + 1;
        memrq(ctx, wta,sizeof(double));

        if (!(wfreq = (double *)calloc((size_t)(ctx->MaxTran + 1),sizeof(double)))) {
            p_err(ctx, -2,1);
            goto CDFin;
        }
        wfreqa = ctx->MaxTran + 1;
        memrq(ctx, wfreqa,sizeof(double));
    }
    err = 0;

    for (i = 0; i < ctx->NOC; ++i) {

        if (ctx->TSVar == -1)           /* get starting time */
            ts = ctx->TSESVal[0];
        else if (ctx->TSVar >= 0)
            ts = get_data(ctx, ctx->TSVar,i);
        else {
            r = v_eval1(ctx, i,ctx->TSESCnt,ctx->TSESTyp,ctx->TSESVal,ctx->ESIdx,&ts,0,0,0,0,0);
            if (r) {
                printf1(ctx, "Error (%d): can't evaluate starting time in case %d.\n",r,i+1);
                prn_emsg2(ctx, r);
                err = -1;
                break;
            }
        }
        if (ctx->TFVar == -1)           /* get ending time */
            tf = ctx->TFESVal[0];
        else if (ctx->TFVar >= 0)
            tf = get_data(ctx, ctx->TFVar,i);
        else {
            r = v_eval1(ctx, i,ctx->TFESCnt,ctx->TFESTyp,ctx->TFESVal,ctx->ESIdx,&tf,0,0,0,0,0);
            if (r) {
                printf1(ctx, "Error (%d): can't evaluate ending time in case %d.\n",r,i+1);
                prn_emsg2(ctx, r);
                err = -1;
                break;
            }
        }

        /* check starting and ending time */

        if (ts < 0.0) {
            printf1(ctx, "Error: found negative starting time (%lg) in case %d.\n",ts,i+1);
            err = -3;
            break;
        }
        if (tf <= ts) {
            printf1(ctx, "Error: found zero or negative duration (%lg) in case %d.\n",tf - ts,i+1);
            err = -4;
            break;
        }

        if (ctx->ORGVar == -1)           /* get origin state */
            org = (int)ctx->ORGESVal[0];
        else if (ctx->ORGVar >= 0)
            org = (int)get_data(ctx, ctx->ORGVar,i);
        else {
            r = v_eval1(ctx, i,ctx->ORGESCnt,ctx->ORGESTyp,ctx->ORGESVal,ctx->ESIdx,&x,0,0,0,0,0);
            if (r) {
                printf1(ctx, "Error (%d): can't evaluate origin state in case %d.\n",r,i+1);
                prn_emsg2(ctx, r);
                err = -1;
                break;
            }
            org = (int)x;
        }

        if (ctx->DESVar == -1)           /* get destination state */
            des = (int)ctx->DESESVal[0];
        else if (ctx->DESVar >= 0)
            des = (int)get_data(ctx, ctx->DESVar,i);
        else {
            r = v_eval1(ctx, i,ctx->DESESCnt,ctx->DESESTyp,ctx->DESESVal,ctx->ESIdx,&x,0,0,0,0,0);
            if (r) {
                printf1(ctx, "Error (%d): can't evaluate destination state in case %d.\n",r,i+1);
                prn_emsg2(ctx, r);
                err = -1;
                break;
            }
            des = (int)x;
        }
        if (org < 0 || des < 0) {
            printf1(ctx, "Error: found negative origin (%d) or destination (%d) state in case %d.\n",org,des,i+1);
            err = -5;
            break;
        }
        idat[i * 3 + 1] = 1;
        idat[i * 3 + 2] = org;
        idat[i * 3 + 3] = des;

        if (ctx->MEFlg) {      /* multiepisode data */

            if (ctx->SNVar == -1)           /* get spell number */
                sn = (int)ctx->SNESVal[0];
            else if (ctx->SNVar >= 0)
                sn = (int)get_data(ctx, ctx->SNVar,i);
            else {
                r = v_eval1(ctx, i,ctx->SNESCnt,ctx->SNESTyp,ctx->SNESVal,ctx->ESIdx,&x,0,0,0,0,0);
                if (r) {
                    printf1(ctx, "Error (%d): can't evaluate destination state in case %d.\n",r,i+1);
                    prn_emsg2(ctx, r);
                    err = -1;
                    break;
                }
                sn = (int)x;
            }
            if (sn < 1) {
                printf1(ctx, "Error: found spell number (%d) less than 1 in case %d.\n",sn,i+1);
                err = -6;
                break;
            }
            idat[i * 3 + 1] = sn;
        }
        if (ctx->WIVar >= 0)                             /* get weights */
            wt[i + 1] = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
    }
    if (err)
        goto CDFin;

    /* make joint frequency table for spell number, origin and destination states */

    r = cfreq(ctx, 3,ctx->NOC,idat,ctx->MaxTran,tran,freq,bf,wflag,wt,wfreq);
    if (r <= 0) {
        printf1(ctx, "Error in checking transitions.\n");
        err = -9;
        if (r == -1) {
            err = -8;
            printf1(ctx, "Number of transitions exceeds the defined maximum (%d).\n",ctx->MaxTran);
        }
        else if (r == -2)
            printf1(ctx, "Stack size is too small.\n");
        else if (r == -3)
            printf1(ctx, "Insufficient memory for frequency distribution.\n");
        else
            printf1(ctx, "Error %d in frequency distribution.\n",r);
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

    ctx->NTran = r;
    ctx->MaxSnn = ctx->MaxOrg = ctx->MaxDes = ctx->NTran1 = 0;

    for (i = 1; i <= ctx->NTran; ++i) {

        j = bf[i];

        sn  = tran[(j - 1) * 3 + 1];
        org = tran[(j - 1) * 3 + 2];
        des = tran[(j - 1) * 3 + 3];

        if (des != org)
            ctx->NTran1++;
        if (ctx->MaxSnn < sn)
            ctx->MaxSnn = sn;
        if (ctx->MaxOrg < org)
            ctx->MaxOrg = org;
        if (ctx->MaxDes < des)
            ctx->MaxDes = des;
    }
    if (ctx->NTran1 == 0) {
        printf1(ctx, "Number of transitions is zero.\n");
        err = -10;
        goto CDFin;
    }
    ctx->MaxOrg1 = ctx->MaxOrg + 1;
    ctx->MaxDes1 = ctx->MaxDes + 1;

    if (alloc_edat(ctx, 1)) {    /* allocate additional data structures */
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

    for (i = 0; i < ctx->NTran; ++i) {

        j = bf[i + 1];
        ctx->SnTran[i]  = tran[(j - 1) * 3 + 1];
        ctx->OrgTran[i] = tran[(j - 1) * 3 + 2];
        ctx->DesTran[i] = tran[(j - 1) * 3 + 3];
        ctx->TranNE[i]  = freq[j];

        if (ctx->WIVar >= 0)
            ctx->TranWE[i] = (float)wfreq[j];
        else
            ctx->TranWE[i] = (float)freq[j];
    }

    /* -------------------------------------------------------------------- */
    /*  create parallel arrays, but only for transitions, in                */
    /*  SnTran1[i]   = spell number                                         */
    /*  OrgTran1[i]  = origin state                                         */
    /*  DesTran1[i]  = destination state                                    */

    k = 0;
    for (i = 0; i < ctx->NTran; ++i) {

        sn  = ctx->SnTran[i];
        org = ctx->OrgTran[i];
        des = ctx->DesTran[i];

        if (des != org) {
            ctx->SnTran1[k] = sn;
            ctx->OrgTran1[k] = org;
            ctx->DesTran1[k] = des;
            k++;
        }
    }
    ctx->SnTran1[k] = ctx->OrgTran1[k] = ctx->DesTran1[k] = -1;


    /********************
    printf1(ctx, "\nTran1...\n");
    for (i = 0; i < NTran1; ++i)
        printf1(ctx, "%2d %2d %2d \n",SnTran1[i],OrgTran1[i],DesTran1[i]);
    ***************/
    /* -------------------------------------------------------------------- */
    /*  create pointer to SnTran1,OrgTran1,DesTran1 in TranPtr.             */
    /*  TranPtr[sn,org] = TranPtr[(sn - 1) * MaxOrg1 + org]                 */
    /*  sn = 1,...,MaxSnn; org = 0,...,MaxOrg. (MaxOrg1 = MaxOrg + 1)       */

    j = ctx->MaxSnn * ctx->MaxOrg1 + 1;
    for (i = 0; i < j; ++i)
        ctx->TranPtr[i] = -1;

    sn = ctx->SnTran1[0] - 1;
    org = ctx->OrgTran1[0] - 1;

    for (i = 0; i < ctx->NTran1; ++i) {

        sn1 = ctx->SnTran1[i];
        org1 = ctx->OrgTran1[i];

        if (sn1 != sn || org1 != org) {
            ctx->TranPtr[(sn1 - 1) * ctx->MaxOrg1 + org1] = i;
            sn = sn1;
            org = org1;
        }
    }
    /***********************
    printf1(ctx, "\nTranPtr\n");
    for (i = 1; i <= MaxSnn; ++i) {
        for (j = 0; j <= MaxOrg; ++j)
            printf1(ctx, "%2d ",TranPtr[(i - 1) * MaxOrg1 + j]);
        printf1(ctx, "\n");
    }
    *********/
    /* -------------------------------------------------------------------- */
    /*  TranMD[i]  = weighted mean duration for transition i                */
    /*  TranTSM[i] = minimum of starting times for transition i             */
    /*  TranTFM[i] = maximum of ending times for transition i               */
    /*  i = 0,...,NTran-1.                                                  */

    x = ctx->DBLMAX;
    for (j = 0; j < ctx->NTran; ++j) {
        ctx->TranTSM[j] = (float)x;
        ctx->TranTFM[j] = -1.0;
    }

    for (i = 0; i < ctx->NOC; ++i) {

        if (get_edat(ctx, i,&sn,&org,&des,&ts,&tf)) {
            p_err(ctx, -39,1);
            err = -7;
            goto CDFin;
        }
        if (ctx->WIVar >= 0)                             /* get weights */
            x = get_data(ctx, ctx->WIVar,i) * ctx->WNorm;
        else
            x = 1.0;

        for (j = 0; j < ctx->NTran; ++j) {
            if (sn == ctx->SnTran[j] && org == ctx->OrgTran[j] && des == ctx->DesTran[j]) {
                ctx->TranMD[j] += (float)(x * (tf - ts));
                if (ctx->TranTSM[j] > (float)ts)
                    ctx->TranTSM[j] = (float)ts;
                if (ctx->TranTFM[j] < (float)tf)
                    ctx->TranTFM[j] = (float)tf;
                break;
            }
        }
    }

    for (j = 0; j < ctx->NTran; ++j)
        ctx->TranMD[j] /= ctx->TranWE[j];

    /* -------------------------------------------------------------------- */
    /*  TWFreq[i]  = weighted number of uncensored episodes                 */
    /*               for transition i = 0,...,NTran1 - 1.                   */
    /*  TWDur[i]   = sum of weighted duration for uncensored episodes       */
    /*               in transition i.                                       */

    for (i = 0; i < ctx->NTran1; ++i) {

        sn = ctx->SnTran1[i];
        org = ctx->OrgTran1[i];
        des = ctx->DesTran1[i];
        w = d = 0.0;
        k = 0;
        for (j = 0; j < ctx->NTran; ++j) {
            if (sn == ctx->SnTran[j] && org == ctx->OrgTran[j]) {
                k++;
                d += (double)((ctx->TranMD[j])) * (double)(ctx->TranWE[j]);
                if (des == ctx->DesTran[j])
                    w += (double)(ctx->TranWE[j]);
            }
            else if (k)
                break;
        }
        ctx->TWDur[i] = (float)(d);
        ctx->TWFreq[i] = (float)(w);
    }
    /***********
    printf1(ctx, "TWFreq TWDur\n");
    for (i = 0; i < NTran1; ++i)
        printf1(ctx, "%lf %lf\n",TWFreq[i],TWDur[i]);
    printf1(ctx, "\n");
    *********/

CDFin:
    if (wfreqa) {
        free((char *)wfreq);
        memrq(ctx, -wfreqa,sizeof(double));
    }
    if (wta) {
        free((char *)wt);
        memrq(ctx, -wta,sizeof(double));
    }
    if (bfa) {
        free((char *)bf);
        memrq(ctx, -bfa,sizeof(int));
    }
    if (trana) {
        free((char *)tran);
        memrq(ctx, -trana,sizeof(int));
    }
    if (freqa) {
        free((char *)freq);
        memrq(ctx, -freqa,sizeof(int));
    }
    if (idata) {
        free((char *)idat);
        memrq(ctx, -idata,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  alloc_edat(mode)    If mode != 0 allocate memory for episode data       */
/*                      structures, else free previously allocated memory.  */
/*                                                                          */
/*                  Return 0 if OK, -1 if error.                            */

int alloc_edat(TDAContext *ctx, int mode)
{
    int err = 0;

    if (mode == 0)
        goto AADFin;

    err = -1;
    if (!(ctx->SnTran = (int *)calloc((size_t)(ctx->NTran),sizeof(int))))
        goto AADFin;
    memrq(ctx, ctx->NTran,sizeof(int));
    ctx->SnTranA = ctx->NTran;

    if (!(ctx->OrgTran = (int *)calloc((size_t)(ctx->NTran),sizeof(int))))
        goto AADFin;
    memrq(ctx, ctx->NTran,sizeof(int));
    ctx->OrgTranA = ctx->NTran;

    if (!(ctx->DesTran = (int *)calloc((size_t)(ctx->NTran),sizeof(int))))
        goto AADFin;
    memrq(ctx, ctx->NTran,sizeof(int));
    ctx->DesTranA = ctx->NTran;

    if (!(ctx->SnTran1 = (int *)calloc((size_t)(ctx->NTran1 + 1),sizeof(int))))
        goto AADFin;
    ctx->SnTran1A = ctx->NTran1 + 1;
    memrq(ctx, ctx->SnTran1A,sizeof(int));

    if (!(ctx->OrgTran1 = (int *)calloc((size_t)(ctx->NTran1 + 1),sizeof(int))))
        goto AADFin;
    ctx->OrgTran1A = ctx->NTran1 + 1;
    memrq(ctx, ctx->OrgTran1A,sizeof(int));

    if (!(ctx->DesTran1 = (int *)calloc((size_t)(ctx->NTran1 + 1),sizeof(int))))
        goto AADFin;
    ctx->DesTran1A = ctx->NTran1 + 1;
    memrq(ctx, ctx->DesTran1A,sizeof(int));

    if (!(ctx->TranPtr = (int *)calloc((size_t)(ctx->MaxSnn) * (size_t)(ctx->MaxOrg1) + 1,sizeof(int))))
        goto AADFin;
    ctx->TranPtrA = ctx->MaxSnn * ctx->MaxOrg1 + 1;
    memrq(ctx, ctx->TranPtrA,sizeof(int));

    if (!(ctx->TranNE = (int *)calloc((size_t)(ctx->NTran),sizeof(int))))
        goto AADFin;
    ctx->TranNEA = ctx->NTran;
    memrq(ctx, ctx->TranNEA,sizeof(int));

    if (!(ctx->TranWE = (float *)calloc((size_t)(ctx->NTran),sizeof(float))))
        goto AADFin;
    ctx->TranWEA = ctx->NTran;
    memrq(ctx, ctx->TranWEA,sizeof(float));

    if (!(ctx->TranMD = (float *)calloc((size_t)(ctx->NTran),sizeof(float))))
        goto AADFin;
    ctx->TranMDA = ctx->NTran;
    memrq(ctx, ctx->TranMDA,sizeof(float));

    if (!(ctx->TranTSM = (float *)calloc((size_t)(ctx->NTran),sizeof(float))))
        goto AADFin;
    ctx->TranTSMA = ctx->NTran;
    memrq(ctx, ctx->TranTSMA,sizeof(float));

    if (!(ctx->TranTFM = (float *)calloc((size_t)(ctx->NTran),sizeof(float))))
        goto AADFin;
    ctx->TranTFMA = ctx->NTran;
    memrq(ctx, ctx->TranTFMA,sizeof(float));

    if (!(ctx->TWFreq = (float *)calloc((size_t)(ctx->NTran1),sizeof(float))))
        goto AADFin;
    ctx->TWFreqA = ctx->NTran1;
    memrq(ctx, ctx->TWFreqA,sizeof(float));

    if (!(ctx->TWDur = (float *)calloc((size_t)(ctx->NTran1),sizeof(float))))
        goto AADFin;
    ctx->TWDurA = ctx->NTran1;
    memrq(ctx, ctx->TWDurA,sizeof(float));

    return(0);

AADFin:
    if (err)
        p_err(ctx, -2,1);

    if (ctx->TWDurA) {
        free((char *)ctx->TWDur);
        memrq(ctx, -ctx->TWDurA,sizeof(float));
        ctx->TWDurA = 0;
    }
    if (ctx->TWFreqA) {
        free((char *)ctx->TWFreq);
        memrq(ctx, -ctx->TWFreqA,sizeof(float));
        ctx->TWFreqA = 0;
    }
    if (ctx->TranTFMA) {
        free((char *)ctx->TranTFM);
        memrq(ctx, -ctx->TranTFMA,sizeof(float));
        ctx->TranTFMA = 0;
    }
    if (ctx->TranTSMA) {
        free((char *)ctx->TranTSM);
        memrq(ctx, -ctx->TranTSMA,sizeof(float));
        ctx->TranTSMA = 0;
    }
    if (ctx->TranMDA) {
        free((char *)ctx->TranMD);
        memrq(ctx, -ctx->TranMDA,sizeof(float));
        ctx->TranMDA = 0;
    }
    if (ctx->TranWEA) {
        free((char *)ctx->TranWE);
        memrq(ctx, -ctx->TranWEA,sizeof(float));
        ctx->TranWEA = 0;
    }
    if (ctx->TranNEA) {
        free((char *)ctx->TranNE);
        memrq(ctx, -ctx->TranNEA,sizeof(int));
        ctx->TranNEA = 0;
    }
    if (ctx->TranPtrA) {
        free((char *)ctx->TranPtr);
        memrq(ctx, -ctx->TranPtrA,sizeof(int));
        ctx->TranPtrA = 0;
    }
    if (ctx->DesTran1A) {
        free((char *)ctx->DesTran1);
        memrq(ctx, -ctx->DesTran1A,sizeof(int));
        ctx->DesTran1A = 0;
    }
    if (ctx->OrgTran1A) {
        free((char *)ctx->OrgTran1);
        memrq(ctx, -ctx->OrgTran1A,sizeof(int));
        ctx->OrgTran1A = 0;
    }
    if (ctx->SnTran1A) {
        free((char *)ctx->SnTran1);
        memrq(ctx, -ctx->SnTran1A,sizeof(int));
        ctx->SnTran1A = 0;
    }
    if (ctx->DesTranA) {
        free((char *)ctx->DesTran);
        memrq(ctx, -ctx->DesTranA,sizeof(int));
        ctx->DesTranA = 0;
    }
    if (ctx->OrgTranA) {
        free((char *)ctx->OrgTran);
        memrq(ctx, -ctx->OrgTranA,sizeof(int));
        ctx->OrgTranA = 0;
    }
    if (ctx->SnTranA) {
        free((char *)ctx->SnTran);
        memrq(ctx, -ctx->SnTranA,sizeof(int));
        ctx->SnTranA = 0;
    }
    ctx->NTran = ctx->NTran1 = ctx->MaxSnn = ctx->MaxOrg = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_edat()      print table with information about episodes.            */

void prn_edat(TDAContext *ctx)
{
    register int i,k;
    int sn,sn1,org,org1,nn,nna;
    double sum,suma;
#ifdef TDA_R_PACKAGE
    /* the Excl column prints "-" or "*"; exported as 0/1 so the flag
       does not have to be read back out of the text */
    static const double zero_one[2] = { 0.0, 1.0 };
#endif

    prnchar(ctx, ' ',41,0);
    printf1(ctx, "Mean\n");
    printf1(ctx, "SN  Org Des   Episodes    Weighted     Duration     TS Min      TF Max  Excl\n");

    /* seed with the same sentinel as org1: the first loop iteration
       always enters the new-group branch through org != org1, and
       SnTran has only NTran entries, so reading SnTran[1] walked past
       the end whenever a dataset has a single transition */
    sn1 = -999;
    nna = nn = k = 0;
    org1 = -999;
    suma = sum = 0.0;

    for (i = 0; i < ctx->NTran; ++i) {

        sn = ctx->SnTran[i];
        org = ctx->OrgTran[i];
        if (sn != sn1 || org != org1) {
            if (nn > 0) {
                printf1(ctx, "Sum         %9d %12.2lf\n",nn,sum);
                nn = 0;
                sum = 0.0;
            }
            prnchar(ctx, '-',76,1);
            sn1 = sn;
            org1 = org;
        }
        printf1(ctx, "%2d %4d %3d %9d %12.2lf %11.2lf %11.2lf %11.2lf    ",
            sn,org,ctx->DesTran[i],ctx->TranNE[i],(double)ctx->TranWE[i],
            (double)ctx->TranMD[i],(double)ctx->TranTSM[i],(double)ctx->TranTFM[i]);
#ifdef TDA_R_PACKAGE
        {
            double erow[8];
            erow[0] = (double)sn;
            erow[1] = (double)org;
            erow[2] = (double)ctx->DesTran[i];
            erow[3] = (double)ctx->TranNE[i];
            erow[4] = (double)ctx->TranWE[i];
            erow[5] = (double)ctx->TranMD[i];
            erow[6] = (double)ctx->TranTSM[i];
            erow[7] = (double)ctx->TranTFM[i];
            tda_export_row(ctx, "episodes.table", erow, 8);
        }
#endif

        if (ctx->TranPtr[(sn - 1) * ctx->MaxOrg1 + org] >= 0) {
            printf1(ctx, "-\n");
#ifdef TDA_R_PACKAGE
            tda_export_row(ctx, "episodes.excluded", &zero_one[0], 1);
#endif
        }
        else {
            printf1(ctx, "*\n");
#ifdef TDA_R_PACKAGE
            tda_export_row(ctx, "episodes.excluded", &zero_one[1], 1);
#endif
            k++;
        }
        nn += ctx->TranNE[i];
        nna += ctx->TranNE[i];
        sum += (double)ctx->TranWE[i];
        suma += (double)ctx->TranWE[i];
    }
    printf1(ctx, "Sum         %9d %12.2lf\n",nn,sum);
    if (nna > nn) {
        prnchar(ctx, '-',76,1);
        printf1(ctx, "Sum         %9d %12.2lf\n",nna,suma);
    }
    if (k) {
        printf1(ctx, "\nWarning: episodes with * in column Excl do not belong to any transition\n");
        printf1(ctx, "and will be excluded from all estimation procedures.\n");
    }
#ifdef TDA_R_PACKAGE
    tda_export_flush(ctx, "episodes.table");
    tda_export_flush(ctx, "episodes.excluded");
#endif
}

/* ------------------------------------------------------------------------ */
/*  get_edat(i,sn,org,des,ts,tf)                                            */
/*                                                                          */
/*  get episode data for case i. Return 0 if OK, -1 if error.               */

int get_edat(TDAContext *ctx, int i,int *sn,int *org,int *des,double *ts,double *tf)
{
    int r;
    double x;

    if (ctx->TSVar == -1)           /* get starting time */
        *ts = ctx->TSESVal[0];
    else if (ctx->TSVar >= 0)
        *ts = get_data(ctx, ctx->TSVar,i);
    else {
        r = v_eval1(ctx, i,ctx->TSESCnt,ctx->TSESTyp,ctx->TSESVal,ctx->ESIdx,ts,0,0,0,0,0);
        if (r)
            return(-1);
    }
    if (ctx->TFVar == -1)           /* get ending time */
        *tf = ctx->TFESVal[0];
    else if (ctx->TFVar >= 0)
        *tf = get_data(ctx, ctx->TFVar,i);
    else {
        r = v_eval1(ctx, i,ctx->TFESCnt,ctx->TFESTyp,ctx->TFESVal,ctx->ESIdx,tf,0,0,0,0,0);
        if (r)
            return(-1);
    }
    if (ctx->ORGVar == -1)           /* get origin state */
        *org = (int)ctx->ORGESVal[0];
    else if (ctx->ORGVar >= 0)
        *org = (int)get_data(ctx, ctx->ORGVar,i);
    else {
        r = v_eval1(ctx, i,ctx->ORGESCnt,ctx->ORGESTyp,ctx->ORGESVal,ctx->ESIdx,&x,0,0,0,0,0);
        if (r)
            return(-1);
        *org = (int)x;
    }
    if (ctx->DESVar == -1)           /* get destination state */
        *des = (int)ctx->DESESVal[0];
    else if (ctx->DESVar >= 0)
        *des = (int)get_data(ctx, ctx->DESVar,i);
    else {
        r = v_eval1(ctx, i,ctx->DESESCnt,ctx->DESESTyp,ctx->DESESVal,ctx->ESIdx,&x,0,0,0,0,0);
        if (r)
            return(-1);
        *des = (int)x;
    }
    if (ctx->MEFlg) {      /* multiepisode data */

        if (ctx->SNVar == -1)           /* get spell number */
            *sn = (int)ctx->SNESVal[0];
        else if (ctx->SNVar >= 0)
            *sn = (int)get_data(ctx, ctx->SNVar,i);
        else {
            r = v_eval1(ctx, i,ctx->SNESCnt,ctx->SNESTyp,ctx->SNESVal,ctx->ESIdx,&x,0,0,0,0,0);
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

int get_id(TDAContext *ctx, int icase,int *err)
{
    int r,id;
    double x;

    *err = 0;
    id = icase + 1;

    if (ctx->MEFlg) {
        if (ctx->IDVar == -1)           /* get spell number */
            id = (int)ctx->IDESVal[0];
        else if (ctx->IDVar >= 0)
            id = (int)get_data(ctx, ctx->IDVar,icase);
        else {
            r = v_eval1(ctx, icase,ctx->IDESCnt,ctx->IDESTyp,ctx->IDESVal,ctx->ESIdx,&x,0,0,0,0,0);
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

int sort_ed(TDAContext *ctx, int opt)
{
    register int i;
    int err,org,des,sn;
    double ts,tf;

    err = 0;
    if (opt == 0)
        goto SEDFin;

    if (!(ctx->TSIdx = (int *)calloc((size_t)(ctx->NOC),sizeof(int))))
        goto SEDFin;
    ctx->TSSortA = ctx->NOC;
    memrq(ctx, ctx->NOC,sizeof(int));

    if (!(ctx->TFIdx = (int *)calloc((size_t)(ctx->NOC),sizeof(int))))
        goto SEDFin;
    ctx->TFSortA = ctx->NOC;
    memrq(ctx, ctx->NOC,sizeof(int));

    for (i = 0; i < ctx->NOC; ++i)
        ctx->TSIdx[i] = i;

    if (alloc_acu(ctx, ctx->NOC))
        goto SEDFin;
                                /* do not sort if ts = constant */

    if (ctx->TSVar == -1 || (ctx->TSVar >= 0 && ctx->VTyp[ctx->TSVar] == 2))
        ;
    else {
        printf1(ctx, "Sorting episodes according to starting times.\n");

        for (i = 0; i < ctx->NOC; ++i) {
            get_edat(ctx, i,&sn,&org,&des,&ts,&tf);
            ctx->AcU[i] = ts;
        }
        tda_qsort_r((char *)ctx->TSIdx,(size_t)(ctx->NOC),sizeof(int), tscomp, ctx);
    }
    for (i = 0; i < ctx->NOC; ++i)
        ctx->TFIdx[i] = i;
                                /* do not sort if tf = constant */

    if (ctx->TFVar == -1 || (ctx->TFVar >= 0 && ctx->VTyp[ctx->TFVar] == 2))
        ;
    else {
        printf1(ctx, "Sorting episodes according to ending times.\n");

        if (alloc_aci(ctx, ctx->NOC))
            goto SEDFin;
        if (alloc_acj(ctx, ctx->NOC))
            goto SEDFin;

        for (i = 0; i < ctx->NOC; ++i) {
            get_edat(ctx, i,&sn,&org,&des,&ts,&tf);
            ctx->AcU[i] = tf;
            ctx->AcI[i] = org;
            ctx->AcJ[i] = des;
        }
        tda_qsort_r((char *)ctx->TFIdx,(size_t)(ctx->NOC),sizeof(int), tfcomp, ctx);
    }
    ctx->ESortFlg = 1;
    return(0);

SEDFin:
    alloc_acu(ctx, 0);
    alloc_acv(ctx, 0);
    if (err)
        p_err(ctx, -2,1);

    if (ctx->TSSortA) {
        free((char *)ctx->TSIdx);
        memrq(ctx, -ctx->TSSortA,sizeof(int));
        ctx->TSSortA = 0;
    }
    if (ctx->TFSortA) {
        free((char *)ctx->TFIdx);
        memrq(ctx, -ctx->TFSortA,sizeof(int));
        ctx->TFSortA = 0;
    }
    ctx->ESortFlg = 0;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  tscomp  Compare function for quicksort. The compare is according to the */
/*          starting times TS.                                              */

int tscomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;
    double ts1,ts2;

    ts1 = ctx->AcU[*(int *)arg1];
    ts2 = ctx->AcU[*(int *)arg2];
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

int tfcomp(const void *arg1, const void *arg2, void *_ctx)
{
    TDAContext *ctx = (TDAContext *)_ctx;
    int org1,org2,des1,des2;
    double tf1,tf2;

    tf1 = ctx->AcU[*(int *)arg1];
    tf2 = ctx->AcU[*(int *)arg2];

    if (tf1 > tf2)
        return(1);
    else if (tf1 < tf2)
        return(-1);
    else {

        org1 = ctx->AcI[*(int *)arg1];
        org2 = ctx->AcI[*(int *)arg2];
        des1 = ctx->AcJ[*(int *)arg1];
        des2 = ctx->AcJ[*(int *)arg2];

        if (org1 != org2)       /* if O states are different ordering   */
            return(0);          /* is not important.                    */
        else {
            int ev1 = (org1 != des1);   /* 1 if event, 0 if censored */
            int ev2 = (org2 != des2);
            if (ev1 != ev2)
                return(ev1 ? -1 : 1);   /* events before censored */
            /* fully tied: use original case index for determinism */
            return (*(int *)arg1) - (*(int *)arg2);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_spell   Get the next spell, or split.                               */
/*                                                                          */

int get_spell(TDAContext *ctx, int init,int *icase,int *sn,int *org,int *des, double *ts, double *tf,int *spl,int *nspl)
{
    double t,ttf,tl,ts1,tf1;
    register int j,k,m;

    if (init) {
        ctx->NSPLEvalErr = ctx->NSplits = 0;
        ctx->s_get_spell_isplit = ctx->s_get_spell_nsplit = 0;
        ctx->s_get_spell_ic = -1;
        return(0);
    }
    if (++ctx->s_get_spell_isplit > ctx->s_get_spell_nsplit) {

        if (++ctx->s_get_spell_ic >= ctx->NOC)         /* no more cases */
            return(0);

        j = 0;
        get_edat(ctx, ctx->s_get_spell_ic,&ctx->s_get_spell_sn1,
                 &ctx->s_get_spell_org1,&ctx->s_get_spell_des1,&ts1,&tf1);

        ctx->SPTTV[++j] = ts1;
        tl = ttf = tf1;
        m = 0;

        while (m >= 0) {
            m = -1;
            for (k = 0; k < ctx->NSP; ++k) {
                t = get_data(ctx, ctx->SPLITVar[k],ctx->s_get_spell_ic);
                if (t > ctx->SPTTV[j] && t < tl) {
                    m = ctx->SPLITVar[k];
                    tl = t;
                }
            }
            if (m >= 0) {
                ctx->SPTTV[++j] = tl;
                ctx->SPIDX[j] = (short)(m);
                tl = ttf;
            }
        }
        ctx->s_get_spell_nsplit = j;
        ctx->s_get_spell_isplit = 1;

        ctx->SPTTV[++j] = ttf;
    }
    *icase = ctx->s_get_spell_ic;

    *sn = ctx->s_get_spell_sn1;
    *ts = ctx->SPTTV[ctx->s_get_spell_isplit];
    *tf = ctx->SPTTV[ctx->s_get_spell_isplit + 1];
    if (ctx->s_get_spell_isplit == 1)
        *org = ctx->s_get_spell_stat = ctx->s_get_spell_org1;
    else
        *org = ctx->s_get_spell_stat;

    if (ctx->s_get_spell_isplit == ctx->s_get_spell_nsplit)
        *des = ctx->s_get_spell_des1;
    else
        *des = ctx->s_get_spell_stat;

    *spl  = ctx->s_get_spell_isplit;
    *nspl = ctx->s_get_spell_nsplit;

    /* evaluate type 5 variables */

    ctx->EDVALSn  = *sn;
    ctx->EDVALOrg = *org;
    ctx->EDVALDes = *des;
    ctx->EDVALTs  = *ts;
    ctx->EDVALTf  = *tf;
    ctx->EDVALTime = *tf;

    if (ctx->NVAR5 > 0) {
        j = ctx->VIFirst;
        while (j >= 0) {

            if (ctx->VTyp[j] == 5) {

                k = v_eval1(ctx, ctx->s_get_spell_ic,ctx->VESCnt[j],ctx->VESTyp[j],ctx->VESVal[j],ctx->ESIdx,&t,0,0,0,0,0);
                if (k) {
                    ctx->NSPLEvalErr++;
                    ctx->NPFlgs[5] += 1;
                    t = 0.0;
                }
                put_data(ctx, t,j,ctx->s_get_spell_ic);
                /*** printf1("put t=%lg into j=%d ctx->s_get_spell_ic=%d\n",t,j,ctx->s_get_spell_ic); **/
            }
            j = ctx->VNxt[j];
        }
    }
    ctx->NSplits++;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_nedat()   check episode data.                                     */
/*                  return 0 if OK, -1 if error.                            */

int check_nedat(TDAContext *ctx)
{
    int nrec,icase,sn,org,des,spl,nspl,id,id1,idn,r;
    double ts,tf;

    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);

    id1 = ctx->INTMAX;
    idn = nrec = 0;
    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        if (ctx->MEFlg) {
            id = get_id(ctx, icase,&r);
            if (r) {
                printf1(ctx, "\nError in evaluating ID variable in case %d\n",icase + 1);
                return(-1);
            }
            if (id != id1) {
                idn++;
                id1 = id;
            }
        }
        nrec++;
    }
    newline(ctx);
    if (ctx->MEFlg)
        printf1(ctx, "Number of individuals: %d\n",idn);
    printf1(ctx, "Number of episodes: %d\n",ctx->NOC);
    if (ctx->NSP > 0)
        printf1(ctx, "Number of splits: %d\n",nrec);

    if (ctx->NVAR5 > 0) {
        printf1(ctx, "Errors in evaluating type 5 variables: %d\n",ctx->NSPLEvalErr);
        if (ctx->NSPLEvalErr > 0)
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

int epdat(TDAContext *ctx)
{
    register int i,j;
    int err,nrec,icase,sn,org,des,spl,nspl,id,r,nr;
    double ts,tf;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    if (ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 5,1,1))
        goto EPDATFin;

    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing episode data to: %s\n",ctx->PMFdName);
    nr = nrec = 0;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 6,2);

    get_spell(ctx, 1,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl);

    while (get_spell(ctx, 0,&icase,&sn,&org,&des,&ts,&tf,&spl,&nspl)) {

        id = get_id(ctx, icase,&r);
        if (r)
            nr++;

        fprintf(ctx->PMFd,"%6d %3d %3d %3d %3d %3d ",id,sn,nspl,spl,org,des);
#ifdef TDA_R_PACKAGE
        /* the episode row epdat writes; its width depends on how many
           v= variables were asked for */
        tda_export_cell(ctx, "epdat.table", (double)id);
        tda_export_cell(ctx, "epdat.table", (double)sn);
        tda_export_cell(ctx, "epdat.table", (double)nspl);
        tda_export_cell(ctx, "epdat.table", (double)spl);
        tda_export_cell(ctx, "epdat.table", (double)org);
        tda_export_cell(ctx, "epdat.table", (double)des);
#endif
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,ts);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,tf);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "epdat.table", ts);
        tda_export_cell(ctx, "epdat.table", tf);
#endif

        for (i = 0; i < ctx->PMNV; ++i) {
            j = ctx->PMVIdx[i];
            rt_fprintf_d(ctx, ctx->PMFd,ctx->VPFmtS[j],get_data(ctx, j,icase));
#ifdef TDA_R_PACKAGE
            tda_export_cell(ctx, "epdat.table", get_data(ctx, j,icase));
#endif
        }
        fprintf(ctx->PMFd,"\n");
#ifdef TDA_R_PACKAGE
        tda_export_endrow(ctx, "epdat.table");
#endif
        nrec++;
        if (ctx->PMNOCFlg && nrec >= ctx->PMNOC)
            break;
    }
    printf1(ctx, "Episode data: %d records written to: %s\n",nrec,ctx->PMFdName);

    if (ctx->PMTDAFDef)
        edtda(ctx, ctx->PMFdName,nrec);

    if (nr > 0)
        printf1(ctx, "Warning: %d errors in evaluating ID variable (ignored).\n",nr);

    err = 0;

EPDATFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  edtda(fname,noc)                                                        */

void edtda(TDAContext *ctx, char *fname,int noc)
{
    register int i,j;

    fprintf(ctx->PMTDAFd,"nvar(\n");
    fprintf(ctx->PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(ctx->PMTDAFd,"  noc = %d,\n",noc);
    fprintf(ctx->PMTDAFd,"  ID   [6.0] = c1, # id number\n");
    fprintf(ctx->PMTDAFd,"  SN   [3.0] = c2, # spell number\n");
    fprintf(ctx->PMTDAFd,"  NSPL [3.0] = c3, # number of splits\n");
    fprintf(ctx->PMTDAFd,"  SPL  [3.0] = c4, # split number\n");
    fprintf(ctx->PMTDAFd,"  ORG  [3.0] = c5, # origin state\n");
    fprintf(ctx->PMTDAFd,"  DES  [3.0] = c6, # destination state\n");
    fprintf(ctx->PMTDAFd,"  TS   [%d.%d] = c7, # starting time\n",ctx->PMFmt1,ctx->PMFmt2);
    fprintf(ctx->PMTDAFd,"  TF   [%d.%d] = c8, # ending time\n",ctx->PMFmt1,ctx->PMFmt2);

    for (i = 0; i < ctx->PMNV; ++i) {
        j = ctx->PMVIdx[i];
        fprintf(ctx->PMTDAFd,"  %s [%d.%d] = c%d,\n",ctx->VName[j],ctx->VPFmt1[j],ctx->VPFmt2[j],i + 9);
    }
    fprintf(ctx->PMTDAFd,");\n");
    printf1(ctx, "TDA description written to: %s\n",ctx->PMTDAFName);
}

/* ------------------------------------------------------------------------ */
/*  epsdat()    epsdat command in CmdBuf: epsdat(...) = fname               */
/*              repuired parameter: t= time points.                         */
/*              write distribution of states to output file.                */
/*              return 0 if OK, -1 if error.                                */

int epsdat(TDAContext *ctx)
{
    register int i,j,k;
    int err,m = 0,n,sn,org,des,r,id,id1,s;
    double ts,tf,t;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    if (ctx->EDAvail == 0) {
        p_err(ctx, -15,1);
        return(0);
    }
    if (parm(ctx, ctx->CmdBuf + 6,1,1))
        goto EPSDATFin;

    if (ctx->PMNTP == 0) {
        p_err(ctx, -17,1);
        goto EPSDATFin;
    }
    if (ctx->SILENTFlg < 2)
        printfe(ctx, "Writing state distributions to: %s\n",ctx->PMFdName);
    printf1(ctx, "Writing state distributions to: %s\nUsing ",ctx->PMFdName);
    if (ctx->MEFlg)
        printf1(ctx, "multi-");
    else
        printf1(ctx, "single ");
    printf1(ctx, "episode data.\n");

    if (alloc_acn(ctx, ctx->MaxOrg + 1))
        goto EPSDATFin;

    if (alloc_aci(ctx, ctx->MaxOrg + 1))
        goto EPSDATFin;

    for (i = 0; i < ctx->NTran; ++i)
        ctx->AcI[ctx->OrgTran[i]] = 1;

    fprintf(ctx->PMFd,"# State distributions.\n");
    fprintf(ctx->PMFd,"# Time    ");
    for (i = 0; i <= ctx->MaxOrg; ++i) {
        if (ctx->AcI[i])
            fprintf(ctx->PMFd," %6d",i);
    }
    fprintf(ctx->PMFd,"    Total  Missing\n");


    for (k = 0; k < ctx->PMNTP; ++k) {        /* for all time points */

        t = ctx->PMTP[k];
        for (j = 0; j <= ctx->MaxOrg; ++j)
            ctx->AcN[j] = 0;

        m = n = 0;
        if (ctx->MEFlg == 0) {
            for (i = 0; i < ctx->NOC; ++i) {

                if (get_edat(ctx, i,&sn,&org,&des,&ts,&tf)) {
                    p_err(ctx, -39,1);
                    goto EPSDATFin;
                }
                if (ts <= t && t < tf) {
                    ctx->AcN[org] += 1;
                    n++;
                }
                m++;
            }
        }
        else {

            s = -1;
            id1 = ctx->INTMAX;
            for (i = 0; i < ctx->NOC; ++i) {

                id = get_id(ctx, i,&r);
                if (id != id1) {
                    if (s >= 0) {
                        ctx->AcN[s] += 1;
                        n++;
                        s = -1;
                    }
                    m++;
                    id1 = id;
                }
                if (r || get_edat(ctx, i,&sn,&org,&des,&ts,&tf)) {
                    p_err(ctx, -39,1);
                    goto EPSDATFin;
                }
                if (s < 0 && ts <= t && t < tf)
                    s = org;
            }
            if (s >= 0) {
                ctx->AcN[s] += 1;
                n++;
            }
        }
        fprintf(ctx->PMFd,"%10.4lf",t);
#ifdef TDA_R_PACKAGE
        /* time, one count per state actually present, then total and
           missing -- the width depends on which states occur */
        tda_export_cell(ctx, "epsdat.table", t);
#endif
        for (j = 0; j <= ctx->MaxOrg; ++j) {
            if (ctx->AcI[j]) {
                fprintf(ctx->PMFd," %6d",ctx->AcN[j]);
#ifdef TDA_R_PACKAGE
                tda_export_cell(ctx, "epsdat.table", (double)ctx->AcN[j]);
#endif
            }
        }
        fprintf(ctx->PMFd," %8d %8d\n",n,m - n);
#ifdef TDA_R_PACKAGE
        tda_export_cell(ctx, "epsdat.table", (double)n);
        tda_export_cell(ctx, "epsdat.table", (double)(m - n));
        tda_export_endrow(ctx, "epsdat.table");
#endif
    }
    if (ctx->MEFlg)
        printf1(ctx, "Number of individuals: %d\n",m);
    else
        printf1(ctx, "Number of single episodes: %d\n",m);

    err = 0;

EPSDATFin:
    p_clean(ctx);
    return(err);
}
