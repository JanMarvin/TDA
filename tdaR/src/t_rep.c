/****************************************************************************/
/*  t_rep                                                                   */
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
#include "t_alloc.h"
#include "t_var.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_matc.h"
#include "t_cmd.h"
#include "t_gf.h"
#include "tda_context.h"

/*  functions in t_rep.c */

int t_exec(TDAContext *ctx, char *cmd);           
int t_repeat(TDAContext *ctx, int typ,int lev);
int t_endrepeat(TDAContext *ctx, int typ,int lev);
int rep_alloc(TDAContext *ctx, int opt,int lev,int typ); 
void rep_free(TDAContext *ctx);
int rep_scmd(TDAContext *ctx, char *cmd);
void rep_ptyp(TDAContext *ctx, int typ,int opt);
void rep_err(TDAContext *ctx, char *s,int opt);
int rep_expr(TDAContext *ctx, int lev,double *res);
int check_break(TDAContext *ctx, char *cmd);   
int check_if(TDAContext *ctx, char *cmd);   
void prn_nex(TDAContext *ctx, char *cmd);   

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

                            /* 0 : while                                    */
                            /* 1 : repeat                                   */



/*--------------------------------------------------------------------------*/
/*  t_repeat(typ,lev)   Repeat command on level lev.                        */
/*                      typ 0 : while                                       */
/*                      typ 1 : repeat                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int t_repeat(TDAContext *ctx, int typ,int lev)
{
    int err,m,idx,row,col,aflag,ivflg;
    register char *p,*q;
    char mname[VNLMax+1];
    double tmp;

    ivflg = aflag = 0;
    err = -1;             
    rep_ptyp(ctx, typ,0);
    printf2(ctx, "[level %d]. Current memory: %d bytes.\n",++lev,ctx->MemReq);

    if (lev >= MaxREP) {
        printf1(ctx, "Error: exceeded maximal level of loops.\n");
        return(-1);
    }
    if (rep_alloc(ctx, 1,lev,typ))
        goto NREPFin;  

    aflag = 1;
    ctx->REPTyp[lev] = (short)(typ);

    if (typ == 0) {                     /* while */

        p = ctx->CmdBuf + 5;
        if (*p++ != '(' || *(q = skip_expr(ctx, p)) != ')') {
            rep_err(ctx, ctx->CmdBuf,0);
            goto NREPFin;
        }
        *q = '\0';          
        m = (int)(strlen(p) + 1);
        if (!(ctx->REPExpr[lev] = (char *)calloc((size_t)(m),sizeof(char *)))) {
            p_err(ctx, -2,1);
            goto NREPFin;  
        }
        memrq(ctx, m,sizeof(char *));
        ctx->REPExprA[lev] = (short)(m);           
        strcpy(ctx->REPExpr[lev],p);

        if (rep_expr(ctx, lev,&tmp))
            goto NREPFin;
            
        printf2(ctx, "expression: %s  Initial value: %lg",ctx->REPExpr[lev],tmp);
        if (tmp < ctx->EPSI1) {
            ctx->REPCond[lev] = 1;
            printf2(ctx, " [false]");
        }
        printf2(ctx, "\n");
    }
    else {                              /* repeat */
        idx = -1;
        p = ctx->CmdBuf;
        if (strncmp(p,"repeat(n=",9)) {
            rep_err(ctx, ctx->CmdBuf,0);
            goto NREPFin;
        }           
        if ((p = n_mexpr(ctx, ctx->CmdBuf + 9,1,0,&row,&col,&ivflg,NULL)) == NULL || row * col != 1) {    
            rep_err(ctx, ctx->CmdBuf,0);
            goto NREPFin;
        }           
        m = (int)ctx->MX[0][1];
        if (m < 1) {
            rep_err(ctx, ctx->CmdBuf,0);
            goto NREPFin;
        }           
        if (*p == ',') {
            if (get_mname(ctx, ++p,mname,1) == NULL)
                goto NREPFin;                     
            if (mat_getidx(ctx, mname,0) >= 0) {
                rep_err(ctx, ctx->CmdBuf,1);
                goto NREPFin;
            }           
            if (m_getmat(ctx, p,1,1,&idx,ctx->CmdBuf,1) == NULL) 
                goto NREPFin;
            ctx->MatVal[idx][1] = 1;
            ctx->REPIdx[lev] = (short)(idx + 1);
        }
        else if (*p != ')') {
            rep_err(ctx, ctx->CmdBuf,0);
            goto NREPFin;
        }
        ctx->REPICNT[lev] = m;
        if (ctx->SILENTFlg < 0) {
            printf2(ctx, "Repeat for %d iterations.",m);
            if (idx >= 0)
                printf2(ctx, " Created matrix: %s",ctx->MatName[idx]);
            printf2(ctx, "\n");
        }
    }
    err = 0;

NREPFin:
    mx_free(ctx);
    if (err && aflag)  
        rep_alloc(ctx, 0,lev,0);
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  t_endrepeat(typ,lev)    Execute loop on level lev.                      */
/*                          typ 0 : while                                   */
/*                          typ 1 : repeat                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int t_endrepeat(TDAContext *ctx, int typ,int lev)
{
    register int i;
    int err,icmd,iter = 0;
    double tmp;

    err = 0;             
    i = ctx->REPTyp[lev];
    if (i > 1)
        i = 1;

    if (lev < 0 || typ != i) {
        printf1(ctx, "Error: unbalanced end");
        rep_ptyp(ctx, typ,1);
        printf1(ctx, "\n");
        err = -1;
        goto EREPFin;
    }
    /**  
    for (i = 0; i < REPCCnt[lev]; ++i) {
        printf1(ctx, "i=%d lev=%d ex: %s\n",i,lev,REPCmd[lev][i]);
    }
    **/  

    if (ctx->BREAKFlg)
        goto EREPFin1;

    /* execute current level */
   
    ctx->REPSFlg[lev] = 2;
    iter = 0;
    while (++iter > 0) {                

        /* check condition for break */

        if (ctx->REPTyp[lev] == 0) {                      /* while */
   
            if (rep_expr(ctx, lev,&tmp)) {
                err = -1;
                goto EREPFin1;
            }
            if (tmp < ctx->EPSI1)
                goto EREPFin1;
        }
        else {                                      /* repeat */
            ctx->REPICNT[lev] -= 1;
            if (ctx->REPICNT[lev] <= 0)
                goto EREPFin1;
            if (ctx->REPIdx[lev] > 0)
                ctx->MatVal[ctx->REPIdx[lev] - 1][1] += 1.0;
        }
        rep_ptyp(ctx, typ,0);
        if (ctx->SILENTFlg < 0) {
            printf2(ctx, "[level %d, iteration %d]. Current memory: %d bytes.\n",
                                                     lev,iter,ctx->MemReq);
            prnchar(ctx, '-',LLEN,1);      
        }
        for (icmd = 0; icmd < ctx->REPCCnt[lev] - 1; ++icmd) {

            strcpy(ctx->CmdBuf,ctx->REPCmd[lev][icmd]);

            /* don't save commands for lower levels */

            for (i = 0; i < lev; ++i)
                ctx->REPSFlg[i] = (short)(-iabs(ctx, (int)ctx->REPSFlg[i]));

            err = t_exec(ctx, ctx->CmdBuf);

            for (i = 0; i < lev; ++i)    
                ctx->REPSFlg[i] = (short)(iabs(ctx, (int)ctx->REPSFlg[i]));

            if (err)
                break;         
        }
        if (err || ctx->BREAKFlg)
            goto EREPFin1;
    }

EREPFin1:
    rep_alloc(ctx, 0,lev,0);
     
    if (err == 0) {
        printf2(ctx, "end");
        rep_ptyp(ctx, typ,0);
        printf2(ctx, "[level %d, iteration %d]. Current memory: %d bytes.\n",
                                                         lev,iter,ctx->MemReq);
    }

EREPFin:
    ctx->BREAKFlg = 0;
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  rep_alloc(opt,lev,typ)      If opt != 0 allocate for level lev,         */
/*                              otherwise free.                             */
/*  Return 0 if OK, -1 if error.                                            */

int rep_alloc(TDAContext *ctx, int opt,int lev,int typ)
{
    register int i,l;

    /**
    printf1(ctx, "repalloc opt=%d lev=%d typ=%d MemReq=%d REPLev=%d\n",opt,lev,typ,MemReq,REPLev); 
    **/ 
   
    if (opt == 0) {
        if (lev < ctx->REPLev) {
            printfe(ctx, "ERROR1 in rep_alloc. opt=%d lev=%d replev=%d\n",opt,lev,ctx->REPLev);
            gerr_exit(ctx, 201);
        }
        if (ctx->REPIdx[lev] > 0) {
            mat_alloc(ctx, ctx->REPIdx[lev] - 1,0);
            ctx->REPIdx[lev] = 0;
        }
        if (ctx->REPExprA[lev] > 0) {         
            free((char *)ctx->REPExpr[lev]);
            memrq(ctx, -ctx->REPExprA[lev],sizeof(char *));
            ctx->REPExprA[lev] = 0;           
        }
        if (ctx->REPCmdA[lev] > 0) {
            for (i = 0; i < ctx->REPCCnt[lev]; ++i) {
                l = (int)(strlen(ctx->REPCmd[lev][i]));
                free((char *)ctx->REPCmd[lev][i]);
                memrq(ctx, -l - 1,sizeof(char *));
            }
            free((char *)ctx->REPCmd[lev]);
            memrq(ctx, -ctx->REPCmdA[lev],sizeof(char **));
            ctx->REPCCnt[lev] = ctx->REPCmdA[lev] = 0;          
            ctx->REPTyp[lev] = ctx->REPSFlg[lev] = 0;
            ctx->REPLev--;
        }
        /**
        printf1(ctx, "ret1 repalloc opt=%d lev=%d typ=%d MemReq=%d \n",opt,lev,typ,MemReq); 
        **/
        return(0);
    }
    if (lev != ctx->REPLev + 1) {
        printfe(ctx, "ERROR2 in rep_alloc. opt=%d lev=%d replev=%d\n",opt,lev,ctx->REPLev);
        gerr_exit(ctx, 202);
    }
    if (!(ctx->REPCmd[lev] = (char **)calloc(MaxREPCMD,sizeof(char **)))) {
        p_err(ctx, -2,1);
        return(-1);    
    }
    memrq(ctx, MaxREPCMD,sizeof(char **));
    ctx->REPCmdA[lev] = MaxREPCMD;
    ctx->REPTyp[lev] = (short)(typ);
    ctx->REPSFlg[lev] = 1;
    ctx->REPCond[lev] = 0;
    ctx->REPLev = lev;
    /**
    printf1(ctx, "ret2 repalloc opt=%d lev=%d typ=%d MemReq=%d \n",opt,lev,typ,MemReq); 
    **/
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  rep_free()      Free all memory for repeat.                             */

void rep_free(TDAContext *ctx)
{
    while (ctx->REPLev >= 0)  
        rep_alloc(ctx, 0,ctx->REPLev,0);
}

/*--------------------------------------------------------------------------*/
/*  rep_scmd(cmd)       Save command for repeat loop.                       */
/*  Return 0 if OK, -1 if error, 1 if command should be skipped.            */

int rep_scmd(TDAContext *ctx, char *cmd)
{
    int l,lev;
  
    if (ctx->REPCond[ctx->REPLev]) {
   
        if (ctx->s_rep_scmd_skiplev < 0)  
            ctx->s_rep_scmd_skiplev++;             
         
        if (!strncmp(cmd,"repeat",6) || !strncmp(cmd,"while(",6)) 
            ctx->s_rep_scmd_skiplev++;

        else if (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8))                              
            ctx->s_rep_scmd_skiplev--;

        if (ctx->s_rep_scmd_skiplev >= 0)
            return(1);
    }
    for (lev = 0; lev <= ctx->REPLev; ++lev) {

        /********
        tda_out("SAVE lev=%d REPLev=%d REPSFlg[lev]=%d cnt=%d cmd: %s\n",lev,REPLev,REPSFlg[lev],REPCCnt[lev],cmd);
        *********/

        if (ctx->REPSFlg[lev] != 1)
            continue;

        /******
        if (lev == REPLev && (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8)))                   
            continue;
        ***/

        if (ctx->REPCCnt[lev] >= MaxREPCMD) {
            printf1(ctx, "Error: exceeded max number of commands in loop.\n");
            return(-1);
        }
        l = (int)(strlen(cmd) + 1);
        if (!(ctx->REPCmd[lev][ctx->REPCCnt[lev]] = (char *)calloc((size_t)(l),sizeof(char *)))) {
            p_err(ctx, -2,1);
            return(-1);    
        }
        memrq(ctx, l,sizeof(char *));
        strcpy(ctx->REPCmd[lev][ctx->REPCCnt[lev]],cmd);
        ctx->REPCCnt[lev] += 1;  
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  rep_ptyp(typ,opt)   print type of loop.                                 */

void rep_ptyp(TDAContext *ctx, int typ,int opt)
{
    (void)opt;        /* unused: the signature is shared */
    if (typ == 0)
        printf2(ctx, "while ");
    else 
        printf2(ctx, "repeat ");
}

/*--------------------------------------------------------------------------*/
/*  rep_err(s,opt)  print error message.                                    */

void rep_err(TDAContext *ctx, char *s,int opt)
{
    printf1(ctx, "Error: %s\n",s);
    if (opt == 1)
        printf1(ctx, "Matrix already exists.\n");
}

/*--------------------------------------------------------------------------*/
/*  rep_expr(lev,res)   Evaluate while expression on level lev.             */
/*                      Return result in res. Return 0 if OK, -1 if error.  */

int rep_expr(TDAContext *ctx, int lev,double *res)
{
    int err,row,col,ivflg;
    double tmp[2];

    err = -1;
    ivflg = 0;
    if (get_mexpr(ctx, ctx->REPExpr[lev],&row,&col,&ivflg))
        goto REPExprErr;

    if (row != 1 || col != 1) {
        printf1(ctx, "Error: while() needs a scalar expression.\n");
        goto REPExprFin;
    }
    if (eval_mexpr(ctx, 0,row,col,tmp,0,tmp,NULL,NULL,-1,0.0))
        goto REPExprErr;

    *res = fabs(tmp[1]);
    err = 0;

REPExprFin:
    return(err);                 

REPExprErr:
    alloc_acx(ctx, 0);
    printf1(ctx, "Error: can't evaluate while expression (level %d).\n",lev);
    printf1(ctx, "Expression: %s\n",ctx->REPExpr[lev]);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  check_break(cmd)    Check for break. Return 1 if continue, else 0.      */
/*                      -1 if error.                                        */

int check_break(TDAContext *ctx, char *cmd)    
{
    if (!strcmp(cmd,"break")) {
        if (ctx->REPLev < 0) {           /* must be in loop */
            printf1(ctx, "Error: misplaced break.\n");
            return(-1);
        }
        ctx->BREAKFlg = 1;
        printf2(ctx, "break\n\n");
        return(0);
    }
    if (ctx->BREAKFlg > 0) {
        if (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8)) {
            ctx->BREAKFlg = 1;
            return(1);
        }
        prn_nex(ctx, ctx->CmdBuf);
        return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_if(cmd)     Check for if/endif. Return 1 if continue, 0 if not.   */
/*                                        If error, return -1.              */

int check_if(TDAContext *ctx, char *cmd)    
{
    int err,row,col,ivflg;
    register char *p,*q;
    double res[2];

    err = -1;
    ivflg = 0;

    if (!strncmp(cmd,"if(",3)) {

        /** printf1("IF%s [level %d]\n",cmd + 2,IFLev + 1); **/

        if (ctx->IFLev >= 0) {
            if (ctx->IFLev >= MaxIFLEV - 1) {
                printf1(ctx, "Error: exceeded maximal if levels.\n\n");
                return(-1);
            }
            if (ctx->IFFLG[ctx->IFLev]) {
                ctx->IFFLG[++ctx->IFLev] = 1;     /* set false */
                return(0);
            }
        }
        ctx->IFLev++;

        /* check expression */

        p = cmd + 3;
        q = skip_expr(ctx, p);              
        if (*q != ')') {
            rep_err(ctx, cmd,0);
            goto CHKIFFin;
        }
        *q = '\0';          
        if (get_mexpr(ctx, p,&row,&col,&ivflg))
            goto CHKIFErr;

        if (row != 1 || col != 1) {
            printf1(ctx, "Error: if() needs a scalar expression.\n");
            goto CHKIFFin;
        }
        if (eval_mexpr(ctx, 0,row,col,res,0,res,NULL,NULL,-1,0.0))
            goto CHKIFErr;

        if (fabs(res[1]) < ctx->EPSI1)
            ctx->IFFLG[ctx->IFLev] = 1;
        else
            ctx->IFFLG[ctx->IFLev] = 0;

        printf2(ctx, "if%s) [level %d]\n\n",cmd + 2,ctx->IFLev);

        err = 0;
        goto CHKIFFin;
    }
    else if (!strcmp(cmd,"else")) {
        /**  printf1("ELSE [level %d]\n",IFLev); **/
        if (ctx->IFLev >= 0) {
            if (ctx->IFLev == 0 || ctx->IFFLG[ctx->IFLev - 1] == 0) {
                if (ctx->IFFLG[ctx->IFLev] == 1)
                    ctx->IFFLG[ctx->IFLev] = 0;
                else
                    ctx->IFFLG[ctx->IFLev] = 1;
                printf2(ctx, "else [level %d]\n\n",ctx->IFLev);
            }
            else
                ctx->IFFLG[ctx->IFLev] = 1;
            return(0);
        }
        else {
            printf1(ctx, "Error: unbalanced else.\n");
            return(-1);
        }
    }
    else if (!strcmp(cmd,"endif")) {
        /**  printf1("endif [level %d]\n",IFLev); **/
        if (ctx->IFLev >= 0) {
            ctx->IFLev--;
            if (ctx->IFLev < 0 || ctx->IFFLG[ctx->IFLev] == 0)
                printf2(ctx, "endif [level %d]\n\n",ctx->IFLev + 1);
            return(0);
        }
        else {
            printf1(ctx, "Error: unbalanced endif.\n");
            return(-1);
        }
    }
    if (ctx->IFLev < 0 || ctx->IFFLG[ctx->IFLev] == 0)
        return(1);
    return(0);

CHKIFFin:
    return(err);

CHKIFErr:
    printf1(ctx, "Error: can't evaluate if expression.\n");
    printf1(ctx, "Expression: %s\n",p);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  prn_nex(cmd)    print: not executed: cmd                                */

void prn_nex(TDAContext *ctx, char *cmd)    
{
    register int i;
    register char *p;

    printf2(ctx, "Not executed: ");
    p = cmd;
    for (i = 0; i < 60; ++i) {
        if (!*p)
            break;
        printf2(ctx, "%c",*p++);
    }
    newline(ctx);
}


