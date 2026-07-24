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

/*  functions in t_rep.c */

int t_exec(char *cmd);           
int t_repeat(int typ,int lev);
int t_endrepeat(int typ,int lev);
int rep_alloc(int opt,int lev,int typ); 
void rep_free(void);
int rep_scmd(char *cmd);
void rep_ptyp(int typ,int opt);
void rep_err(char *s,int opt);
int rep_expr(int lev,double *res);
int check_break(char *cmd);   
int check_if(char *cmd);   
void prn_nex(char *cmd);   

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

int REPLev = -1;            /* current level of repeat                      */
short REPTyp[MaxREP];       /* type of repeat loop                          */
                            /* 0 : while                                    */
                            /* 1 : repeat                                   */
short REPIdx[MaxREP];       /* index of matrix created by repeat            */

short REPSFlg[MaxREP];      /* 1 = save commands, 2 = execute commands      */
char **REPCmd[MaxREP];      /* pointer to commands                          */
short REPCmdA[MaxREP];      /* allocated ...                                */
int REPCCnt[MaxREP];        /* number of commands on level lev              */
int REPICNT[MaxREP];        /* number of iterations for simple repeat       */
char REPCond[MaxREP];       /* True/false result of first evaluation        */
char *REPExpr[MaxREP];      /* saves while expression                       */
short REPExprA[MaxREP];     /* allocated                                    */

int BREAKFlg = 0;           /* set if break command                         */
int IFLev = -1;             /* current if level                             */
char IFFLG[MaxIFLEV];       /* saves truth values of if expressions         */

/*--------------------------------------------------------------------------*/
/*  t_repeat(typ,lev)   Repeat command on level lev.                        */
/*                      typ 0 : while                                       */
/*                      typ 1 : repeat                                      */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int t_repeat(int typ,int lev)
{
    int err,m,idx,row,col,aflag,ivflg;
    register char *p,*q;
    char mname[VNLMax+1];
    double tmp;

    ivflg = aflag = 0;
    err = -1;             
    rep_ptyp(typ,0);
    printf2("[level %d]. Current memory: %d bytes.\n",++lev,MemReq);

    if (lev >= MaxREP) {
        printf1("Error: exceeded maximal level of loops.\n");
        return(-1);
    }
    if (rep_alloc(1,lev,typ))
        goto NREPFin;  

    aflag = 1;
    REPTyp[lev] = typ;

    if (typ == 0) {                     /* while */

        p = CmdBuf + 5;
        if (*p++ != '(' || *(q = skip_expr(p)) != ')') {
            rep_err(CmdBuf,0);
            goto NREPFin;
        }
        *q = '\0';          
        m = strlen(p) + 1;
        if (!(REPExpr[lev] = (char *)calloc(m,sizeof(char *)))) {
            p_err(-2,1);
            goto NREPFin;  
        }
        memrq(m,sizeof(char *));
        REPExprA[lev] = m;           
        strcpy(REPExpr[lev],p);

        if (rep_expr(lev,&tmp))
            goto NREPFin;
            
        printf2("expression: %s  Initial value: %lg",REPExpr[lev],tmp);
        if (tmp < EPSI1) {
            REPCond[lev] = 1;
            printf2(" [false]");
        }
        printf2("\n");
    }
    else {                              /* repeat */
        idx = -1;
        p = CmdBuf;
        if (strncmp(p,"repeat(n=",9)) {
            rep_err(CmdBuf,0);
            goto NREPFin;
        }           
        if ((p = n_mexpr(CmdBuf + 9,1,0,&row,&col,&ivflg,NULL)) == NULL || row * col != 1) {    
            rep_err(CmdBuf,0);
            goto NREPFin;
        }           
        m = (int)MX[0][1];
        if (m < 1) {
            rep_err(CmdBuf,0);
            goto NREPFin;
        }           
        if (*p == ',') {
            if (get_mname(++p,mname,1) == NULL)
                goto NREPFin;                     
            if (mat_getidx(mname,0) >= 0) {
                rep_err(CmdBuf,1);
                goto NREPFin;
            }           
            if (m_getmat(p,1,1,&idx,CmdBuf,1) == NULL) 
                goto NREPFin;
            MatVal[idx][1] = 1;
            REPIdx[lev] = idx + 1;
        }
        else if (*p != ')') {
            rep_err(CmdBuf,0);
            goto NREPFin;
        }
        REPICNT[lev] = m;
        if (SILENTFlg < 0) {
            printf2("Repeat for %d iterations.",m);
            if (idx >= 0)
                printf2(" Created matrix: %s",MatName[idx]);
            printf2("\n");
        }
    }
    err = 0;

NREPFin:
    mx_free();
    if (err && aflag)  
        rep_alloc(0,lev,0);
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  t_endrepeat(typ,lev)    Execute loop on level lev.                      */
/*                          typ 0 : while                                   */
/*                          typ 1 : repeat                                  */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int t_endrepeat(int typ,int lev)
{
    register int i;
    int err,icmd,iter = 0;
    double tmp;

    err = 0;             
    i = REPTyp[lev];
    if (i > 1)
        i = 1;

    if (lev < 0 || typ != i) {
        printf1("Error: unbalanced end");
        rep_ptyp(typ,1);
        printf1("\n");
        err = -1;
        goto EREPFin;
    }
    /**  
    for (i = 0; i < REPCCnt[lev]; ++i) {
        printf1("i=%d lev=%d ex: %s\n",i,lev,REPCmd[lev][i]);
    }
    **/  

    if (BREAKFlg)
        goto EREPFin1;

    /* execute current level */
   
    REPSFlg[lev] = 2;
    iter = 0;
    while (++iter > 0) {                

        /* check condition for break */

        if (REPTyp[lev] == 0) {                      /* while */
   
            if (rep_expr(lev,&tmp)) {
                err = -1;
                goto EREPFin1;
            }
            if (tmp < EPSI1)
                goto EREPFin1;
        }
        else {                                      /* repeat */
            REPICNT[lev] -= 1;
            if (REPICNT[lev] <= 0)
                goto EREPFin1;
            if (REPIdx[lev] > 0)
                MatVal[REPIdx[lev] - 1][1] += 1.0;
        }
        rep_ptyp(typ,0);
        if (SILENTFlg < 0) {
            printf2("[level %d, iteration %d]. Current memory: %d bytes.\n",
                                                     lev,iter,MemReq);
            prnchar('-',LLEN,1);      
        }
        for (icmd = 0; icmd < REPCCnt[lev] - 1; ++icmd) {

            strcpy(CmdBuf,REPCmd[lev][icmd]);

            /* don't save commands for lower levels */

            for (i = 0; i < lev; ++i)
                REPSFlg[i] = -iabs((int)REPSFlg[i]);

            err = t_exec(CmdBuf);

            for (i = 0; i < lev; ++i)    
                REPSFlg[i] = iabs((int)REPSFlg[i]);

            if (err)
                break;         
        }
        if (err || BREAKFlg)
            goto EREPFin1;
    }

EREPFin1:
    rep_alloc(0,lev,0);
     
    if (err == 0) {
        printf2("end");
        rep_ptyp(typ,0);
        printf2("[level %d, iteration %d]. Current memory: %d bytes.\n",
                                                         lev,iter,MemReq);
    }

EREPFin:
    BREAKFlg = 0;
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  rep_alloc(opt,lev,typ)      If opt != 0 allocate for level lev,         */
/*                              otherwise free.                             */
/*  Return 0 if OK, -1 if error.                                            */

int rep_alloc(int opt,int lev,int typ)
{
    register int i,l;

    /**
    printf1("repalloc opt=%d lev=%d typ=%d MemReq=%d REPLev=%d\n",opt,lev,typ,MemReq,REPLev); 
    **/ 
   
    if (opt == 0) {
        if (lev < REPLev) {
            printfe("ERROR1 in rep_alloc. opt=%d lev=%d replev=%d\n",opt,lev,REPLev);
            gerr_exit(201);
        }
        if (REPIdx[lev] > 0) {
            mat_alloc(REPIdx[lev] - 1,0);
            REPIdx[lev] = 0;
        }
        if (REPExprA[lev] > 0) {         
            free((char *)REPExpr[lev]);
            memrq(-REPExprA[lev],sizeof(char *));
            REPExprA[lev] = 0;           
        }
        if (REPCmdA[lev] > 0) {
            for (i = 0; i < REPCCnt[lev]; ++i) {
                l = strlen(REPCmd[lev][i]);
                free((char *)REPCmd[lev][i]);
                memrq(-l - 1,sizeof(char *));
            }
            free((char *)REPCmd[lev]);
            memrq(-REPCmdA[lev],sizeof(char **));
            REPCCnt[lev] = REPCmdA[lev] = 0;          
            REPTyp[lev] = REPSFlg[lev] = 0;
            REPLev--;
        }
        /**
        printf1("ret1 repalloc opt=%d lev=%d typ=%d MemReq=%d \n",opt,lev,typ,MemReq); 
        **/
        return(0);
    }
    if (lev != REPLev + 1) {
        printfe("ERROR2 in rep_alloc. opt=%d lev=%d replev=%d\n",opt,lev,REPLev);
        gerr_exit(202);
    }
    if (!(REPCmd[lev] = (char **)calloc(MaxREPCMD,sizeof(char **)))) {
        p_err(-2,1);
        return(-1);    
    }
    memrq(MaxREPCMD,sizeof(char **));
    REPCmdA[lev] = MaxREPCMD;
    REPTyp[lev] = typ;
    REPSFlg[lev] = 1;
    REPCond[lev] = 0;
    REPLev = lev;
    /**
    printf1("ret2 repalloc opt=%d lev=%d typ=%d MemReq=%d \n",opt,lev,typ,MemReq); 
    **/
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  rep_free()      Free all memory for repeat.                             */

void rep_free(void)
{
    while (REPLev >= 0)  
        rep_alloc(0,REPLev,0);
}

/*--------------------------------------------------------------------------*/
/*  rep_scmd(cmd)       Save command for repeat loop.                       */
/*  Return 0 if OK, -1 if error, 1 if command should be skipped.            */

int rep_scmd(char *cmd)
{
    int l,lev;
    static int skiplev = -1;
  
    if (REPCond[REPLev]) {
   
        if (skiplev < 0)  
            skiplev++;             
         
        if (!strncmp(cmd,"repeat",6) || !strncmp(cmd,"while(",6)) 
            skiplev++;

        else if (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8))                              
            skiplev--;

        if (skiplev >= 0)
            return(1);
    }
    for (lev = 0; lev <= REPLev; ++lev) {

        /********
        printf("SAVE lev=%d REPLev=%d REPSFlg[lev]=%d cnt=%d cmd: %s\n",lev,REPLev,REPSFlg[lev],REPCCnt[lev],cmd);
        *********/

        if (REPSFlg[lev] != 1)
            continue;

        /******
        if (lev == REPLev && (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8)))                   
            continue;
        ***/

        if (REPCCnt[lev] >= MaxREPCMD) {
            printf1("Error: exceeded max number of commands in loop.\n");
            return(-1);
        }
        l = strlen(cmd) + 1;
        if (!(REPCmd[lev][REPCCnt[lev]] = (char *)calloc(l,sizeof(char *)))) {
            p_err(-2,1);
            return(-1);    
        }
        memrq(l,sizeof(char *));
        strcpy(REPCmd[lev][REPCCnt[lev]],cmd);
        REPCCnt[lev] += 1;  
    }
    return(0);
}

/*--------------------------------------------------------------------------*/
/*  rep_ptyp(typ,opt)   print type of loop.                                 */

void rep_ptyp(int typ,int opt)
{
    if (typ == 0)
        printf2("while ");
    else 
        printf2("repeat ");
}

/*--------------------------------------------------------------------------*/
/*  rep_err(s,opt)  print error message.                                    */

void rep_err(char *s,int opt)
{
    printf1("Error: %s\n",s);
    if (opt == 1)
        printf1("Matrix already exists.\n");
}

/*--------------------------------------------------------------------------*/
/*  rep_expr(lev,res)   Evaluate while expression on level lev.             */
/*                      Return result in res. Return 0 if OK, -1 if error.  */

int rep_expr(int lev,double *res)
{
    int err,row,col,ivflg;
    double tmp[2];

    err = -1;
    ivflg = 0;
    if (get_mexpr(REPExpr[lev],&row,&col,&ivflg))
        goto REPExprErr;

    if (row != 1 || col != 1) {
        printf1("Error: while() needs a scalar expression.\n");
        goto REPExprFin;
    }
    if (eval_mexpr(0,row,col,tmp,0,tmp,NULL,NULL,-1,0.0))
        goto REPExprErr;

    *res = fabs(tmp[1]);
    err = 0;

REPExprFin:
    return(err);                 

REPExprErr:
    alloc_acx(0);
    printf1("Error: can't evaluate while expression (level %d).\n",lev);
    printf1("Expression: %s\n",REPExpr[lev]);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  check_break(cmd)    Check for break. Return 1 if continue, else 0.      */
/*                      -1 if error.                                        */

int check_break(char *cmd)    
{
    if (!strcmp(cmd,"break")) {
        if (REPLev < 0) {           /* must be in loop */
            printf1("Error: misplaced break.\n");
            return(-1);
        }
        BREAKFlg = 1;
        printf2("break\n\n");
        return(0);
    }
    if (BREAKFlg > 0) {
        if (!strncmp(cmd,"endrepeat",9) || !strncmp(cmd,"endwhile",8)) {
            BREAKFlg = 1;
            return(1);
        }
        prn_nex(CmdBuf);
        return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_if(cmd)     Check for if/endif. Return 1 if continue, 0 if not.   */
/*                                        If error, return -1.              */

int check_if(char *cmd)    
{
    int err,row,col,ivflg;
    register char *p,*q;
    double res[2];

    err = -1;
    ivflg = 0;

    if (!strncmp(cmd,"if(",3)) {

        /** printf1("IF%s [level %d]\n",cmd + 2,IFLev + 1); **/

        if (IFLev >= 0) {
            if (IFLev >= MaxIFLEV - 1) {
                printf1("Error: exceeded maximal if levels.\n\n");
                return(-1);
            }
            if (IFFLG[IFLev]) {
                IFFLG[++IFLev] = 1;     /* set false */
                return(0);
            }
        }
        IFLev++;

        /* check expression */

        p = cmd + 3;
        q = skip_expr(p);              
        if (*q != ')') {
            rep_err(cmd,0);
            goto CHKIFFin;
        }
        *q = '\0';          
        if (get_mexpr(p,&row,&col,&ivflg))
            goto CHKIFErr;

        if (row != 1 || col != 1) {
            printf1("Error: if() needs a scalar expression.\n");
            goto CHKIFFin;
        }
        if (eval_mexpr(0,row,col,res,0,res,NULL,NULL,-1,0.0))
            goto CHKIFErr;

        if (fabs(res[1]) < EPSI1)
            IFFLG[IFLev] = 1;
        else
            IFFLG[IFLev] = 0;

        printf2("if%s) [level %d]\n\n",cmd + 2,IFLev);

        err = 0;
        goto CHKIFFin;
    }
    else if (!strcmp(cmd,"else")) {
        /**  printf1("ELSE [level %d]\n",IFLev); **/
        if (IFLev >= 0) {
            if (IFLev == 0 || IFFLG[IFLev - 1] == 0) {
                if (IFFLG[IFLev] == 1)
                    IFFLG[IFLev] = 0;
                else
                    IFFLG[IFLev] = 1;
                printf2("else [level %d]\n\n",IFLev);
            }
            else
                IFFLG[IFLev] = 1;
            return(0);
        }
        else {
            printf1("Error: unbalanced else.\n");
            return(-1);
        }
    }
    else if (!strcmp(cmd,"endif")) {
        /**  printf1("endif [level %d]\n",IFLev); **/
        if (IFLev >= 0) {
            IFLev--;
            if (IFLev < 0 || IFFLG[IFLev] == 0)
                printf2("endif [level %d]\n\n",IFLev + 1);
            return(0);
        }
        else {
            printf1("Error: unbalanced endif.\n");
            return(-1);
        }
    }
    if (IFLev < 0 || IFFLG[IFLev] == 0)
        return(1);
    return(0);

CHKIFFin:
    return(err);

CHKIFErr:
    printf1("Error: can't evaluate if expression.\n");
    printf1("Expression: %s\n",p);
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  prn_nex(cmd)    print: not executed: cmd                                */

void prn_nex(char *cmd)    
{
    register int i;
    register char *p;

    printf2("Not executed: ");
    p = cmd;
    for (i = 0; i < 60; ++i) {
        if (!*p)
            break;
        printf2("%c",*p++);
    }
    newline();
}


