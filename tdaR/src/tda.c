/****************************************************************************/
/*  tda                                                                    */
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
#include "t_cmd.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_mat.h"
#include "t_var.h"
#include "t_eval.h"
#include "t_rep.h"
#include "t_gcmd.h"
#include "tda_context.h"

int tda_main(TDAContext *ctx, int argc, char *argv[]);
void prn_head(TDAContext *ctx);
void prn_end(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  MAIN. Start of main program                                             */

int tda_main(TDAContext *ctx, int argc, char *argv[])
{
    register int i,j;
    int fin,len,bflag;
    short *cflag;
    register char *p,*q;

    /* Remember where the executable lives, so help() can find tda.hlp when
       TDA is started from another directory. */
    ctx->ExeDir[0] = '\0';
    if (argc > 0 && argv[0]) {
        size_t n = strlen(argv[0]);
        while (n > 0 && argv[0][n - 1] != PSepC && argv[0][n - 1] != '/')
            --n;
        if (n > 0 && n <= PATHSIZE) {
            memcpy(ctx->ExeDir, argv[0], n);
            ctx->ExeDir[n] = '\0';
        }
    }

    prn_head(ctx);     /* print header message */

    /* if no arguments print message and exit */
    
    if (argc <= 1) {
        tda_out("Syntax: tda commands ...\n");
        exit(0);
    }

    machp(ctx, 0);       /* set basic machine parameters */

    /* allocate flags for command line arguments */

    if (!(cflag = (short *)calloc((size_t)(argc),sizeof(short)))) {
        p_err(ctx, -2,1);
        exit(0);
    }

    /* check for primary commands */

    ctx->MaxNV = ctx->MaxMat = 0;
    ctx->I_Flag = 0;
    for (i = 1; i < argc; ++i) {
        j = 1;
        if (!strcmp(argv[i],"ierr"))
            ctx->IERRFlg = 1;
        else if (!strcmp(argv[i],"i"))
            ctx->I_Flag = 1;
        else if (sscanf(argv[i],"maxnv=%d",&ctx->MaxNV) == 1) {
            if (ctx->MaxNV < 1 || ctx->MaxNV > 32000) {
                tda_out("Error: range of maxnv is 1 ... 32000.\n");
                exit(0);
            }
        }
        else if (sscanf(argv[i],"maxmat=%d",&ctx->MaxMat) == 1) {
            if (ctx->MaxMat < 1 || ctx->MaxMat > 32000) {
                tda_out("Error: range of maxmat is 1 ... 32000.\n");
                exit(0);
            }
        }
        else if (sscanf(argv[i],"cblen=%d",&ctx->CmdBufL) == 1) {
            ;
        }
        else if (sscanf(argv[i],"maxstd=%d",&ctx->MaxSTD) == 1 && ctx->MaxSTD > 0)
            ;
        else
            j = 0;

        cflag[i] = (short)(j); 
    }
    if (ctx->MaxNV <= 0)             /* check max number of variables */
        ctx->MaxNV = MaxNVDef;

    if (ctx->MaxMat <= 0)            /* check max number of matrices */
        ctx->MaxMat = MaxMatDef;

    if (ctx->MaxSTD < 1)             /* max stack size for derivatives */
        ctx->MaxSTD = MaxSTDDef;

    if (ctx->CmdBufL < 1000)  
        ctx->CmdBufL = CmdBufLen;

    /* command buffer */
    if (!(ctx->CmdBuf = (char *)calloc((size_t)(ctx->CmdBufL),sizeof(char)))) {
        p_err(ctx, -2,1);
        exit(0);
    }
    memrq(ctx, ctx->CmdBufL,sizeof(char));
         
    if (alloc_vmax(ctx)) {     /* allocate memory for variables */
        tda_out("Insufficient memory for variables.\n");
        exit(0);
    }
    if (alloc_mat(ctx)) {      /* allocate memory for matrices */
        tda_out("Insufficient memory for matrices.\n");
        exit(0);
    }

    tda_out("Current memory: %d bytes.\n\n",ctx->MxMReq);

    /* begin executing commands */

    for (i = 1; i < argc; ++i) {
        if (cflag[i] == 0) {
            if (t_exec(ctx, argv[i]))
                break;
        }
    }
    if (ctx->I_Flag)  
        tda_out("Starting interactive mode. Max length of command buffer: %d\n",ctx->CmdBufL - 10);

    ctx->SILENTFlg = bflag = 0;
    while (ctx->I_Flag) {    /* interactive mode */

CONT:
        fin = 0;
        p = ctx->CmdBuf;
        len = ctx->CmdBufL - 1;

        while (fin == 0) {

            if (len < 200) {
                p_err(ctx, -43,1);
                goto CONT;
            }
            tda_out(": ");
            fflush(stdin);

            if (fgets(p,len,stdin) == NULL || !*p)      
                continue;
 
            get_sline(ctx, p);
            if (!*p)
                continue;
 
            q = p;
            while (*q) {
                if (*q == '{')
                    bflag++;
                else if (*q == '}')
                    bflag--;

                if (*q == ';' && bflag <= 0) {
                    fin = 1;
                    *q = '\0';
                    break;
                }
                q++;
            }
            if (fin)
                break;
            len -= (int)(q - p);
            p = q;
        }
        if (!strncmp(ctx->CmdBuf,"quit",4) || !strncmp(ctx->CmdBuf,"exit",4)) 
            break;

        t_exec(ctx, ctx->CmdBuf);
    }
    i = 0;
    if (ctx->REPLev >= 0) {
        tda_out("Warning: program ends with incomplete loop.\n");
        rep_free(ctx);         /* free memory for repeat */
        i++;
    }
    if (ctx->BREAKFlg > 0) {
        tda_out("Warning: program ends with unfinished break.\n");
        i++;
    }
    clear_tda(ctx);

    if (i)
        prnchar(ctx, '-',LLEN,1);    
    prn_end(ctx);          /* end of program */
    free(cflag);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_head    print header                                                */

void prn_head(TDAContext *ctx)
{
    tda_out("TDA. Analysis of Transition Data (%3.1fq). ",TDA_Version);
#if TIME_ON
    prn_time(ctx, TDA_CONSOLE);
#else
    newline(ctx);         
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_end()       print end of program                                    */

void prn_end(TDAContext *ctx)
{
    tda_out("Current memory: %d bytes. Max memory used: %d bytes.\n",ctx->MemReq,ctx->MxMReq);
    tda_out("End of program. ");
#if TIME_ON
    prn_time(ctx, TDA_CONSOLE);
#else
    newline(ctx);           
#endif
}

/* The R package links these sources without a main(). */
#ifndef TDA_R_PACKAGE
int main(int argc, char *argv[]) {
    TDAContext *ctx = tda_context_new();
    tda_context_init(ctx);
    int r = tda_main(ctx, argc, argv);
    tda_context_free(ctx);
    return r;
}
#endif
