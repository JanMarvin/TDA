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

int main(int argc, char *argv[]);
void prn_head(void);
void prn_end(void);
int I_Flag = 0;             /* set in interactive mode                      */

/* ------------------------------------------------------------------------ */
/*  MAIN. Start of main program                                             */

int main(int argc, char *argv[])
{
    register int i,j;
    int fin,len,bflag;
    short *cflag;
    register char *p,*q;

    prn_head();     /* print header message */

    /* if no arguments print message and exit */
    
    if (argc <= 1) {
        printf("Syntax: tda commands ...\n");
        exit(0);
    }

    machp(0);       /* set basic machine parameters */

    /* allocate flags for command line arguments */

    if (!(cflag = (short *)calloc(argc,sizeof(short)))) {
        p_err(-2,1);
        exit(0);
    }

    /* check for primary commands */

    MaxNV = MaxMat = 0;
    I_Flag = 0;
    for (i = 1; i < argc; ++i) {
        j = 1;
        if (!strcmp(argv[i],"ierr"))
            IERRFlg = 1;
        else if (!strcmp(argv[i],"i"))
            I_Flag = 1;
        else if (sscanf(argv[i],"maxnv=%d",&MaxNV) == 1) {
            if (MaxNV < 1 || MaxNV > 32000) {
                printf("Error: range of maxnv is 1 ... 32000.\n");
                exit(0);
            }
        }
        else if (sscanf(argv[i],"maxmat=%d",&MaxMat) == 1) {
            if (MaxMat < 1 || MaxMat > 32000) {
                printf("Error: range of maxmat is 1 ... 32000.\n");
                exit(0);
            }
        }
        else if (sscanf(argv[i],"cblen=%d",&CmdBufL) == 1) {
            ;
        }
        else if (sscanf(argv[i],"maxstd=%d",&MaxSTD) == 1 && MaxSTD > 0)
            ;
        else
            j = 0;

        cflag[i] = j; 
    }
    if (MaxNV <= 0)             /* check max number of variables */
        MaxNV = MaxNVDef;

    if (MaxMat <= 0)            /* check max number of matrices */
        MaxMat = MaxMatDef;

    if (MaxSTD < 1)             /* max stack size for derivatives */
        MaxSTD = MaxSTDDef;

    if (CmdBufL < 1000)  
        CmdBufL = CmdBufLen;

    /* command buffer */
    if (!(CmdBuf = (char *)calloc(CmdBufL,sizeof(char)))) {
        p_err(-2,1);
        exit(0);
    }
    memrq(CmdBufL,sizeof(char));
         
    if (alloc_vmax()) {     /* allocate memory for variables */
        printf("Insufficient memory for variables.\n");
        exit(0);
    }
    if (alloc_mat()) {      /* allocate memory for matrices */
        printf("Insufficient memory for matrices.\n");
        exit(0);
    }

    printf("Current memory: %d bytes.\n\n",MxMReq);

    /* begin executing commands */

    for (i = 1; i < argc; ++i) {
        if (cflag[i] == 0) {
            if (t_exec(argv[i]))
                break;
        }
    }
    if (I_Flag)  
        printf("Starting interactive mode. Max length of command buffer: %d\n",CmdBufL - 10);

    SILENTFlg = bflag = 0;
    while (I_Flag) {    /* interactive mode */

CONT:
        fin = 0;
        p = CmdBuf;
        len = CmdBufL - 1;

        while (fin == 0) {

            if (len < 200) {
                p_err(-43,1);
                goto CONT;
            }
            printf(": ");
            fflush(stdin);

            if (fgets(p,len,stdin) == NULL || !*p)      
                continue;
 
            get_sline(p);
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
        if (!strncmp(CmdBuf,"quit",4) || !strncmp(CmdBuf,"exit",4)) 
            break;

        t_exec(CmdBuf);
    }
    i = 0;
    if (REPLev >= 0) {
        printf("Warning: program ends with incomplete loop.\n");
        rep_free();         /* free memory for repeat */
        i++;
    }
    if (BREAKFlg > 0) {
        printf("Warning: program ends with unfinished break.\n");
        i++;
    }
    clear_tda();

    if (i)
        prnchar('-',LLEN,1);    
    prn_end();          /* end of program */
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prn_head    print header                                                */

void prn_head(void)
{
    printf("TDA. Analysis of Transition Data (%3.1fp). ",TDA_Version);
#if TIME_ON
    prn_time(stdout);
#else
    newline();         
#endif
}

/* ------------------------------------------------------------------------ */
/*  prn_end()       print end of program                                    */

void prn_end(void)
{
    printf("Current memory: %d bytes. Max memory used: %d bytes.\n",MemReq,MxMReq);
    printf("End of program. ");
#if TIME_ON
    prn_time(stdout);
#else
    newline();           
#endif
}

