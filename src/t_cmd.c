/****************************************************************************/
/*  t_cmd                                                                  */
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
#include "t_rep.h"
#include "t_var.h"
#include "t_gcmd.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_dmu.h"
#include "t_gdat.h"          
#include "t_mdat.h"
#include "t_pdat.h"
#include "t_rzoo.h"
#include "t_spss.h"
#include "t_stata.h"
#include "t_psf.h"
#include "t_lsreg.h"
#include "t_l1reg.h"
#include "t_qrmod.h"
#include "t_edat.h"
#include "t_ltb.h"
#include "t_ple.h"
#include "t_rate.h"
#include "t_seq.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_dstat.h"
#include "t_ineq.h"
#include "t_eval.h"
#include "t_eval3.h"
#include "t_eval4.h"
#include "t_dplot.h"
#include "t_gdd.h"
#include "t_graph.h"
#include "t_gio.h"
#include "t_plg.h"
#include "t_tree.h"
#include "t_glm.h"
#include "t_smo.h"
#include "t_nlreg.h"
#include "t_xplot.h"
#include "t_help.h"
#include "t_gf.h"
#include "t_brr.h"
#include "t_gdf.h"
#include "t_gg.h"
#include "t_com.h"
#include "t_prox.h"
#include "t_cl.h"
#include "t_mds.h"
#include "t_ass.h"
#include "t_imat.h"
#include "t_ireg.h"
#include "t_loglin.h"
#include "t_bfa.h"
#include "t_mes.h"
#include "t_dem.h"
#include "t_areg.h"
#include "t_plot.h"
#include "t_plot3.h"
#include "t_cplot.h"
#include "t_rod.h"
#include "t_map.h"
#include "t_intp.h"
#include "t_sd.h"
#include "t_sdr.h"
#include "t_sdx.h"
#include "t_e00.h"
#include "t_tri.h"
#include "t_xls.h"
#include "t_clip.h"

#if S_XWIN          /* only included if S_XWIN is set in tda.h */
#include "t_xwin.h"
#endif  

int get_ncmd(int fn);           
int t_exec(char *cmd);           
int t_execute(void);           
void cmd_err(int opt);        
int exec_macro(char *s);   
int check_macro(char *s,int l);
int new_macro(void);
void prn_merr(char *s,int opt);
int mlist(void);       
int mclear(int opt);    
void clear_tda(void);     

#define CFNMAX 10           /* max number of command file levels            */
int CFN = 0;                /* current level of command file                */
FILE *CFFILE[10];           /* command file pointer                         */
int IERRFlg = 0;            /* set by ierr command                          */

/* ------------------------------------------------------------------------ */
/*  clear_tda()     clear all allocated memory                              */

void clear_tda(void)            
{
    clear_a();          /* clear data etc */
    ps_close(0,1,0);    /* close PostScript output file */
    mat_free();         /* clear all matrices */
    mclear(0);          /* clear macros */
    init_hlp(0);        /* free help */
    dblock_alloc(0);    /* free dblock */
    svb_alloc(0);       /* free string buffer */
    gdd_free(0);        /* free gdd (if perm) */
    IERRFlg = 0;
}

/* ------------------------------------------------------------------------ */
/*  get_ncmd(fn)    get next command from command file CFFILE[fn].          */
/*                  if successful put command into CmdBuf and return 1,     */
/*                  if no more command return 0, if error return -1.        */

int get_ncmd(int fn)            
{
    register int fin,cnt;
    register char cflag,*p,*q;
    int bflag,hflag;

    hflag = bflag = 0;
    q = p = CmdBuf;
    fin = cnt = 0;

    while (fin == 0) {

        if (CmdBufL - cnt < 200) {
            p_err(-43,2);
            return(-1);
        }
        if (fgets(p,CmdBufL - cnt,CFFILE[fn]) == NULL)  
            break;

        q = skip_b(p);
        cflag = '\0';
        hflag = 0;

/**     if (*q != '#' && *q != '*' && *q != ';') {          **/
/**     if (*q != '#' && *q != ';') {   **/
        if (*q != '#') {

            while (*q && *q != '\n' && *q != LF && *q != CR) {
            /** if (*q == '"' || *q == '`' || (*q == '\'' && *(q - 1) != '>')) { **/
                if (*q == '"' || *q == '`') {
                    if (cflag == *q)
                        cflag = '\0';
                    else if (!cflag)
                        cflag = *q;
                    else {
                        *p++ = *q;
                        cnt++;
                    }
                }
                else if ((*q != ' ' && *q != '\t') || cflag || hflag) {             
                    if (*q == '{')
                        bflag++;
                    else if (*q == '}')
                        bflag--;

                    if (*q == ';' && !cflag) {
                        if (bflag <= 0)
                            break;
                    }
                    else if (*q == '#' && !cflag)  
                        break;

                    *p++ = *q;
                    cnt++;
                }
                if (cnt == 4 && !strncmp(CmdBuf,"help",4)) {
                    hflag = 1;
                    if (*(q + 1) == '=')    /* ## */
                        q++;
                }
                q++;
            }
            if (*q == ';')
                fin = 1;
            *p = '\0'; 
        }
    }
    if (cnt == 0)  
        return(0);
       
    p = CmdBuf + strlen(CmdBuf) - 1;
    while (p >= CmdBuf && *p && (*p == '\n' || *p == LF || *p == CR)) {
        *p-- = '\0';
        fin = 0;
    }
    if (fin == 0)
        return(-1);

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  t_exec()        execute cmd.                                            */
/*                  return 0 if OK, -1 if error.                            */

int t_exec(char *cmd)            
{
    int n,err = -1;

    if (strncmp(cmd,"cf=",3)) {
        strcpy(CmdBuf,cmd);
        err = t_execute();
    
        if (err == 1)               /* check macro */
            err = exec_macro(CmdBuf);

        if (err > 0)  
            cmd_err(0);             /* unknown command */

        else if (err < 0 && (SILENTFlg == 1 || SILENTFlg == 3)) {
            cmd_err(1);     
            if (IERRFlg == 0)
                SILENTFlg = 0;
        }
        if (err && IERRFlg)         /* ignore errors */
            err = 0;
        goto EXECFin;
    }
    if (CFN >= CFNMAX) {
        p_err(-37,1);       /* reached max level of command files */
        goto EXECFin;
    }   
    if (!(CFFILE[CFN] = fopen(cmd + 3,OPEN_RD))) {
        ps_err(-1,cmd + 3,1);
        goto EXECFin;
    }
    printf1("Reading command file: %s\n",cmd + 3);
    prnchar('=',LLEN,1);

    CFN++;
    while (CFN > 0) {

        n = get_ncmd(CFN - 1);
        if (n == 1) {
            if (t_exec(CmdBuf))  
                goto EXECFin;
        }
        else if (n < 0) {
            printf1("Command file ends with an incomplete command.\n");
            goto EXECFin;
        }
        else {
            CFN--;
            fclose(CFFILE[CFN]);
            break;  
        }
    }
    err = 0;

EXECFin:
    if (err) {
        while (CFN > 0) {
            CFN--;
            fclose(CFFILE[CFN]);
        }
        if (REPLev >= 0)
            rep_free();
        err = -1;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  t_execute()     execute command in CmdBuf.                              */
/*                  return 0 if OK, -1 if error, 1 if cmd not defined.      */

int t_execute(void)         
{
    int n,r,lflg,err;
    char *p = CmdBuf;

    lflg = 0;
    err = 1;

    if (REPLev >= 0) {              /* save command */
        r = rep_scmd(CmdBuf);       /* don't execute if r = 1 */
        if (r) {
            if (r == 1) {
                prn_nex(CmdBuf);
                return(0);
            }
            else   
                return(-1);
        }
    }
    r = check_if(CmdBuf);           /* check for if/endif */
    if (r != 1)  
        return(r);
       
    r = check_break(CmdBuf);        /* check for break */
    if (r != 1)  
        return(r);

    switch (CmdBuf[0]) {

        case 'a':   if (!strncmp(p,"arcd",4))           /* arc */
                        err = arcd();
                    else if (!strcmp(p,"arcc"))         /* arcc */
                        err = arcc();
                    else if (!strncmp(p,"arcvc",5))     /* arcvc */
                        err = arcvc();
                    else if (!strncmp(p,"arcv",4))      /* arcv */
                        err = arcv();
                    else if (!strncmp(p,"atab",4))      /* atab */
                        err = atab();
                    else if (!strncmp(p,"acl",3))       /* acl */
                        err = acl();

                    break;

        case 'b':   if (!strncmp(p,"brr",3))            /* brr */
                        err = brr();
                    else if (!strncmp(p,"bfa",3))       /* bfa */
                        err = bfa();
                    else if (!strncmp(p,"bfc",3))       /* bfc */
                        err = bfc();
                    else if (!strncmp(p,"becl",4))      /* becl */
                        err = becl();
                    break;

        case 'c':   if (!strncmp(p,"ccnt=",5))           /* ccnt */
                        err = dm_ccnt();
                    else if (!strncmp(p,"clearnl",7))   /* clearnl */
                        err = clearnl();
                    else if (!strncmp(p,"clear",5))     /* clear */
                        err = clear();
                    else if (!strncmp(p,"cwt",3))       /* cwt */
                        err = cwt();
                    else if (!strncmp(p,"cov",3))       /* cov */
                        err = pcov(0);
                    else if (!strncmp(p,"corr",4))      /* corr */
                        err = pcov(1);
                    else if (!strncmp(p,"com",3))       /* com */
                        err = com();
                    else if (!strncmp(p,"conj",4))      /* conj */
                        err = conj();
                    else if (!strncmp(p,"cro",3))       /* cro */
                        err = cro();
                    else if (!strncmp(p,"clpyr",5))     /* clpyr */
                        err = clpyr();
                    else if (!strncmp(p,"clp",3))       /* clp */
                        err = clp();
                    else if (!strncmp(p,"clu",3))       /* clu */
                        err = clu();
                    break;

        case 'd':   if (!strncmp(p,"dump",4))           /* dump */
                        err = dm_dump();
                    else if (!strncmp(p,"dsplit",6))    /* dsplit */
                        err = dm_dsplit();
                    else if (!strcmp(p,"data1"))        /* data1 */
                        err = data(1);
                    else if (!strcmp(p,"data"))         /* data */
                        err = data(0);
                    else if (!strncmp(p,"dstat",5))     /* dstat */
                        err = dstat();  
                    else if (!strncmp(p,"dplot",5))     /* dplot */
                        err = dplot();
                    else if (!strncmp(p,"dblock",6))    /* dblock */
                        err = dblock();
                    else if (!strncmp(p,"dma",3))       /* dma */
                        err = dma();
                    else if (!strncmp(p,"dmet1",5))     /* dmet1 */
                        err = dmet1();
                    else if (!strncmp(p,"dmet",4))      /* dmet */
                        err = dmet();
                    else if (!strncmp(p,"dple",4))      /* dple */
                        err = dple();
                    else if (!strncmp(p,"diple",5))     /* diple */
                        err = diple();
                    else if (!strncmp(p,"dltb",4))      /* dltb */
                        err = dltb();
                    break;

        case 'e':   if (!strncmp(p,"edef",4))           /* edef */
                        err = edef();
                    else if (!strcmp(p,"endwhile")) {   /* endwhile */
                        err = t_endrepeat(0,REPLev);
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strcmp(p,"endrepeat")) {  /* endrepeat */
                        err = t_endrepeat(1,REPLev);
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"epdat",5))     /* epdat */
                        err = epdat();
                    else if (!strncmp(p,"epsdat",6))    /* epsdat */
                        err = epsdat();
                    else if (!strncmp(p,"evalfi",6))    /* evalfi */
                        err = evalf(1);
                    else if (!strncmp(p,"evalf",5))     /* evalf */
                        err = evalf(0);
                    else if (!strncmp(p,"esort",5))     /* esort */
                        err = esort();
                    else if (!strncmp(p,"eskip",5))     /* eskip */
                        err = eskip();
                    else if (!strncmp(p,"eselect",7))   /* eselect */
                        err = eselect();
                    else if (!strncmp(p,"emerge",6))    /* emerge */
                        err = emerge();
                    else if (!strncmp(p,"ejoin",5))     /* ejoin */
                        err = ejoin();
                    else if (!strncmp(p,"expm",4))      /* expm */
                        err = e_expm();
                    else if (!strncmp(p,"etest",5))     /* etest */
                        err = etest();
                    break;

        case 'f':   if (!strncmp(p,"fml",3))            /* fml */
                        err = f_min(0);
                    else if (!strncmp(p,"fmin",4))      /* fmin */
                        err = f_min(1);
                    else if (!strncmp(p,"freg",4))      /* freg */
                        err = f_min(2);
                    else if (!strncmp(p,"frml",4))      /* frml */
                        err = f_min(3);
                    else if (!strncmp(p,"freq1",5))     /* freq1 */
                        err = mfreq(1);
                    else if (!strncmp(p,"freq2",5))     /* freq2 */
                        err = mfreq(2);
                    else if (!strncmp(p,"freq",4))      /* freq */
                        err = mfreq(0);
                    break;

        case 'g':   if (!strncmp(p,"glm",3))            /* glm */
                        err = glm();
                    else if (!strncmp(p,"gmin",4))      /* gmin */
                        err = gmin();

                    else if (!strncmp(p,"gdd",3))       /* gdd */
                        err = gdd();
                    else if (!strncmp(p,"gcd",3))       /* gcd */
                        err = gcd();
                    else if (!strncmp(p,"gnc",3))       /* gnc */
                        err = gnc();
                    else if (!strncmp(p,"gcon",4))      /* gcon */
                        err = gcon();
                    else if (!strncmp(p,"gdcyc",5))     /* gdcyc */
                        err = gdcyc();
                    else if (!strncmp(p,"gcyc",4))      /* gcyc */
                        err = gcyc();
                    else if (!strncmp(p,"gdp",3))       /* gdp */
                        err = gdp();
                    else if (!strncmp(p,"gda",3))       /* gda */
                        err = gda();
                    else if (!strncmp(p,"gdu",3))       /* gdu */
                        err = gdu();
                    else if (!strncmp(p,"gde",3))       /* gde */
                        err = gde();
                    else if (!strncmp(p,"gdcset",6))    /* gdcset */
                        err = gdcset();
                    else if (!strncmp(p,"gsym",4))      /* gdot */
                        err = gsym();
                    else if (!strncmp(p,"gdot",4))      /* gdot */
                        err = gdot();
                    else if (!strncmp(p,"gni",3))       /* gni */
                        err = gni();
                    else if (!strncmp(p,"gdln",4))      /* gdln */
                        err = gdln();
                    else if (!strncmp(p,"gio",3))       /* gio */
                        err = gio();
                    else if (!strncmp(p,"gdcon",5))     /* gdcon */
                        err = gdcon();
                    else if (!strncmp(p,"gfcf",4))      /* gfcf */
                        err = gfcf();
                    else if (!strncmp(p,"gfc",3))       /* gfc */
                        err = gfc();
                    else if (!strncmp(p,"gbcf",4))      /* gbcf */
                        err = gbcf();
                    else if (!strncmp(p,"gep",3))       /* gep */
                        err = gep();
                    else if (!strncmp(p,"gtcl",4))      /* gtcl */
                        err = gtcl();
                    else if (!strncmp(p,"gflow",5))     /* gflow */
                        err = gflow();
                    else if (!strncmp(p,"gcliq",5))     /* gcliq */
                        err = gcliq();
                    else if (!strncmp(p,"ggcliq",6))    /* ggcliq */
                        err = ggcliq();
                    else if (!strncmp(p,"gev",3))       /* gev */
                        err = gev();
                    else if (!strncmp(p,"gsort",5))     /* gsort */
                        err = gsort();
                    else if (!strncmp(p,"gst",3))       /* gst */
                        err = gst();
                    else if (!strncmp(p,"gnst",4))      /* gnst */
                        err = gnst();
                    else if (!strncmp(p,"gmst",4))      /* gmst */
                        err = gmst();
                    else if (!strncmp(p,"gcut",4))      /* gcut */
                        err = gcut();
                    else if (!strncmp(p,"gsp",3))       /* gsp */
                        err = gsp();
                    else if (!strncmp(p,"giset",5))     /* giset */
                        err = giset();
                    else if (!strncmp(p,"gpro",4))      /* gpro */
                        err = gpro();
                    else if (!strncmp(p,"gcset",5))     /* gcset */
                        err = gcset();
                    else if (!strncmp(p,"gcni",4))      /* gcni */
                        err = gcni();
                    else if (!strncmp(p,"gap",3))       /* gap */
                        err = gap();
                    else if (!strncmp(p,"gqap",4))      /* gqap */
                        err = gqap();
                    else if (!strncmp(p,"gloc",4))      /* gloc */
                        err = gloc();
                    else if (!strncmp(p,"ghd1",4))      /* ghd1 */
                        err = ghd1();
                    else if (!strncmp(p,"ghd",3))       /* ghd */
                        err = ghd();

                    else if (!strncmp(p,"gdf",3))       /* gdf */
                        err = gdf();
                    else if (!strncmp(p,"gtopo",5))     /* gtopo */
                        err = gtopo();                

                    break;

        case 'h':   if (!strncmp(p,"help=",5))          /* help */
                        err = help(p + 5);
                    else if (!strncmp(p,"help",4))  
                        err = help(p + 4);
                    else if (!strncmp(p,"hclsp",5))     /* hclsp */
                        err = hclsp();
                    else if (!strncmp(p,"hcls",4))      /* hcls */
                        err = hcls();
                    else if (!strncmp(p,"hcld",4))      /* hcld */
                        err = hcld();
                    break;

        case 'i':   if (!strncmp(p,"intp",4))           /* intp */
                        err = intp();
                    else if (!strncmp(p,"int",3))       /* integration */
                        err = t_int();
                    else if (!strcmp(p,"ierr")) {       /* ierr */
                        IERRFlg = 1;
                        err = 0;            
                    }
                    else if (!strncmp(p,"ineq",4))      /* ineq */
                        err = ineq();
                    else if (!strncmp(p,"imean",5))     /* imean */
                        err = imean();
                    else if (!strncmp(p,"ivar1",5))     /* ivar1 */
                        err = ivar1();
                    else if (!strncmp(p,"ivar",4))      /* ivar */
                        err = ivar();
                    else if (!strncmp(p,"igini",5))     /* igini */
                        err = igini();
                    else if (!strncmp(p,"idf(",4))      /* idf */
                        err = idf();
                    else if (!strncmp(p,"iddf(",5))     /* iddf */
                        err = iddf();
                    else if (!strncmp(p,"ilsreg",6))    /* ilsreg */
                        err = ilsreg();
                    else if (!strncmp(p,"inpreg",6))    /* inpreg */
                        err = inpreg();
                    else if (!strncmp(p,"imreg",5))     /* imreg */
                        err = imreg();
                    else if (!strncmp(p,"ivls",4))      /* ivls */
                        err = ivls();
                    else if (!strncmp(p,"ivreg1",6))    /* ivreg1 */
                        err = ivreg1();
                    else if (!strncmp(p,"ivreg2",6))    /* ivreg2 */
                        err = ivreg2();
                    else if (!strncmp(p,"ivreg",5))     /* ivreg */
                        err = ivreg();
                    else if (!strncmp(p,"indep",5))     /* indep */
                        err = indep();
                    break;

        case 'l':   if (!strncmp(p,"lcnt",4))           /* lcnt */
                        err = dm_lcnt();
                    else if (!strncmp(p,"lsreg1",6))    /* lsreg1 */
                        err = lsreg1();
                    else if (!strncmp(p,"lsreg",5))     /* lsreg */
                        err = lsreg();
                    else if (!strncmp(p,"l1reg",5))     /* l1reg */
                        err = l1reg();
                    else if (!strncmp(p,"ltb",3))       /* life table */
                        err = ltb();
                    else if (!strncmp(p,"local(",6))    /* local command */
                        err = alloc_local(p + 6);
                    else if (!strncmp(p,"loglin",6))    /* loglin */
                        err = loglin();
                    break;

        case 'm':   if (!strcmp(p,"mem")) {             /* memory use */
                        err = mem_use();
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strcmp(p,"mach"))         /* machine parameters */
                        err = machp(1);
                    else if (!strncmp(p,"macrodef",8))  /* new macro */
                        err = new_macro();
                    else if (!strncmp(p,"macrolist",9)) /* list macros */
                        err = mlist();
                    else if (!strncmp(p,"macroclear",10))   /* clear macros */
                        err = mclear(1);
                    else if (!strncmp(p,"mdsx",4))      /* mdsx */
                        err = mdsx();
                    else if (!strncmp(p,"mdsm",4))      /* mdsm */
                        err = mdsm();
                    else if (!strncmp(p,"mdsn",4))      /* mdsn */
                        err = mdsn();
                    else if (!strncmp(p,"mdsc",4))      /* mdsc */
                        err = mdsc();
                    else if (!strncmp(p,"mdsr",4))      /* mdsr */
                        err = mdsr();
                    else if (!strncmp(p,"mproc",5))     /* mproc */
                        err = mproc();
                    else if (!strncmp(p,"mreg",4))      /* mreg */
                        err = mreg();
                    else if (!strncmp(p,"mparse",6))    /* parse expression */
                        err = mparse();
                    else {                              /* matrix commands */
                        err = t_mat();   
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    break;

        case 'n':   if (!strncmp(p,"nvar(",5))          /* new variables */
                        err = new_var();
                    else if (!strncmp(p,"nlist",5))     /* nlist */
                        err = nlist();
                    else if (!strncmp(p,"ndvar",5))     /* ndvar */
                        err = ndvar();
                    else if (!strncmp(p,"niset",5))     /* niset */
                        err = niset();
                    else if (!strncmp(p,"nlreg",5))     /* nlreg */
                        err = nlreg();
                    else if (!strncmp(p,"npreg",5))     /* npreg */
                        err = npreg();
                    else if (!strncmp(p,"nncl",4))      /* nncl */
                        err = nncl();
                    else if (!strncmp(p,"nmca",4))      /* nmca */
                        err = nmca();
                    break;

        case 'p':   if (!strncmp(p,"print(",6)) {       /* print text */
                        err = prn_txt();
                        if (err == 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"plotp3",6))    /* plot 3d polygon */
                        err = plotp3(0);
                    else if (!strncmp(p,"plot3",5))     /* plot 3d polygon */
                        err = plotp3(1);
                    else if (!strncmp(p,"plcurv3",7))   /* plot 3d curve */
                        err = plcurv3();
                    else if (!strncmp(p,"plsurf3d",8))  /* plot 3d surface */
                        err = plsurf3d();
                    else if (!strncmp(p,"plsurf3",7))   /* plot 3d surface */
                        err = plsurf3();
                    else if (!strncmp(p,"pltext3",7))   /* plot text */
                        err = pltext3();
                    else if (!strncmp(p,"plcirc3",7))   /* plot circle */
                        err = plcirc3();
                    else if (!strncmp(p,"plglob3",7))   /* plot globe */
                        err = plglob3();

                    else if (!strncmp(p,"parse=",6))    /* parse expression */
                        err = t_parse();
                    else if (!strncmp(p,"pdata",5))     /* pdata */
                        err = pdata();
                    else if (!strncmp(p,"pdatr",5))     /* pdatr */
                        err = pdatr();
                    else if (!strncmp(p,"pdatd",5))     /* pdatd */
                        err = pdatd();
                    else if (!strncmp(p,"ple",3))       /* product-limit */
                        err = ple();
                    else if (!strncmp(p,"pcyc",4))      /* pcyc */
                        err = pcyc();
                    else if (!strncmp(p,"ptree",5))     /* ptree */
                        err = ptree();
                    else if (!strncmp(p,"pltree",6))    /* pltree */
                        err = pltree();

                    else if (!strncmp(p,"psfile=",7))   /* psfile */
                        err = psfile();     
                    else if (!strncmp(p,"psetupg",7))   /* psetupg */
                        err = psetupg();     
                    else if (!strncmp(p,"psetup3",7))   /* psetup3 */
                        err = psetup3();     
                    else if (!strncmp(p,"psetup",6))    /* psetup */
                        err = psetup();     
                    else if (!strcmp(p,"psclose"))      /* psclose */
                        err = ps_close(1,1,0);     
                    else if (!strncmp(p,"plxa",4))      /* plot x axis */
                        err = pl_axis(0,0,0,1,0.2);
                    else if (!strncmp(p,"plya",4))      /* plot y axis */
                        err = pl_axis(1,0,0,1,0.2);
                    else if (!strncmp(p,"plframe",7))   /* plot frame */
                        err = pl_frame();
                    else if (!strncmp(p,"plabel",6))    /* plot label */
                        err = pl_label(0);
                    else if (!strncmp(p,"pxlabel",7))   /* plot label */
                        err = pl_label(1);
                    else if (!strncmp(p,"pylabel",7))   /* plot label */
                        err = pl_label(2);
                    else if (!strncmp(p,"plxgrid",7))   /* grid lines x axis */
                        err = ps_grid(0);
                    else if (!strncmp(p,"plygrid",7))   /* grid lines y axis */
                        err = ps_grid(1);
                    else if (!strncmp(p,"pltext",6))    /* plot text */
                        err = pl_text();
                    else if (!strncmp(p,"plrec",5))     /* plot rectangle */
                        err = pl_rec();

                    else if (!strncmp(p,"plotp",5))     /* plot polygon */
                        err = pl_plotp(0);
                    else if (!strncmp(p,"plotm",5))     /* plot polygon */
                        err = pl_plotm();
                    else if (!strncmp(p,"plotf",5))     /* plot function */
                        err = pl_plotf();
                    else if (!strncmp(p,"ploth",5))     /* plot histogram */
                        err = pl_ploth();
                    else if (!strncmp(p,"plotd",5))     /* plot density */
                        err = pl_plotd();
                    else if (!strncmp(p,"ploto",5))     /* plot circle */
                        err = pl_ploto(0);
                    else if (!strncmp(p,"plote",5))     /* plot ellipse */
                        err = pl_ploto(1);
                    else if (!strncmp(p,"plotk",5))     /* plot arc */
                        err = pl_plotk();
                    else if (!strncmp(p,"plotch",6))    /* convex hull */
                        err = pl_plotch();
                    else if (!strncmp(p,"plotcm",6))    /* contour plot */
                        err = pl_plotcm();
                    else if (!strncmp(p,"plotc",5))     /* contour plot */
                        err = pl_plotc();
                    else if (!strncmp(p,"plotr",5))     /* grey-scaled relief */
                        err = pl_plotr();
                    else if (!strncmp(p,"plotsp",6))    /* smoothing: Akima */
                        err = pl_plots(0);
                    else if (!strncmp(p,"plots",5))     /* smoothing: Akima */
                        err = pl_plots(1);
                    else if (!strncmp(p,"plot",4))      /* plot polygon */
                        err = pl_plotp(1);
                    else if (!strncmp(p,"plg",3))       /* plot graph */
                        err = pl_plg();

                    break;

        case 'q':   if (!strncmp(p,"quant",5))          /* quantiles */
                        err = quant();  
                    else if (!strncmp(p,"qreg",4))      /* qreg */
                        err = qreg();
                    break;

        case 'r':   if (!strncmp(p,"rsys",4))           /* read system file */
                        err = rd_sys();
                    else if (!strncmp(p,"repeat",6)) {  /* repeat */
                        err = t_repeat(1,REPLev);
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"rcorr",5))     /* rcorr */
                        err = pcov(2);
                    else if (!strncmp(p,"recode",6))    /* recode */
                        err = recode();
                    else if (!strncmp(p,"rspss1",6))    /* read spss sav file */
                        err = rd_spss1();
                    else if (!strncmp(p,"rspss",5))     /* read spss por file */
                        err = rd_spss();
                    else if (!strncmp(p,"rstata",6))    /* read stata file */
                        err = rd_stata();
                    else if (!strncmp(p,"rxls",4))      /* read xls file */
                        err = rxls();
                    else if (!strncmp(p,"rucinet",7))   /* read xls file */
                        err = rucinet();
                    else if (!strncmp(p,"rdbf",4))      /* read dbf file */
                        err = rdbf();
                    else if (!strncmp(p,"rplz",4))      /* read plz file */
                        err = rplz();
                    else if (!strncmp(p,"rcsv",4))      /* read csv file */
                        err = rcsv();
/*                  else if (!strncmp(p,"rdxf",4))  */  /* read dxf file */
/*                      err = rdxf();       */         

                    else if (!strncmp(p,"rate",4))      /* rate models */
                        err = rate();
                    else if (!strncmp(p,"rsys",4))      /* read system file */
                        err = rd_sys();
                    else if (!strncmp(p,"repsel",6))    /* repsel */
                        err = repsel();
                    else if (!strncmp(p,"range",5))     /* range */
                        err = range();
                    else if (!strncmp(p,"rmod",4))      /* rmod */
                        err = rmod();
                    else if (!strncmp(p,"rod",3))       /* rod */
                        err = rod();
                    else if (!strncmp(p,"rap",3))       /* rap */
                        err = rap();
                    else if (!strncmp(p,"rfit1",5))      /* rfit1 */
                        err = rfit1();
                    else if (!strncmp(p,"rfit",4))      /* rfit */
                        err = rfit();

                    break;

        case 's':   if (sscanf(p,"silent=%d",&n) == 1) {
                        SILENTFlg = imax(imin(n,3),-1);
                        lflg = 1;
                        err = 0;
                    }
                    else if (!strncmp(p,"silent",6)) {
                        printf1("Current value of silent: %d\n",SILENTFlg);
                        err = 0;
                    }
                    else if (!strncmp(p,"sdnvar",6))    /* sdnvar */
                        err = sdnvar();
                    else if (!strncmp(p,"sdpdata",7))   /* sdpdata */
                        err = sdpdata();
                    else if (!strncmp(p,"sdplot31",8))  /* sdplot31 */
                        err = sdplot31();
                    else if (!strncmp(p,"sdplot32",8))  /* sdplot32 */
                        err = sdplot32();
                    else if (!strncmp(p,"sdplot33",8))  /* sdplot33 */
                        err = sdplot33();
                    else if (!strncmp(p,"sdplot",6))    /* sdplot */
                        err = sdplot();
                    else if (!strncmp(p,"sdrel",5))     /* sdrel */
                        err = sdrel();
                    else if (!strncmp(p,"sdnl",4))      /* sdnl */
                        err = sdnl();
                    else if (!strncmp(p,"sdppol",6))    /* sdppol */
                        err = sdppol();
                    else if (!strncmp(p,"sdlpol",6))    /* sdlpol */
                        err = sdlpol();
                    else if (!strncmp(p,"sdipol",6))    /* sdipol */
                        err = sdipol();
                    else if (!strncmp(p,"sdcpol",6))    /* sdcpol */
                        err = sdcpol();
                    else if (!strncmp(p,"sddcwp",6))    /* sddcwp */
                        err = sddcwp();
                    else if (!strncmp(p,"sdgen",5))     /* sdgen */
                        err = sdgen();
                    else if (!strncmp(p,"sdgshhs",7))   /* sdgshhs */
                        err = sdgshhs();
                    else if (!strncmp(p,"sdshp",5))     /* sdshp */
                        err = sdshp();                 
                    else if (!strncmp(p,"sde00",5))     /* sde00 */
                        err = sde00();                 
                    else if (!strncmp(p,"sdinf",5))     /* sdinf */
                        err = sdinf();                 
                    else if (!strncmp(p,"sdencl",6))    /* sdencl */
                        err = sdencl();                 
                    else if (!strncmp(p,"sdsel",5))     /* sdsel */
                        err = sdsel();                 
                    else if (!strncmp(p,"sdvd",4))      /* sdvd */
                        err = sdvd();
                    else if (!strncmp(p,"sdclip",6))    /* sdclip */
                        err = sdclip();
                    else if (!strncmp(p,"sdpgrat",7))   /* sdpgrat */
                        err = sdpgrat();
                    else if (!strncmp(p,"sdpgeo",6))    /* sdpgeo */
                        err = sdpgeo();
                    else if (!strncmp(p,"sdpmap",6))    /* sdpmap */
                        err = sdpmap();

                    else if (!strncmp(p,"seq",3))       /* sequence data */
                        err = t_seq();
                    else if (!strncmp(p,"segr",4))      /* segr */
                        err = segr();
                    else if (!strncmp(p,"sma",3))       /* sma */
                        err = sma();
                    else if (!strncmp(p,"smd",3))       /* smd */
                        err = smd();
                    else if (!strncmp(p,"spl",3))       /* spl */
                        err = spl();
                    else if (!strncmp(p,"scplot",6))    /* scplot */
                        err = scplot();
                    else if (!strncmp(p,"sddf(",5))     /* sddf */
                        err = sddf();
                    else if (!strncmp(p,"spmod",5))     /* spmod */
                        err = spmod();
                    else if (!strncmp(p,"sga",3))       /* sga */
                        err = sga();
                    else if (!strncmp(p,"subm",4))      /* subm */
                        err = subm();
                    else if (!strncmp(p,"scla",4))      /* scla */
                        err = scla();

                    break;

        case 't':   if (!strcmp(p,"time"))              /* current time */
                        err = prntime();
                    else if (!strncmp(p,"tsel=",5))     /* tsel  */
                        err = tsel();
                    else if (!strncmp(p,"tnet",4))      /* tnet */
                        err = tnet();
                    else if (!strncmp(p,"triang",6))    /* triang */
                        err = triang();
                    break;
        
        case 'u':   if (!strncmp(p,"uds",3))            /* uds */
                        err = uds();
                    else if (!strncmp(p,"unf",3))       /* unf */
                        err = unf();
                    else if (!strncmp(p,"ucl",3))       /* ucl */
                        err = ucl();
                    else
                        break;

        case 'w':   if (!strncmp(p,"wspss1",6))         /* write spss sav file */
                        err = wr_spss1();
                    else if (!strncmp(p,"wspss",5))     /* write spss por file */
                        err = wr_spss();
                    else if (!strncmp(p,"wstata",6))    /* write stata file */
                        err = wr_stata();
                    else if (!strncmp(p,"wsys",4))      /* write system file */
                        err = wr_sys();
                    else if (!strncmp(p,"while",5)) {   /* while */
                        err = t_repeat(0,REPLev);
                        if (err == 0 && SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    break;

        case 'x':   if (!strncmp(p,"xopen",5))          /* xopen */
                        err = xopen();     
                    else if (!strncmp(p,"xdelete",7))   /* xdelete */
                        err = xdelete();    
                    else if (!strcmp(p,"xlog1"))        /* xlog1 */
                        err = xlog1();      
                    else if (!strncmp(p,"xlog",4))      /* xlog */
                        err = xlog();     
                    else if (!strncmp(p,"xconh",5))     /* xconh */
                        err = xconh();     
                    else if (!strncmp(p,"xreg",4))      /* xreg */
                        err = xreg();     
                    else if (!strncmp(p,"xplotf",6))    /* xplotf */
                        err = xplotf();     
                    else if (!strncmp(p,"xplot",5))     /* xplot */
                        err = xplot();     
                    else if (!strncmp(p,"xf",2))        /* xf */
                        err = xfunc();     
                    else if (!strncmp(p,"xdens",5))     /* xdens */
                        err = xdens();     
#if S_XWIN          /* only included if S_XWIN is set in tda.h */
                    else if (!strncmp(p,"xshow",5))     /* open plot window */
                        err = x_show();
#endif
                    break;

        case 'z':   if (!strncmp(p,"zreg1",5))          /* zreg1 */
                        err = zreg1();
                    else if (!strncmp(p,"zreg",4))           /* zreg */
                        err = zreg();
                    break;

        case '$':   err = check_cmd(1);
                    if (err == 0) {
                        fflush(stdout);
                        if (system(CmdBuf+1))  
                            printf1("Note: shell returned error.\n");
                        else
                            printf1("\n> end of shell execution.\n");
                    }        
                    break;

        default:    lflg = 1;
                    break;     

    }
    if (err != 1) {
        if (REPLev >= 0 && err < 0)
            rep_free();
        if (lflg == 0)
            prnchar('-',LLEN,1);    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cmd_err(opt)   Give error message for command in CmdBuf.                */

void cmd_err(int opt)      
{
    register char *p = CmdBuf;

    if (strlen(p) < 1)
        return;

    if (opt == 0)
        printf1("Unknown command: ");
    else
        printf1("Error while executing command: ");

    while (*p && p < CmdBuf + 20)
        printf1("%c",*p++);
    if (*p || opt)
        printf1(" ...");
    printf1("\n");  
    prnchar('-',LLEN,1);    
}

/* ------------------------------------------------------------------------ */
/*  exec_macro(s)   Try to execute macro s.                                 */
/*                                                                          */
/*                  a) if s is known macro expand and execute.              */
/*                  b) if s is expand=mname ...                             */
/*                     show expanded macro on stdout but do not execute.    */  
/*                                                                          */
/*  Return: 1 if unknown, 0 if OK, -1 if error.                             */

int exec_macro(char *s)    
{
    register int i;
    int err,n,m,nlen,aflag,np,l,eflag;
    register char *p,*q,*q1,*qq;
    char c,cflag;
    
    MacroExLevel += 1;
    err = 1; eflag = np = 0;

    if (!strncmp(s,"expand:",7) || !strncmp(s,"expand=",7)) {
        s += 7;
        eflag = 1;
    }
    p = s;
    while (*p) { 
        if (*p == '=' || *p == '(')
            break;
        p++;
    }
    nlen = (int)(p - s);

    /* check for existing macro */

    if (*p == '=' || *p == '(' || !*p) {
        c = *p;
        *p = '\0';
        n = check_macro(s,nlen);
        if (n < 0)
            goto MACROFin;
        *p = c;
        q = s;

        if (SILENTFlg < 0) {
            printf1("Macro: %s\n",s);
            printf1("Execution level: %d\n",MacroExLevel);
            prnchar('-',LLEN,1);    
        }
        err = -1;

        /* get arguments into MacroArgs */

        qq = q = p;
        aflag = 0;
        if (c == '(') {
            q = skip_blev(p);
            if (q > p + 1) {
                if (*(q - 1) != ')') {
                    prn_merr(s,1);
                    goto MACROFin;
                }
                aflag = 1;
                qq = q;
            }
        }
        np = 0;
        if (*q == '=') {
            q++;
            while (*q) {
                if (np >= MaxMacro) {                       
                    prn_merr(s,0);
                    goto MACROFin;
                }
                if (*q == ',') {
                    MacroArgs[np] = NULL;
                }
                else {
                    q1 = skip_nc(q);
                    l = (int)(q1 - q);
                    if (!(MacroArgs[np] = (char *)calloc(l + 1,sizeof(char)))) {
                        p_err(-2,1);
                        goto MACROFin;
                    }
                    memrq(l + 1,sizeof(char));
                    strncpy(MacroArgs[np],q,l);
                    *(MacroArgs[np] + l) = '\0';
                    q = q1;
                }
                np++;
                if (!*q)
                    break;
                q++;
            }      
        }

        if (aflag) {
            while (++p < qq && *p) {
                if (np >= MaxMacro) {                       
                    prn_merr(s,0);
                    goto MACROFin;
                }
                if (*p == ',')
                    MacroArgs[np] = NULL;
                else {
                    q = skip_nc(p);
                    l = (int)(q - p);
                    if (l == 0)
                        MacroArgs[np] = NULL;
                    else {
                        if (!(MacroArgs[np] = (char *)calloc(l + 1,sizeof(char)))) {
                            p_err(-2,1);
                            goto MACROFin;
                        }
                        memrq(l + 1,sizeof(char));
                        strncpy(MacroArgs[np],p,l);
                        *(MacroArgs[np] + l) = '\0';
                    }
                    p = q;
                }
                np++;
            }
        }
        if (np > MacroNP[n])  
            printf1("Warning: macro %s has %d argument(s).\n",MacroName[n],MacroNP[n]);

        /* insert arguments into macro definition and call exec()  */
        /* or show on stdout if eflag = 1 */

        CurMacroNum = n;

        if (eflag)  
            printf1("Expansion of macro: %s\n",s); /*  MacroName[n]);   */
  
        p = MacroDef[n];
        while (*p) {
            q = CmdBuf;
            l = 0;
            cflag = '\0';
            while (*p && *p != ';') {
                if (*p == '`') {
                    if (cflag == *p)
                        cflag = '\0';
                    else 
                        cflag = *p;
                    p++;
                }
                else if (*p == '$' && !cflag) {
                    if (sscanf(p,"$%d",&m) == 1 && m >= 1 && m <= MaxMacro) {
                        p = skip_int(p + 1);
                        m--;
                        if (MacroArgs[m] != NULL) {
                            nlen = strlen(MacroArgs[m]);
                            if (l + nlen > CmdBufL - 10) {
                                p_err(-43,1);
                                goto MACROFin;
                            }     
                            strcpy(q,MacroArgs[m]);
                            q += nlen;
                            l += nlen;
                        }
                    }
                    else {
                        *q++ = *p++;
                        l++;
                    }
                }   
                else {
                    if (l + 1 > CmdBufL - 10) {
                        p_err(-43,1);  
                        goto MACROFin;
                    }     
                    *q++ = *p++;
                    l++;
                }
            }
            *q = '\0';

            if (eflag)  
                printf1("> %s\n",CmdBuf);
            else {
                err = t_exec(CmdBuf);
                if (err)
                    goto MACROFin;
            }
            while (*p == ';')
                p++;
        }
        if (eflag)
            prnchar('-',LLEN,1);    

        else if (SILENTFlg < 0) {
            printf1("End of macro: %s\n",MacroName[n]);
            printf1("New execution level: %d\n",MacroExLevel - 1);
            prnchar('-',LLEN,1);    
        }
        err = 0;
    }

MACROFin:
    for (i = 0; i < np; ++i) {
        if (MacroArgs[i] != NULL) {
            memrq(-strlen(MacroArgs[i]) - 1,sizeof(char));
            free(MacroArgs[i]);
            MacroArgs[i] = NULL;
        }
    }
    free_local();           /* free local matrices */
    MacroExLevel -= 1;
    CurMacroNum = -1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_macro(s,l)  Return number of macro if s is a currently existing   */
/*                    macro name, otherwise return -1.                      */
/*                                                                          */

int check_macro(char *s,int l)
{
    register int i;

    for (i = 0; i < MacroN; ++i) {
        if (!strncmp(s,MacroName[i],l) && strlen(MacroName[i]) == l)
            return(i);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  new_macro()         Get new macro.                                      */
/*                      Syntax: macrodef(name)={...};                       */
/*  ##                  Return 0 if OK, -1 if error.                        */
/*                                                                          */

int new_macro(void)
{
    register int i;
    int err,l,n,m,nlen;
    register char *s,*p,*q;
    short marg[MaxMacro + 1];
    static int first = 1;

    if (first) {
        for (i = 0; i < MaxMacro; ++i) 
            MacroArgs[i] = NULL;
        first = 0;
    }

    err = -1;
    if (check_cmd(2))
        return(-1);

    printf1("New macro definition. Current memory: %d bytes.\n",MemReq);

    if (MacroN >= MaxMacro) {
        printf1("Error: exceeded max number of macro definitions.\n");
        goto GMACROFin;
    }
    s = CmdBuf + 8;
    if (*s++ != '(' || *(p = skip_nc(s)) != ')' || (nlen = (int)(p - s)) < 1     
                                        || *++p != '=' || *++p != '{') {
        prn_merr(CmdBuf,2);        
        goto GMACROFin;
    }
    if (check_macro(s,nlen) >= 0) {
        *(s + nlen) = '\0';
        printf1("Error: macro name %s already used.\n",s);
        goto GMACROFin;
    }
    q = skip_blev(p);
    if (*(q - 1) != '}') {
        prn_merr(CmdBuf,2);
        goto GMACROFin;
    }
    l = strlen(++p);
    if (l < 2) {
        prn_merr(CmdBuf,2);
        goto GMACROFin;
    }
    *--q = '\0';

    for (i = 1; i <= MaxMacro; ++i)
        marg[i] = 0;
    q = p;
    n = 0;

    while (*q) {
        if (sscanf(q,"$%d",&m) == 1 && m >= 0) {
            if (m < 1 || m > MaxMacro) {
                prn_merr(CmdBuf,2);
                printf1("Range of arguments is: $1 - $%d\n",MaxMacro);
                goto GMACROFin;
            }
            marg[m] = 1;
            if (n < m)
                n = m;
            q++;
        }
        q++;
    }
    m = 0;
    for (i = 1; i <= n; ++i) {
        if (marg[i] == 0) {
            prn_merr(CmdBuf,2);
            printf1("Argument numbers should be contiguous.\n");
            goto GMACROFin;
        }
    }
    if (!(MacroName[MacroN] = (char *)calloc(nlen + 1,sizeof(char)))) {
        p_err(-2,1);
        goto GMACROFin;
    }
    if (!(MacroDef[MacroN] = (char *)calloc(l,sizeof(char)))) {
        free(MacroName[MacroN]);
        p_err(-2,1);
        goto GMACROFin;
    }
    memrq(l + nlen + 1,sizeof(char));
    strncpy(MacroName[MacroN],s,nlen);
    strncpy(MacroDef[MacroN],p,l-1);

    MacroNP[MacroN] = n;
    if (MacroNLen < nlen)
        MacroNLen = nlen;

    printf1("Added new macro: %s\n",MacroName[MacroN]);
    MacroN++;
    err = 0;

GMACROFin:
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  prn_merr(s,opt)     print macro error.                                  */

void prn_merr(char *s,int opt)
{
    register int i;
    register char *p;

    if (opt == 0) {
        printf1("Error: exceeded maximal number of macro arguments.\n");
        return;
    }
    else if (opt == 1)
        printf1("Syntax error: ");
    else if (opt == 2)
        printf1("Error in macro definition: ");

    p = s;
    for (i = 0; i < 30; ++i) {
        if (!*p)
            break;
        printf1("%c",*p++);
    }   
    if (*p)
        printf1(" ...\n");
    newline();
}

/* ------------------------------------------------------------------------ */
/*  mlist()         List currently defined macros: macrolist                */

int mlist(void)        
{
    register int i,j,l;
    register char *p;
    char buf[20];

    printf1("Currently defined macros: ");

    if (MacroN == 0) {
        printf1("none\n");
        return(0);
    }
    newline();

    for (i = 0; i < MacroN; ++i) {
        printf1("%3d  %s",i + 1,MacroName[i]);
        l = strlen(MacroName[i]);      
        if (MacroNP[i] > 0) {
            sprintf(buf,"%d",MacroNP[i]);
            printf1("($%s)",buf);
            l += strlen(buf) + 3;
        }
        prnchar(' ',MacroNLen + 5 - l, 0);
        printf1(" [");
        p = MacroDef[i];
        for (j = 0; j < 50; ++j) {
            if (!*p)
                break;
            printf1("%c",*p++);
        }   
        if (*p)
            printf1(" ...");
        printf1("]\n");
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mclear()        Clear currently defined macros: macroclear;             */

int mclear(int opt)        
{
    register int i,l;

    if (MacroN == 0)
        return(0);

    for (i = 0; i < MacroN; ++i) {
        free(MacroName[i]);
        free(MacroDef[i]);
        l = strlen(MacroName[i]) + strlen(MacroDef[i]) + 2;
        memrq(-l,sizeof(char));
    }
    if (opt) 
        printf1("Deleted %d macro definition(s). Current memory: %d bytes.\n",MacroN,MemReq);
    MacroN = MacroNLen = 0;
    return(0);
}


