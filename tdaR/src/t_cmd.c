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
#include <string.h>
#if defined(_WIN32) || defined(_WIN64) || defined(__MINGW32__) || defined(__MINGW64__)
#  include <direct.h>
#  define getcwd _getcwd
#  define chdir  _chdir
#else
#  include <unistd.h>
#endif
#include "t_top.h"
#include "t_seqm.h"
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

#ifdef TDA_R_PACKAGE
/* rdataframe lives in tdaR/src/tda_rdf.c and only the package builds it. */
int rdataframe(TDAContext *ctx);
#endif
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


int get_ncmd(TDAContext *ctx, int fn);           
int t_exec(TDAContext *ctx, char *cmd);           
int t_execute(TDAContext *ctx);           
void cmd_err(TDAContext *ctx, int opt);        
int exec_macro(TDAContext *ctx, char *s);   
int check_macro(TDAContext *ctx, char *s,int l);
int new_macro(TDAContext *ctx);
void prn_merr(TDAContext *ctx, char *s,int opt);
int mlist(TDAContext *ctx);       
int mclear(TDAContext *ctx, int opt);    
void clear_tda(TDAContext *ctx);     

#define CFNMAX 10           /* max number of command file levels            */

/* ------------------------------------------------------------------------ */
/*  clear_tda()     clear all allocated memory                              */

void clear_tda(TDAContext *ctx)            
{
    clear_a(ctx);          /* clear data etc */
    ps_close(ctx, 0,1,0);    /* close PostScript output file */
    mat_free(ctx);         /* clear all matrices */
    mclear(ctx, 0);          /* clear macros */
    init_hlp(ctx, 0);        /* free help */
    dblock_alloc(ctx, 0);    /* free dblock */
    svb_alloc(ctx, 0);       /* free string buffer */
    gdd_free(ctx, 0);        /* free gdd(ctx, if perm) */
    ctx->IERRFlg = 0;
}

/* ------------------------------------------------------------------------ */
/*  get_ncmd(fn)    get next command from command file CFFILE[fn].          */
/*                  if successful put command into CmdBuf and return 1,     */
/*                  if no more command return 0, if error return -1.        */

#ifdef TDA_R_PACKAGE
/* Under R the commands arrive as text, not as a file: tda_cf_set() hands
   the text over before the run, and the top-level "cf=" reads it from
   memory through cf_gets(), line for line as fgets() would.  A nested
   cf= inside the commands still names a real file.  CF_MEM is the
   sentinel that stands in CFFILE[] for the memory source; no stdio call
   ever sees it. */
static const char *cf_mem_text = NULL;
static size_t cf_mem_pos = 0;
static int cf_mem_sentinel;
#define CF_MEM ((FILE *)(void *)&cf_mem_sentinel)

void tda_cf_set(const char *text)
{
    cf_mem_text = text;
    cf_mem_pos = 0;
}

static char *cf_gets(char *buf, int n, FILE *fd)
{
    size_t k = 0;

    if (fd != CF_MEM)
        return fgets(buf, n, fd);
    if (cf_mem_text == NULL || cf_mem_text[cf_mem_pos] == '\0' || n < 2)
        return NULL;
    while (k < (size_t)n - 1 && cf_mem_text[cf_mem_pos] != '\0') {
        buf[k++] = cf_mem_text[cf_mem_pos++];
        if (buf[k - 1] == '\n')
            break;
    }
    buf[k] = '\0';
    return buf;
}

static void cf_close(FILE *fd)
{
    if (fd != CF_MEM)
        fclose(fd);
}
#else
#define cf_gets(buf, n, fd) fgets(buf, n, fd)
#define cf_close(fd) fclose(fd)
#endif

int get_ncmd(TDAContext *ctx, int fn)            
{
    register int fin,cnt;
    register char cflag,*p,*q;
    int bflag,hflag;

    hflag = bflag = 0;
    q = p = ctx->CmdBuf;
    fin = cnt = 0;

    while (fin == 0) {

        if (ctx->CmdBufL - cnt < 200) {
            p_err(ctx, -43,2);
            return(-1);
        }
        if (cf_gets(p,ctx->CmdBufL - cnt,ctx->CFFILE[fn]) == NULL)  
            break;

        q = skip_b(ctx, p);
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
                if (cnt == 4 && !strncmp(ctx->CmdBuf,"help",4)) {
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
       
    p = ctx->CmdBuf + strlen(ctx->CmdBuf) - 1;
    while (p >= ctx->CmdBuf && *p && (*p == '\n' || *p == LF || *p == CR)) {
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

int t_exec(TDAContext *ctx, char *cmd)            
{
    int n,err = -1;

    if (strncmp(cmd,"cf=",3)) {
        if (cmd != ctx->CmdBuf) memmove(ctx->CmdBuf, cmd, strlen(cmd)+1);
        err = t_execute(ctx);
    
        if (err == 1)               /* check macro */
            err = exec_macro(ctx, ctx->CmdBuf);

        if (err > 0)  
            cmd_err(ctx, 0);             /* unknown command */

        else if (err < 0 && (ctx->SILENTFlg == 1 || ctx->SILENTFlg == 3)) {
            cmd_err(ctx, 1);     
            if (ctx->IERRFlg == 0)
                ctx->SILENTFlg = 0;
        }
        if (err && ctx->IERRFlg)         /* ignore errors */
            err = 0;
        goto EXECFin;
    }
    if (ctx->CFN >= CFNMAX) {
        p_err(ctx, -37,1);       /* reached max level of command files */
        goto EXECFin;
    }   
#ifdef TDA_R_PACKAGE
    if (ctx->CFN == 0 && cf_mem_text != NULL)
        ctx->CFFILE[0] = CF_MEM;
    else
#endif
    if (!(ctx->CFFILE[ctx->CFN] = fopen(cmd + 3,OPEN_RD))) {
        ps_err(ctx, -1,cmd + 3,1);
        goto EXECFin;
    }
    {
        /* print only the filename, not the full path,
           so output matches regardless of invocation directory */
        const char *_cf = cmd + 3, *_p;
        for (_p = _cf; *_p; ++_p)
            if (*_p == '/' || *_p == '\\')
                _cf = _p + 1;
        printf1(ctx, "Reading command file: %s\n", _cf);
    }
    prnchar(ctx, '=',LLEN,1);

    {
        const char *cf_path = cmd + 3;
        const char *p, *last_sep = NULL;
        char *orig_dir = NULL;
        int   did_chdir = 0;
        for (p = cf_path; *p; ++p)
            if (*p == '/' || *p == '\\')
                last_sep = p;
        if (last_sep) {
            size_t dlen = (size_t)(last_sep - cf_path);
            char  *cf_dir = (char *)malloc(dlen + 1);
            orig_dir      = (char *)malloc(4096);
            if (cf_dir && orig_dir && getcwd(orig_dir, 4096)) {
                memcpy(cf_dir, cf_path, dlen);
                cf_dir[dlen] = '\0';
                if (chdir(cf_dir) == 0)
                    did_chdir = 1;
            }
            free(cf_dir);
        }
        ctx->CFN++;
        while (ctx->CFN > 0) {
            n = get_ncmd(ctx, ctx->CFN - 1);
            if (n == 1) {
                if (t_exec(ctx, ctx->CmdBuf)) {
                    if (did_chdir && chdir(orig_dir))
                        printf1(ctx, "Warning: could not return to %s.\n",orig_dir);
                    free(orig_dir);
                    goto EXECFin;
                }
            }
            else if (n < 0) {
                printf1(ctx, "Command file ends with an incomplete command.\n");
                if (did_chdir && chdir(orig_dir))
                    printf1(ctx, "Warning: could not return to %s.\n",orig_dir);
                free(orig_dir);
                goto EXECFin;
            }
            else {
                ctx->CFN--;
                cf_close(ctx->CFFILE[ctx->CFN]);
                break;
            }
        }
        if (did_chdir && chdir(orig_dir))
            printf1(ctx, "Warning: could not return to %s.\n",orig_dir);
        free(orig_dir);
    }
    err = 0;

EXECFin:
    if (err) {
        while (ctx->CFN > 0) {
            ctx->CFN--;
            cf_close(ctx->CFFILE[ctx->CFN]);
        }
        if (ctx->REPLev >= 0)
            rep_free(ctx);
        err = -1;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  t_execute()     execute command in CmdBuf.                              */
/*                  return 0 if OK, -1 if error, 1 if cmd not defined.      */

int t_execute(TDAContext *ctx)         
{
    int n,r,lflg,err;
    char *p = ctx->CmdBuf;

    lflg = 0;
    err = 1;

    if (ctx->REPLev >= 0) {              /* save command */
        r = rep_scmd(ctx, ctx->CmdBuf);       /* don't execute if r = 1 */
        if (r) {
            if (r == 1) {
                prn_nex(ctx, ctx->CmdBuf);
                return(0);
            }
            else   
                return(-1);
        }
    }
    r = check_if(ctx, ctx->CmdBuf);           /* check for if/endif */
    if (r != 1)  
        return(r);
       
    r = check_break(ctx, ctx->CmdBuf);        /* check for break */
    if (r != 1)  
        return(r);

    switch (ctx->CmdBuf[0]) {

        case 'a':   if (!strncmp(p,"arcd",4))           /* arc */
                        err = arcd(ctx);
                    else if (!strcmp(p,"arcc"))         /* arcc */
                        err = arcc(ctx);
                    else if (!strncmp(p,"arcvc",5))     /* arcvc */
                        err = arcvc(ctx);
                    else if (!strncmp(p,"arcv",4))      /* arcv */
                        err = arcv(ctx);
                    else if (!strncmp(p,"atab",4))      /* atab */
                        err = atab(ctx);
                    else if (!strncmp(p,"acl",3))       /* acl */
                        err = acl(ctx);

                    break;

        case 'b':   if (!strncmp(p,"brr",3))            /* brr */
                        err = brr(ctx);
                    else if (!strncmp(p,"bfa",3))       /* bfa */
                        err = bfa(ctx);
                    else if (!strncmp(p,"bfc",3))       /* bfc */
                        err = bfc(ctx);
                    else if (!strncmp(p,"becl",4))      /* becl */
                        err = becl(ctx);
                    break;

        case 'c':   if (!strncmp(p,"ccnt=",5))           /* ccnt */
                        err = dm_ccnt(ctx);
                    else if (!strncmp(p,"clearnl",7))   /* clearnl */
                        err = clearnl(ctx);
                    else if (!strncmp(p,"clear",5))     /* clear */
                        err = clear(ctx);
                    else if (!strncmp(p,"cwt",3))       /* cwt */
                        err = cwt(ctx);
                    else if (!strncmp(p,"cov",3))       /* cov */
                        err = pcov(ctx, 0);
                    else if (!strncmp(p,"corr",4))      /* corr */
                        err = pcov(ctx, 1);
                    else if (!strncmp(p,"com",3))       /* com */
                        err = com(ctx);
                    else if (!strncmp(p,"conj",4))      /* conj */
                        err = tda_conj(ctx);
                    else if (!strncmp(p,"cro",3))       /* cro */
                        err = cro(ctx);
                    else if (!strncmp(p,"clpyr",5))     /* clpyr */
                        err = clpyr(ctx);
                    else if (!strncmp(p,"clp",3))       /* clp */
                        err = clp(ctx);
                    else if (!strncmp(p,"clu",3))       /* clu */
                        err = clu(ctx);
                    break;

        case 'd':   if (!strncmp(p,"dump",4))           /* dump */
                        err = dm_dump(ctx);
                    else if (!strncmp(p,"dsplit",6))    /* dsplit */
                        err = dm_dsplit(ctx);
                    else if (!strcmp(p,"data1"))        /* data1 */
                        err = data(ctx, 1);
                    else if (!strcmp(p,"data"))         /* data */
                        err = data(ctx, 0);
                    else if (!strncmp(p,"dstat",5))     /* dstat */
                        err = dstat(ctx);  
                    else if (!strncmp(p,"dplot",5))     /* dplot */
                        err = dplot(ctx);
                    else if (!strncmp(p,"dblock",6))    /* dblock */
                        err = dblock(ctx);
                    else if (!strncmp(p,"dma",3))       /* dma */
                        err = dma(ctx);
                    else if (!strncmp(p,"dmet1",5))     /* dmet1 */
                        err = dmet1(ctx);
                    else if (!strncmp(p,"dmet",4))      /* dmet */
                        err = dmet(ctx);
                    else if (!strncmp(p,"dple",4))      /* dple */
                        err = dple(ctx);
                    else if (!strncmp(p,"diple",5))     /* diple */
                        err = diple(ctx);
                    else if (!strncmp(p,"dltb",4))      /* dltb */
                        err = dltb(ctx);
                    break;

        case 'e':   if (!strncmp(p,"edef",4))           /* edef */
                        err = edef(ctx);
                    else if (!strcmp(p,"endwhile")) {   /* endwhile */
                        err = t_endrepeat(ctx, 0,ctx->REPLev);
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strcmp(p,"endrepeat")) {  /* endrepeat */
                        err = t_endrepeat(ctx, 1,ctx->REPLev);
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"epdat",5))     /* epdat */
                        err = epdat(ctx);
                    else if (!strncmp(p,"epsdat",6))    /* epsdat */
                        err = epsdat(ctx);
                    else if (!strncmp(p,"evalfi",6))    /* evalfi */
                        err = evalf(ctx, 1);
                    else if (!strncmp(p,"evalf",5))     /* evalf */
                        err = evalf(ctx, 0);
                    else if (!strncmp(p,"esort",5))     /* esort */
                        err = esort(ctx);
                    else if (!strncmp(p,"eskip",5))     /* eskip */
                        err = eskip(ctx);
                    else if (!strncmp(p,"eselect",7))   /* eselect */
                        err = eselect(ctx);
                    else if (!strncmp(p,"emerge",6))    /* emerge */
                        err = emerge(ctx);
                    else if (!strncmp(p,"ejoin",5))     /* ejoin */
                        err = ejoin(ctx);
                    else if (!strncmp(p,"expm",4))      /* expm */
                        err = e_expm(ctx);
                    else if (!strncmp(p,"etest",5))     /* etest */
                        err = etest(ctx);
                    break;

        case 'f':   if (!strncmp(p,"fml",3))            /* fml */
                        err = f_min(ctx, 0);
                    else if (!strncmp(p,"fmin",4))      /* fmin */
                        err = f_min(ctx, 1);
                    else if (!strncmp(p,"freg",4))      /* freg */
                        err = f_min(ctx, 2);
                    else if (!strncmp(p,"frml",4))      /* frml */
                        err = f_min(ctx, 3);
                    else if (!strncmp(p,"freq1",5))     /* freq1 */
                        err = mfreq(ctx, 1);
                    else if (!strncmp(p,"freq2",5))     /* freq2 */
                        err = mfreq(ctx, 2);
                    else if (!strncmp(p,"freq",4))      /* freq */
                        err = mfreq(ctx, 0);
                    break;

        case 'g':   if (!strncmp(p,"glm",3))            /* glm */
                        err = glm(ctx);
                    else if (!strncmp(p,"gmin",4))      /* gmin */
                        err = gmin(ctx);

                    else if (!strncmp(p,"gdd",3))       /* gdd */
                        err = gdd(ctx);
                    else if (!strncmp(p,"gcd",3))       /* gcd */
                        err = gcd(ctx);
                    else if (!strncmp(p,"gnc",3))       /* gnc */
                        err = gnc(ctx);
                    else if (!strncmp(p,"gcon",4))      /* gcon */
                        err = gcon(ctx);
                    else if (!strncmp(p,"gdcyc",5))     /* gdcyc */
                        err = gdcyc(ctx);
                    else if (!strncmp(p,"gcyc",4))      /* gcyc */
                        err = gcyc(ctx);
                    else if (!strncmp(p,"gdp",3))       /* gdp */
                        err = gdp(ctx);
                    else if (!strncmp(p,"gda",3))       /* gda */
                        err = gda(ctx);
                    else if (!strncmp(p,"gdu",3))       /* gdu */
                        err = gdu(ctx);
                    else if (!strncmp(p,"gde",3))       /* gde */
                        err = gde(ctx);
                    else if (!strncmp(p,"gdcset",6))    /* gdcset */
                        err = gdcset(ctx);
                    else if (!strncmp(p,"gsym",4))      /* gdot */
                        err = gsym(ctx);
                    else if (!strncmp(p,"gdot",4))      /* gdot */
                        err = gdot(ctx);
                    else if (!strncmp(p,"gni",3))       /* gni */
                        err = gni(ctx);
                    else if (!strncmp(p,"gdln",4))      /* gdln */
                        err = gdln(ctx);
                    else if (!strncmp(p,"gio",3))       /* gio */
                        err = gio(ctx);
                    else if (!strncmp(p,"gdcon",5))     /* gdcon */
                        err = gdcon(ctx);
                    else if (!strncmp(p,"gfcf",4))      /* gfcf */
                        err = gfcf(ctx);
                    else if (!strncmp(p,"gfc",3))       /* gfc */
                        err = gfc(ctx);
                    else if (!strncmp(p,"gbcf",4))      /* gbcf */
                        err = gbcf(ctx);
                    else if (!strncmp(p,"gep",3))       /* gep */
                        err = gep(ctx);
                    else if (!strncmp(p,"gtcl",4))      /* gtcl */
                        err = gtcl(ctx);
                    else if (!strncmp(p,"gflow",5))     /* gflow */
                        err = gflow(ctx);
                    else if (!strncmp(p,"gcliq",5))     /* gcliq */
                        err = gcliq(ctx);
                    else if (!strncmp(p,"ggcliq",6))    /* ggcliq */
                        err = ggcliq(ctx);
                    else if (!strncmp(p,"gev",3))       /* gev */
                        err = gev(ctx);
                    else if (!strncmp(p,"gsort",5))     /* gsort */
                        err = gsort(ctx);
                    else if (!strncmp(p,"gst",3))       /* gst */
                        err = gst(ctx);
                    else if (!strncmp(p,"gnst",4))      /* gnst */
                        err = gnst(ctx);
                    else if (!strncmp(p,"gmst",4))      /* gmst */
                        err = gmst(ctx);
                    else if (!strncmp(p,"gcut",4))      /* gcut */
                        err = gcut(ctx);
                    else if (!strncmp(p,"gsp",3))       /* gsp */
                        err = gsp(ctx);
                    else if (!strncmp(p,"giset",5))     /* giset */
                        err = giset(ctx);
                    else if (!strncmp(p,"gpro",4))      /* gpro */
                        err = gpro(ctx);
                    else if (!strncmp(p,"gcset",5))     /* gcset */
                        err = gcset(ctx);
                    else if (!strncmp(p,"gcni",4))      /* gcni */
                        err = gcni(ctx);
                    else if (!strncmp(p,"gap",3))       /* gap */
                        err = gap(ctx);
                    else if (!strncmp(p,"gqap",4))      /* gqap */
                        err = gqap(ctx);
                    else if (!strncmp(p,"gloc",4))      /* gloc */
                        err = gloc(ctx);
                    else if (!strncmp(p,"ghd1",4))      /* ghd1 */
                        err = ghd1(ctx);
                    else if (!strncmp(p,"ghd",3))       /* ghd */
                        err = ghd(ctx);

                    else if (!strncmp(p,"gdf",3))       /* gdf */
                        err = gdf(ctx);
                    else if (!strncmp(p,"gtopo",5))     /* gtopo */
                        err = gtopo(ctx);                

                    break;

        case 'h':   if (!strncmp(p,"help=",5))          /* help */
                        err = help(ctx, p + 5);
                    else if (!strncmp(p,"help",4))  
                        err = help(ctx, p + 4);
                    else if (!strncmp(p,"hclsp",5))     /* hclsp */
                        err = hclsp(ctx);
                    else if (!strncmp(p,"hcls",4))      /* hcls */
                        err = hcls(ctx);
                    else if (!strncmp(p,"hcld",4))      /* hcld */
                        err = hcld(ctx);
                    break;

        case 'i':   if (!strncmp(p,"intp",4))           /* intp */
                        err = intp(ctx);
                    else if (!strncmp(p,"int",3))       /* integration */
                        err = t_int(ctx);
                    else if (!strcmp(p,"ierr")) {       /* ierr */
                        ctx->IERRFlg = 1;
                        err = 0;            
                    }
                    else if (!strncmp(p,"ineq",4))      /* ineq */
                        err = ineq(ctx);
                    else if (!strncmp(p,"imean",5))     /* imean */
                        err = imean(ctx);
                    else if (!strncmp(p,"icorr",5))     /* icorr */
                        err = icov2(ctx, 5);
                    else if (!strncmp(p,"icov",4))      /* icov */
                        err = icov2(ctx, 4);
                    else if (!strncmp(p,"ivar1",5))     /* ivar1 */
                        err = ivar1(ctx);
                    else if (!strncmp(p,"ivar",4))      /* ivar */
                        err = ivar(ctx);
                    else if (!strncmp(p,"igini",5))     /* igini */
                        err = igini(ctx);
                    else if (!strncmp(p,"idf(",4))      /* idf */
                        err = idf(ctx);
                    else if (!strncmp(p,"iddf(",5))     /* iddf */
                        err = iddf(ctx);
                    else if (!strncmp(p,"ilsreg",6))    /* ilsreg */
                        err = ilsreg(ctx);
                    else if (!strncmp(p,"inpreg",6))    /* inpreg */
                        err = inpreg(ctx);
                    else if (!strncmp(p,"imreg",5))     /* imreg */
                        err = imreg(ctx);
                    else if (!strncmp(p,"ivls",4))      /* ivls */
                        err = ivls(ctx);
                    else if (!strncmp(p,"ivreg1",6))    /* ivreg1 */
                        err = ivreg1(ctx);
                    else if (!strncmp(p,"ivreg2",6))    /* ivreg2 */
                        err = ivreg2(ctx);
                    else if (!strncmp(p,"ivreg",5))     /* ivreg */
                        err = ivreg(ctx);
                    else if (!strncmp(p,"indep",5))     /* indep */
                        err = indep(ctx);
                    break;

        case 'l':   if (!strncmp(p,"lcnt",4))           /* lcnt */
                        err = dm_lcnt(ctx);
                    else if (!strncmp(p,"lsreg1",6))    /* lsreg1 */
                        err = lsreg1(ctx);
                    else if (!strncmp(p,"lsreg",5))     /* lsreg */
                        err = lsreg(ctx);
                    else if (!strncmp(p,"l1reg",5))     /* l1reg */
                        err = l1reg(ctx);
                    else if (!strncmp(p,"ltb",3))       /* life table */
                        err = ltb(ctx);
                    else if (!strncmp(p,"local(",6))    /* local command */
                        err = alloc_local(ctx, p + 6);
                    else if (!strncmp(p,"loglin",6))    /* loglin */
                        err = loglin(ctx);
                    break;

        case 'm':   if (!strcmp(p,"mem")) {             /* memory use */
                        err = mem_use(ctx);
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strcmp(p,"mach"))         /* machine parameters */
                        err = machp(ctx, 1);
                    else if (!strncmp(p,"macrodef",8))  /* new macro */
                        err = new_macro(ctx);
                    else if (!strncmp(p,"macrolist",9)) /* list macros */
                        err = mlist(ctx);
                    else if (!strncmp(p,"macroclear",10))   /* clear macros */
                        err = mclear(ctx, 1);
                    else if (!strncmp(p,"mdsx",4))      /* mdsx */
                        err = mdsx(ctx);
                    else if (!strncmp(p,"mdsm",4))      /* mdsm */
                        err = mdsm(ctx);
                    /*  mdsn1 must be tested BEFORE mdsn: the test is a
                        4-character prefix match, so "mdsn1(...)" matched
                        "mdsn", ran the wrong command, and then reported a
                        syntax error on the leftover "1(...)".  mdsn1 was
                        therefore unreachable -- it showed as 238 lines of
                        dead code with no call site -- although it is a
                        real implementation, the censored counterpart of
                        mdsn in the way lsreg1 is of lsreg.  lsreg1 is
                        ordered correctly a few lines above; this one was
                        simply missing.  */
                    else if (!strncmp(p,"mdsn1",5))     /* mdsn1 */
                        err = mdsn1(ctx);
                    else if (!strncmp(p,"mdsn",4))      /* mdsn */
                        err = mdsn(ctx);
                    else if (!strncmp(p,"mdsc",4))      /* mdsc */
                        err = mdsc(ctx);
                    else if (!strncmp(p,"mdsr",4))      /* mdsr */
                        err = mdsr(ctx);
                    else if (!strncmp(p,"mproc",5))     /* mproc */
                        err = mproc(ctx);
                    else if (!strncmp(p,"mreg",4))      /* mreg */
                        err = mreg(ctx);
                    else if (!strncmp(p,"mparse",6))    /* parse expression */
                        err = mparse(ctx);
                    else {                              /* matrix commands */
                        err = t_mat(ctx);   
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    break;

        case 'n':   if (!strncmp(p,"nvar(",5))          /* new variables */
                        err = new_var(ctx);
                    else if (!strncmp(p,"nlist",5))     /* nlist */
                        err = nlist(ctx);
                    else if (!strncmp(p,"ndvar",5))     /* ndvar */
                        err = ndvar(ctx);
                    else if (!strncmp(p,"niset",5))     /* niset */
                        err = niset(ctx);
                    else if (!strncmp(p,"nlreg",5))     /* nlreg */
                        err = nlreg(ctx);
                    else if (!strncmp(p,"npreg",5))     /* npreg */
                        err = npreg(ctx);
                    else if (!strncmp(p,"nncl",4))      /* nncl */
                        err = nncl(ctx);
                    else if (!strncmp(p,"nmca",4))      /* nmca */
                        err = nmca(ctx);
                    break;

        case 'p':   if (!strncmp(p,"print(",6)) {       /* print text */
                        err = prn_txt(ctx);
                        if (err == 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"plotp3",6))    /* plot 3d polygon */
                        err = plotp3(ctx, 0);
                    else if (!strncmp(p,"plot3",5))     /* plot 3d polygon */
                        err = plotp3(ctx, 1);
                    else if (!strncmp(p,"plcurv3",7))   /* plot 3d curve */
                        err = plcurv3(ctx);
                    else if (!strncmp(p,"plsurf3d",8))  /* plot 3d surface */
                        err = plsurf3d(ctx);
                    else if (!strncmp(p,"plsurf3",7))   /* plot 3d surface */
                        err = plsurf3(ctx);
                    else if (!strncmp(p,"pltext3",7))   /* plot text */
                        err = pltext3(ctx);
                    else if (!strncmp(p,"plcirc3",7))   /* plot circle */
                        err = plcirc3(ctx);
                    else if (!strncmp(p,"plglob3",7))   /* plot globe */
                        err = plglob3(ctx);

                    else if (!strncmp(p,"parse=",6))    /* parse expression */
                        err = t_parse(ctx);
                    else if (!strncmp(p,"pdata",5))     /* pdata */
                        err = pdata(ctx);
                    else if (!strncmp(p,"pdatr",5))     /* pdatr */
                        err = pdatr(ctx);
                    else if (!strncmp(p,"pdatd",5))     /* pdatd */
                        err = pdatd(ctx);
                    else if (!strncmp(p,"ple",3))       /* product-limit */
                        err = ple(ctx);
                    else if (!strncmp(p,"pcyc",4))      /* pcyc */
                        err = pcyc(ctx);
                    else if (!strncmp(p,"ptree",5))     /* ptree */
                        err = ptree(ctx);
                    else if (!strncmp(p,"pltree",6))    /* pltree */
                        err = pltree(ctx);

                    else if (!strncmp(p,"psfile=",7))   /* psfile */
                        err = psfile(ctx);     
                    else if (!strncmp(p,"psetupg",7))   /* psetupg */
                        err = psetupg(ctx);     
                    else if (!strncmp(p,"psetup3",7))   /* psetup3 */
                        err = psetup3(ctx);     
                    else if (!strncmp(p,"psetup",6))    /* psetup */
                        err = psetup(ctx);     
                    else if (!strcmp(p,"psclose"))      /* psclose */
                        err = ps_close(ctx, 1,1,0);     
                    else if (!strncmp(p,"plxa",4))      /* plot x axis */
                        err = pl_axis(ctx, 0,0,0,1,0.2);
                    else if (!strncmp(p,"plya",4))      /* plot y axis */
                        err = pl_axis(ctx, 1,0,0,1,0.2);
                    else if (!strncmp(p,"plframe",7))   /* plot frame */
                        err = pl_frame(ctx);
                    else if (!strncmp(p,"plabel",6))    /* plot label */
                        err = pl_label(ctx, 0);
                    else if (!strncmp(p,"pxlabel",7))   /* plot label */
                        err = pl_label(ctx, 1);
                    else if (!strncmp(p,"pylabel",7))   /* plot label */
                        err = pl_label(ctx, 2);
                    else if (!strncmp(p,"plxgrid",7))   /* grid lines x axis */
                        err = ps_grid(ctx, 0);
                    else if (!strncmp(p,"plygrid",7))   /* grid lines y axis */
                        err = ps_grid(ctx, 1);
                    else if (!strncmp(p,"pltext",6))    /* plot text */
                        err = pl_text(ctx);
                    else if (!strncmp(p,"plrec",5))     /* plot rectangle */
                        err = pl_rec(ctx);

                    else if (!strncmp(p,"plotp",5))     /* plot polygon */
                        err = pl_plotp(ctx, 0);
                    else if (!strncmp(p,"plotm",5))     /* plot polygon */
                        err = pl_plotm(ctx);
                    else if (!strncmp(p,"plotf",5))     /* plot function */
                        err = pl_plotf(ctx);
                    else if (!strncmp(p,"ploth",5))     /* plot histogram */
                        err = pl_ploth(ctx);
                    else if (!strncmp(p,"plotd",5))     /* plot density */
                        err = pl_plotd(ctx);
                    else if (!strncmp(p,"ploto",5))     /* plot circle */
                        err = pl_ploto(ctx, 0);
                    else if (!strncmp(p,"plote",5))     /* plot ellipse */
                        err = pl_ploto(ctx, 1);
                    else if (!strncmp(p,"plotk",5))     /* plot arc */
                        err = pl_plotk(ctx);
                    else if (!strncmp(p,"plotch",6))    /* convex hull */
                        err = pl_plotch(ctx);
                    else if (!strncmp(p,"plotcm",6))    /* contour plot */
                        err = pl_plotcm(ctx);
                    else if (!strncmp(p,"plotc",5))     /* contour plot */
                        err = pl_plotc(ctx);
                    else if (!strncmp(p,"plotr",5))     /* grey-scaled relief */
                        err = pl_plotr(ctx);
                    else if (!strncmp(p,"plotsp",6))    /* smoothing: Akima */
                        err = pl_plots(ctx, 0);
                    else if (!strncmp(p,"plots",5))     /* smoothing: Akima */
                        err = pl_plots(ctx, 1);
                    else if (!strncmp(p,"plot",4))      /* plot polygon */
                        err = pl_plotp(ctx, 1);
                    else if (!strncmp(p,"plg",3))       /* plot graph */
                        err = pl_plg(ctx);

                    break;

        case 'q':   if (!strncmp(p,"quant",5))          /* quantiles */
                        err = quant(ctx);  
                    else if (!strncmp(p,"qreg",4))      /* qreg */
                        err = qreg(ctx);
                    break;

        case 'r':   if (!strncmp(p,"rsys",4))           /* read system file */
                        err = rd_sys(ctx);
                    else if (!strncmp(p,"repeat",6)) {  /* repeat */
                        err = t_repeat(ctx, 1,ctx->REPLev);
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    else if (!strncmp(p,"rcorr",5))     /* rcorr */
                        err = pcov(ctx, 2);
                    else if (!strncmp(p,"recode",6))    /* recode */
                        err = recode(ctx);
                    else if (!strncmp(p,"rspss1",6))    /* read spss sav file */
                        err = rd_spss1(ctx);
                    else if (!strncmp(p,"rspss",5))     /* read spss por file */
                        err = rd_spss(ctx);
                    else if (!strncmp(p,"rstata",6))    /* read stata file */
                        err = rd_stata(ctx);
#ifdef TDA_R_PACKAGE
                    /* Only the R package can be handed an R data frame, so
                       the standalone program never sees this command. */
                    else if (!strncmp(p,"rdataframe",10))
                        err = rdataframe(ctx);
#endif
                    else if (!strncmp(p,"rxls",4))      /* read xls file */
                        err = rxls(ctx);
                    else if (!strncmp(p,"rucinet",7))   /* read xls file */
                        err = rucinet(ctx);
                    else if (!strncmp(p,"rdbf",4))      /* read dbf file */
                        err = rdbf(ctx);
                    else if (!strncmp(p,"rplz",4))      /* read plz file */
                        err = rplz(ctx);
                    else if (!strncmp(p,"rcsv",4))      /* read csv file */
                        err = rcsv(ctx);
/*                  else if (!strncmp(p,"rdxf",4))  */  /* read dxf file */
/*                      err = rdxf();       */         

                    else if (!strncmp(p,"rate",4))      /* rate models */
                        err = rate(ctx);
                    else if (!strncmp(p,"rsys",4))      /* read system file */
                        err = rd_sys(ctx);
                    else if (!strncmp(p,"repsel",6))    /* repsel */
                        err = repsel(ctx);
                    else if (!strncmp(p,"range",5))     /* range */
                        err = range(ctx);
                    else if (!strncmp(p,"rmod",4))      /* rmod */
                        err = rmod(ctx);
                    else if (!strncmp(p,"rod",3))       /* rod */
                        err = rod(ctx);
                    else if (!strncmp(p,"rap",3))       /* rap */
                        err = rap(ctx);
                    else if (!strncmp(p,"rfit1",5))      /* rfit1 */
                        err = rfit1(ctx);
                    else if (!strncmp(p,"rfit",4))      /* rfit */
                        err = rfit(ctx);

                    break;

        case 's':   if (sscanf(p,"silent=%d",&n) == 1) {
                        ctx->SILENTFlg = imax(ctx, imin(ctx, n,3),-1);
                        lflg = 1;
                        err = 0;
                    }
                    else if (!strncmp(p,"silent",6)) {
                        printf1(ctx, "Current value of silent: %d\n",ctx->SILENTFlg);
                        err = 0;
                    }
                    else if (!strncmp(p,"sdnvar",6))    /* sdnvar */
                        err = sdnvar(ctx);
                    else if (!strncmp(p,"sdpdata",7))   /* sdpdata */
                        err = sdpdata(ctx);
                    else if (!strncmp(p,"sdplot31",8))  /* sdplot31 */
                        err = sdplot31(ctx);
                    else if (!strncmp(p,"sdplot32",8))  /* sdplot32 */
                        err = sdplot32(ctx);
                    else if (!strncmp(p,"sdplot33",8))  /* sdplot33 */
                        err = sdplot33(ctx);
                    else if (!strncmp(p,"sdplot",6))    /* sdplot */
                        err = sdplot(ctx);
                    else if (!strncmp(p,"sdrel",5))     /* sdrel */
                        err = sdrel(ctx);
                    else if (!strncmp(p,"sdnl",4))      /* sdnl */
                        err = sdnl(ctx);
                    else if (!strncmp(p,"sdppol",6))    /* sdppol */
                        err = sdppol(ctx);
                    else if (!strncmp(p,"sdlpol",6))    /* sdlpol */
                        err = sdlpol(ctx);
                    else if (!strncmp(p,"sdipol",6))    /* sdipol */
                        err = sdipol(ctx);
                    else if (!strncmp(p,"sdcpol",6))    /* sdcpol */
                        err = sdcpol(ctx);
                    else if (!strncmp(p,"sddcwp",6))    /* sddcwp */
                        err = sddcwp(ctx);
                    else if (!strncmp(p,"sdgen",5))     /* sdgen */
                        err = sdgen(ctx);
                    else if (!strncmp(p,"sdgshhs",7))   /* sdgshhs */
                        err = sdgshhs(ctx);
                    else if (!strncmp(p,"sdshp",5))     /* sdshp */
                        err = sdshp(ctx);                 
                    else if (!strncmp(p,"sde00",5))     /* sde00 */
                        err = sde00(ctx);                 
                    else if (!strncmp(p,"sdinf",5))     /* sdinf */
                        err = sdinf(ctx);                 
                    else if (!strncmp(p,"sdencl",6))    /* sdencl */
                        err = sdencl(ctx);                 
                    else if (!strncmp(p,"sdsel",5))     /* sdsel */
                        err = sdsel(ctx);                 
                    else if (!strncmp(p,"sdvd",4))      /* sdvd */
                        err = sdvd(ctx);
                    else if (!strncmp(p,"sdclip",6))    /* sdclip */
                        err = sdclip(ctx);
                    else if (!strncmp(p,"sdpgrat",7))   /* sdpgrat */
                        err = sdpgrat(ctx);
                    else if (!strncmp(p,"sdpgeo",6))    /* sdpgeo */
                        err = sdpgeo(ctx);
                    else if (!strncmp(p,"sdpmap",6))    /* sdpmap */
                        err = sdpmap(ctx);

                    else if (!strncmp(p,"seq",3))       /* sequence data */
                        err = t_seq(ctx);
                    else if (!strncmp(p,"segr",4))      /* segr */
                        err = segr(ctx);
                    else if (!strncmp(p,"sma",3))       /* sma */
                        err = sma(ctx);
                    else if (!strncmp(p,"smd",3))       /* smd */
                        err = smd(ctx);
                    else if (!strncmp(p,"spl",3))       /* spl */
                        err = spl(ctx);
                    else if (!strncmp(p,"scplot",6))    /* scplot */
                        err = scplot(ctx);
                    else if (!strncmp(p,"sddf(",5))     /* sddf */
                        err = sddf(ctx);
                    else if (!strncmp(p,"spmod",5))     /* spmod */
                        err = spmod(ctx);
                    else if (!strncmp(p,"sga",3))       /* sga */
                        err = sga(ctx);
                    else if (!strncmp(p,"subm",4))      /* subm */
                        err = subm(ctx);
                    else if (!strncmp(p,"scla",4))      /* scla */
                        err = scla(ctx);

                    break;

        case 't':   if (!strcmp(p,"time"))              /* current time */
                        err = prntime(ctx);
                    else if (!strncmp(p,"tsel=",5))     /* tsel  */
                        err = tsel(ctx);
                    else if (!strncmp(p,"tnet",4))      /* tnet */
                        err = tnet(ctx);
                    else if (!strncmp(p,"triang",6))    /* triang */
                        err = triang(ctx);
                    break;
        
        case 'u':   if (!strncmp(p,"uds",3))            /* uds */
                        err = uds(ctx);
                    else if (!strncmp(p,"unf",3))       /* unf */
                        err = unf(ctx);
                    else if (!strncmp(p,"ucl",3))       /* ucl */
                        err = ucl(ctx);
                    else
                        break;

        case 'w':   if (!strncmp(p,"wspss1",6))         /* write spss sav file */
                        err = wr_spss1(ctx);
                    else if (!strncmp(p,"wspss",5))     /* write spss por file */
                        err = wr_spss(ctx);
                    else if (!strncmp(p,"wstata",6))    /* write stata file */
                        err = wr_stata(ctx);
                    else if (!strncmp(p,"wsys",4))      /* write system file */
                        err = wr_sys(ctx);
                    else if (!strncmp(p,"while",5)) {   /* while */
                        err = t_repeat(ctx, 0,ctx->REPLev);
                        if (err == 0 && ctx->SILENTFlg >= 0)
                            lflg = 1;   
                    }
                    break;

        case 'x':   if (!strncmp(p,"xopen",5))          /* xopen */
                        err = xopen(ctx);     
                    else if (!strncmp(p,"xdelete",7))   /* xdelete */
                        err = xdelete(ctx);    
                    else if (!strcmp(p,"xlog1"))        /* xlog1 */
                        err = xlog1(ctx);      
                    else if (!strncmp(p,"xlog",4))      /* xlog */
                        err = xlog(ctx);     
                    else if (!strncmp(p,"xconh",5))     /* xconh */
                        err = xconh(ctx);     
                    else if (!strncmp(p,"xreg",4))      /* xreg */
                        err = xreg(ctx);     
                    else if (!strncmp(p,"xplotf",6))    /* xplotf */
                        err = xplotf(ctx);     
                    else if (!strncmp(p,"xplot",5))     /* xplot */
                        err = xplot(ctx);     
                    else if (!strncmp(p,"xf",2))        /* xf */
                        err = xfunc(ctx);     
                    else if (!strncmp(p,"xdens",5))     /* xdens */
                        err = xdens(ctx);     
                    else if (!strncmp(p,"xshow",5))     /* open plot window */
                        err = x_show(ctx);
                    break;

        case 'z':   if (!strncmp(p,"zreg1",5))          /* zreg1 */
                        err = zreg1(ctx);
                    else if (!strncmp(p,"zreg",4))           /* zreg */
                        err = zreg(ctx);
                    break;

        case '$':   err = check_cmd(ctx, 1);
                    if (err == 0) {
                        tda_out_flush();
                        if (system(ctx->CmdBuf+1))  
                            printf1(ctx, "Note: shell returned error.\n");
                        else
                            printf1(ctx, "\n> end of shell execution.\n");
                    }        
                    break;

        default:    lflg = 1;
                    break;     

    }
    if (err != 1) {
        if (ctx->REPLev >= 0 && err < 0)
            rep_free(ctx);
        if (lflg == 0)
            prnchar(ctx, '-',LLEN,1);    
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cmd_err(opt)   Give error message for command in CmdBuf.                */

void cmd_err(TDAContext *ctx, int opt)      
{
    register char *p = ctx->CmdBuf;

    if (strlen(p) < 1)
        return;

    if (opt == 0)
        printf1(ctx, "Unknown command: ");
    else
        printf1(ctx, "Error while executing command: ");

    while (*p && p < ctx->CmdBuf + 20)
        printf1(ctx, "%c",*p++);
    if (*p || opt)
        printf1(ctx, " ...");
    printf1(ctx, "\n");  
    prnchar(ctx, '-',LLEN,1);    
}

/* ------------------------------------------------------------------------ */
/*  exec_macro(s)   Try to execute macro s.                                 */
/*                                                                          */
/*                  a) if s is known macro expand and execute.              */
/*                  b) if s is expand=mname ...                             */
/*                     show expanded macro on stdout but do not execute.    */  
/*                                                                          */
/*  Return: 1 if unknown, 0 if OK, -1 if error.                             */

int exec_macro(TDAContext *ctx, char *s)    
{
    register int i;
    int err,n,m,nlen,aflag,np,l,eflag;
    register char *p,*q,*q1,*qq;
    char c,cflag;
    
    ctx->MacroExLevel += 1;
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
        n = check_macro(ctx, s,nlen);
        if (n < 0)
            goto MACROFin;
        *p = c;
        q = s;

        if (ctx->SILENTFlg < 0) {
            printf1(ctx, "Macro: %s\n",s);
            printf1(ctx, "Execution level: %d\n",ctx->MacroExLevel);
            prnchar(ctx, '-',LLEN,1);    
        }
        err = -1;

        /* get arguments into MacroArgs */

        qq = q = p;
        aflag = 0;
        if (c == '(') {
            q = skip_blev(ctx, p);
            if (q > p + 1) {
                if (*(q - 1) != ')') {
                    prn_merr(ctx, s,1);
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
                    prn_merr(ctx, s,0);
                    goto MACROFin;
                }
                if (*q == ',') {
                    ctx->MacroArgs[np] = NULL;
                }
                else {
                    q1 = skip_nc(ctx, q);
                    l = (int)(q1 - q);
                    if (!(ctx->MacroArgs[np] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                        p_err(ctx, -2,1);
                        goto MACROFin;
                    }
                    memrq(ctx, l + 1,sizeof(char));
                    strncpy(ctx->MacroArgs[np],q,(size_t)(l));
                    *(ctx->MacroArgs[np] + l) = '\0';
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
                    prn_merr(ctx, s,0);
                    goto MACROFin;
                }
                if (*p == ',')
                    ctx->MacroArgs[np] = NULL;
                else {
                    q = skip_nc(ctx, p);
                    l = (int)(q - p);
                    if (l == 0)
                        ctx->MacroArgs[np] = NULL;
                    else {
                        if (!(ctx->MacroArgs[np] = (char *)calloc((size_t)(l + 1),sizeof(char)))) {
                            p_err(ctx, -2,1);
                            goto MACROFin;
                        }
                        memrq(ctx, l + 1,sizeof(char));
                        strncpy(ctx->MacroArgs[np],p,(size_t)(l));
                        *(ctx->MacroArgs[np] + l) = '\0';
                    }
                    p = q;
                }
                np++;
            }
        }
        if (np > ctx->MacroNP[n])  
            printf1(ctx, "Warning: macro %s has %d argument(s).\n",ctx->MacroName[n],ctx->MacroNP[n]);

        /* insert arguments into macro definition and call exec()  */
        /* or show on stdout if eflag = 1 */

        ctx->CurMacroNum = n;

        if (eflag)  
            printf1(ctx, "Expansion of macro: %s\n",s); /*  MacroName[n]);   */
  
        p = ctx->MacroDef[n];
        while (*p) {
            q = ctx->CmdBuf;
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
                        p = skip_int(ctx, p + 1);
                        m--;
                        if (ctx->MacroArgs[m] != NULL) {
                            nlen = (int)(strlen(ctx->MacroArgs[m]));
                            if (l + nlen > ctx->CmdBufL - 10) {
                                p_err(ctx, -43,1);
                                goto MACROFin;
                            }     
                            strcpy(q,ctx->MacroArgs[m]);
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
                    if (l + 1 > ctx->CmdBufL - 10) {
                        p_err(ctx, -43,1);  
                        goto MACROFin;
                    }     
                    *q++ = *p++;
                    l++;
                }
            }
            *q = '\0';

            if (eflag)  
                printf1(ctx, "> %s\n",ctx->CmdBuf);
            else {
                err = t_exec(ctx, ctx->CmdBuf);
                if (err)
                    goto MACROFin;
            }
            while (*p == ';')
                p++;
        }
        if (eflag)
            prnchar(ctx, '-',LLEN,1);    

        else if (ctx->SILENTFlg < 0) {
            printf1(ctx, "End of macro: %s\n",ctx->MacroName[n]);
            printf1(ctx, "New execution level: %d\n",ctx->MacroExLevel - 1);
            prnchar(ctx, '-',LLEN,1);    
        }
        err = 0;
    }

MACROFin:
    for (i = 0; i < np; ++i) {
        if (ctx->MacroArgs[i] != NULL) {
            memrq(ctx,(int)(-strlen(ctx->MacroArgs[i]) - 1),sizeof(char));
            free(ctx->MacroArgs[i]);
            ctx->MacroArgs[i] = NULL;
        }
    }
    free_local(ctx);           /* free local matrices */
    ctx->MacroExLevel -= 1;
    ctx->CurMacroNum = -1;
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_macro(s,l)  Return number of macro if s is a currently existing   */
/*                    macro name, otherwise return -1.                      */
/*                                                                          */

int check_macro(TDAContext *ctx, char *s,int l)
{
    register int i;

    for (i = 0; i < ctx->MacroN; ++i) {
        if (!strncmp(s,ctx->MacroName[i],(size_t)(l)) && strlen(ctx->MacroName[i]) == (size_t)l)
            return(i);
    }
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  new_macro()         Get new macro.                                      */
/*                      Syntax: macrodef(name)={...};                       */
/*  ##                  Return 0 if OK, -1 if error.                        */
/*                                                                          */

int new_macro(TDAContext *ctx)
{
    register int i;
    int err,l,n,m,nlen;
    register char *s,*p,*q;
    short marg[MaxMacro + 1];

    if (ctx->s_new_macro_first) {
        for (i = 0; i < MaxMacro; ++i) 
            ctx->MacroArgs[i] = NULL;
        ctx->s_new_macro_first = 0;
    }

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);

    printf1(ctx, "New macro definition. Current memory: %d bytes.\n",ctx->MemReq);

    if (ctx->MacroN >= MaxMacro) {
        printf1(ctx, "Error: exceeded max number of macro definitions.\n");
        goto GMACROFin;
    }
    s = ctx->CmdBuf + 8;
    if (*s++ != '(' || *(p = skip_nc(ctx, s)) != ')' || (nlen = (int)(p - s)) < 1     
                                        || *++p != '=' || *++p != '{') {
        prn_merr(ctx, ctx->CmdBuf,2);        
        goto GMACROFin;
    }
    if (check_macro(ctx, s,nlen) >= 0) {
        *(s + nlen) = '\0';
        printf1(ctx, "Error: macro name %s already used.\n",s);
        goto GMACROFin;
    }
    q = skip_blev(ctx, p);
    if (*(q - 1) != '}') {
        prn_merr(ctx, ctx->CmdBuf,2);
        goto GMACROFin;
    }
    l = (int)(strlen(++p));
    if (l < 2) {
        prn_merr(ctx, ctx->CmdBuf,2);
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
                prn_merr(ctx, ctx->CmdBuf,2);
                printf1(ctx, "Range of arguments is: $1 - $%d\n",MaxMacro);
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
            prn_merr(ctx, ctx->CmdBuf,2);
            printf1(ctx, "Argument numbers should be contiguous.\n");
            goto GMACROFin;
        }
    }
    if (!(ctx->MacroName[ctx->MacroN] = (char *)calloc((size_t)(nlen + 1),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto GMACROFin;
    }
    if (!(ctx->MacroDef[ctx->MacroN] = (char *)calloc((size_t)(l),sizeof(char)))) {
        free(ctx->MacroName[ctx->MacroN]);
        p_err(ctx, -2,1);
        goto GMACROFin;
    }
    memrq(ctx, l + nlen + 1,sizeof(char));
    strncpy(ctx->MacroName[ctx->MacroN],s,(size_t)(nlen));
    strncpy(ctx->MacroDef[ctx->MacroN],p,(size_t)(l-1));

    ctx->MacroNP[ctx->MacroN] = (short)(n);
    if (ctx->MacroNLen < nlen)
        ctx->MacroNLen = nlen;

    printf1(ctx, "Added new macro: %s\n",ctx->MacroName[ctx->MacroN]);
    ctx->MacroN++;
    err = 0;

GMACROFin:
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  prn_merr(s,opt)     print macro error.                                  */

void prn_merr(TDAContext *ctx, char *s,int opt)
{
    register int i;
    register char *p;

    if (opt == 0) {
        printf1(ctx, "Error: exceeded maximal number of macro arguments.\n");
        return;
    }
    else if (opt == 1)
        printf1(ctx, "Syntax error: ");
    else if (opt == 2)
        printf1(ctx, "Error in macro definition: ");

    p = s;
    for (i = 0; i < 30; ++i) {
        if (!*p)
            break;
        printf1(ctx, "%c",*p++);
    }   
    if (*p)
        printf1(ctx, " ...\n");
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  mlist()         List currently defined macros: macrolist                */

int mlist(TDAContext *ctx)        
{
    register int i,j,l;
    register char *p;
    char buf[20];

    printf1(ctx, "Currently defined macros: ");

    if (ctx->MacroN == 0) {
        printf1(ctx, "none\n");
        return(0);
    }
    newline(ctx);

    for (i = 0; i < ctx->MacroN; ++i) {
        printf1(ctx, "%3d  %s",i + 1,ctx->MacroName[i]);
        l = (int)(strlen(ctx->MacroName[i]));      
        if (ctx->MacroNP[i] > 0) {
            snprintf(buf,sizeof(buf),"%d",ctx->MacroNP[i]);
            printf1(ctx, "($%s)",buf);
            l += (int)strlen(buf) + 3;
        }
        prnchar(ctx, ' ',ctx->MacroNLen + 5 - l, 0);
        printf1(ctx, " [");
        p = ctx->MacroDef[i];
        for (j = 0; j < 50; ++j) {
            if (!*p)
                break;
            printf1(ctx, "%c",*p++);
        }   
        if (*p)
            printf1(ctx, " ...");
        printf1(ctx, "]\n");
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  mclear()        Clear currently defined macros: macroclear;             */

int mclear(TDAContext *ctx, int opt)        
{
    register int i,l;

    if (ctx->MacroN == 0)
        return(0);

    for (i = 0; i < ctx->MacroN; ++i) {
        /*  The length was measured with strlen() AFTER both strings had
            been freed, so the memory accounting read freed memory --
            heap-use-after-free on macrodef.cf.  Measure first, then
            free.  */
        l = (int)(strlen(ctx->MacroName[i]) + strlen(ctx->MacroDef[i]) + 2);
        free(ctx->MacroName[i]);
        free(ctx->MacroDef[i]);
        memrq(ctx, -l,sizeof(char));
    }
    if (opt) 
        printf1(ctx, "Deleted %d macro definition(s). Current memory: %d bytes.\n",ctx->MacroN,ctx->MemReq);
    ctx->MacroN = ctx->MacroNLen = 0;
    return(0);
}
