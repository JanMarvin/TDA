/****************************************************************************/
/*  t_gmin                                                                  */
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
#include "t_ml.h"  
#include "t_tmin.h"  
#include "t_var.h"  
#include "t_gdat.h"  
#include "t_eval.h"  
#include "t_eval1.h"  
#include "t_eval3.h"  
#include "t_lsei.h"  
#include "t_cdf.h"  
#include "t_min.h"  
#include "t_con.h"  
#include "t_edat.h"  
#include "t_mat.h"  
#include "t_matf.h"  

/*  functions in t_gmin.c */

int f_min(int typ);
void gm_coeff(int n,double *x,int cov,double *diag);
int get_func(char *s,int opt,int *n3,int iflag,char *msg);
char *get_level(char *s,int *err); 
void prn_fexp(char *s,int opt); 
int get_farg(char *s);         
int fcomp(const void *arg1,const void *arg2);
void prn_feval(void);
int get_flval(double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d);
int get_filval(int n,double *x,int deriv,int *ni);
int get_fival(int n,double *x,int deriv,int aflag);
int prn_fsval(int n,double *x,int gmina,double *lb,double *ub,int opt);
void gm_pres(void);          
int get_ifval(double *fl,double *fu,int n,double *xl,double *xu,int deriv,
    double *gl,double *gu);
int get_ifival(int n,double *x,int deriv);
int get_icval(double *fl,double *fu,int n,double *xl,double *xu);          

/* ------------------------------------------------------------------------ */
/*  f_min(typ)      User-defined models/functions.                          */
/*                  typ = 0 : ML                                            */
/*                        1 : general functions                             */
/*                        2 : nonlinear regression                          */
/*                        3 : frml (user-def rate models)                   */
/*                                                                          */
/*          Return 0 if OK, otherwise -1.                                   */

int f_min(int typ)
{
    int err,n,cov,typ1,iflag;

    err = -1;
    iflag = 0;
    if (check_cmd(2))
        return(-1);

    n = 4;
    typ1 = typ;

    if (typ == 0) {
        printf1("ML estimation of user-defined model. ");
        n--;
    }
    else if (typ == 1)
        printf1("Function minimization. ");
    else if (typ == 2)
        printf1("Nonlinear regression. ");
    else if (typ == 3) {
        printf1("ML estimation of user-defined rate model. ");
        typ1 = 0;
    }
    else
        return(-1);

    printf1("Current memory: %d bytes.\n",MemReq);

    if (typ == 3 && EDAvail == 0) {
        p_err(-15,1);
        return(-1);
    }

    set_mldef();                /* set defaults for ML estimation */
                             
    if (parm(CmdBuf + n,12,1))   /* get parameters */
        goto FMLFin;

    if (PMFmtF == 0)            /* default print format */
        pmfmt(12,4);
           
    newline();

    if (FNArgN < 1) {
        printf1("Error: function should contain at least one argument.\n");
        goto FMLFin;
    }
    if (typ != 3 && PMFTYP5 > 0) {
        p_err(-40,1);
        goto FMLFin;
    }
    if (typ == 3)
        FVFlg = 2;

    if (FVFlg) {
        prn_feval();
        if (typ == 3) {
            FVFlg = 2;
            printf1("Using episode data");
            if (NSP > 0)  
                printf1(" (%d episode splits)",NSplits);
            printf1(".\n");
        }
        newline();
    }
    if (PMDOPT < 0)             /* default with derivatives */
        PMDOPT = 2;

    newline();
    set_mlopt();                /* adjust ML options */         
    NParm = FNArgN;             /* number of parameters */
           
    /* covariance matrix only calculated if function refers to data matrix
       variables and there are degrees of freedom. */

    /**************
    if (CCTyp > 2)
        CCTyp = 2;     
    *************/

    cov = 1;
    if (FVFlg == 0 || NParm >= NOC)
        CCTyp = cov = 0;

    if (ml_init(1,typ1,0))  /* init optimization */
        goto FMLFin;

    if (get_dsv(NParm,Par,0,ParLB,ParUB,1))   /* try to get starting values */
        goto FMLFin;

    /* print parameters and starting values */

    if (prn_fsval(NParm,Par + 1,0,ParLB + 1,ParUB + 1,0))
        goto FMLFin;

    if (fnd_alloc(1,2,FNArgN,FNPN,iflag))   /* add memory for derivatives */
        goto FMLFin;
               
    prot_init(3 + typ,UMOD,0);          /* init protocol file */

    if (PMNConS > 0) {                  /* process constraints */
        newline();
        if (con_proc(CmdBuf))
            goto FMLFin;
    }

    if (ffmin(UMOD,typ1,1,1))       /* minimization */
        goto FMLFin;                /* insuff memory or fn error */
                   
    prn_mlres(typ1);                /* print info about min algorithm */

    if (LConv < 0)
        cov = 0;
   
    gm_coeff(NParm,Par,cov,Diag);   /* print coefficients */
             
    if (LConv >= 0) {               /* if convergence */
        newline();
        prn_ml1res(typ1);           /* additional ML results */

        if (PMResFDef && FVFlg)     /* print residuals */
            gm_pres();
    }
    err = 0;

FMLFin:
    if (FVFlg)
        FVFlg = 1;
           
    con_free();
    ml_init(0,0,0);     
    fnd_alloc(0,0,0,0,0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gm_coeff(n,x,cov,diag)                                                  */
/*                                                                          */
/*  Print estimated coefficients x[i], i = 1,...,n.                         */
/*  If cov = 1  diag[] is main diagonal of cov matrix.                      */
/*  If cov = 2  diag[] is full cov matrix.                                  */
   
void gm_coeff(int n,double *x,int cov,double *diag)
{
    register int i,k,l;
    double tmp,tmp1;

    l = 9;
    if (l < FNMLen)
        l = FNMLen;

    printf1("Idx  Parameter ");
    prnchar(' ',l - 9,0);
    prnchar(' ',PMTFmt1 - 4,0); printf1("Value");
    prnchar(' ',PMTFmt1 - 4,0); printf1("Error");
    prnchar(' ',PMTFmt1 - 6,0); printf1("Value/E  Signif\n");
    prnchar('-',14 + l + 3 * (PMTFmt1 + 1),1);
    
    for (i = 1; i <= FNArgN; ++i) {
        k = FNArgSP[i - 1];
        printf1("%3d  %s  ",i,FNArgDef[k]);
        prnchar(' ',l - strlen(FNArgDef[k]),0);
    
        printf1(PMTFmtS,x[i]); 

        if (cov) {
            if (cov == 1)
                tmp = diag[i];
            else
                tmp = diag[(i - 1) * n + i];


            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                printf1(PMTFmtS,tmp);
                if (tmp > 0.0) {
                    tmp1 = x[i] / tmp;
                    tmp = 2.0 * cdnf(fabs(tmp1)) - 1.0;
                    printf1(PMTFmtS,tmp1);
                    printf1("%7.4lf",tmp);
                }
                else {
                    prnchar(' ',PMTFmt1 - 3,0);
                    printf1("---     ---"); 
                }
            }
            else {
                prnchar(' ',PMTFmt1 - 3,0); printf1("--- "); 
                prnchar(' ',PMTFmt1 - 3,0); printf1("---     ---"); 
            }
        }
        newline();           
    }
}

/* ------------------------------------------------------------------------ */
/*  get_func(s,opt,n3,iflag,msg)                                            */
/*                                                                          */
/*                      If opt != 0 get definition of function in s.        */
/*                      If opt == 0 free previously allocated memory.       */
/*                      If iflag != 0 allocate also FNArgVal1[].            */
/*                                                                          */
/*                      FNFlg = 1 if valid function definition.             */
/*                      FCFlg = 1 if valid constraint function definition.  */
/*                      FVFlg = 1 if references to data matrix variables.   */
/*                      FNMLen = max length of parameter strings.           */
/*                                                                          */
/*                      FNArgN = number of arguments. Sorted in FNArgDef,   */
/*                      meaning that FNArgSP[k] is the k.th argument.       */
/*                                                                          */
/*                      n3 is set if function contains type 3 operators     */
/*                      or refers to type 5 variables.                      */
/*                                                                          */
/*                      Introduced optional levels, see get_level().        */  
/*                                                                          */
/*                      Return 0 if OK, -1 if error.                        */

int get_func(char *s,int opt,int *n3,int iflag,char *msg)
{
    register int i;
    int err,np,fin,l,ll,n,nl,nv,nc,n2,fflag,cflag;
    register char c,*p,*q;

    FNNLev = fflag = err = np = 0;
    if (opt == 0)
        goto GFFin;

    printf1("\nFunction definition");
    if (msg != NULL)
        printf1(" (%s):\n",msg);
    else
        printf1(":\n");

    if (!strncmp(s,"level",5)) {    /* get level */
        s =  get_level(s + 5,&err);
        if (err)
            goto GFFin;
    }
    err = -1;
    p = s;
      
    cflag = fin = 0;
    FNMLen = FVFlg = FNFlg = FCFlg = 0;
    *n3 = 0;

    while (*p) {

        s = p;
        if (!strncmp(p,"fn=",3)) {
            fin = 1;
            p += 2;
            *p = '\0';
        }
        else if (!strncmp(p,"con=",4)) {
            fin = -1;
            p += 3;
            *p = '\0';
        }
        else {
            l = 0;
            while (*p && ((*p >= 'a' && *p <= 'z') || isdigit((int)*p))) {
                l++;
                p++;
            }
            if (l == 0 || *p != '=' || !*(p + 1)) {
                q = skip_expr(s);
                if (!*q) {
                    fin = 2;
                    p = s;
                    goto GFCONT;
                }
                goto GFFin;
            }
            *p = '\0';

            if (v_search(s,&nl,&n) >= 0) {
                err = -3;
                goto GFFin;
            }
            if (np >= MaxFNP) {
                err = -4;
                goto GFFin;
            }   
            if (!(FNPDef[np] = (char *)calloc(l + 1,sizeof(char)))) {
                err = -2;
                goto GFFin;   
            }
            FNPLen[np] = l;
            memrq(l + 1,sizeof(char));
            strcpy(FNPDef[np],s);
            FNPN++;
        }
        *p++ = '=';
GFCONT:
        q = skip_expr(p);
        c = *q;
        *q = '\0';

        prn_fexp(s,fin);       /* print expression */

        if (fin >= 0) {
            if (get_farg(p)) { /* check expression for undefined arguments */
                err = -5;
                goto GFFin;
            }
        }

        if ((n = v_parse(p,iflag)) < 0 || ESCnt <= 0) {
            printf1("Syntax error (%d) in expression: %s.\n",n,s);
            if (n < 0)
                prn_emsg1(n);
            err = -5;                       
            goto GFFin;
        }
        if (fin > 0) {
            if (!(FNETyp = (int *)calloc(ESCnt,sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(FNEVal = (double *)calloc(ESCnt,sizeof(double)))) {
                free((char *)FNETyp);
                err = -2;
                goto GFFin;   
            }
            FNECnt = ESCnt;
            memrq(FNECnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ESCnt; ++i) {
                FNETyp[i] = ESTyp[i];
                FNEVal[i] = ESVal[i];
            }

            /* check for unused arguments */

            /*******************************
            for (i = 0; i < FNArgN; ++i) {
                fflag = 0;
                for (j = 0; j < ESCnt; ++j) {
                    if (FNETyp[j] == FAOFFS + i) {
                        fflag = 1;
                        break;
                    }
                }
                if (fflag == 0) {
                    printf1("Error: at least one argument does not occur in function.\n");
                    err = -5;
                    goto GFFin;
                }
            }
            *****************/
            fflag = 1;
        }
        else if (fin < 0) {
            if (fflag == 0) {
                printf1("\nError: using constraints requires previous function definition.\n");
                goto GFFin;
            }
            if (!(FNECTyp = (int *)calloc(ESCnt,sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(FNECVal = (double *)calloc(ESCnt,sizeof(double)))) {
                free((char *)FNECTyp);
                err = -2;
                goto GFFin;   
            }
            FNECCnt = ESCnt;
            memrq(FNECCnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ESCnt; ++i) {
                FNECTyp[i] = ESTyp[i];
                FNECVal[i] = ESVal[i];
            }
            cflag = 1;
        }
        else {
            if (!(FNPETyp[np] = (int *)calloc(ESCnt,sizeof(int)))) {
                err = -2;
                goto GFFin;   
            }
            if (!(FNPEVal[np] = (double *)calloc(ESCnt,sizeof(double)))) {
                free((char *)FNPETyp[np]);
                err = -2;
                goto GFFin;   
            }
            FNPECnt[np] = ESCnt;
            memrq(ESCnt,sizeof(int) + sizeof(double));
            for (i = 0; i < ESCnt; ++i) {
                FNPETyp[np][i] = ESTyp[i];
                FNPEVal[np][i] = ESVal[i];
            }
            np++;
        }
        *q = c;
        if (*q != ',') {
            p = q;
            break;
        }
        p = q + 1;

    }
    if (*p) {        
        printf1("Error: definition must end with fn or con expression.\n");
        err = -5;
        goto GFFin;
    }
    if (fin == 0) {
        printf1("Error: need an fn=... expression.\n");
        err = -5;
        goto GFFin;
    }

    /* allocate memory for parameter values */

    if (FNArgN > 0) {
        if (!(FNArgVal = (double *)calloc(FNArgN,sizeof(double)))) {
            err = -2;
            goto GFFin;   
        }
        FNArgValA = FNArgN;
        memrq(FNArgN,sizeof(double));

        if (iflag) {
            if (!(FNArgVal1 = (double *)calloc(FNArgN,sizeof(double)))) {
                err = -2;
                goto GFFin;   
            }
            FNArgVal1A = FNArgN;
            memrq(FNArgN,sizeof(double));
        }
    }
   
    if (FNArgN > 0) {       /* sorting */

        if (!(FNArgSP = (short *)calloc(FNArgN,sizeof(short)))) { 
            err = -2;     
            goto GFFin;    
        }
        memrq(FNArgN,sizeof(short));
        FNArgSPA = FNArgN;

        if (!(FNArgSPI = (short *)calloc(FNArgN,sizeof(short)))) { 
            err = -2;     
            goto GFFin;    
        }
        memrq(FNArgN,sizeof(short));
        FNArgSPIA = FNArgN;

        for (i = 0; i < FNArgN; ++i)  
            FNArgSP[i] = i;

        qsort((char *)FNArgSP,FNArgN,sizeof(short),fcomp);

        for (i = 0; i < FNArgN; ++i)  
            FNArgSPI[FNArgSP[i]] = i;
    }

    /* check for data matrix variables */

    n2 = nc = nv = 0;
    for (i = 0; i < FNPN; ++i) {

        check_expr(FNPECnt[i],FNPETyp[i],&n,&nl,&l,&ll,0);
        nv += n;
        nc += nl;
        n2 += l;
        *n3 += ll;
    }
    check_expr(FNECnt,FNETyp,&n,&nl,&l,&ll,0);
    nv += n;
    nc += nl;
    n2 += l;
    *n3 += ll;
    if (nc > 0) {
        printf1("\nError: function may not contain ci references.\n");
        err = -5;  
        goto GFFin;
    }
    if (n2 > 0) {
        printf1("\nError: function may not contain type 2 operators or type 4 variables.\n");
        err = -5;  
        goto GFFin;
    }
    if (nv > 0)
        FVFlg = 1;
    else if (FNNLev) {
        printf1("\nError: multilevel function requires reference to data matrix.\n");
        err = -6; 
        goto GFFin;
    }
    FNFlg = 1;

    if (cflag) {        /* check constraint function */

        check_expr(FNECCnt,FNECTyp,&n,&nl,&l,&ll,0);

        if (n > 0 || nl > 0 || l > 0 || ll > 0) {
            printf1("\nError: constraint function may not contain ");
            if (n > 0)  
                printf1("data matrix variables.\n");
            else if (nl > 0)  
                printf1("ci references.\n");
            else if (l > 0)  
                printf1("type 2 operators.\n");
            else if (ll > 0)  
                printf1("type 3 operators.\n");
            err = -5;  
            goto GFFin;
        }
        FCFlg = 1;
    }
    if (FNNLev > 0) {       /* check for correct recognition of levels */
        n = FNNLev;
        for (i = 0; i < FNPN; ++i) {
            if (sscanf(FNPDef[i],"fn%d",&nl) == 1 && nl == n) {
                FNPELev[i] = n;
                n--;
            }
        }
        if (n != 0) {
            printf1("\nError: function definition inconsistent with levels.\n");
            printf1("Need: fn%d",FNNLev);
            for (i = FNNLev - 1; i > 0; --i)
                printf1(", fn%d",i);
            printf1(", and fn expression.\n");   
            err = -6;
            goto GFFin;
        }
    }
    err = 0;

GFFin:
    if (err == -1)
        printf1("Syntax error: %s\n",s);
    else if (err == -2)
        p_err(-2,1);
    else if (err == -3) {
        printf1("Error: %s\n",s);
        printf1("Parameter is reserved, already used, or not allowed.\n");
    }
    else if (err == -4)  
        printf1("Error: exceeded max number of intermediate parameters.\n");

    if (err)
        err = -1;

    if (opt == 0 || err) {

        for (i = 0; i < FNPN; ++i) {
            if (FNPLen[i] > 0) {
                free((char *)FNPDef[i]);
                memrq(-FNPLen[i] - 1,sizeof(char));
                FNPLen[i] = 0;
            }
            if (FNPECnt[i] > 0) {
                free((char *)FNPETyp[i]);
                free((char *)FNPEVal[i]);
                memrq(-FNPECnt[i],sizeof(int) + sizeof(double));
                FNPECnt[i] = 0;
            }
        }
        for (i = 0; i < FNArgN; ++i) {
            if (FNArgLen[i] > 0) {
                free((char *)FNArgDef[i]);
                memrq(-FNArgLen[i] - 1,sizeof(char));
                FNArgLen[i] = 0;
            }
        }
        if (FNECnt > 0) {
            free((char *)FNETyp);
            free((char *)FNEVal);
            memrq(-FNECnt,sizeof(int) + sizeof(double));
            FNECnt = 0;
        }
        if (FNECCnt > 0) {
            free((char *)FNECTyp);
            free((char *)FNECVal);
            memrq(-FNECCnt,sizeof(int) + sizeof(double));
            FNECCnt = 0;
        }
        if (FNArgSPA > 0) {
            free((char *)FNArgSP);
            memrq(-FNArgSPA,sizeof(short));
            FNArgSPA = 0;         
        }
        if (FNArgSPIA > 0) {
            free((char *)FNArgSPI);
            memrq(-FNArgSPIA,sizeof(short));
            FNArgSPIA = 0;         
        }
        if (FNArgValA > 0) {
            free((char *)FNArgVal);
            memrq(-FNArgValA,sizeof(double));
            FNArgValA = 0;         
        }
        if (FNArgVal1A > 0) {
            free((char *)FNArgVal1);
            memrq(-FNArgVal1A,sizeof(double));
            FNArgVal1A = 0;         
        }
        if (FNLVarA > 0) {
            free((char *)FNLVar);
            memrq(-FNLVarA,sizeof(int));
            FNLVarA = 0;         
        }
        FNNLev = FNPN = FNArgN = FNFlg = FCFlg = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_level(s,err)    scan level expression.                              */
/*                                                                          */
/*                      level(1+VName,2+VName,...)                          */
/*                      there might be * instead of +                       */
/*                                                                          */
/*                      FNNLev = number of level                            */
/*                      FNLVar = index of variable                          */
/*                                                                          */
/*  Return err = 0 if OK, -1 if syntax error, -2 if insufficient memory,    */
/*         err = -9 if undefined variable.                                  */
/*                                                                          */
/*  Return pointer to next character in string s.                           */
   
char *get_level(char *s,int *err)  
{
    register int i;
    int n;
    register char *p,*q;
    char vname[VNLMax + 1];

    *err = -1;
    p = s; 
    if (*p++ != '(')
        return(s);
    q = p;
    FNNLev = 1;
    while (*q && *q != ')') {
        if (*q++ == ',')
            FNNLev++; 
    }
    if (*q++ != ')' || *q++ != ',')
        return(s);

    if (!(FNLVar = (int *)calloc(FNNLev,sizeof(int)))) {
        *err = -2;
        return(s);
    }
    memrq(FNNLev,sizeof(int));
    FNLVarA = FNNLev;

    for (i = 1; i <= FNNLev; ++i) {

        if (sscanf(p,"%d",&n) != 1 || n != i)
            return(s);

        p = skip_int(p);

        if (*p++ != '=')
            return(s);

        n = get_vidx1(p,vname);
        if (n < 0) {
            printf1("Undefined variable name in level expression.\n");
            *err = -9;
            return(s);
        }
        FNLVar[i - 1] = n;
        p += strlen(vname);
        if (*p != ',' && *p != ')')
            return(s);
        p++;
    }
    newline();
    for (i = 0; i < FNNLev; ++i)  
        printf1("Level%3d defined by: %s\n",i + 1,VName[FNLVar[i]]);
    newline();
    *err = 0;
    return(q);
}

/* ------------------------------------------------------------------------ */
/*  prn_fexp(s,opt)     Print function expression.                          */

void prn_fexp(char *s,int opt)  
{
    register int l;
    register char *p;

    if (opt == 2)
        printf1("fn      = ");
        
    l = 8;
    p = s;
    while (*p) {
        if (*p == '=') {
            while (l-- > 0)
                printf1(" ");
            printf1("= ");
        }
        else
            printf1("%c",*p);
        if (l > 0)
            l--;
        p++;
    }
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  get_farg(s)     Check expression s for new arguments. Save arguments    */
/*                  in FNArgLen[] and FNArgDef[]. Count FNArgN.             */
/*                  Max number of arguments is MaxP.                        */
/*                  Return 0 if OK, -1 if error (exceeded MaxP, or if       */
/*                  insufficient memory).                                   */

int get_farg(char *s)          
{
    register int l;
    register char c,*p;
    int n,nl,len; 

    while (*s) {

        if (check_vname(s)) {       /* check for variable */
            if (v_search1(s,&len) >= 0)
                s += len;
        }         
        if (!*s)
            break;
        p = s;
        l = 0;

        if (*p >= 'a' && *p <= 'z') {

            while (*p && ((*p >= 'a' && *p <= 'z') || isdigit((int)*p))) {
                l++;
                p++;
            }
            if (*p == '(')  
                l = 0;

            if (l > 0) {
                c = *p;
                *p = '\0';
                if (v_search(s,&nl,&n) < 0) {
                    if (FNArgN >= MaxP) {
                        printf1("Error: exceeded max number of arguments.\n");
                        return(-1);
                    }
                    if (!(FNArgDef[FNArgN] = (char *)calloc(l + 1,sizeof(char)))) {
                        p_err(-2,1);
                        return(-1);
                    }
                    FNArgLen[FNArgN] = l;
                    memrq(l + 1,sizeof(char));
                    strcpy(FNArgDef[FNArgN],s);
                    FNArgN++;
                    if (FNMLen < l)
                        FNMLen = l;
                }
                *p = c;
            }
            else
                p++;
        }
        else
            p++;
        s = p;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  fcomp()    compare function                                             */

int fcomp(const void *arg1,const void *arg2)
{     
    char *p1,*p2;

    p1 = FNArgDef[*(short *)arg1];
    p2 = FNArgDef[*(short *)arg2];
    return(strcmp(p1,p2));

}

/* ------------------------------------------------------------------------ */
/*  prn_feval()     print info about function evaluation.                   */

void prn_feval(void)
{
    register int i,ix,n;
    double tmp,tmp1;

    if (FVFlg) {
        printf1("Function evaluation: ");
        if (FNNLev == 0)
            printf1("sum over %d data matrix cases.\n",NOC);
        else {
            printf1("according to %d level(s) in data matrix.\n",FNNLev);
            n = 0;
            ix = FNLVar[0];
            tmp1 = get_data(ix,0) - 1.0;
            for (i = 0; i < NOC; ++i) {
                tmp = get_data(ix,i);
                if (tmp != tmp1) {
                    n++;
                    tmp1 = tmp;
                }
            }
            printf1("Number of level 1 units: %d\n",n);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  get_flval(f,n,x,deriv,opt,g,h,d)                                        */
/*  ##                                                                      */
/*                      get value from function expression.                 */  
/*                      if FVFlg == 1 sum over all data matrix cases.       */
/*                      if FVFlg == 2 use episode data.                     */
/*                                                                          */
/*                      Function arguments: x[i], i = 0,...,n-1.            */
/*                      Return value in f.                                  */
/*                                                                          */
/*  If deriv >= 1  calculate gradient.      FNGrad[0][...]                  */
/*  If deriv == 2  calculate hessian.       FNHess[0][...]                  */
/*  If CCovFlg ==1 use outer product of gradients instead of hessian        */
/*                 (only if reference to data matrix)                       */
/*  If CGradFlg == 1 calculate gradients in MatVal[MPGradIdx]               */
/*                                                                          */
/*  opt = 0 : do not use g,h,d                                              */
/*  opt = 1 : save gradient in g[1,...,n], full hessian in h[1,...],        */
/*            and do not use d.                                             */
/*  opt = 2 : save gradient in g[1,...,n], lower triangle of hessian in     */
/*            h[1,...] and diagonal of hessian in d[1,...,n].               */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_flval(double *f,int n,double *x,int deriv,int opt,double *g,double *h,double *d)
{
    register int i,j,k,l;
    int err,nn,ii,nii,deriv1,sn,org,des,spl,nspl;
    double tmp,ts,tf;

    NOCUsed = err = 0;
    for (i = 0; i < n; ++i)         /* set values for evaluation */
        FNArgVal[i] = x[i];

    if (FVFlg == 0) {               /* no reference to data matrix, evaluate
                                       expression */
        err = get_fival(0,f,deriv,0);

        if (err == 0 && deriv && opt) {

            for (i = 1; i <= n; ++i) {              /* save gradient */
                g[i] = FNGrad[0][i - 1];

                if (CGradFlg && NOCUsed < MPGradRow)
                    MatVal[MPGradIdx][i] = g[i];        
            }

            if (deriv > 1) {                        /* save hessian */
                if (opt == 1) {
                    k = 0;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = FNHess[0][k++];
                            h[(i - 1) * n + j] = tmp;
                            if (j < i)
                                h[(j - 1) * n + i] = tmp;
                        }
                    }
                }
                else if (opt == 2) {

                    k = 0;
                    l = 1;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            tmp = FNHess[0][k++];
                            if (j < i)
                                h[l++] = tmp;
                            else
                                d[i] = tmp;
                        }
                    }
                }
            }
        }
        NOCUsed = 1;
        return(err);
    }

    /*  reference to data matrix. Sum over all level 1 units */

    *f = 0.0;
    if (deriv && opt) {
        for (i = 1; i <= n; ++i)
            g[i] = 0.0;
        if (deriv > 1) {
            if (opt == 1) {
                nn = n * n;
                for (i = 1; i <= nn; ++i)
                    h[i] = 0.0;
            }
            else if (opt == 2) {
                for (i = 1; i <= n; ++i)
                    d[i] = 0.0;
                nn = n * (n - 1) / 2;
                for (i = 1; i <= nn; ++i)
                    h[i] = 0.0;
            }   
        }
    }
    deriv1 = deriv;
    if (deriv1 == 2 && CCovFlg != 0)
        deriv1 = 1;

    if (FVFlg == 2)
        get_spell(1,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl);
    else
        ii = 0;

    while (1) {
    
        if (FVFlg == 2) {
            if (get_spell(0,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl) == 0) 
                break;
        }
        err = get_filval(ii,&tmp,deriv1,&nii);

        if (err)
            break;

        *f += tmp;

        if (deriv && opt) {

            for (i = 1; i <= n; ++i) {          /* save gradient */
                g[i] += FNGrad[0][i - 1];

                if (CGradFlg && NOCUsed < MPGradRow)
                    MatVal[MPGradIdx][NOCUsed * MPGradCol + i] = FNGrad[0][i - 1];
            }
            if (deriv > 1) {                    /* save hessian */
                if (opt == 1) {
                    k = 0;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {
                            if (CCovFlg == 0)
                                tmp = FNHess[0][k++];
                            else
                                tmp = FNGrad[0][i - 1] * FNGrad[0][j - 1];

                            h[(i - 1) * n + j] += tmp;
                            if (j < i)
                                h[(j - 1) * n + i] += tmp;
                        }
                    }
                }
                else if (opt == 2) {

                    k = 0;
                    l = 1;
                    for (i = 1; i <= n; ++i) {
                        for (j = 1; j <= i; ++j) {

                            if (CCovFlg == 0)
                                tmp = FNHess[0][k++];
                            else
                                tmp = FNGrad[0][i - 1] * FNGrad[0][j - 1];
                            if (j < i)
                                h[l++] += tmp;
                            else
                                d[i] += tmp;
                        }
                    }
                }
            }
        }
        if (CGradFlg && NOCUsed < MPGradRow && PM2NV > 0)
            mp_putvar(NOCUsed,MPGradCol,MatVal[MPGradIdx],PM2NV,PM2VIdx,ii);

        NOCUsed++;
        if (FVFlg != 2) {
            ii = nii;
            if (ii >= NOC)
                break;
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_filval(n,x,deriv,ni)    get value from multilevel function for      */  
/*                              level 1 unit beginning with data            */
/*  ##                          matrix case n.                              */
/*                              Return value in x, return number of next    */
/*                              data matrix case in ni.                     */
/*                                                                          */  
/*  If deriv >= 1  Gradient is in FNGrad[0][...]                            */
/*  If deriv == 2  Hessian is  in FNHess[0][...]                            */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_filval(int n,double *x,int deriv,int *ni)
{
    register int i,j,l;
    int err,last;
    int first[100];
    double ltmp,lcase[100];

    err = 0;
    for (l = 0; l < FNNLev; ++l) {
        first[l] = 1;
        lcase[l] = get_data(FNLVar[l],n);
    }
    l = FNNLev - 1;

FICONT1:
    for (i = 0; i < FNPN; ++i) {

        err = v_eval1(n,FNPECnt[i],FNPETyp[i],FNPEVal[i],ESIdx,x,0,deriv,0,0,0);
        if (err)  
            return(err);

        if (l < 0 || FNPELev[i] != l + 1) {  
            FNPVal[i] = *x;
            if (deriv > 0) {        /* copy gradient and hessian */

                for (j = 0; j < FNGradL; ++j)
                    FNPGrad[i][j] = FNGrad[0][j];

                if (deriv > 1) {    
                    for (j = 0; j < FNHessL; ++j)
                        FNPHess[i][j] = FNHess[0][j];
                }
                FNPEATyp[i] = FNEVATyp;
            }
        }
        else {
            if (first[l]) {
                FNPVal[i] = *x;
                if (deriv > 0) {        /* copy gradient and hessian */

                    for (j = 0; j < FNGradL; ++j)
                        FNPGrad[i][j] = FNGrad[0][j];

                    if (deriv > 1) {    
                        for (j = 0; j < FNHessL; ++j)
                            FNPHess[i][j] = FNHess[0][j];
                    }
                    FNPEATyp[i] = FNEVATyp;
                }
                first[l] = 0;
            }
            else {   
                FNPVal[i] += *x;
                if (deriv > 0) {        /* copy gradient and hessian */

                    for (j = 0; j < FNGradL; ++j)
                        FNPGrad[i][j] += FNGrad[0][j];

                    if (deriv > 1) {    
                        for (j = 0; j < FNHessL; ++j)
                            FNPHess[i][j] += FNHess[0][j];
                    }
                    FNPEATyp[i] = FNEVATyp;
                }
            }
            last = 0;
            if (n + 1 >= NOC)
                last = 1;   
            else {          
                for (j = l; j >= 0; --j) {
                    ltmp = get_data(FNLVar[j],n + 1);
                    if (ltmp != lcase[j]) {
                        first[l] = 1;
                        last = 1;
                        break;
                    }
                }
            }
            if (last == 0) {
                n++;
                for (j = 0; j < FNNLev; ++j) {
                    lcase[j] = get_data(FNLVar[j],n);
                    if (j > l)
                        first[j] = 1;
                }
                l = FNNLev - 1;
                goto FICONT1;
            }
            l--;
        }
    }

    *ni = n + 1;

    /* evaluate final function expression */

    err = v_eval1(n,FNECnt,FNETyp,FNEVal,ESIdx,x,0,deriv,0,0,0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_fival(n,x,deriv,aflag)                                              */
/*                                                                          */
/*  Get value from function expression.                                     */  
/*  if FVFlg != 0 use data matrix case n.                                   */
/*  if aflag != 0 use AVVAL in call of v_eval1().                           */
/*  Return value in x.                                                      */
/*                                                                          */  
/*  If deriv >= 1  Gradient is in FNGrad[0][...]                            */
/*  If deriv == 2  Hessian is  in FNHess[0][...]                            */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_fival(int n,double *x,int deriv,int aflag)
{
    register int i,j;
    int err = 0;

    /* first evaluate intermediate parameters */

    for (i = 0; i < FNPN; ++i) {
        if (FNPECnt[i] > 0) {

            err = v_eval1(n,FNPECnt[i],FNPETyp[i],FNPEVal[i],ESIdx,x,aflag,deriv,0,0,0);
            if (err)  
                return(err);
            FNPVal[i] = *x;

            if (deriv > 0) {        /* copy gradient and hessian */

                for (j = 0; j < FNGradL; ++j)
                    FNPGrad[i][j] = FNGrad[0][j];

                if (deriv > 1) {    
                    for (j = 0; j < FNHessL; ++j)
                        FNPHess[i][j] = FNHess[0][j];
                }
                FNPEATyp[i] = FNEVATyp;
            }
        }
    }

    /* then evaluate function */

    err = v_eval1(n,FNECnt,FNETyp,FNEVal,ESIdx,x,aflag,deriv,0,0,0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  prn_fsval(n,x,gmina,lb,ub,opt)                                          */
/*                                                                          */
/*  Print parameters and starting values: x[i], i = 0,...,n-1               */                    
/*  If gmina != 0 print also bounding box. If opt != 0 do not print         */
/*  parameters.                                                             */  
/*                                                                          */
/*  Return 0 if OK, -1 if error in bounds.                                  */

int prn_fsval(int n,double *x,int gmina,double *lb,double *ub,int opt)
{
    register int i,j;
    int l,err;

    err = 0;
    l = 10;
    if (l < FNMLen)
        l = FNMLen;

    printf1("Idx  Parameter ");
    if (opt == 0) {
        prnchar(' ',l - 9,0);
        printf1(" Starting value");
    }
    if (gmina)
        printf1("      Lower bound     Upper bound");
    newline();
    for (i = 0; i < n; ++i) {
        j = FNArgSP[i];
        printf1("%3d  %s  ",i + 1,FNArgDef[j]);
        prnchar(' ',l - strlen(FNArgDef[j]),0);
        if (opt == 0)
            printf1("%15.8e ",x[i]);
        if (gmina) {
            printf1("%15.8e ",lb[i]);
            printf1("%15.8e ",ub[i]);
            if (opt == 0) {
                if (lb[i] > x[i] + EPSI || ub[i] < x[i] - EPSI)
                    err = -1;
            }
            else {
                if (lb[i] >= ub[i] + EPSI)
                    err = -1;
            }
        }
        printf1("\n");
    }
    newline();
    if (err)  
        printf1("Error in bounding box.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gm_pres()   print residuals to PMResFd.                                 */

void gm_pres(void)           
{
    register int i,j,k;
    int err;
    double tmp;

    if (PMResFDef == 0)
        return;

    for (i = 0; i < FNArgN; ++i)         /* set values for evaluation */
        FNArgVal[i] = Par[i + 1];
  
    for (i = 0; i < NOC; ++i) {
        fprintf(PMResFd,"%6d ",i + 1);
        err = get_fival(i,&tmp,0,0);
        if (err) {
            fprintf(PMResFd,"error\n");
            continue;
        }
        fprintf(PMResFd,PMFmtS,tmp);
        for (j = 0; j < PMNV; ++j) {   /* add variables */
            k = PMVIdx[j];
            tmp = get_data(k,i);
            fprintf(PMResFd,VPFmtS[k],tmp);
        }
        fprintf(PMResFd,"\n");
    }
    printf1("Residuals written to: %s\n",PMResFName);
}

/* ------------------------------------------------------------------------ */
/*  get_ifval(fl,fu,n,xl,xu,deriv,gl,gu)                                    */
/*  ##                                                                      */
/*  Calculates an inclusion function. Function arguments are                */
/*  [xl[i],xu[i]] for i = 0,...,n - 1. Return interval in [fl,fu].          */
/*                                                                          */
/*  If deriv = 1 calculate also gradiend bounds in gl[i],gu[i], i=0,n-1     */
/*                                                                          */
/*  if FVFlg == 1 sum over all data matrix cases.                           */
/*  if FVFlg == 2 use episode data.                                         */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_ifval(double *fl,double *fu,int n,double *xl,double *xu,int deriv,
    double *gl,double *gu) 
{
    register int i;
    int err,ii,sn,org,des,spl,nspl;
    double ts,tf,xx[2];

    err = 0;
    for (i = 0; i < n; ++i) {       /* set values for evaluation */
        FNArgVal[i] = xl[i];
        FNArgVal1[i] = xu[i];
    }
    *fl = *fu = 0.0;

    if (FVFlg == 0) {

        err = get_ifival(0,xx,deriv);
        if (err == 0) {
            *fl = xx[0];
            *fu = xx[1];

            if (deriv) {
                for (i = 0; i < n; ++i) {              /* save gradient */
                    gl[i] = FNGrad[0][i];
                    gu[i] = FNGrad[1][i];
                }
            }
        }
    }
    else {                      /* sum over data matrix cases */
                                /* or use episode data */
        if (deriv) {
            for (i = 0; i < n; ++i)            
                gu[i] = gl[i] = 0.0;           
        }
        if (FVFlg == 2)
            get_spell(1,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl);
        else
            ii = 0;

        while (1) {
    
            if (FVFlg == 2) {
                if (get_spell(0,&ii,&sn,&org,&des,&ts,&tf,&spl,&nspl) == 0) 
                    break;
            }
            err = get_ifival(ii,xx,deriv);
            if (err)
                break;

            *fl += xx[0];
            *fu += xx[1];

            if (deriv) {
                for (i = 0; i < n; ++i) {              /* save gradient */
                    gl[i] += FNGrad[0][i];
                    gu[i] += FNGrad[1][i];
                }
            }
    
            if (FVFlg == 1) {
                if (++ii >= NOC)
                    break;
            }
        }
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_ifival(n,x,deriv)                                                   */
/*                                                                          */
/*  Get value from function expression, using intervals. If FVFlg != 0      */
/*  use data matrix case n. If deriv = 1 calculate also gradient range.     */
/*                                                                          */
/*  Return lower bound: x[0]                                                */
/*         upper bound: x[1]                                                */
/*                                                                          */  
/*  If deriv = 1 then Gradient lower bound in FNGrad[0][...]                */
/*                             upper bound in FNGrad[1][...]                */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_ifival(int n,double *x,int deriv)
{
    register int i,j;
    int err = 0;

    /* first evaluate intermediate parameters */

    for (i = 0; i < FNPN; ++i) {
        if (FNPECnt[i] > 0) {
            err = v_eval1(n,FNPECnt[i],FNPETyp[i],FNPEVal[i],ESIdx,x,0,deriv,0,0,1);
            if (err)  
                return(err);
            FNPVal[i] = x[0];
            FNPVal1[i] = x[1];

            if (deriv) {        /* copy gradient */

                for (j = 0; j < FNGradL; ++j) {
                    FNPGrad[i][j] = FNGrad[0][j];
                    FNPGrad1[i][j] = FNGrad[1][j];
                }
                FNPEATyp[i] = FNEVATyp;
            }
        }
    }
   
    /* then evaluate function */
           
    err = v_eval1(n,FNECnt,FNETyp,FNEVal,ESIdx,x,0,deriv,0,0,1);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_icval(fl,fu,n,xl,xu)                                                */
/*  ##                                                                      */
/*  Calculates an inclusion function for a constraint function.             */  
/*  Function arguments are [ xl[i], xu[i] ], i = 0,...,n - 1.               */
/*  Return interval in [fl,fu].                                             */
/*                                                                          */
/*  Return 0 if ok, else error flag from v_eval()                           */

int get_icval(double *fl,double *fu,int n,double *xl,double *xu)           
{
    register int i;
    int err;
    double xx[2];

    for (i = 0; i < n; ++i) {       /* set values for evaluation */
        FNArgVal[i] = xl[i];
        FNArgVal1[i] = xu[i];
    }

    /* first evaluate intermediate parameters */

    for (i = 0; i < FNPN; ++i) {
        if (FNPECnt[i] > 0) {
            err = v_eval1(0,FNPECnt[i],FNPETyp[i],FNPEVal[i],ESIdx,xx,0,0,0,0,1);
            if (err)  
                return(err);
            FNPVal[i] =  xx[0];
            FNPVal1[i] = xx[1];
        }
    }
   
    /* then evaluate constraint function */
           
    err = v_eval1(0,FNECCnt,FNECTyp,FNECVal,ESIdx,xx,0,0,0,0,1);
    if (err == 0) {
        *fl = xx[0];
        *fu = xx[1];
    }
    return(err);
}

