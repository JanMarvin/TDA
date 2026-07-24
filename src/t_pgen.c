/****************************************************************************/
/*  t_pgen                                                                  */
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
#include "t_parm.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_var.h"
#include "t_cdf.h"
#include "t_qrmod.h"
#include "t_gf.h"

/*  functions in t_pgen.c */

void p_err(int typ,int n);
void ps_err(int typ,char *txt,int n);
void p_warn(int typ,int n);
void p_serr(char *s,int n);
void p_wmsg(int typ,char *fname,int n);
void p_vterr(int typ);
int check_cmd(int opt);  
void p_cmd(char *pcmd,int nw);           
void prn_sve(void);
int eval_sve(int i);
int getd1(double *x,int ix,int ig,int opt);
int getd2(double *x,double *y,int ix,int iy,int opt);
int getdxy(double *x,double *y,int nx,int nif,int opt);
int check_nvar(int opt);
void prn_nwvar(int opt);
void prn_data(int m,int nc,int n,double *x,FILE *fd,char *fmt);
void prn1_data(double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt);
void prn1_coeff(int n,double *x,double *se,int ni,int df,short *vidx,int seflg);
void makefmt(int *fmt1,int *fmt2,char *fmts,int mlen,char sepc,int opt);
void makenfmt(int *fmt,char *fmts,int mlen,char sepc);
void pmfmt(int n,int m);
void pmtfmt(int n,int m);
void prn_cwt(void);        
void prn_sfmt(char *s,int n,char *fmt,double x);
void prn_f1mat(FILE *fd,int prn,char *fmt,int n,double *x);

/* ------------------------------------------------------------------------ */
/*  p_err(typ,n)    Print error message, add n newlines                     */

void p_err(int typ,int n)
{
    switch (typ) {
        case -1:    printf1("Syntax error.");
                    break;
        case -2:    printf1("Insufficient memory.");
                    break;
        case -3:    printf1("Error: no relational/graph data defined.");
                    break;
        case -4:    printf1("Error: syntax error or undefined variables.");
                    break;
        case -5:    printf1("Error: exceeded max number of variables.");
                    break;
        case -6:    printf1("Error: no variables.");
                    break;
        case -7:    printf1("Error: can't open or read input file.");
                    break;
        case -8:    printf1("Error: command needs at least two variables.");
                    break;
        case -9:    printf1("Error: no data archive defined.");
                    break;
        case -10:   printf1("Error: no PostScript output file.");
                    break;
        case -11:   printf1("Error: exceeded max number of sequence data structures.");
                    break;
        case -12:   printf1("Error: need sequence 0 data structure.");
                    break;
        case -13:   printf1("Error: no sequence data defined.");
                    break;
        case -14:   printf1("Error: exceeded max number of model parameters.");
                    break;
        case -15:   printf1("Error: no episode data defined.");
                    break;
        case -16:   printf1("Error: keep and drop options not compatible.");
                    break;
        case -17:   printf1("Error: need time points/periods.");
                    break;
        case -18:   printf1("Error: no data matrix (command ignored).");
                    break;
        case -19:   printf1("Error: no patterns defined.");
                    break;
        case -20:   printf1("Error in pattern definition.");
                    break;
        case -21:   printf1("Pattern definition inconsistent with selected state space.");
                    break;
        case -22:   printf1("Error: exceeded max number of string variables.");
                    break;
        case -23:   printf1("Error: insufficient memory for string variables.");
                    break;
        case -24:   printf1("Error: need an rc parameter.");
                    break;
        case -25:   printf1("Error: need even number of points.");
                    break;
        case -26:   printf1("Error: need at least two points.");
                    break;
        case -27:   printf1("Error: need an output file.");
                    break;
        case -28:   printf1("Error: selection results in zero cases.");
                    break;
        case -29:   printf1("Error: must be only a single wave.");
                    break;
        case -30:   printf1("Error: number of cases is less than number of variables.");
                    break;
        case -31:   printf1("Error: need at least two groups.");
                    break;
        case -32:   printf1("Syntax error in function definition.");
                    break;
        case -33:   printf1("Error: time periods must begin with zero.");
                    break;
        case -34:   printf1("Syntax error in model definition.");
                    break;
        case -35:   printf1("Error: exceeded maximum number of y terms.");
                    break;
        case -36:   printf1("Error: y terms must be recursive.");
                    break;
        case -37:   printf1("Error: reached max level of command files.");
                    break;
        case -38:   printf1("Error: need multi-episode data.");
                    break;
        case -39:   printf1("Error in creating episode data.");
                    break;
        case -40:   printf1("Error: cannot use type 5 variables (type 3 operators).");
                    break;
        case -41:   printf1("Error: reference to cj terms not possible.");
                    break;
        case -42:   printf1("Error: varlist must not contain brackets.");
                    break;
        case -43:   printf1("Error: exceeded max length of command buffer.");
                    break;
        case -44:   printf1("Error: in reading value labels.");
                    break;
        case -45:   printf1("Error: in number of coordinates.");
                    break;
        default:    printf1("Undefined error.");
                    break;
    }
    while (n-- > 0)
        newline();
}

/* ------------------------------------------------------------------------ */
/*  ps_err(typ,txt,n)   Print error message, add n newlines                 */

void ps_err(int typ,char *txt,int n)
{
    switch (typ) {
        case -1:    printf1("Error: can't open: %s",txt);
                    break;

        default:    printf1("Undefined error.");
                    break;
    }
    while (n-- > 0)
        newline();
}

/* ------------------------------------------------------------------------ */
/*  p_serr(s,n)     Print error message, add n newlines                     */

void p_serr(char *s,int n)
{
    printf1("Error: %s",s);
    while (n-- > 0)
        newline();
}

/* ------------------------------------------------------------------------ */
/*  p_warn(typ,n)   Print warning message, add n newlines                   */

void p_warn(int typ,int n)
{
    printf1("Warning: ");
    switch (typ) {
        case -1:    printf1("sel parameter ignored.");
                    break;
        case -2:    printf1("procedure is incompatible with episode splitting.");
                    break;
        case -3:    printf1("varlist contains type 5 variables.");
                    break;
        default:    printf1("Warning ...");
                    break;
    }
    while (n-- > 0)
        newline();
}

/* ------------------------------------------------------------------------ */
/*  p_wmsg(typ,fname,n)   print write message, and add n newlines.          */

void p_wmsg(int typ,char *fname,int n)
{
    switch (typ) {
        case  1:    printf1("Data written to");
                    break;
        case  2:    printf1("Covariance matrix written to");
                    break;
        case  3:    printf1("Estimated values written to");
                    break;
        default:    printf1("... ");
                    break;
    }
    printf1(": %s",fname);
    while (n-- > 0)
        newline();
}

/* ------------------------------------------------------------------------ */
/*  p_vterr(typ)    Print vtyp error.                                       */

void p_vterr(int typ)
{
    printf1("Error: cannot use type %d variables.\n",typ);
}

/* ------------------------------------------------------------------------ */
/*  check_cmd(opt)      Print command in CmdBuf and check syntax.           */
/*  ##                  If opt = 1 print full command, if 0 abbreviated.    */  
/*                      If opt = 2 print cmd()...                           */
/*                      Return 0 if OK, -1 if error.                        */

int check_cmd(int opt)   
{
    int err = -1;
    register char *p = CmdBuf;

    if (opt == 1)
        printf1("%s",p);

    while (*p && *p != '(' && *p != '=') {
        if (opt == 0 || opt == 2)
            printf1("%c",*p);
        p++;
    }        
    if (*p == '(') {
        if (opt != 1)
            printf1("(...");
        p = skip_blev(p);
        if (*(p - 1) != ')')  
            goto PRNCMDFin;
        if (opt != 1)
            printf1(")");
    }
    if (*p) {
        if (opt == 0)
            printf1("%s",p);
        else if (opt == 2)
            printf1("=...");

        if (*p++ != '=' || !*p)
            goto PRNCMDFin;
    }
    printf1("\n");
    err = 0;

PRNCMDFin:
    if (err) {
        printf1("\n");
        printf1("Syntax error.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  p_cmd(pcmd,nw)      print pcmd. if nw > 0 assume right-hand side is a   */
/*                      list of variables for nw panel waves.               */  
/*  ##                                                                      */

void p_cmd(char *pcmd,int nw)            
{
    register int i;
    int l,nf,n,r,nn;
    register char c,*p,*q,*pp,*qq;

    r = nf = 0;
    p = pcmd;
    printf1("Command: ");           
    l = strlen(pcmd);
    if (l < 60)  
        goto PCMDFin;
           
    q = p;
    while (*q) {
        if (*q == '=' || *q == '(')
            break;
        q++;
    }
    if (!*q)
        goto PCMDFin;

    c = *++q;
    *q = '\0';
    printf1("%s\n",p);        
    *q = c;

    while (1) {
        p = q;
        if (!strncmp(p,"xa(",3) || !strncmp(p,"xb(",3) || 
            !strncmp(p,"xc(",3) || !strncmp(p,"xd(",3))   
            q = skip_xa(p);
           
        else if (!strncmp(p,"fn=",3)) 
            q = skip_expr(p + 3);   
        else if (sscanf(p,"y%d",&n) == 1) {
            q = skip_int(p + 1);
            q = skip_expr(q + 1);
        }
        else if (!strncmp(p,"grp=",4)) {
            q = p + 3;
            while (*(q + 1) && (l = get_vnlen(q + 1)) > 0)  
                q += l + 1;
        }
        else {
            pp = p;
            while (*pp) {
                if (nw > 0) {
                    if (sscanf(pp,"nw=%d",&nn) == 1) {    
                        if (nw < nn)
                            nw = nn;
                    }   
                }
                q = skip_nc(pp);

                if (r > 0 && nw > 0) {
                    if (nw == 1) {
                        prnchar(' ',11,0);
                        goto PCMDFin;
                    }
                    for (i = 1; i < nw; ++i)
                        q = skip_nc(q + 1);
                    break;
                }
                qq = q + 1;
                if (isdigit((int)*qq) || *qq == '-' || *qq == '+')
                    pp = q + 1;
   
                else if (*q == ')' && *qq == '=') {
                    q = qq;
                    break;
                }
                else if (*qq == ')' && *++qq == '=') {
                    q = qq;
                    break;
                }
                else
                    break;
            }
        }
        if (*q == ',') {
            *q = '\0';
            printf1("           %s,\n",p);
            *q++ = ',';
            nf++;
        }
        else if (*q == '=' && strlen(q) > 50) {
            *q = '\0';
            printf1("           %s=\n",p);
            *q++ = '=';
            nf++;
            r++;
        }
        else {   
            prnchar(' ',11,0);
            goto PCMDFin;
        }
    }

PCMDFin:
    printf1("%s\n",p);
    /* if (nf)  */
    newline();
}

/* ------------------------------------------------------------------------ */
/*  prn_sve()       print info about case selection.                        */

void prn_sve(void)
{
    if (SVEFlg)  
        printf1("Case selection: %s\n",SVESTR);
}

/* ------------------------------------------------------------------------ */
/*  eval_sve(i)     evaluate sel expression for case i. return 1 if         */
/*                  case to be selected, otherwise zero.                    */

int eval_sve(int i)
{
    int r;
    double x;

    if (SVEFlg == 0)
        return(1);

    r = v_eval1(i,SVECnt,SVETyp,SVEVal,ESIdx,&x,0,0,0,0,0);             

    if (r) {
        printf1("Error (%d): can't evaluate sel expression in case %d (ignored).\n",r,i+1);
        /** prn_emsg2(r); **/
        return(0);
    }
    if (fabs(x) > EPSI1)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  getd1(x,ix,ig,opt)      get data for var ix into x[].                   */  
/*                          recognize select cases option.                  */
/*                          if ig >= 0 take this as a grouping variable     */
/*                          and use only cases where value is != 0.         */
/*                                                                          */
/*                          if opt print cases and error message.           */
/*                                                                          */
/*      Return number of cases, or -1 if error.                             */

int getd1(double *x,int ix,int ig,int opt)
{
    register int i,j;
    double tmp;

    j = 0;
    for (i = 0; i < NOC; ++i) {
        if (eval_sve(i) == 0)   
            continue;

        if (ig >= 0) {
            tmp = get_data(ig,i);
            if (fabs(tmp) < EPSI1)
                continue;
        }
        x[j] = get_data(ix,i);       
        j++;  
    }
    if (opt) {
        if (SVEFlg)  
            printf1("Number of selected cases: %d\n",j);
        if (j == 0)  
            p_err(-28,1);  
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  getd2(x,y,ix,iy,opt)    get data for var ix and iy into x[] and y[].    */  
/*                          if opt print cases and error message.           */

int getd2(double *x,double *y,int ix,int iy,int opt)
{
    register int i,j;

    j = 0;
    for (i = 0; i < NOC; ++i) {
        if (eval_sve(i) == 0)   
            continue;
       
        x[j] = get_data(ix,i);       
        y[j] = get_data(iy,i);       
        j++;  
    }
    if (opt) {
        if (SVEFlg)  
            printf1("Number of selected cases: %d\n",j);
        if (j == 0)  
            p_err(-28,1);  
    }
    return(j);
}

/* ------------------------------------------------------------------------ */
/*  getdxy(x,y,nx,nif,opt)                                                  */
/*                                                                          */  
/*  get data into y[] and x[] based on variables in PMVIdx in standard      */
/*  order. If nif == 0 add intercept. Number of columns in x[] are nx.      */
/*  Recognize sel option.                                                   */
/*  Return number of cases.                                                 */
/*  If opt print cases and error message.                                   */
/*                                                                          */
/*  Note: y[1,...], x(i,j) = x[(i - 1) * nx + j]                            */

int getdxy(double *x,double *y,int nx,int nif,int opt)
{
    register int i,j,k,jj;

    k = 0;
    for (i = 0; i < NOC; ++i) {
        if (eval_sve(i) == 0)   
            continue;
       
        y[k + 1] = get_data(PMVIdx[0],i);       

          
        if (nif == 0) {
            x[k * nx + 1] = 1.0;
            jj = 2;
        }
        else
            jj = 1;

        for (j = 1; j < PMNV; ++j) {
            x[k * nx + jj] = get_data(PMVIdx[j],i);       
            jj++;
        }
        k++;  
    }
    if (opt) {
        if (SVEFlg)  
            printf1("Number of selected cases: %d\n",k);
        if (k == 0)  
            p_err(-28,1);  
    }
    return(k);
}

/* ------------------------------------------------------------------------ */
/*  check_nvar(opt) check number of variables, PMNV, and waves, PMNW.       */
/*                  set NPX = number of x variables.                        */
/*                  if opt != 1 also check z variables and set NPZ.         */
/*                  return number of variables, or 0 if error.              */

int check_nvar(int opt)
{
    int n;

    NPX = NPZ = 0;
    n = PMNV / PMNW;
    if (PMNW * n != PMNV) {
        printf1("Error: number of variables/waves is inconsistent.\n");
        return(0);
    }
    NPX = n - 1;

    if (opt && PMNZ > 0 && PMNQ > 0) {
        NPZ = PMNZ / (PMNQ * PMNW);
        if (NPZ * PMNQ * PMNW != PMNZ) {
            printf1("Error: number of z-variables is inconsistent with nq=%d.\n",PMNQ);
            NPZ = 0;
            return(0);
        }
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  prn_nwvar(opt)  print variables for regression etc model.               */
/*  ##              if opt = 1  print z variables.                          */
/*                  if opt = 2  called by mvreg (YL,YH)                     */

void prn_nwvar(int opt)
{
    register int i,j,k,l;
    int n,iv,len,iz,nw;

    nw = PMNW;
    if (opt == 2)
        nw = 1;
    len = 0;
    if (opt && NPZ > 0)  
        len = 4;

    l = PMNVLEN + 1;
    n = PMNV / nw;
    printf1("Variables");
    k = 9;
    if (nw > 1) {
        printf1(" (%d waves)\n",nw);
        k += 11;
    }
    else if (opt <= 1) {
        printf1(" (cross-section)\n");
        k += 16;
    }
    if (opt <= 1)
        prnchar('-',k,1);
    else
        newline();

    if (opt == 2)
        printf1("YL  ");
    else           
        printf1("Y   ");

    prnchar(' ',len,0);
    printf1(": ");
    k = 0;
    for (i = 0; i < nw; ++i) {
        iv = PMVIdx[k++];
        printf1("%s ",VName[iv]);
        prnchar(' ',l - strlen(VName[iv]),0);
    }
    newline();

    for (j = 1; j < n; ++j) {
        if (opt == 2) {
            if (j == 1)
                printf1("YH  ");
            else
                printf1("X%-3d",j - 1);
        }
        else
            printf1("X%-3d",j);

        prnchar(' ',len,0);
        printf1(": ");
        for (i = 0; i < nw; ++i) {
            iv = PMVIdx[k++];
            printf1("%s ",VName[iv]);
            prnchar(' ',l - strlen(VName[iv]),0);
        }
        newline();
    }
    if (opt == 1 && NPZ > 0) {
                 
        for (iz = 0; iz < NPZ; ++iz) {
            k = iz * nw * PMNQ;
            for (j = 0; j < PMNQ; ++j) {
                if (j == 0)
                    printf1("Z%-3d %2d : ",iz + 1,j + 1);
                else
                    printf1("   %4d : ",j + 1);

                for (i = 0; i < nw; ++i) {
                    iv = PMZIdx[k + i * PMNQ + j];
                    printf1("%s ",VName[iv]);
                    prnchar(' ',l - strlen(VName[iv]),0);
                }
                newline();
            }
        }
    }
    newline();
}

/* ------------------------------------------------------------------------ */
/*  prn_data(m,nc,n,x,fd,fmt)                                               */
/*                                                                          */
/*  print n columns from m,nc matrix x[] to fd with print format fmt.       */
/*  note: x(i,j) = x[(i - 1) * nc + j].                                     */

void prn_data(int m,int nc,int n,double *x,FILE *fd,char *fmt)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= n; ++j)
            fprintf(fd,fmt,x[i * nc + j]);
        fprintf(fd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn1_data(y,x,m,ny,nx,fd,fmt)                                           */
/*                                                                          */
/*  printf first ny columns from y[], then nx columns from x[] to fd        */
/*  with format fmt. m = number of cases.                                   */
/*  note: y(i,j) = y[(i - 1) * ny + j], and same for x[].                   */

void prn1_data(double *y,double *x,int m,int ny,int nx,FILE *fd,char *fmt)
{
    register int i,j;

    for (i = 0; i < m; ++i) {
        for (j = 1; j <= ny; ++j)
            fprintf(fd,fmt,y[i * ny + j]);
        for (j = 1; j <= nx; ++j)
            fprintf(fd,fmt,x[i * nx + j]);
        fprintf(fd,"\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  prn1_coeff(n,x,se,ni,df,vidx,seflg)                                     */
/*                                                                          */
/*  print estimated regression coefficients. x[1,...,n] = coefficients,     */
/*  se[] contains standard errors. if ni = 1 without intercept. df =        */  
/*  degrees of freedom for t distribution. If seflg = 1, se[] contains      */
/*  squared standard errors, otherwise stand. errors.                       */
/*                                                                          */
/*  it is assumed that variables are in standard format in vidx[].          */  
/*                                                                          */
/*  if PMPPFDef print estimated parameters to PMPPFd.                       */

void prn1_coeff(int n,double *x,double *se,int ni,int df,short *vidx,int seflg)
{
    register int i,j,k,l,ii;
    double tmp,e;

    l = PMNVLEN + 1;
    if (l < 9)
        l = 9;

    printf1("\nIdx  Wave  Variable ");
    prnchar(' ',l - 8 + PMTFmt1 - 5,0); printf1("Coeff ");          
    prnchar(' ',PMTFmt1 - 5,0); printf1("Error ");          
    prnchar(' ',PMTFmt1 - 7,0); printf1("Coeff/E  Signif\n");           
    prnchar('-',19 + l + 3 * (PMTFmt1 + 1),1);

    k = 1;
    i = PMNW;

    for (j = 1; j <= n; ++j) {
        if (j == 1 && ni == 0) {
            printf1("%3d     -  Intercept ",k);            
            prnchar(' ',l - 9,0);
        }
        else {
            printf1("%3d     1  %s ",k,VName[vidx[i]]);           
            prnchar(' ',l - strlen(VName[vidx[i]]),0);
        }
        printf1(PMTFmtS,x[j]);            
                              
        if (df > 0 && se[j] > EPSI1) {
            if (seflg)
                e = sqrt(se[j]);
            else
                e = se[j];

            printf1(PMTFmtS,e);              
            tmp = x[j] / e;
            printf1(PMTFmtS,tmp);          
            tmp = 2.0 * cdtf(fabs(tmp),df) - 1.0;
            printf1("%7.4lf",tmp);          
        }
        else {
            prnchar(' ',PMTFmt1 - 3,0); printf1("--- ");          
            prnchar(' ',PMTFmt1 - 3,0); printf1("---     ---");          
        }
        newline();
        if (ni == 1 || j > 1) {
            for (ii = 2; ii <= PMNW; ++ii) {
                printf1("      %3d  %s ",ii,VName[vidx[++i]]);          
                prnchar(' ',l - strlen(VName[vidx[++i]]),0);
                newline();
            }
            i++;
        }
        k++;
    }
    if (PMPPFDef) {         /* print estimated parameters */

        for (j = 1; j <= n; ++j) {
            fprintf(PMPPFd,PMTFmtS,x[j]);            
            if (df > 0 && se[j] > 0.0) {
                if (seflg)
                    e = sqrt(se[j]);
                else
                    e = se[j];

                fprintf(PMPPFd,PMTFmtS,e);              
            }
            fprintf(PMPPFd,"\n");
        }
        printf1("\nParameter estimates written to: %s\n",PMPPFName);
    }
}

/* ------------------------------------------------------------------------ */
/*  makefmt(fmt1,fmt2,fmts,mlen,sepc,opt)                                   */
/*      Make format string with fieldlength fmt1, and precision fmt2.       */
/*      Put string in fmts. mlen is the minimum field length, sepc is       */
/*      the separation character at the end of the format string.           */
/*      If opt == 1 special format for plot labels.                         */
/*      if opt == 2 string variables                                        */
/*      if opt == -1 then format with leading zeros                         */

void makefmt(int *fmt1,int *fmt2,char *fmts,int mlen,char sepc,int opt)
{
    if (opt == 2) {
        sprintf(fmts,"%%%ds%c",imax(iabs(*fmt1),mlen),sepc);
        return;
    }
    if (mlen && *fmt1 > -mlen) {          
        if (*fmt1 < 0)          
            *fmt1 = -mlen;        
        else if (*fmt1 < mlen)
            *fmt1 = mlen;
    }
    if (*fmt2 < 0)
        *fmt2 = 0;
    else if (*fmt2 >= 100)
        *fmt2 = 99;

    if (!*fmt1) {
        if (opt == 0)
            sprintf(fmts,"%%lg%c",sepc);
        else if (opt < 0)
            sprintf(fmts,"%%0lg%c",sepc);
        else
            sprintf(fmts,"(%%g)%c",sepc);
    }
    else if (*fmt1 > 0) {
        if (*fmt1 >= 100)
            *fmt1 = 99;
        if (*fmt2 > *fmt1)
            *fmt2 = *fmt1;

        if (opt == 0)
            sprintf(fmts,"%%%d.%df%c",*fmt1,*fmt2,sepc);
        else if (opt < 0)
            sprintf(fmts,"%%0%d.%df%c",*fmt1,*fmt2,sepc);
        else
            sprintf(fmts,"(%%%d.%df)%c",*fmt1,*fmt2,sepc);
    }
    else {
        *fmt1 = -*fmt1;
        if (*fmt1 >= 100)
            *fmt1 = 99;
        if (*fmt1 - *fmt2 < 8)
            *fmt1 = *fmt2 + 8;
        if (opt == 0)
            sprintf(fmts,"%%%d.%de%c",*fmt1,*fmt2,sepc);
        else if (opt < 0)
            sprintf(fmts,"%%0%d.%de%c",*fmt1,*fmt2,sepc);
        else
            sprintf(fmts,"(%%%d.%de)%c",*fmt1,*fmt2,sepc);
    }
}

/* ------------------------------------------------------------------------ */
/*  makenfmt(fmt,fmts,mlen,sepc)                                            */
/*      Make integer format string with fieldlength fmt.                    */
/*      Put string in fmts. mlen is the minimum field length, sepc is       */
/*      the separation character at the end of the format string.           */

void makenfmt(int *fmt,char *fmts,int mlen,char sepc)
{
    if (*fmt < mlen)
        *fmt = mlen;
    else if (*fmt >= 100)
        *fmt = 99;

    sprintf(fmts,"%%%dd%c",*fmt,sepc);
}
      
/* ------------------------------------------------------------------------ */
/*  pmfmt(n,m)      make a new PMFmtS.                                      */

void pmfmt(int n,int m)
{
    PMFmt1 = n;
    PMFmt2 = m;  
    makefmt(&PMFmt1,&PMFmt2,PMFmtS,0,SEPC,0);
}

/* ------------------------------------------------------------------------ */
/*  pmtfmt(n,m)      make a new PMTFmtS.                                    */

void pmtfmt(int n,int m)
{
    PMTFmt1 = n;
    PMTFmt2 = m;  
    makefmt(&PMTFmt1,&PMTFmt2,PMTFmtS,0,SEPC,0);
}

/* ------------------------------------------------------------------------ */
/*  prn_cwt()   Print case weight information.                              */

void prn_cwt(void)         
{
    if (WIVar >= 0) 
        printf1("Using case weights defined by: %s\n",VName[WIVar]);
}

/* ------------------------------------------------------------------------ */
/*  prn_sfmt(s,n,fmt,x)     print string s, then value x with fmt.          */
/*                          if n > 0 adjust length.                         */

void prn_sfmt(char *s,int n,char *fmt,double x)
{
    int l;

    printf1("%s ",s);
    if (n > 0) {
        l = strlen(s);
        if (n > l)
            prnchar(' ',n - l,0);
    }
    printf1(fmt,x);
    printf1("\n");
}

/* ------------------------------------------------------------------------ */
/*  prn_f1mat(fd,prn,fmt,n,x)   write lower triangle of x matrix to fd.     */
/*                                                                          */
/*                  prn= 0 : lower triangle, including main diagonal        */
/*                  prn= 1 : full square matrix                             */
/*                  prn= 2 : column vector: lower triangle                  */
/*                  prn= 3 : column vector: lower triangle incl. diag.      */
/*                  prn= 4 : column vector: full matrix                     */
/*                                                                          */

void prn_f1mat(FILE *fd,int prn,char *fmt,int n,double *x)
{
    register int i,j,k;

    if (prn >= 1 && prn <= 4) {

        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (prn == 2 && j >= i)
                    break;
                else if (prn == 3 && j > i)
                     break;

                if (j <= i)
                    k = i * (i - 1) / 2 + j - 1;
                else
                    k = j * (j - 1) / 2 + i - 1;

                fprintf(fd,fmt,x[k]);
                if (prn != 1)
                    fprintf(fd,"\n");
            }
            if (prn == 1)
                fprintf(fd,"\n");
        }
    }
    else {
        k = 0;
        for (i = 0; i < n; ++i) {
            for (j = 0; j <= i; ++j)  
                fprintf(fd,fmt,x[k++]); 
            fprintf(fd,"\n");
        }
    }
}



