/****************************************************************************/
/*  t_con                                                                   */
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
#include "t_alloc.h"
#include "t_lsei.h"

/*  functions in t_con.c */

int p_con(char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni);
int get_con(char *pcmd);
int con_proc(char *pcmd);
int con_build(void);
void con_free(void);

/* ------------------------------------------------------------------------ */
int NCon = 0;               /* Number of constraints                        */
int NCon1 = 0;              /* Number of constraints used for minimization  */
double *ConD;               /* Right-hand side of constraints               */
int ConDA = 0;              /* if allocated                                 */
double *ConR;               /* Left-hand side of constraints                */
int ConRA = 0;              /* if allocated                                 */
double *ConQ;               /* Projection matrix (NParm, NParm1)            */
int ConQA = 0;              /* if allocated                                 */
int *CIP;                   /* Pointers to columns of ConQ                  */
int CIPA = 0;               /* if allocated                                 */

/* ------------------------------------------------------------------------ */
/*  p_con(pcmd,nx,nif,mw,nw,w,ne,ni)                                        */
/*                                                                          */
/*  get constraints: lsecon, lsicon from string pcmd.                       */
/*  nx = number of b parameters. if nif = 1 then without intercept.         */
/*  put lsecon constraints at begin, lsicon constraints at the end of       */
/*  data matrix w[] which has mw rows and nw columns.                       */
/*  return *ne = number of lsecon, *ni = number of lsicon constraints.      */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

int p_con(char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni)
{
    int nb,eflag,fflag,err,ir,ctyp;
    register char *p;
    double tmp;

    p = pcmd;         

    *ne = 0;
    *ni = 0;
    err = 0;
    while (*p) {

        ctyp = 0;
        if (!strncmp(p,"lsecon",6)) {     
            printf1("LSECon: ");
            ctyp = 1;
            *ne += 1;
            ir = *ne;
            p += 6;
        }
        else if (!strncmp(p,"lsicon",6)) {     
            printf1("LSICon: ");
            ctyp = 2;
            *ni += 1;
            ir = mw + 1 - *ni;
            p += 6;
        }
        if (ctyp) {
            if (*p++ != '=')
                goto PCONFin;

            err = -1;
            fflag = eflag = 0;
            while (*p) {
                if (sscanf(p,"%lg",&tmp) != 1) {

                    if (*p == 'b')
                        tmp = 1.0;

                    else if (*p == '+' && *(p + 1) == 'b') {
                        tmp = 1.0;
                        p++;
                    }
                    else if (*p == '-' && *(p + 1) == 'b') {
                        tmp = -1.0;
                        p++;
                    }
                    else
                        goto PCONFin;
                }
                else
                    p = skip_dbl(p);

                printf1("%g ",tmp);
                if (eflag) {
                    w[(ir - 1) * nw + nw] = tmp;
                    break;
                }
                if (*p == '*')
                    p++;

                if (sscanf(p,"b%d",&nb) != 1) 
                    goto PCONFin;
                p = skip_int(p + 1);

                printf1("b%d ",nb);
                if (nb < nif || nb > nx)
                    fflag = 1; 
                else {
                    if (nif == 0)
                        nb++;
                    w[(ir - 1) * nw + nb] = tmp;
                }
                if (*p != '+' && *p != '-' && *p != '=')
                    goto PCONFin;

                if (*p != '-')  
                    printf1("%c ",*p);
                if (*p == '=') {
                    p++;
                    eflag = 1;
                }
            }
            if (eflag == 0 || (*p != ',' && *p != ')')) 
                goto PCONFin;
 
            if (fflag) {
                err = -2;
                goto PCONFin;
            }
            newline();
        }
        p++;
    }
    err = 0;  

PCONFin:
    if (err == -1)   
        printf1("\nSyntax error.\n");
    else if (err == -2)  
        printf1("\nError in parameter index.\n");
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  get_con(pcmd)                                                           */
/*  ##                                                                      */
/*  Get constraints: con=... from string pcmd.                              */
/*  There should be PMNConS con=... strings in pcmd (set in parm()).        */
/*  It is assumed that there are NParm parameters, corresponding to         */
/*  b1,b2,... in the definition of constraints.                             */
/*                                                                          */
/*  Save constraints in ConR and ConD.                                      */
/*                                                                          */
/*  return 0 if OK, -1 if insufficient memory, -2 if syntax error.          */

int get_con(char *pcmd)
{
    int nb,eflag,fflag,err,first;
    register char *p,*pp;
    double tmp;

    err = -1;
    if (!(ConR = (double *)calloc(NParm * PMNConS + 1,sizeof(double)))) {
        p_err(-2,1);
        goto GCONFin;
    }
    ConRA = NParm * PMNConS + 1;
    memrq(ConRA,sizeof(double));

    if (!(ConD = (double *)calloc(PMNConS + 1,sizeof(double)))) {
        p_err(-2,1);
        goto GCONFin;
    }
    ConDA = PMNConS + 1;
    memrq(ConDA,sizeof(double));

    p = pcmd;         

    NCon = 0;           /* number of constraints */
    err = -2;
    while (*p) {

        if (!strncmp(p,"con=",4)) {     
            pp = p;
            first = 1;
            p += 4;

            if (++NCon > PMNConS)
                goto GCONFin;

            fflag = eflag = 0;
            while (*p) {
                if (sscanf(p,"%lg",&tmp) != 1) {
                    if (*p == 'b')
                        tmp = 1.0;
                    else if (*p == '+') {
                        tmp = 1.0;
                        p++;
                    }
                    else if (*p == '-') {
                        tmp = -1.0;
                        p++;
                    }
                    else
                        goto GCONFin;
                }
                else
                    p = skip_dbl(p);

                if (first) {
                    printf1("Con: ");
                    first = 0;
                }
                printf1("%g",fabs(tmp));
                if (eflag) {
                    ConD[NCon] = tmp;
                    break;
                }
                if (*p == '*')
                    p++;

                if (sscanf(p,"b%d",&nb) != 1)                          
                    goto GCONFin;

                if (nb < 1 || nb > NParm) {
                    fflag = 1;
                    goto GCONFin;
                }
                p = skip_int(p + 1);

                printf1(" * b%d ",nb);
                ConR[(NCon - 1) * NParm + nb] = tmp;

                if (*p != '+' && *p != '-' && *p != '=')
                    goto GCONFin;

                printf1("%c ",*p);
                if (*p == '=') {
                    p++;
                    eflag = 1;
                }
            }
            if (eflag == 0 || (*p != ',' && *p != ')')) 
                goto GCONFin;
            newline();
        }
        p++;
    }
    err = 0;  

GCONFin:
    if (err)  
        con_free();
    if (err == -2) {
        p = skip_com(pp);
        *p = '\0';
        printf1("\nError: %s\n",pp);
        if (fflag)  
            printf1("Error in parameter index.\n");
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_proc()      process constraints for ML estimation.                  */
/*                  write into protocol file.                               */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

int con_proc(char *pcmd)
{
    register int i,j;
    int err,r,sflag;
    double rnorm,tmp;

    err = -1;
    printf1("Checking constraints.\n");
    if (PMNConS >= NParm) {
        printf1("Error: number of constraints should be less than number of parameters.\n");
        goto CONPROCFin;
    }
    if (get_con(pcmd))
        goto CONPROCFin;

    printf1("\nNumber of constraints: %d\n",NCon);         
    if (NCon == 0)
        goto CONPROCFin;

    if (PMProtFDef) {
        fprintf(PMProtFd,"Constraints (R,D).\n");
        for (i = 1; i <= NCon; ++i) {
            for (j = 1; j <= NParm; ++j)
                fprintf(PMProtFd,PMPFmtS,ConR[(i - 1) * NParm + j]);
            fprintf(PMProtFd,PMPFmtS,ConD[i]);
            fprintf(PMProtFd,"\n");
        }
        fprintf(PMProtFd,"\n");
    }

    /* we need additional parameter vectors. ParS for (new) starting values,
       and ParC for expansion. */

    if (!(ParS = (double *)calloc(NParm + 1,sizeof(double)))) {
        p_err(-2,1);
        goto CONPROCFin;
    }
    ParSA = NParm + 1;
    memrq(ParSA,sizeof(double));

    if (!(ParC = (double *)calloc(NParm + 1,sizeof(double)))) {
        p_err(-2,1);
        goto CONPROCFin;
    }
    ParCA = NParm + 1;
    memrq(ParCA,sizeof(double));

    /* check whether starting values fullfil constraints */

    sflag = 0;
    for (i = 1; i <= NCon; ++i) {
        tmp = 0.0;
        for (j = 1; j <= NParm; ++j)
            tmp += ConR[(i - 1) * NParm + j] * Par[j];
        if (fabs(tmp - ConD[i]) > EPSI1) {
            sflag = 1;
            break;
        }
    }
    if (sflag) {
        printf1("Try to find consistent starting values.\n");         

        if (alloc_acu(ConRA))
            goto CONPROCFin;

        if (alloc_acv(NParm + 1))
            goto CONPROCFin;

        if (alloc_acn(NParm + 1))
            goto CONPROCFin;

        for (i = 1; i < ConRA; ++i)
            AcU[i] = ConR[i];

        for (i = 1; i <= NCon; ++i)
            AcV[i] = ConD[i];

        r = lhhfti(NCon,NParm,NParm,1,AcU,AcV,&rnorm,1,WrkD,WrkH,AcN,EPSI1);

        printf1("Rank of constraints: %d\n",r);          
        if (r != NCon) {
            printf1("Error: should be equal to number of constraints.\n");
            goto CONPROCFin;
        }
        if (rnorm > EPSI1) {
            printf1("Can't find a solution of constraints.\n");
            goto CONPROCFin;
        }
        for (i = 1; i <= NParm; ++i)    /* save new starting values */
            ParS[i] = AcV[i];

        if (PMProtFDef) {
            fprintf(PMProtFd,"Solution of constraints (new starting values).\n");
            for (i = 1; i <= NParm; ++i)  
                fprintf(PMProtFd,PMPFmtS,AcV[i]);
            fprintf(PMProtFd,"\n\n");
        }
        alloc_acu(0);
        alloc_acv(0);
        alloc_acn(0);
    }
    else {
        for (i = 1; i <= NParm; ++i)    /* use original starting values */
            ParS[i] = Par[i];
    }
    NCon1 = NCon;               /* used for estimation */
    NParm1 = NParm - NCon;      /* dimensions of reduced parameter space */
    HSiz1 = NParm1 * (NParm1 - 1) / 2;  

    if (con_build())     
        goto CONPROCFin;

    /* set Par[], used for minimization, to zero */

    for (i = 1; i <= NParm1; ++i) 
        Par[i] = 0.0;

    err = 0;

CONPROCFin:
    if (err)
        con_free();

    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_build.   Build the projection matrix ConQ for function minimization */
/*  ##           in the reduced parameter space. ConQ is a NParm * NParm1   */
/*               matrix. Cf. McCormick 1983, p.265.                         */
/*                                                                          */
/*  Return 0 if OK, -1 if insufficient memory.                              */
/*                  -2 if less than NCon indep constraints.                 */

int con_build(void)
{
    register int i,j,k,l;
    int r,n,n1,err;

    err = -1;
    if (!(ConQ = (double *)calloc(NParm * NParm1 + 1,sizeof(double)))) {
        p_err(-2,1);
        goto BCONFin;
    }
    ConQA = NParm * NParm1 + 1;
    memrq(ConQA,sizeof(double));

    if (!(CIP = (int *)calloc(NParm + 1,sizeof(int)))) {
        p_err(-2,1);
        goto BCONFin;
    }
    CIPA = NParm + 1;
    memrq(CIPA,sizeof(int));

    if (alloc_acu(ConRA))           /* used as copy of ConR */
        goto BCONFin;

    if (alloc_acw(NParm + 1))
        goto BCONFin;

    if (alloc_acn(NParm + 1))
        goto BCONFin;

    /*  Try to find NCon linear independent columns of ConR. */

    n = 0;  /* number of columns already found, indices in CIP    */
    k = 0;  /* CIP is used to record the column interchange. */

    for (j = 1; j <= NParm; ++j)
        CIP[j] = j;

    for (j = 1; j <= NParm; ++j) {
        n1 = n + 1;
        for (i = 1; i <= n; ++i) {
            for (l = 1; l <= NCon; ++l)
                AcU[(l - 1) * n1 + i] = ConR[(l - 1) * NParm + CIP[i]];
        }
        for (l = 1; l <= NCon; ++l)
            AcU[(l - 1) * n1 + n1] = ConR[(l - 1) * NParm + j];

        r = lhhfti(NCon,n1,n1,0,AcU,AcU,AcW,0,WrkD,WrkH,AcN,EPSI1);
        if (r == n1) {
            k = CIP[++n];
            CIP[n] = j;
            CIP[j] = k;
        }
        if (n == NCon)
            break;
    }

    if (PMProtFDef) {
        fprintf(PMProtFd,"Selected columns for projection matrix.\n ");
        for (i = 1; i <= NCon; ++i)  
            fprintf(PMProtFd,"%2d ",CIP[i]);
        fprintf(PMProtFd,"\n\n");
    }
    if (n != NCon) {
        err = -2;
        goto BCONFin;
    }
   
    /*  Now build the first NCon rows of the projection matrix ConQ.    */
    
    for (j = 1; j <= NCon; ++j) {
        for (i = 1; i <= NCon; ++i)
            AcU[(j - 1) * NCon + i] = ConR[(j - 1) * NParm + CIP[i]];

        for (i = 1; i <= NParm1; ++i)  
            ConQ[(j - 1) * NParm1 + i] =
                             -ConR[(j - 1) * NParm + CIP[NCon + i]];
    }
    r = lhhfti(NCon,NCon,NCon,NParm1,AcU,ConQ,AcW,1,WrkD,WrkH,AcN,EPSI1);

    if (r != NCon) {
        n = r;
        err = -2;
        goto BCONFin;
    }
    
    /*  The lower part of ConQ is an identity matrix. */
   
    j = 1;
    for (i = NCon + 1; i <= NParm; ++i) {
        ConQ[(i - 1) * NParm1 + j] = 1.0;
        j++;
    }
    for (i = 1; i <= NParm; ++i) {
        for (j = 1; j <= NParm1; ++j) {
            if (fabs(ConQ[(i - 1) * NParm1 + j]) <= EPSI1)
                ConQ[(i - 1) * NParm1 + j] = 0.0;
        }
    }
    
    /*  Build pointers to the interchanged rows of ConQ in CIP. */
   
    for (i = 1; i <= NParm; ++i)
        AcN[CIP[i]] = i;

    for (i = 1; i <= NParm; ++i)
        CIP[i] = AcN[i];

    if (PMProtFDef) {
        fprintf(PMProtFd,"Parameter interchange for projection.\n ");
        for (i = 1; i <= NParm; ++i)
            fprintf(PMProtFd,"%2d ",CIP[i]);

        fprintf(PMProtFd,"\n\nProjection matrix.\n ");
        for (i = 1; i <= NParm; ++i) {
            for (j = 1; j <= NParm1; ++j)  
                fprintf(PMProtFd," %g",ConQ[(i - 1) * NParm1 + j]);
            fprintf(PMProtFd,"\n ");
        }
        fprintf(PMProtFd,"\n");
    }
    err = 0;

BCONFin:
    if (err == -2)  
        printf1("Error: only %d linear independent constraint(s).\n",n);
    alloc_acu(0);
    alloc_acw(0);
    alloc_acn(0);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  con_free()      free memory used for constraints.                       */
/*                                                                          */
/*  return 0 if OK, -1 if error.                                            */

void con_free(void)
{
    if (ParSA) {
        free((char *)ParS);
        memrq(-ParSA,sizeof(double));
        ParSA = 0;
    }
    if (ParCA) {
        free((char *)ParC);
        memrq(-ParCA,sizeof(double));
        ParCA = 0;
    }
    if (ConRA) {
        free((char *)ConR);
        memrq(-ConRA,sizeof(double));
        ConRA = 0;
    }
    if (ConDA) {
        free((char *)ConD);
        memrq(-ConDA,sizeof(double));
        ConDA = 0;
    }
    if (ConQA > 0) {
        free((char *)ConQ);
        memrq(-ConQA,sizeof(double));
        ConQA = 0;
    }
    if (CIPA > 0) {
        free((char *)CIP);
        memrq(-CIPA,sizeof(int));
        CIPA = 0;
    }
    NParm1 = NParm;
    HSiz1 = HSiz;
    NCon1 = NCon = 0;
}
