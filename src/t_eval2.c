/****************************************************************************/
/*  t_eval2                                                                 */
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
#include "t_var.h"
#include "t_gdat.h"
#include "t_rand.h"
#include "t_gf.h"
#include "t_cdf.h"
#include "t_rzoo.h"
#include "t_edat.h"
#include "t_parm.h"
#include "t_smo.h"
#include "t_alloc.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_mat.h"
#include "t_matf.h"
#include "t_eval.h"
#include "t_eval1.h"
#include "t_eval3.h"
#include "t_mdat.h"
#include "t_sd.h"
#include "t_gm.h"
#include "t_map.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_eval2.c                                                  */

int v_eval2(int vn,short cnt,int *typ,double *val,int *idx,
    double *res,int r1,int r2,int ncol,int bnum);
int edcomp(const void *arg1,const void *arg2);
int es_comp(const void *arg1,const void *arg2);
int ed1comp(const void *arg1,const void *arg2);

/* ------------------------------------------------------------------------ */
/*  check for identical definition in t_eval.c, t_eval2.c                   */

#define DDTYP     400   /* type number for dummy operators                  */
#define CCTYP      99   /* type number for ,, in dummy operators            */
#define COLTYP     98   /* type number for :  in dummy operators            */
#define RVECTYP  1698   /* type number for row vectors <...>'               */
#define CVECTYP  1699   /* type number for column vectors <...>             */

int *BPtr;              /* pointer for sort operators */
float *BVal;            /* values for sort */
short *ZVal;            /* censoring indicator */
char *BPtrSBuf;         /* buffer with values of string variable */
int BPtrSBLen = 0;      /* length of string variable */

/* ------------------------------------------------------------------------ */
/*  v_eval2(int vn,short cnt,int *typ,double *val,int *idx,                 */
/*          double *res,int r1,int r2,int ncol,int bnum)                    */
/*                                                                          */
/*  Evaluate the expression defined in typ and val. Number of               */
/*  stack entries is cnt. The expression may contain type 2 operators.      */
/*                                                                          */
/*  If vn >= 0 this is interpreted as an internal variable number and the   */
/*  result is put into storage with put_data() for all cases:               */
/*  row = r1,...,r2 - 1.                                                    */
/*  Otherwise, the result is returned in the vector res[i],                 */
/*  i = 0,...,r2 - r1.                                                      */
/*                                                                          */
/*  ncol is the current column number, only used for matrix expressions.    */
/*  ncol = 0,1,2,...                                                        */
/*                                                                          */
/*  bnum is the current block number if in block mode.                      */
/*                                                                          */
/*  return 0 if ok, otherwise one of following error indicators:            */
/*                                                                          */
/*  -1  =   undefined X or C or V variable                                  */
/*   1  =   undefined operator                                              */
/*   2  =   division by zero                                                */
/*   3  =   % operator finds non-integers                                   */
/*   4  =   exp out of range                                                */
/*   5  =   argument of log is too small (< LogMin)                         */
/*   6  =   negative argument in sqrt                                       */
/*   7  =   negative or zero argument for lgam                              */
/*   8  =   arguments for icg out of range                                  */
/*   9  =   arguments for icb out of range                                  */
/*  10  =   incorrect range in tr(x,a,b)  (a > b)                           */
/*  11  =   argument for ndi() out ouf range [0,1]                          */
/*  12  =   incorrect arguments for TD                                      */
/*  13  =   incorrect arguments for CHD                                     */
/*  14  =   incorrect arguments for FD                                      */

/*  17  =   undefined intermediate function parameter                       */
/*  18  =   undefined function argument                                     */

/*  32  =   error: cannot calculate incomplete gamma integral               */
/*  33  =   error: cannot calculate binomial coefficient                    */
/*  34  =   error: cannot calculate biv normal distribution                 */
/*  35  =   error: cannot init rdmn generator                               */
/*  36  =   error: wrong index in rdmn operator                             */

/*  42  =   error: in poisson operator                                      */
/*  43  =   error: in negbin operator                                       */
/*  45  =   error: matrix indices must be positive                          */

/*  49  =   error: spatial data required                                    */
/*  50  =   error: cannot read spatial data                                 */
/* ------------------------------------------------------------------------ */
/*  95  =   insufficient memory                                             */
/*  96  =   no cases in data matrix                                         */
/*  97  =   out of stack space                                              */
/*  98  =   if expression contains references to c1,c2,...                  */
/*  99  =   operator that cannot be used in connection with type 2 ops.     */

#define STACKLEVEL 9     /* size of stack */
    
int v_eval2(int vn,short cnt,int *typ,double *val,int *idx,
            double *res,int r1,int r2,int ncol,int bnum)
{
    register int i,j,j1,k,l,l1;
    int jj,j2,j3,j4,t,err,atyp,utyp,narg,noc,n0,n1,n2,nn,row,col,esidx;
    int kk,mrow,mcol,aval,ptr,ptr1;
    double *buf[STACKLEVEL];             
    double tmp,tmp1,tmp2,tmp3,tmp4,e,x,y,z,dval[7];
    float grp,ogrp;
    register char *p,*q;

    err = 0;
    mrow = mcol = 1;
    noc = r2 - r1;
    if (noc < 1)
        return(96);
              
    /* check for type 2 operators, exclude pre (1510) and suc (1509) */

    t = 0;
    for (i = 0; i < cnt; ++i) {
        atyp = iabs(typ[i]);
        if (atyp >= OPT2A && atyp < OPT2B && atyp != 1510 && atyp != 1509) {
            t = 1;
            break;
        }
    }

    /* if the expression does not contain type 2 operators, use v_eval1() */
  
    if (t == 0) {        /* && (vn < 0 || VTyp[vn] < 4)) { */

        E2FNum = r1;
        E2LNum = r2;

        for (l = r1; l < r2; ++l) {

            err = v_eval1(l,cnt,typ,val,idx,&x,0,0,0,0,0);
            if (err)  
                break;           
            if (vn >= 0)
                put_data(x,vn,l);
            else
                res[l - r1] = x;
        }
        E2FNum = 0;
        E2LNum = 0;
        return(err);
    }
  
    /* special procedure for expressions containing type 2 operators */
    /* first allocate memory for noc cases */

    for (i = 0; i < STACKLEVEL; ++i) {
        if (!(buf[i] = (double *)calloc(noc+1,sizeof(double)))) { 
            for (j = 0; j < i; ++j) {
                free((char *)buf[j]);
                memrq(-noc-1,sizeof(double));
            }
            err = 1;
            break;
        }
        else
            memrq(noc+1,sizeof(double));
    }
    if (err) {
        printf1("Error: insufficient memory in v_eval2().\n");
        return(95);
    }
    jj = j = 0;
    for (i = 0; i < cnt; ++i) {

        utyp = typ[i];   
        atyp = iabs(utyp);
        narg = (int)val[i];
        esidx  = idx[i];
        aval = (int)val[i];

        if (esidx >= NLOFFS && esidx < COFFS) {
            if (esidx >= MOFFS && esidx < VOFFS) {
                k = esidx - MOFFS;
                mrow = imax(mrow,MatRow[k]);
                mcol = imax(mcol,MatCol[k]);
            }
            else  
                mrow = imax(mrow,NOC);
        }

        if (atyp >= IFTYP) {        /* ignore if then else */
            ;
        }
        else if (atyp >= COFFS) {
            err = 98;
            goto E2FIN;
        }
        else if (atyp >= VOFFS) {  /* V-variables */
    
            if (jj >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            k = atyp - VOFFS;
    
            if (NOC > 0) {
                for (l = r1; l < r2; ++l) {
                    tmp = get_data(k,l % NOC);

                    if (utyp < 0)  
                        tmp = -tmp;
                    buf[jj][l - r1] =  tmp;                  
                }
            }
            else {
                for (l = r1; l < r2; ++l) {
                    tmp = get_data(k,l);

                    if (utyp < 0)  
                        tmp = -tmp;
                    buf[jj][l - r1] =  tmp;                  
                }
            }
            jj += 1;
        }
        else if (atyp >= MOFFS) {  /* use current column number: ncol */

/* ## */
            if (jj >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            k = atyp - MOFFS;
            row = MatRow[k];      
            col = MatCol[k];
            n0 = (ncol % MatCol[k]) + 1;
            for (l = r1; l < r2; ++l) {
                tmp = MatVal[k][(l % row) * col + n0];
                if (utyp < 0)  
                    tmp = -tmp;
                buf[jj][l - r1] =  tmp; 
            }
            jj += 1;
        }
        else if (atyp == 0) {      /* value */

            if (jj >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            for (l = 0; l < noc; ++l)  
                buf[jj][l] = val[i];
            jj++;      
        }
        else if (atyp < 100) {      /* basic operators */

            if (atyp == CCTYP || atyp == COLTYP) {
                j1 = jj;
                jj++;
            }
            else {
                jj--;     
                j = jj;
                j1 = j - 1;
            }

            for (l = 0; l < noc; ++l) {
       
                switch (atyp) {

                    case CCTYP:                     /* special: CCTYP */
                        buf[j1][l] = (double)INTMAX;
                        break;

                    case COLTYP:                    /* special: COLTYP */
                        buf[j1][l] = (double)(INTMAX - 1);
                        break;

                    case 1:                         /* + */
                      buf[j1][l] += buf[j][l];
                      break;

                    case 2:                         /* - */
                      buf[j1][l] -= buf[j][l];
                      break;

                    case 3:                         /* * */
                      buf[j1][l] *= buf[j][l];
                      break;

                    case 4:                         /* & */
                      x = buf[j1][l] * buf[j][l];
                      if (fabs(x) > EPSI2)
                          buf[j1][l] = 1.0;
                      else
                          buf[j1][l] = 0.0;            
                      break;

                    case 5:                         /* / */
                      if (buf[j][l] == 0.0) {
                          err = 2;
                          goto E2FIN;
                      }
                      else           
                          buf[j1][l] /= buf[j][l];
                      break;

                    case 6:                         /* ^ */
                      x = buf[j1][l];
                      if (fabs(x) < EPSI2)
                          buf[j1][l] = 0.0;
                      else
                          buf[j1][l] = pow(x,buf[j][l]);
                      break;

                    case 7:                         /* % */
                      x = buf[j1][l];
                      y = buf[j][l];
                      if (x != (int)x || y != (int)y) { 
                          err = 3;
                          goto E2FIN;
                      }
                      else               
                          buf[j1][l] = (double)((int)x % (int)y);
                      break;

                    case 8:                         /* | */

                      if (fabs(buf[j1][l]) > EPSI2 || fabs(buf[j][l]) > EPSI2)
                          buf[j1][l] = 1.0;
                      else
                          buf[j1][l] = 0.0;            
                      break;

                    default:
                      err = 99;
                      goto E2FIN;
                }
            }
        }
        else if (atyp == 108 || atyp == 109) {  /* min(...), max(...) */

            for (l = 0; l < noc; ++l) {
                tmp = buf[jj - 1][l];
                for (k = 2; k <= narg; ++k) {
                    if (atyp == 108)
                        tmp = dmin(tmp,buf[jj - k][l]);
                    else
                        tmp = dmax(tmp,buf[jj - k][l]);
                }
                buf[jj - narg][l] = tmp;
            }
            jj -= (narg - 1);
        }
        else if (atyp >= 196 && atyp <= 199) {  /* matrix elements */

            /*    case 196:  X(.,.) */
            /*    case 197:  X(.,j) */
            /*    case 198:  X(i,.) */
            /*    case 199:  X(i,j) */

            n1 = n2 = -1;           /* n1 = row, n2 = col */
            if (atyp == 196)
                j1 = jj;
            else if (atyp == 197) {
                n2 = (int)buf[jj - 1][0] - 1;                    
                j1 = jj - 1;
            }
            else if (atyp == 198) {
                n1 = (int)buf[jj - 1][0] - 1;                    
                j1 = jj - 1;
            }
            else {
                n1 = (int)buf[jj - 2][0] - 1;   
                n2 = (int)buf[jj - 1][0] - 1;   
                j1 = jj - 2;
            }
            kk = iabs(esidx);
            tmp = 0.0;
            if (kk >= COFFS)
                ;            
            else if (kk >= VOFFS) {
                if (n1 >= NOC || n2 >= 1)
                    return(45);
                if (n1 >= 0)
                    tmp = get_data(kk - VOFFS,n1);
                else {
                    for (l = 0; l < NOC; ++l)
                        tmp += get_data(kk - VOFFS,l);
                }
            }
            else if (kk >= MOFFS) {
                n0 = kk - MOFFS;
                row = MatRow[n0];
                col = MatCol[n0];
                if (n1 >= row || n2 >= col)
                    return(45);
                if (n1 >= 0 && n2 >= 0) 
                    tmp = MatVal[n0][n1 * col + n2 + 1];
                else if (n1 >= 0 && n2 < 0) {
                    for (l1 = 1; l1 <= col; ++l1)  
                        tmp += MatVal[n0][n1 * col + l1];
                }
                else if (n1 < 0 && n2 >= 0) {  
                    for (l = 0; l < row; ++l)  
                        tmp += MatVal[n0][l * col + n2 + 1];
                }
                else {  
                    for (l = 0; l < row; ++l) {
                        for (l1 = 1; l1 <= col; ++l1)
                            tmp += MatVal[n0][l * col + l1];
                    }
                }
            }
            else if (kk >= NLOFFS) {
                n0 = kk - NLOFFS;
                col = NLNV[n0];
                if (n1 >= NOC || n2 >= col)
                    return(45);

                if (n1 >= 0 && n2 >= 0) 
                    tmp = get_data(NLVIdx[n0][n2],n1);
                else if (n1 >= 0 && n2 < 0) {
                    for (l1 = 0; l1 < col; ++l1)  
                        tmp += get_data(NLVIdx[n0][l1],n1);
                }
                else if (n1 < 0 && n2 >= 0) {  
                    for (l = 0; l < NOC; ++l)  
                        tmp += get_data(NLVIdx[n0][n2],l);
                }
                else {  
                    for (l = 0; l < NOC; ++l) {
                        for (l1 = 0; l1 < col; ++l1)
                            tmp += get_data(NLVIdx[n0][l1],l);
                    }
                }
            }
            if (esidx < 0)
                tmp = -tmp;

            for (l = 0; l < noc; ++l)
                buf[j1 - 1][l] = tmp;     
            jj = j1;
        }
        else if (atyp < 200) {  /* mathematical operators */
    
            j1 = jj - 1;
            j2 = jj - 2;
            j3 = jj - 3;
            j4 = jj - 4;
            
            for (l = 0; l < noc; ++l) {

                switch (atyp) {

                    case 100:                 /* exp */
                        x = buf[j1][l];
                        if (x < ExpMin || x > ExpMax) {
                            err = 4;
                            goto E2FIN;
                        }
                        buf[j1][l] = exp(x);
                        break;

                    case 101:                 /* log */
                        x = buf[j1][l];
                        if (x <= 0.0) {
                            err = 5;
                            goto E2FIN;
                        }
                        buf[j1][l] = rlog(x);
                        break;

                    case 102:                 /* rnd */
                        x = buf[j1][l];
                        if (x >= 0.0)
                            buf[j1][l] = floor(x + 0.5);
                        else
                            buf[j1][l] = -floor(-x + 0.5);
                        break;

                    case 103:                 /* abs */
                        buf[j1][l] = fabs(buf[j1][l]);
                        break;

                    case 104:                 /* sign */
                        x = buf[j1][l];
                        if (x > 0.0)
                            buf[j1][l] = 1.0;
                        else if (x < 0.0)
                            buf[j1][l] = -1.0;
                        else
                            buf[j1][l] = 0.0;
                        break;

                    case 105:                 /* floor */
                        buf[j1][l] = floor(buf[j1][l]);
                        break;

                    case 106:                 /* ceil */
                        buf[j1][l] = ceil(buf[j1][l]);
                        break;

                    case 107:                 /* tr(x,a,b) */
                        tmp2 = buf[j1][l];
                        tmp1 = buf[j2][l];
                        tmp  = buf[j3][l];

                        if (tmp1 > tmp2) {
                            err = 10;
                            goto E2FIN;
                        }
                        else if (tmp < tmp1)
                            buf[j3][l] = tmp1;
                        else if (tmp > tmp2)
                            buf[j3][l] = tmp2;
                        else
                            buf[j3][l] = tmp;
                        break;

                    case 110:                 /* sqrt */
                        x = buf[j1][l];
                        if (x < 0.0) {
                            err = 6;
                            goto E2FIN;
                        }
                        buf[j1][l] = sqrt(x);
                        break;

                    case 111:                 /* sin */
                        x = buf[j1][l];
                        buf[j1][l] = sin(x);
                        break;

                    case 112:                 /* cos */
                        x = buf[j1][l];
                        buf[j1][l] = cos(x);
                        break;

                    case 113:                 /* lgam */
                        x = buf[j1][l];
                        if (x <= 0.0) {
                            err = 7;
                            goto E2FIN;
                        }
                        buf[j1][l] = loggam(x,&err);
                        if (err) {
                            err = 7;
                            goto E2FIN;
                        }
                        break;

                    case 114:                 /* icg1 */
                        tmp1 = buf[j1][l];
                        tmp  = buf[j2][l];
                        if (tmp <= 0.0 || tmp1 <= 0.0) {
                            err = 8;                
                            goto E2FIN;
                        }
                        buf[j2][l] = icgam(tmp,tmp1,&err);   
                        if (err) {
                            err = 32;
                            goto E2FIN;
                        }
                        break;


                    case 115:                 /* icb */
                        tmp2 = buf[j1][l];
                        tmp1 = buf[j2][l];
                        tmp  = buf[j3][l];
                        if (tmp < 0.0 || tmp > 1.0 || tmp1 <= 0.0 || tmp2 <= 0.0) {
                            err = 9;
                            goto E2FIN;
                        }
                        buf[j3][l] = incbeta(tmp,tmp1,tmp2,&err);
                        if (err) {
                            err = 9;
                            goto E2FIN;
                        }
                        break;

                    case 116:                 /* eexp */
                        x = buf[j1][l];
                        if (x < ExpMin || x > ExpMax) {
                            err = 4;
                            goto E2FIN;
                        }
                        buf[j1][l] = exp(x) / (1.0 + exp(x));
                        break;

                    case 117:                 /* digam */
                        x = buf[j1][l];
                        if (x <= 0.0) {
                            err = 7;
                            goto E2FIN;
                        }
                        buf[j1][l] = digam(x,&err);
                        if (err) {
                            err = 7;
                            goto E2FIN;
                        }
                        break;

                    case 118:                 /* trigam */
                        x = buf[j1][l];
                        if (x <= 0.0) {
                            err = 7;
                            goto E2FIN;
                        }
                        buf[j1][l] = trigam(x,&err);
                        if (err) {
                            err = 7;
                            goto E2FIN;
                        }
                        break;

                    case 120:                 /* icg */
                        tmp1 = buf[j1][l];
                        tmp  = buf[j2][l];
                        if (tmp <= 0.0 || tmp1 <= 0.0) {
                            err = 8;                
                            goto E2FIN;
                        }
                        if (icdgam(tmp,tmp1,dval)) {
                            err = 32;
                            goto E2FIN;
                        }
                        buf[j2][l] = dval[6];                  
                        break;

                    case 121:                 /* bc */
                        tmp1 = buf[j1][l];
                        tmp  = buf[j2][l];
                        tmp = bincoeff((int)tmp,(int)tmp1);
                        if (tmp <= 0.0) {
                            err = 33;                
                            goto E2FIN;
                        }
                        buf[j2][l] = tmp;                      
                        break;

                    case 122:                 /* bivn(x,y,rho) */
                        tmp2 = buf[j1][l];
                        tmp1 = buf[j2][l];
                        tmp  = buf[j3][l];
                        buf[j3][l] = bivn(tmp,tmp1,tmp2,&err);
                        if (err) {
                            err = 34;
                            goto E2FIN;
                        }
                        break;

                    case 124:                 /* poisson(theta,k) */
                        tmp1 = buf[j1][l];
                        tmp  = buf[j2][l];
                        err = l_poisson(buf[j2][l],buf[j1][l],&tmp,&tmp1,&tmp2);
                        if (err) {
                            err = 42;
                            goto E2FIN;
                        }
                        buf[j2][l] = tmp;                      
                        break;

                    case 125:                 /* negbin(alpha,gamma,k) */
                        tmp2 = buf[j1][l];
                        tmp1 = buf[j2][l];
                        tmp  = buf[j3][l];
                        err = l_negbin(tmp,tmp1,tmp2,0,dval);
                        if (err) {
                            err = 43;
                            goto E2FIN;
                        }
                        buf[j3][l] = dval[0];                  
                        break;

                    case 126:                 /* gclen(lon1,lat1,lon2,lat2) */

                        tmp3 = buf[j1][l];
                        tmp2 = buf[j2][l];
                        tmp1 = buf[j3][l];
                        tmp  = buf[j4][l];
                        buf[j4][l] = geod_distance(tmp,tmp1,tmp2,tmp3) * 180.0 / Pi;
                        break;

                    default:
                        err = 99;
                        goto E2FIN;
                }
            }
            if (atyp == 114 || atyp == 120 || atyp == 121 || atyp == 124)
                jj -= 1;
            else if (atyp == 107 || atyp == 115 || atyp == 122 || atyp == 125)
                jj -= 2;
            else if (atyp == 126)
                jj -= 3;
        }
        else if (atyp < 300) {  /* logical operators */

            if (atyp == 200) {      /* not */

                j = jj - 1;
        
                for (l = 0; l < noc; ++l) {

                    tmp = buf[j][l];
                    if (fabs(tmp) > EPSI2)
                        buf[j][l] = 0.0;       
                    else
                        buf[j][l] = 1.0;       
                }
            }
            else if (atyp == 209) {      /* if, three arguments */
                if (jj < 3)  
                    gerr_exit(75);
                jj -= 2;
                j = jj - 1;
                for (l = 0; l < noc; ++l) {
                    if (fabs(buf[j][l]) > EPSI2) 
                        buf[j][l] = buf[jj][l];
                    else
                        buf[j][l] = buf[jj+1][l];
                }
            }
            else {                       /* two arguments  */

                jj--;   
                j = jj - 1;
                j1 = j + 1;

                for (l = 0; l < noc; ++l) {

                    tmp  = buf[j][l];
                    tmp1 = buf[j1][l];

                    switch (atyp) {

                        case 201:                 /* and */
                            if (fabs(tmp) > EPSI2 && fabs(tmp1) > EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 202:                 /* or */
                            if (fabs(tmp) > EPSI2 || fabs(tmp1) > EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 203:                 /* eq */
                            if (fabs(tmp - tmp1) <= EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 204:                 /* ne */
                            if (fabs(tmp - tmp1) > EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 205:                 /* lt */
                            if (tmp + EPSI2 < tmp1)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 206:                 /* le */
                            if (tmp <= tmp1 + EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 207:                 /* gt */
                            if (tmp > tmp1 + EPSI2)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        case 208:                 /* ge */
                            if (tmp + EPSI2 >= tmp1)
                                buf[j][l] = 1.0;       
                            else
                                buf[j][l] = 0.0;       
                            break;

                        default:
                            err = 99;        
                            goto E2FIN;
                    }
                }
            }
        }
        else if (atyp < OPT2S) {  /* some other operators */

            switch (atyp) {

                case 300:                 /* ndf */
                case 301:                 /* nd */
                case 302:                 /* ndi */
                case 303:                 /* mr */
                    j =  jj - 1; 
                    for (l = 0; l < noc; ++l) {
                        tmp = buf[j][l];
                        if (atyp == 300)
                            buf[j][l] = dnf(tmp);  
                        else if (atyp == 301)
                            buf[j][l] = cdnf(tmp);  
                        else if (atyp == 302) {
                            if (tmp <= 0.0 || tmp >= 1.0) {
                                err = 11;
                                goto E2FIN;
                            }
                            buf[j][l] = cdnif1(tmp);  
                        }
                        else if (atyp == 303)
                            buf[j][l] = cdnm(tmp);  
                    }
                    break;

                case 304:                 /* td */
                case 305:                 /* chd */
                case 412:                 /* num(a,d) */
                    jj--;   
                    j = jj - 1;
                    for (l = 0; l < noc; ++l) {
                        tmp1 = buf[jj][l];
                        tmp  = buf[j][l];
                        if (atyp == 304)
                            buf[j][l] = cdtf(tmp,(int)tmp1);
                        else if (atyp == 305) {
                            if (tmp < 0.0 || (int)tmp1 < 1) {
                                err = 13;
                                goto E2FIN;
                            }
                            buf[j][l] = cdchif(tmp,(int)tmp1);
                        }
                        else if (atyp == 412)  
                            buf[j][l] = tmp + tmp1 * (double)l;
                    }
                    break;

                case 306:                 /* fd */

                    jj -= 2;
                    j = jj - 1;
                    j1 = jj + 1;

                    for (l = 0; l < noc; ++l) {
                        tmp2 = buf[j1][l];
                        tmp1 = buf[jj][l];
                        tmp  = buf[j][l];
                        if (atyp == 306) {
                            if (tmp < 0.0 || (int)tmp1 < 1 || (int)tmp2 < 1) {
                                err = 14;
                                goto E2FIN;
                            }
                            buf[j][l] = cdff(tmp,(int)tmp1,(int)tmp2);
                        }
                    }
                    break;

                case 400:                 /* dummy variable operators */

                  jj--;
                  jj -= narg;

                  for (l = 0; l < noc; ++l) {
                    tmp = buf[jj][l];
                    buf[jj][l] = 0.0;
                    
                    tmp1 = 1.0;
                    k = 1;
                    while (k <= narg) {
                        tmp2 = buf[jj + k][l];
                        if (fabs(tmp - tmp2) < EPSI2) {
                            buf[jj][l] = tmp1;
                            break;
                        }
                        if (++k > narg)
                            break;

                        tmp3 = buf[jj + k][l];
                        if ((int)tmp3 == INTMAX - 1) {
                            tmp1 = -1.0;
                            k++;
                        }
                        else if ((int)tmp3 == INTMAX) {
                            k++;
                            if (tmp2 <= tmp && tmp <= buf[jj + k][l]) {
                                 buf[jj][l] = tmp1;
                                 break;
                            }
                        }
                    }
                  }
                  jj++;
                  break;

                case 402:                 /* jul */
                    jj -= 2;
                    j = jj - 1;

                    for (l = 0; l < noc; ++l) {

                        n2 = (int)buf[j+2][l];     /* day */
                        n1 = (int)buf[j+1][l];     /* month */
                        n0 = (int)buf[j][l];       /* year */

                        k  = n2 - 32075 + 1461 *
                             (n0 + 4800 + (n1 - 14) / 12) / 4
                              + 367 * (n1 - 2 - (n1 - 14) / 12 * 12) / 12 
                              - 3 * ((n0 + 4900 + (n1 - 14) / 12) / 100) / 4;

                        buf[j][l] = (double)k;
                    }
                    break;

                case 404:                 /* rd(a,b) */
                    jj -= 1;
                    j = jj - 1;

                    for (l = 0; l < noc; ++l) {
                        tmp  = buf[j][l];
                        tmp1 = buf[j+1][l];
                        buf[j][l] = tmp + random1() * (tmp1 - tmp);
                    }
                    break;

                case 409:                 /* july */
                case 410:                 /* julm */
                case 411:                 /* juld */
                    j = jj - 1;

                    for (l = 0; l < noc; ++l) {
                        k = (int) buf[j][l];    /* julian time */
                        k += 68569;
                        nn = (4 * k) / 146097;
                        k -= (146097 * nn + 3) / 4;
                        n0 = 4000 * (k + 1) / 1461001;
                        k += 31 - (1461 * n0) / 4;
                        n1 = (80 * k) / 2447;
                        n2 = k - (2447 * n1) / 80;
                        k = n1 / 11;
                        n1 += 2 - 12 * k;
                        n0 += 100 * (nn - 49) + k;

                        if (atyp == 409)
                            buf[j][l] = (double)n0;
                        else if (atyp == 410)
                            buf[j][l] = (double)n1;
                        else
                            buf[j][l] = (double)n2;
                    }
                    break;   

                case 415:                 /* row(X) */
                case 416:                 /* col(X) */
                case 417:                 /* begin of row/col expression */
                    j1 = jj - 1;
                    if (atyp == 415) {
                        for (l = 0; l < noc; ++l)  
                            buf[j1][l] = (double)mrow;
                    }
                    else if (atyp == 416) {
                        for (l = 0; l < noc; ++l)  
                            buf[j1][l] = (double)mcol;
                    }
                    mrow = mcol = 1;
                    break;

                case 418:                     /* exists(A) */
                    if (jj >= STACKLEVEL) {
                        err = 97;
                        goto E2FIN;
                    }
                    for (l = 0; l < noc; ++l) {
                        if (val[i] < 0.0)
                            buf[jj][l] = 0.0;
                        else
                            buf[jj][l] = 1.0;
                    }
                    jj++;
                    break;
        
                case 419:                   /* rdp1(z) */

                    j1 = jj - 1;
                    tmp = buf[j1][0];
                    if (tmp <= 0.0 || tmp > ExpMax) {
                        err = 46;
                        goto E2FIN;
                    }
                    for (l = 0; l < noc; ++l)  
                        buf[j1][l] = (double)rdp1(tmp);
                    break;

                case 450:                   /* str(n,m) */

                    j1 = jj - 2;
                    for (l = 0; l < noc; ++l)  
                        buf[j1][l] = 0.0;               
                    jj--;
                    break;

                default:
                    err = 99;
                    goto E2FIN;
            }
        }
        else if (atyp < OPT2A) {  /* ### spatial operators */

            if (SDVarDef == 0)
                return(49);

            j1 = jj - 1;
            j2 = jj - 2;
            
            for (l = 0; l < noc; ++l) {

                switch (atyp) {

                    case 601:                 /* sdp(x,y) */

                        if ((int)get_data(SDVarSDTyp,l) != 3)
                            buf[j2][l] = 0.0;
                        else {
                            if ((nn = sd_getdata(l,0,1,1)) < 1)
                                return(50);

                            tmp1 = buf[j1][l];
                            tmp  = buf[j2][l];
                            buf[j2][l] = (double)g_inpoly(nn - 1,SDVarX,SDVarY,tmp,tmp1);
                        }
                        break;

                    default:
                        err = 99;
                        goto E2FIN;
                }
            }
            if (atyp == 127)
                jj -= 1;
        }
        else if (atyp >= 1880 && atyp < OPT2B) {  /* block mode operators */

            if (atyp == 1880 || atyp == 1881) {
                j = jj - 1;
                y = x = buf[j][0];
                for (l = 1; l < noc; ++l) {
                    tmp = buf[j][l];
                    if (x > tmp)
                        x = tmp;
                    if (y < tmp)
                        y = tmp;
                }
                if (atyp == 1880) {                 /* bmin */
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = x;
                }
                else if (atyp == 1881) {            /* bmax */
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = y;
                }
            }
            else if (atyp == 1882 || atyp == 1883) {    /* bfa, bfr */

                j = jj - 1;
                k = 0;
                for (l = 0; l < noc; ++l) {
                    if (buf[j][l] > EPSI1)
                        k++;
                }   
                for (l = 0; l < noc; ++l)  
                    buf[j][l] = (double)k;
                              
                if (atyp == 1883 && k > 0 && noc > 0) {  
                    for (l = 0; l < noc; ++l)  
                        buf[j][l] /= (double)noc;
                }
            }
            else if (atyp >= 1890 && atyp <= 1894) {

                j = jj;
                if (j >= STACKLEVEL) {
                    err = 97;
                    goto E2FIN;
                }
                for (l = 0; l < noc; ++l) {

                    switch (atyp) {
                        case 1890:  buf[j][l] = (double)noc;        /* bnrec */
                                    break;
                        case 1891:  buf[j][l] = (double)(l + 1);    /* brec */
                                    break;
                        case 1892:                                  /* bfirst */
                        case 1893:  buf[j][l] = 0.0;                /* blast */
                                    break;

                        case 1894:  buf[j][l] = (double)(bnum + 1); /* bnum */
                                    break;
                    }
                }
                if (atyp == 1892)
                    buf[j][0] = 1.0;
                else if (atyp == 1893)
                    buf[j][noc - 1] = 1.0;
                jj++;
            }
            else if (atyp == 1899) {        /* brd */

                j = jj;
                if (j >= STACKLEVEL) {
                    err = 97;
                    goto E2FIN;
                }
                x = random1();
                for (l = 0; l < noc; ++l)  
                    buf[j][l] = x;
                jj++;
            }
        }
        else if (atyp >= OPT2A && atyp < OPT2B) {  /* type 2 operators */
   
            switch (atyp) {

                case 1500:                      /* lag(V,n) */
              
                    j  = jj - 1;
                    j1 = jj - 2;

                    for (l = 0; l < noc; ++l) {
                        l1 = (int)buf[j][l] + l;
                        if (l1 < 0 || l1 >= noc)
                            buf[j][l] = 0.0;
                        else
                            buf[j][l] = buf[j1][l1];
                    }
                    for (l = 0; l < noc; ++l)
                        buf[j1][l] = buf[j][l];
                    jj--;         
                    break;

                case 1501:                      /* sum */

                    j = jj - 1;
                    x = 0.0;
                    for (l = 0; l < noc; ++l)
                        x += buf[j][l];

                    for (l = 0; l < noc; ++l)
                        buf[j][l] = x;
                    break;

                case 1502:                      /* mean */
                    j = jj - 1;
                    x = 0.0;
                    for (l = 0; l < noc; ++l)
                        x += buf[j][l];
                    x /= (double)noc;
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = x;
                    break;

                case 1503:                      /* std */

                    j = jj - 1;
                    x = y = 0.0;
                    for (l = 0; l < noc; ++l) {
                        tmp = buf[j][l];
                        x += tmp;
                        y += tmp * tmp; 
                    }
                    tmp = y - x * x / (double)noc;
                    if (noc > 1) {
                        tmp /= (double)(noc - 1);
                        if (tmp > 0.0)
                            tmp = sqrt(tmp);
                        else
                            tmp = 0.0;
                    }
                    else
                        tmp = 0.0;

                    for (l = 0; l < noc; ++l)
                        buf[j][l] = tmp;
                    break;

                case 1504:                      /* cum */

                    j = jj - 1;
                    x = 0.0;
                    for (l = 0; l < noc; ++l) {
                        x += buf[j][l];
                        buf[j][l] = x;
                    }
                    break;

                case 1505:                      /* cd  */

                    j = jj - 1;
                    x = buf[j][0];
                    j1 = 0;
                    for (l = 1; l < noc; ++l) {
                        y = buf[j][l];
                        if (y > x) {
                            while (j1 < l) {
                                buf[j][j1] = (double)l / (double)noc;
                                j1++;
                            }
                            x = y;
                        }
                    }
                    while (j1 < noc) {
                        buf[j][j1] = 1.0;
                        j1++;
                    }
                    break;

                case 1506:                      /* cdv */

                    j = jj - 1;
                    tmp = x = buf[j][0];
                    j1 = 0;
                    for (l = 1; l < noc; ++l) {
                        y = buf[j][l];
                        if (y > x) {
                            while (j1 < l) {
                                buf[j][j1] = tmp;
                                j1++;
                            }
                            x = y;
                        }
                        tmp += y;
                    }
                    while (j1 < noc) {
                        buf[j][j1] = tmp;
                        j1++;
                    }
                    if (tmp > 0.0) {
                        for (l = 0; l < noc; ++l)
                            buf[j][l] /= tmp;
                    }
                    break;

                case 1507:                     /* vmin */

                    j = jj - 1;
                    x = buf[j][0];
                    for (l = 1; l < noc; ++l) {
                        y = buf[j][l];
                        if (x > y)
                            x = y;
                    }
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = x;
                    break;

                case 1508:                     /* vmax */

                    j = jj - 1;
                    x = buf[j][0];
                    for (l = 1; l < noc; ++l) {
                        y = buf[j][l];
                        if (x < y)
                            x = y;
                    }
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = x;
                    break;

                case 1509:                      /* suc */

                    j = jj;
                    if (j >= STACKLEVEL) {
                        err = 97;
                        goto E2FIN;
                    }
                    k = (int)val[i];
                    for (l = 0; l < noc; ++l) {
                        if (l >= noc - 1)
                            buf[j][l] = 0.0;
                        else
                            buf[j][l] = get_data(k,l + r1 + 1);
                    }
                    jj++;         
                    break;

                case 1510:                      /* pre */

                    j = jj;
                    if (j >= STACKLEVEL) {
                        err = 97;
                        goto E2FIN;
                    }
                    k = (int)val[i];
                    for (l = 0; l < noc; ++l) {
                        if (l == 0)
                            buf[j][l] = 0.0;
                        else
                            buf[j][l] =  get_data(k,l + r1 - 1);
                    }
                    jj++;          
                    break;

                case 1511:                      /* mav(V,n) */
              
                    j  = jj - 1;
                    j1 = jj - 2;
                    l1 = (int)buf[j][0];
                    if (l1 > 0) {
                        for (l = 0; l < noc; ++l) {
                            x = y = 0.0;
                            for (k = l - l1; k <= l + l1; ++k) {
                                if (k >= 0 && k < noc) {
                                    x += buf[j1][k];
                                    y += 1.0;
                                }
                            }
                            buf[j][l] = x / y;
                        }
                        for (l = 0; l < noc; ++l)
                            buf[j1][l] = buf[j][l];
                    }
                    jj--;         
                    break;

                case 1512:                      /* vdif(X,Y) */
              
                    j  = jj - 1;
                    j1 = jj - 2;
                    for (l = 0; l < noc; ++l) {
                        x = buf[j1][l];
                        y = 1.0;
                        for (l1 = 0; l1 < noc; ++l1) {
                            if (fabs(x - buf[j][l1]) <= EPSI1) {
                                y = 0.0;  
                                break;
                            }
                        }   
                        buf[j1][l] = y;             
                    }
                    jj--;         
                    break;

                case 1513:                      /* cmean(a,f) */

                    jj--;
                    j = jj - 1;
                    j1 = jj;

                    tmp = x = 0.0;
                    for (l = noc - 1; l >= 0; --l) {
                        tmp += buf[j1][l];
                        x += buf[j][l] * buf[j1][l];
                        if (tmp > 0.0)
                            buf[j][l] = x / tmp;
                        else
                            buf[j][l] = 0.0;
                    }
                    break;

                case 1514:                      /* cmean1(a,f) */

                    jj--;
                    j = jj - 1;
                    j1 = jj;

                    tmp = x = 0.0;
                    for (l = 0; l < noc; ++l) {
                        tmp += buf[j1][l];
                        x += buf[j][l] * buf[j1][l];
                        if (tmp > 0.0)
                            buf[j][l] = x / tmp;
                        else
                            buf[j][l] = 0.0;
                    }
                    break;

                case 1520:                     /* sort(V) */
                case 1521:                     /* snum(V) */
                case 1522:                     /* rank(V) */
                case 1523:                     /* ndv(V)  */
                case 1524:                     /* ndv1(V,a) */
                case 1525:                     /* ndv2(V,a,b) */

                    j = jj - 1;
                    if (atyp == 1524) {
                        j1 = j;
                        j--;
                        jj--;
                    }
                    else if (atyp == 1525) {
                        j2 = j;
                        j1 = j2 - 1;
                        j = j1 - 1;
                        jj -= 2;
                    }

                    /* allocate a pointer for sorting */

                    if (!(BPtr = (int *)calloc(noc,sizeof(int)))) { 
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(int));

                    if (!(BVal = (float *)calloc(noc,sizeof(float)))) { 
                        free((char *)BPtr);
                        memrq(-noc,sizeof(int));
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(float));

                    for (l = 0; l < noc; ++l) {
                        BPtr[l] = l;
                        BVal[l] = buf[j][l];
                    }
                    qsort((char *)BPtr,noc,sizeof(int),edcomp);

                    if (atyp == 1520) {          /* sort */
                        for (l = 0; l < noc; ++l)
                            buf[j][l] = (double)BVal[BPtr[l]];
                    }
                    else if (atyp == 1521) {     /* snum */
                        k = 0;
                        for (l = 0; l < noc; ++l) {
                            k++;
                            buf[j][BPtr[l]] = (double)k;
                        }
                    }
                    else if (atyp == 1522) {    /* rank */

                        ogrp = BVal[BPtr[0]];
                        k = 0;
                        tmp = x = 0.0;
                        for (l = 0; l < noc; ++l) {
                            grp = BVal[BPtr[l]];
                            if (grp != ogrp) {
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j][BPtr[l1]] = tmp / (double)(l - k);
                                k = l;
                                ogrp = grp;
                                tmp = 0.0;
                            }
                            x += 1;
                            tmp += x;  
                        }
                        l = noc;
                        for (l1 = k; l1 < noc; ++l1)  
                            buf[j][BPtr[l1]] = tmp / (double)(l - k);
                    }
                    else {                       /* ndv, ndv1, ndv2 */
                        k = 0;
                        tmp = (double)BVal[BPtr[0]] - 1.0;
                        for (l = 0; l < noc; ++l) {
                            tmp1 = (double)BVal[BPtr[l]];
                            if (tmp1 != tmp) {
                                if (atyp == 1523)
                                    k++;
                                else if (tmp1 >= buf[j1][l]) {
                                    if (atyp == 1524)
                                        k++;
                                    else if (tmp1 <= buf[j2][l])  
                                        k++;
                                }
                                tmp = tmp1;
                            }
                        }
                        for (l = 0; l < noc; ++l)  
                            buf[j][l] = (double)k;
                    }

                    free((char *)BPtr);
                    memrq(-noc,sizeof(int));
                    free((char *)BVal);
                    memrq(-noc,sizeof(float));
                    break;

                case 1526:                     /* strsp(S) */

                    k = aval - VOFFS;
                    if (k < 0 || k >= NVAR)
                        return(-1);
                    if (VTyp[k] != 1)
                        return(47);
            
                    BPtrSBLen = -VSLen[k];
                    BPtrSBuf  =  VDPtr[k];                

                    j = jj - 1;

                    /* allocate a pointer for sorting */

                    if (!(BPtr = (int *)calloc(noc,sizeof(int)))) { 
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(int));

                    for (l = 0; l < noc; ++l)  
                        BPtr[l] = l;
 
                    qsort((char *)BPtr,noc,sizeof(int),es_comp);
             
                    buf[j][BPtr[0]] = tmp = 1.0;       
                    p = BPtrSBuf + BPtr[0] * BPtrSBLen;
                    for (l = 1; l < noc; ++l) {
                        q = BPtrSBuf + BPtr[l] * BPtrSBLen;
                        if (strncmp(p,q,BPtrSBLen))
                            tmp += 1.0;
                        buf[j][BPtr[l]] = tmp;        
                        p = q;
                    }
                    free((char *)BPtr);
                    memrq(-noc,sizeof(int));
                    break;

            
                case 1527:                     /* quant(X,p) */

                    j1 = jj - 1;
                    j  = jj - 2;
                    jj-= 1;
                    tmp3 = buf[j1][0];

                    if (tmp3 <= 0.0 || tmp3 >= 1.0)  
                        return(15);

                    /* allocate a pointer for sorting */

                    if (!(BPtr = (int *)calloc(noc,sizeof(int)))) { 
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(int));

                    if (!(BVal = (float *)calloc(noc,sizeof(float)))) { 
                        free((char *)BPtr);
                        memrq(-noc,sizeof(int));
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(float));

                    for (l = 0; l < noc; ++l) {
                        BPtr[l] = l;
                        BVal[l] = buf[j][l];
                    }
                    qsort((char *)BPtr,noc,sizeof(int),edcomp);

                    if (noc == 1)
                        tmp = (double)BVal[0];
                    else {
                        tmp1 = tmp3 * (double)(noc - 1);
                        tmp2 = floor(tmp1);
                        tmp = (1.0 - (tmp1 - tmp2)) *
                              (double)BVal[BPtr[(int)tmp2]] +
                              (tmp1 - tmp2) * (double)BVal[BPtr[(int)tmp2 + 1]];  
                    }
                    for (l = 0; l < noc; ++l)
                        buf[j][l] = tmp;                     

                    free((char *)BPtr);
                    memrq(-noc,sizeof(int));
                    free((char *)BVal);
                    memrq(-noc,sizeof(float));
                    break;
/* ## */
                case 1528:                     /* quant1(X,Z,p,m) */

                    j3 = jj - 1;
                    j2 = jj - 2;
                    j1 = jj - 3;
                    j  = jj - 4;
                    jj -= 3;
                    tmp3 = buf[j2][0];
                    tmp4 = buf[j3][0];

                    if (tmp3 <= 0.0 || tmp3 >= 1.0)  
                        return(15);

                    /* allocate a pointer for sorting */

                    if (!(BPtr = (int *)calloc(noc,sizeof(int)))) { 
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(int));

                    if (!(BVal = (float *)calloc(noc,sizeof(float)))) { 
                        free((char *)BPtr);
                        memrq(-noc,sizeof(int));
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(float));

                    if (!(ZVal = (short *)calloc(noc,sizeof(short)))) { 
                        free((char *)BPtr);
                        memrq(-noc,sizeof(int));
                        free((char *)BVal);
                        memrq(-noc,sizeof(float));
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(short));

                    for (l = 0; l < noc; ++l) {
                        BPtr[l] = l;
                        BVal[l] = buf[j][l];
                        if (buf[j1][l] != 0.0)
                            ZVal[l] = 1;           
                    }
                    qsort((char *)BPtr,noc,sizeof(int),ed1comp);

                    tmp1 = tmp = 1.0;
                    y = x = BVal[BPtr[0]];
                    z = (double)noc;
                    l = 0;
                    while (l < noc) {
                        ptr = BPtr[l];
                        if (ZVal[ptr] != 0) {
                            tmp1 = tmp;
                            y = x;
                            x = BVal[ptr];
                            e = 1.0;
                            for (l1 = l + 1; l1 < noc; ++l1) {
                                ptr1 = BPtr[l1];
                                if (BVal[ptr1] == x && ZVal[ptr1] != 0) {
                                    e += 1.0;
                                    l++;
                                }
                                else
                                    break;
                            }
                            tmp *= (1.0 - e / z);
                            if (tmp <= tmp3) {
                                if (tmp1 > tmp)
                                    tmp4 = (tmp1 - tmp3) * (x - y) / (tmp1 - tmp) + y;
                                break;
                            }
                            z -= e;
                        }
                        else 
                            z -= 1.0;
                        l++;
                    }
                    for (l = 0; l < noc; ++l) 
                        buf[j][l] = tmp4;          

                    free((char *)BPtr);
                    memrq(-noc,sizeof(int));
                    free((char *)BVal);
                    memrq(-noc,sizeof(float));
                    free((char *)ZVal);
                    memrq(-noc,sizeof(short));
                    break;


                case 1530:                     /* cnteq */
                case 1531:                     /* cntne */
                case 1532:                     /* cntlt */
                case 1533:                     /* cntgt */
                case 1534:                     /* cntle */
                case 1535:                     /* cntge */

                    j = jj - 1;
                    j1 = j - 1;

                    k = 0;
                    for (l = 0; l < noc; ++l) { 

                        tmp = buf[j1][l] - buf[j][l];

                        switch (atyp) {
                            case 1530:  if (fabs(tmp) < EPSI2) k++; break;
                            case 1531:  if (fabs(tmp) >= EPSI2) k++; break;
                            case 1532:  if (tmp + EPSI2 < 0.0) k++; break;
                            case 1533:  if (tmp > EPSI2) k++; break;
                            case 1534:  if (tmp <= EPSI2) k++; break;
                            case 1535:  if (tmp + EPSI2 >= 0.0) k++; break;
                        }
                    }
                    for (l = 0; l < noc; ++l)  
                        buf[j1][l] = (double)k;
                    jj--;        
                    break;
   
                case 1540:                     /* gcnt  */
                case 1541:                     /* grec  */
                case 1542:                     /* gsn   */
                case 1543:                     /* gfirst*/
                case 1544:                     /* glast */

                case 1550:                     /* gsum  */
                case 1551:                     /* gmean */
                case 1552:                     /* gstd  */
                case 1553:                     /* gmin  */
                case 1554:                     /* gmax  */

                    j = jj - 1;
                    if (atyp >= 1550 && atyp <= 1555)  
                        j1 = j - 1;
                    else               
                        j1 = j;

                    x = buf[j][0];
                    k = t = 0;                            
                    z = tmp = tmp1 = 0.0;
                    if (atyp == 1553 || atyp == 1554)
                        tmp = buf[j1][0];

                    for (l = 0; l <= noc; ++l) { 

                        if (l < noc)  
                            y = buf[j][l];

                        if (fabs(x - y) > EPSI2 || l == noc) {

                            switch (atyp) {

                              case 1540:
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = (double)t;
                                break;

                              case 1541:
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = (double)(l1 + 1 - k);
                                break;

                              case 1542:
                                z += 1.0;
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = z;             
                                break;

                              case 1543:
                                buf[j1][k] = 1.0;
                                for (l1 = k + 1; l1 < l; ++l1)  
                                    buf[j1][l1] = 0.0;
                                break;

                              case 1544:
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = 0.0;
                                buf[j1][l - 1] = 1.0;
                                break;

                              case 1550:
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = tmp;
                                break;

                              case 1551:
                                tmp /= (double)t;
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = tmp;
                                tmp = 0.0;
                                break;

                              case 1552:
                                z = tmp1 - tmp * tmp / (double)t;
                                if (t > 1) {
                                    z /= (double)(t - 1);
                                    if (z > 0.0)
                                        z = sqrt(z);
                                    else
                                        z = 0.0;
                                }
                                else
                                    z = 0.0;

                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = z;
                                tmp = tmp1 = 0.0;
                                break;

                              case 1553:
                              case 1554:
                                for (l1 = k; l1 < l; ++l1)  
                                    buf[j1][l1] = tmp;
                                if (l < noc)
                                    tmp = buf[j1][l];
                                break;
                            }
                            x = y;
                            k = l;
                            t = 0;
                            if (atyp >= 1550 && atyp <= 1552)  
                                tmp = tmp1 = 0.0;
                        }
                        if (l >= noc)
                            break;

                        if (atyp >= 1550 && atyp <= 1554) {
                            z = buf[j1][l];

                            switch (atyp) {
                                case 1550:
                                case 1551: tmp += z; break;
                                case 1552: tmp += z; tmp1 += z * z; break;
                                case 1553: if (tmp > z) tmp = z; break;
                                case 1554: if (tmp < z) tmp = z; break;
                                default: break;
                            }
                        }
                        t++;
                    }  
                    if (atyp >= 1550 && atyp <= 1554)  
                        jj--;          
                    break;
                
                case 1555:                     /* gsort(V,G) */
                case 1556:                     /* gndv(V,G) */
                case 1557:                     /* gndv1(V,a,G) */

                    j = jj - 1;
                    j1 = j - 1;
                    if (atyp == 1557) {
                        j2 = j1;
                        j1--;
                    }

                    /* allocate a pointer for sorting */

                    if (!(BPtr = (int *)calloc(noc,sizeof(int)))) { 
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(int));

                    if (!(BVal = (float *)calloc(noc,sizeof(float)))) { 
                        free((char *)BPtr);
                        memrq(-noc,sizeof(int));
                        err = 95;
                        goto E2FIN;
                    }
                    memrq(noc,sizeof(float));

                    x = buf[j][0];
                    k = t = 0;                            

                    for (l = 0; l <= noc; ++l) { 

                        if (l < noc)  
                            y = buf[j][l];

                        if (fabs(x - y) > EPSI2 || l == noc) {
      
                            qsort((char *)BPtr,t,sizeof(int),edcomp);

                            if (atyp == 1555) {
                                for (l1 = k; l1 < l; ++l1)
                                    buf[j1][l1] = (double)BVal[BPtr[l1 - k]];
                            }
                            else {
                                n0 = 0;
                                tmp = (double)BVal[BPtr[0]] - 1.0;
                                for (l1 = 0; l1 < t; ++l1) {
                                    tmp1 = (double)BVal[BPtr[l1]];
                                    if (tmp1 != tmp) {
                                        if (atyp == 1556)
                                            n0++;
                                        else if (tmp1 >= buf[j2][0])  
                                            n0++;
                                    }
                                    tmp = tmp1;
                                }
                                for (l1 = k; l1 < l; ++l1)
                                    buf[j1][l1] = (double)n0;
                            }
                            x = y;
                            k = l;
                            t = 0;
                        }
                        if (l >= noc)
                            break;
                        BPtr[t] = t;
                        BVal[t] = buf[j1][l];
                        t++;
                    }  
                    jj--;          
                    if (atyp == 1557)
                        jj--;

                    free((char *)BPtr);
                    memrq(-noc,sizeof(int));
                    free((char *)BVal);
                    memrq(-noc,sizeof(float));
                    break;

                case 1600:                      /* change  */
                case 1601:                      /* cntch  */
                case 1602:                      /* ccntch  */
                case 1603:                      /* lagch   */

                    j = jj - 1;
                    k = 1;
                    t = 0;
                    x = buf[j][0];
                    if (atyp == 1603)
                        buf[j][0] = 0.0;
                    else
                        buf[j][0] = 1.0;
                    for (l = 1; l < noc; ++l) {
                        if (fabs(buf[j][l] - x) > EPSI2) {
                            x = buf[j][l];
                            k++;
                            t = 0;
                            if (atyp == 1600)
                                buf[j][l] = 1.0;
                        }
                        else {
                            t++;
                            x = buf[j][l];
                            if (atyp == 1600)
                                buf[j][l] = 0.0;
                        }
                        if (atyp == 1602)
                            buf[j][l] = (double)k;
                        else if (atyp == 1603)
                            buf[j][l] = (double)t;
                    }
                    if (atyp == 1601) {
                        for (l = 0; l < noc; ++l)
                            buf[j][l] = (double)k;
                    }
                    break;
     
                case 1697:                      /* X(i) */

                    j1 = jj - 2;
                    n1 = (int)buf[jj - 1][0];   
                    if (n1 < 1) {
                        err = 45;
                        goto E2FIN;
                    }
                    kk = iabs(esidx);
                    tmp = 0.0;
                    if (kk >= COFFS)
                        ;            
                    else if (kk >= VOFFS) {
                        if (n1 > 1) {
                            err = 45;
                            goto E2FIN;
                        }
                        for (l = 0; l < noc; ++l) {
                            buf[j1][l] = get_data(kk - VOFFS,l % NOC);
                            if (esidx < 0)
                                buf[j1][l] *= -1.0;
                        }
                    }
                    else if (kk >= MOFFS) {
                        n0 = kk - MOFFS;
                        row = MatRow[n0];
                        col = MatCol[n0];
                        if (n1 > col) {
                            err = 45;
                            goto E2FIN;
                        }   
                        for (l = 0; l < noc; ++l) {
                            buf[j1][l] = MatVal[n0][(l % row) * col + n1];
                            if (esidx < 0)
                                buf[j1][l] *= -1.0;
                        }
                    }
                    else if (kk >= NLOFFS) {
                        n0 = kk - NLOFFS;
                        col = NLNV[n0];
                        if (n1 > col) {
                            err = 45;
                            goto E2FIN;
                        }   
                        n1--;
                        for (l = 0; l < noc; ++l) {
                            buf[j1][l] = get_data(NLVIdx[n0][n1],l % NOC);
                            if (esidx < 0)
                                buf[j1][l] *= -1.0;
                        }
                    }
                    jj--;
                    break;

                case 1698:                      /* row vector: <...>'  */
                case 1699:                      /* col vector: <...>   */
                                    /* use current column number: ncol */

                    /****
                    for (k = 0; k < jj; ++k) {
                        printf1("kkk=%2d : bufk ",k);
                        for (l = 0; l < noc; ++l)
                            printf1("%lg ",buf[k][l]);
                        printf1("\n");
                    }
                    ****/

                    k = jj - narg;
                    if (atyp == 1698) {
                        mcol = imax(mcol,narg);
                        tmp = buf[k + (ncol % narg)][0];
                        for (l = 0; l < noc; ++l)
                            buf[k][l] = tmp;
                    }
                    else {
                        mrow = imax(mrow,narg);
                        for (l = r2 - 1; l >= r1; --l)  
                            buf[k][l - r1] = buf[k + (l % narg)][0];
                    }
                    jj -= (narg - 1);
                    break;

                case 1700:                      /* sma[](X) */
                case 1701:                      /* smd[](X) */

                    err = eval_specf(atyp,vn,noc,jj,buf);
                    if (err)
                        goto E2FIN;
                    break;

                default:
                    err = 99;
                    goto E2FIN;

            }
        }
        else if (atyp < IAOFFS) {  /* constante and operators without
                                      arguments */

            j = jj;
            if (j >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            else {
                for (l = 0; l < noc; ++l) {

                    switch (atyp) {

                        case 1900:            /* case */
      
                            if (GDFlg) 
                                buf[j][l] = (double)(l + GDNRec + 1);
                            else
                                buf[j][l] = (double)(l + r1 + 1);
                            break;
            
                        case 1901:            /* nvar */
                            buf[j][l] = (double)NVAR;               
                            break;

                        case 1902:            /* pi */

                            buf[j][l] = Pi;
                            break;
            
                        case 1903:            /* rd */

                            buf[j][l] = random1();
                            break;
            
                        case 1904:            /* rdn */

                            buf[j][l] = normal();
                            break;

                        case 1905:            /* rdn1 */

                            buf[j][l] = normal1();
                            break;

                        case 1906:            /* nocdm */
                            buf[j][l] = (double)NOCDM;               
                            break;
            
                        case 1907:            /* noc */
                            buf[j][l] = (double)NOC;               
                            break;

                        case 1908:            /* bnoc */
                            buf[j][l] = (double)BNOC;               
                            break;

                        case 1910:            /* sdxmin */
                        case 1911:            /* sdxmax */
                        case 1912:            /* sdymin */
                        case 1913:            /* sdymax */
                        case 1914:            /* sdarea */
                        case 1915:            /* sdlen */

                            if (SDVarDef == 0)
                                return(49);

                            if ((nn = sd_getdata(l,0,1,1)) < 1)
                                return(50);
        
                            if (atyp == 1910)
                                buf[j][l] = g_pol_xmin(nn - 1,SDVarX,SDVarY);
                            else if (atyp == 1911)
                                buf[j][l] = g_pol_xmax(nn - 1,SDVarX,SDVarY);
                            else if (atyp == 1912)
                                buf[j][l] = g_pol_ymin(nn - 1,SDVarX,SDVarY);
                            else if (atyp == 1913)
                                buf[j][l] = g_pol_ymax(nn - 1,SDVarX,SDVarY);
                            else if (atyp == 1914) {
                                if ((int)get_data(SDVarSDTyp,l) != 3)
                                    buf[j][l] = 0.0;
                                else
                                    buf[j][l] = g_pol_area(nn - 1,SDVarX,SDVarY);
                            }
                            else if (atyp == 1915) {
                                n0 = (int)get_data(SDVarSDTyp,l);        
                                if (n0 == 2)
                                    buf[j][l] = g_line_len(nn - 1,SDVarX,SDVarY);
                                else if (n0 == 3)
                                    buf[j][l] = g_line_len(nn,SDVarX,SDVarY);
                                else
                                    buf[j][l] = 0.0;
                            }
                            break;

                        default:
                            buf[j][l] = 0.0;         
                            break;
                    }
                }
                jj++;
            }
        }
        else if (atyp < FAOFFS) {  /* intermediate function arguments */

            j = jj;
            if (j >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            else {
                k = atyp - IAOFFS;
                if (k < 0 || k >= FNPN) {
                    err = 17;
                    goto E2FIN;
                }
                tmp = FNPVal[k];
                for (l = 0; l < noc; ++l)  
                    buf[j][l] = tmp;                         
                jj++;
            }
        }
        else if (atyp < FAOFFS + MaxP) {  /* function arguments */

            j = jj;
            if (j >= STACKLEVEL) {
                err = 97;
                goto E2FIN;
            }
            else {
                k = atyp - FAOFFS;
                if (k < 0 || k >= FNPN) {
                    err = 17;
                    goto E2FIN;
                }
                tmp = FNArgVal[FNArgSPI[k]];               
                for (l = 0; l < noc; ++l)  
                    buf[j][l] = tmp;                         
                jj++;
            }
        }
        else {
            err = 99;
            goto E2FIN;
        }  
        if (err)  
            break;
          
        /* check for sign change */

        j = jj - 1;
        if (utyp < 0 && atyp < NLOFFS && j >= 0) {
            for (l = 0; l < noc; ++l)  
                buf[j][l] *= -1.0;                  
        }
    } 
    if (err == 0) {

        if (jj != 1)  
            gerr_exit(53);                 
            
        for (l = 0; l < noc; ++l) {
            if (vn >= 0)
                put_data(buf[0][l],vn,l + r1);
            else  
                res[l] = buf[0][l];
        }
    }
E2FIN:
    for (i = 0; i < STACKLEVEL; ++i) {
        free((char *)buf[i]);
        memrq(-noc-1,sizeof(double));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  edcomp()    compare function for sort operators.                        */

int edcomp(const void *arg1,const void *arg2)
{     
    float x; 

    x = BVal[*(int *)arg1] - BVal[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  es_comp()                                                               */
/*      Compare function for string variables.                              */

int es_comp(const void *arg1,const void *arg2)
{     
    register int i;
    register unsigned char *p1,*p2;

    p1 = (unsigned char *)(BPtrSBuf + *(int *)arg1 * BPtrSBLen);
    p2 = (unsigned char *)(BPtrSBuf + *(int *)arg2 * BPtrSBLen);

    for (i = 0; i < BPtrSBLen; ++i) {
        if (*p1 < *p2)   
            return(-1);
        else if (*p1 > *p2)  
            return(1);
        p1++;
        p2++;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  ed1comp()    compare function for sort operators, with censoring.       */

int ed1comp(const void *arg1,const void *arg2)
{     
    float x; 
    short z;

    x = BVal[*(int *)arg1] - BVal[*(int *)arg2];
    if (x > 0.0)
        return(1);
    else if (x < 0.0)
        return(-1);

    z = ZVal[*(int *)arg1] - ZVal[*(int *)arg2];

    if (z > 0.0)
        return(-1);
    else if (z < 0.0)
        return(1);
    return(0);
}


