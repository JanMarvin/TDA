/****************************************************************************/
/*  t_eval1                                                                 */
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
#include "t_eval3.h"
#include "t_mdat.h"
#include "t_sd.h"
#include "t_gm.h"
#include "t_map.h"

/* ------------------------------------------------------------------------ */
/*  check for identical definition in t_eval.c, t_eval2.c                   */

#define DDTYP     400   /* type number for dummy operators                  */
#define CCTYP      99   /* type number for ,, in dummy operators            */
#define COLTYP     98   /* type number for :  in dummy operators            */
#define RVECTYP  1698   /* type number for row vectors <...>'               */
#define CVECTYP  1699   /* type number for column vectors <...>             */

/* ------------------------------------------------------------------------ */
/*  functions in t_eval1.c                                                  */

int v_eval1(int n,short cnt,int *typ,double *val,int *idx,double *res,int vflag,
    int deriv,int ip,int ivp,int iflag);

/* ------------------------------------------------------------------------ */
/*  v_eval1(int n,short cnt,int *typ,double *val,int idx,double *res,       */
/*  ##      int vflag,int deriv,ip,ivp,iflag)                               */
/*                                                                          */
/*  Evaluate the expression defined in typ and val. Number of               */
/*  stack entries is cnt. If V variables, use n.th case of data matrix.     */
/*  If vflag != 0 use values from AVVAL[] for variables.                    */
/*                                                                          */
/*  If deriv = 1    calculate gradient in FNGrad[0][...]                    */
/*  If deriv = 2    calculate hessian  in FNHess[0][...]                    */
/*                                                                          */
/*  If iflag != 0 perform interval arithmetic. Only with deriv <= 1.        */
/*                                                                          */
/*  Evaluation begins with ip: typ[i],val[i], i = ip,...ip + cnt - 1        */
/*  ivp is pointer to buffer for evaluation.                                */
/*                                                                          */
/*  return 0 if ok, otherwise one of following error indicators:            */
/*                                                                          */
/*  -1  =   undefined X, Y  or C or V variable                              */
/*   1  =   undefined operator                                              */
/*   2  =   division by zero                                                */
/*   3  =   % operator finds non-integers                                   */
/*   4  =   exp out of range                                                */
/*   5  =   argument of log is too small (< LogMin)                         */
/*   6  =   negative argument in sqrt                                       */
/*   7  =   negative or zero argument for lgam, digam                       */
/*   8  =   arguments for icg out of range                                  */
/*   9  =   arguments for icb out of range                                  */
/*  10  =   incorrect range in tr(x,a,b)  (a > b)                           */
/*  11  =   argument for ndi() out ouf range [0,1]                          */
/*  12  =   incorrect arguments for TD                                      */
/*  13  =   incorrect arguments for CHD                                     */
/*  14  =   incorrect arguments for FD                                      */
/*  16  =   rdmn(i,A), i out of range                                       */
/*  17  =   undefined intermediate function parameter                       */
/*  18  =   undefined function argument                                     */
/*                                                                          */
/*  24  =   error: operator cannot be used for derivatives                  */
/*  27  =   error: integration interval not positive                        */
/*  28  =   error: no success in numerical integration                      */
/*  30  =   error: in number of Hermite integration points                  */
/*  31  =   error: misplaced integration variable (t)                       */
/*  32  =   error: cannot calculate incomplete gamma integral               */
/*  33  =   error: cannot calculate binomial coefficient                    */
/*  34  =   error: cannot calculate biv normal distribution                 */
/*  35  =   error: cannot init rdmn generator                               */
/*  36  =   error: wrong index in rdmn operator                             */
/*  37  =   error: need a square matrix                                     */
/*  38  =   error: exceeded max matrix dimensions                           */
/*  39  =   error: in number of arguments                                   */
/*  40  =   error: in mnv evaluation                                        */
/*  41  =   error: in correlation matrix                                    */
/*  42  =   error: in poisson operator                                      */
/*  43  =   error: in negbin operator                                       */
/*  44  =   error: can only calculate first derivative                      */
/*  45  =   error: matrix indices out of range.                             */
/*  46  =   error: argument of rdp1 out of range.                           */
/*  47  =   error: need reference to a string variable                      */
/*  48  =   error: in range of arguments                                    */
/*  49  =   error: spatial data required                                    */
/*  50  =   error: cannot read spatial data                                 */
/*                                                                          */
/*  90  =   syntax error (missing operator?)                                */
/*  92  =   out of stack space for derivatives                              */
/*  97  =   out of stack space for buf[]                                    */

#define E1SSIZ 1000         /* stack size for buf[] */

int v_eval1(int n,short cnt,int *typ,double *val,int *idx,double *res,int vflag,
    int deriv,int ip,int ivp,int iflag)
{
    register int i,j,k,l,ll,kk;
    int err,cnt1,ip1,ip2,r,nn,n0,n1,n2,tflg,narg,utyp,atyp,aval;
    int gn,hn,nhp,jp,htyp,esidx,mrow,mcol,row,col,j1;
    double tmp,tmp1,tmp2,tmp3,ia,ib,iab,ival,iaval,ibval,buf[E1SSIZ];
    double dval[7],*mptr;
    char cc;
          
/*  static double l1val = 0.0;  used for GRD */
    mrow = mcol = 1;
    j = ivp;

    if (deriv) {
        FNGradN = FNHessN = ivp;
        FNEVATyp = 0;
    }
    cnt1 = cnt + ip;
    /****
    printf1("n=%d ip=%d cnt=%d cnt1=%d \n",n,ip,cnt,cnt1);
    for (i = ip; i < cnt1; ++i) {
        printf1("i=%d typ=%d val=%lf idx=%d\n",i,typ[i],val[i],idx[i]);
    }
    ****/

    for (i = ip; i < cnt1; ++i) {
       
        if (j + 3 + iflag >= E1SSIZ)
            return(97);

        utyp = typ[i];   
        atyp = iabs(utyp);
        esidx = idx[i];
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


        if (atyp >= IFTYP) {        /* special for if then else */

            k = iabs((int)val[i]);

            if (k == 1) {           /* A-part */

                if (deriv > 0) {
                    FNGradN--;
                    if (iflag)
                        FNGradN--;

                    if (deriv > 1)
                        FNHessN--;
                }
                if (iflag == 0) {

                    if (fabs(buf[--j]) <= EPSI2) {
                        while (++i < cnt1) {
                            if (typ[i] == atyp && iabs((int)val[i]) == 2)  
                                break;
                        }
                        if (i >= cnt1)
                            gerr_exit(57);
                    }
                }
                else {
                    kk = i;
                    while (++kk < cnt1) {
                        if (typ[kk] == atyp && iabs((int)val[kk]) == 2)  
                            break;
                    }
                    if (kk >= cnt1)
                        gerr_exit(57);

                    tmp1 = buf[--j];
                    tmp2 = buf[--j];
                    if (tmp1 == 0.0 && tmp2 == 0.0)          
                        i = kk;
                    else if (tmp1 != 1 || tmp2 != 1.0)       
                        val[kk] *= -1.0;
                }
            }
            else if (k == 2) {              /* B-part */
                kk = i;
                while (++kk < cnt1) {
                    if (typ[kk] == atyp && iabs((int)val[kk] == 3))
                        break;
                }
                if (kk >= cnt1)
                    gerr_exit(58);

                if (iflag == 0 || (int)val[i] == 2)   
                    i = kk;
                else {   
                    val[i] = fabs(val[i]);
                    val[kk] *= -1.0;
                }
            }
            else if (k == 3 && iflag && (int)val[i] == -3) {
                val[i] = fabs(val[i]);
                if (buf[j - 4] > buf[j - 2])
                    buf[j - 4] = buf[j - 2];
                if (buf[j - 3] < buf[j - 1])
                    buf[j - 3] = buf[j - 1];
                j -= 2;

                if (deriv > 0) {
                    for (kk = 0; kk < FNGradL; ++kk) {
                        if (FNGrad[FNGradN - 4][kk] > FNGrad[FNGradN - 2][kk])
                            FNGrad[FNGradN - 4][kk] = FNGrad[FNGradN - 2][kk];
                        if (FNGrad[FNGradN - 3][kk] < FNGrad[FNGradN - 1][kk])
                            FNGrad[FNGradN - 3][kk] = FNGrad[FNGradN - 1][kk];
                    }
                    FNGradN -= 2;
                }
            }
        }

        else if (atyp > COFFS) {   /* C-variables: i = 1,...,VCJMax */
            k = atyp - COFFS;
            if (k > VCJMax)
                return(-1);
            else if (utyp >= 0) {
                if (deriv) {
                    r = deriv0(0,deriv,j,iflag);
                    if (r)
                        return(r);
                }
                buf[j++] = VCJVal[k];
            }
            else {
                if (deriv) {
                    r = deriv0(0,deriv,j,iflag);
                    if (r)
                        return(r);
                }
                buf[j++] = -VCJVal[k];
            }
            if (iflag) {
                buf[j] = buf[j - 1];
                j++;
            }
        }
        else if (atyp >= VOFFS) {  /* V-variables */

            k = atyp - VOFFS;

            if (vflag && VTypA[k] == 2) {
                if (utyp >= 0) {
                    if (deriv) {
                        r = deriv0(0,deriv,j,iflag);
                        if (r)
                            return(r);
                    }
                    buf[j++] = AVVAL[k];
                }
                else {
                    if (deriv) {
                        r = deriv0(0,deriv,j,iflag);
                        if (r)
                            return(r);
                    }
                    buf[j++] = -AVVAL[k];
                }
            }   
            else {
                if (utyp >= 0) {
                    if (deriv) {
                        r = deriv0(0,deriv,j,iflag);
                        if (r)
                            return(r);
                    }
                    buf[j++] =  get_data(k,n);
                }
                else {
                    if (deriv) {
                        r = deriv0(0,deriv,j,iflag);
                        if (r)
                            return(r);
                    }
                    buf[j++] = -get_data(k,n);
                }
                /** printf1("get var k=%d n=%d val=%lg\n",k,n,buf[j-1]); **/
            }
            if (iflag) {
                buf[j] = buf[j - 1];
                j++;
            }  
        }
        else if (atyp >= MOFFS) {  /* matrices */

            k = atyp - MOFFS;
            buf[j++] = (double)k;                 
            if (iflag) {
                buf[j] = buf[j - 1];
                j++;
            }
        }
        else {

            narg = (int)val[i];     /* number of arguments */

            if (atyp < 100) {       /* basic operators */

                switch (atyp) {

                  case 0:                         /* value */
                    if (deriv) {
                        r = deriv0(0,deriv,j,iflag);
                        if (r)
                            return(r);
                    }
                    buf[j++] = val[i];
                    if (iflag) {
                        buf[j] = buf[j - 1];
                        j++;
                    }
                    break;

                  case CCTYP:                     /* special: CCTYP */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    buf[j++] = (double)INTMAX;
                    break;

                  case COLTYP:                    /* special: COLTYP */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    buf[j++] = (double)(INTMAX - 1);
                    break;

                  case 1:                         /* + */
                    if (deriv) {
                        r = deriv2(1,deriv,j,buf,iflag);
                        if (r)
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        buf[j++] = tmp1 + tmp2;
                    }
                    else  
                        j = ia_op(1,j,buf);

                    break;

                  case 2:                         /* - */
                    if (deriv) {
                        r = deriv2(2,deriv,j,buf,iflag);
                        if (r)
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        buf[j++] = tmp1 - tmp2;
                    }
                    else  
                        j = ia_op(2,j,buf);

                    break;

                  case 3:                         /* * */
                    if (deriv) {
                        r = deriv2(3,deriv,j,buf,iflag);
                        if (r)
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        buf[j++] = tmp1 * tmp2;
                    }
                    else  
                        j = ia_op(4,j,buf);

                    break;

                  case 4:                         /* & */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        if (fabs(tmp1 * tmp2) > EPSI2)
                            buf[j++] = 1.0;
                        else
                            buf[j++] = 0.0;            
                    }
                    else  
                        j = ia_op(60,j,buf);
                    break;

                  case 5:                         /* / */
                    if (deriv) {
                        r = deriv2(4,deriv,j,buf,iflag);
                        if (r)
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        if (tmp1 == 0.0)
                            buf[j++] = 0.0;
                        else if (tmp2 == 0.0)  
                            return(2);
                        else
                            buf[j++] = tmp1 / tmp2;
                    }
                    else {
                        j = ia_op(3,j,buf);
                        if (j < 0)
                            return(2);
                    }
                    break;

                  case 6:                         /* ^ */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        r = deriv2(5,deriv,j,buf,iflag);
                        if (r)
                            return(r);
                    }
                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    if (fabs(tmp1) <= EPSI2)
                        buf[j++] = 0.0;
                    else if (tmp1 < 0.0)
                        buf[j++] = pow(tmp1,floor(tmp2));
                    else
                        buf[j++] = pow(tmp1,tmp2);
                    break;

                  case 7:                         /* % */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    if (tmp1 != (int)tmp1 || tmp2 != (int)tmp2)
                        return(3);

                    buf[j++] = (double)((int)tmp1 % (int)tmp2);
                    break;

                  case 8:                         /* | */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp2 = buf[--j];
                        tmp1 = buf[--j];
                        if (fabs(tmp1) > EPSI2 || fabs(tmp2) > EPSI2)
                            buf[j++] = 1.0;
                        else
                            buf[j++] = 0.0;            
                    }
                    else  
                        j = ia_op(61,j,buf);
                    break;

                  default:
                    gerr_exit(4);
                }
            }
            else if (atyp < 200) {  /* mathematical operators */

                switch (atyp) {

                  case 100:                 /* exp */

                    if (deriv) {
                        if ((r = derivf1(4,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp = buf[--j];
                        /*********************************
                        if (tmp < ExpMin || tmp > ExpMax)
                            return(4);
                        *****************************/
                        buf[j++] = rexp(tmp);
                    }
                    else
                        j = ia_op(5,j,buf);
                    break;

                  case 101:                 /* log */

                    if (deriv) {
                        if ((r = derivf1(6,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp = buf[--j];
                        if (tmp <= 0.0)
                            return(5);
                        buf[j++] = rlog(tmp);
                    }
                    else {
                        j = ia_op(6,j,buf);
                        if (j < 0)
                            return(5);
                    }
                    break;

                  case 102:                 /* rnd */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (tmp >= 0.0)
                        buf[j++] = floor(tmp + 0.5);
                    else
                        buf[j++] = -floor(-tmp + 0.5);
                    break;

                  case 103:                 /* abs */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(1,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    buf[j++] = fabs(tmp);
                    break;

                  case 104:                 /* sign */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (tmp > 0.0)
                        buf[j++] = 1.0;
                    else if (tmp < 0.0)
                        buf[j++] = -1.0;
                    else
                        buf[j++] = 0.0;
                    break;

                  case 105:                 /* floor */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    buf[j++] = floor(tmp);
                    break;

                  case 106:                 /* ceil */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    buf[j++] = ceil(tmp);
                    break;

                  case 107:                 /* tr(x,a,b) */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf3(1,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    tmp  = buf[--j];

                    if (tmp1 > tmp2)
                        return(10);
                    else if (tmp < tmp1)
                        buf[j++] = tmp1;
                    else if (tmp > tmp2)
                        buf[j++] = tmp2;
                    else
                        buf[j++] = tmp;
                    break;

                  case 108:                 /* min */
                    if (deriv)
                        return(24);             

                    tmp = buf[--j];
                    if (iflag)
                        tmp1 = buf[--j];

                    for (k = 2; k <= narg; ++k) {
                        if (tmp > buf[--j])
                            tmp = buf[j];
                        if (iflag) {
                            if (tmp1 > buf[--j])
                                tmp1 = buf[j];
                        }   
                    }
                    if (iflag)
                        buf[j++] = tmp1;
                    buf[j++] = tmp;
                    break;

                  case 109:                 /* max */
                    if (deriv)
                        return(24);             

                    tmp = buf[--j];
                    if (iflag)
                        tmp1 = buf[--j];

                    for (k = 2; k <= narg; ++k) {
                        if (tmp < buf[--j])
                            tmp = buf[j];
                        if (iflag) {
                            if (tmp1 < buf[--j])
                                tmp1 = buf[j];
                        }
                    }
                    if (iflag)
                        buf[j++] = tmp1;
                    buf[j++] = tmp;
                    break;

                  case 110:                 /* sqrt */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(3,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (tmp < 0.0)    
                        return(6);              
                    else
                        buf[j++] = sqrt(tmp);
                    break;

                  case 111:                 /* sin */
                    if (deriv) {
                        if ((r = derivf1(7,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0)
                        buf[j - 1] = sin(buf[j - 1]);
                    else
                        j = ia_op(10,j,buf);

                    break;

                  case 112:                 /* cos */
                    if (deriv) {
                        if ((r = derivf1(8,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp = buf[--j];
                        buf[j++] = cos(tmp);
                    }
                    else
                        j = ia_op(11,j,buf);

                    break;

                  case 113:                 /* lgam */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(9,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (tmp <= 0.0)
                        return(7);
                    buf[j++] = loggam(tmp,&err);
                    if (err)
                        return(7);
                    break;

                  case 114:                 /* icg1 */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp1 = buf[--j];
                    tmp  = buf[--j];
                    if (tmp <= 0.0 || tmp1 <= 0.0)
                        return(8);                
                    buf[j++] = icgam(tmp,tmp1,&err);
                    if (err)
                        return(32); 
                    break;

                  case 115:                 /* icb */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    tmp  = buf[--j];
                    if (tmp < 0.0 || tmp > 1.0 || tmp1 <= 0.0 || tmp2 <= 0.0)
                        return(9);
                    buf[j++] = incbeta(tmp,tmp1,tmp2,&err);
                    if (err)
                        return(9);
                    break;

                  case 116:                 /* eexp */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(5,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (tmp < ExpMin || tmp > ExpMax)
                        return(4);
                    buf[j++] = exp(tmp) / (1.0 + exp(tmp));
                    break;

                  case 117:                 /* digam */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp = buf[--j];
                    if (tmp <= 0.0)
                        return(7);
                    buf[j++] = digam(tmp,&err);
                    if (err)
                        return(7);
                    break;

                  case 118:                 /* trigam */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    tmp = buf[--j];
                    if (tmp <= 0.0)
                        return(7);
                    buf[j++] = trigam(tmp,&err);
                    if (err)
                        return(7);

                    break;

                  case 119:                       /* ## iv */

                    if (deriv) {
                        r = deriv2(1,deriv,j,buf,iflag);
                        if (r)
                            return(r);

                        for (l = 0; l < FNGradL; ++l) {
                            FNGrad[FNGradN - 1][l] *= 0.5;
                            if (iflag)
                                FNGrad[FNGradN - 2][l] *= 0.5;
                        }
                    }
                    if (iflag == 0) {
                        buf[j - 2] = (buf[j - 2] + buf[j - 1]) / 2.0;
                        j--;   
                    }
                    else {
                        ia_iv(buf[j-4],buf[j-3],buf[j-2],buf[j-1],&buf[j-4],&buf[j-3]);
                        j -= 2;
                    }
                    break;

                  case 120:                 /* icg */
                    if (iflag)
                        return(29);

                    tmp1 = buf[j - 1];
                    tmp  = buf[j - 2];
                    if (tmp <= 0.0 || tmp1 <= 0.0)
                        return(8);                
                    if (icdgam(tmp,tmp1,dval))
                        return(32);

                    if (deriv) {
                        r = derivf2(1,deriv,j,val,dval);
                        if (r)
                            return(r);
                    }
                    buf[j - 2] = dval[6];
                    j--;
                    break;

                  case 121:                 /* bc */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))  
                            return(r);
                    }
                    tmp1 = buf[j - 1];
                    tmp  = buf[j - 2];
                    tmp = bincoeff((int)tmp,(int)tmp1);
                    if (tmp < 0.0)   
                        return(33);                
 
                    buf[j - 2] = tmp;        
                    j--;
                    break;

                  case 122:                 /* bivn */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf3(2,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp1 = buf[j - 1];
                    tmp2 = buf[j - 2];
                    tmp3 = buf[j - 3];
                    buf[j - 3] = bivn(tmp3,tmp2,tmp1,&nn);
                    if (nn)
                        return(34);                
                    j -= 2;
                    break;

                  case 123:                 /* mvn(A,k,eps,x1,...,xn) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    n0 = (int)buf[j - narg];
                    n1 = MatRow[n0];
                    n2 = MatCol[n0];

                    if (n1 != n2)               /* need square matrix */
                        return(37);

                    if (n1 > DMVMax)
                        return(38);

                    nn = narg - 3;
                    if (n1 != nn)
                        return(39);

                    ll = j - narg + 1;
                    kk = (int)buf[ll++];
                    tmp1 = buf[ll++];

                    for (l = 1; l <= nn; ++l)  
                        DMVArg[l] = buf[ll++];

                    mptr = MatVal[n0];

                    /********************
                    printf1("kk=%d ", kk);
                    for (l = 1; l <= nn; ++l)
                        printf1("%lg ",DMVArg[l]);
                    printf1("\nmptr:\n");
                    for (l = 1; l <= nn; ++l) {
                        for (ll = 1; ll <= nn; ++ll)
                            printf1("%lg ",mptr[(l - 1) * nn + ll]);
                        newline();
                    }
                    ******************/

                    r = dmv(nn,kk,DMVArg,mptr,&tmp,tmp1,&tmp2);
                    /*** 
                    printf1("r=%d tmp=%lg tmp1=%lg tmp2=%lg\n",r,tmp,tmp1,tmp2);
                    ***/
                    if (r) {
                        if (r == 1)
                            return(40);
                        else
                            return(41);
                    }
                    j -= narg;
                    buf[j++] = tmp;
                    break;

                  case 124:                 /* poisson(theta,k) */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = deriv2(6,deriv,j,buf,iflag)))  
                            return(r);
                    }
                    else {
                        err = l_poisson(buf[j - 2],buf[j - 1],&tmp,&tmp1,&tmp2);
                        if (err)
                            return(42);
                        buf[j - 2] = tmp;
                    }
                    j--;
                    break;

                  case 125:                 /* negbin(alpha,gamma,k) */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = deriv3(1,deriv,j,buf,iflag)))  
                            return(r);
                    }
                    else {
                        err = l_negbin(buf[j-3],buf[j-2],buf[j-1],0,dval);
                        if (err)
                            return(43);
                        buf[j - 3] = dval[0];
                    }
                    j -= 2;
                    break;

                  case 126:                 /* gclen(lon1,lat1,lon2,lat2) */
                    if (iflag)
                        return(29);
                    if (deriv)  
                        return(24);

                    tmp3 = buf[--j];
                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    tmp  = buf[--j];
                    buf[j++] = geod_distance(tmp,tmp1,tmp2,tmp3) * 180.0 / Pi;
                    break;

/* ## */                                    /* numerical integration */
                  case 185:                 /* intn(s,...) */
                  case 188:                 /* inth(s,...) */
                  case 192:                 /* int(a,b,...) */
     
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if (j != FNGradN || (deriv > 1 && FNGradN != FNHessN))
                            gerr_exit(89);
                        if (j >= MaxSTD - 3)
                            return(92);
                    }
                    ip1 = ip2 = i + 1;
                    tflg = 0;
                    for (l = ip; l < cnt1; ++l) {
                        jp = iabs(typ[l]);             
                        if (jp == 1969) {
                            if (l <= i)  
                                return(31);
                            tflg = 1;
                        }
                        if (l > i && (jp == 186 || jp == 189 || jp == 193)) {
                            ip2 = l;
                            break;
                        }
                    }
                    if (atyp == 192) {
                        htyp = nhp = 0;
                        jp = 2;
                        if (ip2 == ip1 || j < jp)
                            gerr_exit(88);

                        ib = buf[j - 1];
                        ia = buf[j - 2];
                        iab = ib - ia;
                        if (iab <= 0.0)
                            return(27);
                    }
                    else {
                        jp = 1;
                        if (ip2 == ip1 || j < jp)
                            gerr_exit(88);

                        nhp = (int)buf[j - 1];
                        if (nhp < 2 || nhp > 12 || nhp == 8 || nhp == 10 || nhp == 11) 
                            return(30);

                        if (atyp == 188)
                            htyp = 0;
                        else
                            htyp = 1;
                    }
                    if (typ[ip2] < 0)
                        utyp = -utyp;

                    gn = FNGradN;           /* save pointer */
                    hn = FNHessN;

                    if (tflg == 0 && atyp == 192) {        /* constant */
                 
                        r = v_eval1(n,ip2 - ip1,typ,val,idx,&tmp,0,deriv,ip1,j,0);
                        if (r)
                            return(r);

                        ival = tmp * iab;           
                        iaval = ibval = tmp;

                        if (deriv) {
                            for (l = 0; l < FNGradL; ++l) {
                                if (deriv > 1)  
                                    FNGrad[MaxSTD - 1][l] = 
                                    FNGrad[MaxSTD - 2][l] = FNGrad[j][l];
                                FNGrad[j][l] *= iab;
                            }
                            if (deriv > 1) {
                                for (l = 0; l < FNHessL; ++l)  
                                    FNHess[j][l] *= iab;
                            }
                        }
                    }
                    else {  

                        if (deriv && atyp == 192) { 

                                       /* calculate function for t=ia, t=ib */

                            NINTTVAR = ia;
                            r = v_eval1(n,ip2 - ip1,typ,val,idx,&iaval,0,deriv,ip1,j,0);
                            if (r)
                                return(r);

                            if (deriv > 1) {
                                for (l = 0; l < FNGradL; ++l)  
                                    FNGrad[MaxSTD - 2][l] = FNGrad[j][l];
                            }
                            NINTTVAR = ib;
                            r = v_eval1(n,ip2 - ip1,typ,val,idx,&ibval,0,deriv,ip1,j,0);
                            if (r)
                                return(r);

                            if (deriv > 1) {
                                for (l = 0; l < FNGradL; ++l)  
                                    FNGrad[MaxSTD - 1][l] = FNGrad[j][l];
                            }
                        }

                        /* evaluate integral */

                        if (atyp == 192)
                            NINTMUsed = NINTMETH;

                        /* copy parser stack to be used by ni_gen() */
    
                        r = 0;
                        for (l = ip1; l < ip2; ++l) {
                            ESTyp[r] = typ[l];
                            ESVal[r++] = val[l];
                        }
                        ESCnt = r;
                        ESCase = n;
                        ESDeriv = 0;
                        ESIVP = j + 1;

                        ival = ni_gen(ia,ib,1,nhp,htyp,&r);
                        if (r)  
                            return(28);

                        if (deriv) {
                            ESDeriv = 1;
                            for (l = 0; l < FNGradL; ++l) {
                                ESDerivI = l;
                                tmp = ni_gen(ia,ib,1,nhp,htyp,&r);
                                if (r)  
                                    return(28);
                                FNGrad[j][l] = tmp;
                            }
                            if (deriv > 1) {
                                ESDeriv = 2;
                                for (l = 0; l < FNHessL; ++l) {
                                    ESDerivI = l;
                                    tmp = ni_gen(ia,ib,1,nhp,htyp,&r);
                                    if (r)  
                                        return(28);
                                    FNHess[j][l] = tmp;
                                }
                            }
                        }
                    }
                    FNGradN = gn;      /* restore pointer */
                    FNHessN = hn;

                    buf[j - jp] = ival;

                    if (deriv) {

                        if (deriv > 1) {
                            kk = 0;
                            for (l = 0; l < FNGradL; ++l) {
                                for (ll = 0; ll <= l; ++ll) {
                                    tmp = 0.0;
                                    if (atyp == 192) {
                                        if (FNATyp[j - 1] == 0) {
                                            if (FNHess[j - 1][kk])
                                                tmp += ibval * FNHess[j - 1][kk];
                                            if (FNGrad[j - 1][l])
                                                tmp += FNGrad[j - 1][l] * FNGrad[MaxSTD - 1][ll];
                                            if (FNGrad[j - 1][ll])
                                                tmp += FNGrad[j - 1][ll] * FNGrad[MaxSTD - 1][l];
                                        }
                                        if (FNATyp[j - 2] == 0) {
                                            if (FNHess[j - 2][kk])
                                                tmp -= iaval * FNHess[j - 2][kk];
                                            if (FNGrad[j - 2][l])
                                                tmp -= FNGrad[j - 2][l] * FNGrad[MaxSTD - 2][ll];
                                            if (FNGrad[j - 2][ll])
                                                tmp -= FNGrad[j - 2][ll] * FNGrad[MaxSTD - 2][l];
                                        }
                                    }
                                    FNHess[j - jp][kk] = tmp + FNHess[j][kk];
                                    kk++;
                                }
                            }
                            if (atyp == 192)
                                FNHessN--;             
                        }
                        for (l = 0; l < FNGradL; ++l) {

                            tmp = FNGrad[j][l];   
                            if (atyp == 192) 
                                    tmp += ibval * FNGrad[j - 1][l] -
                                           iaval * FNGrad[j - 2][l];  

                            FNGrad[j - jp][l] = tmp;             

                        }
                        if (atyp == 192) {
                            FNGradN--;             
                            FNATyp[j - 2] *= FNATyp[j - 1] * FNATyp[j];
                        }
                        else   
                            FNATyp[j - 1] = 0;                     
                    }
                    if (atyp == 192)
                        j--;            
                    i = ip2;
                    break;

                  case 184:
                  case 187:
                  case 190:                 /* begin of integration */
                  case 191:              

                    break;

                  case 196:                 /* X(.,.) */
                  case 197:                 /* X(.,j) */
                  case 198:                 /* X(i,.) */
                  case 199:                 /* X(i,j) */
                    if (iflag)
                        return(29);
                    if (deriv)  
                        return(24);
     
                    n1 = n2 = -1;           /* n1 = row, n2 = col */
                    if (atyp == 196)
                        j1 = j;
                    else if (atyp == 197) {
                        n2 = (int)buf[j - 1] - 1;                    
                        j1 = j - 1;
                    }
                    else if (atyp == 198) {
                        n1 = (int)buf[j - 1] - 1;                    
                        j1 = j - 1;
                    }
                    else {
                        n1 = (int)buf[j - 2] - 1;   
                        n2 = (int)buf[j - 1] - 1;   
                        j1 = j - 2;
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
                            for (ll = 1; ll <= col; ++ll)  
                                tmp += MatVal[n0][n1 * col + ll];
                        }
                        else if (n1 < 0 && n2 >= 0) {  
                            for (l = 0; l < row; ++l)  
                                tmp += MatVal[n0][l * col + n2 + 1];
                        }
                        else {  
                            for (l = 0; l < row; ++l) {
                                for (ll = 1; ll <= col; ++ll)
                                    tmp += MatVal[n0][l * col + ll];
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
                            for (ll = 0; ll < col; ++ll)  
                                tmp += get_data(NLVIdx[n0][ll],n1);
                        }
                        else if (n1 < 0 && n2 >= 0) {  
                            for (l = 0; l < NOC; ++l)  
                                tmp += get_data(NLVIdx[n0][n2],l);
                        }
                        else {  
                            for (l = 0; l < NOC; ++l) {
                                for (ll = 0; ll < col; ++ll)
                                    tmp += get_data(NLVIdx[n0][ll],l);
                            }
                        }
                    }
                    if (esidx < 0)
                        tmp = -tmp;

                    buf[j1 - 1] = tmp;     
                    j = j1;
                    break;

                  default:
                    gerr_exit(5);
                }
            }
            else if (atyp < 300) {  /* logical operators */

                switch (atyp) {
    
                  case 200:                 /* not */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    if (fabs(tmp) > EPSI2)
                        buf[j++] = 0.0;       
                    else
                        buf[j++] = 1.0;       

                    break;

                  case 201:                 /* and */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (fabs(tmp) > EPSI2 && fabs(tmp1) > EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(60,j,buf);
                    break;

                  case 202:                 /* or */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (fabs(tmp) > EPSI2 || fabs(tmp1) > EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(61,j,buf);
                    break;

                  case 203:                 /* eq */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (fabs(tmp - tmp1) <= EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(57,j,buf);

                    break;

                  case 204:                 /* ne */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (fabs(tmp - tmp1) > EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(58,j,buf);
                    break;

                  case 205:                 /* lt */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (tmp + EPSI2 < tmp1)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(54,j,buf);

                    break;

                  case 206:                 /* le */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (tmp <= tmp1 + EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(53,j,buf);
                    break;

                  case 207:                 /* gt */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (tmp > tmp1 + EPSI2)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(56,j,buf);
                    break;

                  case 208:                 /* ge */
                    if (deriv) {
                        if ((r = deriv2(0,deriv,j,buf,iflag)))
                            return(r);
                    }
                    if (iflag == 0) {
                        tmp1 = buf[--j];
                        tmp  = buf[--j];
                        if (tmp + EPSI2 >= tmp1)
                            buf[j++] = 1.0;       
                        else
                            buf[j++] = 0.0;       
                    }
                    else  
                        j = ia_op(55,j,buf);
                    break;

                  case 209:                 /* if */
                    break;

                  default:
                    gerr_exit(6);
                }
            }
            else if (atyp < 400) {  /* density and distribution functions */

                switch (atyp) {

                  case 300:                 /* ndf */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(10,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    buf[j++] = dnf(tmp);  
                    break;

                  case 301:                 /* nd */
                    if (iflag)
                        return(29);

                    if (deriv) {
                        if ((r = derivf1(11,deriv,j,buf,iflag)))
                            return(r);
                    }
                    tmp = buf[--j];
                    buf[j++] = cdnf(tmp);  
                    break;

                  case 302:                 /* ndi */
                    if (iflag)
                        return(29);

                    tmp = buf[--j];
                    if (tmp <= 0.0 || tmp >= 1.0)
                        return(11);
                    buf[j++] = cdnif1(tmp);  
                    if (deriv) {
                        if ((r = derivf1(12,deriv,j,buf,iflag)))
                            return(r);
                    }
                    break;

                  case 303:                 /* mr */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp = buf[--j];
                    buf[j++] = cdnm(tmp);  
                    break;

                  case 304:                 /* td */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp1 = buf[--j];
                    tmp = buf[--j];
                    if ((int)tmp1 < 1)
                        return(12);
                    buf[j++] = cdtf(tmp,(int)tmp1);
                    break;

                  case 305:                 /* chd */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp1 = buf[--j];
                    tmp = buf[--j];
                    if (tmp < 0.0 || (int)tmp1 < 1)
                        return(13);
                    buf[j++] = cdchif(tmp,(int)tmp1);
                    break;

                  case 306:                 /* fd */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp2 = buf[--j];
                    tmp1 = buf[--j];
                    tmp = buf[--j];
                    if (tmp < 0.0 || (int)tmp1 < 1 || (int)tmp2 < 1)
                        return(14);
                    buf[j++] = cdff(tmp,(int)tmp1,(int)tmp2);
                    break;

                  default:
                    gerr_exit(7);
                }
            }
            else if (atyp < OPT2S) {  /* some other operators */

                switch (atyp) {

                  case 400:                 /* dummy variable operators */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    j--;
                    j -= narg;
                    tmp = buf[j];
                    buf[j] = 0.0;
                    
                    tmp1 = 1.0;
                    k = 1;
                    while (k <= narg) {
                        tmp2 = buf[j + k];
                        if (fabs(tmp - tmp2) < EPSI2) {
                            buf[j] = tmp1;
                            break;
                        }
                        if (++k > narg)
                            break;

                        tmp3 = buf[j + k];
                        if ((int)tmp3 == INTMAX - 1) {
                            tmp1 = -1.0;
                            k++;
                        }
                        else if ((int)tmp3 == INTMAX) {
                            k++;
                            if (tmp2 <= tmp && tmp <= buf[j + k]) {
                                buf[j] = tmp1;
                                break;
                            }
                        }
                    }
                    j++;
                    break;

                  case 402:                 /* jul */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    n2 = (int)buf[--j];        /* day */
                    n1 = (int)buf[--j];        /* month */
                    n0 = (int)buf[--j];        /* year */

                    k  = n2 - 32075 + 1461 *
                         (n0 + 4800 + (n1 - 14) / 12) / 4
                             + 367 * (n1 - 2 - (n1 - 14) / 12 * 12) / 12 
                             - 3 * ((n0 + 4900 + (n1 - 14) / 12) / 100) / 4;

                    buf[j++] = (double)k;
                    break;

                  case 403:                 /* jul1 (ACM 199) */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    n2 = (int)buf[--j];        /* day */
                    n1 = (int)buf[--j];        /* month */
                    n0 = (int)buf[--j];        /* year */

                    if (n1 > 2)
                        n1 -= 3;
                    else {
                        n1 += 9;
                        n0--; 
                    }
                    l = n0 / 100;
                    nn = n0 - 100 * l;
                    k = (146097 * l) / 4 + (1461 * nn) / 4 +
                        (153 * n1 + 2) / 5 + n2 + 1721119;

                    buf[j++] = (double)k;
                    break;

                  case 404:                 /* rd(a,b) */
                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             
                    tmp1 = buf[--j];
                    tmp  = buf[--j];
                    buf[j++] = tmp + random1() * (tmp1 - tmp);
                    break;

                  case 405:                 /* rdmn */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    n0 = (int)buf[--j];     /* index to matrix */
                    k  = (int)buf[--j]; 
                                    
                    n1 = MatRow[n0];
                    n2 = MatCol[n0];
                    if (n1 != n2)
                        return(35);
   
                    mptr = MatVal[n0];      /* pointer to matrix */

                    if (RDMNInit == 0) {            /* init */
                        if (rdmn_init(n1,mptr))         
                            return(35);
                    }
                    if (k < 1 || k > RDMNN)
                        return(36);

                    if (k == 1)     /* get next random vector */
                        rdmn(); 
                    buf[j++] = RDMNRd[k];
                    break;

                  case 406:             /* grd(V) */
                    if (iflag)
                        return(29);
           
                    if (deriv)
                        return(24);             
                    tmp = buf[--j];     /* current value of variable V  */
                    if (n == 0) {
                        buf[j++] = val[i] = random1();
                        val[i - 1] = tmp;
                    }
                    else {
                        if (fabs(tmp - val[i - 1]) > EPSI2) {
                            buf[j++] = val[i] = random1();
                            val[i - 1] = tmp;
                        }
                        else  
                            buf[j++] = val[i];
                    }
                    break;

                  case 409:                 /* july */
                  case 410:                 /* julm */
                  case 411:                 /* juld */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    k = (int) buf[--j];     /* julian time */
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
                        buf[j++] = (double)n0;
                    else if (atyp == 410)
                        buf[j++] = (double)n1;
                    else
                        buf[j++] = (double)n2;

                    break;   

                  case 412:             /* num(a,i) */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    tmp = buf[--j];    
                    tmp1 = buf[--j];
                    if (n == 0)  
                        buf[j] = tmp1;
                    else
                        buf[j] = val[i] + tmp;
                    val[i] = buf[j];
                    j++;
                    break;

                  case 413:             /* recode */ 
                  case 414:             /* recode1 */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    if (atyp == 414)
                        tmp1 = buf[--j];

                    n0 = (int)buf[--j];
                    tmp = buf[--j];    
                                    
                    n1 = MatRow[n0];
                    n2 = MatCol[n0];
                    mptr = MatVal[n0] + 1;
           
                    for (l = 0; l < n1; ++l) {
                        if (tmp == *mptr) {
                            buf[j] = *++mptr;
                            break;
                        }
                        mptr += n2;
                    }
                    if (atyp == 414 && l >= n1)
                        buf[j] = tmp1;
                    j++;
                    break;

                  case 415:                 /* row(X) */
                  case 416:                 /* col(X) */
                  case 417:                 /* begin of row/col expression */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             
                    if (atyp == 415)
                        buf[j - 1] = (double)mrow;
                    else if (atyp == 416)
                        buf[j - 1] = (double)mcol;
                    mrow = mcol = 1;
                    break;

                  case 418:                 /* exists(A) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             
                    if (val[i] < 0.0)
                        buf[j++] = 0.0;
                    else
                        buf[j++] = 1.0;
                    break;

                  case 419:                 /* rdp1(x) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             
                    k = rdp1(buf[j - 1]);
                    if (k < 0)
                        return(46);
                    buf[j - 1] = (double)k;
                    break;

                  case 450:                 /* str(n,m) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    j -= 2;
                    buf[j++] = 0.0;
                    break;

                  case 451:                 /* strlen(S) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    k = aval - VOFFS;
                    if (k < 0 || k >= NVAR)
                        return(-1);
                    if (VTyp[k] != 1)
                        return(47);
                    buf[j - 1] = (double)(-VSLen[k]);
                    break;

                  case 452:                 /* strvp(S,n,m) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    k = aval - VOFFS;
                    if (k < 0 || k >= NVAR)
                        return(-1);
                    if (VTyp[k] != 1)
                        return(47);
            
                    n1 = (int)buf[j - 2];
                    n2 = (int)buf[j - 1];  
                    if (n1 < 1 || n2 < n1)
                        return(48);

                    n1--;
                    if (n2 > SVBufLen)
                        n2 = SVBufLen;

                    get_str(SVBuf,k,n);

                    cc = *(SVBuf + n2);
                    *(SVBuf + n2) = '\0';

                    if (sscanf(SVBuf + n1,"%lf",&tmp) != 1)
                        tmp = MBlnkVal;

                    *(SVBuf + n2) = cc;
                    j -= 2;
                    buf[j - 1] = tmp;                    
                    break;


                  case 453:                 /* strv(S) */
                    if (iflag)
                        return(29);
                    if (deriv)
                        return(24);             

                    k = aval - VOFFS;
                    if (k < 0 || k >= NVAR)
                        return(-1);
                    if (VTyp[k] != 1)
                        return(47);

                    get_str(SVBuf,k,n);
                    if (sscanf(SVBuf,"%lf",&tmp) != 1)
                        tmp = MBlnkVal;
                    buf[j - 1] = tmp;                    
                    break;


                  /* operators for sequences */

                  case 501:             /* slen(Y1,,Y2) */
                  case 502:             /* glen(Y1,,Y2) */

                    if (iflag)
                        return(29);

                    if (deriv)
                        return(24);             

                    n0 = 0;
                    n1 = n2 = -1;
                    for (k = narg; k >= 1; --k) {
                        if (buf[j - k] >= 0) {
                            n1 = k;
                            break;
                        }
                    }
                    for (k = 1; k <= narg; ++k) {
                        if (buf[j - k] >= 0) {
                            n2 = k;
                            break;
                        }
                    }
                    if (n1 >= 1) {
                        if (atyp == 501)
                            n0 = n1 - n2 + 1;
                        else {
                            for (k = n2 + 1; k < n1; ++k) {
                                if (buf[j - k] < 0.0)
                                    n0++;
                            }
                        }
                    }
                    j -= narg;
                    buf[j++] = (double)n0;
                    break;

                  default:
                    gerr_exit(8);
                }
            }
            else if (atyp < OPT2A) {  /* ### spatial operators */

                if (iflag)
                    return(29);

                if (deriv)
                    return(24);             

                if (SDVarDef == 0)
                    return(49);

                switch (atyp) {
                  case 601:                 /* sdp(x,y) */

                    tmp2 = buf[--j];
                    tmp1 = buf[--j];

                    if ((int)get_data(SDVarSDTyp,n) != 3)
                        buf[j++] = 0.0;
                    else {
                        if ((nn = sd_getdata(n,0,1,1)) < 1)
                            return(50);

                        buf[j++] = (double)g_inpoly(nn - 1,SDVarX,SDVarY,tmp1,tmp2);
                    }
                    break;

                  default:
                    gerr_exit(8);
                }
            }
            else if (atyp < OPT2B) {  /* special type 2 operators */

                if (iflag)
                    return(29);

                if (deriv)
                    return(24);             

                switch (atyp) {
     
                    case 1510:          /* pre */

                        if (n <= 0 || n == E2FNum)
                            buf[j++] = 0.0;
                        else  
                            buf[j++] =  get_data((int)val[i],n - 1);
                        break;

                    case 1509:          /* suc */

                        if (n >= E2LNum - 1)
                            buf[j++] = 0.0;
                        else  
                            buf[j++] =  get_data((int)val[i],n + 1);
                        break;

                    /********************************
                    case 1890:             bnrec   

                        buf[j++] = (double)(E2LNum - E2FNum);
                        break;
            
                    case 1891:             brec   

                        buf[j++] = (double)(E2FNum + n);
                        break;
                    *****************************************/
                    default:
                        gerr_exit(65);
                }
            }
            else if (atyp < IAOFFS) {  /* constants and operators without
                                          arguments */

                if (iflag)
                    return(29);

                if (deriv) {
                    r = deriv0(0,deriv,j,iflag);
                    if (r)
                        return(r);
                }
   
                switch (atyp) {

                  case 1900:            /* case */
   
                    buf[j++] = (double)(GDNRec + n + 1);
                    break;
            
                  case 1901:            /* nvar */
                    buf[j++] = (double)NVAR;               
                    break;
            
                  case 1902:            /* pi */
                    buf[j++] = Pi;
                    break;
            
                  case 1903:            /* rd */
                    buf[j++] = random1();
                    break;
            
                  case 1904:            /* rdn */
                    buf[j++] = normal();
                    break;

                  case 1905:            /* rdn1 */
                    buf[j++] = normal1();
                    break;

                  case 1906:            /* nocdm */
                    buf[j++] = (double)NOCDM;               
                    break;
            
                  case 1907:            /* noc */
                    buf[j++] = (double)NOC;               
                    break;

                  case 1908:            /* bnoc */
                    buf[j++] = (double)BNOC;               
                    break;

                  case 1910:            /* sdxmin */
                  case 1911:            /* sdxmax */
                  case 1912:            /* sdymin */
                  case 1913:            /* sdymax */
                  case 1914:            /* sdarea */
                  case 1915:            /* sdlen */

                    if (SDVarDef == 0)
                        return(49);
                    if ((nn = sd_getdata(n,0,1,1)) < 1)
                        return(50);
        
                    if (atyp == 1910)
                        buf[j++] = g_pol_xmin(nn - 1,SDVarX,SDVarY);
                    else if (atyp == 1911)
                        buf[j++] = g_pol_xmax(nn - 1,SDVarX,SDVarY);
                    else if (atyp == 1912)
                        buf[j++] = g_pol_ymin(nn - 1,SDVarX,SDVarY);
                    else if (atyp == 1913)
                        buf[j++] = g_pol_ymax(nn - 1,SDVarX,SDVarY);
                    else if (atyp == 1914) {
                        if ((int)get_data(SDVarSDTyp,n) != 3)
                            buf[j++] = 0.0;
                        else
                            buf[j++] = g_pol_area(nn - 1,SDVarX,SDVarY);
                    }
                    else if (atyp == 1915) {
                        l = (int)get_data(SDVarSDTyp,n);        
                        if (l == 2)
                            buf[j++] = g_line_len(nn - 1,SDVarX,SDVarY);
                        else if (l == 3)
                            buf[j++] = g_line_len(nn,SDVarX,SDVarY);
                        else
                            buf[j++] = 0.0;
                    }
                    break;

                  case 1969:            /* t */
                    buf[j++] = NINTTVAR;
                    break; 

                  case 1970:            /* time */
                    buf[j++] = EDVALTime;
                    break;

                  case 1980:            /* sn  */
                    buf[j++] = (double)EDVALSn;
                    break;

                  case 1981:            /* org */
                    buf[j++] = (double)EDVALOrg;
                    break;

                  case 1982:            /* des */
                    buf[j++] = (double)EDVALDes;
                    break;

                  case 1983:            /* ts  */
                    buf[j++] = (double)EDVALTs;
                    break;

                  case 1984:            /* tf  */
                    buf[j++] = (double)EDVALTf;
                    break;

                  default:
                    gerr_exit(9);
                }
                if (iflag) {
                    buf[j] = buf[j - 1];
                    j++;
                }
            }
            else if (atyp < FAOFFS) {  /* intermediate function arguments */

                k = atyp - IAOFFS;
                if (k < 0 || k >= FNPN)
                    return(17);
                /****
                printf1("put intermediate %lg\n",FNPVal[k]);
                ****/

                buf[j++] = FNPVal[k];                   
                if (iflag)
                    buf[j++] = FNPVal1[k];                   

                if (deriv > 0) {            /* copy gradient/ hessian */
       
                    if (FNGradN + iflag >= MaxSTD)
                        return(92);

                    FNATyp[FNGradN] = FNPEATyp[k];

                    for (l = 0; l < FNGradL; ++l)
                        FNGrad[FNGradN][l] = FNPGrad[k][l];
                    FNGradN++;
                         
                    if (iflag) {
                        FNATyp[FNGradN] = FNPEATyp[k];
                        for (l = 0; l < FNGradL; ++l)
                            FNGrad[FNGradN][l] = FNPGrad1[k][l];
                        FNGradN++;
                    }    
                    if (deriv > 1) { 
                        if (FNHessN + iflag >= MaxSTD)
                            return(92);
                        for (l = 0; l < FNHessL; ++l)
                            FNHess[FNHessN][l] = FNPHess[k][l];
                        FNHessN++;
                    }
                }
            }
            else if (atyp < FAOFFS + MaxP) {  /* function arguments */

                k = atyp - FAOFFS;
                if (k < 0 || k >= FNArgN)
                    return(18);

                if (deriv) {
                    if ((r = deriv1(deriv,FNArgSPI[k],j,iflag)))           
                        return(r);
                }
                buf[j++] = FNArgVal[FNArgSPI[k]];                   
 
                /***********
                printf1("k=%d\n",k);
                printf1("##   put argument k===%d    for j=%d %lg\n",k, j-1,buf[j-1]);
                for (k = 0; k < FNArgN; ++k)
                printf1("k=%d fnargsp=%d spi=%d argval=%lg perm=%lg\n",k,FNArgSP[k],
                            FNArgSPI[k], FNArgVal[k],FNArgVal[FNArgSPI[k]]);
                *********************/       

                if (iflag)
                    buf[j++] = FNArgVal1[FNArgSPI[k]];                   
            }
            else            /* undefined operator */
                return(1);

            if (utyp < 0) {        /* change sign */

                if (deriv) {
                    if ((r = deriv0(1,deriv,j,iflag)))
                        return(r);
                }
                if (iflag == 0)
                    buf[j - 1] *= -1.0;
                else {
                    tmp = buf[j - 2];
                    buf[j - 2] = buf[j - 1] * -1.0;
                    buf[j - 1] = -tmp;
                }
            }
        }
    } 
    if (j != ivp + 1 + iflag) {         
        /****
        printfe("Fatal syntax error in expression (j=%d ivp=%d iflag=%d).\n",j,ivp,iflag);
        printfe("Cannot continue.\n");
        gerr_exit(204);         
        ****/
        return(90);
    }
    if (iflag) {
        if (buf[ivp] > buf[ivp + 1]) {
            printfe("j=%d ivp=%d buf=%lg buf=%lg\n",j,ivp,buf[ivp],buf[ivp + 1]);
            gerr_exit(91);
        }
    }
    if (deriv) {
        if (FNGradN != ivp + 1 + iflag || (deriv > 1 && FNHessN != FNGradN)) {
            printfe("\nDERIVATIVES STACK ERROR (%d,%d)\n",FNGradN,FNHessN);
            gerr_exit(205);
        }
        FNEVATyp = FNATyp[ivp];
    }
    res[0] = buf[ivp];
    if (iflag) {
        res[0] = dmin(buf[ivp],buf[ivp + 1]);
        res[1] = dmax(buf[ivp],buf[ivp + 1]);
    }
    return(0);
}
