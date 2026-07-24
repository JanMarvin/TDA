/****************************************************************************/
/*  t_rand                                                                  */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-94 Goetz Rohwer. All rights reserved.           */
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
#include "t_lin.h"
#include "tda_context.h"

/*  functions in t_rand.c */

double random1(TDAContext *ctx);
double random2(TDAContext *ctx);
double random3(TDAContext *ctx);
double normal(TDAContext *ctx);
double normal1(TDAContext *ctx);
int rdmn_init(TDAContext *ctx, int n,double *a);
void rdmn_free(TDAContext *ctx);
void rdmn(TDAContext *ctx);
int rdp1(TDAContext *ctx, double lambda);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define RDInit   13421773   /* default initialization value of random1()    */




/* ------------------------------------------------------------------------ */
/*  random1                                                                 */
/*      Generates pseudo random numbers, equally distributed in the         */
/*      interval (a,b).                                                     */
/*                                                                          */
/*      Initialization at first call with RDInit.                           */
/*                                                                          */
/*      Ref.: M.C. Pike, I.D. Hill, Pseudo-Random Numbers, Algorithm 266    */
/*      Communications of the ACM 8 (1965), No. 10, 605 - 606               */
/*                                                                          */
/*      Note: we use the modificated algorithm for an 32 bit integer        */
/*      arithmetic, cf. M.C. Pike, I.D. Hill, Remark on Algorithm 266,      */
/*      Communications of the ACM 1965, 687.                                */
/*                                                                          */

double random1(TDAContext *ctx)
{
    static TLONG y;

    if (ctx->RD1I) {
        y = ctx->RD1Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        ctx->RD1I = 0;
    }
    ctx->RD1Gen++;    /* count random numbers */
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 5;
    y -= (y / 67108864) * 67108864;
    return((double) y / 67108864.0);
}

/* ------------------------------------------------------------------------ */
/*  random2                                                                 */
/*      Generates pseudo random numbers, equally distributed in the         */
/*      interval (a,b). (Copy of random1)                                   */
/*                                                                          */
/*      Initialization at first call with RDInit.                           */
/*                                                                          */
/*      Ref.: M.C. Pike, I.D. Hill, Pseudo-Random Numbers, Algorithm 266    */
/*      Communications of the ACM 8 (1965), No. 10, 605 - 606               */
/*                                                                          */
/*      Note: we use the modificated algorithm for an 32 bit integer        */
/*      arithmetic, cf. M.C. Pike, I.D. Hill, Remark on Algorithm 266,      */
/*      Communications of the ACM 1965, 687.                                */
/*                                                                          */

double random2(TDAContext *ctx)
{
    static TLONG y;

    if (ctx->RD2I) {
        y = ctx->RD2Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        ctx->RD2I = 0;
    }
    ctx->RD2Gen++;    /* count random numbers */
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 5;
    y -= (y / 67108864) * 67108864;
    return((double) y / 67108864.0);
}

/* ------------------------------------------------------------------------ */
/*  random3                                                                 */
/*      Generates pseudo random numbers, equally distributed in the         */
/*      interval (a,b). (Copy of random1)                                   */
/*                                                                          */
/*      Initialization at first call with RDInit.                           */
/*                                                                          */
/*      Ref.: M.C. Pike, I.D. Hill, Pseudo-Random Numbers, Algorithm 266    */
/*      Communications of the ACM 8 (1965), No. 10, 605 - 606               */
/*                                                                          */
/*      Note: we use the modificated algorithm for an 32 bit integer        */
/*      arithmetic, cf. M.C. Pike, I.D. Hill, Remark on Algorithm 266,      */
/*      Communications of the ACM 1965, 687.                                */
/*                                                                          */

double random3(TDAContext *ctx)
{
    static TLONG y;

    if (ctx->RD3I) {
        y = ctx->RD3Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        ctx->RD3I = 0;
    }
    ctx->RD3Gen++;    /* count random numbers */
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 25;
    y -= (y / 67108864) * 67108864;
    y *= 5;
    y -= (y / 67108864) * 67108864;
    return((double) y / 67108864.0);
}

/* ------------------------------------------------------------------------ */
/*  normal                                                                  */
/*      Generate standard normal distributed pseudo random numbers.         */
/*      Each call of the function returns two random numbers.               */
/*                                                                          */
/*      This function uses random1.                                         */
/*                                                                          */
/*      Ref.: J.R. Bell, Normal Random Deviates, Algorithm 334              */
/*      Communications of the ACM 11 (1968), 498                            */
/*                                                                          */

double normal(TDAContext *ctx)
{
    double s,x,y,xx,yy,l;

    if (ctx->s_normal_dflg) {
        ctx->s_normal_dflg = 0;
        return(ctx->s_normal_dev2);
    }
    else {
        s = 2.0;
        while (s > 1) {
            x = random1(ctx);
            y = 2.0 * random1(ctx) - 1.0;
            xx = x * x;
            yy = y * y;
            s = xx + yy;
        }
        l = sqrt(-2.0 * log(random1(ctx))) / s;
        ctx->s_normal_dev1 = (xx - yy) * l;
        ctx->s_normal_dev2 = 2.0 * x * y * l;
        ctx->s_normal_dflg = 1;
        return(ctx->s_normal_dev1);
    }
}

/* ------------------------------------------------------------------------ */
/*  normal1     Generation of standard normally distributed pseudo          */
/*              random numbers.                                             */
/*              Ref.  ACM Algorithm 488.                                    */

static double DD[61] = {
 0.0,
 0.674489750,0.475859630,0.383771164,0.328611323,0.291142827,0.263684322,                                     
 0.242508452,0.225567444,0.211634166,0.199924267,0.189910758,0.181225181,                                     
 0.173601400,0.166841909,0.160796729,0.155349717,0.150409384,0.145902577,                                     
 0.141770033,0.137963174,0.134441762,0.131172150,0.128125965,0.125279090,                                     
 0.122610883,0.120103560,0.117741707,0.115511892,0.113402349,0.111402720,                                     
 0.109503852,0.107697617,0.105976772,0.104334841,0.102766012,                                    
 0.101265052,0.099827234,0.098448282,0.097124309,0.095851778,0.094627461,                                     
 0.093448407,0.092311909,0.091215482,0.090156838,0.089133867,0.088144619,                                     
 0.087187293,0.086260215,0.085361834,0.084490706,0.083645487,0.082824924,                                     
 0.082027847,0.081253162,0.080499844,0.079766932,0.079053527,0.078358781,                                     
 0.077681899                                                              
};

double normal1(TDAContext *ctx)
{
    register int i;
    double a,v,w;

    a = 0.0;                                                                  
    i = 0;                                                                    
    while (1) {
        ctx->s_normal1_u += ctx->s_normal1_u;                                                                   
        if (ctx->s_normal1_u < 1.0)
            break;
        ctx->s_normal1_u -= 1.0;                                                              
        i++;
        a -= DD[i];                                                             
    }
    while (1) {
        w = DD[i + 1] * ctx->s_normal1_u;                                                             
        v = w * (0.5 * w - a);                                                          
NDLab1:
        ctx->s_normal1_u = random1(ctx);                                                             
        if (v <= ctx->s_normal1_u)
            break;                                                        
        v = random1(ctx);                                                            
        if (ctx->s_normal1_u > v)
            goto NDLab1;                                                  
        ctx->s_normal1_u = (v - ctx->s_normal1_u) / (1.0 - ctx->s_normal1_u);                                                        
    }
    ctx->s_normal1_u = (ctx->s_normal1_u - v) / (1.0 - v);                                                        
    ctx->s_normal1_u += ctx->s_normal1_u;                                                                
    if (ctx->s_normal1_u < 1.0) 
        return(a - w);
    else {
        ctx->s_normal1_u -= 1.0;                                                              
        return(w - a);
    }
}

/* ------------------------------------------------------------------------ */
/*  rdmn_init(n,a)  Init rdmn random number generator.                      */
/*                                                                          */
/*                  Expect an (n,n) correlation matrix a(i,j).              */
/*                  Put Cholesky factor into RDMNFac.                       */
/*                  If successful, RDMNInit = 1.                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int rdmn_init(TDAContext *ctx, int n,double *a)
{
    register int i,j,k;
    int err,nn;

    rdmn_free(ctx);     /* free previously allocated storage */

    err = -1;
    ctx->RDMNN = n;
    nn = n * (n + 1) / 2;
    if (!(ctx->RDMNFac = (double *)calloc((size_t)(nn + 1),sizeof(double))))  
        goto RDMNIFin;
    ctx->RDMNFacA = nn + 1;
    memrq(ctx, ctx->RDMNFacA,sizeof(double));

    if (!(ctx->RDMNRd = (double *)calloc((size_t)(n + 1),sizeof(double))))  
        goto RDMNIFin;
    ctx->RDMNRdA = n + 1;
    memrq(ctx, ctx->RDMNRdA,sizeof(double));

    if (!(ctx->RDMNTmp = (double *)calloc((size_t)(n + 1),sizeof(double))))  
        goto RDMNIFin;
    ctx->RDMNTmpA = n + 1;
    memrq(ctx, ctx->RDMNTmpA,sizeof(double));
              
    k = 1;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= i; ++j)
            ctx->RDMNFac[k++] = a[(i - 1) * n + j];
    }
    err = decomp(ctx, n,ctx->RDMNFac);
    if (err == 0)
        ctx->RDMNInit = 1;

RDMNIFin:
    if (err)
        rdmn_free(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rdmn_free()     Free previously allocated storage for rdmn generator.   */

void rdmn_free(TDAContext *ctx)
{
    if (ctx->RDMNFacA > 0) {
        free((char *)ctx->RDMNFac);
        memrq(ctx, -ctx->RDMNFacA,sizeof(double));
        ctx->RDMNFacA = 0;
    }
    if (ctx->RDMNRdA > 0) {
        free((char *)ctx->RDMNRd);
        memrq(ctx, -ctx->RDMNRdA,sizeof(double));
        ctx->RDMNRdA = 0;
    }
    if (ctx->RDMNTmpA > 0) {
        free((char *)ctx->RDMNTmp);
        memrq(ctx, -ctx->RDMNTmpA,sizeof(double));
        ctx->RDMNTmpA = 0;
    }
    ctx->RDMNInit = 0;
}

/* ------------------------------------------------------------------------ */
/*  rdmn()  Generates multivariate normal random numbers.                   */
/*          The result of one draw is stored in RDMNRd[].                   */
  
void rdmn(TDAContext *ctx)
{
    register int i,j,k;
    double tmp;

    k = 1;
    for (i = 1; i <= ctx->RDMNN; ++i) {
        ctx->RDMNTmp[i] = normal1(ctx);
        tmp = 0.0;
        for (j = 1; j <= i; ++j)
            tmp += ctx->RDMNFac[k++] * ctx->RDMNTmp[j];
        ctx->RDMNRd[i] = tmp;
    }
}

/* ------------------------------------------------------------------------ */
/*  rdp1(lambda)        returns random number distributed according to      */
/*                      a Poisson distribution with mean lambda.            */
/*                      Adapted from ACM 369.                               */
/*  Return -1 if lambda <= 0 or exceeds limits of exp().                    */
/*  Uses a local copy of random().                                          */
  
int rdp1(TDAContext *ctx, double lambda)
{
    register int k;
    double z,t;

    if (lambda <= 0.0 || lambda > ExpMax) 
        return(-1);

    z = exp(-lambda);
    k = 0;
    t = 1.0;
    while ((t *= random3(ctx)) > z)  
        k++;
    return(k);
}

