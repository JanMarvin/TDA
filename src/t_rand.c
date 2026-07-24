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

/*  functions in t_rand.c */

double random1(void);
double random2(void);
double random3(void);
double normal(void);
double normal1(void);
int rdmn_init(int n,double *a);
void rdmn_free(void);
void rdmn(void);
int rdp1(double lambda);

/* ------------------------------------------------------------------------ */
/*  global variables                                                        */

#define RDInit   13421773   /* default initialization value of random1()    */

int RD1Skip = 0;            /* initially skipped random numbers in random1  */
int RD2Skip = 0;            /* initially skipped random numbers in random2  */
int RD3Skip = 0;            /* initially skipped random numbers in random3  */
int RD1Gen = 0;             /* generated random numbers, random1()          */
int RD2Gen = 0;             /* generated random numbers, random2()          */
int RD3Gen = 0;             /* generated random numbers, random3()          */
int RD1Seed = RDInit;       /* Seed for random1(), def. RDInit in tda.h     */
int RD2Seed = RDInit;       /* Seed for random2(), def. RDInit in tda.h     */
int RD3Seed = RDInit;       /* Seed for random3(), def. RDInit in tda.h     */

int RD1I = 1;               /* init flag for random1()                      */    
int RD2I = 1;               /* init flag for random2()                      */
int RD3I = 1;               /* init flag for random3()                      */

int RDMNInit = 0;           /* set for initialization                       */
int RDMNN = 0;              /* number of dimensions                         */
double *RDMNFac;            /* Multivariate normal Cholesky factor          */
int RDMNFacA = 0;           /* if allocated                                 */
double *RDMNRd;             /* storage for MVN random numbers               */
int RDMNRdA = 0;            /* if allocated                                 */
double *RDMNTmp;            /* temporary storage for random numbers         */
int RDMNTmpA = 0;           /* if allocated                                 */

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

double random1(void)
{
    static TLONG y;

    if (RD1I) {
        y = RD1Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        RD1I = 0;
    }
    RD1Gen++;    /* count random numbers */
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

double random2(void)
{
    static TLONG y;

    if (RD2I) {
        y = RD2Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        RD2I = 0;
    }
    RD2Gen++;    /* count random numbers */
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

double random3(void)
{
    static TLONG y;

    if (RD3I) {
        y = RD3Seed;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 25;
        y -= (y / 67108864) * 67108864;
        y *= 5;
        y -= (y / 67108864) * 67108864;
        RD3I = 0;
    }
    RD3Gen++;    /* count random numbers */
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

double normal(void)
{
    static int dflg = 0;
    static double dev1,dev2;
    double s,x,y,xx,yy,l;

    if (dflg) {
        dflg = 0;
        return(dev2);
    }
    else {
        s = 2.0;
        while (s > 1) {
            x = random1();
            y = 2.0 * random1() - 1.0;
            xx = x * x;
            yy = y * y;
            s = xx + yy;
        }
        l = sqrt(-2.0 * log(random1())) / s;
        dev1 = (xx - yy) * l;
        dev2 = 2.0 * x * y * l;
        dflg = 1;
        return(dev1);
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

double normal1(void)
{
    register int i;
    double a,v,w;
    static double u = 0.0;

    a = 0.0;                                                                  
    i = 0;                                                                    
    while (1) {
        u += u;                                                                   
        if (u < 1.0)
            break;
        u -= 1.0;                                                              
        i++;
        a -= DD[i];                                                             
    }
    while (1) {
        w = DD[i + 1] * u;                                                             
        v = w * (0.5 * w - a);                                                          
NDLab1:
        u = random1();                                                             
        if (v <= u)
            break;                                                        
        v = random1();                                                            
        if (u > v)
            goto NDLab1;                                                  
        u = (v - u) / (1.0 - u);                                                        
    }
    u = (u - v) / (1.0 - v);                                                        
    u += u;                                                                
    if (u < 1.0) 
        return(a - w);
    else {
        u -= 1.0;                                                              
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

int rdmn_init(int n,double *a)
{
    register int i,j,k;
    int err,nn;

    rdmn_free();     /* free previously allocated storage */

    err = -1;
    RDMNN = n;
    nn = n * (n + 1) / 2;
    if (!(RDMNFac = (double *)calloc(nn + 1,sizeof(double))))  
        goto RDMNIFin;
    RDMNFacA = nn + 1;
    memrq(RDMNFacA,sizeof(double));

    if (!(RDMNRd = (double *)calloc(n + 1,sizeof(double))))  
        goto RDMNIFin;
    RDMNRdA = n + 1;
    memrq(RDMNRdA,sizeof(double));

    if (!(RDMNTmp = (double *)calloc(n + 1,sizeof(double))))  
        goto RDMNIFin;
    RDMNTmpA = n + 1;
    memrq(RDMNTmpA,sizeof(double));
              
    k = 1;
    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= i; ++j)
            RDMNFac[k++] = a[(i - 1) * n + j];
    }
    err = decomp(n,RDMNFac);
    if (err == 0)
        RDMNInit = 1;

RDMNIFin:
    if (err)
        rdmn_free();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  rdmn_free()     Free previously allocated storage for rdmn generator.   */

void rdmn_free(void)
{
    if (RDMNFacA > 0) {
        free((char *)RDMNFac);
        memrq(-RDMNFacA,sizeof(double));
        RDMNFacA = 0;
    }
    if (RDMNRdA > 0) {
        free((char *)RDMNRd);
        memrq(-RDMNRdA,sizeof(double));
        RDMNRdA = 0;
    }
    if (RDMNTmpA > 0) {
        free((char *)RDMNTmp);
        memrq(-RDMNTmpA,sizeof(double));
        RDMNTmpA = 0;
    }
    RDMNInit = 0;
}

/* ------------------------------------------------------------------------ */
/*  rdmn()  Generates multivariate normal random numbers.                   */
/*          The result of one draw is stored in RDMNRd[].                   */
  
void rdmn(void)
{
    register int i,j,k;
    double tmp;

    k = 1;
    for (i = 1; i <= RDMNN; ++i) {
        RDMNTmp[i] = normal1();
        tmp = 0.0;
        for (j = 1; j <= i; ++j)
            tmp += RDMNFac[k++] * RDMNTmp[j];
        RDMNRd[i] = tmp;
    }
}

/* ------------------------------------------------------------------------ */
/*  rdp1(lambda)        returns random number distributed according to      */
/*                      a Poisson distribution with mean lambda.            */
/*                      Adapted from ACM 369.                               */
/*  Return -1 if lambda <= 0 or exceeds limits of exp().                    */
/*  Uses a local copy of random().                                          */
  
int rdp1(double lambda)
{
    register int k;
    double z,t;

    if (lambda <= 0.0 || lambda > ExpMax) 
        return(-1);

    z = exp(-lambda);
    k = 0;
    t = 1.0;
    while ((t *= random3()) > z)  
        k++;
    return(k);
}

