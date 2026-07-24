#include "tda_rhooks.h"
/****************************************************************************/
/*  t_lp                                                                    */
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
#include "t_gf.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_lp.c                                                     */

int lpf1(TDAContext *ctx, int m,int n,int p,double *tab,double *x,double *y,double *z,double tdd);
int phase1(TDAContext *ctx, int p);
void progamma(TDAContext *ctx);
int prophi(TDAContext *ctx);
int rowtrans(TDAContext *ctx, int ii,int jj);
int lpi(TDAContext *ctx, int n,int m,int nest,int *a0,int *b0,int *opts,int *itmp,int *vmax);
int lpia(TDAContext *ctx, int n,int m,int nest,int *an,int *an0,int **av,int **ar, int *b0,int *opts,int *itmp,int *vmax);
int lpia_get1j(TDAContext *ctx, int j,int *an,int **av,int **ar);
int lpia_get1j0(TDAContext *ctx, int j,int *an0,int **av,int **ar);
int lpia_getij(TDAContext *ctx, int i,int j,int *an,int **av,int **ar);
int lpia_getij0(TDAContext *ctx, int i,int j,int *an0,int **av,int **ar);
void lpia_null(TDAContext *ctx, int j,int *an);
void lpia_addb(TDAContext *ctx, int j,int *b,int *an,int **av,int **ar);
void lpia_addb0(TDAContext *ctx, int j,int *b,int *an0,int **av,int **ar);
void lpia_subb(TDAContext *ctx, int j,int *b,int *an,int **av,int **ar);
void lpia_subb0(TDAContext *ctx, int j,int *b,int *an0,int **av,int **ar);
void lpia_copy(TDAContext *ctx, int j,int *an,int *an0,int **av,int **ar);
void lpia_prn(TDAContext *ctx, int m,int n,int *an,int *an0,int **av,int **ar);


/* ------------------------------------------------------------------------ */
/*      lpf1()                                                              */
/*                                                                          */
/*      Lineare Programmierung, optional mit Equality Constraints           */
/*                                                                          */
/*      Ref. R.C. Salazar, S.K. Sen, Minit Algorithm for Linear Programming */
/*      Algorithm 333, Comm. ACM 11 (1968), 437 - 440.                      */
/*                                                                          */
/*      Es sei a eine (m,n)-Matrix, b ein (m)-Vektor, c ein (n)-Vektor      */
/*      Dann werden ein (n)-Vektor x und ein (m)-Vektor y berechnet, so     */
/*      dass:                                                               */
/*              c'x ---> max, mit x >= 0 und  ax <= b (primal)              */
/*              b'y ---> min, mit y >= 0 und a'y >= c (dual)                */
/*                                                                          */
/*      Die Funktion erwartet die Inputdaten in Form einer Tableau-Matrix   */
/*      die folgendermassen aufgebaut ist:                                  */
/*                                                                          */
/*                  c(1)      ...   c(n)      0                             */
/*                                                                          */
/*                  a(1, 1)  ...   a(1,n)   b(1)                            */
/*        tab[] =   a(2, 1)  ...   a(2,n)   b(2)                            */
/*                  ...                                                     */
/*                  a(m, 1)  ...   a(m,n)   b(m)                            */
/*                                                                          */
/*      Beachte: Es ist moeglich, einige der Restriktionen als Gleichheits- */
/*      Relationen zu definieren. Wenn dies der Fall sein soll, muss es     */
/*      sich um die letzten Zeilen der Tableau-Matrix handeln.              */
/*                                                                          */
/*      Beachte: wenn Equality Constraints definiert werden, werden fuer    */
/*      die entsprechenden Constraints keine Werte des dualen Loesungs-     */
/*      vektors berechnet. Sie werden dann = 0 gesetzt.                     */
/*                                                                          */
/*                                                                          */
/*      Aufruf der Funktion                                                 */
/*                                                                          */
/*      int lpf1(m,n,p,tab,x,y,z,tdd)                                       */
/*                                                                          */
/*      int m       input:  Anzahl der Constraints                          */
/*      int n       input:  Anzahl der Variablen                            */
/*      int p       input:  Anzahl der Equality-Constraints.                */
/*      double *tab input:  die Tableau-Matrix, wie oben beschrieben; d.h.  */
/*                          eine (m + 1,n + 1)-Matrix als ein eindimens.    */
/*                          Feld: tab(i,j) = tab[(i - 1) * (m + 1) + j]     */
/*                          wobei i = 1,m + 1; j = 1,n + 1.                 */
/*                  output: unveraendert                                    */
/*      double *x   input:  eindimensionales Feld x(i) (i = 1,n)            */
/*                  output: Loesungsvektor fuer das primale Problem         */
/*      double *y   input:  eindimensionales Feld y(i) (i = 1,m)            */
/*                  output: Loesungsvektor fuer das duale Problem           */
/*      double *z   output: Wert der Zielfunktion                           */
/*      double tdd  input:  eine Konstante; es wird angenommen dass ein     */
/*                          Element der Tableau-Matrix = 0 ist, wenn es     */
/*                          groesser als -tdd ist. Der Wert sollte der      */
/*                          Groessenordnung der Elemente der Tableau-Matrix */
/*                          entsprechen; etwa 1.e-8                         */
/*                                                                          */
/*      Return:     0  wenn erfolgreich                                     */
/*                 -1  wenn Fehler in den Inputparametern                   */
/*                 -2  wenn keine Loesung existiert (moeglicherweise        */
/*                     linear abhaengige Equality Constraints)              */
/*                 -3  wenn das primale Problem keine zulaessige Loesung    */
/*                     hat, weil die duale Zielfunktion unbeschraenkt ist   */
/*                 -4  wenn das duale Problem keine zulaessige Loesung      */
/*                     hat, weil die primale Zielfunktion unbeschraenkt ist */
/*                 -5  wenn insufficient memory.                            */
/*                                                                          */

/*  global variables */

static int k,l,m1,n1,lcol,im,iimax,jm,jmin;
static int *iimin,*jmax,*ind,*ind1,*chk;
static double ggmin,phimax,*e,*delmax,*thmin,td,big,big1;

int lpf1(TDAContext *ctx, int m,int n,int p,double *tab,double *x,double *y,double *z,double tdd)
{
    int i,j,r;

    if (m < 1 || n < 1 || p < 0 || p > m || tdd < 0.0)
        return(-1);

    r = 0;
    td = tdd;
    big = 1.e100;
    big1 = 1.e102;

    m1 = m + 1;
    n1 = n + 1;
    lcol = m + n - p + 1;

    r = -5;
    if (!(e = (double *) calloc ((size_t)(m1 * lcol + 1),sizeof(double))))  
        goto LPFin;
    memrq(ctx, m1 * lcol + 1,sizeof(double));

    if (!(delmax = (double *) calloc ((size_t)(m1 + 1),sizeof(double))))  
        goto LPFin1;
    memrq(ctx, m1 + 1,sizeof(double));

    if (!(thmin = (double *) calloc ((size_t)(lcol + 1),sizeof(double))))  
        goto LPFin2;
    memrq(ctx, lcol + 1,sizeof(double));

    if (!(iimin = (int *) calloc ((size_t)(lcol + 1),sizeof(int))))  
        goto LPFin3;
    memrq(ctx, lcol + 1,sizeof(int));

    if (!(jmax = (int *) calloc ((size_t)(m1 + 1),sizeof(int))))   
        goto LPFin4;
    memrq(ctx, m1 + 1,sizeof(int));

    if (!(ind = (int *) calloc ((size_t)(lcol + 1),sizeof(int))))  
        goto LPFin5;
    memrq(ctx, lcol + 1,sizeof(int));

    if (!(ind1 = (int *) calloc ((size_t)(m1 + 1),sizeof(int))))  
        goto LPFin6;
    memrq(ctx, m1 + 1,sizeof(int));

    if (!(chk = (int *) calloc ((size_t)(m1 + 1),sizeof(int))))
        goto LPFin7;
    memrq(ctx, m1 + 1,sizeof(int));

    r = 0;             

    /* Konstruktion der neuen Tableau-Matrix e */

    for (j = 1; j <= n; ++j)
        e[j] = -tab[j];
    for (j = n1; j <= lcol; ++j)
        e[j] = 0.0;
    for (i = 2; i <= m1; ++i) {
        for (j = 1; j <= n; ++j)
            e[(i - 1) * lcol + j] = tab[(i - 1) * n1 + j];
        for (j = n1; j < lcol; ++j) {
            if (j - n + 1 == i)
                e[(i - 1) * lcol + j] = 1.0;
            else
                e[(i - 1) * lcol + j] = 0.0;
        }
        e[i * lcol] = tab[i * n1];
    }

    /* hier beginnt der Algorithmus MINIT entsprechend der oben angegebenen */
    /* Quelle. Verwendung findet die Tableau-Matrix e                       */

    for (i = 2; i <= m1; ++i)  
        chk[i] = 0;

    if (p && phase1(ctx, p)) {
        r = -2;             /* no solution */
        goto LPFin;
    }
RCS:
    l = k = 1;
    for (j = 1; j < lcol; ++j) {
        if (e[j] < -td)  
            ind[l++] = j;
    }
    for (i = 2; i <= m1; ++i) {
        if (e[i * lcol] < -td)  
            ind1[k++] = i;
    }
    if (l == 1) {
        if (k == 1) {   /* found solution */

            *z = e[lcol];

            for (i = 1; i <= n; ++i)
                x[i] = 0.0;

            for (i = 2; i <= m1; ++i) {
                if (chk[i] > n)
                    chk[i] = 0;
                if (chk[i] > 0)
                    x[chk[i]] = e[i * lcol];
            }
            for (j = 1; j <= m; ++j)
                y[j] = 0.0;

            for (j = n1; j < lcol; ++j)
                y[j - n] = e[j];

            r = 0;
            goto LPFin;
        }
        else if (k == 2) {
            for (j = 1; j < lcol; ++j)   
                if (e[(ind1[1] - 1) * lcol + j] < 0.0)
                    goto R;

            /* primal problem has no feasible solution */     
            /* dual objective function is unbounded    */

            r = -3;
            goto LPFin;
        }
        else
            goto R;
    }
    else {
        if (l == 2) {
            if (k == 1) {
                for (i = 2; i <= m1; ++i)  
                    if (e[(i - 1) * lcol + ind[1]] > 0.0)
                        goto C;

                /* primal objective function is unbounded */
                /* dual problem has no feasible solution  */

                r = -4;
                goto LPFin;
            }
            else
                goto S;
        }
        if (k == 1)
            goto C;
        else 
            goto S;
    }
R:
    prophi(ctx);
    if (rowtrans(ctx, iimax,jm)) {
        r = -2;             /* no solution */
        goto LPFin;
    }
    goto RCS;
C:
    progamma(ctx);
    if (rowtrans(ctx, im,jmin)) {
        r = -2;             /* no solution */
        goto LPFin;
    }
    goto RCS;
S:
    progamma(ctx);
    prophi(ctx);
    if (ggmin == big) {
        if (rowtrans(ctx, iimax,jm)) {
            r = -2;             /* no solution */
            goto LPFin;
        }
        goto RCS;
    }
    if (phimax == -big) {
        if (rowtrans(ctx, im,jmin)) {
            r = -2;             /* no solution */
            goto LPFin;
        }
        goto RCS;
    }
    if (fabs(phimax) > fabs(ggmin)) {
        if (rowtrans(ctx, iimax,jm)) {
            r = -2;             /* no solution */
            goto LPFin;
        }
    }
    else {
        if (rowtrans(ctx, im,jmin)) {
            r = -2;             /* no solution */
            goto LPFin;
        }
    }
    goto RCS;

LPFin:
    free((char *)chk);
    memrq(ctx, -m1 - 1,sizeof(int));
LPFin7:
    free((char *)ind1);
    memrq(ctx, -m1 - 1,sizeof(int));
LPFin6:
    free((char *)ind);
    memrq(ctx, -lcol - 1,sizeof(int));
LPFin5:
    free((char *)jmax);
    memrq(ctx, -m1 - 1,sizeof(int));
LPFin4:
    free((char *)iimin);
    memrq(ctx, -lcol - 1,sizeof(int));
LPFin3:
    free((char *)thmin);
    memrq(ctx, -lcol - 1,sizeof(double));
LPFin2:
    free((char *)delmax);
    memrq(ctx, -m1 - 1,sizeof(double));
LPFin1:
    free((char *)e);
    memrq(ctx, -m1 * lcol - 1,sizeof(double));
    return(r);
}

/* -##--------------------------------------------------------------------- */
/*  phase1                                                                  */
/*      Applied only to equality constraints if any; number of equality     */
/*      Constraints is p.                                                   */
/*      Returns: 0 if no error; -2 if no solution                           */

int phase1(TDAContext *ctx, int p)
{
    int i,j,r,first,im1,jmin1;
    double gamma,theta;

    im1 = jmin1 = first = 0;

    /* da der Algorithmus voraussetzt, dass die Elemente des b-Vektors  */
    /* nicht-negativ sein duerfen, wenn es sich um equality constraints */
    /* handelt, wird hier eine entsprechende Korrektur vorgenommen:     */

    for (i = m1 - p + 1; i <= m1; ++i) {
        if (e[i * lcol] < 0.0) {
            for (j = 1; j <= lcol; ++j)
                e[(i - 1) * lcol + j] = -e[(i - 1) * lcol + j];
        }
    }
    for (r = 1; r <= p; ++r) {
        ggmin = big;
        l = 1;
        jmin = 0; 
        first = 1;
        for (j = 1; j < n1; ++j) {
            thmin[j] = big;
            if (e[j] < -td)
                ind[l++] = j;
        }
L1:     if (l == 1) {
            for (j = 1; j < n1; ++j)
                ind[j] = j;
            l = n1;
        }
        for (k = 1; k < l; ++k) {
            for (i = m1 - p + 1; i <= m1; ++i) {
                if (chk[i] == 0) {
                    if (e[(i - 1) * lcol + ind[k]] > td) {
                        theta = e[i * lcol] / e[(i - 1) * lcol + ind[k]];
                        if (theta < thmin[ind[k]]) {
                            thmin[ind[k]] = theta;
                            iimin[ind[k]] = i;
                        }
                    }
                }
            }
            if (thmin[ind[k]] >= big)
                gamma = big1;
            else {
                gamma = thmin[ind[k]] * e[ind[k]];

                if (gamma < ggmin && thmin[ind[k]] < big) {
                    ggmin = gamma;
                    jmin = ind[k];
                }
            }
        }
        if (jmin == 0) {
            if (first) {
                first = 0;
                l = 1;
                goto L1;
            }
            else
                im = 0;
        }
        else
            im = iimin[jmin];

        if (im == im1 && jmin == jmin1) {
            l = 1;
            goto L1;
        }
        if (rowtrans(ctx, im,jmin)) {
            return(-2);         /* no solution */
        }
        im1 = im;
        jmin1 = jmin;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  progamma()                                                              */
/*      Performs calculations over columns to determine the pivot element   */
/*      Returns: nothing                                                    */

void progamma(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,l1;
    double gamma,theta;

    ggmin = big;         /* any "large" number */
    jmin = 0;

    for (l1 = 1; l1 < l; ++l1) {
        iimin[ind[l1]] = 0;
        thmin[ind[l1]] = big;

        for (i = 2; i <= m1; ++i) {
            if (e[(i - 1) * lcol + ind[l1]] > td && e[i * lcol] >= -td) {
                theta = e[i * lcol] / e[(i - 1) * lcol + ind[l1]];
                if (theta < thmin[ind[l1]]) {
                    thmin[ind[l1]] = theta;
                    iimin[ind[l1]] = i;
                }
            }
        }
        if (thmin[ind[l1]] == big)
            gamma = big1;
        else
            gamma = thmin[ind[l1]] * e[ind[l1]];
        if (gamma < ggmin) {
            ggmin = gamma;
            jmin = ind[l1];
        }
    }
    if (jmin > 0)
        im = iimin[jmin];
}

/* ------------------------------------------------------------------------ */
/*  prophi()                                                                */
/*      Performs calculations over rows to determine the pivot element      */
/*      Returns: nothing                                                    */

int prophi(TDAContext *ctx)
{
    (void)ctx;        /* unused: the signature is shared */
    int j,k1;
    double delta,phi;

    phimax = -big;
    iimax = 0;

    for (k1 = 1; k1 < k; ++k1) {
        jmax[ind1[k1]] = 0;
        delmax[ind1[k1]] = -big;

        for (j = 1; j < lcol; ++j) {
            if (e[(ind1[k1] - 1) * lcol + j] < -td && e[j] >= -td) {
                delta = e[j] / e[(ind1[k1] - 1) * lcol + j];
                if (delta > delmax[ind1[k1]]) {
                    delmax[ind1[k1]] = delta;
                    jmax[ind1[k1]] = j;
                }
            }
        }
        if (delmax[ind1[k1]] == -big)
            phi = -big1;
        else
            phi = delmax[ind1[k1]] * e[ind1[k1] * lcol];

        if (phi > phimax) {
            phimax = phi;
            iimax = ind1[k1];    
        }
    }
    if (iimax > 0)
        jm = jmax[iimax];
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  rowtrans(ii,jj)                                                         */
/*     Performs the tableau transformations; (ii,jj) is the pivotal element */
/*     Returns: 0 if no error; -2 if no solution                            */

int rowtrans(TDAContext *ctx, int ii,int jj)
{
    (void)ctx;        /* unused: the signature is shared */
    int i,j;
    double tmp;  

    if (!ii || !jj)
        return(-2);   /* no solution */

    tmp = e[(ii - 1) * lcol + jj];

    for (j = 1; j <= lcol; ++j)
        e[(ii - 1) * lcol + j] /= tmp;  

    for (i = 1; i <= m1; ++i) {
        if (i != ii) {
            if (e[(i - 1) * lcol + jj] != 0.0) {
                tmp = e[(i - 1) * lcol + jj];
                for (j = 1; j <= lcol; ++j)
                    e[(i - 1) * lcol + j] -= e[(ii - 1) * lcol + j] * tmp;  
            }
        }
    }
    chk[ii] = jj;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*      lpi()                                                               */
/*                                                                          */
/*      Linear programming in 0-1 variables.                                */
/*      Adapted from CACM 449.                                              */
/*                                                                          */
/*      a0[m,n]                                                             */
/*      b0[m]                                                               */
/*      opts[nest,n]                                                        */
/*                                                                          */
/*      nest = upper limit to the number of solutions in opts[].            */
/*                                                                          */
/*      itmp is auxiliary integer storage of length                         */
/*      4 * m + 6 * n + m * n + 2                                           */
/*                                                                          */
/*      Return: -1 if constraints are inconsistent (no solution)            */
/*              -2 if nest is too small                                     */
/*              -3 wrong input arguments                                    */ 
/*              >0 number of solutions in opts[]                            */  
/*                                                                          */

int lpi(TDAContext *ctx, int n,int m,int nest,int *a0,int *b0,int *opts,int *itmp,int *vmax)
{
    register int i,j,k1,ij,j1;
    int in,nopt,ns,ni,nat,v,vneg,itest,l1,newv,is,it1,it2,ip;
    int *b,*b1,*s1,*ind_,*c,*s,*t,*s0,*bc,*a,*x;



    if (nest < 1)
        return(-3);

    b = itmp;  
    b1 = b + m;
    s1 = b1 + m; 
    ind_ = s1 + m;
    c = ind_ + m;
    x = c + n;
    s = x + n;
    t = s + n;
    s0 = t + n;
    bc = s0 + n;
    a = bc + n;
        
    nopt = ns = ni = nat = 0;
    *vmax = 0;

    for (j = 1; j <= n; ++j) 
        t[j] = 0;

    for (i = 1; i <= m; ++i) {
        b[i] = b0[i];
        in = (i - 1) * n;
        for (j = 1; j <= n; ++j)
            a[in + j] = a0[in + j];
    }
    vneg = -1;
    for (j = 1; j <= n; ++j) {
        x[j] = 2;
        if (a[j] < 0)
            vneg += a[j];
    }
    b[1] = vneg;
    v = vneg;
L50:
    for (i = 1; i <= m; ++i)
        ind_[i] = 0;
L70:
    for (i = 1; i <= m; ++i)
        b1[i] = b[i];
    ni++;
    itest = 1;
    for (i = 1; i <= m; ++i) {
        in = (i - 1) * n;
        s1[i] = 0;
        if (ind_[i] == 1)
            continue;

        for (j = 1; j <= n; ++j) {
            if (a[in + j] < 0)
                b1[i] -= a[in + j]; 
            s1[i] += iabs(ctx, a[in + j]);
        }
        if (b1[i] <= 0)
            ind_[i] = 1;
        else
            itest = 0;
    }
    if (itest == 1)
        goto L420;

    for (i = 1; i <= m; ++i) {
        if (ind_[i] == 1)
            continue;
        if (s1[i] < b1[i])   
            goto L560;
    }
    i = 1;
L130:
    if (ind_[i] == 1)
        goto L360;
    if (s1[i] > b1[i])
        goto L200; 

L140:
    for (j = 1; j <= n; ++j) {
        if (a[(i-1)*n + j] == 0)
            continue;
        ns++;
        bc[ns] = 1;
        if (a[(i-1)*n + j] < 0) {
            s[ns] = -j;
            x[j] = 0;
        }
        else {
            s[ns] = j;
            x[j] = 1;
            for (ij = 1; ij <= m; ++ij)  
                b[ij] -= a[(ij-1)*n + j];
        }
        for (ij = 1; ij <= m; ++ij)
            a[(ij-1)*n + j] = 0;
    }

    goto L70;

L200:
    for (j = 1; j <= n; ++j) 
        c[j] = iabs(ctx, a[(i-1)*n + j]);
    l1 = 1;
L220:
    j = l1 + 1;
L230:
    if (c[l1] < c[j]) {
        ip = c[l1];
        c[l1] = c[j];
        c[j] = ip;
    }
    j++;
    if (j <= n)
        goto L230;

    l1++;
    if (l1 < n)
        goto L220;
L260:
    if (c[l1] > 0)
        goto L270;
    l1--;
    goto L260;
L270:
    if (s1[i] - c[l1] < b1[i])
        goto L140;
    if (s1[i] - c[1] - b1[i] >= 0)
        goto L360;
    ns++;
    bc[ns] = 1;
L280:
    for (j = 1; j <= n; ++j) {
        if (iabs(ctx, a[(i - 1) * n + j]) == c[1])
            break;          
    }
    if (a[(i-1)*n + j] < 0)
        goto L330;
L310:
    s[ns] = j; 
    x[j] = 1;
    for (ij = 1; ij <= m; ++ij)
        b[ij] -= a[(ij - 1) * n + j];
    goto L340;
L330:
    s[ns] = -j;
    x[j] = 0;
L340:
    for (ij = 1; ij <= m; ++ij)
        a[(ij - 1) * n + j] = 0;
    goto L70;

L360:
    i++;
    if (i <= m)
        goto L130;
    if (ns == n)
        goto L480;

    for (j = 1; j <= n; ++j)
        c[j] = iabs(ctx, a[j]);

    for (j = 2; j <= n; ++j) {
        if (c[1] >= c[j])
            continue;  
        c[1] = c[j];
    }
    if (c[1] == 0)
        goto L390;
    ns++;
    bc[ns] = 0;
    i = 1;
    goto L280;
L390:
    for (j = 1; j <= n; ++j) {
        for (j1 = 1; j1 <= ns; ++j1) {
            if (j == iabs(ctx, s[j1]))
                goto L410;
        }
        ns++;
        bc[ns] = 0;
        goto L310;
L410: ;
    }
L420:
    for (j = 1; j <= n; ++j) {
        if (ns == n)
            goto L480;
        if ((x[j] != 2) || (a[j] == 0))
            continue;  
        ns++;
        bc[ns] = 1;
        if (a[j] < 0)
            goto L440;
        s[ns] = j;
        x[j] = 1;
        for (i = 1; i <= m; ++i)
            b[i] -= a[(i - 1) * n + j];
        goto L450;
L440:
        s[ns] = -j;
        x[j] = 0;
L450:
        for (i = 1; i <= m; ++i)
            a[(i - 1) * n + j] = 0;
    }
L480:
    newv = 0;
    for (j = 1; j <= n; ++j) 
        newv += x[j] * a0[j];
    for (j = 1; j <= ns; ++j) {
        k1 = ns + 1 - j;
        if (bc[k1] == 0)
            t[k1] = 1;
    }
    if (newv > v)
        goto L510;
    nopt++;
    if (nopt <= nest)
        goto L540;
    return(-2);

L510:
    nopt = 1;
    v = newv;
    b[1] = v;
    for (j = 1; j <= n; ++j) {
        if (x[j] != 1)
            continue;  
        b[1] -= a0[j];
    }
    for (j = 1; j <= n; ++j)
        s0[j] = s[j];
L540:
    for (j = 1; j <= n; ++j)  
        opts[(nopt-1) * n + j] = x[j];
         
L560:
    if (ns == 0)
        goto L580;
    is = 0;
    for (j = 1; j <= ns; ++j)  
        is += bc[j];
          
    if (is < ns)
        goto L600;
    if (v > vneg)
        goto L590;
L580:
    return(-1);           

L590:
    v += b0[1];
    *vmax = v;
    return(nopt);             

L600:
    for (j1 = 1; j1 <= ns; ++j1) {
        k1 = ns + 1 - j1;
        if (bc[k1] == 0)
            goto L620;
    }
L620:
    if (t[k1] == 1)
        goto L750;
L630:
    for (j1 = k1; j1 <= ns; ++j1) {
        for (j = 1; j <= n; ++j) {
            if (j == iabs(ctx, s[j1]))
                goto L650;
        }
L650:
        if (k1 == j1)
            goto L700;
        if (x[j] == 1)
            goto L670;
        for (i = 1; i <= m; ++i) {
            a[(i-1)*n + j] = a0[(i-1)*n + j];
        }
        goto L690;
L670:
        for (i = 1; i <= m; ++i) {
            a[(i-1)*n + j] = a0[(i-1)*n + j];
            b[i] += a[(i-1)*n + j];
        }
L690:
        x[j] = 2;
        goto L740;
L700:
        s[k1] = -s[k1];
        bc[k1] = 1;
        x[j] = 1 - x[j];
        if (x[j] == 0)
            goto L720;

        for (i = 1; i <= m; ++i)  
            b[i] -= a0[(i-1)*n + j];

        goto L740;

L720:
        for (i = 1; i <= m; ++i)
            b[i] += a0[(i-1)*n + j];
L740: ;
    }
    ns = k1;
    goto L50;

L750:
    t[k1] = 0;
    it1 = it2 = 0;
    for (j1 = k1; j1 <= n; ++j1) {
        for (j = 1; j <= n; ++j) {
            if (j == iabs(ctx, s0[j1]))
                goto L770;
        }
L770:
        if (k1 == j1)
            goto L780;
        if (((x[j] == 0) && (a0[j] > 0)) ||
            ((x[j] == 1) && (a0[j] < 0)))
                it2 += iabs(ctx, a0[j]);
        goto L790;
L780:
        it1 = iabs(ctx, a0[j]);
L790: ;
    }
    if (it1 <= it2)
        goto L630;
    bc[k1] = 1;
    nat++;
    goto L560;
}

/* -###-------------------------------------------------------------------- */
/*      lpia()                                                              */
/*                                                                          */
/*      Linear programming in 0-1 variables.                                */
/*      Adapted from CACM 449.                                              */
/*                                                                          */
/*      a0[m,n]                                                             */
/*      b0[m]                                                               */
/*      opts[nest,n]                                                        */
/*                                                                          */
/*      nest = upper limit to the number of solutions in opts[].            */
/*                                                                          */
/*      itmp is auxiliary integer storage of length                         */
/*      4 * m + 6 * n + m * n + 2                                           */
/*                                                                          */
/*      Return: -1 if constraints are inconsistent (no solution)            */
/*              -2 if nest is too small                                     */
/*              -3 wrong input arguments                                    */ 
/*              >0 number of solutions in opts[]                            */  
/*                                                                          */

int lpia(TDAContext *ctx, int n,int m,int nest,int *an,int *an0,int **av,int **ar, int *b0,int *opts,int *itmp,int *vmax)
{
    register int i,j,k1,j1;
    int nopt,ns,ni,nat,v,vneg,itest,l1,newv,is,it1,it2,ip,atmp = 0,a0j;
    int *b,*b1,*s1,*ind_,*c,*s,*t,*s0,*bc,*x;

    if (nest < 1)
        return(-3);

    b = itmp;  
    b1 = b + m;
    s1 = b1 + m; 
    ind_ = s1 + m;
    c = ind_ + m;
    x = c + n;
    s = x + n;
    t = s + n;
    s0 = t + n;
    bc = s0 + n;
        
    nopt = ns = ni = nat = 0;
    *vmax = 0;

    for (j = 1; j <= n; ++j) 
        t[j] = 0;

         
    for (i = 1; i <= m; ++i)    
        b[i] = b0[i];
         
    for (j = 1; j <= n; ++j)   
        lpia_copy(ctx, j,an,an0,av,ar);

    vneg = -1;
    for (j = 1; j <= n; ++j) {
        x[j] = 2;
        if ((atmp = lpia_get1j(ctx, j,an,av,ar)) < 0)
            vneg += atmp;
    }
    b[1] = vneg;
    v = vneg;
L50:
    tda_err("L50\n");         
    for (i = 1; i <= m; ++i)
        ind_[i] = 0;
L70:
/*  tda_err("  L70\n");      */
    for (i = 1; i <= m; ++i)
        b1[i] = b[i];
    ni++;
    itest = 1;
    for (i = 1; i <= m; ++i) {
        s1[i] = 0;
        if (ind_[i] == 1)
            continue;

        for (j = 1; j <= n; ++j) {
            atmp = lpia_getij(ctx, i,j,an,av,ar);
            if (atmp < 0)
                b1[i] -= atmp; 
            s1[i] += iabs(ctx, atmp);
        }
        if (b1[i] <= 0)
            ind_[i] = 1;
        else
            itest = 0;
    }
    if (itest == 1)
        goto L420;

    for (i = 1; i <= m; ++i) {
        if (ind_[i] == 1)
            continue;
        if (s1[i] < b1[i])   
            goto L560;
    }
    i = 1;
L130:
    if (ind_[i] == 1)
        goto L360;
    if (s1[i] > b1[i])
        goto L200; 

L140:
    for (j = 1; j <= n; ++j) {

        atmp = lpia_getij(ctx, i,j,an,av,ar);

        /***************************
        if (a[(i-1)*n + j] == 0)
            continue;
        **********************/

        if (atmp == 0)
            continue;

        ns++;
        bc[ns] = 1;
/*****  if (a[(i-1)*n + j] < 0) {   ****/

        if (atmp < 0) {
            s[ns] = -j;
            x[j] = 0;
        }
        else {
            s[ns] = j;
            x[j] = 1;

            lpia_subb(ctx, j,b,an,av,ar);
            /*****************************
            for (ij = 1; ij <= m; ++ij)  
                b[ij] -= a[(ij-1)*n + j];
            *****************************/
        }
        lpia_null(ctx, j,an);
        /********************
        for (ij = 1; ij <= m; ++ij)
            a[(ij-1)*n + j] = 0;
        ****************************/
    }

    goto L70;

L200:
    for (j = 1; j <= n; ++j) {
        /**  c[j] = iabs(a[(i-1)*n + j]); **/
        atmp = lpia_getij(ctx, i,j,an,av,ar);
        c[j] = iabs(ctx, atmp);
    }
    l1 = 1;
L220:
    j = l1 + 1;
L230:
    if (c[l1] < c[j]) {
        ip = c[l1];
        c[l1] = c[j];
        c[j] = ip;
    }
    j++;
    if (j <= n)
        goto L230;

    l1++;
    if (l1 < n)
        goto L220;
L260:
    if (c[l1] > 0)
        goto L270;
    l1--;
    goto L260;
L270:
    if (s1[i] - c[l1] < b1[i])
        goto L140;
    if (s1[i] - c[1] - b1[i] >= 0)
        goto L360;
    ns++;
    bc[ns] = 1;
L280:
    for (j = 1; j <= n; ++j) {
        atmp = lpia_getij(ctx, i,j,an,av,ar);
        if (iabs(ctx, atmp) == c[1])
            break;          
    }
    if (atmp < 0)
        goto L330;
L310:
    s[ns] = j; 
    x[j] = 1;

    lpia_subb(ctx, j,b,an,av,ar);
    /********************************
    for (ij = 1; ij <= m; ++ij)
        b[ij] -= a[(ij - 1) * n + j];
    *********************************/

    goto L340;
L330:
    s[ns] = -j;
    x[j] = 0;
L340:
    lpia_null(ctx, j,an);
    /****************************
    for (ij = 1; ij <= m; ++ij)
        a[(ij - 1) * n + j] = 0;
    ****************************/
    goto L70;

L360:
    i++;
    if (i <= m)
        goto L130;
    if (ns == n)
        goto L480;

    for (j = 1; j <= n; ++j)
        c[j] = iabs(ctx, lpia_get1j(ctx, j,an,av,ar));

    for (j = 2; j <= n; ++j) {
        if (c[1] >= c[j])
            continue;  
        c[1] = c[j];
    }
    if (c[1] == 0)
        goto L390;
    ns++;
    bc[ns] = 0;
    i = 1;
    goto L280;
L390:
    for (j = 1; j <= n; ++j) {
        for (j1 = 1; j1 <= ns; ++j1) {
            if (j == iabs(ctx, s[j1]))
                goto L410;
        }
        ns++;
        bc[ns] = 0;
        goto L310;
L410: ;
    }
L420:
    for (j = 1; j <= n; ++j) {
        if (ns == n)
            goto L480;
        if ((x[j] != 2) || ((atmp = lpia_get1j(ctx, j,an,av,ar)) == 0))
            continue;  
        ns++;
        bc[ns] = 1;
        if ((atmp = lpia_get1j(ctx, j,an,av,ar)) < 0)
            goto L440;
        s[ns] = j;
        x[j] = 1;

        lpia_subb(ctx, j,b,an,av,ar);
        /*******************************
        for (i = 1; i <= m; ++i)
            b[i] -= a[(i - 1) * n + j];
        *******************************/
        goto L450;
L440:
        s[ns] = -j;
        x[j] = 0;
L450:
        lpia_null(ctx, j,an);
        /****************************
        for (i = 1; i <= m; ++i)
            a[(i - 1) * n + j] = 0;
        *****************************/
    }
L480:
    newv = 0;
    for (j = 1; j <= n; ++j) {
        newv += x[j] * lpia_get1j0(ctx, j,an0,av,ar);

/*      newv += x[j] * a0[j];           */ 

    }

    for (j = 1; j <= ns; ++j) {
        k1 = ns + 1 - j;
        if (bc[k1] == 0)
            t[k1] = 1;
    }
    if (newv > v)
        goto L510;
    nopt++;
    if (nopt <= nest)
        goto L540;
    return(-2);

L510:
    nopt = 1;
    v = newv;
    b[1] = v;
    for (j = 1; j <= n; ++j) {
        if (x[j] != 1)
            continue;  
        b[1] -= lpia_get1j0(ctx, j,an0,av,ar);

/**     b[1] -= a0[j];          ***/
    }
    for (j = 1; j <= n; ++j)
        s0[j] = s[j];
L540:
    for (j = 1; j <= n; ++j)  
        opts[(nopt-1) * n + j] = x[j];
         
L560:
    if (ns == 0)
        goto L580;
    is = 0;

    for (j = 1; j <= ns; ++j) {
        is += bc[j];
    }    
    if (is < ns)
        goto L600;
    if (v > vneg)
        goto L590;
L580:
    return(-1);           

L590:
    v += b0[1];
    *vmax = v;

    return(nopt);             

L600:
    tda_err("  L600\n");     
    for (j1 = 1; j1 <= ns; ++j1) {
        k1 = ns + 1 - j1;
        if (bc[k1] == 0)
            goto L620;
    }
L620:
    if (t[k1] == 1)
        goto L750;
L630:
    for (j1 = k1; j1 <= ns; ++j1) {
        for (j = 1; j <= n; ++j) {
            if (j == iabs(ctx, s[j1]))
                goto L650;
        }
L650:
        if (k1 == j1)
            goto L700;
        if (x[j] == 1)
            goto L670;


        lpia_copy(ctx, j,an,an0,av,ar); 
        /***************************
        for (i = 1; i <= m; ++i) {
            a[(i-1)*n + j] = a0[(i-1)*n + j];
        }
        ************************************/

        goto L690;
L670:
        lpia_copy(ctx, j,an,an0,av,ar); 
        lpia_addb(ctx, j,b,an,av,ar);
        /************************************
        for (i = 1; i <= m; ++i) {
            a[(i-1)*n + j] = a0[(i-1)*n + j];
            b[i] += a[(i-1)*n + j];
        }
        *************************************/

L690:
        x[j] = 2;
        goto L740;
L700:
        s[k1] = -s[k1];
        bc[k1] = 1;
        x[j] = 1 - x[j];
        if (x[j] == 0)
            goto L720;

        lpia_subb0(ctx, j,b,an0,av,ar);
          
/******
        for (i = 1; i <= m; ++i)
            b[i] -= lpia_getij0(ctx, i,j,an0,av,ar);
          
**/
        /***************************
        for (i = 1; i <= m; ++i)  
            b[i] -= a0[(i-1)*n + j];
        *****************************/

        goto L740;
L720:
        lpia_addb0(ctx, j,b,an0,av,ar);                   
        /****************************
        for (i = 1; i <= m; ++i)
            b[i] += a0[(i-1)*n + j];
        *****************************/

L740: ;
    }

    ns = k1;
    goto L50;

L750:
    t[k1] = 0;
    it1 = it2 = 0;
    for (j1 = k1; j1 <= n; ++j1) {
        for (j = 1; j <= n; ++j) {
            if (j == iabs(ctx, s0[j1]))
                goto L770;
        }
L770:
        a0j = lpia_get1j0(ctx, j,an0,av,ar);
        if (k1 == j1)
            goto L780;

        if (((x[j] == 0) && (a0j > 0)) ||
            ((x[j] == 1) && (a0j < 0)))
                it2 += iabs(ctx, a0j);
        goto L790;
L780:
        it1 = iabs(ctx, a0j);
L790: ;
    }
    if (it1 <= it2)
        goto L630;
    bc[k1] = 1;
    nat++;
    goto L560;
}

/* -###-------------------------------------------------------------------- */
/*  lpia_get1j(j,an,av,ar)   return a[1,j]                                  */
/*                                                                          */

int lpia_get1j(TDAContext *ctx, int j,int *an,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    if (an[j] > 0 && ar[j][1] == 1)
        return(av[j][1]);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  lpia_get1j0(j,an0,av,ar)  return a0[1,j]                                */
/*                                                                          */

int lpia_get1j0(TDAContext *ctx, int j,int *an0,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    if (an0[j] > 0 && ar[j][an0[j] + 1] == 1)
        return(av[j][an0[j] + 1]);
    else
        return(0);
}

/* ------------------------------------------------------------------------ */
/*  lpia_getij(i,j,an,av,ar)   return a[i,j]                                */
/*                                                                          */

int lpia_getij(TDAContext *ctx, int i,int j,int *an,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k1;

    for (k1 = 1; k1 <= an[j]; ++k1) {
        if (ar[j][k1] == i)
            return(av[j][k1]);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  lpia_getij0(i,j,an,av,ar)   return a0[i,j]                              */
/*                                                                          */

int lpia_getij0(TDAContext *ctx, int i,int j,int *an0,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k1,kk;

    for (k1 = 1; k1 <= an0[j]; ++k1) {
        kk = k1 + an0[j];
        if (ar[j][kk] == i)
            return(av[j][kk]);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  lpia_null(j,an,av,ar)   set column j to zero                            */
/*                                                                          */

void lpia_null(TDAContext *ctx, int j,int *an)
{
    (void)ctx;        /* unused: the signature is shared */
    an[j] = 0;
}

/* ------------------------------------------------------------------------ */
/*  lpia_addb(j,b,an,av,ar)  add column j of a[] to b[]                     */
/*                                                                          */

void lpia_addb(TDAContext *ctx, int j,int *b,int *an,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    for (i = 1; i <= an[j]; ++i)  
        b[ar[j][i]] += av[j][i];
}

/* ------------------------------------------------------------------------ */
/*  lpia_addb0(j,b,an0,av,ar)  add column j of a0[] to b[]                  */
/*                                                                          */

void lpia_addb0(TDAContext *ctx, int j,int *b,int *an0,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,ii;

    for (i = 1; i <= an0[j]; ++i) {
        ii = an0[j] + i;
        b[ar[j][ii]] += av[j][ii];
    }
}

/* ------------------------------------------------------------------------ */
/*  lpia_subb(j,b,an,av,ar)  subtract column j of a[] to b[]                */
/*                                                                          */

void lpia_subb(TDAContext *ctx, int j,int *b,int *an,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;

    for (i = 1; i <= an[j]; ++i)  
        b[ar[j][i]] -= av[j][i];
}

/* ------------------------------------------------------------------------ */
/*  lpia_subb0(j,b,an0,av,ar)  subtract column j of a0[] to b[]             */
/*                                                                          */

void lpia_subb0(TDAContext *ctx, int j,int *b,int *an0,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,ii;

    for (i = 1; i <= an0[j]; ++i) {
        ii = an0[j] + i;
        b[ar[j][ii]] -= av[j][ii];
    }
}

/* ------------------------------------------------------------------------ */
/*  lpia_copy(j,an,an0,av,ar)  copy column j from a0[] to a[].              */
/*                                                                          */

void lpia_copy(TDAContext *ctx, int j,int *an,int *an0,int **av,int **ar)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i,ii;

    for (i = 1; i <= an0[j]; ++i) {
        ii = an0[j] + i;
        av[j][i] = av[j][ii];
        ar[j][i] = ar[j][ii];
    }
    an[j] = an0[j];
}

void tda_reset_t_lp(void)
{
    k = l = m1 = n1 = lcol = im = iimax = jm = jmin = 0;
    iimin = jmax = ind = ind1 = chk = NULL;
    e = delmax = thmin = NULL;
    ggmin = phimax = td = big = big1 = 0.0;
}
