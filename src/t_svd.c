/****************************************************************************/
/*  t_svd                                                                   */
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
#include "t_graph.h"
#include "t_rand.h"

/* ------------------------------------------------------------------------ */
/*  functions in t_svd.c                                                    */

int eigen(int n,double *a,double *evr,double *evi,double *vec,int bal);
void balance(int n,double *a,double *d,int *low,int *hi);
void balbak(int n,double *a,double *d,int low,int hi);
void orthes(int n,double *a,int k,int l,double *d);
void ortrans(int n,double *h,double *d,int low,int upp,double *v);
int hqr2(int n,double *h,int low,int upp,double *vecs,double *wr,double *wi);
void cdiv(double xr,double xi,double yr,double yi,double *zr,double *zi);
int svdecomp(int m,int n,double *a,double *q,double *u,double *v,int mode);
int ginv(int m,int n,double *a);
int ginv1(int m,int n,double *a,double *sv);
int evecf(int n, double *d, double *z);
void tridf1(int n, double *d, double *e, double *z);
int evectf(int n, double *d, double *e, double *z);
int simitz(int n,int p,int km,double eps,int em,double *x,double *d,int *nit);
double g_ranf(void);
double g_ip(int n,double *z,double *w);
void g_op(int n,double *z,double *w);
int g_imtql2(int nm,int n,double *d,double *e,double *z,int mxit);
void g_tred2(int nm,int n,double *a,double *d,double *e,double *z);
int eigen2(int n,double *a,double *evr,double *evi,double *vec,int mode,int itmax);
int eigen1(int n,double *a,double *evr,double *evi,
    double *vecr,double *veci,int *indic);
double scale(int n,double *a,double *veci,double *prfact);
double hesqr(int n,double *a,double *veci,double *evr,double *evi,
    double *subdia,int *indic,double ex);
int realve(int n,int m,int ivec,double *a,double *vecr,double *evr,double *evi,
    int *iwork,double *work,int *indic,double eps,double ex);
int compve(int n,int m,int ivec,double *a,double *vecr,double *veci,double *evr,
    double *evi,int *indic,int *iwork,double *subdia,double *work1,
    double *work2,double *work,double eps,double ex);
void prxxx(char *txt,int n,int m,double *x);


/* ------------------------------------------------------------------------ */
/*  Function: eigen()                                                       */
/*                                                                          */
/*      Berechnung der Eigenwerte und Eigenvektoren einer allgemeinen       */
/*      reellen Matrix. Die Berechnung erfolgt in folgenden Schritten:      */
/*                                                                          */
/*      (1) Die Matrix wird mit der Funktion balance "balanciert", so       */
/*          dass eine moeglichst geringe Norm erreicht wird. (Dieser        */
/*          Schritt ist optional, vgl. den Parameter bal.)                  */
/*      (2) Die Matrix wird mit der Funktion orthes auf eine Hessenberg-    */
/*          Form gebracht.                                                  */
/*      (3) Es wird die Matrix der akkumulierten Transformationen be-       */
/*          rechnet, mit der Funktion: ortrans.                             */
/*      (4) Mit der Funktion hqr2 werden die Eigenwerte und Eigenvektoren   */
/*          der Hessenberg-Matrix berechnet.                                */
/*      (5) Wenn die Inputmatrix mit balance "balanciert" worden ist,       */
/*          werden die Eigenvektoren fuer die urspruengliche Matrix mit     */
/*          der Funktion balbak zurueckgewonnen.                            */
/*      (6) Die Eigenvektoren werden normiert. Reelle Eigenwerte werden     */
/*          so normiert, dass ihre euklidische Norm = 1 ist. Komplexe       */
/*          Eigenwerte werden so normiert, dass die Komponente mit dem      */
/*          groessten Absolutbetrag den Realteil 1 und den Imaginaerteil 0  */
/*          erhaelt. Dies entspricht der Normierungsmethode des Algorithmus */
/*          343 in Comm. ACM 11 (1963), 820ff.                              */
/*                                                                          */
/*          Ref.: Peters, G., Wilkinson, J.H. (1971). Eigenvectors of       */
/*          Real and Complex Matrices by LR and QR Triangularizations.      */
/*          In: J.\,H. Wilkinson, C. Reinsch, eds.,                         */
/*          Linear Algebra, pp. 372 -- 395. New York: Springer.             */
/*                                                                          */
/*      int eigen(n,a,evr,evi,vec,bal)                                      */
/*                                                                          */
/*      int n         input:  Ordnung der Inputmatrix                       */
/*      double *a     input:  die (n,n)-Inputmatrix als ein eindimens.      */
/*                            Feld a(i,j) = a[(i - 1) * n + j]              */
/*                    output: undefiniert                                   */
/*      double *evr   input:  eindimensionales Feld evr(i) (i = 1,n)        */
/*                    output: die Realteile der Eigenwerte                  */
/*      double *evi   input:  eindimensionales Feld evi(i) (i = 1,n)        */
/*                    output: die Imaginaerteile der Eigenwerte             */
/*      double *vec   input:  Speicherplatz fuer eine (n,n)-Matrix          */
/*                    output: die Eigenvektoren, spaltenweise im Feld       */
/*                            vec(i,j) = vec[(i - 1) * n + j], und zwar:    */
/*                            wenn der i.te Eigenwert reell ist, dann       */
/*                            enthaelt die i.te Spalte von vec den zu-      */
/*                            hoerigen Eigenvektor; wenn der i.te Eigen-    */
/*                            wert komplex ist, dann enthaelt die i.te      */
/*                            Spalte von vec den Realteil und die i + 1.te  */
/*                            Spalte von vec den Imaginaerteil des Eigen-   */
/*                            vektors (und entsprechend fuer den konjugiert */
/*                            komplexen Eigenwert/vektor).                  */
/*      int bal       input:  wenn ungleich 0, wird die Inputmatrix         */
/*                            mit der Funktion balance "balanciert".        */
/*                                                                          */
/*      Return:       0  if successful.                                     */
/*                    1  if exceeded max number of iterations in hqr2.      */
/*                   -1  insufficient memory.                               */

int eigen(int n,double *a,double *evr,double *evi,double *vec,int bal)
{
    int i,j,jj,jn,r,low,hi;
    double *d,tmp,tmp1,vr,vi,vri;                       

    low = 1; hi = n;

    if (!(d = (double *) calloc (n + 1,sizeof(double))))
        return(-1);
    memrq(n + 1,sizeof(double));

    if (bal)    
        balance(n,a,d,&low,&hi);

    orthes(n,a,low,hi,d);   /* reduction to Hessenberg form */
   
    ortrans(n,a,d,low,hi,vec);  /* calculate transformation matrix */

    r = hqr2(n,a,low,hi,vec,evr,evi);   /* calc eigenvalues and vectors */

    if (bal)                        /* back transformation */
        balbak(n,a,d,low,hi);

    /* normalization of eigenvectors */

    i = 1;
    while (i <= n) {
        if (fabs(evi[i]) < EPSI) {    /* if eigenvalue real */
            tmp = 0.0;
            for (j = 1; j <= n; ++j) {
                tmp1 = vec[(j - 1) * n + i];
                tmp += tmp1 * tmp1;
            }
            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                for (j = 1; j <= n; ++j)  
                    vec[(j - 1) * n + i] /= tmp;
            }
            i++;
        }
        else {                      /* complex eigenvalues */
            jj = 0;
            vri = 0.0;
            for (j = 1; j <= n; ++j) {
                jn = (j - 1) * n + i;
                vr = vec[jn];
                vi = vec[jn + 1];
                tmp  = vr * vr + vi * vi;
                if (vri < tmp) {
                    vri = tmp;
                    jj = j;
                }
            }
            if (vri > 0.0) {
                jn = (jj - 1) * n + i;
                vr = vec[jn] / vri;
                vi = vec[jn + 1] / vri;
                for (j = 1; j <= n; ++j) {
                    jn = (j - 1) * n + i;
                    tmp = vec[jn];
                    tmp1 = vec[jn + 1];
                    vec[jn] = tmp * vr + tmp1 * vi;
                    vec[jn + 1] = tmp1 * vr - tmp * vi;
                }
            }
            i += 2;
        }
    }
    free((char *)d);
    memrq(-n - 1,sizeof(double));
    return(r);
}

/* ------------------------------------------------------------------------ */
/*      balance()                                                           */
/*                                                                          */
/*      Eine allgemeine reelle quadratische Matrix wird durch Aehnlich-     */
/*      keitstransformationen so veraendert, dass ihre Norm moeglichst      */
/*      gering ist. Dies ist eine zweckmaessige Behandlung, bevor zum       */
/*      Beispiel Eigenwerte und Eigenvektoren berechnet werden. Zur         */
/*      Ruecktransformation vgl. die Funktion: balbak                       */
/*                                                                          */
/*      Ref.: B.N. Parlett, C. Reinsch, Balancing a Matrix for Calculation  */
/*      of Eigenvalues and Eigenvectors, in: J.H. Wilkinson, C. Reinsch,    */
/*      (eds.), Linear Algebra, Handbook for Automatic Computation, Vol II  */
/*      Springer-Verlag 1971, 315 - 326.                                    */
/*                                                                          */
/*      balance(n,a,d,low,hi)                                               */
/*                                                                          */
/*      int n       input:  dimension of input matrix                       */
/*      double *a   input:  the (n,n) input matrix as a one-dimensional     */
/*                          array a(i,j) = a[(i - 1) * n + j]               */
/*                  output: the transformed matrix                          */
/*      double *d   input:  one-dimensional array d(i) (i = 1,n)            */
/*                  output: information for back transformation.            */
/*      int *low    output: information for back transformation             */
/*      int *hi     output: information for back transformation             */
/*                                                                          */

void balance(int n,double *a,double *d,int *low,int *hi)
{
    register int i,j,k,l;
    int in,jn,kn,ln,noconv;
    double b2,c,g,f,r,s;
    double base = 10.0;     /* base of floating point system */

    b2 = base * base;
    l = 1;
    k = n;
L1:
    kn = (k - 1) * n;

    for (j = k; j >= 1; --j) {
        r = 0.0;
        jn = (j - 1) * n;

        for (i = 1; i <= k; ++i) {
            if (i != j)  
                r += fabs(a[jn + i]);
        }
        if (r == 0.0) {

            d[k] = j; 
            if (j != k) {
                for (i = 1; i <= k; ++i) {
                    in = (i - 1) * n;
                    f = a[in + j];
                    a[in + j] = a[in + k];
                    a[in + k] = f;
                }
                for (i = l; i <= n; ++i) {
                    f = a[jn + i];
                    a[jn + i] = a[kn + i];
                    a[kn + i] = f;
                }
            }
            k--;
            goto L1;
        }
    }

L2:
    ln = (l - 1) * n;

    for (j = l; j <= k; ++j) {
        jn = (j - 1) * n;
        c = 0.0;
        for (i = l; i <= k; ++i) {
            if (i != j)
                c += fabs(a[(i - 1) * n + j]);
        }
        if (c == 0.0) {

            d[l] = j; 
            if (j != l) {
                for (i = 1; i <= k; ++i) {
                    in = (i - 1) * n;
                    f = a[in + j];
                    a[in + j] = a[in + l];
                    a[in + l] = f;
                }
                for (i = l; i <= n; ++i) {
                    f = a[jn + i];
                    a[jn + i] = a[ln + i];
                    a[ln + i] = f;
                }
            }
            l++;
            goto L2;
        }
    }
    *low = l; 
    *hi  = k;

    for (i = l; i <= k; ++i)
        d[i] = 1.0;
    
Iteration:
    noconv = 0;
    for (i = l; i <= k; ++i) {
        in = (i - 1) * n;
        c = r = 0.0;
        for (j = l; j <= k; ++j) {
            if (j != i) {
                c += fabs(a[(j - 1) * n + i]);
                r += fabs(a[in + j]);
            }
        }
        g = r / base;
        f = 1.0;
        s = c + r;
L3:
        if (c < g) {
            f *= base;
            c *= b2;
            goto L3;
        }
        g = r * base;
L4:
        if (c >= g) {
            f /= base;
            c /= b2;
            goto L4;
        }
        if ((c + r) / f < 0.95 * s) {
            g = 1.0 / f;
            d[i] *= f;
            noconv = 1;
            for (j = l; j <= n; ++j)
                a[in + j] *= g;
            for (j = 1; j <= k; ++j)
                a[(j - 1) * n + i] *= f;
        }
    }
    if (noconv)
        goto Iteration;
}

/* ------------------------------------------------------------------------ */
/*      balbak()                                                            */
/*                                                                          */
/*      Ruecktransformation der Rechtseigenvektoren einer Matrix, die       */
/*      bevor ihre Eigenvektoren berechnet worden sind, mit der Funktion    */
/*      balance transformiert worden ist. Beachte: Wenn Linkseigenvektoren  */
/*      ruecktransformiert werden sollen, muss die Funktion veraendert      */
/*      werden; vgl. die im folg. angegebene Quelle.                        */
/*                                                                          */
/*      Ref.: B.N. Parlett, C. Reinsch, Balancing a Matrix for Calculation  */
/*      of Eigenvalues and Eigenvectors, in: J.H. Wilkinson, C. Reinsch,    */
/*      (eds.), Linear Algebra, Handbook for Automatic Computation, Vol II  */
/*      Springer-Verlag 1971, 315 - 326.                                    */
/*                                                                          */
/*      balbak(n,a,d,low,hi)                                                */
/*                                                                          */
/*      int n       input:  dimension of input matrix                       */
/*      double *a   input:  the (n,n) input matrix as a one-dimensional     */
/*                          array a(i,j) = a[(i - 1) * n + j]. the columns  */
/*                          should contain the eigenvectors of a matrix     */
/*                          that was transformed by balance().              */
/*                  output: backtransformed eigenvectors.                   */
/*      double *d   input:  one-dimensional array d(i) (i=1,n) as produced  */
/*                          by the function balance().                      */
/*                  output: not changed.                                    */
/*      int low,hi  input:  values as given by function balance().          */

void balbak(int n,double *a,double *d,int low,int hi)
{
    register int i,j;
    int k,in,kn;
    double tmp;

    for (i = low; i <= hi; ++i) {
        in = (i - 1) * n;
        for (j = 1; j <= n; ++j) 
            a[in + j] *= d[i];
    }
    for (i = 1; i <= n; ++i) {
        if (i < low || i > hi) {
            k = (int)d[i];
            if (k != i) {
                in = (i - 1) * n;
                kn = (k - 1) * n;
                for (j = 1; j <= n; ++j) {
                    tmp = a[in + j];
                    a[in + j] = a[kn + j];
                    a[kn + j] = tmp;
                }
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*      orthes()                                                            */
/*                                                                          */
/*      Eine allgemeine reelle Matrix wird durch orthogonale Transforma-    */
/*      tionen in eine Hessenberg-Form ueberfuehrt.                         */
/*                                                                          */
/*      Ref.: R.S. Martin, J.H. Wilkinson, Similarity Reduction of a        */
/*      General Matrix to Hessenberg-Form, in: J.H. Wilkinson, C. Reinsch,  */
/*      (eds.), Linear Algebra, Handbook for Automatic Computation, Vol II  */
/*      Springer-Verlag 1971, 339 - 358.                                    */
/*                                                                          */
/*      orthes(n,a,k,l,d)                                                   */
/*                                                                          */
/*      int n       input:  dimension of input matrix.                      */
/*      double *a   input:  the (n,n) input matrix as a one-dimensional     */
/*                          array a(i,j) = a[(i - 1) * n + j]               */
/*                  output: der obere Teil des Feldes enthaelt den relev.   */
/*                          Teil der Hessenberg-Matrix, der untere Teil     */
/*                          enthaelt Informationen ueber die verwendeten    */
/*                          orthogonalen Transformationen, wie sie von      */
/*                          der Ruecktransformationsfunktion ortbak be-     */
/*                          noetigt werden.                                 */
/*      int k,l     input:  die Reduktion auf Hessenberg-Form wird nur      */
/*                          fuer die quadr. Submatrix durchgefuehrt, die    */
/*                          bei a(k,k) beginnt und bei a(l,l) endet.        */
/*      double *d   input:  eindimensionales Feld d(i) (i = 1,n)            */
/*                  output: enthaelt weitere Informationen uber die ver-    */
/*                          wendeten orthogonalen Transformationen, wie     */
/*                          sie von der Ruecktransformationsfunktion        */
/*                          ortbak benoetigt werden.                        */

void orthes(int n,double *a,int k,int l,double *d)
{
    register int i,j,m,in;
    double f,g,h;
    double tol;

    tol = DBLMIN / EPSI;  /* DBLMIN smallest number representable   */
                          /* EPSI is machine's epsilon            */

    for (m = k + 1; m < l; ++m) {
        h = 0.0;
        for (i = l; i >= m; --i) {
            f = d[i] = a[(i - 1) * n + m - 1];
            h += f * f;
        }
        if (h <= tol)
            g = 0.0;
        else {
            if (f >= 0.0)
                g = -sqrt(h);
            else
                g = sqrt(h);
            h -= f * g;
            d[m] = f - g;
            for (j = m; j <= n; ++j) {
                f = 0.0;
                for (i = l; i >= m; --i)  
                    f += d[i] * a[(i - 1) * n + j];
                f /= h;
                for (i = m; i <= l; ++i)
                    a[(i - 1) * n + j] -= f * d[i];
            }
            for (i = 1; i <= l; ++i) {
                in = (i - 1) * n;
                f = 0.0;
                for (j = l; j >= m; --j)  
                    f += d[j] * a[in + j];
                f /= h;
                for (j = m; j <= l; ++j)
                    a[in + j] -= f * d[j];
            }
        }
        a[(m - 1) * n + m - 1] = g;
    }
}   

/* ------------------------------------------------------------------------ */
/*      ortrans()                                                           */
/*                                                                          */
/*      Aufbereitung einer durch orthes in Hessenberg-Form gebrachten       */
/*      Matrix fuer die Berechnung der Eigenwerte und Eigenvektoren durch   */
/*      die Funktion hqr2.                                                  */
/*                                                                          */
/*      Ref. G. Peters, J.H. Wilkinson, Eigenvectors of Real and Complex    */
/*      Matrices by LR und QR Triangularizations, in: J.H. Wilkinson,       */
/*      C. Reinsch (eds.), Linear Algebra, Handbook for Automatic           */
/*      Computation, Vol. II, Springer-Verlag 1971, 372 - 395               */
/*                                                                          */
/*      ortrans(n,h,d,low,upp,v)                                            */
/*                                                                          */
/*      int n       input:  Ordnung der Inputmatrix                         */
/*      double *h   input:  die (n,n)-Inputmatrix als eindimensionales      */
/*                          Feld h(i,j) = h[(i - 1) * n + j]. Es muss sich  */
/*                          um eine durch orthes erzeugte Matrix in der     */
/*                          Form einer oberen Hessenbergmatrix handeln.     */
/*                  output: unveraendert                                    */
/*      int low,upp input:  Integerwerte, die von der Funktion balance      */
/*                          erzeugt worden sind, wenn diese vor Anwendung   */
/*                          der Funktion orthes verwendet worden ist;       */
/*                          andernfalls low = 1, upp = n.                   */
/*      double *d   input:  eindimensionales Feld d(i) (i = 1,n), das von   */
/*                          der Funktion orthes erzeugt worden ist.         */
/*                  output: undefiniert                                     */
/*      double *v   input:  Speicherplatz fuer eine (n,n)-Matrix v(i,j) =   */
/*                          v[(i - 1) * n + j]                              */
/*                  output: die Aehnlichkeitstransformationsmatrix, durch   */
/*                          die von orthes die Ursprungsmatrix in die       */
/*                          obere Hessenbergform gebracht worden ist.       */
/*                                                                          */

void ortrans(int n,double *h,double *d,int low,int upp,double *v)
{
    register int i,j,k,m;
    double x,y;

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j)
            v[(i - 1) * n + j] = 0.0;
        v[(i - 1) * n + i] = 1.0;
    }
    for (k = upp - 2; k >= low; --k) {
        m = k + 1;
        y = h[(m - 1) * n + k];
        if (y != 0.0) {
            y *= d[m];
            for (i = k + 2; i <= upp; ++i)
                d[i] = h[(i - 1) * n + k];
            for (j = m; j <= upp; ++j) {
                x = 0.0;
                for (i = m; i <= upp; ++i)  
                    x += d[i] * v[(i - 1) * n + j];
                x /= y;
                for (i = m; i <= upp; ++i) 
                    v[(i - 1) * n + j] += x * d[i];
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*      hqr2()                                                              */
/*                                                                          */
/*      Berechnung der Eigenwerte und Eigenvektoren einer allgemeinen       */
/*      reellen Matrix, die die Form einer oberen Hessenberg-Matrx hat.     */
/*                                                                          */
/*      Ref. G. Peters, J.H. Wilkinson, Eigenvectors of Real and Complex    */
/*      Matrices by LR und QR Triangularizations, in: J.H. Wilkinson,       */
/*      C. Reinsch (eds.), Linear Algebra, Handbook for Automatic           */
/*      Computation, Vol. II, Springer-Verlag 1971, 372 - 395               */
/*                                                                          */
/*      int hqr2(n,h,low,upp,vecs,wr,wi)                                    */
/*                                                                          */
/*      int n           input:  Ordnung der Matrix                          */
/*      double *h       input:  die (n,n)-Matrix h als ein eindimension.    */
/*                              Feld h(i,j) = h[(i - 1) * n + j].           */
/*                      output: zerstoert (enthaelt die Eigenwerte in der   */
/*                              Diagonalen bzw. der Subdiagonalen)          */
/*      int low,upp     input:  wenn die Inputmatrix zuvor mit der Funktion */
/*                              balance transformiert worden ist, werden    */
/*                              diese Integerwerte von balance erzeugt;     */
/*                              andernfalls low = 1 und upp = n.            */
/*      double *vecs    input:  ein eindim. Feld fuer eine (n,n)-Matrix     */
/*                              in der Form vecs(i,j) = vecs[(i-1) * n + j] */
/*                      output: die Eigenvektoren als Spaltenvektoren. Wenn */
/*                              der i.te Eigenwert reell, dann befindet     */
/*                              sich der zugehoerige Eigenvektor in der     */
/*                              i.ten Spalte; wenn der i.te Eigenwert ein   */
/*                              konj. komplexes Paar, dann befindet sich    */
/*                              der Realteil des Eigenvektors in der i.ten, */
/*                              der Imaginaerteil in der i + 1.ten Spalte   */
/*      double *wr      input:  eindim. Feld wr(i) (i = 1,n)                */
/*                      output: die Realteile der Eigenwerte                */
/*      double *wi      input:  eindim. Feld wi(i) (i = 1,n)                */
/*                      output: die Imaginaerteile der Eigenwerte           */
/*                                                                          */
/*      Return:         0  if successfull                                   */
/*                      1  if exceeded the max number of iterations.        */

int hqr2(int n,double *h,int low,int upp,double *vecs,double *wr,double *wi)
{
    register int i,j,k,l;
    int m,en,na,its;
    double norm,p,q,r,ra,s,sa,t,vr,vi,w,x,y,z,tmp1,tmp2;
    int itmax = 50;   /* max number of iterations */

    for (i = 1; i < low; ++i) {
        wr[i] = h[(i - 1) * n + i];
        wi[i] = 0.0;

        /**********
        cnt[i] = 0;
        **********/
    }
    for (i = upp + 1; i <= n; ++i) {
        wr[i] = h[(i - 1) * n + i];
        wi[i] = 0.0;

        /***********
        cnt[i] = 0;
        ***********/
    }
    en = upp;
    t = 0.0;

Nextw:
    if (en < low)
        goto Fin;

    its = 0;
    na = en - 1;

NextIt:
    for (l = en; l > low; --l) {
        tmp1 = fabs(h[(l - 1) * n + l - 1]);
        tmp2 = fabs(h[(l - 2) * n + l - 1]) + fabs(h[(l - 1) * n + l]);
        tmp2 *= EPSI;
        if (tmp1 <= tmp2)
            goto Cont1;
    }
    l = low;

Cont1:
    x = h[(en - 1) * n + en];
    if (l == en)
        goto Onew;
    y = h[(na - 1) * n + na];
    w = h[(en - 1) * n + na] * h[(na - 1) * n + en];
    if (l == na)
        goto Twow;

    if (its >= itmax) {

        /**************
        cnt[en] = 31;
        **************/

        return(1);
    }
    if (its == 10 || its == 20) {
        t += x;
        for (i = low; i <= en; ++i) 
            h[(i - 1) * n + i] -= x;
        s = fabs(h[(en - 1) * n + na]) + fabs(h[(na - 1) * n + en - 2]);
        x = y = 0.75 * s;
        w = -0.4375 * s * s;
    }
    its++;

    for (m = en - 2; m >= l; --m) {

        z = h[(m - 1) * n + m];
        r = x - z;
        s = y - z;
        p = (r * s - w) / h[m * n + m] + h[(m - 1) * n + m + 1];
        q = h[m * n + m + 1] - z - r - s;
        r = h[(m + 1) * n + m + 1];
        s = fabs(p) + fabs(q) + fabs(r);
        p /= s;
        q /= s;
        r /= s;
        if (m == l)
            goto Cont2;
        tmp1 = fabs(h[(m - 1) * n + m - 1]) * (fabs(q) + fabs(r));
        tmp2 = fabs(h[(m - 2) * n + m - 1]) + fabs(z) + fabs(h[m * n + m + 1]);
        tmp2 *= fabs(p) * EPSI;
        if (tmp1 <= tmp2)
            goto Cont2;
    }

Cont2:
    for (i = m + 2; i <= en; ++i)
        h[(i - 1) * n + i - 2] = 0.0;
    for (i = m + 3; i <= en; ++i)
        h[(i - 1) * n + i - 3] = 0.0;

    for (k = m; k <= na; ++k) {
        if (k != m) {
            p = h[(k - 1) * n + k - 1];
            q = h[k * n + k - 1];
            if (k != na)
                r = h[(k + 1) * n + k - 1];
            else
                r = 0.0;
            x = fabs(p) + fabs(q) + fabs(r);
            if (x == 0.0)
                goto Cont3;
            p /= x;   
            q /= x;
            r /= x;
        }
        s = sqrt(p * p + q * q + r * r);
        if (p < 0.0)
            s = -s;
        if (k != m)
            h[(k - 1) * n + k - 1] = -s * x;
        else if (l != m)
            h[(k - 1) * n + k - 1] = -h[(k - 1) * n + k - 1];
        p += s;
        x = p / s;
        y = q / s;
        z = r / s;
        q /= p;
        r /= p;

        for (j = k; j <= n; ++j) {
            p = h[(k - 1) * n + j] + q * h[k * n + j];
            if (k != na) {
                p += r * h[(k + 1) * n + j];
                h[(k + 1) * n + j] -= p * z;
            }
            h[k * n + j] -= p * y;
            h[(k - 1) * n + j] -= p * x;
        }
        j = k + 3;
        if (j > en)
            j = en;
        for (i = 1; i <= j; ++i) {
            p = x * h[(i - 1) * n + k] + y * h[(i - 1) * n + k + 1];
            if (k != na) {
                p += z * h[(i - 1) * n + k + 2];
                h[(i - 1) * n + k + 2] -= p * r;
            }
            h[(i - 1) * n + k + 1] -= p * q;
            h[(i - 1) * n + k] -= p;
        }
        for (i = low; i <= upp; ++i) {
            p = x * vecs[(i - 1) * n + k] + y * vecs[(i - 1) * n + k + 1];
            if (k != na) {
                p += z * vecs[(i - 1) * n + k + 2];
                vecs[(i - 1) * n + k + 2] -= p * r;
            }
            vecs[(i - 1) * n + k + 1] -= p * q;
            vecs[(i - 1) * n + k] -= p;
        }

Cont3:  ;
    }
    goto NextIt;

Onew:
    wr[en] = h[(en - 1) * n + en] = x + t;
    wi[en] = 0.0;

    /**************
    cnt[en] = its;
    **************/

    en = na;
    goto Nextw;

Twow:
    p = (y - x) / 2.0;
    q = p * p + w;
    z = sqrt(fabs(q));
    x = h[(en - 1) * n + en] = x + t;
    h[(na - 1) * n + na] = y + t;

    /***************
    cnt[en] = -its;
    cnt[na] = its;
    ***************/

    if (q > 0.0) {
        if (p < 0.0)
            z = p - z;
        else
            z = p + z;
        wr[na] = x + z;
        wr[en] = s = x - w / z;
        wi[na] = wi[en] = 0.0;
        x = h[(en - 1) * n + na];
        r = sqrt(x * x + z * z);
        p = x / r;
        q = z / r;
        for (j = na; j <= n; ++j) {
            z = h[(na - 1) * n + j];
            h[(na - 1) * n + j] = q * z + p * h[(en - 1) * n + j];
            h[(en - 1) * n + j] = q * h[(en - 1) * n + j] - p * z;
        }
        for (i = 1; i <= en; ++i) {
            z = h[(i - 1) * n + na];
            h[(i - 1) * n + na] = q * z + p * h[(i - 1) * n + en];
            h[(i - 1) * n + en] = q * h[(i - 1) * n + en] - p * z;
        }
        for (i = low; i <= upp; ++i) {
            z = vecs[(i - 1) * n + na];
            vecs[(i - 1) * n + na] = q * z + p * vecs[(i - 1) * n + en];
            vecs[(i - 1) * n + en] = q * vecs[(i - 1) * n + en] - p * z;
        }
    }
    else {
        wr[na] = wr[en] = x + p;
        wi[na] = z;
        wi[en] = -z;
    }
    en -= 2;
    goto Nextw;

Fin:
    norm = 0.0;
    k = 1;
    for (i = 1; i <= n; ++i) {
        for (j = k; j <= n; ++j)
            norm += fabs(h[(i - 1) * n + j]);
        k = i;
    }
    for (en = n; en >= 1; --en) {
        p = wr[en];
        q = wi[en];
        na = en - 1;
        if (q == 0.0) {
            m = en;
            h[(en - 1) * n + en] = 1.0;
            for (i = na; i >= 1; --i) {
                w = h[(i - 1) * n + i] - p;
                r = h[(i - 1) * n + en];
                for (j = m; j <= na; ++j)  
                    r += h[(i - 1) * n + j] * h[(j - 1) * n + en];
                if (wi[i] < 0.0) {
                    z = w;
                    s = r;
                }
                else {
                    m = i;
                    if (wi[i] == 0.0) {
                        if (w != 0.0)
                            h[(i - 1) * n + en] = -r / w;
                        else
                            h[(i - 1) * n + en] = -r / (norm * EPSI);
                    }
                    else {
                        x = h[(i - 1) * n + i + 1];
                        y = h[i * n + i];
                        tmp1 = wr[i] - p;
                        q = tmp1 * tmp1 + wi[i] * wi[i];
                        h[(i - 1) * n + en] = t = (x * s - z * r) / q;
                        if (fabs(x) > fabs(z))
                            h[i * n + en] = (-r - w * t) / x;
                        else
                            h[i * n + en] = (-s - y * t) / z;
                    }
                }
            }
        }
        else {
            if (q < 0.0) {
                m = na;
                if (fabs(h[(en - 1) * n + na]) > fabs(h[(na - 1) * n + en])) {
                    h[(na - 1) * n + na] = -(h[(en - 1) * n + en] - p) /
                                                         h[(en - 1) * n + na];
                    h[(na - 1) * n + en] = -q / h[(en - 1) * n + na];
                }
                else {
                    cdiv(-h[(na - 1) * n + en],0.0,h[(na - 1) * n + na] - p,q,
                                &h[(na - 1) * n + na],&h[(na - 1) * n + en]);
                }
                h[(en - 1) * n + na] = 1.0;
                h[(en - 1) * n + en] = 0.0;

                for (i = na - 1; i >= 1; --i) {
                    w = h[(i - 1) * n + i] - p;
                    ra = h[(i - 1) * n + en];
                    sa = 0.0;
                    for (j = m; j <= na; ++j) {
                        ra += h[(i - 1) * n + j] * h[(j - 1) * n + na];
                        sa += h[(i - 1) * n + j] * h[(j - 1) * n + en];
                    }
                    if (wi[i] < 0.0) {
                        z = w;
                        r = ra;
                        s = sa;
                    }
                    else {
                        m = i;
                        if (wi[i] == 0.0) {
                            cdiv(-ra,-sa,w,q,&h[(i - 1) * n + na],
                                             &h[(i - 1) * n + en]);
                        }      
                        else {
                            x = h[(i - 1) * n + i + 1];
                            y = h[i * n + i];
                            tmp1 = wr[i] - p;
                            vr = tmp1 * tmp1 + wi[i] * wi[i] - q * q;
                            vi = tmp1 * 2.0 * q;
                            if (vr == 0.0 && vi == 0.0) {
                                vr = fabs(w) + fabs(q) + fabs(x) +
                                     fabs(y) + fabs(z);
                                vr *= norm * EPSI;
                            }
                            cdiv(x * r - z * ra + q * sa,x * s - z * sa - q * ra,
                              vr,vi,&h[(i - 1) * n + na],&h[(i - 1) * n + en]);

                            if (fabs(x) > fabs(z) + fabs(q)) {
                                h[i * n + na] = (-ra - w * h[(i - 1) * n + na] +
                                                  q * h[(i - 1) * n + en]) / x;
                                h[i * n + en] = (-sa - w * h[(i - 1) * n + en] -
                                                  q * h[(i - 1) * n + na]) / x;
                            }
                            else {
                                cdiv(-r - y * h[(i - 1) * n + na],
                                     -s - y * h[(i - 1) * n + en],z,q,
                                           &h[i * n + na],&h[i * n + en]);
                            }
                        }
                    }
                }
            }
        }
    }

    for (i = 1; i < low; ++i) {
        for (j = i + 1; j <= n; ++j)
            vecs[(i - 1) * n + j] = h[(i - 1) * n + j];
    }
    for (i = upp + 1; i <= n; ++i) {
        for (j = i + 1; j <= n; ++j)
            vecs[(i - 1) * n + j] = h[(i - 1) * n + j];
    }
    for (j = n; j >= low; --j) {
        m = j;
        if (m > upp)
            m = upp;
        l = j - 1;
        if (wi[j] < 0.0) {
            for (i = low; i <= upp; ++i) {
                y = z = 0.0;
                for (k = low; k <= m; ++k) {
                    y += vecs[(i - 1) * n + k] * h[(k - 1) * n + l];
                    z += vecs[(i - 1) * n + k] * h[(k - 1) * n + j];
                }
                vecs[(i - 1) * n + l] = y;
                vecs[(i - 1) * n + j] = z;
            }
        }
        else {
            if (wi[j] == 0.0) {
                for (i = low; i <= upp; ++i) {
                    z = 0.0;
                    for (k = low; k <= m; ++k)  
                        z += vecs[(i - 1) * n + k] * h[(k - 1) * n + j];
                    vecs[(i - 1) * n + j] = z;
                }
            }
        }
    }
    return(0);
}

/* Function: cdiv (complex division) */

void cdiv(double xr,double xi,double yr,double yi,double *zr,double *zi)
{
    double h;          

    if (fabs(yr) > fabs(yi)) {
        h = yi / yr;
        yr = h * yi + yr;
        *zr = (xr + h * xi) / yr;
        *zi = (xi - h * xr) / yr;
    }
    else {
        h = yr / yi;
        yi = h * yr + yi;
        *zr = (h * xr + xi) / yi;
        *zi = (h * xi - xr) / yi;
    }
}

/* ------------------------------------------------------------------------ */
/*      svdecomp()                                                          */
/*                                                                          */
/*      Singular Value Decomposition einer (m,n)-Matrix a (m >= n)          */
/*                                                                          */
/*      Diese Funktion berechnet:                                           */
/*                                                                          */
/*      eine orthogonale (m,n)-Matrix U  (U'U = In)                         */
/*      eine orthogonale (n,n)-Matrix V  (V'V = VV' = In)                   */
/*      eine Diagonalmatrix Q = { q1,...,qn }                               */
/*                                                                          */
/*      so dass: A = U * Q * V' (singular value decomposition).             */
/*                                                                          */
/*      Die Werte q1,...,qn sind die Singular Values von A. Es gilt, dass   */
/*      q1 * q1,..., qn * qn die Eigenwerte der symmetrischen positiv       */
/*      definiten (n,n)-Matrix A'A und zugleich die n groessten Eigenwerte  */
/*      der symmetrischen pos. definiten (m,m)-Matrix AA' sind.             */
/*                                                                          */
/*      Die n Spaltenvektoren der Matrix V sind die zugehoerigen Eigen-     */
/*      vektoren von A'A. Die n Spaltenvektoren der Matrix U sind die zu    */
/*      den n groessten Eigenwerten von AA' zugehoerigen Eigenvektoren.     */
/*                                                                          */
/*      Es gilt: Wenn r = rang(A), genau dann sind genau r der n Singular   */
/*      Values von A ungleich 0. Die Anzahl der Singular Values = 0 zeigt   */
/*      also den Defekt von A an.                                           */
/*                                                                          */
/*      Ref.: G.H. Golub, C. Reinsch, Singular Value Decomposition and      */
/*      Least Squares Solutions, in: J.H. Wilkinson, C. Reinsch, Linear     */
/*      Algebra (= Handbook for Automatic Computation, Vol. II), Berlin     */
/*      (Springer) 1971.                                                    */
/*                                                                          */
/*          svdecomp(m,n,a,q,u,v,mode)                                      */
/*                                                                          */
/*          int m       input:  Anzahl Zeilen der Matrix a                  */
/*          int n       input:  Anzahl Spalten der Matrix a (m >= n)        */
/*          double *a   input:  Die Matrix a als eindimensionales Feld      */
/*          double *q   input:  Eindimensionales Feld: 1,...,n              */
/*                      output: Die Singular Values der Matrix a            */
/*          double *u   input:  Eindimensionales Feld: 1,...,mn             */
/*                      output: Die (m,n)-Matrix u                          */
/*          double *v   input:  Eindimensionales Feld: 1,...,nn             */
/*                      output: Die (n,n)-Matrix v                          */
/*          int mode    input:  0 Berechnung nur der Singular Values        */
/*                              1 Ausserdem wird die u-Matrix berechnet     */
/*                              2 Ausserdem wird die v-Matrix berechnet     */
/*                              3 Es werden v- und u-Matrix berechnet       */
/*                                                                          */
/*          return: 0 if successful                                         */
/*                 -1 if insufficient memory                                */
/*                 -2 if m < n                                              */

int svdecomp(int m,int n,double *a,double *q,double *u,double *v,int mode)
{
    register int i,j,k,l,l1,l2;
    double c,f,f1,g,h,s,x,y,z,eps,tol;
    double *e;

    if (m < n) 
        return(-2);

    /*********************************************************************
    Bestimmung der Konstanten eps und tol. eps darf nicht kleiner sein als
    die kleine Zahl eta, fuer die gilt 1. + eta > 1. Sie wird als 2. * EPSI
    definiert.

    tol ist definert als DBLMIN / EPSI, wobei DBLMIN kleinste positive Zahl
    ist, die dargestellt werden kann.
    *********************************************************************/

    eps = 2.0 * EPSI;
    tol = DBLMIN / EPSI;

    if (!(e = (double *) calloc (n + 1,sizeof(double))))
        return(-1);
    memrq(n + 1,sizeof(double));

    for (i = 1; i <= m; ++i)
        for (j = 1; j <= n; ++j) 
            u[(i - 1) * n + j] = a[(i - 1) * n + j];

    x = g = 0.0;
    for (i = 1; i <= n; ++i) {
        e[i] = g; s = 0.0; l = i + 1;
        for (j = i; j <= m; ++j)
            s += u[(j - 1) * n + i] * u[(j - 1) * n + i];
        g = 0.0;
        if (s >= tol) {
            f = u[(i - 1) * n + i];
            g = sqrt(s);
            if (f >= 0.0) g = -g;
            h = f * g - s;
            u[(i - 1) * n + i] = f - g;
            for (j = l; j <= n; ++j) {
                s = 0.0;
                for (k = i; k <= m; ++k)
                    s += u[(k - 1) * n + i] * u[(k - 1) * n + j];
                if (h)
                    f = s / h;
                for (k = i; k <= m; ++k)
                    u[(k - 1) * n + j] += f * u[(k - 1) * n + i];
            }
        }
        q[i] = g; s = 0.0;
        for (j = l; j <= n; ++j) 
            s += u[(i - 1) * n + j] * u[(i - 1) * n + j];
        g = 0.0;
        if (s >= tol) {
            if (i != n) {
                f = u[(i - 1) * n + i + 1]; g = sqrt(s); if (f >= 0.0) g = -g;
                h = f * g - s;
                u[(i - 1) * n + i + 1] = f - g;
                if (h) {
                    for (j = l; j <= n; ++j) 
                        e[j] = u[(i - 1) * n + j] / h;
                }
            }
            for (j = l; j <= m; ++j) {
                s = 0.0;
                for (k = l; k <= n; ++k)
                    s += u[(j - 1) * n + k] * u[(i - 1) * n + k];
                for (k = l; k <= n; ++k)
                    u[(j - 1) * n + k] += s * e[k];
            }
        }
        y = fabs(q[i]) + fabs(e[i]); if (y > x) x = y;
    }
    if (mode == 2 || mode == 3) {
        i = n;
        while (1) {
            v[(i - 1) * n + i] = 1.0; g = e[i]; l = i;
            if (--i == 0) break;
            if (g != 0.0) {
                h = u[(i - 1) * n + i + 1] * g;
                if (h) {
                    for (j = l; j <= n; ++j)
                        v[(j - 1) * n + i] = u[(i - 1) * n + j] / h;
                }
                for (j = l; j <= n; ++j) {
                    s = 0.0;
                    for (k = l; k <= n; ++k) 
                        s += u[(i - 1) * n + k] * v[(k - 1) * n + j];
                    for (k = l; k <= n; ++k) 
                        v[(k - 1) * n + j] += s * v[(k - 1) * n + i];
                }
            }
            for (j = l; j <= n; ++j)  
                v[(i - 1) * n + j] = v[(j - 1) * n + i] = 0.0;
        }
    }
    if (mode == 1 || mode == 3) {
        i = n;
        while (1) {
            l = i + 1; g = q[i];
            for (j = l; j <= n; ++j)
                u[(i - 1) * n + j] = 0.0;
            if (g != 0.0) {
                h = u[(i - 1) * n + i] * g;
                for (j = l; j <= n; ++j) {
                    s = 0.0;
                    for (k = l; k <= m; ++k) 
                        s += u[(k - 1) * n + i] * u[(k - 1) * n + j];
                    if (h)
                        f = s / h;
                    for (k = i; k <= m; ++k) 
                        u[(k - 1) * n + j] += f * u[(k - 1) * n + i];
                }
                if (g) {
                    for (j = i; j <= m; ++j)  
                        u[(j - 1) * n + i] /= g;
                }
            }
            else {
                for (j = i; j <= m; ++j)
                    u[(j - 1) * n + i] = 0.0;
            }
            u[(i - 1) * n + i] += 1.0;
            if (--i == 0) break;
        }
    }
    eps *= x; k = n;

    while (1) {                                          /* main loop    */

        for (l = k; l >= 1; --l) {
            if (fabs(e[l]) <= eps)
                goto lab37;
            if (l > 1 && fabs(q[l - 1]) <= eps )
                break;
        }

        c = 0.0; s = 1.0; l1 = l - 1;
        for (i = l; i <= k; ++i) {
            f = s * e[i]; e[i] *= c;
            if (fabs(f) <= eps) break;
            g = q[i]; q[i] = sqrt(f * f + g * g);
            h = q[i];
            if (h) {
                c = g / h;
                s = -f / h;
            }
            if (mode == 1 || mode == 3) {
                for (j = 1; j <= m; ++j) {
                    y = 0.0;
                    if (l1 >= 1) y = u[(j - 1) * n + l1];
                    z = u[(j - 1) * n + i];
                    if (l1 >= 1) u[(j - 1) * n + l1] = y * c + z * s;
                    u[(j - 1) * n + i] = -y * s + z * c;
                }
            }
        }
lab37:
        z = q[k];
        if (l == k) {
            if (z < 0.0) {
                q[k] = -z;
                if (mode == 2 || mode == 3) {
                    for (j = 1; j <= n; ++j)
                        v[(j - 1) * n + k] = -v[(j - 1) * n + k];
                }
            }
            if (--k == 0) break;                   /* break main loop      */
        }
        else {
            x = q[l];
            y = q[k - 1];
            g = e[k - 1];
            h = e[k];
            if (h && y)
                f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = sqrt(f * f + 1.0);
            f1 = f - g;
            if (f >= 0.0) f1 = f + g;
            if (x && f1)
                f = ((x - z) * (x + z) + h * (y / f1 - h)) / x;
            c = s = 1.0; l2 = l + 1;
            for (i = l2; i <= k; ++i) {
                g = e[i]; y = q[i]; h = s * g; g = c * g;
                z = sqrt(f * f + h * h);
                e[i - 1] = z;
                if (z) {
                    c = f / z;
                    s = h / z;
                }
                f = x * c + g * s; g = -x * s + g * c;
                h = y * s; y = y * c;

                if (mode == 2 || mode == 3) {
                    for (j = 1; j <= n; ++j) {
                        x = v[(j - 1) * n + i - 1];
                        z = v[(j - 1) * n + i];
                        v[(j - 1) * n + i - 1] = x * c + z * s;
                        v[(j - 1) * n + i] = -x * s + z * c;
                    }  
                }     
                z = sqrt(f * f + h * h);
                q[i - 1] = z;
                if (z) {
                    c = f / z;
                    s = h / z;
                }
                f = c * g + s * y;
                x = -s * g + c * y;
      
                if (mode == 1 || mode == 3) {
                    for (j = 1; j <= m; ++j) {
                        y = u[(j - 1) * n + i - 1];
                        z = u[(j - 1) * n + i];
                        u[(j - 1) * n + i - 1] = y * c + z * s;
                        u[(j - 1) * n + i] = -y * s + z * c;
                    }  
                }
            }
            e[l] = 0.0; e[k] = f; q[k] = x;
        }
    }
    free((char *)e);
    memrq(-n - 1,sizeof(double));

    /* sort in descending order */

    for (i = 1; i < n; ++i) {
        c = q[i];        
        k = i;
        for (j = i + 1; j <= n; ++j) {
            if (q[j] > c) {
                k = j;
                c = q[j];
            }
        }
        if (k != i) {
            q[k] = q[i];
            q[i] = c;
            if (mode == 1 || mode == 3) {
                for (j = 0; j < m; ++j) {
                    x = u[j * n + i];
                    u[j * n + i] = u[j * n + k];
                    u[j * n + k] = x;
                }
            }
            if (mode == 2 || mode == 3) {
                for (j = 0; j < n; ++j) {
                    x = v[j * n + i];
                    v[j * n + i] = v[j * n + k];
                    v[j * n + k] = x;
                }
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*      ginv()                                                              */
/*                                                                          */
/*      Berechnung der verallgemeinerten Inversen einer (m,n)-Matrix unter  */
/*      Verwendung der Funktion svdecomp(). Dabei wird vorausgesetzt: m>=n  */
/*                                                                          */
/*      Sei a eine (m,n)-Matrix (m >= n). Dann gibt es genau eine (n,m)-    */
/*      Matrix b, so dass gilt:                                             */
/*                                                                          */
/*      (1)  a b a = a                                                      */
/*      (2)  b a b = b                                                      */
/*      (3)  (ab)' = a b                                                    */
/*      (4)  (ba)' = b a                                                    */
/*                                                                          */
/*      Die Matrix b heisst verallgemeinerte Inverse von a. Wenn a eine     */
/*      regulaere quadratische Matrix ist, ist b mit der gewoehnlichen      */
/*      Inversen von a identisch.                                           */
/*                                                                          */
/*      Zur Berechnung von b wird eine Singular Value Decomposition von a   */
/*      benutzt: a = u q v'  (vgl. die Funktion svdecomp). Es gilt:         */
/*      b = v p u', wobei p die verallgemeinerte Inverse von q ist. Die     */
/*      verallgemeinerte Inverse von q ist folgendermassen definiert:       */
/*                                                                          */
/*      p(i,i) := 1./q(i,i) wenn q(i,i) != 0, 0 sonst: 0                    */
/*                                                                          */
/*      Ref.: G.H. Golub, C. Reinsch, Singular Value Decomposition and      */
/*      Least Squares Solutions, in: J.H. Wilkinson, C. Reinsch, Linear     */
/*      Algebra (= Handbook for Automatic Computation, Vol. II), Berlin     */
/*      (Springer) 1971.                                                    */
/*                                                                          */
/*          ginv(m,n,a)                                                     */
/*                                                                          */
/*          int m       input:  Anzahl Zeilen der Matrix (m >= n)           */
/*          int n       input:  Anzahl Spalten der Matrix                   */
/*          double *a   input:  Eindimens. Feld mit der (m,n)-Matrix a      */
/*                              a(i,j) = a[ (i - 1) * n + j ]               */
/*                      output: Die verallgemeinerte Inverse (n,m)-Matrix   */
/*                                                                          */
/*          Return: wenn >= 0 den Pseudorank von a                          */
/*                  andernfalls einen Fehlercode:                           */
/*                  -1 if insufficient memory.                              */
/*                  -2 if m < n                                             */
/*                                                                          */

int ginv(int m,int n,double *a)
{
    register int i,j,k,r;
    double tmp,*q,*u,*v;

    if (m < n)
        return (-2);

    if (!(q = (double *) calloc (n + 1,sizeof(double)))) {
        return(-1);
    }
    if (!(u = (double *) calloc (m * n + 1,sizeof(double)))) {
        free((char *)q);
        return(-1);
    }
    if (!(v = (double *) calloc (n * n + 1,sizeof(double)))) {
        free((char *)q);
        free((char *)u);
        return(-1);
    }
    memrq(n * (1 + n + m) + 3,sizeof(double));

    if ((r = svdecomp(m,n,a,q,u,v,3)) < 0)    
        goto GINV_Fin;
           
    r = 0;
    for (i = 1; i <= n; ++i) {
        if (fabs(q[i]) > EPSI1)
            r++;
        for (j = 1; j <= m; ++j) {
            tmp = 0.0;
            for (k = 1; k <= n; ++k) {
                if (fabs(q[k]) > EPSI1)
                    tmp += v[(i - 1) * n + k] * u[(j - 1) * n + k] / q[k];
            }
            a[(i - 1) * m + j] = tmp;
        }
    }

GINV_Fin:
    free((char *)v);
    free((char *)u);
    free((char *)q);
    memrq(-n * (1 + n + m) - 3,sizeof(double));
    return (r);
}
   
/* ------------------------------------------------------------------------ */
/*      ginv1()                                                             */
/*                                                                          */
/*      Identical to ginv() but, in addition, returns singular values.      */
/*                                                                          */
/*          ginv1(m,n,a,sv)                                                 */
/*                                                                          */
/*          int m       input:  Anzahl Zeilen der Matrix (m >= n)           */
/*          int n       input:  Anzahl Spalten der Matrix                   */
/*          double *a   input:  Eindimens. Feld mit der (m,n)-Matrix a      */
/*                              a(i,j) = a[ (i - 1) * n + j ]               */
/*                      output: Die verallgemeinerte Inverse (n,m)-Matrix   */
/*          double *sv  input:  array of dimension n                        */
/*                      output: singular values                             */
/*                                                                          */
/*          Return: wenn >= 0 den Pseudorank von a                          */
/*                  andernfalls einen Fehlercode:                           */
/*                  -1 if insufficient memory.                              */
/*                  -2 if m < n                                             */
/*                                                                          */

int ginv1(int m,int n,double *a,double *sv)
{
    register int i,j,k,r;
    double tmp,*u,*v;

    if (m < n)
        return (-2);

    if (!(u = (double *) calloc (m * n + 1,sizeof(double)))) {
        return(-1);
    }
    if (!(v = (double *) calloc (n * n + 1,sizeof(double)))) {
        free((char *)u);
        return(-1);
    }
    memrq(n * (n + m) + 2,sizeof(double));

    if ((r = svdecomp(m,n,a,sv,u,v,3)) < 0)    
        goto GINV1Fin;
           
    r = 0;
    for (i = 1; i <= n; ++i) {
        if (sv[i] != 0.0)
            r++;
        for (j = 1; j <= m; ++j) {
            tmp = 0.0;
            for (k = 1; k <= n; ++k) {
                if (sv[k] != 0.0)
                    tmp += v[(i - 1) * n + k] * u[(j - 1) * n + k] / sv[k];
            }
            a[(i - 1) * m + j] = tmp;
        }
    }

GINV1Fin:
    free((char *)v);
    free((char *)u);
    memrq(-n * (n + m) - 2,sizeof(double));
    return (r);
}

/* ------------------------------------------------------------------------ */
/*  evecf                                                                   */
/*      Calculate eigenvalues and eigenvectors of a real symmetric matrix.  */
/*      Code is adapted from:                                               */
/*                                                                          */
/*      D.N. Sparks, A.D. Todd, Latent Roots of a Symmetric Matrix,         */
/*      AS 60, Applied Statistics 22 (1973), 260 - 265. Comment in:         */
/*      Applied Statistics 23 (1974), 101 - 102.                            */
/*                                                                          */
/*      Call: int evecf(n,d,z)                                              */
/*                                                                          */
/*      n is the dimension of the input matrix z, d is a n-vector, z is     */
/*      a (n,n)-matrix. If successful at return d contains the eigenvalues, */
/*      and z contains the associated eigenvectors. The eigenvector are     */
/*      normalized to an Euklidean norm of 1. Note that only the lower      */
/*      triangle of the input matrix z is needed for input.                 */
/*                                                                          */
/*      The eigenvalues are in descending order. The resulting matrix z     */
/*      of eigenvectors is scaled so that it is orthogonal: z z' = I.       */
/*                                                                          */
/*      Return -1 if not successful, that is the maximum number of          */
/*                iterations is reached.                                    */
/*             -2 if insufficient memory.                                   */

int evecf(int n, double *d, double *z)
{
    int r; 
    double *e;

    if (!(e = (double *) calloc (n + 1,sizeof(double))))  
        return(-2);

    tridf1(n,d,e,z);
    r = evectf(n,d,e,z);

    free((char *) e);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  tridf1  Tridiagonalization of a real symmetric matrix, used by evecf.   */

void tridf1(int n, double *d, double *e, double *z)
{
    register int i,i1,j,k;
    int l,in,j2;
    double f,g,h,hh,tol;

    tol = DBLMIN / EPSI;

    i = n;
    for (i1 = 2; i1 <= n; ++i1) {
        in = (i - 1) * n;
        l = i - 2;
        f = z[in + i - 1];
        g = 0.0;
        for (k = 1; k <= l; ++k) 
            g += z[in + k] * z[in + k];
        h = g + f * f;
        if (g <= tol) {
            e[i] = f;
            d[i] = 0.0;
        }
        else {
            l++;
            g = sqrt(h);
            if (f >= 0.0)
                g = -g;
            e[i] = g;
            h -= f * g;
            z[in + i - 1] = f - g;
            f = 0.0;
            for (j = 1; j <= l; ++j) {
                j2 = (j - 1) * n;
                z[j2 + i] = z[in + j] / h;
                g = 0.0;
                for (k = 1; k <= j; ++k)
                    g += z[j2 + k] * z[in + k];
                for (k = j + 1; k <= l; ++k)
                    g += z[(k - 1) * n + j] * z[in + k];
                e[j] = g / h;
                f += g * z[j2 + i];
            }
            hh = f / (h + h);
            for (j = 1; j <= l; ++j) {
                f = z[in + j];
                g = e[j] - hh * f;
                e[j] = g;
                for (k = 1; k <= j; ++k) 
                    z[(j - 1) * n + k] -= (f * e[k] + g * z[in + k]);
            }
            d[i] = h;
        }
        i--;
    }
    d[1] = e[1] = 0.0;

    for (i = 1; i <= n; ++i) {
        in = (i - 1) * n;
        l = i - 1;
        if (d[i] != 0.0 && l > 0) {
            for (j = 1; j <= l; ++j) {
                g = 0.0;
                for (k = 1; k <= l; ++k)  
                    g += z[in + k] * z[(k - 1) * n + j];
                for (k = 1; k <= l; ++k)
                    z[(k - 1) * n + j] -= g * z[(k - 1) * n + i];
            }
        }
        d[i] = z[in + i];
        z[in + i] = 1.0;
        for (j = 1; j <= l; ++j)
            z[in + j] = z[(j - 1) * n + i] = 0.0;
    }   
}

/* ------------------------------------------------------------------------ */
/*  evectf  Calculate eigenvalues. Used by evecf.                           */

int evectf(int n, double *d, double *e, double *z)
{
    register int i,j,k,l;
    int i1,jj,m,m1,kn,j2;
    double b,c,f,g,h,p,pr,r,s;
    int maxits = 50; /* maximum number of iterations */

    for (i = 2; i <= n; ++i)
        e[i - 1] = e[i];
    e[n] = b = f = 0.0;
    for (l = 1; l <= n; ++l) {
        jj = 0;
        h = EPSI * (fabs(d[l]) + fabs(e[l]));
        if (b < h)
            b = h;
        for (m1 = l; m1 <= n; ++m1) {
            m = m1;
            if (fabs(e[m]) <= b)
                break;
        }
        if (m != l) {
            while (++jj <= maxits) {
                p = (d[l + 1] - d[l]) / (2.0 * e[l]);
                r = sqrt(p * p + 1.0);
                if (p < 0.0)
                    pr = p - r;
                else
                    pr = p + r;
                h = d[l] - e[l] / pr;
                for (i = l; i <= n; ++i)
                    d[i] -= h;
                f += h;
                p = d[m];
                c = 1.0;
                s = 0.0;
                m1 = m - 1;
                i = m;
                for (i1 = l; i1 <= m1; ++i1) {
                    j = i;
                    i--;
                    g = c * e[i];
                    h = c * p;
                    if (fabs(p) < fabs(e[i])) {
                        c = p / e[i];
                        r = sqrt(c * c + 1.0);
                        e[j] = s * e[i] * r;
                        s = 1.0 / r;
                        c /= r;
                    }
                    else {
                        c = e[i] / p;
                        r = sqrt(c * c + 1.0);
                        e[j] = s * p * r;
                        s = c / r;
                        c = 1.0 / r;
                    }
                    p = c * d[i] - s * g;
                    d[j] = h + s * (c * g + s * d[i]);
                    for (k = 1; k <= n; ++k) {
                        kn = (k - 1) * n;
                        h = z[kn + j];
                        z[kn + j] = s * z[kn + i] + c * h;
                        z[kn + i] = c * z[kn + i] - s * h;
                    }
                }
                e[l] = s * p;
                d[l] = c * p;
                if (fabs(e[l]) <= b)
                    break;
            }
            if (jj > maxits)
                return(-1);
        }   
        d[l] += f;
    }
    for (i = 1; i < n; ++i) {
        k = i;
        p = d[i];

        for (j = i + 1; j <= n; ++j) {
            if (d[j] > p) {                /* use descending order */
                k = j;
                p = d[j];
            }
        }
        if (k != i) {
            d[k] = d[i];
            d[i] = p;
            for (j = 1; j <= n; ++j) {
                j2 = (j - 1) * n;
                p = z[j2 + i];
                z[j2 + i] = z[j2 + k];
                z[j2 + k] = p;
            }
        }
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  simitz                                                                  */
/*  ##                                                                      */
/*  Iterative computation of eigenvalues largest in magnitude and           */
/*  corresponding eigenvectors of a real symmetric matrix.                  */
/*                                                                          */
/*  Adapted from: Paul J. Nikolai, CACM 538.                                */ 
/*  Transactions in math. software 5 (1979), pp. 118 - 125.                 */
/*                                                                          */
/*  n = order of C matrix.                                                  */
/*  p = 1 plus number of required eigenvalues/vectors                       */
/*  km = max number of iterations                                           */
/*  eps = desired accuracy                                                  */
/*  x = (p,n) matrix            at return contains eigenvectors             */
/*  d = (p)-vector              at return contains eigenvalues              */
/*  em = # of requested eigenvalues, must be 1 <= em < p <= n.              */ 
/*  at return: *nit contains number of iterations performed.                */
/*                                                                          */
/*  Return: # of eigenvalues calculated, or -1 if error (insuff. memory)    */
             
int simitz(int n,int p,int km,double eps,int em,double *x,double *d,int *nit)
{
    register int i,j,l,k;
    int ik,jk,g,ig,h,z1,z2,ks,m,m1,lf,l1,jp,orig,gflag;
    int err,cxa,fa,rqa,ua,ba,mxit;
    double e,e1,e2,s,t,tol;
    double *xp,*xp1,*cx,*f,*rq,*u,*b;

    mxit = km;
    *nit = 0;
    tol = DBLMIN / EPSI;  /* DBLMIN smallest number representable   */
                          /* EPSI is machine's epsilon            */

    /*  need working arrays   

        cx[i], i = 1,...,p
        f[i],  i = 1,...,p
        rq[i], i = 1,...,p
        u[i],  i = 1,...,2 * n
        b(i,j) i,j = 1,...,p
    */

    err = -1;
    cxa = fa = rqa = ua = ba = 0;

    if (!(cx = (double *) calloc (p + 1,sizeof(double))))
        goto SIMITZFin;
    cxa = p + 1;
    memrq(cxa,sizeof(double));

    if (!(f = (double *) calloc (p + 1,sizeof(double))))
        goto SIMITZFin;
    fa = p + 1;
    memrq(fa,sizeof(double));

    if (!(rq = (double *) calloc (p + 1,sizeof(double))))
        goto SIMITZFin;
    rqa = p + 1;
    memrq(rqa,sizeof(double));

    if (!(u = (double *) calloc (2 * n + 1,sizeof(double))))
        goto SIMITZFin;
    ua = 2 * n + 1;
    memrq(ua,sizeof(double));

    if (!(b = (double *) calloc (p * p + 1,sizeof(double))))
        goto SIMITZFin;
    ba = p * p + 1;
    memrq(ba,sizeof(double));

    if (SILENTFlg < 2)
        printfe("\n  Iter  Eigenvalues  Eigenvectors  Criterion\n");

    e = 0.0;  
    jp = g =  h = z1 = z2 = ks = 0;
    m = ig = 1;

    cx[p] = 0.0;
    for (l = 1; l <= p; ++l) {
        f[l] = 4.0;
        rq[l] = 0.0;
    }
    if (km <= 0)
        goto L50;

    for (l = 1; l <= p; ++l) {
        xp = x + (l - 1) * n;
        for (j = 1; j <= n; ++j) {
            xp[j] = 2.0 * g_ranf() - 1.0;    
        }
    }

L50:
    km = iabs(km);

    gflag = 60;
    lf = ig;
    l1 = p;
    goto L990;

    /* RAYLEIGH-RITZ STEP */

L60:
    for (k = ig; k <= p; ++k) {
        g_op(n,x + (k - 1) * n,u);
        xp = x + (k - 1) * n;
        for (j = 1; j <= n; ++j)  
            xp[j] = u[j];
    }
    gflag = 90;
    lf = ig;
    l1 = p;
    goto L990;

L90:
    if (ks != 0)
        goto L150;

    /* MEASURES AGAINST UNHAPPY CHOICE OF INITIAL VECTORS */

    for (k = 1; k <= p; ++k) {
        if (b[(k-1)*p + k] != 0.0) 
            continue;
        xp = x + (k - 1) * n;
        for (i = 1; i <= n; ++i)  
            xp[i] = 2.0 * g_ranf() - 1.0;   
        ks = 1;
    }
    if (ks != 1)
        goto L150;  

L140:
    gflag = 60;
    lf = 1;
    l1 = p;
    goto L990;

L150:
    for (k = ig; k <= p; ++k) {

        for (l = k; l <= p; ++l) {
            s = 0.0;    
            for (i = l; i <= p; ++i)  
                s += b[(k - 1) * p + i] * b[(l - 1) * p + i];
            b[(k - 1) * p + l] = -s;

        }
    }
    g_tred2(p,p - g,b + (ig-1)*p + ig-1,d + ig-1,u,b + (ig-1)*p + ig-1);
    l = g_imtql2(p,p - g,d + ig-1,u,b + (ig-1)*p + ig-1,mxit);
   
    f[ig] = dmax(f[ig],1000.0 * (double)l);

    for (k = ig; k <= p; ++k) 
        d[k] = sqrt(dmax(-d[k],0.0));

    /* REORDERING EIGENVALUES AND EIGENVECTORS ACCORDING TO SIZE OF
       THE FORMER IS ACCOMPLISHED IN SUBROUTINE IMTQL2. */

    for (j = 1; j <= n; ++j) {
        for (k = ig; k <= p; ++k) {
            s = 0.0;    
                 
            for (l = ig; l <= p; ++l)
                s += x[(l - 1) * n + j] * b[(k - 1) * p + l];/* WK(L,K,1);  */
            u[k] = s;
        }
        for (k = ig; k <= p; ++k)   
            x[(k - 1) * n + j] = u[k];       
    }
    ks++;          
    e = dmax(d[p],e);

    /* RANDOMIZATION */

    if (3 < z1)
        goto L260;

    xp = x + (p - 1) * n;
    for (j = 1; j <= n; ++j)  
        xp[j] = 2.0 * g_ranf() - 1.0;    
        
    jp = p - 1;
    gflag = 260;
    lf = l1 = p;
    goto L990;

    /* COMPUTE CONTROL QUANTITIES CX(I). */

L260:
    for (k = ig; k <= jp; ++k) {
        s = (d[k] - e) * (d[k] + e);
        if (s > 0.0)
            goto L280;

        cx[k] = 0.0;
        continue; 
L280:
        if (e != 0.0)
            goto L300;

        cx[k] = 1000.0 + log(d[k]);
        continue; 

L300:
        cx[k] = log((d[k] + sqrt(s)) / e);
    }

    /* ACCEPTANCE TEST FOR EIGENVALUES INCLUDING ADJUSTMENT OF EM AND
       H SUCH THAT D(EM) .GT. E, D(H) .GT. E AND D(EM) DOES NOT
       OSCILLATE STRONGLY */

    i = z1 - 1;
    k = g;
L320:
    k++;      
    if (em < k) 
        goto L370;

    if (d[k] <= e) 
        goto L360;

    if (i <= 0) 
        goto L320;

    if (d[k] > 0.999 * rq[k])                     
        goto L320;

L360:
    em = k - 1;

L370:
    if (em == 0) 
        goto L1130;

L380:
    k = h;
    s = 1.0 + 0.1 * eps;

L390:
    k++;       

    if (d[k] == 0.0)                    
        goto L410;

    if (d[k] <= s * rq[k]) 
        goto L390;

L410:   
    h = k - 1;
    k = em;

L420:
    if (++k > h)
        goto L450;

    if (d[k] > e) 
        goto L420;
                
    h = k - 1;

    /* ACCEPTANCE TEST FOR EIGENVECTORS */

L450:
    l = g;
    e2 = 0.0;     
    for (k = ig; k <= jp; ++k) {
        if (k != l + 1)                     
            goto L510;

        /* CHECK FOR NESTED EIGENVALUES */

        l = l1 = k;
        s = 0.5 / (double)ks;  
        t = 1.0 / (double)(ks * m);

L470:
        l++;       

        if (l > jp) 
            goto L490;

        if (cx[l] * (cx[l] + s) + t - cx[l-1] * cx[l-1] > 0.0)
            goto L470;
L490:
        l--;      
        if (l <= h) 
            goto L510;

        l = l1 - 1;
        goto L600;

L510:
        g_op(n,x + (k - 1) * n,u);
        s = 0.0;     
        for (j = 1; j <= l; ++j) {
            if (fabs(d[j] - d[k]) >= 0.01 * d[k]) 
                goto L540;
   
            t = g_ip(n,u,x + (j - 1) * n);
            xp = x + (j - 1) * n;
            for (i = 1; i <= n; ++i)  
                u[i] -= t * xp[i];

            s += t * t;
L540:       ;
        }
        t = 1.0;    
        if (s != 0.0)  
            t = g_ip(n,u,u);
   
        e2 = dmax(e2,sqrt(t / (s + t)));

        if (k != l) 
            goto L590;

        /* TEST FOR ACCEPTANCE OF GROUP OF EIGENVECTORS */
   
        if (l >= em && d[em] * f[em] < eps * (d[em] - e))
            g = em;

        if (e2 - f[l] >= 0.0) 
            goto L580;

        for (j = l1; j <= l; ++j)
            f[j] = e2;

L580:
        if (l <= em && d[l] * f[l] < eps * (d[l] - e))
            g = l;

L590:   ;

    }

L600:
    ig = g + 1;
    if (e > 0.04 * d[1]) 
        goto L620;

    m = k = 1;
    goto L630;

L620:
    e2 = 2.0 / e;
    e1 = 0.51 * e2;
    k = 2 * (int)(4.0 / dmin(cx[1],4.0));
    m = imin(m,k); 

    /* REDUCE EM IF CONVERGENCE WOULD BE TOO SLOW. */

L630:
    if (f[em] == 0.0)                    
        goto L690;

    if ((double)ks >= 0.9 * (double)km)                  
        goto L690;

    s = (double)k * cx[em];    

    if (s >= 0.05)                  
        goto L670;

    t = 0.5 * s * cx[em];    
    goto L680;

L670:
    t = cx[em] + log(0.5 + 0.5 * exp(-2.0 * s)) / (double)k;    

L680:
    s = log(d[em] * f[em] / (eps * (d[em] - e)));

    if (s * (double)ks > t * (double)((km - ks) * km) )
        em--;       

L690:
    for (k = ig; k <= jp; ++k)
        rq[k] = d[k];

    if (SILENTFlg < 2)
        printfe("%6d  %8d     %8d      %17.10e\n",ks,h,g,f[1]);

    if (g >= em || ks >= km)               
        goto L1130;

L710:
    if (ks + m <= km)                   
        goto L730;

    z2 = -1;

    if (m > 1)
        m = 2 * ((km - ks + 1) / 2);

L730:
    m1 = m;

    /* SHORTCUT LAST INTERMEDIATE BLOCK IF ALL F(I) ARE SUFFICIENTLY SMALL */

    if (l < em)                  
        goto L780;

    s = d[em] * f[em] / (eps * (d[em] - e));
    t = s * s - 1.0;  

    if (t <= 0.0)             
        goto L60;

    s = log(s + sqrt(t)) / (cx[em] - cx[h + 1]);
    m1 = 2 * (int)(0.5 * s + 1.01);

    if (m1 <= m) 
        goto L770;

    m1 = m;
    goto L780;

L770:
    z2 = -1;

    /* CHEBYSHEV ITERATION */

L780:
    if (m < 1)
        goto L900;
    else if (m > 1)
        goto L820;

L790:
    for (k = ig; k <= p; ++k) {
        g_op(n,x + (k - 1) * n,u);
        xp = x + (k - 1) * n;
        for (i = 1; i <= n; ++i)
            xp[i] = u[i];     
    }
    goto L900;

L820:
    l1 = m1 - 4;
    for (k = ig; k <= p; ++k) {
        g_op(n,x + (k - 1) * n,u);
        for (i = 1; i <= n; ++i) {
            ik = i + n;
            u[ik] = e1 * u[i];        
        }
        g_op(n,u + n,u);
        xp = x + (k - 1) * n;
        for (i = 1; i <= n; ++i) 
            xp[i] = e2 * u[i] - xp[i];     

        if (l1 < 0)                
            goto L890; 

        for (j = 4; j <= m1; j += 2) {
            g_op(n,x + (k - 1) * n,u);
            for (i = 1; i <= n; ++i) {
                ik = i + n;
                u[ik] = e2 * u[i] - u[ik];     
            }
            g_op(n,u + n,u);

            xp = x + (k - 1) * n;
            for (i = 1; i <= n; ++i) 
                xp[i] = e2 * u[i] - xp[i]; 

        }
L890: ;

    }

L900:
    gflag = 910;
    lf = ig;
    l1 = p;
    goto L990;

    /* DISCOUNTING THE ERROR QUANTITIES F */

L910:
    if (g >= h)                   
        goto L970;

    if (m != 1)                 
        goto L950;

    for (k = ig; k <= h; ++k)
        f[k] *= (d[h + 1] / d[k]);

    goto L970;

L950:
    t = exp(-(double)m1 * cx[h + 1]);      

    for (k = ig; k <= h; ++k) {
        s = exp(-(double)m1 * (cx[k] - cx[h + 1]));
        f[k] *= s * (1.0 + t * t) / (1.0 + (s * t) * (s * t));
    }

L970:         
    ks += m1;   
    z2 -= m1;

    /* POSSIBLE REPETITION OF INTERMEDIATE STEPS */

    if (z2 >= 0)                     
        goto L710;

    z1++;        
    z2 = 2 * z1;
    m *= 2;    
    goto L60;

    /* PERFORMS ORTHONORMALIZATION OF COLUMNS 1 THROUGH L1 OF ARRAY
        X ASSUMING THAT COLUMNS 1 THROUGH LF - 1 ARE ALREADY ORTHONORMAL */

L990:
    for (k = lf; k <= l1; ++k) {
        orig = 1;     
L1000:
        t = 0.0;  
        jk = k - 1;

        if (jk <= 0)
            goto L1040;

        for (i = 1; i <= jk; ++i) {
            s = g_ip(n,x + (i - 1) * n,x + (k - 1) * n);
            if (orig)
                b[(i - 1) * p + k] = s;

            t += s * s;
            xp = x + (k - 1) * n;
            xp1 = x + (i - 1) * n;
            for (j = 1; j <= n; ++j) 
                xp[j] -= s * xp1[j];             
        }

L1040:
        s = g_ip(n,x + (k - 1) * n, x + (k - 1) * n);
        t += s;   
        if (s <= 0.01 * t)
            goto L1060;

        if (t > tol)
            goto L1080;

L1060:
        orig = 0;       
        if (s > tol) 
            goto L1000;

        s = 0.0;   

L1080:
        s = sqrt(s);
        b[(k - 1) * p + k] = s;
        if (s == 0.0) 
            goto L1100;

        s = 1.0 / s;

L1100:
        xp = x + (k - 1) * n;
        for (j = 1; j <= n; ++j) 
            xp[j] *= s;                   
    } 
    if (gflag == 60)
        goto L60;
    else if (gflag == 90)
        goto L90;
    else if (gflag == 260)
        goto L260;
    else if (gflag == 910)
        goto L910;
    else
        goto L1140;

L1130:
    em = g;

    /* SOLVE EIGENVALUE PROBLEM OF PROJECTION OF MATRIX C. */

    gflag = 1140;
    lf = 1;
    l1 = jp;
    goto L990;

L1140:
    for (k = 1; k <= jp; ++k) {
        g_op(n,x + (k - 1) * n,x + (p - 1) * n);
        for (i = 1; i <= k; ++i)  
            b[(i - 1) * p + k] = -g_ip(n,x + (i - 1) * n,x + (p - 1) * n);
    }
    g_tred2(p,jp,b,d,u,b);
    l = g_imtql2(p,jp,d,u,b,mxit);
    
    f[ig] = dmax(f[ig],1000.0 * (double)l);

    /* ARRANGE EIGENVALUES IN ORDER OF DECREASING ABSOLUTE VALUE. */

    for (j = 1; j <= jp; ++j) {
        k = j;
        for (i = j; i <= jp; ++i) {
            if (fabs(d[i]) > fabs(d[k]))
                k = i;
        }
        if (k <= j)
            goto L1200;

        t = d[k];
        d[k] = d[j];
        d[j] = t;

        for (i = 1; i <= jp; ++i) {
            t = b[(k - 1) * p + i];
            b[(k - 1) * p + i] = b[(j - 1) * p + i];
            b[(j - 1) * p + i] = t;
        }
L1200:
        d[j] = -d[j];
    }
    for (j = 1; j <= n; ++j) {
        for (i = 1; i <= jp; ++i) {
            s = 0.0;    
            for (k = 1; k <= jp; ++k)
                s += x[(k - 1) * n + j] * b[(i - 1) * p + k];
            u[i] = s;
        }
        for (i = 1; i <= jp; ++i)
            x[(i - 1) * n + j] = u[i];
    }
    km = ks;
    *nit = km;
    d[p] = e;
    err = em;

SIMITZFin:
    if (cxa > 0) {
        free((char *)cx);
        memrq(-cxa,sizeof(double));
    }
    if (fa > 0) {
        free((char *)f);
        memrq(-fa,sizeof(double));
    }
    if (rqa > 0) {
        free((char *)rq);
        memrq(-rqa,sizeof(double));
    }
    if (ua > 0) {
        free((char *)u);
        memrq(-ua,sizeof(double));
    }
    if (ba > 0) {
        free((char *)b);
        memrq(-ba,sizeof(double));
    }
    return(err);
}            

/* ------------------------------------------------------------------------ */
/*  g_ranf()    returns a random number in (0,1).                           */
/*                                                                          */

double g_ranf(void)
{
    static int k = 0; 
    static int init = 0;
    static double seed = 0.1;

    if (init == 0) {
        if (seed > 1.0)
            seed = 1.0 / seed;
        k = 2 * (int)(16384.0 * seed) - 1;        
        init = 1;
    }
    k = (3125 * k) % 65536;
    return(fabs((double)(k - 32768) / 32768.0));  
}

/* ------------------------------------------------------------------------ */
/*  g_ip(n,z,w) returns inner product of z and w.                           */
/*                                                                          */

double g_ip(int n,double *z,double *w)
{
    register int i;
    double tmp = 0.0;

    for (i = 1; i <= n; ++i)
        tmp += z[i] * w[i];
    return(tmp);
}

/* ------------------------------------------------------------------------ */
/*  g_op(n,z,w) returns w = Cz where C is the n x n matrix to be used       */
/*              for the calculation of eigenvalues/vectors.                 */

void g_op(int n,double *z,double *w)
{
    g_op1(n,z,w);
}

/* ------------------------------------------------------------------------ */
/*  g_imtql2(nm,n,d,e,z,mxit)                                               */
/*                                                                          */
/*  This function corresponds to subroutine IMTQL2 in the SIMITZ            */
/*  package. Originally a translation of the ALGOL procedure IMTQL2,        */
/*  NUM. MATH. 12, 377-383(1968) BY MARTIN AND WILKINSON,                   */
/*  AS MODIFIED IN NUM. MATH. 15, 450(1970) BY DUBRULLE.                    */
/*  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).         */
/*                                                                          */
/*  nm must be column dimension of z.                                       */  
/*                                                                          */
/*  Return:  0 if OK, j if the jth eigenvalue has not been determined       */
/*           after 30 iterations.                                           */

int g_imtql2(int nm,int n,double *d,double *e,double *z,int mxit)
{
    register int i,l,j,m,ii,k;
    int err,mml;
    double p,g,r,s,c,f,b,tmp;
              
    err = 0;
    if (n == 1) 
        return(0);   
 
    for (i = 2; i <= n; ++i)  
        e[i-1] = e[i];
 
    e[n] = 0.0;
 
    for (l = 1; l <= n; ++l) {
        j = 0;  

        /* LOOK FOR SMALL SUB-DIAGONAL ELEMENT */            
L105:
        for (m = l; m <= n; ++m) {
            if (m == n)              
                break;
            if (fabs(e[m]) <= EPSI * (fabs(d[m]) + fabs(d[m + 1])))
                break;    
        }
        p = d[l];
        if (m == l)
            goto L240;

        if (j >= mxit)
            return(l);

        j++;       

        /* FORM SHIFT */          

        g = (d[l + 1] - p) / (2.0 * e[l]);
        r = sqrt(g * g + 1.0);

        if (g >= 0.0)
            tmp = r;
        else
            tmp = -r;

        g = d[m] - p + e[l] / (g + tmp);
        c = s = 1.0;
        p = 0.0;
        mml = m - l;

        /* FOR I=M-1 STEP -1 UNTIL L DO -- */          

        for (ii = 1; ii <= mml; ++ii) {
            i = m - ii;
            f = s * e[i];
            b = c * e[i];
            if (fabs(f) < fabs(g))
                goto L150;

            c = g / f;
            r = sqrt(c * c + 1.0);
            e[i + 1] = f * r;
            s = 1.0 / r;
            c *= s;    
            goto L160;
L150:
            s = f / g;
            r = sqrt(s * s + 1.0);
            e[i + 1] = g * r;
            c = 1.0 / r;
            s *= c;    
L160:
            g = d[i + 1] - p;
            r = (d[i] - g) * s + 2.0 * c * b;
            p = s * r;
            d[i + 1] = g + p;
            g = c * r - b;
 
            /* FORM VECTOR */            

            for (k = 1; k <= n; ++k) {
                f = z[i * nm + k];
                z[i * nm + k] = s * z[(i-1) * nm + k] + c * f;
                z[(i - 1) * nm + k] = c * z[(i - 1) * nm + k] - s * f;
            }
        }
        d[l] -= p;
        e[l] = g;
        e[m] = 0.0;
        goto L105;
L240: ;
    }

    /* ORDER EIGENVALUES AND EIGENVECTORS */             
       
    for (ii = 2; ii <= n; ++ii) {
        i = ii - 1; 
        k = i;
        p = d[i];
 
        for (j = ii; j <= n; ++j) {
            if (d[j] >= p)
                continue;  
            k = j;
            p = d[j];
        }
        if (k == i)
            goto L300;
        d[k] = d[i];
        d[i] = p;
 
        for (j = 1; j <= n; ++j) {
            p = z[(i - 1) * nm + j];
            z[(i - 1) * nm + j] = z[(k - 1) * nm + j];
            z[(k - 1) * nm + j] = p;
        }
  
L300:    ;     
    }
    return(0);    
}

/* ------------------------------------------------------------------------ */
/*  g_tred2(nm,n,a,d,e,z)                                                   */
/*                                                                          */
/*  This function corresponds to subroutine TRED2 in the SIMITZ             */
/*  package. Originally a translation of the ALGOL procedure TRED2,         */
/*  NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.         */
/*  HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).         */
/*                                                                          */
  
void g_tred2(int nm,int n,double *a,double *d,double *e,double *z) 
{
    register int i,j,ii,k,l;
    int jp1;
    double scale,h,f,g,hh;

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= i; ++j)  
            z[(j - 1) * nm + i] = a[(j - 1) * nm + i];
    }
    if (n == 1)
        goto L320;

    /* FOR I=N STEP -1 UNTIL 2 DO -- */

    for (ii = 2; ii <= n; ++ii) {
        i = n + 2 - ii;
        l = i - 1;
        scale = h = 0.0;
        if (l < 2)
            goto L130;

        /* SCALE ROW (ALGOL TOL THEN NOT NEEDED) */

        for (k = 1; k <= l; ++k) 
            scale += fabs(z[(k - 1) * nm + i]);
  
        if (scale != 0.0)
            goto L140;

L130:
        e[i] = z[(l - 1) * nm + i];
        goto L290;
  
L140:
        for (k = 1; k <= l; ++k) {
            z[(k - 1) * nm + i] /= scale;
            h += z[(k - 1) * nm + i] * z[(k - 1) * nm + i];
        }
        f = z[(l - 1) * nm + i];   
        if (f >= 0.0)
            g = -sqrt(h); 
        else
            g = sqrt(h); 

        e[i] = scale * g;
        h -= f * g;
        z[(l - 1) * nm + i] = f - g;
        f = 0.0;
  
        for (j = 1; j <= l; ++j) {
            z[(i - 1) * nm + j] = z[(j - 1) * nm + i] / h;
            g = 0.0;

            /* FORM ELEMENT OF A*U */           

            for (k = 1; k <= j; ++k)  
                g += z[(k - 1) * nm + j] * z[(k - 1) * nm + i];
   
            jp1 = j + 1;
            if (l < jp1)
                goto L220;
   
            for (k = jp1; k <= l; ++k)
                g += z[(j - 1) * nm + k] * z[(k - 1) * nm + i];

            /* FORM ELEMENT OF P */
L220:
            e[j] = g / h;
            f += e[j] * z[(j - 1) * nm + i];      
        }
        hh = f / (h + h);

        /* FORM REDUCED A */            

        for (j = 1; j <= l; ++j) {
            f = z[(j - 1) * nm + i];    
            g = e[j] - hh * f;
            e[j] = g;
            for (k = 1; k <= j; ++k)
                z[(k - 1) * nm + j] = z[(k - 1) * nm + j] - f * e[k] -
                                        g * z[(k - 1) * nm + i];
        }
L290:
        d[i] = h;
    }
   
L320:
    d[1] = e[1] = 0.0;

    /* ACCUMULATION OF TRANSFORMATION MATRICES */

    for (i = 1; i <= n; ++i) {
        l = i - 1;
        if (d[i] == 0.0)
            goto L380;
  
        for (j = 1; j <= l; ++j) {
            g = 0.0;
            for (k = 1; k <= l; ++k)  
                g += z[(k - 1) * nm + i] * z[(j - 1) * nm + k];
               
            for (k = 1; k <= l; ++k)
                z[(j - 1) * nm + k] -= g * z[(i - 1) * nm + k];
        }
L380:
        d[i] = z[(i - 1) * nm + i];
        z[(i - 1) * nm + i] = 1.0;

        if (l < 1)
            continue;               
  
        for (j = 1; j <= l; ++j) {
            z[(j - 1) * nm + i] = 0.0;
            z[(i - 1) * nm + j] = 0.0;
        }
    }
}

/* ------------------------------------------------------------------------ */
/*                                                                          */
/*  Function: eigen2()                                                      */
/*                                                                          */
/*      Berechnung der Eigenwerte und Eigenvektoren einer allgemeinen       */
/*      reellen Matrix.                                                     */
/*                                                                          */
/*      Ref.: P.J. Eberlein, J. Boothroyd, Solution to the Eigenproblem     */
/*      by a Norm Reducing Jacobi Type Method, in: Handbook for Automatic   */
/*      Computation, Vol. II, = J.H. Wilkinson, C. Reinsch (eds.), Linear   */
/*      Algebra, Springer-Verlag 1971, S. 327 - 338                         */
/*                                                                          */
/*      Beachte: Die Eigenvektoren werden normiert. Reelle Eigenwerte       */
/*      werden so normiert, dass ihre euklidische Norm = 1 ist. Komplexe    */
/*      Eigenwerte werden so normiert, dass die Komponente mit dem          */
/*      groessten Absolutbetrag den Realteil 1 und den Imaginaerteil 0      */
/*      erhaelt. Dies entspricht der Normierungsmethode des Algorithmus     */
/*      343 in Comm. ACM 11 (1963), 820ff.                                  */
/*                                                                          */
/*                                                                          */
/*      Aufruf der Funktion                                                 */
/*                                                                          */
/*      int eiegen2(n,a,evr,evi,vec,mode,itmax)                             */
/*                                                                          */
/*      int n       input:  Ordnung der Inputmatrix a                       */
/*      double *a   input:  die (n,n)-Inputmatrix a als ein eindimens.      */
/*                          Feld a(i,j) = a[(i - 1) * n + j]                */
/*                  output: zerstoert                                       */
/*      double *evr input:  eindimensionales Feld evr(i) (i = 1,n)          */
/*                  output: die Realteile der Eigenwerte                    */
/*      double *evi input:  eindimensionales Feld evr(i) (i = 1,n)          */
/*                  output: die Imaginaerteile der Eigenwerte               */
/*      double *vec input:  ein eindimensionales Feld fuer eine (n,n)-      */
/*                          Matrix                                          */
/*                  output: die Eigenvektoren (in komprimierter Form)       */
/*      int itmax   input:  maximale Anzahl der Iterationen (etwa 50)       */
/*      int mode    input:   0 dann werden nur Eigenwerte berechnet         */
/*                          +1 dann werden Rechtseigenvektoren berechnet    */
/*                          -1 dann werden Linkseigenvektoren berechnet     */
/*                                                                          */
/*      Return:     > 0  Anzahl der Iterationen (erfolgreiche Rechnung)     */
/*                   -2  keine Berechnung moeglich                          */
/*                   -3  keine Konvergenz erreicht                          */
/*                                                                          */

int eigen2(int n,double *a,double *evr,double *evi,double *vec,int mode,int itmax)
{
    register int i,j,k,m;
    int jj,jn,it,mark;
    double aii,aij,aji,h,g,hj,yh,aik,aim,te,tee,aki,ami,tep,tem;
    double d,akm,amk,c,e,cx,sx,cot2x,sig,cotx,cos2x,sin2x,den,tanhy,chy;
    double shy,c1,c2,s1,s2,tik,tim,tki,tmi,ep,eps,vr,vi,vri,tmp,tmp1;

    ep = 0.5 * EPSI;
    eps = 1.1 * sqrt(ep);

    mark = 0;
    if (mode) {
        k = 1;
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (i == j)
                    vec[k++] = 1.0;
                else
                    vec[k++] = 0.0;
            }
        }
    }

    for (it = 1; it <= itmax; ++it) {

        if (mark)
            return(-2);
        for (i = 1; i < n; ++i) {
            aii = a[(i - 1) * n + i];
            for (j = i + 1; j <= n; ++j) {
                aij = a[(i - 1) * n + j];
                aji = a[(j - 1) * n + i];
                if (fabs(aij + aji) > eps ||
                   (fabs(aij - aji) > eps &&
                    fabs(aii - a[(j - 1) * n + j]) > eps))
                    goto Cont;
            }
        }
        goto Fin;
Cont:
        mark = 1;
        for (k = 1; k < n; ++k) {
            for (m = k + 1; m <= n; ++m) {
                h = g = hj = yh = 0.0;
                for (i = 1; i <= n; ++i) {
                    aik = a[(i - 1) * n + k];
                    aim = a[(i - 1) * n + m];
                    te  = aik * aik;
                    tee = aim * aim;
                    yh += te - tee;
                    if (i != k && i != m) {
                        aki = a[(k - 1) * n + i];
                        ami = a[(m - 1) * n + i];
                        h   = h + aki * ami - aik * aim;
                        tep = te + ami * ami;
                        tem = tee + aki * aki;
                        g  += tep + tem;
                        hj += tem - tep;
                    }
                }
                h  += h;
                d   = a[(k - 1) * n + k] - a[(m - 1) * n + m];
                akm = a[(k - 1) * n + m];
                amk = a[(m - 1) * n + k];
                c   = akm + amk;
                e   = akm - amk;
                if (fabs(c) <= ep) {
                    cx = 1.0;
                    sx = 0.0;
                }
                else {
                    if ((cot2x = d / c) < 0.0)
                        sig = -1.0;
                    else
                        sig = 1.0;
                    cotx = cot2x + (sig * sqrt(1.0 + cot2x * cot2x));
                    sx = sig / sqrt(1.0 + cotx * cotx);
                    cx = sx * cotx;
                }
                if (yh < 0.0) {
                    tem = cx;
                    cx  = sx; 
                    sx  = -tem;
                }
                cos2x = cx * cx - sx * sx;
                sin2x = 2.0 * sx * cx;
                d     = d * cos2x + c * sin2x;
                h     = h * cos2x - hj * sin2x;
                den   = g + 2.0 * (e * e + d * d);
                tanhy = (e * d - h / 2.0) / den;
                if (fabs(tanhy) <= ep) {
                    chy = 1.0;
                    shy = 0.0;
                }
                else {
                    chy = 1.0 / sqrt(1.0 - tanhy * tanhy);
                    shy = chy * tanhy;
                }
                c1 =  chy * cx - shy * sx;
                c2 =  chy * cx + shy * sx;
                s1 =  chy * sx + shy * cx;
                s2 = -chy * sx + shy * cx;
                if (fabs(s1) > ep || fabs(s2) > ep) {
                    mark = 0;
                    for (i = 1; i <= n; ++i) {
                        aki = a[(k - 1) * n + i];
                        ami = a[(m - 1) * n + i];
                        a[(k - 1) * n + i] = c1 * aki + s1 * ami;
                        a[(m - 1) * n + i] = s2 * aki + c2 * ami;
                        if (mode < 0) {
                            tki = vec[(k - 1) * n + i];
                            tmi = vec[(m - 1) * n + i];
                            vec[(k - 1) * n + i] = c1 * tki + s1 * tmi;
                            vec[(m - 1) * n + i] = s2 * tki + c2 * tmi;
                        }
                    }
                    for (i = 1; i <= n; ++i) {
                        aik = a[(i - 1) * n + k];
                        aim = a[(i - 1) * n + m];
                        a[(i - 1) * n + k] = c2 * aik - s2 * aim;
                        a[(i - 1) * n + m] = -s1 * aik + c1 * aim;
                        if (mode > 0) {
                            tik = vec[(i - 1) * n + k];
                            tim = vec[(i - 1) * n + m];
                            vec[(i - 1) * n + k] = c2 * tik - s2 * tim;
                            vec[(i - 1) * n + m] = -s1 * tik + c1 * tim;
                        }
                    }
                }
            }
        }
    }
    return(-3);
 
Fin:                            /* Uebertragen der Eigenwerte aus a in  */
                                /* die Felder evr und evi               */
    i = 0;
    while (++i <= n) {
        tmp = a[(i - 1) * n + i];
        tmp1 = a[(i - 1) * n + i + 1];
        if (i < n && fabs(tmp1) > eps) {
            evr[i] = tmp;
            evi[i] = tmp1;
            evr[i + 1] = tmp;
            evi[i + 1] = -tmp1;
            i++;
        }
        else {
            evr[i] = tmp;
            evi[i] = 0.0;
        }
    }
    if (!mode)
        return(it);
                                    /* Normalisierung der Eigenvektoren */
    i = 1;
    while (i <= n) {
        if (fabs(evi[i]) < ep) {    /* wenn Eigenwert reell */
            tmp = 0.0;
            for (j = 1; j <= n; ++j) {
                tmp1 = vec[(j - 1) * n + i];
                tmp += tmp1 * tmp1;
            }
            if (tmp > 0.0) {
                tmp = sqrt(tmp);
                for (j = 1; j <= n; ++j)  
                    vec[(j - 1) * n + i] /= tmp;
            }
            i++;
        }
        else {                      /* wenn konj. komplexe Eigenwerte */
            jj = 0;
            vri = 0.0;
            for (j = 1; j <= n; ++j) {
                jn = (j - 1) * n + i;
                vr = vec[jn];
                vi = vec[jn + 1];
                tmp  = vr * vr + vi * vi;
                if (vri < tmp) {
                    vri = tmp;
                    jj = j;
                }
            }
            if (vri > 0.0) {
                jn = (jj - 1) * n + i;
                vr = vec[jn] / vri;
                vi = vec[jn + 1] / vri;
                for (j = 1; j <= n; ++j) {
                    jn = (j - 1) * n + i;
                    tmp = vec[jn];
                    tmp1 = vec[jn + 1];
                    vec[jn] = tmp * vr + tmp1 * vi;
                    vec[jn + 1] = tmp1 * vr - tmp * vi;
                }
            }
            i += 2;
        }
    }
    return(it);
}

/* -###-------------------------------------------------------------------- */
/*                                                                          */
/*  Function: eigen1()                                                      */
/*                                                                          */
/*      Berechnung der Eigenwerte und Eigenvektoren einer allgemeinen       */
/*      reellen quadratischen Matrix                                        */
/*                                                                          */
/*      Ref. J. Grad, M.A. Brebner, Eigenvalues and Eigenvectors of a       */
/*      Real General Matrix, Algorithm 343, Communications of the ACM       */
/*      Vol. 11 (1963), 820 - 826. Comment in: Comm. ACM 13 (1970), 122     */
/*                                                                          */
/*      eigen1(n,a,evr,evi,vecr,veci,indic)                                 */
/*                                                                          */
/*      int n         input:  Ordnung der (n,n)-Matrix a                    */
/*      double *a     input:  Die Matrix a als ein eindimensionales Feld    */
/*      double *evr   input:  Eindimensionales Feld evr[i] (i = 1,n)        */
/*                    output: Realteile der Eigenwerte (i = 1,n)            */
/*      double *evi   input:  Eindimensionales Feld evi[i] (i = 1,n)        */
/*                    output: Imaginaerteile der Eigenwerte (i = 1,n)       */
/*      double *vecr  input:  Eine (n,n)-Matrix als eindimensionales Feld   */
/*                    output: Matrix der Eigenvektoren (real)               */
/*      double *veci  input:  Eine (n,n)-Matrix als eindimensionales Feld   */
/*                    output: Matrix der Eigenvektoren (imaginaer)          */
/*      int *indic    input:  Ein eindimensionales Feld indic[i] (i = 1,n)  */
/*                    output: Indikatoren mit folgender Bedeutung:          */
/*                            2 = Eigenwert und Eigenvektor gefunden        */
/*                            1 = Nur der Eigenwert wurde gefunden          */
/*                            0 = Weder Eigenwert noch Eigenvektor gefunden */
/*                                                                          */
/*      Return:  0 wenn indic[i] == 2 fuer alle i = 1,n                     */
/*               1 andernfalls                                              */
/*             < 0 wenn ein Fehlercode (AllocErr)                           */
/*                                                                          */
/*      Alle Vektoren (a,vecr,veci) werden als eindimensionale Felder       */
/*      behandelt, zum Beispiel:                                            */
/*                                                                          */
/*          a(i,j) = a[(i - 1) * n + j]                                     */
/*                                                                          */
/*      Die Eigenvektoren werden als Spaltenvektoren in vecr und veci       */
/*      uebergeben.                                                         */
/*                                                                          */
/*      Die Normierung der Eigenvektoren erfolgt folgendermassen: Reelle    */
/*      Eigenvektoren werden so normiert, dass ihre euklische Norm (Summe   */
/*      der Quadrate ihrer Komponenten) = 1 ist. Komplexe Eigenvektoren     */
/*      werden so normiert, dass ihre groesste Komponente = 1 und ihr       */
/*      zugehoeriger Imaginaerteil = 0 ist.                                 */
/*                                                                          */
/*      Literatur: Algorithmus 343 aus: Communications of the ACM, Vol 11   */
/*      No. 12, Dezember 1963, S. 820 - 826. Comment in Comm. ACM, Vol 13   */
/*      No. 2, Februar 1970, S. 122 - 124.                                  */
/*                                                                          */
/*      Methode                                                             */
/*                                                                          */
/*      First in the subroutine scale the matrix is scaled so that the      */
/*      corresponding rows and columns are approximately balanced and then  */
/*      the matrix is normalised so that the value of the euclidian norm    */
/*      of the matrix is equal to one.                                      */
/*                                                                          */
/*      the eigenvalues are computed by the QR double step method in the    */
/*      subroutine HESQR                                                    */
/*                                                                          */
/*      the eigenvectors are computed by inverse iteration in the           */
/*      subroutine REALVE, for the real eigenvalues, or in the subroutine   */
/*      COMPVE for the complex eigenvalues.                                 */
/*                                                                          */

int eigen1(int n,double *a,double *evr,double *evi,
    double *vecr,double *veci,int *indic)
{
    int i,j,k,k1;
    int err,l,l1,m,kon,ivec,iworka,locala,prfacta,subdiaa,worka,work1a,work2a;
    int *iwork,*local;
    double d1,d2,d3,r,r1,eps,enorm,ex;
    double *prfact,*subdia,*work,*work1,*work2;

    if (n == 1) {
        evr[1] = a[1]; evi[1] = 0.0; vecr[1] = 1.0; veci[1]  = 0.0;
        indic[1] = 2;
        err = 0;
        goto fin;
    }
    err = -1;

    /* temporary storage allocation */

    iworka = locala = prfacta = subdiaa = worka = work1a = work2a = 0;

    if (!(iwork = (int *)calloc(n + 1,sizeof(int))))
        goto fin;   
    memrq(n + 1,sizeof(int));
    iworka = n + 1;

    if (!(local = (int *)calloc(n + 1,sizeof(int))))
        goto fin;   
    memrq(n + 1,sizeof(int));
    locala = n + 1;

    if (!(prfact = (double *)calloc(n + 1,sizeof(double))))
        goto fin;   
    memrq(n + 1,sizeof(double));
    prfacta = n + 1;

    if (!(subdia = (double *)calloc(n + 1,sizeof(double))))
        goto fin;  
    memrq(n + 1,sizeof(double));
    subdiaa = n + 1;

    if (!(work = (double *)calloc(n + 1,sizeof(double))))
        goto fin;  
    memrq(n + 1,sizeof(double));
    worka = n + 1;

    if (!(work1 = (double *)calloc(n + 1,sizeof(double))))
        goto fin;  
    memrq(n + 1,sizeof(double));
    work1a = n + 1;

    if (!(work2 = (double *)calloc(n + 1,sizeof(double))))
        goto fin;   
    memrq(n + 1,sizeof(double));
    work2a = n + 1;

    err = 0;

    enorm = scale(n,a,veci,prfact);

    /*  the computation of the eigenvalues of the normalised matrix  */

    /****
    t = 53.0;     Anzahl signifikante Bits     
    ex = exp(-t * log(2.0));        
    printf1("\nex=%24.18lf",ex);
    ****/

    ex = EPSI;

    eps = hesqr(n,a,veci,evr,evi,subdia,indic,ex);

    /* the possible decomposition of the upper hessenberg matrix into the
    the submatrices of lower order is indicated in the array local. the
    decomposition occurs when some subdiagonal elements are in modulus
    less then a small positive number eps defined in the subroutine
    hesqr. the amount of work in the eigenvector problem may be
    diminished in this way. */
 
    j = n;
    i = 1;
    local[1] = 1;
    while (j > 1) {
        if (fabs(subdia[j - 1]) <= eps) {
            ++i; local[i] = 0;
        }
        local[i] = local[i] + 1;
        --j;
    }
 
    /* the eigenvector problem */
 
    k = 1;
    kon = 0;
    l = local[1];
    m = n;
    for (i = 1; i <= n; ++i) {
        ivec = n - i + 1;
        if (i > l) {
            ++k; m = n - l; l += local[k];
        }
        if (indic[ivec] != 0) {
            if (evi[ivec] == 0.0) {    
 
                /* transfer of an upper-hessenberg matrix of the order m from the
                arrays veci and subdia into the array a. */
 
                for (k1 = 1; k1 <= m; ++k1) {
                    for (l1 = k1; l1 <= m; ++l1) {
                        a[(k1 - 1) * n + l1] = veci[(k1 - 1) * n + l1];
                    }
                    if (k1 != 1)
                        a[(k1 - 1) * n + k1 - 1] = subdia[k1 - 1];
                }
 
                /* the computation of the real eigenvector ivec of the upper-
                hessenberg matrix corresponding to the real eigenvalue evr[ivec] */

                realve(n,m,ivec,a,vecr,evr,evi,iwork,work,indic,eps,ex);
            }
            else {

            /* the computation of the complex eigenvector ivec of the upper-
            hessenberg matrix corresponding to the complex eigenvalue
            evr[ivec]+i*evi[ivec]. if the value of kon is not equal to zero
            then this complex eigenvector has already been found from its
            conjugate. */

                if (kon != 0)
                    kon = 0;
                else {
                    kon = 1;
                    compve(n,m,ivec,a,vecr,veci,evr,evi,indic,iwork,subdia,
                          work1,work2,work,eps,ex);     
                }
            }
        }
    }

    /* the reconstruction of the matrix used in the reduction of matrix a
    to an upper hessenberg form by householder method. */
 
    for (i = 1; i <= n; ++i) {
        for (j = i; j <= n; ++j)  
            a[(i - 1) * n + j] = a[(j - 1) * n + i] = 0.0;
        a[(i - 1) * n + i] = 1.0;
    }
    if (n > 2) {
        m = n - 2;
        for (k = 1; k <= m; ++k) {
            l = k + 1;
            for (j = 2; j <= n; ++j) {
                d1 = 0.0;
                for (i = l; i <= n; ++i) {
                   d2 = veci[(i - 1) * n + k]; d1 += d2 * a[(j - 1) * n + i];
                }
                for (i = l; i <= n; ++i)   
                    a[(j - 1) * n + i] -= veci[(i - 1) * n + k] * d1;
            }
        }
    }
 
    /* the computation of the eigenvectors of the original nonscaled matrix */
 
    kon = 1;
    for (i = 1; i <= n; ++i) {
        l = 0;
        if (evi[i] == 0.0)
            goto lab1016;     

        l = 1;
        if (kon == 0) goto lab1016;
        kon = 0;
        goto lab1024;
lab1016:
        for (j = 1; j <= n; ++j) {
            d1 = d2 = 0.0;
            for (k = 1; k <= n; ++k) {
                d3 = a[(j - 1) * n + k];
                d1 = d1 + d3 * vecr[(k - 1) * n + i];
                if (l != 0) 
                    d2 += d3 * vecr[(k - 1) * n + i - 1];
            }
            if (prfact[j] == 0.0) {
                printfe("\nerror in eigen1: prfact(j) = %lg\n",prfact[j]);
                gerr_exit(216);
            }
            work[j] = d1 / prfact[j];

            if (l != 0)
                subdia[j] = d2 / prfact[j];
        }

        /* the normalisation of the eigenvectors and the computation of the
        eigenvalues of the original non-scaled matrix. */

        if (l == 1)
            goto lab1021;
        d1 = 0.0;
        for (m = 1; m <= n; ++m)   
            d1 += work[m] * work[m];

        if (d1 == 0.0) {
            printfe("\nerror in eigen1: d1 = %24.18e m=%d\n",d1,m);
            gerr_exit(217);
        }
        d1 = sqrt(d1);

        for (m = 1; m <= n; ++m) {
            veci[(m - 1) * n + i] = 0.0;
            vecr[(m - 1) * n + i] = work[m] / d1;
        }
        evr[i] *= enorm;
        goto lab1024;
lab1021:
        kon = 1;
        evr[i] *= enorm; evr[i - 1] =  evr[i];
        evi[i] *= enorm; evi[i - 1] = -evi[i];
        r = 0.0;
        for (j = 1; j <= n; ++j) {
            r1 = work[j] * work[j] + subdia[j] * subdia[j];
            if (r < r1) {
                r = r1; l = j;
            }
        }
        d3 = work[l];
        r1 = subdia[l];

        for (j = 1; j <= n; ++j) {
            if (r == 0) {
                printfe("\nerror in eigen1: r = %20.16f\n",r);
                gerr_exit(218);
            }
            d1 = work[j];
            d2 = subdia[j];
            vecr[(j - 1) * n + i] = (d1 * d3 + d2 * r1) / r;
            veci[(j - 1) * n + i] = (d2 * d3 - d1 * r1) / r;
            vecr[(j - 1) * n + i - 1] =  vecr[(j - 1) * n + i];
            veci[(j - 1) * n + i - 1] = -veci[(j - 1) * n + i];
        }
lab1024: ;
    }

fin:
    if (iworka > 0) {
        free((char *)iwork);
        memrq(-iworka,sizeof(int));
    }
    if (locala > 0) {
        free((char *)local);
        memrq(-locala,sizeof(int));
    }
    if (prfacta > 0) {
        free((char *)prfact);
        memrq(-prfacta,sizeof(double));
    }
    if (subdiaa > 0) {
        free((char *)subdia);
        memrq(-subdiaa,sizeof(double));
    }
    if (worka > 0) {
        free((char *)work);
        memrq(-worka,sizeof(double));
    }
    if (work1a > 0) {
        free((char *)work1);
        memrq(-work1a,sizeof(double));
    }
    if (work2a > 0) {
        free((char *)work2);
        memrq(-work2a,sizeof(double));
    }
    if (err)
        return(err);
   
    for (i = 1; i <= n; ++i) 
        if (indic[i] != 2) 
            return(1);
    return(0);
}
 
/*  scale

    this subroutine stores the matrix a of the order n from the
    array a into the array h. afterward the matrix in the array a
    is scaled so that the quotient of the absolute sum of the
    off-diagonal elements of column i and the absolute sum of the
    off-diagonal elements of row i lies within the values of bound1
    and bound2.
 
    the component i of the eigenvector obtained by using the scaled
    matrix must be divided by the value found in the prfact[i] of 
    the array prfact. in this way the eigenvector of the non-
    scaled matrix is obtained.
    
    after the matrix is scaled it is normalised so that the value
    of the euclidian norm is equal to one. if the process of scaling
    was not succesful the original matrix from the array h would be
    stored back into a and the eigenprobelem would be solved by
    using this matrix.
 
    the eigenvalues of the normalised matrix must be multiplied by
    the scalar enorm in order that they cecome the eigenvalues of
    the non-normalised matrix. */
     
double scale(int n,double *a,double *veci,double *prfact)
{
    int i,iter,j,ncount;
    double enorm,bound1,bound2,column,row,factor,fnorm,q,fabs(),sqrt();

    for (i = 1; i <= n; ++i) {
        for (j = 1; j <= n; ++j)
            veci[(i - 1) * n + j] = a[(i - 1) * n + j];
        prfact[i] = 1.0;
    }
    bound1 = 0.75; bound2 = 1.33; iter = ncount = 0;

    while (ncount < n) {
        for (i = 1; i <= n; ++i) {
            column = row = 0.0;
            for (j = 1; j <= n; ++j) {
                if (i != j) {
                    column += fabs(a[(j - 1) * n + i]);
                    row += fabs(a[(i - 1) * n + j]);
                }
            }
            if ((column == 0.0) || (row == 0.0)) ++ncount;
            else {
                q = column / row;
                if ((q >= bound1) && (q <= bound2)) ++ncount;
                else {

                    /***************************************************/
                    if (q == 0.0) {
                        printfe("\nerror in eigen1 (scale): q = %20.18lf\n",q);
                        gerr_exit(219);
                    }
                    /***************************************************/

                    factor = sqrt(q);
                    for (j = 1; j <= n; ++j) {
                        if (i != j) {
                            a[(i - 1) * n + j] *= factor;
                            a[(j - 1) * n + i] /= factor;
                        }
                    }
                    prfact[i] *= factor;
                }
            }
        }
        if (++iter > 30)
            break;
    }
    if (iter <= 30) {
        fnorm = 0.0;

        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) {
                if (i != j)
                    a[(i - 1) * n + j] = veci[(i - 1) * n + j] * prfact[i] /
                                                                 prfact[j];
                q = a[(i - 1) * n + j];
                fnorm += q * q;
            }
        }

        /**********************************************************/
        if (fnorm == 0.0) {
            printfe("\nerror in eigen1 (scale): fnorm = %20.18lf\n",fnorm);
            gerr_exit(220);
        }
        /**********************************************************/

        fnorm = sqrt(fnorm);
        for (i = 1; i <= n; ++i) {
            for (j = 1; j <= n; ++j) 
                a[(i - 1) * n + j] /= fnorm;
        }
        enorm = fnorm;
    }
    else {

        /*******************************************
        printf1("\neigen1 (scale): iter = %d\n",iter);
        *******************************************/

        for (i = 1; i <= n; ++i) {
            prfact[i] = 1.0;
            for (j = 1; j <= n; ++j)
                a[(i - 1) * n + j] = veci[(i - 1) * n + j];
        }
        enorm = 1.0;
    }
    return(enorm);
}
 
 
/*  hesqr

    this subroutine finds all eigenvalues of a real general matrix.
    the original matrix a of order n is reduced to the upper 
    hessenberg form h by means of similarity transformations
    (householder method). the matrix h is preserved in the upper
    half of the array h and in the array subdia. the special
    vectors used in the definition of the householder transformation
    matrices are stored in the lower part of the array h.

    the real parts of the n eigenvalues will be found in the first
    n places of the array evr and the imaginary parts in the first
    n places of the array evi. the array indic indicates the succes
    of the routine as follows:

    value of indic[i]    eigenvalue i
       0                   not found
       1                   found
  
    eps is a small positive number that numerically represents zero
    in the program. eps = (euclidian norm of h) * ex, where ex =
    2 ** (-t). t is the number of binary digits in the mantissa of
    a floating point number. */
      
double hesqr(int n,double *a,double *veci,double *evr,double *evi,
    double *subdia,int *indic,double ex)
{
    int i,j,k,l,m,m1,maxst,ns;
    double eps,r,s,sr,sr2,shift,t,x,y,z;

    /* reduction of the matrix a to an upper hessenberg form h. */
 
    if (n == 2) {
        subdia[1] = a[n + 1];
    }
    else if (n > 2) {
        m = n - 2;
        for (k = 1; k <= m; ++k) {
            l = k + 1; s = 0.0;
            for (i = l; i <= n; ++i) {
                veci[(i - 1) * n + k] = a[(i - 1) * n + k];
                s += fabs(a[(i - 1) * n + k]);
            }
            if (fabs(s - fabs(a[k * n + k])) == 0.0) {
                subdia[k] = a[k * n + k]; veci[k * n + k] = 0.0;
            }
            else {
                sr2 = 0.0;
                for (i = l; i <= n; ++i) {
                    a[(i - 1) * n + k] = sr = a[(i - 1) * n + k] / s;
                    sr2 = sr2 + sr * sr;
                }
                sr = sqrt(sr2);
                if (a[(l - 1) * n + k] >= 0.0) sr = -sr;
                sr2 = sr2 - sr * a[(l - 1) * n + k];
                a[(l - 1) * n + k] -= sr;
                veci[(l - 1) * n + k] -= sr * s;
                subdia[k] = sr * s;

                /**********************************************************/
                if (sr2 == 0.0) {
                    printfe("\nerror in eigen1 (hesqr): sr2 = %20.18lf\n",sr2);
                    gerr_exit(221);
                }
                /**********************************************************/

                x = s * sqrt(sr2);
                for (i = l; i <= n; ++i) {
                    veci[(i - 1) * n + k] /= x;
                    subdia[i] = a[(i - 1) * n + k] / sr2;
                }

                /* premultiplication by the matrix pr. */
 
                for (j = l; j <= n; ++j) {
                    sr = 0.0;
                    for (i = l; i <= n; ++i)
                        sr += a[(i - 1) * n + k] * a[(i - 1) * n + j];
                    for (i = l; i <= n; ++i)
                        a[(i - 1) * n + j] -= subdia[i] * sr;
                }
 
                /* postmultiplication by the matrix pr. */
 
                for (j = 1; j <= n; ++j) {
                    sr = 0.0;
                    for (i = l; i <= n; ++i)
                        sr += a[(j - 1) * n + i] * a[(i - 1) * n + k];
                    for (i = l; i <= n; ++i)
                        a[(j - 1) * n + i] -= subdia[i] * sr;
                }
            }
        }
        for (k=1; k<=m; ++k)
            a[k * n + k] = subdia[k];
 
        /* transfer of the upper half of the matrix a into the array h and the
        calculation of the small positive number eps. */
 
        subdia[n-1] = a[(n - 1) * n + n - 1];
    }
    eps = 0.0;
    for (k = 1; k <= n; ++k) {
        indic[k] = 0;
        if (k != n)
            eps = eps + subdia[k] * subdia[k];
        for (i = k; i <= n; ++i) {
            veci[(k - 1) * n + i] = a[(k - 1) * n + i];
            eps += a[(k - 1) * n + i] * a[(k - 1) * n + i];
        }
    }

    /********************************************************** 
    if (eps == 0.0) {
        printfe("\nerror in eigen1 (hesqr): eps = %20.18lf\n",eps);
        gerr_exit(222);
    }
    **********************************************************/

    eps = ex * sqrt(eps);
 
    /* the QR iterative process. the upper hessenberg matrix h is
    reduced to the upper modified triangular form.
 
    determination of the shift of origin for the first step of the
    QR iterative process. */
 
    shift = a[(n - 1) * n + n - 1];

    if (n <= 2 || a[(n - 1) * n + n] != 0.0 ||
                  a[(n - 2) * n + n] != 0.0 ||
                  a[(n - 2) * n + n - 1] != 0.0)
        shift = 0.0;

    m = n; ns = 0; maxst = n * 10;
 
    /* testing if the upper half of the matrix is equal to zero. if it is equal
    to zero the upper QR process is not necessary. */   
 
    for (i = 2; i <= n; ++i)
        for (k = i; k <= n; ++k)
            if (a[(i - 2) * n + k] != 0.0)
                goto lab3018;

    for (i = 1; i <= n; ++i) {
        indic[i] = 1;
        evr[i]   = a[(i - 1) * n + i];
        evi[i]   = 0.0;
    }
    return(eps);
 
    /* start the main loop of the QR process */
 
lab3018:
    k = m - 1; i = m1 = k;
 
    /* find any decompositions of the matrix. j*/
 
    if (k < 0) return(eps);

    if (k == 0) goto lab3034;

    /* compute the last eigenvalue */

    if (fabs(a[(m - 1) * n + k]) <= eps) goto lab3034;
    if (m == 2) goto lab3035;

    while (k > 0) {
        i--;
        if (fabs(a[(k - 1) * n + i]) <= eps) break;
        k = i;
    }
lab3021:
    if (k == m1) goto lab3035;
 
    /* transformation of the matrix of the order greater than two. */
 
    s = a[(m - 1) * n + m] + a[(m1 - 1) * n + m1] + shift;
    sr = a[(m - 1) * n + m ] * a[(m1 - 1) * n + m1] -
         a[(m - 1) * n + m1] * a[(m1 - 1) * n + m ] + 0.25 * shift * shift;
    a[(k + 1) * n + k] = 0.0;
 
    /* calculate x1,y1,z1 for the submatrix obtained by the decomposition */
 
    x = a[(k - 1) * n + k] * (a[(k - 1) * n + k] - s) +
        a[(k - 1) * n + k + 1] * a[k * n + k] + sr;
    y = a[k * n + k] * (a[(k - 1) * n + k] + a[k * n + k + 1] - s);
    r = fabs(x) + fabs(y);

    if (r / eps - eps / ex <= 0.0) {
        if (shift - a[(m - 1) * n + m - 1] != 0.0) {
            shift = a[(m - 1) * n + m - 1];
            goto lab3021;  
        }
    }
    z = a[(k + 1) * n + k + 1] * a[k * n + k];

    /***
    if (r == 0.0) {
        shift = a[(m - 1) * n + m - 1];
        goto lab3021;  
    }
    z = a[(k + 1) * n + k + 1] * a[k * n + k];
    ***/

    shift = 0.0;
    ns++;

    /* the loop for one step of the QR process. */
 
    for (i=k; i<=m1; ++i) {
        if (i != k) {
            x = a[(i - 1) * n + i - 1];
            y = a[ i * n + i - 1];
            z = 0.0;
            if (i + 2 <= m)
                z = a[(i + 1) * n + i - 1];
        }     
        sr2 = fabs(x) + fabs(y) + fabs(z);

        if (sr2 != 0.0) {
            x = x / sr2;
            y = y / sr2;
            z = z / sr2;
        }
        s = sqrt(x*x + y*y + z*z);   
        if (x >= 0.0)
            s = -s;
        if (i != k)
            a[(i - 1) * n + i - 1] = s * sr2;

        if (sr2 == 0.0) {  
            if (i + 3 > m)
                goto lab3033;
            else
                goto lab3032;
        }
        sr = 1.0 - x / s;
        s  = x - s;
        x  = y / s;
        y  = z / s;

        /* premultiplication by the matrix pr. */
 
        for (j = i; j <= m; ++j) {
            s = a[(i - 1) * n + j] + a[i * n + j] * x;
            if (i + 2 <= m)
                s += a[(i + 1) * n + j] * y;
            s = s * sr;
            a[(i - 1) * n + j] -= s;
            a[i * n + j] -= s * x;
            if (i + 2 <= m)
                a[(i + 1) * n + j] -= s * y;
        } 
 
        /* postmultiplication by the matrix pr. */
 
        l = i + 2;
        if (i >= m1)
            l = m;
        for (j = k; j <= l; ++j) {
            s = a[(j - 1) * n + i] + a[(j - 1) * n + i + 1] * x;
            if (i + 2 <= m)
                s += a[(j - 1) * n + i + 2] * y;
            s = s * sr;
            a[(j - 1) * n + i] -= s;
            a[(j - 1) * n + i + 1] -= s * x;
            if (i + 2 <= m)
                a[(j - 1) * n + i + 2] -= s * y;
        }
        if (i + 3 > m)   
            goto lab3033;
        s = -a[(i + 2) * n + i + 2] * y * sr;
lab3032:
        a[(i + 2) * n + i] = s;
        a[(i + 2) * n + i + 1] = s * x;
        a[(i + 2) * n + i + 2] = s * y + a[(i + 2) * n + i + 2];
lab3033: ;
    }  

    if (ns <= maxst) goto lab3018;
    else
        return(eps);

lab3034:
    evr[m] = a[(m - 1) * n + m];
    evi[m] = 0.0;
    indic[m] = 1;
    m = k;
    goto lab3018;
 
lab3035:

    /* compute the eigenvalues of the last 2 x 2-matrix obtained by the
    decomposition */
 
    r = 0.5 * (a[(k - 1) * n + k] + a[(m - 1) * n + m]);
    s = 0.5 * (a[(m - 1) * n + m] - a[(k - 1) * n + k]);
    s = s * s + a[(k - 1) * n + m] * a[(m - 1) * n + k];
    indic[k] = 1;
    indic[m] = 1;
    if (s >= 0.0) {
        t = sqrt(s);
        evr[k] = r - t; evr[m] = r + t; evi[k] = 0.0; evi[m] = 0.0;
    }  
    else {
        t = sqrt(-s);
        evr[k] = r; evi[k] = t; evr[m] = r; evi[m] = -t;
    }
    m -= 2;
    goto lab3018;
}


/*  realve

    this subroutine finds the real eigenvectors of the real upper-
    hessenberg matrix in the array a, corresponding to the real
    eigenvalues stored in evr[ivec]. the inverse iteration method
    is used.

    the matrix in a is destroyed by the subroutine.
 
    n is the order of the upper hessenberg matrix. m is the order
    of the submatrix obtained by a suitable decomposition of the
    upper hessenberg matrix if some subdiagonal elements are equal
    to zero. the value of m is chosen so that the last n-m
    components of the eigenvector are zero.
 
    ivec gives the position of the eigenvalue in the array evr
    for which the corresponding eigenvector is computed.

    the array evi would contain the imaginary parts of the n eigen-
    values if they existed.
 
    the m components of the computed real eigenvector will be found
    in the first m places of the column ivec of the two dimensional
    array vecr.

    iwork and work are the working stores used during the gaussian
    elimination and backsubstitution process.
  
    the array indic indicates the success of the routine as follows
         value of indic[i]    eigenvector i
              1                  not found
              2                  found

    eps is a small positive number that numerically represents
    zero in the program. eps = (euclidian norm of a) * ex, where
    ex = 2 ** (-t). t is the number of binary digits in the mantissa
    of a floating point number. */
     
int realve(int n,int m,int ivec,double *a,double *vecr,double *evr,double *evi,
    int *iwork,double *work,int *indic,double eps,double ex)
{
    int i,j,k,l,iter,ns;
    double evalue,bound,r,r1,s,sr,t,previs;
      
    vecr[ivec] = 1.0;
    if (m == 1) {
        indic[ivec] = 2;
        for (i = m + 1; i <= n; ++i)
            vecr[(i - 1) * n + ivec] = 0.0;
        return(0);
    }

    /* small pertubation of equal eigenvalues to obtain a full set of
    eigenvectors */

    evalue = evr[ivec];

    if (ivec != m) {
        k = ivec + 1;
        r = 0.0;
        for (i = k; i <= m; ++i) {

            if (fabs(evalue - evr[i]) == 0.0 && fabs(evi[i]) == 0.0)
                r += 3.0;
        }
        evalue += r * ex;
    }
    for (k = 1; k <= m; ++k)
        a[(k - 1) * n + k] -= evalue;
    
    /* gaussian elimination of the upper hessenberg matrix a. all row
    interchanges are indicated in the array iwork. all the multipliers
    are stored as the subdiagonal elements of a. */
 
    k = m - 1;
    for (i = 1; i <= k; ++i) {
        l = i + 1;
        iwork[i] = 0;

        if (a[i * n + i] != 0.0) {
            if (fabs(a[(i - 1) * n + i]) < fabs(a[i * n + i])) {
                iwork[i] = 1;
                for (j = i; j <= m; ++j) {
                    r = a[(i - 1) * n + j];
                    a[(i - 1) * n + j] = a[i * n + j];
                    a[i * n + j] = r;
                }
            }
            r = -a[i * n + i] / a[(i - 1) * n + i];
            a[i * n + i] = r;
            for (j = l; j <= m; ++j)  
                a[i * n + j] += r * a[(i - 1) * n + j];
        }
        else {
            if (a[(i - 1) * n + i] == 0.0) {
                a[(i - 1) * n + i] = eps;
            }
        }
    }
    if (a[(m - 1) * n + m] == 0.0) 
        a[(m - 1) * n + m] = eps;

    /* the vector (1,1,...,1) is stored in the place of the right hand
    side column vector */
 
    for (i=1; i<=n; ++i) {
        if (i > m) work[i] = 0.0; else work[i] = 1.0;
    }
 
    /* the inverse iteration is performed on the matrix until the infinite
    norm of the right hand side vector is greater than the bound 
    defined as 0.01 / (n * ex). */
      
    bound = 0.01 / (ex * (double)n);
    ns = 0;
    iter = 1;
 
    /* the backsubstitution */
 
    while (1) {
        r = 0.0;
        for (i = 1; i <= m; ++i) {
            j = m - i + 1;
            s = work[j];
            if (j != m) {
                l = j + 1;
                for (k = l; k <= m; ++k) {
                    sr = work[k];
                    s -= sr * a[(j - 1) * n + k];
                }
            }
            work[j] = s / a[(j - 1) * n + j];
            t = fabs(work[j]);
            if (r < t)   
                r = t;
        }
   
        /* the computation of the right hand side vector for the new
        iteration step */

        for (i=1; i<=m; ++i)
            work[i] /= r;   
   
        /* the computation of the residuals and comparison of the residuals
        of the two successive steps of the inverse iteration. if the infinite
        norm of the residual vector is greater than the infinite norm of
        the previous residual vector the computed eigenvector of the 
        previous step is taken as the final eigenvector. */
 
        r1 = 0.0;
        for (i = 1; i <= m; ++i) {
            t = 0.0;
            for (j = i; j <= m; ++j)
                t += a[(i - 1) * n + j] * work[j];
            t = fabs(t);
            if (r1 < t)
                r1 = t;
        }
        if (iter != 1) {
            if (previs <= r1)
                break;
        }
        for (i = 1; i <= m; ++i)
            vecr[(i - 1) * n + ivec] = work[i];
        previs = r1;
        if (ns == 1)
            break;
        if (iter > 6)
            goto lab4025;
        ++iter;

        if (r >= bound)
            ns = 1;
    
        /* gaussian elimination of the right hand side vector */
 
        k = m - 1;
        for (i = 1; i <= k; ++i) {
            r = work[i + 1];
            if (iwork[i] == 0)
                work[i + 1] += work[i] * a[i * n + i];
            else {
                work[i + 1] = work[i] + work[i + 1] * a[i * n + i];
                work[i] = r;
            }
        }
    }
    indic[ivec] = 2;
lab4025:
    for (i = m + 1; i <= n; ++i)
        vecr[(i - 1) * n + ivec] = 0.0;
    return(0);
}
 
/*  compve
  
    this subroutine finds the complex eigenvector of the real upper
    hessenberg matrix of order n corresponding to the complex
    eigenvalue with the real part in evr[ivec] and the imaginary
    part in evi[ivec]. the inverse iteration method is used,
    modified to avoid the use of complex arithmetic.
  
    the matrix on which the inverse iteration is performed is built
    up in the array a by using the upper hessenberg matrix 
    preserved in the upper half of the array h and in the array subdia.
  
    m is the order of the submatrix obtained by a suitable 
    decomposition of the upper hessenberg matrix if some subdiagonal
    elements are equal to zero. the value of m is chosen so that
    the last n-m components of the complex eigenvector are zero.
  
    the real parts of the first m components of the computed
    complex eigenvector will be found in the first m places of
    the column whose top element is vecr[1,ivec] and the corresponding
    imaginary parts of the first m components of the complex eigen-
    vector will be found in the first m places of the column whose
    top element is vecr[i,ivec-1].
  
    the array indic indicates the success of the routine as follows
          value of indic[i]    eigenvector i
                1                not found
                2                found
  
    the arrays iwork, work1, work2, work are the working stores used
    during the inverse iteration process. eps and ex are defined
    as in realve. */
    
int compve(int n,int m,int ivec,double *a,double *vecr,double *veci,double *evr,
    double *evi,int *indic,int *iwork,double *subdia,double *work1,
    double *work2,double *work,double eps,double ex)
{
    int i,i1,i2,j,k,l,iter,ns;
    double  b,bound,d,d1,eta1,fksi,r,s,u,v,previs;

    fksi = evr[ivec];
    eta1  = evi[ivec];
 
    /* the modification of the eigenvalue (fksi + i * eta1) if more eigenvalues
    are equal */
 
    if (ivec != m) {
        k = ivec + 1;
        r = 0.0;
        for (i = k; i <= m; ++i) {
            if (fabs(fksi - evr[i]) == 0.0 &&
                fabs(fabs(eta1) - fabs(evi[i])) == 0.0)
                r += 3.0;
        }
        r *= ex; fksi += r; eta1 += r;
    }

    /* the matrix ((h - fksi * i) * (h - fksi * i) + (eta1 * eta1) * i)
    is stored in a */
 
    r = fksi * fksi + eta1 * eta1; s = 2.0 * fksi; l = m - 1;

    for (i = 1; i <= m; ++i) {
        for (j = i; j <= m; ++j) {
            a[(j - 1) * n + i] = d = 0.0;
            for (k = i; k <= j; ++k)
                d += veci[(i - 1) * n + k] * veci[(k - 1) * n + j];
            a[(i - 1) * n + j] = d - s * veci[(i - 1) * n + j];
        }
        a[(i - 1) * n + i] += r;
    }
    for (i = 1; i <= l; ++i) {
        r = subdia[i];
        a[i * n + i] = -s * r;
        i1 = i + 1;
        for (j = 1; j <= i1; ++j)
            a[(j - 1) * n + i] += r * veci[(j - 1) * n + i + 1];
        if (i != 1)
            a[i * n + i - 1] = r * subdia[i - 1];
        for (j = i; j <= m; ++j)
            a[i * n + j] += r * veci[(i - 1) * n + j];
    }
 
    /* the gaussian elimination of the matrix ((h - fksi * I) * (h - fksi * i)+
    (eta1 * eta1) * i) in the array a. the row interchanges that occur
    are indicated in the array iwork. all the multipliers are stored
    in the first and in the second subdiagonal of the array a. */
 
    k = m - 1;
    for (i = 1; i <= k; ++i) {
        i1 = i + 1;
        i2 = i + 2;
        iwork[i] = 0;
        if (i == k) goto lab5010;
        if (a[(i + 1) * n + i] != 0.0) goto lab5011;

lab5010:
        if (a[i * n + i] != 0.0) goto lab5011;
        if (a[(i - 1) * n + i] != 0.0) goto lab5018;
        a[(i - 1) * n + i] = eps;
        goto lab5018;

lab5011:
        if (i == k) goto lab5012;
        if (fabs(a[i * n + i]) - fabs(a[(i + 1) * n + i]) >= 0.0)
            goto lab5012;

        if (fabs(a[(i - 1) * n + i]) - fabs(a[(i + 1) * n + i]) >= 0.0)
            goto lab5016;

        l = i + 2;
        iwork[i] = 2;
        goto lab5013;

lab5012:
        if (fabs(a[(i - 1) * n + i]) - fabs(a[i * n + i]) >= 0.0)
            goto lab5015;

        l = i + 1;
        iwork[i] = 1;

lab5013:
        for (j = i; j <= m; ++j) {
            r = a[(i - 1) * n + j];
            a[(i - 1) * n + j] = a[(l - 1) * n + j];
            a[(l - 1) * n + j] = r;
        }
lab5015:
        if (i == k) i2 = i1;
lab5016:
        for (l = i1; l <= i2; ++l) {

            /*******************************************************************/
            r = -a[(l - 1) * n + i] / a[(i - 1) * n + i];
            a[(l - 1) * n + i] = r;
            for (j = i1; j <= m; ++j)
                a[(l - 1) * n + j] += r * a[(i - 1) * n + j];
        }
lab5018: ;
    }

    if (a[(m - 1) * n + m] == 0.0)  
        a[(m - 1) * n + m] = eps;
 
    /* the vector (1,1...,1) is stored into the right hand side vectors
    vecr[ ,ivec] and vecr[ ,ivec-1] representing the complex right
    hand side vector */
 
    for (i = 1; i <= n; ++i) {
        if (i > m)  
            vecr[(i - 1) * n + ivec] = vecr[(i - 1) * n + ivec - 1] = 0.0;
        else  
            vecr[(i - 1) * n + ivec] = vecr[(i - 1) * n + ivec - 1] = 1.0;
    }
 
    /* the inverse iteration is performed on the matrix until the infinite
    norm of the right hand side vector is greater than the bound
    defined as 0.01 / (n * ex). */
 
    bound = 0.01 / (ex * (double)n); 
    ns = 0;
    iter = 1;
    for (i = 1; i <= m; ++i) 
        work[i] = veci[(i - 1) * n + i] - fksi;

    while (1) {
        for (i = 1; i <= m; ++i) {
            d = work[i] * vecr[(i - 1) * n + ivec];
            if (i != 1)
                d += subdia[i-1] * vecr[(i - 2) * n + ivec];
            l = i + 1;
            if (l <= m) {
                for (k = l; k <= m; ++k)
                    d += veci[(i - 1) * n + k] * vecr[(k - 1) * n + ivec];
            }
            vecr[(i - 1) * n + ivec - 1] =
                                      d - eta1 * vecr[(i - 1) * n + ivec - 1];
        }
    
        /* gaussian elimination of the right hand side vector */
 
        k = m - 1;
        for (i = 1; i <= k; ++i) {
            l = i + iwork[i];
            r = vecr[(l - 1) * n + ivec - 1];
            vecr[(l - 1) * n + ivec - 1] = vecr[(i - 1) * n + ivec - 1];
            vecr[(i - 1) * n + ivec - 1] = r;
            vecr[i * n + ivec - 1] += a[i * n + i] * r;
            if (i != k)
                vecr[(i + 1) * n + ivec - 1] += a[(i + 1) * n + i] * r;
        }
 
        /* the computation of the real part u(s+1) of the complex vector
        w(s+1). the vector u(s+1) is obtained after the backsubstitution */
 
        for (i = 1; i <= m; ++i) {
            j = m - i + 1;
            d = vecr[(j - 1) * n + ivec - 1];
            if (j != m) {
                l = j + 1;
                for (k = l; k <= m; ++k) {
                    d1 = a[(j - 1) * n + k];
                    d -= d1 * vecr[(k - 1) * n + ivec - 1];
                }
            }

            /***************************************************************/

            vecr[(j - 1) * n + ivec - 1] = d / a[(j - 1) * n + j];
        }
    
        /* the computation of the imaginary part v(s+1) of the vector
        w(s+1), where v(s+1) = (p(s) - (a - fksi * i) * u(s+1)) / eta1. */
 
        for (i = 1; i <= m; ++i) {
            d = work[i] * vecr[(i - 1) * n + ivec - 1];
            if (i != 1)
                d += subdia[i - 1] * vecr[(i - 2) * n + ivec - 1];
            l = i + 1;
            if (l <= m) {
                for (k = l; k <= m; ++k)
                    d += veci[(i - 1) * n + k] * vecr[(k - 1) * n + ivec - 1];
            }

            /***************************************************************/

            vecr[(i - 1) * n + ivec] = (vecr[(i - 1) * n + ivec] - d) / eta1;
        }

        /* the computation of (infin.norm of w(s+1))**2 */
 
        l = 1;
        s = 0.0;
        for (i = 1; i <= m; ++i) {
            r  = vecr[(i - 1) * n + ivec] * vecr[(i - 1) * n + ivec];
            r += vecr[(i - 1) * n + ivec - 1] * vecr[(i - 1) * n + ivec - 1];
            if (r - s > 0.0) {
                s = r; l = i;
            }
        }
    
        /* the computation of the vector z(s+1), where z(s+1) = w(s+1)/
        (component of w(s+1) with the largest absolute value) */
 
        u = vecr[(l - 1) * n + ivec - 1];
        v = vecr[(l - 1) * n + ivec];
        for (i = 1; i <= m; ++i) {
            b = vecr[(i - 1) * n + ivec];
            r = vecr[(i - 1) * n + ivec - 1];

            /***************************************************************/

            vecr[(i - 1) * n + ivec] = (r * u + b * v) / s;
            vecr[(i - 1) * n + ivec - 1] = (b * u - r * v) / s;
        }
    
        /* computation of the residuals and comparison of the residuals
        of the two successive steps of the inverse iteration. if the
        infinite norm of the residual vector is greater than the infinite
        norm of the previous residual vector the computed vector of the
        previous step is taken as the computed approximation to the
        eigenvector */
 
        b = 0.0;
        for (i = 1; i <= m; ++i) {
            r = work[i] * vecr[(i - 1) * n + ivec - 1]
                                         - eta1 * vecr[(i - 1) * n + ivec];
            u = work[i] * vecr[(i - 1) * n + ivec]
                                         + eta1 * vecr[(i - 1) * n + ivec - 1];
            if (i != 1) {
                r += subdia[i - 1] * vecr[(i - 2) * n + ivec - 1];
                u += subdia[i - 1] * vecr[(i - 2) * n + ivec];
            }
            l = i + 1;
            if (l <= m) {
                for (j = l; j <= m; ++j) {
                    r += veci[(i - 1) * n + j] * vecr[(j - 1) * n + ivec - 1];
                    u += veci[(i - 1) * n + j] * vecr[(j - 1) * n + ivec];
                }
            }
            u = r * r + u * u;
            if (b - u < 0.0)
                b = u;
        }
        if (iter == 1 || previs > b) {
            for (i = 1; i <= n; ++i) {
                work1[i] = vecr[(i - 1) * n + ivec];
                work2[i] = vecr[(i - 1) * n + ivec - 1];
            }
            previs = b;
            if (ns == 1 || iter > 6) break;
            ++iter;
            if (bound - sqrt(s) <= 0.0)
                ns = 1;
        }
        else {
            for (i = 1; i <= n; ++i) {
                vecr[(i - 1) * n + ivec] = work1[i];
                vecr[(i - 1) * n + ivec - 1] = work2[i];
            }
            break;
        }
    }
    if (iter <= 6)  
        indic[ivec - 1] = indic[ivec] = 2;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  prxxx(txt,n,m,x)                                                        */
/*      Prints a matrix with n,m elements. The string txt is printed first. */
/*      Used for protocol file.                                             */

void prxxx(char *txt,int n,int m,double *x)
{
    register int i,j;

    printf1("%s\n",txt); 
    for (i = 0; i < n; ++i) {
        for (j = 1; j <= m; ++j) 
            printf1("%16.12lf ",x[i * m + j]);
        printf1("\n");
    }
}

