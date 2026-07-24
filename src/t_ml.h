/* t_ml.h */

#ifndef _TML_H
#define _TML_H

/*  functions in t_ml.c */

void set_mldef(void);
void set_mlopt(void);
int get_dsv(int n,double *par,int gmina,double *lb,double *ub,int nl);
int ml_init(int opt,int typ,int gmina);
int syminv1(int n, double *d, double *h);
int syminv2(int n, double *x);
void fn(double *x, int mod, int typ, int scal, int opt,int *err);
int checkov(void);
void prn_mlres(int typ);
void prn_ml1res(int typ);
void prvec(char *txt,int n,double *x);
void prmat(char *txt,int n,int m,double *x);
void prhess(char *txt,int n);
void prval(char *txt,double x);
void prot_init(int typ,int mod,int gmina);

/*--------------------------------------------------------------------------*/
extern double  OFMax;       /* Values used for overflow and underflows      */
extern double  OFMin;       /* checks in checkov().                         */

extern int PCovFlg;         /* set by pcov command to print cov matrix      */
extern FILE *PCovFd;        /* file descriptor for pcov command.            */

extern double FNTStat;      /* Test statistic calculated by fn()            */
extern int FNTRank;         /* Rank of hessian with FNTStat calculation     */
extern int GCTestFlg;       /* Set if gctest command                        */
extern double GCTest;       /* GC test statistic                            */
extern int GCRank;          /* Rank of hessian with GC calculation          */
extern int LMTestFlg;       /* Set if lmtest command                        */
extern double LMTest;       /* LM test statistic                            */
extern int LMRank;          /* Rank of hessian with LM calculation          */
extern int NWTest;          /* Number of wtest commands                     */
extern short **WTest;       /* Storage for WTest[i][j]                      */
extern short *WTestN;       /* Number of entries in WTest                   */
/*--------------------------------------------------------------------------*/
extern int MINA;            /* Type of minimization algorithm               */
extern int CCTyp;           /* Type of covariance matrix calculation        */
extern int CCovFlg;         /* Type set during cov matrix calculation       */
extern int MCovFlg;         /* Set during covariance matrix calculation     */
extern int CGradFlg;        /* Set when calculating gradients into matrix   */

extern int ParTyp;          /* Type of model parametrization                */
extern int NOCUsed;         /* number of cases used                         */

extern int LConv;           /* Return Code of minimization algorithms       */
extern int NumF;            /* Number of function calls                     */
extern int NumFG;           /* Number of function calls with gradient       */
extern int NumFH;           /* Number of function calls with hessian        */

extern int NParm;           /* Total number of model parameter              */
extern int NParm1;          /* Same in reduced space                        */
extern int NParmA;          /* Used to alloc Par[]                          */
extern int NParmX;          /* Used temporarily for 2S estimation           */
extern int HSiz;            /* Size of lower hessian triangle               */
extern int HSiz1;           /* Same in reduced space                        */
extern int HSizA;           /* Used to alloc Hess[]                         */
/*--------------------------------------------------------------------------*/
extern int Iter;            /* Counter for iterations                       */
extern int MxIter;          /* Max number of iterations                     */
extern int MxItFlg;         /* Set if mxit command used                     */
extern int MxIt1;           /* Max number of iterations for linear search   */
extern int MxItR;           /* Max number of iterations for random search   */

extern int Crite;           /* Convergence Criterion (1 - 6)                */
extern int CritFlg;         /* Set if user defined criterion                */
extern double TOLG;         /* Tolerance for gradient                       */
extern double TOLF;         /* Tolerance for function value changes         */
extern double TOLP;         /* Tolerance for parameter changes              */
extern double TOLV;         /* Tolerance for Simplex algorithm              */
extern double TOLS;         /* Tolerance for direct search algorithm        */
extern double TOLI;         /* Tolerance for Monte Carlo integration        */
extern double TOLSG;        /* Tolerance for scaled gradient                */
extern double TOLSP;        /* Tolerance for scaled parameter change        */
extern double TOLSGF;       /* Final scaled gradient                        */
extern double TOLSPF;       /* Final scaled parameter change                */
extern double TOLA;         /* A-convergence tolerance, used for rsearch    */

extern int CCheck;          /* Check convergence interval in Simplex alg.   */
extern int STFlg;           /* Set if function evaluation without           */

extern double AMue;         /* Mue for Armijo conditions in fmin            */
extern double SMin;         /* Minimal step size                            */
extern double SLen;         /* Step length in direct search algorithm       */
extern double SRed;         /* Step reduction in direct search algorithm    */
extern double SRLen;        /* Step length for random search algorithm      */
extern double SRRed;        /* Step length reduction for random search      */

extern double CValG;        /* Value for convergence criterion 1            */
extern double CValF;        /* Value for convergence criterion 2            */
extern double CValP;        /* Value for convergence criterion 3            */
extern double CValV;        /* Value for convergence criterion 4            */

extern int SVFlg;           /* Set if singular value calculation for final  */
                            /* hessian requested                            */
/*--------------------------------------------------------------------------*/
extern int LLTN;            /* Set for log likelihood sensitivity test      */
extern double LLTA;
/*--------------------------------------------------------------------------*/
extern double FMax0;        /* Max log-likelihood (exponential null model)  */
extern double FMax;         /* Max log-likelihood (estimated model)         */
extern double FMax1;        /* Log-likelihood (at starting values)          */
extern double FMin;         /* Value of function after minimization         */
extern double FTmp;         /* Same, but temporary                          */

extern double *Par;         /* Parameter vector (NParm)                     */
extern double *Par1;        /* Parameter vector (NParm)                     */
extern double *ParLB;       /* lower bound                                  */
extern double *ParUB;       /* upper bound                                  */

extern double *ParS;        /* Parameter vector with starting values        */
extern int ParSA;           /* if allocated                                 */
extern double *ParC;        /* Parameter vector used with constraints       */
extern int ParCA;           /* if allocated                                 */
extern double *TPar;        /* Temporary pointer for Par                    */
extern double *Grad;        /* Gradient vector                              */
extern double *GTmp;        /* Temporary for outer product                  */
extern double *Srch;        /* Search vektor                                */
extern double *Diag;        /* Diagonal of Hessian                          */
extern double *Diag1;       /* Diagonal of Hessian                          */
extern int Diag1A;          /* if allocated                                 */
extern double *WrkD;        /* Working area with NParm dimension            */
extern double *WrkH;        /* Working area with HSiz dimension             */
extern double *Hess;        /* Strict lower triangle of Hessian (HSiz)      */
extern double *Hess1;       /* Strict lower triangle of Hessian (HSiz)      */
extern int Hess1A;          /* if allocated                                 */
extern double *TSHess;      /* Full Hessian for MATyp 7 and 8.              */
extern double *TSGrad;      /* Gradient for MATyp 7 and 8.                  */
extern double *WrkE;

extern double DScal;        /* Scaling factor                               */
extern int DScalFlg;        /* Set if user defined scaling factor           */

extern int LFunc;           /* Flag: function value calculation             */
extern int LGrad;           /* Flag: calculation of gradient                */
extern int LSec;            /* Flag: calculation of Hessian                 */
extern int FN_First;        /* set for first call of fn()                   */
extern int DSVFlg;          /* set if externally supplied starting values   */

/* ------------------------------------------------------------------------ */

#endif /* _TML_H */

/*  end of t_ml.h */
