/* t_imat.h */

#ifndef _TIMAT_H
#define _TIMAT_H

/*  functions in t_imat.c */

int m_midf(char *cmd);
int m_midf1(char *cmd);
int m_midf2(char *cmd);
int m_midf3(char *cmd);
int gmin(void);
int range(void);
int idf(void);
int sddf(void);
int iddf(void);
int imean(void);
int ivar(void);
int igini(void);
int ivar1(void);
int alloc_par(int n,int opt);
int alloc_list(int n,int m);
int igmin(int opt,int narg,int nbmax,double *par,double *lb,double *ub,int mxit,
    double tolbw,double tolfd,double tolfe,int gc,int typ);
int igmin_res(double tolfe,int opt);
void igmin_cpar(int n,double *x,double *par);

void i_add(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_sub(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_mul(double xl,double xh,double yl,double yh,double *rl,double *rh);
int  i_div(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_abs(double xl,double xh,double *rl,double *rh);
void i_max(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_min(double xl,double xh,double yl,double yh,double *rl,double *rh);
void i_square(double xl,double xh,double *rl,double *rh);
void i_sqrt(double xl,double xh,double *rl,double *rh);
double ivarf(int n,double x,double *xl,double *xh);


extern double TOLBW;        /* tolerance for box length, branch and bound   */
extern double TOLBC;        /* tolerance for constraint function boxes      */
extern double TOLFD;        /* tolerance for function range                 */
extern double TOLFE;        /* tolerance for global minimum                 */

extern double *RPar;        /* parameters                                   */
extern double *RParL;       /* lower bounds                                 */
extern double *RParU;       /* upper bounds                                 */

extern double *GO_LB;       /* box list                                     */
extern double *GO_UB; 
extern double *GO_LBF; 
extern double *GO_UBF; 
extern short  *GO_FLG;

extern int GO_IT;           /* number of iterations performed               */
extern int GO_FN;           /* number of function evaluations               */
extern int GO_IFN;          /* number of inclusion function evaluations     */
extern int GO_NB;           /* number of boxes in final list                */
extern int GO_NA;           /* number of accepted boxes                     */
extern int GO_NBU;          /* number of boxes used                         */
extern double GO_FMIN;      /* global function value                        */
extern double IGMINFUNVal;  /* argument in igmin_fun                        */

#endif /* _TIMAT_H */

/*  end of t_imat.h */















