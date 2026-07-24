/* t_matf.h */

#ifndef _TMATF_H
#define _TMATF_H

/*  functions in t_matf.c */

void mp_info(void);                        
int mp_alloc(int t,int m,int n);                 
void mp_putvar(int i,int n,double *mat,int nv,short *vidx,int icase); 
int mp_putlog(double x);
int mp_putpar(int n,double *x);
int mp_putmpar(int m,int n,double *x);
int mp_putcov(int n,double *x);
int get_mexpr(char *exp,int *row,int *col,int *ivflg);
int eval_mexpr(int off,int row,int col,double *mat,int ivflg,
    double *mat1,int *rsel,int *csel,int idx,double dval);

extern int MPLogIdx;        /* index of PMMPLog  matrix                     */
extern int MPParIdx;        /* index of PMMPPar  matrix                     */
extern int MPCovIdx;        /* index of PMMPCov  matrix                     */
extern int MPGradIdx;       /* index of PMMPGrad matrix                     */
extern int MPResIdx;        /* index of PMMPRes  matrix                     */

extern int MPGradRow;       /* allocated rows in PMMPGrad                   */
extern int MPGradCol;       /* allocated columns in PMMPGrad                */
extern int MPResRow;        /* allocated rows in PMMPRes                    */
extern int MPResCol;        /* allocated columns in PMMPRes                 */

#endif /* _TMATF_H */

/*  end of t_matf.h */















