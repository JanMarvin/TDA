/* t_gdd.h */

#ifndef _TGDD_H
#define _TGDD_H

/*  functions in t_gdd.c */

int gdd(void);
int gdd_alloc(int n); 
void gdd_free(int opt);
void gdd_setgt(int typ,int gt,int gtt);
int gdd_setnd(void);
int gdd_setptr(int opt);
int gdd_setap(void);
void gdd_prot(void);
int gdd_fndni(int k);
int check_gdj(int j);
double gdd_adj(int i,int j,int gn);
int g_suc(int i,int gn);
int g_pre(int i,int gn);
int g_loop(int i,int gn);
int g_getnd(int n,int *nd);
int g_getne(int n,int *ni,int *nj);
int gdd_check(int gn,int opt);
int gdd_check2(int gn,int gn1,int opt);
int gdd_tcheck(int typ,int dir,int val);
int gdd_ei(int k);             
int gdd_ej(int k);             
double gdd_ev(int k,int gn);      
int gdd_node(int i);             
void n_info(int i,int n);
void n_info_e(void);
int gcd(void);

extern int GD_TYP;      /* type of data structure (= opt in gdd)            */
extern int GD_GT;       /* type of graph                                    */
extern int GD_GTT;      /* handling of bidirectional edges                  */
extern int GD_NG;       /* number of graphs                                 */
extern int GD_NP;       /* number of nodes                                  */
extern int GD_NDMAX;    /* max node number                                  */
extern int *GD_ND;      /* list of node numbers                             */
extern int GD_EI;       /* index of variable for starting node              */
extern int GD_EJ;       /* index of variable for ending node                */
extern int *GD_NE;      /* number of edges                                  */
extern int *GD_NL;      /* number of loops                                  */
extern int *GD_EV;      /* list of variable numbers for edge values         */
extern double *GD_EVMax;    /* max value of edges                           */
extern int *GD_APtr;    /* pointer to adjacency matrix (option 2)           */
extern int *GD_FPN;     /* number of forward links                          */
extern int **GD_FPI;    /* forward pointer: internal node numbers           */
extern int **GD_FPK;    /* forward pointer: data matrix rows                */
extern int **GD_FPK1;   /* backward pointer: undirected graph               */
extern int *GD_BPN;     /* number of backward links                         */
extern int **GD_BPI;    /* backward pointer: internal node numbers          */
extern int **GD_BPK;    /* backward pointer: data matrix rows               */
extern int GD_FPMax;    /* max value of GD_FPN[]                            */   
extern int GD_BPMax;    /* max value of GD_BPN[]                            */   
extern int GD_PERM;     /* set if gdd1 command is used                      */
extern int GD_NOC;      /* number of cases for graph data                   */
extern int *GD_EIP;     /* values of EI                                     */
extern int *GD_EJP;     /* values of EJ                                     */
extern float *GD_EP;    /* edge values                                      */

#endif /* _TGDD_H */

/* end of t_gdd.h */


