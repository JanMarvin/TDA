/* t_com.h */

#ifndef _TCOM_H
#define _TCOM_H

/*  functions in t_com.c */

int com(void);
int comb_n(int n,int *com,int first); 
int comb_nm(int n,int m,int *com,int first);
int comb_nm1(int n,int m,int *com,int first);
int perm(int n,int *com,int *p,int *d,int first);
int partit(int *k,int *p,int last);
int partit1(int n,int k,int *c,int *d,int last);
int pcyc(void);
int indep(void);

#endif /* _TCOM_H */

/* end of t_com.h */


