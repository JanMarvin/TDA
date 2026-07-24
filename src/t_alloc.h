/* t_alloc.h */

#ifndef _TALLOC_H
#define _TALLOC_H

/*  functions in t_alloc.c */

void a_clean(void);
int alloc_actmp(int n);
int alloc_actmp1(int n);
int alloc_acx(int n);
int alloc_acy(int n);
int alloc_acy1(int n);
int alloc_acz(int n);
int alloc_act(int n);
int alloc_acu(int n);
int alloc_acv(int n);
int alloc_acw(int n);
int alloc_acxf(int n);
int alloc_acyf(int n);
int alloc_aczf(int n);
int alloc_acuf(int n);
int alloc_acvf(int n);
int alloc_acm(int n);
int alloc_acn(int n);
int alloc_aci(int n);
int alloc_acj(int n);
int alloc_acr(int n);
int alloc_acs(int n);
int alloc_ack(int n);
int alloc_acl(int n);
int alloc_acc(int n);
int alloc_acd(int n);
int alloc_ace(int n);
int alloc_acf(int n);
int alloc_acns(int n);
int alloc_acms(int n);
int alloc_acptr(int n);
int alloc_aciptr(int n);
int alloc_acjptr(int n);

extern int AcTmpN;
extern double *AcTmp;
extern int AcTmp1N;
extern double *AcTmp1;
extern int AcXN;
extern double *AcX;
extern int AcYN;
extern double *AcY;
extern int AcY1N;
extern double *AcY1;
extern int AcZN;
extern double *AcZ;
extern int AcTN;
extern double *AcT;
extern int AcUN;
extern double *AcU;
extern int AcVN;
extern double *AcV;
extern int AcWN;
extern double *AcW;
extern int AcXFN;
extern float *AcXF;
extern int AcYFN;
extern float *AcYF;
extern int AcZFN;
extern float *AcZF;
extern int AcUFN;
extern float *AcUF;
extern int AcVFN;
extern float *AcVF;
extern int AcNM;
extern int *AcM;
extern int AcNN;
extern int *AcN;
extern int AcIN;
extern int *AcI;
extern int AcJN;
extern int *AcJ;
extern int AcRN;
extern int *AcR;
extern int AcSN;
extern int *AcS;
extern int AcKN;
extern int *AcK;
extern int AcLN;
extern int *AcL;
extern int AcCN;
extern char *AcC;
extern int AcDN;
extern char *AcD;
extern int AcEN;
extern char *AcE;
extern int AcFN;
extern char *AcF;
extern int AcNNS;
extern short *AcNS;
extern int AcNMS;
extern short *AcMS;
extern int AcNPtr;
extern char **AcPtr;
extern int AcNIPtr;
extern int **AcIPtr;
extern int AcNJPtr;
extern int **AcJPtr;

#endif /* _TALLOC_H */

/* end of t_alloc.h */


