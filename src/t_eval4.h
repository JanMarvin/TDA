/* t_eval4.h */

#ifndef _TEVAL4_H
#define _TEVAL4_H

/*  functions in t_eval4.c                                                  */

int alloc_mex(int n,int opt);
int mparse(void);          
int alloc_mex_eval(int opt,int n,int m);
int mex_eval(char *s);

extern int MEXMaxLen;   /* length of stack                                  */
extern int *MEXTyp;     /* type of stack entry                              */
extern double *MEXVal;  /* value of operand                                 */
extern int MEXCnt;      /* number of entries in stack                       */
extern double *MEXOP[];
extern int MEXOPRow[];             
extern int MEXOPCol[];  
extern int MEXOPMaxL;  
extern int MEXOPMaxM;  

#endif /* _TEVAL_H */

/*  end of t_eval.h */




