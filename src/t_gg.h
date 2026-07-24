/* t_gg.h */

#ifndef _TGG_H
#define _TGG_H

/*  functions in t_gg.c */

int giset(void);
int gcni(void);
int gcset(void);
double g_mineval(int n1,int *nodes1,int n2,int *nodes2,int opt);
int gcliq(void); 
int ggcliq(void);
int g_gclptr_i(int max);
int g_gclptr_c(int n,int *nodes);

extern int **GCL_PTR;
extern int *GCL_NPTR;
extern int GCL_PN;

#endif /* _TGG_H */

/* end of t_gg.h */


