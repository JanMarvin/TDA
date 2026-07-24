/* t_con.h */

#ifndef _TCON_H
#define _TCON_H

/*  functions in t_con.c */

int p_con(char *pcmd,int nx,int nif,int mw,int nw,double *w,int *ne,int *ni);
int con_proc(char *pcmd);
void con_free(void);

extern int NCon;            /* Number of constraints                        */
extern int NCon1;           /* Number of constraints used for minimization  */
extern double *ConQ;        /* Projection matrix (NParm, NParm1)            */
extern int *CIP;            /* Pointers to columns of ConQ                  */

#endif /* _TCON_H */

/* end of t_con.h */


