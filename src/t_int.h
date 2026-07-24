/* t_int.h */

#ifndef _TINT_H
#define _TINT_H

/*  functions in t_int.c */

int niset(void);
void ni_info(void);
int t_int(void);
int prn_int(int m,double rerr);
double ni_gen(double a,double b,int ftyp,int nhp,int htyp,int *err);
int qsniff(double a,double b,double rerr,double *val,int *nf,int *nfu,int ftyp);

extern int NINTMETH;        /* default method of integration (QNG)          */
extern double NINTRERR;     /* default relative error                       */
extern double NINTAERR;     /* default absolute error                       */
extern double NINTTVAR;     /* integration variable t, used in v_eval1()    */
extern int NINTMUsed;       /* method actually used in v_eval1()            */
extern double HIW[][7];     /* parameter for Hermite integration            */
extern double HIZ[][7];

#endif /* _TINT_H */

/* end of t_int.h */


