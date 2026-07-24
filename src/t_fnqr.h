/* t_fnqr.h */

#ifndef _TFNQR_H
#define _TFNQR_H

/*  functions in t_fnqr.c */

int fn_logit1(void);
int fn_logit2(void);
int fn_logit3(void);
int fn_logit4(void);
int fn_probit1(void);
int fn_probit2(void);
int fn_probit3(void);
int fn_prob3(int icase,int wave,int cat,double *par,double *prob);
int fn_probit4(void);
int fn_rmod1(void);

extern int NFLUsed;     /* number of cases used in fn_logit4() */

#endif /* _TFNQR_H */

/*  end of t_fnqr.h */


