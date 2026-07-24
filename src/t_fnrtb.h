/* t_fnrtb.h */

#ifndef _TFNRTB_H
#define _TFNRTB_H

/*  functions in t_fnrtb.c */

int fn_exp(int resid);
int fn_exp1(int resid);
int fn_exp2(int resid);
int fn_pol(int resid);
int fn_pol1(int resid);
double get_mpol1(double t,int *err); 
double get_mpol2(double t,int *err); 
double get_mpol3(double t,int *err); 
int fn_gm(int resid);
int fn_wei(int resid);
int fn_sic(int resid);
int fn_ll(int resid);
int fn_ll2(int resid);
int fn_ln(int resid);
int fn_ig(int resid);
int fn_gam(int resid);

#endif /* _TFNRTB_H */

/*  end of t_fnrtb.h */




