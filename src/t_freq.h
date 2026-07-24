/* t_freq.h */

#ifndef _TFREQ_H
#define _TFREQ_H

/*  functions in t_freq.c */

int cfreq(int n, int m, int *x, int mcel, int *value, int *freq, int *bf,
    int wflg, double *wt, double *wfreq);
int cfreq1(int n,short *idx,int mcel,int *value,int *freq,int *bf,
    int widx,double *wfreq);

#endif /* _TFREQ_H */

/*  end of t_freq.h */

