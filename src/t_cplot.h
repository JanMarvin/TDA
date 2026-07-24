/* t_cplot.h */

#ifndef _TCPLOT_H
#define _TCPLOT_H

/*  functions in t_cplot.c */

int pl_plotc(void);
void cont_plot(int nlev,double *flev,int opt);
int pl_plotcm(void);
int pl_plotr(void);

extern int PCNX;        /* number of x axis intervals                       */
extern int PCNY;        /* number of y axis intervals                       */
extern double PCDX;     /* width of x axis intervals                        */
extern double PCDY;     /* width of y axis intervals                        */
extern float *PCFV;     /* array with function values                       */
extern char *PCBitM;    /* array used in cont_plot()                        */

#endif /* _TCPLOT_H */

/* end of t_cplot.h */









