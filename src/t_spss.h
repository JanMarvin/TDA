/* t_spss.h */

#ifndef _TSPSS_H
#define _TSPSS_H

/*  functions in t_spss.c */

int rd_spss(void); 
int wr_spss(void);
int rd_spss1(void); 
int wr_spss1(void);
int get_dvarp(int *np,char *vname);
void make_arcd(int fn,char *fname,int len,int nrec,int vn);

#endif /* _TSPSS_H */

/* end of t_spss.h */


