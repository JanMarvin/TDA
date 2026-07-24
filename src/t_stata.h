/* t_stata.h */

#ifndef _TSTATA_H
#define _TSTATA_H

/*  functions in t_stata.c */

int rd_stata(void); 
int wr_stata(void);
int st_gets(char *p);
int st_geti(char *p);
double st_getsf(char *p,int l); 
double st_getf(char *p);
double st_getd(char *p);

extern int HILO;        /* 1 if hilo byte order, 2 if lohi byte order       */

#endif /* _TSTATA_H */

/* end of t_stata.h */


