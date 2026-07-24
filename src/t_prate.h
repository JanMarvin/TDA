/* t_prate.h */

#ifndef _TPRATE_H
#define _TPRATE_H

/*  functions in t_prate.c */

int prate(char *pcmd,int mod);

extern FILE *PRTFd;     /* output file name                                 */
extern int PRTIdx;      /* index for table                                  */

#endif /* _TPRATE_H */

/* end of t_prate.h */


