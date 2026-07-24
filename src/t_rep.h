/* t_rep.h */

#ifndef _TREP_H
#define _TREP_H

/*  functions in t_rep.c */

int t_repeat(int typ,int lev);
int t_endrepeat(int typ,int lev);
void rep_free(void);
int rep_scmd(char *cmd);
int check_break(char *cmd);   
int check_if(char *cmd);   
void prn_nex(char *cmd);   

extern int REPLev;          /* current level of repeat                      */
extern int BREAKFlg;        /* set if break command                         */

#endif /* _TREP_H */

/* end of t_rep.h */


