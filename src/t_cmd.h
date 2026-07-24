/* t_cmd.h */

#ifndef _TCMD_H
#define _TCMD_H

/*  functions in t_cmd.c */

int get_ncmd(int fn);           
int t_exec(char *cmd);           
int t_execute(void);           
void cmd_err(int opt);        
int exec_macro(char *s);   
int check_macro(char *s,int l);
int new_macro(void);
void prn_merr(char *s,int opt);
int mlist(void);       
int mclear(int opt);      
void clear_tda(void);   

extern int IERRFlg;            /* set by ierr command */

#endif /* _TCMD_H */

/* end of t_cmd.h */


