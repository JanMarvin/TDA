/* t_cmd.h */

#ifndef _TCMD_H
#define _TCMD_H
#include "tda_ctx_fwd.h"

/*  functions in t_cmd.c */

int get_ncmd(TDAContext *ctx, int fn);           
int t_exec(TDAContext *ctx, char *cmd);           
int t_execute(TDAContext *ctx);           
void cmd_err(TDAContext *ctx, int opt);        
int exec_macro(TDAContext *ctx, char *s);   
int check_macro(TDAContext *ctx, char *s,int l);
int new_macro(TDAContext *ctx);
void prn_merr(TDAContext *ctx, char *s,int opt);
int mlist(TDAContext *ctx);       
int mclear(TDAContext *ctx, int opt);      
void clear_tda(TDAContext *ctx);   


#endif /* _TCMD_H */

/* end of t_cmd.h */


