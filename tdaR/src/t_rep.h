/* t_rep.h */

#ifndef _TREP_H
#define _TREP_H
#include "tda_ctx_fwd.h"

/*  functions in t_rep.c */

int t_repeat(TDAContext *ctx, int typ,int lev);
int t_endrepeat(TDAContext *ctx, int typ,int lev);
void rep_free(TDAContext *ctx);
int rep_scmd(TDAContext *ctx, char *cmd);
int check_break(TDAContext *ctx, char *cmd);   
int check_if(TDAContext *ctx, char *cmd);   
void prn_nex(TDAContext *ctx, char *cmd);   


#endif /* _TREP_H */

/* end of t_rep.h */


