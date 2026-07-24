/* t_gio.h */

#ifndef _TGIO_H
#define _TGIO_H
#include "tda_ctx_fwd.h"

/*  functions in t_gio.c */

int gnc(TDAContext *ctx);
int gcon(TDAContext *ctx);
int gdcon(TDAContext *ctx);
int gst(TDAContext *ctx);
int gmst(TDAContext *ctx);
int g_mst(TDAContext *ctx, int gn,int n,int *nodes,int *aci,int *acj,int mflag);
int g_mst1(TDAContext *ctx, int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag);
int gcut(TDAContext *ctx);
int gnst(TDAContext *ctx);
void put_bit(TDAContext *ctx, int i,char *bf,int v);
int get_bit(TDAContext *ctx, int i,char *bf);
void put_bit2(TDAContext *ctx, int i,int j,int n,char *bf,int v);
int get_bit2(TDAContext *ctx, int i,int j,int n,char *bf);
int gcyc(TDAContext *ctx);
int gio(TDAContext *ctx);
int gfcf(TDAContext *ctx);
int gbcf(TDAContext *ctx);
int gep(TDAContext *ctx);
int gflow(TDAContext *ctx);
int gfc(TDAContext *ctx);
void u2_visit(TDAContext *ctx, int k,int gn);
void u2_visit1(TDAContext *ctx, int k,int gn);


#endif /* _TGIO_H */

/* end of t_gio.h */


