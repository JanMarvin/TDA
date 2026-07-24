/* t_gio.h */

#ifndef _TGIO_H
#define _TGIO_H

/*  functions in t_gio.c */

int gnc(void);
int gcon(void);
int gdcon(void);
int gst(void);
int gmst(void);
int g_mst(int gn,int n,int *nodes,int *aci,int *acj,int mflag);
int g_mst1(int gn,int gtyp,int n,int *nodes,int *aci,int *acj,int mflag);
int gcut(void);
int gnst(void);
void put_bit(int i,char *bf,int v);
int get_bit(int i,char *bf);
void put_bit2(int i,int j,int n,char *bf,int v);
int get_bit2(int i,int j,int n,char *bf);
int gcyc(void);
int gio(void);
int gfcf(void);
int gbcf(void);
int gep(void);
int gflow(void);
int gfc(void);
void u2_visit(int k,int gn);
void u2_visit1(int k,int gn);

extern int GSCON_ID;
extern int GSCON_NN;
extern int GSCON_NP;
extern int GSCON_NSP;
extern int GSCON_CNT;
extern int GSCON_LEN;
extern int GSCON_MINLEN;
extern int GSCON_MAXLEN;
extern double GSCON_VAL;
extern double GSCON_MINVAL;
extern double GSCON_MAXVAL;

#endif /* _TGIO_H */

/* end of t_gio.h */


