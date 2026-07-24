/* t_var.h */

#ifndef _TVAR_H
#define _TVAR_H

/*  functions in t_var.c */

char get_vnchar(char c);
int check_vname(char *p);
int get_vnlen(char *p);
int get_vnl(int n,short *vinum);
int alloc_vmax(void);
int get_nidx(void);
void clear_vidx(int i);
int get_vidx(char *p,char *name);
int get_vidx1(char *p,char *name);
int get_nvidx(char *p,short *vinum);
int save_var(char *vd,int opt);
void var_err(char *vd);
int alloc_vdat(int idx,int opt);
int alloc_vsdat(int n,short *ivar);
int get_slen(int j,int noc);
void prn_var(int idx);
void prn_vname(int i);
void prn_vlabel(int i);
void prn_hvar(void);
void prn_hlabel(void);
void clear_dm(void);
void clear_var(int idx);
void clear_avar(int idx);
void free_var(int i);
int nlist(void);
int nl_check(char *nl);
void nl_free(int j);
int alloc_vl(int n);
char *get_nvi(char *s,int *n,int nv,short *vidx,int *nb);
char *get_nvia(char *s,int *n,int opt,int *nb);
int recode(void);
int ndvar(void);
int get_mxvlen(int n,short *vidx);
void prn_vlist(int n,short *vidx);
int svb_alloc(int n);

extern int MaxNV;       /* max number of variables (tda.cfg)                */
extern int NVAR;        /* number of variables                              */

extern short VIFirst;   /* index of first variable                          */
extern short VILast;    /* index of last variable                           */

extern char *VAlloc;    /* flags for allocated data structures              */
extern short *VNxt;     /* pointer to next variables                        */

extern char **VName;    /* names of variables                               */
extern char **VDef;     /* definitions of variables                         */
extern char **VLabel;   /* labels of variables                              */

extern char **VDPtr;    /* pointer to data fields                           */

extern short *VSLen;    /* storage length                                   */
extern char *VTyp;      /* type of variable                                 */
extern char *VTypA;     /* set if archive variable                          */
extern short *VPFmt1;   /* print format of variables                        */
extern short *VPFmt2;
extern char **VPFmtS;

extern int **VESTyp;    /* parser stack for variables                       */
extern double **VESVal; /* parser value for variables                       */
extern short *VESCnt;   /* size of parser arrays                            */

extern int VNameLen;    /* max length of variable names                     */
extern int VLabelLen;   /* max length of variables labels, also used to     */
                        /* flag if at least one var label present.          */

extern int NVArc;       /* number of archive variables                      */
extern int NVArc1;      /* number of new archive variables, set in save_var */
extern short *AVIdx;    /* internal variable number                         */
extern short *AVOff;    /* offset of variable                               */
extern short *AVLen;    /* length of variable                               */
extern short *AVFmt1;   /* format of variable                               */
extern short *AVFmt2;       
extern int *AVMBlnk;    /* missing values: blank                            */
extern int *AVMStar;    /* missing values: star                             */
extern int *AVMPnt ;    /* missing values: point                            */
extern int *AVMGen ;    /* missing values: general                          */
extern int AVDFN;       /* index to data file for archive variables         */
extern int AVIdxA;      /* flags for mem allocation                         */
extern int AVOffA;
extern int AVLenA;
extern int AVFmt1A;
extern int AVFmt2A;
extern int AVMBlnkA;
extern int AVMStarA;
extern int AVMPntA;
extern int AVMGenA;

extern short *VStrN;    /* first column of string variable                  */

extern int NNL;         /* number of namelists                              */
extern char NLName[][VNLMax + 1];   /* names                                */
extern short NLNV[];           /* number of variables                       */
extern short *NLVIdx[];        /* indices of variables                      */
extern short VLTyp[];          /* counting variable types                   */

extern int VLNV;        /* number of variables in VLVIdx[]                  */
extern short *VLVIdx;

extern char *SVBuf;     /* buffer for longest string                        */
extern int SVBufA;      /* allocated                                        */
extern int SVBufLen;    /* current maximal string length                    */

#endif /* _TVAR_H */

/* end of t_var.h */


