/* t_gdat.h */

#ifndef _TGDAT_H
#define _TGDAT_H

/*  functions in t_gdat.c */

int new_var(void);
double get_data(int j,int i);
void put_data(double x,int j,int i);
void put_str(char *buf,int blen,int j,int i,int bflag);
void get_str(char *buf,int j,int i);
double dscan(char *p,int len,int *mval);
char *skip_sep(char *p);
char *skip_dval(char *p);
void make_vfmt(int nidx);
int sdnvar(void);
int check_nvdef(void);
void sdnvar_close(void);
int get_nsdxy(char *buf,double *x,double *y);
int sdnvar_alloc(int opt,int n);
int check_sd(int opt);

/*--------------------------------------------------------------------------*/
extern int SDVarDef;        /* set to 1 if spatial data defined             */
extern int SDVarFDef;       /* set to 1 if SDVarFd open                     */
extern FILE *SDVarFd;       /* file handle                                  */
extern int SDVarSDID;       /* index of SDID variable                       */
extern int SDVarSDPtr;      /* index of SDPtr variable                      */
extern int SDVarSDTyp;      /* index of SDTyp variable                      */
extern int SDVarSDN;        /* index of SDN variable                        */
extern int SDVarMax;        /* max number of points in object               */
extern int SDVarNT;         /* total number of points                       */
extern double SDVarXMin;    /* bounding box                                 */
extern double SDVarXMax;
extern double SDVarYMin;
extern double SDVarYMax;

extern double *SDVarX;      /* standard array for spatial objects           */
extern double *SDVarY;
extern int SDVarN;          /* number of allocated entries                  */

/*--------------------------------------------------------------------------*/
extern int DMDef;       /* set to 1 if data matrix is available.            */
extern int NOCMaxA;     /* max number of cases.                             */
extern int NOCDM;       /* actual number of cases in data matrix.           */
extern int NOC;         /* number of data matrix rows to be used, this      */
                        /* depends on the actual tsel command.              */

extern int WIVar;       /* internal variable number of weights              */
extern double WSum;     /* sum of weights                                   */
extern double WSumS;    /* factor defined by wnorm option                   */
extern double WNorm;    /* sum W(i) * WNOrm = (required) WSum               */
extern int WNormFlag;   /* 1 if wnorm used in cwt command                   */

extern int VCJMax;      /* Highest Cj index                                 */
extern double *VCJVal;  /* array with values of required Cj                 */

extern double MBlnkVal; /* missing value code: blank                        */
extern char BMsk[];     /* Bitmasks for bit-wise stored variables           */
extern double *AVVAL;   /* temporary storage of archive variables           */

extern int TSelFlg;         /* set if tsel command active                   */
extern int *TSelect;        /* array with select indices                    */
extern int TSelectA;        /* if allocated                                 */

extern int GDFlg;           /* set by df parameter.                         */
extern int GDNRec;          /* number of records written to GDFd            */

#endif /* _TGDAT_H */

/*  end of t_gdat.h */











