/* t_rzoo.h */

#ifndef _TRZOO_H
#define _TRZOO_H

/*  functions in t_rzoo.c */

void arcd_off(void);
int arcd(void);
void prn_afiles(void);
int get_drec(int fn, char *buf, int nmax);
int get_record(int fn, char *buf, int nmax);
int alloc_avar(int idx,int n,int opt);
int check_avar(int n,int adic);
int get_avar(void);
void get_astr(char *buf,int i);
int arcc(void);
int arcv(void);
int arcvc(void);
int get_afmt(int *w1,int *w2);

/* ------------------------------------------------------------------------ */
extern int ARCDef;          /* set if archive sucessfully opened            */
extern char ZADNam[];       /* name of archive description file             */
extern char ZOONam[];       /* name of zoo data archive                     */
extern FILE *ZOOFd;         /* file handle for zoo archive                  */
extern int ZANF;            /* number of files in description file          */
extern char **ZAFNam;       /* names of files                               */
extern int *ZAFNRec;        /* number of records                            */
extern short *ZAFNum;       /* logical file number                          */
extern short *ZAFTyp;       /* type of file                                 */
extern short *ZAFRLen;      /* record length                                */
extern short *ZAFNVar;      /* number of variables                          */
extern int ZABLen;          /* max record length of files                   */
extern char *ZABuf;         /* Record buffer  for files                     */
extern int FDefLen;         /* max length of file names                     */
extern char *ZAFReq;        /* flags for required files                     */
extern int VFN;             /* Number of variable description file          */
/* ------------------------------------------------------------------------ */
extern int *ZAZOfs;         /* Direntry offset in archive                   */
extern int *ZAZSiz;         /* Original size of file in bytes               */
extern char *ZAZTyp;        /* Packing method (1 - 2)                       */
extern char ZAInit;         /* Set if decompression initialized             */
extern int ZAAlloc;         /* Set if memory allocated                      */
extern char ZAEOF;          /* Set if EOF reached                           */
extern long FPtr;           /* File pointer for the files                   */
/* ------------------------------------------------------------------------ */
/*  For each file we allocate a read buffer where the decompression is      */
/*  done, and a line buffer.                                                */

extern char *Out_Buf_Adr;   /* Output buffer used for decompression         */
extern char *In_Buf_Adr;    /* Input buffer used for decompression          */
extern char *Out_Buf_Ptr;   /* Pointer to output buffer                     */
extern int Out_Buf_Cnt;     /* Character count of output buffer             */
/* ------------------------------------------------------------------------ */


#endif /* _TRZOO_H */

/*  end of t_rzoo.h */












