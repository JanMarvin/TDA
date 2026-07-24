/* t_gen.h */

#ifndef _TGEN_H
#define _TGEN_H

/*  functions in t_gen.c */

void gerr_exit(int n);
void memrq(int n, int m);
void prn_mem(void);      
void newline(void);           
char *check_comment(char *p);
int check_drec(char *buf);
void prn_message(int rec,int opt,int wflag);
int get_sline(char *p); 
void prnchar(unsigned char c, int n, int mode);
void fprnchar(FILE *fd, unsigned char c, int n, int mode);
int mpr(int m,int n,double *a,int ivflg,double *b,char *fmt,char *fn,char *s,int aflag);
int mprf(int m,int n,float *a,char *fmt,FILE *fd);
int mprd(int m,int n,double *a,char *fmt,FILE *fd);
char *skip_b(char *p);
char *skip_c(char *p);
char *skip_cb(char *p);
char *skip_int(char *p);
char *skip_dbl(char *p);
char *skip_expr(char *p);
char *skip_blev(char *p);
char *skip_nc(char *p);
char *skip_com(char *p);
char *skip_xa(char *p);
int get_fname(char *p,char *q);
int get_float(char *p,int nmax,float *x);
void printf1(const char *fmt, ...);
void printf2(const char *fmt, ...);
void printfe(const char *fmt, ...); 
void fflushe(void);
#if TIME_ON
void prn_time(FILE *fd);      
void prn1_time(char *buf); 
#endif
  
/*--------------------------------------------------------------------------*/
extern int SILENTFlg;       /* controls output printing                     */

extern int MemReq;          /* Number of bytes of requested memory          */
extern int MxMReq;          /* Max memory request (bytes)                   */

extern int CmdBufL;         /* actual command buffer length                 */
extern char *CmdBuf;        /* command buffer                               */

extern int INTMAX;          /* largest integer                              */
extern double EPSI;         /* machine's epsilon                            */
extern double EPSI1;        /* sqrt(EPSI)                                   */
extern double EPSI2;        /* 1000 * EPSI                                  */
extern double DBLMAX;       /* largest number in machine                    */
extern double DBLMIN;       /* smallest number in machine                   */
extern float MINFLOAT;
/*--------------------------------------------------------------------------*/

#endif










                                                

