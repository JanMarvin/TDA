/* t_mat.h */

#ifndef _TMAT_H
#define _TMAT_H

/*  functions in t_mat.c */

int t_mat(void);
void mdefcpy(char *d,char *s);
int check_local(char *name);
int m_cmdmsg(void);
int alloc_local(char *cmd);
void free_local(void);
int alloc_mat(void);
void mat_free(void);
int mat_getidx(char *mname,int opt);   
char *get_mname(char *p,char *mname,int opt);
int mat_newmat(char *mname,int m,int n);
int mat_ncheck(char *p,int opt);
int mat_alloc(int idx,int opt);
int mat_newidx(char *mname,int m,int n);
void mat_err(int opt);   
int mat_ncopy(int row,int col,double *x,char *name);
int mat_info(void);     

extern int MaxMat;          /* maximum number of matrices                   */
extern int *MatAlloc;       /* memory allocation for matrices               */
extern double **MatVal;     /* matrix data                                  */
extern char **MatName;      /* names of matrices                            */
extern char **MatDef;       /* definition of matrices                       */
extern int *MatRow;         /* number of rows                               */
extern int *MatCol;         /* number of columns                            */

extern int PMATFmt1;        /* print format: mfmt=                          */
extern int PMATFmt2;
extern int PMATFmtF;
extern char PMATFmtS[];        
extern char *MatCmdBuf;     /* saves command for error messages             */

#endif /* _TMAT_H */

/*  end of t_mat.h */















