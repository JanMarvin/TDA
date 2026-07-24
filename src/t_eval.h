/* t_eval.h */

#ifndef _TEVAL_H
#define _TEVAL_H

/*  functions in t_eval.c                                                  */

int alloc_est(int n,int opt);
void parse(short cnt,int *typ,double *val,int *idx);
int v_parse(char *s,int iflag);
int v_search(char *p,int *len,int *narg); 
int v_search1(char *p,int *len);
char *skip_const(char *p);                   
void prn_emsg1(int n);
void prn_emsg2(int n);

/* ------------------------------------------------------------------------ */
extern int IVEXPRFlg;   /* set if expression contains intervals             */
extern int ESErr;       /* global error flag                                */
extern int *ESTyp;      /* arrays for result of parser                      */
extern int *ESIdx;      /* number of matrix, variable, or namelist          */
extern char *ESFlg;
extern double *ESVal;
extern short ESCnt;
extern int ESCase;  
extern int ESCase;      /* case number for temporary use                    */
extern int ESIVP;       /* pointer to val[]                                 */
extern int ESDeriv;     /* type of derivative                               */
extern int ESDerivI;    /* index of argument                                */

extern int EDVALSn;
extern int EDVALOrg;
extern int EDVALDes;
extern double EDVALTs;
extern double EDVALTf;
extern double EDVALTime; 

extern int IAOFFS;      /* Offset intermediate arguments                    */  
extern int FAOFFS;      /* Offset function arguments ... FAOFFS + MaxP      */  
extern int NLOFFS;      /* Offset namelists                                 */
extern int MOFFS;       /* Offset matrices                                  */
extern int VOFFS;       /* Offset V variables                               */
extern int COFFS;       /* Offset C variables                               */
extern int COFFMAX;     /* Max offset C variables                           */
extern int IFTYP;       /* special types for if then else operators         */

extern int IFLevel;     /* level of if then else                            */
extern int E2FNum;      /* number for first record in block, set in eval2() */
extern int E2LNum;      /* number for last record in block, set in eval2()  */

extern int CIdxMax;     /* max index of ci terms                            */
extern int XIdxMax;     /* max index of xi terms                            */

extern int OPT2S;       /* begin of spatial operators                       */
extern int OPT2A;       /* begin of type 2 operators                        */
extern int OPT2C;
extern int OPT2B;       /* end of type 2 operators                          */

extern char *PREVName;  /* pointer to vname for pre(vname)                  */
extern int PREVNum;     /* vnum of vname. set by save_var().                */

extern int FNFlg;       /* set if valid function definition                 */
extern int FCFlg;       /* set if valid constraint function definition      */
extern int FNPN;        /* number of parameters                             */
extern int FNMLen;      /* max length of parameter strings                  */
extern int FVFlg;       /* set if function refers to data matrix variables  */

extern char *FNPDef[];  /* array with parameter strings                     */
extern short FNPLen[];  /* array with parameter length                      */
extern double FNPVal[]; /* array with parameter values                      */
extern double FNPVal1[];           

extern double *FNPGrad[];   /* gradients                                    */
extern double *FNPGrad1[];  /* gradients, upper bounds                      */
extern double *FNPHess[];   /* hessians                                     */
extern int FNPGradA;        /* allocated                                    */
extern int FNPGrad1A;       /* allocated                                    */
extern int FNPHessA;        /* allocated                                    */

extern short FNPEATyp[];    /* flag for dependency on arguments             */
extern short FNPECnt[];     /* array with parser stacks                     */
extern int *FNPETyp[];
extern double *FNPEVal[];
extern short FNPELev[];     /* level of intermediate expression             */

extern short FNECnt;        /* parser stack for function expression         */
extern int *FNETyp;
extern double *FNEVal;

extern short FNECCnt;       /* parser stack for constraint function         */
extern int *FNECTyp;
extern double *FNECVal;

extern int FNArgN;          /* number of parameters                         */
extern char *FNArgDef[];    /* array with parameter strings                 */
extern short FNArgLen[];    /* array with parameter length                  */
extern short *FNArgSP;      /* pointer for sorting arguments                */
extern int FNArgSPA;        /* if allocated                                 */
extern short *FNArgSPI;     /* inverse pointer                              */
extern int FNArgSPIA;       /* if allocated                                 */

extern int FNGradAA;
extern int FNHessAA;
extern short *FNATyp;       /* type of argument                             */
extern int FNATypA;
extern int FNEVATyp;        /* global for type of arguments                 */  

extern int FNNLev;          /* number of level in function definition       */
extern int *FNLVar;         /* index of variable used to define level       */
extern int FNLVarA;

extern double *FNArgVal;
extern int FNArgValA;
extern double *FNArgVal1;
extern int FNArgVal1A;

extern int MaxSTD;      /* max size of stacks                               */
extern int FNGradN;     /* stack counter gradient                           */
extern int FNGradA;     /* max size of gradient stack                       */
extern int FNGradL;     /* length of gradient vector                        */
extern int FNHessN;     /* stack counter hessian                            */
extern int FNHessA;     /* max size of hessian stack                        */
extern int FNHessL;     /* length of hessian                                */
extern double **FNGrad; /* stack for gradient                               */
extern double **FNHess; /* stack for hessian                                */
extern int FNEVATyp;    /* global for type of arguments                     */  

extern int MatExprFlg;  /* set while parsing matrix expressions             */
extern int MatExprRow;
extern int MatExprCol;

#endif /* _TEVAL_H */

/*  end of t_eval.h */




