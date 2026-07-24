/* t_rand.h */

#ifndef _TRAND_H
#define _TRAND_H

/*  functions in t_rand.c */

double random1(void);
double random2(void);
double normal(void);
double normal1(void);
int rdmn_init(int n,double *a);
void rdmn_free(void);
void rdmn(void);
int rdp1(double lambda);

/* ------------------------------------------------------------------------ */
extern int RD1Skip;         /* initially skipped random numbers in random1  */
extern int RD2Skip;         /* initially skipped random numbers in random2  */
extern int RD1Gen;          /* generated random numbers, random1()          */
extern int RD2Gen;          /* generated random numbers, random2()          */
extern int RD1Seed;         /* Seed for random1(), def. RDInit in tda.h     */
extern int RD2Seed;         /* Seed for random2(), def. RDInit in tda.h     */

extern int RD1I;            /* init flag for random1()                      */    
extern int RD2I;            /* init flag for random2()                      */

extern int RDMNInit;        /* set for initialization, by gdat()            */
extern int RDMNN;           /* number of dimensions                         */
extern double *RDMNRd;      /* storage for RDMN random numbers              */

#endif /* _TRAND_H */

/*  end of t_rand.h */

