/* t_seq.h */

#ifndef _TSEQ_H
#define _TSEQ_H

/*  functions in t_seq.c */

int t_seq(void);
int check_sdj(int j);
int seq_sget(int i,int t,int n);
int seq_getsn(int sn,int opt);
int seq_scheck(int i,int j);
void seq_afree(int opt);                   
void prn_ssel(int n);
int seq_csel(void);
void seq_dtda(int typ,int sn,char *fname,int noc,int nv,short *vidx);

extern int SeqDN;           /* number of sequences                          */  
extern int SeqDNH;          /* highest sequence number                      */
extern char SeqDT[];        /* type of sequence                             */
extern short SeqDNV[];      /* number of variables                          */
extern short *SeqDV[];      /* indices of variables                         */

extern short SeqSTN[];      /* number of different states                   */
extern short *SeqSTNI[];    /* array with SeqSTN[] state numbers            */
extern short *SeqSTNII[];   /* array with SeqSTNH[] inverse state numbers   */
extern short SeqSTNH[];     /* highest state number                         */
extern int SeqSTNHMax;      /* highest state number in all sequences        */
extern int SeqTMin[];       /* minimum of time points                       */
extern int SeqTMax[];       /* maximum of time points                       */
extern int SeqTCMin;        /* minimum of time points for all sequences     */
extern int SeqTCMax;        /* maximum of time points for all sequences     */

#endif /* _TSEQ_H */

/* end of t_seq.h */


