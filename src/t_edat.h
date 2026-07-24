/* t_edat.h */

#ifndef _TEDAT_H
#define _TEDAT_H

/*  functions in t_edat.c */

int check_edj(int j);
int edef(void);
int def_edat(int mode);
void edat_off(int opt);
void prn_edef(void);
void edef_info(void);
int check_edat(void);
void prn_edat(void);
int get_edat(int i,int *sn,int *org,int *des,double *ts,double *tf);
int sort_ed(int opt); 
int get_spell(int init,int *icase,int *sn,int *org,int *des,
    double *ts,double *tf,int *spl,int *nspl);
int epdat(void);
int epsdat(void);

/*--------------------------------------------------------------------------*/
extern int EDAvail;     /* 1 if duration data allocated and available       */
extern int MEFlg;       /* set for multiepisode data                        */

extern int ID_Idx;      /* index in CmdDef of id command                    */
extern int SN_Idx;      /* index in CmdDef of sn command                    */
extern int TS_Idx;      /* index in CmdDef of ts command                    */
extern int TF_Idx;      /* index in CmdDef of tf command                    */
extern int ORG_Idx;     /* index in CmdDef of org command                   */
extern int DES_Idx;     /* index in CmdDef of des command                   */

extern int MaxTran;     /* max number of (sn,org,des) combinations          */

extern int MaxSnn;      /* highest spell number                             */
extern int MaxOrg;      /* highest origin state number                      */
extern int MaxOrg1;     /* MaxOrg1 = MaxOrg + 1                             */
extern int MaxDes;      /* highest destination state number                 */
extern int MaxDes1;     /* MaxDes1 = MaxDes + 1                             */
extern int NTran;       /* number of (sn,org,des) combinations              */
extern int *SnTran;     /* this is an array of lengt NTran for all          */
extern int *OrgTran;    /* combinations of (sn,org,des) in the input data,  */
extern int *DesTran;    /* sorted according to sn,org,des.                  */
extern int NTran1;      /* number of (sn,org,des) transitions               */
extern int *SnTran1;    /* this is an array of lengt NTran1 for all         */
extern int *OrgTran1;   /* transitions (sn,org,des) in the input data,      */
extern int *DesTran1;   /* sorted according to sn,org,des.                  */

extern int *TranPtr;    /* pointer to entries in SnTran1,OrgTran1,DesTran1  */

extern int *TranNE;     /* number of episodes                               */
extern float *TranWE;   /* weighted number of episodes                      */
extern float *TranMD;   /* weighted mean duration                           */
extern float *TranTSM;  /* minimum of starting times                        */
extern float *TranTFM;  /* maximum of ending times                          */

extern float *TWFreq;   /* weighted number of uncensored episodes           */
extern float *TWDur;    /* summed durations                                 */

extern int ESortFlg;    /* 1 if sorted                                      */
extern int *TSIdx;      /* pointer for sorting according to ts              */
extern int *TFIdx;      /* pointer for sorting according to tf              */
/*--------------------------------------------------------------------------*/
extern int NSP;             /* number of split variables                    */
extern short SPLITVar[];    /* indices of split variables                   */
extern int NSplits;         /* number of splits, set in get_spell()         */
extern int NSPLEvalErr;     /* number of errors when evaluating type 5 vars */

#endif /* _TEDAT_H */

/* end of t_edat.h */


