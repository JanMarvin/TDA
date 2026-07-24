/* t_mdat.h */

#ifndef _TMDAT_H
#define _TMDAT_H

/*  functions in t_mdat.c */

int clear(void);
void clear_a(void);
int clearnl(void);
int tsel(void);       
void tsel_off(int opt);       
int cwt(void);       
int wr_sys(void);
int rd_sys(void);
int dblock(void);          
int dblock_alloc(int n);   
int repsel(void);          
int repsel_off(void);   
int repsel_alloc(int n);   

extern int *DBlckPtr;       /* maps cases to blocks                         */
extern int BNOC;            /* number of blocks                             */ 
extern int *REPSelect;      /* indices for repsel                           */
extern int REPSelFlg;       /* set if repsel active                         */

#endif /* _TMDAT_H */

/*  end of t_mdat.h */











