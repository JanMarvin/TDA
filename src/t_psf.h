/* t_psf.h */

#ifndef _TPSF_H
#define _TPSF_H

/*  functions in t_psf.c */

int psfile(void);
int psetup(void);
void setup_2dsys(void);
int ps_close(int opt,int scx,int kflag);
int pl_axis(int typ,int opt,int opt1,int lt,double lw);
int psetup3(void);
void setup_3dsys(void);
void setup_3proj(double lon,double lat);

#endif /* _TPSF_H */

/*  end of t_psf.h */











