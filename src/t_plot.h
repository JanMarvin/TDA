/* t_plot.h */

#ifndef _TPLOT_H
#define _TPLOT_H

/*  functions in t_plot.c */

int check_pcmd(int opt,int dim);              
int check_ps(int dim);              
void upd_bbox(int opt,double x,double y); 
void set_clip(void); 
void ps_ltyp(int typ);
void ps_lwidth(double lw);
void ps_fill(double g); 
void ps_2dplot(double x,double y, int opt);
double ps_2dx(double x);
double ps_2dy(double y);
void plot_str(double x,double y,char *s, double siz, int opt,int adj,
    int r,int xopt,int cflag);
void ps_nlab(int n,double x,double y,double fs,int opt);
void ps_sym(int typ, double px, double py,double siz);
int ps_grid(int typ);
int pl_label(int typ);
int pl_text(void);
int pl_frame(void);
int pl_rec(void);
int pl_plotp(int typ);
int pl_plotm(void);
int pl_plotf(void);
int pl_ploth(void);
int pl_plotd(void);
int pl_ploto(int typ);
int pl_plotk(void); 
int pl_plotch(void);
int pl_plots(int typ);
int pl_plg(void);


/*--------------------------------------------------------------------------*/
extern int PSFFlg;          /* set if PostScript output file defined        */
extern int PS3DFlg;         /* set by psetup3().                            */
extern int PSPROJ;          /* set by psetupg()                             */

extern char PSFName[];      /* name of PostScript output file               */
extern FILE *PSFd;          /* file handle for PostScript output file       */
extern int PSHdFlg;         /* set if header is written                     */
extern int PSIFlg;          /* set if basic initialization is done          */
extern int BBLX,BBLY,BBUX,BBUY;    /* bounding box                          */
extern int PSONUM;          /* plot object number                           */ 

extern double PtMM;         /* points per mm                                */

extern double LWCSMax;      /* max line width used for axes                 */
extern double LWDef;        /* default line width                           */
extern double LW1Def;       /* default second line width (grids etc.)       */
extern double FSDef;        /* default font size                            */
extern double PTLen;        /* tick length in mm                            */

extern int PSXOrg;          /* origin of PostScript figure                  */
extern int PSYOrg;
extern int XOrg;            /* origin of PostScript figure, redefined       */
extern int YOrg;            /* by psorg commands.                           */
extern double SCALXFac;     /* PostScript x scaling factor                  */
extern double SCALYFac;     /* PostScript y scaling factor                  */
extern double ROTFac;       /* PostScript rotation factor                   */
                         
extern double PXLen;        /* length of X axis in mm                       */
extern double PYLen;        /* length of Y axis in mm                       */
extern double PSXLen;       /* lenght of X axis in points                   */
extern double PSYLen;       /* lenght of Y axis in points                   */

extern int PSLog[];         /* flags for logarithmic axes                   */
extern double PA1[];        /* logical start of axis                        */
extern double PA2[];        /* logical end of axis                          */

extern double UXLen;        /* length of x axis in user units               */
extern double UYLen;        /* length of y axis in user units               */
extern double UZLen;        /* length of z axis in user units               */

extern double PX3D[];       /* logical x axis                               */
extern double PY3D[];       /* logical y axis                               */
extern double PZ3D[];       /* logical z axis                               */
extern double PX3C;         /* center of x axis                             */
extern double PY3C;         /* center of y axis                             */
extern double PZ3C;         /* center of z axis                             */

extern double PSLon;        /* direction of projection, longitude           */
extern double PSLat;        /* direction of projection, latitude            */

extern double PSPR11,PSPR12,PSPR13;
extern double PSPR21,PSPR22,PSPR23;
extern double PSPR31,PSPR32,PSPR33;

extern double PSPI11,PSPI12,PSPI13;
extern double PSPI21,PSPI22,PSPI23;
extern double PSPI31,PSPI32,PSPI33;

#endif /* _TPLOT_H */

/* end of t_plot.h */









