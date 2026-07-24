/****************************************************************************/
/*  t_xwin                                                                  */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-97 Goetz Rohwer. All rights reserved.           */
/*                                                                          */
/*  This file is part of TDA.                                               */
/*                                                                          */
/*  TDA is free software; you can redistribute it and/or modify             */
/*  it under the terms of the GNU General Public License as published by    */
/*  the Free Software Foundation; either version 2 of the License, or       */
/*  (at your option) any later version.                                     */
/*                                                                          */
/*  TDA is distributed in the hope that it will be useful,                  */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/*  GNU General Public License for more details.                            */
/*                                                                          */
/*  You should have received a copy of the GNU General Public License       */
/*  along with this program (it should be in a file named COPYING),         */
/*  if not, write to the Free Software Foundation, Inc.,                    */
/*  675 Mass Ave, Cambridge, MA 02139, USA.                                 */
/*                                                                          */

#include "tda.h"

#if S_XWIN          /* only included if S_XWIN is set in tda.h */

#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_plot.h"
#include "t_ml.h"
#include "t_gf.h"

#if S_UNIX
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <X11/Xatom.h>
#include <X11/cursorfont.h>
#include <X11/keysym.h>
#endif
#if S_DOS
#include "t_xnt.h"
#endif

/* ------------------------------------------------------------------------ */
/*  functions in t_xwin                                                     */

int x_show(void);
int x_graph(void);
int x_getx(double x,int mode);
int x_gety(double y,int mode);
void x_draw(int n,int mode);
void x_fill(int n,int mode,double g);
void x_arc(double r,double a1,double a2);
void x_text(char *s);
void x_symbol(char *s);
void x_setclip(int n);
void x_offclip(void);
void x_dfonts(char *s);

/*--------------------------------------------------------------------------*/
/*  Global variables and defines                                            */

#define MAX_XPOINTS 1000        /* set in xlib.c, this is the maximum       */

Display      *display;
Window       win,root;
int          screen;
unsigned int xwidth = 400, xheight = 200;
long int     black, white;
XEvent       event;
GC           gc;
XGCValues    gc_values;

#define XFontDef "7x13"     /* default font                                 */
XFontStruct  *font_info;    
char *XFontD = NULL;        /* actual font                                  */
XCharStruct min_b,max_b;

Visual *default_visual;
Colormap default_cmap;
XColor exact_def;

#define MAX_COLORS 8
int x_ncolors = 0;
int x_colors[MAX_COLORS];
char *ColName[] = {
    "red",        /*  255,   0,   0,   */
    "blue",       /*    0,   0, 255,   */
    "green",      /*    0, 255,   0,   */
/*  "orange",  */ /*  255, 165,   0,   */
    "magenta",    /*  255,   0, 255,   */
    "brown",      /*  165,  42,  42,   */
/*  "yellow",*/   /*  255, 255,   0,   */
    "cyan",       /*    0, 255, 255,   */
    "plum",       /*  221, 160, 221,   */ 
/*  "gold",  */   /*  255, 215,   0,   */ 
/*  "pink",  */   /*  255, 192, 203,   */
    "purple"      /*  160,  32, 240,   */ 
};
                                                                                   
int XSHOWFIRST = 1;                                                          

XPoint *xpoints;    
int XPAllocN = 0;

XRectangle xrec[1];   
int XBLX,XBLY,XBUX,XBUY;
double XPPnt = 0.0;         /* X axis pixels per PostScript points          */
double YPPnt = 0.0;         /* Y axis pixels per PostScript points          */
double XPPMM = 0.0;         /* X axis pixels per mm                         */
double YPPMM = 0.0;         /* Y axis pixels per mm                         */
double XLEN = 0.0;          /* length of plot in mm                         */
double YLEN = 0.0;          /* height of plot in mm                         */
int XNXOrg = 0;
int XNYOrg = 0;

int XLWID = 1;              /* line width                                   */
int XLTYP = 1;              /* line type                                    */

int XTRANSL = 0;
int YTRANSL = 0;

int XClipSet = 0;           /* if clipping set                              */
int XArcFlg = 0;            /* if arc                                       */
int XScalFlg = 0;           /* set for scaling                              */
double XXScal,XYScal;       /* scaling factors                              */

#if S_UNIX
#define icon_bitmap_width 40
#define icon_bitmap_height 40
static char icon_bitmap_bits[] = {
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
   0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00};
XIconSize *size_list;
Pixmap icon_pixmap;
XSizeHints *size_hints;
XClassHint *class_hints;
XWMHints *wm_hints;
XTextProperty windowName, iconName;
char *progname  = "TDA plot window";
char *icon_name = "TDA";
#endif

/* ------------------------------------------------------------------------ */
/*  x_show()    create a plot window and show current plot.                 */
/*                                                                          */
/*              xshow(                                                      */
/*                  dscal=...,      scaling factor, def. 1.0                */
/*                  font=...,       font, def. 7x13 [not implemented,yet]   */
/*                                                                          */
/*              ) [= PostScript_File];                                      */
/*                                                                          */
/*              Note: the basic ideas are taken from examples in            */
/*              Arian Nye, XLib Programming Manual, O'Reilly & Ass.         */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int x_show(void)
{
    int err,n,fnd,x,y,sw,sh,swm,shm,done,oflag;
    unsigned int bw,dr;
    /** char errt[101]; **/

    err = -1;                
    XClipSet = XTRANSL = YTRANSL = oflag = 0;
    XLWID = XLTYP = 1;

    if (check_cmd(1))
        return(-1);

    printf1("Open plot window. Current memory: %d bytes.\n",MemReq);


    if (parm(CmdBuf + 5,2,0))     /* get parameters */
        goto XSHOWFin;

    if (DScalFlg == 0 || DScal <= EPSI1)
        DScal = 1.0;
    else
        DScal = sqrt(DScal);
   
    if (alloc_acc(RLMaxDef))
        goto XSHOWFin;

    if (PMFDef == 0) {
        if (PSFFlg != 2) {
            printf1("Error: no current plot file.\n");
            goto XSHOWFin;
        }
        printf1("Using current plot file: %s\n",PSFName);
        fflush(PSFd);
        if (!(PMFd = fopen(PSFName,OPEN_RD))) {   
            printf1("Error: can't reopen this file.\n");
            goto XSHOWFin;
        }
        oflag = 1;

        XBLX = BBLX - 5;
        XBLY = BBLY - 5;
        XBUX = BBUX + 5;
        XBUY = BBUY + 5;
        fnd = 2;
    }
    else {
        printf1("Using PostScript file: %s\n",PMFdName);
        fnd = 0;
        while (fgets(AcC,RLMaxDef,PMFd)) {
            if (!strncmp(AcC,"%%Creator: TDA",14))  
                fnd = 1;

            if (fnd == 1 && sscanf(AcC + 2,"BoundingBox: %d %d %d %d",&XBLX,&XBLY,&XBUX,&XBUY) == 4) {
                fnd = 2;
                break;
            }
        }
    }
    if (fnd != 2) {
        printf1("Error: probably not a TDA PostScript file.\n");
        goto XSHOWFin;
    }
    XLEN = (double)(XBUX - XBLX) / PtMM;
    YLEN = (double)(XBUY - XBLY) / PtMM;

    if (XLEN <= EPSI1 || YLEN <= EPSI1) {
        printf1("Error: no valid bounding box.\n");
        goto XSHOWFin;
    }
    printf1("PostScript bounding box: %d %d %d %d\n",XBLX,XBLY,XBUX,XBUY);

    /* allocate memory for xpoints structure */

    n = MAX_XPOINTS;    
                        
    /************************************************************
     n = XMaxRequestSize(display);     
     n = (n - 3) / 2;        
     if (n < 10) {
         printf1("Error: cannot allocate XPoint structure.\n");
         goto XSHOWFin;
     }
    ************************************************************/

    if (!(xpoints = (XPoint *) calloc(n,sizeof(XPoint)))) {
        p_err(-2,1);
        goto XSHOWFin;
    }
    XPAllocN = n;
    memrq(n,sizeof(XPoint));

    if (XSHOWFIRST) {

        if ( (display=XOpenDisplay("")) == NULL) {
            printf1("Error: cannot open display.\n");
            goto XSHOWFin;
        }
        screen = DefaultScreen(display);
        root   = RootWindow(display, screen);

        black = BlackPixel(display, screen);
        white = WhitePixel(display, screen);
   
        sw = DisplayWidth(display,screen);
        sh = DisplayHeight(display,screen);
        swm = DisplayWidthMM(display,screen);
        shm = DisplayHeightMM(display,screen);

        XPPMM = (double)sw / (double)swm;
        YPPMM = (double)sh / (double)shm;

        xwidth = (int)(DScal * XLEN * XPPMM + 0.5);
        xheight = (int)(DScal * YLEN * YPPMM + 0.5);


        win = XCreateSimpleWindow(display, root,
                  0, 0, xwidth, xheight, 2, white,white );

        XSelectInput(display, win, ExposureMask | KeyPressMask |
                                ButtonPressMask | StructureNotifyMask);
  
#if S_UNIX
        if (!(size_hints = XAllocSizeHints()) ||           
            !(wm_hints = XAllocWMHints())     ||
            !(class_hints = XAllocClassHint())) ;
          
        else {

            XStringListToTextProperty(&progname,1,&windowName);        
            XStringListToTextProperty(&icon_name,1,&iconName);         

            icon_pixmap = XCreateBitmapFromData(display,win,icon_bitmap_bits, 
                                      icon_bitmap_width, icon_bitmap_height);

            size_hints->flags = PPosition | PSize | PMinSize;
            size_hints->min_width = xwidth / 5;              
            size_hints->min_height = xheight / 5;          
   
            wm_hints->initial_state = NormalState;
            wm_hints->input = True;
            wm_hints->icon_pixmap = icon_pixmap;    
            wm_hints->flags = StateHint | IconPixmapHint | InputHint;      

            class_hints->res_name = progname;
            class_hints->res_class = progname;
           
            XSetWMProperties(display,win,&windowName,&iconName, 
                         &progname,0,size_hints, wm_hints,class_hints);
        }

        if (XFontD == NULL) {
            if (!(XFontD = (char *)calloc(strlen(XFontDef)+1,sizeof(char)))) {           
                p_err(-2,1);
                goto XSHOWFin;
            }
            strcpy(XFontD,XFontDef);
        }
        if ((font_info = XLoadQueryFont(display,XFontD)) == NULL) {
            printf1("Error: cannot open font: %s\n",XFontD);
            goto XSHOWFin;
        }
#endif

/*      XMapWindow(display,win);  */
        gc = XCreateGC(display, win, (long)0, &gc_values);

#if S_UNIX
        default_visual = XDefaultVisual(display,screen);
        default_cmap = XDefaultColormap(display,screen);
#endif
        x_ncolors = 0;
        for (n = 0; n < MAX_COLORS; ++n) {
            if (!XParseColor(display,default_cmap,ColName[n],&exact_def)) {
                continue;                  
            }
            if (!XAllocColor(display,default_cmap,&exact_def)) {
                continue;
            }
            x_colors[x_ncolors++] = exact_def.pixel;
        }
    }
    else {
        xwidth = (int)(DScal * XLEN * XPPMM + 0.5);
        xheight = (int)(DScal * YLEN * YPPMM + 0.5);

        XResizeWindow(display,win,xwidth,xheight);
        /*****
        XClipSet = 1;
            x_offclip();
        ***/
    }
    XMapWindow(display,win);
    XSetForeground(display,gc,black);     
    /* XClearWindow(display,win); */ 

    XFlush(display);       

    done = 0;
    while(!done) {
        XNextEvent(display, &event);

        switch(event.type) {
             
            case Expose:
                if (event.xexpose.count != 0)
                    break;

                XGetGeometry(display,win,&root,&x,&y,&xwidth,&xheight,&bw,&dr);

                if (x_graph())          /* call plot function */
                    goto XSHOWFin1;

                XSetInputFocus(display,win,RevertToParent,CurrentTime);
                XFlush(display);        
                break;

            case ConfigureNotify:
                break;

            case ButtonPress:
                break;

            case KeyPress:
                XUnmapWindow(display,win);  
                XFlush(display);        
                done = 1;
                break;

            default:
                break;
        }
    }
XSHOWFin1:
    XUnmapWindow(display,win);  
    XFlush(display);        
    XSHOWFIRST = 0;

XSHOWFin:
    if (oflag)
        fclose(PMFd);

    if (XPAllocN > 0) {
        free((char *)xpoints);
        memrq(-XPAllocN,sizeof(XPoint));
        XPAllocN = 0;
    }      
    p_clean();         
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  x_graph()   plot current state of PostScript file.                      */

/*                                                                          */
/*  Return 0 if OK, -1 if error                                             */

int x_graph(void)
{
    int n,w,h,xpnt,ypnt,xorg,yorg,newp,mode;
    double x,y,r,a1,a2;
    register char *p;

    XScalFlg = XArcFlg = XClipSet = XTRANSL = YTRANSL = 0;
    XLWID = XLTYP = 1;
    xpnt = XBUX - XBLX;
    ypnt = XBUY - XBLY;

    XClipSet = 1;
        x_offclip();

    fseek(PMFd,0,0);
    xorg = yorg = -1;
    while (fgets(AcC,RLMaxDef,PMFd)) {

        if (xorg < 0) {
            if (sscanf(AcC,"/xorg %d",&xorg) == 1) ;
        }
        if (yorg < 0) {
            if (sscanf(AcC,"/yorg %d",&yorg) == 1)  
            break;
        }
    }
    if (xorg < 0 || yorg < 0) {
        printf1("Error: can't find origin.\n");
        return(-1);
    }
    XPPnt = (double)xwidth / (double)xpnt;
    YPPnt = (double)xheight / (double)ypnt;

    XNXOrg = (int)((double)(xorg - XBLX) * XPPnt);
    XNYOrg = xheight - (int)((double)(yorg - XBLY) * YPPnt);

    mode = CoordModeOrigin;
    n = 0;
    newp = 1;
    while (fgets(AcC,RLMaxDef,PMFd)) {
        p = skip_b(AcC);   
        if (sscanf(p,"%lf %lf",&x,&y) == 2) {
            p = skip_dbl(p);
            p = skip_b(p);
            p = skip_dbl(p);
            p = skip_b(p);
            if (*p == 'm') {
                n = newp = 1;
                xpoints[0].x = x_getx(x,1);
                xpoints[0].y = x_gety(y,1);
            }
            else if (*p == 'l') {
                if (newp >= 1) {
                    if (n < XPAllocN) {
                        xpoints[n].x = x_getx(x,1);
                        xpoints[n].y = x_gety(y,1);
                        n++;
                    }
                    if (newp == 1) {
                        newp= 2;
                        mode = CoordModeOrigin;
                    }
                }
            }
            else if (*p == 'r' && *(p + 1) == 'l') {
                if (newp >= 1) {
                    if (n < XPAllocN) {
                        xpoints[n].x = x_getx(x,0);
                        xpoints[n].y = x_gety(y,0);
                        n++;
                    }
                    if (newp == 1) {
                        newp= 2;
                        mode = CoordModePrevious;
                    }
                }
            }
            else if (!strncmp(p,"translate",9)) {
                XTRANSL = YTRANSL = 0;
                XTRANSL = x_getx(x,1);
                YTRANSL = x_gety(y,1);
            }
            else if (!strncmp(p,"scale",5)) {
                XXScal = x;            
                XYScal = y;            
                XScalFlg = 1;
            }
            else if (sscanf(p,"%lf %lf %lf",&r,&a1,&a2) == 3) {
                p = skip_dbl(p);
                p = skip_b(p);
                p = skip_dbl(p);
                p = skip_b(p);
                p = skip_dbl(p);
                p = skip_b(p);
                if (!strncmp(p,"arc",3)) {
                    xpoints[0].x = x_getx(x,1);
                    xpoints[0].y = x_gety(y,1);
                    XArcFlg = 1;
                }
            }
        }
        else if (sscanf(p,"%lf",&x) == 1) {
            p = skip_dbl(p);
            p = skip_b(p);
            if (!strncmp(p,"setlinewidth",12)) {
                XLWID = (int)(DScal * XPPMM * x / PtMM + 0.5);
                if (XLWID < 1)
                    XLWID = 1;
            }
            else if (!strncmp(p,"setgray fill",12)) {
                x_fill(n,mode,x);
            }
        }
        else if (*p == '[') {                       /* line type */
            XLTYP = 1;
            if (sscanf(p,"[%d %d]",&w,&h) == 2) {
                if (w == 1 && h == 2)  
                    XLTYP = 2;
                else if (w == 1 && h == 3)  
                    XLTYP = 3;
                else if (w == 1 && h == 5)  
                    XLTYP = 4;
                else if (w == 2 && h == 2)  
                    XLTYP = 5;
                else if (w == 4 && h == 2)  
                    XLTYP = 6;
                else if (w == 6 && h == 4)  
                    XLTYP = 7;
                else if (w == 6 && h == 2)  
                    XLTYP = 8;
                else if (w == 9 && h == 4)  
                    XLTYP = 9;
            }
        }
        else if (!strncmp(p,"%#text:",7)) {
            x_text(p + 1);
        }
        else if (!strncmp(p,"%#symbol:",9)) {
            x_symbol(p + 1);
        }
        else if (!strncmp(p,"stroke",6)) {
            if (XArcFlg) {
                x_arc(r,a1,a2);
            }
            else if (newp == 2) {
                x_draw(n,mode);
                if (XClipSet)
                    x_offclip();
            }
            XScalFlg = XArcFlg = n = 0;
    
        }
        else if (!strncmp(p,"%#stroke",8)) {
           if (newp == 2) {
               x_draw(n,mode);
               n = 0;
           }
        }
        else if (!strncmp(p,"show",4))  
            XScalFlg = XTRANSL = YTRANSL = 0;
        else if (!strncmp(p,"grestore",8))  
            XScalFlg = XTRANSL = YTRANSL = n = 0;
        else if (!strncmp(p,"clip",4))  {
            x_setclip(n);
            n = 0;
        }
    }
    if (XClipSet)      
        x_offclip();

    XScalFlg = XArcFlg = XClipSet = XTRANSL = YTRANSL = 0;
    XLWID = XLTYP = 1;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  x_getx(x,mode)   get XWin coordinate for x.                             */

int x_getx(double x,int mode)
{
    int n;

    n = (int)((x * XPPnt) + 0.5);
    if (mode > 0) {
        if (XTRANSL != 0)
            n += XTRANSL;
        else if (mode)                   
            n += XNXOrg;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  x_gety(y,mode)   get XWin coordinate for y.                             */

int x_gety(double y,int mode)
{
    int n = 0;

    if (mode > 0) {
        if (YTRANSL != 0)
            n = YTRANSL;
        else if (mode)
            n = XNYOrg;
    }
    n -= (int)((y * YPPnt) + 0.5);
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  x_draw(n,mode)   plot line with n points.                               */

void x_draw(int n,int mode)
{
    int cn,ls = LineSolid;
    char dlist[4];

    if (n < 2)
        return;

    /******************************
    if (XLTYP > 1) { 
        if (XLTYP <= 5)
            ls = LineOnOffDash;
        else              
            ls = LineDoubleDash;
    }
    *********************************/

    XSetLineAttributes(display,gc,XLWID,ls,CapButt,JoinRound);

    if (XLTYP > 1) { 
        cn = XLTYP - 2;
        if (cn < x_ncolors)
            XSetForeground(display,gc,x_colors[cn]);
        else {
          switch (XLTYP) {
            case 2: dlist[0] =  3; dlist[1] = 3; break;
            case 3: dlist[0] =  4; dlist[1] = 4; break;
            case 4: dlist[0] =  4; dlist[1] = 6; break;
            case 5: dlist[0] =  5; dlist[1] = 3; break;
            case 6: dlist[0] =  7; dlist[1] = 4; break;
            case 7: dlist[0] =  9; dlist[1] = 6; break;
            case 8: dlist[0] = 10; dlist[1] = 4; break;
            case 9: dlist[0] = 12; dlist[1] = 5; break;
          }
          XSetDashes(display,gc,0,dlist,2);
        }
    }
    XDrawLines(display,win,gc,xpoints,n,mode);           
    if (XLTYP > 1)
        XSetForeground(display,gc,black);
}

/* ------------------------------------------------------------------------ */
/*  x_fill(n,mode,g)    fill with g.                                        */

void x_fill(int n,int mode,double g)
{
    if (n < 3)
        return;
    
    if (g > 0.5) {
        XSetForeground(display,gc,white);
#if S_UNIX
        XSetBackground(display,gc,black);                
#endif
    }
        
    
    XSetLineAttributes(display,gc,1,LineSolid,CapButt,JoinRound);

    XFillPolygon(display,win,gc,xpoints,n,Nonconvex,mode);  
              
    if (g > 0.5) {
        XSetForeground(display,gc,black);
#if S_UNIX
        XSetBackground(display,gc,white);                
#endif
    }
}

/* ------------------------------------------------------------------------ */
/*  x_arc(n,mode)   draw an arc.                                            */

void x_arc(double r,double a1,double a2)
{
    int irx,iry,ix,iy,ia1,ia2;
    double rx,ry;

    rx = ry = r;
    if (XScalFlg) {
        rx *= XXScal;
        ry *= XYScal;
    }
    irx = iabs(x_getx(rx,-1));      
    iry = iabs(x_gety(ry,-1));      
    ix = iabs(xpoints[0].x - irx);          
    iy = iabs(xpoints[0].y - iry);               
    ia1 = (int)a1 * 64;
    ia2 = (int)(a2 - a1) * 64;

    XDrawArc(display,win,gc,ix,iy,2 * irx,2 * iry,ia1,ia2);
}

/* ------------------------------------------------------------------------ */
/*  x_text(s)   plot text.                                                  */

void x_text(char *s)
{
    int adj,r,ix,iy,slen,w;
    double x,y,siz;
    register char *p;

    if (sscanf(s,"#text: %lf %lf %lf %d %d",&x,&y,&siz,&adj,&r) != 5)
        return;

    if (r != 0)         /* cannot rotate strings */
        return;

    s = skip_b(s + 6);
    s = skip_dbl(s);
    s = skip_b(s);
    s = skip_dbl(s);
    s = skip_b(s);
    s = skip_dbl(s);
    s = skip_b(s);
    s = skip_int(s);
    s = skip_b(s);
    s = skip_int(s);
    s = skip_b(s);
    p = s + strlen(s);
    while (--p > s) {
        if (*p == CR || *p == LF)
            *p = '\0';
        else
            break;
    }
    ix = x_getx(x,1);
    iy = x_gety(y,1);
    slen = strlen(s);
    if (adj > 0) {  

#if S_UNIX
        w = XTextWidth(font_info,s,slen);                   
#endif
#if S_DOS                   /* use fixed character width */
        w = 7 * slen;
#endif  
        if (adj == 1)                       /* center */
            ix -= w / 2;
        else if (adj == 2)                  /* right justified */
            ix -= w;
    }
#if S_UNIX
    XSetBackground(display,gc,white);                
#endif
    XDrawImageString(display,win,gc,ix,iy,s,slen);
}

/* -##--------------------------------------------------------------------- */
/*  x_symbol(s)   plot symbol.                                              */

void x_symbol(char *s)
{
    int n,cn,ix,iy,ixx,iyy,iw,iww;
    double x,y,siz;

    if (sscanf(s,"#symbol: %d %lf %lf %lf",&n,&x,&y,&siz) != 4 || siz <= 0.0)
        return;

    ix = x_getx(x,1);
    iy = x_gety(y,1);
    iww = iabs(x_getx(siz,-1));
    iw = iabs(x_getx(2.0 * siz,-1));
    ixx = ix - iww;
    iyy = iy - iww;
           
    XSetForeground(display,gc,white);
    
    XFillRectangle(display,win,gc,ixx,iyy,iw,iw);  
    
    cn = -1;
    if (n > 1) 
        cn = (n - 2) % 8;

    if (cn >= 0 && cn < x_ncolors)
        XSetForeground(display,gc,x_colors[cn]);
    else {
        XSetForeground(display,gc,black);
        cn = -1;
    }

    if (n == 1 || n == 3 || n == 6 || n == 14 || n == 17) {
        xpoints[0].x = ixx;
        xpoints[0].y = iy;
        xpoints[1].x = ixx + iw;
        xpoints[1].y = iy;
        XDrawLines(display,win,gc,xpoints,2,CoordModeOrigin);
        xpoints[0].x = ix;
        xpoints[0].y = iyy;
        xpoints[1].x = ix;
        xpoints[1].y = iyy + iw;
        XDrawLines(display,win,gc,xpoints,2,CoordModeOrigin);
    }
    if (n == 2 || n == 3 || n == 7) {
        xpoints[0].x = ixx;
        xpoints[0].y = iyy;
        xpoints[1].x = ixx + iw;
        xpoints[1].y = iyy + iw;
        XDrawLines(display,win,gc,xpoints,2,CoordModeOrigin);
        xpoints[0].x = ixx;
        xpoints[0].y = iyy + iw;
        xpoints[1].x = ixx + iw;
        xpoints[1].y = iyy;
        XDrawLines(display,win,gc,xpoints,2,CoordModeOrigin);
    }
    if (n == 4 || n == 5 || n == 6 || n == 7) {
        XDrawArc(display,win,gc,ixx,iyy,iw,iw,0,360 * 64);
        if (n == 5) {
        /*  XSetForeground(display,gc,black); */
            XFillArc(display,win,gc,ixx,iyy,iw,iw,0,360 * 64);
        }
    }
    if (n == 8 || n == 14) 
        XDrawRectangle(display,win,gc,ixx,iyy,iw,iw);

    if (n == 9)
        XFillRectangle(display,win,gc,ixx,iyy,iw,iw);  

    if (n == 10 || n == 11) {
        xpoints[0].x = ixx;
        xpoints[0].y = iyy + iw;
        xpoints[1].x = ix;
        xpoints[1].y = iyy;
        xpoints[2].x = ixx + iw;
        xpoints[2].y = iyy + iw;
        xpoints[3].x = ixx;
        xpoints[3].y = iyy + iw;
        XDrawLines(display,win,gc,xpoints,4,CoordModeOrigin);
        if (n == 11)  
            XFillPolygon(display,win,gc,xpoints,4,Convex,CoordModeOrigin);
    }
    if (n == 12 || n == 13) {
        xpoints[0].x = ixx;
        xpoints[0].y = iyy;
        xpoints[1].x = ixx + iw;
        xpoints[1].y = iyy;
        xpoints[2].x = ix;
        xpoints[2].y = iyy + iw;
        xpoints[3].x = ixx;
        xpoints[3].y = iyy;
        XDrawLines(display,win,gc,xpoints,4,CoordModeOrigin);
        if (n == 13)  
            XFillPolygon(display,win,gc,xpoints,4,Convex,CoordModeOrigin);
    }
    if (n == 15 || n == 16 || n == 17) {
        xpoints[0].x = ixx;
        xpoints[0].y = iy;
        xpoints[1].x = ix;
        xpoints[1].y = iyy;
        xpoints[2].x = ixx + iw;
        xpoints[2].y = iy;
        xpoints[3].x = ix;
        xpoints[3].y = iyy + iw;
        xpoints[4].x = ixx;
        xpoints[4].y = iy;
        XDrawLines(display,win,gc,xpoints,5,CoordModeOrigin);
        if (n == 16)  
            XFillPolygon(display,win,gc,xpoints,5,Convex,CoordModeOrigin);
    }
    if (cn >= 0)
        XSetForeground(display,gc,black);
}
   
/* ------------------------------------------------------------------------ */
/*  x_setclip    set clipping region.                                       */

void x_setclip(int n)
{
         
    if (n < 4)   
        return;
 
    xrec[0].x = xpoints[1].x;
    xrec[0].y = xpoints[1].y;
    xrec[0].width = iabs(xpoints[2].x - xpoints[1].x);
    xrec[0].height = iabs(xpoints[0].y - xpoints[1].y);

    XSetClipRectangles(display,gc,0,0,xrec,1,Unsorted);
    XClipSet = 1;
}

/* ------------------------------------------------------------------------ */
/*  x_offclip    reset clipping to full area.                               */

void x_offclip(void)
{
       
    if (XClipSet == 0)
        return;
 
    xrec[0].x = 0;
    xrec[0].y = 0;
    xrec[0].width = xwidth;                              
    xrec[0].height = xheight;

    XSetClipRectangles(display,gc,0,0,xrec,1,Unsorted);
    XClipSet = 0;
}

/*****************************************************************************
    x_dfonts(s)  display available fonts (matching string s)                  

void x_dfonts(char *s)
{
    register int i;
    int fh,fw;
    int n = 0;                number of fonts    
    char **fnames;            array with font names   
    XFontStruct *font_i;
    XCharStruct mx_b;
    char *p; 

    fnames = XListFonts(display,s,XFontNMax,&n);
    if (n <= 0) {
       	printf1("No fonts found.\n");
       	return;
    }
    i = 0;
    while (i < n) {
       	p = fnames[i];
       	printf1("Font %4d of %4d: %s\n",i + 1,n,p);

       	if ((font_i = XLoadQueryFont(display,p)) != NULL) {
       	    mx_b = font_i->max_bounds;
       	    fw = mx_b.width;
       	    fh = mx_b.ascent + mx_b.descent;
       	    printf1("Font width: %d  height: %d\n",fw,fh);
        }
       	i++;
    }
}
*****************************************************************/

#endif  /* end of #if S_XWIN */


