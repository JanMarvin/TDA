/****************************************************************************/
/*  t_plg                                                                   */
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
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_gf.h"
#include "t_plot.h"
/*   #include "t_graph.h"   */

/*  functions in t_plg.c */

int pl_plg(void);
void plg_lab(double ic,double fs,double xi,double yi); 
void plg_nodes(int n);
void plg_arrow(double x,double y,double sa,double sb,int dir,double rot);
void plg_init(void);
void plg_fin(void);
void plg_line(double xa,double ya,double xb,double yb);
void plg_arc(double x,double y,int na,int nb,double r,int s);

/*--------------------------------------------------------------------------*/
/*  pl_plg()    Plot graph.                                                 */
/*                                                                          */
/*              plg(                                                        */
/*                node(gt=,n=,str=,rd=,lt=,lw=,gs=,fs=) = x,y,              */
/*                edge(a=,ic=,fs=,lt=,lw=,dir=,rd=) = i,j,                  */
/*              );                                                          */
/*                                                                          */
/*              Return 0 if OK, -1 if error.                                */

int pl_plg(void)
{
    register int k;
    int err,eerr,ii,nmax,nn,i,j,ptyp,aflag,gti,gtj;
    register char c,*p,*q;
    double x,y,xi,yi,xj,yj,alpha,tmp,tmp1,tmp2;
    double a,rot,dx,dy,r,rai,raj,rbi,rbj,sa,sb;

    nmax = 100;        /* default max number of nodes and edges */
           
    err = -1;
    if (check_pcmd(0,2))
        return(-1);

    printf1("Plot graph. Current memory: %d bytes.\n",MemReq);

    p = CmdBuf;
    if (sscanf(p + 4,"nmax=%d",&nn) == 1 && nn > 0) {
        nmax = nn;
        p = skip_int(p + 9);
    }
    printf1("Max number of nodes: %d\n",nmax);

    if (alloc_acns(nmax + 1))   /* type of node */
        goto PLGFin;

    if (alloc_acx(nmax))    /* x coordinate of node */
        goto PLGFin;

    if (alloc_acy(nmax))    /* y coordinate of node */
        goto PLGFin;

    if (alloc_acxf(nmax))   /* rda */
        goto PLGFin;

    if (alloc_acyf(nmax))   /* rdb */
        goto PLGFin;

    if (alloc_acv(nmax))    /* font size */
        goto PLGFin;

    if (alloc_acn(nmax))    /* node number */
        goto PLGFin;

    if (alloc_acm(nmax))    /* line type */
        goto PLGFin;

    if (alloc_acw(nmax))    /* line width */
        goto PLGFin;

    if (alloc_actmp(nmax))  /* grey value */
        goto PLGFin;

    if (alloc_acptr(nmax))  /* label */
        goto PLGFin;

    fprintf(PSFd,"\n%%#%d: graph\n",++PSONUM);

    nn = 0;
    for (ii = 0; ii < 2; ++ii) {        /* for nodes and edges */

        p = CmdBuf;     
        while (*p) {
            eerr = 0;
            if ((ii == 0 && !strncmp(p,"node",4)) || (ii == 1 && !strncmp(p,"edge",4))) {

                q = p + 4;
                if (*q != '(') {
                    p_err(-1,1);
                    goto PLGFin;
                }
                while (*++q && *q != ')') ;

                if (*q++ != ')' || *q++ != '=') {
                    p_err(-1,1);
                    goto PLGFin;
                }
                q = skip_dbl(q);
                q = skip_dbl(q + 1);
                if (*q != ',' && *q++ != ')') {
                    p_err(-1,1);
                    goto PLGFin;
                }
                c = *q;
                *q = '\0';
                printf1("> %s\n",p);
                if (parm(p + 4,5,1))    /* get parameters */
                    goto PLGFin;

                if (nn >= nmax) { 
                    printf1("Error: exceeded max number of nodes.\n");
                    goto PLGFin;
                }
                if (ii == 0) {      /* nodes */

                    if (PMGT < 1 || PMGT > 4)
                        PMGT = 1;

                    if (PMSTRA > 0) {
                        k = strlen(PMSTR) + 1;
                        if (!(AcPtr[nn]  = (char *)calloc(k,sizeof(char)))) {
                            p_err(-2,1);
                            goto PLGFin;
                        }
                        memrq(k,sizeof(char));
                        strcpy(AcPtr[nn],PMSTR);
                    }
                    if (PMRDA < EPSI1)
                        PMRDA = 4.0;

                    if (PMGT <= 2)
                        PMRDB = PMRDA;
                    else if (PMRDB < EPSI1)
                        PMRDB = 2.0 * PMRDA;
                         
                    AcNS[nn] = PMGT;
                    AcX[nn] = PMRHSA;
                    AcY[nn] = PMRHSB;
                    AcXF[nn] = PMRDA / 2.0;
                    AcYF[nn] = PMRDB / 2.0;
                    AcN[nn] = PMN;
                    AcM[nn] = PMLT;
                    AcW[nn] = PMLW;
                    if (PMFSFlg)
                        AcV[nn] = PMFS;
                    else
                        AcV[nn] = PMRDA / 2.0;

                    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0)
                        AcTmp[nn] = PMGS;
                    else
                        AcTmp[nn] = 1.0;
                    nn++;
                }
                else {          /* plot edge */

                    eerr = 1;
                    i = j = -1;
                    for (k = 0; k < nn; ++k) {
                        if (AcN[k] == (int)PMRHSA)
                            i = k;
                        if (AcN[k] == (int)PMRHSB)
                            j = k;
                    }
                    if (i < 0 || j < 0)  
                        goto PLGCONT;

                    gti = AcNS[i];
                    gtj = AcNS[j];
                    xi = ps_2dx(AcX[i]);
                    yi = ps_2dy(AcY[i]);
                    xj = ps_2dx(AcX[j]);
                    yj = ps_2dy(AcY[j]);
                    rai = AcXF[i] * PtMM;
                    rbi = AcYF[i] * PtMM;
                    raj = AcXF[j] * PtMM;
                    rbj = AcYF[j] * PtMM;

                    if (PMAFlg && PMA1 > 0.0 && PMA2 > 0.0) {
                        aflag = 1;
                        a = PMA1 * PtMM;
                        sa = PMA1 * PtMM / 7.0;
                        sb = PMA2 * PtMM / 4.0;
                    }
                    else {
                        aflag = 0;
                        a = 0.0;
                    }
                    if (i == j) {       /*** ## loop ***/

                        if (PMRDA <= 0.0)
                            PMRDA = 1.5 * AcXF[i];

                        r = PMRDA * PtMM;
                        tmp = sqrt(0.5);

                        if (PMDIR == 1) {
                            xi -= (rbi - rai);
                            plg_init();
                            plg_arc(xi - r,yi + r,0,270,r,1);
                            plg_line(xi,yi,xi - r,yi);
                            plg_line(xi,yi,xi,yi + r);
                            plg_fin();
                            if (aflag)
                                plg_arrow(xi,yi + rai,sa,sb,-1,270.0);
                            if (PMICFlg)                           
                                plg_lab(PMIC,PMFS,xi - 2.0 * r,yi + r);
                        }
                        else if (PMDIR == 2) {
                            xi -= (rbi - rai);
                            plg_init();
                            plg_arc(xi - r,yi - r,90,360,r,1);
                            plg_line(xi,yi,xi - r,yi);
                            plg_line(xi,yi,xi,yi - r);
                            plg_fin();
                            if (aflag)
                                plg_arrow(xi,yi - rai,sa,sb,-1,90.0);
                            if (PMICFlg)                           
                                plg_lab(PMIC,PMFS,xi - 2.0 * r,yi - r);
                        }
                        else if (PMDIR == 3) {
                            xi += (rbi - rai);
                            plg_init();
                            plg_arc(xi + r,yi - r,180,90,r,1);
                            plg_line(xi,yi,xi + r,yi);
                            plg_line(xi,yi,xi,yi - r);
                            plg_fin();
                            if (aflag)
                                plg_arrow(xi,yi - rai,sa,sb,-1,90.0);
                            if (PMICFlg)                           
                                plg_lab(PMIC,PMFS,xi + 2.0 * r,yi - r);
                        }
                        else {              
                            xi += (rbi - rai);
                            plg_init();
                            plg_arc(xi + r,yi + r,270,180,r,1);
                            plg_line(xi,yi,xi + r,yi);
                            plg_line(xi,yi,xi,yi + r);
                            plg_fin();
                            if (aflag)
                                plg_arrow(xi,yi + rai,sa,sb,-1,270.0);
                            if (PMICFlg)                           
                                plg_lab(PMIC,PMFS,xi + 2.0 * r,yi + r);
                        }
                    }
                    else if (PMRDA == 0.0) {  /*** ## line ***/
     
                        if (aflag) {
    
                            if (gtj == 3) {
                                if (xi < xj - rbj)
                                    xj -= rbj;
                                else if (xi > xj + rbj)
                                    xj += rbj;
                                else
                                    xj = xi;
                            }
                            dx = xj - xi;
                            dy = yj - yi;

                            if (dx != 0.0) {
                                alpha = atan(fabs(dy / dx));
                                rot = alpha * 360.0 / (2.0 * Pi);
                            }
                            else
                                alpha = rot = 0.0;

                            if (dx > 0.0) {
                                if (dy > 0.0) {
                                    if (gtj == 1 || gtj == 3) {
                                        yj -= sin(alpha) * raj;
                                        xj -= cos(alpha) * raj;
                                    }
                                    else {
                                        x = xj - rbj;
                                        y = yj - raj;
                                        if (x != xi)
                                            tmp = atan((y - yi) / (x - xi));
                                        else
                                            tmp = 0.0;

                                        if (alpha >= tmp) {
                                            xj -= rbj;
                                            yj -= tan(alpha) * rbj;
                                        }
                                        else {
                                            yj -= raj;
                                            xj -= tan(Pi / 2.0 - alpha) * raj;
                                        }
                                    }
                                }   
                                else if (dy < 0.0) {
                                    if (gtj == 1 || gtj == 3) {
                                        yj += sin(alpha) * raj;
                                        xj -= cos(alpha) * raj;
                                    }
                                    else {
                                        x = xj - rbj;
                                        y = yj + raj;
                                        if (x != xi)
                                            tmp = atan(fabs((y - yi) / (x - xi)));
                                        else
                                            tmp = 0.0;

                                        if (alpha >= tmp) {
                                            xj -= rbj;
                                            yj += tan(alpha) * rbj;
                                        }
                                        else {
                                            yj += raj;
                                            xj -= tan(Pi / 2.0 - alpha) * raj;
                                        }
                                    }
                                    rot = 360.0 - rot;
                                }   
                                else {
                                    xi += rbi;
                                    if (gtj == 3)
                                        xj -= raj;
                                    else
                                        xj -= rbj;
                                    rot = 0.0;
                                }
                            }                            
                            else if (dx < 0.0) {
                                if (dy > 0.0) {
                                    if (gtj == 1 || gtj == 3) {
                                        yj -= sin(alpha) * raj;
                                        xj += cos(alpha) * raj;
                                    }
                                    else {
                                        x = xj + rbj;
                                        y = yj - raj;
                                        if (x != xi)
                                            tmp = atan(fabs((y - yi) / (x - xi)));
                                        else
                                            tmp = 0.0;

                                        if (alpha >= tmp) {
                                            xj += rbj;
                                            yj -= tan(alpha) * rbj;
                                        }
                                        else {
                                            yj -= raj;
                                            xj += tan(Pi / 2.0 - alpha) * raj;
                                        }
                                    }
                                    rot = 180.0 - rot;
                                }   
                                else if (dy < 0.0) {
                                    if (gtj == 1 || gtj == 3) {
                                        yj += sin(alpha) * raj;
                                        xj += cos(alpha) * raj;
                                    }
                                    else {
                                        x = xj + rbj;
                                        y = yj + raj;
                                        if (x != xi)
                                            tmp = atan(fabs((y - yi) / (x - xi)));
                                        else
                                            tmp = 0.0;

                                        if (alpha >= tmp) {
                                            xj += rbj;
                                            yj += tan(alpha) * rbj;
                                        }
                                        else {
                                            yj += raj;
                                            xj += tan(Pi / 2.0 - alpha) * raj;
                                        }
                                    }
                                    rot = 180.0 + rot;
                                }   
                                else {
                                    xi -= rbi;
                                    if (gtj == 3)
                                        xj += raj;
                                    else
                                        xj += rbj;
                                    rot = 180.0;
                                }
                            }                            
                            else {
                                if (dy > 0.0) {
                                    yj -= raj;
                                    yi += rai;
                                    rot = 90.0;
                                }
                                else if (dy < 0.0) {
                                    yj += raj;
                                    yi -= rai;
                                    rot = 270.0;
                                }
                            }

                        }
                        plg_init();
                        plg_line(xi,yi,xj,yj);
                        plg_fin();

                        if (aflag)
                            plg_arrow(xj,yj,sa,sb,-1,rot);

                        if (PMICFlg)                           
                            plg_lab(PMIC,PMFS,(xi + xj) / 2.0,(yi + yj) / 2.0);
                    }

                    else if (PMRDA != 0.0) {    /*** ##  arc  ***/
     
                        r = PMRDA * PtMM;
                        ptyp = 0;
                        if (r < 0) {
                            r = -r;
                            ptyp = 1;
                        }
                        if (ptyp == 0) {
                            tmp1 = rbi;
                            tmp2 = rbj;
                            if (gti == 3)  
                                tmp1 += rai;
                            if (gtj == 3)
                                tmp2 += raj;
                             
                            if (xi + tmp1 + r <= xj && yi + r + raj <= yj) {
                                xi += tmp1;
                                yj -= raj;
                                plg_init();
                                plg_arc(xj - r,yi + r,270,360,r,0);
                                plg_line(xi,yi,xj - r,yi);
                                plg_line(xj,yi + r,xj,yj - a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,90.0);
                                if (PMICFlg) {                         
                                    if (yj - yi > xj - xi)
                                        plg_lab(PMIC,PMFS,xj,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yi);
                                }
                            }
                            else if (yi + 2.0 * r <= yj && xj + tmp2 > xi) {
                                xi += tmp1;
                                xj += tmp2;
                                x = dmax(xi,xj) + a;
                                plg_init();
                                plg_arc(x,yi + r,270,360,r,1);
                                plg_arc(x,yj - r,0,90,r,0);
                                plg_line(xi,yi,x,yi);
                                plg_line(x + r,yi + r,x + r,yj - r);
                                plg_line(xj + a / 2.0,yj,x,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,180.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,x + r,(yi + yj) / 2.0);

                            }
                            else if (xi <= xj - tmp2 - r && yi - rai - r >= yj) {
                                yi -= rai;
                                xj -= tmp2;
                                plg_init();
                                plg_arc(xi + r,yj + r,180,270,r,0);
                                plg_line(xi,yi,xi,yj + r);
                                plg_line(xi + r,yj,xj - a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,0.0);
                                if (PMICFlg) {                         
                                    if (yi - yj > xj - xi)
                                        plg_lab(PMIC,PMFS,xi,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yj);
                                }
                            }
                            else if (xi + 2.0 * r <= xj) {
                                yi -= rai;
                                yj -= raj;
                                if (yi >= yj)
                                    y = yj - a;
                                else
                                    y = yi + a;
                                plg_init();
                                plg_arc(xi + r,y,180,270,r,1);
                                plg_arc(xj - r,y,270,360,r,0);
                                plg_line(xi,yi,xi,y);
                                plg_line(xi + r,y - r,xj - r,y - r);
                                plg_line(xj,yj - a / 2.0,xj,y);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,90.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,(xi + xj) / 2.0,y - r);

                            }
                            else if (xi - tmp1 >= xj + r && yi >= yj + raj + r) {
                                xi -= tmp1;
                                yj += raj;
                                plg_init();
                                plg_arc(xj + r,yi - r,90,180,r,0);
                                plg_line(xi,yi,xj + r,yi);
                                plg_line(xj,yi - r,xj,yj + a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,270.0);
                                if (PMICFlg) {                         
                                    if (yi - yj > xi - xj)
                                        plg_lab(PMIC,PMFS,xj,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yi);
                                }
                            }
                            else if (yj + 2.0 * r <= yi) {
                                xi -= tmp1;
                                xj -= tmp2;
                                x = dmin(xi,xj) - a;
                                plg_init();
                                plg_arc(x,yi - r,90,180,r,1);
                                plg_arc(x,yj + r,180,270,r,0);
                                plg_line(xi,yi,x,yi);
                                plg_line(x - r,yi - r,x - r,yj + r);
                                plg_line(x,yj,xj - a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,0.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,x - r,(yi + yj) / 2.0);
                            }
                            else if (yi + rai + r <= yj && xj + tmp2 <= xi) {
                                yi += rai;
                                xj += tmp2;
                                plg_init();
                                plg_arc(xi - r,yj - r,0,90,r,0);
                                plg_line(xi,yi,xi,yj - r);
                                plg_line(xi - r,yj,xj + a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,180.0);
                                if (PMICFlg) {                         
                                    if (yj - yi > xi - xj)
                                        plg_lab(PMIC,PMFS,xi,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yj);
                                }   
                            }
                            else if (xj + 2.0 * r <= xi) {
                                yi += rai;
                                yj += raj;
                                y = dmax(yi,yj) + a;
                                plg_init();
                                plg_arc(xi - r,y,0,90,r,1);
                                plg_arc(xj + r,y,90,180,r,0);
                                plg_line(xi,yi,xi,y);
                                plg_line(xi - r,y + r,xj + r,y + r);
                                plg_line(xj,y,xj,yj + a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,270.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,(xi + xj) / 2.0, y + r);
                            }
                            else     
                                goto PLGCONT;
                        }
                        else {                      /* ## ptyp = 1 */
                            tmp1 = rbi;
                            tmp2 = rbj;
                            if (gti == 3)  
                                tmp1 += rai;
                            if (gtj == 3)
                                tmp2 += raj;

                            if (yi + rai + r <= yj && xi <= xj - tmp2 - r) {
                                yi += rai;
                                xj -= tmp2;
                                plg_init();
                                plg_arc(xi + r,yj - r,90,180,r,0);
                                plg_line(xi,yi,xi,yj - r);
                                plg_line(xi + r,yj,xj - a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,0.0);
                                if (PMICFlg) {                         
                                    if (yj - yi > xj - xi)
                                        plg_lab(PMIC,PMFS,xi,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yj);
                                }
                            }
                            else if (xi + tmp1 + r <= xj && yj + rbj + r <= yi) {
                                xi += tmp1;
                                yj += raj; 
                                plg_init();
                                plg_arc(xj - r,yi - r,0,90,r,0);
                                plg_line(xi,yi,xj - r,yi);
                                plg_line(xj,yi - r,xj,yj + a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,270.0);
                                if (PMICFlg) {                         
                                    if (yi - yj > xj - xi)
                                        plg_lab(PMIC,PMFS,xj,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yi);
                                }
                            }
                            else if (xj + tmp2 + r <= xi && yj + r <= yi - rai) {
                                yi -= rai;
                                xj += tmp2;
                                plg_init();
                                plg_arc(xi - r,yj + r,270,360,r,0);
                                plg_line(xi,yi,xi,yj + r);
                                plg_line(xi - r,yj,xj  + a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,180.0);
                                if (PMICFlg) {                         
                                    if (yi - yj > xi - xj)
                                        plg_lab(PMIC,PMFS,xi,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yj);
                                }
                            }
                            else if (xj <= xi - r - tmp1 && yi <= yj - r - raj) {
                                xi -= tmp1;
                                yj -= raj;
                                plg_init();
                                plg_arc(xj + r,yi + r,180,270,r,0);
                                plg_line(xi,yi,xj + r,yi);
                                plg_line(xj,yi + r,xj,yj - a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,90.0);
                                if (PMICFlg) {                         
                                    if (yj - yi > xi - xj)
                                        plg_lab(PMIC,PMFS,xj,(yi + yj) / 2.0);
                                    else
                                        plg_lab(PMIC,PMFS,(xi + xj) / 2.0,yi);
                                }   
                            }

                            else if (yj >= yi + 2.0 * r) {
                                xi -= tmp1;
                                xj -= tmp2;
                                x = dmin(xi,xj) - a;
                                plg_init();
                                plg_arc(x,yi + r,180,270,r,1);
                                plg_arc(x,yj - r,90,180,r,0);
                                plg_line(xi,yi,x,yi);
                                plg_line(x - r,yi + r,x - r,yj - r);
                                plg_line(x,yj,xj - a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,0.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,x - r,(yi + yj) / 2.0);
                            }
                            else if (xi + 2.0 * r <= xj) {
                                yi += rai;
                                yj += raj;
                                y = dmax(yi,yj) + a;
                                plg_init();
                                plg_arc(xi + r,y,90,180,r,1);
                                plg_arc(xj - r,y,0,90,r,0);
                                plg_line(xi,yi,xi,y);
                                plg_line(xi + r,y + r,xj - r,y + r);
                                plg_line(xj,y,xj,yj + a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,270.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,(xi + xj) / 2.0,y + r);

                            }
                            else if (yj + 2.0 * r <= yi) {
                                xi += tmp1;
                                xj += tmp2;
                                x = dmax(xi,xj) + a;
                                plg_init();
                                plg_arc(x,yi - r,0,90,r,1);
                                plg_arc(x,yj + r,270,360,r,0);
                                plg_line(xi,yi,x,yi);
                                plg_line(x + r,yi - r,x + r,yj + r);
                                plg_line(x,yj,xj + a / 2.0,yj);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,180.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,x + r,(yi + yj) / 2.0);
                            }
                            else if (xj + 2.0 * r <= xi) {
                                yi -= rai;
                                yj -= raj;
                                y = dmin(yi,yj) - a;
                                plg_init();
                                plg_arc(xi - r,y,270,360,r,1);
                                plg_arc(xj + r,y,180,270,r,0);
                                plg_line(xi,yi,xi,y);
                                plg_line(xi - r,y - r,xj + r,y - r);
                                plg_line(xj,y,xj,yj - a / 2.0);
                                plg_fin();
                                if (aflag)
                                    plg_arrow(xj,yj,sa,sb,-1,90.0);
                                if (PMICFlg)                           
                                    plg_lab(PMIC,PMFS,(xi + xj) / 2.0, y - r);
                            }
                            else     
                                goto PLGCONT;
                        }
                    }
                }           
                eerr = 0;
PLGCONT:
                if (eerr)
                    printf1("Cannot draw this edge.\n");

                *q = c;
                p = q + 1;

            }
            else
                p++;
        }
    }
    plg_nodes(nn);       /* plot nodes */
    err = 0;

PLGFin:
    for (k = 0; k < nn; ++k) {
        if (AcPtr[k] != NULL) {
            j = strlen(AcPtr[k]) + 1;
            free((char *)AcPtr[k]);
            memrq(-j,sizeof(char));
        }
    }
    p_clean();
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_lab()    Plot label for pl_plg().                                    */
/*              xi,yi in Postscript coordinates.                            */

void plg_lab(double ic,double fs,double xi,double yi)
{
    int l;
    char buf[200];
    double h,w,h2,w2;

    sprintf(buf,"%g",ic);
    l = strlen(buf);
    h2 = fs * PtMM;
    h = 2.0 * h2;
    w2 = (double)l * h2 / 2.0;
    w = 2.0 * w2;

    fprintf(PSFd,"gsave\n%5.2f %5.2f m\n",xi,yi);
    fprintf(PSFd,"%5.2f %5.2f rm\n",-w2,-h2);
    fprintf(PSFd,"0 %5.2f rl\n",h);
    fprintf(PSFd,"%5.2f 0 rl\n",w);
    fprintf(PSFd,"0 %5.2f rl\n",-h);
    fprintf(PSFd,"%5.2f 0 rl\n",-w);
    fprintf(PSFd,"closepath\ngsave\n1 setgray\nfill\ngrestore\n");
    fprintf(PSFd,"gsave\n%5.2f %5.2f rm\n",w2,h2 / 2.0);
    fprintf(PSFd,"/fsiz %5.2f def FT (%g) center show\ngrestore\n",1.5 * fs * PtMM,ic);
    fprintf(PSFd,"grestore\n"); 
}

/*--------------------------------------------------------------------------*/
/*  plg_nodes()  Plot nodes.                                                */

void plg_nodes(int n)
{
    register int k;
    double ra,rb,x,y;

    for (k = 0; k < n; ++k) {   

        ra = AcXF[k] * PtMM;
        rb = AcYF[k] * PtMM;

        fprintf(PSFd,"gsave\n");
        ps_ltyp(AcM[k]);
        ps_lwidth(AcW[k]);
        x = ps_2dx(AcX[k]);     /* translate to PostScript coordinates */
        y = ps_2dy(AcY[k]);

        switch (AcNS[k]) {
            case 1:                                 /* circle */

                fprintf(PSFd,"%5.2f %5.2f %5.2f 0 360 arc\n",x,y,ra);
                break;

            case 2:                                         /* square */
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);       
                fprintf(PSFd,"%5.2f %5.2f rm\n",-ra,-ra);
                fprintf(PSFd,"0 %5.2f rl\n",2.0 * ra);
                fprintf(PSFd,"%5.2f 0 rl\n",2.0 * ra);
                fprintf(PSFd,"0 %5.2f rl\n",-2.0 * ra);
                fprintf(PSFd,"closepath\n");
                break;

            case 3:                                         /* oval */
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);       
                fprintf(PSFd,"%5.2f %5.2f rm\n",-rb,-ra);
                fprintf(PSFd,"%5.2f 0 rl\n",2.0 * rb);
                fprintf(PSFd,"0 %5.2f rm\n",2.0 * ra);
                fprintf(PSFd,"%5.2f 0 rl\n",-2.0 * rb);
                fprintf(PSFd,"%5.2f %5.2f %5.2f 90 270 arc\n",x - rb,y,ra);
                fprintf(PSFd,"%5.2f %5.2f %5.2f 270 90 arc\n",x + rb,y,ra);
                break;

            case 4:                                         /* rectangle */
                fprintf(PSFd,"%5.2f %5.2f m\n",x,y);       
                fprintf(PSFd,"%5.2f %5.2f rm\n",-rb,-ra);
                fprintf(PSFd,"0 %5.2f rl\n",2.0 * ra);
                fprintf(PSFd,"%5.2f 0 rl\n",2.0 * rb);
                fprintf(PSFd,"0 %5.2f rl\n",-2.0 * ra);
                fprintf(PSFd,"closepath\n");
                break;
        }
    
        fprintf(PSFd,"gsave\n");
        ps_fill(AcTmp[k]);
        fprintf(PSFd,"grestore\n");
    
        fprintf(PSFd,"stroke\ngrestore\n");

        rb = AcV[k] * PtMM;                 /* plot label */
        if (AcV[k] > 0.0) {
            if (AcPtr[k] != NULL)
                plot_str(x,y - rb / 2.0,AcPtr[k],1.5 * rb,0,1,0,1,0);
            else if (AcN[k] > 0)
                ps_nlab(AcN[k],AcX[k],AcY[k],AcV[k],0);
        }
    }
}

/*--------------------------------------------------------------------------*/
/*  plg_arrow(x,y,dir,rot)   Plot arrow                                     */

void plg_arrow(double x,double y,double sa,double sb,int dir,double rot)
{
    fprintf(PSFd,"gsave\n%5.2f %5.2f m\n",x,y);          
    if (rot != 0.0)
        fprintf(PSFd,"%5.2f rotate\n",rot);
    fprintf(PSFd,"%5.2f %5.2f scale\n",sa,sb);
    if (dir == 0)
        fprintf(PSFd,"-1 2 rl\n7 -2 rl\n-7 -2 rl\n");
    else if (dir == 1)
        fprintf(PSFd,"1 2 rl\n-7 -2 rl\n7 -2 rl\n");
    else if (dir == -1)
        fprintf(PSFd,"-7 2 rl\n1 -2 rl\n-1 -2 rl\n");
    fprintf(PSFd,"closepath\nfill\nstroke\ngrestore\n");
}   

/*--------------------------------------------------------------------------*/
/*  plg_init()                                                              */

void plg_init(void)
{
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);
}

/*--------------------------------------------------------------------------*/
/*  plg_fin()                                                               */

void plg_fin(void)
{
    fprintf(PSFd,"stroke\ngrestore\n");
}

/*--------------------------------------------------------------------------*/
/*  plg_line()                                                              */

void plg_line(double xa,double ya,double xb,double yb)
{
    fprintf(PSFd,"%5.2f %5.2f m\n%5.2f %5.2f l\n",xa,ya,xb,yb);
}

/*--------------------------------------------------------------------------*/
/*  plg_arc()                                                              */

void plg_arc(double x,double y,int na,int nb,double r,int s)
{
    fprintf(PSFd,"%5.2f %5.2f %5.2f %d %d arc\n",x,y,r,na,nb);
    if (s)
        fprintf(PSFd,"stroke\n");
}




