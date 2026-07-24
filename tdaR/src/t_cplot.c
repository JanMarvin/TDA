/****************************************************************************/
/*  t_cplot                                                                 */
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

/* open_memstream() is a POSIX extension (used below to fix a real
   contour-fill draw-order bug), not part of strict C99 -- this has to
   be defined before any system header is included (via tda.h, next),
   or glibc hides its own declaration and the compiler silently
   assumes it returns int, corrupting the pointer on a 64-bit build. */
#define _POSIX_C_SOURCE 200809L

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_parm.h"
#include "t_plot.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_eval.h"
#include "t_gf.h"
#include "t_gmin.h"
#include "t_int.h"
#include "t_mat.h"
#include "tda_context.h"

/*  functions in t_cplot.c */

int pl_plotc(TDAContext *ctx);
void cont_plot(TDAContext *ctx, int nlev,double *flev,int opt);
void cont_draw(TDAContext *ctx, double x, double y, int n, int level,int nlev,double *flev, int opt);
int pl_plotcm(TDAContext *ctx);
int pl_plotr(TDAContext *ctx);
void pl_plotr_c(TDAContext *ctx, int im,int row,int col,double lev,double dx,double dy);

/*--------------------------------------------------------------------------*/
/*  global variables                                                        */

#define PCMAX 1.e6  /* max function value in contour plots                  */



/*--------------------------------------------------------------------------*/
/*  pl_plotc()          Contour plot.                                       */
/*                      CmdBuf: plotc(n,lt,lw,nc,x) = function;             */
/*                      Return 0 if OK, -1 if error.                        */

int pl_plotc(TDAContext *ctx)
{
    register int i,j;
    int nn,err,r,pcbita,pcfva;
    double tmp;
           
    pcbita = pcfva = 0;
    err = -1;

    if (check_pcmd(ctx, 2,2))
        return(-1);

    printf1(ctx, "Contour plot.\n");

    if (parm(ctx, ctx->CmdBuf + 5,9,1))    /* get parameters */
        goto PLCFin;

    if (ctx->FNFlg == 0 || ctx->PMNTP < 1) {
        p_err(ctx, -1,1);
        goto PLCFin;
    }
    if (ctx->FNArgN != 2) {
        printf1(ctx, "Error: there must be exactly two function arguments.\n");
        goto PLCFin;
    }
    if (ctx->PMFTYP5) {
        p_err(ctx, -40,1);
        goto PLCFin;
    }
    if (ctx->FVFlg)
        printf1(ctx, "Sum over %d data matrix cases.\n",ctx->NOC);

    ctx->NINTMUsed = -1;

    if (ctx->PMNNFlg == 0 || ctx->PMNN1 < 1 || ctx->PMNN2 < 1)
        ctx->PMNN1 = ctx->PMNN2 = 10;
      
    ctx->PCDX = ctx->UXLen / (double) ctx->PMNN1;
    ctx->PCDY = ctx->UYLen / (double) ctx->PMNN2;
    ctx->PCNX = ctx->PMNN1 + 1;
    ctx->PCNY = ctx->PMNN2 + 1;

    if (!(ctx->PCFV = (float *) calloc((size_t)(ctx->PCNX) * (size_t)(ctx->PCNY) + 1,sizeof(float)))) {
        p_err(ctx, -2,1);
        goto PLCFin;
    }
    memrq(ctx, ctx->PCNX * ctx->PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * ctx->PCNX * ctx->PCNY * ctx->PMNTP / 8 + 1;
    if (!(ctx->PCBitM = (char *)calloc((size_t)(nn),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto PLCFin;
    }
    memrq(ctx, nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= ctx->PCNX; ++i) {
        for (j = 1; j <= ctx->PCNY; ++j) {

            ctx->FNArgVal[0] = ctx->PA1[0] + (double)(i - 1) * ctx->PCDX;
            ctx->FNArgVal[1] = ctx->PA1[1] + (double)(j - 1) * ctx->PCDY;

            r = get_flval(ctx, &tmp,ctx->FNArgN,ctx->FNArgVal,0,0,&tmp,&tmp,&tmp);
            if (r) {      /* r from v_eval1(ctx) */
                printf1(ctx, "Can't evaluate function.\n");
                prn_emsg2(ctx, r);
                goto PLCFin;
            }
            ctx->PCFV[(i - 1) * ctx->PCNY + j] = (float)tmp;
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    {
        /* Set up per-fragment buffering (see cont_draw()'s own
           comment for the full reasoning) only when fills are
           actually requested and there is more than one level to
           possibly misorder; nothing to reorder otherwise. */
        int buffering = ctx->PMGSFlg && ctx->PMNTP > 1;
        int fi;

        if (buffering) {
            ctx->PCFragCap = 16;
            ctx->PCFragCount = 0;
            ctx->PCFragList = (char **)calloc((size_t)(ctx->PCFragCap), sizeof(char *));
            ctx->PCFragAreaList = (double *)calloc((size_t)(ctx->PCFragCap), sizeof(double));
            if (!ctx->PCFragList || !ctx->PCFragAreaList) {
                free(ctx->PCFragList); ctx->PCFragList = NULL;
                free(ctx->PCFragAreaList); ctx->PCFragAreaList = NULL;
                ctx->PCFragCap = 0;
                buffering = 0;
            }
        }

        cont_plot(ctx, ctx->PMNTP,ctx->PMTP,0);

        if (buffering && ctx->PCFragList) {
            /* Sort by descending area (largest painted first), with
               un-filled fragments (area sentinel -1, stroke only, see
               cont_draw()) always last so their own line stays on top
               and visible rather than getting painted over. A plain
               insertion sort: PCFragCount is the number of fragments
               in one contour command, never large enough for this to
               matter. */
            for (fi = 1; fi < ctx->PCFragCount; fi++) {
                double a = ctx->PCFragAreaList[fi];
                char *m = ctx->PCFragList[fi];
                int fj = fi - 1;
                while (fj >= 0 && ctx->PCFragAreaList[fj] < a) {
                    ctx->PCFragAreaList[fj + 1] = ctx->PCFragAreaList[fj];
                    ctx->PCFragList[fj + 1] = ctx->PCFragList[fj];
                    fj--;
                }
                ctx->PCFragAreaList[fj + 1] = a;
                ctx->PCFragList[fj + 1] = m;
            }
            for (fi = 0; fi < ctx->PCFragCount; fi++) {
                if (ctx->PCFragList[fi]) {
                    fputs(ctx->PCFragList[fi], ctx->PSFd);
                    free(ctx->PCFragList[fi]);
                }
            }
            free(ctx->PCFragList);
            free(ctx->PCFragAreaList);
            ctx->PCFragList = NULL;
            ctx->PCFragAreaList = NULL;
            ctx->PCFragCount = 0;
            ctx->PCFragCap = 0;
        }
    }
    err = 0;

PLCFin:
    if (pcfva) {
        free((char *)ctx->PCFV);
        memrq(ctx, -ctx->PCNX * ctx->PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)ctx->PCBitM);
        memrq(ctx, -nn,sizeof(char));
    }
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  cont_plot   Contour plot routine.                                       */
/*                                                                          */
/*    The code is adopted from                                              */
/*    Snyder, ACM algorithm 531, Trans. in Math. Software 4(3), 290 - 294   */

void cont_plot(TDAContext *ctx, int nlev,double *flev,int opt)
{
    register int k,l,ii,jj,q;
    int ix,icv,ni,ks,jump;
    int idir,nxidir,ibkey,iflag,icur,jcur,iedge;
    int ij[3],l1[5],l2[5],i1[3],i2[3],i3[7];
    double tmp,cval,dmax,z1,z2,zz,xint[5],xy[3];

    ks = 0;
    l1[1] = ctx->PCNX;
    l1[2] = ctx->PCNY;
    l1[3] = l1[4] = -1;
    i1[1] = i2[1] =  1;
    i1[2] = 0;
    i2[2] = -1;
    i3[1] = i3[4] = i3[5] = 1;
    i3[2] = i3[3] = i3[6] = 0;

    dmax = PCMAX;

    /* set current pen position, default corresponds to z(1,1). */
    
    xy[1] = xy[2] = 1.0;
    
    icur = (int)xy[1];
    if (icur > ctx->PCNX) icur = ctx->PCNX;
    if (icur < 1)    icur = 1;
    jcur = (int)xy[2];
    if (jcur > ctx->PCNY) jcur = ctx->PCNY;
    if (jcur < 1)    jcur = 1;
    ibkey = 0;

PCL10:
    ij[1] = icur;
    ij[2] = jcur;
PCL20:
    l2[1] = ij[1];
    l2[3] = -ij[1];
    l2[2] = ij[2];
    l2[4] = -ij[2];
    idir = 0;

PCL30:
    nxidir = idir + 1;
    k = nxidir;
    if (nxidir > 3)
        nxidir = 0;
PCL40:
    if (ij[1] < 0) ij[1] = -ij[1];
    if (ij[2] < 0) ij[2] = -ij[2];

    if ((double)(ctx->PCFV[(ij[1] - 1) * ctx->PCNY + ij[2]]) > dmax)
        goto PCL140;

    l = 1;
    while (ij[l] >= l1[l]) {
        if (++l > 2)
            goto PCL140;
    }
    while (1) {
        ii = ij[1] + i1[l];
        jj = ij[2] + i1[3 - l];

        while ((double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]) > dmax) {
            if (++l > 2)
                goto PCL140;
            else {
                while (ij[l] >= l1[l]) {
                    if (++l > 2)
                        goto PCL140;
                }
                ii = ij[1] + i1[l];
                jj = ij[2] + i1[3 - l];
            }
        }
        jump = 1;

PCL60:
        ix = 1;
        if (ij[3 - l] != 1) {
            ii = ij[1] - i1[3 - l];
            jj = ij[2] - i1[l];
    
            if ((double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]) <= dmax) {
                ii = ij[1] + i2[l];
                jj = ij[2] + i2[3 - l];

                if ((double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]) < dmax)
                    ix = 0;
            }
            if (ij[3 - l] >= l1[3 - l])
                goto PCL90;
        }
        ii = ij[1] + i1[3 - l];
        jj = ij[2] + i1[l];

        if ((double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]) <= dmax &&
                                     (double)(ctx->PCFV[ij[1] * ctx->PCNY + ij[2] + 1]) < dmax) {
            if (jump)
                goto PCL100;
            else {
                if (ix != 0)
                    iflag = 4;
                iedge = ks + 2;
                if (iedge > 4)
                    iedge -= 4;
                xint[iedge] = xint[ks];
                goto PCL200;
            }
        }
PCL90:
        ix += 2;
        if (!jump) {
            if (ix != 0)
                iflag = 4;
            iedge = ks + 2;
            if (iedge > 4)
                iedge -= 4;
            xint[iedge] = xint[ks];
            goto PCL200;
        }

PCL100:
        if (ix != 3 && (ix + ibkey) != 0) {

            ii = ij[1] + i1[l];
            jj = ij[2] + i1[3 - l];
            z1 = (double)(ctx->PCFV[(ij[1] - 1) * ctx->PCNY + ij[2]]);
            z2 = (double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]);

            for (icv = 1; icv <= nlev; ++icv) {

                q = 2 * (ctx->PCNX * (ctx->PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
                if (!(*(ctx->PCBitM + q / 8) & ctx->PBMsk[q % 8])) {
                    tmp = z1;
                    if (tmp > z2)
                        tmp = z2;
                    if (flev[icv - 1] > tmp) {
                        tmp = z1;
                        if (tmp < z2)
                            tmp = z2;
                        if (flev[icv - 1] <= tmp) {
                            iedge = l;
                            cval = flev[icv - 1];
                            if (ix != 1)
                                iedge += 2;
                            iflag = 2 + ibkey;
                            xint[iedge] = (cval - z1) / (z2 - z1);
                            goto PCL200;
                        }
                    }
                    q = 2 * (ctx->PCNX * (ctx->PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
                    *(ctx->PCBitM + q / 8) |= ctx->PBMsk[q % 8];
                }
            }
        }
        if (++l > 2)  
           break;

        while (ij[l] >= l1[l]) {
            if (++l > 2)
                goto PCL140;
        }
    }

PCL140:
    l = (idir % 2) + 1;
    if (l1[k] < 0)
        ij[l] = -ij[l];

    while (1) {
        if (ij[l] < l1[k]) {
            ij[l] += 1;
            if (ij[l] <= l2[k])
                goto PCL40;

            l2[k] = ij[l];
            idir = nxidir;
            goto PCL30;
        }
        if (idir != nxidir) {
            nxidir++;
            ij[l] = l1[k];
            k = nxidir;
            l = 3 - l;
            ij[l] = l2[k];
            if (nxidir > 3)
                nxidir = 0;
        }
        else if (ibkey != 0)  
            return;
        else 
            break;
    }
    ibkey = 1;
    goto PCL10;

PCL200:
    while (1) {
    
        xy[l] = (double)ij[l] + xint[iedge];
        xy[3 - l] = (double) ij[3 - l];
        q = 2 * (ctx->PCNX * (ctx->PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
        *(ctx->PCBitM + q / 8) |= ctx->PBMsk[q % 8];

        cont_draw(ctx, xy[1],xy[2],iflag,icv,nlev,flev,opt);
    
        if (iflag >= 4) {
            icur = ij[1];
            jcur = ij[2];
            goto PCL20;
        }
        ni = 1;
        if (iedge >= 3) {
            ij[1] -= i3[iedge];
            ij[2] -= i3[iedge + 2];
        }
        for (k = 1; k <= 4; ++k) {
    
            if (k != iedge) {
                ii = ij[1] + i3[k];
                jj = ij[2] + i3[k + 1];
                z1 = (double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]);
                ii = ij[1] + i3[k + 1];
                jj = ij[2] + i3[k + 2];
                z2 = (double)(ctx->PCFV[(ii - 1) * ctx->PCNY + jj]);
    
                tmp = z1;
                if (tmp > z2)
                    tmp = z2;
                if (cval > tmp) {
                    tmp = z1;
                    if (tmp < z2)
                        tmp = z2;
                    if (cval <= tmp) {
                        if (k == 1 || k == 4) {
                            zz = z1;
                            z1 = z2;
                            z2 = zz;
                        }
                        xint[k] = (cval - z1) / (z2 - z1);
                        ni++;
                        ks = k;
                    }
                }
            }
        }
        if (ni != 2) {
            ks = 5 - iedge;
            if (xint[3] >= xint[1]) {
                ks = 3 - iedge;
                if (ks <= 0)
                    ks += 4;
            }
        }
        l = ks;
        iflag = 1;
        jump = 0;
        if (ks >= 3) {
            ij[1] += i3[ks];
            ij[2] += i3[ks + 2];
            l = ks - 2;
        }
        q = 2 * (ctx->PCNX * (ctx->PCNY * (icv - 1) + ij[2] - 1) + ij[1] - 1) + l - 1;
        if (!(*(ctx->PCBitM + q / 8) & ctx->PBMsk[q % 8]))  
            goto PCL60;
    
        iflag = 5;
        iedge = ks + 2;
        if (iedge > 4)
            iedge -= 4;
    
        xint[iedge] = xint[ks];
    }
}   

/* ------------------------------------------------------------------------ */
/*  cont_draw   Draw function to be used with cont_plot.                    */

void cont_draw(TDAContext *ctx, double x, double y, int n, int level,int nlev,double *flev, int opt)
{
    register int i,j,fflag;
    double g,xx,yy,px,py,cx,cy;
    static double eps = 1.e-3;

    fflag = 0;
    cx = cy = -1.0;

    i = (int)x;
    j = (int)y;
    xx = ctx->PA1[0] + (double)(i - 1) * ctx->PCDX;
    yy = ctx->PA1[1] + (double)(j - 1) * ctx->PCDY;
    xx += (x - (double)i) * ctx->PCDX;
    yy += (y - (double)j) * ctx->PCDY;

    if (opt) {
        xx += ctx->PCDX / 2.0;
        yy += ctx->PCDY / 2.0;
    }


    px = ctx->PtMM * (xx - ctx->PA1[0]) * ctx->PXLen / ctx->UXLen;
    py = ctx->PtMM * (yy - ctx->PA1[1]) * ctx->PYLen / ctx->UYLen;

    switch (n) {
        case  1:  if (ctx->PCFragList && ctx->PCFragBuf) {
                      ctx->PCFragArea += ctx->PCFragXPrev * py - px * ctx->PCFragYPrev;
                      ctx->PCFragXPrev = px;
                      ctx->PCFragYPrev = py;
                  }
                  fprintf(ctx->PSFd,"%4.2f %4.2f l\n",px,py);
                  break;
        case  2:
        case  3:  /* Contour fill draw-order fix: buffer this fragment
                      (one gsave..grestore block) in its own memstream
                      instead of writing straight to the real
                      PostScript file, so pl_plotc() can flush every
                      fragment traced this call afterward sorted by
                      measured filled area (largest first) rather than
                      in whatever order the underlying grid-tracing
                      algorithm happened to discover them.

                      That order is not the same thing as ascending
                      level number, and assuming otherwise was a real
                      mistake caught before landing: checked directly
                      against a real run's own polygon data, a fringe
                      contour fragment near a plot corner can end up
                      geometrically larger than a *lower*-level
                      fragment nearby (the corner-closing straight line
                      below has to sweep back through more of the plot
                      to close a curve that exits further from the
                      corner, inflating that fragment's own closed
                      area independently of its level), so sorting by
                      level number alone got the paint order backwards
                      in exactly the cases it was meant to fix. Sorting
                      by each fragment's own actual, measured area
                      (via the shoelace formula, accumulated below as
                      each point is written) is correct regardless of
                      why a fragment ended up whatever size it is. */
                  if (ctx->PCFragList) {
                      ctx->PCFragMem = NULL;
                      ctx->PCFragMemSize = 0;
#if defined(_WIN32) || defined(__CYGWIN__)
                      ctx->PCFragBuf = tmpfile();
#else
                      ctx->PCFragBuf = open_memstream(&ctx->PCFragMem, &ctx->PCFragMemSize);
#endif
                      if (ctx->PCFragBuf) {
                          ctx->PCFragRealPSFd = ctx->PSFd;
                          ctx->PSFd = ctx->PCFragBuf;
                      }
                      ctx->PCFragArea = 0.0;
                      ctx->PCFragX0 = ctx->PCFragXPrev = px;
                      ctx->PCFragY0 = ctx->PCFragYPrev = py;
                  }

                  fprintf(ctx->PSFd,"gsave %% level(%d) = %g\n",level,flev[level - 1]);

                  ps_lwidth(ctx, ctx->PMLW);
                  ps_ltyp(ctx, ctx->PMLT);

                  fprintf(ctx->PSFd,"%4.2f %4.2f m\n",px,py);
                  ctx->s_cont_draw_xf = px;
                  ctx->s_cont_draw_yf = py; 
                  break;
        case  4:
        case  5:  if (ctx->PCFragList && ctx->PCFragBuf) {
                      ctx->PCFragArea += ctx->PCFragXPrev * py - px * ctx->PCFragYPrev;
                      ctx->PCFragXPrev = px;
                      ctx->PCFragYPrev = py;
                  }
                  fprintf(ctx->PSFd,"%4.2f %4.2f l\n",px,py);

                  /* fill only curves which can be closed */

                  if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {

                      fflag = 0;
                      if (fabs(px - ctx->s_cont_draw_xf) <= eps || fabs(py - ctx->s_cont_draw_yf) <= eps)  
                          fflag = 1;
                      else if ((fabs(px) <= eps && fabs(ctx->s_cont_draw_yf) <= eps) ||
                          (fabs(py) <= eps && fabs(ctx->s_cont_draw_xf) <= eps)) {
                          fprintf(ctx->PSFd,"0 0 l %% aaa\n");
                          cx = 0.0; cy = 0.0;
                          fflag = 1;
                      }
                      else if ((fabs(px) <= eps && fabs(ctx->s_cont_draw_yf - ctx->PSYLen) <= eps) ||
                               (fabs(py - ctx->PSYLen) <= eps && fabs(ctx->s_cont_draw_xf) <= eps)) {
                          fprintf(ctx->PSFd,"0 %5.2f l %% bbb\n",ctx->PSYLen);
                          cx = 0.0; cy = ctx->PSYLen;
                          fflag = 1;
                      }
                      else if ((fabs(px - ctx->PSXLen) <= eps && fabs(ctx->s_cont_draw_yf) <= eps) ||
                               (fabs(py) <= eps && fabs(ctx->s_cont_draw_xf - ctx->PSXLen) <= eps)) {
                          fprintf(ctx->PSFd,"%5.2f 0 l %% ccc\n",ctx->PSXLen);
                          cx = ctx->PSXLen; cy = 0.0;
                          fflag = 1;
                      }
                      else if ((fabs(px - ctx->PSXLen) <= eps && fabs(ctx->s_cont_draw_yf - ctx->PSYLen) <= eps) ||
                             (fabs(py - ctx->PSYLen) <= eps && fabs(ctx->s_cont_draw_xf - ctx->PSYLen) <= eps)) {
                          fprintf(ctx->PSFd,"%5.2f %5.2f l %% ddd\n",ctx->PSXLen,ctx->PSYLen);
                          cx = ctx->PSXLen; cy = ctx->PSYLen;
                          fflag = 1;
                      }
                      if (fflag) {
                          if (ctx->PCFragList && ctx->PCFragBuf) {
                              if (cx >= 0.0) {
                                  ctx->PCFragArea += ctx->PCFragXPrev * cy - cx * ctx->PCFragYPrev;
                                  ctx->PCFragXPrev = cx;
                                  ctx->PCFragYPrev = cy;
                              }
                              ctx->PCFragArea += ctx->PCFragXPrev * ctx->PCFragY0 - ctx->PCFragX0 * ctx->PCFragYPrev;
                          }

                          if (ctx->PMLW > 0.0)
                              fprintf(ctx->PSFd,"gsave\n");

                          if (ctx->PMGSFlg == 1)
                              g = ctx->PMGS * (1.0 + (double)(1 - level) / (double)nlev);
                          else if (nlev > 1)  
                              g = ctx->PMGS - (double)(level - 1) * (ctx->PMGS - ctx->PMGS1) / (double)(nlev - 1);
                          else
                              g = ctx->PMGS;

                          ps_fill(ctx, g);
                          if (ctx->PMLW > 0.0)
                              fprintf(ctx->PSFd,"grestore\n");
                      }
                  }
                  fprintf(ctx->PSFd,"stroke\n");
                  fprintf(ctx->PSFd,"grestore\n");

                  if (ctx->PCFragList && ctx->PCFragBuf) {
                      double area;
#if defined(_WIN32) || defined(__CYGWIN__)
                      if (ctx->PCFragBuf) {
                          fseek(ctx->PCFragBuf, 0, SEEK_END);
                          long size = ftell(ctx->PCFragBuf);
                          rewind(ctx->PCFragBuf);

                          ctx->PCFragMem = malloc((size_t)size + 1);
                          ctx->PCFragMemSize = (size_t)fread(ctx->PCFragMem, 1, (size_t)size, ctx->PCFragBuf);
                          ctx->PCFragMem[ctx->PCFragMemSize] = '\0';
                          fclose(ctx->PCFragBuf);
                      }
#else
                      fclose(ctx->PCFragBuf);
#endif
                      area = fflag ? fabs(ctx->PCFragArea) / 2.0 : -1.0;
                      if (ctx->PCFragCount >= ctx->PCFragCap) {
                          int newcap = ctx->PCFragCap > 0 ? ctx->PCFragCap * 2 : 16;
                          char **nl = (char **)realloc(ctx->PCFragList, (size_t)newcap * sizeof(char *));
                          double *na = (double *)realloc(ctx->PCFragAreaList, (size_t)newcap * sizeof(double));
                          if (nl) ctx->PCFragList = nl;
                          if (na) ctx->PCFragAreaList = na;
                          if (nl && na) ctx->PCFragCap = newcap;
                      }
                      if (ctx->PCFragCount < ctx->PCFragCap) {
                          ctx->PCFragList[ctx->PCFragCount] = ctx->PCFragMem;
                          ctx->PCFragAreaList[ctx->PCFragCount] = area;
                          ctx->PCFragCount++;
                      } else {
                          free(ctx->PCFragMem);
                      }
                      ctx->PSFd = ctx->PCFragRealPSFd;
                      ctx->PCFragBuf = NULL;
                      ctx->PCFragMem = NULL;
                  }
                  break;
        default:  fprintf(ctx->PSFd,"%% cont_draw n=%d level=%d\n",n,level);
                  break;
    }
}


/*--------------------------------------------------------------------------*/
/*  pl_plotcm   Contour plot.                                               */
/*                                                                          */
/*              plotcm(                                                     */
/*                  opt=...,    option, def. 1                              */
/*                              1 use all matrix elements                   */
/*                              2 use only non-negative matrix elements     */
/*                  dopt=...,   option for matrix rows, def. 1              */
/*                              1 use matrix rows in standard order         */
/*                              2 use matrix rows in reverse order          */
/*                  lt=...,     line type, def. 1                           */
/*                  lw=...,     line width, def.                            */
/*                  x=...,      level                                       */
/*                  gs=...,     grey-scale option                           */
/*              ) = matrix_name;                                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int pl_plotcm(TDAContext *ctx)
{
    register int i,j,k;
    int nn,err,pcbita,pcfva,im;
    double tmp;
           
    pcbita = pcfva = 0;
    err = -1;

    if (check_pcmd(ctx, 2,2))
        return(-1);

    printf1(ctx, "Contour plot.\n");
    if (parm(ctx, ctx->CmdBuf + 6,14,1))    /* get parameters */
        goto PLCMFin;

    if (ctx->PMNVTyp != 1 || ctx->PMNV != 1) {
        printf1(ctx, "Error: need exactly one matrix on right-hand side.\n");
        goto PLCMFin;
    }
    if (ctx->PMOPT != 2)
        ctx->PMOPT = 1;
    if (ctx->PMDOPT != 2)
        ctx->PMDOPT = 1;

    im = ctx->PMVIdx[0];
    ctx->PCNY = ctx->MatRow[im];
    ctx->PCNX = ctx->MatCol[im];

    ctx->PCDX = ctx->UXLen / (double)ctx->PCNX;
    ctx->PCDY = ctx->UYLen / (double)ctx->PCNY;

    if (!(ctx->PCFV = (float *) calloc((size_t)(ctx->PCNX) * (size_t)(ctx->PCNY) + 1,sizeof(float)))) {
        p_err(ctx, -2,1);
        goto PLCMFin;
    }
    memrq(ctx, ctx->PCNX * ctx->PCNY + 1,sizeof(float));
    pcfva = 1;

    nn = 2 * ctx->PCNX * ctx->PCNY * ctx->PMNTP / 8 + 1;
    if (!(ctx->PCBitM = (char *)calloc((size_t)(nn),sizeof(char)))) {
        p_err(ctx, -2,1);
        goto PLCMFin;
    }
    memrq(ctx, nn,sizeof(char));
    pcbita = 1;

    /* get function values in PCFV */
      
    for (i = 1; i <= ctx->PCNX; ++i) {
        for (j = 0; j < ctx->PCNY; ++j) {
            tmp = (double)((float)ctx->MatVal[im][j * ctx->PCNX + i]);       
            if (ctx->PMOPT == 2 && tmp < 0.0)
                tmp = PCMAX - 1.0;  

            if (ctx->PMDOPT == 2)
                k = j + 1;
            else
                k = ctx->PCNY - j;

            ctx->PCFV[(i - 1) * ctx->PCNY + k] = (float)tmp;
        }
    }
    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    cont_plot(ctx, ctx->PMNTP,ctx->PMTP,1);
    err = 0;

PLCMFin:
    if (pcfva) {
        free((char *)ctx->PCFV);
        memrq(ctx, -ctx->PCNX * ctx->PCNY - 1,sizeof(float));
    }
    if (pcbita) {
        free((char *)ctx->PCBitM);
        memrq(ctx, -nn,sizeof(char));
    }
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotr    Plot grey-scaled relief.                                    */
/*                                                                          */
/*              plotr(                                                      */
/*                  dopt=...,   option for matrix rows, def. 1              */
/*                              1 use matrix rows in standard order         */
/*                              2 use matrix rows in reverse order          */
/*                  gs = min,max, def. 0,1                                  */
/*              ) = matrix_name;                                            */

int pl_plotr(TDAContext *ctx)
{
    register int i,j,k;
    int err,im,row,col,nm;
    double a,b,d,dx,dy,x,y,tmp,vmin,vmax,g;
           
    err = -1;
    if (check_pcmd(ctx, 1,2))
        return(-1);

    if (parm(ctx, ctx->CmdBuf + 5,14,1))    /* get parameters */
        goto PLRFin;

    if (ctx->PMNVTyp != 1 || ctx->PMNV != 1) {
        printf1(ctx, "Error: need exactly one matrix on right-hand side.\n");
        goto PLRFin;
    }
    if (ctx->PMDOPT != 2)
        ctx->PMDOPT = 1;

    a = 0.0;
    b = 1.0;
    if (ctx->PMGSFlg == 2) {
        a = ctx->PMGS;
        b = ctx->PMGS1;
    }
    im = ctx->PMVIdx[0];
    row = ctx->MatRow[im];
    col = ctx->MatCol[im];

    dx = ctx->UXLen / (double)col;
    dy = ctx->UYLen / (double)row;

    vmin = ctx->DBLMAX;
    vmax = -1.0;      
    nm = 0;
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = ctx->MatVal[im][i * col + j];
            if (tmp < 0.0)
                nm++;
            else {
                vmin = dmin(ctx, vmin,tmp);
                vmax = dmax(ctx, vmax,tmp);
            }
        }
    }
    printf1(ctx, "Number of cells: %d x %d = %d\n",row,col,row * col);
    printf1(ctx, "Number of empty cells: %d\n",nm);
    if (nm == row * col) {
        err = 0;
        goto PLRFin;
    }
    tda_out("Minimum: %g\n",vmin);
    tda_out("Maximum: %g\n",vmax);
    d = vmax - vmin;
    if (d < ctx->EPSI1) {
        printf1(ctx, "Error: values almost identical.\n");
        goto PLRFin;
    }

    fprintf(ctx->PSFd,"\n%%#%d: %s\n",++ctx->PSONUM,ctx->CmdBuf);

    fprintf(ctx->PSFd,"gsave\n");
    if (ctx->PMNC != 1)
        set_clip(ctx);
    ps_lwidth(ctx, ctx->PMLW);
    ps_ltyp(ctx, ctx->PMLT);
   
    for (i = row - 1; i >= 0; --i) {

        if (ctx->PMDOPT == 2)
            k = row - i - 1;
        else
            k = i;

        for (j = 1; j <= col; ++j) {

            tmp = ctx->MatVal[im][k * col + j];
            if (tmp < 0.0)
                continue;

            g = (vmax - tmp) / d;
            g = a + g * (b - a);
            x = ctx->PA1[0] + dx * (double)(j - 1);
            y = ctx->PA1[1] + dy * (double)(row - i - 1);

            fprintf(ctx->PSFd,"gsave %4.2f setgray\n",g);    
            ps_2dplot(ctx, x,y,0);    
            ps_2dplot(ctx, x,y + dy,1);    
            ps_2dplot(ctx, x + dx,y + dy,1);    
            ps_2dplot(ctx, x + dx,y,1);    
            ps_2dplot(ctx, x,y,1);    
            fprintf(ctx->PSFd,"fill\ngrestore\n");
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
/**
    pl_plotr_c(ctx, im,row,col,30.0,dx,dy);
    pl_plotr_c(ctx, im,row,col,50.0,dx,dy);
    pl_plotr_c(ctx, im,row,col,80.0,dx,dy);
    pl_plotr_c(ctx, im,row,col,100.0,dx,dy);
    pl_plotr_c(ctx, im,row,col,120.0,dx,dy);
**/ 

    err = 0;

PLRFin:
    p_clean(ctx);
    return(err);
}

/*--------------------------------------------------------------------------*/
/*  pl_plotr_c(                                                             */

void pl_plotr_c(TDAContext *ctx, int im,int row,int col,double lev,double dx,double dy)
{
    register int i,j;
    double x,y,a,b;                                

    fprintf(ctx->PSFd,"gsave\n");    
    for (i = row - 1; i >= 0; --i) {
        for (j = 1; j <= col; ++j) {
            a = ctx->MatVal[im][i * col + j];
            x = ctx->PA1[0] + dx * (double)(j - 1);
            y = ctx->PA1[1] + dy * (double)(row - i - 1);

            if (j < col) {
                b = ctx->MatVal[im][i * col + j + 1];
                if ((a < lev && b >= lev) || (b < lev && a >= lev)) {
                    ps_2dplot(ctx, x + dx,y,0);    
                    ps_2dplot(ctx, x + dx,y + dy,1);    
                    fprintf(ctx->PSFd,"stroke\n");
                }
            }
            if (i > 0) {
                b = ctx->MatVal[im][(i - 1) * col + j];
                if ((a < lev && b >= lev) || (b < lev && a >= lev)) {
                    ps_2dplot(ctx, x,y + dy,0);    
                    ps_2dplot(ctx, x + dx,y + dy,1);    
                    fprintf(ctx->PSFd,"stroke\n");
                }
            }
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
}
