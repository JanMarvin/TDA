/****************************************************************************/
/*  t_prox                                                                  */
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
#include "t_var.h"
#include "t_gdd.h"
#include "t_alloc.h"
#include "t_gdat.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_gio.h"
#include "t_mat.h"
#include "t_matc.h"
#include "t_svd.h"
#include "tda_context.h"

/*  functions in t_prox.c */

int mproc(TDAContext *ctx);

/* ------------------------------------------------------------------------ */
/*  mproc       Procrustes rotation. Matrix command: mproc(X,Y,Z)           */  
/*              Performs procrustes rotation for X,Y and creates            */
/*              the new configuration in matrix Z.                          */
/*              Return: 0 if OK, -1 if error.                               */

int mproc(TDAContext *ctx)
{
    register int i,j,k;
    int err,r,row,col,row1,col1,idx,ivflg;
    register char *p;
    double x,y,ytr,tmp,scal;

    err = -1;
    printf2(ctx, "%s\n",ctx->CmdBuf);

    ivflg = 0;

    p = ctx->CmdBuf + 6;
    if ((p = n_mexpr(ctx, p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MPROCFin;
    if ((p = m_check1(ctx, p)) == NULL)
        goto MPROCFin;
    if ((p = n_mexpr(ctx, p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MPROCFin;
    if (row1 != row || col1 != col) {
        mat_err(ctx, 2);
        goto MPROCFin;
    }
    if ((p = m_check1(ctx, p)) == NULL)      /* get output matrix */
        goto MPROCFin;
    if ((p = m_getmat(ctx, p,row,col,&idx,ctx->CmdBuf,1)) == NULL)
        goto MPROCFin;
            
    /* calculate mean values */

    if (alloc_act(ctx, col + 1))   
        goto MPROCFin;
    if (alloc_acw(ctx, col + 1))   
        goto MPROCFin;
    if (alloc_actmp(ctx, col + 1))   
        goto MPROCFin;

    for (j = 1; j <= col; ++j) {
        x = y = 0.0;
        for (i = 0; i < row; ++i) {
            x += ctx->MX[0][i * col + j];
            y += ctx->MX[1][i * col + j];
        }
        ctx->AcT[j] = x / (double)row;
        ctx->AcW[j] = y / (double)row;
    }
                      
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            ctx->MX[0][i * col + j] -= ctx->AcT[j];
            ctx->MX[1][i * col + j] -= ctx->AcW[j];
        }
    }
                  

    /* create z and calculate traces */

    if (alloc_acz(ctx, col * col + 1))   
        goto MPROCFin;
    if (alloc_acu(ctx, col * col + 1))   
        goto MPROCFin;
    if (alloc_acv(ctx, col * col + 1))   
        goto MPROCFin;
    if (alloc_acx(ctx, col + 1))   
        goto MPROCFin;
              
    ytr = 0.0;
    for (i = 1; i <= col; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (k = 0; k < row; ++k) {
                x = ctx->MX[0][k * col + j];
                y = ctx->MX[1][k * col + i];
                tmp += x * y;
                if (j == 1)  
                    ytr += y * y;
            }
            ctx->AcZ[(i - 1) * col + j] = tmp;
        }
    }
     
    r = svdecomp(ctx, col,col,ctx->AcZ,ctx->AcX,ctx->AcU,ctx->AcV,3); 
    if (r) {
        p_err(ctx, -2,1);
        goto MPROCFin;
    }

    /* calculate scaling factor */
     
    tmp = 0.0; 
    for (i = 1; i <= col; ++i)  
        tmp += ctx->AcX[i];

    if (ytr > 0.0 && fabs(ytr - tmp) > ctx->EPSI1)  
        scal = tmp / ytr;
    else
        scal = 1.0;

    printf2(ctx, "Scaling factor: ");
    rt_printf2_d(ctx, ctx->PMATFmtS,scal);
    printf2(ctx, " (ytr=%lg tmp=%lg)\nRotation matrix\n",ytr,tmp);
    
    /* calculate rotation matrix */
     
    for (i = 0; i < col; ++i) {
        for (j = 0; j < col; ++j) {
            tmp = 0.0;
            for (k = 1; k <= col; ++k)
                tmp += ctx->AcU[i * col + k] * ctx->AcV[j * col + k];
            ctx->AcZ[i * col + j + 1] = tmp;
            rt_printf2_d(ctx, ctx->PMATFmtS,tmp);
        }
        printf2(ctx, "\n");
    }
    printf2(ctx, "Translation vector\n");
    for (j = 1; j <= col; ++j) {
        tmp = 0.0;
        for (i = 0; i < col; ++i)  
            tmp += ctx->AcW[i + 1] * ctx->AcZ[i * col + j];
        ctx->AcTmp[j] = ctx->AcT[j] - tmp * scal;
        rt_printf2_d(ctx, ctx->PMATFmtS,ctx->AcTmp[j]);
    }
    printf2(ctx, "\n");

    /* calculate new configuration */
     
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (k = 0; k < col; ++k)  
                tmp += ctx->MX[1][i * col + k + 1] * ctx->AcZ[k * col + j];
            tmp *= scal;
            tmp += ctx->AcT[j];
            ctx->MatVal[idx][i * col + j] = tmp;               
        }
    }
    tmp = 0.0; /* calculate norm */
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            x = ctx->MX[0][i * col + j] + ctx->AcT[j];
            y = ctx->MatVal[idx][i * col + j];
            tmp += (x - y) * (x - y);
        }
    }
    tmp = sqrt(tmp);
    printf2(ctx, "Norm of difference: ");
    rt_printf2_d(ctx, ctx->PMATFmtS,tmp);
    printf2(ctx, "\n");
    err = 0;

MPROCFin:
    mx_free(ctx);         
    p_clean(ctx);    
    return(err);
}


