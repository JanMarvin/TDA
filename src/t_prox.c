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

/*  functions in t_prox.c */

int mproc(void);

/* ------------------------------------------------------------------------ */
/*  mproc       Procrustes rotation. Matrix command: mproc(X,Y,Z)           */  
/*              Performs procrustes rotation for X,Y and creates            */
/*              the new configuration in matrix Z.                          */
/*              Return: 0 if OK, -1 if error.                               */

int mproc(void)
{
    register int i,j,k;
    int err,r,row,col,row1,col1,idx,ivflg;
    register char *p;
    double x,y,ytr,tmp,tmp1,scal;

    err = -1;
    printf2("%s\n",CmdBuf);

    ivflg = 0;

    p = CmdBuf + 6;
    if ((p = n_mexpr(p,1,0,&row,&col,&ivflg,NULL)) == NULL) 
        goto MPROCFin;
    if ((p = m_check1(p)) == NULL)
        goto MPROCFin;
    if ((p = n_mexpr(p,1,1,&row1,&col1,&ivflg,NULL)) == NULL) 
        goto MPROCFin;
    if (row1 != row || col1 != col) {
        mat_err(2);
        goto MPROCFin;
    }
    if ((p = m_check1(p)) == NULL)      /* get output matrix */
        goto MPROCFin;
    if ((p = m_getmat(p,row,col,&idx,CmdBuf,1)) == NULL)
        goto MPROCFin;
            
    /* calculate mean values */

    if (alloc_act(col + 1))   
        goto MPROCFin;
    if (alloc_acw(col + 1))   
        goto MPROCFin;
    if (alloc_actmp(col + 1))   
        goto MPROCFin;

    for (j = 1; j <= col; ++j) {
        x = y = 0.0;
        for (i = 0; i < row; ++i) {
            x += MX[0][i * col + j];
            y += MX[1][i * col + j];
        }
        AcT[j] = x / (double)row;
        AcW[j] = y / (double)row;
    }
                      
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            MX[0][i * col + j] -= AcT[j];
            MX[1][i * col + j] -= AcW[j];
        }
    }
                  

    /* create z and calculate traces */

    if (alloc_acz(col * col + 1))   
        goto MPROCFin;
    if (alloc_acu(col * col + 1))   
        goto MPROCFin;
    if (alloc_acv(col * col + 1))   
        goto MPROCFin;
    if (alloc_acx(col + 1))   
        goto MPROCFin;
              
    ytr = 0.0;
    for (i = 1; i <= col; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (k = 0; k < row; ++k) {
                x = MX[0][k * col + j];
                y = MX[1][k * col + i];
                tmp += x * y;
                if (j == 1)  
                    ytr += y * y;
            }
            AcZ[(i - 1) * col + j] = tmp;
        }
    }
     
    r = svdecomp(col,col,AcZ,AcX,AcU,AcV,3); 
    if (r) {
        p_err(-2,1);
        goto MPROCFin;
    }

    /* calculate scaling factor */
     
    tmp = 0.0; 
    for (i = 1; i <= col; ++i)  
        tmp += AcX[i];

    if (ytr > 0.0 && fabs(ytr - tmp) > EPSI1)  
        scal = tmp / ytr;
    else
        scal = 1.0;

    printf2("Scaling factor: ");
    printf2(PMATFmtS,scal);
    printf2(" (ytr=%lg tmp=%lg)\nRotation matrix\n",ytr,tmp);
    
    /* calculate rotation matrix */
     
    for (i = 0; i < col; ++i) {
        for (j = 0; j < col; ++j) {
            tmp = 0.0;
            for (k = 1; k <= col; ++k)
                tmp += AcU[i * col + k] * AcV[j * col + k];
            AcZ[i * col + j + 1] = tmp;
            printf2(PMATFmtS,tmp);
        }
        printf2("\n");
    }
    printf2("Translation vector\n");
    for (j = 1; j <= col; ++j) {
        tmp = 0.0;
        for (i = 0; i < col; ++i)  
            tmp += AcW[i + 1] * AcZ[i * col + j];
        AcTmp[j] = AcT[j] - tmp * scal;
        printf2(PMATFmtS,AcTmp[j]);
    }
    printf2("\n");

    /* calculate new configuration */
     
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            tmp = 0.0;
            for (k = 0; k < col; ++k)  
                tmp += MX[1][i * col + k + 1] * AcZ[k * col + j];
            tmp *= scal;
            tmp += AcT[j];
            MatVal[idx][i * col + j] = tmp;               
        }
    }
    tmp = 0.0; /* calculate norm */
    for (i = 0; i < row; ++i) {
        for (j = 1; j <= col; ++j) {
            x = MX[0][i * col + j] + AcT[j];
            y = MatVal[idx][i * col + j];
            tmp += (x - y) * (x - y);
        }
    }
    tmp = sqrt(tmp);
    printf2("Norm of difference: ");
    printf2(PMATFmtS,tmp);
    printf2("\n");
    err = 0;

MPROCFin:
    mx_free();         
    p_clean();    
    return(err);
}


