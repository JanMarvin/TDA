/****************************************************************************/
/*  t_pdat                                                                  */
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
#include "t_gdat.h"  
#include "t_var.h"  
#include "t_gf.h"  
#include "t_alloc.h"  
#include "t_sort.h"  
#include "t_sd.h"  

/*  functions in t_pdat.c */

int pdata(void);
void dtda(char *fname,int noc,int nv,short *vidx,int nq);
void dspss(char *fname,int nv,short *vidx,int nq);
int pdatr(void);
int pdatd(void);
void pdatdr(int n,double *x);
double pdatd1(int opt,int n,double *x,double *y,int nw,double *w);


/* ------------------------------------------------------------------------ */
/*  pdata()     pdata command.                                              */
/*                                                                          */
/*          pdata (                                                         */
/*              keep=varlist,                                               */
/*              drop=varlist,                                               */
/*              sort=varlist,                                               */
/*              noc=...,        only first noc cases                        */
/*              nn = n1,n2      range of cases                              */
/*              transp          transpose data                              */
/*              nc = ...        line feed after nc variables                */
/*              nq = ...        line feed after nq records, def. 1          */
/*              ap=1,           append                                      */
/*              prn=...,        output options, def. 0                      */
/*                              1 triangle                                  */
/*                              2 square matrix                             */
/*                              3 write data as single column               */
/*              dtda=...,       tda description file                        */
/*              dspss=...,      SPSS description file                       */
/*              sepc=' ',       separation character, def ' '               */
/*                              also possible: sepc=none                    */
/*              l0=...,         1 if leading zeros, def. 0                  */
/*                                                                          */
/*              sd=...,         if 1 write s-data records, def. 0           */
/*                              (only with prn=0)                           */
/*              sdfmt=...,      print format for s-data records, def. 12.4  */
/*          ) = fname;          output file                                 */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */
 
int pdata(void)
{
    register int i,j,k,ii;
    int err,nrec,n,nn,nv,na,nb,sflag,lcnt,w,d,r,sdrec;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    if (parm(CmdBuf + 5,1,1)) {
        p_clean();
        return(-1);
    }
    if (PMSD) {
        if (check_sd(0)) {
            printf1("Error: no spatial data defined.\n");
            goto PDATFin;
        }
        if (PMPRNO != 0) {
            printf1("Error: sd option requires prn=0.\n");
            goto PDATFin;
        }
    }
    if (PMKeep && PMDrop) {
        p_err(-16,1);
        goto PDATFin;
    }
    if (PMNQ < 1)
        PMNQ = 1;

    /* get list of variables in AcNS[] */

    if (alloc_acns(imax(NVAR,PMNV)))
        goto PDATFin;

    nv = 0;
    if (PMKeep) {
        for (k = 0; k < PMNV; ++k)           
            AcNS[nv++] = PMVIdx[k];              
    }
    else {
        j = VIFirst;
        while (j >= 0) {
            i = 1;
            if (PMDrop) {
                for (k = 0; k < PMNV; ++k) {
                    if (PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                AcNS[nv++] = j;
            j = VNxt[j];
        }
    }
       
    if (nv == 0) {
        printf1("No variables selected.\n");
        goto PDATFin;
    }
    if (PMPRNO == 1 || PMPRNO == 2) {
        if (nv != 1) {
            printf1("Error: prn options can only be used with a single variable.\n");
            goto PDATFin;
        }
    }

    sflag = 0;
    if (PM1NV > 0) {        /* sort */
        err = vsort(PM1NV,PM1VIdx,1,0,1);
        if (err)
            goto PDATFin;
        sflag = 1;
    }
    na = 0; 
    nb = NOC;
    if (PMNNFlg) {
        na = PMNN1 - 1;
        nb = PMNN2;
    }
    else if (PMNOCFlg)
        nb = PMNOC;
    if (na < 0)
        na = 0;
    if (nb > NOC)
        nb = NOC;

    nn = 0;
    if (PMPRNO == 2) {
        nn = (int)sqrt((double)(nb - na));
        if (nn * nn != nb - na) {
            printf1("Error: number of cases not compatible with prn=2 option.\n");
            goto PDATFin;
        }
    }
    if (PML0 != 0)
        PML0 = -1;

    if (PML0 != 0 || XSEPC != SEPC) {        /* change print format */
        for (k = 0; k < nv; ++k) {
            j = AcNS[k];
            if (VTyp[j] == 1)
                continue;
            w = (int)VPFmt1[j];             
            d = (int)VPFmt2[j];             
            makefmt(&w,&d,VPFmtS[j],0,XSEPC,PML0);
        }
    }

    if (SILENTFlg < 2)
        printfe("Writing: %s\n",PMFdName);

    sdrec = nrec = 0;

    n = 0;
    if (PMTransp == 0) {
        lcnt = 0;
        for (i = na; i < nb; ++i) {
            ii = i;
            if (sflag)
                ii = VSORTPtr[i];

            if (PMPRNO == 1 || PMPRNO == 2) {
                j = AcNS[0];       
                if (VTyp[j] == 1)
                    continue;
                tmp = get_data(j,ii);
                fprintf(PMFd,VPFmtS[j],tmp);
                if (PMPRNO == 1) {
                    if (++n > nn) {
                        fprintf(PMFd,"\n");
                        nrec++;
                        prn_message(nrec,0,1);
                        n = 0;
                        nn++;
                    }
                }
                else if (++n >= nn) {
                    fprintf(PMFd,"\n");
                    nrec++;
                    prn_message(nrec,0,1);
                    n = 0;
                }
            }
            else {              /* standard printing */
                n = 0;
                for (k = 0; k < nv; ++k) {
                    j = AcNS[k];
                    if (VTyp[j] != 1) {
                        tmp = get_data(j,ii);                  
                        fprintf(PMFd,VPFmtS[j],tmp);
                    }
                    else {          /* string variables */
                        get_str(SVBuf,j,ii);
                        fprintf(PMFd,"%s",SVBuf);
                        if (XSEPC)
                            fprintf(PMFd,"%c",XSEPC);
                    }
                    if (k == nv - 1)
                        break;

                    if (PMPRNO == 3 || (PMNC > 0 && ++n >= PMNC)) {
                        fprintf(PMFd,"\n");
                        nrec++;
                        prn_message(nrec,0,1);
                        n = 0;
                    }
                }
                if (++lcnt >= PMNQ) {
                    fprintf(PMFd,"\n");
                    nrec++;
                    prn_message(nrec,0,1);
                    lcnt = 0;
                }
                if (PMSD) {
                    if ((r = sd_getdata(ii,0,0,0)) < 1) {
                        printf1("Error: cannot read s-data records.\n");
                        goto PDATFin;
                    }
                    for (j = 0; j < r; ++j) {
                        fprintf(PMFd,PMSDFmtS,SDVarX[j]);
                        fprintf(PMFd,PMSDFmtS,SDVarY[j]);
                        fprintf(PMFd,"\n");
                        sdrec++;
                    }
                }
            }
        }
        if (PMNQ > 1 && lcnt > 0 && lcnt < PMNQ) {
            fprintf(PMFd,"\n");
            nrec++;
            prn_message(nrec,0,1);
        }
    }
    else {
        for (k = 0; k < nv; ++k) {
            j = AcNS[k];
            if (VTyp[j] == 1)               /* string variables */
                continue;

            for (i = na; i < nb; ++i) {
                ii = i;
                if (sflag)
                    ii = VSORTPtr[i];
                tmp = get_data(j,ii);
                fprintf(PMFd,VPFmtS[j],tmp);
            }
            fprintf(PMFd,"\n");
            nrec++;
            prn_message(nrec,0,1);
        }
        nv = nb - na;
    }
    prn_message(nrec,1,1);

    if (PMPRNO == 3)
        printf1("%d records written to: %s\n",nrec,PMFdName);
    else {
        printf1("%d records with %d variables written to: %s\n",nrec,nv,PMFdName);
        if (PMSD)
            printf1("Added %d s-data records.\n",sdrec);
    }
    if (PMTransp == 0 && PMPRNO == 0) {

        if (PMTDAFDef)
            dtda(PMFdName,nrec,nv,AcNS,PMNQ);

        if (PMSPSSFDef)
            dspss(PMFdName,nv,AcNS,PMNQ);
    }
    if (PML0 != 0 || XSEPC != SEPC) {        /* reset format */
        for (k = 0; k < nv; ++k) {
            j = AcNS[k];
            w = (int)VPFmt1[j];             
            d = (int)VPFmt2[j];             
            makefmt(&w,&d,VPFmtS[j],0,SEPC,0);
        }
    }
    err = 0;

PDATFin:
    if (PM1NV > 0)                     
        vsort(0,PM1VIdx,0,0,1);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  dtda(fname,noc,nv,vidx,nq)                                              */
 
void dtda(char *fname,int noc,int nv,short *vidx,int nq)
{
    register int i,j,n,ni;

    fprintf(PMTDAFd,"nvar(\n");
    fprintf(PMTDAFd,"  dfile = %s,\n",fname);
    fprintf(PMTDAFd,"  noc = %d,\n",noc);
    ni = 1;
    for (n = 1; n <= nq; ++n) {
        for (i = 0; i < nv; ++i) {
            j = vidx[i];
            fprintf(PMTDAFd,"  %s",VName[j]);
            if (nq > 1)
                fprintf(PMTDAFd,"%d",n);                                                 

            fprnchar(PMTDAFd,' ',VNameLen - strlen(VName[j]),0);
            fprintf(PMTDAFd,"<%d>[%d.%d]",VSLen[j],VPFmt1[j],VPFmt2[j]);

            if (VLabel[j] != NULL)                                                      
                fprintf(PMTDAFd,"(%s)",VLabel[j]);                                           
            fprintf(PMTDAFd," = c%d,\n",ni++);
        }
    }
    fprintf(PMTDAFd,");\n");
    
    printf1("TDA description written to: %s\n",PMTDAFName);
}

/* ------------------------------------------------------------------------ */
/*  dspss(fname,noc,nv,vidx,nq)                                             */
 
void dspss(char *fname,int nv,short *vidx,int nq)
{
    register int i,j,n,nn,ii;

    fprintf(PMSPSSFd,"DATA LIST FILE='%s' FREE/\n",fname);                                              
    ii = 0;
    for (nn = 1; nn <= nq; ++nn) {
        for (i = 0; i < nv; ++i) {                                                      
            j = vidx[i];                                                                
            fprintf(PMSPSSFd," %s",VName[j]);                                                 
            if (nq > 1)
                fprintf(PMSPSSFd,"%d",nn);                                                 
            ii++;
            n = ii / 8;
            if (8 * n == ii) {
                fprintf(PMSPSSFd,"\n");
                ii = 0;
            }     
        }
    }
    fprintf(PMSPSSFd,".\nVARIABLE LABELS\n");
    for (nn = 1; nn <= nq; ++nn) {
        for (i = 0; i < nv; ++i) {
            j = vidx[i];
            fprintf(PMSPSSFd," %s",VName[j]);
            if (nq > 1)
                fprintf(PMSPSSFd,"%d",nn);                                                 
            fprnchar(PMSPSSFd,' ',VNameLen - strlen(VName[j]),0);
            if (VLabel[j] != NULL)                                                      
                fprintf(PMSPSSFd," '%s'",VLabel[j]);                                           
            else
                fprintf(PMSPSSFd," '%s'",VName[j]);                                           
            if (i == nv - 1 && nn == nq)
                fprintf(PMSPSSFd,".\n");
            else
                fprintf(PMSPSSFd,"/\n");
        }
    }                                                                               
    printf1("SPSS description written to: %s\n",PMSPSSFName);
} 

/* -##--------------------------------------------------------------------- */
/*  pdatr() print data for all pairs (i,j) of cases.                        */
/*                                                                          */
/*          pdatr (                                                         */
/*              keep = varlist,                                             */
/*              drop = varlist,                                             */
/*              nfmt = ...,     print format for i,j, def. 4                */  
/*          ) = fname;          output file                                 */
/*                                                                          */
/*              return 0 if OK, -1 if error.                                */
 
int pdatr(void)
{
    register int i,j,k,jv;
    int err,nrec,nv;
    double tmp;

    err = -1;
    if (check_cmd(0))
        return(-1);

    if (DMDef == 0) {
        p_err(-18,1);
        return(0);
    }
    if (parm(CmdBuf + 5,1,1)) {
        p_clean();
        return(-1);
    }
    if (PMKeep && PMDrop) {
        p_err(-16,1);
        goto PDATRFin;
    }

    /* get list of variables in AcNS[] */

    if (alloc_acns(imax(NVAR,PMNV)))
        goto PDATRFin;

    nv = 0;
    if (PMKeep) {
        for (k = 0; k < PMNV; ++k)           
            AcNS[nv++] = PMVIdx[k];              
    }
    else {
        j = VIFirst;
        while (j >= 0) {
            i = 1;
            if (PMDrop) {
                for (k = 0; k < PMNV; ++k) {
                    if (PMVIdx[k] == j) {
                        i = 0;
                        break;
                    }   
                }
            }
            if (i) 
                AcNS[nv++] = j;
            j = VNxt[j];
        }
    }
       
    if (nv == 0) {
        printf1("No variables selected.\n");
        goto PDATRFin;
    }
    if (SILENTFlg < 2)
        printfe("Writing: %s\n",PMFdName);

    nrec = 0;
    for (i = 0; i < NOC; ++i) {
        for (j = 0; j < NOC; ++j) {
            fprintf(PMFd,PMNFmtS,i+1);  
            fprintf(PMFd,PMNFmtS,j+1);  

            for (k = 0; k < nv; ++k) {
                jv = AcNS[k];
                if (VTyp[jv] != 1) {
                    tmp = get_data(jv,i);                  
                    fprintf(PMFd,VPFmtS[jv],tmp);
                }
                else {          /* string variables */
                    get_str(SVBuf,jv,i);
                    fprintf(PMFd,"%s",SVBuf);
                    if (XSEPC)
                        fprintf(PMFd,"%c",XSEPC);
                }
            }
            for (k = 0; k < nv; ++k) {
                jv = AcNS[k];
                if (VTyp[jv] != 1) {
                    tmp = get_data(jv,j);                  
                    fprintf(PMFd,VPFmtS[jv],tmp);
                }
                else {          /* string variables */
                    get_str(SVBuf,jv,j);
                    fprintf(PMFd,"%s",SVBuf);
                    if (XSEPC)
                        fprintf(PMFd,"%c",XSEPC);
                }
            }
            fprintf(PMFd,"\n");
            nrec++;
        }
    }
    prn_message(nrec,1,1);
    printf1("%d records written to: %s\n",nrec,PMFdName);
    err = 0;

PDATRFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  pdatd() Create a distance matrix.                                       */
/*                                                                          */
/*          pdatd (                                                         */
/*              v=...,      optional varlist, def. all variables            */
/*              opt=...,    option, def. 1                                  */
/*                           1 euclidean distances                          */
/*                           2 city-block distance                          */
/*                           3 number of different values                   */
/*                           4 dissimilarity index                          */

/*              wt=...,      optional weights, def. none                    */
/*              fmt = ...,   print format, def. 10.4                        */  
/*          ) = fname;          output file                                 */
/*                                                                          */
/*          Return 0 if OK, -1 if error.                                    */
 
int pdatd(void)
{
    register int i,j,k;
    int err,nv;
    double d;

    err = -1;
    if (check_cmd(0))
        return(-1);

    printf1("Create a distance matrix. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 5,1,1)) {
        p_clean();
        return(-1);
    }
    if (PMFmtF == 0)
        pmfmt(10,4);

    if (PMOPT == 2)  
        printf1("City-block distance.\n");
    else if (PMOPT == 3)  
        printf1("Number of different values.\n");
    else if (PMOPT == 4)  
        printf1("Dissimilarity index.\n");
    else {
        PMOPT = 1;    
        printf1("Euclidean distance.\n");
    }   
    if (PMNV > 0)
        nv = PMNV;
    else
        nv = NVAR;

    if (alloc_aci(nv))
        goto PDATDFin;

    if (PMNV > 0) {
        for (j = 0; j < PMNV; ++j)
            AcI[j] = PMVIdx[j];
    }
    else { 
        j = 0;
        i = VIFirst;
        while (i >= 0) {
            AcI[j++] = i;
            i = VNxt[i];
        }
    }
    printf1("Variables: %s",VName[AcI[0]]);
    for (k = 1; k < nv; ++k)  
        printf1(", %s",VName[AcI[k]]);
    newline();
    newline();

    if (PMNTP > 0) {
        if (PMNTP != nv) {
            printf1("Error: would need %d weights.\n",nv);
            goto PDATDFin;
        }
    }
    if (alloc_acx(NOC * NOC + 1))
        goto PDATDFin;

    if (alloc_acu(nv + 1))
        goto PDATDFin;
    if (alloc_acv(nv + 1))
        goto PDATDFin;

    for (i = 1; i < NOC; ++i) {
        for (k = 0; k < nv; ++k)  
            AcU[k] = get_data(AcI[k],i);

        if (PMOPT == 4)         /* change into relative values */
            pdatdr(nv,AcU);

        for (j = 0; j < i; ++j) {
            for (k = 0; k < nv; ++k)  
                AcV[k] = get_data(AcI[k],j);

            if (PMOPT == 4)                                               
                pdatdr(nv,AcV);

            AcX[i * NOC + j] = pdatd1(PMOPT,nv,AcU,AcV,PMNTP,PMTP);
        }
    }

    for (i = 0; i < NOC; ++i) {
        for (j = 0; j < NOC; ++j) {
            if (j <= i)
                d = AcX[i * NOC + j];
            else           
                d = AcX[j * NOC + i];
            fprintf(PMFd,PMFmtS,d);  
        }
        fprintf(PMFd,"\n");  
    }
    printf1("%d records written to: %s\n",NOC,PMFdName);
    err = 0;

PDATDFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  pdatdr(n,x)     Change x into relative values                           */

void pdatdr(int n,double *x)
{
    register int i;
    double s = 0.0;

    for (i = 0; i < n; ++i)
        s += x[i];
    if (s != 0.0) {
        for (i = 0; i < n; ++i)
            x[i] /= s;
    }
}

/* -##--------------------------------------------------------------------- */
/*  pdatd1(opt,n,x,y,nw,w)                                                  */
/*                                                                          */
/*  Return distance of n-vectors x and y, depending on opt.                 */
/*  opt = 1 : euclidean distance                                            */
/*        2 : city-block distance                                           */
/*        3 : number of different values                                    */
/*        4 : dissimilarity index (x and y already contain rel. values)     */
/*                                                                          */
/*  If nw > 0 use weights w[k], k = 0,1,...,n-1                             */
  
double pdatd1(int opt,int n,double *x,double *y,int nw,double *w)
{
    register int i;
    double tmp,d,wt;

    wt = 1.0;
    d = 0.0;
    for (i = 0; i < n; ++i) {
        if (nw > 0)
            wt = w[i];

        tmp = x[i] - y[i];

        if (opt == 1)  
            d += wt * tmp * tmp;
        else if (opt == 2 || opt == 4)
            d += wt * fabs(tmp);
        else if (opt == 3) {
            if (fabs(tmp) > EPSI1)
                d += wt;
        }
    }
    if (opt == 1)
        d = sqrt(d);
    else if (opt == 4)
        d /= 2.0;       

    return(d);
}







