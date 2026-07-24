/****************************************************************************/
/*  t_gg                                                                    */
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
#include "t_ml.h"
#include "t_gf.h"
#include "t_gio.h"
#include "t_com.h"
#include "tda_context.h"

/*  functions in t_gg.c */

int giset(TDAContext *ctx);
void gis_gmis(TDAContext *ctx, int *adj,int *adjn,int *adjn0,double alpha,double beta, int *biter,int *bkset,int *btup,int *deg,int *deg0,int *errcnd,int *fdmiq, int *fdmq,int *fwdmap,int *incdij,int *iniad,int *iniad0,int *invmap, int *iq,int *kset,int *kset0,int look4,int lstype, int lswhen,int max2a,int max2n2,int maxa,int maxitr,int maxn,int maxn2, int maxt,int maxtd,int narcs0,int *nkset0,int nlowdg,int nnode0,int *node1, int *node2,int *ptrn,int *ptrn0,int *q,int *rcl,int *seed, int tupdim,int *tupind,int *tuple,int tuplmt,int *tuplst,int *vtup,int *gitr);
void gis_mkds(TDAContext *ctx, int *adjn0,int *deg0,int *incdij,int *iniad0,int max2a,int maxa, int maxn,int maxn2,int narcs0,int nnode0,int *node1,int *node2,int *ptrn0);
void gis_mktup(TDAContext *ctx, int *adj,int *adjn0,double beta,int *deg0,int *errcnd,int *fdmiq, int *fdmq,int *incdij,int *iniad0,int *iq,int *lowdeg,int max2a,int max2n2, int maxn,int maxn2,int maxt,int maxtd,int nlowdg,int nnode0,int *ntuple, int *ptrn0,int *q,int *tupdim,int *tupind,int *tuple,int *tuplst,int *vtup);
void gis_tupfdm(TDAContext *ctx, int *adj,int *adjn0,int *errcnd,int *fdmiq,int *fdmq, int *incdij,int *iniad0,int *lowdeg,int max2a,int max2n2,int maxn2, int maxn,int maxt,int maxtd,int *minfdm,int *nfdmq,int nnode0,int *ntuple, int *ptrn0,int tupdim,int *tupind,int *tuple,int *vtup);
void gis_mkrgrf(TDAContext *ctx, int *adjn,int *adjn0,int *deg,int *deg0,int *fwdmap,int *iniad, int *iniad0,int *invmap,int max2a,int max2n2,int maxn,int maxtd,int *narcs, int *nnode,int nnode0,int *ptrn,int *ptrn0,int tupdim,int *vthere,int *vtup);
void gis_grspss(TDAContext *ctx, int *adjn,double alpha,int *biter0,int *bkset,int *deg, int *degv,int *gitr,int *incdij,int *iniad,int *invmap, int *kset,int look4,int lstype,int lswhen,int max2a,int max2n2,int maxn, int maxn2,int maxitr,int *nkset,int nkset0,int nnode,int nnode0, int *ptrn,int *rcl,int *seed,int tpl,int tupdim,int *vset);
void gis_iprn(TDAContext *ctx, int tpl,int iter,int size,int best);
void gis_sprn(TDAContext *ctx, int n,int *kset);                       
void gis_build(TDAContext *ctx, int *adjn,double alpha,int *bkset,int *deg,int *iniad, int max2a,int maxn,int *nbkset,int nnode,int *ptrn,int *rcl,int *seed);
void gis_gls12(TDAContext *ctx, int *adjn,int *bkset,int *incdij,int *iniad,int *invmap, int max2a,int maxn,int maxn2,int *nbkset,int nnode, int nnode0,int *ptrn,int *vset);
void gis_gls23(TDAContext *ctx, int *adjn,int *bkset,int *incdij,int *iniad,int *invmap, int max2a,int maxn,int maxn2,int *nbkset,int nnode,int nnode0, int *ptrn,int *vset);
void gis_insrtq(TDAContext *ctx, int dimq,int *iq,int iv,int *q,int *sizeq,int v);
void gis_removq(TDAContext *ctx, int dimq,int *iq,int *iv,int *q,int *sizeq,int *v);
void gis_cpi4(TDAContext *ctx, int dimn,int n,int *x,int *y);
double gis_randp(TDAContext *ctx, int *ix);

int gcni(TDAContext *ctx);
int g_cni_n(TDAContext *ctx, int i,int d,int *nodes,char *nflag,int gn);
void g_cni_prn(TDAContext *ctx, int i,int ne,int n,int *nodes);
int g_cni_ne(TDAContext *ctx, int n,int *nodes,int gn);
int g_cni_c(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int ne1,int opt,int gn);
int g_cni_ci(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int gn);
int g_cni_ci1(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int gn);
int gcset(TDAContext *ctx);
int g_cset(TDAContext *ctx, int m,int *in,int *jn,int cn,int opt,int *nr,int dflag);
int g_cset_df(TDAContext *ctx, int n,int *nodes);
double g_mineval(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int opt);
int gcliq(TDAContext *ctx);
void c_visit(TDAContext *ctx, int k,int gn);
void c_visit1(TDAContext *ctx, int k,int gn);
int g_cliq(TDAContext *ctx, int nc,int n,int *nodes,int min);
int g_cliq_1(TDAContext *ctx, int *old,int ne,int ce,int *nmap,int nc,int min);
int g_cliq_2(TDAContext *ctx, int nc,int n,int *nodes,int min);
int g_cliq_2getd(TDAContext *ctx, int i,int nmax,char *nflag,char *nflag1,int n,int *d,int *rs);
int ggcliq(TDAContext *ctx);
int g_gclfind(TDAContext *ctx, int in,int n,int *nodes,int gn,double s,int *nflag, int level,int min,int *itmp);
int g_gcladd(TDAContext *ctx, int i,int n,int *nodes,int gn,double s);
int g_gclptr_i(TDAContext *ctx, int max);
int g_gclptr_a(TDAContext *ctx, int n,int *nodes);
int g_gclptr_aa(TDAContext *ctx, int n,int *nodes,int *itmp);
int g_gclptr_c(TDAContext *ctx, int n,int *nodes);
int g_gclptr_cc(TDAContext *ctx, int n,int *nodes,int is);
void g_sp1(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec,double s);
void g_sp2(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec,double s,int nn,int *nodes);
int g_nfind(TDAContext *ctx, int i,int n,int *nodes);
int g_gclcheck(TDAContext *ctx, int n,int *nodes,int gn,double s);
int g_gclcheck1(TDAContext *ctx, int n,int *nodes,int gn,double s);
int g_gclclubs(TDAContext *ctx, int n,int *nodes,int gn,double s,int *idx);
void g_gclprn(TDAContext *ctx, int n,int *nodes);
int g_gclfind1(TDAContext *ctx, int in,int gn,double s,int min,int *idx,char *iflag,int *itmp);
void xxx(TDAContext *ctx);
int comb_xx(TDAContext *ctx, int n,int m,int *com,int first,int *idx);




/* ------------------------------------------------------------------------ */
/*  giset       Independent sets.                                           */
/*              Requires an undirected graph.                               */
/*                                                                          */
/*              giset(                                                      */
/*                  mxit=...,       def. 100                                */
/*                  min=...,        min requested size                      */
/*                  seed=...,       def. 270001                             */
/*                  dim=...,        tupdim, def. 0                          */
/*                  gn  = ...,      graph number, def. 1                    */
/*                  cg  = ...,      1 to use compl. graph, def. 0           */
/*                  nfmt=...,       integer format, def. 4                  */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int giset(TDAContext *ctx)
{
    register int i,j,k;
    int err,n,in,jn,lstype,lswhen,gitr,biter,btup;
    int nnode0,narcs0,maxt,maxtd;
    int seed,nkset0,maxa,max2a,maxn,maxn2,max2n2,errcnd,nlowdg,tupdim,tuplmt;
    int *node1,*node2,*adj,*adjn,*adjn0,*bkset,*deg,*deg0,*fdmiq,*fdmq;
    int *fwdmap,*incdij,*iniad,*iniad0,*invmap,*iq,*kset,*kset0,*ptrn,*ptrn0;
    int *q,*rcl,*tupind,*tuple,*tuplst,*vtup;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);

    printf1(ctx, "Independent node sets. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GISFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GISFin;
    if (gdd_tcheck(ctx, 0,1,0))          /* undirected graph */
        goto GISFin;

    if (ctx->MxIter < 1 || ctx->MxItFlg == 0)
        ctx->MxIter = 100;
    if (ctx->PMIDFA <= 0.0 || ctx->PMIDFA > 1.0)
        ctx->PMIDFA = 0.1;
    if (ctx->PMIDFB <= 0.0 || ctx->PMIDFB > 1.0)
        ctx->PMIDFB = 0.1;

    narcs0 = 0;                     /* calculate number of edges */
    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = i + 1; j < ctx->GD_NP; ++j) {
            if (gdd_adj(ctx, i,j,ctx->PMGN) >= 0.0) {
                if (ctx->PMCG != 1)
                    narcs0++;
            }
            else if (ctx->PMCG == 1)
                narcs0++;
        }
    }
    nnode0 = ctx->GD_NP;             /* number of nodes */

    if (ctx->PMMin < 1 || ctx->PMMin > nnode0)
        ctx->PMMin = nnode0;
    if (ctx->PMSEED >= 1 && ctx->PMSEED < ctx->INTMAX)
        seed = ctx->PMSEED;
    else
        seed = 270001;

    if (ctx->PMDIM < 0 || ctx->PMDIM >= nnode0)       /* dimension of tuple vector */
        ctx->PMDIM = 0;
    tupdim = ctx->PMDIM;
    
    /*  local search parameters 
        If lstype = 0, there is no local search.
                  = 1, local search is (1,2)-exchange.
                  = 2, local search is (2,3)-exchange.
 
        If lswhen = 0, local search is done only if phase 1 solution
                       is better than the average phase 1 solution.
                  = 1, local search is no matter what phase 1 
                       solution was. */ 

    lstype = 1; 
    lswhen = 0; 


    gitr = 0;               /* number of iterations */
    biter = 0;              /* number of iterations for best solution */
    btup = 0;               /* number of tuples to find best set (output) */
    nkset0 = 0;             /* number of nodes in ind set */
    errcnd = 0;             /* error condition */
    nlowdg = 50;            /* number of low degree nodes to make tuples */
    tuplmt = nnode0;        /* max number of tuples */
    maxt = nnode0;          /* max number of tuples in freedom queue */
    maxtd = imax(ctx, 1,tupdim); /* maxtd >= 1 and maxtd >= tupdim required */
    maxn = nnode0;
    maxn2 = maxn * maxn;
    max2n2 = 2 * maxn2;
    maxa = narcs0;
    max2a = 2 * maxa;

    n = 2 * maxa + 4 * maxn2 + 4 * max2a + 10 * maxn + 2 * maxt + 2 * maxtd + 2 * max2n2 + 2;
    if (alloc_acn(ctx, n))
        goto GISFin;

    node1  = ctx->AcN;
    node2  = node1 + maxa;
    adj    = node2 + maxa;
    adjn   = adj + maxn2;
    adjn0  = adjn + max2a;
    bkset  = adjn0 + max2a;
    deg    = bkset + maxn;
    deg0   = deg + maxn;
    fdmiq  = deg0 + maxn;
    fdmq   = fdmiq + maxt;
    fwdmap = fdmq + maxt;
    incdij = fwdmap + maxn;
    iniad  = incdij + maxn2;
    iniad0 = iniad + maxn;
    invmap = iniad0 + maxn;
    iq     = invmap + maxn;
    kset   = iq + maxn2;
    ptrn   = kset + maxn;
    ptrn0  = ptrn + max2a;
    q      = ptrn0 + max2a;
    rcl    = q + maxn2;
    tupind = rcl + maxn;
    tuple  = tupind + maxtd;
    tuplst = tuple + max2n2;
    vtup   = tuplst + max2n2;                 
    kset0  = vtup + maxtd;      /* maxn */

    /*  create edge list in node1[],node2[] */ 

    k = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {
        for (j = i + 1; j < ctx->GD_NP; ++j) {
            if (gdd_adj(ctx, i,j,ctx->PMGN) >= 0.0) {
                if (ctx->PMCG != 1) {
                    k++;
                    node1[k] = i + 1;
                    node2[k] = j + 1;
                }
            }   
            else if (ctx->PMCG == 1) {
                k++;
                node1[k] = i + 1;
                node2[k] = j + 1;
            }
        }
    }
    if (ctx->PMCG == 1)  
        printf1(ctx, "Using complementary graph with %d edges (%d).\n",narcs0,k);

    printf1(ctx, "Max number of iterations: %d (alpha %g, beta %g)\n",
        ctx->MxIter,ctx->PMIDFA,ctx->PMIDFB);
    printf1(ctx, "Seed of random number generator: %d\n",seed);
    printf1(ctx, "Dimension of tuple vector: %d\n",tupdim);
    if (ctx->PMMin < nnode0)
        printf1(ctx, "Requested minimal size: %d\n",ctx->PMMin);
    newline(ctx);

    
    gis_gmis(ctx, adj,adjn,adjn0,ctx->PMIDFA,ctx->PMIDFB,&biter,bkset,&btup,deg,deg0,
        &errcnd,fdmiq,fdmq,fwdmap,incdij,iniad,iniad0,invmap,iq,kset,
        kset0,ctx->PMMin,lstype,lswhen,max2a,max2n2,maxa,ctx->MxIter,maxn,maxn2,
        maxt,maxtd,narcs0,&nkset0,nlowdg,nnode0,node1,node2,ptrn,ptrn0,q,
        rcl,&seed,tupdim,tupind,tuple,tuplmt,tuplst,vtup,&gitr);
     
    /*****
    if (errcnd != 0)  
        printf1(ctx, "errcnd=%d\n",errcnd);
    ***/

    printf1(ctx, "Number of iterations performed: %d\n",gitr);
    printf1(ctx, "Found maximal independent set with %d elements in iteration %d.\n",nkset0,biter);
    printf1(ctx, "Found on tuple: %d\n",btup);

    /* check for ind. set or clique */

    k = 0;
    for (i = 1; i < nkset0; ++i) {
        in = kset0[i] - 1;
        for (j = i + 1; j <= nkset0; ++j) {
            jn = kset0[j] - 1;
            if (gdd_adj(ctx, in,jn,ctx->PMGN) >= 0.0) {
                if (ctx->PMCG != 1) {
                    k = 1;
                    break;
                }
            }
            else if (ctx->PMCG == 1) {
                k = 1;
                break;
            }
        }
        if (k) {
            printf1(ctx, "Check failed (in=%d,jn=%d).\n",gdd_node(ctx, in),gdd_node(ctx, jn));
            break;
        }
    }

    for (i = 1; i <= nkset0; ++i) {
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,kset0[i]);
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, kset0[i] - 1));
        fprintf(ctx->PMFd,"\n");
    }
    printf1(ctx, "%d records written to: %s\n",nkset0,ctx->PMFdName);
    err = 0;

GISFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  gis_gmis()      finds a large independent set.                          */
/*                                                                          */
/*  The code is adapted from a FORTRAN program developted by                */
/*  Mauricio G.C. Resende, Thomas A. Feo, and Stuart H. Smith.              */
/*                                                                          */
/*  See the paper:                                                          */
/*                                                                          */
/*  Thomas A. Feo, Mauricio G.C. Resende, and S.H. Smith, "A greedy         */
/*  randomized adaptive search procedure for maximum independent            */
/*  set," Operations Research, vol. 42, pp. 860--878, 1994.                 */
/*                                                                          */
/*  Passed input scalars:                                                   */
/*  alpha  - GRASP RCL parameter                                            */
/*  beta   - tuple list size restriction parameter                          */
/*  look4  - size of independent set sought                                 */
/*  lstype - type of local search                                           */
/*  lswhen - indicates when local search is done                            */
/*  max2a  - array dimension = 2*maxa                                       */
/*  max2n2 - array dimension = 2*maxn2                                      */
/*  maxa   - array dimension                                                */
/*  maxitr - maximum number of GRASP iterations                             */
/*  maxn   - array dimension                                                */
/*  maxn2  - array dimension                                                */
/*  maxt   - maximum number of tuples in freedom queue                      */
/*  maxtd  - array dimension                                                */
/*  narcs0 - number of arcs in original graph                               */
/*  nlowdg - number of low degree vertices used to make tuples              */
/*  nnode0 - number of nodes in original graph                              */
/*  tupdim - dimension of tuple vector                                      */
/*  tuplmt - maximum number of tuples allowed                               */
/*                                                                          */
/*  Passed input/output scalars:                                            */
/*  errcnd - error condition                                                */
/*  seed   - seed for pseudo random number generator                        */
/*                                                                          */
/*  Passed working scalars:                                                 */
/*  gitr   - number of grasp iterations taken                               */
/*  narcs  - number of arcs in reduced graph                                */
/*  nkset  - number of vertices in independent set of reduced               */
/*           graph                                                          */
/*  nnode  - number of vertices in reduced graph                            */
/*  ntuple - number of tuples produced by subroutine mktup                  */
/*                                                                          */
/*  Passed output scalars:                                                  */
/*  biter  - number of iterations to find best set                          */
/*  btup   - number of tuples to find best set                              */
/*  nkset0 - number of vertices in independent set                          */
/*                                                                          */
/*  Passed input arrays:                                                    */
/*  node1  - vertex 1 of egde                                               */
/*  node2  - vertex 2 of edge                                               */
/*                                                                          */
/*  Passed working arrays:                                                  */
/*  adj    - indicator array of adjacent nodes                              */
/*  adjn   - vertex in adjacency list (copy)                                */
/*  adjn0  - vertex in adjacency list (original)                            */
/*  bkset  - set of vertices in independent set                             */
/*  deg    - degree of vertex (copy)                                        */
/*  deg0   - degree of vertex (original)                                    */
/*  fdmiq  - freedom heap indices                                           */
/*  fdmq   - freedom heap values                                            */
/*  fwdmap - node i in original graph is fwdmap(i) in reduced               */
/*           graph                                                          */
/*  incdij - graph incidence matrix                                         */
/*  iniad  - pointer to start of adjacency list (copy)                      */
/*  iniad0 - pointer to start of adjacency list (original)                  */
/*  invmap - node i in reduced graph is invmap(i) in original               */
/*           graph, passed to mkrgrf, grspss                                */
/*  iq     - heap index array                                               */
/*  kset   - set of independent nodes                                       */
/*  ptrn   - pointer to next element of adjacency list (copy)               */
/*  ptrn0  - pointer to next element of adjacency list (original)           */
/*  q      - heap value array                                               */
/*  rcl    - restricted candidate list                                      */
/*  tupind - index of tuple                                                 */
/*  tuple  - list of tuples                                                 */
/*  tuplst - list of tuples                                                 */
/*  vtup   - tuple vertices                                                 */
/*                                                                          */
/*  Passed output array:                                                    */
/*  kset0  - best set of independent nodes, in original graph               */

void gis_gmis(TDAContext *ctx, int *adj,int *adjn,int *adjn0,double alpha,double beta, int *biter,int *bkset,int *btup,int *deg,int *deg0,int *errcnd,int *fdmiq, int *fdmq,int *fwdmap,int *incdij,int *iniad,int *iniad0,int *invmap, int *iq,int *kset,int *kset0,int look4,int lstype, int lswhen,int max2a,int max2n2,int maxa,int maxitr,int maxn,int maxn2, int maxt,int maxtd,int narcs0,int *nkset0,int nlowdg,int nnode0,int *node1, int *node2,int *ptrn,int *ptrn0,int *q,int *rcl,int *seed, int tupdim,int *tupind,int *tuple,int tuplmt,int *tuplst,int *vtup,int *gitr)
{
    register int i;
    int tpl,k,narcs,nkset,nnode,ntuple,biter0,maxtup,nq;

    ntuple = nkset = 0;
    *errcnd = 0;

    if (ctx->SILENTFlg < 0)
        printfe(ctx, "\nTuple  Iteration   Size   Best\n");

    /* Create adjacency data structure. */

    gis_mkds(ctx, adjn0,deg0,incdij,iniad0,max2a,maxa,maxn,maxn2,narcs0,
         nnode0,node1,node2,ptrn0);

    /* Create list of tuples for processing reduced graph    
       G - VTUP - adjv(VTUP). 
       Array invmap is passed as a work array named lowdeg in mktup. */

    if (tupdim != 0) {     
        gis_mktup(ctx, adj,adjn0,beta,deg0,errcnd,fdmiq,fdmq,incdij,
              iniad0,iq,invmap,max2a,max2n2,maxn,maxn2,maxt,maxtd,
              nlowdg,nnode0,&ntuple,ptrn0,q,&tupdim,tupind,tuple,
              tuplst,vtup);
    }      

    /* Initialize nkset: number of vertices in max indep set.
       Limit number of tuples to be examined to at most tuplmt. */

    *nkset0 = 0;
    if (tupdim > 0)      
        maxtup = imin(ctx, tuplmt,ntuple);
    else
        maxtup = 1;

    /* For each tuple (vtup), run GRASP on reduced graph
       G - vtup - adjv(vtup). */

    maxitr = maxitr / maxtup;

    for (tpl = 1; tpl <= maxtup; ++tpl) {

        /* Identify vertices (vtup) of tuple. */

        for (k = 1; k <= tupdim; ++k) {
            vtup[k] = tuplst[tupdim * (tpl - 1) + k];
        }

        /* Make the reduced graph G - vtup - adjv(vtup)
           The graph is placed in data structure (nnode,narcs,iniad,
           ptrn,adjn,deg).
 
           Node i in reduced  graph is invmap(i) in original graph. 
           Node i in original graph is fwdmap(i) in reduced  graph. 
 
           Array tuple is passed as a work array named vthere in mkrgrf. */
    
        gis_mkrgrf(ctx, adjn,adjn0,deg,deg0,fwdmap,iniad,iniad0,invmap,
               max2a,max2n2,maxn,maxtd,&narcs,&nnode,nnode0,ptrn,
               ptrn0,tupdim,tuple,vtup);
      
        /* Run maxitr GRASP iterations on reduced graph.
           Largest independent set found is of size nkset and its
           vertices (of reduced graph) are returned in kset.
 
           Array tuple is passed as a work array that is called 
           degv in grspss.
 
           Array q is passed as a work array that is called 
           vset in grspss. */
          
        gis_grspss(ctx, adjn,alpha,&biter0,bkset,deg,tuple,gitr,incdij,
               iniad,invmap,kset,look4,lstype,lswhen,
               max2a,max2n2,maxn,maxn2,maxitr,&nkset,*nkset0,
               nnode,nnode0,ptrn,rcl,seed,tpl,tupdim,q);
    
        /* If the independent set is the largest found so far,
           unpack it and save it in array kset0.
           The size of the set is nkset0. */

        if (nkset + tupdim > *nkset0) {     
            for (i = 1; i <= nkset; ++i)  
                kset[i] = invmap[kset[i]];
               
            for (k = 1; k <= tupdim; ++k)  
                kset[nkset + k] = vtup[k];

            nkset = nkset + tupdim;
            *nkset0 = nkset;
            *btup = tpl;
            *biter = biter0;
            gis_cpi4(ctx, maxn,nkset,kset,kset0);

            /* If the set found is of the size (look4) sought or
               greater, return. */

            if (*nkset0 >= look4) {
                if (ctx->PMF1Def) 
                    gis_sprn(ctx, *nkset0,kset0);
                else
                    break;     
            }
        }
    }

    /* Sort the nkset0 indepedent vertices in kset0 to return. */

    nq = 0;
    for (i = 1; i <= *nkset0; ++i) {
        gis_insrtq(ctx, maxn2,iq,i,q,&nq,kset0[i]);
    }
    for (i = 1; i <= *nkset0; ++i) {
        gis_removq(ctx, maxn2,iq,&k,q,&nq,&kset0[i]);
    }
}

/* ------------------------------------------------------------------------ */
/*  gis_mkds()                                                              */
/*                                                                          */
/*  Builds adjacency data structure for graph, computes degree.             */

void gis_mkds(TDAContext *ctx, int *adjn0,int *deg0,int *incdij,int *iniad0,int max2a,int maxa, int maxn,int maxn2,int narcs0,int nnode0,int *node1,int *node2,int *ptrn0)
{
    (void)ctx; (void)max2a; (void)maxa; (void)maxn; (void)maxn2;        /* unused: the signature is shared */
    register int i,index;
    int n1,n2,ptradj;

    /* Initialize pointers to list of adjacent nodes and degrees. */

    ptradj = 0;
    for (i = 1; i <= nnode0; ++i)  
        iniad0[i] = deg0[i] = 0;

    /* Initialize adjacency matrix. */

    for (index = 1; index <= nnode0 * nnode0; ++index)  
        incdij[index] = 0;

    /* Build adjacency list and adjacency matrix and compute degrees */

    for (i = 1; i <= narcs0; ++i) {
        n1 = node1[i];
        n2 = node2[i];
        incdij[(n1 - 1) * nnode0 + n2] = 1;
        incdij[(n2 - 1) * nnode0 + n1] = 1;
        ptradj = ptradj + 1;
        adjn0[ptradj] = n2;
        ptrn0[ptradj] = iniad0[n1];
        iniad0[n1] = ptradj;
        deg0[n1] = deg0[n1] + 1;
        ptradj = ptradj + 1;
        adjn0[ptradj] = n1;
        ptrn0[ptradj] = iniad0[n2];
        iniad0[n2] = ptradj;
        deg0[n2] = deg0[n2] + 1;
    }
}
  
/* ------------------------------------------------------------------------ */
/*  gis_mktup()                                                             */
/*                                                                          */
/*  Computes the list of tuples having freedom between                      */
/*  minfdm and minfdm + beta * (maxfdm-minfdm), where:                      */
/*                                                                          */
/*  minfdm = min freedom over all pairs of nodes                            */
/*  maxfdm = max freedom over all pairs of nodes.                           */

void gis_mktup(TDAContext *ctx, int *adj,int *adjn0,double beta,int *deg0,int *errcnd,int *fdmiq, int *fdmq,int *incdij,int *iniad0,int *iq,int *lowdeg,int max2a,int max2n2, int maxn,int maxn2,int maxt,int maxtd,int nlowdg,int nnode0,int *ntuple, int *ptrn0,int *q,int *tupdim,int *tupind,int *tuple,int *tuplst,int *vtup)
{
    register int i,k;
    int index,sizeq,cutoff,dv,freedm,iv,maxfdm,minfdm,nfdmq,tpl;

    /* Sort nodes, using heap sort, in increasing order of degrees.
       Retrieve the nlowdg nodes with smallest degrees.  Place them in
       array lowdeg. */

    nfdmq = sizeq = *ntuple = 0;
    if (nlowdg > nnode0)
        nlowdg = nnode0;

    for (i = 1; i <= nnode0; ++i) {
        gis_insrtq(ctx, maxn2,iq,i,q,&sizeq,deg0[i]);
    }
    for (i = 1; i <= nlowdg; ++i) {
        gis_removq(ctx, maxn2,iq,&iv,q,&sizeq,&dv);
        lowdeg[i] = iv;
    }
    minfdm = nnode0;
    if (*tupdim > nlowdg)
        *tupdim = nlowdg;

    /* Generate all tuples of lowdeg nodes of size tupdim and compute 
       freedom for all those that are made up of nodes that are pairwise 
       nonadjacent.  Sort freedoms, using heap sort, in decreasing order.
 
       Generate first tuple = (1,2,...,tupdim). */

    sizeq = 0;
    for (i = 1; i <= *tupdim; ++i) {
        tupind[i] = i;
        sizeq++;            
        q[sizeq] = i;
    }

    gis_tupfdm(ctx, adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,lowdeg,
           max2a,max2n2,maxn2,maxn,maxt,maxtd,&minfdm,&nfdmq,
           nnode0,ntuple,ptrn0,*tupdim,tupind,tuple,vtup);

    /* begin one level deeper */

L40:
    if (sizeq == 0)
        goto L80;

    i = q[sizeq];
    tupind[i] += 1;      

    if (tupind[i] != nlowdg - *tupdim + i) {      

        for (k = i + 1; k <= *tupdim; ++k) {
            tupind[k] = tupind[i] + k - i;
            sizeq++;      
            q[sizeq] = k;
        }
    }
    else {
        for (k = i + 1; k <= *tupdim; ++k) {
            tupind[k] = tupind[i] + k - i;
        }

        gis_tupfdm(ctx, adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,
               lowdeg,max2a,max2n2,maxn2,maxn,
               maxt,maxtd,&minfdm,&nfdmq,nnode0,ntuple,ptrn0,
               *tupdim,tupind,tuple,vtup);

        if (sizeq == 0)       
            i = 0;
        else {
            i = q[sizeq];
            sizeq--;       
        }        
        goto L40;
    }
L70:
    gis_tupfdm(ctx, adj,adjn0,errcnd,fdmiq,fdmq,incdij,iniad0,lowdeg,
           max2a,max2n2,maxn2,maxn,maxt,maxtd,&minfdm,&nfdmq,
           nnode0,ntuple,ptrn0,*tupdim,tupind,tuple,vtup);

    if (tupind[i] == nlowdg - *tupdim + i) {       

        if (sizeq == 0)       
            i = 0;
        else {
            i = q[sizeq];
            sizeq--;       
        }           
        goto L40;
    }
    else {
        i = q[sizeq];
        tupind[i] += 1;      
        goto L70;
    }
L80:

    /* Setup the sorted tuple list.   
       Get first tuple (with largest freedom) from heap. */

    gis_removq(ctx, maxt,fdmiq,&tpl,fdmq,&nfdmq,&freedm);
    index = *tupdim * (tpl - 1);
    for (i = 1; i <= *tupdim; ++i)  
        tuplst[i] = tuple[index + i];

    /* Compute cutoff for inclusion in tuple list.
       All but one of the accepted tuples will have freedom between 
       minfdm and minfdm + beta*(maxfdm-minfdm). One tuple will have
       the largest freedom less than minfdm + beta*(maxfdm-minfdm). */

    maxfdm = -freedm;
    cutoff = (int)(minfdm + beta * (maxfdm - minfdm));

    /* Get all tuples with freedom >= cutoff.
       Get one tuple with largest freedom < cutoff.

       Node pairs making up tuple list are placed in array tuplst. */

    for (k = 2; k <= *ntuple; ++k) {

        gis_removq(ctx, maxt,fdmiq,&tpl,fdmq,&nfdmq,&freedm);

        for (i = 1; i <= *tupdim; ++i) {
            tuplst[*tupdim * (k - 1) + i] = tuple[*tupdim * (tpl - 1) + i];
        }
        if(-freedm < cutoff)
            goto L120;
    }

    /* Reset number of tuples. */

L120:
    *ntuple = k;
}
   
/* ------------------------------------------------------------------------ */
/*  gis_tupfdm()                                                            */
/*                                                                          */
/*  The freedom of a set of vertices is the number of vertices              */
/*  in the graph that are not adjacent to any member vertex                 */
/*  of the set.                                                             */
/*                                                                          */
/*  Compute the freedom of set define by lowdeg (tupind).                   */

void gis_tupfdm(TDAContext *ctx, int *adj,int *adjn0,int *errcnd,int *fdmiq,int *fdmq, int *incdij,int *iniad0,int *lowdeg,int max2a,int max2n2,int maxn2, int maxn,int maxt,int maxtd,int *minfdm,int *nfdmq,int nnode0,int *ntuple, int *ptrn0,int tupdim,int *tupind,int *tuple,int *vtup)
{
    (void)max2a; (void)max2n2; (void)maxn; (void)maxn2; (void)maxtd;        /* unused: the signature is shared */
    register int i,j,k;
    int freedm,ptr,sumadj;

    for (i = 1; i <= tupdim; ++i)  
        vtup[i] = lowdeg[tupind[i]];

    /* Check if vertices in set vtup are pairwise nonadjacent.
       Get next tuple if they are (return). */

    for (i = 1; i < tupdim; ++i) {
        for (j = i + 1; j <= tupdim; ++j) {
            if (incdij[(vtup[i] - 1) * nnode0 + vtup[j]] == 1)
                return;
        }   
    }

    /* Vertices in set vtup are pairwise nonadjacent. 
       Setup indicator arrays adj of adjacent nodes:
  
       adj(nnode0*(vtup(i)-1)+u) = 1 if node u is adjacent 
                                     to vtup(i)
                                 = 0 otherwise
 
       This array is used to compute freedom of the set vtup. */

    for (i = 1; i <= tupdim; ++i) {
        for (k = 1; k <= nnode0; ++k)  
            adj[nnode0 * (i - 1) + k] = 0;
    }
    for (i = 1; i <= tupdim; ++i) {

        ptr = iniad0[vtup[i]];
L60:
        if(ptr == 0)
            goto L70;
        adj[nnode0 * (vtup[i] - 1) + adjn0[ptr]] = 1;
        ptr = ptrn0[ptr];
        goto L60;
L70: ;                   
    }

    /* Compute freedom of tuple (vtup) and insert into heap 
       for sorting.
 
       Find minimum value of freedom as freedoms are computed. */

    freedm = 0;

    for (k = 1; k <= nnode0; ++k) {
        for (i = 1; i <= tupdim; ++i) {
            if (k == vtup[i])
                goto L110;
        }
        sumadj = 0;
        for (i = 1; i <= tupdim; ++i) {
            sumadj += adj[nnode0 * (vtup[i] - 1) + k];
        }
        if(sumadj == 0)
            freedm++;       
L110: ;            
    }

    *ntuple += 1;     

    for (i = 1; i <= tupdim; ++i)  
        tuple[tupdim * (*ntuple - 1) + i] = vtup[i];
      
    if (*ntuple > maxt)     
        *errcnd = 100;
    else {
        gis_insrtq(ctx, maxt,fdmiq,*ntuple,fdmq,nfdmq,-freedm);
        if (freedm < *minfdm)
            *minfdm = freedm;
    }                      
}

/* ------------------------------------------------------------------------ */
/*  gis_mkrgrf()                                                            */
/*                                                                          */
/*  Deletes nodes in vtup from the original graph generating                */
/*  the reduced graph for the GRASP iterations.                             */
/*                                                                          */
/*  Reduced graph is returned in data structure                             */
/*  (nnode,narcs,iniad,ptrn,adjn,deg).                                      */
/*                                                                          */
/*  Arrays (fwdmap and invmap) that map corresponding                       */
/*  vertices in original and reduced graphs are set up.                     */
/*                                                                          */
/*      + Node i in reduced  graph is invmap(i)                             */
/*        in original graph.                                                */
/*                                                                          */
/*      + Node i in original graph is fwdmap(i)                             */
/*      in reduced  graph.                                                  */
/*                                                                          */
  
void gis_mkrgrf(TDAContext *ctx, int *adjn,int *adjn0,int *deg,int *deg0,int *fwdmap,int *iniad, int *iniad0,int *invmap,int max2a,int max2n2,int maxn,int maxtd,int *narcs, int *nnode,int nnode0,int *ptrn,int *ptrn0,int tupdim,int *vthere,int *vtup)
{
    (void)ctx; (void)deg0; (void)max2a; (void)max2n2; (void)maxn; (void)maxtd;        /* unused: the signature is shared */
    register int i,k;
    int adjn0p,nz,ptr,ptrad,ri;

    /* Prepare forward (fwdmap) and inverse (invmap) maps. */

    for (i = 1; i <= nnode0; ++i)  
        vthere[i] = 1;
      
    if (tupdim > 0) {     

        for (k = 1; k <= tupdim; ++k)  
            vthere[vtup[k]] = 0;
         
        for (k = 1; k <= tupdim; ++k) {

            ptr = iniad0[vtup[k]];
L30:
            if (ptr == 0)
                goto L40;
            vthere[adjn0[ptr]] = 0;
            ptr = ptrn0[ptr];
            goto L30;
L40:        ;
        }
        nz = 0;
        for (i = 1; i <= nnode0; ++i) {
            if (vthere[i] == 0) {     
                nz++;    
                fwdmap[i] = 0;
            }
            else {
                fwdmap[i] = i - nz;
                invmap[i - nz] = i;
            }            
        }
    }
    else {
        nz = 0;
        for (i = 1; i <= nnode0; ++i) {
            fwdmap[i] = i;
            invmap[i] = i;
        }
    }

    /* Compress adjacency list. */

    *nnode = nnode0 - nz;

    for (i = 1; i <= *nnode; ++i)  
        iniad[i] = deg[i] = 0;
      
    *narcs = ptrad = 0;

    for (i = 1; i <= nnode0; ++i) {

        if (vthere[i] == 1) {     

            ptr = iniad0[i]; 
L90:
            if (ptr == 0)
                goto L100;

            ri = fwdmap[i];
            adjn0p = adjn0[ptr];

            if (vthere[adjn0p] == 1) {     
                ptrad++;      
                ptrn[ptrad] = iniad[ri];
                iniad[ri] = ptrad;
                adjn[ptrad] = fwdmap[adjn0p];
                deg[ri] = deg[ri] + 1;
                *narcs += 1;    
            }      
            ptr = ptrn0[ptr];
            goto L90;
L100:       ;                  
        }         
    }
    *narcs /= 2;    
}
   
/* ------------------------------------------------------------------------ */
/*  gis_grspss()                                                            */
/*                                                                          */
/*  Applies GRASP to construct an independent set in graph                  */
/*  G=(iniad,ptrn,adjn,deg).                                                */
/*                                                                          */
/*  Set of size nkset is returned in array kset.                            */

void gis_grspss(TDAContext *ctx, int *adjn,double alpha,int *biter0,int *bkset,int *deg, int *degv,int *gitr,int *incdij,int *iniad,int *invmap, int *kset,int look4,int lstype,int lswhen,int max2a,int max2n2,int maxn, int maxn2,int maxitr,int *nkset,int nkset0,int nnode,int nnode0, int *ptrn,int *rcl,int *seed,int tpl,int tupdim,int *vset)
{
    (void)max2n2;        /* unused: the signature is shared */
    int avgph1,nbkset;

    /* Initialize independent set size and average size of phase 1
       independent sets. */

    *nkset = nkset0 - tupdim;
    avgph1 = 0.0;  

    /* Save vertex degrees in degv.  Array deg is over-written and must
       be recovered to start each GRASP iteration. */

    gis_cpi4(ctx, maxn,nnode,deg,degv);

    /* Do maxitr GRASP iterations. */

    for (*gitr = 1; *gitr <= maxitr; *gitr += 1) {

        /* GRASP construction phase (phase 1). */

        gis_build(ctx, adjn,alpha,bkset,deg,iniad,max2a,maxn,&nbkset,
              nnode,ptrn,rcl,seed);

        /* Store best solution, if found. */

        if (nbkset > *nkset) {      
            *nkset = nbkset;
            *biter0 = *gitr;
            gis_cpi4(ctx, maxn,*nkset,bkset,kset);
            gis_iprn(ctx, tpl,*gitr,nbkset + tupdim,*nkset + tupdim);

            if (*nkset + tupdim >= look4) {
                if (ctx->PMF1Def) 
                    gis_sprn(ctx, *nkset,kset);                       
                else
                    return;    
            }
        }
        else  
            gis_iprn(ctx, tpl,*gitr,nbkset + tupdim,*nkset + tupdim);

        /* Save value of phase 1 solution and update average of phase 1
           solutions. */

        avgph1 = (avgph1 * (*gitr - 1) + nbkset) / *gitr;

        /* GRASP local search phase.
 
           If lstype = 0, there is no local search.
                     = 1, local search is (1,2)-exchange.
                     = 2, local search is (2,3)-exchange.
 
           If lswhen = 0, local search is done only if phase 1 solution
                          is better than the average phase 1 solution.
                     = 1, local search is no matter what phase 1 
                          solution was. */

        if (lstype != 0) {     

            if (lswhen == 0) {     

                if (nbkset >= avgph1) {     

                    if (lstype == 1) {     

                        /* (1,2)-exchange local search */

                        gis_gls12(ctx, adjn,bkset,incdij,iniad,
                                  invmap,max2a,maxn,maxn2,
                                  &nbkset,nnode,nnode0,ptrn,vset);
                    }
                    else {

                        /* (2,3)-exchange local search */

                        gis_gls23(ctx, adjn,bkset,incdij,iniad,
                                  invmap,max2a,maxn,maxn2,
                                  &nbkset,nnode,nnode0,ptrn,vset);
                    }        
                }     
            }
            else {
                if (lstype == 1) {     

                    /* (1,2)-exchange local search */

                    gis_gls12(ctx, adjn,bkset,incdij,iniad,
                              invmap,max2a,maxn,maxn2,
                              &nbkset,nnode,nnode0,ptrn,vset);
                }
                else {

                    /* (2,3)-exchange local search */

                    gis_gls23(ctx, adjn,bkset,incdij,iniad,invmap,
                              max2a,maxn,maxn2,&nbkset,nnode,
                              nnode0,ptrn,vset);
                }      
            }     
        }            

        /* Store best solution, if found. */

        if (nbkset > *nkset) {    
            *nkset = nbkset;
            gis_cpi4(ctx, maxn,*nkset,bkset,kset);
            /** gis_iprn(tpl,*gitr,nbkset + tupdim,*nkset + tupdim); */

            if (*nkset + tupdim >= look4) {
                if (ctx->PMF1Def) 
                    gis_sprn(ctx, *nkset,kset);                       
                else
                    return;    
            }
        }
        /**
        else  
            gis_iprn(ctx, tpl,*gitr,nbkset + tupdim,*nkset + tupdim);
        **/

        /* Restore array deg to its original state. */ 

        gis_cpi4(ctx, maxn,nnode,degv,deg);

    }
    *gitr = maxitr;   
} 

/* ------------------------------------------------------------------------ */
/*  gis_iprn()      print iteration info.                                   */

void gis_iprn(TDAContext *ctx, int tpl,int iter,int size,int best)
{
    if (ctx->SILENTFlg < 0)
        printfe(ctx, "%5d  %9d %6d %6d\n",tpl,iter,size,best);
}

/* ------------------------------------------------------------------------ */
/*  gis_sprn()  print current set to PMFd1                                  */

void gis_sprn(TDAContext *ctx, int n,int *kset)                        
{
    register int i;

    if (ctx->PMF1Def) {
        for (i = 1; i <= n; ++i)
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,kset[i]);
        fprintf(ctx->PMF1d,"\n");
    }
}
     
/* ------------------------------------------------------------------------ */
/*  gis_build()                                                             */
/*                                                                          */
/*  GRASP construction phase (phase 1)                                      */

void gis_build(TDAContext *ctx, int *adjn,double alpha,int *bkset,int *deg,int *iniad, int max2a,int maxn,int *nbkset,int nnode,int *ptrn,int *rcl,int *seed)
{
    (void)max2a; (void)maxn;        /* unused: the signature is shared */
    register int i;
    int cutoff,degi,maxdeg,mindeg,nrcl,nselct,ptr,ptrv,u,v,vfound,vrtx;

    /* Initialize size of indepedendent set */

    *nbkset = 0;

    /* Independent set is constructed, one vertex at a time. */

L10:
 
    /* Set indicator that candidate vertex exists to false. */

    vfound = 0;

    /* Find minimum and maximum degrees of all candidate nodes.
       In process, determine if construction phase is done, i.e.
       if there no more candidates to be included in the indep set. */

    mindeg = nnode;
    maxdeg = 0;

    for (i = 1; i <= nnode; ++i) {

        degi = deg[i];
        if (degi > -1) {     

            /* At least one vertex is a candidate to be in the
               independent set. */

            vfound = 1;
            if (degi < mindeg)      
                mindeg = degi;

            if (degi > maxdeg)        
                maxdeg = degi;
        }     
    }

    /* Test if construction phase is done. */

    if (vfound == 0)
        goto L80;

    /* Construction phase not done yet.
 
       Compute cut-off for restricted candidate list (RCL).
 
       A vertex will be put in the RCL if its degree with respect 
       to vertices no in the independent set is between mindeg and 
       mindeg + alpha*(maxdeg-mindeg). */

    cutoff = (int)(mindeg + alpha * (maxdeg - mindeg));

    /* Build the RCL.  Place the nrcl vertices in array rcl. */

    nrcl = 0;
    for (i = 1; i <= nnode; ++i) {
        degi = deg[i];
        if (degi > -1 && degi <= cutoff) {     
            nrcl++;     
            rcl[nrcl] = i;
        }                  
    }

    /* Select a candidate at random from RCL to be put in the
       independent set.
 
       There are nrcl elements in array rcl.   The value nselct
       is a pseudo random integer between 1 and nrcl, inclusive.
       Function randp returns a seed between 1 and 2147483647. */

    if (nrcl < 1)
        nrcl = 1;
    gis_randp(ctx, seed);
    nselct = 1 + *seed / (2147483647 / nrcl);

    /* Increment size of independent set and put vertex
       vrtx = rcl(nselct) in independent set array bkset. */

    *nbkset += 1;      
    vrtx = rcl[nselct];
    bkset[*nbkset] = vrtx;

    /* Adaptive component of GRASP. Adjust degrees of remaining 
       nodes, thus changing the greedy function.
 
       By convention, a vertex with degree < 0 is not a candidate.
       Either it is in the independent set or is adjacent to a
       vertex that is in the independent set. */

    deg[vrtx] = -1;

    /* Scan adjacent vertices of independent vertex vrtx. */

    ptr = iniad[vrtx];
L40:
    if (ptr == 0)
        goto L70;

    v = adjn[ptr];

    /* If adjacent vertex (v) is not in the independent set,
       scan its neighbors and decrease their degrees by one. */

    if (deg[v] >= 0) {     
        deg[v] = -1;
        ptrv = iniad[v];
L50:
        if (ptrv == 0)
            goto L60;
        u = adjn[ptrv];
        deg[u] -= 1;        
        ptrv = ptrn[ptrv];
        goto L50;
L60: ;
    }                  
    ptr = ptrn[ptr];
    goto L40;
L70: ;

    /* Begin a new construction phase iteration. */

    goto L10;

L80:    ;
}
  
/* ------------------------------------------------------------------------ */
/*  gis_gls12()                                                             */
/*                                                                          */
/*  GRASP local search phase: (1,2)-exchange                                */
/*                                                                          */
/*  For each vertex v in independent set S, examine all pairs               */
/*  of nonadjacent vertices (u,w) not in the independent and                */
/*  check if S - v + u + w is an independent set.                           */
/*                                                                          */
/*  If so, remove v from S and add u and w to S, thus                       */
/*  increasing the size of the independent set by one.                      */
/*                                                                          */
/*  Repeat until local optimum (with respect to this                        */
/*  neighborhood structure) is found.                                       */

void gis_gls12(TDAContext *ctx, int *adjn,int *bkset,int *incdij,int *iniad,int *invmap, int max2a,int maxn,int maxn2,int *nbkset,int nnode, int nnode0,int *ptrn,int *vset)
{
    (void)ctx; (void)max2a; (void)maxn; (void)maxn2;        /* unused: the signature is shared */
    register int i,k,k1,k2;
    int bksetk,index,invu1,invu2,nvset,ptr,u1,u2,v1,w;

    /* Local search of neighborhood N(S) of solution S starts. */

L10:
 
    /* For all vertices v1 in current indep set (bkset), check
       if there exist nonadjacent vertices (u1,u2) such that u1
       and u2 are not adjacent to bkset-v1. */

    for (i = 1; i <= *nbkset; ++i) {

        v1 = bkset[i];

        /* Initialize array that indicates whether a node is
           adjacent to set bkset-v1. */

        for (k = 1; k <= nnode; ++k)  
            vset[k] = 0;

        /* Mark node if it is adjacent to node in bkset-v1  */

        for (k = 1; k <= *nbkset; ++k) {

            /* If node is v1, get next node in bkset. */

            if (k != i) {        
                bksetk = bkset[k];

                /* Mark bkset node as adjacent and scan its 
                   neighbors, marking them as adjacent. */

                vset[bksetk] = 1;
                ptr = iniad[bksetk];
L30:
                if (ptr == 0)
                    goto L40;
                w = adjn[ptr];
                vset[w] = 1;
                ptr = ptrn[ptr];
                goto L30;
L40:            ;              
            }      
        }

        /* Put all nvset nonadjacent nodes in vset. */

        nvset = 0;
        for (k = 1; k <= nnode; ++k) {
            if (k != v1 && vset[k] == 0) {    
                nvset++;     
                vset[nvset] = k;
            }      
        }

        /* If there are enough vertices to improve independent set
           (i.e. nvset .ge. 2), see if improvement occurs. */

        if (nvset >= 2) {     

            /* For all pairs (u1,u2) in vset, check if u1 and u2
               are nonadjacent.  
               If so, remove v1 from S and add
               u1 and u2 to S, increasing size of independent set. */

            for (k1 = 1; k1 < nvset; ++k1) {

                u1 = vset[k1];
                invu1 = invmap[u1];
                    
                for (k2 = k1 + 1; k2 <= nvset; ++k2) {

                    u2 = vset[k2];

                    /* Check if u1 and u2 are nonadjacent. */ 

                    invu2 = invmap[u2];
                    index = (invu1 - 1) * nnode0 + invu2;
                    if (incdij[index] != 1) {     

                        /* Nodes u1 and u2 are nonadjacent:
                           remove v1 from the independent set 
                           and add u1 and u2 to set, thus 
                           increasing size of independent set 
                           by one. */

                        bkset[i] = u1;
                        *nbkset += 1;    
                        bkset[*nbkset] = u2;

                        /* Restart local search in the
                           neighborhood of new solution. */

                        goto L10;
                    }        
                }
            }
        }     
    }
}            
  
/* ------------------------------------------------------------------------ */
/*  gis_gls23()                                                             */
/*                                                                          */
/*   GRASP local search phase: (2,3)-exchange                               */
/*                                                                          */
/*   For each pair of verteices (v,t) in independent set S,                 */
/*   examine all triples of parwise nonadjacent vertices (u,w,y)            */
/*   not in the independent and check if S - v - t + u + w + y              */
/*   is an independent set.                                                 */
/*                                                                          */
/*   If so, remove v and t from S and add u, w, adn y to S, thus            */
/*   increasing the size of the independent set by one.                     */
/*                                                                          */
/*   Repeat until local optimum (with respect to this                       */
/*   neighborhood structure) is found.                                      */
 
void gis_gls23(TDAContext *ctx, int *adjn,int *bkset,int *incdij,int *iniad,int *invmap, int max2a,int maxn,int maxn2,int *nbkset,int nnode,int nnode0, int *ptrn,int *vset)
{
    (void)ctx; (void)max2a; (void)maxn; (void)maxn2;        /* unused: the signature is shared */
    register int i,j,k,k1,k2,k3;
    int bksetk,index,invu1,invu2,invu3,nvset,ptr,u1,u2,u3,w;

    /* Local search of neighborhood N(S) of solution S starts. */

L10:            

    /* For all pairs of vertices (v1,v2) in current indep set 
       (bkset), check if there exist pairwise nonadjacent vertices 
       (u1,u2,u3) such that u1, u2, and and u3 are not adjacent to 
       bkset-v1-v2. */

    for (i = 1; i < *nbkset; ++i) {


        for (j = i + 1; j <= *nbkset; ++j) {


            /* Initialize array that indicates whether a node is
               adjacent to set bkset-v1-v2. */

            for (k = 1; k <= nnode; ++k)  
                vset[k] = 0;

            /* Mark node if it is adjacent to node in bkset-v1-v2. */

            for (k = 1; k <= *nbkset; ++k) {

                /* Mark bkset node as adjacent. */

                bksetk = bkset[k];
                vset[bksetk] = 1;

                /* If node is not v1 nor v2, scan neighbors of 
                   bkset node, marking them as adjacent. */

                if (k != i && k != j) {     
                    ptr = iniad[bksetk];
L30:
                    if (ptr == 0)
                        goto L40;
                    w = adjn[ptr];
                    vset[w] = 1;
                    ptr = ptrn[ptr];
                    goto L30;
L40:                ;
                }     

            }

            /* Put all nvset nonadjacent nodes in vset. */

            nvset = 0;
            for (k = 1; k <= nnode; ++k) {
                if (vset[k] == 0) {     
                    nvset++;     
                    vset[nvset] = k;
                }                 
            }

            /* If there are enough vertices to improve independent
               set (i.e. nvset .ge. 3), see if improvement occurs. */

            if (nvset >= 3) {     

                for (k1 = 1; k1 < nvset - 1; ++k1) {

                    /* For all triples (u1,u2,u3) in vset, check
                       if u1, u2, and u3 are pairwise nonadjacent.  
                       If so, remove (v1,v2) from S and add
                       (u1,u2,u3) to S, thus increasing size of 
                       the independent set by one. */

                    u1 = vset[k1];

                    for (k2 = k1 + 1; k2 < nvset; ++k2) {
                            
                        u2 = vset[k2];

                        for (k3 = k2 + 1; k3 <= nvset; ++k3) {
                                
                            u3 = vset[k3];

                            /* Check if u1, u2, and u3 are  pairwise
                               nonadjacent. */

                            invu1 = invmap[u1];
                            invu2 = invmap[u2];
                            invu3 = invmap[u3];
                            index = (invu1 - 1) * nnode0 + invu2;
                            if (incdij[index] == 1)
                                goto L70;
                            index = (invu1 - 1) * nnode0 + invu3;
                            if (incdij[index] == 1)
                                goto L70;
                            index = (invu2 - 1) * nnode0 + invu3;
                            if (incdij[index] == 1)
                                goto L70;

                            /* Nodes u1, u2, and u3 are pairwise
                               nonadjacent.     
                               Remove v1, v2 from independent set 
                               and add u1, u2, and u3, thus 
                               increasing size of independent set 
                               by one. */

                            bkset[i] = u1;
                            bkset[j] = u2;
                            *nbkset += 1;    
                            bkset[*nbkset] = u3;

                            /* Restart local search in neighborhood 
                               of new solution. */

                            goto L10;
L70:    ;
                        }         
                    }           
                }
            }              
        }
    }
}         
   
/* ------------------------------------------------------------------------ */
/*  gis_insrtq(dimq,iq,iv,q,sizeq,v)                                        */
/*                                                                          */
/*  Insert an element (v,iv) into a queue (q,iq).                           */

void gis_insrtq(TDAContext *ctx, int dimq,int *iq,int iv,int *q,int *sizeq,int v)
{
    (void)ctx; (void)dimq;        /* unused: the signature is shared */
    int sq,tsz;

    /* Insert element into heap. */

    *sizeq += 1;   
    q[*sizeq] = v;
    iq[*sizeq] = iv;

    /* Update heap to proper order. */

    sq = *sizeq;
    v = q[sq];  
    iv = iq[sq];
L10:
    tsz = sq / 2;
    if (tsz != 0) {     
        if (q[tsz] >= v) {     
            q[sq] = q[tsz];
            iq[sq] = iq[tsz];
            sq = tsz;
            goto L10;
        }       
    }       
    q[sq] = v;
    iq[sq] = iv;
}

/* ------------------------------------------------------------------------ */
/*  gis_removq(dimq,iq,iv,q,sizeq,v)                                        */
/*                                                                          */
/*  Remove smallest element (v,iv) from a priority queue (q,iq).            */

void gis_removq(TDAContext *ctx, int dimq,int *iq,int *iv,int *q,int *sizeq,int *v)
{
    (void)ctx; (void)dimq;        /* unused: the signature is shared */
    register int j,k;
    int ivtmp,szqd2,vtmp;

    /* Remove element from heap. */

    *v = q[1];
    *iv = iq[1];
    q[1] = q[*sizeq];
    iq[1] = iq[*sizeq];
    *sizeq -= 1;     

    /* Update heap to proper order. */

    k = 1;
    vtmp = q[k];
    ivtmp = iq[k];
    szqd2 = *sizeq / 2;
L10:
    if (k <= szqd2) {     
        j = k + k;
        if (j < *sizeq) {     
            if (q[j] > q[j+1])
                j++;  
        }        
        if (vtmp > q[j]) {     
            q[k] = q[j];
            iq[k] = iq[j];
            k = j;
            goto L10;
        }       
    }      
    q[k] = vtmp;
    iq[k] = ivtmp;
}

/* ------------------------------------------------------------------------ */
/*  gis_cpi4(dimn,n,x,y)                                                    */
/*                                                                          */
/*  Copy integer array x of size n into y.  Dimension of array is dimn.     */

void gis_cpi4(TDAContext *ctx, int dimn,int n,int *x,int *y)
{
    (void)ctx; (void)dimn;        /* unused: the signature is shared */
    register int i;

    for (i = 1; i <= n; ++i)
        y[i] = x[i];
}

/* ------------------------------------------------------------------------ */    
/*  (double)gis_randp(ix)                                                   */
/*                                                                          */
/*  Portable pseudo-random number generator.                                */
/*  Reference: L. Schrage, "A More Portable Fortran                         */
/*  Random Number Generator", ACM Transactions on                           */
/*  Mathematical Software, Vol. 2, No. 2, (June, 1979).                     */
/*                                                                          */

double gis_randp(TDAContext *ctx, int *ix)
{
    (void)ctx;        /* unused: the signature is shared */
    int xhi,xalo,leftlo,fhi,k;
    double randp;

    static int a = 16807;
    static int b15 = 32768;
    static int b16 = 65536;
    static int p = 2147483647;

    xhi = *ix / b16;
    xalo = (*ix - xhi * b16) * a;
    leftlo = xalo / b16;
    fhi = xhi * a + leftlo;
    k = fhi / b15;
    *ix = (((xalo - leftlo * b16) - p) + (fhi - k * b15) * b16) + k;
    if (*ix < 0)
        *ix += p;
    randp=(double)(*ix) * 4.656612875e-10;
    return(randp);
}

/* ------------------------------------------------------------------------ */
/*  gcni        Comparing Neighborhoods.                                    */
/*              Graph must be undirected unvalued, defined with GD_TYP 1    */
/*                                                                          */
/*              gcni(                                                       */
/*                  opt=...,        option, def. 1                          */
/*                                  1 = compares size of neighborhoods      */
/*                                  2 = compares size and n of edges        */
/*                                  3 = check for isomorphic subgraphs      */
/*                                      i,j fixed                           */
/*                                  4 = check for isomorphic subgraphs      */
/*                                      i,j not fixed                       */
/*                  prn=...,        format of output file, def. 1           */
/*                                  1 = one record for each class           */  
/*                                  2 = one record for each node            */  
/*                  n=...,          depth of neighborhood, def. 1           */  
/*                  gn =            graph number,                           */
/*                  df = ...        list of neighborhoods                   */
/*                  nfmt = ...,     integer format, def. 4                  */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcni(TDAContext *ctx)
{
    register int i,j,k;
    int err,nr,nr1,nr2,n1,n2,r,lev,nc,ne;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Comparing neighborhoods. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto GCNIFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCNIFin;

    if (gdd_tcheck(ctx, 1,1,1))      /* must be undirected, unvalued, GD_TYP 1 */
        goto GCNIFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMN < 1)
        ctx->PMN = 1;

    if (ctx->PMPRNO != 2)
        ctx->PMPRNO = 1;
    if (ctx->PMOPT > 4)
        ctx->PMOPT = 4;

    printf1(ctx, "Option %d: ",ctx->PMOPT);
    if (ctx->PMOPT == 1)
        printf1(ctx, "check size of neighborhoods.\n");
    else if (ctx->PMOPT == 2)
        printf1(ctx, "check size of neighborhoods and number of edges.\n");
    else if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
        printf1(ctx, "check for isomorphic neighborhoods (starting points ");
        if (ctx->PMOPT == 4)
            printf1(ctx, "not ");
        printf1(ctx, "fixed).\n");
    }
    printf1(ctx, "Depth of neighborhood: %d\n\n",ctx->PMN);

    if (alloc_acn(ctx, ctx->GD_NP))           
        goto GCNIFin;
    if (alloc_acm(ctx, ctx->GD_NP))           
        goto GCNIFin;
    if (alloc_ack(ctx, ctx->GD_NP))       /* level, equiv. classes */
        goto GCNIFin;
    if (alloc_acr(ctx, ctx->GD_NP + 1))   /* number of elements in classes */
        goto GCNIFin;
    if (alloc_acs(ctx, ctx->GD_NP + 1))   /* size of neighborhood */
        goto GCNIFin;
    if (alloc_aci(ctx, ctx->GD_NP + 1))   /* number of edges */
        goto GCNIFin;

    lev = nc = nr = nr1 = nr2 = 0;
    for (i = 0; i < ctx->GD_NP; ++i) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 50) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }

        if (alloc_acc(ctx, ctx->GD_NP))   
            goto GCNIFin;

        n1 = g_cni_n(ctx, i,ctx->PMN,ctx->AcN,ctx->AcC,ctx->PMGN);
        ne = g_cni_ne(ctx, n1,ctx->AcN,ctx->PMGN);

        if (ctx->PMProtFDef) {      /* write neighborhoods */
            g_cni_prn(ctx, i,ne,n1,ctx->AcN);
            nr2++;
        }
        r = 0;
        for (j = 0; j < i; ++j) {

            if (alloc_acc(ctx, ctx->GD_NP))   
                goto GCNIFin;

            n2 = g_cni_n(ctx, j,ctx->PMN,ctx->AcM,ctx->AcC,ctx->PMGN);

            r = g_cni_c(ctx, n1,ctx->AcN,n2,ctx->AcM,ne,ctx->PMOPT,ctx->PMGN);
            if (r) {
                if (r < 0)
                    goto GCNIFin;
                break;
            }
        }
        if (r == 0) {
            ctx->AcK[i] = ++lev;
            ctx->AcR[lev] = 1;
            ctx->AcS[lev] = n1;
            ctx->AcI[lev] = ne;
            nc++; 
        }
        else {
            ctx->AcK[i] = ctx->AcK[j];
            ctx->AcR[ctx->AcK[j]] += 1;
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 50) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }

    printf1(ctx, "Number of classes: %d\n",nc);

    for (i = 1; i <= lev; ++i) {
        k = 0;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcK[j] != i)
                continue;

            if (k == 0 || ctx->PMPRNO == 2) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,i);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcS[i]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcI[i]);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcR[i]);
            }
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));
            if (ctx->PMPRNO == 2) {
                fprintf(ctx->PMFd,"\n");
                nr++;
            }
            k++;
        }
        if (ctx->PMPRNO == 1) {
            fprintf(ctx->PMFd,"\n");
            nr++;
        }
    }
    printf1(ctx, "%d records written to: %s\n",nr,ctx->PMFdName);
    if (nr1 > 0)
        printf1(ctx, "%d records written to: %s\n",nr1,ctx->PMF1dName);
    if (nr2 > 0)
        printf1(ctx, "%d records written to: %s\n",nr2,ctx->PMProtFName);
    err = 0;

GCNIFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cni_n(i,d,nodes)  Return d-neighborhood of i in nodes[].              */
/*                      Return number of nodes in nodes{}.                  */

int g_cni_n(TDAContext *ctx, int i,int d,int *nodes,char *nflag,int gn)
{
    register int k,j,l,jn,jnn;
    int n,n1,n2,m,j1,j2;

    nodes[0] = i;
    nflag[i] = 1;
    n1 = 0;
    n = n2 = 1;
    for (k = 1; k <= d; ++k) {
        for (j = n1; j < n2; ++j) {
            jn = nodes[j];
            m = ctx->GD_FPN[jn];
            for (l = 0; l < m; ++l) {
                jnn = ctx->GD_FPI[jn][l];             /* edge from jn to jnn */
                if (ctx->AcC[jnn] == 0) {
                    j1 = ctx->GD_FPK[jn][l];
                    j2 = ctx->GD_FPK1[jn][l];
                    if ((j1 >= 0 && gdd_ev(ctx, j1,gn) >= 0.0) ||
                        (j2 >= 0 && gdd_ev(ctx, j2,gn) >= 0.0)) {  
                        nodes[n++] = jnn;
                        ctx->AcC[jnn] = 1;
                    }
                }
            }
        }
        n1 = n2;
        n2 = n;
    }
    return(n);
}

/* ------------------------------------------------------------------------ */
/*  g_cni_prn(i,ne,n,nodes)    pint neighborhood of i to PMProtFd           */

void g_cni_prn(TDAContext *ctx, int i,int ne,int n,int *nodes)
{
    register int j;

    rt_fprintf_i(ctx, ctx->PMProtFd,ctx->PMNFmtS,i + 1);
    rt_fprintf_i(ctx, ctx->PMProtFd,ctx->PMNFmtS,gdd_node(ctx, i));
    rt_fprintf_i(ctx, ctx->PMProtFd,ctx->PMNFmtS,ne);
    rt_fprintf_i(ctx, ctx->PMProtFd,ctx->PMNFmtS,n);
    for (j = 0; j < n; ++j)
        rt_fprintf_i(ctx, ctx->PMProtFd,ctx->PMNFmtS,gdd_node(ctx, nodes[j]));
    fprintf(ctx->PMProtFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  g_cni_ne(n,nodes)   Return number of edges in subgraph nodes[].         */

int g_cni_ne(TDAContext *ctx, int n,int *nodes,int gn)
{
    register int i,j,in,jn;
    int ne = 0;

    for (i = 1; i < n; ++i) {
        in = nodes[i];
        for (j = 0; j < i; ++j) {
            jn = nodes[j];
            if (gdd_adj(ctx, in,jn,gn) >= 0.0)   
                ne++;
        }
    }
    return(ne);
}

/* ------------------------------------------------------------------------ */
/*  g_cni_c(n1,nodes1,n2,nodes2,ne1,opt,gn)                                 */
/*                                                                          */
/*  Compare node sets nodes1[] and nodes2[].                                */
/*  ne1 is number of edges in first subgraph.                               */
/*                                                                          */
/*  opt = 1:    return 1 if same size, otherwise 0.                         */
/*  opt = 2:    return 1 if same size and same number of edges, else 0      */
/*              take nodes1[0] and nodes2[0] as fixed.                      */
/*  opt = 3:    return 1 if nodes1[] and nodes2[] isomorphic.               */
/*              take nodes1[0] and nodes2[0] as fixed.                      */
/*  opt = 4:    return 1 if nodes1[] and nodes2[] isomorphic.               */
/*              take nodes1[0] and nodes2[0] NOT as fixed.                  */
/*                                                                          */
/*  Return -1 if insufficient memory.                                       */

int g_cni_c(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int ne1,int opt,int gn)
{
    int ne2,r;

    /*********
    printf1(ctx, "\nn1=%d nodes1: ",n1);
    for (i = 0; i < n1; ++i)
        printf1(ctx, "%d ",nodes1[i]);
    newline(ctx);
       
    printf1(ctx, "n2=%d nodes2: ",n2);
    for (i = 0; i < n2; ++i)
        printf1(ctx, "%d ",nodes2[i]);
    newline(ctx);
    ****/

    if (n1 != n2)           /* compare number of nodes */
        return(0);
    if (opt == 1)
        return(1);

    /* compare number of edges */

    ne2 = g_cni_ne(ctx, n2,nodes2,gn);

    if (ne1 != ne2)
        return(0);
    if (opt == 2)
        return(1);

    r = 0;
    if (opt == 3)
        r = g_cni_ci(ctx, n1,nodes1,n2,nodes2,gn);
    else if (opt == 4)
        r = g_cni_ci1(ctx, n1,nodes1,n2,nodes2,gn);
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  g_cni_ci(n1,nodes1,n2,nodes2,gn)                                        */
/*                                                                          */
/*  Compare node sets nodes1[] and nodes2[].                                */
/*  Return 1 if isomorphic, otherwise 0.                                    */
/*  Take nodes1[0] and nodes2[0] as fixed.                                  */
/*  Use simple enumerative procedure.                                       */
/*  Return -1 if insufficient memory.                                       */

int g_cni_ci(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int gn)
{
    register int i,j,in,jn,inn,jnn;
    int n,first,r,is;
    int *com,*p,*d;

    if (n1 != n2)
        return(0);

    n = n1 - 1;
    if (n < 1)
        return(1);

    if (!(com = (int *)calloc((size_t)(n),sizeof(int))))   
        return(-1);
    if (!(p = (int *)calloc((size_t)(n),sizeof(int)))) { 
        free((char *)com);
        return(-1);
    }
    if (!(d = (int *)calloc((size_t)(n),sizeof(int)))) { 
        free((char *)com);
        free((char *)p);        /* p is the live allocation here; the
                                   original freed d, which is the one
                                   that just failed (cppcheck) */
        return(-1);
    }
    memrq(ctx, 3 * n,sizeof(int));

    r = 0;
    first = 1;
    while (perm(ctx, n,com,p,d,first)) {
        first = 0;
        is = 1;
        for (i = 1; i < n1; ++i) {
            in = nodes1[i];
            inn = nodes2[com[i - 1] + 1];

            for (j = 0; j < i; ++j) {
                jn = nodes1[j];
                if (j == 0)
                    jnn = nodes2[0];
                else
                    jnn = nodes2[com[j - 1] + 1];

                if (gdd_adj(ctx, in,jn,gn) >= 0.0) {
                    if (gdd_adj(ctx, inn,jnn,gn) < 0.0) {
                        is = 0;
                        break;
                    }
                }
                else {
                    if (gdd_adj(ctx, inn,jnn,gn) >= 0.0) {
                        is = 0;
                        break;
                    }
                }
            }
            if (is == 0)
                break;
        }
        if (is) {
            r = 1;
            break;
        }
    }
    free((char *)com);
    free((char *)p);
    free((char *)d);
    memrq(ctx, -3 * n,sizeof(int));
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  g_cni_ci1(n1,nodes1,n2,nodes2,gn)                                       */
/*                                                                          */
/*  Compare node sets nodes1[] and nodes2[].                                */
/*  Return 1 if isomorphic, otherwise 0.                                    */
/*  Take nodes1[0] and nodes2[0] NOT as fixed.                              */
/*  Use simple enumerative procedure.                                       */
/*  Return -1 if insufficient memory.                                       */

int g_cni_ci1(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int gn)
{
    register int i,j,in,jn,inn,jnn;
    int n,first,r,is;
    int *com,*p,*d;

    if (n1 != n2)
        return(0);

    n = n1;
    if (n <= 1)
        return(1);

    if (!(com = (int *)calloc((size_t)(n),sizeof(int))))   
        return(-1);
    if (!(p = (int *)calloc((size_t)(n),sizeof(int)))) { 
        free((char *)com);
        return(-1);
    }
    if (!(d = (int *)calloc((size_t)(n),sizeof(int)))) { 
        free((char *)com);
        free((char *)p);        /* p is the live allocation here; the
                                   original freed d, which is the one
                                   that just failed (cppcheck) */
        return(-1);
    }
    memrq(ctx, 3 * n,sizeof(int));

    r = 0;
    first = 1;
    while (perm(ctx, n,com,p,d,first)) {
        first = 0;
        is = 1;
        for (i = 1; i < n1; ++i) {
            in = nodes1[i];
            inn = nodes2[com[i]];

            for (j = 0; j < i; ++j) {
                jn = nodes1[j];
                jnn = nodes2[com[j]];

                if (gdd_adj(ctx, in,jn,gn) >= 0.0) {
                    if (gdd_adj(ctx, inn,jnn,gn) < 0.0) {
                        is = 0;
                        break;
                    }
                }
                else {
                    if (gdd_adj(ctx, inn,jnn,gn) >= 0.0) {
                        is = 0;
                        break;
                    }
                }
            }
            if (is == 0)
                break;
        }
        if (is) {
            r = 1;
            break;
        }
    }
    free((char *)com);
    free((char *)p);
    free((char *)d);
    memrq(ctx, -3 * n,sizeof(int));
    return(r);
}


/* -##--------------------------------------------------------------------- */
/*  gcset       Compact sets.       Graph is always undirected/valued.      */
/*                                  Loops are ignored.                      */
/*                                  Always uses Kruskal alg. for MST        */

/*              gcset(                                                      */
/*                  opt=...,        output option, def. 1                   */
/*                                  1 = one record for each compact set     */
/*                                  2 = one record for each node            */
/*                  gn =            graph number,                           */
/*                  df = ...        pseudo image graph as edge list         */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  fmt = ...,      print format for values, def. 10.4      */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int gcset(TDAContext *ctx)
{
    register int i,j;
    int err,r,ns,nr,drec,nrec,ne,nt,nte;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Compact sets. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GCSETFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCSETFin;

    if (gdd_tcheck(ctx, 0,1,2))      /* must be undirected */
        goto GCSETFin;

    if (ctx->PMFmtF == 0)            /* default print format */
        pmfmt(ctx, 10,4);

    if (ctx->PMOPT > 2)
        ctx->PMOPT = 2;

    if (alloc_acn(ctx, ctx->GD_NP))                 
        goto GCSETFin;
    if (alloc_acm(ctx, ctx->GD_NP))                 
        goto GCSETFin;
    if (alloc_aci(ctx, ctx->GD_NP))                 
        goto GCSETFin;
    if (alloc_acj(ctx, ctx->GD_NP))                 
        goto GCSETFin;

    if (ctx->PMF1Def) {      /* extra storage for pseudo image graph */
        if (g_gclptr_i(ctx, ctx->GD_NP + 1))
            goto GCSETFin;
    }

    ctx->GSCON_ID = 1;
    drec = ns = nte = nt = i = nrec = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
            printfe(ctx, "Component: %7d     %c",ctx->GSCON_ID,CR);
            fflushe(ctx);
        }
        ctx->GSCON_NN = 0;
        if (ctx->GD_TYP == 1)
            u2_visit(ctx, i,ctx->PMGN);
        else
            u2_visit1(ctx, i,ctx->PMGN);
     
        if (ctx->GSCON_NN == 0)
            break;

        if (ctx->GSCON_NN == 1) {
            ctx->AcI[0] = ctx->AcJ[0] = ctx->AcM[0];
            ne = 1;
            nte++;
        }
        else {

            /* note, this algorithm allocates AcX,AcR,AcS */

            ne = g_mst1(ctx, ctx->PMGN,ctx->GD_TYP,ctx->GSCON_NN,ctx->AcM,ctx->AcI,ctx->AcJ,0);
            if (ne <= 0) {
                if (ne == 0)  
                    printf1(ctx, "ERROR: gmst.\n");
                goto GCSETFin;
            }   
        }
        r = g_cset(ctx, ne,ctx->AcI,ctx->AcJ,ctx->GSCON_ID,ctx->PMOPT,&nr,ctx->PMF1Def);
        if (r < 0)
            goto GCSETFin;

        ns += r;
        nrec += nr;
        nt++;
        ctx->GSCON_ID++;  

        if (ctx->PMF1Def) {      /* create pseudo image graph */
            r = g_cset_df(ctx, ctx->GSCON_NN,ctx->AcM);    
            if (r < 0)
                goto GCSETFin;
            drec += r;
        }

        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == 0) {
                i = j;
                break;
            }
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP > 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "\nNumber of components: %d",nt);
    if (nte > 0)
        printf1(ctx, " [includes %d isolated nodes]",nte);
    printf1(ctx, "\nNumber of compact sets: %d\n",ns);
    printf1(ctx, "%d records written to: %s\n",nrec,ctx->PMFdName);
    if (drec > 0)
        printf1(ctx, "%d records written to: %s\n",drec,ctx->PMF1dName);
    err = 0;

GCSETFin:       
    g_gclptr_i(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  g_cset(m,in,jn,cn,opt,nr)                                               */
/*                                                                          */
/*  Find compact sets for current tree. For the algorithm see               */
/*  Chiou-Kuo Liang, An O(n2) algorithm for finding the compact sets        */
/*  of a graph, BIT 33 (1993), 390 - 395.                                   */
/*                                                                          */
/*  m = number of edges.                                                    */
/*  edges are (in[j],jn[j]), j = 0,...,m - 1.                               */
/*  cn is number of current conn. component                                 */
/*  opt controls writing to output file                                     */  
/*  return in nr the number of records written to output file.              */
/*  if dflag != 0 add csets to global list.                                 */
/*                                                                          */  
/*  Return: number of compact sets, or -1 if error.                         */

int g_cset(TDAContext *ctx, int m,int *in,int *jn,int cn,int opt,int *nr,int dflag)
{
    (void)cn;        /* unused: the signature is shared */
    register int i,j,k,l,ii;
    int err,n,i0,j0,lev,li,lj,ns;
    double tmp,dia,dm;

    err = -1;
    ns = 0;
    *nr = 0;

    if (alloc_ack(ctx, m))
        goto G_CSETFin;
    if (alloc_acx(ctx, m))
        goto G_CSETFin;
    if (alloc_acr(ctx, ctx->GD_NP + 1))
        goto G_CSETFin;
    if (alloc_acs(ctx, ctx->GD_NP + 1))
        goto G_CSETFin;
    if (alloc_acu(ctx, ctx->GD_NP + 1))
        goto G_CSETFin;

    for (j = 0; j < m; ++j)  
        ctx->AcX[j] = gdd_adj(ctx, in[j],jn[j],ctx->PMGN);

    if (sortdp(ctx, m,ctx->AcX,ctx->AcK))
        goto G_CSETFin;

    lev = 0;
    for (j = 0; j < m; ++j) {
        k = ctx->AcK[j];
        i0 = in[k] + 1;
        j0 = jn[k] + 1;
        li = ctx->AcS[i0];
        lj = ctx->AcS[j0];

        if (li == 0 && lj == 0) {       /* new level */
            l = ++lev;

            ctx->AcS[i0] = ctx->AcS[j0] = l;
            dia = ctx->AcU[l] = ctx->AcX[k];          
        }
        else if (li == 0) {
            l = lj;
            ctx->AcS[i0] = l;      

            /* update diameter */
            dia = dmax(ctx, ctx->AcU[l],ctx->AcX[k]);            

            for (i = 1; i <= ctx->GD_NP; ++i) {
                if (i == in[k] + 1)
                    continue;

                if (ctx->AcS[i] == l && (tmp = gdd_adj(ctx, in[k],i - 1,ctx->PMGN)) > dia)  
                    dia = tmp;
            }
            ctx->AcU[l] = dia;                            
        }
        else if (lj == 0) {
            l = li;
            ctx->AcS[j0] = l;      

            /* update diameter */
            dia = dmax(ctx, ctx->AcU[l],ctx->AcX[k]);            

            for (i = 1; i <= ctx->GD_NP; ++i) {
                if (i == jn[k] + 1)
                    continue;

                if (ctx->AcS[i] == l && (tmp = gdd_adj(ctx, jn[k],i - 1,ctx->PMGN)) > dia)  
                    dia = tmp;
            }
            ctx->AcU[l] = dia;                            
        }
        else {
            l = imax(ctx, li,lj);

            /* update diameter */
            dia = dmax(ctx, ctx->AcU[l],ctx->AcX[k]);            

            for (i = 1; i <= ctx->GD_NP; ++i) {
                if (ctx->AcS[i] != li)
                    continue;
                ctx->AcS[i] = l;                          

                for (ii = 1; ii <= ctx->GD_NP; ++ii) {
                    if (ctx->AcS[ii] != lj || ii == i)
                        continue;
                    ctx->AcS[ii] = l;

                    if ((tmp = gdd_adj(ctx, i - 1,ii - 1,ctx->PMGN)) > dia)  
                        dia = tmp;
                }
            }
            ctx->AcU[l] = dia;                            
        }

        /* find min distance to group complement */

        dm = ctx->DBLMAX;
        n = 0;
        for (i = 0; i < m; ++i) {
            i0 = in[i] + 1;
            j0 = jn[i] + 1;
            if ((ctx->AcS[i0] == l && ctx->AcS[j0] != l) || (ctx->AcS[i0] != l && ctx->AcS[j0] == l)) {    
                dm = dmin(ctx, dm,ctx->AcX[i]);
                n++;
            }
        }
        if (n > 0 && dia < dm) {        /* if compact set */
            ns++;
            n = 0;
            for (i = 1; i <= ctx->GD_NP; ++i) {
                if (ctx->AcS[i] == l)
                    ctx->AcR[n++] = i - 1;
            }   

            if (dflag) {                    /* add to global list */
                if (g_gclptr_a(ctx, n,ctx->AcR))
                    goto G_CSETFin;
            }
            if (opt == 1) {
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ns);
                rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dia);
                rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dm);
                for (i = 0; i < n; ++i) 
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcR[i]));
                fprintf(ctx->PMFd,"\n");
                *nr += 1;
            }
            else if (opt == 2) {
                for (i = 0; i < n; ++i) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GSCON_ID);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ns);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dia);
                    rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,dm);
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, ctx->AcR[i]));
                    fprintf(ctx->PMFd,"\n");
                    *nr += 1;
                }
            }
        }
    }
    err = ns;

G_CSETFin:
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  g_cset_df(n,nodes)                                                      */
/*                                                                          */

int g_cset_df(TDAContext *ctx, int n,int *nodes)
{
    register int i,j;
    int nrec,ns,np,n1,n2,in,jn,it,jt,*lptr;
    double val;

    nrec = 0;
    ns = 0;     /* create list of maximal csets */

tda_out("GCL_PN=%d\n",ctx->GCL_PN);
tda_out("n=%d nodes: ",n);
for (i = 0; i < n; ++i)
tda_out("%d ",nodes[i]);
newline(ctx);


    if (alloc_acr(ctx, ctx->GCL_PN))                 
        return(-1);     

    for (i = 0; i < ctx->GCL_PN; ++i) {

        np = ctx->GCL_NPTR[i];

        if (g_gclptr_cc(ctx, np,ctx->GCL_PTR[i],i) == 0)  
            ctx->AcR[ns++] = i;
    }

    /* create AcS[i] = 0 for isolated node, ptr+1 for starting node of cset,
       -1 for other members of cset */

    if (alloc_acs(ctx, ctx->GD_NP))                 
        return(-1);     

    for (i = 0; i < ns; ++i) {
        np = ctx->GCL_NPTR[ctx->AcR[i]];
        lptr = ctx->GCL_PTR[ctx->AcR[i]];
        ctx->AcS[lptr[0]] = ctx->AcR[i] + 1;
        for (j = 1; j < np; ++j)
            ctx->AcS[lptr[j]] = -1;
    }

    /* create pseudo image graph */

    for (i = 0; i < n; ++i) {
        in = nodes[i];
        it = ctx->AcS[in];
        if (it < 0)
            continue;

        for (j = i + 1; j < n; ++j) {
            jn = nodes[j];
            jt = ctx->AcS[jn];
            if (jt < 0)
                continue;

tda_out("j=%d n=%d\n",j,n);

            /* calculate new edge values */
     
            n1 = n2 = 1;
            if (it == 0 && jt == 0)  
                val = gdd_adj(ctx, in,jn,ctx->PMGN);
            else if (it == 0) {
                n2 = ctx->GCL_NPTR[jt - 1];
                val = g_mineval(ctx, 1,&in,n2,ctx->GCL_PTR[jt-1],0);
            }
            else if (jt == 0) {
                n1 = ctx->GCL_NPTR[it - 1];
                val = g_mineval(ctx, n1,ctx->GCL_PTR[it-1],1,&jn,0);
            }
            else {              
                n1 = ctx->GCL_NPTR[it - 1];
                n2 = ctx->GCL_NPTR[jt - 1];
                val = g_mineval(ctx, n1,ctx->GCL_PTR[it-1],n2,ctx->GCL_PTR[jt-1],0);
            }
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,in + 1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,jn + 1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, in));
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, jn));
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,n1);
            rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,n2);
            rt_fprintf_d(ctx, ctx->PMF1d,ctx->PMFmtS,val);
            fprintf(ctx->PMF1d,"\n");
            nrec++;
        }
    }
    return(nrec);
}

/* ------------------------------------------------------------------------ */
/*  g_mineval(n1,nodes1,n2,nodes,opt)                                       */
/*                                                                          */
/*  If opt = 0, return min edge values for edges v(i,j), in in nodes1,      */
/*  j in nodes2. If opt != 0 return max edge value.                         */
/*  Skip loops.                                                             */

double g_mineval(TDAContext *ctx, int n1,int *nodes1,int n2,int *nodes2,int opt)
{
    register int i,j,in,jn;
    double val,tmp;

    if (opt == 0)
        val = ctx->DBLMAX;
    else
        val = -1.0;

    for (i = 0; i < n1; ++i) {
        in = nodes1[i];
        for (j = 0; j < n2; ++j) {
            jn = nodes2[j];
            if (in == jn)
                continue;
            tmp = gdd_adj(ctx, in,jn,ctx->PMGN);
            if (tmp >= 0.0) {
                if (opt == 0)
                    val = dmin(ctx, val,tmp);
                else
                    val = dmax(ctx, val,tmp);
            }
        }
    }
    if (val >= ctx->DBLMAX)
        val = -1.0;
    return(val);
} 

/* ------------------------------------------------------------------------ */
/*  gcliq       Cliques. Graph is always undirected and unvalued.           */  
/*                       Opt 1 requires gdd option 1.                       */
/*                                                                          */
/*              gcliq(                                                      */
/*                  alg=...,        1 = Bron-Kerbosch algorithm             */
/*                                  2 = Harary-Ross algorithm               */
/*                  min=...,        min elements of cliques, def. 3         */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  sort,           sort node numbers in cliques            */  
/*                  df=...,         second output file                      */  
/*              ) = fname;          output file (required)                  */
/*                                                                          */
/*  The df option can be used to request a second output file organized     */
/*  as a node list. Each record will contain three entries:                 */
/*  column 1 : node number                                                  */
/*  column 2 : number of component the node belongs to                      */
/*  column 3 : number of clique the node belongs to                         */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int gcliq(TDAContext *ctx)
{
    register int i,j;
    int err,nc,ni,ncm,ncl;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Cliques. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 5,1,1))       /* get parameters */
        goto GCLFin;

    if (ctx->PMALG != 2)
        ctx->PMALG = 1;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GCLFin;

    printf1(ctx, "Algorithm: ");
    if (ctx->PMALG == 1) {
        printf1(ctx, "Bron-Kerbosch.\n");
        if (gdd_tcheck(ctx, 0,1,0))      /* must be undirected */
            goto GCLFin;
        if (ctx->PMMin < 2)
            ctx->PMMin = 3;
    }
    else {
        printf1(ctx, "Harary-Ross.\n");
        if (gdd_tcheck(ctx, 1,1,0))      /* must be undirected, gdd option 1 */
            goto GCLFin;
        if (ctx->PMMin < 3)
            ctx->PMMin = 3;
    }
    if (alloc_acn(ctx, ctx->GD_NP))       /* used in c_visit */
        goto GCLFin;
    if (alloc_acm(ctx, ctx->GD_NP))       /* used in c_visit */
        goto GCLFin;

    ctx->GCL_CNT1 = 0;               /* unique clique id across components */
    ctx->GCL_M = 0;                  /* max clique size */
    ctx->GG_REC = 0;                 /* number of records for output file */
    ctx->GG_ID = 1;
    ncl = ni = ncm = nc = i = 0;

    while (i >= 0) {

        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
            printfe(ctx, "Component: %7d     %c",ctx->GG_ID,CR);
            fflushe(ctx);
        }
        ctx->GG_CNT = 0;      /* number of nodes in current component */
        if (ctx->GD_TYP == 1)
            c_visit(ctx, i,ctx->PMGN);
        else
            c_visit1(ctx, i,ctx->PMGN);
     
        nc++;                       /* number of components */

        if (ctx->GG_CNT < ctx->PMMin) {
            if (ctx->GG_CNT == 1)     
                ni++;
        }
        else {
            ncm++;          /* components with min size 3 */
      
            ctx->GCL_CNT = 0;  /* number of cliques in current component */

            if (ctx->PMALG == 1) {
                if (g_cliq(ctx, nc,ctx->GG_CNT,ctx->AcM,ctx->PMMin))
                    goto GCLFin;
            }
            else {
                if (g_cliq_2(ctx, nc,ctx->GG_CNT,ctx->AcM,ctx->PMMin))
                    goto GCLFin;
            }
            ncl += ctx->GCL_CNT;
        }
        ctx->GG_ID++;  

        i = -1;
        for (j = 0; j < ctx->GD_NP; ++j) {
            if (ctx->AcN[j] == 0) {
                i = j;
                break;
            }
        }
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 500) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    printf1(ctx, "Number of components: %d",nc);
    if (ni > 0)
        printf1(ctx, " [includes %d isolated nodes]",ni);
    printf1(ctx, "\nComponents with at least %d members: %d\n",ctx->PMMin,ncm);
    printf1(ctx, "Cliques with at least %d members: %d\n",ctx->PMMin,ncl);
    printf1(ctx, "Maximal clique size: %d\n",ctx->GCL_M);
                
    printf1(ctx, "%d records written to: %s\n",ctx->GG_REC,ctx->PMFdName);
    err = 0;

GCLFin:       
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  c_visit(k,gn)    called by gcon. k is node number, gn is graph number.   */

void c_visit(TDAContext *ctx, int k,int gn)
{
    register int j,l,k1,k2;
    int n;

    ctx->AcN[k] = ctx->GG_ID;
    ctx->AcM[ctx->GG_CNT++] = k;

    n = ctx->GD_FPN[k];
    if (n <= 0)
        return;          

    for (l = 0; l < n; ++l) {

        j = ctx->GD_FPI[k][l];                    /* edge from k to j */
        if (ctx->AcN[j] == 0) {
            k1 = ctx->GD_FPK[k][l];
            k2 = ctx->GD_FPK1[k][l];
            if ((k1 >= 0 && gdd_ev(ctx, k1,gn) >= 0.0) ||
                (k2 >= 0 && gdd_ev(ctx, k2,gn) >= 0.0)) {  
                c_visit(ctx, j,gn);
            }
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  c_visit1(k,gn)   called by gcon. k is node number, gn is graph number.   */

void c_visit1(TDAContext *ctx, int k,int gn)
{
    register int j;

    ctx->AcN[k] = ctx->GG_ID;
    ctx->AcM[ctx->GG_CNT++] = k;

    for (j = 0; j < ctx->GD_NP; ++j) {

        if (ctx->AcN[j] == 0) {
            if (gdd_adj(ctx, k,j,gn) >= 0.0)  
                c_visit1(ctx, j,gn);
        }
    }
}

/* ------------------------------------------------------------------------ */
/*  g_cliq(nc,n,nodes,min)                                                  */
/*                                                                          */
/*  Calculate cliques of an undirected unvalued graph. Algorithm adapted    */
/*  from: C. Bron, J. Kerbosch, Algorithm 457 -- Finding all cliques of     */
/*  an undirected graph, Comm. ACM 16, 1973, 575 - 577.                     */
/*                                                                          */
/*  nodes[i] (i = 0,...,n-1) contains the node numbers of current           */    
/*  component. nc is number of current component.                           */
/*  min is min number of elements in cliques.                               */
/*                                                                          */
/*  Record maximal clique size in GCL_M.                                    */
/*  Count number of cliques in GCL_CNT and GCL_CNT1                         */
/*  Return 0 if OK, or -1 if insuff. memory.                                */

int g_cliq(TDAContext *ctx, int nc,int n,int *nodes,int min)
{
    register int i;
    int err;
         
    err = -1;
    if (alloc_ack(ctx, n + 1))       /* working space in g_cliq */
        return(-1);     
    if (alloc_aci(ctx, n + 1))       /* records node number in g_cliq */
        return(-1);     
    if (alloc_acj(ctx, n + 1))       /* mapping of node numbers for g_cliq */
        return(-1);     

    for (i = 1; i <= n; ++i) {
        ctx->AcK[i] = i;
        ctx->AcJ[i] = nodes[i - 1];
    }
    ctx->GCL_N   = 0;                /* number of nodes in cliques */

    if (g_cliq_1(ctx, ctx->AcK,0,n,ctx->AcJ,nc,min))
        goto G_CLIQFin;

    err = 0;        

G_CLIQFin:
    return(err);
}

/* -------------------------------------------------------------------- */
/*  The following function is recursively callled by g_clique()         */
/*                                                                      */
/*  old[i], i = ne + 1,...,ne contains node number                      */
/*  mmap[] maps local to internal node numbers.                         */
/*  AcI[] is used to record nodes in cliques, counter is GCL_N.         */
/*  Cliques are  written to PMFd. GG_REC is record count.               */
/*  Counter for cliques is: GCL_CNT (and GCL_CNT1)                      */
/*  nc is number of current component                                   */
/*  Print and count cliques with at least min members.                  */
/*  Record maximal clique size in GCL_M.                                */
/*                                                                      */
/*  Return 0 if OK, -1 if error (insuff memory).                        */

int g_cliq_1(TDAContext *ctx, int *old,int ne,int ce,int *nmap,int nc,int min)
{
    register int i = 0,j = 0,nod = 0,p = 0;
    int r = 0,minnod = 0,count = 0,fixp = 0,s = 0,sel = 0,pos = 0,newce = 0,newne = 0;
    int *nnew = NULL;

    if (!(nnew = (int *)calloc((size_t)(ce + 1),sizeof(int))))   
        return(-1);
    memrq(ctx, ce + 1,sizeof(int));

    minnod = ce;
    r = i = nod = 0;                                                    

    while (++i <= ce && minnod != 0) {
        p = old[i];
        count = 0;
        j = ne;

        while (++j <= ce && count < minnod) {
            if (p != old[j] &&
                gdd_adj(ctx, nmap[p],nmap[old[j]],ctx->PMGN) < 0.0) {
                count++;
                pos = j;
            }
        }
        if (count < minnod) {
            fixp = p;
            minnod = count;
            if (i <= ne)
                s = pos;
            else {
                s = i;
                nod = 1;
            }
        }
    }
    for (nod = minnod + nod; nod >= 1; --nod) {
        p = old[s];
        old[s] = old[ne + 1];
        sel = old[ne + 1] = p;

        newne = i = 0;
        while (++i <= ne) {
            if (sel == old[i] ||
                gdd_adj(ctx, nmap[sel],nmap[old[i]],ctx->PMGN) >= 0.0) {
                newne++;
                nnew[newne] = old[i];
            }
        }
        newce = newne;
        i = ne + 1;
        while (++i <= ce) {
            if (sel == old[i] ||
                gdd_adj(ctx, nmap[sel],nmap[old[i]],ctx->PMGN) >= 0.0) {
                newce++;
                nnew[newce] = old[i];
            }
        }
        ctx->AcI[ctx->GCL_N++] = sel;

        if (newce == 0 && ctx->GCL_N >= min) {
            ctx->GCL_CNT++;
            ctx->GCL_CNT1++;
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GCL_CNT);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GCL_N);

            if (ctx->PMSORTFlg) {
                if (alloc_acs(ctx, ctx->GCL_N)) {
                    r = -1;
                    goto GCL1Fin;     
                }   
                for (i = 0; i < ctx->GCL_N; ++i)  
                    ctx->AcS[i] = gdd_node(ctx, nmap[ctx->AcI[i]]);

                if (sorti(ctx, ctx->GCL_N,ctx->AcS,0)) {
                    r = -1;
                    goto GCL1Fin;     
                }   
                for (i = 0; i < ctx->GCL_N; ++i) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcS[i]);

                    if (ctx->PMF1Def) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcS[i]);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nc);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT1);
                        fprintf(ctx->PMF1d,"\n");
                    }
                }
            }
            else {
                for (i = 0; i < ctx->GCL_N; ++i) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nmap[ctx->AcI[i]]));

                    if (ctx->PMF1Def) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, nmap[ctx->AcI[i]]));
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nc);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT1);
                        fprintf(ctx->PMF1d,"\n");
                    }
                }
            }
            fprintf(ctx->PMFd,"\n");
            ctx->GG_REC++;
            if (ctx->GCL_M < ctx->GCL_N)
                ctx->GCL_M = ctx->GCL_N;
        }
        else {
            if (newne < newce) {
                r = g_cliq_1(ctx, nnew,newne,newce,nmap,nc,min);
                if (r)
                    goto GCL1Fin;
            }
        }
        ctx->GCL_N--;             
        ne++;
        if (nod > 1) {
            s = ne;
            while (s++ &&  old[s]) {
                if (fixp != old[s] &&
                    gdd_adj(ctx, nmap[fixp],nmap[old[s]],ctx->PMGN) < 0.0)
                    break;
            }
        }
    }
GCL1Fin:
    free((char *)nnew);
    memrq(ctx, -ce - 1,sizeof(int));
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  g_cliq_2(nc,n,nodes,min)     [recursive function]                       */
/*                                                                          */
/*  Find cliques with Harary-Ross algorithm.                                */
/*  nodes[i] (i = 0,...,n-1) contains the nodes of a connected component.   */
/*  nc is number of component.                                              */
/*  min is min number of elements in cliques.                               */
/*  Record maximal clique size in GCL_M.                                    */
/*  Count number of cliques in GCL_CNT and GCL_CNT1                         */
/*  Write cliques to PMFd, count records in GG_REC.                         */
/*  Return 0 if OK, or -1 if insufficient memory.                           */

int g_cliq_2(TDAContext *ctx, int nc,int n,int *nodes,int min)
{
    register int i,j,ii;
    int err,rs,ne,nu,nn,nmax,nflaga,nflag1a,nnda,im,rmin;
    int *nnd;
    char *nflag,*nflag1;
  
    err = -1;
    nn = n;
    nnda = nflaga = nflag1a = nmax = 0;
    for (i = 0; i < n; ++i)  
        nmax = imax(ctx, nmax,nodes[i]);

    if (!(nflag = (char *)calloc((size_t)(nmax + 1),sizeof(char))))   
        goto GCL2Fin;
    nflaga = nmax + 1;
    memrq(ctx, nflaga,sizeof(char));

    for (i = 0; i < n; ++i)   
        nflag[nodes[i]] = 2;

    if (n < 1)
        goto GCL2Fin;
    if (!(nnd = (int *)calloc((size_t)(n),sizeof(int))))   
        goto GCL2Fin;
    nnda = n;
    memrq(ctx, nnda,sizeof(int));

    while (1) {
        nu = 0;     /* set if we find unicliqual node */
        im = -1;
        for (i = 0; i < n; ++i) {   
            ii = nodes[i];
            if (nflag[ii] == 0)  
                continue;

            ne = g_cliq_2getd(ctx, ii,nmax,nflag,NULL,n,nnd,&rs);
            if (ne < 0)
                goto GCL2Fin;

            if (ne == 0) {                  /* not in any clique */
                if (nflag[ii]) {
                    nflag[ii] = 0;
                    nn--;
                }
                continue;
            }
            if (im < 0 || rmin > rs) {
                im = ii,
                rmin = rs;
            }   
            if (rs == ne * (ne - 1)) {      /*  unicliqual point */
                nu = 1;
                break;
            }
        }
        if (nu == 0)        /* no more unicliqual nodes */
            break;

        /* found unicliqual node */

        if (ne + 1 >= min) {

            ctx->GCL_CNT++;
            ctx->GCL_CNT1++;
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,nc);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GCL_CNT);
            rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ne + 1);

            if (ctx->PMSORTFlg) {
                if (alloc_acs(ctx, ne + 1))  
                    goto GCL2Fin;     

                for (i = 0; i < ne; ++i)  
                    ctx->AcS[i] = gdd_node(ctx, nnd[i]);
                ctx->AcS[ne] = gdd_node(ctx, ii);

                if (sorti(ctx, ne + 1,ctx->AcS,0))  
                    goto GCL2Fin;     

                for (i = 0; i <= ne; ++i) {
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->AcS[i]);

                    if (ctx->PMF1Def) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->AcS[i]);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nc);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT1);
                        fprintf(ctx->PMF1d,"\n");
                    }
                }
            }
            else {
                for (i = -1; i < ne; ++i) {
                    if (i < 0) j = ii;
                    else       j = nnd[i];
                    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, j));

                    if (ctx->PMF1Def) {
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,gdd_node(ctx, j));
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,nc);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT);
                        rt_fprintf_i(ctx, ctx->PMF1d,ctx->PMNFmtS,ctx->GCL_CNT1);
                        fprintf(ctx->PMF1d,"\n");
                    }
                }
            }
            fprintf(ctx->PMFd,"\n");
            ctx->GG_REC++;
            if (ctx->GCL_M < ne + 1)
                ctx->GCL_M = ne + 1;
        }
        if (nn == ne + 1) {         /* all found */
            err = 0;
            goto GCL2Fin;
        }

        /* remove all unicliqual nodes contained in current clique */

        for (i = 0; i < ne; ++i) {
            j = g_cliq_2getd(ctx, nnd[i],nmax,nflag,NULL,n,NULL,&rs);
            if (j < 0)
                goto GCL2Fin;

            if (rs == j * (j - 1))  
                nflag[nnd[i]] = 3;          
        }
        for (i = 0; i < ne; ++i) {
            if (nflag[nnd[i]] == 3) {
                nflag[nnd[i]] = 0;   
                nn--;
            }
        }
        nflag[ii] = 0;
        nn--;
    }
    if (im < 0) {
        err = 0;
        goto GCL2Fin;
    }
                
    /* create first group: all nodes cocliqual with im */

    ne = g_cliq_2getd(ctx, im,nmax,nflag,NULL,n,nnd,&rs);
    nnd[ne++] = im;

    /* create second group */
   
    for (i = 0; i < ne; ++i)  
        nflag[nnd[i]] = 1;
 
    if (!(nflag1 = (char *)calloc((size_t)(nmax + 1),sizeof(char))))   
        goto GCL2Fin;
    nflag1a = nmax + 1;     
    memrq(ctx, nflag1a,sizeof(char));
    
    for (i = 0; i <= nmax; ++i) {
        if (nflag[i] == 2) {
            nflag1[i] = 2;
            if (g_cliq_2getd(ctx, i,nmax,nflag,nflag1,n,NULL,&rs) < 0)
                goto GCL2Fin;
        }
    }
     
    /* process first group */

    if (ne >= min) {       
        if (g_cliq_2(ctx, nc,ne,nnd,min))
            goto GCL2Fin;
    }      

    /* process second group */

    ne = 0;
    for (i = 0; i <= nmax; ++i) {
        if (nflag1[i])       
            nnd[ne++] = i;
    }
    if (ne >= min) {      
        if (g_cliq_2(ctx, nc,ne,nnd,min))
            goto GCL2Fin;
    }      
    err = 0;

GCL2Fin:
    if (nflaga > 0) {
        free((char *)nflag);
        memrq(ctx, -nflaga,sizeof(char));
    }
    if (nflag1a > 0) {
        free((char *)nflag1);
        memrq(ctx, -nflag1a,sizeof(char));
    }
    if (nnda > 0) {
        free((char *)nnd);
        memrq(ctx, -nnda,sizeof(int));
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_cliq_2getd()      Get non-zero elements of i.th row of AA x A         */
/*                      into  d[j], j = 0,1,...                             */
/*                                                                          */
/*  Assume: gdd option 1.                                                   */
/*  Return row sum in rs.                                                   */
/*  if d != NULL: put node  number of nonzero elements of rows into d.      */
/*  if nflag1 != NULL: flag nonzero elements of rows in nflag1.             */
/*  Return number of nonzero elements in row.                               */

int g_cliq_2getd(TDAContext *ctx, int i,int nmax,char *nflag,char *nflag1,int n,int *d,int *rs)
{
    (void)n;        /* unused: the signature is shared */
    register int l,ll,j,k;
    int nn,cnt,ne;

    *rs = ne = 0;
    nn = ctx->GD_FPN[i];
    if (nn <= 0)  
        return(0);       
      
    for (l = 0; l < nn; ++l) {

        j = ctx->GD_FPI[i][l];                       /* edge from i to j */
        if (j == i || j > nmax || nflag[j] == 0)
            continue;

        if (gdd_adj(ctx, i,j,ctx->PMGN) >= 0.0) {

            cnt = 0;
            for (ll = 0; ll < nn; ++ll) {

                k = ctx->GD_FPI[i][ll];              /* edge from i to k */
                if (k == i || k == j || k > nmax || nflag[k] == 0)
                    continue;

                if (gdd_adj(ctx, i,k,ctx->PMGN) >= 0.0) {
                    if (gdd_adj(ctx, k,j,ctx->PMGN) >= 0.0)
                        cnt++;
                }
            }
            if (cnt > 0) {
                if (d != NULL)   
                    d[ne] = j;

                if (nflag1 != NULL)  
                    nflag1[j] = 2;

                ne++;
                *rs += cnt;
            }
        }
    }
    return(ne);
}

/* ------------------------------------------------------------------------ */
/*  ggcliq      Generalized cliques.  Only with gdd option 1, undirected.   */
/*                                                                          */  
/*              ggcliq(                                                     */
/*                  opt=...,        1 = cliques                             */
/*                                  2 = clans                               */
/*                                  3 = clubs                               */
/*                  alg=...,        1 = agglomerative                       */
/*                                  2 = enumerative                         */
/*                  sc=...,         max length, def. 1                      */
/*                  min=...,        min number of elements, def. 3          */
/*                  gn = ...,       graph number,                           */
/*                  nfmt = ...,     integer format, def. 4                  */
/*                  max = ...,      memory for max cliques, def. 1000       */
/*              ) = fname;                                                  */
/*                                                                          */
/*              Return: 0 if OK, -1 if error.                               */

int ggcliq(TDAContext *ctx)
{
    register int i,j;
    int err,r;

    err = -1;
    if (check_cmd(ctx, 1))
        return(-1);

    printf1(ctx, "Generalized cliques. Current memory: %d bytes.\n",ctx->MemReq);
        
    if (parm(ctx, ctx->CmdBuf + 6,1,1))       /* get parameters */
        goto GGCLFin;

    if (gdd_check(ctx, ctx->PMGN,0))
        goto GGCLFin;
    if (gdd_tcheck(ctx, 1,1,0))          /* must be undirected, gdd option 1 */
        goto GGCLFin;

    if (ctx->PMOPT == 1)
        printf1(ctx, "Option 1: cliques (");
    else if (ctx->PMOPT == 2)
        printf1(ctx, "Option 2: clans (");
    else {          
        printf1(ctx, "Option 3: clubs (");
        ctx->PMOPT = 3;
    }
    if (ctx->PMSCFlg == 0 || ctx->PMSC < 0.0)
        ctx->PMSC = 1.0;

    if (ctx->PMMin < 3)
        ctx->PMMin = 3;

    printf1(ctx, "max path length: %g, ",ctx->PMSC);
    printf1(ctx, "min size: %d).\n",ctx->PMMin);


/**
    xxx(ctx);
exit(0);
**/





    if (ctx->PMALG != 2)
        ctx->PMALG = 1;
    printf1(ctx, "Algorithm: ");
    if (ctx->PMALG == 1)  
        printf1(ctx, "agglomerative.\n");
    else  
        printf1(ctx, "enumerative.\n");

    /* allocate memory for maximal max cliques */

    if (ctx->PMMax < 1)
        ctx->PMMax = 1000;

    if (g_gclptr_i(ctx, ctx->PMMax))
        goto GGCLFin;

    if (alloc_acn(ctx, ctx->GD_NP))   
        goto GGCLFin;
    if (alloc_acm(ctx, ctx->GD_NP))   
        goto GGCLFin;
    if (alloc_acx(ctx, ctx->GD_NP))   
        goto GGCLFin;
    if (alloc_ack(ctx, ctx->GD_NP + 1))   
        goto GGCLFin;
    if (alloc_acs(ctx, ctx->GD_NP))   
        goto GGCLFin;

    if (ctx->PMALG == 2) {
        if (alloc_ace(ctx, ctx->GD_NP))   
            goto GGCLFin;
    }
    ctx->GCL_CNT = 0;            /* count number of cliques */
    ctx->GG_REC = 0;             /* number of records */

    for (i = 0; i < ctx->GD_NP; ++i) {

        if (ctx->PMALG == 2 && ctx->AcE[i] != 0)
            continue;

        if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 50) {
            printfe(ctx, "Node: %7d     %c",i + 1,CR);
            fflushe(ctx);
        }
        if (ctx->PMALG == 1) {
            ctx->AcN[0] = i;
            ctx->AcM[i] = 1;
            r = g_gclfind(ctx, i,1,ctx->AcN,ctx->PMGN,ctx->PMSC,ctx->AcM,1,ctx->PMMin,ctx->AcS);
        }   
        else
            r = g_gclfind1(ctx, i,ctx->PMGN,ctx->PMSC,ctx->PMMin,ctx->AcM,ctx->AcE,ctx->AcS);

        if (r < 0)
            goto GGCLFin;        
     
    }
    if (ctx->SILENTFlg < 2 && ctx->GD_NP >= 50) {
        printfe(ctx, "\n");
        fflushe(ctx);
    }
    
   
    printf1(ctx, "END pn=%d\n",ctx->GCL_PN);
    for (i = 0; i < ctx->GCL_PN; ++i) {
        printf1(ctx, "nptr=%2d : ",ctx->GCL_NPTR[i]);
        for (j = 0; j < ctx->GCL_NPTR[i]; ++j) 
            printf1(ctx, "%d ",ctx->GCL_PTR[i][j]);
        newline(ctx);
    }
      
    /* print, depending on option */

    for (i = 0; i < ctx->GCL_PN; ++i) {

        r = ctx->GCL_NPTR[i];

        if (ctx->PMOPT == 3) {       /* separate procedure for clubs */
                 
            if (g_gclclubs(ctx, r,ctx->GCL_PTR[i],ctx->PMGN,ctx->PMSC,ctx->AcN))
                goto GGCLFin;
        }
        else {
            if (ctx->PMOPT == 2 && g_gclcheck(ctx, r,ctx->GCL_PTR[i],ctx->PMGN,ctx->PMSC) == 0) 
                continue;

            g_gclprn(ctx, r,ctx->GCL_PTR[i]);
        }
    }
    printf1(ctx, "%d records written to: %s\n",ctx->GG_REC,ctx->PMFdName);
    err = 0;

GGCLFin:       
    g_gclptr_i(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_gclfind(in,n,nodes,gn,s,nflag,level,min,itmp)                         */
/*                                                                          */
/*  nodes[i] (i = 0,...,n-1) is a list of nodes. Find all subgraphs         */
/*  which contain this set of nodes by adding nodes. If at least one        */
/*  node can be added to nodes[], return 1, otherwise 0. If insufficient    */
/*  memory, return -1. in is the starting node from which the function      */
/*  is called for the first time. This implies that the function is         */
/*  already called with all i < in. min is minimal number of elements.      */
/*                                                                          */  

int g_gclfind(TDAContext *ctx, int in,int n,int *nodes,int gn,double s,int *nflag, int level,int min,int *itmp)
{
    register int i,j,l;
    int err,ip,m,r,fnd,add;

    err = -1;
/**
    printf1(ctx, "n=%d nodes: ",n);
    for (i = 0; i < n; ++i)
        printf1(ctx, "%d ",nodes[i]);
    newline(ctx);
**/

    for (i = 0; i < n; ++i) {
        if (nodes[i] < in) {
/*          printf1("RETURN 111111111111\n");
*/
            return(1);
        }
    }
    fnd = 0;
    for (ip = 0; ip < n; ++ip) {

        i = nodes[ip];
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];         /* edge from i to j */

            if (j == i || nflag[j])
                continue;

            /* check whether j can be added to nodes[] */

            if (gdd_adj(ctx, i,j,ctx->PMGN) >= 0.0) {    

                nflag[j] = level + 1;

                add = g_gcladd(ctx, j,n,nodes,gn,s);
                if (add < 0)  
                    goto GGCLFFin;

                else if (add > 0) {     /* j can be added */

                    if (j < in) {
/*                      printf1("RETURN 999999999999\n");
*/
                        err = 1;
                        goto GGCLFFin;
                    }
                    fnd++;
                    nodes[n++] = j;

                    r = g_gclfind(ctx, in,n,nodes,gn,s,nflag,level + 1,min,itmp);

                    if (r < 0)   
                        goto GGCLFFin;

                    else if (r == 0 && n >= min) {
                            
/**
                        printf1(ctx, "NEW group n=%d nodes: ",n);
                        for (k = 0; k < n; ++k) {
                            printf1(ctx, "%d ",nodes[k]);
                        }
                        newline(ctx);
**/

                        /* if not already found add to pointer list */

                        if (g_gclptr_c(ctx, n,nodes) == 0) {
                            if (g_gclptr_aa(ctx, n,nodes,itmp))
                                goto GGCLFFin;
                        }
                    }
                    n--;              
                }
            }
        }

    }
    err = fnd;

GGCLFFin:
    for (i = 0; i < ctx->GD_NP; ++i) {
        if (nflag[i] == level + 1)
            nflag[i] = 0;
    }
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  g_gcladd(i,n,nodes,gn,s)                                                */
/*                                                                          */
/*  Check whether node i can be added to nodes[j] (j = 0,...,n-1).          */
/*  Return 1 if i can be added, 0 if i cannot be added, -1 if error.        */
/*                                                                          */
/*  Note: this function requires that AcX and AcK are allocated for         */
/*  GD_NP nodes.                                                            */

int g_gcladd(TDAContext *ctx, int i,int n,int *nodes,int gn,double s)
{
    register int j;

    g_sp1(ctx, i,ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s);      

    for (j = 0; j < n; ++j) {
        if (ctx->AcX[nodes[j]] >= ctx->DBLMAX)  
            return(0);
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  g_gclptr_i(max)     create pointer for max subgraphs. If max = 0        */
/*                      free previously allocated memory.                   */
/*                      Return 0 if OK, -1 if error.                        */

int g_gclptr_i(TDAContext *ctx, int max)
{
    register int i;

    if (max > 0) {
        if (!(ctx->GCL_NPTR = (int *)calloc((size_t)(max),sizeof(int)))) { 
            p_err(ctx, -2,1);
            return(-1);          
        }
        if (!(ctx->GCL_PTR = (int **)calloc((size_t)(max),sizeof(int *)))) { 
            p_err(ctx, -2,1);
            free((char *)ctx->GCL_NPTR);
            return(-1);
        }
        ctx->GCL_PTRA = max;        
        memrq(ctx, max,sizeof(int) + sizeof(int *));
        return(0);
    }
    if (ctx->GCL_PTRA == 0)
        return(0);

    for (i = 0; i < ctx->GCL_PN; ++i) {
        if (ctx->GCL_NPTR[i] > 0) {
            free((char *)ctx->GCL_PTR[i]);
            memrq(ctx, -ctx->GCL_NPTR[i],sizeof(int));
        }
    }
    free((char *)ctx->GCL_PTR);
    memrq(ctx, -ctx->GCL_PTRA,sizeof(int *));
    free((char *)ctx->GCL_NPTR);
    memrq(ctx, -ctx->GCL_PTRA,sizeof(int));
    ctx->GCL_PN = ctx->GCL_PTRA = 0;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclptr_a(n,nodes)     Add node list to clique pointers.               */
/*                          Return 0 if OK, -1 if error.                    */

int g_gclptr_a(TDAContext *ctx, int n,int *nodes)
{
    register int i;

    if (ctx->GCL_PN >= ctx->GCL_PTRA) {
        printf1(ctx, "Error: insufficient storage with max = %d.\n",ctx->GCL_PTRA);
        return(-1);
    }
    if (!(ctx->GCL_PTR[ctx->GCL_PN] = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);          
    }
    memrq(ctx, n,sizeof(int));
    ctx->GCL_NPTR[ctx->GCL_PN] = n;

    for (i = 0; i < n; ++i)  
        ctx->GCL_PTR[ctx->GCL_PN][i] = nodes[i];
    ctx->GCL_PN++;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclptr_aa(n,nodes)    Add node list to clique pointers.               */
/*                          Sort before adding.                             */
/*                          Return 0 if OK, -1 if error.                    */

int g_gclptr_aa(TDAContext *ctx, int n,int *nodes,int *itmp)
{
    register int i;

    if (ctx->GCL_PN >= ctx->GCL_PTRA) {
        printf1(ctx, "Error: insufficient storage with max = %d.\n",ctx->GCL_PTRA);
        return(-1);
    }
    if (sortdpi(ctx, n,nodes,itmp))
        return(-1);

    if (!(ctx->GCL_PTR[ctx->GCL_PN] = (int *)calloc((size_t)(n),sizeof(int)))) { 
        p_err(ctx, -2,1);
        return(-1);          
    }
    memrq(ctx, n,sizeof(int));
    ctx->GCL_NPTR[ctx->GCL_PN] = n;

    for (i = 0; i < n; ++i)  
        ctx->GCL_PTR[ctx->GCL_PN][i] = nodes[itmp[i]];
    ctx->GCL_PN++;
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclptr_c(n,nodes)     Return 1 if nodes[i] (i = 0,...,n-1) is         */
/*                          contained in one already stored subgraph.       */
/*                          Otherwise return 0.                             */

int g_gclptr_c(TDAContext *ctx, int n,int *nodes)
{
    register int i,j,k;
    int ni,m,fnd,fnd1;

    for (i = 0; i < ctx->GCL_PN; ++i) {
        m = ctx->GCL_NPTR[i];
        if (n > m)
            continue;

        fnd = 1;
        for (j = 0; j < n; ++j) {
            ni = nodes[j];
            if (ni < ctx->GCL_PTR[i][0] || ni > ctx->GCL_PTR[i][m - 1]) {
                fnd = 0;
                break;
            }
            fnd1 = 0;
            for (k = 0; k < m; ++k) {
                if (ni == ctx->GCL_PTR[i][k]) {
                    fnd1 = 1;
                    break;
                }
            }
            if (fnd1 == 0) {
                fnd = 0;
                break;
            }
        }
        if (fnd)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclptr_cc(n,nodes,is)   Return 1 if nodes[i] (i = 0,...,n-1) is       */
/*                            contained in one already stored subgraph.     */
/*                            Otherwise return 0. Skip entry is.            */

int g_gclptr_cc(TDAContext *ctx, int n,int *nodes,int is)
{
    register int i,j,k;
    int ni,m,fnd,fnd1;

    for (i = 0; i < ctx->GCL_PN; ++i) {
        if (i == is)
            continue;

        m = ctx->GCL_NPTR[i];
        if (n > m)
            continue;

        fnd = 1;
        for (j = 0; j < n; ++j) {
            ni = nodes[j];
            if (ni < ctx->GCL_PTR[i][0] || ni > ctx->GCL_PTR[i][m - 1]) {
                fnd = 0;
                break;
            }
            fnd1 = 0;
            for (k = 0; k < m; ++k) {
                if (ni == ctx->GCL_PTR[i][k]) {
                    fnd1 = 1;
                    break;
                }
            }
            if (fnd1 == 0) {
                fnd = 0;
                break;
            }
        }
        if (fnd)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_sp1(in,n,gn,d,w,nflg,s)                                               */
/*                                                                          */
/*  Find all nodes that are reachable from node in with a path that has     */
/*  maximal length s. (Modified version of CALGO 562).                      */
/*                                                                          */
/*  in = node number                                                        */
/*  n  = number of nodes                                                    */
/*  gn = graph number                                                       */
/*  d[i] (i = 0,...,GD_NP-1)    list of shortest path values                */
/*  ndec[i] array used for double-sided deque.                              */
/*                                                                          */

void g_sp1(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec,double s)
{
    (void)gn;        /* unused: the signature is shared */
    register int i,j,l,np;   
    int m,max;
    double dj,tmp;   

    for (i = 0; i < n; ++i) {
        d[i] = ctx->DBLMAX;
        ndec[i] = 0;
    }
    max = n + 2;
    d[in] = 0.0;
    ndec[in] = max;  
    np = i = in;

    while (1) {
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];                   
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0) {
                dj = d[i] + tmp;
                if (dj < d[j] && dj <= s) {
                    d[j] = dj;
                    /** w[j] = i; **/
                    if (ndec[j] == 0) {
                        ndec[np] = j + 1;
                        np = j;
                        ndec[j] = max;
                    }
                    else if (ndec[j] < 0) {
                        ndec[j] = ndec[i];
                        ndec[i] = j + 1;
                        if (np == i)
                            np = j;
                    }
                }        
            }
        }
        j = ndec[i];
        if (j >= max)
            break;
        ndec[i] = -j;        
        i = j - 1;
        if (i < 0)
            gerr_exit(ctx, 104);
    }
}

/* ------------------------------------------------------------------------ */
/*  g_sp2(in,n,gn,d,w,nflg,s,nn,nodes)                                      */
/*                                                                          */
/*  Find all nodes that are reachable from node in with a path that has     */
/*  maximal length s. Use only path in nodes[j] (j = 0,...,nn-1).           */
/*                                                                          */
/*  in = node number                                                        */
/*  n  = number of nodes                                                    */
/*  gn = graph number                                                       */
/*  d[i] (i = 0,...,GD_NP-1)    list of shortest path values                */
/*  ndec[i] array used for double-sided deque.                              */
/*                                                                          */

void g_sp2(TDAContext *ctx, int in,int n,int gn,double *d,int *ndec,double s,int nn,int *nodes)
{
    (void)gn;        /* unused: the signature is shared */
    register int i,j,l,np;   
    int m,max;
    double dj,tmp;   

    for (i = 0; i < n; ++i) {
        d[i] = ctx->DBLMAX;
        ndec[i] = 0;
    }
    max = n + 2;
    d[in] = 0.0;
    ndec[in] = max;  
    np = i = in;

    while (1) {
        m = ctx->GD_FPN[i];
        for (l = 0; l < m; ++l) {
            j = ctx->GD_FPI[i][l];                   
            tmp = gdd_adj(ctx, i,j,ctx->PMGN);
            if (tmp >= 0.0) {
                dj = d[i] + tmp;
                if (dj < d[j] && dj <= s) {
                    if (g_nfind(ctx, j,nn,nodes)) {
                        d[j] = dj;
                        if (ndec[j] == 0) {
                            ndec[np] = j + 1;
                            np = j;
                            ndec[j] = max;
                        }
                        else if (ndec[j] < 0) {
                            ndec[j] = ndec[i];
                            ndec[i] = j + 1;
                            if (np == i)
                                np = j;
                        }
                    }
                }        
            }
        }
        j = ndec[i];
        if (j >= max)
            break;
        ndec[i] = -j;        
        i = j - 1;
        if (i < 0)
            gerr_exit(ctx, 105);
    }
}

/* ------------------------------------------------------------------------ */
/*  g_nfind(i,n,nodes)  Return 1 if i is contained in nodes[j] (j=0,n-1)    */
/*                      Otherwise return 0.                                 */

int g_nfind(TDAContext *ctx, int i,int n,int *nodes)
{
    (void)ctx;        /* unused: the signature is shared */
    register int j;

    for (j = 0; j < n; ++j) {
        if (nodes[j] == i)
            return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclcheck(n,nodes,gn,s)                                                */
/*                                                                          */
/*  Return 1 if path length between all pair of nodes in nodes[] is <= s,   */
/*  otherwise return 0.                                                     */

int g_gclcheck(TDAContext *ctx, int n,int *nodes,int gn,double s)
{
    register int i,j;
     
    for (i = 1; i < n; ++i) {
        g_sp2(ctx, nodes[i],ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s,n,nodes);
        for (j = 0; j < i; ++j) {
            if (ctx->AcX[nodes[j]] >= ctx->DBLMAX)  
                return(0);
        }
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  g_gclcheck1(n,nodes,gn,s)                                               */
/*                                                                          */
/*  Return 1 if path length between all pair of nodes in nodes[] is <= s,   */
/*  otherwise return 0. Use full graph.                                     */

int g_gclcheck1(TDAContext *ctx, int n,int *nodes,int gn,double s)
{
    register int i,j;
     
    for (i = 1; i < n; ++i) {
        g_sp1(ctx, nodes[i],ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s);
        for (j = 0; j < i; ++j) {
            if (ctx->AcX[nodes[j]] >= ctx->DBLMAX)  
                return(0);
        }
    }
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  g_gclclubs. Find clubs in generalized cliques.                          */
/*  Return 0 if OK, -1 if insufficient memory.                              */

int g_gclclubs(TDAContext *ctx, int n,int *nodes,int gn,double s,int *idx)            
{
    register int i,j;
    int m,ne,nn,f,fnd;

    if (n <= 4) {
        g_gclprn(ctx, n,nodes);
        return(0);
    }

    if (alloc_acc(ctx, n))   
        return(-1);   
    if (alloc_aci(ctx, n))   
        return(-1);   
    if (alloc_acj(ctx, n))   
        return(-1);   

    ne = 0;
    for (i = 1; i < n; ++i) {
        g_sp2(ctx, nodes[i],ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s,n,nodes);
        for (j = 0; j < i; ++j) {
            if (ctx->AcX[nodes[j]] >= ctx->DBLMAX) {
                if (ctx->AcC[i] == 0) {
                    ctx->AcC[i] = 1;
                    ctx->AcI[ne++] = i;
                }
                if (ctx->AcC[j] == 0) {
                    ctx->AcC[j] = 1;
                    ctx->AcI[ne++] = j;
                }
            }
        }
    }
    if (ne == 0) {
        g_gclprn(ctx, n,nodes);
        return(0);
    }

/***
printf1(ctx, "acc ne=%d : ",ne);
for (i = 0; i < n; ++i)
    printf1(ctx, "%d ",AcC[i]);
newline(ctx);

printf1(ctx, "aci: ");
for (i = 0; i < ne; ++i)
    printf1(ctx, "%d ",AcI[i]);
newline(ctx);

**/



    m = 1;
    idx[0] = 0;

    while (idx[0] < ne) {
        fnd = 0;

        printf1(ctx, "m=%d : ",m);
        for (i = 0; i < m; ++i)
            printf1(ctx, "%d ",idx[i]);
        newline(ctx);

        nn = 0;
        for (i = 0; i < n; ++i) {
            f = 1;
            for (j = 0; j < m; ++j) {
                if (ctx->AcI[idx[j]] == nodes[i]) {
                    f = 0;
                    break;
                }
            }
            if (f)
                ctx->AcJ[nn++] = nodes[i];
        }
/**
printf1(ctx, "check nodes : ");
for (i = 0; i < nn; ++i)
    printf1(ctx, "%d ",AcJ[i]);
newline(ctx);
**/
        if (g_gclcheck(ctx, nn,ctx->AcJ,gn,s)) {
            g_gclprn(ctx, nn,ctx->AcJ);
            fnd = 1;
        }
/**
tda_out("fnd=%d\n",fnd);
**/
        if (m == 1 && idx[0] == 2)
            fnd = 1;

        if (m < ne && idx[m - 1] + 1 < ne) {
            if (fnd)  
                idx[m - 1] += 1;
            else {
                idx[m] = idx[m - 1] + 1;
                m++;
            }
        }
        else if (--m > 0 && idx[m - 1] < ne) {
            idx[m - 1] += 1;
        }
        else
            break;
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  g_gclprn(n,nodes)       print to output file.                           */

void g_gclprn(TDAContext *ctx, int n,int *nodes)
{
    register int i;

    ctx->GG_REC++;
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,ctx->GG_REC);
    rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,n);
    for (i = 0; i < n; ++i)
        rt_fprintf_i(ctx, ctx->PMFd,ctx->PMNFmtS,gdd_node(ctx, nodes[i]));
    fprintf(ctx->PMFd,"\n");
}

/* ------------------------------------------------------------------------ */
/*  g_gclfind1(in,gn,s,min,idx,iflag)                                       */
/*                                                                          */
/*  Assume that AcX and AcK already allocated.                              */

int g_gclfind1(TDAContext *ctx, int in,int gn,double s,int min,int *idx,char *iflag,int *itmp)
{
    register int i,j;
    int err,n,ne,m,f,nn,first,r,mm;
    double tmp;
    int fnd[100];

    for (i = 0; i < 100; ++i)
        fnd[i] = 0;
    
    /*******
    printf1(ctx, "in=%d\n",in);
    ****/
     
    err = -1;
    if (alloc_acn(ctx, ctx->GD_NP))   
        goto GCLF1Fin;

    g_sp1(ctx, in,ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s);      
/*  
    printf1(ctx, "in=%d AcX: ",in);
    for (j = 0; j < GD_NP; ++j)
        printf1(ctx, "%g ",AcX[j]);
    newline(ctx);
*/     
    tmp = s / 2.0;
    n = 0;
    for (j = 0; j < ctx->GD_NP; ++j) {
        if (ctx->AcX[j] < ctx->DBLMAX) {
            ctx->AcN[n++] = j;
            if (ctx->AcX[j] <= tmp)
                iflag[j] = 1;
        }   
    }
    if (n < min)
        return(0);



/*
    printf1(ctx, "iflag: ");
    for (i = 0; i < GD_NP; ++i)
        printf1(ctx, "%d ",iflag[i]);
    newline(ctx);
*/ 
    printf1(ctx, "n=%d AcN: ",n);
    for (i = 0; i < n; ++i)
        printf1(ctx, "%d ",ctx->AcN[i]);
    newline(ctx);
    
    if (alloc_acc(ctx, n))   
        goto GCLF1Fin;
    if (alloc_aci(ctx, n))   
        goto GCLF1Fin;

    ne = 0;
    for (i = 1; i < n; ++i) {
        g_sp1(ctx, ctx->AcN[i],ctx->GD_NP,gn,ctx->AcX,ctx->AcK,s);
        for (j = 0; j < i; ++j) {
            if (ctx->AcX[ctx->AcN[j]] >= ctx->DBLMAX) {
                if (ctx->AcC[i] == 0) {
                    ctx->AcC[i] = 1;
                    ctx->AcI[ne++] = i;
                }
                if (ctx->AcC[j] == 0) {
                    ctx->AcC[j] = 1;
                    ctx->AcI[ne++] = j;
                }
            }
        }
    }
  
    printf1(ctx, "ne=%d AcI: ",ne);
    for (i = 0; i < ne; ++i)
        printf1(ctx, "%d ",ctx->AcI[i]);
    newline(ctx);
   
    if (ne == 0) {  /* if not already found add to pointer list */

        if (g_gclptr_c(ctx, n,ctx->AcN) == 0) {
            if (g_gclptr_aa(ctx, n,ctx->AcN,itmp))
                goto GCLF1Fin;
        }
    }
    if (alloc_acj(ctx, ctx->GD_NP))   
        goto GCLF1Fin;

tda_out("n=%d ne=%d\n",n,ne);

    mm = imin(ctx, ne,n - min);
tda_out("mm=%d\n",mm);



    for (m = 1; m <= mm; ++m) {

        first = 1;
        while ((r = comb_xx(ctx, ne,m,idx,first,fnd))) {
            first = 0;


            printf1(ctx, "fnd : ");
            for (i = 0; i < ne; ++i)
                printf1(ctx, "%d ",fnd[i]);
            newline(ctx);

            printf1(ctx, "m=%d r=%2d : ",m,r);
            for (i = 0; i < m; ++i)
                printf1(ctx, "%d ",idx[i]);
            newline(ctx);


            nn = 0;
            for (i = 0; i < n; ++i) {
                f = 1;
                for (j = 0; j < m; ++j) {
                    if (ctx->AcI[idx[j]] == i) {
                        f = 0;
                        break;
                    }
                }
                if (f)
                    ctx->AcJ[nn++] = ctx->AcN[i];
            }
     
printf1(ctx, "check nodes : ");
for (i = 0; i < nn; ++i)
    printf1(ctx, "%d ",ctx->AcJ[i]);
newline(ctx);
     
            if (g_gclcheck1(ctx, nn,ctx->AcJ,gn,s)) {
   printf1(ctx, "FOUND: ");
for (i = 0; i < nn; ++i)
    printf1(ctx, "%d ",ctx->AcJ[i]);
newline(ctx);

                for (j = 0; j < m; ++j) {
                    fnd[idx[j]] = 1;
                }
  
  
                if (g_gclptr_c(ctx, nn,ctx->AcJ) == 0) {
        printf1(ctx, "added\n");
  

                    if (g_gclptr_aa(ctx, nn,ctx->AcJ,itmp))
                        goto GCLF1Fin;
                }
            }
        }
    }
    err = 0;

GCLF1Fin:
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  g_gclprn(n,nodes)       print to output file.                           */

void xxx(TDAContext *ctx)
{
    register int i;
    int m,first,r;
    int idx[100],fnd[100];

    for (i = 0; i < 20; ++i)
        fnd[i] = 0;


    for (m = 1; m <= 3; ++m) {

        first = 1;
        while ((r = comb_xx(ctx, 4,m,idx,first,fnd))) {
            first = 0;


            printf1(ctx, "m=%d r=%2d : ",m,r);
            for (i = 0; i < m; ++i)
                printf1(ctx, "%d ",idx[i]);
            newline(ctx);
/*  
            if (m == 1 && idx[0] == 2)
                fnd[2] = 1;
*/
            if (m == 2 && idx[0] == 0 && idx[1] == 3) {
                fnd[0] = 1;
                fnd[3] = 1;
            }

        }
    }
}

/* ------------------------------------------------------------------------ */
/*  comb_xx(n,m,com,first)                                                  */
/*                                                                          */  
/*  Generate all subsets of {0,...,n-1} with m elements.                    */
/*  Return 1 if new subset found, otherwise return 0.                       */

int comb_xx(TDAContext *ctx, int n,int m,int *com,int first,int *idx)
{
    (void)ctx;        /* unused: the signature is shared */
    register int k,j,l;
    int r;

    if (first) {
        for (k = 0; k < m; ++k)
            com[k] = n;
    }
CONT:
    r = 1;
    for (k = 1; k <= m; ++k) {
        j = m - k;
        if (com[j] < n - k) {
            l = com[j];
            for (k = j; k < m; ++k) {
                com[k] = ++l;
                if (idx[l])
                    r = -1;        
            }
            if (r < 0)
                goto CONT;
            return(r);
        }
    }
    for (k = 0; k < m; ++k) {
        com[k] = k;
        if (idx[k])
            r = -1;         
    }
    if (first) {
        first = 0;
        if (r < 0)
            goto CONT;
        return(r);
    }
    else
        return(0);
}

