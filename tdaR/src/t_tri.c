#include "tda_rhooks.h"
/****************************************************************************/
/*  t_tri                                                                   */
/*                                                                          */
/*  TDA. Program for Transition Data Analysis, written by Goetz Rohwer.     */
/*  Copyright (C) 1989,1991-2002 Goetz Rohwer. All rights reserved.         */
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
/*  NOTE: this file contains source code from a C program developed by      */
/*  Steve J. Fortune: an algorithm to calculate Voronoi diagrams and        */
/*  Delaunay triangulations. The copyright is, of course, with Fortune.     */

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_gm.h"
#include "t_sd.h"
#include "tda_context.h"

/* ------------------------------------------------------------------------ */
/* Global variables and definitions                                         */

#define F_DELETED -2
#define F_le 0
#define F_re 1
           
/* Free-list nodes are recycled Site/Edge/Halfedge objects.  The link is
   stored in the first sizeof(void *) bytes of the node with memcpy, so the
   node is never accessed through a type other than its own. */
struct  Freelist    {
    void        *head;
    int         nodesize;
};
struct Point    {
    double x,y;
};
/* structure used both for sites and for nodes */
struct Site {
    struct  Point   coord;
    int     sitenbr;
    int     refcnt;
};
struct Edge {
    double      a,b,c;
    struct  Site    *ep[2];
    struct  Site    *reg[2];
    int     edgenbr;
};
struct Halfedge {
    struct Halfedge *ELleft, *ELright;
    struct Edge *ELedge;
    int     ELrefcnt;
    char        ELpm;
    struct  Site    *vertex;
    double      ystar;
    struct  Halfedge *PQnext;
};
struct Freelist sfl;
struct  Freelist efl;
struct Freelist hfl;
struct  Site    *sites;
struct  Site    *bottomsite;
struct  Halfedge *ELleftend, *ELrightend;
struct  Halfedge **ELhash;
struct  Halfedge *PQhash;


 

float F_xmin, F_xmax, F_ymin, F_ymax, F_deltax, F_deltay;
float pxmin, pxmax, pymin, pymax, cradius;

/*  functions in t_tri.c */

int sdvd(TDAContext *ctx);
int reg_memory(TDAContext *ctx, int opt);

/* the following functions are taken from Fortune's C program */

int voronoi(TDAContext *ctx, struct Site *(*nextsite)(TDAContext *));
struct Site *nextone(TDAContext *ctx);
int ELinitialize(TDAContext *ctx);
struct Halfedge *HEcreate(TDAContext *ctx, struct Edge *e,int pm,int *err);
void ELinsert(TDAContext *ctx, struct Halfedge *lb,struct Halfedge *new);
struct Halfedge *ELgethash(TDAContext *ctx, int b);
struct Halfedge *ELleftbnd(TDAContext *ctx, struct Point *p);
void ELdelete(TDAContext *ctx, struct Halfedge *he);
struct Halfedge *ELright(TDAContext *ctx, struct Halfedge *he);
struct Halfedge *ELleft(TDAContext *ctx, struct Halfedge *he);
struct Site *leftreg(TDAContext *ctx, struct Halfedge *he);
struct Site *rightreg(TDAContext *ctx, struct Halfedge *he);
void geominit(TDAContext *ctx);
struct Edge *bisect(TDAContext *ctx, struct Site *s1,struct Site *s2,int *err);
struct Site *intersect(TDAContext *ctx, struct Halfedge *el1,struct Halfedge *el2,int *err);
int right_of(TDAContext *ctx, struct Halfedge *el,struct Point *p);
void endpoint(TDAContext *ctx, struct Edge *e,int lr,struct Site *s);
float dist(TDAContext *ctx, struct Site *s,struct Site *t);
void makevertex(TDAContext *ctx, struct Site *v);
void deref(TDAContext *ctx, struct   Site *v);
void ref(TDAContext *ctx, struct Site *v);
void PQinsert(TDAContext *ctx, struct Halfedge *he,struct Site *v,float  offset);
void PQdelete(TDAContext *ctx, struct Halfedge *he);
int PQbucket(TDAContext *ctx, struct Halfedge *he);
int PQempty(TDAContext *ctx);
struct Point PQ_min(TDAContext *ctx);
struct Halfedge *PQextractmin(TDAContext *ctx);
int PQinitialize(TDAContext *ctx);
void freeinit(TDAContext *ctx, struct    Freelist *fl,int    size);
char *getfree(TDAContext *ctx, struct    Freelist *fl,int *err);
void makefree(TDAContext *ctx, void *curr,struct Freelist *fl);
void out_bisector(TDAContext *ctx, struct Edge *e);
void out_ep(TDAContext *ctx, struct Edge *e);
void out_vertex(TDAContext *ctx, struct Site *v);
void out_triple(TDAContext *ctx, struct Site *s1,struct Site *s2,struct Site *s3);
char *myalloc(TDAContext *ctx, unsigned n,int *err);
void myalloc_free(TDAContext *ctx);
/*   char *malloc();    */    

/* ------------------------------------------------------------------------ */
/*  sdvd        Create and optionally polot Voronoi diagrams and            */
/*              Delaunay triangulations.                                    */
/*                                                                          */  
/*              sdvd(                                                       */
/*                  xyv=...,  two variables with coordinates                */
/*                  opt=...,  option, def. 1                                */
/*                            1 = Voronoi diagram (nodes and edges)         */
/*                            2 = Voronoi diagram (edge list)               */
/*                            3 = Delaunay triangulation (triangles)        */
/*                            4 = Delaunay triangulation (edge list)        */
/*                  rec=...,  bounding box for opt 1 (xmin,ymin,xmax,ymax)  */
/*                  m=...,    maximal number of memory requests, def. 10000 */
/*                  fmt=...,  print format for coordinates, def. 10.4       */
/*              ) output_file;                                              */
/*                                                                          */
/*  The output file is required. If the xyv parameter specifies two         */
/*  variables, these are taken for the x and y coordinates, respectively.   */
/*  By default, the command uses all type 1 objects (points) from the       */
/*  current spatial data structure.                                         */
/*                                                                          */
/*  The command uses an algorithm developed by S. J. Fortune.               */
/*                                                                          */
/*  Return: 0 if OK, -1 if error.                                           */

int sdvd(TDAContext *ctx)
{
    register int i,k;
    int err,r,nn,dflag,ix,iy;
    double dx,dy;

    err = -1;
    printf1(ctx, "Voronoi diagram and Delaunay triangulation. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 4,1,1))       /* get parameters */
        goto SDVDFin;

    if (ctx->PMFmtF == 0)
        pmfmt(ctx, 10,4);

    dflag = 0;
    if (ctx->PMNV == 2) {
        dflag = 1;
    }
    else if (check_sd(ctx, 0) == 0)    
        dflag = 2;

    if (dflag == 0) {
        printf1(ctx, "Error: need a spatial data matrix or two variables on left-hand side.\n");
        goto SDVDFin;
    }
    if (ctx->PMOPT > 4)
        ctx->PMOPT = 1;
    printf1(ctx, "Option %d. ",ctx->PMOPT);
    if (ctx->PMOPT == 3 || ctx->PMOPT == 4) {
        printf1(ctx, "Delaunay triangulation.\n");
        ctx->F_triangulate = 1;
    }
    else {
        printf1(ctx, "Voronoi diagram.\n");
        ctx->F_triangulate = 0;
    }
    if (ctx->PMMFlg == 0)
        ctx->F_MaxAlloc = 10000;
    else
        ctx->F_MaxAlloc = ctx->PMM;

    if (ctx->F_MaxAlloc < 10)
        ctx->F_MaxAlloc = 10;

    if (alloc_acxf(ctx, ctx->NOC + 1))
        goto SDVDFin; 
    if (alloc_acyf(ctx, ctx->NOC + 1))
        goto SDVDFin; 
    if (alloc_ack(ctx, ctx->NOC + 1))
        goto SDVDFin; 

    if (dflag == 1) {                   /* variables from data matrix */
        ix = ctx->PMVIdx[0];
        iy = ctx->PMVIdx[1];

        for (i = 0; i < ctx->NOC; ++i) {
            ctx->AcXF[i] = (float)get_data(ctx, ix,i);
            ctx->AcYF[i] = (float)get_data(ctx, iy,i);
        }
        nn = ctx->NOC;
    }
    else {                              /* spatial data */
        nn = 0;
        for (i = 0; i < ctx->NOC; ++i) {
            if ((int)get_data(ctx, ctx->SDVarSDTyp,i) != 1)
                continue;

            if (sd_getdata(ctx, i,0,0,1) < 1)
                goto SDVDFin;
             
            ctx->AcXF[nn] = (float)ctx->SDVarX[0];
            ctx->AcYF[nn] = (float)ctx->SDVarY[0];
            nn++;
        }
    }
    if (nn < 3) {
        printf1(ctx, "Error: need at least three points.\n");
        goto SDVDFin;
    }

    /* calculate Voronoi diagram or Delaunay triangulation */
    /* first get sorted input data into sites */

    ctx->F_NMReq = 0;
    if (reg_memory(ctx, ctx->F_MaxAlloc))
        goto SDVDFin;

    freeinit(ctx, &sfl, sizeof *sites);

    sites = (struct Site *)(void *)myalloc(ctx,(unsigned int)((size_t)nn * sizeof *sites),&r);
    if (r)
        goto SDVDFin;
    ctx->F_nsites = nn;
        
    if (sortdp2f(ctx, nn,ctx->AcYF,ctx->AcXF,ctx->AcK))
        goto SDVDFin;

    for (i = 0; i < nn; ++i) {
        k = ctx->AcK[i];
        sites[i].sitenbr = k;
        sites[i].refcnt = 0;
        sites[i].coord.x = (double)(ctx->AcXF[k]);
        sites[i].coord.y = (double)(ctx->AcYF[k]);
    }
    F_xmin = (float)(sites[0].coord.x); 
    F_xmax = (float)(sites[0].coord.x);
    for(i = 1; i < ctx->F_nsites; ++i) {
        if((double)((sites[i].coord.x)) < (double)(F_xmin))
            F_xmin = (float)(sites[i].coord.x);
        if((double)((sites[i].coord.x)) > (double)(F_xmax))
            F_xmax = (float)(sites[i].coord.x);
    } 
    F_ymin = (float)(sites[0].coord.y);
    F_ymax = (float)(sites[ctx->F_nsites - 1].coord.y);

    printf1(ctx, "\nBounding box of input data: %lg, %lg, %lg, %lg\n",
        (double)(F_xmin),(double)(F_xmax),(double)(F_ymin),(double)(F_ymax));

    dx = (double)((F_xmax - F_xmin)) / 10.0;
    dy = (double)((F_ymax - F_ymin)) / 10.0;
        
    if (dx <= ctx->EPSI1 || dy <= ctx->EPSI1) {
        printf1(ctx, "Error: all points are (almost) collinear.\n");
        goto SDVDFin;
    }
    if (ctx->PMRECFlg == 0) {
        ctx->PMRECXMin = (double)(F_xmin) - dx;
        ctx->PMRECYMin = (double)(F_ymin) - dy;
        ctx->PMRECXMax = (double)(F_xmax) + dx;
        ctx->PMRECYMax = (double)(F_ymax) + dy;
    }
    else {
        ctx->PMRECXMin = (double)(dmin(ctx, ctx->PMRECXMin,(double)F_xmin));
        ctx->PMRECYMin = (double)(dmin(ctx, ctx->PMRECYMin,(double)F_ymin));
        ctx->PMRECXMax = (double)(dmax(ctx, ctx->PMRECXMax,(double)F_xmax));
        ctx->PMRECYMax = (double)(dmax(ctx, ctx->PMRECYMax,(double)F_ymax));
    }
    if (ctx->F_triangulate == 0) {
        printf1(ctx, "Bounding box for Voronoi diagram: %lg, %lg, %lg, %lg\n",
                 ctx->PMRECXMin,ctx->PMRECYMin,ctx->PMRECXMax,ctx->PMRECYMax);
    } 
    ctx->F_NTriangles = ctx->F_NNodes = ctx->F_NNEdges = 0;
    ctx->F_sdid = ctx->F_nrec = ctx->F_siteidx = 0;

    geominit(ctx);
    if (voronoi(ctx, nextone))
        goto SDVDFin;

    printf1(ctx, "\nNumber of data points: %d\n",ctx->F_nsites);
    printf1(ctx, "Number of Voronoi nodes: %d\n",ctx->F_NNodes);
    printf1(ctx, "Number of Voronoi edges: %d\n",ctx->F_NNEdges);
    printf1(ctx, "Number of Delaunay triangles: %d\n",ctx->F_NTriangles);
    printf1(ctx, "%d records written to: %s\n",ctx->F_nrec,ctx->PMFdName);

    printf1(ctx, "\nNumber of memory requests: %d\n",ctx->F_NMReq);
    err = 0;

SDVDFin:
    myalloc_free(ctx);
    reg_memory(ctx, 0);
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  reg_memory(opt)                                                         */
/*                                                                          */  
/*  if opt = 1 allocate F_ASize and F_APtr for maximal F_MaxAlloc           */
/*  memory requests (by myalloc), otherwise free previously allocated       */
/*  memory. Return 0 if OK, -1 if error.                                    */

int reg_memory(TDAContext *ctx, int opt)
{
    if (ctx->F_ASizeA > 0) {
        free((char *)ctx->F_ASize);
        memrq(ctx, -ctx->F_ASizeA,sizeof(unsigned int));
        ctx->F_ASizeA = 0;
    }
    if (ctx->F_APtrA > 0) {
        free((char *)ctx->F_APtr);
        memrq(ctx, -ctx->F_APtrA,sizeof(char *));
        ctx->F_APtrA = 0;
    }
    if (opt == 0)
        return(0);

    if (!(ctx->F_ASize = (unsigned int *)calloc((size_t)(ctx->F_MaxAlloc),sizeof(unsigned int)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->F_ASizeA = ctx->F_MaxAlloc;
    memrq(ctx, ctx->F_ASizeA,sizeof(unsigned int));

    if (!(ctx->F_APtr = (char **)calloc((size_t)(ctx->F_MaxAlloc),sizeof(char *)))) {
        p_err(ctx, -2,1);
        return(-1);
    }
    ctx->F_APtrA = ctx->F_MaxAlloc;
    memrq(ctx, ctx->F_APtrA,sizeof(char *));
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  voronoi(next)                                                           */
/*                                                                          */
/*  If F_triangulate = 0, this function calculates a Voronoi diagram,       */
/*  otherwise a Delaunay triangulation.                                     */
/*                                                                          */
/*  This function and all function called by this function are taken from   */
/*  the source code written by S.J. Fortune.                                */
/*                                                                          */
/*  Return 0 if successfull, -1 if insufficient memory.                     */

int voronoi(TDAContext *ctx, struct Site *(*nextsite)(TDAContext *))
{
    struct Site *newsite, *bot, *top, *temp, *p;
    struct Site *v;
    struct Point newintstar = {0.0,0.0};
    int pm,err;
    struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
    struct Edge *e;

    if (PQinitialize(ctx))
        return(-1);

    bottomsite = (*nextsite)(ctx);
    /** out_site(bottomsite); **/
    if (ELinitialize(ctx))
        return(-1);

    newsite = (*nextsite)(ctx);
    while(1) {
 
        if (!PQempty(ctx))
            newintstar = PQ_min(ctx);

        if (newsite != (struct Site *)NULL 
               && (PQempty(ctx) 
                 || newsite -> coord.y < newintstar.y
                 || (newsite->coord.y == newintstar.y 
                 && newsite->coord.x < newintstar.x))) {
                                                     /* new site is smallest */
            /** out_site(newsite); **/
                lbnd = ELleftbnd(ctx, &(newsite->coord));
                rbnd = ELright(ctx, lbnd);
                bot = rightreg(ctx, lbnd);
                e = bisect(ctx, bot, newsite,&err);
                if (err)
                    return(-1);

                bisector = HEcreate(ctx, e,F_le,&err);
                if (err)
                    return(-1);

                ELinsert(ctx, lbnd, bisector);

                p = intersect(ctx, lbnd, bisector,&err);
                if (err)
                    return(-1);

                if (p != (struct Site *) NULL) {
                PQdelete(ctx, lbnd);
                        PQinsert(ctx, lbnd, p, dist(ctx, p,newsite));
            } 
                lbnd = bisector;
                bisector = HEcreate(ctx, e,F_re,&err);
                if (err)
                    return(-1);

                ELinsert(ctx, lbnd, bisector);
                p = intersect(ctx, bisector, rbnd,&err);
                if (err)
                    return(-1);

                if (p != (struct Site *) NULL) {
                PQinsert(ctx, bisector, p, dist(ctx, p,newsite)); 
                } 
                newsite = (*nextsite)(ctx);    
        }
        else if (!PQempty(ctx)) {
                                        /* intersection is smallest */
            lbnd = PQextractmin(ctx);
                llbnd = ELleft(ctx, lbnd);
                rbnd = ELright(ctx, lbnd);
                rrbnd = ELright(ctx, rbnd);
                bot = leftreg(ctx, lbnd);
                top = rightreg(ctx, rbnd);
                out_triple(ctx, bot, top, rightreg(ctx, lbnd));
                v = lbnd->vertex;
                makevertex(ctx, v);
                endpoint(ctx, lbnd->ELedge,lbnd->ELpm,v);
                endpoint(ctx, rbnd->ELedge,rbnd->ELpm,v);
                ELdelete(ctx, lbnd); 
                PQdelete(ctx, rbnd);
                ELdelete(ctx, rbnd); 
                pm = F_le;
                if (bot->coord.y > top->coord.y) {
                    temp = bot; bot = top; top = temp; pm = F_re;
            }
                e = bisect(ctx, bot,top,&err);
                if (err)
                    return(-1);

                bisector = HEcreate(ctx, e,pm,&err);
                if (err)
                    return(-1);

                ELinsert(ctx, llbnd, bisector);
                endpoint(ctx, e, F_re - pm, v);
                deref(ctx, v);
                p = intersect(ctx, llbnd, bisector,&err);
                if (err)
                    return(-1);

                if (p != (struct Site *) NULL) {
                PQdelete(ctx, llbnd);
                        PQinsert(ctx, llbnd, p, dist(ctx, p,bot));
            } 
                p = intersect(ctx, bisector, rrbnd,&err);
                if (err)
                    return(-1);

                if (p != (struct Site *) NULL) {
                PQinsert(ctx, bisector, p, dist(ctx, p,bot));
                } 
        }
        else
            break;
    } 
    for(lbnd=ELright(ctx, ELleftend); lbnd != ELrightend; lbnd=ELright(ctx, lbnd)) {
        e = lbnd -> ELedge;
            out_ep(ctx, e);
    } 
    return(0);
}

/* return a single in-storage site */

struct Site *nextone(TDAContext *ctx)
{
    for ( ; ctx->F_siteidx < ctx->F_nsites; ctx->F_siteidx += 1) {
        if (ctx->F_siteidx == 0 || sites[ctx->F_siteidx].coord.x != sites[ctx->F_siteidx-1].coord.x 
                               || sites[ctx->F_siteidx].coord.y != sites[ctx->F_siteidx-1].coord.y) {
            ctx->F_siteidx += 1;
                return (&sites[ctx->F_siteidx - 1]);
        } 
    } 
    return((struct Site *)NULL);
}


/* ======================================================================== */
/* code from edgelist.c                                                     */
/* ======================================================================== */

int ELinitialize(TDAContext *ctx)
{
    int i,err;

    freeinit(ctx, &hfl, sizeof **ELhash);
    ctx->ELhashsize = 2 * ctx->F_sqrt_nsites;
    ELhash = (struct Halfedge **)(void *)myalloc(ctx,(unsigned int)(sizeof *ELhash * (size_t)ctx->ELhashsize),&err);
    if (err)
        return(-1);

    for(i=0; i<ctx->ELhashsize; i +=1)
        ELhash[i] = (struct Halfedge *)NULL;
    ELleftend = HEcreate(ctx, (struct Edge *)NULL,0,&err);
    if (err)
        return(-1);

    ELrightend = HEcreate(ctx, (struct Edge *)NULL,0,&err);
    if (err)
        return(-1);

    ELleftend -> ELleft = (struct Halfedge *)NULL;
    ELleftend -> ELright = ELrightend;
    ELrightend -> ELleft = ELleftend;
    ELrightend -> ELright = (struct Halfedge *)NULL;
    ELhash[0] = ELleftend;
    ELhash[ctx->ELhashsize-1] = ELrightend;
    return(0);
}

struct Halfedge *HEcreate(TDAContext *ctx, struct Edge *e,int pm,int *err)
{
    struct Halfedge *answer;
    answer = (struct Halfedge *)(void *)getfree(ctx, &hfl,err);
    if (*err)
        return(answer);

    answer -> ELedge = e;
    answer -> ELpm = (char)(pm);
    answer -> PQnext = (struct Halfedge *) NULL;
    answer -> vertex = (struct Site *) NULL;
    answer -> ELrefcnt = 0;
    return(answer);
}

void ELinsert(TDAContext *ctx, struct Halfedge *lb,struct Halfedge *new)
{
    (void)ctx;        /* unused: the signature is shared */
    new -> ELleft = lb;
    new -> ELright = lb -> ELright;
    (lb -> ELright) -> ELleft = new;
    lb -> ELright = new;
}

/* Get entry from hash table, pruning any deleted nodes */

struct Halfedge *ELgethash(TDAContext *ctx, int b)
{
    struct Halfedge *he;

    if (b<0 || b>=ctx->ELhashsize)
        return((struct Halfedge *) NULL);
    he = ELhash[b]; 
    if (he == (struct Halfedge *) NULL || 
        he -> ELedge != (struct Edge *)F_DELETED )
        return (he);

    /* Hash table points to deleted half edge.  Patch as necessary. */

    ELhash[b] = (struct Halfedge *) NULL;
    if ((he -> ELrefcnt -= 1) == 0)                  
        makefree(ctx, he, &hfl);
    return ((struct Halfedge *) NULL);
}   

struct Halfedge *ELleftbnd(TDAContext *ctx, struct Point *p)
{
    int i, bucket;
    struct Halfedge *he;

    /* Use hash table to get close to desired halfedge */

    bucket = (int)((double)((p->x - (double)(F_xmin))) / (double)(F_deltax) * ctx->ELhashsize);
    if (bucket<0)
        bucket =0;
    if(bucket>=ctx->ELhashsize)
        bucket = ctx->ELhashsize - 1;
    he = ELgethash(ctx, bucket);
    if(he == (struct Halfedge *) NULL) {
        for(i=1; 1 ; i += 1) {
            if ((he=ELgethash(ctx, bucket-i)) != (struct Halfedge *) NULL)
                break;
            if ((he=ELgethash(ctx, bucket+i)) != (struct Halfedge *) NULL)
                break;
        }
        ctx->F_totalsearch += i;
    } 
    ctx->F_ntry += 1;

    /* Now search linear list of halfedges for the corect one */

    if (he==ELleftend  || (he != ELrightend && right_of(ctx, he,p))) {
        do {he = he -> ELright;} while (he!=ELrightend && right_of(ctx, he,p));
         he = he -> ELleft;
    }
    else 
        do {he = he -> ELleft;} while (he!=ELleftend && !right_of(ctx, he,p));

    /* Update hash table and reference counts */

    if(bucket > 0 && bucket <ctx->ELhashsize-1) {
        if(ELhash[bucket] != (struct Halfedge *) NULL) 
                    ELhash[bucket] -> ELrefcnt -= 1;
            ELhash[bucket] = he;
            ELhash[bucket] -> ELrefcnt += 1;
    } 
    return (he);
}
    
/* This delete routine can't reclaim node, since pointers from hash
   table may be present. */

void ELdelete(TDAContext *ctx, struct Halfedge *he)
{
    (void)ctx;        /* unused: the signature is shared */
    (he -> ELleft) -> ELright = he -> ELright;
    (he -> ELright) -> ELleft = he -> ELleft;
    he -> ELedge = (struct Edge *)F_DELETED;
}

struct Halfedge *ELright(TDAContext *ctx, struct Halfedge *he)
{
    (void)ctx;        /* unused: the signature is shared */
    return (he -> ELright);
}

struct Halfedge *ELleft(TDAContext *ctx, struct Halfedge *he)
{
    (void)ctx;        /* unused: the signature is shared */
    return (he -> ELleft);
}

struct Site *leftreg(TDAContext *ctx, struct Halfedge *he)
{
    (void)ctx;        /* unused: the signature is shared */
    if(he -> ELedge == (struct Edge *)NULL)
        return(bottomsite);
    return( he -> ELpm == F_le ? 
                he -> ELedge -> reg[F_le] : he -> ELedge -> reg[F_re]);
}

struct Site *rightreg(TDAContext *ctx, struct Halfedge *he)
{
    (void)ctx;        /* unused: the signature is shared */
    if(he -> ELedge == (struct Edge *)NULL)
        return(bottomsite);
    return( he -> ELpm == F_le ? 
                he -> ELedge -> reg[F_re] : he -> ELedge -> reg[F_le]);
}

/* ======================================================================== */
/* code from geometry.c                                                     */
/* ======================================================================== */

void geominit(TDAContext *ctx)
{
    struct Edge e;
    double sn;

    freeinit(ctx, &efl, sizeof e);
    ctx->F_nvertices = 0;
    ctx->F_nedges = 0;
    sn = (double)((float)(ctx->F_nsites + 4));
    ctx->F_sqrt_nsites = (int)sqrt(sn);
    F_deltay = F_ymax - F_ymin;
    F_deltax = F_xmax - F_xmin;
}

struct Edge *bisect(TDAContext *ctx, struct Site *s1,struct Site *s2,int *err)
{
    double dx,dy,adx,ady;
    struct Edge *newedge;

    newedge = (struct Edge *)(void *)getfree(ctx, &efl,err);
    if (*err)
        return(newedge);

    newedge -> reg[F_le] = s1;
    newedge -> reg[F_re] = s2;
    ref(ctx, s1); 
    ref(ctx, s2);
    newedge -> ep[F_le] = (struct Site *) NULL;
    newedge -> ep[F_re] = (struct Site *) NULL;

    dx = s2->coord.x - s1->coord.x;
    dy = s2->coord.y - s1->coord.y;
    adx = dx>0 ? dx : -dx;
    ady = dy>0 ? dy : -dy;
    newedge -> c = s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5;
    if (adx>ady) {
        newedge -> a = 1.0; newedge -> b = dy/dx; newedge -> c /= dx;
    }
    else {
        newedge -> b = 1.0; newedge -> a = dx/dy; newedge -> c /= dy;
    } 
    newedge -> edgenbr = ctx->F_nedges;
    out_bisector(ctx, newedge);
    ctx->F_nedges += 1;
    return(newedge);
}

struct Site *intersect(TDAContext *ctx, struct Halfedge *el1,struct Halfedge *el2,int *err)
{
    struct  Edge *e1,*e2, *e;
    struct  Halfedge *el;
    double d, xint, yint;
    int right_of_site;
    struct Site *v;

    e1 = el1 -> ELedge;
    e2 = el2 -> ELedge;
    if(e1 == (struct Edge*)NULL || e2 == (struct Edge*)NULL) 
            return ((struct Site *) NULL);
    if (e1->reg[F_re] == e2->reg[F_re])
        return ((struct Site *) NULL);

    d = e1->a * e2->b - e1->b * e2->a;

    /* tda_out("intersect: d=%g\n", d); */

    if (-1.0e-17 < d && d < 1.0e-17) {
        return ((struct Site *) NULL);
    } 
    xint = (e1->c*e2->b - e2->c*e1->b)/d;
    yint = (e2->c*e1->a - e1->c*e2->a)/d;

    if ((e1->reg[F_re]->coord.y < e2->reg[F_re]->coord.y) ||
        (e1->reg[F_re]->coord.y == e2->reg[F_re]->coord.y &&
            e1->reg[F_re]->coord.x < e2->reg[F_re]->coord.x)) {
        el = el1; e = e1;
    }
    else {
        el = el2; e = e2;
    } 
    right_of_site = xint >= e -> reg[F_re] -> coord.x;
    if ((right_of_site && el -> ELpm == F_le) ||
           (!right_of_site && el -> ELpm == F_re))
        return ((struct Site *) NULL);

    v = (struct Site *)(void *)getfree(ctx, &sfl,err);
    if (*err)
        return(v);

    v -> refcnt = 0;
    v -> coord.x = xint;
    v -> coord.y = yint;
    return(v);
}

/* returns 1 if p is to right of halfedge e */

int right_of(TDAContext *ctx, struct Halfedge *el,struct Point *p)
{
    (void)ctx;        /* unused: the signature is shared */
    struct Edge *e;
    struct Site *topsite;
    int right_of_site, above, fast;
    double dxp, dyp, dxs, t1, t2, t3, yl;

    e = el -> ELedge;
    topsite = e -> reg[F_re];
    right_of_site = p -> x > topsite -> coord.x;
    if(right_of_site && el -> ELpm == F_le)
        return(F_re);
    if(!right_of_site && el -> ELpm == F_re)
        return (F_le);

    if (e->a == 1.0) {
        dyp = p->y - topsite->coord.y;
        dxp = p->x - topsite->coord.x;
        fast = 0;
        if ((!right_of_site & (e->b < 0.0)) | (right_of_site & (e->b >= 0.0))) {
            above = dyp>= e->b*dxp; 
                fast = above;
        }
        else {
            above = p->x + p->y * e->b > e-> c;
                if(e->b<0.0)
                above = !above;
                if (!above)
                fast = 1;
        } 
        if (!fast) {
            dxs = topsite->coord.x - (e->reg[F_le])->coord.x;
                above = e->b * (dxp*dxp - dyp*dyp) <
                dxs * dyp * (1.0 + 2.0 * dxp/ dxs + e->b*e->b);
                if(e->b < 0.0)
                above = !above;
        } 
    }
    else {    /*e->b==1.0 */
        yl = e->c - e->a*p->x;
        t1 = p->y - yl;
        t2 = p->x - topsite->coord.x;
        t3 = yl - topsite->coord.y;
        above = t1*t1 > t2*t2 + t3*t3;
    } 
    return (el->ELpm == F_le ? above : !above);
}

void endpoint(TDAContext *ctx, struct Edge *e,int lr,struct Site *s)
{
    e -> ep[lr] = s;
    ref(ctx, s);
    if(e -> ep[F_re - lr]== (struct Site *) NULL)
        return;
    out_ep(ctx, e);
    deref(ctx, e->reg[F_le]);
    deref(ctx, e->reg[F_re]);
    makefree(ctx, e,&efl);
}

float dist(TDAContext *ctx, struct Site *s,struct Site *t)
{
    (void)ctx;        /* unused: the signature is shared */
    double dx,dy;
    dx = s->coord.x - t->coord.x;
    dy = s->coord.y - t->coord.y;
    return ((float)(sqrt(dx*dx + dy*dy)));
}

void makevertex(TDAContext *ctx, struct Site *v)
{
    v -> sitenbr = ctx->F_nvertices;
    ctx->F_nvertices += 1;
    out_vertex(ctx, v);
}

void deref(TDAContext *ctx, struct Site *v)
{
    v -> refcnt -= 1;
    if (v -> refcnt == 0 )  {
        makefree(ctx, v,&sfl);
    }
}

void ref(TDAContext *ctx, struct Site *v)
{
    (void)ctx;        /* unused: the signature is shared */
    v -> refcnt += 1;
}

/* ======================================================================== */
/* code from heap.c                                                         */
/* ======================================================================== */

void PQinsert(TDAContext *ctx, struct Halfedge *he,struct Site *v,float  offset)
{
    struct Halfedge *last, *next;
    he -> vertex = v;
    ref(ctx, v);
    he -> ystar = (double)((v -> coord.y)) + (double)(offset);
    last = &PQhash[PQbucket(ctx, he)];
    while ((next = last -> PQnext) != (struct Halfedge *) NULL &&
           (he -> ystar  > next -> ystar  ||
           (he -> ystar == next -> ystar && v -> coord.x > next->vertex->coord.x))) {
        last = next;
    } 
    he -> PQnext = last -> PQnext; 
    last -> PQnext = he;
    ctx->PQcount += 1;
}

void PQdelete(TDAContext *ctx, struct Halfedge *he)
{
    struct Halfedge *last;

    if(he ->  vertex != (struct Site *) NULL) {
        last = &PQhash[PQbucket(ctx, he)];
        while (last -> PQnext != he)
            last = last -> PQnext;
        last -> PQnext = he -> PQnext;
        ctx->PQcount -= 1;
        deref(ctx, he -> vertex);
        he -> vertex = (struct Site *) NULL;
    } 
}

int PQbucket(TDAContext *ctx, struct Halfedge *he)
{
    int bucket;

    if  ((double)((he->ystar)) < (double)(F_ymin))
        bucket = 0;
    else if ((double)((he->ystar)) >= (double)(F_ymax))
        bucket = ctx->PQhashsize-1;
    else
               bucket = (int)((double)((((he->ystar - (double)(F_ymin))))) / (double)(F_deltay) * ctx->PQhashsize);

    if (bucket<0)
        bucket = 0;
    if (bucket>=ctx->PQhashsize)
        bucket = ctx->PQhashsize-1 ;
    if (bucket < ctx->PQmin)
        ctx->PQmin = bucket;
    return(bucket);
}

int PQempty(TDAContext *ctx)
{
    return(ctx->PQcount==0);
}

struct Point PQ_min(TDAContext *ctx)
{
    struct Point answer;

    while(PQhash[ctx->PQmin].PQnext == (struct Halfedge *)NULL) {
        ctx->PQmin += 1;
    } 
    answer.x = PQhash[ctx->PQmin].PQnext -> vertex -> coord.x;
    answer.y = PQhash[ctx->PQmin].PQnext -> ystar;
    return (answer);
}

struct Halfedge *PQextractmin(TDAContext *ctx)
{
    struct Halfedge *curr;

    curr = PQhash[ctx->PQmin].PQnext;
    PQhash[ctx->PQmin].PQnext = curr -> PQnext;
    ctx->PQcount -= 1;
    return(curr);
}

int PQinitialize(TDAContext *ctx)
{
    int i,err; struct Point;

    ctx->PQcount = 0;
    ctx->PQmin = 0;
    ctx->PQhashsize = 4 * ctx->F_sqrt_nsites;
    PQhash = (struct Halfedge *)(void *)myalloc(ctx,(unsigned int)((size_t)ctx->PQhashsize * sizeof *PQhash),&err);
    if (err)
        return(-1);

    for(i=0; i<ctx->PQhashsize; i+=1)
        PQhash[i].PQnext = (struct Halfedge *)NULL;
    return(0);
}

/* ======================================================================== */
/* code from memory.c                                                       */
/* ======================================================================== */

void freeinit(TDAContext *ctx, struct Freelist *fl,int size)
{
    (void)ctx;        /* unused: the signature is shared */
    fl -> head = NULL;
    fl -> nodesize = size;
}

char *getfree(TDAContext *ctx, struct Freelist *fl,int *err)
{
    int i; char *t;

    if(fl->head == NULL) {
        t = myalloc(ctx,(unsigned int)(ctx->F_sqrt_nsites * fl->nodesize),err);
        if (*err)
            return(t);

        for (i = 0; i < ctx->F_sqrt_nsites; i += 1)
                makefree(ctx, t+i*fl->nodesize, fl);
    }
    t = fl -> head;
    memcpy(&fl->head, t, sizeof fl->head);
    return(t);
}

void makefree(TDAContext *ctx, void *curr,struct Freelist *fl)
{
    (void)ctx;        /* unused: the signature is shared */
    memcpy(curr, &fl->head, sizeof fl->head);
    fl -> head = curr;
}

char *myalloc(TDAContext *ctx, unsigned n,int *err)
{
    char *t = NULL;

    ctx->F_NMReq += 1;
    if (ctx->F_ACnt >= ctx->F_MaxAlloc) {
        printf1(ctx, "Error: insufficient memory (due to counting limits).\n");
        *err = 1;
        return(t);
    }
    if ((t = malloc(n)) == (char *)0) {
        printf1(ctx, "Error: insufficient memory (allocated: %d).\n",ctx->F_total_alloc);
        *err = 1;
        return(t);
    } 
    ctx->F_ASize[ctx->F_ACnt] = n;
    ctx->F_APtr[ctx->F_ACnt] = t;      
    ctx->F_ACnt++;
    memrq(ctx,(int)(n),sizeof(char));
    ctx->F_total_alloc = ctx->F_total_alloc + (int)n;
    *err = 0;
    return(t);
}

void myalloc_free(TDAContext *ctx)
{
    int i;

    for (i = 0; i < ctx->F_ACnt; ++i) {
        free(ctx->F_APtr[i]);
        memrq(ctx,(int)(-ctx->F_ASize[i]),sizeof(char));
    }
    ctx->F_ACnt = 0;
}

/* =###==================================================================== */
/* modified code from output.c                                              */
/* ======================================================================== */

void out_bisector(TDAContext *ctx, struct Edge *e)
{
    (void)ctx; (void)e;        /* unused: the signature is shared */

    /********************************************
    if(!F_triangulate) {
        tda_out("l %f %f %f\n", e->a, e->b, e->c);
    ********************************************/

    /*****************************************************************
    if(debug)
            tda_out("line(%d) %gx+%gy=%g, bisecting %d %d\n", e->edgenbr,
        e->a, e->b, e->c, e->reg[F_le]->sitenbr, e->reg[F_re]->sitenbr);
    *****************************************************************/
}

void out_ep(TDAContext *ctx, struct Edge *e)
{
    int i,j;
    double x1,y1,x2,y2;
    struct Site *s1,*s2;

    ctx->F_NNEdges += 1;

    if (ctx->F_triangulate == 0) {
        /**************************
        tda_out("e %d", e->edgenbr);
        tda_out(" %d ", e->ep[F_le] != (struct Site *)NULL ? e->ep[F_le]->sitenbr : -1);
        tda_out("%d\n", e->ep[F_re] != (struct Site *)NULL ? e->ep[F_re]->sitenbr : -1);
        **************************************************************************/

        if (e -> a == 1.0 && e ->b >= 0.0) {
            s1 = e -> ep[F_re];
            s2 = e -> ep[F_le];
        }
        else {
            s1 = e -> ep[F_le];
            s2 = e -> ep[F_re];
        } 
        if (s1 != (struct Site *)NULL)  
            i = 1 + s1->sitenbr;
        else
            i = -1;

        if (s2 != (struct Site *)NULL)  
            j = 1 + s2->sitenbr;
        else
            j = -1;
     
        if(e -> a == 1.0) {
            y1 = ctx->PMRECYMin;
            if (s1 != (struct Site *)NULL && s1->coord.y > ctx->PMRECYMin)
                y1 = s1->coord.y;
            if (y1 > ctx->PMRECYMax)
                return;

            x1 = e -> c - e -> b * y1;
            y2 = ctx->PMRECYMax;
            if (s2 != (struct Site *)NULL && s2->coord.y < ctx->PMRECYMax) 
                y2 = s2->coord.y;
            if (y2 < ctx->PMRECYMin)
                return;
            x2 = e -> c - e -> b * y2;
            if ((x1 > ctx->PMRECXMax && x2 > ctx->PMRECXMax) ||
                (x1 < ctx->PMRECXMin && x2 < ctx->PMRECXMin))
                return;

            if (x1 > ctx->PMRECXMax) {
                x1 = ctx->PMRECXMax;
                y1 = (e -> c - x1) / e -> b;
            } 
            if (x1 < ctx->PMRECXMin) {
                x1 = ctx->PMRECXMin;
                y1 = (e -> c - x1) / e -> b;
            } 
            if (x2 > ctx->PMRECXMax) {
                x2 = ctx->PMRECXMax;
                y2 = (e -> c - x2) / e -> b;
            } 
            if (x2 < ctx->PMRECXMin) {
                x2 = ctx->PMRECXMin;
                y2 = (e -> c - x2) / e -> b;
            } 
        }
        else {
            x1 = ctx->PMRECXMin;
            if (s1 != (struct Site *)NULL && s1->coord.x > ctx->PMRECXMin) 
                x1 = s1->coord.x;
            if (x1 > ctx->PMRECXMax)
                return;
            y1 = e -> c - e -> a * x1;
            x2 = ctx->PMRECXMax;
            if (s2 != (struct Site *)NULL && s2->coord.x < ctx->PMRECXMax) 
                x2 = s2->coord.x;
            if (x2 < ctx->PMRECXMin)
                return;
            y2 = e -> c - e -> a * x2;
            if ((y1 > ctx->PMRECYMax && y2 > ctx->PMRECYMax) ||
                (y1 < ctx->PMRECYMin && y2 < ctx->PMRECYMin))
                return;
            if (y1 > ctx->PMRECYMax) {
                y1 = ctx->PMRECYMax;
                x1 = (e -> c - y1) / e -> a;
            } 
            if (y1 < ctx->PMRECYMin) {
                y1 = ctx->PMRECYMin;
                x1 = (e -> c - y1) / e -> a;
            } 
            if (y2 > ctx->PMRECYMax) {
                y2 = ctx->PMRECYMax;
                x2 = (e -> c - y2) / e -> a;
            } 
            if (y2 < ctx->PMRECYMin) {
                y2 = ctx->PMRECYMin;
                x2 = (e -> c - y2) / e -> a;
            } 
        } 
        if (fabs(x1 - x2) <= ctx->EPSI1 && fabs(y1 - y2) <= ctx->EPSI1)
            return;

        fprintf(ctx->PMFd,"%8d 2 %6d %6d %6d\n",++ctx->F_sdid,2,i,j);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,x1);               
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y1);               
        fprintf(ctx->PMFd,"\n");
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,x2);               
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,y2);               
        fprintf(ctx->PMFd,"\n");
        ctx->F_nrec += 3;
    } 
}

void out_vertex(TDAContext *ctx, struct Site *v)
{
    int n;

    ctx->F_NNodes += 1;

    if (!ctx->F_triangulate) {
        /*******************************************
        printf ("v %f %f\n", v->coord.x, v->coord.y);
        *******************************************/

        n = 1 + v->sitenbr;
        fprintf(ctx->PMFd,"%8d 1 %6d %6d %6d\n",++ctx->F_sdid,1,n,n);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,v->coord.x);
        rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,v->coord.y);
        fprintf(ctx->PMFd,"\n");
        ctx->F_nrec += 2;
    }

    /***********************************************************************
    if(debug)
        tda_out("vertex(%d) at %f %f\n", v->sitenbr, v->coord.x, v->coord.y);
    ***********************************************************************/
}

void out_triple(TDAContext *ctx, struct Site *s1,struct Site *s2,struct Site *s3)
{
    int n1,n2,n3;

    ctx->F_NTriangles += 1;

    if (ctx->F_triangulate) {
        /********************************************************* 
        tda_out("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
        **********************************************************/
    
        n1 = 1 + (s1->sitenbr);
        n2 = 1 + (s2->sitenbr);
        n3 = 1 + (s3->sitenbr);

        if (ctx->PMOPT == 3) {
            fprintf(ctx->PMFd,"%8d 3 %6d %6d %6d %6d\n",++ctx->F_sdid,3,n1,n2,n3);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s1->coord.x);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s1->coord.y);
            fprintf(ctx->PMFd,"\n");
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s2->coord.x);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s2->coord.y);
            fprintf(ctx->PMFd,"\n");
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s3->coord.x);
            rt_fprintf_d(ctx, ctx->PMFd,ctx->PMFmtS,s3->coord.y);
            fprintf(ctx->PMFd,"\n");
            ctx->F_nrec += 4;
        }
        else if (ctx->PMOPT == 4) {
            fprintf(ctx->PMFd,"%6d %6d\n",n1,n2);
            fprintf(ctx->PMFd,"%6d %6d\n",n2,n3);
            fprintf(ctx->PMFd,"%6d %6d\n",n3,n1);
            ctx->F_nrec += 3;
        }
    }

    /**************************************************************
    if (debug)
        tda_out("circle through left=%d right=%d bottom=%d\n", 
                            s1->sitenbr, s2->sitenbr, s3->sitenbr);
    **************************************************************/
}








void tda_reset_t_tri(void)
{
    F_xmin = F_xmax = F_ymin = F_ymax = F_deltax = F_deltay = 0.0f;
    pxmin = pxmax = pymin = pymax = cradius = 0.0f;
}
