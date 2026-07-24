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

/* ------------------------------------------------------------------------ */
/* Global variables and definitions                                         */

#define F_DELETED -2
#define F_le 0
#define F_re 1
           
struct  Freenode    {
    struct  Freenode    *nextfree;
};
struct  Freelist    {
    struct  Freenode    *head;
    int         nodesize;
};
struct Point    {
    float x,y;
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

int F_nrec = 0;                 /* records written to output file */
int F_sdid = 0;                 /* ID for spatial output file */

int F_MaxAlloc = 0;             /* maximal allocations */
unsigned int *F_ASize;
char **F_APtr;
int F_ASizeA = 0;
int F_APtrA = 0;
int F_ACnt = 0;
int F_total_alloc = 0;
int F_NMReq = 0;
 
int F_NTriangles = 0;
int F_NNodes = 0;
int F_NNEdges = 0;

int F_nsites = 0;
int F_siteidx = 0;
int F_sqrt_nsites = 0;
int F_nvertices = 0;
int F_nedges = 0;
int ELhashsize = 0;
int F_ntry = 0;
int F_totalsearch = 0;
int F_triangulate = 0;
int PQhashsize = 0;
int PQcount = 0;
int PQmin = 0;
float F_xmin, F_xmax, F_ymin, F_ymax, F_deltax, F_deltay;
float pxmin, pxmax, pymin, pymax, cradius;

/*  functions in t_tri.c */

int sdvd(void);
int reg_memory(int opt);

/* the following functions are taken from Fortune's C program */

int voronoi(struct Site *(*nextsite)());
struct Site *nextone();
int ELinitialize(void);
struct Halfedge *HEcreate(struct Edge *e,int pm,int *err);
void ELinsert(struct Halfedge *lb,struct Halfedge *new);
struct Halfedge *ELgethash(int b);
struct Halfedge *ELleftbnd(struct Point *p);
void ELdelete(struct Halfedge *he);
struct Halfedge *ELright(struct Halfedge *he);
struct Halfedge *ELleft(struct Halfedge *he);
struct Site *leftreg(struct Halfedge *he);
struct Site *rightreg(struct Halfedge *he);
void geominit(void);
struct Edge *bisect(struct Site *s1,struct Site *s2,int *err);
struct Site *intersect(struct Halfedge *el1,struct Halfedge *el2,int *err);
int right_of(struct Halfedge *el,struct Point *p);
void endpoint(struct Edge *e,int lr,struct Site *s);
float dist(struct Site *s,struct Site *t);
void makevertex(struct Site *v);
void deref(struct   Site *v);
void ref(struct Site *v);
void PQinsert(struct Halfedge *he,struct Site *v,float  offset);
void PQdelete(struct Halfedge *he);
int PQbucket(struct Halfedge *he);
int PQempty(void);
struct Point PQ_min(void);
struct Halfedge *PQextractmin(void);
int PQinitialize(void);
void freeinit(struct    Freelist *fl,int    size);
char *getfree(struct    Freelist *fl,int *err);
void makefree(struct Freenode *curr,struct Freelist *fl);
void out_bisector(struct Edge *e);
void out_ep(struct Edge *e);
void out_vertex(struct Site *v);
void out_triple(struct Site *s1,struct Site *s2,struct Site *s3);
char *myalloc(unsigned n,int *err);
void myalloc_free(void);
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

int sdvd(void)
{
    register int i,k;
    int err,r,nn,dflag,ix,iy;
    double dx,dy;

    err = -1;
    printf1("Voronoi diagram and Delaunay triangulation. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 4,1,1))       /* get parameters */
        goto SDVDFin;

    if (PMFmtF == 0)
        pmfmt(10,4);

    dflag = 0;
    if (PMNV == 2) {
        dflag = 1;
    }
    else if (check_sd(0) == 0)    
        dflag = 2;

    if (dflag == 0) {
        printf1("Error: need a spatial data matrix or two variables on left-hand side.\n");
        goto SDVDFin;
    }
    if (PMOPT > 4)
        PMOPT = 1;
    printf1("Option %d. ",PMOPT);
    if (PMOPT == 3 || PMOPT == 4) {
        printf1("Delaunay triangulation.\n");
        F_triangulate = 1;
    }
    else {
        printf1("Voronoi diagram.\n");
        F_triangulate = 0;
    }
    if (PMMFlg == 0)
        F_MaxAlloc = 10000;
    else
        F_MaxAlloc = PMM;

    if (F_MaxAlloc < 10)
        F_MaxAlloc = 10;

    if (alloc_acxf(NOC + 1))
        goto SDVDFin; 
    if (alloc_acyf(NOC + 1))
        goto SDVDFin; 
    if (alloc_ack(NOC + 1))
        goto SDVDFin; 

    if (dflag == 1) {                   /* variables from data matrix */
        ix = PMVIdx[0];
        iy = PMVIdx[1];

        for (i = 0; i < NOC; ++i) {
            AcXF[i] = (float)get_data(ix,i);
            AcYF[i] = (float)get_data(iy,i);
        }
        nn = NOC;
    }
    else {                              /* spatial data */
        nn = 0;
        for (i = 0; i < NOC; ++i) {
            if ((int)get_data(SDVarSDTyp,i) != 1)
                continue;

            if (sd_getdata(i,0,0,1) < 1)
                goto SDVDFin;
             
            AcXF[nn] = (float)SDVarX[0];
            AcYF[nn] = (float)SDVarY[0];
            nn++;
        }
    }
    if (nn < 3) {
        printf1("Error: need at least three points.\n");
        goto SDVDFin;
    }

    /* calculate Voronoi diagram or Delaunay triangulation */
    /* first get sorted input data into sites */

    F_NMReq = 0;
    if (reg_memory(F_MaxAlloc))
        goto SDVDFin;

    freeinit(&sfl, sizeof *sites);

    sites = (struct Site *)myalloc(nn * sizeof *sites,&r);
    if (r)
        goto SDVDFin;
    F_nsites = nn;
        
    if (sortdp2f(nn,AcYF,AcXF,AcK))
        goto SDVDFin;

    for (i = 0; i < nn; ++i) {
        k = AcK[i];
        sites[i].sitenbr = k;
        sites[i].refcnt = 0;
        sites[i].coord.x = AcXF[k];
        sites[i].coord.y = AcYF[k];
    }
    F_xmin = sites[0].coord.x; 
    F_xmax = sites[0].coord.x;
    for(i = 1; i < F_nsites; ++i) {
        if(sites[i].coord.x < F_xmin)
            F_xmin = sites[i].coord.x;
        if(sites[i].coord.x > F_xmax)
            F_xmax = sites[i].coord.x;
    } 
    F_ymin = sites[0].coord.y;
    F_ymax = sites[F_nsites - 1].coord.y;

    printf1("\nBounding box of input data: %lg, %lg, %lg, %lg\n",
        F_xmin,F_xmax,F_ymin,F_ymax);

    dx = (F_xmax - F_xmin) / 10.0;
    dy = (F_ymax - F_ymin) / 10.0;
        
    if (dx <= EPSI1 || dy <= EPSI1) {
        printf1("Error: all points are (almost) collinear.\n");
        goto SDVDFin;
    }
    if (PMRECFlg == 0) {
        PMRECXMin = F_xmin - dx;
        PMRECYMin = F_ymin - dy;
        PMRECXMax = F_xmax + dx;
        PMRECYMax = F_ymax + dy;
    }
    else {
        PMRECXMin = dmin(PMRECXMin,F_xmin);
        PMRECYMin = dmin(PMRECYMin,F_ymin);
        PMRECXMax = dmax(PMRECXMax,F_xmax);
        PMRECYMax = dmax(PMRECYMax,F_ymax);
    }
    if (F_triangulate == 0) {
        printf1("Bounding box for Voronoi diagram: %lg, %lg, %lg, %lg\n",
                 PMRECXMin,PMRECYMin,PMRECXMax,PMRECYMax);
    } 
    F_NTriangles = F_NNodes = F_NNEdges = 0;
    F_sdid = F_nrec = F_siteidx = 0;

    geominit();
    if (voronoi(nextone))
        goto SDVDFin;

    printf1("\nNumber of data points: %d\n",F_nsites);
    printf1("Number of Voronoi nodes: %d\n",F_NNodes);
    printf1("Number of Voronoi edges: %d\n",F_NNEdges);
    printf1("Number of Delaunay triangles: %d\n",F_NTriangles);
    printf1("%d records written to: %s\n",F_nrec,PMFdName);

    printf1("\nNumber of memory requests: %d\n",F_NMReq);
    err = 0;

SDVDFin:
    myalloc_free();
    reg_memory(0);
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  reg_memory(opt)                                                         */
/*                                                                          */  
/*  if opt = 1 allocate F_ASize and F_APtr for maximal F_MaxAlloc           */
/*  memory requests (by myalloc), otherwise free previously allocated       */
/*  memory. Return 0 if OK, -1 if error.                                    */

int reg_memory(int opt)
{
    if (F_ASizeA > 0) {
        free((char *)F_ASize);
        memrq(-F_ASizeA,sizeof(unsigned int));
        F_ASizeA = 0;
    }
    if (F_APtrA > 0) {
        free((char *)F_APtr);
        memrq(-F_APtrA,sizeof(char *));
        F_APtrA = 0;
    }
    if (opt == 0)
        return(0);

    if (!(F_ASize = (unsigned int *)calloc(F_MaxAlloc,sizeof(unsigned int)))) {
        p_err(-2,1);
        return(-1);
    }
    F_ASizeA = F_MaxAlloc;
    memrq(F_ASizeA,sizeof(unsigned int));

    if (!(F_APtr = (char **)calloc(F_MaxAlloc,sizeof(char *)))) {
        p_err(-2,1);
        return(-1);
    }
    F_APtrA = F_MaxAlloc;
    memrq(F_APtrA,sizeof(char *));
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

int voronoi(struct Site *(*nextsite)())
{
    struct Site *newsite, *bot, *top, *temp, *p;
    struct Site *v;
    struct Point newintstar;
    int pm,err;
    struct Halfedge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
    struct Edge *e;

    if (PQinitialize())
        return(-1);

    bottomsite = (*nextsite)();
    /** out_site(bottomsite); **/
    if (ELinitialize())
        return(-1);

    newsite = (*nextsite)();
    while(1) {
 
        if (!PQempty())
            newintstar = PQ_min();

        if (newsite != (struct Site *)NULL 
               && (PQempty() 
                 || newsite -> coord.y < newintstar.y
                 || (newsite->coord.y == newintstar.y 
                 && newsite->coord.x < newintstar.x))) {
                                                     /* new site is smallest */
            /** out_site(newsite); **/
                lbnd = ELleftbnd(&(newsite->coord));
                rbnd = ELright(lbnd);
                bot = rightreg(lbnd);
                e = bisect(bot, newsite,&err);
            if (err)
                return(-1);

                bisector = HEcreate(e,F_le,&err);
            if (err)
                return(-1);

                ELinsert(lbnd, bisector);

                p = intersect(lbnd, bisector,&err);
            if (err)
                return(-1);

                if (p != (struct Site *) NULL) {
                PQdelete(lbnd);
                        PQinsert(lbnd, p, dist(p,newsite));
            } 
                lbnd = bisector;
                bisector = HEcreate(e,F_re,&err);
            if (err)
                return(-1);

                ELinsert(lbnd, bisector);
                p = intersect(bisector, rbnd,&err);
            if (err)
                return(-1);

                if (p != (struct Site *) NULL) {
                PQinsert(bisector, p, dist(p,newsite)); 
                } 
                newsite = (*nextsite)();    
        }
        else if (!PQempty()) {
                                        /* intersection is smallest */
            lbnd = PQextractmin();
                llbnd = ELleft(lbnd);
                rbnd = ELright(lbnd);
                rrbnd = ELright(rbnd);
                bot = leftreg(lbnd);
                top = rightreg(rbnd);
                out_triple(bot, top, rightreg(lbnd));
                v = lbnd->vertex;
                makevertex(v);
                endpoint(lbnd->ELedge,lbnd->ELpm,v);
                endpoint(rbnd->ELedge,rbnd->ELpm,v);
                ELdelete(lbnd); 
                PQdelete(rbnd);
                ELdelete(rbnd); 
                pm = F_le;
                if (bot->coord.y > top->coord.y) {
                    temp = bot; bot = top; top = temp; pm = F_re;
            }
                e = bisect(bot,top,&err);
            if (err)
                return(-1);

                bisector = HEcreate(e,pm,&err);
            if (err)
                return(-1);

                ELinsert(llbnd, bisector);
                endpoint(e, F_re - pm, v);
                deref(v);
                p = intersect(llbnd, bisector,&err);
            if (err)
                return(-1);

                if (p != (struct Site *) NULL) {
                PQdelete(llbnd);
                        PQinsert(llbnd, p, dist(p,bot));
            } 
                p = intersect(bisector, rrbnd,&err);
            if (err)
                return(-1);

                if (p != (struct Site *) NULL) {
                PQinsert(bisector, p, dist(p,bot));
                } 
        }
        else
            break;
    } 
    for(lbnd=ELright(ELleftend); lbnd != ELrightend; lbnd=ELright(lbnd)) {
        e = lbnd -> ELedge;
            out_ep(e);
    } 
    return(0);
}

/* return a single in-storage site */

struct Site *nextone()
{
    struct Site *s;
    for ( ; F_siteidx < F_nsites; F_siteidx += 1) {
        if (F_siteidx == 0 || sites[F_siteidx].coord.x != sites[F_siteidx-1].coord.x 
                               || sites[F_siteidx].coord.y != sites[F_siteidx-1].coord.y) {
            F_siteidx += 1;
                return (&sites[F_siteidx - 1]);
        } 
    } 
    return((struct Site *)NULL);
}


/* ======================================================================== */
/* code from edgelist.c                                                     */
/* ======================================================================== */

int ELinitialize(void)
{
    int i,err;

    freeinit(&hfl, sizeof **ELhash);
    ELhashsize = 2 * F_sqrt_nsites;
    ELhash = (struct Halfedge **)myalloc(sizeof *ELhash * ELhashsize,&err);
    if (err)
        return(-1);

    for(i=0; i<ELhashsize; i +=1)
        ELhash[i] = (struct Halfedge *)NULL;
    ELleftend = HEcreate((struct Edge *)NULL,0,&err);
    if (err)
        return(-1);

    ELrightend = HEcreate((struct Edge *)NULL,0,&err);
    if (err)
        return(-1);

    ELleftend -> ELleft = (struct Halfedge *)NULL;
    ELleftend -> ELright = ELrightend;
    ELrightend -> ELleft = ELleftend;
    ELrightend -> ELright = (struct Halfedge *)NULL;
    ELhash[0] = ELleftend;
    ELhash[ELhashsize-1] = ELrightend;
    return(0);
}

struct Halfedge *HEcreate(struct Edge *e,int pm,int *err)
{
    struct Halfedge *answer;
    answer = (struct Halfedge *)getfree(&hfl,err);
    if (*err)
        return(answer);

    answer -> ELedge = e;
    answer -> ELpm = pm;
    answer -> PQnext = (struct Halfedge *) NULL;
    answer -> vertex = (struct Site *) NULL;
    answer -> ELrefcnt = 0;
    return(answer);
}

void ELinsert(struct Halfedge *lb,struct Halfedge *new)
{
    new -> ELleft = lb;
    new -> ELright = lb -> ELright;
    (lb -> ELright) -> ELleft = new;
    lb -> ELright = new;
}

/* Get entry from hash table, pruning any deleted nodes */

struct Halfedge *ELgethash(int b)
{
    struct Halfedge *he;

    if (b<0 || b>=ELhashsize)
        return((struct Halfedge *) NULL);
    he = ELhash[b]; 
    if (he == (struct Halfedge *) NULL || 
        he -> ELedge != (struct Edge *)F_DELETED )
        return (he);

    /* Hash table points to deleted half edge.  Patch as necessary. */

    ELhash[b] = (struct Halfedge *) NULL;
    if ((he -> ELrefcnt -= 1) == 0)                  
        makefree((struct Freenode *)he, &hfl);
    return ((struct Halfedge *) NULL);
}   

struct Halfedge *ELleftbnd(struct Point *p)
{
    int i, bucket;
    struct Halfedge *he;

    /* Use hash table to get close to desired halfedge */

    bucket = (p->x - F_xmin) / F_deltax * ELhashsize;
    if (bucket<0)
        bucket =0;
    if(bucket>=ELhashsize)
        bucket = ELhashsize - 1;
    he = ELgethash(bucket);
    if(he == (struct Halfedge *) NULL) {
        for(i=1; 1 ; i += 1) {
            if ((he=ELgethash(bucket-i)) != (struct Halfedge *) NULL)
                break;
                if ((he=ELgethash(bucket+i)) != (struct Halfedge *) NULL)
                break;
            } 
            F_totalsearch += i;
    } 
    F_ntry += 1;

    /* Now search linear list of halfedges for the corect one */

    if (he==ELleftend  || (he != ELrightend && right_of(he,p))) {
        do {he = he -> ELright;} while (he!=ELrightend && right_of(he,p));
         he = he -> ELleft;
    }
    else 
        do {he = he -> ELleft;} while (he!=ELleftend && !right_of(he,p));

    /* Update hash table and reference counts */

    if(bucket > 0 && bucket <ELhashsize-1) {
        if(ELhash[bucket] != (struct Halfedge *) NULL) 
                    ELhash[bucket] -> ELrefcnt -= 1;
            ELhash[bucket] = he;
            ELhash[bucket] -> ELrefcnt += 1;
    } 
    return (he);
}
    
/* This delete routine can't reclaim node, since pointers from hash
   table may be present. */

void ELdelete(struct Halfedge *he)
{
    (he -> ELleft) -> ELright = he -> ELright;
    (he -> ELright) -> ELleft = he -> ELleft;
    he -> ELedge = (struct Edge *)F_DELETED;
}

struct Halfedge *ELright(struct Halfedge *he)
{
    return (he -> ELright);
}

struct Halfedge *ELleft(struct Halfedge *he)
{
    return (he -> ELleft);
}

struct Site *leftreg(struct Halfedge *he)
{
    if(he -> ELedge == (struct Edge *)NULL)
        return(bottomsite);
    return( he -> ELpm == F_le ? 
                he -> ELedge -> reg[F_le] : he -> ELedge -> reg[F_re]);
}

struct Site *rightreg(struct Halfedge *he)
{
    if(he -> ELedge == (struct Edge *)NULL)
        return(bottomsite);
    return( he -> ELpm == F_le ? 
                he -> ELedge -> reg[F_re] : he -> ELedge -> reg[F_le]);
}

/* ======================================================================== */
/* code from geometry.c                                                     */
/* ======================================================================== */

void geominit(void)
{
    struct Edge e;
    float sn;

    freeinit(&efl, sizeof e);
    F_nvertices = 0;
    F_nedges = 0;
    sn = (float)(F_nsites + 4);
    F_sqrt_nsites = (int)sqrt(sn);
    F_deltay = F_ymax - F_ymin;
    F_deltax = F_xmax - F_xmin;
}

struct Edge *bisect(struct Site *s1,struct Site *s2,int *err)
{
    double dx,dy,adx,ady;
    struct Edge *newedge;

    newedge = (struct Edge *)getfree(&efl,err);
    if (*err)
        return(newedge);

    newedge -> reg[F_le] = s1;
    newedge -> reg[F_re] = s2;
    ref(s1); 
    ref(s2);
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
    newedge -> edgenbr = F_nedges;
    out_bisector(newedge);
    F_nedges += 1;
    return(newedge);
}

struct Site *intersect(struct Halfedge *el1,struct Halfedge *el2,int *err)
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

    /* printf("intersect: d=%g\n", d); */

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

    v = (struct Site *)getfree(&sfl,err);
    if (*err)
        return(v);

    v -> refcnt = 0;
    v -> coord.x = xint;
    v -> coord.y = yint;
    return(v);
}

/* returns 1 if p is to right of halfedge e */

int right_of(struct Halfedge *el,struct Point *p)
{
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

void endpoint(struct Edge *e,int lr,struct Site *s)
{
    e -> ep[lr] = s;
    ref(s);
    if(e -> ep[F_re - lr]== (struct Site *) NULL)
        return;
    out_ep(e);
    deref(e->reg[F_le]);
    deref(e->reg[F_re]);
    makefree((struct Freenode *)e,&efl);
}

float dist(struct Site *s,struct Site *t)
{
    float dx,dy;
    dx = s->coord.x - t->coord.x;
    dy = s->coord.y - t->coord.y;
    return(sqrt(dx*dx + dy*dy));
}

void makevertex(struct Site *v)
{
    v -> sitenbr = F_nvertices;
    F_nvertices += 1;
    out_vertex(v);
}

void deref(struct Site *v)
{
    v -> refcnt -= 1;
    if (v -> refcnt == 0 )  {
        makefree((struct Freenode *)v,&sfl);
    }
}

void ref(struct Site *v)
{
    v -> refcnt += 1;
}

/* ======================================================================== */
/* code from heap.c                                                         */
/* ======================================================================== */

void PQinsert(struct Halfedge *he,struct Site *v,float  offset)
{
    struct Halfedge *last, *next;
    he -> vertex = v;
    ref(v);
    he -> ystar = v -> coord.y + offset;
    last = &PQhash[PQbucket(he)];
    while ((next = last -> PQnext) != (struct Halfedge *) NULL &&
           (he -> ystar  > next -> ystar  ||
           (he -> ystar == next -> ystar && v -> coord.x > next->vertex->coord.x))) {
        last = next;
    } 
    he -> PQnext = last -> PQnext; 
    last -> PQnext = he;
    PQcount += 1;
}

void PQdelete(struct Halfedge *he)
{
    struct Halfedge *last;

    if(he ->  vertex != (struct Site *) NULL) {
        last = &PQhash[PQbucket(he)];
        while (last -> PQnext != he)
            last = last -> PQnext;
        last -> PQnext = he -> PQnext;
        PQcount -= 1;
        deref(he -> vertex);
        he -> vertex = (struct Site *) NULL;
    } 
}

int PQbucket(struct Halfedge *he)
{
    int bucket;

    if  (he->ystar < F_ymin)
        bucket = 0;
    else if (he->ystar >= F_ymax)
        bucket = PQhashsize-1;
    else
               bucket = (he->ystar - F_ymin) / F_deltay * PQhashsize;

    if (bucket<0)
        bucket = 0;
    if (bucket>=PQhashsize)
        bucket = PQhashsize-1 ;
    if (bucket < PQmin)
        PQmin = bucket;
    return(bucket);
}

int PQempty(void)
{
    return(PQcount==0);
}

struct Point PQ_min(void)
{
    struct Point answer;

    while(PQhash[PQmin].PQnext == (struct Halfedge *)NULL) {
        PQmin += 1;
    } 
    answer.x = PQhash[PQmin].PQnext -> vertex -> coord.x;
    answer.y = PQhash[PQmin].PQnext -> ystar;
    return (answer);
}

struct Halfedge *PQextractmin(void)
{
    struct Halfedge *curr;

    curr = PQhash[PQmin].PQnext;
    PQhash[PQmin].PQnext = curr -> PQnext;
    PQcount -= 1;
    return(curr);
}

int PQinitialize(void)
{
    int i,err; struct Point *s;

    PQcount = 0;
    PQmin = 0;
    PQhashsize = 4 * F_sqrt_nsites;
    PQhash = (struct Halfedge *)myalloc(PQhashsize * sizeof *PQhash,&err);
    if (err)
        return(-1);

    for(i=0; i<PQhashsize; i+=1)
        PQhash[i].PQnext = (struct Halfedge *)NULL;
    return(0);
}

/* ======================================================================== */
/* code from memory.c                                                       */
/* ======================================================================== */

void freeinit(struct Freelist *fl,int size)
{
    fl -> head = (struct Freenode *) NULL;
    fl -> nodesize = size;
}

char *getfree(struct Freelist *fl,int *err)
{
    int i; struct Freenode *t;

    if(fl->head == (struct Freenode *) NULL) {
        t =  (struct Freenode *)myalloc(F_sqrt_nsites * fl->nodesize,err);
        if (*err)
            return((char *)t);

        for (i = 0; i < F_sqrt_nsites; i += 1)  
                makefree((struct Freenode *)((char *)t+i*fl->nodesize), fl);
    } 
    t = fl -> head;
    fl -> head = (fl -> head) -> nextfree;
    return((char *)t);
}

void makefree(struct Freenode *curr,struct Freelist *fl)
{
    curr -> nextfree = fl -> head;
    fl -> head = curr;
}

char *myalloc(unsigned n,int *err)
{
    char *t;

    F_NMReq += 1;
    if (F_ACnt >= F_MaxAlloc) {
        printf1("Error: insufficient memory (due to counting limits).\n");
        *err = 1;
        return(t);
    }
    if ((t = malloc(n)) == (char *)0) {
        printf1("Error: insufficient memory (allocated: %d).\n",F_total_alloc);
        *err = 1;
        return(t);
    } 
    F_ASize[F_ACnt] = n;
    F_APtr[F_ACnt] = t;      
    F_ACnt++;
    memrq(n,sizeof(char));
    F_total_alloc += n;
    *err = 0;
    return(t);
}

void myalloc_free(void)
{
    int i;

    for (i = 0; i < F_ACnt; ++i) {
        free(F_APtr[i]);
        memrq(-F_ASize[i],sizeof(char));
    }
    F_ACnt = 0;
}

/* =###==================================================================== */
/* modified code from output.c                                              */
/* ======================================================================== */

void out_bisector(struct Edge *e)
{

    /********************************************
    if(!F_triangulate) {
        printf("l %f %f %f\n", e->a, e->b, e->c);
    ********************************************/

    /*****************************************************************
    if(debug)
            printf("line(%d) %gx+%gy=%g, bisecting %d %d\n", e->edgenbr,
        e->a, e->b, e->c, e->reg[F_le]->sitenbr, e->reg[F_re]->sitenbr);
    *****************************************************************/
}

void out_ep(struct Edge *e)
{
    int i,j;
    float x1,y1,x2,y2;
    struct Site *s1,*s2;

    F_NNEdges += 1;

    if (F_triangulate == 0) {
        /**************************
        printf("e %d", e->edgenbr);
        printf(" %d ", e->ep[F_le] != (struct Site *)NULL ? e->ep[F_le]->sitenbr : -1);
        printf("%d\n", e->ep[F_re] != (struct Site *)NULL ? e->ep[F_re]->sitenbr : -1);
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
            y1 = PMRECYMin;
            if (s1 != (struct Site *)NULL && s1->coord.y > PMRECYMin)
                y1 = s1->coord.y;
            if (y1 > PMRECYMax)
                return;

            x1 = e -> c - e -> b * y1;
            y2 = PMRECYMax;
            if (s2 != (struct Site *)NULL && s2->coord.y < PMRECYMax) 
                y2 = s2->coord.y;
            if (y2 < PMRECYMin)
                return;
            x2 = e -> c - e -> b * y2;
            if ((x1 > PMRECXMax && x2 > PMRECXMax) ||
                (x1 < PMRECXMin && x2 < PMRECXMin))
                return;

            if (x1 > PMRECXMax) {
                x1 = PMRECXMax;
                y1 = (e -> c - x1) / e -> b;
            } 
            if (x1 < PMRECXMin) {
                x1 = PMRECXMin;
                y1 = (e -> c - x1) / e -> b;
            } 
            if (x2 > PMRECXMax) {
                x2 = PMRECXMax;
                y2 = (e -> c - x2) / e -> b;
            } 
            if (x2 < PMRECXMin) {
                x2 = PMRECXMin;
                y2 = (e -> c - x2) / e -> b;
            } 
        }
        else {
            x1 = PMRECXMin;
            if (s1 != (struct Site *)NULL && s1->coord.x > PMRECXMin) 
                x1 = s1->coord.x;
            if (x1 > PMRECXMax)
                return;
            y1 = e -> c - e -> a * x1;
            x2 = PMRECXMax;
            if (s2 != (struct Site *)NULL && s2->coord.x < PMRECXMax) 
                x2 = s2->coord.x;
            if (x2 < PMRECXMin)
                return;
            y2 = e -> c - e -> a * x2;
            if ((y1 > PMRECYMax && y2 > PMRECYMax) ||
                (y1 < PMRECYMin && y2 < PMRECYMin))
                return;
            if (y1 > PMRECYMax) {
                y1 = PMRECYMax;
                x1 = (e -> c - y1) / e -> a;
            } 
            if (y1 < PMRECYMin) {
                y1 = PMRECYMin;
                x1 = (e -> c - y1) / e -> a;
            } 
            if (y2 > PMRECYMax) {
                y2 = PMRECYMax;
                x2 = (e -> c - y2) / e -> a;
            } 
            if (y2 < PMRECYMin) {
                y2 = PMRECYMin;
                x2 = (e -> c - y2) / e -> a;
            } 
        } 
        if (fabs(x1 - x2) <= EPSI1 && fabs(y1 - y2) <= EPSI1)
            return;

        fprintf(PMFd,"%8d 2 %6d %6d %6d\n",++F_sdid,2,i,j);
        fprintf(PMFd,PMFmtS,x1);               
        fprintf(PMFd,PMFmtS,y1);               
        fprintf(PMFd,"\n");
        fprintf(PMFd,PMFmtS,x2);               
        fprintf(PMFd,PMFmtS,y2);               
        fprintf(PMFd,"\n");
        F_nrec += 3;
    } 
}

void out_vertex(struct Site *v)
{
    int n;

    F_NNodes += 1;

    if (!F_triangulate) {
        /*******************************************
        printf ("v %f %f\n", v->coord.x, v->coord.y);
        *******************************************/

        n = 1 + v->sitenbr;
        fprintf(PMFd,"%8d 1 %6d %6d %6d\n",++F_sdid,1,n,n);
        fprintf(PMFd,PMFmtS,v->coord.x);
        fprintf(PMFd,PMFmtS,v->coord.y);
        fprintf(PMFd,"\n");
        F_nrec += 2;
    }

    /***********************************************************************
    if(debug)
        printf("vertex(%d) at %f %f\n", v->sitenbr, v->coord.x, v->coord.y);
    ***********************************************************************/
}

void out_triple(struct Site *s1,struct Site *s2,struct Site *s3)
{
    int n1,n2,n3;

    F_NTriangles += 1;

    if (F_triangulate) {
        /********************************************************* 
        printf("%d %d %d\n", s1->sitenbr, s2->sitenbr, s3->sitenbr);
        **********************************************************/
    
        n1 = 1 + (s1->sitenbr);
        n2 = 1 + (s2->sitenbr);
        n3 = 1 + (s3->sitenbr);

        if (PMOPT == 3) {
            fprintf(PMFd,"%8d 3 %6d %6d %6d %6d\n",++F_sdid,3,n1,n2,n3);
            fprintf(PMFd,PMFmtS,s1->coord.x);
            fprintf(PMFd,PMFmtS,s1->coord.y);
            fprintf(PMFd,"\n");
            fprintf(PMFd,PMFmtS,s2->coord.x);
            fprintf(PMFd,PMFmtS,s2->coord.y);
            fprintf(PMFd,"\n");
            fprintf(PMFd,PMFmtS,s3->coord.x);
            fprintf(PMFd,PMFmtS,s3->coord.y);
            fprintf(PMFd,"\n");
            F_nrec += 4;
        }
        else if (PMOPT == 4) {
            fprintf(PMFd,"%6d %6d\n",n1,n2);
            fprintf(PMFd,"%6d %6d\n",n2,n3);
            fprintf(PMFd,"%6d %6d\n",n3,n1);
            F_nrec += 3;
        }
    }

    /**************************************************************
    if (debug)
        printf("circle through left=%d right=%d bottom=%d\n", 
                            s1->sitenbr, s2->sitenbr, s3->sitenbr);
    **************************************************************/
}






