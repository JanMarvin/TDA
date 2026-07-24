/****************************************************************************/
/*  t_map                                                                   */
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

#include "tda.h"
#include "t_gen.h"
#include "t_pgen.h"
#include "t_var.h"
#include "t_gdat.h"
#include "t_parm.h"
#include "t_alloc.h"
#include "t_sort.h"
#include "t_gf.h"
#include "t_psf.h"
#include "t_plot.h"
#include "t_plot3.h"
#include "t_sd.h"
#include "t_clip.h"
#include "t_gm.h"
#include "tda_context.h"

/*  functions in t_map.c */

int psetupg(TDAContext *ctx); 
int check_geo(TDAContext *ctx, int opt);
void gd_view(TDAContext *ctx, int opt);
int gd_check_inverse(TDAContext *ctx);
double gd_adlon(TDAContext *ctx, double lon);
double gd_adlat(TDAContext *ctx, double lat);
int gd_proj_p(TDAContext *ctx, double *x,double *y);
int gd_proj_inv(TDAContext *ctx, double *x,double *y);
int check_x_in_region(TDAContext *ctx, double x);
int check_y_in_region(TDAContext *ctx, double y);
int check_inside(TDAContext *ctx, double lon,double lat);
int check_clip(TDAContext *ctx, double x,double y);
double gd_proj_gety(TDAContext *ctx, double lat);
int gd_proj_getxy(TDAContext *ctx, double *x,double *y);
double gd_adjust_lon(TDAContext *ctx, double lon);

int sdpgrat(TDAContext *ctx); 
void mapgrid_bounds(TDAContext *ctx);
void mapgrid_bounds_r(TDAContext *ctx);
void mapgrid_arc(TDAContext *ctx, double x,double y,double a,double b,double r,int gflag);
void mapgrid_lon(TDAContext *ctx, double sc);
void mapgrid_line(TDAContext *ctx, double xa,double ya,double xb,double yb);
void mapgrid_label(TDAContext *ctx, int opt,double lab,double x,double y,double fs,int m);
void mapgrid_lon_azi(TDAContext *ctx, double lon);
void mapgrid_lon_azi_lab0(TDAContext *ctx, double lab,double lon,double fs);
void mapgrid_lon_azi_lab1(TDAContext *ctx, double lab,double lon,double fs);
void mapgrid_lon_azi1(TDAContext *ctx, double lon,double sc);
void mapgrid_lat(TDAContext *ctx);
void mapgrid_lat_a0(TDAContext *ctx, double lat);
void mapgrid_lat_a3(TDAContext *ctx, double lat);

int sdpgeo(TDAContext *ctx); 
int sdpmap(TDAContext *ctx); 
int map_point(TDAContext *ctx, double x,double y,int nc);
int map_point_cyl(TDAContext *ctx, double x,double y);
int map_point_azi(TDAContext *ctx, double x,double y,int nc);

int map_line(TDAContext *ctx, int n,double *x,double *y,int gflag);
int map_line_cyl_s(TDAContext *ctx, int n,double *x,double *y);
void map_find_n(TDAContext *ctx, double x,double y,double *xa,double *ya);
int map_line_cyl_g(TDAContext *ctx, int n,double *x,double *y);
int map_line_cyl_geod(TDAContext *ctx, double lon1,double lat1,double lon2,double lat2);
void map_line_clip(TDAContext *ctx, double xa,double ya,double xb,double yb);
int map_line_azi(TDAContext *ctx, int n,double *x,double *y,int gflag);
int map_line_straight(TDAContext *ctx, double lona,double lata,double lonb,double latb);
int map_geodesic_a(TDAContext *ctx, double lona,double lata,double lonb,double latb);
void err_geodesic(TDAContext *ctx, int n,double lona,double lata,double lonb,double latb);

int map_polygon(TDAContext *ctx, int n,double *x,double *y,int sdid);
void map_polygon_p2(TDAContext *ctx, int n,double *xx,double *yy);
void plsd_arc(TDAContext *ctx, double x0,double y0,double x1,double y1);
int map_polygon_cyl(TDAContext *ctx, int n,double *x,double *y,int sdid);
int map_polygon_cyl_clip(TDAContext *ctx, int n,double *x,double *y,int sdid);
void map_polygon_cyl_pbox(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax);
void map_polygon_cyl_ppol(TDAContext *ctx, int n,double *x,double *y);
int map_polygon_cyl_plist(TDAContext *ctx, int n,double *x,double *y,int *nf,int *d,int sdid);
int map_check_pole(TDAContext *ctx, int n,double *x,double *y);

double geod_distance(TDAContext *ctx, double lon1,double lat1,double lon2,double lat2);

/*--------------------------------------------------------------------------*/
/*  Global variables                                                        */

#define GDPROJ10  10        /* cylindrical (equidistant)                    */
#define GDPROJ11  11        /* cylindrical (equal-area)                     */
#define GDPROJ12  12        /* cylindrical (Mercator)                       */

#define GDPROJ20  20        /* azimuthal (orthographic)                     */


#define GDPROJ1 1           /* parallel projection (rectangular region)     */
#define GDPROJ2 2           /* parallel projection (hemisphere)             */
#define GDPROJ3 3           /* cylindrical projection (equidistant)         */
#define GDPROJ4 4           /* cylindrical projection (Mercator)            */






/* -##--------------------------------------------------------------------- */
/*  psetupg     Set up a geographical coordinate system and projection.     */
/*                                                                          */
/*  psetupg(                                                                */
/*      proj=...,   typ of projection, def. 10                              */
/*                  10  cylindrical (equidistant)                           */
/*                  11  cylindrical (equal-area)                            */
/*                  12  cylindrical (Mercator)                              */
/*                  20  azimuthal (orthographic)                            */
/*                                                                          */
/*      view=...,   direction: lon,lat, default 0,0                         */
/*      region=..., dlon,dlat (no default)                                  */
/*      rhem=...,   size of hemisphere, def. 90                             */
/*      pxlen=...,  width of plot in mm, def. 120                           */
/*      pylen=...,  height of plot in mm, def. depends on proj              */
/*      psorg=...,  origin of Postscript coordinate system, def. 100,100    */
/*      psrot=...,  rotation, def. 0                                        */
/*  );                                                                      */
/*                                                                          */
/*  The range of lon and lat is always restricted in the following way:     */
/*                                                                          */
/*      -180 <= lon <= +180, -90 <= lat <= +90                              */
/*                                                                          */
/*  Additional restrictions depend on the type of projection.               */
/*                                                                          */
/*  The range of dlon and dlat has the following restrictions: Both dlon    */
/*  and dlat must be greater than zero. Additional remarks apply            */
/*  separately for different projections.                                   */ 
/*                                                                          */
/*  PROJ10  cylindrical (equidistant)                                       */
/*                                                                          */
/*      The plot area is rectangular. The allowed range of dlon is:         */
/*                                                                          */
/*          -180 <= dlon <= +180                                            */
/*                                                                          */
/*      The allowed range of dlat is                                        */
/*                                                                          */
/*          -90 <= dlat <= +90                                              */
/*                                                                          */
/*      The plot is restricted to the range of latitudes: -90 to +90.       */
/*      The user coordinate system has its left lower point at              */
/*      lon - dlon  and  dmax(lat - dlat,-90).                              */
/*                                                                          */
/*  PROJ11  cylindrical (equal-area)                                        */
/*                                                                          */  
/*      lon and lat can be specified in the same way as for projection 10.  */
/*      The same restrictions apply.                                        */
/*                                                                          */  
/*  PROJ12  cylindrical (Mercator)                                          */
/*                                                                          */
/*      Additional restriction: -85 <= dlat <= +85.                         */
/*                                                                          */
/*      lon and lat can be specified in the same way as for projection 10.  */
/*      but the plot area is restricted to +/- 89.                          */  
/*                                                                          */
/*  PROJ20  azimuthal (orthographic)                                        */
/*                                                                          */
/*      a) If the region parameter is not used, the projection area will    */
/*         be a circle. By default, it will be a full hemisphere. In two    */
/*         special cases, view = lat,90 and view = lat,-90, the size of the */
/*         hemisphere can be specified with the rhem parameter. Allowed     */
/*         range is                                                         */
/*                                                                          */
/*              0 < rhem <= 90                                              */
/*                                                                          */
/*         The outside parallel of latitude is then taken to be             */
/*         90 - rhem, if view = lon,90, or rhem - 90, if view = lon,-90.    */  
/*                                                                          */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int psetupg(TDAContext *ctx)  
{
   	int	err;
    double x,y,ymin,ymax;
      
    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Set up geographical coordinate system. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,3,0)) {     /* get parameters */
        goto PSGFin;
    }   
    if (ctx->PMPROJ == 0)                /* default projection */
        ctx->PMPROJ = 10;

    if (ctx->PSFFlg == 0) {
        printf1(ctx, "Error: need a PostScript file.\n");
        return(-1);
    }
    ctx->PSPROJ = 0;                 /* make current system invalid */
    ctx->GDLon = ctx->PMViewLon;          /* center of projection */
    if (ctx->GDLon <= -180.0)
        ctx->GDLon = 180.0;
    ctx->GDLat = ctx->PMViewLat;
    ctx->GDPROJA = 0;
           
    switch (ctx->PMPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:

            if (ctx->PMRegionFlg == 0) {
                if (ctx->GDLat == 0) {
                    ctx->PMRegion1 = 180.0;
                    ctx->PMRegion2 =  90.0;
                }   
                else {
                    printf1(ctx, "Error: need region parameter.\n");
                    goto PSGFin;
                }
            }
            if (ctx->PMRegion1 < ctx->EPSI1 || ctx->PMRegion1 > 180.0 || ctx->PMRegion2 < ctx->EPSI1 || ctx->PMRegion2 > 90.0) {   
                printf1(ctx, "Error in region parameter.\n");
                goto PSGFin;
            }
            ctx->GDLonD = ctx->PMRegion1;
            ctx->GDLatD = ctx->PMRegion2;

            if (ctx->PMPROJ == GDPROJ10) {
                ctx->GDPROJ = GDPROJ10;
                printf1(ctx, "Cylindrical projection (equidistant).\n");
                gd_view(ctx, 0);
                ymin = -90.0;
                ymax =  90.0;
            }
            else if (ctx->PMPROJ == GDPROJ11) {
                ctx->GDPROJ = GDPROJ11;
                printf1(ctx, "Cylindrical projection (equal-area).\n");
                gd_view(ctx, 0);
                ymin = -90.0;
                ymax =  90.0;
            }
            if (ctx->PMPROJ == GDPROJ12) {
                ctx->GDPROJ = GDPROJ12;
                printf1(ctx, "Cylindrical projection (Mercator).\n");
                gd_view(ctx, 0);

                if (ctx->GDLat < -85.0 || ctx->GDLat > 85.0) {
                    printf1(ctx, "Error in view (lat) parameter.\n");
                    goto PSGFin;
                }
                ymin = -89.0;
                ymax =  89.0;
            }
            ctx->GDLonMin = ctx->GDLonA = gd_adlon(ctx, ctx->GDLon - ctx->GDLonD);
            ctx->GDLonMax = ctx->GDLonB = gd_adlon(ctx, ctx->GDLon + ctx->GDLonD);

            ctx->GDLatMin = ctx->GDLatA = dmax(ctx, ctx->GDLat - ctx->GDLatD,ymin);
            ctx->GDLatMax = ctx->GDLatB = dmin(ctx, ctx->GDLat + ctx->GDLatD,ymax);
            ctx->GDLatD = (ctx->GDLatB - ctx->GDLatA) / 2.0;

            ctx->GDPX = ctx->GDLon;
            ctx->GDPXA = gd_adlon(ctx, ctx->GDPX - ctx->GDLonD);
            ctx->GDPXB = ctx->GDPXA + 2.0 * ctx->GDLonD;

            ctx->GDPY = gd_proj_gety(ctx, ctx->GDLat);
            ctx->GDPYA = gd_proj_gety(ctx, ctx->GDLatA);
            ctx->GDPYB = gd_proj_gety(ctx, ctx->GDLatB);

            ctx->PXLen = ctx->PMXLen;
            if (ctx->PMYLenFlg)
                ctx->PYLen = ctx->PMYLen;
            else {
                if (ctx->PMPROJ == GDPROJ12)
                    ctx->PYLen = ctx->PXLen;
                else
                    ctx->PYLen = ctx->PXLen * (ctx->GDPYB - ctx->GDPYA) / (ctx->GDPXB - ctx->GDPXA);
            }
            printf1(ctx, "Selected region longitude: %9.4f (%9.4f) %9.4f\n",ctx->GDLonA,ctx->GDLon,ctx->GDLonB);
            printf1(ctx, "Selected region latitude:  %9.4f (%9.4f) %9.4f\n",ctx->GDLatA,ctx->GDLat,ctx->GDLatB);
            break;                          

        case GDPROJ20:
            ctx->GDPROJ = GDPROJ20;
            printf1(ctx, "Azimuthal projection (orthographic).\n");

            if (ctx->PMRegionFlg != 0) {
                ctx->GDPROJA = 0;
                if (ctx->PMRegion1 < ctx->EPSI1 || ctx->PMRegion1 > 90.0 || ctx->PMRegion2 < ctx->EPSI1 || ctx->PMRegion2 > 90.0) {   
                    printf1(ctx, "Error in region parameter.\n");
                    goto PSGFin;
                }
                ctx->GDLonD = ctx->PMRegion1;
                ctx->GDLatD = ctx->PMRegion2;

                printf1(ctx, "View: %g,%g. Rectangular region (%g,%g).\n",ctx->GDLon,ctx->GDLat,ctx->GDLonD,ctx->GDLatD);

                ctx->GDLonA = gd_adlon(ctx, ctx->GDLon - ctx->GDLonD);
                ctx->GDLonB = gd_adlon(ctx, ctx->GDLon + ctx->GDLonD);
                ctx->GDLatA = ctx->GDLat - ctx->GDLatD;
                ctx->GDLatB = ctx->GDLat + ctx->GDLatD;
                if (ctx->GDLatA < -90.0 || ctx->GDLatB > 90.0) {
                    printf1(ctx, "The selected region is not supported.\n");
                    goto PSGFin;
                }
                setup_3proj(ctx, ctx->GDLon,ctx->GDLat);   /* set up projection matrix */

                ctx->GDPX = ctx->GDPY = 0.0;
                ctx->GDPXA = ctx->GDLon - ctx->GDLonD; y = ctx->GDLat;
                gd_proj_p(ctx, &ctx->GDPXA,&y);
                ctx->GDPXB = ctx->GDLon + ctx->GDLonD; y = ctx->GDLat;
                gd_proj_p(ctx, &ctx->GDPXB,&y);

                ctx->GDPYA = ctx->GDLat - ctx->GDLatD; x = ctx->GDLon;
                gd_proj_p(ctx, &x,&ctx->GDPYA);
                ctx->GDPYB = ctx->GDLat + ctx->GDLatD; x = ctx->GDLon;
                gd_proj_p(ctx, &x,&ctx->GDPYB);

/**               
                tda_out("GDPX=%g GDPY=%g\n",GDPX,GDPY);
                tda_out("xa=%g xb=%g\n",GDPXA,GDPXB);
                tda_out("ya=%g yb=%g\n",GDPYA,GDPYB);
**/              
                /* calculate coordinates of corner points and check
                   visibility. */ 

                ctx->GDLonXAYA = ctx->GDPXA;
                ctx->GDLatXAYA = ctx->GDPYA;
                ctx->GDLonXAYB = ctx->GDPXA;
                ctx->GDLatXAYB = ctx->GDPYB;
                ctx->GDLonXBYA = ctx->GDPXB;
                ctx->GDLatXBYA = ctx->GDPYA;
                ctx->GDLonXBYB = ctx->GDPXB;
                ctx->GDLatXBYB = ctx->GDPYB;

                if (gd_proj_inv(ctx, &ctx->GDLonXAYA,&ctx->GDLatXAYA) == 0 ||
                    gd_proj_inv(ctx, &ctx->GDLonXAYB,&ctx->GDLatXAYB) == 0 ||
                    gd_proj_inv(ctx, &ctx->GDLonXBYA,&ctx->GDLatXBYA) == 0 ||
                    gd_proj_inv(ctx, &ctx->GDLonXBYB,&ctx->GDLatXBYB) == 0) {
                    printf1(ctx, "Error: selected region not completely visible.\n");
                    goto PSGFin;
                }                  
                printf1(ctx, "\nSelected region (lon,lat)\n");
                printf1(ctx, "(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    ctx->GDLonXAYB,ctx->GDLatXAYB,ctx->GDLon,ctx->GDLatB,ctx->GDLonXBYB,ctx->GDLatXBYB);
                printf1(ctx, "(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    ctx->GDLonA,ctx->GDLat,ctx->GDLon,ctx->GDLat,ctx->GDLonB,ctx->GDLat);
                printf1(ctx, "(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    ctx->GDLonXAYA,ctx->GDLatXAYA,ctx->GDLon,ctx->GDLatA,ctx->GDLonXBYA,ctx->GDLatXBYA);

                if (ctx->GDLat >= 0.0) {
                    ctx->GDLatMin = ctx->GDLatXAYA;
                    ctx->GDLatMax = ctx->GDLatB;
                    ctx->GDLonMin = ctx->GDLonXAYB;
                    ctx->GDLonMax = ctx->GDLonXBYB;
                }
                else {
                    ctx->GDLatMin = ctx->GDLatA;
                    ctx->GDLatMax = ctx->GDLatXAYB;
                    ctx->GDLonMin = ctx->GDLonXAYA;
                    ctx->GDLonMax = ctx->GDLonXBYA;
                }
                tda_out("GDlatmin=%g %g\n",ctx->GDLatMin,ctx->GDLatMax);
                tda_out("GDlonmin=%g %g\n",ctx->GDLonMin,ctx->GDLonMax);





/***
                x = GDLonXAYA;
                y = GDLatXAYA;
                tda_out("x=%g y=%g ",x,y);
                gd_proj_p(ctx, &x,&y);
                tda_out(" -> %g %g\n",x,y);

                x = GDLonXBYA;
                y = GDLatXBYA;
                tda_out("x=%g y=%g ",x,y);
                gd_proj_p(ctx, &x,&y);
                tda_out(" -> %g %g\n",x,y);

                x = GDLonXAYB;
                y = GDLatXAYB;
                tda_out("x=%g y=%g ",x,y);
                gd_proj_p(ctx, &x,&y);
                tda_out(" -> %g %g\n",x,y);

                x = GDLonXBYB;
                y = GDLatXBYB;
                tda_out("x=%g y=%g ",x,y);
                gd_proj_p(ctx, &x,&y);
                tda_out(" -> %g %g\n",x,y);
**/
               

                x = fabs(ctx->GDPXB - ctx->GDPXA);
                y = fabs(ctx->GDPYB - ctx->GDPYA);

                ctx->PXLen = ctx->PMXLen;
                if (ctx->PMYLenFlg)
                    ctx->PYLen = ctx->PMYLen;
                else 
                    ctx->PYLen = ctx->PXLen * y / x;                   


            }
            else {
                setup_3proj(ctx, ctx->GDLon,ctx->GDLat);   /* set up projection matrix */

                ctx->GDRHem = ctx->PMRHem;
                if (ctx->GDRHem > 90.0)
                    ctx->GDRHem = 90.0;

                if (fabs(ctx->GDLat - 90.0) <= ctx->EPSI1) {  /* view from north pole */
                    ctx->GDPROJA = 1;
                    ctx->GDLat = 90.0;
                    ctx->GDLonMin = -180.0;
                    ctx->GDLonMax =  180.0;
                }
                else if (fabs(ctx->GDLat + 90.0) <= ctx->EPSI1) {  /* view from south pole */
                    ctx->GDPROJA = 2;
                    ctx->GDLat = -90.0;
                    ctx->GDLonMin = -180.0;
                    ctx->GDLonMax =  180.0;
                }
                else {
                    ctx->GDPROJA = 3;
                    ctx->GDRHem = 90.0;
                }
                printf1(ctx, "View: %g,%g. Circular area. Size (rhem) = %g\n\n",ctx->GDLon,ctx->GDLat,ctx->GDRHem);

                ctx->GDLonD = 90.0;
                ctx->GDLonMin = ctx->GDLonA = gd_adlon(ctx, ctx->GDLon - ctx->GDLonD);
                ctx->GDLonMax = ctx->GDLonB = gd_adlon(ctx, ctx->GDLon + ctx->GDLonD);

tda_out("GDLonD=%g GDLonMin=%g GDLonMax=%g\n",ctx->GDLonD,ctx->GDLonMin,ctx->GDLonMax);



                x = degree_to_arc(ctx, 90.0 - ctx->GDRHem);
                ctx->GDRHemP = cos(x);

                ctx->GDPX = ctx->GDPY = 0.0;
                ctx->GDPXA = ctx->GDPYA = -ctx->GDRHemP;
                ctx->GDPXB = ctx->GDPYB =  ctx->GDRHemP;
                ctx->GDLonD = ctx->GDLatD = ctx->GDRHem;

                ctx->PXLen = ctx->PMXLen;
                if (ctx->PMYLenFlg)
                    ctx->PYLen = ctx->PMYLen;
                else 
                    ctx->PYLen = ctx->PXLen;
            }

            break;




        case GDPROJ1:
            ctx->GDPROJ = GDPROJ1;
            printf1(ctx, "Parallel projection (rectangular region).\n");
            gd_view(ctx, 0);

            setup_3proj(ctx, ctx->GDLon,ctx->GDLat);   /* set up projection matrix */

            ctx->GDLonMin = ctx->GDLonA = gd_adlon(ctx, ctx->GDLon - ctx->GDLonD);
            ctx->GDLonMax = ctx->GDLonB = gd_adlon(ctx, ctx->GDLon + ctx->GDLonD);

            ctx->GDLatMin = dmax(ctx, ctx->GDLat - ctx->GDLatD,-90.0);
            ctx->GDLatMax = dmin(ctx, ctx->GDLat + ctx->GDLatD, 90.0);

            if (ctx->GDLat + ctx->GDLatD > 90.0 || ctx->GDLat - ctx->GDLatD < -90.0) {
                ctx->GDLonMin = -180.0;
                ctx->GDLonMax =  180.0;
            }
            ctx->GDLatA = gd_adlat(ctx, ctx->GDLat - ctx->GDLatD);
            ctx->GDLatB = gd_adlat(ctx, ctx->GDLat + ctx->GDLatD);

            ctx->GDPX = ctx->GDPY = 0.0;

            ctx->GDPXA = ctx->GDLon - ctx->GDLonD; y = ctx->GDLat;
            gd_proj_p(ctx, &ctx->GDPXA,&y);

            ctx->GDPXB = ctx->GDLon + ctx->GDLonD; y = ctx->GDLat;
            gd_proj_p(ctx, &ctx->GDPXB,&y);

            ctx->GDPYA = ctx->GDLat - ctx->GDLatD; x = ctx->GDLon;
            gd_proj_p(ctx, &x,&ctx->GDPYA);

            ctx->GDPYB = ctx->GDLat + ctx->GDLatD; x = ctx->GDLon;
            gd_proj_p(ctx, &x,&ctx->GDPYB);

            ctx->PXLen = ctx->PMXLen;
            if (ctx->PMYLenFlg)
                ctx->PYLen = ctx->PMYLen;
            else 
                ctx->PYLen = ctx->PXLen * (ctx->GDPYB - ctx->GDPYA) / (ctx->GDPXB - ctx->GDPXA);

            if (gd_check_inverse(ctx) == 0) {
                printf1(ctx, "Error: rectangle exceeds limits.\n");
                goto PSGFin;
            }
            break;                          

        case GDPROJ2:
            ctx->GDPROJ = GDPROJ2;
            ctx->GDLonD = ctx->GDLatD = 90.0;
            printf1(ctx, "Parallel projection (hemisphere).\n");
            gd_view(ctx, 1);
        
            setup_3proj(ctx, ctx->GDLon,ctx->GDLat);   /* set up projection matrix */

            ctx->GDLonA = gd_adlon(ctx, ctx->GDLon - ctx->GDLonD);
            ctx->GDLonB = gd_adlon(ctx, ctx->GDLon + ctx->GDLonD);

            ctx->GDLonMin = -180.0;
            ctx->GDLonMax =  180.0;

            ctx->GDLatMin = dmax(ctx, ctx->GDLat - ctx->GDLatD,-90.0);
            ctx->GDLatMax = dmin(ctx, ctx->GDLat + ctx->GDLatD, 90.0);

            ctx->GDLatA = gd_adlat(ctx, ctx->GDLat - ctx->GDLatD);
            ctx->GDLatB = gd_adlat(ctx, ctx->GDLat + ctx->GDLatD);

            ctx->GDPX = ctx->GDPY = 0.0;
            ctx->GDPXA = ctx->GDPYA = -1.0;
            ctx->GDPXB = ctx->GDPYB =  1.0;

            ctx->PXLen = ctx->PMXLen;
            if (ctx->PMYLenFlg)
                ctx->PYLen = ctx->PMYLen;
            else 
                ctx->PYLen = ctx->PXLen;
            break;


        default:                            
            printf1(ctx, "Error: unknown projection %d.\n",ctx->PMPROJ);
            goto PSGFin;
    }

    /* create 2d coordinate system */

    newline(ctx);
    ctx->PSLog[0] = ctx->PSLog[1] = 0;    /* no logarithmic axes */
    ctx->XOrg = ctx->PMXOrg;              /* origin of PostScript figure */
    ctx->YOrg = ctx->PMYOrg;
    ctx->SCALXFac = 1.0;             /* PostScript x scaling factor */
    ctx->SCALYFac = 1.0;             /* PostScript y scaling factor */
    ctx->ROTFac = ctx->PMPSRot;           /* PostScript rotation factor */
    ctx->PA1[0] = ctx->GDPXA;            
    ctx->PA2[0] = ctx->GDPXB;            
    ctx->PA1[1] = ctx->GDPYA;              
    ctx->PA2[1] = ctx->GDPYB;              
    setup_2dsys(ctx);

    ctx->GDPXPS = ps_2dx(ctx, ctx->GDPX);      /* center in PostScript coordinates */
    ctx->GDPYPS = ps_2dy(ctx, ctx->GDPY);
    ctx->GDPSRadius = ps_2dx(ctx, ctx->GDPX) - ps_2dx(ctx, ctx->GDPXA);

    ctx->PSPROJ = ctx->GDPROJ;        /* make valid projection globally known */
    err = 0;

PSGFin:
    p_clean(ctx);
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_geo(opt)  Return 0 if a valid geographical coordinate system      */
/*                  exists, otherwise return -1. If opt != 0 error message. */

int check_geo(TDAContext *ctx, int opt)
{
    if (ctx->PSPROJ)
        return(0);

    if (opt)
        printf1(ctx, "Error: no geographical coordinate system.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  gd_view(opt)   Print view and region.                                   */

void gd_view(TDAContext *ctx, int opt)
{
    printf1(ctx, "View: %g,%g",ctx->GDLon,ctx->GDLat);
    if (opt == 0)
        printf1(ctx, ". Region: %g,%g",ctx->GDLonD,ctx->GDLatD);
    newline(ctx);
    newline(ctx);
}

/* ------------------------------------------------------------------------ */
/*  gd_check_inverse(). Return 1 if rectangle (GDPROJ1) is completely in    */
/*                      projection of globe, otherwise return 0.            */

int gd_check_inverse(TDAContext *ctx)
{
    double x,y;

    x = ctx->GDPXA;
    y = ctx->GDPYA;
    if (gd_proj_inv(ctx, &x,&y) == 0)
        return(0);

    x = ctx->GDPXA;
    y = ctx->GDPYB;
    if (gd_proj_inv(ctx, &x,&y) == 0)
        return(0);

    x = ctx->GDPXB;
    y = ctx->GDPYA;
    if (gd_proj_inv(ctx, &x,&y) == 0)
        return(0);

    x = ctx->GDPXB;
    y = ctx->GDPYB;
    if (gd_proj_inv(ctx, &x,&y) == 0)
        return(0);

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  gd_adlon(lon)   Return adjusted longitude.                              */

double gd_adlon(TDAContext *ctx, double lon)
{
    (void)ctx;        /* unused: the signature is shared */
    if (lon > 180.0)
        return(lon - 360.0);
    else if (lon < -180.0)
        return(lon + 360.0);
    return(lon);
}

/* ------------------------------------------------------------------------ */
/*  gd_adlat(lat)   Return adjusted latitude.                               */

double gd_adlat(TDAContext *ctx, double lat)
{
    (void)ctx;        /* unused: the signature is shared */
    if (lat < -90.0)
        return(-180.0 - lat);
    else if (lat > 90.0)
        return(180.0 - lat);    
    return(lat);
}

/* ------------------------------------------------------------------------ */
/*  gd_proj_p(x,y)  Assume x and y are geographical coordinates. Return     */
/*                  projected coordinates based on parallel projection.     */
/*                  Return 0 if point is visible, otherwise return 1.       */

int gd_proj_p(TDAContext *ctx, double *x,double *y)
{
    double lon,lat,r,z;

    lon = *x;
    lat = *y;
    r = 1.0;
    ps_3dgeo(ctx, &lon,&lat,&r);
    ps_3dprj1(ctx, lon,lat,r,x,y,&z);
    if (z < 0.0)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gd_proj_inv(x,y)                                                        */
/*                                                                          */
/*  Given user coordinates x and y, return geographical coordinates based   */
/*  on parallel projection. Return 1 if the inverse projection exists,      */  
/*  otherwise return 0.                                                     */

int gd_proj_inv(TDAContext *ctx, double *x,double *y)
{
    double lon,lat,r,z;

    z = 1.0 - *x * *x - *y * *y;
    if (z <= 0.0) {
        if (z < -ctx->EPSI1)   
            return(0);
        z = 0.0;
    }
    else
        z = sqrt(z);

    ps_3dprj_inv(ctx, *x,*y,z,&lon,&lat,&r);
    ps_3dgeo_inv(ctx, &lon,&lat,&r);
    *x = lon;
    *y = lat;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_x_in_region(x)                                                    */
/*                                                                          */
/*  Given a longitude x, return 1 if x is in the selected region,           */ 
/*  otherwise return 0.                                                     */

int check_x_in_region(TDAContext *ctx, double x)
{
    if (ctx->GDLonMin < ctx->GDLonMax - ctx->EPSI1) {
        if (x >= ctx->GDLonMin && x <= ctx->GDLonMax)
            return(1);
        return(0);
    }
    if (x >= ctx->GDLonMin || x <= ctx->GDLonMax)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_y_in_region(y)                                                    */
/*                                                                          */
/*  Given a latitude y, return 1 if y is in the selected region,            */ 
/*  otherwise return 0.                                                     */

int check_y_in_region(TDAContext *ctx, double y)
{
    if (y >= ctx->GDLatMin && y <= ctx->GDLatMax)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_inside(lon,lat)   Return 1 if (lon,lat) is in selected region.    */

int check_inside(TDAContext *ctx, double lon,double lat)
{
    if (check_x_in_region(ctx, lon) && check_y_in_region(ctx, lat))
        return(1);
    return(0);
} 

/* ------------------------------------------------------------------------ */
/*  check_clip(x,y) Assume x and y are user coordinates. Return 0 if (x,y)  */
/*                  is in current plot region, otherwise return 1.          */ 

int check_clip(TDAContext *ctx, double x,double y)
{
    if (x >= ctx->GDPXA && x <= ctx->GDPXB && y >= ctx->GDPYA && y <= ctx->GDPYB)
        return(0);
    return(1);
} 

/* ------------------------------------------------------------------------ */
/*  gd_proj_gety(lat)                                                       */
/*                                                                          */ 
/*  Given a latitude lon, this function returns the projected coordinate    */
/*  depending on GDPROJ.                                                    */

double gd_proj_gety(TDAContext *ctx, double lat)
{
    double y;

    switch (ctx->GDPROJ) {

        case GDPROJ10:
            y = lat;         
            break;

        case GDPROJ11:
            y = sin(degree_to_arc(ctx, lat)) * 90.0;
            break;

        case GDPROJ12:
            if (lat > 89.0)
                lat = 89.0;
            else if (lat < -89.0)
                lat = -89.0;

            y = degree_to_arc(ctx, 45.0 + lat / 2.0);              
            y = log(tan(y));
            break;
/************************
        case GDPROJ20:

            r = gd_proj_p(ctx, x,y);
            if (r)
                return(-1);

            if (GDPROJA == 1 || GDPROJA == 2) {
                *x *= GDRHem;
                *y *= GDRHem;
            }
**************/








            break;

        default:
            printf1(ctx, "gd_proj_gety: unknown projection %d\n",ctx->GDPROJ);
            exit(0); 
    }
    return(y);
}

/* ------------------------------------------------------------------------ */
/*  gd_proj_getxy(x,y)                                                      */
/*                                                                          */ 
/*  Assume geographical coordinates x and y. Return the projection of x     */
/*  and y depending on GDPROJ.                                              */
/*                                                                          */
/*  Return   0  if visible and inside clip region                           */
/*           1  if not in clip region                                       */
/*          -1  if not visible                                              */

int gd_proj_getxy(TDAContext *ctx, double *x,double *y)
{
    int r;
    double lon,lat;

    lon = *x;
    lat = *y;

    switch (ctx->GDPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            *x = gd_adjust_lon(ctx, lon);
            *y = gd_proj_gety(ctx, lat);
            r = check_clip(ctx, *x,*y);
            return(r);

        case GDPROJ20:
            if (ctx->GDPROJA == 1) {
                if (90.0 - *y >= ctx->GDRHem)
                    return(-1);
                r = gd_proj_p(ctx, x,y);
                if (r)
                    return(-1);
                return(0);
            }
            else if (ctx->GDPROJA == 2) {
                if (90.0 + *y >= ctx->GDRHem)
                    return(-1);
                r = gd_proj_p(ctx, x,y);
                if (r)
                    return(-1);
                return(0);
            }
            else {
                tda_out("lon=%g lat=%g\n",lon,lat);


                if (check_inside(ctx, lon,lat) == 0) {
                    tda_out("not inside\n"); 
                    return(1);
                }
                r = gd_proj_p(ctx, x,y);
                if (r) {
                    tda_out("not visible\n");
                    return(-1);
                }

                /* check whether the projected point is inside clip region */

                tda_out("projected: %g %g\n",*x,*y);
   
                if (check_clip(ctx, *x,*y) != 0) {
                    tda_out("not in clip region\n");
                    return(1);
                }
                tda_out("OK\n");


                return(0);
            }

        default:
            printf1(ctx, "Error: unknown projection %d\n",ctx->GDPROJ);
            exit(0);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gd_adjust_lon(lon)                                                      */
/*                                                                          */
/*  Given longitude lon, return rotated value that is nearest to the        */
/*  current clip area.                                                      */
                      
double gd_adjust_lon(TDAContext *ctx, double lon)
{
    if (lon < ctx->GDPXA && lon + 360.0 <= ctx->GDPXB)
        return(lon + 360.0);
    else if (lon > ctx->GDPXB && lon - 360.0 >= ctx->GDPXA)
        return(lon - 360.0);
    return(lon);
}

/* -##--------------------------------------------------------------------- */
/*  sdpgrat   Meridians and parallels.                                      */
/*                                                                          */
/*  sdpgrat(                                                                */
/*      cont=...,   if 1 draw contour of selected region, def. 0            */
/*                  (depends on type of projection)                         */
/*      lt=...,     line typ for contour, def. 1                            */
/*      lw=...,     line width for contour, def. 0.2 (mm)                   */
/*      gs=...,     greyscale value for cont option, def. 1 (white)         */
/*      lon=...,    sequence of longitudes for grid lines                   */
/*      lat=...,    sequence of latitudes for grid lines                    */
/*      sc=...,     degrees around pole, def. 10                            */ 
/*      lt1=...,    line typ for grid lines, def. 1                         */
/*      lw1=...,    line width for grid lines, def. 0.05 (mm)               */
/*      fsx=...,    font size for labels in x direction, def. 0 (mm)        */
/*      fsy=...,    font size for lables in y direction, def. 0 (mm)        */
/*      fmt=...,    print format for labels, def. 0.0                       */  
/*  );                                                                      */
/*                                                                          */
/*  Supported projections:                                                  */
/*      GDPROJ10:   cylindrical (equidistant)                               */
/*      GDPROJ11:   cylindrical (equal-area)                                */
/*                                                                          */
/*  Labels are only drawn if fsx / fsy > 0. Also, this depends on the       */  
/*  type of projection.                                                     */
/*                                                                          */
/*  GDPROJ10: Labels are drawn for all x and y grid lines.                  */
/*  GDPROJ11: Labels are drawn for all x and y grid lines. Note: this       */
/*      can result in overlapping labels, then use separate sdpgrat         */
/*      commands, or add labels with the pltext command.                    */
/*  GDPROJ12: Labels are drawn for all x and y grid lines. Note: this       */
/*      can result in overlapping labels, then use separate sdpgrat         */
/*      commands, or add labels with the pltext command.                    */
/*                                                                          */  
/*  GDPROJ20:                                                               */
/*      GDPROJA = 1,2,3     Boundary area is a circle.                      */
/*      GDPROJA = 1,2       lon defines straight lines, lat defines         */
/*                          concentric circles. Labels are only drawn for   */
/*                          lon parameter.                                  */
/*      GDPROJA = 3         no labels are drawn.                            */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdpgrat(TDAContext *ctx)  
{
   	int	err;              

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Graticule for geographical coordinate system. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 7,3,0)) {     /* get parameters */
        goto SDPGRATFin;
    }   
    if (ctx->PSPROJ == 0) {
        printf1(ctx, "Error: need a geographical coordinate system.\n");
        return(-1);
    }
    if (ctx->PMCONT) {
        printf1(ctx, "Drawing bounds of selected region.\n");
        mapgrid_bounds(ctx);
    }
    if (ctx->PMNTP > 0) {
        printf1(ctx, "Drawing grid lines (meridians).\n");
        if (ctx->PMSCFlg == 0)
            ctx->PMSC = 10.0;
        else {
            if (ctx->PMSC < 0.0)
                ctx->PMSC = 0.0;
            else if (ctx->PMSC > 60.0)
                ctx->PMSC = 60.0;
        }
        mapgrid_lon(ctx, ctx->PMSC);
    }
    if (ctx->PMNTP1 > 0) {
        printf1(ctx, "Drawing grid lines (parallels).\n");
        mapgrid_lat(ctx);
    }
    err = 0;

SDPGRATFin:
    p_clean(ctx);
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_bounds    Draw bounds of selected region.                       */

void mapgrid_bounds(TDAContext *ctx)
{
    double r;

    fprintf(ctx->PSFd,"\n%%#%d: plgeogrid (bounds)\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);
    ps_lwidth(ctx, ctx->PMLW);

    switch (ctx->GDPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            mapgrid_bounds_r(ctx);     
            break;

        case GDPROJ20:
            if (ctx->GDPROJA == 0)  
                mapgrid_bounds_r(ctx);     
            else {
                r = degree_to_arc(ctx, 90.0 - ctx->GDRHem);
                r = cos(r);
                mapgrid_arc(ctx, ctx->GDPX,ctx->GDPY,0.0,360.0,r,1);
            }
            break;

        default: printf1(ctx, "Unknown projection %d.\n",ctx->GDPROJ);
            return;
    }
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n");
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_bounds_r    draw rectangular bounds.                            */

void mapgrid_bounds_r(TDAContext *ctx)
{
    ps_2dplot(ctx, ctx->GDPXA,ctx->GDPYA,0);
    ps_2dplot(ctx, ctx->GDPXA,ctx->GDPYB,1);
    ps_2dplot(ctx, ctx->GDPXB,ctx->GDPYB,1);
    ps_2dplot(ctx, ctx->GDPXB,ctx->GDPYA,1);
    ps_2dplot(ctx, ctx->GDPXA,ctx->GDPYA,1);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_arc(x,y,a,b,r,gflag)                                            */
/*                                                                          */
/*  Draw arc from a to b with radius r. Center is given by (x,y) in         */
/*  user coordinates for the projection.                                    */
/*  If rflag != 0 fill according to PMGSFlg.                                */

void mapgrid_arc(TDAContext *ctx, double x,double y,double a,double b,double r,int gflag)
{
    double rp;

    /* translate to PostScript coordinates */

    rp = (ps_2dx(ctx, x + r) - ps_2dx(ctx, x - r)) / 2.0;
    x = ps_2dx(ctx, x);  
    y = ps_2dy(ctx, y);

    if (a < 0.0)
        a = 0.0;
    if (b > 360.0)
        b = 360.0;
            
    fprintf(ctx->PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f arc\n",x,y,rp,a,b);

    if (gflag != 0 && ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        if (a > 0.0 || b < 360.0)  
            fprintf(ctx->PSFd,"%5.2f %5.2f l\nclosepath\n",x,y);
        fprintf(ctx->PSFd,"gsave\n");
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon(sc)   Draw grid lines: meridians. If PMFSX > 0 also plot    */
/*                    labels. Depends on type of projection. sc is degrees  */
/*                    around pole for azimuthal projections.                */

void mapgrid_lon(TDAContext *ctx, double sc)
{
    int i;
    double x,xa,ya,yb;

    if (ctx->PMLW <= 0.0)
        return;

    fprintf(ctx->PSFd,"\n%%#%d: plgeogrid (meridians)\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT1);
    ps_lwidth(ctx, ctx->PMLW1);

    if (ctx->GDPROJ == GDPROJ10 || ctx->GDPROJ == GDPROJ11 || ctx->GDPROJ == GDPROJ12) {

        for (i = 0; i < ctx->PMNTP; ++i) {
            x = ctx->PMTP[i];                          
            if (x < -180.0 || x > 180.0)
                continue;
     
            if (check_x_in_region(ctx, x)) {
                xa = x;
                if (xa < ctx->GDPXA) {
                    if (fabs(xa + 180.0) <= ctx->EPSI1)
                        continue;
                    xa += 360.0;
                }
                if (ctx->PMLW > 0.0)
                    mapgrid_line(ctx, xa,ctx->GDPYA,xa,ctx->GDPYB);        
                if (ctx->PMFSX > 0.0) {
                    mapgrid_label(ctx, 0,x,xa,ctx->GDPYA,ctx->PMFSX,2);
                    if (fabs(xa - ctx->GDPXA) <= ctx->EPSI1 && fabs(ctx->GDPXB - ctx->GDPXA - 360.0) <= ctx->EPSI1 
                        && fabs(xa + 180.0) > ctx->EPSI1)
                        mapgrid_label(ctx, 0,x,ctx->GDPXB,ctx->GDPYA,ctx->PMFSX,2);
                }
            }
        }
    }
    if (ctx->GDPROJ == GDPROJ20 && ctx->GDPROJA != 0) {
                  
        for (i = 0; i < ctx->PMNTP; ++i) {
            x = ctx->PMTP[i];                          
            if (x < -180.0 || x > 180.0)
                continue;

            if (ctx->GDPROJA == 1 || ctx->GDPROJA == 2) {
                if (ctx->PMLW > 0.0)
                    mapgrid_lon_azi(ctx, x);
                if (ctx->PMFSX > 0.0)
                    mapgrid_lon_azi_lab1(ctx, x,x,ctx->PMFSX);
            }
            else if (ctx->GDPROJA == 3) {
                if (ctx->PMLW > 0.0)
                    mapgrid_lon_azi1(ctx, x,sc);
            }
        }
    }
    if (ctx->GDPROJ == GDPROJ20 && ctx->GDPROJA == 0) {

        if (ctx->PMLW > 0.0) {
            set_clip(ctx);     

            for (i = 0; i < ctx->PMNTP; ++i) {
                x = ctx->PMTP[i];                          
                if (x < -180.0 || x > 180.0)
                    continue;

                if (check_x_in_region(ctx, x)) {
                    ya = dmax(ctx, ctx->GDLatMin,sc - 90.0);
                    yb = dmin(ctx, ctx->GDLatMax,90.0 - sc);
                    if (ya < yb)
                        map_geodesic_a(ctx, x,ya,x,yb);
                }
            }
        }
        if (ctx->PMFSX > 0.0) {
            fprintf(ctx->PSFd,"grestore\n");
            fprintf(ctx->PSFd,"gsave\n");

            for (i = 0; i < ctx->PMNTP; ++i) {
                x = ctx->PMTP[i];                          
                if (x < -180.0 || x > 180.0)
                    continue;

                if (check_x_in_region(ctx, x))  
                    mapgrid_lon_azi_lab0(ctx, x,x,ctx->PMFSX);          
            }
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
}
                  
/* -##--------------------------------------------------------------------- */
/*  mapgrid_line(xa,ya,xb,yb)                                               */
/*                                                                          */
/*  Draw 2d line from (xa,ya) to (xb,yb). Assume that values are given in   */
/*  user coordinates.                                                       */

void mapgrid_line(TDAContext *ctx, double xa,double ya,double xb,double yb)
{
    fprintf(ctx->PSFd,"gsave\n");
    ps_2dplot(ctx, xa,ya,0);
    ps_2dplot(ctx, xb,yb,1);
    fprintf(ctx->PSFd,"stroke\ngrestore\n");
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_label(opt,lab,x,y,fs,m)                                         */
/*                                                                          */
/*  Plot labels. If opt = 0 then for the X axis, otherwise for the Y axis.  */
/*  Assume that x and y are already given in user coordinates.              */

void mapgrid_label(TDAContext *ctx, int opt,double lab,double x,double y,double fs,int m)
{
    int l;
    double px,py,size;
    char buf[30];

    rt_snprintf_d(buf,sizeof(buf),ctx->PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    px = (x - ctx->GDPXA) / ctx->UXLen;
    py = (y - ctx->GDPYA) / ctx->UYLen;
    px *= ctx->PSXLen;
    py *= ctx->PSYLen;

    size = 1.5 * fs * ctx->PtMM;           
    l = (int)(strlen(buf));

    if (opt == 0) {
        py -= 1.5 * size;
        plot_str(ctx, px,py,buf,size,0,1,0,1,0);
    }
    else {
        px += m * size;                     
        py -= 0.3 * size;
        plot_str(ctx, px,py,buf,size,0,2,0,1,0);
    } 
    upd_bbox(ctx, 1,px + (double)l * size,py);     
    upd_bbox(ctx, 1,px - (double)l * size,py);     
    upd_bbox(ctx, 1,px,py + 1.3 * fs * ctx->PtMM);
    upd_bbox(ctx, 1,px,py - 0.5 * fs * ctx->PtMM);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lon_azi(lon)                                                    */
/*                                                                          */
/*  Draw meridian as a straight line, beginning at projection center.       */
/*  Called for GDPROJ20 and GDPROJA = 1 or 2.                               */

void mapgrid_lon_azi(TDAContext *ctx, double lon)
{
    double a,x,y;

    a = degree_to_arc(ctx, lon) - Pi / 2.0;
    x = cos(a) * ctx->GDRHemP;
    y = sin(a) * ctx->GDRHemP;

    if (ctx->GDPROJA == 2)  
        y = -y;
        
    ps_2dplot(ctx, x * 0.05,y * 0.05,0);
    ps_2dplot(ctx, x,y,1);
    fprintf(ctx->PSFd,"stroke\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon_azi_lab0(lab,lon,fs)                                        */
/*                                                                          */
/*  Plot label lab at longitude lon assuming GDPROJA = 0.                   */

void mapgrid_lon_azi_lab0(TDAContext *ctx, double lab,double lon,double fs)
{
    int l;
    double x,y,px,py,size;
    char buf[30];

    x = lon;
    if (ctx->GDLat >= 0.0)
        y = ctx->GDLatMin;
    else
        y = ctx->GDLatMax;

    if (gd_proj_p(ctx, &x,&y) != 0 || x < ctx->GDPXA || x > ctx->GDPXB)  
        return;

    if (fabs(lab + 180.0) < ctx->EPSI1)  /* don't use -180.0 */
        return;

    rt_snprintf_d(buf,sizeof(buf),ctx->PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    size = 1.5 * fs * ctx->PtMM;           
    px = ps_2dx(ctx, x);

    if (ctx->GDLat >= 0.0)  
        py = ps_2dy(ctx, ctx->GDPYA) - 1.5 * size;
    else
        py = ps_2dy(ctx, ctx->GDPYB) + size;
                         
    l = (int)(strlen(buf));
    plot_str(ctx, px,py,buf,size,0,1,0,1,0);

    upd_bbox(ctx, 1,px + (double)l * size,py);     
    upd_bbox(ctx, 1,px - (double)l * size,py);     
    upd_bbox(ctx, 1,px,py + 1.3 * fs * ctx->PtMM);
    upd_bbox(ctx, 1,px,py - 0.5 * fs * ctx->PtMM);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lon_azi_lab1(lab,lon,fs)                                        */
/*                                                                          */
/*  Plot label lab at longitude lon assuming GDPROJA = 1 or 2.              */

void mapgrid_lon_azi_lab1(TDAContext *ctx, double lab,double lon,double fs)
{
    int l,cflag;
    double a,x,y,px,py,size;
    char buf[30];

    if (fabs(lab + 180.0) < ctx->EPSI1)  /* don't use -180.0 */
        return;

    rt_snprintf_d(buf,sizeof(buf),ctx->PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    a = degree_to_arc(ctx, lon) - Pi / 2.0;
    x = cos(a) * ctx->GDRHemP * 1.05;
    y = sin(a) * ctx->GDRHemP * 1.05;

    if (ctx->GDPROJA == 2)  
        y = -y;
      
    px = ps_2dx(ctx, x);
    py = ps_2dy(ctx, y);

    size = fs * ctx->PtMM;           
    py -= size / 2.0;
    size *= 1.5;

    l = (int)(strlen(buf));

    if (lon < 0.0)
        lon += 360.0;

    if (lon <= 30.0)
        cflag = 1;
    else if (lon <= 150.0)
        cflag = 0;
    else if (lon <= 210.0)      
        cflag = 1;
    else if (lon <= 330.0)
        cflag = 2;
    else
        cflag = 1;

    plot_str(ctx, px,py,buf,size,0,cflag,0,1,0);

    upd_bbox(ctx, 1,px + (double)l * size,py);     
    upd_bbox(ctx, 1,px - (double)l * size,py);     
    upd_bbox(ctx, 1,px,py + 1.3 * fs * ctx->PtMM);
    upd_bbox(ctx, 1,px,py - 0.5 * fs * ctx->PtMM);
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon_azi1(lat,sc)                                                */
/*                                                                          */
/*  Draw meridians. Called for GDPROJ20, GDPROJA = 3. sc is degrees around  */
/*  pole.                                                                   */

void mapgrid_lon_azi1(TDAContext *ctx, double lon,double sc)
{
    double a,b,dlon;

    dlon = gd_adlon(ctx, lon + 180.0);
        
    if (ctx->GDLat >= 0.0) {
        b = 90.0 - sc;
        a = dmin(ctx, b,ctx->GDLat);
        if (a < b) {
            map_geodesic_a(ctx, lon,a,lon,b);
            map_geodesic_a(ctx, dlon,a,dlon,b);
        }
        b = sc - 90.0;
        if (a > b)   
            map_geodesic_a(ctx, lon,a,lon,b);
    }
    else {
        b = sc - 90.0;
        a = dmax(ctx, b,ctx->GDLat);
        if (a > b) {
            map_geodesic_a(ctx, lon,a,lon,b);
            map_geodesic_a(ctx, dlon,a,dlon,b);
        }
        b = 90.0 - sc;
        if (a < b)   
            map_geodesic_a(ctx, lon,a,lon,b);
    }
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lat    Draw grid lines: parallels. If PMFSY > 0 also plot       */
/*                 labels. Depends on type of projection.                   */

void mapgrid_lat(TDAContext *ctx)
{
    int i;
    double y,yp;

    if (ctx->PMLW <= 0.0)
        return;

    fprintf(ctx->PSFd,"\n%%#%d: plgeogrid (parallels)\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT1);
    ps_lwidth(ctx, ctx->PMLW1);

    if (ctx->GDPROJ == GDPROJ20 && ctx->GDPROJA == 0) 
        set_clip(ctx);     
          
    for (i = 0; i < ctx->PMNTP1; ++i) {

        y = ctx->PMTP1[i];                          
        if (y < -90.0 || y > 90.0)
            continue;
     
        switch (ctx->GDPROJ) {

            case GDPROJ10:
            case GDPROJ11:
            case GDPROJ12:
                if (check_y_in_region(ctx, y)) {
                    yp = gd_proj_gety(ctx, y);
                    if (ctx->PMLW > 0.0)
                        mapgrid_line(ctx, ctx->GDPXA,yp,ctx->GDPXB,yp);
                    if (ctx->PMFSY > 0.0)  
                        mapgrid_label(ctx, 1,y,ctx->GDPXB,yp,ctx->PMFSY,2);
                }
                break;  

            case GDPROJ20:
                if (ctx->GDPROJA == 1) {
                    if (y != 90.0 && 90.0 - y <= ctx->GDRHem) {
                        y = degree_to_arc(ctx, y);
                        y = cos(y);
                        mapgrid_arc(ctx, ctx->GDPX,ctx->GDPY,0.0,360.0,y,0);
                    }
                }
                else if (ctx->GDPROJA == 2) {
                    if (y != -90.0 && 90.0 + y <= ctx->GDRHem) {        
                        y = degree_to_arc(ctx, y);
                        y = cos(y);
                        mapgrid_arc(ctx, ctx->GDPX,ctx->GDPY,0.0,360.0,y,0);
                    }
                }
                else if (ctx->GDPROJA == 3) {    
                    mapgrid_lat_a3(ctx, y);                     
                }
                else if (ctx->GDPROJA == 0) {            /* rectangular region */
                    mapgrid_lat_a0(ctx, y);          
                }
                break; 

            default:
                printf1(ctx, "Unknown projection %d\n",ctx->GDPROJ);
                break;  
        }
    }
    fprintf(ctx->PSFd,"grestore\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lat_a0(lat)                                                     */
/*                                                                          */
/*  Draw parallel. Called for GDPROJ20, GDPROJA = 0 (rectangular region).   */

void mapgrid_lat_a0(TDAContext *ctx, double lat)
{
    double x,y,xc,yc,r,beta,sinb,cosb,tlat,r1,r2,s,t,phi;
                                
    if (ctx->GDLat >= 90.0 - ctx->EPSI1 || ctx->GDLat <= ctx->EPSI1 - 90.0)
        return;

    if (lat <= ctx->GDLatMin || lat >= ctx->GDLatMax)
        return;
    
    beta = degree_to_arc(ctx, lat);
    sinb = sin(beta);
    cosb = cos(beta);

    if (fabs(ctx->GDLat) <= ctx->EPSI1) {     /* straight lines */ 
        fprintf(ctx->PSFd,"gsave\n");
        x = ctx->GDPXA;
        y = ctx->GDPY + sinb;
        ps_2dplot(ctx, x,y,0);
        x = ctx->GDPXB;
        ps_2dplot(ctx, x,y,1);
        fprintf(ctx->PSFd,"stroke\n");
        fprintf(ctx->PSFd,"grestore\n");
        return;
    }
    ps_3dprj(ctx, 0.0,0.0,sinb,&xc,&yc);              /* xc should be zero */

    s = cos(degree_to_arc(ctx, 90.0 - fabs(ctx->GDLat)));
    tlat = tan(degree_to_arc(ctx, ctx->GDLat));
    t = -tlat * sinb / cosb;

    if (t <= -1.0 + ctx->EPSI1 || t >= 1.0 - ctx->EPSI1) {
        r1 = 0.0;
        r2 = 360.0;
    }
    else {  
        phi = acos(t) * 180.0 / Pi;
        if (ctx->GDLat >= 0.0) {
            r1 = 270.0 - phi;
            r2 = phi - 90.0;
        }
        else {
            r1 = 90.0 - phi;
            r2 = 90.0 + phi;
        }
    }
    xc = ps_2dx(ctx, xc);
    yc = ps_2dy(ctx, yc);
    r = cosb * (ps_2dx(ctx, 1.0) - ps_2dx(ctx, -1.0)) / 2.0;

    fprintf(ctx->PSFd,"gsave\n");
    fprintf(ctx->PSFd,"%7.3f %7.3f translate\n",xc,yc);
    fprintf(ctx->PSFd,"1.0 %7.4f scale\n",s);   
    fprintf(ctx->PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",r,r1,r2);
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lat_a3(lat)                                                     */
/*                                                                          */
/*  Draw parallel. Called for GDPROJ20, GDPROJA = 3.                        */

void mapgrid_lat_a3(TDAContext *ctx, double lat)
{
    double x,y,xc,yc,r,beta,sinb,cosb,tlat,r1,r2,s,t,phi;
                                
    if (ctx->GDLat >= 90.0 - ctx->EPSI1 || ctx->GDLat <= ctx->EPSI1 - 90.0)
        return;

    if (lat <= ctx->GDLat - 90.0 + ctx->EPSI1 ||
        lat >= ctx->GDLat + 90.0 - ctx->EPSI1) 
        return;
    
    beta = degree_to_arc(ctx, lat);
    sinb = sin(beta);
    cosb = cos(beta);

    if (fabs(ctx->GDLat) <= ctx->EPSI1) {     /* straight lines */ 
        fprintf(ctx->PSFd,"gsave\n");
        x = ctx->GDPX - cosb * (ctx->GDPX - ctx->GDPXA);
        y = ctx->GDPY + sinb * (ctx->GDPY - ctx->GDPYA);
        ps_2dplot(ctx, x,y,0);
        x = ctx->GDPX + cosb * (ctx->GDPX - ctx->GDPXA);
        ps_2dplot(ctx, x,y,1);
        fprintf(ctx->PSFd,"stroke\n");
        fprintf(ctx->PSFd,"grestore\n");
        return;
    }
    ps_3dprj(ctx, 0.0,0.0,sinb,&xc,&yc);              /* xc should be zero */
    s = cos(degree_to_arc(ctx, 90.0 - fabs(ctx->GDLat)));
    tlat = tan(degree_to_arc(ctx, ctx->GDLat));
    t = -tlat * sinb / cosb;

    if (t <= -1.0 + ctx->EPSI1 || t >= 1.0 - ctx->EPSI1) {
        r1 = 0.0;
        r2 = 360.0;
    }
    else {  
        phi = acos(t) * 180.0 / Pi;
        if (ctx->GDLat >= 0.0) {
            r1 = 270.0 - phi;
            r2 = phi - 90.0;
        }
        else {
            r1 = 90.0 - phi;
            r2 = 90.0 + phi;
        }
    }
    xc = ps_2dx(ctx, xc);
    yc = ps_2dy(ctx, yc);
    r = cosb * (ps_2dx(ctx, ctx->GDPXB) - ps_2dx(ctx, ctx->GDPXA)) / 2.0;

    fprintf(ctx->PSFd,"gsave\n");
    fprintf(ctx->PSFd,"%7.3f %7.3f translate\n",xc,yc);
    fprintf(ctx->PSFd,"1.0 %7.4f scale\n",s);   
    fprintf(ctx->PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",r,r1,r2);
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
}

/* -###-------------------------------------------------------------------- */
/*  sdpgeo  Plot points/lines in geographical coordinate system.            */
/*                                                                          */
/*  sdpgeo(                                                                 */
/*      ct=...,     type of connection, def. 0                              */
/*                  0 = connect by straight line in projection plane        */  
/*                  1 = connect by geodesic line                            */
/*      lt=...,     line typ, def. 1 (solid line)                           */
/*      lw=...,     line width, def. 0.05 (mm)                              */
/*      s=...,      number of marker symbol, def. 0                         */
/*      fs=...,     size of marker symbol, def. 2 mm                        */
/*                                                                          */
/*  ) = x1,y1,x2,y2,...;                                                    */
/*                                                                          */
/*                                                                          */  
/*  Return 0 if OK, -1 if error.                                            */

int sdpgeo(TDAContext *ctx)  
{
    register int i,j;
   	int	err,n;
    double x,y;

    err = -1;
    if (check_cmd(ctx, 2))
        return(-1);
                    
    printf1(ctx, "Plot in geographical coordinate system. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,7,1)) {     /* get parameters */
        goto SDPGEOFin;
    }   
    if (ctx->PSPROJ == 0) {
        printf1(ctx, "Need a geographical coordinate system.\n");
        return(-1);
    }
    if (ctx->PMLWFlg == 0)
        ctx->PMLW = 0.05;

    if (ctx->PMCT != 1)
        ctx->PMCT = 0;

    if (ctx->PMS < 1 || ctx->PMS > 17)        /* check for valid marker symbol */
        ctx->PMS = 0;           

    n = ctx->PMRHSN / 2;
    if (2 * n != ctx->PMRHSN) {
        printf1(ctx, "Error: right-hand side must define pairs of coordinates.\n");
        goto SDPGEOFin;
    }
    if (n == 1 && ctx->PMS == 0) {
        printf1(ctx, "Single points require marker symbol.\n");
        err = 0;
        goto SDPGEOFin;
    }
    if (alloc_acx(ctx, n + 1))
        return(-1);              
    if (alloc_acy(ctx, n + 1))
        return(-1);                  

    j = 0;
    for (i = 0; i < n; ++i) {
        x = ctx->PMRHSX[j++];
        y = ctx->PMRHSX[j++];
        if (x < -180.0 || x > 180.0 || y < -90.0 || y > 90.0) {
            printf1(ctx, "Error: right-hand side violates range of geographical coordinates.\n");
            goto SDPGEOFin;
        }
        ctx->AcX[i] = x;
        ctx->AcY[i] = y;
    }
    fprintf(ctx->PSFd,"\n%%#%d: plgeod\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMS != 0)
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);

    /* set_clip(); not used in this function */

    if (ctx->PMLW > 0.0)
        map_line(ctx, n,ctx->AcX,ctx->AcY,ctx->PMCT);

    if (ctx->PMS != 0) {
        for (i = 0; i < n; ++i)  
            map_point(ctx, ctx->AcX[i],ctx->AcY[i],0);
    }
    fprintf(ctx->PSFd,"grestore\n");
    printf1(ctx, "PostScript output written to: %s\n",ctx->PSFName);
    err = 0;

SDPGEOFin:
    p_clean(ctx);
    return(err);
}

/* -###-------------------------------------------------------------------- */
/*  sdpmap   Plot spatial objects in geographical coordinate system.        */
/*                                                                          */
/*  sdpmap(                                                                 */
/*      lt=...,     line typ for grid lines, def. 1                         */
/*      lw=...,     line width for grid lines, def. 0.05 (mm)               */
/*      s=...,      number of marker symbol, def. 0                         */
/*      fs=...,     size of marker symbol, def. 2 mm                        */
/*  );                                                                      */
/*                                                                          */
/*  This command requires the definition of a geographical coordinate       */
/*  system and a projection with the psetupg command, and a spatial data    */
/*  structure defined with the sdnvar command.                              */
/*                                                                          */
/*  The command plots all currently selected spatial objects based on the   */
/*  specified coordinate system and projection.                             */
/*                                                                          */
/*  It is assumed that all polygons, viewed on the surface of a sphere,     */
/*  are simple.                                                             */
/*                                                                          */
/*  Return 0 if OK, -1 if error.                                            */

int sdpmap(TDAContext *ctx)  
{
    register int i,j;
   	int	err,n1,n2,n3,n1p,n2p,n3p,n,cflag,typ,r,sdid;

    err = -1;
    if (check_cmd(ctx, 0))
        return(-1);
                    
    printf1(ctx, "Plot in geographical coordinate system. Current memory: %d bytes.\n",ctx->MemReq);

    if (parm(ctx, ctx->CmdBuf + 6,3,0)) {     /* get parameters */
        goto SDPMAPFin;
    }   
    if (check_sd(ctx, 0)) {
        printf1(ctx, "No spatial data defined.\n");
        return(-1);
    }
    if (ctx->PSPROJ == 0) {
        printf1(ctx, "Need a geographical coordinate system.\n");
        return(-1);
    }
    if (ctx->PMLWFlg == 0)
        ctx->PMLW = 0.05;

    if (ctx->PMS < 1 || ctx->PMS > 17)        /* check for valid marker symbol */
        ctx->PMS = 0;           

    fprintf(ctx->PSFd,"\n%%#%d: plgeosd\n",++ctx->PSONUM);
    fprintf(ctx->PSFd,"gsave\n");
    ps_ltyp(ctx, ctx->PMLT);  
    ps_lwidth(ctx, ctx->PMLW);
    if (ctx->PMS != 0)
        fprintf(ctx->PSFd,"/ssiz %4.2f def\n",ctx->PtMM * ctx->PMFS / 2.0);
    /* set_clip(); */

    n1p = n1 = 0;     /* number of points */
    n2p = n2 = 0;     /* number of lines */
    n3p = n3 = 0;     /* number of polygons */

    for (i = 0; i < ctx->NOC; ++i) {

        if ((typ  = (int)get_data(ctx, ctx->SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        sdid = (int)get_data(ctx, ctx->SDVarSDID,i);                      

        if (typ == 3)  
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(ctx, i,0,cflag,1)) < 1)
            goto SDPMAPFin;
             
        if (typ == 1) {
            n1++;
            if (ctx->PMS) {
                if (map_point(ctx, ctx->SDVarX[0],ctx->SDVarY[0],0))
                    n1p++;
            }
        }
        else if (typ == 2) {
            n2++;
            r = map_line(ctx, n,ctx->SDVarX,ctx->SDVarY,0);
            if (r < 0)
                goto SDPMAPFin;
            n2p += r;

            if (ctx->PMS) {
                for (j = 0; j < n; ++j)  
                    map_point(ctx, ctx->SDVarX[j],ctx->SDVarY[j],0);
            }
        }
        else if (typ == 3) {
            n3++;
            r = map_polygon(ctx, n - 1,ctx->SDVarX,ctx->SDVarY,sdid);   
            if (r < 0)  
                r = 0;
            else
                n3p += r;

            if (ctx->PMS) {
                if ((n = sd_getdata(ctx, i,0,1,1)) < 1)
                    goto SDPMAPFin;
                for (j = 1; j < n; ++j)  
                    map_point(ctx, ctx->SDVarX[j],ctx->SDVarY[j],0);
            }
        }
    }
    fprintf(ctx->PSFd,"grestore\n");

    printf1(ctx, "\nNumber of type 1 objects: %d",n1);
    if (n1p > 0)
        printf1(ctx, ". Plotted: %d",n1p);
    printf1(ctx, "\nNumber of type 2 objects: %d",n2);
    if (n2p > 0)
        printf1(ctx, ". Plotted: %d (parts)",n2p);
    printf1(ctx, "\nNumber of type 3 objects: %d",n3);
    if (n3p > 0)
        printf1(ctx, ". Plotted: %d (parts)",n3p);
    newline(ctx);
    printf1(ctx, "PostScript output written to: %s\n",ctx->PSFName);
    err = 0;

SDPMAPFin:
    p_clean(ctx);
    return(err);
}


/* ------------------------------------------------------------------------ */
/*  map_point(x,y,nc)    Plot point with geographical coordinates.          */
/*                                                                          */
/*  x and y are the geographical coordinates of a point. Plot the point     */
/*  with the current marker symbol PMS and size PMFS. If nc = 0 plot only   */
/*  if the point is in current clip region.                                 */
/*                                                                          */
/*  Return 1 if the point is plotted, otherwise return 0.                   */
    
int map_point(TDAContext *ctx, double x,double y,int nc)
{
    int r = 0;

    if (ctx->PMS == 0)
        return(0);

    switch (ctx->PSPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            r = map_point_cyl(ctx, x,y);  
            break;

        case GDPROJ20:
            r = map_point_azi(ctx, x,y,nc);
            break;

        default:
            tda_out("projection not supported.\n");
            break;  
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  map_point_cyl(x,y)  Called by map_point for cylindrical projections.    */

int map_point_cyl(TDAContext *ctx, double x,double y)
{
    double d;

    if (gd_proj_getxy(ctx, &x,&y) == 0) {    
        x = ps_2dx(ctx, x);
        y = ps_2dy(ctx, y);
        d = ctx->PtMM * ctx->PMFS / 2.0;
        ps_sym(ctx, ctx->PMS,x,y,d);
        return(1);
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  map_point_azi(x,y,nc)  Called by map_point for azimuthal projections.   */

int map_point_azi(TDAContext *ctx, double x,double y,int nc)
{
    double d;
            
tda_out("map_point_azi lon=%g %g\n",x,y);

    if (nc == 0) {
        if (gd_proj_getxy(ctx, &x,&y) == 0) {        
            x = ps_2dx(ctx, x);
            y = ps_2dy(ctx, y);
            d = ctx->PtMM * ctx->PMFS / 2.0;
            ps_sym(ctx, ctx->PMS,x,y,d);
            return(1);
        }
    }
    else {
        if (gd_proj_p(ctx, &x,&y)) {
            tda_out(" point  not visible\n");
            return(0);
        }
        upd_bbox(ctx, 0,x,y);
        x = ps_2dx(ctx, x);
        y = ps_2dy(ctx, y);
        d = ctx->PtMM * ctx->PMFS / 2.0;
        ps_sym(ctx, ctx->PMS,x,y,d);
        return(1);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  map_line(n,x,y,gflag)   Plot line with geographical coordinates.        */
/*                                                                          */
/*  x[i], y[i] (i=0,n-1) contains the points of a line. This function       */
/*  draws the line assuming geographical coordinates. If gflag = 0 connect  */
/*  points by straight lines, otherwise by geodesics.                       */
/*                                                                          */
/*  Return the number of (parts of) lines plotted.                          */
    
int map_line(TDAContext *ctx, int n,double *x,double *y,int gflag)
{
    int r = 0;

    switch (ctx->PSPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            if (gflag)
                r = map_line_cyl_g(ctx, n,x,y);
            else
                r = map_line_cyl_s(ctx, n,x,y);  
            break;
   
        case GDPROJ20:
            r = map_line_azi(ctx, n,x,y,gflag);  
            break;

        default:
            tda_out("projection not supported.\n");
            break;  
    }
    return(r);
}
           
/* -###-------------------------------------------------------------------- */
/*  map_line_cyl_s(n,x,y)                                                   */
/*                                                                          */
/*  Called by map_line to plot a line for cylindrical projections.          */
/*  Points are connected by straight lines.                                 */
/*  Return number of lines plotted.                                         */

int map_line_cyl_s(TDAContext *ctx, int n,double *x,double *y)
{
    register int i; 
    int rt,r1,r2,ca,cb,ra,rb,dir;
    double xa,ya,xb,yb,sax1,say1,sbx1,sby1,sax2,say2,sbx2,sby2;
    double yp;

    rt = 0;
    for (i = 1; i < n; ++i) {
        xa = x[i - 1];
        ya = y[i - 1];
        xb = x[i];
        yb = y[i];
        xa = gd_adjust_lon(ctx, xa);
        xb = gd_adjust_lon(ctx, xb);
    
        map_find_n(ctx, xa,ya,&xb,&yb);
        r1 = clip_line(ctx, xa,ya,xb,yb,ctx->GDPXA,ctx->GDLatA,ctx->GDPXB,ctx->GDLatB,
             &sax1,&say1,&sbx1,&sby1,&ca,&cb,&ra,&rb,&dir);

        xa = x[i - 1];
        ya = y[i - 1];
        xb = x[i];
        yb = y[i];
        xa = gd_adjust_lon(ctx, xa);
        xb = gd_adjust_lon(ctx, xb);

        map_find_n(ctx, xb,yb,&xa,&ya);
        r2 = clip_line(ctx, xa,ya,xb,yb,ctx->GDPXA,ctx->GDLatA,ctx->GDPXB,ctx->GDLatB,
             &sax2,&say2,&sbx2,&sby2,&ca,&cb,&ra,&rb,&dir);
    
        if (r1) {
            yp = gd_proj_gety(ctx, say1);
            ps_2dplot(ctx, sax1,yp,0);
            yp = gd_proj_gety(ctx, sby1);
            ps_2dplot(ctx, sbx1,yp,1);
            fprintf(ctx->PSFd,"stroke\n");    
            rt++;    

            if (r2) {
                if (fabs(sax1 - sax2) > ctx->EPSI1 || fabs(say1 - say2) > ctx->EPSI1 ||
                    fabs(sbx1 - sbx2) > ctx->EPSI1 || fabs(sby1 - sby2) > ctx->EPSI1) {  
                    yp = gd_proj_gety(ctx, say2);
                    ps_2dplot(ctx, sax2,yp,0);
                    yp = gd_proj_gety(ctx, sby2);
                    ps_2dplot(ctx, sbx2,yp,1);
                    fprintf(ctx->PSFd,"stroke\n");    
                    rt++;
                }
            }
        }
        else if (r2) {
            yp = gd_proj_gety(ctx, say2);
            ps_2dplot(ctx, sax2,yp,0);
            yp = gd_proj_gety(ctx, sby2);
            ps_2dplot(ctx, sbx2,yp,1);
            fprintf(ctx->PSFd,"stroke\n");    
            rt++;     
        }
    }
    return(rt);
}

/* ------------------------------------------------------------------------ */
/*  map_find_n(x,y,xa,ya)       Find nearest point.                         */

void map_find_n(TDAContext *ctx, double x,double y,double *xa,double *ya)
{
    (void)ctx; (void)y; (void)ya;        /* unused: the signature is shared */
    double d,xx;

    d = fabs(x - *xa);
    if (*xa <= x) {
        xx = *xa + 360.0;
        if (fabs(x - xx) < d)
            *xa = xx;
    }
    else {
        xx = *xa - 360.0;
        if (fabs(x - xx) < d)
            *xa = xx;
    }
}

/* -####------------------------------------------------------------------- */
/*  map_line_cyl_g(n,x,y)                                                   */
/*                                                                          */
/*  Called by map_line to plot a line for cylindrical projections.          */
/*  Points are connected by geodesics.                                      */
/*  Return number of lines plotted.                                         */

int map_line_cyl_g(TDAContext *ctx, int n,double *x,double *y)
{
    register int i; 
    int rt;
    double xa,ya,xb,yb;

    rt = 0;
    for (i = 1; i < n; ++i) {
        xa = x[i - 1];
        ya = y[i - 1];
        xb = x[i];
        yb = y[i];

        if (xa < ctx->GDPXA)
            xa += 360.0;
        else if (xa > ctx->GDPXB)
            xa -= 360.0;

        map_find_n(ctx, xa,ya,&xb,&yb); 
        rt = map_line_cyl_geod(ctx, xa,ya,xb,yb);
    }
    return(rt);
}

/* -####------------------------------------------------------------------- */
/*  map_line_cyl_geod(lon1,lat1,lon2,lat2)                                  */
/*                                                                          */
/*  Draw a geodesic line from (lon1,lat1) to (lon2,lat2). Assume            */
/*  cylindrical projection. Return 1 if at least some part of the geodesic  */
/*  is drawn, otherwise return 0.                                           */

#define NPNT_GEOD 100

int map_line_cyl_geod(TDAContext *ctx, double lon1,double lat1,double lon2,double lat2)
{
    int i,j,first,rt,ca,cb,ra,rb,dir;
    double d,sind,slon1,slon2,slat1,slat2,clon1,clon2,clat1,clat2;
    double f,x,y,z,a,b,t,rm,x0,xl,yl,delta,sax,say,sbx,sby;
    double xp[NPNT_GEOD + 1],yp[NPNT_GEOD + 1];

    /* check for antipodal points */

    if (fabs(lat1 + lat2) <= ctx->EPSI1 && fabs(fabs(lon1 - lon2) - 180.0) <= ctx->EPSI1)
        return(0);

    if (lon1 > lon2) {
        t = lon1;
        lon1 = lon2;
        lon2 = t;
        t = lat1;
        lat1 = lat2;
        lat2 = t;
    }
    rt = 0;
    d = geod_distance(ctx, lon1,lat1,lon2,lat2);
    if (d <= ctx->EPSI1)
        return(0);

    if (fabs(d - Pi) <= ctx->EPSI1)  
        return(0);            

    sind = sin(d);

    lon1 = degree_to_arc(ctx, lon1);
    lat1 = degree_to_arc(ctx, lat1);
    lon2 = degree_to_arc(ctx, lon2);
    lat2 = degree_to_arc(ctx, lat2);

    slon1 = sin(lon1);
    slon2 = sin(lon2);
    slat1 = sin(lat1);
    slat2 = sin(lat2);
    clon1 = cos(lon1);
    clat1 = cos(lat1);
    clon2 = cos(lon2);
    clat2 = cos(lat2);

    rm = 180.0 / Pi;
    delta = 1.0 / (double)NPNT_GEOD;

    for (i = 0; i <= NPNT_GEOD; ++i) {
        f = (double)i * delta;
        a = sin((1.0 - f) * d) / sind;
        b = sin(f * d) / sind;
        x = a * clat1 * clon1 + b * clat2 * clon2;
        y = a * clat1 * slon1 + b * clat2 * slon2;
        z = a * slat1 + b * slat2;
        t = sqrt(x * x + y * y);

        if (x == 0.0 && y == 0.0)
            return(0);
        xp[i] = atan2(y,x) * rm;

        if (z == 0.0 && t == 0.0)
            return(0);
        yp[i] = atan2(z,t) * rm;
    }
    xl = yl = x0 = 0.0;
    first = 0;
    j = -1;
    for (i = 0; i <= NPNT_GEOD; ++i) {
        x = xp[i];
        y = yp[i];
        if (gd_proj_getxy(ctx, &x,&y) != 0) {
            if (first) {
                xl = xp[i - 1];
                yl = yp[i - 1];
                gd_proj_getxy(ctx, &xl,&yl);             

                if (clip_line(ctx, xl,yl,x,y,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,
                    &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) && 
                    fabs(xl - sax) <= ctx->EPSI1 && fabs(yl - say) <= ctx->EPSI1) {
                    ps_2dplot(ctx, sbx,sby,first);
                }
/*              map_line_clip(xl,yl,x,y);           */
                fprintf(ctx->PSFd,"stroke\n");    
                rt = 1;
                first = 0;
            }
            continue;
        }
        if (j < 0)  
            x0 = x;
        else if (x < x0 - ctx->EPSI1)
            break;

        j = i;

        if (first == 0 && i > 0) {
            xl = xp[i - 1];
            yl = yp[i - 1];
            gd_proj_getxy(ctx, &xl,&yl);             

            if (clip_line(ctx, xl,yl,x,y,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,
                 &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) &&  
                 fabs(sbx - x) <= ctx->EPSI1 && fabs(sby - y) <= ctx->EPSI1) {
                 ps_2dplot(ctx, sax,say,first);
                 first = 1;
            }
        }
        ps_2dplot(ctx, x,y,first);
        first = 1;
    }
    if (j < 0)
        return(0);

    if (first) {
        fprintf(ctx->PSFd,"stroke\n");    
        first = 0;
        rt = 1;
    }
    if (j != NPNT_GEOD) {
        first = 0;
        for (i = NPNT_GEOD; i >= 0; --i) {
            x = xp[i];
            y = yp[i];
            if (gd_proj_getxy(ctx, &x,&y) != 0) {
                if (first) {

                    xl = xp[i + 1];
                    yl = yp[i + 1];
                    gd_proj_getxy(ctx, &xl,&yl);             

                    if (clip_line(ctx, xl,yl,x,y,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,
                        &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) && 
                        fabs(xl - sax) <= ctx->EPSI1 && fabs(yl - say) <= ctx->EPSI1) {
                        ps_2dplot(ctx, sbx,sby,first);
                    }
/*                  map_line_clip(xl,yl,x,y);       */
                    fprintf(ctx->PSFd,"stroke\n");    
                    rt = 1;
                    first = 0;
                }
                continue;
            }
            if (i <= j || x >= x0 - ctx->EPSI1)
                break;

            if (first == 0 && i < NPNT_GEOD) {
                xl = xp[i + 1];
                yl = yp[i + 1];
                gd_proj_getxy(ctx, &xl,&yl);             

                if (clip_line(ctx, xl,yl,x,y,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,
                     &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) &&  
                     fabs(sbx - x) <= ctx->EPSI1 && fabs(sby - y) <= ctx->EPSI1) {
                     ps_2dplot(ctx, sax,say,first);
                     first = 1;
                }
            }
            ps_2dplot(ctx, x,y,first);
            first = 1;
        }
    }
    if (first) {
        fprintf(ctx->PSFd,"stroke\n");    
        rt = 1;
    }
    return(rt);
}

/* -###-------------------------------------------------------------------- */
/*  map_line_clip(xa,ya,xb,yb)                                              */
/*                                                                          */
/*  (xa,ya) and (xb,yb) are points in the projection plane. also (xa,ya)    */
/*  is contained in the clip region. Connect by a straight line.            */

void map_line_clip(TDAContext *ctx, double xa,double ya,double xb,double yb)
{                 
    int ca,cb,ra,rb,dir;
    double sax,say,sbx,sby;

tda_out("maplineclip %g %g %g %g\n",xa,ya,xb,yb);


    if (clip_line(ctx, xa,ya,xb,yb,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,
             &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir)) {
tda_out("sax %g %g %g %g\n",sax,say,sbx,sby);
        if (fabs(xa - sax) <= ctx->EPSI1 && fabs(ya - say) <= ctx->EPSI1) {
tda_out("plotted sbx=%g %g\n",sbx,sby);
            ps_2dplot(ctx, sbx,sby,1);
        }
    }        
}

/* -###-------------------------------------------------------------------- */
/*  map_line_azi(n,x,y,gflag)                                               */
/*                                                                          */
/*  Called by map_line to plot a line for azimuthal projections. If         */
/*  gflag = 0 connect points by straight lines, otherwise by geodesics.     */
/*                                                                          */  
/*  Return number of lines plotted.                                         */

int map_line_azi(TDAContext *ctx, int n,double *x,double *y,int gflag)
{
    register int i;
    int r,nl = 0;

    tda_out("map_line_azi gflag=%d  \n",gflag);

    for (i = 1; i < n; ++i) {
        if (gflag == 0)
            r = map_line_straight(ctx, x[i - 1],y[i - 1],x[i],y[i]);
        else  
            r = map_geodesic_a(ctx, x[i - 1],y[i - 1],x[i],y[i]);

        if (r == 1)
            nl++;
    }
    return(nl);
}

/* -###-------------------------------------------------------------------- */
/*  map_line_straight(lona,lata,lonb,latb)                                  */
/*                                                                          */
/*  Draw straight line from (lona,lata) to (lonb,latb), restricted to the   */
/*  current clip region. Return 1 if at least some part of the line is      */
/*  drawn, otherwise return 0.                                              */
/*                                                                          */
/*  Note: a line is only drawn if both points are visible in the current    */
/*  hemisphere.                                                             */

int map_line_straight(TDAContext *ctx, double lona,double lata,double lonb,double latb)
{
    int ra,rb,ca,cb,dir;
    double ax,ay,bx,by,sax,say,sbx,sby;

    ax = lona;
    ay = lata;
    if (gd_proj_p(ctx, &ax,&ay)) {
        tda_out("first point  not visible\n");
        return(0);
    }
    bx = lonb;
    by = latb;
    if (gd_proj_p(ctx, &bx,&by)) {
        tda_out("second point  not visible\n");
        return(0);
    }
    if (clip_line(ctx, ax,ay,bx,by,ctx->GDPXA,ctx->GDPYA,ctx->GDPXB,ctx->GDPYB,&sax,&say,&sbx,&sby, 
        &ca,&cb,&ra,&rb,&dir) == 0)  
        return(0);
           
    ps_2dplot(ctx, sax,say,0);
    ps_2dplot(ctx, sbx,sby,1);
    fprintf(ctx->PSFd,"stroke\n");
    return(1);
}

/* -###-------------------------------------------------------------------- */
/*  map_geodesic_a(lona,lata,lonb,latb)                                     */
/*                                                                          */  
/*  Draw geodesic from (lona,lata) to (lonb,latb).                          */
/*  Return 1 if some part is successfully drawn, 0 if nothing is drawn,     */
/*  -1 if error.                                                            */

int map_geodesic_a(TDAContext *ctx, double lona,double lata,double lonb,double latb)
{
    int la,lb;
    double x,y,ax,ay,az,bx,by,bz,ra,rb,rd,rm,r1 = 0.0,r2 = 0.0;
    double apx,apy,apz,bpx,bpy,bpz,cpx,cpy,cpz;
    double b,d,rx,ry,rz,cx,cy,cz,nx,ny,nz,cosphi,psi;

    rm = 180.0 / Pi;
      
   
    if (fabs(lona - lonb) <= ctx->EPSI1 && fabs(lata - latb) <= ctx->EPSI1)
        return(0);
         
    ax = lona;
    ay = lata;
    az = 1.0;
    ps_3dgeo(ctx, &ax,&ay,&az);

    /** tda_out("ax=%g %g %g\n",ax,ay,az); **/ 

    ps_3dprj1(ctx, ax,ay,az,&apx,&apy,&apz);           

    bx = lonb;
    by = latb;
    bz = 1.0;
    ps_3dgeo(ctx, &bx,&by,&bz);
    ps_3dprj1(ctx, bx,by,bz,&bpx,&bpy,&bpz);           

    /**  tda_out("bx=%g %g %g bpx=%g %g %g\n",bx,by,bz,bpx,bpy,bpz); **/

    if (apz < 0.0 && bpz < 0.0)    /* both points not visible */
        return(0);

    rx = ctx->GDLon;
    ry = ctx->GDLat;
    rz = 1.0;
    ps_3dgeo(ctx, &rx,&ry,&rz);
    
    d = rx * bx + ry * by + rz * bz;

           
    if (fabs(d) <= ctx->EPSI1) {
        cx = bx;
        cy = by;
        cz = bz;
    }
    else {
        b = (rx * ax + ry * ay + rz * az);
        if (fabs(b - d) <= ctx->EPSI1)
            d = 1.0;
        else
            d = b / d;

        cx = ax - d * bx;
        cy = ay - d * by;
        cz = az - d * bz;
        d = sqrt(cx * cx + cy * cy + cz * cz);
        if (fabs(d) <= ctx->EPSI1) {
            err_geodesic(ctx, -1,lona,lata,lonb,latb);
            return(-1);
        }
        cx /= d;
        cy /= d;
        cz /= d;
    }
    /** tda_out("cx=%g %g %g\n",cx,cy,cz); **/

    /** g_spat(rx,ry,rz,cx,cy,cz,&dx,&dy,&dz); **/

    /** tda_out("dx=%g %g %g\n",dx,dy,dz); **/ 

    g_spat(ctx, bx,by,bz,ax,ay,az,&nx,&ny,&nz);


    d = sqrt(nx * nx + ny * ny + nz * nz);
    if (fabs(d) <= ctx->EPSI1) {
        err_geodesic(ctx, -3,lona,lata,lonb,latb);
        return(-3);
    }
    cosphi = (nx * rx + ny * ry + nz * rz) / d;
   
        tda_out("cosphi=%g\n",cosphi); 

    /* length of minor semi-axis */
                                                              
    b = fabs(cosphi); 

     tda_out("b=%g\n",b);

    ps_3dprj1(ctx, cx,cy,cz,&cpx,&cpy,&cpz);           
    if (fabs(cpz) > ctx->EPSI1 || (fabs(cpx) <= ctx->EPSI1 && fabs(cpy) <= ctx->EPSI1)) {
        err_geodesic(ctx, -2,lona,lata,lonb,latb);
        return(-2);
    }
    /** tda_out("cpx=%g cpy=%g cpz=%g \n",cpx,cpy,cpz); */

    psi = rm * atan2(cpy,cpx);
    /** tda_out("proj c... x=%g %g psi=%g\n",cpx,cpy,psi); **/

    la = g_left(ctx, cpx,cpy,-cpx,-cpy,apx,apy);
    lb = g_left(ctx, cpx,cpy,-cpx,-cpy,bpx,bpy);

    tda_out("la=%d lb=%d\n",la,lb); 

    /* calculate arc between a and c, and b and c */

    ra = (cx * ax + cy * ay + cz * az) / sqrt(ax * ax + ay * ay + az * az);
    rb = (cx * bx + cy * by + cz * bz) / sqrt(bx * bx + by * by + bz * bz);

    if (ra < -1.0)
        ra = -1.0;
    else if (ra > 1.0)
        ra = 1.0;
    if (rb < -1.0)
        rb = -1.0;
    else if (rb > 1.0)
        rb = 1.0;

    ra = acos(ra) * rm;
    rb = acos(rb) * rm;

    tda_out("\nra=%g rb=%g\n",ra,rb);

    /* check for visibility */ 


    if (apz >= 0.0 && bpz >= 0.0) {         /* both points visible */
        if (la == 0 && lb == 0) {           /* both to the right */
            r1 = dmin(ctx, ra,rb);
            r2 = dmax(ctx, ra,rb);
        }
        else {                              /* both to the left */
            ra = 360.0 - ra;
            rb = 360.0 - rb;
            r1 = dmin(ctx, ra,rb);
            r2 = dmax(ctx, ra,rb);
        }
    }          
    if (apz >= 0.0 && bpz < 0.0) {          /* only point a is visible */
        if (la == 0) {
            r1 = 0.0;
            r2 = ra;
        }
        else {
            r1 = 360.0 - ra;
            r2 = 360.0;
        }
    }
    if (bpz >= 0.0 && apz < 0.0) {                  /* only point b is visible */
        if (lb == 0) {
            r1 = 0.0;
            r2 = rb;
        }
        else {
            r1 = 360.0 - rb;
            r2 = 360.0;
        }
    }

    tda_out("\nr1=%g r2=%g b=%g\n",r1,r2,b);

    if (fabs(r1 - r2) <= ctx->EPSI1)  
        return(0);
          
/*  rd = (ps_2dx(GDPXB) - ps_2dx(GDPXA)) / 2.0;     */


    rd = (ps_2dx(ctx, 1.0) - ps_2dx(ctx, -1.0)) / 2.0;
    x = ps_2dx(ctx, ctx->GDPX);   
    y = ps_2dy(ctx, ctx->GDPY);

    fprintf(ctx->PSFd,"gsave\n");
    fprintf(ctx->PSFd,"%7.3f %7.3f translate\n",x,y);
    fprintf(ctx->PSFd,"%7.4f rotate\n",psi);
    fprintf(ctx->PSFd,"1.0 %7.4f scale\n",b);   
    fprintf(ctx->PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",rd,r1,r2);
    fprintf(ctx->PSFd,"stroke\n");
    fprintf(ctx->PSFd,"grestore\n");
    return(1);
}

/* -###-------------------------------------------------------------------- */
/*  err_geodesic(n,lona,lata,lonb,latb)                                     */

void err_geodesic(TDAContext *ctx, int n,double lona,double lata,double lonb,double latb)
{
    printf1(ctx, "Error (%d) in calculation of geodesic: (%g,%g), (%g,%g)\n",
        n,lona,lata,lonb,latb);
}

/* ------------------------------------------------------------------------ */
/*  map_polygon(n,x,y,sdid)   Plot polygon with geographical coordinates.   */
/*                                                                          */
/*  x[i], y[i] (i=0,n-1) contains the points.                               */
    
int map_polygon(TDAContext *ctx, int n,double *x,double *y,int sdid)
{
    int np = 0;

    switch (ctx->PSPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            np = map_polygon_cyl(ctx, n,x,y,sdid);
            break;



        case GDPROJ2:
               map_polygon_p2(ctx, n,x,y);
            break;

        default:
            tda_out("Projection not supported.\n");
            break;  
    }
    return(np);
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_p2(n,x,y)   Called by map_polygon for projection type 2.    */
    
void map_polygon_p2(TDAContext *ctx, int n,double *xx,double *yy)
{
    register int i;
    int first,r,vflag,nflag,npi;
    double x,y,xf,yf,xl,yl,xn,yn;

    i = 0;
    while (i < n) {         /* find first visible point */ 
        x = xx[i];
        y = yy[i];
        r = gd_proj_p(ctx, &x,&y);
        if (r == 0)
            break;
        i++;
    }
    if (i >= n - 1)
        return;

    xf = x;
    yf = y;
    npi = nflag = vflag = first = 0;
    for ( ; i < n; ++i) {
        x = xx[i];
        y = yy[i];
        r = gd_proj_p(ctx, &x,&y);
        if (r) {
            if (nflag == 0) {
                xn = xl;
                yn = yl;
                nflag = 1;
            }
            continue;
        }
        if (nflag) {
            if (npi >= 2)
                plsd_arc(ctx, xn,yn,x,y);           
            npi = nflag = 0;
        }
        else {
            ps_2dplot(ctx, x,y,first);
            first = 1;
            npi++;
        }
        xl = x;
        yl = y; 

    }
    if (first) {
        if (nflag)   
            plsd_arc(ctx, xn,yn,xf,yf);           

        if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
            fprintf(ctx->PSFd,"closepath\ngsave\n%4.2f setgray\n",ctx->PMGS);    
            fprintf(ctx->PSFd,"fill\ngrestore\n");
        }
        if (ctx->PMLW > 0.0)
            fprintf(ctx->PSFd,"stroke\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  plsd_arc(x0,y0,x1,y1)   Draw arc from (x0,y0) to (x1,y1). Assume that   */
/*                          points are given in user coordinates.           */
    
void plsd_arc(TDAContext *ctx, double x0,double y0,double x1,double y1)
{
    double a,b,d;

    if ((a = atan2(y0,x0) * 180.0 / Pi) < 0.0) {
        /** fprintf(PSFd,"%%a=%g ",a); **/
        a += 360.0;
    }
    /** fprintf(PSFd,"%%a=%g\n",a); **/

    if ((b = atan2(y1,x1) * 180.0 / Pi) < 0.0)  {
        /** fprintf(PSFd,"%%b=%g ",b); **/
        b += 360.0;
    }
    /** fprintf(PSFd,"%%b=%g\n",b); **/

    d = fabs(a - b);
    if (d <= ctx->EPSI1)
        return;

    fprintf(ctx->PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f ",ctx->GDPXPS,ctx->GDPYPS,ctx->GDPSRadius,a,b);
    if (d <= 180.0) {
        if (a <= b)
            fprintf(ctx->PSFd,"arc\n");
        else
            fprintf(ctx->PSFd,"arcn\n");
    }
    else {
        if (a <= b)
            fprintf(ctx->PSFd,"arcn\n");
        else
            fprintf(ctx->PSFd,"arc\n");
    }
}
 
/* -###-------------------------------------------------------------------- */
/*  map_polygon_cyl(n,x,y,sdid)                                             */
/*                                                                          */
/*  Plot polygon with geographical coordinates. Called by map_polygon for   */
/*  cylindrical projections.                                                */
/*                                                                          */
/*  sdid is the SDID of the polygon (used for error messages).              */
/*                                                                          */
/*  Return number of polygons plotted, or -1 if error.                      */

int map_polygon_cyl(TDAContext *ctx, int n,double *x,double *y,int sdid)
{
    register int i,j;
    int rt,r,txflag,pflag,imin,imax;
    double x0,y0,x1,y1,xmin,xmax,ymin,ymax;

    pflag = map_check_pole(ctx, n,x,y);


/**   
    tda_out("\nNew polygon\n");
    for (i = 0; i < n; ++i)
        tda_out("i=%4d x=%f y=%f\n",i,x[i],y[i]);
**/   

    /* transformation into a planar polygon */

    j = -1;         /* find one point inside plot region */

    for (i = 0; i < n; ++i) {
        if (check_inside(ctx, x[i],y[i])) {
            j = i;
            break;
        }
    }
    if (j < 0)
        return(0);

    xmin = xmax = x0 = x[j];
    ymin = ymax = y0 = y[j];
    imin = imax = j;


    for (i = j + 1; i < n; ++i) {
        x1 = x[i];
        y1 = y[i];
        map_find_n(ctx, x0,y0,&x1,&y1);
        x0 = x[i] = x1;
        y0 = y[i] = y1;
        if (x0 < xmin) {
            xmin = x0;
            imin = i;
        }
        if (x0 > xmax) {
            xmax = x0;
            imax = i;
        }
        ymin = dmin(ctx, y0,ymin);
        ymax = dmax(ctx, y0,ymax);
    }        
    for (i = j - 1; i >= 0; --i) {
        x1 = x[i];
        y1 = y[i];
        map_find_n(ctx, x0,y0,&x1,&y1);
        x0 = x[i] = x1;
        y0 = y[i] = y1;
        if (x0 < xmin) {
            xmin = x0;
            imin = i;
        }
        if (x0 > xmax) {
            xmax = x0;
            imax = i;
        }
        ymin = dmin(ctx, y0,ymin);
        ymax = dmax(ctx, y0,ymax);
    }        





    if (pflag != 0) {               /* add points for polar regions */
        if (pflag == 1)
            y0 = 90.0;
        else
            y0 = -90.0;

        if (imin == 0 && imax == n - 1) {
            if (xmax < xmin + 360.0 - ctx->EPSI1) {
                x[n] = xmin + 360.0;
                y[n] = y[0];
                n++;
                xmax = xmin + 360.0;
            }
            x[n] = x[n - 1];
            y[n] = y0;        
            n++;
            x[n] = x[0];
            y[n] = y0;
            n++;
        }
        else if (imin == n - 1 && imax == 0) {
            if (xmin > xmax - 360.0 + ctx->EPSI1) {
                x[n] = xmax - 360.0;
                y[n] = y[0];
                n++;
                xmin = xmax - 360.0;
            }
            x[n] = x[n - 1];
            y[n] = y0;        
            n++;
            x[n] = x[0];
            y[n] = y0;
            n++;
        }
        else {
            printf1(ctx, "Error: cannot plot polygon with ID %d containing a pole (%d).\n",sdid,pflag);
            printf1(ctx, "Will be skipped.\n");
            return(0);
        }
        ymin = dmin(ctx, ymin,y0);
        ymax = dmax(ctx, ymax,y0);
    }
/**   
    tda_out("\nNew polygon nach Umwandlung imin=%d imax=%d\n",imin,imax);
    for (i = 0; i < n; ++i)
        tda_out("i=%4d x=%f y=%f\n",i,x[i],y[i]);
      
    tda_out("xmin=%g xmax=%g\n",xmin,xmax);
    tda_out("ymin=%g ymax=%g\n",ymin,ymax);
**/         
    if (ymin >= ctx->GDLatB || ymax <= ctx->GDLatA)
        return(0);

    if (xmax <= ctx->GDPXA) {
        if (xmin + 360.0 >= ctx->GDPXB)
            return(0);

        for (i = 0; i < n; ++i)
            x[i] += 360.0;

        xmin += 360.0;
        xmax += 360.0;
    }
    if (xmin >= ctx->GDPXB) {
        if (xmax - 360.0 <= ctx->GDPXA)
            return(0);

        for (i = 0; i < n; ++i)
            x[i] -= 360.0;

        xmin -= 360.0;
        xmax -= 360.0;
    }

    txflag = 0;
    if (xmax > ctx->GDPXB && xmax - 360.0 > ctx->GDPXA) {
        txflag = -1;
    }
    else if (xmin < ctx->GDPXA && xmin + 360.0 < ctx->GDPXB) {
        txflag = 1;
    }


    /* upstream debug prints removed here (per-vertex dump of the
       shifted polygon and its transformation flag); they wrote into
       the protocol on every clipped polygon */
    rt = 0;
    r = map_polygon_cyl_clip(ctx, n,x,y,sdid);
    if (r < 0)
        return(r);
    rt += r;

    if (txflag) {
tda_out("\nHier begin zweiter Teil txflag=%d\n",txflag);

        x0 = (double)txflag * 360.0;
        for (i = 0; i < n; ++i)  
            x[i] += x0;

        r = map_polygon_cyl_clip(ctx, n,x,y,sdid);
        if (r < 0)
            return(r);
        rt += r;
    }
    /* one more upstream debug print (fragment count) removed */


    return(rt);
}

/* -###-------------------------------------------------------------------- */
/*  map_polygon_cyl_clip(n,x,y,sdid)                                        */
/*                                                                          */
/*  Clip polygon defined by x[i], y[i] (i = 0,n-1) to rectangular clip      */
/*  region defined by GDPXA,GDLatA,GDPXB,GDLatB. Then plot the clipped      */
/*  polygon. sdid is the SDID of the polygon (used for error messages).     */
/*                                                                          */
/*  Return -1 if error, or number of polygons plotted.                      */

int map_polygon_cyl_clip(TDAContext *ctx, int n,double *x,double *y,int sdid)
{
    int m,r;  

    if (alloc_acu(ctx, 2 * n + 10))
        return(-1);              
    if (alloc_acv(ctx, 2 * n + 10))
        return(-1);                  
    if (alloc_acn(ctx, 2 * n + 10))
        return(-1);                      
    if (alloc_acr(ctx, 2 * n + 10))
        return(-1);                      
      
    m = sdclip_poly_rec(ctx, n,x,y,ctx->AcU,ctx->AcV,ctx->AcN,ctx->AcR,ctx->GDPXA,ctx->GDLatA,ctx->GDPXB,ctx->GDLatB);

    /* an upstream debug print sat here (m===) and would have polluted
       the protocol on the first partially-clipped polygon */

    if (m <= 0)
        return(m);
     
    if (m == 1) {                     /* polygon completely in clip region */
        map_polygon_cyl_ppol(ctx, n,x,y);
        return(1);
    }
    else if (m == 2) {                /* clip region completely in polygon */
        return(0);
        /***********
        map_polygon_cyl_pbox(ctx, GDPXA,GDLatA,GDPXB,GDLatB);
        return(1);
        *************/
    }
    else {
        r = map_polygon_cyl_plist(ctx, m,ctx->AcU,ctx->AcV,ctx->AcN,ctx->AcR,sdid);
        if (r < 0) {
            sdclip_err_msg(ctx, 1,sdid);
            r = 0;
        }
        return(r);
    }
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_pbox(xmin,ymin,xmax,ymax)                               */
/*                                                                          */
/*  Plot rectangular clip region.                                           */

void map_polygon_cyl_pbox(TDAContext *ctx, double xmin,double ymin,double xmax,double ymax)
{
    ps_2dplot(ctx, xmin,ymin,0);
    ps_2dplot(ctx, xmax,ymin,1);
    ps_2dplot(ctx, xmax,ymax,1);
    ps_2dplot(ctx, xmin,ymax,1);
    ps_2dplot(ctx, xmin,ymin,1);
    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n");
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_ppol(n,x,y)                                             */
/*                                                                          */
/*  Plot polygon x[i], y[i] (i=0,ln-1).                                     */

void map_polygon_cyl_ppol(TDAContext *ctx, int n,double *x,double *y)
{
    register int i;
    double yp;

    yp = gd_proj_gety(ctx, y[0]);  
    ps_2dplot(ctx, x[0],yp,0);

    for (i = 1; i < n; ++i) {
        yp = gd_proj_gety(ctx, y[i]);  
        ps_2dplot(ctx, x[i],yp,1);
    }
    yp = gd_proj_gety(ctx, y[0]);  
    ps_2dplot(ctx, x[0],yp,1);

    if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
        fprintf(ctx->PSFd,"gsave\n");
        ps_fill(ctx, ctx->PMGS);
        fprintf(ctx->PSFd,"grestore\n");
    }
    if (ctx->PMLW > 0.0)
        fprintf(ctx->PSFd,"stroke\n");
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_plist(n,x,y,nf,d,sdid)                                  */
/*                                                                          */
/*  This function gets an edge list created by sdclip_poly_rec and plots    */
/*  all polygons. Return number of polygons plotted, or -1 if error.        */

int map_polygon_cyl_plist(TDAContext *ctx, int n,double *x,double *y,int *nf,int *d,int sdid)
{
    (void)sdid;        /* unused: the signature is shared */
    register int i,j,k;
    int m,np,first;
    double yp;

    /************* 
    tda_out("map_polygon_plist\n");
    for (i = 0; i < n; ++i) {
        tda_out("i=%3d x=%f y=%f nf=%3d d=%d\n",i,x[i],y[i],nf[i],d[i]);
    }
    ********/

    np = 0;
    for (i = 0; i < n; ++i) {
        if (d[i] == 2)
            d[i] = -2;
        else
            d[i] = -1;
    }

    /* first check whether the edge list represents a set of valid
       simple polygons. */

    while (1) {
        j = -1;
        for (i = 0; i < n; ++i) {
            if (d[i] == -1) {
                j = i;
                break;
            }
        }
        if (j < 0)
            break;

        m = 1;
        k = j;
        d[k] = j;

        while (1) {
            k = nf[k];
            if (k == j)
                break;

            if (++m > n || d[k] == j) {
                m = 0;
                break;
            }
            d[k] = j;
        }
        if (m == 0)             /* no valid simple polygon */
            return(-1);
    }

    /* now write to output file */

    for (i = 0; i < n; ++i) {
        if (d[i] != -2)
            d[i] = -1;
    }

    while (1) {
        j = -1;
        for (i = 0; i < n; ++i) {
            if (d[i] == -1) {
                j = i;
                break;
            }
        }
        if (j < 0)
            break;

        m = 1;
        k = j;
        d[k] = j;

        while (1) {
            k = nf[k];
            if (k == j)
                break;
            if (++m > n || d[k] == j) {
                m = 0;
                break;
            }
            d[k] = j;
        }
        if (m == 0)         
            return(-1);

        if (m < 3)
            continue;

        first = 0;
        k = j;
        while (1) {
            yp = gd_proj_gety(ctx, y[k]);  
            ps_2dplot(ctx, x[k],yp,first);
            first = 1;
            k = nf[k];
            if (k == j)
                break;
        }
        if (first) {
            yp = gd_proj_gety(ctx, y[j]);  
            ps_2dplot(ctx, x[j],yp,first);

            if (ctx->PMGSFlg && ctx->PMGS >= 0.0 && ctx->PMGS <= 1.0) {
                fprintf(ctx->PSFd,"gsave\n");
                ps_fill(ctx, ctx->PMGS);
                fprintf(ctx->PSFd,"grestore\n");
            }
            if (ctx->PMLW > 0.0)
                fprintf(ctx->PSFd,"stroke\n");
            np++;
        }
    }
    return(np);
}


/* ------------------------------------------------------------------------ */
/*  map_check_pole(n,x,y)                                                   */
/*                                                                          */
/*  Return 1 if north pole inside polygon, -1 if south pole inside polygon, */
/*  otherwise 0. Method: use parallel projection of polygon, then standard  */ 
/*  inside check.                                                           */

int map_check_pole(TDAContext *ctx, int n,double *x,double *y)
{
    (void)ctx;        /* unused: the signature is shared */
    register int i;
    int c,m;
    double p,li,lj,bi,bj,xi,xj,yi,yj;

    p = Pi / 180.0;
    m = c = 0;

    lj = x[n - 1] * p;
    bj = y[n - 1] * p;
    xj = sin(lj) * cos(bj);
    yj = cos(lj) * cos(bj);

    for (i = 0; i < n; ++i) {
        if (y[i] > 0.0)
            m++;
        li = x[i] * p;
        bi = y[i] * p;
        xi = sin(li) * cos(bi);
        yi = cos(li) * cos(bi);

        if ((((yi <= 0.0) && (0.0 < yj)) ||
             ((yj <= 0.0) && (0.0 < yi))) &&
             (0.0 < (xj - xi) * (0.0 - yi) / (yj - yi) + xi))  
            c = !c;

        xj = xi;
        yj = yi;
    }

    if (c) {        /* decide whether north or south pole */
        if (m > n - m)
            c = 1;
        else
            c = -1;
    }
    return(c);
}


/* ------------------------------------------------------------------------ */
/*  geod_distance(lon1,lat1,lon2,lat2)                                      */
/*                                                                          */
/*  Return the great circle distance between (lon1,lat1) and (lon2,lat2).   */

double geod_distance(TDAContext *ctx, double lon1,double lat1,double lon2,double lat2)
{
    double d,slon,slat,t;

    lon1 = degree_to_arc(ctx, lon1);
    lat1 = degree_to_arc(ctx, lat1);
    lon2 = degree_to_arc(ctx, lon2);
    lat2 = degree_to_arc(ctx, lat2);
    slon = sin((lon1 - lon2) / 2.0);
    slat = sin((lat1 - lat2) / 2.0);
    t = sqrt(slat * slat + cos(lat1) * cos(lat2) * slon * slon);
    if (t > 1.0)
        t = 1.0;
    else if (t < -1.0)
        t = -1.0;
    d = 2.0 * asin(t);
    return(d);
}



