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

/*  functions in t_map.c */

int psetupg(void); 
int check_geo(int opt);
void gd_view(int opt);
int gd_check_inverse(void);
double gd_adlon(double lon);
double gd_adlat(double lat);
int gd_proj_p(double *x,double *y);
int gd_proj_inv(double *x,double *y);
int check_x_in_region(double x);
int check_y_in_region(double y);
int check_inside(double lon,double lat);
int check_clip(double x,double y);
double gd_proj_gety(double lat);
int gd_proj_getxy(double *x,double *y);
double gd_adjust_lon(double lon);

int sdpgrat(void); 
void mapgrid_bounds(void);
void mapgrid_bounds_r(void);
void mapgrid_arc(double x,double y,double a,double b,double r,int gflag);
void mapgrid_lon(double sc);
void mapgrid_line(double xa,double ya,double xb,double yb);
void mapgrid_label(int opt,double lab,double x,double y,double fs,int m);
void mapgrid_lon_azi(double lon);
void mapgrid_lon_azi_lab0(double lab,double lon,double fs);
void mapgrid_lon_azi_lab1(double lab,double lon,double fs);
void mapgrid_lon_azi1(double lon,double sc);
void mapgrid_lat(void);
void mapgrid_lat_a0(double lat);
void mapgrid_lat_a3(double lat);

int sdpgeo(void); 
int sdpmap(void); 
int map_point(double x,double y,int nc);
int map_point_cyl(double x,double y);
int map_point_azi(double x,double y,int nc);

int map_line(int n,double *x,double *y,int gflag);
int map_line_cyl_s(int n,double *x,double *y);
void map_find_n(double x,double y,double *xa,double *ya);
int map_line_cyl_g(int n,double *x,double *y);
int map_line_cyl_geod(double lon1,double lat1,double lon2,double lat2);
void map_line_clip(double xa,double ya,double xb,double yb);
int map_line_azi(int n,double *x,double *y,int gflag);
int map_line_straight(double lona,double lata,double lonb,double latb);
int map_geodesic_a(double lona,double lata,double lonb,double latb);
void err_geodesic(int n,double lona,double lata,double lonb,double latb);

int map_polygon(int n,double *x,double *y,int sdid);
void map_polygon_p2(int n,double *xx,double *yy);
void plsd_arc(double x0,double y0,double x1,double y1);
int map_polygon_cyl(int n,double *x,double *y,int sdid);
int map_polygon_cyl_clip(int n,double *x,double *y,int sdid);
void map_polygon_cyl_pbox(double xmin,double ymin,double xmax,double ymax);
void map_polygon_cyl_ppol(int n,double *x,double *y);
int map_polygon_cyl_plist(int n,double *x,double *y,int *nf,int *d,int sdid);
int map_check_pole(int n,double *x,double *y);

double geod_distance(double lon1,double lat1,double lon2,double lat2);

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

int GDPROJ = 0;             /* type of projection (duplicated in PSPROJ)    */
int GDPROJA = 0;            /* records special cases for azimuthal proj.    */
double GDRHem = 0.0;        /* size of hemisphere for azimuthal proj.       */
double GDRHemP = 0.0;       /* cos(90 - GDRHem)                             */

double GDLon = 0.0;         /* defined by view parameter                    */
double GDLat = 0.0;
double GDLonD = 0.0;        /* defined by region parameter                  */
double GDLatD = 0.0;         
double GDLonA = 0.0;
double GDLonB = 0.0;
double GDLatA = 0.0;
double GDLatB = 0.0;
double GDLonXAYA = 0.0;
double GDLatXAYA = 0.0;
double GDLonXBYA = 0.0;
double GDLatXBYA = 0.0;
double GDLonXAYB = 0.0;
double GDLatXAYB = 0.0;
double GDLonXBYB = 0.0;
double GDLatXBYB = 0.0;

double GDLonMin = 0.0;      /* extension of selected region: longitudes     */
double GDLonMax = 0.0;      /* used for inside checks.                      */
double GDLatMin = 0.0;      /* extension of selected region: latitudes      */
double GDLatMax = 0.0;      /* restricted to -90 and +90.                   */

double GDPX  = 0.0;         /* projected coordinates of selected region     */
double GDPY  = 0.0;
double GDPXA = 0.0;     
double GDPYA = 0.0;
double GDPXB = 0.0;     
double GDPYB = 0.0;

double GDPXPS = 0.0;        /* PostScript coordinates of center             */
double GDPYPS = 0.0;
double GDPSRadius = 1.0;    /* radius in PostScript coordinates             */

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

int psetupg(void)  
{
   	int	err,n,r;            
    double x,y,ymin,ymax;
      
    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Set up geographical coordinate system. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,3,0)) {     /* get parameters */
        goto PSGFin;
    }   
    if (PMPROJ == 0)                /* default projection */
        PMPROJ = 10;

    if (PSFFlg == 0) {
        printf1("Error: need a PostScript file.\n");
        return(-1);
    }
    PSPROJ = 0;                 /* make current system invalid */
    GDLon = PMViewLon;          /* center of projection */
    if (GDLon <= -180.0)
        GDLon = 180.0;
    GDLat = PMViewLat;
    GDPROJA = 0;
           
    switch (PMPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:

            if (PMRegionFlg == 0) {
                if (GDLat == 0) {
                    PMRegion1 = 180.0;
                    PMRegion2 =  90.0;
                }   
                else {
                    printf1("Error: need region parameter.\n");
                    goto PSGFin;
                }
            }
            if (PMRegion1 < EPSI1 || PMRegion1 > 180.0 || PMRegion2 < EPSI1 || PMRegion2 > 90.0) {   
                printf1("Error in region parameter.\n");
                goto PSGFin;
            }
            GDLonD = PMRegion1;
            GDLatD = PMRegion2;

            if (PMPROJ == GDPROJ10) {
                GDPROJ = GDPROJ10;
                printf1("Cylindrical projection (equidistant).\n");
                gd_view(0);
                ymin = -90.0;
                ymax =  90.0;
            }
            else if (PMPROJ == GDPROJ11) {
                GDPROJ = GDPROJ11;
                printf1("Cylindrical projection (equal-area).\n");
                gd_view(0);
                ymin = -90.0;
                ymax =  90.0;
            }
            if (PMPROJ == GDPROJ12) {
                GDPROJ = GDPROJ12;
                printf1("Cylindrical projection (Mercator).\n");
                gd_view(0);

                if (GDLat < -85.0 || GDLat > 85.0) {
                    printf1("Error in view (lat) parameter.\n");
                    goto PSGFin;
                }
                ymin = -89.0;
                ymax =  89.0;
            }
            GDLonMin = GDLonA = gd_adlon(GDLon - GDLonD);
            GDLonMax = GDLonB = gd_adlon(GDLon + GDLonD);

            GDLatMin = GDLatA = dmax(GDLat - GDLatD,ymin);
            GDLatMax = GDLatB = dmin(GDLat + GDLatD,ymax);
            GDLatD = (GDLatB - GDLatA) / 2.0;

            GDPX = GDLon;
            GDPXA = gd_adlon(GDPX - GDLonD);
            GDPXB = GDPXA + 2.0 * GDLonD;

            GDPY = gd_proj_gety(GDLat);
            GDPYA = gd_proj_gety(GDLatA);
            GDPYB = gd_proj_gety(GDLatB);

            PXLen = PMXLen;
            if (PMYLenFlg)
                PYLen = PMYLen;
            else {
                if (PMPROJ == GDPROJ12)
                    PYLen = PXLen;
                else
                    PYLen = PXLen * (GDPYB - GDPYA) / (GDPXB - GDPXA);
            }
            printf1("Selected region longitude: %9.4f (%9.4f) %9.4f\n",GDLonA,GDLon,GDLonB);
            printf1("Selected region latitude:  %9.4f (%9.4f) %9.4f\n",GDLatA,GDLat,GDLatB);
            break;                          

        case GDPROJ20:
            GDPROJ = GDPROJ20;
            printf1("Azimuthal projection (orthographic).\n");

            if (PMRegionFlg != 0) {
                GDPROJA = 0;
                if (PMRegion1 < EPSI1 || PMRegion1 > 180.0 || PMRegion2 < EPSI1 || PMRegion2 > 90.0) {   
                    printf1("Error in region parameter.\n");
                    goto PSGFin;
                }
                GDLonD = PMRegion1;
                GDLatD = PMRegion2;

                printf1("View: %g,%g. Rectangular region (%g,%g).\n",GDLon,GDLat,GDLonD,GDLatD);

                GDLonA = gd_adlon(GDLon - GDLonD);
                GDLonB = gd_adlon(GDLon + GDLonD);
                GDLatA = GDLat - GDLatD;
                GDLatB = GDLat + GDLatD;
                if (GDLatA < -90.0 || GDLatB > 90.0) {
                    printf1("The selected region is not supported.\n");
                    goto PSGFin;
                }
                setup_3proj(GDLon,GDLat);   /* set up projection matrix */

                GDPX = GDPY = 0.0;
                GDPXA = GDLon - GDLonD; y = GDLat;
                gd_proj_p(&GDPXA,&y);
                GDPXB = GDLon + GDLonD; y = GDLat;
                gd_proj_p(&GDPXB,&y);

                GDPYA = GDLat - GDLatD; x = GDLon;
                gd_proj_p(&x,&GDPYA);
                GDPYB = GDLat + GDLatD; x = GDLon;
                gd_proj_p(&x,&GDPYB);

/**               
                printf("GDPX=%g GDPY=%g\n",GDPX,GDPY);
                printf("xa=%g xb=%g\n",GDPXA,GDPXB);
                printf("ya=%g yb=%g\n",GDPYA,GDPYB);
**/              
                /* calculate coordinates of corner points and check
                   visibility. */ 

                GDLonXAYA = GDPXA;
                GDLatXAYA = GDPYA;
                GDLonXAYB = GDPXA;
                GDLatXAYB = GDPYB;
                GDLonXBYA = GDPXB;
                GDLatXBYA = GDPYA;
                GDLonXBYB = GDPXB;
                GDLatXBYB = GDPYB;

                if (gd_proj_inv(&GDLonXAYA,&GDLatXAYA) == 0 ||
                    gd_proj_inv(&GDLonXAYB,&GDLatXAYB) == 0 ||
                    gd_proj_inv(&GDLonXBYA,&GDLatXBYA) == 0 ||
                    gd_proj_inv(&GDLonXBYB,&GDLatXBYB) == 0) {
                    printf1("Error: selected region not completely visible.\n");
                    goto PSGFin;
                }                  
                printf1("\nSelected region (lon,lat)\n");
                printf1("(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    GDLonXAYB,GDLatXAYB,GDLon,GDLatB,GDLonXBYB,GDLatXBYB);
                printf1("(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    GDLonA,GDLat,GDLon,GDLat,GDLonB,GDLat);
                printf1("(%9.4f,%9.4f) (%9.4f,%9.4f) (%9.4f,%9.4f)\n",
                    GDLonXAYA,GDLatXAYA,GDLon,GDLatA,GDLonXBYA,GDLatXBYA);

                if (GDLat >= 0.0) {
                    GDLatMin = GDLatXAYA;
                    GDLatMax = GDLatB;
                    GDLonMin = GDLonXAYB;
                    GDLonMax = GDLonXBYB;
                }
                else {
                    GDLatMin = GDLatA;
                    GDLatMax = GDLatXAYB;
                    GDLonMin = GDLonXAYA;
                    GDLonMax = GDLonXBYA;
                }
                printf("GDlatmin=%g %g\n",GDLatMin,GDLatMax);
                printf("GDlonmin=%g %g\n",GDLonMin,GDLonMax);





/***
                x = GDLonXAYA;
                y = GDLatXAYA;
                printf("x=%g y=%g ",x,y);
                gd_proj_p(&x,&y);
                printf(" -> %g %g\n",x,y);

                x = GDLonXBYA;
                y = GDLatXBYA;
                printf("x=%g y=%g ",x,y);
                gd_proj_p(&x,&y);
                printf(" -> %g %g\n",x,y);

                x = GDLonXAYB;
                y = GDLatXAYB;
                printf("x=%g y=%g ",x,y);
                gd_proj_p(&x,&y);
                printf(" -> %g %g\n",x,y);

                x = GDLonXBYB;
                y = GDLatXBYB;
                printf("x=%g y=%g ",x,y);
                gd_proj_p(&x,&y);
                printf(" -> %g %g\n",x,y);
**/
               

                x = fabs(GDPXB - GDPXA);
                y = fabs(GDPYB - GDPYA);

                PXLen = PMXLen;
                if (PMYLenFlg)
                    PYLen = PMYLen;
                else 
                    PYLen = PXLen * y / x;                   


            }
            else {
                setup_3proj(GDLon,GDLat);   /* set up projection matrix */

                GDRHem = PMRHem;
                if (GDRHem > 90.0)
                    GDRHem = 90.0;

                if (fabs(GDLat - 90.0) <= EPSI1) {  /* view from north pole */
                    GDPROJA = 1;
                    GDLat = 90.0;
                    GDLonMin = -180.0;
                    GDLonMax =  180.0;
                }
                else if (fabs(GDLat + 90.0) <= EPSI1) {  /* view from south pole */
                    GDPROJA = 2;
                    GDLat = -90.0;
                    GDLonMin = -180.0;
                    GDLonMax =  180.0;
                }
                else {
                    GDPROJA = 3;
                    GDRHem = 90.0;
                }
                printf1("View: %g,%g. Circular area. Size (rhem) = %g\n\n",GDLon,GDLat,GDRHem);

                GDLonD = 90.0;
                GDLonMin = GDLonA = gd_adlon(GDLon - GDLonD);
                GDLonMax = GDLonB = gd_adlon(GDLon + GDLonD);

printf("GDLonD=%g GDLonMin=%g GDLonMax=%g\n",GDLonD,GDLonMin,GDLonMax);



                x = degree_to_arc(90.0 - GDRHem);
                GDRHemP = cos(x);

                GDPX = GDPY = 0.0;
                GDPXA = GDPYA = -GDRHemP;
                GDPXB = GDPYB =  GDRHemP;
                GDLonD = GDLatD = GDRHem;

                PXLen = PMXLen;
                if (PMYLenFlg)
                    PYLen = PMYLen;
                else 
                    PYLen = PXLen;
            }

            break;




        case GDPROJ1:
            GDPROJ = GDPROJ1;
            printf1("Parallel projection (rectangular region).\n");
            gd_view(0);

            setup_3proj(GDLon,GDLat);   /* set up projection matrix */

            GDLonMin = GDLonA = gd_adlon(GDLon - GDLonD);
            GDLonMax = GDLonB = gd_adlon(GDLon + GDLonD);

            GDLatMin = dmax(GDLat - GDLatD,-90.0);
            GDLatMax = dmin(GDLat + GDLatD, 90.0);

            if (GDLat + GDLatD > 90.0 || GDLat - GDLatD < -90.0) {
                GDLonMin = -180.0;
                GDLonMax =  180.0;
            }
            GDLatA = gd_adlat(GDLat - GDLatD);
            GDLatB = gd_adlat(GDLat + GDLatD);

            GDPX = GDPY = 0.0;

            GDPXA = GDLon - GDLonD; y = GDLat;
            gd_proj_p(&GDPXA,&y);

            GDPXB = GDLon + GDLonD; y = GDLat;
            gd_proj_p(&GDPXB,&y);

            GDPYA = GDLat - GDLatD; x = GDLon;
            gd_proj_p(&x,&GDPYA);

            GDPYB = GDLat + GDLatD; x = GDLon;
            gd_proj_p(&x,&GDPYB);

            PXLen = PMXLen;
            if (PMYLenFlg)
                PYLen = PMYLen;
            else 
                PYLen = PXLen * (GDPYB - GDPYA) / (GDPXB - GDPXA);

            if (gd_check_inverse() == 0) {
                printf1("Error: rectangle exceeds limits.\n");
                goto PSGFin;
            }
            break;                          

        case GDPROJ2:
            GDPROJ = GDPROJ2;
            GDLonD = GDLatD = 90.0;
            printf1("Parallel projection (hemisphere).\n");
            gd_view(1);
        
            setup_3proj(GDLon,GDLat);   /* set up projection matrix */

            GDLonA = gd_adlon(GDLon - GDLonD);
            GDLonB = gd_adlon(GDLon + GDLonD);

            GDLonMin = -180.0;
            GDLonMax =  180.0;

            GDLatMin = dmax(GDLat - GDLatD,-90.0);
            GDLatMax = dmin(GDLat + GDLatD, 90.0);

            GDLatA = gd_adlat(GDLat - GDLatD);
            GDLatB = gd_adlat(GDLat + GDLatD);

            GDPX = GDPY = 0.0;
            GDPXA = GDPYA = -1.0;
            GDPXB = GDPYB =  1.0;

            PXLen = PMXLen;
            if (PMYLenFlg)
                PYLen = PMYLen;
            else 
                PYLen = PXLen;
            break;


        default:                            
            printf1("Error: unknown projection %d.\n",PMPROJ);
            goto PSGFin;
    }

    /* create 2d coordinate system */

    newline();
    PSLog[0] = PSLog[1] = 0;    /* no logarithmic axes */
    XOrg = PMXOrg;              /* origin of PostScript figure */
    YOrg = PMYOrg;
    SCALXFac = 1.0;             /* PostScript x scaling factor */
    SCALYFac = 1.0;             /* PostScript y scaling factor */
    ROTFac = PMPSRot;           /* PostScript rotation factor */
    PA1[0] = GDPXA;            
    PA2[0] = GDPXB;            
    PA1[1] = GDPYA;              
    PA2[1] = GDPYB;              
    setup_2dsys();

    GDPXPS = ps_2dx(GDPX);      /* center in PostScript coordinates */
    GDPYPS = ps_2dy(GDPY);
    GDPSRadius = ps_2dx(GDPX) - ps_2dx(GDPXA);

    PSPROJ = GDPROJ;        /* make valid projection globally known */
    err = 0;

PSGFin:
    p_clean();
    return(err);
}

/* ------------------------------------------------------------------------ */
/*  check_geo(opt)  Return 0 if a valid geographical coordinate system      */
/*                  exists, otherwise return -1. If opt != 0 error message. */

int check_geo(int opt)
{
    if (PSPROJ)
        return(0);

    if (opt)
        printf1("Error: no geographical coordinate system.\n");
    return(-1);
}

/* ------------------------------------------------------------------------ */
/*  gd_view(opt)   Print view and region.                                   */

void gd_view(int opt)
{
    printf1("View: %g,%g",GDLon,GDLat);
    if (opt == 0)
        printf1(". Region: %g,%g",GDLonD,GDLatD);
    newline();
    newline();
}

/* ------------------------------------------------------------------------ */
/*  gd_check_inverse(). Return 1 if rectangle (GDPROJ1) is completely in    */
/*                      projection of globe, otherwise return 0.            */

int gd_check_inverse(void)
{
    double x,y;

    x = GDPXA;
    y = GDPYA;
    if (gd_proj_inv(&x,&y) == 0)
        return(0);

    x = GDPXA;
    y = GDPYB;
    if (gd_proj_inv(&x,&y) == 0)
        return(0);

    x = GDPXB;
    y = GDPYA;
    if (gd_proj_inv(&x,&y) == 0)
        return(0);

    x = GDPXB;
    y = GDPYB;
    if (gd_proj_inv(&x,&y) == 0)
        return(0);

    return(1);
}

/* ------------------------------------------------------------------------ */
/*  gd_adlon(lon)   Return adjusted longitude.                              */

double gd_adlon(double lon)
{
    if (lon > 180.0)
        return(lon - 360.0);
    else if (lon < -180.0)
        return(lon + 360.0);
    return(lon);
}

/* ------------------------------------------------------------------------ */
/*  gd_adlat(lat)   Return adjusted latitude.                               */

double gd_adlat(double lat)
{
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

int gd_proj_p(double *x,double *y)
{
    double lon,lat,r,z;

    lon = *x;
    lat = *y;
    r = 1.0;
    ps_3dgeo(&lon,&lat,&r);
    ps_3dprj1(lon,lat,r,x,y,&z);
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

int gd_proj_inv(double *x,double *y)
{
    double lon,lat,r,z;

    z = 1.0 - *x * *x - *y * *y;
    if (z <= 0.0) {
        if (z < -EPSI1)   
            return(0);
        z = 0.0;
    }
    else
        z = sqrt(z);

    ps_3dprj_inv(*x,*y,z,&lon,&lat,&r);
    ps_3dgeo_inv(&lon,&lat,&r);
    *x = lon;
    *y = lat;
    return(1);
}

/* ------------------------------------------------------------------------ */
/*  check_x_in_region(x)                                                    */
/*                                                                          */
/*  Given a longitude x, return 1 if x is in the selected region,           */ 
/*  otherwise return 0.                                                     */

int check_x_in_region(double x)
{
    if (GDLonMin < GDLonMax - EPSI1) {
        if (x >= GDLonMin && x <= GDLonMax)
            return(1);
        return(0);
    }
    if (x >= GDLonMin || x <= GDLonMax)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_y_in_region(y)                                                    */
/*                                                                          */
/*  Given a latitude y, return 1 if y is in the selected region,            */ 
/*  otherwise return 0.                                                     */

int check_y_in_region(double y)
{
    if (y >= GDLatMin && y <= GDLatMax)
        return(1);
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  check_inside(lon,lat)   Return 1 if (lon,lat) is in selected region.    */

int check_inside(double lon,double lat)
{
    if (check_x_in_region(lon) && check_y_in_region(lat))
        return(1);
    return(0);
} 

/* ------------------------------------------------------------------------ */
/*  check_clip(x,y) Assume x and y are user coordinates. Return 0 if (x,y)  */
/*                  is in current plot region, otherwise return 1.          */ 

int check_clip(double x,double y)
{
    if (x >= GDPXA && x <= GDPXB && y >= GDPYA && y <= GDPYB)
        return(0);
    return(1);
} 

/* ------------------------------------------------------------------------ */
/*  gd_proj_gety(lat)                                                       */
/*                                                                          */ 
/*  Given a latitude lon, this function returns the projected coordinate    */
/*  depending on GDPROJ.                                                    */

double gd_proj_gety(double lat)
{
    double y;
    double m12 = 18.77103004;   /* log(tan(45 + 89/2)) */

    switch (GDPROJ) {

        case GDPROJ10:
            y = lat;         
            break;

        case GDPROJ11:
            y = sin(degree_to_arc(lat)) * 90.0;
            break;

        case GDPROJ12:
            if (lat > 89.0)
                lat = 89.0;
            else if (lat < -89.0)
                lat = -89.0;

            y = degree_to_arc(45.0 + lat / 2.0);              
            y = log(tan(y));
            break;
/************************
        case GDPROJ20:

            r = gd_proj_p(x,y);
            if (r)
                return(-1);

            if (GDPROJA == 1 || GDPROJA == 2) {
                *x *= GDRHem;
                *y *= GDRHem;
            }
**************/








            break;

        default:
            printf1("gd_proj_gety: unknown projection %d\n",GDPROJ);
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

int gd_proj_getxy(double *x,double *y)
{
    int r;
    double lon,lat;

    lon = *x;
    lat = *y;

    switch (GDPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            *x = gd_adjust_lon(lon);
            *y = gd_proj_gety(lat);
            r = check_clip(*x,*y);
            return(r);

        case GDPROJ20:
            if (GDPROJA == 1) {
                if (90.0 - *y >= GDRHem)
                    return(-1);
                r = gd_proj_p(x,y);
                if (r)
                    return(-1);
                return(0);
            }
            else if (GDPROJA == 2) {
                if (90.0 + *y >= GDRHem)
                    return(-1);
                r = gd_proj_p(x,y);
                if (r)
                    return(-1);
                return(0);
            }
            else {
                printf("lon=%g lat=%g\n",lon,lat);


                if (check_inside(lon,lat) == 0) {
                    printf("not inside\n"); 
                    return(1);
                }
                r = gd_proj_p(x,y);
                if (r) {
                    printf("not visible\n");
                    return(-1);
                }

                /* check whether the projected point is inside clip region */

                printf("projected: %g %g\n",*x,*y);
   
                if (check_clip(*x,*y) != 0) {
                    printf("not in clip region\n");
                    return(1);
                }
                printf("OK\n");


                return(0);
            }

        default:
            printf1("Error: unknown projection %d\n",GDPROJ);
            exit(0);
    }
    return(0);
}

/* ------------------------------------------------------------------------ */
/*  gd_adjust_lon(lon)                                                      */
/*                                                                          */
/*  Given longitude lon, return rotated value that is nearest to the        */
/*  current clip area.                                                      */
                      
double gd_adjust_lon(double lon)
{
    if (lon < GDPXA && lon + 360.0 <= GDPXB)
        return(lon + 360.0);
    else if (lon > GDPXB && lon - 360.0 >= GDPXA)
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

int sdpgrat(void)  
{
   	int	err;              

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Graticule for geographical coordinate system. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 7,3,0)) {     /* get parameters */
        goto SDPGRATFin;
    }   
    if (PSPROJ == 0) {
        printf1("Error: need a geographical coordinate system.\n");
        return(-1);
    }
    if (PMCONT) {
        printf1("Drawing bounds of selected region.\n");
        mapgrid_bounds();
    }
    if (PMNTP > 0) {
        printf1("Drawing grid lines (meridians).\n");
        if (PMSCFlg == 0)
            PMSC = 10.0;
        else {
            if (PMSC < 0.0)
                PMSC = 0.0;
            else if (PMSC > 60.0)
                PMSC = 60.0;
        }
        mapgrid_lon(PMSC);
    }
    if (PMNTP1 > 0) {
        printf1("Drawing grid lines (parallels).\n");
        mapgrid_lat();
    }
    err = 0;

SDPGRATFin:
    p_clean();
    return(err);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_bounds    Draw bounds of selected region.                       */

void mapgrid_bounds(void)
{
    double r;

    fprintf(PSFd,"\n%%#%d: plgeogrid (bounds)\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);
    ps_lwidth(PMLW);

    switch (GDPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            mapgrid_bounds_r();     
            break;

        case GDPROJ20:
            if (GDPROJA == 0)  
                mapgrid_bounds_r();     
            else {
                r = degree_to_arc(90.0 - GDRHem);
                r = cos(r);
                mapgrid_arc(GDPX,GDPY,0.0,360.0,r,1);
            }
            break;

        default: printf1("Unknown projection %d.\n",GDPROJ);
            return;
    }
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_bounds_r    draw rectangular bounds.                            */

void mapgrid_bounds_r(void)
{
    ps_2dplot(GDPXA,GDPYA,0);
    ps_2dplot(GDPXA,GDPYB,1);
    ps_2dplot(GDPXB,GDPYB,1);
    ps_2dplot(GDPXB,GDPYA,1);
    ps_2dplot(GDPXA,GDPYA,1);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_arc(x,y,a,b,r,gflag)                                            */
/*                                                                          */
/*  Draw arc from a to b with radius r. Center is given by (x,y) in         */
/*  user coordinates for the projection.                                    */
/*  If rflag != 0 fill according to PMGSFlg.                                */

void mapgrid_arc(double x,double y,double a,double b,double r,int gflag)
{
    double rp;

    /* translate to PostScript coordinates */

    rp = (ps_2dx(x + r) - ps_2dx(x - r)) / 2.0;
    x = ps_2dx(x);  
    y = ps_2dy(y);

    if (a < 0.0)
        a = 0.0;
    if (b > 360.0)
        b = 360.0;
            
    fprintf(PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f arc\n",x,y,rp,a,b);

    if (gflag != 0 && PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        if (a > 0.0 || b < 360.0)  
            fprintf(PSFd,"%5.2f %5.2f l\nclosepath\n",x,y);
        fprintf(PSFd,"gsave\n");
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon(sc)   Draw grid lines: meridians. If PMFSX > 0 also plot    */
/*                    labels. Depends on type of projection. sc is degrees  */
/*                    around pole for azimuthal projections.                */

void mapgrid_lon(double sc)
{
    int i;
    double x,xa,ya,yb;

    if (PMLW <= 0.0)
        return;

    fprintf(PSFd,"\n%%#%d: plgeogrid (meridians)\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT1);
    ps_lwidth(PMLW1);

    if (GDPROJ == GDPROJ10 || GDPROJ == GDPROJ11 || GDPROJ == GDPROJ12) {

        for (i = 0; i < PMNTP; ++i) {
            x = PMTP[i];                          
            if (x < -180.0 || x > 180.0)
                continue;
     
            if (check_x_in_region(x)) {
                xa = x;
                if (xa < GDPXA) {
                    if (fabs(xa + 180.0) <= EPSI1)
                        continue;
                    xa += 360.0;
                }
                if (PMLW > 0.0)
                    mapgrid_line(xa,GDPYA,xa,GDPYB);        
                if (PMFSX > 0.0) {
                    mapgrid_label(0,x,xa,GDPYA,PMFSX,2);
                    if (fabs(xa - GDPXA) <= EPSI1 && fabs(GDPXB - GDPXA - 360.0) <= EPSI1 
                        && fabs(xa + 180.0) > EPSI1)
                        mapgrid_label(0,x,GDPXB,GDPYA,PMFSX,2);
                }
            }
        }
    }
    if (GDPROJ == GDPROJ20 && GDPROJA != 0) {
                  
        for (i = 0; i < PMNTP; ++i) {
            x = PMTP[i];                          
            if (x < -180.0 || x > 180.0)
                continue;

            if (GDPROJA == 1 || GDPROJA == 2) {
                if (PMLW > 0.0)
                    mapgrid_lon_azi(x);
                if (PMFSX > 0.0)
                    mapgrid_lon_azi_lab1(x,x,PMFSX);
            }
            else if (GDPROJA == 3) {
                if (PMLW > 0.0)
                    mapgrid_lon_azi1(x,sc);
            }
        }
    }
    if (GDPROJ == GDPROJ20 && GDPROJA == 0) {

        if (PMLW > 0.0) {
            set_clip();     

            for (i = 0; i < PMNTP; ++i) {
                x = PMTP[i];                          
                if (x < -180.0 || x > 180.0)
                    continue;

                if (check_x_in_region(x)) {
                    ya = dmax(GDLatMin,sc - 90.0);
                    yb = dmin(GDLatMax,90.0 - sc);
                    if (ya < yb)
                        map_geodesic_a(x,ya,x,yb);
                }
            }
        }
        if (PMFSX > 0.0) {
            fprintf(PSFd,"grestore\n");
            fprintf(PSFd,"gsave\n");

            for (i = 0; i < PMNTP; ++i) {
                x = PMTP[i];                          
                if (x < -180.0 || x > 180.0)
                    continue;

                if (check_x_in_region(x))  
                    mapgrid_lon_azi_lab0(x,x,PMFSX);          
            }
        }
    }
    fprintf(PSFd,"grestore\n");
}
                  
/* -##--------------------------------------------------------------------- */
/*  mapgrid_line(xa,ya,xb,yb)                                               */
/*                                                                          */
/*  Draw 2d line from (xa,ya) to (xb,yb). Assume that values are given in   */
/*  user coordinates.                                                       */

void mapgrid_line(double xa,double ya,double xb,double yb)
{
    fprintf(PSFd,"gsave\n");
    ps_2dplot(xa,ya,0);
    ps_2dplot(xb,yb,1);
    fprintf(PSFd,"stroke\ngrestore\n");
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_label(opt,lab,x,y,fs,m)                                         */
/*                                                                          */
/*  Plot labels. If opt = 0 then for the X axis, otherwise for the Y axis.  */
/*  Assume that x and y are already given in user coordinates.              */

void mapgrid_label(int opt,double lab,double x,double y,double fs,int m)
{
    int l;
    double px,py,size;
    char buf[30];

    sprintf(buf,PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    px = (x - GDPXA) / UXLen;
    py = (y - GDPYA) / UYLen;
    px *= PSXLen;
    py *= PSYLen;

    size = 1.5 * fs * PtMM;           
    l = strlen(buf);

    if (opt == 0) {
        py -= 1.5 * size;
        plot_str(px,py,buf,size,0,1,0,1,0);
    }
    else {
        px += m * size;                     
        py -= 0.3 * size;
        plot_str(px,py,buf,size,0,2,0,1,0);
    } 
    upd_bbox(1,px + (double)l * size,py);     
    upd_bbox(1,px - (double)l * size,py);     
    upd_bbox(1,px,py + 1.3 * fs * PtMM);
    upd_bbox(1,px,py - 0.5 * fs * PtMM);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lon_azi(lon)                                                    */
/*                                                                          */
/*  Draw meridian as a straight line, beginning at projection center.       */
/*  Called for GDPROJ20 and GDPROJA = 1 or 2.                               */

void mapgrid_lon_azi(double lon)
{
    double a,x,y;

    a = degree_to_arc(lon) - Pi / 2.0;
    x = cos(a) * GDRHemP;
    y = sin(a) * GDRHemP;

    if (GDPROJA == 2)  
        y = -y;
        
    ps_2dplot(x * 0.05,y * 0.05,0);
    ps_2dplot(x,y,1);
    fprintf(PSFd,"stroke\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon_azi_lab0(lab,lon,fs)                                        */
/*                                                                          */
/*  Plot label lab at longitude lon assuming GDPROJA = 0.                   */

void mapgrid_lon_azi_lab0(double lab,double lon,double fs)
{
    int l;
    double x,y,px,py,size;
    char buf[30];

    x = lon;
    if (GDLat >= 0.0)
        y = GDLatMin;
    else
        y = GDLatMax;

    if (gd_proj_p(&x,&y) != 0 || x < GDPXA || x > GDPXB)  
        return;

    if (fabs(lab + 180.0) < EPSI1)  /* don't use -180.0 */
        return;

    sprintf(buf,PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    size = 1.5 * fs * PtMM;           
    px = ps_2dx(x);

    if (GDLat >= 0.0)  
        py = ps_2dy(GDPYA) - 1.5 * size;
    else
        py = ps_2dy(GDPYB) + size;
                         
    l = strlen(buf);
    plot_str(px,py,buf,size,0,1,0,1,0);

    upd_bbox(1,px + (double)l * size,py);     
    upd_bbox(1,px - (double)l * size,py);     
    upd_bbox(1,px,py + 1.3 * fs * PtMM);
    upd_bbox(1,px,py - 0.5 * fs * PtMM);
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lon_azi_lab1(lab,lon,fs)                                        */
/*                                                                          */
/*  Plot label lab at longitude lon assuming GDPROJA = 1 or 2.              */

void mapgrid_lon_azi_lab1(double lab,double lon,double fs)
{
    int l,cflag;
    double a,x,y,px,py,size;
    char buf[30];

    if (fabs(lab + 180.0) < EPSI1)  /* don't use -180.0 */
        return;

    sprintf(buf,PMFmtS,lab);
    *(buf + strlen(buf) - 1) = '\0';    /* drop last blank */

    a = degree_to_arc(lon) - Pi / 2.0;
    x = cos(a) * GDRHemP * 1.05;
    y = sin(a) * GDRHemP * 1.05;

    if (GDPROJA == 2)  
        y = -y;
      
    px = ps_2dx(x);
    py = ps_2dy(y);

    size = fs * PtMM;           
    py -= size / 2.0;
    size *= 1.5;

    l = strlen(buf);

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

    plot_str(px,py,buf,size,0,cflag,0,1,0);

    upd_bbox(1,px + (double)l * size,py);     
    upd_bbox(1,px - (double)l * size,py);     
    upd_bbox(1,px,py + 1.3 * fs * PtMM);
    upd_bbox(1,px,py - 0.5 * fs * PtMM);
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lon_azi1(lat,sc)                                                */
/*                                                                          */
/*  Draw meridians. Called for GDPROJ20, GDPROJA = 3. sc is degrees around  */
/*  pole.                                                                   */

void mapgrid_lon_azi1(double lon,double sc)
{
    double a,b,dlon;

    dlon = gd_adlon(lon + 180.0);
        
    if (GDLat >= 0.0) {
        b = 90.0 - sc;
        a = dmin(b,GDLat);
        if (a < b) {
            map_geodesic_a(lon,a,lon,b);
            map_geodesic_a(dlon,a,dlon,b);
        }
        b = sc - 90.0;
        if (a > b)   
            map_geodesic_a(lon,a,lon,b);
    }
    else {
        b = sc - 90.0;
        a = dmax(b,GDLat);
        if (a > b) {
            map_geodesic_a(lon,a,lon,b);
            map_geodesic_a(dlon,a,dlon,b);
        }
        b = 90.0 - sc;
        if (a < b)   
            map_geodesic_a(lon,a,lon,b);
    }
}

/* -##--------------------------------------------------------------------- */
/*  mapgrid_lat    Draw grid lines: parallels. If PMFSY > 0 also plot       */
/*                 labels. Depends on type of projection.                   */

void mapgrid_lat(void)
{
    int i;
    double y,yp;

printf("mapgrid_lat: GDPROJ=%d GDPROJA=%d\n",GDPROJ,GDPROJA);

    if (PMLW <= 0.0)
        return;

    fprintf(PSFd,"\n%%#%d: plgeogrid (parallels)\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT1);
    ps_lwidth(PMLW1);

    if (GDPROJ == GDPROJ20 && GDPROJA == 0) 
        set_clip();     
          
    for (i = 0; i < PMNTP1; ++i) {

        y = PMTP1[i];                          
        if (y < -90.0 || y > 90.0)
            continue;
     
        switch (GDPROJ) {

            case GDPROJ10:
            case GDPROJ11:
            case GDPROJ12:
                if (check_y_in_region(y)) {
                    yp = gd_proj_gety(y);
                    if (PMLW > 0.0)
                        mapgrid_line(GDPXA,yp,GDPXB,yp);
                    if (PMFSY > 0.0)  
                        mapgrid_label(1,y,GDPXB,yp,PMFSY,2);
                }
                break;  

            case GDPROJ20:
                if (GDPROJA == 1) {
                    if (y != 90.0 && 90.0 - y <= GDRHem) {
                        y = degree_to_arc(y);
                        y = cos(y);
                        mapgrid_arc(GDPX,GDPY,0.0,360.0,y,0);
                    }
                }
                else if (GDPROJA == 2) {
                    if (y != -90.0 && 90.0 + y <= GDRHem) {        
                        y = degree_to_arc(y);
                        y = cos(y);
                        mapgrid_arc(GDPX,GDPY,0.0,360.0,y,0);
                    }
                }
                else if (GDPROJA == 3) {    
                    mapgrid_lat_a3(y);                     
                }
                else if (GDPROJA == 0) {            /* rectangular region */
                    mapgrid_lat_a0(y);          
                }
                break; 

            default:
                printf1("Unknown projection %d\n",GDPROJ);
                break;  
        }
    }
    fprintf(PSFd,"grestore\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lat_a0(lat)                                                     */
/*                                                                          */
/*  Draw parallel. Called for GDPROJ20, GDPROJA = 0 (rectangular region).   */

void mapgrid_lat_a0(double lat)
{
    double x,y,xc,yc,r,a,b,c,beta,sinb,cosb,tlat,r1,r2,s,t,w,phi;
                                
    if (GDLat >= 90.0 - EPSI1 || GDLat <= EPSI1 - 90.0)
        return;

    if (lat <= GDLatMin || lat >= GDLatMax)
        return;
    
    beta = degree_to_arc(lat);
    sinb = sin(beta);
    cosb = cos(beta);

    if (fabs(GDLat) <= EPSI1) {     /* straight lines */ 
        fprintf(PSFd,"gsave\n");
        x = GDPXA;
        y = GDPY + sinb;
        ps_2dplot(x,y,0);
        x = GDPXB;
        ps_2dplot(x,y,1);
        fprintf(PSFd,"stroke\n");
        fprintf(PSFd,"grestore\n");
        return;
    }
    ps_3dprj(0.0,0.0,sinb,&xc,&yc);              /* xc should be zero */

    s = cos(degree_to_arc(90.0 - fabs(GDLat)));
    tlat = tan(degree_to_arc(GDLat));
    t = -tlat * sinb / cosb;

    if (t <= -1.0 + EPSI1 || t >= 1.0 - EPSI1) {
        r1 = 0.0;
        r2 = 360.0;
    }
    else {  
        phi = acos(t) * 180.0 / Pi;
        if (GDLat >= 0.0) {
            r1 = 270.0 - phi;
            r2 = phi - 90.0;
        }
        else {
            r1 = 90.0 - phi;
            r2 = 90.0 + phi;
        }
    }
    xc = ps_2dx(xc);
    yc = ps_2dy(yc);
    r = cosb * (ps_2dx(1.0) - ps_2dx(-1.0)) / 2.0;

    fprintf(PSFd,"gsave\n");
    fprintf(PSFd,"%7.3f %7.3f translate\n",xc,yc);
    fprintf(PSFd,"1.0 %7.4f scale\n",s);   
    fprintf(PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",r,r1,r2);
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
}

/* -###-------------------------------------------------------------------- */
/*  mapgrid_lat_a3(lat)                                                     */
/*                                                                          */
/*  Draw parallel. Called for GDPROJ20, GDPROJA = 3.                        */

void mapgrid_lat_a3(double lat)
{
    double x,y,xc,yc,r,a,b,c,beta,sinb,cosb,tlat,r1,r2,s,t,w,phi;
                                
    if (GDLat >= 90.0 - EPSI1 || GDLat <= EPSI1 - 90.0)
        return;

    if (lat <= GDLat - 90.0 + EPSI1 ||
        lat >= GDLat + 90.0 - EPSI1) 
        return;
    
    beta = degree_to_arc(lat);
    sinb = sin(beta);
    cosb = cos(beta);

    if (fabs(GDLat) <= EPSI1) {     /* straight lines */ 
        fprintf(PSFd,"gsave\n");
        x = GDPX - cosb * (GDPX - GDPXA);
        y = GDPY + sinb * (GDPY - GDPYA);
        ps_2dplot(x,y,0);
        x = GDPX + cosb * (GDPX - GDPXA);
        ps_2dplot(x,y,1);
        fprintf(PSFd,"stroke\n");
        fprintf(PSFd,"grestore\n");
        return;
    }
    ps_3dprj(0.0,0.0,sinb,&xc,&yc);              /* xc should be zero */
    s = cos(degree_to_arc(90.0 - fabs(GDLat)));
    tlat = tan(degree_to_arc(GDLat));
    t = -tlat * sinb / cosb;

    if (t <= -1.0 + EPSI1 || t >= 1.0 - EPSI1) {
        r1 = 0.0;
        r2 = 360.0;
    }
    else {  
        phi = acos(t) * 180.0 / Pi;
        if (GDLat >= 0.0) {
            r1 = 270.0 - phi;
            r2 = phi - 90.0;
        }
        else {
            r1 = 90.0 - phi;
            r2 = 90.0 + phi;
        }
    }
    xc = ps_2dx(xc);
    yc = ps_2dy(yc);
    r = cosb * (ps_2dx(GDPXB) - ps_2dx(GDPXA)) / 2.0;

    fprintf(PSFd,"gsave\n");
    fprintf(PSFd,"%7.3f %7.3f translate\n",xc,yc);
    fprintf(PSFd,"1.0 %7.4f scale\n",s);   
    fprintf(PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",r,r1,r2);
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
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

int sdpgeo(void)  
{
    register int i,j;
   	int	err,n;
    double x,y;

    err = -1;
    if (check_cmd(2))
        return(-1);
                    
    printf1("Plot in geographical coordinate system. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,7,1)) {     /* get parameters */
        goto SDPGEOFin;
    }   
    if (PSPROJ == 0) {
        printf1("Need a geographical coordinate system.\n");
        return(-1);
    }
    if (PMLWFlg == 0)
        PMLW = 0.05;

    if (PMCT != 1)
        PMCT = 0;

    if (PMS < 1 || PMS > 17)        /* check for valid marker symbol */
        PMS = 0;           

    n = PMRHSN / 2;
    if (2 * n != PMRHSN) {
        printf1("Error: right-hand side must define pairs of coordinates.\n");
        goto SDPGEOFin;
    }
    if (n == 1 && PMS == 0) {
        printf1("Single points require marker symbol.\n");
        err = 0;
        goto SDPGEOFin;
    }
    if (alloc_acx(n + 1))
        return(-1);              
    if (alloc_acy(n + 1))
        return(-1);                  

    j = 0;
    for (i = 0; i < n; ++i) {
        x = PMRHSX[j++];
        y = PMRHSX[j++];
        if (x < -180.0 || x > 180.0 || y < -90.0 || y > 90.0) {
            printf1("Error: right-hand side violates range of geographical coordinates.\n");
            goto SDPGEOFin;
        }
        AcX[i] = x;
        AcY[i] = y;
    }
    fprintf(PSFd,"\n%%#%d: plgeod\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);
    if (PMS != 0)
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);

    /* set_clip(); not used in this function */

    if (PMLW > 0.0)
        map_line(n,AcX,AcY,PMCT);

    if (PMS != 0) {
        for (i = 0; i < n; ++i)  
            map_point(AcX[i],AcY[i],0);
    }
    fprintf(PSFd,"grestore\n");
    printf1("PostScript output written to: %s\n",PSFName);
    err = 0;

SDPGEOFin:
    p_clean();
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

int sdpmap(void)  
{
    register int i,j;
   	int	err,n1,n2,n3,n1p,n2p,n3p,n,cflag,typ,r,sdid;
    double x,y,d;

    err = -1;
    if (check_cmd(0))
        return(-1);
                    
    printf1("Plot in geographical coordinate system. Current memory: %d bytes.\n",MemReq);

    if (parm(CmdBuf + 6,3,0)) {     /* get parameters */
        goto SDPMAPFin;
    }   
    if (check_sd(0)) {
        printf1("No spatial data defined.\n");
        return(-1);
    }
    if (PSPROJ == 0) {
        printf1("Need a geographical coordinate system.\n");
        return(-1);
    }
    if (PMLWFlg == 0)
        PMLW = 0.05;

    if (PMS < 1 || PMS > 17)        /* check for valid marker symbol */
        PMS = 0;           

    fprintf(PSFd,"\n%%#%d: plgeosd\n",++PSONUM);
    fprintf(PSFd,"gsave\n");
    ps_ltyp(PMLT);  
    ps_lwidth(PMLW);
    if (PMS != 0)
        fprintf(PSFd,"/ssiz %4.2f def\n",PtMM * PMFS / 2.0);
    /* set_clip(); */

    n1p = n1 = 0;     /* number of points */
    n2p = n2 = 0;     /* number of lines */
    n3p = n3 = 0;     /* number of polygons */

    for (i = 0; i < NOC; ++i) {

        if ((typ  = (int)get_data(SDVarSDTyp,i)) < 1 || typ > 3)
            continue;

        sdid = (int)get_data(SDVarSDID,i);                      

        if (typ == 3)  
            cflag = 1;
        else
            cflag = 0;

        if ((n = sd_getdata(i,0,cflag,1)) < 1)
            goto SDPMAPFin;
             
        if (typ == 1) {
            n1++;
            if (PMS) {
                if (map_point(SDVarX[0],SDVarY[0],0))
                    n1p++;
            }
        }
        else if (typ == 2) {
            n2++;
            r = map_line(n,SDVarX,SDVarY,0);
            if (r < 0)
                goto SDPMAPFin;
            n2p += r;

            if (PMS) {
                for (j = 0; j < n; ++j)  
                    map_point(SDVarX[j],SDVarY[j],0);
            }
        }
        else if (typ == 3) {
            n3++;
            r = map_polygon(n - 1,SDVarX,SDVarY,sdid);   
            if (r < 0)  
                r = 0;
            else
                n3p += r;

            if (PMS) {
                if ((n = sd_getdata(i,0,1,1)) < 1)
                    goto SDPMAPFin;
                for (j = 1; j < n; ++j)  
                    map_point(SDVarX[j],SDVarY[j],0);
            }
        }
    }
    fprintf(PSFd,"grestore\n");

    printf1("\nNumber of type 1 objects: %d",n1);
    if (n1p > 0)
        printf1(". Plotted: %d",n1p);
    printf1("\nNumber of type 2 objects: %d",n2);
    if (n2p > 0)
        printf1(". Plotted: %d (parts)",n2p);
    printf1("\nNumber of type 3 objects: %d",n3);
    if (n3p > 0)
        printf1(". Plotted: %d (parts)",n3p);
    newline();
    printf1("PostScript output written to: %s\n",PSFName);
    err = 0;

SDPMAPFin:
    p_clean();
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
    
int map_point(double x,double y,int nc)
{
    int r = 0;

    if (PMS == 0)
        return(0);

    switch (PSPROJ) {

        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            r = map_point_cyl(x,y);  
            break;

        case GDPROJ20:
            r = map_point_azi(x,y,nc);
            break;

        default:
            printf("projection not supported.\n");
            break;  
    }
    return(r);
}

/* ------------------------------------------------------------------------ */
/*  map_point_cyl(x,y)  Called by map_point for cylindrical projections.    */

int map_point_cyl(double x,double y)
{
    double d;

    if (gd_proj_getxy(&x,&y) == 0) {    
        x = ps_2dx(x);
        y = ps_2dy(y);
        d = PtMM * PMFS / 2.0;
        ps_sym(PMS,x,y,d);
        return(1);
    }
    return(0);
}

/* -###-------------------------------------------------------------------- */
/*  map_point_azi(x,y,nc)  Called by map_point for azimuthal projections.   */

int map_point_azi(double x,double y,int nc)
{
    int r;
    double d;
            
printf("map_point_azi lon=%g %g\n",x,y);

    if (nc == 0) {
        if (gd_proj_getxy(&x,&y) == 0) {        
            x = ps_2dx(x);
            y = ps_2dy(y);
            d = PtMM * PMFS / 2.0;
            ps_sym(PMS,x,y,d);
            return(1);
        }
    }
    else {
        if (gd_proj_p(&x,&y)) {
            printf(" point  not visible\n");
            return(0);
        }
        upd_bbox(0,x,y);
        x = ps_2dx(x);
        y = ps_2dy(y);
        d = PtMM * PMFS / 2.0;
        ps_sym(PMS,x,y,d);
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
    
int map_line(int n,double *x,double *y,int gflag)
{
    int r = 0;

    switch (PSPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            if (gflag)
                r = map_line_cyl_g(n,x,y);
            else
                r = map_line_cyl_s(n,x,y);  
            break;
   
        case GDPROJ20:
            r = map_line_azi(n,x,y,gflag);  
            break;

        default:
            printf("projection not supported.\n");
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

int map_line_cyl_s(int n,double *x,double *y)
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
        xa = gd_adjust_lon(xa);
        xb = gd_adjust_lon(xb);
    
        map_find_n(xa,ya,&xb,&yb);
        r1 = clip_line(xa,ya,xb,yb,GDPXA,GDLatA,GDPXB,GDLatB,
             &sax1,&say1,&sbx1,&sby1,&ca,&cb,&ra,&rb,&dir);

        xa = x[i - 1];
        ya = y[i - 1];
        xb = x[i];
        yb = y[i];
        xa = gd_adjust_lon(xa);
        xb = gd_adjust_lon(xb);

        map_find_n(xb,yb,&xa,&ya);
        r2 = clip_line(xa,ya,xb,yb,GDPXA,GDLatA,GDPXB,GDLatB,
             &sax2,&say2,&sbx2,&sby2,&ca,&cb,&ra,&rb,&dir);
    
        if (r1) {
            yp = gd_proj_gety(say1);
            ps_2dplot(sax1,yp,0);
            yp = gd_proj_gety(sby1);
            ps_2dplot(sbx1,yp,1);
            fprintf(PSFd,"stroke\n");    
            rt++;    

            if (r2) {
                if (fabs(sax1 - sax2) > EPSI1 || fabs(say1 - say2) > EPSI1 ||
                    fabs(sbx1 - sbx2) > EPSI1 || fabs(sby1 - sby2) > EPSI1) {  
                    yp = gd_proj_gety(say2);
                    ps_2dplot(sax2,yp,0);
                    yp = gd_proj_gety(sby2);
                    ps_2dplot(sbx2,yp,1);
                    fprintf(PSFd,"stroke\n");    
                    rt++;
                }
            }
        }
        else if (r2) {
            yp = gd_proj_gety(say2);
            ps_2dplot(sax2,yp,0);
            yp = gd_proj_gety(sby2);
            ps_2dplot(sbx2,yp,1);
            fprintf(PSFd,"stroke\n");    
            rt++;     
        }
    }
    return(rt);
}

/* ------------------------------------------------------------------------ */
/*  map_find_n(x,y,xa,ya)       Find nearest point.                         */

void map_find_n(double x,double y,double *xa,double *ya)
{
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

int map_line_cyl_g(int n,double *x,double *y)
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

        if (xa < GDPXA)
            xa += 360.0;
        else if (xa > GDPXB)
            xa -= 360.0;

        map_find_n(xa,ya,&xb,&yb); 
        rt = map_line_cyl_geod(xa,ya,xb,yb);
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

int map_line_cyl_geod(double lon1,double lat1,double lon2,double lat2)
{
    int i,j,first,rt,ca,cb,ra,rb,dir;
    double d,sind,slon1,slon2,slat1,slat2,clon1,clon2,clat1,clat2;
    double f,x,y,z,a,b,t,rm,x0,xl,yl,delta,sax,say,sbx,sby;
    double xp[NPNT_GEOD + 1],yp[NPNT_GEOD + 1];

    /* check for antipodal points */

    if (fabs(lat1 + lat2) <= EPSI1 && fabs(fabs(lon1 - lon2) - 180.0) <= EPSI1)
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
    d = geod_distance(lon1,lat1,lon2,lat2);
    if (d <= EPSI1)
        return(0);

    if (fabs(d - Pi) <= EPSI1)  
        return(0);            

    sind = sin(d);

    lon1 = degree_to_arc(lon1);
    lat1 = degree_to_arc(lat1);
    lon2 = degree_to_arc(lon2);
    lat2 = degree_to_arc(lat2);

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
        if (gd_proj_getxy(&x,&y) != 0) {
            if (first) {
                xl = xp[i - 1];
                yl = yp[i - 1];
                gd_proj_getxy(&xl,&yl);             

                if (clip_line(xl,yl,x,y,GDPXA,GDPYA,GDPXB,GDPYB,
                    &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) && 
                    fabs(xl - sax) <= EPSI1 && fabs(yl - say) <= EPSI1) {
                    ps_2dplot(sbx,sby,first);
                }
/*              map_line_clip(xl,yl,x,y);           */
                fprintf(PSFd,"stroke\n");    
                rt = 1;
                first = 0;
            }
            continue;
        }
        if (j < 0)  
            x0 = x;
        else if (x < x0 - EPSI1)
            break;

        j = i;

        if (first == 0 && i > 0) {
            xl = xp[i - 1];
            yl = yp[i - 1];
            gd_proj_getxy(&xl,&yl);             

            if (clip_line(xl,yl,x,y,GDPXA,GDPYA,GDPXB,GDPYB,
                 &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) &&  
                 fabs(sbx - x) <= EPSI1 && fabs(sby - y) <= EPSI1) {
                 ps_2dplot(sax,say,first);
                 first = 1;
            }
        }
        ps_2dplot(x,y,first);
        first = 1;
    }
    if (j < 0)
        return(0);

    if (first) {
        fprintf(PSFd,"stroke\n");    
        first = 0;
        rt = 1;
    }
    if (j != NPNT_GEOD) {
        first = 0;
        for (i = NPNT_GEOD; i >= 0; --i) {
            x = xp[i];
            y = yp[i];
            if (gd_proj_getxy(&x,&y) != 0) {
                if (first) {

                    xl = xp[i + 1];
                    yl = yp[i + 1];
                    gd_proj_getxy(&xl,&yl);             

                    if (clip_line(xl,yl,x,y,GDPXA,GDPYA,GDPXB,GDPYB,
                        &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) && 
                        fabs(xl - sax) <= EPSI1 && fabs(yl - say) <= EPSI1) {
                        ps_2dplot(sbx,sby,first);
                    }
/*                  map_line_clip(xl,yl,x,y);       */
                    fprintf(PSFd,"stroke\n");    
                    rt = 1;
                    first = 0;
                }
                continue;
            }
            if (i <= j || x >= x0 - EPSI1)
                break;

            if (first == 0 && i < NPNT_GEOD) {
                xl = xp[i + 1];
                yl = yp[i + 1];
                gd_proj_getxy(&xl,&yl);             

                if (clip_line(xl,yl,x,y,GDPXA,GDPYA,GDPXB,GDPYB,
                     &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir) &&  
                     fabs(sbx - x) <= EPSI1 && fabs(sby - y) <= EPSI1) {
                     ps_2dplot(sax,say,first);
                     first = 1;
                }
            }
            ps_2dplot(x,y,first);
            first = 1;
        }
    }
    if (first) {
        fprintf(PSFd,"stroke\n");    
        rt = 1;
    }
    return(rt);
}

/* -###-------------------------------------------------------------------- */
/*  map_line_clip(xa,ya,xb,yb)                                              */
/*                                                                          */
/*  (xa,ya) and (xb,yb) are points in the projection plane. also (xa,ya)    */
/*  is contained in the clip region. Connect by a straight line.            */

void map_line_clip(double xa,double ya,double xb,double yb)
{                 
    int ca,cb,ra,rb,dir;
    double sax,say,sbx,sby;

printf("maplineclip %g %g %g %g\n",xa,ya,xb,yb);


    if (clip_line(xa,ya,xb,yb,GDPXA,GDPYA,GDPXB,GDPYB,
             &sax,&say,&sbx,&sby,&ca,&cb,&ra,&rb,&dir)) {
printf("sax %g %g %g %g\n",sax,say,sbx,sby);
        if (fabs(xa - sax) <= EPSI1 && fabs(ya - say) <= EPSI1) {
printf("plotted sbx=%g %g\n",sbx,sby);
            ps_2dplot(sbx,sby,1);
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

int map_line_azi(int n,double *x,double *y,int gflag)
{
    register int i;
    int r,nl = 0;

    printf("map_line_azi gflag=%d  \n",gflag);

    for (i = 1; i < n; ++i) {
        if (gflag == 0)
            r = map_line_straight(x[i - 1],y[i - 1],x[i],y[i]);
        else  
            r = map_geodesic_a(x[i - 1],y[i - 1],x[i],y[i]);

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

int map_line_straight(double lona,double lata,double lonb,double latb)
{
    int ra,rb,ca,cb,dir;
    double ax,ay,bx,by,sax,say,sbx,sby;

    ax = lona;
    ay = lata;
    if (gd_proj_p(&ax,&ay)) {
        printf("first point  not visible\n");
        return(0);
    }
    bx = lonb;
    by = latb;
    if (gd_proj_p(&bx,&by)) {
        printf("second point  not visible\n");
        return(0);
    }
    if (clip_line(ax,ay,bx,by,GDPXA,GDPYA,GDPXB,GDPYB,&sax,&say,&sbx,&sby, 
        &ca,&cb,&ra,&rb,&dir) == 0)  
        return(0);
           
    ps_2dplot(sax,say,0);
    ps_2dplot(sbx,sby,1);
    fprintf(PSFd,"stroke\n");
    return(1);
}

/* -###-------------------------------------------------------------------- */
/*  map_geodesic_a(lona,lata,lonb,latb)                                     */
/*                                                                          */  
/*  Draw geodesic from (lona,lata) to (lonb,latb).                          */
/*  Return 1 if some part is successfully drawn, 0 if nothing is drawn,     */
/*  -1 if error.                                                            */

int map_geodesic_a(double lona,double lata,double lonb,double latb)
{
    int la,lb;
    double x,y,ax,ay,az,bx,by,bz,ra,rb,rd,rm,r1,r2;
    double apx,apy,apz,bpx,bpy,bpz,cpx,cpy,cpz;
    double b,d,rx,ry,rz,cx,cy,cz,nx,ny,nz,cosphi,psi;

    rm = 180.0 / Pi;
      
    printf("\ngeodesic: lona=%g %g lonb=%g %g\n",lona,lata,lonb,latb); 
   
    if (fabs(lona - lonb) <= EPSI1 && fabs(lata - latb) <= EPSI1)
        return(0);
         
    ax = lona;
    ay = lata;
    az = 1.0;
    ps_3dgeo(&ax,&ay,&az);

    /** printf("ax=%g %g %g\n",ax,ay,az); **/ 

    ps_3dprj1(ax,ay,az,&apx,&apy,&apz);           
        printf("ax=%g %g %g apx=%g %g %g\n",ax,ay,az,apx,apy,apz); 

    bx = lonb;
    by = latb;
    bz = 1.0;
    ps_3dgeo(&bx,&by,&bz);
    ps_3dprj1(bx,by,bz,&bpx,&bpy,&bpz);           

    /**  printf("bx=%g %g %g bpx=%g %g %g\n",bx,by,bz,bpx,bpy,bpz); **/

    if (apz < 0.0 && bpz < 0.0)    /* both points not visible */
        return(0);

    rx = GDLon;
    ry = GDLat;
    rz = 1.0;
    ps_3dgeo(&rx,&ry,&rz);
    printf("rx=%g %g %g\n",rx,ry,rz); 
    
    d = rx * bx + ry * by + rz * bz;

           
    if (fabs(d) <= EPSI1) {
        cx = bx;
        cy = by;
        cz = bz;
    }
    else {
        b = (rx * ax + ry * ay + rz * az);
        if (fabs(b - d) <= EPSI1)
            d = 1.0;
        else
            d = b / d;

        cx = ax - d * bx;
        cy = ay - d * by;
        cz = az - d * bz;
        d = sqrt(cx * cx + cy * cy + cz * cz);
        if (fabs(d) <= EPSI1) {
            err_geodesic(-1,lona,lata,lonb,latb);
            return(-1);
        }
        cx /= d;
        cy /= d;
        cz /= d;
    }
    /** printf("cx=%g %g %g\n",cx,cy,cz); **/

    /** g_spat(rx,ry,rz,cx,cy,cz,&dx,&dy,&dz); **/

    /** printf("dx=%g %g %g\n",dx,dy,dz); **/ 

    g_spat(bx,by,bz,ax,ay,az,&nx,&ny,&nz);

    printf("nx=%g %g %g\n",nx,ny,nz);       

    d = sqrt(nx * nx + ny * ny + nz * nz);
    if (fabs(d) <= EPSI1) {
        err_geodesic(-3,lona,lata,lonb,latb);
        return(-3);
    }
    cosphi = (nx * rx + ny * ry + nz * rz) / d;
   
        printf("cosphi=%g\n",cosphi); 

    /* length of minor semi-axis */
                                                              
    b = fabs(cosphi); 

     printf("b=%g\n",b);

    ps_3dprj1(cx,cy,cz,&cpx,&cpy,&cpz);           
    if (fabs(cpz) > EPSI1 || (fabs(cpx) <= EPSI1 && fabs(cpy) <= EPSI1)) {
        err_geodesic(-2,lona,lata,lonb,latb);
        return(-2);
    }
    /** printf("cpx=%g cpy=%g cpz=%g \n",cpx,cpy,cpz); */

    psi = rm * atan2(cpy,cpx);
    /** printf("proj c... x=%g %g psi=%g\n",cpx,cpy,psi); **/

    la = g_left(cpx,cpy,-cpx,-cpy,apx,apy);
    lb = g_left(cpx,cpy,-cpx,-cpy,bpx,bpy);

    printf("la=%d lb=%d\n",la,lb); 

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

    printf("\nra=%g rb=%g\n",ra,rb);

    /* check for visibility */ 

    printf("apz=%g bpz=%g\n",apz,bpz);

    if (apz >= 0.0 && bpz >= 0.0) {         /* both points visible */
        if (la == 0 && lb == 0) {           /* both to the right */
            r1 = dmin(ra,rb);
            r2 = dmax(ra,rb);
        }
        else {                              /* both to the left */
            ra = 360.0 - ra;
            rb = 360.0 - rb;
            r1 = dmin(ra,rb);
            r2 = dmax(ra,rb);
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

    printf("\nr1=%g r2=%g b=%g\n",r1,r2,b);

    if (fabs(r1 - r2) <= EPSI1)  
        return(0);
          
/*  rd = (ps_2dx(GDPXB) - ps_2dx(GDPXA)) / 2.0;     */


    rd = (ps_2dx(1.0) - ps_2dx(-1.0)) / 2.0;
    x = ps_2dx(GDPX);   
    y = ps_2dy(GDPY);

    fprintf(PSFd,"gsave\n");
    fprintf(PSFd,"%7.3f %7.3f translate\n",x,y);
    fprintf(PSFd,"%7.4f rotate\n",psi);
    fprintf(PSFd,"1.0 %7.4f scale\n",b);   
    fprintf(PSFd,"0.0 0.0 %7.3f %7.3f %7.3f arc\n",rd,r1,r2);
    fprintf(PSFd,"stroke\n");
    fprintf(PSFd,"grestore\n");
    return(1);
}

/* -###-------------------------------------------------------------------- */
/*  err_geodesic(n,lona,lata,lonb,latb)                                     */

void err_geodesic(int n,double lona,double lata,double lonb,double latb)
{
    printf1("Error (%d) in calculation of geodesic: (%g,%g), (%g,%g)\n",
        n,lona,lata,lonb,latb);
}

/* ------------------------------------------------------------------------ */
/*  map_polygon(n,x,y,sdid)   Plot polygon with geographical coordinates.   */
/*                                                                          */
/*  x[i], y[i] (i=0,n-1) contains the points.                               */
    
int map_polygon(int n,double *x,double *y,int sdid)
{
    int np = 0;

    switch (PSPROJ) {
        case GDPROJ10:
        case GDPROJ11:
        case GDPROJ12:
            np = map_polygon_cyl(n,x,y,sdid);
            break;



        case GDPROJ2:
               map_polygon_p2(n,x,y);
            break;

        default:
            printf("Projection not supported.\n");
            break;  
    }
    return(np);
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_p2(n,x,y)   Called by map_polygon for projection type 2.    */
    
void map_polygon_p2(int n,double *xx,double *yy)
{
    register int i;
    int first,r,vflag,nflag,npi,np;   
    double x,y,xf,yf,xl,yl,xn,yn;

    np = i = 0;
    while (i < n) {         /* find first visible point */ 
        x = xx[i];
        y = yy[i];
        r = gd_proj_p(&x,&y);
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
        r = gd_proj_p(&x,&y);
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
                plsd_arc(xn,yn,x,y);           
            npi = nflag = 0;
        }
        else {
            ps_2dplot(x,y,first);
            first = 1;
            np++;
            npi++;
        }
        xl = x;
        yl = y; 

    }
    if (first) {
        if (nflag)   
            plsd_arc(xn,yn,xf,yf);           

        if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
            fprintf(PSFd,"closepath\ngsave\n%4.2f setgray\n",PMGS);    
            fprintf(PSFd,"fill\ngrestore\n");
        }
        if (PMLW > 0.0)
            fprintf(PSFd,"stroke\n");
    }
}

/* ------------------------------------------------------------------------ */
/*  plsd_arc(x0,y0,x1,y1)   Draw arc from (x0,y0) to (x1,y1). Assume that   */
/*                          points are given in user coordinates.           */
    
void plsd_arc(double x0,double y0,double x1,double y1)
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
    if (d <= EPSI1)
        return;

    fprintf(PSFd,"%5.2f %5.2f %5.2f %5.2f %5.2f ",GDPXPS,GDPYPS,GDPSRadius,a,b);
    if (d <= 180.0) {
        if (a <= b)
            fprintf(PSFd,"arc\n");
        else
            fprintf(PSFd,"arcn\n");
    }
    else {
        if (a <= b)
            fprintf(PSFd,"arcn\n");
        else
            fprintf(PSFd,"arc\n");
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

int map_polygon_cyl(int n,double *x,double *y,int sdid)
{
    register int i,j;
    int rt,r,ns,txflag,pflag,imin,imax;
    double x0,y0,x1,y1,xmin,xmax,ymin,ymax,d;

    if ((pflag = map_check_pole(n,x,y))) {
        printf("POLE r=%d  n=%d sdid=%d\n",r, n,sdid);
                      
    }


/**   
    printf("\nNew polygon\n");
    for (i = 0; i < n; ++i)
        printf("i=%4d x=%f y=%f\n",i,x[i],y[i]);
**/   

    /* transformation into a planar polygon */

    j = -1;         /* find one point inside plot region */

    for (i = 0; i < n; ++i) {
        if (check_inside(x[i],y[i])) {
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
        map_find_n(x0,y0,&x1,&y1);
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
        ymin = dmin(y0,ymin);
        ymax = dmax(y0,ymax);
    }        
    for (i = j - 1; i >= 0; --i) {
        x1 = x[i];
        y1 = y[i];
        map_find_n(x0,y0,&x1,&y1);
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
        ymin = dmin(y0,ymin);
        ymax = dmax(y0,ymax);
    }        

    printf("imin=%d imax=%d\n",imin,imax);
    printf("xmin=%g xmax=%g\n",xmin,xmax);
    printf("ymin=%g ymax=%g\n",ymin,ymax);




    if (pflag != 0) {               /* add points for polar regions */
        if (pflag == 1)
            y0 = 90.0;
        else
            y0 = -90.0;

        if (imin == 0 && imax == n - 1) {
            if (xmax < xmin + 360.0 - EPSI1) {
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
            if (xmin > xmax - 360.0 + EPSI1) {
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
            printf1("Error: cannot plot polygon with ID %d containing a pole (%d).\n",sdid,pflag);
            printf1("Will be skipped.\n");
            return(0);
        }
        ymin = dmin(ymin,y0);
        ymax = dmax(ymax,y0);
    }
/**   
    printf("\nNew polygon nach Umwandlung imin=%d imax=%d\n",imin,imax);
    for (i = 0; i < n; ++i)
        printf("i=%4d x=%f y=%f\n",i,x[i],y[i]);
      
    printf("xmin=%g xmax=%g\n",xmin,xmax);
    printf("ymin=%g ymax=%g\n",ymin,ymax);
**/         
    if (ymin >= GDLatB || ymax <= GDLatA)
        return(0);

    if (xmax <= GDPXA) {
        if (xmin + 360.0 >= GDPXB)
            return(0);

        for (i = 0; i < n; ++i)
            x[i] += 360.0;

        xmin += 360.0;
        xmax += 360.0;
    }
    if (xmin >= GDPXB) {
        if (xmax - 360.0 <= GDPXA)
            return(0);

        for (i = 0; i < n; ++i)
            x[i] -= 360.0;

        xmin -= 360.0;
        xmax -= 360.0;
    }

    txflag = 0;
    if (xmax > GDPXB && xmax - 360.0 > GDPXA) {
        txflag = -1;
    }
    else if (xmin < GDPXA && xmin + 360.0 < GDPXB) {
        txflag = 1;
    }


    printf("\nNew polygon nach Verschiebung\n");
    for (i = 0; i < n; ++i)
        printf("i=%4d x=%f y=%f\n",i,x[i],y[i]);


    printf("TXFLAG=%d\n",txflag);

    rt = 0;
    r = map_polygon_cyl_clip(n,x,y,sdid);
    if (r < 0)
        return(r);
    rt += r;

    if (txflag) {
printf("\nHier begin zweiter Teil txflag=%d\n",txflag);

        x0 = (double)txflag * 360.0;
        for (i = 0; i < n; ++i)  
            x[i] += x0;

        r = map_polygon_cyl_clip(n,x,y,sdid);
        if (r < 0)
            return(r);
        rt += r;
    }
printf("rtttt=%d\n",rt);


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

int map_polygon_cyl_clip(int n,double *x,double *y,int sdid)
{
    int m,r;  

    if (alloc_acu(2 * n + 10))
        return(-1);              
    if (alloc_acv(2 * n + 10))
        return(-1);                  
    if (alloc_acn(2 * n + 10))
        return(-1);                      
    if (alloc_acr(2 * n + 10))
        return(-1);                      
      
    m = sdclip_poly_rec(n,x,y,AcU,AcV,AcN,AcR,GDPXA,GDLatA,GDPXB,GDLatB);

printf("m===%d\n",m);

    if (m <= 0)
        return(m);
     
    if (m == 1) {                     /* polygon completely in clip region */
        map_polygon_cyl_ppol(n,x,y);
        return(1);
    }
    else if (m == 2) {                /* clip region completely in polygon */
        return(0);
        /***********
        map_polygon_cyl_pbox(GDPXA,GDLatA,GDPXB,GDLatB);
        return(1);
        *************/
    }
    else {
        r = map_polygon_cyl_plist(m,AcU,AcV,AcN,AcR,sdid);
        if (r < 0) {
            sdclip_err_msg(1,sdid);
            r = 0;
        }
        return(r);
    }
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_pbox(xmin,ymin,xmax,ymax)                               */
/*                                                                          */
/*  Plot rectangular clip region.                                           */

void map_polygon_cyl_pbox(double xmin,double ymin,double xmax,double ymax)
{
    ps_2dplot(xmin,ymin,0);
    ps_2dplot(xmax,ymin,1);
    ps_2dplot(xmax,ymax,1);
    ps_2dplot(xmin,ymax,1);
    ps_2dplot(xmin,ymin,1);
    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_ppol(n,x,y)                                             */
/*                                                                          */
/*  Plot polygon x[i], y[i] (i=0,ln-1).                                     */

void map_polygon_cyl_ppol(int n,double *x,double *y)
{
    register int i;
    double yp;

    yp = gd_proj_gety(y[0]);  
    ps_2dplot(x[0],yp,0);

    for (i = 1; i < n; ++i) {
        yp = gd_proj_gety(y[i]);  
        ps_2dplot(x[i],yp,1);
    }
    yp = gd_proj_gety(y[0]);  
    ps_2dplot(x[0],yp,1);

    if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
        fprintf(PSFd,"gsave\n");
        ps_fill(PMGS);
        fprintf(PSFd,"grestore\n");
    }
    if (PMLW > 0.0)
        fprintf(PSFd,"stroke\n");
}

/* ------------------------------------------------------------------------ */
/*  map_polygon_cyl_plist(n,x,y,nf,d,sdid)                                  */
/*                                                                          */
/*  This function gets an edge list created by sdclip_poly_rec and plots    */
/*  all polygons. Return number of polygons plotted, or -1 if error.        */

int map_polygon_cyl_plist(int n,double *x,double *y,int *nf,int *d,int sdid)
{
    register int i,j,k;
    int m,np,first;
    double yp;

    /************* 
    printf("map_polygon_plist\n");
    for (i = 0; i < n; ++i) {
        printf("i=%3d x=%f y=%f nf=%3d d=%d\n",i,x[i],y[i],nf[i],d[i]);
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
            yp = gd_proj_gety(y[k]);  
            ps_2dplot(x[k],yp,first);
            first = 1;
            k = nf[k];
            if (k == j)
                break;
        }
        if (first) {
            yp = gd_proj_gety(y[j]);  
            ps_2dplot(x[j],yp,first);

            if (PMGSFlg && PMGS >= 0.0 && PMGS <= 1.0) {
                fprintf(PSFd,"gsave\n");
                ps_fill(PMGS);
                fprintf(PSFd,"grestore\n");
            }
            if (PMLW > 0.0)
                fprintf(PSFd,"stroke\n");
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

int map_check_pole(int n,double *x,double *y)
{
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

double geod_distance(double lon1,double lat1,double lon2,double lat2)
{
    double d,slon,slat,t;

    lon1 = degree_to_arc(lon1);
    lat1 = degree_to_arc(lat1);
    lon2 = degree_to_arc(lon2);
    lat2 = degree_to_arc(lat2);
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



