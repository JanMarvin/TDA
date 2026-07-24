/* t_map.h */

#ifndef _TMAP_H
#define _TMAP_H

/*  functions in t_map.c */

int psetupg(void); 
int sdpgrat(void);
int sdpgeo(void); 
int sdpmap(void); 
double geod_distance(double lon1,double lat1,double lon2,double lat2);

#endif /* _TMAP_H */

/* end of t_map.h */


