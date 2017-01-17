#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <math.h>


double get_dist(double,double,double,double);

main(int argc, char *argv[]) {
/***********************************************************************
  cell_area.c           Keith Cherkauer             

  This program computes the area being routed using the basin, or
  sub-basin flow direction file, and the drainage fraction file.  

***********************************************************************/

  
  double cell_lat, cell_lng;
  double ll_lat, ll_lng;
  double tmplng, tmplat;
  double cellsize;
  double fract;
  double tmpsum, areasum;
  double routarea, tmprout;
  double NODATA;
  double ddist;

  if(argc!=5) {
    fprintf(stderr,"Usage: %s <longitude> <latitude> <longitude> <latitude>\n",argv[0]);
    exit(0);
  }


  cell_lat = atof(argv[2]);
  cell_lng = atof(argv[1]);
  tmplat = atof(argv[4]);
  tmplng = atof(argv[3]);
  
  ddist = get_dist(cell_lat,cell_lng,tmplat,tmplng);


  fprintf(stdout,"%lf\n",ddist);

}
