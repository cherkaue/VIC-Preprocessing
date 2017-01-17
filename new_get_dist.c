#include <stdio.h>
#include <math.h>

/***************************************************************************
  Function: double distance(double lat1, double long1, double lat2, double long2)
  Returns : distance between two locations
****************************************************************************/

#ifndef _E_RADIUS
#define E_RADIUS 6371.0         /* average radius of the earth */
#endif

#ifndef _PI
#define PI 3.1415
#endif

double get_dist(double lat1, double long1, double lat2, double long2)
{
  double theta1;
  double phi1;
  double theta2;
  double phi2;
  double dtor;
  double term1;
  double term2;
  double term3;
  double temp;
  double dist;

  dtor = 2.0*PI/360.0;
  theta1 = dtor*long1;
  phi1 = dtor*lat1;
  theta2 = dtor*long2;
  phi2 = dtor*lat2;
  term1 = cos(phi1)*cos(theta1)*cos(phi2)*cos(theta2);
  term2 = cos(phi1)*sin(theta1)*cos(phi2)*sin(theta2);
  term3 = sin(phi1)*sin(phi2);
  temp = term1+term2+term3;
  temp = (double) (1.0 < temp) ? 1.0 : temp;
  dist = E_RADIUS*acos(temp);

  return dist;
}  

#undef E_RADIUS
#undef PI
