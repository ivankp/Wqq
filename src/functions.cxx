#include "functions.h"

#include <cmath>
#include <functional>

#include "inverse.h"

// theoretical functions ********************************************
double opening_angle(double tilt_angle, double gamma)
{
  // beta gamma sin
  double bgs = sqrt(gamma*gamma - 1.)*sin(tilt_angle);
  return atan2( 2.*bgs, bgs*bgs - 1. );
}
double min_opening_angle(double gamma)
{
  //return opening_angle(M_PI/2.,gamma);
  return atan2( 1.,sqrt(gamma*gamma - 1.) )*2.;
}
double opening_angle_over_min(double tilt_angle, double gamma)
{
  return opening_angle(tilt_angle,gamma)/min_opening_angle(gamma);
}

double tilt_angle_prob(double tilt_angle, int helicity)
{
  if (helicity==1) {
    return 0.375*pow(1.+cos(tilt_angle),2)*sin(tilt_angle);
  } else if (helicity==0) {
    return 0.75*pow(sin(tilt_angle),3);
  } else if (helicity==-1) {
    return 0.375*pow(1.-cos(tilt_angle),2)*sin(tilt_angle);
  } else return 0.;
}

double tilt_angle_cos_prob(double tilt_angle_cos, int helicity)
{
  if (helicity==1) {
    return 0.375*pow(1.+tilt_angle_cos,2);
  } else if (helicity==0) {
    return 0.75*( 1.-pow(tilt_angle_cos,2) );
  } else if (helicity==-1) {
    return 0.375*pow(1.-tilt_angle_cos,2);
  } else return 0.;
}

double tilt_angle_abs_cos_prob(double tilt_angle_cos, int helicity)
{
  return tilt_angle_cos_prob( tilt_angle_cos,helicity)
       + tilt_angle_cos_prob(-tilt_angle_cos,helicity);
}

double opening_angle_prob(double angle, int helicity, const inverse& tilt_from_op)
{
  double derivative;
  double tilt_angle = tilt_from_op(angle,derivative);
  derivative = std::abs(derivative);
  if (tilt_angle==0.) return 0;
  else return tilt_angle_prob(      tilt_angle ,helicity)*derivative +
              tilt_angle_prob((M_PI-tilt_angle),helicity)*derivative;
}

// ******************************************************************

double tilt_angle(double angle, double gamma)
{
  return asin( 1./( tan(angle/2.) * sqrt(gamma*gamma-1.) ) );
}
double tilt_angle(double angle, double gamma, double& derivative)
{
  double half_angle = angle/2.;
  double a = sqrt(gamma*gamma-1.);
  double b = 1./( tan(half_angle) * a );
  double c = sin(half_angle);
  derivative = -1./( 2.*c*c*a*sqrt(1.-b) );
  return asin(b);
}
double opening_angle_prob(double angle, int helicity, double gamma)
{
  double derivative;
  double tilt = tilt_angle(angle,gamma,derivative);
  derivative = std::abs(derivative);
  if (tilt==0.) return 0;
  else return tilt_angle_prob(      tilt ,helicity)*derivative +
              tilt_angle_prob((M_PI-tilt),helicity)*derivative;
}

// ******************************************************************

double cut_angle(double gamma, int num_tilt_op, int num_opt)
{
  std::function<double (double)> f = bind(opening_angle,std::placeholders::_1,gamma);
  inverse tilt_from_op(f,num_tilt_op,0,M_PI/2.); // turns out doesn't have to be so precise

  const unsigned n = num_opt;

  double x[n], y[2][n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[0][i] = opening_angle_prob(x[i],0,tilt_from_op);
    y[1][i] = opening_angle_prob(x[i],1,tilt_from_op);
  }

  // ************************
  // Integrate and compute **
  for (unsigned i=n-1;i!=0;i--) {
    for (unsigned j=0;j<i;j++)
      for (unsigned k=0;k<2;k++) y[k][i] += y[k][j];
    if (y[1][i]>0.) y[0][i] /= sqrt(y[1][i]);
  }
  // ************************
  // ************************

  // find maximum ***********
  // ************************
  unsigned maxi = 2./(gamma*step); // start from derived minimum angle
  while (y[0][maxi]<y[0][n-1]) maxi++;
  while (y[0][maxi+1] > y[0][maxi]) maxi++;
  // ************************
  // ************************

  return x[maxi];
}
double cut_angle(double gamma) { return cut_angle(gamma,100,200); }
