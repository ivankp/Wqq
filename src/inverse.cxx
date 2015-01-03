#include "inverse.h"

#include <algorithm>

inverse::inverse(const ddfcn& f, int n, double xmin, double xmax)
: n(n), x(new double[n]), y(new double[n])
{
  double step = (xmax-xmin)/(n-1);
  if (f(xmin) < f(xmax)) {
    for (int i=0;i<n;i++) {
      x[i] = xmin+step*i; // inverse value
      y[i] = f(x[i]);     // inverse coordinate
    }
  } else {
    for (int i=n-1;i>=0;i--) {
      x[i] = xmax-step*i; // inverse value
      y[i] = f(x[i]);     // inverse coordinate
    }
  }
}
inverse::inverse() : n(0), x(NULL), y(NULL) { }
inverse::inverse(const inverse& other) : n(other.n), x(new double[n]), y(new double[n])
{
  std::copy(other.x, other.x+n, x);
  std::copy(other.y, other.y+n, y);
}
inverse::~inverse() {
  if (x) delete[] x;
  if (y) delete[] y;  
}

double inverse::eval(double val) const {
  // return 0 for undefined values
  if (val<y[0] || y[n-1]<val) return 0.;

  // find encompasing indices
  int i1, i2;

  // initial guess : linear in y
  int i = (val-y[0])*(n-1)/(y[n-1]-y[0]);

  int dir = 1;
  double di = val-y[i];
  if (di==0.) return x[i];
  else {
    if (di<0.) { dir = -1; } // determine direction

    while (dir*val > dir*y[i]) i += dir;

    if (dir==1) {i1=i-1; i2=i  ;}
    else        {i1=i  ; i2=i+1;}
  }
  
  // linear interpolation
  derivative = (x[i2]-x[i1])/(y[i2]-y[i1]);
  double b = x[i1] - derivative*y[i1];
  return derivative*val + b;
}

double inverse::operator()(double val) const {
  return eval(val);
}

double inverse::operator()(double val, double& der) const {
  double ret = eval(val);
  if (val<y[0] || y[n-1]<val) der = 0.;
  else der = derivative;
  return ret;
}
