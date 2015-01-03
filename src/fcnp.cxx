#include "fcnp.h"

void fcnp::init(std::function<double (double)>& f,
                      int n, double xmin, double xmax) {
  this->n = n;
  x = new double[n];
  y = new double[n];
  if (xmin < xmax) {
    this->xmin = xmin;
    this->xmax = xmax;
  } else {
    this->xmin = xmax;
    this->xmax = xmin;
  }
  double step = (this->xmax-this->xmin)/(n-1);
  for (int i=0;i<n;i++) {
    x[i] = this->xmin+step*i;
    y[i] = f(x[i]);
  }
}
fcnp::fcnp(std::function<double (double)>& f,
                       int n_, double xmin, double xmax) {
  init(f,n,xmin,xmax);
}
fcnp::fcnp(double (*f)(double),
                       int n, double xmin, double xmax) {
  std::function<double (double)> fcn = f;
  init(fcn,n,xmin,xmax);
}
fcnp::fcnp() {
  x=0; y=0;
}
fcnp::~fcnp() {
  if (x) delete[] x;
  if (y) delete[] y;  
}

double fcnp::operator()(double val) const {
  // return 0 for out of range
  if (val<xmin || xmax<val) return 0.;

  // compute the lower index
  int i = (val-xmin)*(n-1)/(xmax-xmin);

  // linear interpolation
  double m = (y[i+1]-y[i])/(x[i+1]-x[i]);
  double b = y[i] - m*x[i];
  return m*val + b;
}
