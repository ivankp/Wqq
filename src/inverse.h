#ifndef inverse_h
#define inverse_h

#include <functional>

typedef std::function<double (double)> ddfcn;

class inverse : public ddfcn {
  public:
    // constructor
    inverse(const ddfcn& f, int n, double xmin, double xmax);
    // NOTE: f must be monotonic

    // copy constructor
    inverse(const inverse& other);

    inverse(); // dummy constructor
    virtual ~inverse(); // destructor

    // return inverse
    virtual double operator()(double val) const;
    virtual double operator()(double val, double& der) const; // also return derivative

  protected:
    const int n; // size of arrays
    // value arrays
    double * const x;
    double * const y;

    mutable double derivative;

    virtual double eval(double val) const;
};

#endif
