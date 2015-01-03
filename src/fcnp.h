#ifndef fcnp_h
#define fcnp_h

#include <functional>

class fcnp : public std::function<double (double)> {
  public:
    // constructor
    fcnp(std::function<double (double)>& f, int n, double xmin, double xmax);
    fcnp(double (*f)(double), int n, double xmin, double xmax);

    fcnp(); // dummy constructor
    virtual ~fcnp(); // destructor

    // return fcnp
    virtual double operator()(double val) const;

  protected:
    int n; // size of arrays
    double *x, *y; // value arrays
    double xmin, xmax; // range

    // common code for constructors
    void init(std::function<double (double)>& f, int n, double xmin, double xmax);
};

#endif
