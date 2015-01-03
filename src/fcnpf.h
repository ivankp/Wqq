#ifndef fcnpf_h
#define fcnpf_h

#include <vector>

class fcnpf {
  public:
    // constructor
    fcnpf(const char* file);

    fcnpf() { } // dummy constructor
    virtual ~fcnpf() { } // destructor

    // return fcnpf
    virtual double operator()(int col, double val) const;

  protected:
    std::vector< std::vector<double> > x; // value arrays

    int last;
};

#endif
