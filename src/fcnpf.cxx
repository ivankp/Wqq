#include "fcnpf.h"

#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <cctype>

using namespace std;

class fcnpfexpt: public exception {
  private:
    string msg;
  public:
    fcnpfexpt(const char* msg, const char* file) {
      this->msg = msg;
      this->msg += file;
    }
    virtual ~fcnpfexpt() throw() { }
    virtual const char* what() const throw() {
      return msg.c_str();
    }
};

fcnpf::fcnpf(const char* file) {
  ifstream dat(file);
  if (dat.is_open()) {

    int ncols = 0;
    char c[2];
    while (c[0]!='\n') {
      dat.get(c[1]);
      if ( isspace(c[1]) && (!isspace(c[0])) ) ncols++;
      c[0] = c[1];
    }
    dat.seekg(0,dat.beg); // rewind

    // push_back vectors
    for (int i=0;i<ncols;i++) x.push_back(vector<double>());

    double val;
    bool ok = true;
    while (true) {
      for (int i=0;i<ncols;i++) {
        if (dat >> val) x[i].push_back(val);
        else { ok = false; break; }
      }
      if (!ok) break;
    }
    dat.close();

    ok = true;
    for (int i=1;i<ncols;i++) {
      ok = (x[0].size() == x[i].size());
      if (!ok ) throw fcnpfexpt("Unequal length columns in ",file);
    }

    last = x[0].size() - 1;

  } else throw fcnpfexpt("Unable to open file ",file);
}

double fcnpf::operator()(int col, double val) const {
  // return 0 for out of range
  if (val<x[0][0] || x[0][last]<val) return 0.;

  // compute the lower index
  int i = (val-x[0][0])*last/(x[0][last]-x[0][0]);

  // linear interpolation
  double m = (x[col][i+1]-x[col][i])/(x[0][i+1]-x[0][i]);
  double b = x[col][i] - m*x[0][i];
  return m*val + b;
}
