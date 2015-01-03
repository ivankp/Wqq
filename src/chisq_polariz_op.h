#ifndef chisq_polariz_op_h
#define chisq_polariz_op_h

#include "chisq_polariz_base.h"

class fcnpf;

class chisq_polariz_op: public chisq_polariz_base {
protected:
  const fcnpf* dist;
  double norm_const[2];

  virtual Double_t fit_fnc(Double_t x, Double_t *par) const;

  void norm_dist();

public:
  chisq_polariz_op(const TH1F* data, const char* dist_file, unsigned nbins_fit=1000);
  chisq_polariz_op(const TH1F* data, const char* dist_file, double min, double max, unsigned nbins_fit=1000);
  virtual ~chisq_polariz_op();

  virtual void Fit();
};

#endif
