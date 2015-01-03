#ifndef chisq_polariz_ediff_abs_h
#define chisq_polariz_ediff_abs_h

#include "chisq_polariz_base.h"

class chisq_polariz_ediff_abs: public chisq_polariz_base {
protected:
  virtual Double_t fit_fnc(Double_t x, Double_t *par) const;

public:
  chisq_polariz_ediff_abs(const TH1F* data, unsigned nbins_fit=1000);
  chisq_polariz_ediff_abs(const TH1F* data, double min, double max, unsigned nbins_fit=1000);
  virtual ~chisq_polariz_ediff_abs();

  virtual void Fit();
};

#endif
