#ifndef chisq_polariz_tilt_h
#define chisq_polariz_tilt_h

#include "chisq_polariz_base.h"

class chisq_polariz_tilt: public chisq_polariz_base {
protected:
  virtual Double_t fit_fnc(Double_t x, Double_t *par) const;

public:
  chisq_polariz_tilt(const TH1F* data, unsigned nbins_fit=1000);
  chisq_polariz_tilt(const TH1F* data, double min, double max, unsigned nbins_fit=1000);
  virtual ~chisq_polariz_tilt();

  virtual void Fit();
};

#endif
