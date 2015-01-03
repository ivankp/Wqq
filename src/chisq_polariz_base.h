#ifndef chisq_polariz_base_h
#define chisq_polariz_base_h

#define compare_integrals 0

#include "TMinuitObj.h"

#include <utility>

class TH1F;
class TLegend;
class TPaveText;

class chisq_polariz_base: public MnFcnBase {
protected:
  const std::pair<double,double> range; // variable range
  const unsigned nbins_data, nbins_fit;
  TH1F * const data;
  TH1F * fit;
  TH1F * chisq;

  const double data_integral;

  unsigned nbinavg;

  // Must use functions normalized to unity
  virtual Double_t fit_fnc(Double_t x, Double_t *par) const =0;

  static TLegend* leg;
  static TPaveText *frac_text, *nbins_text;

  virtual void ChiSq();

public:
  chisq_polariz_base(const TH1F* data, unsigned nbins_fit=1000);
  chisq_polariz_base(const TH1F* data, double min, double max, unsigned nbins_fit=1000);
  virtual ~chisq_polariz_base();

  virtual void operator()
  (Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) const;

  static void Draw(chisq_polariz_base* fit_p, chisq_polariz_base* fit_m);
  static void DrawChiSq(chisq_polariz_base* fit_p, chisq_polariz_base* fit_m);

  virtual void Fit();

  virtual double get_frac(unsigned i) const; // 1:+  0:0  2:-
  virtual double get_unc (unsigned i) const;

  virtual double get_chisq() const;
  virtual double get_chisq(double x) const;

  virtual void set_nbinavg(unsigned n);

  static const int  colors[2];
  static const char sign  [2];
};

#endif
