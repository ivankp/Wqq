#include "chisq_polariz_ediff_abs.h"

#include "functions.h"

chisq_polariz_ediff_abs::chisq_polariz_ediff_abs(const TH1F* data, unsigned nbins_fit)
: chisq_polariz_base(data,nbins_fit)
{

}
chisq_polariz_ediff_abs::chisq_polariz_ediff_abs(const TH1F* data, double min, double max, unsigned nbins_fit)
: chisq_polariz_base(data,min,max,nbins_fit)
{

}
chisq_polariz_ediff_abs::~chisq_polariz_ediff_abs() { }

Double_t chisq_polariz_ediff_abs::fit_fnc(Double_t x, Double_t *par) const
{
    return ( par[0] *tilt_angle_abs_cos_prob(x, 0) + // h = 0
        (1.-par[0]) *tilt_angle_abs_cos_prob(x, 1) ) // h = 1
    * data_integral;
}

void chisq_polariz_ediff_abs::Fit()
{
  TMinuitObj minuit(1); // arg - max num params
  minuit.SetPrintLevel(-1); // -1 for nothing
  minuit.SetFCN(this);

  minuit.DefineParameter(
    0,       // parameter number
    "h0",    // parameter name
    0.5,     // start value
    0.01,    // step size
    0.,      // mininum
    1.       // maximum
  );

  minuit.Migrad(); // run Migrad minimization algorithm

  chisq_polariz_base::Fit();
}
