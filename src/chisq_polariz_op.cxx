#include "chisq_polariz_op.h"

#include "fcnpf.h"
#include "functions.h"

chisq_polariz_op::chisq_polariz_op(const TH1F* data, const char* dist_file, unsigned nbins_fit)
: chisq_polariz_base(data,nbins_fit), dist(new fcnpf(dist_file))
{
  norm_dist();
}
chisq_polariz_op::chisq_polariz_op(const TH1F* data, const char* dist_file,
  double min, double max, unsigned nbins_fit)
: chisq_polariz_base(data,min,max,nbins_fit), dist(new fcnpf(dist_file))
{
  norm_dist();
}
chisq_polariz_op::~chisq_polariz_op()
{
  delete dist;
}

Double_t chisq_polariz_op::fit_fnc(Double_t x, Double_t *par) const
{
  return ((1.-par[0])*norm_const[1]*(*dist)(2,x) + // h = 1 (col 2)
              par[0] *norm_const[1]*(*dist)(1,x)); // h = 0 (col 1)
}

void chisq_polariz_op::Fit()
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

void chisq_polariz_op::norm_dist()
{
  const double step = (range.second-range.first)/nbins_data;
  for (int h=0;h<=1;h++) {
    double integ = 0.;
    for (unsigned i=0;i<nbins_data;i++) integ += (*dist)(h+1,step*(i+0.5));
    norm_const[h] = data_integral/(integ*step);
  }
}
