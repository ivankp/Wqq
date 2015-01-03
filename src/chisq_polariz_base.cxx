#include "chisq_polariz_base.h"

#include <cmath>
#include <iostream>

#include "TH1.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TPaveText.h"

chisq_polariz_base::chisq_polariz_base(const TH1F* data, unsigned nbins_fit)
: range(data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax()),
  nbins_data(data->GetNbinsX()), nbins_fit(nbins_fit),
  data(new TH1F(*data)), fit(NULL), chisq(NULL),
  data_integral( data->Integral(
      data->FindFixBin(range.first), data->FindFixBin(range.second)-1, "width"
  ) ),
  nbinavg(25)
{ }

chisq_polariz_base::chisq_polariz_base(const TH1F* data, double min, double max, unsigned nbins_fit)
: range(min,max),
  nbins_data(data->GetNbinsX()), nbins_fit(nbins_fit),
  data(new TH1F(*data)), fit(NULL), chisq(NULL),
  data_integral( data->Integral(
      data->FindFixBin(range.first), data->FindFixBin(range.second)-1, "width"
  ) ),
  nbinavg(25)
{ }

chisq_polariz_base::~chisq_polariz_base()
{
  delete data;
  if (fit) delete fit;
  if (chisq) delete chisq;
}

void chisq_polariz_base::operator()
(Int_t npar, Double_t *grad, Double_t &fval, Double_t *par, Int_t flag) const
{
  //std::cout << nbinavg << std::endl;
  double chisq=0., delta, y, f, err;
  for (int i=1;i<=data->GetNbinsX();i++)
  {
    // data
    y = data->GetBinContent(i);
    y==0. ? err = 1 : err = sqrt(y);

    // average over theory
    if (nbinavg>0) {
      double step = data->GetBinWidth(i)/nbinavg;
      f = 0.;
      for (double x=data->GetBinLowEdge(i);x<data->GetBinLowEdge(i+1);x+=step)
        f += fit_fnc(x,par);
      f /= nbinavg;
    }
    // or take the bin center
    else f = fit_fnc(data->GetBinCenter(i),par);

    // accumulate chi-squared
    delta = (y-f)/err;
    chisq += delta*delta;
  }

  // chi-squared per degree of freedom
  fval = chisq/(data->GetNbinsX()-npar); // return
}

double chisq_polariz_base::get_chisq() const
{
  Double_t fval;
  operator()(npar,NULL,fval,par,0);
  return fval;
}
double chisq_polariz_base::get_chisq(double x) const
{
  if (npar==1) {
    Double_t fval;
    operator()(npar,NULL,fval,&x,0);
    return fval;
  } else return 0.;
}

double chisq_polariz_base::get_frac(unsigned i) const // 1:+  0:0  2:-
{
  if (npar==2) {
    switch (i) {
      case 0: return par[0]; break;
      case 1: return par[1]; break;
      case 2: return 1.-par[0]-par[1]; break;
      default: return 0.; break;
    }
  } else if (npar==1) {
    switch (i) {
      case 0: return    par[0]; break;
      case 1: return 1.-par[0]; break;
      default: return 0.;
    }
  } else return 0.;
}
double chisq_polariz_base::get_unc (unsigned i) const
{
  if (npar==2) {
    switch (i) {
      case 0: return err[0]; break;
      case 1: return err[1]; break;
      case 2: return sqrt( err[0]*err[0] + err[1]*err[1] ); break;
      default: return 0.;
    }
  } else if (npar==1) return err[0];
    else return 0.;
}

void chisq_polariz_base::Fit()
{
  fit = new TH1F(Form("%s_%f_%f",data->GetName(),range.first,range.second),"",
                 nbins_fit,range.first,range.second);
  for (unsigned i=0;i<nbins_fit;i++)
    fit->SetBinContent(i,fit_fnc(fit->GetXaxis()->GetBinCenter(i),par));

  //data->Scale(1./data_integral);
  //fit ->Scale(1./data_integral);
}

void chisq_polariz_base::ChiSq()
{
  chisq = new TH1F(Form("%s_chisq",fit->GetName()),"",1000,0.,1.);
  for (unsigned i=0;i<nbins_fit;i++)
    chisq->SetBinContent(i,get_chisq(chisq->GetXaxis()->GetBinCenter(i)));
}

void chisq_polariz_base::Draw(chisq_polariz_base* fit_p, chisq_polariz_base* fit_m)
{
  if (leg)        delete[] leg;
  if (frac_text)  delete[] frac_text;
  if (nbins_text) delete   nbins_text;

  chisq_polariz_base* const chisq[2] = { fit_p, fit_m };
  TH1* const data [2] = { fit_p->data, fit_m->data };
  TH1* const fit  [2] = { fit_p->fit,  fit_m->fit  };

  leg = new TLegend[2] {
    TLegend(0.15,0.84,0.30,0.905),
    TLegend(0.73,0.84,0.88,0.905)
  };

  frac_text = new TPaveText[2] {
    TPaveText(0.15,0.73,0.30,0.84,"NDCbr"),
    TPaveText(0.73,0.73,0.88,0.84,"NDCbr")
  };

  nbins_text = new TPaveText(0.91,0.75,0.99,0.9,"NDCbr");
  nbins_text->SetBorderSize(1);
  nbins_text->SetFillColor(0);
  nbins_text->SetTextFont(42);
  nbins_text->AddText(Form("nbins: %d",data[0]->GetNbinsX()));

  //data[0]->SetMinimum(0.);
  data[0]->SetMaximum(1.050*std::max(data[0]->GetMaximum(),data[1]->GetMaximum()));
  data[0]->SetMinimum(0.975*std::min(data[0]->GetMinimum(1),data[1]->GetMinimum(1)));

  for (unsigned char i=0;i<2;i++) {
    data[i]->SetLineColor(colors[i]);
    data[i]->SetMarkerColor(colors[i]);
    data[i]->SetLineWidth(2);
    data[i]->SetStats(0);
    if (i==0) {
      data[i]->Draw();
    } else data[i]->Draw("same");
    fit[i]->SetLineColor(colors[i]);
    fit[i]->SetMarkerColor(colors[i]);
    fit[i]->Draw("same");

    leg[i].AddEntry(data[i],Form("W%c [%d]",sign[i],(int)data[i]->GetEntries()));
    leg[i].SetFillColor(0);
    leg[i].SetFillStyle(0);
    leg[i].SetBorderSize(0);
    leg[i].Draw();

    frac_text[i].SetBorderSize(1);
    frac_text[i].SetFillColor(0);
    frac_text[i].SetTextFont(42);
    if (fit_p->npar==2) {

      frac_text[i].AddText(Form("h_{+} = %.3f#pm%.3f",chisq[i]->par[1],chisq[i]->err[1]));
      frac_text[i].AddText(Form("h_{0} = %.3f#pm%.3f",chisq[i]->par[0],chisq[i]->err[0]));
      frac_text[i].AddText(Form("h_{#font[122]{-}} = %.3f#pm%.3f",
        1.-chisq[i]->par[0]-chisq[i]->par[1],
        sqrt( chisq[i]->err[0]*chisq[i]->err[0] + chisq[i]->err[1]*chisq[i]->err[1] )
      ));

    } else if (fit_p->npar==1) {

      frac_text[i].AddText(Form(  "h_{0} = %.3f#pm%.3f",   chisq[i]->par[0],chisq[i]->err[0]));
      frac_text[i].AddText(Form("h_{#pm} = %.3f#pm%.3f",1.-chisq[i]->par[0],chisq[i]->err[0]));

      nbins_text->AddText(Form("#chi^{2}_{%c}=%.2f",sign[i],chisq[i]->get_chisq()));

    }
    frac_text[i].Draw();

    // compare integrals
    #if compare_integrals
      printf("∫ data[%d] : %.3f\n",i,data[i]->Integral(
        data[i]->FindFixBin(chisq[i]->range.first),
        data[i]->FindFixBin(chisq[i]->range.second)-1,
        "width"
      ));
      printf("∫ fit [%d] : %.3f\n",i,fit [i]->Integral(
        fit[i]->FindFixBin(chisq[i]->range.first),
        fit[i]->FindFixBin(chisq[i]->range.second)-1,
        "width"
      ));
    #endif
  }

  nbins_text->Draw();

}

void chisq_polariz_base::DrawChiSq(chisq_polariz_base* fit_p, chisq_polariz_base* fit_m)
{
  if (leg) delete[] leg;

  chisq_polariz_base* const chisq[2] = { fit_p, fit_m };

  for (unsigned char i=0;i<2;i++) if (!chisq[i]->chisq) chisq[i]->ChiSq();
  TH1* const hist[2] = { fit_p->chisq, fit_m->chisq };

  leg = new TLegend[2] {
    TLegend(0.15,0.84,0.30,0.905),
    TLegend(0.73,0.84,0.88,0.905)
  };
/*
  hist[0]->SetMinimum(0.90*std::min(
    hist[0]->GetMinimum(0.1),
    hist[1]->GetMinimum(0.1)
  ));
  hist[0]->SetMaximum(1.05*std::max(
    hist[0]->GetMaximum(),
    hist[1]->GetMaximum()
  ));
*/

  hist[0]->SetMinimum(0);
  hist[0]->SetMaximum(5);

  hist[0]->SetTitle("#chi^{2}");
  hist[0]->SetXTitle("Parameter h_{0}");

  for (unsigned char i=0;i<2;i++) {
    hist[i]->SetLineColor(colors[i]);
    hist[i]->SetMarkerColor(colors[i]);
    hist[i]->SetLineWidth(2);
    hist[i]->SetStats(0);
    if (i==0) {
      hist[i]->Draw();
    } else hist[i]->Draw("same");

    leg[i].AddEntry(hist[i],Form("W%c",sign[i]));
    leg[i].SetFillColor(0);
    leg[i].SetFillStyle(0);
    leg[i].SetBorderSize(0);
    leg[i].Draw();
  }
}

void chisq_polariz_base::set_nbinavg(unsigned n) { nbinavg = n; }

const int  chisq_polariz_base::colors[2] = {46, 602};
const char chisq_polariz_base::sign  [2] = {'+','-'};
TLegend* chisq_polariz_base::leg = NULL;
TPaveText* chisq_polariz_base::frac_text = NULL;
TPaveText* chisq_polariz_base::nbins_text = NULL;
