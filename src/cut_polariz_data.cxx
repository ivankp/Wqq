#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "LHEFReader.h"

#include "TStyle.h"
#include "TH1.h"
#include "TAxis.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGaxis.h"
#include "TString.h"

#include "inverse.h"
#include "functions.h"
#include "fcnp.h"

using namespace std;

int colors[] = {51,56,65,78,91,99}; // color scheme
TCanvas canv("canv","",600,400); // Default canvas

const double Wmass = 80.4; // GeV

bool first_hist = true;
void add_hist(TH1F& hist,int i,TLegend& leg,const char* opt="") {
  hist.SetLineColor(colors[i]);
  hist.SetMarkerColor(colors[i]);
  hist.SetLineWidth(2);
  hist.SetStats(0);
  leg.AddEntry(&hist,hist.GetName());
  if (first_hist) {
    hist.Draw(opt);
    first_hist = false;
  } else hist.Draw(Form("same%s",opt));
}

int main(int argc, char* argv[])
{
  string file_name;

  if (argc!=2) {
    cout <<"Usage "<<argv[0]<<" file.lhe"<<endl;
    return 1;
  } else {
    file_name = argv[1];
    if (file_name.substr(file_name.size()-4,4).compare(".lhe")) {
      cout <<file_name<<" doesn't have .lhe extension"<<endl;
      return 1;
    }
  }

  // Create LHEF Readed
  LHEFReader data(file_name,true);

  // Create histogram for W masses
  const int nbins = 200;
  const double min_Wene = 100., max_Wene = 1200.;
  TH1F hist_total("All","Opening angle cut for W polarization",
      nbins,min_Wene,max_Wene);
  TH1F hist_below_lowest("below min allowed","below min allowed",
     nbins,min_Wene,max_Wene);
  TH1F hist_below_cut("below cut (mostly h = 0)","below cut",
     nbins,min_Wene,max_Wene);
  TH1F hist_above_cut("above cut (mostly h = #pm1)","above cut",
     nbins,min_Wene,max_Wene);
  TH1F hist_mean_cut("mean cut angle","mean_cut",
      nbins,min_Wene,max_Wene);

  // precompute cut angle function
  fcnp cut_angle_pre(cut_angle,200,min_Wene/Wmass,max_Wene/Wmass);

  // LOOP **************************
  while ( data.readEvent() ) { // loop over events

    double W_energy, W_op_angle, op_angle_cut, op_angle_min, gamma;
    for (unsigned i=0; i<data.size(); i++) {
      int id = data[i]->id();
      if (id==24 || id==-24) {
        W_energy   = data[i]->energy();
        W_op_angle = data[i]->opening_angle();
        gamma = data[i]->gamma();
        op_angle_min = opening_angle(M_PI/2.,gamma);
        op_angle_cut = cut_angle_pre(gamma);
        hist_total.Fill(W_energy);
        if (op_angle_min-W_op_angle>1e-8) {
          hist_below_lowest.Fill(W_energy);
          printf("\n%0.10f  %0.10f\n",op_angle_min,W_op_angle);
        } else if (W_op_angle<op_angle_cut) {
          hist_below_cut.Fill(W_energy);
        } else {
          hist_above_cut.Fill(W_energy);
        }
      } // end if W
    }

  } // end while
  data.done_msg();

  // fill mean cut plot
  for (Int_t i=0;i<hist_mean_cut.GetNbinsX();i++) {
    hist_mean_cut.SetBinContent(i,
      cut_angle_pre( hist_mean_cut.GetXaxis()->GetBinCenter(i)/Wmass )
    );
  }
  // Have to do it this way because on TH1 Yaxis doesn't have bins
  const double cut_hist_max_y =
    hist_mean_cut.GetMaximum()*(gStyle->GetHistTopMargin()+1.);
  const double tot_hist_max_y =
    hist_total.GetMaximum()*(gStyle->GetHistTopMargin()+1.);
  hist_mean_cut.Scale(tot_hist_max_y/cut_hist_max_y);

  // show if there are unexpectedly low opening angles
  const int below_lowest = hist_below_lowest.GetEntries();
  if (below_lowest) {
    printf("\033[31mUnexpectedly low opening angles: %d\033[0m\n",below_lowest);
    hist_below_lowest.SetName(Form("%s (n = %d)",
      hist_below_lowest.GetName(), below_lowest
    ));
  }

  // Draw the histograms **********************
  canv.Clear();

  TLegend leg(0.5,0.6,0.84,0.85);

  TPaveText frac_text(0.67,0.45,0.84,0.59,"NDCbr");
  frac_text.SetBorderSize(1);
  frac_text.SetFillColor(0);
  frac_text.SetTextFont(42);
  frac_text.AddText("Polarization fractions");
  frac_text.AddLine(0,0.7,0,0.7);
  frac_text.AddText(Form("Long : %5.2f%%",hist_below_cut.GetEntries()*100./hist_total.GetEntries()));
  frac_text.AddText(Form("Trans: %5.2f%%",hist_above_cut.GetEntries()*100./hist_total.GetEntries()));

  hist_total.SetXTitle("W energy, GeV");
  hist_total.SetYTitle("Num counts");
  add_hist(hist_total,0,leg);
  if (below_lowest) add_hist(hist_below_lowest,5,leg);
  add_hist(hist_below_cut,3,leg);
  add_hist(hist_above_cut,2,leg);
  add_hist(hist_mean_cut,4,leg,"C");

  // second axis
  TGaxis axis2(
    hist_total.GetXaxis()->GetXmax(),0.,
    hist_total.GetXaxis()->GetXmax(),tot_hist_max_y,
    0.,cut_hist_max_y,510,"+L");
  axis2.SetTitle("Cut angle, rad");
  axis2.SetLabelFont(hist_total.GetYaxis()->GetLabelFont());
  axis2.SetTitleFont(hist_total.GetYaxis()->GetTitleFont());
  axis2.SetLabelSize(hist_total.GetYaxis()->GetLabelSize());
  axis2.SetTitleSize(hist_total.GetYaxis()->GetTitleSize());
  axis2.Draw();

  leg.SetFillColor(0);
  leg.Draw();
  frac_text.Draw();

  canv.SaveAs("polariz_cut.pdf");

  return 0;
}
