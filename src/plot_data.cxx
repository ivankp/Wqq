#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>

#include "LHEF2Wqq.h"

#include "TString.h"
#include "TMath.h"
#include "TH1.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;

typedef unsigned char small;

// default canvas
TCanvas canv("canv","",600,400);
string file_name, base_name;

class TH1F2 {
private:
  const char* var;
  const bool logy;
  TH1F * const p, * const m;
  static small N, Nhist;

public:
  TH1F2(const char* var, const char* title, const char* xtitle,
        Int_t nbinsx, Double_t xlow, Double_t xup, bool logy=false)
  : var(var), logy(logy),
    p(new TH1F(Form("W+%s",var),title,nbinsx,xlow,xup)),
    m(new TH1F(Form("W-%s",var),title,nbinsx,xlow,xup))
  {
    p->SetLineColor(46);
    p->SetMarkerColor(46);
    m->SetLineColor(602);
    m->SetMarkerColor(602);
    p->SetStats(false);
    m->SetStats(false);
    p->SetXTitle(xtitle);

    Nhist++;
  }
  ~TH1F2()
  {
    delete p;
    delete m;
  }
  void Fill(const int& pid, const double& x)
  {
    if (pid>0) p->Fill(x);
    else       m->Fill(x);
  }
  void Draw()
  {
    canv.Clear();
    canv.SetLogy(logy);

    double maxy = max(p->GetMaximum(),m->GetMaximum());
    if (logy) {
      p->SetMaximum(pow(maxy,1.05));
      //p->SetMinimum(pow(min(p->GetMinimum(),m->GetMinimum()),0.9));
    } else {
      p->SetMaximum(1.05*maxy);
      //p->SetMinimum(0.);
    }

    p->Draw();
    m->Draw("same");

    if (p->GetEntries()!=m->GetEntries())
      printf("Warning!!!: Unequal number of W+ and W-\n");

    TPaveText stat(0.75,0.80,0.95,0.95,"NDCbr");
    stat.SetBorderSize(1);
    stat.SetFillColor(0);
    stat.SetTextFont(42);
    stat.AddText(Form("Num events: %.0f",p->GetEntries()));
    stat.AddText(Form("Num bins: %d",p->GetNbinsX()));
    stat.AddText(Form("W+ mean: %.2f",p->GetMean()));
    stat.AddText(Form("W- mean: %.2f",m->GetMean()));
    stat.Draw();

    int p_u = p->GetBinContent(0);
    int p_o = p->GetBinContent(p->GetNbinsX()+1);
    int m_u = m->GetBinContent(0);
    int m_o = m->GetBinContent(m->GetNbinsX()+1);

    TPaveText overflow(0.91,0.64,0.99,0.79,"NDCbr");
    overflow.SetBorderSize(1);
    overflow.SetFillColor(0);
    overflow.SetTextFont(82);
    overflow.AddText(Form("+O %d",  p_o ));
    overflow.AddText(Form("-O %d",  m_o ));
    overflow.AddText(Form("+U %d", p_u ));
    overflow.AddText(Form("-U %d", m_u ));
    overflow.SetTextAlign(12);
    overflow.Draw();

    TLegend leg(0.8,0.71,0.89,0.79);
    leg.SetFillColor(0);
    leg.AddEntry(p,"W+");
    leg.AddEntry(m,"W-");
    leg.Draw();

    string save_name = base_name;
    save_name += "pdf";
    if (N==0) save_name += '(';
    else if (N==Nhist-1) save_name += ')';

    canv.SaveAs(save_name.c_str());

    N++;
  }

  static small num() { return Nhist; }
};
small TH1F2::N = 0;
small TH1F2::Nhist = 0;

int main(int argc, char* argv[])
{
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

  const size_t bn0 = file_name.rfind('/')+1;
  base_name = file_name.substr(bn0,file_name.size()-bn0-3);

  // Create LHEF Readed
  LHEFReader data(file_name,true);

  // Create histogram for W masses
  const int nbins = 100;
  TH1F2 hist[] = {
    TH1F2("mass","W mass","Mass, GeV",nbins,60.,100.),
    TH1F2("qmass","W #rightarrow qq  quarks mass^{2}","Mass^{2}, MeV^{2}",nbins,-.5,.5),
    TH1F2("energy","W energy","Energy, GeV",nbins,0.,1200.,true),
    TH1F2("Pt","W Pt","Pt, GeV",nbins,0.,1000.,true),
    TH1F2("opening_angle","W #rightarrow qq  opening angle",
          "Opening angle, rad",nbins,0.,M_PI),
    TH1F2("tilt_angle","W #rightarrow qq  tilt angle",
          "Tilt angle, rad",nbins,0.,M_PI),
    TH1F2("prod_angle","W production angle","Production angle, rad",nbins,0.,M_PI),
    TH1F2("cos_prod_angle","W cos(production angle)",
          "cos(production angle)",nbins,-1.,1.),
    TH1F2("polar_angle","W polar angle","Polar angle, rad",nbins,0.,M_PI,true),
    TH1F2("dE","W #rightarrow qq  #DeltaE","Energy, GeV",nbins,-1200.,1200.,true),
    TH1F2("dPt","W #rightarrow qq  #DeltaPt","Pt, GeV",nbins,-800.,800.,true),
    TH1F2("tilt_angle_cos","W #rightarrow qq  cos(tilt angle)",
          "cos(tilt angle)",nbins,-1.,1.),
    TH1F2("dE/p_W","W #rightarrow qq  #DeltaE/P_{W}","#DeltaE/P_{W}",nbins,-1.,1.),
    TH1F2("dPt/p_W","W #rightarrow qq  #DeltaPt/P_{W}","#DeltaPt/P_{W}",nbins,-1.,1.),
    TH1F2("dPz/p_W","W #rightarrow qq  #DeltaPz/P_{W}","#DeltaPz/P_{W}",10,-1.,1.),
    TH1F2("sPt_q/Pt_W","W #rightarrow qq  #SigmaPt_{q}/Pt_{W}","#SigmaPt_{q}/Pt_{W}",nbins,0.9,5.)
  };

  double tilt_angle, prod_angle;
  while ( data.readEvent() ) { // loop over events

    for (unsigned i=0; i<data.size(); i++) {
      int id = data[i]->id();
      int j=0;
      if (id==24 || id==-24) {
        tilt_angle = data[i]->tilt_angle();
        prod_angle = data[i]->prod_angle();

        //if (abs(cos(data[i]->polar_angle()))>0.1) continue;

        hist[j].Fill(id,data[i]->mass()); j++;
        hist[j].Fill(id,data[i]->prod(0)->mass2()*1e6);
        hist[j].Fill(id,data[i]->prod(1)->mass2()*1e6); j++;
        hist[j].Fill(id,data[i]->energy()); j++;
        hist[j].Fill(id,data[i]->Pt()); j++;
        hist[j].Fill(id,data[i]->opening_angle()); j++;
        hist[j].Fill(id,tilt_angle); j++;
        hist[j].Fill(id,prod_angle); j++;
        hist[j].Fill(id,cos(prod_angle)); j++;
        hist[j].Fill(id,data[i]->polar_angle()); j++;
        hist[j].Fill(id,data[i]->dE()); j++;
        hist[j].Fill(id,data[i]->dPt()); j++;
        hist[j].Fill(id,cos(tilt_angle)); j++;
        hist[j].Fill(id,data[i]->dE()/data[i]->momvecmag()); j++;
        hist[j].Fill(id,data[i]->dPt()/data[i]->momvecmag()); j++;
        hist[j].Fill(id,data[i]->dPz()/data[i]->momvecmag()); j++;
        hist[j].Fill(id,(data[i]->prod(0)->Pt()+data[i]->prod(1)->Pt())/data[i]->Pt());
      }
    }

  } // end while
  data.done_msg();

  for (small i=0;i<TH1F2::num();i++) hist[i].Draw();

  return 0;
}
