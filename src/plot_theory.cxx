#include <cmath>
#include <cstdio>

#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGaxis.h"
#include "TMarker.h"
#include "TLine.h"

#include "inverse.h"
#include "functions.h"

#define defnpts 1000

using namespace std;

int colors[] = {51,56,65,78,91,99}; // color scheme
TCanvas canv("canv","",600,400); // Default canvas

// Graphs ***********************************************************
TGraph* opening_angle_graph(double gamma, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = opening_angle(x[i],gamma)/M_PI;
    x[i] /= M_PI;
  }
  return new TGraph(n,x,y);
}
TGraph* min_opening_angle_graph(unsigned n=defnpts)
{
  double x[n], y[n];
  double step = (100.-1.)/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i+1.;
    y[i] = min_opening_angle(x[i])/M_PI;
  }
  return new TGraph(n,x,y);
}
TGraph* tilt_angle_prob_graph(int helicity, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = tilt_angle_prob(x[i],helicity);
    x[i] /= M_PI;
  }
  return new TGraph(n,x,y);
}
TGraph* tilt_angle_cos_prob_graph(int helicity, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = 2./(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i-1.;
    y[i] = tilt_angle_cos_prob(x[i],helicity);
  }
  return new TGraph(n,x,y);
}
TGraph* tilt_angle_abs_cos_prob_graph(int helicity, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = 1./(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = tilt_angle_abs_cos_prob(x[i],helicity);
  }
  return new TGraph(n,x,y);
}

/*TGraph* tilt_angle_graph(double gamma, unsigned n=defnpts)
{
  using namespace std::placeholders;
  function<double (double)> f = bind(opening_angle,_1,gamma);
  inverse g(f,200,0,M_PI/2.);
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = g(x[i])/M_PI;
    x[i] /= M_PI;
  }
  return new TGraph(n,x,y);
}*/
TGraph* tilt_angle_graph(double gamma, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = tilt_angle(x[i],gamma)/M_PI;
    x[i] /= M_PI;
  }
  return new TGraph(n,x,y);
}
TGraph* opening_angle_prob_graph(int helicity, double gamma, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[i] = opening_angle_prob(x[i],helicity,gamma);

    if ( !isnormal(y[i]) ) y[i] = 0.;
    x[i] /= M_PI;
  }

  // Normalization is needed to correct numeric error
  #if true
  // Integral ***************
  double integral = 0.;
  for (unsigned i=0;i<n;i++) integral += y[i];
  integral *= step;

  // Normalize **************
  for (unsigned i=0;i<n;i++) y[i] /= integral;
  #endif

  return new TGraph(n,x,y);
}
TGraph* opening_angle_prob_int_graph(int helicity, double gamma, unsigned n=defnpts)
{
  TGraph* gr = opening_angle_prob_graph(helicity,gamma,n);
  double * const y = gr->GetY();
  const double step = gr->GetX()[1]*M_PI;

  for (unsigned i=0;i<n;i++) y[i] *= step;

  // Integrate **************
  for (int i=n-1;i!=0;i--) for (int j=0;j<i;j++) y[i] += y[j];

  return gr;
}

TGraph* tilt_angle_deriv_graph(double gamma, unsigned n=defnpts)
{
  double x[n], y[n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    tilt_angle(x[i],gamma,y[i]);
    x[i] /= M_PI;
  }
  return new TGraph(n,x,y);
}

// Plots ************************************************************
void draw_opening_angle()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 6;
  double gamma[ngr] = {1.,1.2,2.,5.,10.,100};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = opening_angle_graph(gamma[i]);
    gr[i]->SetLineColor(colors[i]);
    gr[i]->SetMarkerColor(colors[i]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    TString name("#gamma = ");
    name += gamma[i];
    leg.AddEntry(gr[i],name);
  }

  mg->Draw("ac");
  mg->SetTitle("Opening angle");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetYaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("Tilt angle/#pi");
  mg->GetYaxis()->SetTitle("Opening angle/#pi");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("opening_angle.pdf");

  delete mg;
}
void draw_min_opening_angle()
{
  canv.Clear();
  TGraph* gr;

  gr = min_opening_angle_graph();
  gr->SetLineColor(602);
  gr->SetMarkerColor(602);
  gr->SetFillColor(0);
  gr->SetLineWidth(2);

  gr->Draw("ac");
  gr->SetTitle("Minimum opening angle");
  gPad->Update();

  gr->GetXaxis()->SetRangeUser(1.,100.);
  gr->GetYaxis()->SetRangeUser(0.,1.);
  gr->GetXaxis()->SetTitle("#gamma");
  gr->GetYaxis()->SetTitle("Opening angle/#pi");

  canv.SaveAs("min_opening_angle.pdf");

  delete gr;
}

void draw_tilt_angle_prob()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 3;
  int helicity[ngr] = {1,0,-1};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = tilt_angle_prob_graph(helicity[i]);
    gr[i]->SetLineColor(colors[i+2]);
    gr[i]->SetMarkerColor(colors[i+2]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    leg.AddEntry(gr[i],Form("h = %d",helicity[i]));
  }

  mg->Draw("ac");
  mg->SetTitle("Tilt angle probability");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetYaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("Tilt angle/#pi");
  mg->GetYaxis()->SetTitle("Probability");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("tilt_angle_prob.pdf");

  delete mg;
}

void draw_tilt_angle_cos_prob()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 3;
  int helicity[ngr] = {1,0,-1};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = tilt_angle_cos_prob_graph(helicity[i]);
    gr[i]->SetLineColor(colors[i+2]);
    gr[i]->SetMarkerColor(colors[i+2]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    leg.AddEntry(gr[i],Form("h = %d",helicity[i]));
  }

  mg->Draw("ac");
  mg->SetTitle("Tilt angle cosine probability");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(-1.,1.);
  //mg->GetYaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("cos(Tilt angle)");
  mg->GetYaxis()->SetTitle("Probability");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("tilt_angle_cos_prob.pdf");

  delete mg;
}

void draw_tilt_angle_abs_cos_prob()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 3;
  int helicity[ngr] = {1,0,-1};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = tilt_angle_abs_cos_prob_graph(helicity[i]);
    gr[i]->SetLineColor(colors[i+2]);
    gr[i]->SetMarkerColor(colors[i+2]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    leg.AddEntry(gr[i],Form("h = %d",helicity[i]));
  }

  mg->Draw("ac");
  mg->SetTitle("Folded tilt angle cosine probability");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("|cos(Tilt angle)|");
  mg->GetYaxis()->SetTitle("Probability");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("tilt_angle_abs_cos_prob.pdf");

  delete mg;
}

void draw_tilt_angle()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 5;
  double gamma[ngr] = {1.2,2.,5.,10.,100};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = tilt_angle_graph(gamma[i]);
    gr[i]->SetLineColor(colors[i+1]);
    gr[i]->SetMarkerColor(colors[i+1]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    TString name("#gamma = ");
    name += gamma[i];
    leg.AddEntry(gr[i],name);
  }

  mg->Draw("ap");
  mg->SetTitle("Tilt angle");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetYaxis()->SetRangeUser(0.,0.5);
  mg->GetXaxis()->SetTitle("Opening angle/#pi");
  mg->GetYaxis()->SetTitle("Tilt angle/#pi");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("tilt_angle.pdf");

  delete mg;
}

void draw_opening_angle_prob(double gamma)
{
  canv.Clear();
  canv.SetLogy(1);

  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 3;
  int helicity[ngr] = {1,0,-1};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);
  leg.SetHeader(Form("#gamma = %.1f",gamma));

  int draw_order[] = {2,1,0};
  for (unsigned i, j=0;j<ngr;j++) {
    i=draw_order[j];
    gr[i] = opening_angle_prob_graph(helicity[i],gamma);
    gr[i]->SetLineColor(colors[i+2]);
    gr[i]->SetMarkerColor(colors[i+2]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    leg.AddEntry(gr[i],Form("h = %d",helicity[i]));
  }

  mg->Draw("ap");
  mg->SetTitle("Opening angle probability");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->SetMinimum(1e-6);
  mg->GetYaxis()->SetTitleOffset(1.4);
  mg->GetXaxis()->SetTitle("Opening angle/#pi");
  mg->GetYaxis()->SetTitle("Probability density");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("opening_angle_prob.pdf");

  delete mg;
  canv.SetLogy(0);
}

void draw_opening_angle_prob_int(double gamma)
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 2;
  int helicity[ngr] = {1,0};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);
  leg.SetHeader(Form("#gamma = %.1f",gamma));

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = opening_angle_prob_int_graph(helicity[i],gamma);
    gr[i]->SetLineColor(colors[i+2]);
    gr[i]->SetMarkerColor(colors[i+2]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    leg.AddEntry(gr[i],Form("h = %d",helicity[i]));
  }

  mg->Draw("ac");
  mg->SetTitle("Opening angle probability integral");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("Opening angle/#pi");
  mg->GetYaxis()->SetTitle("Probability integral");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("opening_angle_prob_int.pdf");

  delete mg;
}

void draw_cut_acceptance(double gamma, double C, bool reverse = false)
{
  const unsigned n = 10000;
  const unsigned ngr = 2;

  TGraph* integral[] = {
    opening_angle_prob_graph(0,gamma,n),
    opening_angle_prob_graph(1,gamma,n)
  };
  double* x  = integral[0]->GetX();
  double* y[] = { integral[0]->GetY(), integral[1]->GetY() };

  // Integrate **************
  for (unsigned k=0;k<ngr;k++)
    if (!reverse) {
      for (unsigned i=n-1;i!=0;i--)
        for (unsigned j=0;j<i;j++) y[k][i] += y[k][j];
    } else {
      for (unsigned i=0;i<n;i++)
        for (unsigned j=n-1;j>i;j--) y[k][i] += y[k][j];
    }
  // ************************

  // Normalize **************
  if (!reverse) {
    for (unsigned i=0;i<n;i++) {
      y[0][i] *= (C/y[0][n-1]);
      y[1][i] *= ((1.-C)/y[1][n-1]);
    }
  } else {
    for (unsigned i=0;i<n;i++) {
      y[0][i] *= (C/y[0][0]);
      y[1][i] *= ((1.-C)/y[1][0]);
    }
  }
  // ************************

  TGraph* gr[] = { new TGraph(n,x,y[0]), new TGraph(n,x,y[1]) };

  for (unsigned k=0;k<ngr;k++)
    for (unsigned i=0;i<n;i++)
      if (y[!k][i]>0.) gr[k]->GetY()[i] = y[k][i]/sqrt(y[!k][i]);

  // find maximum ***********
  // ************************
  unsigned maxi = 2./(gamma*M_PI/(n-1)); // start from derived minimum angle
  while (gr[0]->GetY()[maxi]<gr[0]->GetY()[n-1]) maxi++;
  while (gr[0]->GetY()[maxi+1] > gr[0]->GetY()[maxi]) maxi++;
  // ************************
  // ************************

  // Draw *******************
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  //TLegend leg(0.91,0.8,0.99,0.9);
  TLegend leg(0.15,0.5,0.23,0.3);
  //leg.SetHeader(Form("#gamma = %.1f",gamma));

  for (unsigned i=0;i<ngr;i++) {
    gr[i]->SetLineColor(colors[3-i]);
    gr[i]->SetMarkerColor(colors[3-i]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
  }
  for (unsigned i=0;i<ngr;i++) {
    integral[i]->SetLineColor(colors[3-i]);
    integral[i]->SetMarkerColor(colors[3-i]);
    integral[i]->SetFillColor(0);
    integral[i]->SetLineWidth(2);
    integral[i]->SetLineStyle(3);
    mg->Add(integral[i]);
  }

  leg.AddEntry(gr[0],"#scale[0.8]{#epsilon_{0} /#sqrt{#epsilon_{1}}}");
  leg.AddEntry(gr[1],"#scale[0.8]{#epsilon_{1} /#sqrt{#epsilon_{0}}}");
  leg.AddEntry(integral[0],"#epsilon_{0}");
  leg.AddEntry(integral[1],"#epsilon_{1}");

  mg->Draw("ac");
  mg->SetTitle("Cut acceptance curves");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetXaxis()->SetTitle("Opening angle/#pi");
  mg->GetYaxis()->SetTitle("#epsilon_{s} /#sqrt{#epsilon_{b}}");

  // maximum marker
	TLine line(gr[0]->GetX()[maxi],0.,gr[0]->GetX()[maxi],mg->GetYaxis()->GetXmax());
  line.SetLineColor(2);
  line.Draw();

/*
  TMarker m1(gr[0]->GetX()[maxi],gr[0]->GetY()[maxi],2);
  m1.SetMarkerColor(2); // red
  m1.SetMarkerSize(2);
  m1.Draw();
*/

  // second axis
  TGaxis *axis2 = new TGaxis(1.,0.,1.,mg->GetYaxis()->GetXmax(),0.,1.,510,"+L");
  axis2->SetTitle("Acceptance, #epsilon");
  axis2->SetLabelFont(mg->GetYaxis()->GetLabelFont());
  axis2->SetTitleFont(mg->GetYaxis()->GetTitleFont());
  axis2->SetLabelSize(mg->GetYaxis()->GetLabelSize());
  axis2->SetTitleSize(mg->GetYaxis()->GetTitleSize());
  axis2->Draw();

  // Scale integrals for drawing
  for (unsigned i=0;i<n;i++) {
    y[0][i] *= mg->GetYaxis()->GetXmax()/y[0][n-1];
    y[1][i] *= mg->GetYaxis()->GetXmax()/y[1][n-1];
  }
  // ************************

  leg.SetFillColor(0);
  leg.Draw();

  TPaveText text(0.15,0.875,0.35,0.52,"NDCbr");
  text.SetBorderSize(1);
  text.SetFillColor(0);
  text.SetTextFont(42);
  text.AddText("C: Fraction of h = 0");
  text.AddText("");
  text.AddText("#epsilon_{h}(#theta) = #int_{0}^{#theta} P_{h}(#theta) #scale[3]{/}"
    "#int_{0}^{#pi} P_{h}(#theta)");
  text.AddText("");
  text.AddText("");
  text.AddText(Form("#gamma = %.3f",gamma));
  text.AddText(Form("C = %.3f",C));
  text.AddText(Form("#theta_{cut} = %.3f #pi",gr[0]->GetX()[maxi]));
  text.AddText(Form("h=0 accepted: %.0f%%",y[0][maxi]*100./y[0][n-1]));
  text.AddText(Form("h=1 accepted: %.0f%%",y[1][maxi]*100./y[1][n-1]));
  text.AddText(Form("h=1 fraction: %.0f%%",y[1][maxi]*100./(y[0][maxi]*y[0][n-1]/y[1][n-1]+y[1][maxi])));
  text.AddLine(.0,.55,1.,.55);
  text.Draw();

  canv.SaveAs("cut_acceptance.pdf");

  delete mg;

}

void draw_tilt_angle_deriv()
{
  canv.Clear();
  TMultiGraph *mg = new TMultiGraph();
  const unsigned ngr = 5;
  double gamma[ngr] = {1.2,2.,5.,10.,100};
  TGraph* gr[ngr];
  TLegend leg(0.91,0.7,0.99,0.9);

  for (unsigned i=0;i<ngr;i++) {
    gr[i] = tilt_angle_deriv_graph(gamma[i]);
    gr[i]->SetLineColor(colors[i+1]);
    gr[i]->SetMarkerColor(colors[i+1]);
    gr[i]->SetFillColor(0);
    gr[i]->SetLineWidth(2);
    mg->Add(gr[i]);
    TString name("#gamma = ");
    name += gamma[i];
    leg.AddEntry(gr[i],name);
  }

  mg->Draw("ac");
  mg->SetTitle("Tilt angle derivative");
  gPad->Update();

  mg->GetXaxis()->SetRangeUser(0.,1.);
  mg->GetYaxis()->SetRangeUser(-2,0);
  mg->GetXaxis()->SetTitle("Opening angle/#pi");
  mg->GetYaxis()->SetTitle("d#theta_{tilt}/d#theta_{op}");

  leg.SetFillColor(0);
  leg.Draw();

  canv.SaveAs("tilt_angle_deriv.pdf");

  delete mg;
}

// MAIN *************************************************************
int main(int argc, char* argv[])
{
  if ( argc>1 && (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-h")) ) {
    printf("Usage: %s [option]\n"
      "  --opening\n\t"
        "Draw opening angle as a function of tilt angle\n"
      "  --opening-min\n\t"
        "Draw minimum opening angle as a function of gamma\n"
      "  --tilt-prob\n\t"
        "Draw tilt angle probability distribution\n"
      "  --tilt-cos-prob\n\t"
        "Draw cos(tilt angle) probability distribution\n"
      "  --tilt-abs-cos-prob\n\t"
        "Draw |cos(tilt angle)| probability distribution\n"
      "  --tilt\n\t"
        "Draw tilt angle derivative as a function of opening angle\n"
      "  --tilt-deriv\n\t"
        "Draw tilt angle as a function of opening angle\n"
      "  --opening-prob [gamma=2]\n\t"
        "Draw opening angle probability distribution\n"
      "  --opening-prob-int [gamma=2]\n\t"
        "Draw opening angle probability distribution integral\n"
      "  --cut-accep [gamma=2] [longitutinal fraction=0.5]\n\t"
        "Draw acceptance curves\n",
      argv[0]
    );
  } else {
    bool draw_opening_angle_check = false;
    bool draw_opening_angle_min_check = false;
    bool draw_tilt_angle_prob_check = false;
    bool draw_tilt_angle_cos_prob_check = false;
    bool draw_tilt_angle_abs_cos_prob_check = false;
    bool draw_tilt_angle_check = false;
    bool draw_tilt_angle_deriv_check = false;
    bool draw_opening_angle_prob_check = false;
    bool draw_opening_angle_prob_int_check = false;
    bool draw_cut_acceptance_check = false;

    // defaults
    double gamma = 2;
    double Cfraction = 0.5;

    if (argc==1 || (argc==2 && !strcmp(argv[1],"--all"))) {
      draw_opening_angle_check = true;
      draw_opening_angle_min_check = true;
      draw_tilt_angle_prob_check = true;
      draw_tilt_angle_cos_prob_check = true;
      draw_tilt_angle_abs_cos_prob_check = true;
      draw_tilt_angle_check = true;
      draw_tilt_angle_deriv_check = true;
      draw_opening_angle_prob_check = true;
      draw_opening_angle_prob_int_check = true;
    } else {
      if (!strcmp(argv[1],"--opening")) draw_opening_angle_check = true;
      if (!strcmp(argv[1],"--opening-min")) draw_opening_angle_min_check = true;
      if (!strcmp(argv[1],"--tilt-prob")) draw_tilt_angle_prob_check = true;
      if (!strcmp(argv[1],"--tilt-cos-prob")) draw_tilt_angle_cos_prob_check = true;
      if (!strcmp(argv[1],"--tilt-abs-cos-prob")) draw_tilt_angle_abs_cos_prob_check = true;
      if (!strcmp(argv[1],"--tilt")) draw_tilt_angle_check = true;
      if (!strcmp(argv[1],"--tilt-deriv")) draw_tilt_angle_deriv_check = true;
      if (!strcmp(argv[1],"--opening-prob")) {
        draw_opening_angle_prob_check = true;
        if (argc>2) gamma = atof(argv[2]);
      }
      if (!strcmp(argv[1],"--opening-prob-int")) {
        draw_opening_angle_prob_int_check = true;
        if (argc>2) gamma = atof(argv[2]);
      }
      if (!strcmp(argv[1],"--cut-accep")) {
        draw_cut_acceptance_check = true;
        if (argc>2) gamma = atof(argv[2]);
        if (argc>3) Cfraction = atof(argv[3]);
      }
    }

    if (draw_opening_angle_check) draw_opening_angle();
    if (draw_opening_angle_min_check) draw_min_opening_angle();
    if (draw_tilt_angle_check) draw_tilt_angle();
    if (draw_tilt_angle_deriv_check) draw_tilt_angle_deriv();
    if (draw_tilt_angle_prob_check) draw_tilt_angle_prob();
    if (draw_tilt_angle_cos_prob_check) draw_tilt_angle_cos_prob();
    if (draw_tilt_angle_abs_cos_prob_check) draw_tilt_angle_abs_cos_prob();
    if (draw_opening_angle_prob_check) draw_opening_angle_prob(gamma);
    if (draw_opening_angle_prob_int_check) draw_opening_angle_prob_int(gamma);
    if (draw_cut_acceptance_check) draw_cut_acceptance(gamma,Cfraction);
  }
  return 0;
}
