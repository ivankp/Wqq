#include <cmath>
#include <cstdio>

#include "TMarker.h"
#include "TString.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TF1.h"

#include "inverse.h"

#include "functions.h"

using namespace std;

int colors[] = {51,56,65,78,91,99}; // color scheme
TCanvas canv("canv","",600,400); // Default canvas

// Plots ************************************************************
double cut_angle2(double gamma, bool draw=false, bool units_of_pi=false)
{
  using namespace std::placeholders;
  function<double (double)> f = bind(opening_angle,_1,gamma);
  inverse g(f,100,0,M_PI/2.); // turns out doesn't have to be so precise

  const unsigned n = 200;

  double x[n], y[2][n];
  double step = M_PI/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = step*i;
    y[0][i] = opening_angle_prob(x[i],0,g);
    y[1][i] = opening_angle_prob(x[i],1,g);
  }

  // ************************
  // Integrate and compute **
  for (unsigned i=n-1;i!=0;i--) {
    for (unsigned j=0;j<i;j++)
      for (unsigned k=0;k<2;k++) y[k][i] += y[k][j];
    if (y[1][i]>0.) y[0][i] /= sqrt(y[1][i]);
  }
  // ************************
  // ************************

  // find maximum ***********
  // ************************
  unsigned maxi = 2./(gamma*step); // start from derived minimum angle
  while (y[0][maxi]<y[0][n-1]) maxi++;
  while (y[0][maxi+1] > y[0][maxi]) maxi++;
  // ************************
  // ************************

  if (units_of_pi) // in units of pi
    for (unsigned i=0;i<n;i++) x[i] /= M_PI;

  // Draw *******************
  if (draw) {

    canv.Clear();
    TGraph *gr = new TGraph(n,x,y[0]);
    TLegend leg(0.91,0.7,0.99,0.9);
    leg.SetHeader(Form("#gamma = %.1f",gamma));

    gr->SetLineColor(colors[3]);
    gr->SetMarkerColor(colors[3]);
    gr->SetFillColor(0);
    gr->SetLineWidth(2);
    leg.AddEntry(gr,"#epsilon_{0} /#sqrt{#epsilon_{1}}");

    gr->Draw("ac");
    gr->SetTitle("Acceptance curve");
    gPad->Update();

    gr->GetXaxis()->SetRangeUser(0.,M_PI);
    gr->GetXaxis()->SetTitle("Opening angle");
    gr->GetYaxis()->SetTitle("Acceptance");

    leg.SetFillColor(0);
    leg.Draw();

    TMarker m1(x[maxi],y[0][maxi],2);
    m1.SetMarkerColor(4); // blue
    m1.SetMarkerSize(2);
    m1.Draw();

    canv.SaveAs("accep.pdf");

    delete gr;
  }

  return x[maxi];
}

void cut_vs_gamma_plot()
{
  const unsigned n = 100;

  double x[n], y[n];
  double g1 = 1.01, g2 = 50.;
  double step = (g2-g1)/(n-1);
  for (unsigned i=0;i<n;i++) {
    x[i] = g1+step*i;
    y[i] = cut_angle2(x[i]);
    //printf("%3u %7.3f %7.3f\n",i,x[i],y[i]);
  }

  canv.Clear();
  canv.SetLogy();
  TGraph *gr = new TGraph(n,x,y);

  //TF1 *fcn = new TF1("fitfcn","[0]*exp(-[2]*x)+[1]*exp(-[3]*x)",g1,g2);
  //gr->Fit("fitfcn");

  gr->SetLineColor(96);
  gr->SetMarkerColor(96);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.5);
  gr->SetFillColor(0);
  gr->SetLineWidth(2);

  gr->Draw("ap");
  gr->SetTitle("Cut dependence on #gamma");
  gPad->Update();

  gr->GetXaxis()->SetRangeUser(1.,g2);
  gr->GetXaxis()->SetTitle("#gamma");
  gr->GetYaxis()->SetTitle("Opening angle cut");

  canv.SaveAs("cut_gamma.pdf");

  delete gr;
  //delete fcn;
}



// MAIN *************************************************************
int main(int argc, char* argv[])
{
  if (argc==2) {
    printf("%.5f\n",cut_angle2(atof(argv[1]),true));
  } else {
    cut_vs_gamma_plot();
    
  }

  return 0;
}
