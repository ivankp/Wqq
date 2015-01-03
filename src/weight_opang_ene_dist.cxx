#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <boost/program_options.hpp>

#include "TStyle.h"
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

#include "LHEF2Wqq.h"

#include "inverse.h"
#include "functions.h"

using namespace std;
namespace po = boost::program_options;

TCanvas canv("canv","",600,400); // Default canvas
string base_name;

class dist_base {

public:
  dist_base(float Emin, float Emax, unsigned nvals, unsigned npoints)
  : Emin(Emin), Emax(Emax), nvals(nvals), npoints(npoints)
  {
    dist = new double*[nvals];
    for (unsigned i=0;i<nvals;i++)
      dist[i] = new double[npoints];
    x   = new double[npoints];
    tmp = new double[npoints];

    normalized = false;
  }
  virtual ~dist_base() {
    for (unsigned i=0;i<nvals;i++) delete [] dist[i];
    delete [] dist;
    delete [] x;
    delete [] tmp;
  }

  float getEmin() const { return Emin; }
  float getEmax() const { return Emax; }
  bool inrange(float E) const { return (Emin<E && E<Emax); }

  void save_dat() {
    normalize();

    stringstream ss;
    ss << base_name << "_opang_ew_" << Emin << '-' << Emax << ".dat";
    cout << "Writing data file: " << ss.str() << endl;

    ofstream dat(ss.str());
    dat.precision(15);
    if (dat.is_open()) {
      for (unsigned i=0;i<npoints;i++) {
        dat << x[i];
        for (unsigned j=0;j<nvals;j++) dat << '\t' << dist[j][i];
        dat << endl;
      }
      dat.close();
    } else {
      cout << "Unable to open file " << ss.str() << endl;
      exit(1);
    }
  }

  virtual void save_plot() =0;

protected:
  const float Emin, Emax;
  const unsigned nvals, npoints;
  double* x;
  double** dist;
  double* tmp;
  bool normalized;

  void normalize() {
    if (!normalized) {
      for (unsigned j=0;j<nvals;j++) {
        double integral = 0.;
        for (unsigned i=0;i<npoints;i++) integral += dist[j][i];
        for (unsigned i=0;i<npoints;i++) dist[j][i] /= integral;
      }
    }
    normalized = true;
  }

};

class dist_opang: public dist_base {

public:
  dist_opang(float Emin, float Emax, unsigned npoints=1000)
  : dist_base(Emin,Emax,2,npoints)
  {
    const double step = M_PI/(npoints-1);
    for (unsigned i=0;i<npoints;i++) x[i] = step*i;
  }
  virtual ~dist_opang() { }

  virtual void add(int helicity, double gamma)
  {
    function<double (double)> f = bind(opening_angle,std::placeholders::_1,gamma);
    inverse g(f,200,0,M_PI/2.);

    double integral = 0.;
    for (unsigned i=0;i<npoints;i++) {
      tmp[i] = opening_angle_prob(x[i],helicity,g);
      integral += tmp[i];
    }
    for (unsigned i=0;i<npoints;i++)
      dist[helicity][i] += tmp[i]/integral;

    normalized = false;
  }

  virtual void save_plot() {
    normalize();

    canv.Clear();
    TMultiGraph *mg = new TMultiGraph();
    TGraph* gr[2];
    TLegend leg(0.91,0.7,0.99,0.9);
    int colors[] = {78,65};

    for (unsigned h=0;h<2;h++) {
      gr[h] = new TGraph(npoints,x,dist[h]);
      gr[h]->SetFillColor(0);
      gr[h]->SetLineWidth(2);
      gr[h]->SetLineColor(colors[h]);
      gr[h]->SetMarkerColor(colors[h]);
      mg->Add(gr[h]);
      leg.AddEntry(gr[h],Form("h = %d",h));
    }

    mg->Draw("ac");
    mg->SetTitle("Energy weighted Opening angle distribution");
    gPad->Update();

    mg->GetXaxis()->SetRangeUser(0.,M_PI);
    mg->SetMinimum(0.);
    mg->GetYaxis()->SetTitleOffset(1.4);
    mg->GetXaxis()->SetTitle("Opening angle");
    mg->GetYaxis()->SetTitle("Probability density");

    leg.SetFillColor(0);
    leg.Draw();

    stringstream ss;
    ss << base_name << "_opang_ew_" << Emin << '-' << Emax << ".pdf";

    canv.SaveAs(ss.str().c_str());

    delete mg;
  }

};

int main (int argc, char* argv[])
{
  // START OPTIONS **************************************************
  string input_file, data_dir, cut_str;
  vector<string> ene_str;
  bool cut_applied = false;
  unsigned npoints;
  bool make_plots;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("lhe,i", po::value<string>(&input_file),
       "input lhe file")
      ("dir,o", po::value<string>(&data_dir),
       "output data directory")
      ("energy,e", po::value< vector<string> >(&ene_str),
       "add energy range")
      ("cut,c", po::value<string>(&cut_str),
       "apply a data cut")
      ("num,n", po::value<unsigned>(&npoints)->default_value(1000),
       "set number of computation points")
      ("plot,p", po::bool_switch(&make_plots),
       "draw plots")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    po::notify(vm);

    // Options Properties ---------------------------------
    if (argc == 1 || vm.count("help")) {
      cout << all_opt << endl;
      exit(0);
    } else if (vm.count("cut")) {
      cut_applied = true;
    }

    // Necessary options ----------------------------------
    if (!vm.count("lhe"))
    { cerr << "no input lhe file provided" << endl; exit(1); }
    if (!vm.count("dir"))
    { cerr << "no output data provided" << endl; exit(1); }
    if (!vm.count("energy"))
    { cerr << "no energy ranges provided" << endl; exit(1); }
  }
  catch(exception& e) {
    cerr << "Error: " <<  e.what() << endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Book distributions ***********************************
  const size_t bn0 = input_file.rfind('/')+1;

  base_name = data_dir;
  base_name += '/';
  base_name += input_file.substr(bn0,input_file.size()-bn0-4);

  const unsigned ne = ene_str.size();
  dist_opang* dist[ne];

  for (unsigned e=0; e<ne; e++) { // energy ranges
    const size_t colon = ene_str[e].find(':');
    float E[2] = {
      (float)atof(ene_str[e].substr(0,colon).c_str()),
      (float)atof(ene_str[e].substr(colon+1).c_str())
    };
    if (E[0]>E[1]) swap(E[0],E[1]);
    dist[e] = new dist_opang(E[0],E[1]);
  }

  // Create LHEF Readed ***********************************
  LHEF2Wqq data(input_file,true);

  if (cut_applied) {
    data.set_cut(cut_str);
    cout << "Applying cut: " << cut_str << endl;
  }

  while ( data.readEvent() ) { // loop over events

    for (unsigned i=0; i<data.size(); i++) {
      int pid = data[i]->id();
      if (pid==24 || pid==-24)
        for (unsigned e=0;e<ne;e++)
          if (dist[e]->inrange(data[i]->energy())) {
            for (int h=0;h<2;h++) dist[e]->add(h,data[i]->gamma());
            break;
          }
    }

  } // end while
  data.done_msg();

  // Save output ******************************************
  ofstream dat;
  for (unsigned e=0; e<ne; e++) {
    // write data file
    dist[e]->save_dat();

    // make plot
    dist[e]->save_plot();

    delete dist[e];
  }

  return 0;
}
