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

template <class T>
struct two {
  T x1, x2;
  two(const T& x1, const T& x2) : x1(x1), x2(x2) {}
  //two(const T& x) : x1(x), x2(x) {}
  T& operator[](unsigned i) {
    if (i) return x2;
    else   return x1;
  }
};

TCanvas canv("canv","",900,600); // Default canvas

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

  void save_dat(const string& prefix) {
    normalize();

    stringstream ss;
    ss << prefix << Emin << '-' << Emax << ".dat";
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

  virtual void save_plot(const string& prefix) =0;

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

bool div_by_min;
double (*opening_angle_fcn)(double,double); // <-- work here

class dist_opang: public dist_base {

public:
  dist_opang(float Emin, float Emax, unsigned npoints=1000)
  : dist_base(Emin,Emax,2,npoints)
  {
    double step, min;
    if (div_by_min) {
      min  = 1.;
      step = (10.-min)/(npoints-1);
    } else {
      min  = 0.;
      step = M_PI/(npoints-1);
    }
    for (unsigned i=0;i<npoints;i++) x[i] = step*i+min;
  }
  virtual ~dist_opang() { }

  virtual void add(int helicity, double gamma, bool flip=false)
  {
    const ddfcn f = bind(opening_angle_fcn,std::placeholders::_1,gamma);

    const inverse* tilt_from_op;
    if (!flip) tilt_from_op = new inverse(f,200,0.,M_PI/2.);
    else       tilt_from_op = new inverse(f,200,M_PI/2.,M_PI);

    double integral = 0.;
    for (unsigned i=0;i<npoints;i++) {
      tmp[i] = opening_angle_prob(x[i],helicity,*tilt_from_op);
      integral += tmp[i];
    }
    for (unsigned i=0;i<npoints;i++)
      dist[helicity][i] += tmp[i]/integral;

    normalized = false;

    delete tilt_from_op;
  }

  virtual void save_plot(const string& prefix) {
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
      gr[h]->SetMarkerStyle(20);
      gr[h]->SetMarkerSize(0.3);
      mg->Add(gr[h]);
      leg.AddEntry(gr[h],Form("h = %d",h));
    }

    mg->Draw("ap");
    mg->SetTitle("Energy weighted opening angle distribution");
    gPad->Update();

    mg->SetMinimum(0.);
    mg->GetYaxis()->SetTitleOffset(1.4);
    if (div_by_min) {
      mg->GetXaxis()->SetRangeUser(1.,10.);
      mg->GetXaxis()->SetTitle("Opening angle/min");
    } else {
      mg->GetXaxis()->SetRangeUser(0.,M_PI);
      mg->GetXaxis()->SetTitle("Opening angle, rad");
    }
    mg->GetYaxis()->SetTitle("Probability density");

    leg.SetFillColor(0);
    leg.Draw();

    stringstream ss;
    ss << prefix << Emin << '-' << Emax << ".png";

    canv.SaveAs(ss.str().c_str());

    delete mg;
  }

};

int main (int argc, char* argv[])
{
  // START OPTIONS **************************************************
  string input_file, data_dir;
  string bin_str, var_bin_str;
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

      ("bins,b", po::value<string>(&bin_str),
       "min:max:step")
      ("var-bin,v", po::value<string>(&var_bin_str),
       "binning variable")
      ("div_by_min", po::bool_switch(&div_by_min),
       "divide opening angle by minimum possible")

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
    }

    // Necessary options ----------------------------------
    if (!vm.count("lhe"))
    { cerr << "no input lhe file provided" << endl; exit(1); }
    if (!vm.count("dir"))
    { cerr << "no output data provided" << endl; exit(1); }
    if (!vm.count("bins"))
    { cerr << "no binning ranges provided" << endl; exit(1); }
    if (!vm.count("var-bin"))
    { cerr << "no binning variable provided" << endl; exit(1); }
  }
  catch(exception& e) {
    cerr << "Error: " <<  e.what() << endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Binning variable *************************************
  double (particle::* fcn_bin)() const;

  if      (!var_bin_str.compare("energy-lab")) {
    fcn_bin = &particle::energy;
  }
  else if (!var_bin_str.compare("energy-cm")) {
    fcn_bin = &particle::energy_cm;
  }
  else {
    cerr << "Undefined binning variable: " << var_bin_str << endl;
    exit(1);
  }

  if (div_by_min) opening_angle_fcn = &opening_angle_over_min;
  else            opening_angle_fcn = &opening_angle;

  // Book distributions ***********************************
  vector< two<dist_opang*> > dist;
  double bin_min=0., bin_max=0., bin_step=0.;

  // find colons (:)
  size_t colon[2];
  unsigned ncol = 0;
  for (size_t i=0;i<bin_str.size();i++) {
    if (bin_str[i]==':') {
      if (ncol>1) {
        cerr << "Badly formated bins argument " << bin_str << endl;
        exit(1);
      } else {
        colon[ncol] = i; ncol++;
      }
    }
  }

  bin_min  = atof(bin_str.substr(0,colon[0]).c_str());
  bin_max  = atof(bin_str.substr(colon[0]+1,colon[1]-colon[0]-1).c_str());
  bin_step = atof(bin_str.substr(colon[1]+1).c_str());

  if ( bin_min>bin_max || bin_max-bin_min<bin_step ) {
    cerr << "Badly formated bins argument " << bin_str << endl;
    exit(1);
  }

  for (double b=bin_min;b+bin_step<=bin_max;b+=bin_step)
    dist.push_back( two<dist_opang*>( new dist_opang(b,b+bin_step), new dist_opang(b,b+bin_step) ) );

  // Create LHEF Readed ***********************************
  LHEF2Wqq data(input_file,true);

  while ( data.readEvent() ) { // loop over events

    const particle* W[2] = { data.Wp(), data.Wm() };

    for (unsigned char i=0;i<2;i++)
      for (auto d=dist.begin(); d!=dist.end(); ++d)
        if ((*d)[i]->inrange((W[i]->*fcn_bin)())) {

          //bool flip = bool( W[i]->tilt_angle() > (M_PI/2.) ); // CHEAT with trueth info

          for (int h=0;h<2;h++) (*d)[i]->add(h,W[i]->gamma());

          break;
        }

  } // end while
  data.done_msg();

  // Save output ******************************************
  ofstream dat;
  for (auto d=dist.begin(); d!=dist.end(); ++d)
    for (unsigned char i=0;i<2;i++) {
      stringstream ss;
      ss << data_dir;
      if (*(data_dir.end()-1) != '/') ss << '/';
      if (i) ss << "Wm";
      else   ss << "Wp";

      (*d)[i]->save_dat(ss.str());
      if (make_plots) (*d)[i]->save_plot(ss.str());

      delete (*d)[i];
    }

  return 0;
}
