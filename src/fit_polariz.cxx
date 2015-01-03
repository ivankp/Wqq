#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <utility>

#include <boost/program_options.hpp>

#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TText.h"
#include "TPaveText.h"

#include "LHEF2Wqq.h"

#include "chisq_polariz_op.h"
#include "chisq_polariz_tilt.h"
#include "chisq_polariz_ediff.h"
#include "chisq_polariz_ediff_abs.h"

using namespace std;
namespace po = boost::program_options;

enum class fittype : char { op, tilt, ediff, ediff_abs, op_norm } fit_type;
enum class bintype : char { none, energy_lab, energy_cm, cos_prod } bin_type;

TCanvas canv("canv","",600,400); // Default canvas

class ftos {
private:
  char* str;
public:
  ftos (const double& f) {
    double fracpart, intpart;
    fracpart = modf(f,&intpart);

    const unsigned lenint  = log10(abs(intpart )) + 1;
          unsigned lenfrac = abs(log10(fracpart)) + 1;
    if (lenfrac>4) lenfrac=0;
    const unsigned len = lenint+lenfrac+3; // minus, period, null
    char fmt[5];
    str = new char[len];

    sprintf(fmt,"%%.%df",lenfrac);
    sprintf(str,fmt,f);
  }
  ~ftos() { delete[] str; }
  const char* c_str() const { return str; }
};

// Histogram wrapper ************************************************

class hist {
public:
  static int nbins;
  virtual const TH1F* operator[](unsigned i) const {
    if (i==0) return &p;
    else      return &m;
  }
  virtual ~hist() { }
protected:
  hist(const char* name, const char* title, const char* xtitle, double min, double max)
  : p(TH1F(Form("%s+",name),title,nbins,min,max)),
    m(TH1F(Form("%s-",name),title,nbins,min,max))
  {
    p.SetXTitle(xtitle);
  }
  TH1F p, m;
};
int hist::nbins=1;

class hist_all : public hist {
public:
  hist_all(const char* title, const char* xtitle, double min, double max)
  : hist("hist_all",Form("%s   All",title),xtitle,min,max)
  { }
  virtual void Fill(const double& pval, const double& mval) {
    p.Fill(pval);
    m.Fill(mval);
  }
};
class hist_bin : public hist {
public:
  const double b1, b2;
  hist_bin(const char* title, const char* xtitle, double min, double max, double b1, double b2)
  : hist(Form("hist_%f:%f",b1,b2),
         Form("%s   Bin [%s,%s]",title,ftos(b1).c_str(),ftos(b2).c_str()),xtitle,min,max),
    b1(b1), b2(b2)
  { }
  virtual void Fill(const double& pval, const double& mval,
                    const double& pbin, const double& mbin) {
    if (b1<pbin && pbin<b2) p.Fill(pval);
    if (b1<mbin && mbin<b2) m.Fill(mval);
  }
};

// MAIN *************************************************************
int main(int argc, char* argv[])
{
  // START OPTIONS **************************************************
  string input_file, output_plot, fit_type_str, bin_str, ew_dir;
  unsigned nbinavg;
  bool draw_chisq;
  bin_type = bintype::none;

  try {
    // General Options ------------------------------------
    po::options_description all_opt("Options");
    all_opt.add_options()
      ("help,h", "produce help message")
      ("lhe,i", po::value<string>(&input_file),
       "input lhe file")
      ("output,o", po::value<string>(&output_plot)->
       default_value("fit.pdf"),
       "output plot file")

      ("fit,f", po::value<string>(&fit_type_str),
       "type of fit")

      ("bins-energy-lab", po::value<string>(&bin_str),
       "energy bins min:max:step GeV")
      ("bins-energy-cm", po::value<string>(&bin_str),
       "energy bins min:max:step GeV")
      ("bins-cos-prod", po::value<string>(&bin_str),
       "|cos(θprod)| bins min:max:step")

      ("nbins,n", po::value<int>(&hist::nbins)->default_value(100),
       "number of bins for data histograms")
      ("nbinavg,a", po::value<unsigned>(&nbinavg)->default_value(25),
       "number of point for averaging theoretical distr within data bin")
      ("weighted,d", po::value<string>(&ew_dir),
       "directory with the energy weighted distributions")

      ("draw-chisq", po::bool_switch(&draw_chisq),
       "draw the χ² distribution after fit")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, all_opt), vm);
    po::notify(vm);

    // Options Properties ---------------------------------
    if (argc==1 || vm.count("help")) {
      cout << all_opt << endl;
      exit(0);
    }

    int check_num_bins_args = 0;
    if (vm.count("bins-energy-lab")) { bin_type = bintype::energy_lab; check_num_bins_args++; }
    if (vm.count("bins-energy-cm"))  { bin_type = bintype::energy_cm;  check_num_bins_args++; }
    if (vm.count("bins-cos-prod"))   { bin_type = bintype::cos_prod;   check_num_bins_args++; }

    if (check_num_bins_args>1)
    { cerr << "Only one type of binning can be specified" << endl;
      exit(1); }

    // Necessary options ----------------------------------
    if (!vm.count("lhe"))
    { cerr << "No input lhe file provided." << endl; exit(1); }
    if (!vm.count("fit"))
    { cerr << "No fit type provided" << endl;
      cerr << "Options are: op, tilt, ediff." << endl;
      exit(1); }

    if (!fit_type_str.compare("op")) {
      if (!vm.count("weighted")) {
        cerr << "No directory with the energy weighted distributions provided." << endl;
        exit(1);
      }
      if (!check_num_bins_args) {
        cerr << "No binning ranges provided." << endl;
        exit(1);
      }
    }

  } catch(exception& e) {
    cerr << "Arguments error: " <<  e.what() << endl;
    exit(1);
  }
  // END OPTIONS ****************************************************

  // Select fit type **************************************
  double (particle::* fcn_fill)() const;

  if      (!fit_type_str.compare("op")) {
    fit_type = fittype::op;
    fcn_fill = &particle::opening_angle;
  }
  else if (!fit_type_str.compare("tilt")) {
    fit_type = fittype::tilt;
    fcn_fill = &particle::tilt_angle;
  }
  else if (!fit_type_str.compare("ediff")) {
    fit_type = fittype::ediff;
    fcn_fill = &particle::ediff;
  }
  else if (!fit_type_str.compare("ediff_abs")) {
    fit_type = fittype::ediff_abs;
    fcn_fill = &particle::ediff_abs;
  }
  else if (!fit_type_str.compare("op_norm")) {
    fit_type = fittype::op_norm;
    fcn_fill = &particle::op_over_min;
  }
  else {
    cerr << "Undefined fit type: " << fit_type_str << endl;
    exit(1);
  }

  // Select bin type **************************************
  double (particle::* fcn_bin)() const;

  switch (bin_type) {
    case bintype::energy_lab: fcn_bin = &particle::energy;       break;
    case bintype::energy_cm:  fcn_bin = &particle::energy_cm;    break;
    case bintype::cos_prod:   fcn_bin = &particle::abs_cos_prod; break;
    default: fcn_bin = nullptr; break;
  }

  // Created histograms ***********************************
  vector<hist*> hdata;
  double bin_min=0., bin_max=0., bin_step=0.;

  if (bin_type != bintype::none) {
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

    // Book histograms
    for (double b=bin_min;b+bin_step<=bin_max;b+=bin_step) {
      switch (fit_type) {
        case fittype::op:  hdata.push_back(
          new hist_bin("W #rightarrow qq  opening angle","Opening angle, rad",0.,M_PI,b,b+bin_step)
        ); break;
        case fittype::tilt:  hdata.push_back(
          new hist_bin("W #rightarrow qq  tilt angle","Tilt angle, rad",0.,M_PI,b,b+bin_step)
        ); break;
        case fittype::ediff: hdata.push_back(
          new hist_bin("W #rightarrow qq  #DeltaE/p_{W}","#DeltaE/p_{W}",-1.,1.,b,b+bin_step)
        ); break;
        case fittype::ediff_abs: hdata.push_back(
          new hist_bin("W #rightarrow qq  |#DeltaE/p_{W}|","|#DeltaE/p_{W}|",0.,1.,b,b+bin_step)
        ); break;
        case fittype::op_norm:  hdata.push_back(
          new hist_bin("W #rightarrow qq  opening angle","Opening angle/min",1.,4.,b,b+bin_step)
        ); break;
      }
    }
  } else {
    switch (fit_type) {
      case fittype::tilt:  hdata.push_back(
        new hist_all("W #rightarrow qq  tilt angle","Tilt angle, rad",0.,M_PI)
      ); break;
      case fittype::ediff: hdata.push_back(
        new hist_all("W #rightarrow qq  #DeltaE/p_{W}","#DeltaE/p_{W}",-1.,1.)
      ); break;
      case fittype::ediff_abs: hdata.push_back(
        new hist_all("W #rightarrow qq  |#DeltaE/p_{W}|","|#DeltaE/p_{W}|",0.,1.)
      ); break;
      default:
        cerr << "This fit type can only be done for a specific bin" << endl;
        exit(1);
    }
  }

  // Create LHEF Readed ***********************************
  LHEF2Wqq data(input_file,true);

  if (bin_type == bintype::none) {
    while ( data.readEvent() ) { // loop over events

        ((hist_all*)*hdata.begin())->Fill(
          (data.Wp()->*fcn_fill)(), (data.Wm()->*fcn_fill)()
        );

    } // end while
  } else {
    while ( data.readEvent() ) { // loop over events

      for (auto it=hdata.begin();it!=hdata.end();++it)
        ((hist_bin*)*it)->Fill(
          (data.Wp()->*fcn_fill)(), (data.Wm()->*fcn_fill)(),
          (data.Wp()->*fcn_bin )(), (data.Wm()->*fcn_bin )()
        );

    } // end while
  }
  data.done_msg();

  // FIT **************************************************
  vector<chisq_polariz_base**> fits;
  for (auto h=hdata.begin();h!=hdata.end();++h) {
    chisq_polariz_base** fit = new chisq_polariz_base*[2];

    for (unsigned i=0;i<2;i++) {
      stringstream ss;
      switch (fit_type) {
        case fittype::op_norm:
        case fittype::op:
          ss << ew_dir;
          if ( *(ew_dir.end()-1)!='/' ) ss << '/';
          if (i) ss << "Wm"; else   ss << "Wp";
          ss << ((hist_bin*)*h)->b1 << '-' << ((hist_bin*)*h)->b2 << ".dat";
          fit[i] = new chisq_polariz_op( (**h)[i], ss.str().c_str() );
          break;
        case fittype::tilt:  fit[i] = new chisq_polariz_tilt ((**h)[i]); break;
        case fittype::ediff: fit[i] = new chisq_polariz_ediff((**h)[i]); break;
        case fittype::ediff_abs: fit[i] = new chisq_polariz_ediff_abs((**h)[i]); break;
      }
      fit[i]->set_nbinavg(nbinavg);
      fit[i]->Fit();
    }

    fits.push_back(fit);
  }

  // Draw *************************************************
  canv.SaveAs(Form("%s[",output_plot.c_str()));

  // Draw summary
  TH1F*** summary = nullptr;
  unsigned npolfrac = 0; // number of polarization fractions
  if (bin_type != bintype::none) {

    switch (fit_type) {
      case fittype::op:        npolfrac=2; break;
      case fittype::tilt:      npolfrac=3; break;
      case fittype::ediff:     npolfrac=3; break;
      case fittype::ediff_abs: npolfrac=2; break;
      case fittype::op_norm:   npolfrac=2; break;
    }

    summary = new TH1F**[2];
    for (unsigned i=0;i<2;i++) {    // W+, W-
      summary[i] = new TH1F*[npolfrac];
      for (unsigned j=0;j<npolfrac;j++) {  // helicity // 0:0  1:+  2:-
        // Book histogram
        summary[i][j] = new TH1F(Form("summary_%u_%u",i,j),"",
          hdata.size(),bin_min,bin_min+bin_step*hdata.size());

        // Fill histogram
        for (unsigned k=0;k<fits.size();k++) {    // fit
          summary[i][j]->SetBinContent(k+1,fits[k][i]->get_frac(j));
          summary[i][j]->SetBinError  (k+1,fits[k][i]->get_unc (j));
        }
      }
    }

    const int* const color = chisq_polariz_base::colors;
    unsigned jj[3] = {1,0,2};
    for (unsigned j, j2=0; j2<npolfrac; j2++) {  // helicity // 1:+  0:0  2:-
      j = jj[j2];
      canv.Clear();

      for (unsigned i=0;i<2;i++) {
        TH1F* const h = summary[i][j];
        h->SetStats(0);
        h->SetLineColor(color[i]);
        h->SetMarkerColor(color[i]);
        if (i) h->Draw("same");
        else   {
          h->SetMinimum(0.);
          h->SetMaximum(1.);
          const char* title = "dependence of polarization fraction";
          if (bin_type==bintype::energy_lab) {
            title = Form("Lab Energy %s",title);
            h->SetXTitle("W energy, GeV");
          } else if (bin_type==bintype::energy_cm) {
            title = Form("CM Energy %s",title);
            h->SetXTitle("W energy, GeV");
          } else if (bin_type==bintype::cos_prod) {
            title = Form("|cos(#theta_{prod})| %s",title);
            h->SetXTitle("|cos(#theta_{prod})|");
          }
          if (j==0) h->SetTitle(Form("%s h_{0}",title));
          else if (j==1) {
            if      (npolfrac==3) h->SetTitle(Form("%s h_{+}",  title));
            else if (npolfrac==2) h->SetTitle(Form("%s h_{#pm}",title));
          }
          else if (j==2) h->SetTitle(Form("%s h_{#font[122]{-}}",title));
          h->Draw();
        }
      }
      canv.SaveAs(output_plot.c_str());
    }

  }

  if (fit_type==fittype::op_norm //||
      //fit_type==fittype::op
  ) canv.SetLogy();

  // Draw fits
  for (auto fit=fits.begin();fit!=fits.end();++fit) {
    canv.Clear();

    chisq_polariz_base::Draw((*fit)[0],(*fit)[1]);

    canv.SaveAs(output_plot.c_str());

    if (draw_chisq) {
      chisq_polariz_base::DrawChiSq((*fit)[0],(*fit)[1]);
      canv.SaveAs(output_plot.c_str());
    }
  }

  canv.SaveAs(Form("%s]",output_plot.c_str()));

  // Collect garbage **************************************
  for (auto h=hdata.begin();h!=hdata.end();++h) delete *h;
  for (auto fit=fits.begin();fit!=fits.end();++fit) {
    for (unsigned i=0;i<2;i++) delete (*fit)[i];
    delete[] *fit;
  }
  if (summary) {
    for (unsigned i=0;i<2;i++) {
      for (unsigned j=0;j<npolfrac;j++) delete summary[i][j];
      delete[] summary[i];
    }
    delete[] summary;
  }

  return 0;
}
