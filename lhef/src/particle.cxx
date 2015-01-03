#include "particle.h"

#include <cmath>
#include <cstring>
#include <exception>

#include <iostream>
#include <iomanip>

using namespace std;

ostream& operator<< (ostream& out, const TLorentzVector& x) {
  out << '('
      << setw(7) << setprecision(4) << x.Px() << ", "
      << setw(7) << setprecision(4) << x.Py() << ", "
      << setw(7) << setprecision(4) << x.Pz() << ", "
      << setw(7) << setprecision(4) << x.E () << ')';
  return out;
}

// ******************************************************************

class motherex: public exception {
public:
  motherex(unsigned i) {
    sprintf(str,"Invalid particle mother index %u",i);
  }
  virtual const char* what() const throw() {
    return str;
  }
private:
  char str[40];
};

// ******************************************************************
// constructor
particle::particle(int num, int id, int stat, const vector<double>& pup) {
  n = num;
  pid = id;
  status = stat;
  P = TLorentzVector(pup[0],pup[1],pup[2],pup[3]);
  m = pup[4];
  mother[0] = 0;
  mother[1] = 0;
}
particle::~particle() {  } // destructor

// ******************************************************************

int particle::id() const { return pid; }
int particle::get_status() const { return status; }

void particle::add(const particle* p) { // add product
  products.push_back(p);
}
void particle::set_mother(unsigned i,const particle* p) {
  if (i==0 || i==1) mother[i] = p;
  else throw motherex(i);
}
const particle* particle::get_mother(unsigned i) const {
  if (i==0 || i==1) return mother[i];
  else throw motherex(i);
}

// ******************************************************************

unsigned particle::nprod() const { // get number of products
  return products.size();
}

// get product from particle
const particle* particle::prod(unsigned i) const {
  return products[i];
}
// get product from particle
const particle* particle::operator[](unsigned i) const {
  return products[i];
}

// ******************************************************************

double particle::mass() const { // get particle mass
  return m;
}
double particle::mass2() const { // compute particle mass
  return pow(energy(),2) - pow(momvecmag(),2);
}
double particle::energy() const { // get particle energy
  return P.E();
}
double particle::energy_cm() const { // get particle energy
  if (!( mother[0]==0 || mother[1]==0 )) {

    TLorentzVector boost = mother[0]->P + mother[1]->P;
    boost *= (-1./boost.E());

    TLorentzVector this_particle(P);
    this_particle.Boost(boost.Px(),boost.Py(),boost.Pz());

    return this_particle.E(); // energy in CM frame

  } else return 0.;
}
double particle::momvecmag() const { // get particle energy
  return P.Vect().Mag();
}
double particle::gamma() const { // get particle Lorentz factor
  return P.E()/m;
}
double particle::Pt() const { // get particle transverse momentum
  return P.Pt();
}

 // get the angle between the products' momenta
double particle::opening_angle() const {
  if (nprod()==2)
    return products[0]->P.Angle(products[1]->P.Vect());
  else return 0.;
}

// The left-handed quark is the same sign quark

// get the tilt angle of the left-handed quark
double particle::tilt_angle() const {
  if (nprod()==2) {
    TLorentzVector* prod;
//    if ( (products[0]->id()>0) == (pid>0) )
    if ( products[0]->id()<0 )
      prod = new TLorentzVector(products[0]->P);
    else
      prod = new TLorentzVector(products[1]->P);

    // Boost the product to the parent rest frame
    prod->Boost(-P.Px()/P.E(),-P.Py()/P.E(),-P.Pz()/P.E());

    double angle = prod->Angle(P.Vect()); // tilt angle

    // if (pid>0) angle = M_PI - angle;

    delete prod;
    return angle;

  } else return 0.;
}
double particle::polar_angle() const { // get the angle wrt beam axis (z)
  return P.Angle(TVector3(0.,0.,1.));
}
double particle::prod_angle() const { // production angle
  if (!( mother[0]==0 || mother[1]==0 )) {

    TLorentzVector boost = mother[0]->P + mother[1]->P;
    boost *= (-1./boost.E());

    TLorentzVector this_particle(P);
    this_particle.Boost(boost.Px(),boost.Py(),boost.Pz());

    return this_particle.Angle(TVector3(0.,0.,1.)); // production angle

  } else return 0.;
}

double particle::dE() const {
  if (nprod()==2) {
    if (products[0]->pid < 0)
      return products[0]->P.E() - products[1]->P.E();
    else
      return products[1]->P.E() - products[0]->P.E();
  } else return 0.;
}
double particle::abs_dE() const {
  if (nprod()==2)
    return abs(products[0]->P.E() - products[1]->P.E());
  else return 0.;
}
double particle::dPt() const {
  if (nprod()==2) {
    if (products[0]->pid < 0)
      return products[0]->P.Pt() - products[1]->P.Pt();
    else
      return products[1]->P.Pt() - products[0]->P.Pt();
  } else return 0.;
}
double particle::abs_dPt() const {
  if (nprod()==2)
    return abs(products[0]->P.Pt() - products[1]->P.Pt());
  else return 0.;
}
double particle::dPz() const {
  if (nprod()==2) {
    if (products[0]->pid < 0)
      return products[0]->P.Pz() - products[1]->P.Pz();
    else
      return products[1]->P.Pz() - products[0]->P.Pz();
  } else return 0.;
}

// *****************************

double particle::ediff() const {
  return dE()/momvecmag();
}
double particle::ediff_abs() const {
  return abs((products[0]->P.E() - products[1]->P.E())/momvecmag());
}
double particle::abs_cos_prod() const {
  return abs(cos(prod_angle()));
}

double particle::op_over_min() const {
  double op  = opening_angle();
  double g   = gamma();
  double min = atan2( 1.,sqrt(g*g - 1) )*2.;
  return op/min;
}

// ******************************************************************

bool particle::is_quark() const {
  for (int i=1;i<=6;i++) if (pid==i) return true;
  return false;
}
bool particle::is_lepton() const {
  for (int i=11;i<=16;i++) if (pid==i) return true;
  return false;
}

// ******************************************************************

void particle::print4vector() const { // print particle's 4 momentum
  printf("(%f,%f,%f,%f)\n",P.E(),P.Px(),P.Py(),P.Pz());
}
