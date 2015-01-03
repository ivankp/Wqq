#include "LHEF2Wqq.h"

#include <iostream>
#include <iomanip>
#include <exception>

using namespace std;

// ******************************************************************
class lhefex: public exception {
public:
  lhefex(int type) { this->type = type; str=""; }
  virtual const char* what() const throw() { return str; }
  bool is_warning() { return (type==1); }
protected:
  const char* str;
  int type;
};
class lhef_warning: public lhefex {
public:
  lhef_warning(const char* msg) : lhefex(2) { str = msg; }
};
class lhef_err: public lhefex {
public:
  lhef_err(const char* msg) : lhefex(1) { str = msg; }
};

// ******************************************************************

LHEF2Wqq::LHEF2Wqq(std::istream& is)
: LHEFReader(is), fEventsSelected(0)
{ }

LHEF2Wqq::LHEF2Wqq(const std::string& filename, bool msg)
: LHEFReader(filename,msg), fEventsSelected(0)
{ }

bool LHEF2Wqq::vars_set = false;

// ******************************************************************

bool LHEF2Wqq::readEvent() {
  bool status;

  status = LHEFReader::readEvent();
  if (status) {
    IdentifyParticles();
    fEventsSelected++;
  }

  return status;
}

const unsigned long& LHEF2Wqq::nEventsSelected() const {
  return fEventsSelected;
}

void LHEF2Wqq::done_msg() {
  LHEFReader::done_msg();
  cout << "Events selected : " << setw(ndig) << fEventsSelected << endl;
}

void LHEF2Wqq::IdentifyParticles() {
  int isq_i=0, fsq_i=0, Wp_i=0, Wm_i=0;
  for (unsigned i=0; i<eventSize; i++) {
    if (particles[i].get_status()==-1) {
      isq_p[isq_i] = &particles[i]; isq_i++;
      if (isq_i>2) throw lhef_err("Too many initial state particles");
    } else if (particles[i].get_status()==1) {
      fsq_p[fsq_i] = &particles[i]; fsq_i++;
      if (fsq_i>4) throw lhef_err("Too many final state particles");
    } else if (particles[i].get_status()==2 && particles[i].id()==24) {
      Wp_p = &particles[i]; Wp_i++;
      if (Wp_i>1) throw lhef_err("More then one W+");
    } else if (particles[i].get_status()==2 && particles[i].id()==-24) {
      Wm_p = &particles[i]; Wm_i++;
      if (Wm_i>1) throw lhef_err("More then one W-");
    }
  }

  if (isq_i!=2) throw lhef_warning("Not 2 initial state particles");
  if (fsq_i!=4) throw lhef_warning("Not 4 initial state particles");
  if (Wp_i==0)  throw lhef_warning("No W+");
  if (Wm_i==0)  throw lhef_warning("No W-");
}

// ******************************************************************

const particle* LHEF2Wqq::Wp() const { return Wp_p; }
const particle* LHEF2Wqq::Wm() const { return Wm_p; }
const particle* LHEF2Wqq::isq(int i) const {
  if (i>2) throw lhef_err("Requested initial state quark out of index");
  return isq_p[i];
}
const particle* LHEF2Wqq::fsq(int i) const {
  if (i>4) throw lhef_err("Requested final state quark out of index");
  return fsq_p[i];
}

