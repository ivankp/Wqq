#ifndef LHEF2Wqq_h
#define LHEF2Wqq_h

#include "LHEFReader.h"

namespace parser {
  class parse;
}

class LHEF2Wqq : public LHEFReader {
public:
  LHEF2Wqq(std::istream& is);
  LHEF2Wqq(const std::string& filename, bool msg = false);

  virtual bool readEvent();
  virtual const unsigned long& nEventsSelected() const;
  virtual void done_msg();

  const particle* Wp() const;
  const particle* Wm() const;
  const particle* isq(int i) const;
  const particle* fsq(int i) const;

protected:
  static bool vars_set;
  unsigned long fEventsSelected;

  const particle *Wp_p, *Wm_p, *isq_p[2], *fsq_p[4];

  virtual void IdentifyParticles();
};

#endif
