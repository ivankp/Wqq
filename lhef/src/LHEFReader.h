#ifndef LHEFReader_h
#define LHEFReader_h

#include <vector>

#include "particle.h"

namespace LHEF {
  class Reader;
}

class LHEFReader {
public:
  LHEFReader(std::istream& is);
  LHEFReader(const std::string& filename, bool msg = false);
  ~LHEFReader();

  virtual bool readEvent();
  virtual const unsigned long& nEventsRead() const;
  virtual void done_msg();

  virtual unsigned size() const; // event size
  virtual const particle* at(unsigned i) const; // get particle in the event
  virtual const particle* operator[](unsigned i) const; // get particle in the event

protected:
  LHEF::Reader* reader;
  std::vector<particle> particles;
  unsigned eventSize;
  unsigned long fEventsRead;
  int ndig;
  clock_t last_time;

};

#endif
