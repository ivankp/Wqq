#include "LHEFReader.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>

#include "LHEF.h"

using namespace std;

// Constructors *******************************************

LHEFReader::LHEFReader(istream& is) : fEventsRead(0)
{ reader = new LHEF::Reader(is);
  last_time = clock();
}

LHEFReader::LHEFReader(const string& filename, bool msg) : fEventsRead(0)
{ reader = new LHEF::Reader(filename);
  if (msg) cout<<"Opened file: "<<filename<<endl;
  last_time = clock();
}

LHEFReader::~LHEFReader()
{ delete reader; }

// ********************************************************

bool LHEFReader::readEvent() {
  bool status = reader->readEvent();

  if (status) {
    // number of particles in the event
    eventSize = reader->hepeup.IDUP.size();

    // reset particles vector
    particles.clear();

    // loop over particles to put them into vector
    for (unsigned i=0; i<eventSize; i++) {
      particles.push_back( particle(i+1,
        reader->hepeup.IDUP[i],
        reader->hepeup.ISTUP[i],
        reader->hepeup.PUP[i]
      ) );
    }
    // assign particles' mothers and products
    for (unsigned i=0; i<eventSize; i++) {
      int mother1 = reader->hepeup.MOTHUP[i].first;
      int mother2 = reader->hepeup.MOTHUP[i].second;
      if (mother1) {
        particles[i].set_mother(0,&particles[mother1-1]);
        particles[mother1-1].add(&particles[i]);
      }
      if (mother2 && (mother2!=mother1)) {
        particles[i].set_mother(1,&particles[mother2-1]);
        particles[mother2-1].add(&particles[i]);
      }
    }

    fEventsRead++;

    if ( (clock()-last_time)/CLOCKS_PER_SEC > 1 ) {
      cout << setw(10) << fEventsRead;
      cout.flush();
      for (char i=0;i<10;i++) cout << '\b';
      last_time = clock();
    }

  }

  return status;
}

const unsigned long& LHEFReader::nEventsRead() const {
  return fEventsRead;
}

void LHEFReader::done_msg() {
  ndig = 0;
  if (fEventsRead>0) ndig = log10(fEventsRead)+1;
  cout << "Events processed: " << setw(ndig) << fEventsRead << endl;
}

// ********************************************************

unsigned LHEFReader::size() const {
  return eventSize;
}
const particle* LHEFReader::at(unsigned i) const {
  return &particles.at(i);
}
const particle* LHEFReader::operator[](unsigned i) const {
  return &particles.at(i);
}
