#ifndef particle_h
#define particle_h

#include <vector>

#include "TLorentzVector.h"

class particle {
public:
  particle(int num, int id, int stat, const std::vector<double>& pup);
  ~particle();

  int      id() const;
  int      get_status() const;
  void     add(const particle* p);
  void     set_mother(unsigned i,const particle* p);
  const particle* get_mother(unsigned i) const;

  unsigned nprod() const;
  const    particle* prod(unsigned i) const;
  const    particle* operator[](unsigned i) const;

  double   mass() const;
  double   mass2() const;
  double   energy() const;
  double   energy_cm() const;
  double   momvecmag() const;
  double   gamma() const;
  double   Pt() const;
  double   opening_angle() const;
  double   tilt_angle() const;
  double   polar_angle() const;
  double   prod_angle() const;
  double   dE() const;
  double   abs_dE() const;
  double   dPt() const;
  double   abs_dPt() const;
  double   dPz() const;

  double   ediff() const;
  double   ediff_abs() const;
  double   abs_cos_prod() const;
  double   op_over_min() const;

  bool     is_quark() const;
  bool     is_lepton() const;

  void     print4vector() const;

private:
  int n, pid, status;
	TLorentzVector P;
  double m;
	const particle* mother[2];
	std::vector<const particle*> products;
};

#endif
