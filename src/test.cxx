#include <iostream>
#include <vector>

#include "particle.h"
#include "parser.h"

using namespace std;
using namespace parser;
/*
class particle_theta_op: public unitary_op {
  double x;
  public:
  particle_theta_op() : unitary_op("double") { }
  void* operator()(void* r) {
    const particle* p = *static_cast<const particle**>(r);
    x = p->theta();
    cout << p << ' ';
    cout << p->id() << ' ' << p->theta() << endl;
    return &x;
  }
};*/

int main(int argc, char* argv[])
{
  /*vector<double> pup {0.10132882063E+03, -0.12179056536E+03,  0.18479194345E+03,  0.25620766518E+03,  0.79961450919E+02};
  particle* p = new particle(1,24,2,pup);

  parse::add_var("p","particle",&p);
  unitary_op::newop("theta","particle",new particle_theta_op());

  parse cut(argv[1]); // "(abs cos theta @p) > 0.8"

  cout << "cut = " << cut.eval<int>() << endl;
*/
  return 0;
}
