#include "variables.hpp"
#include <random>

using namespace std;

// init velocity
// such that...
// Norm: V0
// Direction: random [0,2pi]
void Variables::set_init_vel(const double V0){
  mt19937 mt(2);
  uniform_real_distribution<double> ud(0.0,1.0);

  double avx,avy;
  avx=avy=0.0;

  for(auto &a : atoms){
    double phi = 2.0 * ud(mt) * M_PI; // return [0, 2pi]
    double vx = V0 * cos(phi);
    double vy = V0 * sin(phi);

    a.px = vx;
    a.py = vy;
    avx += vx;
    avy += vy;
  }

  const int pn = atoms.size();
  avx /= static_cast<double>(pn);
  avy /= static_cast<double>(pn);

  // averaged velocity subtracted
  for(auto &a : atoms){
    a.px -= avx;
    a.py -= avy;
  }
}

void Variables::add_atoms(double x, double y){
  Atom a;
  a.qx = x;
  a.qy = y;
  a.px = 0.0;
  a.py = 0.0;
  atoms.push_back(a);
}

void Variables::add_atoms(double x, double y, double px, double py){
  Atom a;
  a.qx = x;
  a.qy = y;
  a.px = px;
  a.py = py;
  atoms.push_back(a);
}
