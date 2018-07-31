#pragma once

#include <vector>

struct Atom {
  double qx, qy;
  double px, py;
};


class Variables {
  public:
    std::vector<Atom> atoms;
    double time;
    Variables(void){time = 0.0;}
    void add_atoms(double x, double y);
    void add_atoms(double x, double y, double px, double py);
    int num_of_atoms(void){ return static_cast<int>(atoms.size()); }
    void set_init_vel(const double);
};
