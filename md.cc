#include "md.hpp"
#include "variables.hpp"
#include "observer.hpp"
#include "systemparam.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdio>

using namespace std;

MD::MD(void){
  vars = new Variables();
  obs = new Observer();
}

MD::MD(int N, double* mass, double* posx, double* posy, double* velx, double* vely){
  vars = new Variables();
  obs = new Observer();

  for(int i = 0; i < N; ++i){
    vars->add_atoms(posx[i], posy[i], velx[i], vely[i]);
  }

}

MD::~MD(void){
  delete vars;
  delete obs;
}

void MD::write2pointers(int N, double* mass, double* posx, double* posy, double* velx, double* vely){
  int i = 0;
  for(auto &a : vars->atoms){
    posx[i] = a.qx;
    posy[i] = a.qy;
    velx[i] = a.px;
    vely[i] = a.py;
    ++i;
  }

}

void MD::makeconf(void){
  const double density = 0.50;
//  const double s = 1.0/pow(density*0.25,1.0/3.0);
  const double s = 1.0/pow(density*0.25,1.0/2.0);
  const double hs = s * 0.5;
  const int is = static_cast<int>(L/s);

  // init configuration of postition
  for(int iy = 0; iy < is; ++iy){
    for(int ix = 0; ix < is; ++ix){
      vars->add_atoms(ix * s, iy * s);
      //vars->add_atoms(ix * s + hs, iy * s);
    }
  }

  vars->set_init_vel(1.0);
}

void MD::update_pos(){
  const double dt2 = dt * 0.5;

  for(auto &a : vars->atoms){
    a.qx += a.px * dt2;
    a.qy += a.py * dt2;
  }

}

void MD::calc_force(){
  const int pn = vars->num_of_atoms();
  Atom *atoms = vars->atoms.data();

  for(int i = 0; i < pn-1; ++i){
    for(int j = i+1; j < pn; ++j){
      double dx,dy;
      dx = atoms[j].qx - atoms[i].qx;
      dy = atoms[j].qy - atoms[i].qy;

      adjust_periodic(dx,dy);
      double r2 = dx*dx+dy*dy;

      if(r2 > CL2)
        continue;

      double r6 = r2*r2*r2;
      double df = (24.0 * r6 - 48.0) / (r6*r6*r2)*dt;

      atoms[i].px += df*dx;
      atoms[i].py += df*dy;

      atoms[j].px -= df*dx;
      atoms[j].py -= df*dy;
    }
  }

}

void MD::periodic(){
  for(auto &a : vars->atoms){
    if(a.qx < 0.0) a.qx += L;
    if(a.qy < 0.0) a.qy += L;
    if(a.qx > L) a.qx -= L;
    if(a.qy > L) a.qy -= L;

    assert(a.qx < L);
    assert(a.qy < L);
  }
}

void MD::calc(){
  update_pos();
  calc_force();
  update_pos();
  periodic();
  vars->time += dt;
}

void MD::run(){
  makeconf();
  const int STEPS = 10000;
  const int OBSERVE = 100;

  for(int i = 0; i < STEPS; ++i){
    if( (i % OBSERVE) == 0){
      double k = obs->kinetic_energy(vars);
      double v = obs->potential_energy(vars);
      cout << vars->time << " ";
      cout << k << " ";
      cout << v << " ";
      cout << k+v << endl;
    }
    calc();
  }
}

void MD::step(int dstep){
  for(int i = 0; i < dstep; ++i){
    calc();
  }
}
