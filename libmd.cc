#include "md.hpp"
#include <cstdio>

extern "C" void calc_step_md(int dstep, int N, double* mass, double* posx, double* posy, double* velx, double* vely){
  MD *md = new MD(N,mass,posx,posy,velx,vely);

  md->step(dstep);

  md->write2pointers(N,mass,posx,posy,velx,vely);
  delete md;
}
