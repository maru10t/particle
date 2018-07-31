#pragma once
#include "variables.hpp"
#include "observer.hpp"

class MD{
  private:
    Variables *vars;
    Observer *obs;
    void makeconf(void);
    void update_pos(void);
    void calc_force(void);
    void periodic(void);
    void calc(void);

  public:
    MD(void);
    MD(int N, double* mass, double* posx, double* posy, double* velx, double* vely);
    ~MD(void);

    void run(void);
    void step(int dstep);
    void write2pointers(int N, double* mass, double* posx, double* posy, double* velx, double* vely);
};
