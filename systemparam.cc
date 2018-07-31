#include "systemparam.hpp"

void adjust_periodic(double &dx, double &dy){
  const double LH = L * 0.5;
  if(dx < -LH) dx += L;
  if(dy < -LH) dy += L;
  if(dx > +LH) dx -= L;
  if(dy > +LH) dy -= L;
}
