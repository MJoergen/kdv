#include <vector>
#include <iostream>

#include "deriv.h"

const double twelvth = 1.0/12.0;

// This calculates the spatial first derivative,
// assuming periodic boundary conditions.
void calc_deriv_1(dv& x1, const dv& x)
{
   const unsigned int len = x.size();

   x1[0] = (-x[2] + 8.0*x[1] - 8.0*x[len-1] + x[len-2])*twelvth;
   x1[1] = (-x[3] + 8.0*x[2] - 8.0*x[0]     + x[len-1])*twelvth;
   for (unsigned int i=2; i<len-2; ++i)
   {
      x1[i] = (-x[i+2] + 8.0*x[i+1] - 8.0*x[i-1] + x[i-2])*twelvth;
   }
   x1[len-2] = (-x[0] + 8.0*x[len-1] - 8.0*x[len-3] + x[len-4])*twelvth;
   x1[len-1] = (-x[1] + 8.0*x[0]     - 8.0*x[len-2] + x[len-3])*twelvth;
} // end of calc_deriv_1


// This calculates the spatial second derivative,
// assuming periodic boundary conditions.
void calc_deriv_2(dv& x2, const dv& x)
{
   const unsigned int len = x.size();

   x2[0] = (-x[2] + 16.0*x[1] - 30*x[0] + 16.0*x[len-1] - x[len-2])*twelvth;
   x2[1] = (-x[3] + 16.0*x[2] - 30*x[1] + 16.0*x[0]     - x[len-1])*twelvth;
   for (unsigned int i=2; i<len-2; ++i)
   {
      x2[i] = (-x[i+2] + 16.0*x[i+1] - 30*x[i] + 16.0*x[i-1] - x[i-2])*twelvth;
   }
   x2[len-2] = (-x[0] + 16.0*x[len-1] - 30*x[len-2] + 16.0*x[len-3] - x[len-4])*twelvth;
   x2[len-1] = (-x[1] + 16.0*x[0]     - 30*x[len-1] + 16.0*x[len-2] - x[len-3])*twelvth;
} // end of calc_deriv_2


// This calculates the spatial third derivative,
// assuming periodic boundary conditions.
void calc_deriv_3(dv& x3, const dv& x)
{
   const unsigned int len = x.size();

   x3[0] = (x[2] - 2.0*x[1] + 2.0*x[len-1] - x[len-2])*0.5;
   x3[1] = (x[3] - 2.0*x[2] + 2.0*x[0]     - x[len-1])*0.5;
   for (unsigned int i=2; i<len-2; ++i)
   {
      x3[i] = (x[i+2] - 2.0*x[i+1] + 2.0*x[i-1] - x[i-2])*0.5;
   }
   x3[len-2] = (x[0] - 2.0*x[len-1] + 2.0*x[len-3] - x[len-4])*0.5;
   x3[len-1] = (x[1] - 2.0*x[0]     + 2.0*x[len-2] - x[len-3])*0.5;
} // end of calc_deriv_3


