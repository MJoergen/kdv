#include <vector>
#include <iostream>

#include "deriv.h"

// This calculates the spatial first derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_deriv_1(dv& x1, const dv& x)
{
   const unsigned int len = x.size();
   x1[0] = (-3.0*x[4] + 16.0*x[3] - 36.0*x[2] + 48.0*x[1] - 25.0*x[0])/12.0;
   x1[1] = (     x[4] -  6.0*x[3] + 18.0*x[2] - 10.0*x[1] -  3.0*x[0])/12.0;

   x1[len-2] = (   -x[len-5] +  6.0*x[len-4] - 18.0*x[len-3] + 10.0*x[len-2] +  3.0*x[len-1])/12.0;
   x1[len-1] = (3.0*x[len-5] - 16.0*x[len-4] + 36.0*x[len-3] - 48.0*x[len-2] + 25.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x1[i] = (-x[i+2] + 8.0*x[i+1] - 8.0*x[i-1] + x[i-2])/12.0;
   }
} // end of calc_deriv_1


// This calculates the spatial second derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_deriv_2(dv& x2, const dv& x)
{
   const unsigned int len = x.size();
   x2[0] = (11.0*x[4] - 56.0*x[3] + 114.0*x[2] - 104.0*x[1] + 35.0*x[0])/12.0;
   x2[1] = (    -x[4] +  4.0*x[3] +   6.0*x[2] -  20.0*x[1] + 11.0*x[0])/12.0;

   x2[len-2] = (    -x[len-5] +  4.0*x[len-4] +   6.0*x[len-3] -  20.0*x[len-2] + 11.0*x[len-1])/12.0;
   x2[len-1] = (11.0*x[len-5] - 56.0*x[len-4] + 114.0*x[len-3] - 104.0*x[len-2] + 35.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x2[i] = (-x[i+2] + 16.0*x[i+1] - 30*x[i] + 16.0*x[i-1] - x[i-2])/12.0;
   }
} // end of calc_deriv_2


// This calculates the spatial third derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_deriv_3(dv& x3, const dv& x)
{
   const unsigned int len = x.size();
   x3[0] = (-3.0*x[4] + 14.0*x[3] - 24.0*x[2] + 18.0*x[1] - 5.0*x[0])*0.5;
   x3[1] = (    -x[4] +  6.0*x[3] - 12.0*x[2] + 10.0*x[1] - 3.0*x[0])*0.5;

   x3[len-2] = (    x[len-5] -  6.0*x[len-4] + 12.0*x[len-3] - 10.0*x[len-2] + 3.0*x[len-1])*0.5;
   x3[len-1] = (3.0*x[len-5] - 14.0*x[len-4] + 24.0*x[len-3] - 18.0*x[len-2] + 5.0*x[len-1])*0.5;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x3[i] = (x[i+2] - 2.0*x[i+1] + 2.0*x[i-1] - x[i-2])*0.5;
   }
} // end of calc_deriv_3


