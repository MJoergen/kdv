#include <vector>
#include <iostream>
#include <math.h>

#include "deriv.h"


int main()
{
   //
   // CUSTOMIZABLE AREA BEGINS HERE vvv
   //

   // The values below determine the region of interest.
   const double xmax = 20.0;
   const double tmax = 20.0;

   // The values below control the accuracy of the calculations.
   const double dx = 0.01;
   const double dt = 0.01;

   // This function defines the initial condition.
   auto u0 = [&] (double x)
   {
      double c = cosh((x-(xmax/2.0)) / 2.0);
      return -0.5/(c*c);
   };

   // This function defines the partial differential equation
   auto calc_ut = [&] (dv& ut, const dv& u)
   {
      dv u1(u.size());
      calc_deriv_1(u1, u);

      calc_deriv_3(ut, u);
      for (unsigned int n=0; n<u.size(); ++n)
      {
         ut[n] = -ut[n]/(dx*dx*dx) + 6.0*u[n]*u1[n]/dx;
      }
   }; // calc_ut

   //
   // CUSTOMIZABLE AREA ENDS HERE ^^^
   //


   const unsigned int nmax = xmax/dx;
   const unsigned int mmax = tmax/dx;

   dv u(nmax+1);
   for (unsigned int n=0; n<=nmax; ++n)
   {
      double x = n*dx;
      u[n] = u0(x);
   }

   dv ut(nmax+1);
   for (unsigned int m=0; m<=mmax; ++m)
   {
      calc_ut(ut, u);
      for (unsigned int n=0; n<=nmax; ++n)
      {
         double x = n*dx;
         std::cout << "x=" << x << " , u=" << u[n] << " , ut=" << ut[n] << std::endl;
         u[n] += ut[n]*dt;
      }
      break;
   }


   return 0;
} // main

