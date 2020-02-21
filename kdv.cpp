#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

#include "deriv.h"
#include "integrator.h"


int main()
{
   //
   // CUSTOMIZABLE AREA BEGINS HERE vvv
   //

   // The values below determine the region of interest.
   const double xmax = 100.0;
   const double tmax = 10.0;

   // The values below control the accuracy of the calculations.
   const double dx = 0.05;
   const double dt = 0.000001;

   // This function defines the initial condition.
   auto u0 = [&] (double x)
   {
      const double xmid = xmax/2.0;

      // Temporary variable.
      const double c = cosh((x-xmid) / 2.0);

      return -0.5/(c*c);
   };

   // Output file name prefix
   const std::string file_name_base = "kdv_out";

   //
   // CUSTOMIZABLE AREA ENDS HERE ^^^
   //


   const unsigned int nmax = xmax/dx + 0.5;  // Add 0.5 for rounding.
   const unsigned int mmax = tmax/dt + 0.5;  // Add 0.5 for rounding.

   // Prepare initial condition.
   // Make sure initial values are (approximately) periodic, even though the
   // function u0() is not.
   dv u(nmax+1);
   for (unsigned int n=0; n<=nmax; ++n)
   {
      const double x = n*dx;
      u[n] = u0(x) + u0(x-xmax) + u0(x+xmax);
   }

   Integrator integrator(mmax, nmax, dt, dx, file_name_base);

   for (unsigned int m=0; m<=mmax; ++m)
   {
      const double t = m*dt;
      integrator.log(t, u);
      integrator.single_step(u);
   }

   return 0;
} // main

