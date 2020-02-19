#include <vector>
#include <iostream>
#include <math.h>

// Just a convenient shorthand.
typedef std::vector<double> dv;

int main()
{
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
   //
   // CUSTOMIZABLE AREA ENDS HERE ^^^


   const unsigned int nmax = xmax/dx;
   const unsigned int mmax = tmax/dx;

   dv u(nmax+1);

   for (unsigned int n=0; n<=nmax; ++n)
   {
      double x = n*dx;
      u[n] = u0(x);
   }

   return 0;
} // main

