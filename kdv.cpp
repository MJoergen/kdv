#include <iostream>
#include <iomanip>
#include <math.h>

#include "deriv.h"


int main()
{
   //
   // CUSTOMIZABLE AREA BEGINS HERE vvv
   //

   // The values below determine the region of interest.
   const double xmax = 40.0;
   const double tmax = 10.0;

   // The values below control the accuracy of the calculations.
   const double dx = 0.1;
   const double dt = 0.00001;

   // This function defines the initial condition.
   auto u0 = [&] (double x)
   {
      // Boundary values must be handled separately.
      if (x<=0.0 || x>=xmax)
         return 0.0;

      // This transforms linearly from [0 .. xmax] to [-1 .. 1].
      const double s = 2.0*x/xmax-1.0;

      // This determines the amount of stretching.
      const unsigned int STRETCH = 8;

      // This transforms nonlinearly from [-1 .. 1] to [-infinity .. infinity].
      const double t = s/(1-pow(s,STRETCH))*xmax*0.5;

      // Temporary variable.
      const double c = cosh(t / 2.0);

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
   const unsigned int mmax = tmax/dt;

   dv u(nmax+1);
   for (unsigned int n=0; n<=nmax; ++n)
   {
      const double x = n*dx;
      u[n] = u0(x);
   }

   dv ut(nmax+1);
   // Simple Euler method.
   for (unsigned int m=0; m<=mmax; ++m)
   {
      calc_ut(ut, u);

      if ((m%(mmax/1000)) == 0)
      {
         const double t = m*dt;
         for (unsigned int n=0; n<=nmax; ++n)
         {
            const double x = n*dx;
            std::cout <<    "t="  << std::setw(8) << std::setprecision(4) << t
               << " , x="  << std::setw(8) << std::setprecision(4) << x
               << " , u="  << std::setw(8) << u[n]
               << " , ut=" << std::setw(8) << ut[n]
               << std::endl;
         }
         std::cout << std::endl;
      }

      for (unsigned int n=0; n<=nmax; ++n)
      {
         u[n] += ut[n]*dt;
      }
   }

   return 0;
} // main

