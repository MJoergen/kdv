#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

#include "deriv.h"


int main()
{
   //
   // CUSTOMIZABLE AREA BEGINS HERE vvv
   //

   // The values below determine the region of interest.
   const double xmax = 40.0;
   const double tmax = 100.0;

   // The values below control the accuracy of the calculations.
   const double dx = 0.1;
   const double dt = 0.00001;

   // This function defines the initial condition.
   auto u0 = [&] (double x)
   {
      const double xmid = xmax/2.0;

      // Temporary variable.
      const double c = cosh((x-xmid) / 2.0);

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

   // This function calculates the mass
   auto calc_mass = [&] (const dv& u)
   {
      double sum = 0.0;
      for (unsigned int n=0; n<u.size(); ++n)
      {
         sum += u[n];
      }
      return sum;
   }; // calc_mass

   // This function calculates the momentum
   auto calc_momentum = [&] (const dv& u)
   {
      double sum = 0.0;
      for (unsigned int n=0; n<u.size(); ++n)
      {
         sum += u[n]*u[n];
      }
      return sum;
   }; // calc_momentum

   // Initial values of mass and momentum
   const double exp_mass     = -20.0;
   const double exp_momentum = 20.0/3.0;

   //
   // CUSTOMIZABLE AREA ENDS HERE ^^^
   //


   const unsigned int nmax = xmax/dx + 0.5;  // Add 0.5 for rounding.
   const unsigned int mmax = tmax/dt + 0.5;  // Add 0.5 for rounding.

   const std::string file_name_base = "kdv_out";

   // Prepare initial condition.
   // Make sure initial values are (approximately) periodic, even though the
   // function u0() is not.
   dv u(nmax+1);
   for (unsigned int n=0; n<=nmax; ++n)
   {
      const double x = n*dx;
      u[n] = u0(x) + u0(x-xmax) + u0(x+xmax);
   }

   dv f(nmax+1);
   dv utmp(nmax+1);
   dv ftmp(nmax+1);
   std::ofstream oflog(file_name_base + ".txt");
   for (unsigned int m=0; m<=mmax; ++m)
   {
      const double t = m*dt;
      if ((m%(mmax/1000)) == 0)
      {
         // Output error in conserved quantities.
         oflog <<          std::setw(19) << std::setprecision(4) << t
               << " , " << std::setw(19) << std::setprecision(4) << calc_mass(u)     - exp_mass
               << " , " << std::setw(19) << std::setprecision(4) << calc_momentum(u) - exp_momentum
               << std::endl;

         std::ofstream ofs(file_name_base + std::to_string(t) + ".txt");
         for (unsigned int n=0; n<=nmax; ++n)
         {
            const double x = n*dx;
            ofs <<          std::setw(19) << std::setprecision(4) << x
                << " , " << std::setw(19) << std::setprecision(4) << u[n] << std::endl;
         }
      }

      // Heun's method.
      calc_ut(f, u);
      for (unsigned int n=0; n<=nmax; ++n)
      {
         utmp[n] = u[n] + dt*f[n];
      }
      calc_ut(ftmp, utmp);
      for (unsigned int n=0; n<=nmax; ++n)
      {
         u[n] = u[n] + dt/2.0*(f[n]+ftmp[n]);
      }
   }

   return 0;
} // main

