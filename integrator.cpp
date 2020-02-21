#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

#include "deriv.h"
#include "integrator.h"

void IntegratorBase::single_step(dv& u)
{
   // Heun's method.
   calc_ut(m_f, u);
   for (unsigned int n=0; n<=m_nmax; ++n)
   {
      m_utmp[n] = u[n] + m_dt*m_f[n];
   }
   calc_ut(m_ftmp, m_utmp);
   for (unsigned int n=0; n<=m_nmax; ++n)
   {
      u[n] = u[n] + m_dt*0.5*(m_f[n]+m_ftmp[n]);
   }
} // single_step


void Integrator::calc_ut(dv& ut, const dv& u)
{
   calc_deriv_1(m_u1, u);
   calc_deriv_3(ut, u);
   for (unsigned int n=0; n<u.size(); ++n)
   {
      ut[n] = -ut[n]/(m_dx*m_dx*m_dx) + 6.0*u[n]*m_u1[n]/m_dx;
   }
} // calc_ut

double Integrator::calc_mass(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n];
   }
   return sum*m_dx;
}; // calc_mass

double Integrator::calc_momentum(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n]*u[n];
   }
   return sum*m_dx;
}; // calc_momentum

void Integrator::log(double t, const dv& u)
{
   if (m_first)
   {
      m_mass     = calc_mass(u);
      m_momentum = calc_momentum(u);
      m_first    = false;
   }

   if (t >= m_next_t)
   {
      m_next_t += m_dt*m_mmax/100.0;

      // Output error in conserved quantities.
      m_oflog <<          std::setw(19) << std::setprecision(4) << t
              << " , " << std::setw(19) << std::setprecision(4) << calc_mass(u)     - m_mass
              << " , " << std::setw(19) << std::setprecision(4) << calc_momentum(u) - m_momentum
              << std::endl;

      std::ofstream ofs(m_file_name_base + std::to_string(t) + ".txt");
      for (unsigned int n=0; n<=m_nmax; ++n)
      {
         const double x = n*m_dx;
         ofs <<          std::setw(19) << std::setprecision(4) << x
             << " , " << std::setw(19) << std::setprecision(4) << u[n] << std::endl;
      }
   }

}
