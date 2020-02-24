#include <iomanip>
#include "kdvintegrator.h"
#include "deriv.h"

void KdvIntegrator::calc_ut(dv& ut, const dv& u)
{
   calc_deriv_1(m_u1, u);
   calc_deriv_3(ut, u);
   for (unsigned int n=0; n<u.size(); ++n)
   {
      ut[n] = -ut[n]/(m_dx*m_dx*m_dx) + 6.0*u[n]*m_u1[n]/m_dx;
   }
} // calc_ut

void KdvIntegrator::single_step(double t, dv& u)
{
   log(t, u);
   IntegratorBase::single_step(u);
} // single_step

double KdvIntegrator::calc_mass(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n];
   }
   return sum*m_dx;
}; // calc_mass

double KdvIntegrator::calc_momentum(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n]*u[n];
   }
   return sum*m_dx;
}; // calc_momentum

void KdvIntegrator::log(double t, const dv& u)
{
   if (m_first)
   {
      m_mass     = calc_mass(u);
      m_momentum = calc_momentum(u);
      m_first    = false;
   }

   if (t >= m_next_t)
   {
      m_next_t += m_dt*m_mmax/10.0;

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
} // log

