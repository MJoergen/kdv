#include <iomanip>
#include <math.h>
#include "kdvintegrator.h"
#include "deriv.h"

// This calculates the right-hand side of the partial differential equation,
// i.e. u_t = -u_xxx + 6uu_x.
void KdvIntegrator::calc_ut(dv& ut, const dv& u)
{
   calc_deriv_1(m_u1, u);  // Calculate u_x
   calc_deriv_3(ut, u);    // Calculate u_xxx
   for (unsigned int n=0; n<u.size(); ++n)
   {
      ut[n] = -ut[n]/(m_dx*m_dx*m_dx) + 6.0*u[n]*m_u1[n]/m_dx;
   }
} // calc_ut

// This calculates the initial condition,
// i.e. -1/2*sech^2(x/2).
double KdvIntegrator::u0(double x) const
{
   const double xmid = (m_nmax/2)*m_dx;
   const double c = cosh((x-xmid) / 2.0);
   return -0.5/(c*c);
} // u0

// This populates the vector with the initial condition.
void KdvIntegrator::initialize(dv& u)
{
   const double xmax = m_nmax*m_dx;
   for (unsigned int n=0; n<=m_nmax; ++n)
   {
      const double x = n*m_dx;
      u[n] = u0(x) + u0(x-xmax) + u0(x+xmax);
   }
} // initialize

// This advances the solution with one time-step.
void KdvIntegrator::single_step(double t, dv& u)
{
   log(t, u);
   IntegratorBase::single_step(u);
} // single_step

// This calculates the mass,
// i.e. integral of u.
double KdvIntegrator::calc_mass(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n];
   }
   return sum*m_dx;
}; // calc_mass

// This calculates the momentum,
// i.e. integral of u^2.
double KdvIntegrator::calc_momentum(const dv& u) const
{
   double sum = 0.0;
   for (unsigned int n=0; n<u.size(); ++n)
   {
      sum += u[n]*u[n];
   }
   return sum*m_dx;
}; // calc_momentum

// This dumps results to files.
void KdvIntegrator::log(double t, const dv& u)
{
   // Calculate initial values of mass and momentum.
   if (m_first)
   {
      m_mass     = calc_mass(u);
      m_momentum = calc_momentum(u);
      m_first    = false;
   }

   if (t >= m_next_t)
   {
      m_next_t += m_dt*m_mmax/10.0;

      // Dump current state to a new file.
      std::ofstream ofs(m_file_name_base + std::to_string(t) + ".txt");
      for (unsigned int n=0; n<=m_nmax; ++n)
      {
         const double x = n*m_dx;
         ofs <<          std::setw(19) << std::setprecision(4) << x
             << " , " << std::setw(19) << std::setprecision(4) << u[n] << std::endl;
      }

      // Output error in conserved quantities.
      m_oflog <<          std::setw(19) << std::setprecision(4) << t
              << " , " << std::setw(19) << std::setprecision(4) << calc_mass(u)     - m_mass
              << " , " << std::setw(19) << std::setprecision(4) << calc_momentum(u) - m_momentum
              << std::endl;
   }
} // log

