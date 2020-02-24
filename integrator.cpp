#include "integrator.h"
#include "deriv.h"

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

