#include "deriv.h"

// This base class encapsulates Heuns
class IntegratorBase
{
   public:
      IntegratorBase(unsigned int nmax, double dt, double dx)
         : m_nmax(nmax),
           m_dt  (dt),
           m_dx  (dx),
           m_f   (dv(nmax+1)),
           m_utmp(dv(nmax+1)),
           m_ftmp(dv(nmax+1))
      {}

      virtual void single_step(dv& u);

   protected:
      virtual void calc_ut(dv& ut, const dv& u) = 0;

      unsigned int m_nmax;
      double       m_dt;
      double       m_dx;
      dv           m_f;
      dv           m_utmp;
      dv           m_ftmp;

}; // end of class IntegratorBase

