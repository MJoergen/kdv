#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <math.h>

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

      void single_step(dv& u);

   protected:
      virtual void calc_ut(dv& ut, const dv& u) = 0;

      unsigned int m_nmax;
      double       m_dt;
      double       m_dx;
      dv           m_f;
      dv           m_utmp;
      dv           m_ftmp;

}; // end of class IntegratorBase

class Integrator : public IntegratorBase
{
   public:
      Integrator(unsigned int mmax, unsigned int nmax, double dt, double dx, std::string file_name_base) :
         IntegratorBase(nmax, dt, dx),
         m_mmax(mmax),
         m_u1(dv(nmax+1)),
         m_oflog(file_name_base+".txt"),
         m_file_name_base(file_name_base),
         m_next_t(0.0),
         m_first(true),
         m_mass(0.0),
         m_momentum(0.0)
      {}

      virtual void calc_ut(dv& ut, const dv& u);

      double calc_mass(const dv& u) const;
      double calc_momentum(const dv& u) const;
      void log(double t, const dv& u);

   private:
      unsigned int  m_mmax;
      dv            m_u1;
      std::ofstream m_oflog;
      std::string   m_file_name_base;
      double        m_next_t;
      bool          m_first;
      double        m_mass;
      double        m_momentum;

}; // end of class Integrator

