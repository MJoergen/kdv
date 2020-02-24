#include <fstream>
#include <string>

#include "integrator.h"

class KdvIntegrator : public IntegratorBase
{
   public:
      KdvIntegrator(unsigned int mmax, unsigned int nmax, double dt, double dx, std::string file_name_base) :
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

      void calc_ut(dv& ut, const dv& u);
      void single_step(double t, dv& u);
      void initialize(dv& u);

   private:
      double u0(double x) const;
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

}; // end of class KdvIntegrator

