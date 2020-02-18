#include <vector>
#include <iostream>

typedef std::vector<double> dv;

// This calculates the spatial first derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_x1(dv& x1, const dv& x)
{
   const unsigned int len = x.size();
   x1[0] = (-3.0*x[4] + 16.0*x[3] - 36.0*x[2] + 48.0*x[1] - 25.0*x[0])/12.0;

   x1[1] = (x[4] - 6.0*x[3] + 18.0*x[2] - 10.0*x[1] - 3.0*x[0])/12.0;

   x1[len-2] = (-x[len-5] + 6.0*x[len-4] - 18.0*x[len-3] + 10.0*x[len-2] + 3.0*x[len-1])/12.0;

   x1[len-1] = (3.0*x[len-5] - 16.0*x[len-4] + 36.0*x[len-3] - 48.0*x[len-2] + 25.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x1[i] = (-x[i+2] + 8.0*x[i+1] - 8.0*x[i-1] + x[i-2])/12.0;
   }
} // end of calc_x1

// This calculates the spatial second derivative.
// The calculation has a relative error of second order, even at the edge.
// That means x2[i] = x''[i] + O(x''''[i]).
// The coefficients are chosen such that an arbitrary third degree polynomial
// gives the exact result.
void calc_x2(dv& x2, const dv& x)
{
   const unsigned int len = x.size();
   x2[0] = -x[3]+4.0*x[2]-5.0*x[1]+2.0*x[0];
   x2[len-1] = 2.0*x[len-1]-5.0*x[len-2]+4.0*x[len-3]-x[len-4];
   for (unsigned int i=1; i<len-1; ++i)
   {
      x2[i] = x[i+1]-2.0*x[i]+x[i-1];
   }
} // end of calc_x2

int main()
{
   const unsigned int len = 10;
   dv x(len);
   dv xt(len);

   for (unsigned int i=0; i<len; ++i)
   {
      x[i] = 6.0*i*i*i+3.0*i*i+4.0*i+5.0; // arbitrary third degree polynomial
   }
   calc_x1(xt,x);
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << xt[i] << ", err=" << xt[i] - (18.0*i*i+6.0*i+4.0) << std::endl; 
   }
   std::cout << std::endl;

   calc_x2(xt,x);
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << xt[i] << ", err=" << xt[i] - (36*i+6.0) << std::endl; 
   }
}
