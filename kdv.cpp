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
   x1[1] = (     x[4] -  6.0*x[3] + 18.0*x[2] - 10.0*x[1] -  3.0*x[0])/12.0;

   x1[len-2] = (-x[len-5] + 6.0*x[len-4] - 18.0*x[len-3] + 10.0*x[len-2] + 3.0*x[len-1])/12.0;
   x1[len-1] = (3.0*x[len-5] - 16.0*x[len-4] + 36.0*x[len-3] - 48.0*x[len-2] + 25.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x1[i] = (-x[i+2] + 8.0*x[i+1] - 8.0*x[i-1] + x[i-2])/12.0;
   }
} // end of calc_x1

// This calculates the spatial second derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_x2(dv& x2, const dv& x)
{
   const unsigned int len = x.size();
   x2[0] = (11.0*x[4] - 56.0*x[3] + 114.0*x[2] - 104.0*x[1] + 35.0*x[0])/12.0;
   x2[1] = (    -x[4] +  4.0*x[3] +   6.0*x[2] -  20.0*x[1] + 11.0*x[0])/12.0;

   x2[len-2] = (    -x[len-5] +  4.0*x[len-4] +   6.0*x[len-3] -  20.0*x[len-2] + 11.0*x[len-1])/12.0;
   x2[len-1] = (11.0*x[len-5] - 56.0*x[len-4] + 114.0*x[len-3] - 104.0*x[len-2] + 35.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x2[i] = (-x[i+2] + 16.0*x[i+1] - 30*x[i] + 16.0*x[i-1] - x[i-2])/12.0;
   }
} // end of calc_x2

// This calculates the spatial third derivative.
// The coefficients are chosen such that an arbitrary fourth degree polynomial
// gives the exact result.
void calc_x3(dv& x2, const dv& x)
{
   const unsigned int len = x.size();
   x2[0] = (-18.0*x[4] + 84.0*x[3] - 144.0*x[2] + 108.0*x[1] - 30.0*x[0])/12.0;
   x2[1] = ( -6.0*x[4] + 36.0*x[3] -  72.0*x[2] +  60.0*x[1] - 18.0*x[0])/12.0;

   x2[len-2] = (  6.0*x[len-5] - 36.0*x[len-4] +  72.0*x[len-3] -  60.0*x[len-2] + 18.0*x[len-1])/12.0;
   x2[len-1] = ( 18.0*x[len-5] - 84.0*x[len-4] + 144.0*x[len-3] - 108.0*x[len-2] + 30.0*x[len-1])/12.0;

   for (unsigned int i=2; i<len-2; ++i)
   {
      x2[i] = (6.0*x[i+2] - 12.0*x[i+1] + 12.0*x[i-1] - 6.0*x[i-2])/12.0;
   }
} // end of calc_x3

int main()
{
   const unsigned int len = 10;
   dv x(len);
   dv xt(len);

   for (unsigned int i=0; i<len; ++i)
   {
      x[i] = 7.0*i*i*i*i + 6.0*i*i*i + 3.0*i*i + 4.0*i + 5.0; // arbitrary fourth degree polynomial
   }

   calc_x1(xt,x);
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << xt[i] << ", err=" << xt[i] - (28.0*i*i*i+18.0*i*i+6.0*i+4.0) << std::endl; 
   }
   std::cout << std::endl;

   calc_x2(xt,x);
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << xt[i] << ", err=" << xt[i] - (84.0*i*i+36*i+6.0) << std::endl; 
   }
   std::cout << std::endl;

   calc_x3(xt,x);
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << xt[i] << ", err=" << xt[i] - (168.0*i+36) << std::endl; 
   }
   std::cout << std::endl;
}
