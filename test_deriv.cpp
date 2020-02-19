#include <vector>
#include <iostream>

#include "deriv.h"

int main()
{
   const unsigned int len = 10;
   dv x(len);        // Input vector
   dv x1(len);       // Calculated first derivative
   dv x2(len);       // Calculated second derivative
   dv x3(len);       // Calculated third derivative
   dv exp_x1(len);   // Expected first derivative
   dv exp_x2(len);   // Expected second derivative
   dv exp_x3(len);   // Expected third derivative

   // Define the parameters of an arbitrary fourth degree polynomial
   double a = 7.0;
   double b = 6.0;
   double c = 5.0;
   double d = 4.0;
   double e = 3.0;

   // Prepare for test
   for (unsigned int i=0; i<len; ++i)
   {
      x[i]      =    a*i*i*i*i +   b*i*i*i +   c*i*i + d*i + e;  // arbitrary fourth degree polynomial
      exp_x1[i] =  4*a*i*i*i   + 3*b*i*i   + 2*c*i   + d;        // First derivative
      exp_x2[i] = 12*a*i*i     + 6*b*i     + 2*c;                // Second derivative
      exp_x3[i] = 24*a*i       + 6*b;                            // Third derivative
   }
   calc_deriv_1(x1,x);
   calc_deriv_2(x2,x);
   calc_deriv_3(x3,x);

   std::cout << "Verifying first derivative:" << std::endl;
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << x1[i] << ", err=" << x1[i] - exp_x1[i] << std::endl; 
   }
   std::cout << std::endl;

   std::cout << "Verifying second derivative:" << std::endl;
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << x2[i] << ", err=" << x2[i] - exp_x2[i] << std::endl; 
   }
   std::cout << std::endl;

   std::cout << "Verifying third derivative:" << std::endl;
   for (unsigned int i=0; i<len; ++i)
   {
      std::cout << x3[i] << ", err=" << x3[i] - exp_x3[i] << std::endl; 
   }
   std::cout << std::endl;
}
