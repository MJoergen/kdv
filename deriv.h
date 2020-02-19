#include <vector>
#include <iostream>

// Just a convenient shorthand.
typedef std::vector<double> dv;

// This calculates the spatial first derivative.
void calc_deriv_1(dv& x1, const dv& x);

// This calculates the spatial second derivative.
void calc_deriv_2(dv& x2, const dv& x);

// This calculates the spatial third derivative.
void calc_deriv_3(dv& x3, const dv& x);

