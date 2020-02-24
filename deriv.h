#include <vector>

// The functions in this block calculate the first, second, and third
// derivative of a function sampled at equidistant points.
// The function is assumed to be periodic, meaning that x[len] = x[0], where
// len is the number of elements in the vector x.

// The functions all make use of a five point stencil, as explained here:
// https://en.wikipedia.org/wiki/Five-point_stencil


// Just a convenient shorthand.
typedef std::vector<double> dv;


// This calculates the spatial first derivative.
void calc_deriv_1(dv& x1, const dv& x);

// This calculates the spatial second derivative.
void calc_deriv_2(dv& x2, const dv& x);

// This calculates the spatial third derivative.
void calc_deriv_3(dv& x3, const dv& x);

