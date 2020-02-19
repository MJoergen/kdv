# Approximating derivatives in x.
To approximate the first spatial derivative, we will use an equation of the
following form:

* dx\*u\_x = a\*u[n+2] + b\*u[n+1] + c\*u[n] + d\*u[n-1] + e\*u[n-2]

for suitably chosen constants a, ..., e. The aim in the following is to
determine values for these constants, so that the above approximation
is exact for all polynomials up to fourth order.

Inserting therefore the polynomials u(x)=1, u(x)=x, u(x)=x^2, u(x)=x^3, and
u(x)=x^4 gives the five equations:
* 0 = a+b+c+d+e
* 1 = (a+b+c+d+e)\*x + 2a+b-d-2e
* 2x = (a+b+c+d+e)\*x^2 + 2(2a+b-d-2e)\*x + 4a+b+d+4e
* 3x^2 = (a+b+c+d+e)\*x^3 + 3(2a+b-d-2e)\*x^2 + 3(4a+b+d+4e)x + 8a+b-d-8e
* 4x^3 = (a+b+c+d+e)\*x^4 + 4(2a+b-d-2e)\*x^3 + 6(4a+b+d+4e)x^2 + 4(8a+b-d-8e)x + 16a+b+d+16e

which conveniently reduce to:
* 0 = a+b+c+d+e
* 1 = 2a+b-d-2e
* 0 = 4a+b+d+4e
* 0 = 8a+b-d-8e
* 0 = 16a+b+d+16e

whose solution is a=-1/12, b=2/3, c=0, d=-2/3, and e=1/12.  So the first order
derivative can be approximated by:

```
dx*u_x = (-u[n+2] + 8*u[n+1] - 8*u[n-1] + u[n-2])/12
```

However, the above equation only works for points well inside the region of
interest. At the boundary, we can only access values to one side of the point.
So we will also consider an equations of the forms:
* dx\*u\_x = a u[n+3] + b u[n+2] + c u[n+1] + d u[n] + e u[n-1]
* dx\*u\_x = a u[n+4] + b u[n+3] + c u[n+2] + d u[n+1] + e u[n]

By following a similar procedure we arrive at the following results:
 
```
dx*u_x = (   u[n+3] -  6*u[n+2] + 18*u[n+1] - 10*u[n]   -  3*u[n-1])/12
dx*u_x = (-3*u[n+4] + 16*u[n+3] - 36*u[n+2] + 48*u[n+1] - 25*u[n])/12
```

To approximate the second spatial derivative, we follow the same procedure
and arrive at the following set of equations:

```
dx^2*u_xx = (  -u[n+2] + 16*u[n+1] -  30*u[n]   +  16*u[n-1] -    u[n-2])/12
dx^2*u_xx = (  -u[n+3] +  4*u[n+2] +   6*u[n+1] -  20*u[n]   + 11*u[n-1])/12
dx^2*u_xx = (11*u[n+4] - 56*u[n+3] + 114*u[n+2] - 104*u[n+1] + 35*u[n])/12
```

To approximate the third spatial derivative, we follow the same procedure
and arrive at the following set of equations:

```
dx^3*u_xxx = (   u[n+2] -  2*u[n+1] +              2*u[n-1] -   u[n-2])/2
dx^3*u_xxx = (  -u[n+3] +  6*u[n+2] - 12*u[n+1] + 10*u[n]   - 3*u[n-1])/2
dx^3*u_xxx = (-3*u[n+4] + 14*u[n+3] - 24*u[n+2] + 18*u[n+1] - 5*u[n])/2
```

The above equations are all verified in the file test\_deriv.cpp.
