# Simulation of the Korteweg-de-Vries equation

In this project we will perform numerical simulations of the
[Korteweg-de-Vries
equation](https://en.wikipedia.org/wiki/Korteweg%E2%80%93de_Vries_equation).


## What is the Korteweg-de-Vries equation?
The Korteweg-de-Vries equation is the partial differential equation given by:
```
u_t = -u_xxx + 6*u*u_x.
```
Here the function u(x,t) is defined for all x in R, and for t >= 0. The
boundary condition is:

* u(x,t) -> 0 for x -> infinity, at all times t >= 0.

The initial condition is:

* u(x,t) = u0(x), for t=0.

where of course the function u0(x) satisfies the boundary condition too.


## Overview
The numerical simulations we will perform here are based on a given
initial function u0(x), and solving the partial differential equation
numerically to determine how the solution function u(x,t) evolves as time
t increases.

In order to perform numerical simulations, we will discretize (sample) the
function u(x,t) at regular intervals dx and dt, where the values dx and dt are
small but non-zero. The values chosen should be small compared to the length
and time scales of interest. Since we will mainly be interested in length and
time scales of the order of 1, we will choose values of dx and dt both around
0.01, at least initially.

The first task is to approximate the continuous partial differential equation
into a finite difference equation, discrete in both x and t, and then
to talk about how to actually solve the resulting finite difference.

Next we need to discuss how to verify the error introduced by the finite
difference equation, and how to choose appropriate values of x and t.


## Discretizing
The original partial differential equation involves the function u(x,t) where x
and t are real numbers. In a numerical simulation we need to replace x and t by
discrete values x = n\*dx and t = m\*dt, where n and m are integers belonging
to a finite and discrete range.

This finiteness of n and m imply a finiteness of x and t. So we must choose a
finite range of values for x and t. For simplicity, lets choose x in [0 ..
x\_max] and t in [0 .. t\_max], where the maximum values are somewhat
arbitrary. For now, let's put x\_max = 20 and t\_max = 20.

The corresponding ranges for n and m then appear when dividing by the sampling
interval, i.e. m=0, 1, 2, ..., M and n=0, 1, 2, ..., N, where M = t\_max/dt and
N = x\_max/dx.

So instead of working with the function u(x,t) for x in R and t >=0 we now work
with a matrix of values u[n,m] where n=0, 1, 2, ..., N and m=0, 1, 2, ..., M.

Choosing a finite and discrete range of values for x and t also affects the
boundary condition and the initial condition. Starting with the initial
condition, we now have a vector of values u0[n] for n=0, 1, 2, ..., N.

The boundary condition is a bit more tricky, because setting e.g. u[0,m] = 0
will introduce spurious boundary effects. Instead, we will choose a finite
difference equation where the boundary condition does not explicitly appear.

## Approximating derivatives in x.
The first order derivative can be approximated by one of the following equations:
```
dx*u_x = (  -u[n+2] +  8*u[n+1] -              8*u[n-1] +    u[n-2])/12
dx*u_x = (   u[n+3] -  6*u[n+2] + 18*u[n+1] - 10*u[n]   -  3*u[n-1])/12
dx*u_x = (-3*u[n+4] + 16*u[n+3] - 36*u[n+2] + 48*u[n+1] - 25*u[n])/12
```
Notice how the first equation is centered around u[n]; this is suitable for
calculations within the interval of interest.
The other two equations are skewed, and are suitable for calculations at or near
the boundary.

The second order derivative can be approximated by one of the following equations:
```
dx^2*u_xx = (  -u[n+2] + 16*u[n+1] -  30*u[n]   +  16*u[n-1] -    u[n-2])/12
dx^2*u_xx = (  -u[n+3] +  4*u[n+2] +   6*u[n+1] -  20*u[n]   + 11*u[n-1])/12
dx^2*u_xx = (11*u[n+4] - 56*u[n+3] + 114*u[n+2] - 104*u[n+1] + 35*u[n])/12
```

The third order derivative can be approximated by one of the following equations:
```
dx^3*u_xxx = (   u[n+2] -  2*u[n+1] +              2*u[n-1] -   u[n-2])/2
dx^3*u_xxx = (  -u[n+3] +  6*u[n+2] - 12*u[n+1] + 10*u[n]   - 3*u[n-1])/2
dx^3*u_xxx = (-3*u[n+4] + 14*u[n+3] - 24*u[n+2] + 18*u[n+1] - 5*u[n])/2
```

The derivation of these equations can be found [here](derivatives.md).

