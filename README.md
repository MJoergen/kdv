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

* u(x,t) -> 0 for x -> infinity at all times t >= 0.

The initial condition is:

* u(x,t) = u0(x) for t=0.

where of course the function u0(x) satisfies the boundary condition too.

## Overview
The numerical simulations we will perform here are based on a given
initial function u0(x), and solving the partial differential equation
numerically to determine how the solution function u(x,t) evolves as time
t increases.

In order to perform numerical simulations, we will discretize (sample) the
function u(x,t) at regular intervals dx and dt, where the values dx and dt
are small but non-zero.

The first task is to approximate the continuous partial differential equation
into a finite difference equation, discrete in both x and t


## Approximating derivatives in x.
The first order derivative can be approximated by one of the following equations:
```
dx*u_x = (  -u[n+2] +  8*u[n+1] -              8*u[n-1] +    u[n-2])/12
dx*u_x = (   u[n+3] -  6*u[n+2] + 18*u[n+1] - 10*u[n]   -  3*u[n-1])/12
dx*u_x = (-3*u[n+4] + 16*u[n+3] - 36*u[n+2] + 48*u[n+1] - 25*u[n])/12
```

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

