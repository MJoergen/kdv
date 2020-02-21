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
t increases from t=0.

In order to perform numerical simulations, we will discretize (sample) the
function u(x,t) at regular intervals dx and dt, where the values dx and dt are
small but non-zero. The values chosen should be small compared to the length
and time scales of interest.

The first task is to approximate the continuous partial differential equation
into a finite difference equation, discrete in both x and t, and then
to talk about how to actually solve the resulting finite difference.

Next we will discuss how to verify the error introduced by the finite
difference equation, and how to choose appropriate values of dx and dt.


## Discretizing
The original partial differential equation involves the function u(x,t) where x
and t are real numbers. In the numerical simulations we will replace x and t by
discrete values: x = n\*dx and t = m\*dt, where n and m are integers belonging
to a finite and discrete range.

This finiteness of n and m imply a finiteness of x and t. So we must choose a
finite range of values for x and t. For simplicity, lets choose x in [0 ..
x\_max] and t in [0 .. t\_max], where the maximum values are somewhat
arbitrary.

The corresponding ranges for n and m then appear when dividing by the sampling
interval, i.e. m=0, 1, 2, ..., M and n=0, 1, 2, ..., N, where M = t\_max/dt and
N = x\_max/dx.

So instead of working with the function u(x,t) for x in R and t >=0 we now work
with a matrix of values u[n,m] where n=0, 1, 2, ..., N and m=0, 1, 2, ..., M.

Choosing a finite and discrete range of values for x and t also affects the
boundary condition and the initial condition. Starting with the initial
condition, we now have a vector of values u0[n] for n=0, 1, 2, ..., N.

The boundary condition is a bit more tricky, because setting e.g. u[0,m] = 0
will introduce spurious boundary effects. Instead, we will choose periodic
boundary conditions, i.e. u0[0] = u0[N+1].

## Approximating derivatives in x.
The derivatives in space will be approximated using a [five point stencil](https://en.wikipedia.org/wiki/Five-point_stencil)
The first order derivative can be approximated by one of the following equations:
```
dx*u_x     = (-u[n+2] +  8*u[n+1] -            8*u[n-1] + u[n-2])/12
dx^2*u_xx  = (-u[n+2] + 16*u[n+1] - 30*u[n] + 16*u[n-1] - u[n-2])/12
dx^3*u_xxx = ( u[n+2] -  2*u[n+1] +            2*u[n-1] - u[n-2])/2
```

## Conserved quantities
Certain definite integrals over x involving u(x,t) turn out to be
independent of time. Two such are the integral of u and the integral of u^2.

Due to the discrete sampling in x, we will instead consider the following
two sums:

```
mass     = sum_{n=0}^N u[n,m]
momentum = sum_{n=0}^N u[n,m] * u[n,m]
```

These two sums can be evaluated for every value m. According to the above, the
value of the sums should not change. This is a way to monitor the correctness
of the calculation, including monitoring the errors introduced by the discrete
sampling in x and t.

## Numerical integration of the partial differential equation.
It is well known that the simple Euler's method has problems with
numerical stability, so instead we will make use of [Heun's
method](https://en.wikipedia.org/wiki/Heun%27s_method).

