# LSODA.jl: the LSODA algorithm for solving ordinary differential equations

```@contents
```

## Introduction

**LSODA.jl** is a Julia package that interfaces to the [liblsoda](https://github.com/sdwfrost/liblsoda) library, developped by [Simon Frost](http://www.vet.cam.ac.uk/directory/sdf22@cam.ac.uk) ([@sdwfrost](http://github.com/sdwfrost)), thereby providing a way to use the LSODA algorithm from Linda Petzold and Alan Hindmarsh from [Julia](http://julialang.org/). **[Clang.jl](https://github.com/ihnorton/Clang.jl)** has been used to write the library and **[Sundials.jl](https://github.com/JuliaDiffEq/Sundials.jl)** was a inspiring source.

## Installation

To install this package, run the command `Pkg.clone("https://github.com/rveltz/LSODA.jl.git")`

## Example

We first need to load the library.

```example 1
using LSODA
```

We next need to define a function that provides the derivatives `ydot` given the time, `t`, initial conditions, `y`, and optionally some additional data, `data`. Note that this function modifies `ydot` in-place and returns `nothing`.

```example 1
function rhs!(t, y, ydot, data)
	ydot[1]=1.0E4 * y[2] * y[3] - .04E0 * y[1]
	ydot[3]=3.0E7 * y[2] * y[2]
	ydot[2]=-ydot[1] - ydot[3]
  nothing
end
```

The model can be solved by providing an initial condition for the state variables, and a time span over which to simulate.

```example 1
y0 = [1.,0.,0.]
tspan = [0., 0.4]
res =  lsoda(rhs!, y0, tspan, reltol= 1e-4, abstol = Vector([1.e-6,1.e-10,1.e-6]))
```

This should give the following.

```example 1
at t =   4.0000e-01 y=   9.851712e-01   3.386380e-05   1.479493e-02
at t =   4.0000e+00 y=   9.055333e-01   2.240655e-05   9.444430e-02
at t =   4.0000e+01 y=   7.158403e-01   9.186334e-06   2.841505e-01
at t =   4.0000e+02 y=   4.505250e-01   3.222964e-06   5.494717e-01
at t =   4.0000e+03 y=   1.831976e-01   8.941774e-07   8.168016e-01
at t =   4.0000e+04 y=   3.898729e-02   1.621940e-07   9.610125e-01
at t =   4.0000e+05 y=   4.936362e-03   1.984221e-08   9.950636e-01
at t =   4.0000e+06 y=   5.161832e-04   2.065786e-09   9.994838e-01
at t =   4.0000e+07 y=   5.179811e-05   2.072030e-10   9.999482e-01
at t =   4.0000e+08 y=   5.283524e-06   2.113420e-11   9.999947e-01
at t =   4.0000e+09 y=   4.658945e-07   1.863579e-12   9.999995e-01
at t =   4.0000e+10 y=   1.423392e-08   5.693574e-14   1.000000e+00
```

## Application programming interface

### Functions

```@docs
lsoda
```

```@docs
lsoda_evolve!
```
