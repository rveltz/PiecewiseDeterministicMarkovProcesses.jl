# PiecewiseDeterministicMarkovProcesses.jl 

| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://rveltz.github.io/PiecewiseDeterministicMarkovProcesses.jl/dev) | [![Build status](https://github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl/workflows/CI/badge.svg)](https://github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl/actions) [![codecov](https://codecov.io/gh/rveltz/PiecewiseDeterministicMarkovProcesses.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/rveltz/PiecewiseDeterministicMarkovProcesses.jl) |

PiecewiseDeterministicMarkovProcesses.jl is a Julia package that allows simulation of *Piecewise Deterministic Markov Processes* (PDMP); these encompass hybrid systems and jump processes, comprised of continuous and discrete components, as well as processes with time-varying rates. The aim of the package is to provide methods for the simulation of these processes that are "exact" up to the ODE integrator. A lot of care has been devoted to reduce allocations as much as possible.

To install this package, run the command 

```julia
add PiecewiseDeterministicMarkovProcesses
```


Please, have a look at the [documention](https://rveltz.github.io/PiecewiseDeterministicMarkovProcesses.jl/latest).

# Authors

This is a joint work of [Romain Veltz](https://romainveltz.pythonanywhere.com/) ([@rveltz](http://github.com/rveltz)) and [Simon Frost](http://www.vet.cam.ac.uk/directory/sdf22@cam.ac.uk) ([@sdwfrost](http://github.com/sdwfrost)).
