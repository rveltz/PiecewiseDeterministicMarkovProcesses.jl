# PiecewiseDeterministicMarkovProcesses.jl 

[![Build Status](https://travis-ci.org/rveltz/PDMP.jl.svg?branch=master)](https://travis-ci.org/rveltz/PDMP.jl)
<a href='https://coveralls.io/github/rveltz/PDMP.jl?branch=master'><img src='https://coveralls.io/repos/github/rveltz/PDMP.jl/badge.svg?branch=master' alt='Coverage Status' /></a>
<!--[![Build status](https://ci.appveyor.com/api/projects/status/github/rveltz/PDMP.jl?svg=true&branch=master)](https://ci.appveyor.com/project/rveltz/pdmp-jl/branch/master)
--><!--[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://rveltz.github.io/PDMP.jl/stable)-->
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://rveltz.github.io/PDMP.jl/latest) 

PiecewiseDeterministicMarkovProcesses.jl is a Julia package that allows simulation of *Piecewise Deterministic Markov Processes* (PDMP); these encompass hybrid systems and jump processes, comprised of continuous and discrete components, as well as processes with time-varying rates. The aim of the package is to provide methods for the simulation of these processes that are "exact" up to the ODE integrator. A lot of care has been devoted to reduce allocations as much as possible.

To install this (unregistered) package, run the command 

```julia
add https://github.com/rveltz/PDMP.jl.git
```


Please, have a look at the [documention](https://rveltz.github.io/PiecewiseDeterministicMarkovProcesses.jl/latest).

# Authors

This is a joint work of [Romain Veltz](https://romainveltz.pythonanywhere.com/) ([@rveltz](http://github.com/rveltz)) and [Simon Frost](http://www.vet.cam.ac.uk/directory/sdf22@cam.ac.uk) ([@sdwfrost](http://github.com/sdwfrost)).
