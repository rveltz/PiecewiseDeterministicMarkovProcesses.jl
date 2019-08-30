# PDMP Solvers

`solve(prob::PDMPProblem,alg;kwargs)`

Solves the PDMP defined by `prob` using the algorithm `alg`.

## Simulation methods

We provide several methods for the simulation:

- a relatively recent trick, called **CHV**, explained in [paper-2015](http://arxiv.org/abs/1504.06873) which allows to implement the **True Jump Method** without the need to use event detection schemes for the ODE integrator. These event detections can be quite numerically unstable as explained in [paper-2015](http://arxiv.org/abs/1504.06873) and CHV provide a solution to this problem.
- **rejection methods** for which the user is asked to provide a bound on the **total** reaction rate. These last methods are the most "exact" but not the fastest if the reaction rate bound is not tight. In case the flow is known **analytically**, a method is also provided.


These methods require solving stiff ODEs (for CHV) in an efficient manner. [```Sundials.jl```](https://github.com/JuliaLang/Sundials.jl) and [```LSODA.jl```](https://github.com/rveltz/LSODA.jl) are great, but other solvers can also be considered (see [stiff ode solvers](http://lh3lh3.users.sourceforge.net/solveode.shtml) and also the [solvers](http://docs.juliadiffeq.org/stable/solvers/ode_solve.html) from [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl)). Hence, the current package allows the use of all solvers in `DifferentialEquations.jl` thereby giving access to a wide range of solvers. In particular, we can test different solvers to see how precise they are. Here is an example from `examples/pdmpStiff.jl` for which an analytical expression is available which allows computation of the errors

```julia
Comparison of solvers
--> norm difference = 0.00019114008823351014  - solver = cvode
--> norm difference = 0.00014770067837588385  - solver = lsoda
--> norm difference = 0.00018404736432131585  - solver = CVODEBDF
--> norm difference = 6.939603217404056e-5    - solver = CVODEAdams
--> norm difference = 2.216652299580346e-5    - solver = tsit5
--> norm difference = 2.2758951345736023e-6   - solver = rodas4P-noAutoDiff
--> norm difference = 2.496987313804766e-6    - solver = rodas4P-AutoDiff
--> norm difference = 0.0004373003700521849   - solver = RS23
--> norm difference = 2.216652299580346e-5    - solver = AutoTsit5RS23
```

!!! note "ODE Solvers"
    A lot of [care](https://discourse.julialang.org/t/help-reduce-large-gc-time/17215) have been taken to be sure that the algorithms do not allocate and hence are fast. This is based on an iterator interface of `DifferentialEquations`. If you chose `save_positions = (false, false)`, the allocations should be independent from the requested jump number. However, the iterator solution is not yet available for `LSODA` in `DifferentialEquations`. Hence, you can pass `ode = :lsoda` to access an old version of the algorithm (which allocates), or any other solver like `ode = Tsit5()` to access the **new** solver.

## How to chose an algorithm?

The choice of the method `CHV` vs `Rejection` only depends on how much you know about the system. 

More precisely, if the total rate function does not vary much in between jumps, use the rejection method. For example, if the rate is $R(x_c(t)) = 1+0.1\cos(t)$,  then $1+0.1$ will provide a tight bound to use for the rejection method and almost no (fictitious) jumps will be rejected. 

In all other cases, one should try the CHV method where no a priori knowledge of the rate function is needed.

!!! warning "CHV Method"
    A strong requirement for the CHV method is that the total rate (*i.e.* `sum(rate)`) must be positive. This can be easily achieved by adding a dummy Poisson process with very low intensity (see examples).

## Common Solver Options

To simulate a PDMP, one uses `solve(prob::PDMPProblem,alg;kwargs)`. The field are as follows

- `alg` can be `CHV(ode)` (for the [CHV algorithm](https://arxiv.org/abs/1504.06873)), `Rejection(ode)` for the Rejection algorithm and `RejectionExact()` for the rejection algorithm in case the flow in between jumps is known analytically. In this latter case, `prob.F` is used for the specification of the Flow. The ODE solver `ode` can be any solver of [DifferentialEquations.jl](https://github.com/JuliaDiffEq/DifferentialEquations.jl) like `Tsit5()` for example or anyone of the list `[:cvode, :lsoda, :adams, :bdf, :euler]`. Indeed, the package implement an iterator interface which does not work yet with `ode = LSODA()`. In order to have access to the ODE solver `LSODA()`, one should use `ode = :lsoda`.
- `n_jumps = 10`: requires the solver to only compute the first 10 jumps.
- `save_position = (true, false)`: (output control) requires the solver to save the pre-jump and the post-jump states `xc, xd`.
- `verbose = true`: requires the solver to print information concerning the simulation of the PDMP
- `reltol`: relative tolerance used in the ODE solver
- `abstol`: absolute tolerance used in the ODE solver
- `ind_save_c`: which indices of `xc` should be saved
- `ind_save_d`: which indices of `xd` should be saved
- `save_rate = true`: requires the solver to solve the total rate. Can be useful when estimating the rate bounds to use the Rejection algorithm.
-  `X_extended = zeros(Tc, 1 + 1)`: (advanced use) options used to provide the shape of the extended array in the [CHV algorithm](https://arxiv.org/abs/1504.06873). Can be useful to use `StaticArrays.jl`.