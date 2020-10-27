## Mathematical Specification of a PDMP Problem

### Vector field

To define a PDMP Problem, you first need to give the function $F$ and the initial condition $x_{c,0}$ which define an ODE:

```math
\frac{dx_c}{dt} = F(x_c(t),x_d(t),p,t)
```

where $F$ should be specified in-place as `F(dxc,xc,xd,p,t)`, and `xc0` should be an `AbstractArray` (or number) whose geometry matches the desired geometry of `xc`. Note that we are not limited to numbers or vectors for `xc0`; one is allowed to provide u₀ as arbitrary matrices / higher dimension tensors as well.

### Jumps

Jumps are defined as a Jump process which changes states at some rate $R$ which is a scalar function of the type 

```math
R(x_c(t),x_d(t),p,t).
```

Note, that in between jumps, $x_d(t)$ is constant but $x_c(t)$ is allowed to evolve.
$R$ should be specified in-place as `R(rate,xc,xd,p,t,issum::Bool)` where it mutates `rate`. Note that a boolean `issum` is provided and the behavior of `R` should be as follows

- if `issum == true`, we only require `R` to return the total rate, *e.g.* `sum(rate)`. We use this formalism because sometimes you can compute the `sum` without mutating `rate`.
- if `issum == false`, `R` must populate `rate` with the updated rates

We then need to provide the way the jumps affect the state variable. There are two possible ways here:

- either give a transition matrix `nu`: it will only affect the discrete component `xd` and leave `xc` unaffected.
- give a function to implement jumps `Delta(xc, xd, parms, t, ind_reaction::Int64)` where you can mutate `xc,xd` or `parms`. The argument `ind_reaction` is the index of the reaction at which the jump occurs. See `examples/pdmp_example_eva.jl` for an example.

## Problem Type

### Constructors

- `PDMPProblem(F,R,Delta,nu,xc0,xd0,p,tspan)`
- `PDMPProblem(F,R,nu,xc0,xd0,p,tspan)` when ones does not want to provide the function `Delta`
- `PDMPProblem(F,R,Delta,reaction_number::Int64,xc0,xd0,p,tspan)` when ones does not want to provide the transition matrix `nu`. The length `reaction_number` of the rate vector must then be provided.

We also provide a wrapper to [DiffEqJump.jl](https://github.com/JuliaDiffEq/DiffEqJump.jl). This is quite similar to how a `JumpProblem` would be created.

- `PDMPProblem(prob, jumps...)` where `prob` can be an `ODEProblem`. For an example, please consider `example/examplediffeqjumpwrapper.jl`.

### Fields
- `F`: the function of the ODE
- `R`: the function to compute the transition rates
- `Delta` [Optional]: the function to effect the jumps
- `nu` [Optional]: the transition matrix
- `xc0`: the initial condition of the continuous part
- `xd0`: the initial condition of the discrete part
- `p`: the parameters to be provided to the functions `F, R, Delta`
- `tspan`: The timespan for the problem.
