# PiecewiseDeterministicMarkovProcesses.jl 


PiecewiseDeterministicMarkovProcesses.jl is a Julia package that allows simulation of *Piecewise Deterministic Markov Processes* (PDMP); these encompass hybrid systems and jump processes, comprised of continuous and discrete components, as well as processes with time-varying rates. The aim of the package is to provide methods for the simulation of these processes that are "statistically exact" up to the ODE integrator.

*If you find this package useful, please star the repo. If you use it in your work, please cite this code and send us an email so that we can cite your work here.*

## Definition of the Jump process

We briefly recall facts about a simple class of PDMPs. They are described by a couple $(x_c, x_d)$ where $x_c$ is solution of the differential equation $$\frac{dx_c(t)}{dt} = F(x_c(t),x_d(t),p,t).$$ The second component $x_d$ is a piecewise constant array with type `Int` and `p` are some parameters. The jumps occur at rates $R(x_c(t),x_d(t),p,t)$. At each jump, $x_d$ or $x_c$ can be affected.


## Related projects

- [Gillespie.jl](https://github.com/sdwfrost/Gillespie.jl): a package for simulation of pure Jump processes, *i.e.* without the continuous part $x_c$.
- [DiffEqJump.jl](https://github.com/JuliaDiffEq/DiffEqJump.jl): similar to our setting with different sampling algorithm
- [PDSampler.jl](https://github.com/alan-turing-institute/PDSampler.jl)
- [ConstrainedPDMP.jl](https://github.com/tlienart/ConstrainedPDMP.jl)

## Installation

To install this package, run the command 

```julia
add PiecewiseDeterministicMarkovProcesses
```

## References
- R. Veltz, [A new twist for the simulation of hybrid systems using the true jump method](https://arxiv.org/abs/1504.06873), arXiv preprint, 2015.
- A. Drogoul and R. Veltz [Hopf bifurcation in a nonlocal nonlinear transport equation stemming from stochastic neural dynamics](https://aip.scitation.org/doi/abs/10.1063/1.4976510), Chaos: An Interdisciplinary Journal of Nonlinear Science, 27(2), 2017
- Aymard, Campillo, and Veltz, [Mean-Field Limit of Interacting 2D Nonlinear Stochastic Spiking Neurons](https://arxiv.org/abs/1906.10232), arXiv preprint, 2019.






