# PDMP.jl 


PDMP.jl is a Julia package that allows simulation of *Piecewise Deterministic Markov Processes* (PDMP); these encompass hybrid systems and jump processes, comprised of continuous and discrete components, as well as processes with time-varying rates. The aim of the package is to provide methods for the simulation of these processes that are "exact" up to the ODE integrator.

We provide several methods for the simulation:

- a recent trick, called **CHV**, explained in [paper-2015](http://arxiv.org/abs/1504.06873) which allows to implement the **True Jump Method** without the need to use event detection schemes for the ODE integrator. These event detections can be quite unstable as explained in [paper-2015](http://arxiv.org/abs/1504.06873) and CHV provide a solution to this problem.
- **rejection methods** for which the user is asked to provide a bound on the reaction rates. These last methods are the most "exact" but not the fastest if the reaction rate bound is not tight. In case the flow is known analytically, a method is also provided.


These methods require solving stiff ODEs (for CHV ) in an efficient manner. [```Sundials.jl```](https://github.com/JuliaLang/Sundials.jl) and [```LSODA.jl```](https://github.com/rveltz/LSODA.jl) are used, but other solvers could be easily added. (See [stiff ode solvers](http://lh3lh3.users.sourceforge.net/solveode.shtml)).

We briefly recall facts about a simple class of PDMPs. They are described by a couple $(x_c,x_d)$ where $x_c$ is solution of the differential equation $\frac{dx_c}{dt} = F(x_c,x_d,t)$. The second component $x_d$ is a jump process with rates $R(x_c,x_d,t)$. At each jump of $x_d$, a jump can also be added to the continuous variable $x_c$.


## Installation

To install this (unregistered) package, run the command 

```julia
Pkg.clone("https://github.com/rveltz/PDMP.jl.git")
```

## Basic example

See the [examples directory](https://github.com/rveltz/PDMP.jl/tree/master/examples).

A simple example of a TCP process is given below.More precisely, we look at the following process of switching dynamics where X(t) = $(x_c(t), x_d(t)) \in\mathbb R\times\lbrace-1,1\rbrace$. In between jumps, $x_c$ evolves according to $\dot x_c(t) = x_d(t)$. 

We first need to load the library.
```julia
using PDMP
```
We then define a function that encodes the dynamics in between jumps. We need to provide the vector field of the ODE with a function. Hence, we need to define a function that given continuous state $x_c$ and discrete state $x_d$ at time $t$ return the vector field. In addition some parameters can be passed with the variable `parms`.

```julia
function F_tcp(xc, xd, t, parms)
  # vector field used for the continuous variable
  return xd[1]
end
```

Let's consider a stochastic process with following transitions.


| Transition | Rate |
|---|---|---|
|$x_d\to x_d-2$ if $x_d>0$ | 1 |
|$x_d\to x_d+2$ if $x_d<0$ | 1 |

This is encoded in the following function


```julia
function R_tcp(xc, xd, t, parms, sum_rate::Bool)
  # rate function for each transition
  # in this case,  the transitions are xd->xd+2 or xd->xd-2
  # sum_rate is a boolean which tells R_tcp the type which must be returned:
  # i.e. the sum of the rates or the vector of the rates
  if sum_rate==false
      if xd[1] > 0
          return [0,1]
      else
          return [1,0]
      end
  else
    return 1
  end
end

# initial conditions for the continuous/discrete variables
# initial conditions for the continuous/discrete variables
xc0 = vec([0.05])
xd0 = vec([1])

# matrix of jumps for the discrete variables, analogous to chemical reactions
const nu = reshape([[2];[-2]],2,1)


# parameters
parms = [0.]
tf = 2000.

# compile the program:
dummy =  PDMP.pdmp(2,xc0,xd0,F_tcp,R_tcp,nu,parms,0.0,tf,false)

# compute a trajectory
result =  @time PDMP.pdmp(100,xc0,xd0,F_tcp,R_tcp,nu,parms,0.0,tf,false)

# plotting
using Plots
Plots.plot(result.time, result.xd[1,:],line=:step,title = string("#Jumps = ",length(result.time)),label="Xd")
Plots.plot(result.time, result.xc',title = string("#Jumps = ",length(result.time)),label="Xc")
```

# Application programming interface

## Functions


```@docs
chv!
```


```@docs
rejection_exact
```
