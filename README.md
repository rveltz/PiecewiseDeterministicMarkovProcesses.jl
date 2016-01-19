# PDMP 

[![Build Status](https://travis-ci.org/sdwfrost/PDMP.jl.svg?branch=master)](https://travis-ci.org/sdwfrost/PDMP.jl)

This is a joint work of (Veltz Romain, Frost Simon).

PDMP.jl is a Julia package that allows simulation of Piecewise Deterministic Markov Processes (PDMP) also called hybrid systems.

It is based on an implementation of the [True Jump Method method](http://arxiv.org/abs/1504.06873) for performing stochastic simulations of PDMP. (A rejection method will be added soon.)

We briefly recall facts about a simple class of PDMPs. They are decribed by a couple $(x_c,x_d)$ where $x_c$ is solution of the differential equation $\frac{dx_c}{dt} = F(x_c,x_d,t)$. The second component $x_d$ is a jump process with rates $R(x_c,x_d,t)$. At each jump of $x_d$, a jump can be added to the continuous variable $x_c$ too.

##Installation
To install this (unregistered) package, run the command 	```Pkg.clone(https://github.com/sdwfrost/PDMP.jl.git)```

##Examples
See the [examples directory](https://github.com/sdwfrost/PDMP.jl/tree/master/examples).

A simple example of a TCP process is given below:

```julia
using PDMP

function F_tcp(xcdot, xc, xd, t, parms )
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = xc[1]
  else
    xcdot[1] = -xc[1]
  end
  nothing
end

function R_tcp(xc, xd, t, parms, sum_rate::Bool)
  # rate function
  if sum_rate==false
    return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]]/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1)
  else
    return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1
  end
end

function Delta_xc_tcp(xc, xd, t, parms, ind_reaction::Int64)
	# jump on the continuous variable
  return true #in this example, no jump
end


xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = [0.]
tf = 2000.

dummy =  PDMP.chv(2,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)
result =  @time PDMP.chv(2000,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)

println("#jumps = ", length(result.time))

using GR
GR.inline()
ind = find(result.time.<49)
GR.plot(result.time[ind],result.xc[1,:][ind],"k",result.time[ind],result.xd[1,:][ind],"r",title = string("#Jumps = ",length(result.time)))
```

![SIR](examples/tcp.png)

Passing functions as arguments in Julia (currently) incurs a performance penalty. One can circumvent this by passing an immutable object, with ```call``` overloaded. An example of this approach is given [here](https://).
