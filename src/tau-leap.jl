using Distributions

"""
tauleap


This function performs a simulation using the midpoint tau-leap method.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. The returned vector has components. If sum_rate is `False`, one must return rate_vector, bound_ where bound_ is a bound on the total rate vector. In the case sum_rate is `True`, one must return total_rate,bound_ where total_rate is a `Float64` that is the sum of the rates.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
- **dt**: stepsize for the tau-leap method
"""
function tauleap(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Function,R::Function,DX::Function,nu::AbstractArray{Int64},parms::Vector{T},ti::Float64, tf::Float64;verbose::Bool = false,ode = :cvode,dt = 0.1,algo=:tauleap) where T
    # it is faster to pre-allocate arrays and fill it at run time
    @assert algo in [:tauleap,:binomial_tauleap]
    n_max += 1 #to hold initial vector
    nsteps = 1
    # warn("-->Jump on the continuous variable not considered!")

    step! = tau_leap_step!
    if algo==:binomial_tauleap
        step! = binomial_tau_leap_step!
    end


    # Args
    args = pdmpArgs(xc0,xd0,F,R,F,nu,parms,tf)
    if verbose println("--> Args saved!") end

    # Set up initial variables
    t::Float64 = ti
    xc0 = reshape(xc0,1,length(xc0))
    X0  = vec(xc0)
    # to hold the vector field
    dX0  = zeros(X0)
    xd0 = reshape(xd0,1,length(xd0))
    Xd  = deepcopy(xd0)
    deltaxc = copy(nu[1,:]) #declare this variable
    deltaxd = copy(nu[1,:]) # declare this variable, variable to hold discrete jump
    numpf   = size(nu,1)    # number of reactions
    rate    = zeros(numpf)  #vector of rates

    # arrays for storing history, pre-allocate storage
    t_hist  = Array{Float64}(n_max)
    xc_hist = Array{Float64}(length(xc0), n_max)
    xd_hist = Array{Int64}(length(xd0), n_max)
    res_ode = Array{Float64,2}

    # initialise arrays
    t_hist[nsteps] = t
    xc_hist[:,nsteps] = copy(xc0)
    xd_hist[:,nsteps] = copy(xd0)

    # Main loop
    termination_status = "finaltime"

    tp = [0.,0.]
    R(rate,X0,Xd,t,parms,false)
    while (t < tf) && (nsteps < n_max)
        # println("#################\nX  = $X0\nXd = $Xd\n,$nsteps,t=$t")
        @assert sum(Xd.<0)==0 "You have negative discrete component, $Xd, $t, $nsteps"

        verbose && println("--> step : ",nsteps," / ",n_max )
        t += dt
        nsteps += 1

        # Poisson approximation for discrete part
        R(rate,X0,Xd,t,parms,false)
        step!(t,dt,rate,nu,X0,Xd,DX,parms)

        # Euler scheme for ODE part
        F(dX0,X0,Xd,t,parms)
        @. X0 .= X0 + dt * dX0

        t_hist[nsteps] = t
        xc_hist[:,nsteps] = copy(X0)
        xd_hist[:,nsteps] = copy(Xd)
    end

    # nsteps -=1

    verbose && println("-->Done")
    stats = pdmpStats(termination_status,nsteps)

    verbose && println("--> xc = ",xd_hist[:,1:nsteps])
    result = pdmpResult(t_hist[1:nsteps],xc_hist[:,1:nsteps],xd_hist[:,1:nsteps],stats,args)

    return(result)
end

function tau_leap_step!(t,dt,rates,nu,Xc,Xd,DX,p)
    nb_mol = 0
    for i=1:length(rates)
        if rates[i]>0
            nb_mol = rand(Poisson(dt * rates[i])) # tau-leap method
            LinearAlgebra.BLAS.axpy!(nb_mol, nu[i,:], Xd)
            DX(Xc,Xd,t,p,i)
        end
    end
end

function binomial_tau_leap_step!(dt,rates,nu,Xc,Xd,DX,p)
    nb_mol = 0
    for i=1:length(ppf)
        if rates[i]>0
            # nb_mol = rand(Poisson(dt * ppf[i])) # tau-leap method
            ind = findall(nu[i,:] .== -1)
            if length(ind) == 0
                nb_mol = rand(Poisson(dt * ppf[i])) # tau-leap method
            else
                # @show ind,nsteps,t,nu[i,:],i, dt * ppf[i],Xd[ind[1]]
                println("--> dt = $dt, $(dt * ppf[i]/Xd[ind[1]])")
                nb_mol = rand(Distributions.Binomial(Xd[ind[1]],dt * ppf[i]/Xd[ind[1]]))
            end
            LinearAlgebra.BLAS.axpy!(nb_mol, nu[i,:], Xd)
			DX(Xc,Xd,t,p,i)
        end
    end
end
