function Phi_dummy(out::Array{Float64,2}, xc::Vector{Float64},xd,t::Array{Float64},parms)
    # vector field used for the continuous variable
    # trivial dynamics
    out[1,:] .= xc
    out[2,:] .= xc
    nothing
end


"""
This function performs a simulation using the rejection method.
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
"""
function rejection{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Function,R::Function,DX::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode = :cvode,save_rejected=false)
    @assert ode in [:cvode,:lsoda]
    # it is faster to pre-allocate arrays and fill it at run time
    n_max += 1 #to hold initial vector
    nsteps = 1
    npoints = 2 # number of points for ODE integration
	njumps = 1

    # Args
    args = pdmpArgs(xc0,xd0,F,R,DX,nu,parms,tf)
    if verbose println("--> Args saved!") end

    # Set up initial variables
    t::Float64 = ti
    xc0 = reshape(xc0,1,length(xc0))
    X0  = vec(xc0)
    xd0 = reshape(xd0,1,length(xd0))
    Xd  = deepcopy(xd0)
    deltaxc = copy(nu[1,:]) #declare this variable

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

    reject = true
    lambda_star = 0.0 # this is the bound for the rejection method
    tp = [0.,0.]
    ppf = R(X0,Xd,t,parms,true)
    while (t < tf) && (nsteps < n_max)
        if verbose println("--> step : ",nsteps," / ",n_max ) end
        reject = true
        while (reject) && (nsteps < n_max)
            tp = [t, min(tf, t - log(rand())/ppf[2]) ] #mettre un lambda_star?
            if verbose println("----> tspan : ",tp ) end
            if ode==:cvode
                res_ode = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd,t,parms), X0, tp, abstol = 1e-9, reltol = 1e-7)
            elseif ode==:lsoda
                res_ode = LSODA.lsoda((t,x,xdot,data)->F(xdot,x,Xd,t,parms), X0, tp, abstol = 1e-9, reltol = 1e-7)
            end
            X0 = vec(res_ode[end,:])
            t = tp[end]
            ppf = R(X0,Xd,t,parms,true)
            @assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf, t=$t"
            if t == tf
                reject = false
            else
                reject = rand() < (1. - ppf[1] / ppf[2])
            end
            if save_rejected
                nsteps += 1
                t_hist[nsteps] = t
                xc_hist[:,nsteps] = copy(X0[1:end])
                xd_hist[:,nsteps] = copy(Xd)
            end
			njumps +=1

        end

        if save_rejected==false
            nsteps += 1
            t_hist[nsteps] = t
            xc_hist[:,nsteps] = copy(X0[1:end])
            xd_hist[:,nsteps] = copy(Xd)
        end
        # there is a jump!
        ppf = R(X0,Xd,t,parms,false)
        pf = StatsBase.Weights(convert(Array{Float64,1},ppf[1])) #this is to ease sampling

        if (t < tf)
            # make a jump
            ev = Distributions.sample(pf)
            deltaxd = nu[ev,:]

            # Xd = Xd .+ deltaxd
            Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

            # Xc = Xc .+ deltaxc
            DX(X0,Xd,X0[end],parms,ev)
        end
    end
    println("njumps = ",nsteps," / rejections = ", njumps, ", lambda_star = ",lambda_star)
    if verbose println("-->Done") end
    stats = pdmpStats(termination_status,nsteps)
    if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
    result = pdmpResult(t_hist[1:nsteps],xc_hist[:,1:nsteps],xd_hist[:,1:nsteps],stats,args)
    return(result)
end

"""

rejection_exact

This function performs a simulation using the rejection method when the flow **is known analytically**.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **Phi!** : a `Function` or a callable type, which itself takes 6 arguments to represent the vector field; rate a `Vector` of `Float64` representing the **flow** of the vector which needs to be filled with values of the rates, xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system, sum_of_rate a `Bool` stating if the function must return the total rate.
- **R!** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. The returned vector has components. If sum_rate is `False`, one must return rate_vector, bound_ where bound_ is a bound on the total rate vector. In the case sum_rate is `True`, one must return total_rate,bound_ where total_rate is a `Float64` that is the sum of the rates. In any case, the function must return a couple (total_rates, bound) where bound is a bound for the total rate.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"""
function rejection_exact{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},Phi::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false, xd_jump::Bool=true)
    # it is faster to pre-allocate arrays and fill it at run time
    n_max += 1 #to hold initial vector
    const nsteps = 1
    const npoints = 2 # number of points for ODE integration
    const njumps = 1

    # Args
    args = pdmpArgs(xc0,xd0,Phi,R,DX,nu,parms,tf)

    # Set up initial variables
    t::Float64 = ti
    xc0 = reshape(xc0,1,length(xc0))
    X0  = vec(xc0)
    xd0 = reshape(xd0,1,length(xd0))
    Xd  = deepcopy(xd0)
    deltaxc = copy(nu[1,:]) #declare this variable

    # arrays for storing history, pre-allocate storage
    t_hist  = Array{Float64}( n_max)
    xc_hist = Array{Float64}(length(xc0), n_max)
    xd_hist = Array{Int64}(length(xd0), n_max)
    res_ode = Array{Float64}(2,length(xc0))
    rate_vector = zeros(length(nu[:,1]))


    # initialise arrays
    t_hist[nsteps] = t
    xc_hist[:,nsteps] = copy(xc0)
    xd_hist[:,nsteps] = copy(xd0)

    # Main loop
    termination_status = "finaltime"

    const reject = true
    const nb_rejet::Int = 0
    const lambda_star = 0.0 # this is the bound for the rejection method
    tp = [0.,0.]
    lambda_star = R(rate_vector,X0,Xd,t,parms,true)[2]

    t_hist[njumps] = t
    xc_hist[:,njumps] = copy(X0[1:end])
    xd_hist[:,njumps] = copy(Xd)


    while (t < tf) && (njumps < n_max)
        verbose && println("--> step : $njumps, / $n_max, #reject = $nsteps" , ", lambda_star = ",lambda_star, ", t=$t")
        reject = true
        nsteps = 1
        while (reject) && (nsteps < 10^6) && (t < tf)
            tp = [t, t - log(rand())/lambda_star ]		# mettre un lambda_star?
            Phi(res_ode, X0, Xd, tp, parms) 				# we evolve the flow inplace
            X0 = vec(res_ode[end,:])
            t = tp[end]
            ppf = R(rate_vector, X0, Xd, t, parms, true) 	# we don't want the full rate vector, just the sum of rates
            # @show ppf[1]
			verbose && @show X0, tp,ppf
            @assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong at t = $t"
            reject = rand() <  (1. - ppf[1] / ppf[2])
            nsteps += 1
        end
        # keep track of nb of rejections
        nb_rejet += nsteps

        @assert(nsteps <= 10^6,"Error, too many rejections!!")
        njumps += 1
        t_hist[njumps] = t
        xc_hist[:,njumps] = copy(X0[1:end])
        xd_hist[:,njumps] = copy(Xd)

        # there is a jump!
        lambda_star = R(rate_vector,X0,Xd,t,parms,false)[2]
        verbose && println("----> rate = $rate_vector,Xc = $X0" )
        pf = StatsBase.Weights(convert(Array{Float64,1},rate_vector)) #this is to ease sampling
        @assert(pf.sum>0,"Error, rate vector is null for some reason")
        if (t < tf)
            # make a jump
            ev = Distributions.sample(pf)
            verbose && println("----> reaction = $ev" )

            if xd_jump
                deltaxd = nu[ev,:]
                # Xd = Xd .+ deltaxd
                Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)
            end
            # Xc = Xc .+ deltaxc
            DX(X0,Xd,t,parms,ev)
        end
    end
    println("njumps = ",njumps," / rejections = ", nb_rejet, ", lambda_star = ",lambda_star)
    verbose && println("-->Done")
    stats = pdmpStats(termination_status,nsteps)
    # if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
    result = pdmpResult(t_hist[1:njumps],xc_hist[:,1:njumps],xd_hist[:,1:njumps],stats,args)
    return(result)
end


rejection_exact{T}(n_max::Int64,xd0::Array{Int64,1},R::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false, xd_jump::Bool=true) = PDMP.rejection_exact(n_max,[0.],xd0,Phi_dummy,R,Delta_dummy,nu,parms,ti, tf,verbose, xd_jump)
