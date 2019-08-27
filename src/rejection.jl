struct RejectionExact <: AbstractRejectionExact end

"""
This function performs a simulation using the rejection method.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R!** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. The returned vector has components. If sum_rate is `False`, one must return rate_vector, bound_ where bound_ is a bound on the total rate vector. In the case sum_rate is `True`, one must return total_rate,bound_ where total_rate is a `Float64` that is the sum of the rates. `R!(rate,xc,xd,t,parms,sum_rate)` returns `Float64,Float64`. In case `sum_rate = true`, you are not allowed to modify the first argument e.g. `rate`
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : data for the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
"""
function solve(problem::PDMPProblem, algo::Rejection{Tode}; verbose::Bool = false, save_rejected=false, ind_save_d=-1:1, ind_save_c=-1:1, n_jumps = Inf64, reltol = 1e-7, abstol = 1e-9) where {Tode <: Symbol}
	verbose && println("#"^30)
	ode = algo.ode
	@assert ode in [:cvode, :lsoda, :adams, :bdf]
	verbose && printstyled(color=:red,"--> Start rejection method\n")

	# define the ODE flow
	if ode == :cvode || ode == :bdf
		Flow = (X0_,Xd_,tp_)->Sundials.cvode(  (tt,x,xdot)->problem.caract.F(xdot,x,Xd,tt,problem.caract.parms), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :BDF)
	elseif	ode == :adams
		Flow = (X0_,Xd_,tp_)->Sundials.cvode(  (tt,x,xdot)->problem.caract.F(xdot,x,Xd,tt,problem.caract.parms), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :Adams)
	elseif ode == :lsoda
		Flow = (X0_,Xd_,tp_)->LSODA.lsoda((tt,x,xdot,data)->problem.caract.F(xdot,x,Xd,tt,problem.caract.parms), X0_, tp_, abstol = abstol, reltol = reltol)
	end

	ti, tf = problem.interval
	# it is faster to pre-allocate arrays and fill it at run time
	n_jumps  += 1 #to hold initial vector
	nsteps  = 1
	npoints = 2 # number of points for ODE integration

	xc0 = problem.caract.xc
	xd0 = problem.caract.xd

	# Set up initial variables
	t = ti
	X0, _, Xd, t_hist, xc_hist, xd_hist, res_ode = allocate_arrays(ti, xc0, xd0, n_jumps, true)

	deltaxd = copy(problem.caract.pdmpjump.nu[1, :]) # declare this variable
	numpf   = size(problem.caract.pdmpjump.nu,1)	 # number of reactions
	rate	= zeros(numpf)  # vector of rates
	tp = [ti, tf]		   # vector to hold the time interval over which to integrate the flow

	#variables for rejection algorithm
	reject = true
	lambda_star = 0.0 # this is the bound for the rejection method
	ppf = problem.caract.R(rate, X0, Xd, t, problem.caract.parms, true)

	# @assert ppf[2] == R(rate,X0+0.1265987*cumsum(ones(length(X0))),Xd,t+0.124686489,parms,true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"
	# rate *= 0;ppf = R(rate,X0,Xd,t,parms,true)
	# @assert sum(rate) == 0 "You cannot modify the first argument of your rate function when sum_rate = true"

	δt = problem.simjptimes.tstop_extended

	while (t < tf) && (nsteps < n_jumps)
		verbose && println("--> step : ",nsteps," / ", n_jumps)
		reject = true
		while reject && (nsteps < n_jumps)
			tp .= [t, min(tf, t + δt / ppf[2]) ] #mettre un lambda_star?
			res_ode .= Flow(X0, Xd, tp)

			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end, ii]
			end
			verbose && println("----> δt = ", δt, ", t∈", tp, ", dt = ", tp[2]-tp[1], ", xc = ", X0)

			t = tp[end]
			ppf = problem.caract.R(rate, X0, Xd, t, problem.caract.parms, true)
			@assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf"
			if t == tf
				reject = false
			else
				reject = rand() < 1 - ppf[1] / ppf[2]
			end
			δt = -log(rand())
		end

		# there is a jump!
		ppf = problem.caract.R(rate,X0,Xd,t,problem.caract.parms,false)

		if (t < tf)
			verbose && println("----> Jump!, ratio = ", ppf[1] / ppf[2], ", xd = ", Xd)
			# make a jump
			ev = pfsample(rate, sum(rate), numpf)
			deltaxd .= problem.caract.pdmpjump.nu[ev,:]

			# Xd = Xd .+ deltaxd
			LinearAlgebra.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			problem.caract.pdmpjump.Delta(X0, Xd, X0[end], problem.caract.parms, ev)
		end

		nsteps += 1
		t_hist[nsteps] = t
		@inbounds for ii in eachindex(X0)
			xc_hist[ii,nsteps] = X0[ii]
		end
		@inbounds for ii in eachindex(Xd)
			xd_hist[ii,nsteps] = Xd[ii]
		end
	end
	if verbose println("-->Done") end
	if verbose println("--> xd = ",xd_hist[:,1:nsteps]) end
	result = PDMPResult(t_hist[1:nsteps], xc_hist[:,1:nsteps], xd_hist[:,1:nsteps], Float64[])
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
- **parms** : data for the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"""
function solve(problem::PDMPProblem, algo::Talgo; verbose::Bool = false, save_rejected = false, ind_save_d=-1:1, ind_save_c=-1:1, n_jumps = Inf64, xd_jump::Bool=true) where {Talgo <: AbstractRejectionExact}
	ti, tf = problem.interval
	# it is faster to pre-allocate arrays and fill it at run time
	n_jumps += 1 #to hold initial vector
	nsteps = 1
	npoints = 2 # number of points for ODE integration

	# Set up initial variables
	t = ti
	X0, _, Xd, t_hist, xc_hist, xd_hist, res_ode = allocate_arrays(ti, xc0, xd0, n_jumps, true)

	deltaxd = copy(problem.caract.pdmpjump.nu[1, :]) # declare this variable
	numpf   = size(problem.caract.pdmpjump.nu,1)	 # number of reactions
	rate_vector = zeros(numpf)#vector of rates
	tp = [0., 1.]

	reject = true
	nb_rejet = 0
	lambda_star = 0.0 # this is the bound for the rejection method
	tp = [0.,0.]
	lambda_star = problem.caract.R(rate_vector, X0, Xd, t, problem.caract.parms, true)[2]

	@assert lambda_star == problem.caract.R(rate_vector,X0,Xd,t+rand(),problem.caract.parms,true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"

	δt = problem.simjptimes.tstop_extended

	while (t < tf) && (nsteps < n_jumps)
		if verbose println("--> step : $njumps, / $n_jumps, #reject = $nsteps" ) end
		reject = true
		nsteps = 1
		while (reject) && (nsteps < 10^6) && (t < tf)
			tp = [t, min(tf, t + δt / lambda_star)] 		# mettre un lambda_star?
			problem.caract.F(res_ode, X0, Xd, tp, parms) 	# we evolve the flow inplace

			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end, ii]
			end

			t = tp[end]
			ppf = problem.caract.R(rate_vector, X0, Xd, t, problem.caract.parms, true)
			@assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf"
			if t == tf
				reject = false
			else
				reject = rand() < 1 - ppf[1] / ppf[2]
			end
			δt = -log(rand())
			nb_rejet += 1
		end

		# there is a jump!
		ppf = problem.caract.R(rate_vector, X0, Xd, t, problem.caract.parms, false)

		if (t < tf)
			verbose && println("----> Jump!, ratio = ", ppf[1] / ppf[2], ", xd = ", Xd)
			# make a jump
			ev = pfsample(rate_vector, sum(rate_vector), numpf)
			deltaxd .= problem.caract.pdmpjump.nu[ev,:]

			# Xd = Xd .+ deltaxd
			LinearAlgebra.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			problem.caract.pdmpjump.Delta(X0, Xd, t, problem.caract.parms, ev)
		end

		nsteps += 1
		t_hist[nsteps] = t
		@inbounds for ii in eachindex(X0)
			xc_hist[ii,nsteps] = X0[ii]
		end
		@inbounds for ii in eachindex(Xd)
			xd_hist[ii,nsteps] = Xd[ii]
		end
	end
	println("njumps = ",nsteps, " / rejections = ", nb_rejet)
	if verbose println("-->Done") end

	# if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
	result = PDMPResult(t_hist[1:njumps],xc_hist[:,1:njumps],xd_hist[:,1:njumps],Float64[])
	return(result)
end
