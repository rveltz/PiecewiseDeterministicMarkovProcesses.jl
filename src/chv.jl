# """
# This is a wrapper implementing the change of variable method to simulate the PDMP.
# see https://arxiv.org/abs/1504.06873
# """

"
Function copied from Gillespie.jl and StatsBase

This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:
- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.
"
function pfsample(w::Array{Float64,1},s::Float64,n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end

"""
This is a wrapper implementing the change of variable method to simulate the PDMP.
This wrapper is meant to be called by Sundials.CVode
see https://arxiv.org/abs/1504.06873
"""
function cvode_ode_wrapper(t, x_nv, xdot_nv, user_data)
	# Reminder: user_data = [F R Xd params]
	x    = convert(Vector, x_nv)
	xdot = convert(Vector, xdot_nv)

	tau = x[end]
	# the first x is a dummy variable, it will be seen as the rate vector but it
	# must not be modified
	sr = user_data[2](user_data[5], x, user_data[3], tau, user_data[4], true)[1]::Float64
	@assert sr > 0.0 "Total rate must be positive"

	isr = min(1.0e9,1.0 / sr)
	user_data[1](xdot, x, user_data[3], tau, user_data[4])
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] * isr
	end
	xdot[end] = isr
	return Sundials.CV_SUCCESS
end

function f_CHV!(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64}, xd::Vector{Int64}, parms,rate_::Array{Float64})
	# used for the exact method
	# we put [1] to use it in the case of the rejection method as well
	tau = x[end]
	sr = R(rate_,x,xd,tau,parms,true)[1]
	@assert sr > 0.0 "Total rate must be positive"
	isr = min(1.0e9,1.0 / sr)
	F(xdot,x,xd,tau,parms)
	xdot[end] = 1.0
	ly = length(xdot)
	scale!(xdot, isr)
	nothing
end

function f_CHV_Wrap!(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64}, p::DataForODE{T} where T)
	f_CHV!(F,R,t, x, xdot, p.xd, p.parms,p.rate)
end


function Flow_Wrap!(X0,Xd,dt,r,prob::ODEProblem)
	(prob.p).rate .= r
	(prob.p).xd .= Xd
	prob.u0[:] .= X0
	println("flow wrap done")
	sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),save_start=true,save_end=true,save_everystep = false).u
	# res[1,:] .= sol[1]
	# res[2,:] .= sol[2]
    return hcat(sol[1],sol[2])'
	return nothing
end
"""

chv!

This function performs a pdmp simulation using the Change of Variable (CHV) method see https://arxiv.org/abs/1504.06873.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F!** : an inplace `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : an inplace `Function` or a callable type, which itself takes six arguments to represent the rate functions associated to the jumps;rate `Vector` of `Float64` holding the different reaction rates, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **DX** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : data for the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
"""

function chv!(n_max::Int64,xc0::AbstractVector{Float64},xd0::AbstractVector{Int64},F::Function,R::Function,DX::Function,nu::AbstractArray{Int64},parms,ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,ind_save_d=-1:1,ind_save_c=-1:1)
	@assert ode in [:cvode,:lsoda,:OrdinaryDiffEq]
	# it is faster to pre-allocate arrays and fill it at run time
	n_max += 1 #to hold initial vector
	nsteps = 1 #index for the current jump number
	npoints = 2 # number of points for ODE integration

	# permutation to choose randomly a given number of data Args
	args = pdmpArgs(xc0,xd0,F,R,DX,nu,parms,tf)
	if verbose println("--> Args saved!") end

	# Set up initial variables
	t = ti         # initial simulation time
	X0, Xc, Xd, t_hist, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c = allocate_arrays(ti,xc0,xd0,n_max,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	nsteps += 1

	deltaxd = copy(nu[1,:]) # declare this variable, variable to hold discrete jump
	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)  #vector of rates

	# define the ODE flow, this leads to big memory saving
	if ode==:cvode
		Flow = (X0_,Xd_,dt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., dt], abstol = 1e-9, reltol = 1e-7)
	elseif ode==:lsoda
		Flow = (X0_,Xd_,dt,r_)->LSODA.lsoda((tt,x,xdot,data)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., dt], abstol = 1e-9, reltol = 1e-7)
	elseif ode==:OrdinaryDiffEq
		data_ode = DataForODE(parms,Xd,rate)
		prob_CHV = ODEProblem((xdot,x,data,tt)->f_CHV_Wrap!(F,R,tt,x,xdot,data),X0,(0.,1.),data_ode)
		Flow = (X0_,Xd_,dt,r_)-> Flow_Wrap!(X0_,Xd_,dt,r_,prob_CHV)
	end

	# Main loop
	termination_status = "finaltime"
	while (t < tf) && (nsteps < n_max)

		# dt = -log(rand())
		dt = - log(rand())
		verbose && println("**********************\n--> t = ",t," - dt = ",dt, ",nstep =  ",nsteps)

        @show Flow(X0,Xd,dt,rate)
		@show res_ode

		res_ode .= Flow(X0,Xd,dt,rate)

		verbose && println("--> ode solve is done!")

		@inbounds for ii in eachindex(X0)
			X0[ii] = res_ode[end,ii]
		end
		t = res_ode[end,end]

		R(rate,Xc,Xd,t,parms, false)
		# jump time:
		if (t < tf)
			# Update event
			ev = pfsample(rate,sum(rate),numpf)
			deltaxd .= nu[ev,:]
			# Xd = Xd .+ deltaxd
			Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			DX(Xc,Xd,t,parms,ev) #requires allocation!!

			verbose && println("--> Which reaction? => ",ev)
			# save state
			t_hist[nsteps] = t
			save_data(nsteps,X0,Xd,xc_hist,xd_hist,ind_save_d, ind_save_c)

		else
			if ode==:cvode
				res_ode_last =   Sundials.cvode((tt,x,xdot)->F(xdot,x,Xd,tt,parms), X0[1:end-1], [t_hist[nsteps-1], tf], abstol = 1e-9, reltol = 1e-7)
			elseif ode==:lsoda
				res_ode_last = LSODA.lsoda((tt,x,xdot,data)->F(xdot,x,Xd,tt,parms), X0[1:end-1], [t_hist[nsteps-1], tf], abstol = 1e-9, reltol = 1e-7)
			end
			t = tf

			# save state
			t_hist[nsteps] = tf
			@inbounds for ii in eachindex(ind_save_c)
				xc_hist[ii,nsteps] = res_ode_last[end,ind_save_c[ii]]
		    end
		    @inbounds for ii in eachindex(ind_save_d)
				xd_hist[ii,nsteps] = Xd[ind_save_d[ii]]
		    end
		end
		nsteps += 1
	end
	verbose && println("-->Done")
	stats = pdmpStats(termination_status,nsteps)
	verbose && println("--> xc = ",xd_hist[:,1:nsteps-1])
	return pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
end


"""

chv_optim!

This function performs a pdmp simulation using the Change of Variable (CHV) method, see https://arxiv.org/abs/1504.06873. Its use of Sundials solver is optimized in term of memory consumption. It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : data for the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
"""
function chv_optim!(n_max::Int64,xc0::AbstractVector{Float64},xd0::AbstractVector{Int64}, F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::AbstractArray{Int64}, parms,ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,ind_save_d=-1:1,ind_save_c=-1:1)
	@assert ode in [:cvode] string("Sorry, ",ode," is not available for chv_optim yet")
	# it is faster to pre-allocate arrays and fill it at run time
	n_max  += 1 #to hold initial vector
	nsteps  = 1
	npoints = 2 # number of points for ODE integration

	# Args
	args = pdmpArgs(xc0,xd0,F,R,DX,nu,parms,tf)
	if verbose println("--> Args saved!") end

	# Set up initial variables
	t::Float64 = ti
	X0, Xc, Xd, t_hist, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c = allocate_arrays(ti,xc0,xd0,n_max,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	nsteps += 1

	deltaxd = copy(nu[1,:]) # declare this variable
	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)#vector of rates
	nsteps += 1

	# Main loop
	termination_status = "finaltime"

	# save ODE context, reduces allocation of memory
	if ode==:cvode
		ctx = cvode_ctx(F,R,Xd,parms,rate, X0, [0.0, 1.0], abstol = 1e-9, reltol = 1e-7)
	else
		# ctx = LSODA.lsoda_context_t()
		dt_lsoda = 0.
	end
	#   prgs = Progress(n_max, 1)
	while (t < tf) && (nsteps<n_max)
		#     update!(prgs, nsteps)
		dt = -log(rand())
		if verbose println("--> t = ",t," - dt = ",dt) end

		if ode==:cvode
			# println(" --> CVODE solve #",nsteps,", X0 = ", X0)
			cvode_evolve!(res_ode, ctx[1],F,R,Xd,parms, rate, X0, [0.0, dt])
			# println(" ----> res_ode = ", res_ode)
			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end,ii]
			end
		else
			if nsteps == 2
				println(" --> LSODA solve #",nsteps,", X0 = ", X0)
				res_ode = LSODA.lsoda((t,x,xdot,data)->f_CHV(F,R,t,x,xdot,Xd,parms,rate), X0, [0.0, dt], abstol = 1e-9, reltol = 1e-7)
				X0 = vec(res_ode[end,:])
				dt_lsoda += dt
				println(" ----> res_ode = ", res_ode, ", neq = ",ctx)
			else
				println(" --> lsoda_evolve #",nsteps,", X0 = ",X0,", res_ode = ",res_ode,",dt = ", [dt_lsoda, dt_lsoda + dt])
				LSODA.lsoda_evolve!(ctx[1], X0, [dt_lsoda, dt_lsoda + dt])
				dt_lsoda += dt
			end
		end
		if verbose println(" --> ode solve is done!") end

		R(rate,Xc,Xd,t,parms, false)

		# Update time
		t = X0[end] #t = res_ode[end,end]
		# @assert t == X0[end]
		# Update event
		if (t < tf)
			ev = pfsample(rate,sum(rate),numpf)
			deltaxd .= nu[ev,:]
			# Xd = Xd .+ deltaxd
			Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			DX(Xc,Xd,t,parms,ev) #requires allocation!!

			if verbose println(" --> Which reaction? => ",ev) end

			# save state
			t_hist[nsteps] = t

			# copy cols: faster, cf. performance tips in JuliaLang
            save_data(nsteps,X0,Xd,xc_hist,xd_hist,ind_save_d, ind_save_c)
		else
			if ode==:cvode
				res_ode_last = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd ,t,parms), X0[1:end-1], [t_hist[nsteps-1], tf], abstol = 1e-8, reltol = 1e-7)
			end
			t = tf

			# save state
			t_hist[nsteps] = tf
			@inbounds for ii in eachindex(ind_save_c)
				xc_hist[ii,nsteps] = res_ode_last[end,ind_save_c[ii]]
		    end
		    @inbounds for ii in eachindex(ind_save_d)
				xd_hist[ii,nsteps] = Xd[ind_save_d[ii]]
		    end
		end
		nsteps += 1
	end

	if ode==:cvode
		# Sundials.CVodeFree(Ref([ctx]))
		Sundials.empty!(ctx[1])
		Sundials.empty!(ctx[2])
		Sundials.empty!(ctx[3])
	end
	# collect the data
	if verbose println("-->Done") end
	stats = pdmpStats(termination_status,nsteps)
	if verbose println("--> xc = ",xd_hist[:,1:nsteps-1]) end
	if verbose println("--> time = ",t_hist[1:nsteps-1]) end
	if verbose println("--> chv_optim, #jumps = ",length(t_hist[1:nsteps-1])) end
	result = pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
	return(result)
end
