# """
# This is a wrapper implementing the change of variable method to simulate the PDMP.
# see https://arxiv.org/abs/1504.06873
# """

"
Function copied from Gillespie.jl

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

	# the first x is a dummy variable
	const sr = user_data[2](x, x, user_data[3], t, user_data[4], true)::Float64
	@assert sr > 0.0 "Total rate must be positive"

	const isr = min(1.0e9,1.0 / sr)
	user_data[1](xdot, x, user_data[3], t, user_data[4])
	const ly = length(xdot)
	for i in 1:ly
		xdot[i] = xdot[i] * isr
	end
	xdot[end] = isr
	return Sundials.CV_SUCCESS
end

function f_CHV!{T}(F::Base.Callable,R::Base.Callable,t::Float64, x::Vector{Float64}, xdot::Vector{Float64}, xd::Array{Int64,2}, parms::Vector{T})
	# used for the exact method
	const sr = R(x,x,xd,t,parms,true)
	@assert sr > 0.0 "Total rate must be positive"
	const isr = min(1.0e9,1.0 / sr)
	F(xdot,x,xd,t,parms)
	xdot[end] = 1.0
	const ly = length(xdot)
	scale!(xdot, isr)
	nothing
end

"""
This function performs a pdmp simulation using the Change of Variable (CHV) method see https://arxiv.org/abs/1504.06873.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F!** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
"""

function chv!{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,ind_save_d=-1:1,ind_save_c=-1:1)
	@assert ode in [:cvode,:lsoda]
	# it is faster to pre-allocate arrays and fill it at run time
	n_max += 1 #to hold initial vector
	nsteps = 1
	npoints = 2 # number of points for ODE integration

	# booleans to know if we save data
	save_c = true
	save_d = true
	# permutation to choose randomly a given number of data

	# if ind_save_d[1]==0
# 		save_d = false
# 	elseif ind_save_d[1]==-1
# 		ind_save_d = 1:length(xd0)
# 	end
#
# 	if ind_save_c[1]==0
# 		save_c = false
# 	elseif ind_save_c[1]==-1
# 		ind_save_c = 1:length(xc0)
# 	end

	# Args
	args = pdmpArgs(xc0,xd0,F,R,DX,nu,parms,tf)
	if verbose println("--> Args saved!") end

	# Set up initial variables
	t::Float64 = ti
	xc0     = reshape(xc0,1,length(xc0))
	X0      = vec([xc0 t])
	xd0     = reshape(xd0,1,length(xd0))
	Xd      = deepcopy(xd0)
	deltaxc = copy(nu[1,:]) # declare this variable
	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)#vector of rates

	# arrays for storing history, pre-allocate storage
	t_hist  = Array{Float64}(n_max)
	xc_hist = Array{Float64}(length(xc0), n_max)
	xd_hist = Array{Int64}(length(xd0), n_max)
	res_ode = Array{Float64,2}


	# initialise arrays
	t_hist[nsteps] = t
    xc_hist[:,nsteps] = copy(xc0)
    xd_hist[:,nsteps] = copy(Xd)
	nsteps += 1

	# Main loop
	termination_status = "finaltime"
	while (t < tf) && (nsteps < n_max)

		dt = -log(rand())
		if verbose println("--> t = ",t," - dt = ",dt, ",nstep =  ",nsteps) end
		if ode==:cvode
			res_ode = Sundials.cvode((t,x,xdot)->f_CHV!(F,R,t,x,xdot,Xd,parms), X0, [0.0, dt], abstol = 1e-9, reltol = 1e-7)
		elseif ode==:lsoda
			res_ode = LSODA.lsoda((t,x,xdot,data)->f_CHV!(F,R,t,x,xdot,Xd,parms), X0, [0.0, dt], abstol = 1e-9, reltol = 1e-7)
		end
		if verbose println("--> ode solve is done!") end
		X0 .= res_ode[end,:]
		R(rate,X0[1:end-1],Xd,X0[end],parms, false)
		# jump time:
		t = res_ode[end,end]
		if (t < tf)
			# Update event
			ev = pfsample(rate,sum(rate),numpf)
			deltaxd = nu[ev,:]
			# Xd = Xd .+ deltaxd
			Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			DX(X0,Xd,X0[end],parms,ev)

			if verbose println("--> Which reaction? => ",ev) end
			# save state
			t_hist[nsteps] = t

            @inbounds for ii in eachindex(xc0)
				xc_hist[ii,nsteps] = X0[ii]
            end
            @inbounds for ii in eachindex(Xd)
				xd_hist[ii,nsteps] = Xd[ii]
            end
			# save_c && (xc_hist[:,nsteps] .= X0[ind_save_c])
			# save_d && (xd_hist[:,nsteps] .= Xd[ind_save_d])
		else
			if ode==:cvode
				res_ode = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd,t,parms), X0[1:end-1], [t_hist[end-1], tf], abstol = 1e-9, reltol = 1e-7)
			elseif ode==:lsoda
				res_ode = LSODA.lsoda((t,x,xdot,data)->F(xdot,x,Xd,t,parms), X0[1:end-1], [t_hist[end-1], tf], abstol = 1e-9, reltol = 1e-7)
			end
			t = tf

			# save state
			t_hist[nsteps] = t
			# xc_hist[:,nsteps] = copy(vec(res_ode[end,:]))
			# xd_hist[:,nsteps] = copy(Xd)
			# save_c && (xc_hist[:,nsteps] .= X0[ind_save_c])
			# save_d && (xd_hist[:,nsteps] .= Xd[ind_save_d])
            @inbounds for ii in eachindex(xc0)
				xc_hist[ii,nsteps] = X0[ii]
            end
            @inbounds for ii in eachindex(Xd)
				xd_hist[ii,nsteps] = Xd[ii]
            end
		end
		nsteps += 1
	end
	if verbose println("-->Done") end
	stats = pdmpStats(termination_status,nsteps)
	if verbose println("--> xc = ",xd_hist[:,1:nsteps-1]) end
	return pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
end


"""
This function performs a pdmp simulation using the Change of Variable (CHV) method, see https://arxiv.org/abs/1504.06873. Its use of Sundials solver is optimized in term of memory consumption. It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda
"""
function chv_optim!{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1}, F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64}, parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode)
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
	xc0     = reshape(xc0,1,length(xc0))
	X0      = vec([xc0 t])
	xd0     = reshape(xd0,1,length(xd0))
	Xd      = deepcopy(xd0)
	deltaxc = copy(nu[1,:]) #declare this variable
	numpf   = size(nu,1) #number of reactions
	rate    = zeros(numpf)#vector of rates

	# arrays for storing history, pre-allocate storage
	t_hist  = Array{Float64}( n_max)
	xc_hist = Array{Float64}(length(xc0), n_max)
	xd_hist = Array{Int64}(   length(xd0), n_max)
	res_ode = zeros(2, length(X0))

	# initialise arrays
	t_hist[nsteps]    = t
	xc_hist[:,nsteps] = copy(xc0)
	xd_hist[:,nsteps] = copy(xd0)
	nsteps += 1

	# Main loop
	termination_status = "finaltime"

	# save ODE context, reduces allocation of memory
	if ode==:cvode
		ctx = cvode_ctx(F,R,Xd,parms, X0, [0.0, 1.0], abstol = 1e-9, reltol = 1e-7)
	else
		
		if nsteps == 2
			println(" --> LSODA solve #",nsteps,", X0 = ", X0)
			ctx, res_ode = LSODA.lsoda((t,x,xdot,data)->f_CHV(F,R,t,x,xdot,Xd,parms), X0, [0.0, dt], abstol = 1e-9, reltol = 1e-7)
			X0 = vec(res_ode[end,:])
			dt_lsoda += dt
			println(" ----> res_ode = ", res_ode, ", neq = ",ctx)
		else
			@assert 1==0
			if nsteps == 2
				println(" --> LSODA solve #",nsteps,", X0 = ", X0)
				res_ode = LSODA.lsoda((t,x,xdot,data)->f_CHV!(F,R,t,x,xdot,Xd,parms), X0, [0.0, dt], abstol = 1e-9, reltol = 1e-7)
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

		R(rate,X0[1:end-1],Xd,X0[end],parms, false)

		# Update time
		t = X0[end] #t = res_ode[end,end]
		# @assert t == X0[end]
		# Update event
		if (t < tf)
			ev = pfsample(rate,sum(rate),numpf)
			deltaxd = nu[ev,:]

			# Xd = Xd .+ deltaxd
			Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

			# Xc = Xc .+ deltaxc
			DX(X0,Xd,X0[end],parms,ev)

			if verbose println(" --> Which reaction? => ",ev) end

			# save state
			t_hist[nsteps] = t

			# copy cols: faster, cf. performance tips in JuliaLang
			xc_hist[:,nsteps] = X0[1:end-1]
			xd_hist[:,nsteps] = Xd
		else
			if ode==:cvode
				res_ode = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd ,t,parms), X0[1:end-1], [t_hist[end-1], tf], abstol = 1e-8, reltol = 1e-7)
			end
			t = tf

			# save state
			t_hist[nsteps] = t
			xc_hist[:,nsteps] = copy(vec(res_ode[end,:]))
			xd_hist[:,nsteps] = copy(Xd)
		end
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
	return pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
end
