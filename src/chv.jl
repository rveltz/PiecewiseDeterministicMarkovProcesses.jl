function f_CHV!(F::Function,R::Function,t::Float64, x, xdot, xd, parms,rate)
	# used for the exact method
	# we put [1] to use it in the case of the rejection method as well
	tau = x[end]
	sr = R(rate,x,xd,tau,parms,true)[1]
	@assert sr > 0.0 "Total rate must be positive"
	isr = min(1.0e9,1.0 / sr)
	F(xdot, x, xd, tau, parms)
	xdot[end] = 1.0
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] * isr
	end
	nothing
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
- **ode**: ode time stepper, must be one of those: [:cvode,:lsoda,:Adams,:BDF]
- **save_at**: array of ordered time at which the solution is required
"""
function chv!(xc0::vecc,xd0::vecd,
				F::Function,R::Function,DX::Function,
				nu::Tnu,parms,
				ti::Float64, tf::Float64,
				verbose::Bool = false;
				ode=:cvode,ind_save_d=-1:1,ind_save_c=-1:1,dt=0.001,n_max::Int64 = Inf64) where {vecc <: AbstractVector{Float64}, vecd <: AbstractVector{Int64}, Tnu <:AbstractArray{Int64}}

	@assert ode in [:cvode,:lsoda,:adams,:bdf,:euler]

	n_max  += 1 # to hold initial vector
	nsteps  = 1 # index for the current jump number

	# Set up initial simulation time
	t = ti

	X0 = similar(xc0,length(xc0)+1)
	for ii in eachindex(xc0)
		X0[ii] = xc0[ii]
	end
	X0[end] = ti

	t_hist  = [ti]
	Xd     = copy(xd0)
	if ind_save_c[1] == -1
		ind_save_c = 1:length(xc0)
	end

	if ind_save_d[1] == -1
		ind_save_d = 1:length(xd0)
	end
	xc_hist = VectorOfArray([copy(xc0)[ind_save_c]])
	xd_hist = VectorOfArray([copy(xd0)[ind_save_d]])
	res_ode = zeros(2,length(X0))

	nsteps += 1

	deltaxd = copy(nu[1,:]) # declare this variable, variable to hold discrete jump
	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)  # vector of rates

	# define the ODE flow, this leads to big memory saving
	if ode==:cvode
		Flow = (X0_,Xd_,Δt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., Δt], abstol = 1e-9, reltol = 1e-7, integrator = :BDF)
	elseif	ode==:bdf
		Flow = (X0_,Xd_,Δt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., Δt], abstol = 1e-9, reltol = 1e-7, integrator = :BDF)
	elseif	ode==:adams
		Flow = (X0_,Xd_,Δt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., Δt], abstol = 1e-9, reltol = 1e-7, integrator = :Adams)
	elseif ode==:lsoda
		Flow = (X0_,Xd_,Δt,r_)->LSODA.lsoda((tt,x,xdot,data)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, [0., Δt], abstol = 1e-9, reltol = 1e-7)
	elseif ode==:euler
		Flow = (X0_,Xd_,Δt,r_)->euler( (tt,x,xdot)->f_CHV!(F,R,tt,x,xdot,Xd_,parms,r_), X0_, dt, 0., Δt)
	end
	δt = - log(rand())

	# Main loop
	while (t < tf) && (nsteps < n_max)

		verbose && println("--> t = ",t," - δt = ",δt, ",nstep =  ",nsteps)

		res_ode .= Flow(X0,Xd,δt,rate)

		verbose && println("--> ode solve is done!")

		# this holds the new state of the continuous component
		@inbounds for ii in eachindex(X0)
			X0[ii] = res_ode[end,ii]
		end

		# this is the next jump time
		t = res_ode[end,end]

		R(rate,X0,Xd,t,parms, false)

		# jump time:
		if (t < tf) && nsteps < n_max
			# Update event
			ev = pfsample(rate,sum(rate),numpf)
			deltaxd .= nu[ev,:]

			# Xd = Xd .+ deltaxd
			@inbounds for ii in eachindex(Xd)
				Xd[ii] += deltaxd[ii]
			end

			# Xc = Xc .+ deltaxc
			DX(X0,Xd,t,parms,ev)

			verbose && println("--> Which reaction? => ",ev)

			# save state, post-jump
			push!(t_hist,t)
			push!(xc_hist, X0[ind_save_c])
			push!(xd_hist, Xd[ind_save_d])

			δt = - log(rand())

		else
			if ode in [:cvode,:bdf,:adams]
				res_ode_last =   Sundials.cvode((tt,x,xdot)->F(xdot,x,Xd,tt,parms), xc_hist[end], [t_hist[end], tf], abstol = 1e-9, reltol = 1e-7)
			else#if ode==:lsoda
				res_ode_last = LSODA.lsoda((tt,x,xdot,data)->F(xdot,x,Xd,tt,parms), xc_hist[end], [t_hist[end], tf], abstol = 1e-9, reltol = 1e-7)
			end
			t = tf

			# save state
			push!(t_hist,tf)
			push!(xc_hist, res_ode_last[end,ind_save_c])
			push!(xd_hist, Xd[ind_save_d])
		end
		nsteps += 1
	end
	verbose && println("-->Done")
	verbose && println("--> xc = ",xd_hist[:,1:nsteps-1])
	return PDMPResult(t_hist,xc_hist,xd_hist,Float64[])
end
