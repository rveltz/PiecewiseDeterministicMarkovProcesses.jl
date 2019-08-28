include("chvdiffeq.jl")

### WARNING This is an old ODE solver which is not based on an iterator implementation. We keep it until LSODA has an iterator implementation

function f_CHV!(F, R, t::Float64, x, xdot, xd, parms, rate)
	# used for the exact method
	# we put [1] to use it in the case of the rejection method as well
	tau = x[end]
	sr = R(rate, x, xd, tau, parms, true)[1]
	@assert sr > 0.0 "Total rate must be positive"
	isr = min(1.0e9, 1.0 / sr)
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
function solve(problem::PDMPProblem, algo::CHV{Tode}; verbose::Bool = false, ind_save_d=-1:1, ind_save_c=-1:1, dt=0.001, n_jumps = Inf64, reltol = 1e-7, abstol = 1e-9, save_positions = (false,true),) where {Tode <: Symbol}
	verbose && println("#"^30)
	ode = algo.ode
	@assert ode in [:cvode, :lsoda, :adams, :bdf, :euler]

	ti, tf = problem.interval
	n_jumps  += 1 # to hold initial vector
	nsteps  = 1 # index for the current jump number

	xc0 = problem.caract.xc
	xd0 = problem.caract.xd

	# Set up initial simulation time
	t = ti

	X_extended = similar(xc0, length(xc0) + 1)
	for ii in eachindex(xc0)
		X_extended[ii] = xc0[ii]
	end
	X_extended[end] = ti

	t_hist  = [ti]
	Xd = copy(xd0)
	if ind_save_c[1] == -1
		ind_save_c = 1:length(xc0)
	end

	if ind_save_d[1] == -1
		ind_save_d = 1:length(xd0)
	end
	xc_hist = VectorOfArray([copy(xc0)[ind_save_c]])
	xd_hist = VectorOfArray([copy(xd0)[ind_save_d]])
	res_ode = zeros(2, length(X_extended))

	nsteps += 1

	deltaxd = copy(problem.caract.pdmpjump.nu[1,:]) # declare this variable, variable to hold discrete jump
	numpf   = size(problem.caract.pdmpjump.nu,1)    # number of reactions
	rate    = zeros(numpf)  # vector of rates

	# define the ODE flow, this leads to big memory saving
	if ode == :cvode || ode == :bdf
		Flow = (X0_,Xd_,Δt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(problem.caract.F,problem.caract.R,tt,x,xdot,Xd_,problem.caract.parms,r_), X0_, [0., Δt], abstol = abstol, reltol = reltol, integrator = :BDF)
	elseif	ode==:adams
		Flow = (X0_,Xd_,Δt,r_)->Sundials.cvode(  (tt,x,xdot)->f_CHV!(problem.caract.F,problem.caract.R,tt,x,xdot,Xd_,problem.caract.parms,r_), X0_, [0., Δt], abstol = abstol, reltol = reltol, integrator = :Adams)
	elseif ode==:lsoda
		Flow = (X0_,Xd_,Δt,r_)->LSODA.lsoda((tt,x,xdot,data)->f_CHV!(problem.caract.F,problem.caract.R,tt,x,xdot,Xd_,problem.caract.parms,r_), X0_, [0., Δt], abstol = abstol, reltol = reltol)
	elseif ode==:euler
		Flow = (X0_,Xd_,Δt,r_)->euler( (tt,x,xdot)->f_CHV!(problem.caract.F,problem.caract.R,tt,x,xdot,Xd_,problem.caract.parms,r_), X0_, dt, 0., Δt)
	end

	# we use the first time interval from the one generated by the constructor PDMPProblem
	δt = problem.simjptimes.tstop_extended

	# Main loop
	while (t < tf) && (nsteps < n_jumps)

		verbose && println("--> t = ", t," - δt = ", δt, ",nstep =  ", nsteps)

		res_ode .= Flow(X_extended, Xd, δt, rate)

		verbose && println("--> ode solve is done!")

		# this holds the new state of the continuous component
		@inbounds for ii in eachindex(X_extended)
			X_extended[ii] = res_ode[end, ii]
		end

		# this is the next jump time
		t = res_ode[end, end]

		problem.caract.R(rate, X_extended, Xd, t, problem.caract.parms, false)

		# jump time:
		if (t < tf) && nsteps < n_jumps
			# Update event
			ev = pfsample(rate, sum(rate), numpf)
			deltaxd .= problem.caract.pdmpjump.nu[ev,:]

			# Xd = Xd .+ deltaxd
			@inbounds for ii in eachindex(Xd)
				Xd[ii] += deltaxd[ii]
			end

			# Xc = Xc .+ deltaxc
			problem.caract.pdmpjump.Delta(X_extended, Xd, t, problem.caract.parms, ev)

			verbose && println("--> Which reaction? => ", ev)
			verbose && println("--> xd = ", Xd)

			# save state, post-jump
			push!(t_hist, t)
			push!(xc_hist, X_extended[ind_save_c])
			push!(xd_hist, Xd[ind_save_d])

			δt = - log(rand())

		else
			if ode in [:cvode, :bdf, :adams]
				res_ode_last = Sundials.cvode((tt, x, xdot)->problem.caract.F(xdot,x,Xd,tt,problem.caract.parms), xc_hist[end], [t_hist[end], tf], abstol = 1e-9, reltol = 1e-7)
			else#if ode==:lsoda
				res_ode_last = LSODA.lsoda((tt, x, xdot, data)->problem.caract.F(xdot,x,Xd,tt,problem.caract.parms), xc_hist[end], [t_hist[end], tf], abstol = 1e-9, reltol = 1e-7)
			end
			t = tf

			# save state
			push!(t_hist, tf)
			push!(xc_hist, res_ode_last[end,ind_save_c])
			push!(xd_hist, Xd[ind_save_d])
		end
		nsteps += 1
	end
	verbose && println("-->Done")
	verbose && println("--> xc = ", xd_hist[:,1:nsteps-1])
	return PDMPResult(t_hist, xc_hist, xd_hist,Float64[])

end
