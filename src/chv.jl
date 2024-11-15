### WARNING This is an old ODE solver which is not based on an iterator implementation. We keep it until LSODA has an iterator implementation

"""
Same as the `solve` for `CHV(::DiffEqBase.DEAlgorithm)` but for `CHV(::Symbol)`. This is an old implementation of the CHV algorithm which can be used with `:lsoda`. For all other solvers, use the the new solver.
"""
function solve(problem::PDMPProblem, algo::CHV{Tode};
				verbose::Bool = false,
				ind_save_d = -1:1,
				ind_save_c = -1:1,
				n_jumps = Inf64,
				reltol = 1e-7,
				abstol = 1e-9,
				save_positions = (false,
				true),
				save_rate = false,
				finalizer = finalize_dummy,
				kwargs...) where {Tode <: Symbol}
	verbose && println("#"^30)
	ode = algo.ode
	@assert ode in [:cvode, :lsoda, :adams, :BDF]
	verbose && printstyled(color=:red,"--> Start CHV method (algo::Symbol)\n")

	# table to use DiffEqBase
	odeTable = Dict(:lsoda => lsoda(),
					:BDF => CVODE_BDF(),
					:adams => CVODE_Adams(),
					:cvode => CVODE_BDF())

	# initialise the problem. If I call twice this solve function, it should give the same result...
	init!(problem)

	# we declare the characteristics for convenience
	caract = problem.caract
	ratecache = caract.ratecache

	ti, tf = problem.tspan
	n_jumps += 1 # to hold initial vector
	nsteps   = 1 # index for the current jump number

	xc0 = caract.xc0
	xd0 = caract.xd0

	# Set up initial simulation time
	t = ti

	X_extended = similar(xc0, length(xc0) + 1)
	for ii in eachindex(xc0)
		X_extended[ii] = xc0[ii]
	end
	X_extended[end] = ti

	#useful to use the same array, as it can be used in CHV(ode)
	Xd = caract.xd
	if ind_save_c[1] == -1
		ind_save_c = 1:length(xc0)
	end

	if ind_save_d[1] == -1
		ind_save_d = 1:length(xd0)
	end
	xc_hist = VectorOfArray([copy(xc0)[ind_save_c]])
	xd_hist = VectorOfArray([copy(xd0)[ind_save_d]])

	res_ode = zeros(length(X_extended))

	nsteps += 1

	probExtLsoda = ODEProblem((du, u, p, _t) -> algo(du, u, caract, _t), copy(X_extended), (ti, tf))

	function Flow(_X0, _Xd, Δt, _r; _alg = odeTable[ode])
		prob = DiffEqBase.remake(probExtLsoda; tspan = (0, Δt))
		prob.u0 .= _X0
		sol = solve(prob, _alg; abstol = abstol, reltol = reltol, save_everystep = false)
		return sol.u[end]
	end

	# we use the first time interval from the one generated by the constructor PDMPProblem
	δt = problem.simjptimes.tstop_extended

	# Main loop
	while (t < tf) && (nsteps < n_jumps)

		verbose && println("├─── t = $t, -log(U) = $δt, nstep =  $nsteps")

		res_ode .= Flow(X_extended, Xd, δt, get_tmp(ratecache, X_extended))

		verbose && println("│    ode solve has been performed!")

		if (res_ode[end] < tf) && nsteps < n_jumps
			verbose && println("│    Δt = ", res_ode[end] - t)
			# this is the next jump time
			t = res_ode[end]

			# this holds the new state of the continuous component
			@inbounds for ii in eachindex(X_extended)
				X_extended[ii] = res_ode[ii]
			end

			caract.R(get_tmp(ratecache, X_extended), X_extended, Xd, caract.parms, t, false)

			# Update event
			ev = pfsample(get_tmp(ratecache, X_extended))

			# we perform the jump, it changes Xd and (possibly) X_extended
			affect!(caract.pdmpjump, ev, X_extended, Xd, caract.parms, t)

			verbose && println("│    reaction = ", ev)
			# verbose && println("--> xd = ", Xd)

			# save state, post-jump
			if save_positions[2] || (nsteps == n_jumps - 1)
				pushTime!(problem, t)
				push!(xc_hist, copy(X_extended[ind_save_c]))
				push!(xd_hist, copy(Xd[ind_save_d]))
			end

			save_rate && push!(problem.rate_hist, caract.R(get_tmp(ratecache, X_extended), X_extended, Xd, caract.parms, t, true)[1])

			finalizer(get_tmp(ratecache, X_extended), caract.xc, caract.xd, caract.parms, t)

			δt = - log(rand())

		else
			probLast = ODEProblem((du, u, p, _t) -> caract.F(du, u, Xd, caract.parms, _t), X_extended[1:end-1], (t, tf))
			res_ode_last = solve(probLast, odeTable[ode]; abstol = 1e-9, reltol = 1e-7, save_everystep = false)

			t = tf

			# save state
			pushTime!(problem, tf)
			push!(xc_hist, copy(res_ode_last[end][ind_save_c]))
			push!(xd_hist, copy(Xd[ind_save_d]))
		end
		nsteps += 1
	end
	verbose && println("--> Done")
	if verbose && save_positions[2]
		println("--> xc = ", xd_hist[:, 1:nsteps-1])
	end
	return PDMPResult(copy(problem.time), xc_hist, xd_hist, problem.rate_hist, save_positions, length(problem.time), 0)
end
