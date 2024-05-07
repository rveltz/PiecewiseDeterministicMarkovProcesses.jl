struct RejectionExact <: AbstractRejectionExact end

function solve(problem::PDMPProblem, Flow::Function; verbose::Bool = false, save_rejected = false, ind_save_d = -1:1, ind_save_c = -1:1, n_jumps = Inf64, save_positions = (false, true), save_rate = false, finalizer = finalize_dummy)
	verbose && println("#"^30)
	verbose && printstyled(color=:red,"--> Start Rejection method\n")

	# initialise the problem. If I call twice this function, it should give the same result...
	init!(problem)

	# we declare the characteristics for convenience
	caract = problem.caract
	ratecache = caract.ratecache
	simjptimes = problem.simjptimes

	ti, tf = problem.tspan
	# it is faster to pre-allocate arrays and fill it at run time
	n_jumps  += 0 # to hold initial vector
	n_reject  = 0  # to hold the number of rejects
	nsteps  = 1
	npoints = 2 # number of points for ODE integration

	xc0 = caract.xc
	xd0 = caract.xd

	# Set up initial variables
	t = ti

	X0 = copy(xc0)
	Xd = copy(xd0)
	res_ode = zeros(2, length(X0))

	X0, _, Xd, _, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c = allocate_arrays(ti, xc0, xd0, n_jumps; rejection = true)

	tp = [ti, tf]		   # vector to hold the time interval over which to integrate the flow

	#variables for rejection algorithm
	reject = true
	lambda_star = 0.0 # this is the bound for the rejection method
	ppf = caract.R(get_tmp(ratecache, X0), X0, Xd, caract.parms, t, true)

	δt = simjptimes.tstop_extended

	while (t < tf) && (nsteps < n_jumps)
		verbose && println("--> step : ",nsteps," / ", n_jumps)
		reject = true
		while reject && (nsteps < n_jumps)
			tp .= [t, min(tf, t + δt / ppf[2]) ] #mettre un lambda_star?
			Flow(res_ode, X0, Xd, tp)

			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end, ii]
			end
			verbose && println("----> δt = ", δt, ", t∈", tp, ", dt = ", tp[2]-tp[1], ", xc = ", X0)

			t = tp[end]
			ppf = caract.R(get_tmp(ratecache, X0), X0, Xd, caract.parms, t, true)
			@assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf"
			if t == tf
				reject = false
			else
				reject = rand() < 1 - ppf[1] / ppf[2]
			end
			δt = -log(rand())
			if reject
				n_reject += 1
			end
		end

		# there is a jump!
		ppf = caract.R(get_tmp(ratecache, X0), X0, Xd, caract.parms, t, false)

		if (t < tf)
			verbose && println("----> Jump!, ratio = ", ppf[1] / ppf[2], ", xd = ", Xd)
			# make a jump
			ev = pfsample(get_tmp(ratecache, X0))

			# we perform the jump
			affect!(caract.pdmpjump, ev, X0, Xd, caract.parms, t)
		end

		nsteps += 1
		pushTime!(problem, t)
		push!(xc_hist, X0[ind_save_c])
		push!(xd_hist, Xd[ind_save_d])
		save_rate && push!(problem.rate_hist, sum(get_tmp(ratecache, X0)))
		finalizer(get_tmp(ratecache, X0), caract.xc, caract.xd, caract.parms, t)
	end
	if verbose println("--> Done") end
	if verbose println("--> xd = ",xd_hist[:,1:nsteps]) end
	return PDMPResult(problem.time, xc_hist, xd_hist, problem.rate_hist, save_positions, nsteps, n_reject)
end

function solve(problem::PDMPProblem, algo::Rejection{Tode}; reltol = 1e-7, abstol = 1e-9, kwargs...) where {Tode <: Symbol}
	ode = algo.ode
	@assert ode in [:cvode, :lsoda, :adams, :bdf]

	caract = problem.caract

	# define the ODE flow
	if ode == :cvode || ode == :bdf
		Flow0 = (X0_,Xd,tp_) -> Sundials.cvode(  (tt,x,xdot) -> caract.F(xdot,x,Xd,caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :BDF)
	elseif	ode == :adams
		Flow0 = (X0_,Xd,tp_) -> Sundials.cvode(  (tt,x,xdot) -> caract.F(xdot,x,Xd,caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :Adams)
	elseif ode == :lsoda
		Flow0 = (X0_,Xd,tp_) -> LSODA.lsoda((tt,x,xdot,data) -> caract.F(xdot,x,Xd,caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol)
	end

	Flow = (out,X0_,Xd,tp_) -> (out .= Flow0(X0_,Xd,tp_))

	return solve(problem, Flow; kwargs...)
end

function solve(problem::PDMPProblem, algo::Talgo; kwargs...) where {Talgo <: AbstractRejectionExact}
	Flow = (res_ode, X0, Xd, tp) -> problem.caract.F(res_ode, X0, Xd, problem.caract.parms, tp)
	solve(problem, Flow; kwargs...)
end
