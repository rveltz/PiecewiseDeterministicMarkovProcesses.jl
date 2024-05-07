struct Rejection{Tode} <: AbstractCHVIterator
	ode::Tode	# ODE solver to use for the flow in between jumps
end

function (rej::Rejection{Tode})(xdot, x, prob::Tpb, t) where {Tode, Tpb <: PDMPCaracteristics}
	prob.F(xdot, x, prob.xd, prob.parms, t)
	return nothing
end
###################################################################################################
### implementation of the rejection algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function rejectionjump(integrator, prob::PDMPProblem, save_pre_jump, save_rate, verbose)
	# we declare the characteristics for convenience
	caract = prob.caract
	ratecache = caract.ratecache
	simjptimes = prob.simjptimes

	# final simulation time
	tf = prob.tspan[2]

	# find the next jump time
	t = integrator.t
	simjptimes.lastjumptime = t

	verbose && printstyled(color=:red,"--> REJECTION JUMP, t = $t\n")
	verbose && printstyled(color=:red,"----> xc = $(integrator.u)\n")

	verbose && printstyled(color=:green,"--> Fictitous jump at t = $t, # = ",simjptimes.fictitous_jumps," !!\n")

	simjptimes.ppf .= caract.R(get_tmp(ratecache, integrator.u), integrator.u, caract.xd, caract.parms, t, true)
	@assert simjptimes.ppf[1] < simjptimes.ppf[2] "Error, your bound on the rates is not high enough!, $(simjptimes.ppf)"

	simjptimes.reject = rand() < 1 - simjptimes.ppf[1] / simjptimes.ppf[2]
	δt = -log(rand()) / simjptimes.ppf[2]

	verbose && printstyled(color=:green,"----> xc = ",caract.xc ,", xd = ",caract.xd,", reject = ",simjptimes.reject,", rates = ",simjptimes.ppf,"\n")

	# execute the jump
	if t < tf && simjptimes.reject == false
		verbose && println("----> Jump!, ratio = ", simjptimes.ppf[1] / simjptimes.ppf[2], ", xd = ", caract.xd)
		simjptimes.ppf .= caract.R(get_tmp(ratecache, integrator.u), integrator.u, caract.xd, caract.parms, t, false)

		if (save_pre_jump) && (t <= tf)
			verbose && printstyled(color=:green,"----> save pre-jump\n")
			pushXc!(prob, integrator.u)
			pushXd!(prob, copy(caract.xd))
			pushTime!(prob, t)
			#save rates for debugging
			save_rate && push!(prob.rate_hist, sum(ratecache.rate))
		end

		# Update event
		ev = pfsample(get_tmp(ratecache, integrator.u))

		# we perform the jump
		affect!(caract.pdmpjump, ev, integrator.u, caract.xd, caract.parms, t)
		u_modified!(integrator, true)

		@inbounds for ii in eachindex(caract.xc)
			caract.xc[ii] = integrator.u[ii]
		end
		simjptimes.njumps += 1
	else
		simjptimes.fictitous_jumps += 1
	end
	verbose && printstyled(color=:green,"----> jump effectued, xd = ",caract.xd,"\n")

	# we register the next time interval to solve the extended ode
	simjptimes.tstop_extended += δt
	add_tstop!(integrator, simjptimes.tstop_extended)
	verbose && printstyled(color=:green,"--> End jump\n\n")
end

function rejection_diffeq!(problem::PDMPProblem,
				ti::Tc, tf::Tc, verbose = false; ode = Tsit5(),
				save_positions = (false, true), n_jumps::Td = Inf64, reltol=1e-7, abstol=1e-9, save_rate = false, finalizer = finalize_dummy) where {Tc, Td}
	verbose && println("#"^30)
	verbose && printstyled(color=:red,"Entry in rejection_diffeq\n")
	ti, tf = problem.tspan
	algopdmp = Rejection(ode)

#ISSUE HERE, IF USING A PROBLEM p MAKE SURE THE TIMES in p.sim ARE WELL SET
	# set up the current time as the initial time
	t = ti
	# previous jump time, needed because problem.simjptimes.lastjumptime contains next jump time even if above tf
	tprev = t

	# initialise the problem. If I call twice this function, it should give the same result...
	init!(problem)

	# we declare the characteristics for convenience
	caract = problem.caract
	ratecache = caract.ratecache
	simjptimes = problem.simjptimes

	# vector to hold the state space
	X0 = copy(caract.xc)

	# current jump number
	# njumps = 0
	simjptimes.njumps = 1
	simjptimes.lambda_star = 0.0	# this is the bound for the rejection method
	simjptimes.ppf .= caract.R(get_tmp(ratecache, X0), X0, caract.xd, caract.parms, t, true)

	simjptimes.tstop_extended = simjptimes.tstop_extended / simjptimes.ppf[2] + ti

	simjptimes.reject = true

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, integrator -> rejectionjump(integrator, problem, save_positions[1], save_rate, verbose), save_positions = (false, false))

	# define the ODE flow, this leads to big memory saving
	prob_REJ = ODEProblem((xdot, x, data, tt) -> algopdmp(xdot, x, caract, tt), X0, (ti, 1e9))
	integrator = init(prob_REJ, ode, tstops = simjptimes.tstop_extended, callback = cb, save_everystep = false, reltol = reltol, abstol = abstol, advance_to_tstop = true)

	while (t < tf) && simjptimes.njumps < n_jumps #&& simjptimes.fictitous_jumps < 10
		verbose && println("--> n = $(simjptimes.njumps), t = $t -> ", simjptimes.tstop_extended)
		step!(integrator)
		if t >= simjptimes.lastjumptime
			@warn "Could not compute next jump time $(simjptimes.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(simjptimes.lastjumptime),\n solver = $ode"
			return PDMPResult(problem, save_positions)
		end
		t, tprev = simjptimes.lastjumptime, t

		# the previous step was a jump!
		if save_positions[2] && (t <= tf) && simjptimes.reject == false
			verbose && println("----> save post-jump, xd = ", problem.Xd)
			pushTime!(problem ,t)
			pushXc!(problem, copy(caract.xc))
			pushXd!(problem, copy(caract.xd))
			#save rates for debugging
			save_rate && push!(problem.rate_hist, sum(get_tmp(ratecache, X0)))
			verbose && println("----> end save post-jump, ")
			#put the flag for rejection
			simjptimes.reject = true
		end
		finalizer(get_tmp(ratecache, X0), caract.xc, caract.xd, caract.parms, t)
	end
	# we check whether the last bit [t_last_jump, tf] is missing
	if t>tf
		verbose && println("----> LAST BIT!!, xc = ", caract.xc[end], ", xd = ",caract.xd)
		prob_last_bit = ODEProblem((xdot,x,data,tt) -> caract.F(xdot, x, caract.xd, caract.parms, tt), copy(caract.xc), (tprev, tf))
		sol = DiffEqBase.solve(prob_last_bit, ode)
		verbose && println("-------> xc[end] = ", sol.u[end])
		pushXc!(problem, sol.u[end])
		pushXd!(problem, copy(caract.xd))
		pushTime!(problem, sol.t[end])
	end
	return PDMPResult(problem, save_positions)
end


function solve(problem::PDMPProblem{Tc, Td, vectype_xc, vectype_xd, Tcar}, algo::Rejection{Tode}; verbose = false, n_jumps = Inf64, save_positions = (false, true), reltol = 1e-7, abstol = 1e-9, save_rate = true, finalizer = finalize_dummy) where {Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar, Tode <: DiffEqBase.DEAlgorithm}

	return rejection_diffeq!(problem, problem.tspan[1], problem.tspan[2], verbose; ode = algo.ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol, save_rate = save_rate, finalizer = finalizer )
end
