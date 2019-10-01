###################################################################################################
include("jumps.jl")
include("utils.jl")

struct CHV{Tode} <: AbstractCHVIterator
	ode::Tode	# ODE solver to use for the flow in between jumps
end

function (chv::CHV{Tode})(xdot, x, prob::Tpb, t) where {Tode, Tpb <: PDMPCaracteristics}
	tau = x[end]
	rate = get_rate(prob.ratecache, x)
	sr = prob.R(rate, x, prob.xd, prob.parms, tau, true)[1]
	prob.F(xdot, x, prob.xd, prob.parms, tau)
	xdot[end] = 1.0
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] / sr
	end
	return nothing
end

###################################################################################################
### implementation of the CHV algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function chvjump(integrator, prob::PDMPProblem, save_pre_jump, save_rate, verbose)
	# we declare the characteristics for convenience
	caract = prob.caract
	ratecache = caract.ratecache
	simjptimes = prob.simjptimes

	# final simulation time
	tf = prob.tspan[2]

	# find the next jump time
	t = integrator.u[end]
	simjptimes.lastjumptime = t

	verbose && printstyled(color=:green, "--> Jump detected at t = $t !!\n")
	verbose && printstyled(color=:green, "--> jump not yet performed, xd = ", caract.xd,"\n")

	if (save_pre_jump) && (t <= tf)
		verbose && printstyled(color=:green, "----> saving pre-jump\n")
		pushXc!(prob, (integrator.u[1:end-1]))
		pushXd!(prob, copy(caract.xd))
		pushTime!(prob, t)
		#save rates for debugging
		save_rate && push!(prob.rate_hist, sum(ratecache.rate))
	end

	# execute the jump
	caract.R(get_rate(ratecache, integrator.u), integrator.u, caract.xd, caract.parms, t, false)
	if (t < tf)
		#save rates for debugging
		save_rate && push!(prob.rate_hist, sum(ratecache.rate))

		# Update event
		ev = pfsample(ratecache.rate, sum(ratecache.rate), length(ratecache.rate))

		# we perform the jump
		affect!(caract.pdmpjump, ev, integrator.u, caract.xd, caract.parms, t)

		u_modified!(integrator, true)

		@inbounds for ii in eachindex(caract.xc)
			caract.xc[ii] = integrator.u[ii]
		end
	end
	verbose && printstyled(color=:green,"--> jump computed, xd = ",caract.xd,"\n")
	# we register the next time interval to solve the extended ode
	simjptimes.njumps += 1
	simjptimes.tstop_extended += -log(rand())
	add_tstop!(integrator, simjptimes.tstop_extended)
	verbose && printstyled(color=:green,"--> End jump\n\n")
end

function chv_diffeq!(problem::PDMPProblem,
			ti::Tc, tf::Tc, X_extended::vece,
			verbose = false; ode = Tsit5(), save_positions = (false, true), n_jumps::Td = Inf64, reltol=1e-7, abstol=1e-9, save_rate = false) where {Tc, Td, vece}
	verbose && println("#"^30)
	verbose && printstyled(color=:red,"Entry in chv_diffeq\n")

	ti, tf = problem.tspan
	algopdmp = CHV(ode)

	# initialise the problem. If I call twice this solve function, it should give the same result...
	init!(problem)

	# we declare the characteristics for convenience
	caract = problem.caract
	ratecache = caract.ratecache
	simjptimes = problem.simjptimes

#ISSUE HERE, IF USING A PROBLEM p MAKE SURE THE TIMES in p.sim ARE WELL SET
	# set up the current time as the initial time
	t = ti
	# previous jump time, needed because problem.simjptimes.lastjumptime contains next jump time even if above tf
	tprev = t

	# vector to hold the state space for the extended system
	# X_extended = similar(problem.xc, length(problem.xc) + 1)
	# @show typeof(X_extended) vece

	for ii in eachindex(caract.xc)
		X_extended[ii] = caract.xc[ii]
	end
	X_extended[end] = ti

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, integrator -> chvjump(integrator, problem, save_positions[1], save_rate, verbose), save_positions = (false, false))

	# define the ODE flow, this leads to big memory saving
	# prob_CHV = ODEProblem((xdot,x,data,tt) -> problem(xdot, x, data, tt), X_extended, (0.0, 1e9))
	prob_CHV = ODEProblem((xdot, x, data, tt) -> algopdmp(xdot, x, caract, tt), X_extended, (0.0, 1e9))
	integrator = init(prob_CHV, ode, tstops = simjptimes.tstop_extended, callback = cb, save_everystep = false, reltol = reltol, abstol = abstol, advance_to_tstop = true)

	# current jump number
	njumps = 0

	while (t < tf) && simjptimes.njumps < n_jumps-1
		verbose && println("--> n = $(problem.simjptimes.njumps), t = $t, δt = ", simjptimes.tstop_extended)
		step!(integrator)

		@assert( t < simjptimes.lastjumptime, "Could not compute next jump time $(simjptimes.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(simjptimes.lastjumptime),\n solver = $ode. dt = $(t - simjptimes.lastjumptime)")
		t, tprev = simjptimes.lastjumptime, t

		# the previous step was a jump! should we save it?
		if njumps < simjptimes.njumps && save_positions[2] && (t <= tf)
			verbose && println("----> save post-jump, xd = ",problem.Xd)
			pushXc!(problem, copy(caract.xc))
			pushXd!(problem, copy(caract.xd))
			pushTime!(problem, t)
			njumps +=1
			verbose && println("----> end save post-jump, ")
		end
	end
	# we check that the last bit [t_last_jump, tf] is not missing
	if t>tf
		verbose && println("----> LAST BIT!!, xc = ", caract.xc[end], ", xd = ", caract.xd, ", t = ", problem.time[end])
		prob_last_bit = ODEProblem((xdot,x,data,tt) -> caract.F(xdot, x, caract.xd, caract.parms, tt), copy(caract.xc), (tprev, tf))
		sol = DiffEqBase.solve(prob_last_bit, ode)
		verbose && println("-------> xc[end] = ",sol.u[end])
		pushXc!(problem, sol.u[end])
		pushXd!(problem, copy(caract.xd))
		pushTime!(problem, sol.t[end])
	end
	return PDMPResult(problem, save_positions)
end

function solve(problem::PDMPProblem{Tc, Td, vectype_xc, vectype_xd, Tcar}, algo::CHV{Tode}, X_extended; verbose = false, n_jumps = Inf64, save_positions = (false, true), reltol = 1e-7, abstol = 1e-9, save_rate = false) where {Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar, Tode <: DiffEqBase.DEAlgorithm}

	return chv_diffeq!(problem, problem.tspan[1], problem.tspan[2], X_extended, verbose; ode = algo.ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol, save_rate = save_rate)
end

function solve(problem::PDMPProblem{Tc, Td, vectype_xc, vectype_xd, Tcar}, algo::CHV{Tode}; verbose = false, n_jumps = Inf64, save_positions = (false, true), reltol = 1e-7, abstol = 1e-9, save_rate = false) where {Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar, Tode <: DiffEqBase.DEAlgorithm}

	# resize the extended vector to the proper dimension
	X_extended = zeros(Tc, length(problem.caract.xc) + 1)

	return chv_diffeq!(problem, problem.tspan[1], problem.tspan[2], X_extended, verbose; ode = algo.ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol, save_rate = save_rate)
end
