###################################################################################################
include("jumps.jl")
include("utils.jl")

struct CHV{Tode} <: AbstractCHVIterator
	ode::Tode	# ODE solver to use for the flow in between jumps
end

function (chv::CHV{Tode})(xdot, x, prob::Tpb, t) where {Tode, Tpb <: PDMPCaracteristics}
	tau = x[end]
	sr = prob.R(prob.ratecache, x, prob.xd, tau, prob.parms, true)[1]
	prob.F(xdot, x, prob.xd, tau, prob.parms)
	xdot[end] = 1.0
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] / sr
	end
	return nothing
end


###################################################################################################
###################################################################################################
### implementation of the CHV algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function chvjump(integrator, prob::PDMPProblem, save_pre_jump)
	# final simulation time
	tf = prob.interval[2]
	
	# find the next jump time
	t = integrator.u[end]
	prob.simjptimes.lastjumptime = t

	# we declare the characteristics for convenience
	caract = prob.caract

	prob.verbose && printstyled(color=:green, "--> Jump detected at t = $t !!\n")
	prob.verbose && printstyled(color=:green, "--> jump not yet performed, xd = ", caract.xd,"\n")

	if (save_pre_jump) && (t <= prob.tf)
		prob.verbose && printstyled(color=:green, "----> saving pre-jump\n")
		push!(prob.Xc, (integrator.u[1:end-1]))
		push!(prob.Xd, copy(caract.xd))
		push!(prob.time,t)
	end

	# execute the jump
	caract.R(caract.ratecache, integrator.u, caract.xd, t, caract.parms, false)
	if (t < tf)
		#save rates for debugging
		prob.save_rate && push!(prob.rate_hist, sum(prob.rateCache.rate))

		# Update event
		ev = pfsample(caract.ratecache, sum(caract.ratecache), length(caract.ratecache))

		deltaxd = view(caract.pmdpjump.nu, ev, :)

		# Xd = Xd .+ deltaxd
		# LinearAlgebra.BLAS.axpy!(1.0, deltaxd, prob.xd)
		@inbounds for ii in eachindex(caract.xd)
			caract.xd[ii] += deltaxd[ii]
		end

		# Xc = Xc .+ deltaxc
		caract.pmdpjump.Delta(integrator.u, caract.xd, t, caract.parms, ev)
		u_modified!(integrator, true)

		@inbounds for ii in eachindex(caract.xc)
			caract.xc[ii] = integrator.u[ii]
		end
	end
	prob.verbose && printstyled(color=:green,"--> jump computed, xd = ",prob.xd,"\n")
	# we register the next time interval to solve the extended ode
	prob.simjptimes.njumps += 1
	prob.simjptimes.tstop_extended += -log(rand())
	add_tstop!(integrator, prob.simjptimes.tstop_extended)
	prob.verbose && printstyled(color=:green,"--> End jump\n\n")
end

"""
Implementation of the CHV method to sample a PDMP using the package `DifferentialEquations`. The advantage of doing so is to lower the number of calls to `solve` using an `integrator` method. The reason why we can pass `rate` vector and `xc0_extended` is to allow the use of `StaticArrays.jl` for which the constructs differ from Base.Array
"""
function chv_diffeq!(xc0::vecc, xd0::vecd,
		F::TF, R::TR, DX::TD,
		nu::Tnu, parms::Tp,
		ti::Tc, tf::Tc,
		verbose::Bool = false;
		ode = Tsit5(), n_jumps::Td = Inf64,
		save_positions		= (false, true),
		saverate			= false,
		rate::vecrate		= zeros(Tc, size(nu, 1)),
		xc0_extended::vece	= zeros(Tc, length(xc0) + 1),
		reltol=1e-7, abstol=1e-9,
		alg = Tsit5() ) where {Tc,Td,Tnu <: AbstractArray{Td}, Tp, TF ,TR ,TD,
		vecc <: AbstractVector{Tc},
		vecd <: AbstractVector{Td},
		vecrate <: AbstractVector{Tc},
		vece <: AbstractVector{Tc}}

	# custom type to collect all parameters in one structure
	problem  = PDMPProblem(xc0,xd0,rate,F,R,DX,nu,parms,(ti,tf),save_positions[1],verbose,alg,saverate)

	return	chv_diffeq!(problem, ti, tf, copy(xc0_extended), verbose; ode = ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol)
end


function chv_diffeq!(problem::PDMPProblem,
			ti::Tc, tf::Tc, X_extended::vece,
			verbose = false; ode = Tsit5(), save_positions = (false, true), n_jumps::Td = Inf64, reltol=1e-7, abstol=1e-9) where {Tc, Td, vece}
	problem.verbose && printstyled(color=:red,"Entry in chv_diffeq\n")

	algopdmp = CHV(ode)

	# we declare the characteristics for convenience
	caract = problem.caract

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
	cb = DiscreteCallback(problem, integrator -> chvjump(integrator, problem, save_positions[1]), save_positions = (false, false))

	# define the ODE flow, this leads to big memory saving
	# prob_CHV = ODEProblem((xdot,x,data,tt) -> problem(xdot, x, data, tt), X_extended, (0.0, 1e9))
	prob_CHV = ODEProblem((xdot,x,data,tt) -> algopdmp(xdot, x, caract, tt), X_extended, (0.0, 1e9))
	integrator = init(prob_CHV, ode, tstops = problem.simjptimes.tstop_extended, callback = cb, save_everystep = false, reltol = reltol, abstol = abstol, advance_to_tstop = true)

	# current jump number
	njumps = 0

	while (t < tf) && problem.simjptimes.njumps < n_jumps-1
		problem.verbose && println("--> n = $(problem.simjptimes.njumps), t = $t, δt = ",problem.simjptimes.tstop_extended)
		step!(integrator)
		@assert( t < problem.simjptimes.lastjumptime, "Could not compute next jump time $(problem.simjptimes.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(problem.simjptimes.lastjumptime),\n solver = $ode. dt = $(t - problem.simjptimes.lastjumptime)")
		t, tprev = problem.simjptimes.lastjumptime, t

		# the previous step was a jump! should we save it?
		if njumps < problem.simjptimes.njumps && save_positions[2] && (t <= tf)
			problem.verbose && println("----> save post-jump, xd = ",problem.Xd)
			push!(problem.Xc, copy(caract.xc))
			push!(problem.Xd, copy(caract.xd))
			push!(problem.time, t)
			njumps +=1
			problem.verbose && println("----> end save post-jump, ")
		end
	end
	# we check that the last bit [t_last_jump, tf] is not missing
	if t>tf
		problem.verbose && println("----> LAST BIT!!, xc = ",caract.xc[end], ", xd = ",caract.xd, ", t = ", problem.time[end])
		prob_last_bit = ODEProblem((xdot,x,data,tt) -> caract.F(xdot, x, caract.xd, tt, caract.parms), copy(caract.xc),(tprev,tf))
		sol = DiffEqBase.solve(prob_last_bit, ode)
		problem.verbose && println("-------> xc[end] = ",sol.u[end])
		push!(problem.Xc, sol.u[end])
		push!(problem.Xd, copy(caract.xd))
		push!(problem.time, sol.t[end])
	end
	return PDMPResult(problem.time, problem.Xc, problem.Xd, problem.rate_hist)
end


function solve(problem::PDMPProblem{Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar}, algo::CHV{Tode}; verbose = false, n_jumps = Inf64, X_extended = zeros(Tc, 1 + 1), save_positions = (false, true), reltol = 1e-7, abstol = 1e-9) where {Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar, Tode <: DiffEqBase.DEAlgorithm}
	# hack to resize the extended vector to the proper dimension
	resize!(X_extended, length(problem.caract.xc) + 1)

	return chv_diffeq!(problem, problem.simjptimes.lastjumptime, problem.tf, X_extended, verbose; ode = algo.ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol )
end
