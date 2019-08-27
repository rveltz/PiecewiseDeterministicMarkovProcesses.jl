###################################################################################################
include("jumps.jl")
include("utils.jl")

struct Rejection{Tode} <: AbstractCHVIterator
	ode::Tode	# ODE solver to use for the flow in between jumps
end

function (rej::Rejection{Tode})(xdot, x, prob::Tpb, t) where {Tode, Tpb <: PDMPCaracteristics}
	prob.F(xdot, x, prob.xd, t, prob.parms)
	return nothing
end
###################################################################################################
### implementation of the rejection algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function rejectionjump(integrator, prob::PDMPProblem, save_pre_jump, verbose)
	verbose && printstyled(color=:red,"--> REJECTION JUMP\n")
	# final simulation time
	tf = prob.interval[2]

	# find the next jump time
	t = integrator.t
	prob.simjptimes.lastjumptime = t

	# we declare the characteristics for convenience
	caract = prob.caract

	verbose && printstyled(color=:green,"--> Fictitous jump at t = $t, # = ",prob.simjptimes.fictitous_jumps," !!\n")

	prob.simjptimes.ppf .= caract.R(caract.ratecache, integrator.u, caract.xd, t, caract.parms, true)
	@assert prob.simjptimes.ppf[1] < prob.simjptimes.ppf[2] "Error, your bound on the rates is not high enough!, $(prob.simjptimes.ppf)"
	prob.simjptimes.reject = rand() < (1 - prob.simjptimes.ppf[1] / prob.simjptimes.ppf[2])

	verbose && printstyled(color=:green,"----> xc = ",caract.xc ,", xd = ",caract.xd,", reject = ",prob.simjptimes.reject,", rates = ",prob.simjptimes.ppf,"\n")

	# execute the jump
	if t < tf && prob.simjptimes.reject == false
		verbose && printstyled(color=:blue,"--> TRUE jump!\n")
		prob.simjptimes.ppf .= caract.R(caract.ratecache, caract.xc, caract.xd, t, caract.parms, false)

		if (save_pre_jump) && (t <= tf)
			verbose && printstyled(color=:green,"----> save pre-jump\n")
			push!(prob.Xc, (integrator.u[1:end-1]))
			push!(prob.Xd, copy(caract.xd))
			push!(prob.time,t)
		end

		#save rates for debugging
		prob.save_rate && push!(prob.rate_hist, sum(caract.ratecache))

		# Update event
		ev = pfsample(caract.ratecache, sum(caract.ratecache), length(caract.ratecache))

		deltaxd = view(caract.pdmpjump.nu,ev,:)

		# Xd = Xd .+ deltaxd
		# LinearAlgebra.BLAS.axpy!(1.0, deltaxd, prob.xd)
		@inbounds for ii in eachindex(caract.xd)
			caract.xd[ii] += deltaxd[ii]
		end

		# Xc = Xc .+ deltaxc
		caract.pdmpjump.Delta(integrator.u, caract.xd, t, caract.parms, ev)
		u_modified!(integrator, true)

		@inbounds for ii in eachindex(caract.xc)
			caract.xc[ii] = integrator.u[ii]
		end

		prob.simjptimes.njumps += 1
	else
		prob.simjptimes.fictitous_jumps += 1
	end
	verbose && printstyled(color=:green,"----> jump effectued, xd = ",caract.xd,"\n")
	# we register the next time interval to solve the extended ode
	prob.simjptimes.tstop_extended += -log(rand()) / prob.simjptimes.ppf[2]
	add_tstop!(integrator, prob.simjptimes.tstop_extended)
	verbose && printstyled(color=:green,"--> End jump\n\n")
end

"""
Implementation of the rejection method to sample a PDMP using the package `DifferentialEquations`. The advantage of doing so is to lower the number of calls to `solve` using an `integrator` method.
"""
#TODO The two following functions should be merged with solve
function rejection_diffeq!(xc0::vecc, xd0::vecd,
		F::TF, R::TR, DX::TD,
		nu::Tnu, parms::Tp,
		ti::Tc, tf::Tc,
		verbose::Bool = false;
		ode = Tsit5(), save_positions=(false,true), n_jumps::Int64 = Inf64, saverate = false, rate::vecrate = zeros(Tc, size(nu,1))) where {Tc,Td,Tnu <: AbstractArray{Td}, Tp, TF ,TR ,TD,
		vecc <: AbstractVector{Tc},
		vecd <:  AbstractVector{Td},
		vecrate <: AbstractVector{Tc}}
	@assert 1==0 "WIP, function to be removed"

	# custom type to collect all parameters in one structure
	problem  = PDMPProblem{Tc,Td,vecc,vecd,vecrate,Tnu,Tp,TF,TR,TD}(false,xc0,xd0,rate,F,R,DX,nu,parms,ti,tf,save_positions[1],verbose,saverate)

	rejection_diffeq!(problem,ti,tf;ode = ode, save_positions = save_positions,n_jumps = n_jumps)
end

function rejection_diffeq!(problem::PDMPProblem,
				ti::Tc, tf::Tc, verbose = false; ode = Tsit5(),
				save_positions = (false,true), n_jumps::Td = Inf64, reltol=1e-7, abstol=1e-9) where {Tc, Td}
	verbose && printstyled(color=:red,"Entry in rejection_diffeq\n")
	ti, tf = problem.interval
	algopdmp = CHV(ode)

#ISSUE HERE, IF USING A PROBLEM p MAKE SURE THE TIMES in p.sim ARE WELL SET
	# set up the current time as the initial time
	t = ti
	# previous jump time, needed because problem.simjptimes.lastjumptime contains next jump time even if above tf
	tprev = t

	# we declare the characteristics for convenience
	caract = problem.caract

	# vector to hold the state space for the extended system
	X0 = similar(caract.xc,length(caract.xc))
	for ii in eachindex(caract.xc)
		X0[ii] = caract.xc[ii]
	end

	# current jump number
	njumps = 0
	problem.simjptimes.lambda_star = 0	# this is the bound for the rejection method
	problem.simjptimes.ppf .= caract.R(caract.ratecache, X0, caract.xd, t, caract.parms, true)

	# @assert problem.simjptimes.ppf[2] == problem.pdmpFunc.R(problem.rate,X0+0.1265987*cumsum(ones(length(X0))),problem.xd,t+0.124686489,problem.parms,true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"
	# problem.rate .*= 0;problem.simjptimes.ppf .= problem.pdmpFunc.R(problem.rate,X0,problem.xd,t,problem.parms,true)
	# @assert sum(problem.rate) == 0 "You cannot modify the first argument of your rate function when sum_rate = true"

	problem.simjptimes.tstop_extended = problem.simjptimes.tstop_extended / problem.simjptimes.ppf[2] + ti
	problem.simjptimes.reject = true

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, integrator -> rejectionjump(integrator, problem, save_positions[1], verbose), save_positions = (false, false))

	# define the ODE flow, this leads to big memory saving
	prob_REJ = ODEProblem((xdot,x,data,tt) -> algopdmp(xdot, x, caract, tt), X0, (ti, 1e9))
	integrator = init(prob_REJ, ode, tstops = problem.simjptimes.tstop_extended, callback = cb, save_everystep = false, reltol = reltol, abstol = abstol, advance_to_tstop = true)


	while (t < tf) && problem.simjptimes.njumps < n_jumps-1 #&& problem.simjptimes.fictitous_jumps < 10
		verbose && println("--> n = $(problem.simjptimes.njumps), t = $t, δt = ",integrator.t)
		step!(integrator)
		@assert( t < problem.simjptimes.lastjumptime, "Could not compute next jump time $(problem.simjptimes.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(problem.simjptimes.lastjumptime),\n solver = $ode")
		t, tprev = problem.simjptimes.lastjumptime, t

		# the previous step was a jump!
		if njumps < problem.simjptimes.njumps && save_positions[2] && (t <= tf) && problem.simjptimes.reject == false
			verbose && println("----> save post-jump, xd = ",problem.Xd)
			push!(problem.Xc, copy(caract.xc))
			push!(problem.Xd, copy(caract.xd))
			push!(problem.time,t)
			njumps +=1
			verbose && println("----> end save post-jump, ")
			#put the flag fpr rejection
			problem.simjptimes.reject = true
		end
	end
	# we check that the last bit [t_last_jump, tf] is not missing
	if t>tf
		verbose && println("----> LAST BIT!!, xc = ",caract.xc[end], ", xd = ",caract.xd)
		prob_last_bit = ODEProblem((xdot,x,data,tt) -> caract.F(xdot, x, caract.xd, tt, caract.parms), copy(caract.xc), (tprev,tf))
		sol = DiffEqBase.solve(prob_last_bit, ode)
		verbose && println("-------> xc[end] = ",sol.u[end])
		push!(problem.Xc, sol.u[end])
		push!(problem.Xd, copy(caract.xd))
		push!(problem.time, sol.t[end])
	end
	return PDMPResult(problem.time,problem.Xc,problem.Xd,problem.rate_hist)
end


function solve(problem::PDMPProblem{Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar}, algo::Rejection{Tode}; verbose = false, n_jumps = Inf64, save_positions = (false, true), reltol = 1e-7, abstol = 1e-9) where {Tc, Td, vectype_xc, vectype_xd, vectype_rate, Tnu, Tp, TF, TR, Tcar, Tode <: DiffEqBase.DEAlgorithm}

	return rejection_diffeq!(problem, problem.interval[1], problem.interval[2], verbose; ode = algo.ode, save_positions = save_positions, n_jumps = n_jumps, reltol = reltol, abstol = abstol )
end
