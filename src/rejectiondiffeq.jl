###################################################################################################
include("jumps.jl")
include("utils.jl")

struct Rejection{Tode} <: AbstractCHVIterator
	ode::Tode	# ODE solver to use for the flow in between jumps
end
###################################################################################################
### implementation of the rejection algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function rejectionjump(integrator,prob::PDMPProblem)
	prob.verbose && printstyled(color=:red,"--> REJECTION JUMP\n")
	# find the next jump time
	t = integrator.t
	prob.sim.lastjumptime = t
	prob.verbose && printstyled(color=:green,"--> Fictitous jump at t = $t, # = ",prob.sim.fictitous_jumps," !!\n")

	prob.sim.ppf .= prob.pdmpFunc.R(prob.rateCache,integrator.u,prob.xd,t,prob.parms,true)
	@assert prob.sim.ppf[1] < prob.sim.ppf[2] "Error, your bound on the rates is not high enough!, $(prob.sim.ppf)"
	prob.sim.reject = rand() < (1 - prob.sim.ppf[1] / prob.sim.ppf[2])

	prob.verbose && printstyled(color=:green,"----> xc = ",prob.xc ,", xd = ",prob.xd,", reject = ",prob.sim.reject,", rates = ",prob.sim.ppf,"\n")

	# execute the jump
	if t < prob.tf && prob.sim.reject == false
		prob.verbose && printstyled(color=:blue,"--> TRUE jump!\n")
		prob.sim.ppf .= prob.pdmpFunc.R(prob.rateCache,prob.xc,prob.xd,t,prob.parms, false)

		if (prob.save_pre_jump) && (t <= prob.tf)
			prob.verbose && printstyled(color=:green,"----> save pre-jump\n")
			push!(prob.Xc, (integrator.u[1:end-1]))
			push!(prob.Xd, copy(prob.xd))
			push!(prob.time,t)
		end

		#save rates for debugging
		prob.save_rate && push!(prob.rate_hist, sum(prob.rateCache))

		# Update event
		ev = pfsample(prob.rateCache, sum(prob.rateCache), length(prob.rateCache))

		deltaxd = view(prob.nu,ev,:)

		# Xd = Xd .+ deltaxd
		# LinearAlgebra.BLAS.axpy!(1.0, deltaxd, prob.xd)
		@inbounds for ii in eachindex(prob.xd)
			prob.xd[ii] += deltaxd[ii]
		end

		# Xc = Xc .+ deltaxc
		prob.pdmpFunc.Delta(integrator.u,prob.xd,t,prob.parms,ev)
		u_modified!(integrator,true)
		@inbounds for ii in eachindex(prob.xc)
			prob.xc[ii] = integrator.u[ii]
		end

		prob.sim.njumps += 1
	else
		prob.sim.fictitous_jumps += 1
	end
	prob.verbose && printstyled(color=:green,"----> jump effectued, xd = ",prob.xd,"\n")
	# we register the next time interval to solve the extended ode
	prob.sim.tstop_extended += -log(rand()) / prob.sim.ppf[2]
	add_tstop!(integrator, prob.sim.tstop_extended)
	prob.verbose && printstyled(color=:green,"--> End jump\n\n")
end

"""
Implementation of the rejection method to sample a PDMP using the package `DifferentialEquations`. The advantage of doing so is to lower the number of calls to `solve` using an `integrator` method.
"""
function rejection_diffeq!(xc0::vecc, xd0::vecd,
				F::TF, R::TR, DX::TD,
				nu::Tnu, parms::Tp,
				ti::Tc, tf::Tc,
				verbose::Bool = false;
				ode = Tsit5(), save_positions=(false,true), n_jumps::Int64 = Inf64, saverate = false, rate::vecrate = zeros(Tc, size(nu,1))) where {Tc,Td,Tnu <: AbstractArray{Td}, Tp, TF ,TR ,TD,
				vecc <: AbstractVector{Tc},
				vecd <:  AbstractVector{Td},
				vecrate <: AbstractVector{Tc}}

				# custom type to collect all parameters in one structure
				problem  = PDMPProblem{Tc,Td,vecc,vecd,vecrate,Tnu,Tp,TF,TR,TD}(false,xc0,xd0,rate,F,R,DX,nu,parms,ti,tf,save_positions[1],verbose,saverate)

				rejection_diffeq!(problem,ti,tf;ode = ode, save_positions = save_positions,n_jumps = n_jumps)
end

function rejection_diffeq!(problem::PDMPProblem,
				ti::Tc, tf::Tc, verbose = false;ode=Tsit5(),
				save_positions = (false,true), n_jumps::Td = Inf64) where {Tc,Td}
	problem.verbose && printstyled(color=:red,"Entry in rejection_diffeq\n")

#ISSUE HERE, IF USING A PROBLEM p MAKE SURE THE TIMES in p.sim ARE WELL SET
	t = ti
	# previous jump time, needed because problem.sim.lastjumptime contains next jump time even if above tf
	tprev = t

	# vector to hold the state space for the extended system
	X0 = similar(problem.xc,length(problem.xc))
	for ii in eachindex(problem.xc)
		X0[ii] = problem.xc[ii]
	end

	# current jump number
	njumps = 0
	problem.sim.lambda_star = 0	# this is the bound for the rejection method
	problem.sim.ppf .= problem.pdmpFunc.R(problem.rateCache,X0,problem.xd,t,problem.parms,true)

	# @assert problem.sim.ppf[2] == problem.pdmpFunc.R(problem.rate,X0+0.1265987*cumsum(ones(length(X0))),problem.xd,t+0.124686489,problem.parms,true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"
	# problem.rate .*= 0;problem.sim.ppf .= problem.pdmpFunc.R(problem.rate,X0,problem.xd,t,problem.parms,true)
	# @assert sum(problem.rate) == 0 "You cannot modify the first argument of your rate function when sum_rate = true"

	problem.sim.tstop_extended = problem.sim.tstop_extended / problem.sim.ppf[2] + ti
	problem.sim.reject = true

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, integrator -> rejectionjump(integrator,problem), save_positions = (false,false))

	# define the ODE flow, this leads to big memory saving
	prob_REJ = ODEProblem((xdot,x,data,tt)->problem(xdot,x,data,tt),X0,(ti,1_000_000_000))
	integrator = init(prob_REJ, ode, tstops = problem.sim.tstop_extended, callback=cb, save_everystep = false, reltol=1e-7, abstol=1e-9, advance_to_tstop=true)


	while (t < tf) && problem.sim.njumps < n_jumps-1 #&& problem.sim.fictitous_jumps < 10
		problem.verbose && println("--> n = $(problem.sim.njumps), t = $t, δt = ",integrator.t)
		step!(integrator)
		@assert( t < problem.sim.lastjumptime, "Could not compute next jump time $(problem.sim.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(problem.sim.lastjumptime),\n solver = $ode")
		t, tprev = problem.sim.lastjumptime, t

		# the previous step was a jump!
		if njumps < problem.sim.njumps && save_positions[2] && (t <= problem.tf) && problem.sim.reject == false
			problem.verbose && println("----> save post-jump, xd = ",problem.Xd)
			push!(problem.Xc,copy(problem.xc))
			push!(problem.Xd,copy(problem.xd))
			push!(problem.time,t)
			njumps +=1
			problem.verbose && println("----> end save post-jump, ")
			#put the flag fpr rejection
			problem.sim.reject = true
		end
	end
	# we check that the last bit [t_last_jump, tf] is not missing
	if t>tf
		problem.verbose && println("----> LAST BIT!!, xc = ",problem.xc[end], ", xd = ",problem.xd)
		prob_last_bit = ODEProblem((xdot,x,data,tt)->problem.pdmpFunc.F(xdot,x,problem.xd,tt,problem.parms),
					copy(problem.xc),(tprev,tf))
		sol = solve(prob_last_bit, ode)
		problem.verbose && println("-------> xc[end] = ",sol.u[end])
		push!(problem.Xc,sol.u[end])
		push!(problem.Xd,copy(problem.xd))
		push!(problem.time,sol.t[end])
	end
	return PDMPResult(problem.time,problem.Xc,problem.Xd,problem.rate_hist)
end
