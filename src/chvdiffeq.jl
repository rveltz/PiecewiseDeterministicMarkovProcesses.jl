###################################################################################################
###################################################################################################
###################################################################################################
### implementation of the CHV algo using DiffEq

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
# callable struct
function chvjump(integrator, prob::PDMPProblem)#(prob::PDMPProblem)(integrator)
	# find the next jump time
	t = integrator.u[end]
	prob.sim.lastjumptime = t
	prob.verbose && printstyled(color=:green,"--> Jump detected at t = $t !!\n")

	prob.verbose && printstyled(color=:green,"--> jump not yet performed, xd = ",prob.xd,"\n")

	if (prob.save_pre_jump) && (t <= prob.tf)
		prob.verbose && printstyled(color=:green,"----> save pre-jump\n")
		push!(prob.Xc, copy(prob.xc))
		push!(prob.Xd, copy(prob.xd))
		push!(prob.time,t)
	end

	# execute the jump
	prob.pdmpFunc.R(prob.rate, integrator.u, prob.xd, t, prob.parms, false)
	if (t < prob.tf)
		#save rates for debugging
		prob.save_rate && push!(prob.rate_hist, sum(prob.rate))

		# Update event
		ev = pfsample(prob.rate,sum(prob.rate),length(prob.rate))

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
	end
	prob.verbose && printstyled(color=:green,"--> jump effectued, xd = ",prob.xd,"\n")
	# we register the next time interval to solve the extended ode
	prob.sim.njumps += 1
	prob.sim.tstop_extended += -log(rand())
	add_tstop!(integrator, prob.sim.tstop_extended)
	prob.verbose && printstyled(color=:green,"--> End jump\n\n")
end

"""
Implementation of the CHV method to sample a PDMP using the package `DifferentialEquations`. The advantage of doing so is to lower the number of calls to `solve` using an `integrator` method.
"""
function chv_diffeq!(xc0::vecc, xd0::vecd,
		F::TF, R::TR, DX::TD,
		nu::Tnu, parms::Tp,
		ti::Tc, tf::Tc,
		verbose::Bool = false;
		ode = Tsit5(), n_jumps::Int64 = Inf64, save_at::vecc = [], save_positions = (false, true), saverate = false) where {Tc,Td,Tnu <: AbstractArray{Td}, Tp, TF ,TR ,TD,
		vecc <: AbstractVector{Tc},
		vecd <: AbstractVector{Td}}

	# custom type to collect all parameters in one structure
	problem  = PDMPProblem{Tc,Td,vecc,vecd,Tnu,Tp,TF,TR,TD}(xc0,xd0,F,R,DX,nu,parms,ti,tf,save_positions[1],verbose,saverate)

	chv_diffeq!(problem, ti, tf, save_at, verbose; ode = ode, save_positions = save_positions, n_jumps = n_jumps)
end


function chv_diffeq!(problem::PDMPProblem,
				ti::Tc, tf::Tc, verbose = false; ode = Tsit5(),
				save_positions = (false, true), n_jumps::Td = Inf64) where {Tc, Td}
	problem.verbose && printstyled(color=:red,"Entry in chv_diffeq\n")

	# set up the current time as the initial time
	t = ti

	# vector to hold the state space for the extended system
	X_extended = similar(problem.xc,length(problem.xc)+1)
	for ii in eachindex(problem.xc)
		X_extended[ii] = problem.xc[ii]
	end
	X_extended[end] = ti

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, integrator -> chvjump(integrator, problem), save_positions = (false,false))

	# define the ODE flow, this leads to big memory saving
	prob_CHV = ODEProblem((xdot,x,data,tt)->problem(xdot,x,data,tt),X_extended,(0.0,1000_000_000.0))
	integrator = init(prob_CHV, ode, tstops = problem.sim.tstop_extended, callback = cb, save_everystep = false, reltol=1e-7, abstol=1e-9, advance_to_tstop = true)

	# current jump number
	njumps = 0

	while (t < tf) && problem.sim.njumps < n_jumps-1
		problem.verbose && println("--> n = $(problem.sim.njumps), t = $t, δt = ",problem.sim.tstop_extended)
		step!(integrator)
		@assert( t < problem.sim.lastjumptime, "Could not compute next jump time $(problem.sim.njumps).\nReturn code = $(integrator.sol.retcode)\n $t < $(problem.sim.lastjumptime),\n solver = $ode. dt = $(t - problem.sim.lastjumptime)")
		t = problem.sim.lastjumptime

		# the previous step was a jump! should we save it?
		if njumps < problem.sim.njumps && save_positions[2] && (t <= problem.tf)
			problem.verbose && println("----> save post-jump, xd = ",problem.Xd)
			push!(problem.Xc, copy(problem.xc))
			push!(problem.Xd, copy(problem.xd))
			push!(problem.time,t)
			njumps +=1
			problem.verbose && println("----> end save post-jump, ")
		end
	end
	# we check that the last bit [t_last_jump, tf] is not missing
	if t>tf
		problem.verbose && println("----> LAST BIT!!, xc = ",problem.Xc[end], ", xd = ",problem.xd)
		prob_last_bit = ODEProblem((xdot,x,data,tt) -> problem.pdmpFunc.F(xdot, x, problem.xd, tt, problem.parms), problem.Xc[end],(problem.time[end],tf))
		sol = solve(prob_last_bit, ode)
		push!(problem.Xc,sol.u[end])
		push!(problem.Xd,copy(problem.xd))
		push!(problem.time,sol.t[end])
	end
	return PDMPResult(problem.time, problem.Xc, problem.Xd, problem.rate_hist)
end

function chv_diffeq!(problem::PDMPProblem,
				ti::Tc, tf::Tc, save_at::vecc, verbose = false; ode = Tsit5(),
				save_positions = (false,true), n_jumps::Td = Inf64) where {Tc, Td, vecc <: AbstractVector{Tc}}
	# if isempty(save_at)
		return chv_diffeq!(problem, ti, tf, verbose; ode = ode, save_positions = save_positions, n_jumps = n_jumps)
	# end

	filter_saveat!(save_at, ti, tf)
	t = ti

	for tnext in save_at
		@show (t, tnext)
		res = chv_diffeq!(problem, t, tnext, verbose, ode = ode,	save_positions = save_positions, n_jumps = n_jumps)
		problem.xc .= res.xc[:,end]
		problem.xd .= res.xd[:,end]
		problem.sim.lastjumptime = tnext
		problem.sim.tstop_extended = log(rand())
		t = tnext
	end
	return PDMPResult(problem.time, problem.Xc, problem.Xd, problem.rate_hist)
end
