###################################################################################################
###################################################################################################
###################################################################################################
### WIP implementation of the CHV algo using DiffEq

@inline function (prob::PDMPProblem)(u,t,integrator)
    t == prob.tstop_extended
end

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
function (prob::PDMPProblem)(integrator)
	# find the next jump time
	t = integrator.u[end]

	# state of the continuous variable right before the jump
# HOW TO USE VIEW FOR THIS?
	X0 = integrator.u[1:end-1]

	# execute the jump
	prob.R(prob.rate,X0,prob.xd,t,prob.parms, false)
	if (t < prob.tf)
		# Update event
		ev = pfsample(prob.rate,sum(prob.rate),length(prob.rate))

		deltaxd = view(prob.nu,ev,:)

		# Xd = Xd .+ deltaxd
		LinearAlgebra.BLAS.axpy!(1.0, deltaxd, prob.xd)

		# Xc = Xc .+ deltaxc
		prob.Delta(X0,prob.xd,t,prob.parms,ev)
	end
	# we register the next time interval to solve ode
	prob.njumps += 1
	prob.tstop_extended += -log(rand())
	add_tstop!(integrator, prob.tstop_extended)
end

function (prob::PDMPProblem)(xdot,x,p,t)
	# used for the exact method
	# we put [1] to use it in the case of the rejection method as well
	tau = x[end]
	sr = prob.R(prob.rate,x,prob.xd,tau,prob.parms,true)[1]
	@assert sr > 0.0 "Total rate must be positive"
	isr = min(1.0e9,1.0 / sr)
	prob.F(xdot,x,prob.xd,tau,prob.parms)
	xdot[end] = 1.0
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] * isr
	end
	nothing
end

"""
Implementation of the CHV method to sample a PDMP using the package `DifferentialEquations`. The advantage of doing so is to lower the number of calls to `solve` using an `integrator` method.
"""
function chv_diffeq!(n_max::Int64,xc0::AbstractVector{Float64},xd0::AbstractVector{Int64},
				F::Function,R::Function,DX::Function,
				nu::AbstractArray{Int64},parms,
				ti::Float64, tf::Float64,
				verbose::Bool = false;
				ode=:euler,ind_save_d=-1:1,ind_save_c=-1:1,save_at=[])
	# @warn "Modifying parameters is fine if the Interpolant is not lazy like a Vern method\n"

	# Set up initial simulation time
	t = ti

	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)  # vector of rates

	# vector to hold the state space for the extended system
	X_extended = copy(xc0); push!(X_extended,0.)

	# custom type to collect all parameters in one structure
	problem  = PDMPProblem(copy(xc0),copy(xd0),F,R,DX,nu,parms,tf,rate,-log(rand()),0)

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, problem)


	X_hist = copy(xc0)
	t_hist = [ti]
	# define the ODE flow, this leads to big memory saving
	prob_CHV = ODEProblem((xdot,x,data,tt)->problem(xdot,x,data,tt),X_extended,(ti,tf))
	integrator = init(prob_CHV, Tsit5(), tstops = problem.tstop_extended, callback=cb, save_everystep = false, save_positions=(false,false))

	njumps = 0
# I DONT NEED TO TAKE CARE OF THE FACT THAT the last jump might be after tf	
	while (t < tf) && problem.njumps < n_max
		verbose && println("--> n = $(problem.njumps), t = $t")

		step!(integrator)
		t = integrator.u[end]

		# the previous step was a jump!
		if njumps < problem.njumps
			# need to find a good way to solve the jumps, not done YET
			# append!(X_hist,integrator.u[1:end-1])
			# push!(t_hist,t)
			njumps +=1
		end
	end
	return t_hist,X_hist,integrator.sol
end
