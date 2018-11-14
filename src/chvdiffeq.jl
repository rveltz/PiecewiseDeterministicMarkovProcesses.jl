###################################################################################################
###################################################################################################
###################################################################################################
### WIP implementation of the CHV algo using DiffEq

"""
a type to hold data
"""
struct DataForODE{T}
	parms::T
	xd::Vector{Int64}
	rate::Array{Float64}
end

# function f_CHV_Wrap!(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64}, p::DataForODE{Ty}) where Ty
# 	f_CHV!(F,R,t, x, xdot, p.xd, p.parms,p.rate)
# end
#
# function Flow_Wrap!(X0,Xd,dt,r,prob::ODEProblem)
# 	(prob.p).rate .= r
# 	(prob.p).xd   .= Xd
# 	prob.u0[:]    .= X0
# 	println("flow wrap done")
# 	sol = DifferentialEquations.solve(prob, DifferentialEquations.Tsit5(),save_start=false,save_end=true,save_everystep = false).u
#     return sol[end]
# end

@inline function (prob::PDMPProblem)(u,t,integrator)
    t == prob.tstop_extended
end

# The following function is a callback to discrete jump. Its role is to perform the jump on the solution given by the ODE solver
function (prob::PDMPProblem)(integrator)
	# print("----> Jump detected!, enter callback to deal with it")
	# find the next jump time
	t = integrator.u[end]
	# println(", tjump = $t")

	# state of the continuous variable right before the jump
	X0 = integrator.u[1:end-1]

	# execute the jump
	prob.R(prob.rate,X0,prob.xd,t,prob.parms, false)
	if (t < prob.tf)
		# Update event
		ev = pfsample(prob.rate,sum(prob.rate),length(prob.rate))

		deltaxd = prob.nu[ev,:]
		# Xd = Xd .+ deltaxd
		LinearAlgebra.BLAS.axpy!(1.0, deltaxd, prob.xd)

		# Xc = Xc .+ deltaxc
		prob.Delta(X0,prob.xd,t,prob.parms,ev)
	end
	# println("----> x = ",X0)
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
	@warn "Modifying parameters is fine if the Interpolant is not lazy like a Vern method\n"

	# Set up initial simulation time
	t = ti

	deltaxd = copy(nu[1,:]) # declare this variable, variable to hold discrete jump
	numpf   = size(nu,1)    # number of reactions
	rate    = zeros(numpf)  # vector of rates
	# vector to hold the state space for the extended system
	X_extended = copy(xc0); push!(X_extended,0.)

	# custom type to collect all parameters in one structure
	data_ode = DataForODE(parms,copy(xd0),rate)
	problem  = PDMPProblem(copy(xc0),copy(xd0),F,R,DX,nu,parms,tf,rate,-log(rand()),0)

	# definition of the callback structure passed to DiffEq
	cb = DiscreteCallback(problem, problem)


	X_hist = copy(xc0)
	t_hist = [ti]
	# define the ODE flow, this leads to big memory saving
	prob_CHV = ODEProblem((xdot,x,data,tt)->problem(xdot,x,data,tt),X_extended,(0.,1.))
	integrator = init(prob_CHV, Tsit5(), tstops = problem.tstop_extended, callback=cb)

	njumps = 0
	println("--> n = 0, t = $t, tf = $tf")
	while (t < tf) && problem.njumps < n_max
		verbose && println("--> n = $(problem.njumps), t = $t")
		# we compute a jump
		step!(integrator)

		t = integrator.u[end]

		if njumps < problem.njumps
			# need to find a good way to solve the jumps, not done YET
			append!(X_hist,integrator.u[1:end-1])
			push!(t_hist,t)
			njumps +=1
		end
	end
	return t_hist,X_hist
end
