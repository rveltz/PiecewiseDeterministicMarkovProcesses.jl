# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
	const PDMP = PiecewiseDeterministicMarkovProcesses
using JumpProcesses

rate = (u,p,t) -> .1 + u[1]
affect! = (integrator) -> (integrator.u[1] = integrator.u[1]/2; integrator.u[2] +=1)
jump = JumpProcesses.VariableRateJump(rate, affect!, interp_points = 1000)
jumpprint = JumpProcesses.VariableRateJump((u,p,t) -> 10.0, x -> x.u[3] +=1, interp_points = 1000)

vf = function (du,u,p,t)
	if mod(u[2],2) == 0
	  du[1] = u[1]
	else
		du[1] = -u[1]
	end
	du[2] = 0.
	du[3] = 0.
	nothing
end

prob = ODEProblem(vf, [0.2, 0.0, 0.0], (0.0,10.0))
jump_prob = JumpProcesses.JumpProblem(prob, Direct(), jump, jumpprint)

# let us solve the PDMD with JumpProcesses
Random.seed!(123)
	soldj = @time JumpProcesses.solve(jump_prob, Tsit5())
# plot(soldj,ylims=(0, 2))

# wrapper to PDMP
pb = PDMP.PDMPProblem(prob, jump, jumpprint)
Random.seed!(123)
	solwp = @time PDMP.solve(pb, CHV(Tsit5()); save_positions = (false, true))
# plot(solwp.time, solwp.xc[1,:])
# 	plot!(solwp.time, solwp.xc[2,:], line=:step, marker=:d)
