using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, Sundials

function R_sir!(rate,xc,xd,parms,t::Float64,issum::Bool)
	(S,I,R,~) = xd
	(beta,mu) = parms
	infection = beta*S*I
	recovery = mu*I
	rate_display = 0.01
	if issum == false
			rate[1] = infection
			rate[2] = recovery
			rate[3] = rate_display
			return 0.
	else
		return infection+recovery + rate_display
	end
end

function F_sir!(xdot,xc,xd,parms,t::Float64)
	# vector field used for the continuous variable
	xdot[1] = 0.0
	nothing
end

xc0 = [0.0]
xd0 = [99,10,0,0]
nu = [[-1 1 0 0];[0 -1 1 0];[0 0 0 1]]
parms = [0.1/100.0,0.01]
tf = 150.0

Random.seed!(1234)
problem = PDMP.PDMPProblem(F_sir!,R_sir!,nu, xc0, xd0, parms, (0.0, tf))
result = PDMP.solve(problem, CHV(Tsit5()); n_jumps = 1000)
result = PDMP.solve(problem, CHV(:cvode); n_jumps = 1000)
result = PDMP.solve(problem, CHV(CVODE_BDF()); n_jumps = 1000)
