using PiecewiseDeterministicMarkovProcesses, Random, DifferentialEquations
	const PDMP = PiecewiseDeterministicMarkovProcesses

function F_tcp!(ẋ, xc, xd, parms, t)
	if mod(xd[1],2)==0
		 ẋ[1] =  1.0
		 ẋ[2] = -1.0 * xc[2]
	else
		 ẋ[1] = -1.0 * xc[1]
		 ẋ[2] =  1.0
	end
	nothing
end

R(x) = x

function R_tcp!(rate, xc, xd, parms, t, issum::Bool)
		rate[1] = R(xc[1]) +  R(xc[2])
		rate[2] = parms[1] * xd[1] * xc[1]
	if issum == false
		return 0.
	else
		return sum(rate)
	end
end

xc0 = [0.05, 0.075]
xd0 = [0, 1]

nu_tcp = [1 0;0 -1]
parms = [0.1]
tf = 10000.0
nj = 1000


Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
	result1 = @time PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))

Random.seed!(1234)
	result2 = PDMP.solve(problem, CHV(:cvode); n_jumps = nj, save_positions = (false, false))

Random.seed!(1234)
	result3 = PDMP.solve(problem, CHV(:lsoda); n_jumps = 2, save_positions = (false, false))

#test auto-differentiation
# result4 = @time PDMP.solve(problem, CHV(Rodas5()))

####################################################################################################
# DEBUG DEBUG
#
# algo = CHV(Tsit5())
# xd1 = zeros(Float64, length(xc0)+1)
# xdd1 = similar(xd1)
#
# J = zeros(3,3)
# algo(xdd1,xd1,problem.caract,0.)
#
# vf = (dx,x) -> algo(dx,x,problem.caract,0.)
# #it works!
# vf(xdd1,xd1)
#
# ForwardDiff.jacobian!(J, vf, xdd1, xd1)
