# using Revise
using PiecewiseDeterministicMarkovProcesses,DifferentialEquations, LinearAlgebra, Random
const PDMP = PiecewiseDeterministicMarkovProcesses

function F_eva!(xcdot, xc, xd, parms::Vector{Float64}, t::Float64)
	# vector field used for the continuous variable
	xcdot[1] = -(xc[1] - 1.5)
	nothing
end

function R(x)
	return x^4
end

function R_eva(rate, xc, xd, parms, t::Float64, issum::Bool)
	# rate function
	rate_print = parms[1]
	if issum == false
		if xd[1] == 0
			rate[1] = R(xc[1])
			rate[2] = 0.0
			rate[3] = rate_print
			return 0.0, 4.95 #transition 0->1
		else
			rate[1] = 0.0
			rate[2] = 1.0
			rate[3] = rate_print
			return 0.0, 4.95 #transition 1->0
		end
	else
		if xd[1] == 0
			return R(xc[1]) + rate_print, 5. #transition 0->1
		else
			return 1.0 + rate_print, 5. #transition 1->0
		end
	end
end

function Delta_xc_eva(xc, xd, parms, t::Float64, ind_reaction::Int64)
	# this function return the jump in the continuous component
	if ind_reaction == 2
		xc[1] = 0.0
	end
	return true
end

xc0 = [0.0]
xd0 = [0, 1]

nu_eva = [1 0;-1 0;0 1]
parms = [1.]
tf = 100.

println("--> Case simple chv:")
	Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_eva!,R_eva,Delta_xc_eva,nu_eva, xc0, xd0, parms, (0.0, tf))
	dummy_t =  @time PDMP.solve(problem, CHV(Tsit5()); n_jumps = 200)


println("--> For simulations rejection (Tsit5):")
	Random.seed!(123)
	problem = PDMP.PDMPProblem(F_eva!,R_eva,Delta_xc_eva,nu_eva, xc0, xd0, parms, (0.0, tf))
	result1 =  @time PDMP.solve(problem, Rejection(:lsoda); n_jumps = 200)

println("--> Simulation using save_at to see sampling behaviour")
	nj = 51
	parms = [10.0]
	Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_eva!,R_eva,Delta_xc_eva,nu_eva, xc0, xd0, parms, (0.0, tf))
	result3 =  @time PDMP.solve(problem, CHV(Tsit5()); n_jumps = 200, save_positions = (false,true))
