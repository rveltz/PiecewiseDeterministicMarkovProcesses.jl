# using Revise
using PiecewiseDeterministicMarkovProcesses,DifferentialEquations, LinearAlgebra, Random
const PDMP = PiecewiseDeterministicMarkovProcesses

function F_eva!(xcdot, xc, xd, t::Float64, parms::Vector{Float64})
	# vector field used for the continuous variable
	xcdot[1] = -(xc[1] - 1.5)
	nothing
end

function R(x)
	return x^4
end

function R_eva(rate,xc, xd, t::Float64, parms, sum_rate::Bool)
	# rate function
	rate_print = parms[1]
	if sum_rate == false
		if xd[1] == 0
			rate[1] = R(xc[1])
			rate[2] = 0.0
			rate[3] = rate_print
			return 0.0, 4.95 #transition 0->1
		else
			rate[1] = 0.0
			rate[2] = 1.0
			rate[3] = rate_print
			return 0.0,4.95 #transition 1->0
		end
	else
		if xd[1] == 0
			return R(xc[1]) + rate_print,4.95 #transition 0->1
		else
			return 1.0 + rate_print,4.95 #transition 1->0
		end
	end
end

function Delta_xc_eva(xc, xd, t::Float64, parms, ind_reaction::Int64)
	# this function return the jump in the continuous component
	if ind_reaction == 2
		xc[1] = 0.0
	end
	return true
end

xc0 = vec([0.0])
xd0 = vec([0, 1])

nu_eva = [[1 0];[-1 0];[0 1]]
parms = [1.]
tf = 100.

println("--> Case simple chv:")
	dummy_t =  PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,n_jumps=1)
	Random.seed!(1234)
	dummy_t =  @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,n_jumps=200000)

println("For simulations (lsoda):")
	result = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=1)
	Random.seed!(1234)
	result = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=200000)

println("--> Case tauleap:")
	resultt = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=1,algo=:tauleap)
	Random.seed!(1234)
	resultt = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=20000,algo=:tauleap,dt=0.01)

println("--> For simulations (Tsit5):")
	result1 = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=Tsit5(),n_jumps=1)
	Random.seed!(1234)
	result1 = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=Tsit5(),n_jumps=200)

println("--> For simulations rejection (Tsit5):")
	result1 = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=Tsit5(),n_jumps=1,algo=:rejection)
	Random.seed!(1234)
	result1 = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=Tsit5(),n_jumps=200,saverate=true,algo=:rejection)

println("--> Simulation using save_at to see sampling behaviour")
	nj = 51
	parms = [10.0]
	Random.seed!(1234)
	result3 = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.4,11.,ode=Tsit5(),n_jumps=nj, save_positions = (false,true))

	Random.seed!(1234)
	result4 = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.4,11.,ode=Tsit5(),n_jumps=nj, save_positions = (true,false))

	Random.seed!(1234)
	result5 = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.4,11.,ode=Tsit5(),n_jumps=nj, save_positions = (true,true))
