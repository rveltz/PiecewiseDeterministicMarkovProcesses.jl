# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
const PDMP = PiecewiseDeterministicMarkovProcesses

function F_tcp!(ẋ, xc, xd, parms, t)
	# vector field used for the continuous variable
	if mod(xd[1],2) == 0
		 ẋ[1] = 1.0
	else
		 ẋ[1] = -10.0 * xc[1]
	end
	nothing
end

rate_tcp(x) = 1/(1+exp(-x))

function R_tcp!(rate, xc, xd, parms, t, issum::Bool)
	if issum == false
		rate[1] = rate_tcp(xc[1])
		rate[2] = 0.0
		return rate_tcp(xc[1]), 1.0
	else
		return rate_tcp(xc[1]), 1.0
	end
end

function AnalyticalSample(xc0,xd0,ti,nj::Int64; verbose = false)
	verbose && printstyled(color=:red,"--> Start analytical method\n")
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	t = ti
	xc = xc0[1]
	njumps = 1
	rt = zeros(2)
	lambda_star = R_tcp!(rt,xc0,xd0,ti,Float64[],true)[2]
	rate		= R_tcp!(rt,xc0,xd0,ti,Float64[],true)[1]
	S = -log(rand()) / lambda_star
	while njumps < nj
		xd  = [xdh[end] ,1]

		t += S
		if mod(xd[1],2) == 0
			xc = xc + S
		else
			xc = xc * exp(-10S)
		end
		verbose && println("--> S = $S, t = $t, xc = $xc, xd = $(xd[1]), λ_* = ", lambda_star)
		#reject?
		lambda_star = R_tcp!(rt,[xc],xd,ti,Float64[],true)[2]
		rate		= R_tcp!(rt,[xc],xd,ti,Float64[],true)[1]
		reject = rand() < (1 - rate / lambda_star)
		S = -log(rand()) / lambda_star
		if ~reject
			verbose && println("----> Jump!, ratio = ",rate / lambda_star)
			push!(th,t)
			push!(xch,xc)
			push!(xdh,xdh[end] + 1)
			njumps += 1
			# dummy call to rand to emulate sampling pfsample
			dum = -log(rand())
		end
	end
	return th, xch, xdh
end

xc0 = [ 0.0 ]
xd0 = [0, 0]

nu_tcp = [1 0;0 -1]
parms = [0.0]
tf = 100000.
nj = 50
errors = Float64[]

Random.seed!(1234)
res_a = AnalyticalSample(xc0, xd0, 0.0, nj, verbose=false)

println("\n\nComparison of solvers")
for ode in [(:cvode,"cvode"),
			(:lsoda,"lsoda"),
			(CVODE_BDF(),"CVODEBDF"),
			(CVODE_Adams(),"CVODEAdams"),
			(Tsit5(),"tsit5"),
			(Rodas4P(autodiff=true),"rodas4P-AutoDiff"),
			(Rodas4P(),"rodas4P-AutoDiff"),
			(Rosenbrock23(),"RS23"),
			(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
	Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
	res =  PDMP.solve(problem, Rejection(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time[1:nj] - res_a[1],Inf64), "  - solver = ", ode[2])
	push!(errors, norm(res.xc[1,1:nj] - res_a[2], Inf64))
end

println("test for allocations, should not depend on")
Random.seed!(1234)
problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
res =  PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj, save_positions = (false, false))
Random.seed!(1234)
res =  @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj, save_positions = (false, false))
Random.seed!(1234)
res =  @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 2nj, save_positions = (false, false))
Random.seed!(1234)
res =  @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 3nj, save_positions = (false, false))

println("test for multiple calls, the result should not depend on")
Random.seed!(1234)
problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
res1 =  PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj)
res2 =  PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj)
@assert res1.time != res2.time
Random.seed!(1234)
res1 =  PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj)
Random.seed!(1234)
res2 =  PDMP.solve(problem, Rejection(Tsit5()); n_jumps = nj)
@assert res1.time == res2.time

# Random.seed!(1234)
# 	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
# 	alloc1 =  @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 2nj, save_positions = (false, false))
#
# Random.seed!(1234)
# 	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
# 	alloc2 =  @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 4nj, save_positions = (false, false))
#
# Random.seed!(1234)
# 	PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
# 	res = @time PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 4nj, save_positions = (false, false))
#
# Random.seed!(1234)
# 	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
# 	res =  PDMP.solve(problem, Rejection(:lsoda); n_jumps = nj, save_positions = (false, false), save_rate = true)


# test the number of rejected jumps
Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
	res1 =  PDMP.solve(problem, Rejection(:cvode); n_jumps = nj)
Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
	res2 =  PDMP.solve(problem, Rejection(CVODE_BDF()); n_jumps = nj)
