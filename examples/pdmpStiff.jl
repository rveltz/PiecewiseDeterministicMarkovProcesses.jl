using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
const PDMP = PiecewiseDeterministicMarkovProcesses

function AnalyticalSample(xc0, xd0, ti, nj::Int64)
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	t = ti
	while length(th)<nj
		xc = xch[end]
		xd = xdh[end]
		S = -log(rand())
		# @show xd,S
		if mod(xd,2) == 0
			t += 1/10*log(1+10*S/xc)
			push!(xch,xc + 10 * S )
		else
			t += 1/(3xc)*(exp(3S)-1)
			push!(xch,xc * exp(-3S) )
		end
		push!(xdh,xd + 1 )
		push!(th,t)

		S = -log(rand())
	end
	return th, xch, xdh
end

function F!(ẋ, xc, xd, parms, t)
	if mod(xd[1],2)==0
		ẋ[1] = 10xc[1]
	else
		ẋ[1] = -3xc[1]^2
	end
end

R(x) = x

function R!(rate, xc, xd, parms, t, issum::Bool)
	# rate function
	if issum == false
		rate[1] = R(xc[1])
		rate[2] = parms[1]
		return 0.
	else
		return R(xc[1]) + parms[1]
	end
end

xc0 = [1.0]
	xd0 = [0, 0]

	nu = [1 0;0 -1]
	parms = [.0]
	ti = 0.322156
	tf = 100000.
	nj = 50

errors = Float64[]

Random.seed!(8)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

println("\n\nComparison of solvers")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23(autodiff=true)),"AutoTsit5-RS23")]
	Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time - res_a[1], Inf64), "  - solver = ",ode[2])
	push!(errors,norm(res.time - res_a[1],Inf64))
end


# here, we write the jump problem with a function
using SparseArrays
nusp = spzeros(Int64, 2, 2)

function Delta!(xc, xd, t, parms, ind_reaction::Int64)
	if ind_reaction == 1
		xd[1] += 1
	else
		xd[2] -= 1
	end
	nothing
end

println("\n\nComparison of solvers, with function Delta")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
	Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!,  Delta!, nusp, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time - res_a[1],Inf64), "  - solver = ", ode[2])
	push!(errors, norm(res.time - res_a[1], Inf64))
end


Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(:lsoda); n_jumps = nj, verbose = false)
