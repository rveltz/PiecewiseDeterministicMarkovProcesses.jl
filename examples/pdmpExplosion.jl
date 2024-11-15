# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, OrdinaryDiffEq, Sundials
const r = 10.

function AnalyticalSample(xc0,xd0,ti,nj::Int64)
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	list_rng = Float64[]
	t = ti
	while length(th)<nj
		xc = xch[end]
		xd = xdh[end]
		push!(list_rng, rand())
		S = -log(list_rng[end])
		a = -r * (2mod(xd,2)-1)
		dt = log(a*S/xc+1)/a
		t += dt
		push!(th, t)
		push!(xch,xc + a * S )
		push!(xdh,xd .+ 1 )
		push!(list_rng, rand())
		S = -log(list_rng[end])
	end
	return th, xch, xdh, list_rng
end


function F!(ẋ, xc, xd, parms, t)
	ẋ[1] = -r * (2mod(xd[1],2)-1) * xc[1]
end

R(x) = x

function R!(rate, xc, xd, parms, t, issum::Bool)
	# rate fonction
	if issum == false
		rate[1] = R(xc[1])
		rate[2] = 0.0
		return R(xc[1]), 40.
	else
		return R(xc[1]), 40.
	end
end

xc0 = [1.0]
xd0 = [0, 0]

nu = [[1 0];[0 -1]]
parms = [0.0]
ti = 0.332
tf = 100000.
nj = 50

Random.seed!(18)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

errors = Float64[]
# state of the random generator
rnd_state = 0.

println("\n\nComparison of solvers")
	for ode in [(:lsoda,"lsoda"),
				(:cvode,"cvode"),
				(CVODE_BDF(),"CVODEBDF"),
				(CVODE_Adams(),"CVODEAdams"),
				(Tsit5(),"tsit5"),
				(Rodas4P(autodiff=true),"rodas4P-AutoDiff"),
				(Rodas5(),"rodas5"),
				(Rosenbrock23(),"RS23"),
				(AutoTsit5(Rosenbrock23()),"AutoTsit5-RS23")]
		Random.seed!(18)

		problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
		res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
		# this is to check the state of the random generator at the end of the simulation
		if ode[1] == :lsoda
			global rnd_state = rand()
		else
			@assert rnd_state == rand()
		end

		println("--> norm difference = ", norm(res.time - res_a[1],Inf64), "  - solver = ",ode[2])
		push!(errors, norm(res.time - res_a[1], Inf64))
	end
