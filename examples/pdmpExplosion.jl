# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
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


function F!(ẋ, xc, xd, t, parms)
	ẋ[1] = -r * (2mod(xd[1],2)-1) * xc[1]
end

R(x) = x

function R!(rate, xc, xd, t, parms, sum_rate::Bool)
	# rate fonction
	if sum_rate==false
		rate[1] = R(xc[1])
		rate[2] = 0.0
		return R(xc[1]),40.
	else
		return R(xc[1]),40.
	end
end

xc0 = [1.0]
xd0 = [0,0]

nu = [[1 0];[0 -1]]
parms = [0.0]
ti = 0.332
tf = 100000.
nj = 50

Random.seed!(8)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

errors = Float64[]
# state of the random generator
rnd_state = 0.

println("\n\nComparison of solvers")
	for ode in [(:lsoda,"lsoda"),(:cvode,"cvode"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas5(),"rodas5"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23()),"AutoTsit5-RS23")]
		Random.seed!(8)
		res =  PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F!, R!, nu, parms, ti, tf, n_jumps = nj, ode = ode[1], verbose = false)
		# this is to check the state of the random generator at the end of the simulation
		if ode[1] == :lsoda
			global rnd_state = rand()
		else
			@assert rnd_state == rand()
		end

		println("--> norm difference = ", norm(res.time - res_a[1],Inf64), "  - solver = ",ode[2])
		push!(errors,norm(res.time - res_a[1],Inf64))
	end

#
# println("\n\nComparison of solvers for splitting the computation")
# # we compute up to time 0.5 and then up to the end
# println("\n\nComparison of solvers")
# 	for ode in [(:lsoda,"lsoda"),(:cvode,"cvode"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas5(),"rodas5"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
# 		@show ode
# 		Random.seed!(8)
# 		res1 =  PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F!, R!, nu, parms, ti, tf, n_jumps = 4, ode = ode[1], verbose = false)
#
# 		println(length(res1.time), " ", res_a[end][9])
# 		@show res1.time
#
# 		# Random.GLOBAL_RNG.idxF +=-0
#
# 		@show rand()
#
# 		res2 =  PiecewiseDeterministicMarkovProcesses.pdmp!(res1.xc[:,end], res1.xd[:,end], F!, R!, nu, parms, 0.8, tf, n_jumps = nj - 4 , ode = ode[1], verbose = false)
#
# 		tt = vcat(res1.time[1:end],res2.time[1:end])
#
# 		@show tt - res_a[1][1:length(tt)]
#
#
# 		println("--> norm difference = ", norm(tt - res_a[1],Inf64), "  - solver = ", ode[2])
# 		push!(errors,norm(tt - res_a[1],Inf64))
# 	end
