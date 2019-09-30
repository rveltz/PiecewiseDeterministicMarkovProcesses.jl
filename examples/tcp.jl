# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
const PDMP = PiecewiseDeterministicMarkovProcesses

function AnalyticalSample(xc0,xd0,ti,nj::Int64)
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	t = ti
	while length(th)<nj
		xc = xch[end]
		xd = xdh[end]
		S = -log(rand())
		a = mod(xd[1],2)==0 ? -1 : 1
		dt = (exp(a*S)-1)*exp(-a*S)/(a*xc)
		t += dt
		push!(th, t)
		push!(xch, xc * exp(a*S) )
		push!(xdh, xd .+ 1 )
		S = -log(rand())
	end
	return th,xch,xdh
end

function F_tcp!(ẋ, xc, xd, parms, t)
	# vector field used for the continuous variable
	if mod(xd[1], 2)==0
		 ẋ[1] = 1.
	else
		 ẋ[1] = -1.
	end
	nothing
end

rate_tcp(x) = 1/x

function R_tcp!(rate, xc, xd, parms, t, issum::Bool)
	if issum==false
		rate[1] = rate_tcp(xc[1])
		rate[2] = 0.0
		return 0., 100.
	else
		return rate_tcp(xc[1]), 100.
	end
end

xc0 = [ 1.0 ]
xd0 = [0, 1]

nu_tcp = [1 0;0 -1]
parms = [0.0]
tf = 100000.
nj = 100

Random.seed!(1234)
	res_a = AnalyticalSample(xc0,xd0,0.,nj)

errors = Float64[]

println("\n\nComparison of solvers")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Rosenbrock23(),"RS23"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas5(),"rodas5"),(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
	Random.seed!(1234)
	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
	res =  @time PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	printstyled(color=:green,"--> norm difference = ", norm(res.time - res_a[1],Inf64), "  - solver = ",ode[2],"\n\n")
	push!(errors,norm(res.time - res_a[1],Inf64))
end

# plot!(res.time,res.xc[:,1])
# case with no allocations  0.000721 seconds (330 allocations: 26.266 KiB)
# Random.seed!(1234)
# 	problem = PDMP.PDMPProblem(F_tcp!, R_tcp!, nu_tcp, xc0, xd0, parms, (0.0, tf))
# 	res =  @time PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))

# res =  @timed PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
# res[end].poolalloc

# # Random.seed!(1234)
# #	 using PiecewiseDeterministicMarkovProcesses
# #	 println("\n"*"+"^40)
# #	 res =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = 10,   ode =Tsit5(), algo=:rejection, verbose=true)
# #
# #	 res.time |> println
