# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, OrdinaryDiffEq, Test
const PDMP = PiecewiseDeterministicMarkovProcesses

function AnalyticalSampleCHV(xc0, xd0, ti, nj::Int64)
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	t = ti
	while length(th)<nj
		xc = xch[end]
		xd = xdh[end]
		S = -log(rand())
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

function AnalyticalSampleRejection(xc0,xd0,ti,nj::Int64; verbose = false)
	verbose && printstyled(color=:red,"--> Start analytical method\n")
	xch = [xc0[1]]
	xdh = [xd0[1]]
	th  = [ti]
	t = ti
	xc = xc0[1]
	njumps = 1
	rt = zeros(2)
	lambda_star = R!(rt,xc0,xd0,ti,Float64[],true)[2]
	rate		= R!(rt,xc0,xd0,ti,Float64[],true)[1]
	S = -log(rand()) / lambda_star
	while njumps < nj
		xd  = [xdh[end] ,1]

		t += S
		if mod(xd[1],2) == 0
			xc = xc * exp(10 * S)#xc + S
		else
			xc = xc / (3 * S * xc + 1)#xc * exp(-10S)
		end
		verbose && println("--> S = $S, t = $t, xc = $xc, xd = $(xd[1]), λ_* = ", lambda_star)
		#reject?
		lambda_star = R!(rt,[xc],xd,ti,Float64[],true)[2]
		rate		= R!(rt,[xc],xd,ti,Float64[],true)[1]
		@assert rate <= lambda_star "Wrong bound"
		reject = rand() < (1 - rate / lambda_star)
		S = -log(rand()) / lambda_star
		if ~reject
			verbose && println("----> Jump!, ratio = ", rate / lambda_star)
			@assert rate <= lambda_star "Wrong bound"
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

function F!(ẋ, xc, xd, parms, t)
	if mod(xd[1], 2)==0
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
		return 0., parms[1] + 100.
	else
		return R(xc[1]) + parms[1], parms[1] + 100.
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
res_a_chv = AnalyticalSampleCHV(xc0,xd0,ti,nj)
Random.seed!(8)
res_a_rej = AnalyticalSampleRejection(xc0,xd0,ti,nj)

algos = [(:cvode,"cvode"),
			(:lsoda,"lsoda"),
			(CVODE_BDF(),"CVODEBDF"),
			(CVODE_Adams(),"CVODEAdams"),
			(Tsit5(),"tsit5"),
			(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),
			(Rodas4P(),"rodas4P-AutoDiff"),
			(AutoTsit5(Rosenbrock23(autodiff=true)),"AutoTsit5-RS23")]

problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
println("\n\nComparison of solvers - CHV")
for ode in algos
	Random.seed!(8)
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj, reltol = 1e-8, abstol = 1e-11)
	println("--> norm difference = ", norm(res.time - res_a_chv[1], Inf64), "  - solver = ",ode[2])
		# compare jump times
		@test norm(res.time - res_a_chv[1], Inf64) < 3e-3
		# compare xc end values
		@test norm(res.xc[end][1] - res_a_chv[2][end], Inf64) < 4e-6
end

println("\n\nComparison of solvers - CHV (without saving solution)")
	for ode in algos
	Random.seed!(8)
	res1 =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	Random.seed!(8)
	res2 =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj, save_positions = (true, false))
		@test norm(res1.time[end] - res2.time[end]) ≈ 0
	@test norm(res1.xc[end] - res2.xc[end]) ≈ 0
	if ode[1] isa Symbol
		@test norm(res1.xd[end] - res2.xd[end]) ≈ 0
	end
end

println("\n\nComparison of solvers - CHV (limited by simulation time)")
	problem.tspan[2] = 4.0
	jumpsana = res_a_chv[1][res_a_chv[1] .< problem.tspan[2]]
	for ode in algos
	Random.seed!(8)
	res1 =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
		# same without recording the intermediate jumps
	Random.seed!(8)
	res2 =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj, save_positions = (true, false))
		@test norm(res1.time[1:end-1] .- jumpsana, Inf) < 2e-5
		@test norm(res1.time[end] - res2.time[end]) ≈ 0
	@test norm(res1.xc[end] - res2.xc[end]) ≈ 0
	if ode[1] isa Symbol
		@test norm(res1.xd[end] - res2.xd[end]) ≈ 0
	end
end

# idem as above but with tf limited simulation
prob2 = deepcopy(problem)
prob2.tspan[2] = 4.0
Random.seed!(8)
res3 =  @time PDMP.solve(prob2, CHV(:lsoda); n_jumps = 200)
Random.seed!(8)
res4 =  @time PDMP.solve(prob2, CHV(:lsoda); n_jumps = 20, save_positions = (true, false) )
@test res3.time[end] ≈ res4.time[end]
@test res3.xc[end] ≈ res3.xc[end]

# using Plots
# plot(res1.time, res1.xc[:,:]')

problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
println("\n\nComparison of solvers - rejection")
	for ode in algos
	Random.seed!(8)
	res =  PDMP.solve(problem, Rejection(ode[1]); n_jumps = 4, verbose = false)
	println("--> norm difference = ", norm(res.time - res_a_rej[1][1:4], Inf64), "  - solver = ",ode[2])
	@test norm(res.time - res_a_rej[1][1:4], Inf64) < 0.0043
end

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res1 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj)
	Random.seed!(8)
	res2 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj)

println("Alternate between calls CHV - Rej")
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	res_chv = PDMP.solve(problem, CHV(Rodas4()); n_jumps = 50)
	Random.seed!(8)
	res_rej = PDMP.solve(problem, Rejection(Rodas4()); n_jumps = 4)
	println("--> norm diff (CHV) = ", norm(res_chv.time - res_a_chv[1], Inf64))
	println("--> norm diff (Rej) = ", norm(res_rej.time - res_a_rej[1][1:4], Inf64))
	@test norm(res_chv.time - res_a_chv[1], Inf64) < 1e-3
	@test norm(res_rej.time - res_a_rej[1][1:4], Inf64) < 0.0043

# here, we write the jump problem with a function
function Delta!(xc, xd, t, parms, ind_reaction::Int64)
	if ind_reaction == 1
		xd[1] += 1
	else
		xd[2] -= 1
	end
	nothing
end

problem = PDMP.PDMPProblem(F!, R!,  Delta!, 2, xc0, xd0, parms, (ti, tf))
println("\n\nComparison of solvers, with function Delta")
	for ode in algos
	Random.seed!(8)
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time - res_a_chv[1],Inf64), "  - solver = ", ode[2])
	push!(errors, norm(res.time - res_a_chv[1], Inf64))
end

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(:lsoda); n_jumps = nj)

# test for allocations, should not depend on the requested number of jumps
Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, 1e9))
	res =  PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
	alloc1 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	alloc1 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	alloc2 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = 2nj, save_positions = (false, false))
	Random.seed!(8)
	alloc3 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = 3nj, save_positions = (false, false))
	println("--> allocations = ", (alloc1, alloc2, alloc3)) #--> allocations = (58736, 13024)

# test for many calls to solve, the trajectories should be the same
problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	res = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, true))
	restime1 = res.time
	Random.seed!(8)
	res12 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, true))
