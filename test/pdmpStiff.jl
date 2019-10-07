using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials, Test
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
		return 0., parms[1] + 50.
	else
		return R(xc[1]) + parms[1], parms[1] + 50.
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

problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
println("\n\nComparison of solvers")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23(autodiff=true)),"AutoTsit5-RS23")]
	Random.seed!(8)
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time - res_a_chv[1], Inf64), "  - solver = ",ode[2])
	@test norm(res.time - res_a_chv[1],Inf64) < 1e3
end

println("\n\nComparison of solvers - rejection")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23(autodiff=true)),"AutoTsit5-RS23")]
	Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, Rejection(ode[1]); n_jumps = 4, verbose = false)
	println("--> norm difference = ", norm(res.time - res_a_rej[1][1:4], Inf64), "  - solver = ",ode[2])
end

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res1 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj)
	Random.seed!(8)
	res2 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj)


println("Alternate between calls CHV - Rej")
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	res_chv = PDMP.solve(problem, CHV(Tsit5()); n_jumps = 50)
	Random.seed!(8)
	res_rej = PDMP.solve(problem, Rejection(Tsit5()); n_jumps = 4)
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

println("\n\nComparison of solvers, with function Delta")
	for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
	Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!,  Delta!, 2, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(ode[1]); n_jumps = nj)
	println("--> norm difference = ", norm(res.time - res_a_chv[1],Inf64), "  - solver = ", ode[2])
	push!(errors, norm(res.time - res_a_chv[1], Inf64))
end

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res =  PDMP.solve(problem, CHV(:lsoda); n_jumps = nj)

# test for allocations, should not depend on the requested number of jumps
Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	alloc1 =  PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
	alloc1 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	alloc2 =  @allocated PDMP.solve(problem, CHV(Tsit5()); n_jumps = 2nj, save_positions = (false, false))
	println("--> allocations = ", (alloc1, alloc2))

# test for many calls to solve, the trajectories should be the same
problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	res = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, true))
	restime1 = res.time
	Random.seed!(8)
	res12 = PDMP.solve(problem, CHV(Tsit5()); n_jumps = nj, save_positions = (false, true))
