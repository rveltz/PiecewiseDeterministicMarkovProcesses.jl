using Revise
	using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations
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
	# rate fonction
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

	nu = [1 0; 0 -1]
	parms = [0.0]
	ti = 0.322156
	tf = 100000.
	nj = 50

	errors = Float64[]

Random.seed!(8)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

Random.seed!(8) #0.000637 seconds (327 allocations: 26.219 KiB)
	pb = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res = @time PDMP.solve(pb, CHV(Tsit5()), n_jumps = nj, save_positions = (false, false))

##########################################
# build a PDMP Problem
Random.seed!(8)
	pb = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res = @time PDMP.solve(pb, CHV(Tsit5()), n_jumps = nj, save_positions = (false, true))
	norm(res.time - res_a[1],Inf64)

Random.seed!(8)
	pb = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res = PDMP.solve(pb, CHV(:lsoda), n_jumps = nj)
	norm(res.time - res_a[1],Inf64)
