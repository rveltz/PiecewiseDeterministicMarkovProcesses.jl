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

function F!(ẋ, xc, xd, t, parms)
	if mod(xd[1],2)==0
		ẋ[1] = 10xc[1]
	else
		ẋ[1] = -3xc[1]^2
	end
end

R(x) = x

function R!(rate, xc, xd, t, parms, sum_rate::Bool)
	# rate fonction
	if sum_rate == false
		rate[1] = R(xc[1])
		rate[2] = parms[1]
		return 0.
	else
		return R(xc[1]) + parms[1]
	end
end

xc0 = [1.0]
	xd0 = [0, 0]

	nu = [[1 0]; [0 -1]]
	parms = [0.0]
	ti = 0.322156
	tf = 100000.
	nj = 50

	errors = Float64[]

Random.seed!(8)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

Random.seed!(8) #0.000521 seconds (364 allocations: 28.313 KiB)
	res = @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F!, R!, nu, parms, ti, tf; n_jumps = nj, ode = Tsit5(), save_positions = (false, false))

res_a[1]

res.time
##########################################
# build a PDMP Problem

pb = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms)

#
# Random.seed!(8)
# 	sol = @time PDMP.solve(prob, CHV(Tsit5()), xc0, xd0, parms, (ti, tf), false; save_positions = (false, true), n_jumps = nj)
# 	norm(sol.time - res_a[1], Inf64)


const algo = PDMP.CHV(Tsit5())
# const probcache = PDMP.solve(prob, algo, xc0, xd0, parms, (0., 1.), false;return_pb = true)

using BenchmarkTools
xdot0 = similar(xc0)
rate0 = zeros(2)

xe = zeros(2)
xed = zeros(2)

# This looks good, it does not allocate
@btime (algo)($xed, $xe, $pb, 0.)
