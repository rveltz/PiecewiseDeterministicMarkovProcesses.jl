# using Revise, Test
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
	const PDMP = PiecewiseDeterministicMarkovProcesses

function F!(ẋ, xc, xd, parms, t)
	if mod(xd[1], 2) == 0
		ẋ[1] = 10xc[1]
	else
		ẋ[1] = -3xc[1]^2
	end
end

R(x) = 10.0 + x

function R!(rate, xc, xd, parms, t, issum::Bool)
	if issum == false
		rate[1] = R(xd[1])
		rate[2] = parms[1]
		return 0., 0.
	else
		r = R(xd[1]) + parms[1]
		return r, r
	end
end

xc0 = [1.0]
	xd0 = [0, 0]

	nu = [1 0;0 -1]
	parms = [.0]
	ti = 0.0
	tf = 100000.
	nj = 14

# test the different way to write rates
rate0 = zeros(2)
Rvar = PDMP.VariableRate(R!)
Rcst = PDMP.ConstantRate(R!)

out1 = R!(rate0, [1.0], [2,0], parms, 0., true)
	outv = Rvar(rate0, [1.0], [2,0], parms, 0., true)
	@test out1 == outv
	outc = Rcst(rate0, [10.0], [2,0], parms, 0., true)
	@test out1 == outc
	outc = Rcst(rate0, [10.0], [3,0], parms, 0., true)
	out1 = R!(rate0, [1.0], [3,0], parms, 0., true)
	@test out1 != outc
	outc = Rcst(rate0, [10.0], [3,0], parms, 0., false)
	outc = Rcst(rate0, [10.0], [3,0], parms, 0., true)
	@test out1 == outc

algo = CHV(Tsit5())

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)

# using Plots
# plot(res0.time, res0.xc[1,:], marker = :d)

# here the rate function is constant in between jumps
pbcst = PDMP.PDMPProblem(F!, Rcst, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescst = PDMP.solve(pbcst, algo; n_jumps = nj)
	@test rescst.time == res0.time
	# plot!(rescst.time, rescst.xc[1,:])

# here the rate function is constant in between jumps
pbvar = PDMP.PDMPProblem(F!, Rvar, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	resvar = PDMP.solve(pbvar, algo; n_jumps = nj)
	@test resvar.time == res0.time
	# plot!(resvar.time, resvar.xc[1,:])

# using PrettyTables
# datat = hcat(res0.time, resvar.time, rescst.time)
# datax = hcat(res0.xc', resvar.xc', rescst.xc')
#
# PrettyTables.pretty_table(datat)
# 	PrettyTables.pretty_table(datax)


# test for allocations
Random.seed!(8)
	res0 = PDMP.solve(problem, algo; n_jumps = 3nj, save_positions = (false, false))
	resvar = PDMP.solve(pbvar, algo; n_jumps = 3nj, save_positions = (false, false))
	rescst = PDMP.solve(pbcst, algo; n_jumps = 3nj, save_positions = (false, false))

	Random.seed!(8)
	res0 = @timed PDMP.solve(problem, algo; n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	resvar = @timed PDMP.solve(pbvar, algo; n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	rescst = @timed PDMP.solve(pbcst, algo; n_jumps = nj, save_positions = (false, false))
	alloc0 = res0[5].poolalloc
	allocvar = resvar[5].poolalloc
	alloccst = rescst[5].poolalloc
	Random.seed!(8)
	res0 = @timed PDMP.solve(problem, algo; n_jumps = 3nj, save_positions = (false, false))
	Random.seed!(8)
	resvar = @timed PDMP.solve(pbvar, algo; n_jumps = 3nj, save_positions = (false, false))
	Random.seed!(8)
	rescst = @timed PDMP.solve(pbcst, algo; n_jumps = 3nj, save_positions = (false, false))
	# @test res0[5].poolalloc == alloc0
	# @test resvar[5].poolalloc == allocvar
	# @test rescst[5].poolalloc == alloccst

# test for the end time
tf = 0.6
Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)
	# plot(res0.time, res0.xc[1,:], marker = :d)

pbcst = PDMP.PDMPProblem(F!, Rcst, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescst = PDMP.solve(pbcst, algo; n_jumps = nj)
	@test rescst.time == res0.time
	# plot!(rescst.time, res0.xc[1,:])

pbvar = PDMP.PDMPProblem(F!, Rvar, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	resvar = PDMP.solve(pbvar, algo; n_jumps = nj)
	@test resvar.time == res0.time
	# plot!(resvar.time, resvar.xc[1,:])

Random.seed!(8)
	res0 =  PDMP.solve(problem, CHV(:lsoda); n_jumps = nj)
	Random.seed!(8)
	rescst = PDMP.solve(pbcst, CHV(:lsoda); n_jumps = nj)
	@test rescst.time == res0.time
	Random.seed!(8)
	resvar = PDMP.solve(pbvar, CHV(:lsoda); n_jumps = nj)
	@test resvar.time == res0.time
