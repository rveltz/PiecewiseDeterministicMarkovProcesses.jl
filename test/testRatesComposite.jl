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
		rate[1] = R(xc[1])
		rate[2] = xd[2]
		return 0., 0.
	else
		r = R(xc[1]) + xd[2]
		return r, r
	end
end

function Rcst!(rate, xc, xd, parms, t, issum::Bool)
	if issum == false
		# rate[1] = R(xc[1])
		rate[2] = xd[2]
		return 0., 0.
	else
		r = xd[2]
		return r, r
	end
end

function Rvar!(rate, xc, xd, parms, t, issum::Bool)
	if issum == false
		rate[1] = R(xc[1])
		# rate[2] = xd[2]
		return 0., 0.
	else
		r = R(xc[1])
		return r, r
	end
end

xc0 = [1.0]
	xd0 = [0, 0]

	nu = [1 0;0 -1]
	parms = [.0]
	ti = 0.0
	tf = 100000.
	nj = 20

# test the different way to write rates
rate0 = zeros(2)
Rvar = PDMP.VariableRate(R!)
Rcmp = PDMP.CompositeRate(Rcst!, Rvar!)

# using BenchmarkTools
# @btime R!($rate0, $[10.0], $[3,0], $parms, 0., false)
# @btime Rvar!($rate0, $[10.0], $[3,0], $parms, 0., false)
# @btime Rcst!($rate0, $[10.0], $[3,0], $parms, 0., false)
# @btime $Rvar($rate0, $[10.0], $[3,0], $parms, 0., false)
# @btime $Rcmp($rate0, $[10.0], $[3,0], $parms, 0., false)
#
# @btime $Rvar($rate0, $[10.0], $[3,0], $parms, 0., true)
# @btime $Rcmp($rate0, $[10.0], $[3,0], $parms, 0., true)

out1 =   R!(rate0, [1.0], [2,0], parms, 0., true)
outv = Rvar(rate0, [1.0], [2,0], parms, 0., true)
@test out1 == outv
outc = Rcmp(rate0, [1.0], [2,0], parms, 0., true)
@test out1 == outc
outc = Rcmp(rate0, [10.0], [3,0], parms, 0., true)
out1 =   R!(rate0,  [1.0], [3,0], parms, 0., true)
@test out1 != outc
outc = Rcmp(rate0, [10.0], [3,0], parms, 0., false)
outc = Rcmp(rate0, [10.0], [3,0], parms, 0., true)
out1 =   R!(rate0, [10.0], [3,0], parms, 0., true)
@test out1 == outc

algo = CHV(Tsit5())

# using Plots

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)
	# plot(res0.time, res0.xc[1,:], marker = :d)

# here the rate function is constant in between jumps
pbvar = PDMP.PDMPProblem(F!, Rvar, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	resvar = PDMP.solve(pbvar, algo; n_jumps = nj)
	@test resvar.time == res0.time
	# plot!(resvar.time, resvar.xc[1,:], label = "Variable")

# here the rate function is constant in between jumps
pbcmp = PDMP.PDMPProblem(F!, Rcmp, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescmp = PDMP.solve(pbcmp, algo; n_jumps = nj)
	@test rescmp.time == res0.time
	# plot!(rescmp.time, rescmp.xc[1,:], label = "Composite")

# using PrettyTables
# datat = hcat(res0.time, resvar.time, rescst.time)
# datax = hcat(res0.xc', resvar.xc', rescst.xc')
#
# PrettyTables.pretty_table(datat)
# 	PrettyTables.pretty_table(datax)

# test for the end time
tf = 0.6
Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)
	# plot(res0.time, res0.xc[1,:], marker = :d)

pbvar = PDMP.PDMPProblem(F!, Rvar, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	resvar = PDMP.solve(pbvar, algo; n_jumps = nj)
	@test resvar.time == res0.time
	# plot!(rescmp.time, resvar.xc[1,:])

pbcmp = PDMP.PDMPProblem(F!, Rcmp, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescmp = PDMP.solve(pbcmp, algo; n_jumps = nj)
	@test rescmp.time == res0.time
	# plot!(rescmp.time, rescmp.xc[1,:])

# test for allocations
Random.seed!(8)
	res0 = PDMP.solve(problem, algo; n_jumps = 3nj, save_positions = (false, false))
	resvar = PDMP.solve(pbvar, algo; n_jumps = 3nj, save_positions = (false, false))
	rescmp = PDMP.solve(pbcmp, algo; n_jumps = 3nj, save_positions = (false, false))

	Random.seed!(8)
	res0 = @timed PDMP.solve(problem, algo; n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	resvar = @timed PDMP.solve(pbvar, algo; n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	rescmp = @timed PDMP.solve(pbcmp, algo; n_jumps = nj, save_positions = (false, false))
	alloc0 = res0[5].poolalloc
	allocvar = resvar[5].poolalloc
	alloccmp = rescmp[5].poolalloc
	Random.seed!(8)
	res0 = @timed PDMP.solve(problem, algo; n_jumps = 3nj, save_positions = (false, false))
	Random.seed!(8)
	resvar = @timed PDMP.solve(pbvar, algo; n_jumps = 3nj, save_positions = (false, false))
	Random.seed!(8)
	rescmp = @timed PDMP.solve(pbcmp, algo; n_jumps = 3nj, save_positions = (false, false))
	# @test res0[5].poolalloc != alloc0
	# @test resvar[5].poolalloc != allocvar
	# @test rescmp[5].poolalloc != alloccmp

println("--> Allocations with Composite struct")
	Random.seed!(8)
	pbcmp = PDMP.PDMPProblem(F!, Rcmp, nu, xc0, xd0, parms, (ti, 100000.))
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = nj, save_positions = (false, false))
	# plot(rescmp.time, rescmp.xc[1,:], label = "Composite", marker=:d)
	Random.seed!(8)
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = 2nj, save_positions = (false, false))
	# plot!(rescmp.time, rescmp.xc[1,:], label = "Composite", marker=:c)
	Random.seed!(8)
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = 3nj, save_positions = (false, false))
	# plot!(rescmp.time, rescmp.xc[1,:], label = "Composite")


# test with different CompositeStruct made of all variables
tf = 100000.
problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)

	Rcmp = CompositeRate(VariableRate(Rcst!), VariableRate(Rvar!))
	pbcmp = PDMP.PDMPProblem(F!, Rcmp, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescmp = PDMP.solve(pbcmp, algo; n_jumps = nj)
	@test rescmp.time == res0.time

println("--> Allocations with Composite struct made of VariableRates")
	Random.seed!(8)
	pbcmp = PDMP.PDMPProblem(F!, Rcmp, nu, xc0, xd0, parms, (ti, 100000.))
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = nj, save_positions = (false, false))
	Random.seed!(8)
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = 2nj, save_positions = (false, false))
	Random.seed!(8)
	rescmp = @time PDMP.solve(pbcmp, algo; n_jumps = 3nj, save_positions = (false, false))
