using Revise, Test
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
	const PDMP = PiecewiseDeterministicMarkovProcesses

function F!(ẋ, xc, xd, parms, t)
	if mod(xd[1] ,2)==0
		ẋ[1] = 10xc[1]
	else
		ẋ[1] = -3xc[1]^2
	end
end

R(x) = 10.0 + x

function R!(rate, xc, xd, parms, t, issum::Bool)
	if issum == false
		rate[1] = R(mod(xd[1], 2))
		rate[2] = parms[1]
		return 0., 0.
	else
		r = R(mod(xd[1], 2)) + parms[1]
		return r, r
	end
end

xc0 = [1.0]
	xd0 = [0, 0]

	nu = [1 0;0 -1]
	parms = [.0]
	ti = 0.0
	tf = 100000.
	nj = 3

# test the different way to write rates
rate0 = zeros(2)
Rvar = PDMP.VariableRate(R!)
Rcst = PDMP.ConstantRate(R!)

out1 = R!(rate0, [1.0], [2,0], parms, 0., true)
outv = Rvar(rate0, [1.0], [2,0], parms, 0., true)
@test out1 == outv
outc = Rcst(rate0, [10.0], [2,0], parms, 0., true)
outc = Rcst(rate0, [10.0], [2,0], parms, 0., true)


algo = CHV(Tsit5())

Random.seed!(8)
	problem = PDMP.PDMPProblem(F!, R!, nu, xc0, xd0, parms, (ti, tf))
	res0 =  PDMP.solve(problem, algo; n_jumps = nj)

using Plots
plot(res0.time, res0.xc[1,:], marker = :d)

# here the rate function is constant in between jumps
pbcst = PDMP.PDMPProblem(F!, Rcst, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	rescst =  PDMP.solve(pbcst, algo; n_jumps = nj)
	plot!(rescst.time, rescst.xc[1,:])


# here the rate function is constant in between jumps
pbvar = PDMP.PDMPProblem(F!, Rvar, nu, xc0, xd0, parms, (ti, tf))
	Random.seed!(8)
	resvar =  PDMP.solve(pbvar, algo; n_jumps = nj)
	plot!(resvar.time, resvar.xc[1,:])
	@assert resvar.time == res0.time

using PrettyTables
datat = hcat(res0.time, resvar.time, rescst.time)
datax = hcat(res0.xc', resvar.xc', rescst.xc')

PrettyTables.pretty_table(datat)

PrettyTables.pretty_table(datax)

rate0 = zeros(2)
R!(rate0, [1.0], [2,0], parms, 0., true)
Rvar(rate0, [1.0], [2,0], parms, 0., true)
Rcst(rate0, [10.0], [2,0], parms, 0., true)
