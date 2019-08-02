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
xd0 = [0,0]

nu = [[1 0];[0 -1]]
parms = [0.0]
ti = 0.322156
tf = 100000.
nj = 50

errors = Float64[]

Random.seed!(8)
	res_a = AnalyticalSample(xc0,xd0,ti,nj)

# build a PDMP Problem
prob = PDMP.PDMPProblem(F!, R!, nu)

Random.seed!(8)
	sol = @time PDMP.solve(prob, CHV(Tsit5()), xc0, xd0, parms, (ti, tf), false; save_positions = (false, true), n_jumps = nj)
	norm(sol.time - res_a[1], Inf64)

# This should give the timings 0.000511 seconds (411 allocations: 29.859 KiB)
# however, we are far from this, more like 0.000985 seconds (6.47 k allocations: 123.891 KiB)
Random.seed!(1234)
	sol = @time PDMP.solve(prob, CHV(Tsit5()), xc0, xd0, parms, (ti, tf), false; save_positions = (false, false), n_jumps = nj)


##########################################
# The above issue comes from the call step!(integrator) in the while loop of `solve` in chvdiffeq.jl
const algo = CHV(Tsit5())
const probcache = PDMP.solve(prob, algo, xc0, xd0, parms, (0., 1.), false;return_pb = true)

using BenchmarkTools
xdot0 = similar(xc0)
rate0 = zeros(2)

xe = zeros(2)
xed = zeros(2)

# This looks good, it does not allocate
@btime (algo)($xed, $xe, $probcache, 0.)
@btime (probcache)($xed, $xe, $parms, 0.)

# seems type stable:
@code_warntype algo(xed, xe, probcache, 0.)

# the vector field
@btime $(probcache)($xed, $xe, $parms, 0.)

# in the `solve` in chvdiffeq.jl, I define this vector field that I pass to ODEProblem
vf = (xdot, x, data, tt) -> algo(xdot, x, data, tt)
# this gives 39.222 ns (1 allocation: 16 bytes) which I don't understand since the line L85 above does not allocate
@btime vf($xed, $xe, $probcache, 0.)
# this seems type stable
@code_warntype vf(xed, xe, parms, 0.)

# trial with an iterator
X_extended = rand(2)
cb = DiscreteCallback(probcache, integrator -> PDMP.chv_callback(integrator, probcache), save_positions = (false, false))


prob_CHV = ODEProblem(probcache, xe, (0.0, 1e9))
integrator = init(prob_CHV, Tsit5(), tstops = probcache.algocache.tstop_extended, save_everystep = false, reltol = 1e-7, abstol = 1e-9, advance_to_tstop = true, callback = cb)

function f(integ, n)
	for ii=1:n
		step!(integ)
		# @show integra
	end
end

@time f(integrator, 100)

####################################################################################################
# I have an example where the iterator solution does not allocate
function f(dx,x,p,t)
	dx[1] = -0.001 * x[1]
	dx[2] = -0.001 * x[2]
end

X_extended = rand(2)
prob_CHV = ODEProblem(f, X_extended, (0.0, 1e9))
integrator = init(prob_CHV, Tsit5(), tstops = 1:1:1000, save_everystep = false, advance_to_tstop = true)

function flow(integ, n)
	for ii=1:n
		step!(integrator)
	end
end

@time flow(integrator,10)
	@show integrator

@time flow(integrator,10)
	@show integrator

@time flow(integrator,50)
	@show integrator
