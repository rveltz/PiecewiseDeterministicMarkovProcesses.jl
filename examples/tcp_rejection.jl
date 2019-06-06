# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials

function F_tcp!(ẋ, xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
         ẋ[1] = 1
    else
         ẋ[1] = -10xc[1]
    end
    nothing
end

rate_tcp(x) = 1/(1+exp(-x))

function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)
    if sum_rate==false
        rate[1] = rate_tcp(xc[1])
        rate[2] = 0.0
        return rate_tcp(xc[1]), 1.
    else
        return rate_tcp(xc[1]), 1.
    end
end

function AnalyticalSample(xc0,xd0,ti,nj::Int64; verbose = false)
    verbose && printstyled(color=:red,"--> Start analytical method\n")
    xch = [xc0[1]]
    xdh = [xd0[1]]
    th  = [ti]
    t = ti
    xc = xc0[1]
    njumps = 1
    rt = zeros(2)
    lambda_star = R_tcp!(rt,xc0,xd0,ti,Float64[],true)[2]
    rate        = R_tcp!(rt,xc0,xd0,ti,Float64[],true)[1]
    while njumps<nj
        xd  = [xdh[end] ,1]
        S = -log(rand()) / lambda_star

        t += S
        if mod(xd[1],2) == 0
            xc = xc + S
        else
            xc = xc * exp(-10S)
        end
        verbose && println("--> S = $S, t = $t, xc = $xc, xd = $(xd[1]), λ = ",lambda_star)
        #reject?
        lambda_star = R_tcp!(rt,[xc],xd,ti,Float64[],true)[2]
        rate        = R_tcp!(rt,[xc],xd,ti,Float64[],true)[1]
        if rand() > (1 - rate / lambda_star)
            verbose && println("----> Jump!, ratio = ",rate / lambda_star)
            push!(th,t)
            push!(xch,xc)
            push!(xdh,xdh[end] + 1)
            njumps += 1
            # dummy call to rand to emulate sampling pfsample
            dum = -log(rand())
        end
    end
    return th,xch,xdh
end

xc0 = [ 0.0 ]
xd0 = [1, 1]

nu_tcp = [[1 0];[0 -1]]
parms = [0.0]
tf = 100000.
nj = 500
errors = Float64[]

Random.seed!(1234)
    res_a = AnalyticalSample(xc0,xd0,0.0,nj,verbose=false)

println("\n\nComparison of solvers")
    for ode in [(:cvode,"cvode"),(:lsoda,"lsoda"),(CVODE_BDF(),"CVODEBDF"),(CVODE_Adams(),"CVODEAdams"),(Tsit5(),"tsit5"),(Rodas4P(autodiff=false),"rodas4P-noAutoDiff"),(Rodas4P(),"rodas4P-AutoDiff"),(Rosenbrock23(),"RS23"),(AutoTsit5(Rosenbrock23()),"AutoTsit5RS23")]
    Random.seed!(1234)
    res =  PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = nj,   ode = ode[1], verbose = false, save_positions=(false, true), algo=:rejection)
    println("--> norm difference = ", norm(res.time[1:nj] - res_a[1],Inf64), "  - solver = ",ode[2])
    push!(errors,norm(res.xc[1,1:nj] - res_a[2],Inf64))
    end


# using Plots
# plot(res_old.time,res_old.xc',label="Old-rejection")
#     plot!(res.time,[x[1] for x in res.xc],label="Iterator")
#
