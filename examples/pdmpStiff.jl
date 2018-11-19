using Revise, PDMP, LinearAlgebra, Random, DifferentialEquations

function AnalyticalSample(xc0,xd0,ti,nj::Int64)
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
    return th,xch,xdh
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

xc0 = vec([1.0])
xd0 = vec([0,0])

nu = [[1 0];[0 -1]]
parms = vec([0.0])
tf = 100000.
nj = 100

Random.seed!(8)
    res_a = AnalyticalSample(xc0,xd0,0.,nj)

println("\n\nComparison of solvers")
    for ode in [(:cvode,:cvode),(:lsoda,:lsoda),(Tsit5(),:tsit5),(AutoTsit5(Rosenbrock23()),:tsit5RS23),(Rodas4P(autodiff=false),:rodas4p)]#,(Rodas5(autodiff=false),:rodas5)]:dp5)]
    Random.seed!(8)
    res =  PDMP.pdmp!(xc0, xd0, F!, R!, nu, parms, 0.0, tf, n_jumps = nj, ode = ode[1], verbose = false)
    println("--> norm difference = ", res.time - res_a[1] |> norm, "  - solver = ",ode[2])
end

# using Plots
# plot(res.time,res.xc[1,:],label = "CHV",marker=:d)
#     plot!(res_a[1],res_a[2])
# plot(res.time - res_a[1])
