using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations

function F_tcp!(ẋ, xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
         ẋ[1] = 1.
    else
         ẋ[1] = -xc[1]
    end
    nothing
end

rate_tcp(x) = 1/(1+exp(-x))

function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)
    if sum_rate==false
        rate[1] = rate_tcp(xc[1])
        rate[2] = 0.0
        return 0., 1.
    else
        return rate_tcp(xc[1]), 1.
    end
end

xc0 = [ 1.0 ]
xd0 = [0, 1]

nu_tcp = [[1 0];[0 -1]]
parms = [0.0]
tf = 100000.
nj = 200

Random.seed!(1234)
    res_chv =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = nj,   ode =Tsit5(), save_positions=(false, true))

Random.seed!(1234)
    res_old =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = nj,   ode =:lsoda, save_positions=(false, true), algo=:rejection)
    # println(res_old.time[1:7])

Random.seed!(1234)
    res =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = nj,   ode =Tsit5(), save_positions=(false, true), algo=:rejection,verbose=false)
    # println(res_old.time[1:7])

println("--> Error with new iterator scheme = ",norm(res_old.time[1:nj] - res.time,Inf64))


using Plots
plot(res_old.time,res_old.xc',label="Old-rejection")
    plot!(res.time,[x[1] for x in res.xc],label="Iterator")
