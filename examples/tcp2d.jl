using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations
const PDMP = PiecewiseDeterministicMarkovProcesses
function F_tcp!(ẋ, xc, xd, t, parms)
    if mod(xd[1],2)==0
         ẋ[1] =  1.0
         ẋ[2] = -1.0
    else
         ẋ[1] = -1.0
         ẋ[2] =  1.0
    end
    nothing
end

R(x) = 5.0/(1.0 + exp(-x/1.0 + 5.0))

function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)
    if sum_rate==false
        rate[1] = R(xc[1])
        rate[2] = parms[1]
        return 0.
    else
        return R(xc[1]) + parms[1]
    end
end

xc0 = vec([0.05,0.075])
xd0 = vec([0, 1])

nu_tcp = [[1 0];[0 -1]]
parms = vec([0.1])
tf = 100000.

println("--> inplace implementation,\n ----> cvode")
Random.seed!(1234)
result2 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = 2,   ode = :cvode)
println(result2.time)
Random.seed!(1234)
result2 =  @time PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = 1000, ode = :cvode)

Random.seed!(1234)
println(" ----> lsoda")
result3 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode=:lsoda, n_jumps = 2)
println(result3.time)
Random.seed!(1234)
result3 =  @time PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode=:lsoda, n_jumps = 10000)

Random.seed!(1234)
println(" ----> DiffEq")

result4 =       PDMP.chv_diffeq!(xc0,xd0,
                F_tcp!,R_tcp!,PDMP.Delta_dummy,
                nu_tcp,parms,0.0,tf,false, n_jumps = 2,ode = Tsit5(), save_positions = (false,true))

Random.seed!(1234)
result4 =  @time PDMP.chv_diffeq!(xc0,xd0,
                F_tcp!,R_tcp!,PDMP.Delta_dummy,
                nu_tcp,parms,0.0,tf,false, n_jumps = 10000,ode = Tsit5(), save_positions = (false,false))

Random.seed!(1234)
    result4 =  @time PDMP.chv_diffeq!(xc0,xd0,
                F_tcp!,R_tcp!,PDMP.Delta_dummy,
                nu_tcp,parms,0.0,tf,true, n_jumps = 100,ode = Rodas5(), save_positions = (false,false))

# using StaticArrays
# sxc0 = @MVector [x for x in xc0]
# sxd0 = @MVector [x for x in xd0]
# spb = PDMP.PDMPPb(sxc0,sxd0,
#                 F_tcp!,R_tcp!,PDMP.Delta_dummy,
#                 nu_tcp,parms,0.1,tf,false)
#
# sdxc0 = @MVector rand(3)
# spb(sdxc0,sxc0,[1.],0.0)
# @show sdxc0, sxc0
# result5 =  @time PDMP.chv_diffeq!(spb,0.1,tf,true,ode = RK4())


println("--> stopping time == tf? (not more) ",maximum(result2.time) == tf)
println("#jumps = ", length(result2.time))

# println("--> check for time stop at tf")
# tf = 40.5
# nj = 5
# Random.seed!(1234)
# result2 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = nj,   ode = :cvode)
# println(result2)
#
# Random.seed!(1234)
# result3 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode=:lsoda, n_jumps = nj)
# println(result3)
#
# Random.seed!(1234)
#     result4 =  PDMP.chv_diffeq!(xc0,xd0,F_tcp!,R_tcp!,PDMP.Delta_dummy,
#                     nu_tcp,parms,0.0,tf,false, n_jumps = nj, save_positions = (false,true))
#     println(result4)
