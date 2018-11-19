using PDMP, LinearAlgebra, Random, DifferentialEquations


function F_tcp!(ẋ, xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
         ẋ[1] = 1.#xc[1]
    else
         ẋ[1] = -1.#xc[1]
    end
    nothing
end

rate_tcp(x) = 5.0/(1.0 + exp(-x*1 + 5.0)) + 0.1

function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)
    # rate fonction
    if sum_rate==false
        rate[1] = rate_tcp(xc[1])
        rate[2] = parms[1]
        return 0.
    else
        return rate_tcp(xc[1]) + parms[1]
    end
end

xc0 = vec([0.0])
xd0 = vec([0, 1])

nu_tcp = [[1 0];[0 -1]]
parms = vec([0.]) # sampling rate
tf = 100000.

println("\n\n--> inplace implementation,\n ----> cvode")
# more efficient way, inplace modification
Random.seed!(1234)
result2 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = 2,   ode = :cvode)
println(result2.time)
Random.seed!(1234)
result2 =  @time PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, n_jumps = 10000, ode = :cvode)

Random.seed!(1234)
println(" ----> lsoda")
result3 =        PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode=:lsoda, n_jumps = 2)
println(result3.time)
Random.seed!(1234)
result3 =  @time PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode=:lsoda, n_jumps = 10000)

Random.seed!(1234)
println(" ----> DiffEq")
ode = Tsit5()
ode = AutoTsit5(Rosenbrock23())
result4 =  PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode = ode, n_jumps = 2, save_positions = (false,false))
println(result4.time)

Random.seed!(1234)
result4 =  @time PDMP.pdmp!(xc0, xd0, F_tcp!, R_tcp!, nu_tcp, parms, 0.0, tf, ode = ode, n_jumps = 10000,save_positions = (false,false))


# sxc0 = @SVector [x for x in xc0]
# sxd0 = @SVector [x for x in xd0]
# result5 =  @time PDMP.chv_diffeq!(sxc0,sxd0,
#                 F_tcp!,R_tcp!,PDMP.Delta_dummy,
#                 nu_tcp,parms,0.0,tf,false, n_jumps = 10000,ode = Tsit5(), save_positions = (false,true))

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
