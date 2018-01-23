# push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia")
using PDMP

function F_tcp!(xcdot::Vector, xc::Vector, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = xc[1]
  else
    xcdot[1] = -xc[1]
  end
  nothing
end

function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)
    # rate fonction
    if sum_rate==false
        rate[1] = 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1
        rate[2] = parms[1]
        return 0.
    else
        return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]
    end
end

xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = vec([0.1]) # sampling rate
tf = 200.

println("--> inplace implementation,\n ----> cvode")
# more efficient way, inplace modification
srand(1234)
result2=        PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false,n_jumps = 2)
result2=  @time PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false,n_jumps = 100)
srand(1234)
println(" ----> lsoda")
result2=        PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false,ode=:lsoda,n_jumps = 2)
result2=  @time PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false,ode=:lsoda,n_jumps = 100)

println("--> Case optimised:")
srand(1234)
dummy_t =        PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false, algo=:chv_optim,n_jumps = 2)
dummy_t =  @time PDMP.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu_tcp,parms,0.0,tf,false, algo=:chv_optim,n_jumps = 100)

println("--> stopping time == tf? (not more) ",maximum(result2.time) == tf)
println("#jumps = ", length(result2.time))
