using PDMP

function F_tcp!(ẋ, xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
         ẋ[1] = xc[1]
    else
         ẋ[1] = -xc[1]
    end
    nothing
end

function F_tcp(xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
        return vec([xc[1]])
    else
        return vec([-xc[1]])
    end
end

function R_tcp(xc, xd, t, parms, sum_rate::Bool)
    # rate fonction
    if sum_rate==false
        return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]])
    else
        return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]
    end
end

xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = vec([0.1]) # sampling rate
tf = 200.

println("--> basic implementation")
srand(1234)
result =  PDMP.pdmp(2,        xc0,xd0,F_tcp,R_tcp,nu_tcp,parms,0.0,tf,false)
result =  @time PDMP.pdmp(200,xc0,xd0,F_tcp,R_tcp,nu_tcp,parms,0.0,tf,false)

println("--> inplace implementation")
# more efficient way, inplace modification
srand(1234)
result2=        PDMP.pdmp(2,  xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false)
result2=  @time PDMP.pdmp(200,xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false)
result2=        PDMP.pdmp(2,  xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false,ode=:lsoda)
result2=  @time PDMP.pdmp(200,xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false,ode=:lsoda)

println("--> Case optimised:")
srand(1234)
dummy_t =  PDMP.pdmp(2,xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false, algo=:chv_optim)
dummy_t =  @time PDMP.pdmp(200,xc0,xd0,F_tcp!,R_tcp,nu_tcp,parms,0.0,tf,false, algo=:chv_optim)

println("--> stopping time == tf? (not more) ",maximum(result.time) == tf)
println("#jumps = ", length(result.time))
