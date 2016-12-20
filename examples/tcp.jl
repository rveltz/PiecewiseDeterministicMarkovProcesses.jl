# push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia")
using PDMP

function F_tcp(xcdot, xc, xd, t, parms )
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = xc[1]
  else
    xcdot[1] = -xc[1]
  end
  nothing
end

function R_tcp(xc, xd, t, parms, sum_rate::Bool)
  # rate fonction
  if sum_rate==false
    return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]])
  else
    return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]
  end
end

function Delta_xc_tcp(xc, xd, t, parms, ind_reaction::Int64)
  # inplace implementation of the jump of the continuous variable xc
  return true
end


xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = [0.1]
tf = 200.

srand(1234)
result =  PDMP.chv(2,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false,algo=:cvode)
result =  @time PDMP.chv(200,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false,algo=:cvode)

srand(1234)
result =  PDMP.chv(2,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false,algo=:lsoda)
result =  @time PDMP.chv(200,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false,algo=:lsoda)

println("--> stopping time == tf? (not more) ",maximum(result.time) == tf)
println("#jumps = ", length(result.time))