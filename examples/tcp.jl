# push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia")
using PDMP

function F_tcp(xcdot::Vector, xc::Vector, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = xc[1]
  else
    xcdot[1] = -xc[1]
  end
  nothing
end

function R_tcp(xc::Vector, xd::Array, t::Float64, parms::Vector, sum_rate::Bool)
  # rate fonction
  if sum_rate==false
    return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]])
  else
    return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]
  end
end

function Delta_xc_tcp(xc::Array{Float64,1}, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
  # inplace implementation of the jump of the continuous variable xc
  return true
end


xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = vec([0.1]) # sampling rate
tf = 200.

srand(1234)
result =  PDMP.chv(2,        xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)
result =  @time PDMP.chv(200,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)

println("--> stopping time == tf? (not more) ",maximum(result.time) == tf)
println("#jumps = ", length(result.time))