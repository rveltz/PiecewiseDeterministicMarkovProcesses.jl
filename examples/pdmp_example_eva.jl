# push!(LOAD_PATH,"/Users/rveltz/work/prog_gd/julia")
using PDMP
using Compat

function F_eva(xcdot::Vector{Float64}, xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  xcdot[1] = -xc[1]+1.
  nothing
end

function R_eva(xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, sum_rate::Bool)
  # rate function
  rate_print = 1.
  if sum_rate == false
    if xd[1] == 0
      return vec([1.0,0.0,rate_print]) #transition 0->1
    else
      return vec([0.0,1.0, rate_print]) #transition 1->0
    end
  else
    if xd[1] == 0
      return 1.0 + rate_print #transition 0->1
    else
      return 1.0 + rate_print #transition 1->0
    end
  end
end

function Delta_xc_eva(xc::Array{Float64,1}, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
  # this function return the jump in the continuous component
  if ind_reaction==2
    xc[1] = 0.0
  end
  return true
end

xc0 = vec([0.0])
xd0 = vec([0, 1])

const nu_eva = [[1 0];[-1 0];[0 1]]
parms = [0.1,0.01]
tf = 100.

println("--> Case simple chv:")
dummy_t =  PDMP.sample(2,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.sample(200000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)

println("--> Case chv optimised:")
dummy_t =  PDMP.sample(20,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,algo=:chv_optim)
srand(1234)
dummy_t =  @time PDMP.sample(200000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,algo=:chv_optim)

println("--> Case with types optimised:")
dummy_t =  PDMP.sample(20,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,algo=:chv_optim)
srand(1234)
dummy_t =  @time PDMP.sample(200000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,algo=:chv_optim)

println("--> #jumps = ", length(dummy_t.time))

println("For simulations (lsoda):")
result = PDMP.sample(2,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,ode=:lsoda)
srand(1234)
result = @time PDMP.sample(10000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false,ode=:lsoda)