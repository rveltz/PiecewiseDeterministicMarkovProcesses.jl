# using PDMP, Plots
# push!(LOAD_PATH,"/Users/rveltz/work/prog_gd/julia")
# workspace()
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

immutable F_type_eva; end
@compat (::Type{F_type_eva})(xcd, xc, xd, t, parms) = F_eva(xcd, xc, xd, t, parms)

immutable R_type_eva; end
@compat (::Type{R_type_eva})(xc, xd, t, parms, sr) = R_eva(xc, xd, t, parms, sr)

immutable DX_type_eva; end
@compat (::Type{DX_type_eva})(xc, xd, t, parms, ind_reaction) = Delta_xc_eva(xc, xd, t, parms, ind_reaction)

xc0 = vec([0.0])
xd0 = vec([0, 1])

const nu_eva = [[1 0];[-1 0];[0 1]]
parms = [0.1,0.01]
tf = 100.

println("--> Case with types:")
dummy_t =  PDMP.chv(2,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv(200000,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)

println("--> Case optimised:")
dummy_t =  PDMP.chv_optim(20,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv_optim(200000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)

println("--> Case with types optimised:")
dummy_t =  PDMP.chv_optim(20,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv_optim(200000,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)
dummy_t2 =  @time PDMP.chv_optim(200000,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)

println("--> Case with functions:")
dummy_f =  PDMP.chv(2,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)
srand(1234)
dummy_f =  @time PDMP.chv(200000,xc0,xd0,F_eva,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,false)

println("--> #jumps = ", length(dummy_f.time))
println("--> xc_f-xc_t = ",norm(dummy_f.xc-dummy_t.xc))
println("--> xd_f-xd_t = ",norm(dummy_f.xd-dummy_t.xd))

println("For simulations:")
srand(1234)
result = @time PDMP.chv(10000,xc0,xd0,F_type_eva,R_type_eva,DX_type_eva,nu_eva,parms,0.0,tf,false)

# println(size(result.time))
# ind = find(result.time.<134)
# Plots.plotlyjs()
# Plots.plot(result.time[ind],result.xc[1,ind],title = string("#Jumps = ",length(dummy_f.time)))
