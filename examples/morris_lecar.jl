push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia/")
cd("/Users/rveltz/work/prog_gd/julia/PDMP.jl/examples")
using PDMP, JSON, GR
GR.inline()
const p  = convert(Dict{AbstractString,Float64}, JSON.parsefile("ml.json")["type II"])
const p1  = ( JSON.parsefile("ml.json"))

function F_ml(xcdot::Vector{Float64}, xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector)
  # vector field used for the continuous variable
  #compute the current, v = xc[1]
#   const fna = p["g_Na"] * (p["v_Na"] - xc[1])
#   const fk  = p["g_K"]  * (p["v_K"]  - xc[1])
#   const fl  = p["g_L"]  * (p["v_L"]  - xc[1])

#   xcdot[1] = xd[2] / p["N"] * fna + xd[4] / p["M"] * fk  + fl + p["I_app"]
  xcdot[1] = xd[2] / p["N"] * (p["g_Na"] * (p["v_Na"] - xc[1])) + xd[4] / p["M"] * (p["g_K"]  * (p["v_K"]  - xc[1]))  + (p["g_L"]  * (p["v_L"]  - xc[1])) + p["I_app"]
  nothing
end

function R_ml(xc::Vector{Float64},xd::Array{Int64},t::Float64, parms::Vector, sum_rate::Bool)
  if sum_rate==false
    return vec([p["beta_na"] * exp(4.0 * p["gamma_na"] * xc[1] + 4.0 * p["k_na"]) * xd[1],
                p["beta_na"] * xd[2],
                p["beta_k"] * exp(p["gamma_k"] * xc[1] + p["k_k"]) * xd[3],
                p["beta_k"] * exp(-p["gamma_k"] * xc[1]  -p["k_k"]) * xd[4]])
  else
    return p["beta_na"] * exp(4.0 * p["gamma_na"] * xc[1] + 4.0 * p["k_na"]) * xd[1] +
           p["beta_na"] * xd[2] +
           p["beta_k"] * exp( p["gamma_k"] * xc[1] + p["k_k"]) * xd[3] +
           p["beta_k"] * exp(-p["gamma_k"] * xc[1] - p["k_k"]) * xd[4]
  end
end

function Delta_ml(xc::Array{Float64},xd::Array{Int64},t::Float64,parms::Vector,ind_reaction::Int64)
  # this function return the jump in the continuous component
  return vec([0.])
end

immutable F_type; end
call(::Type{F_type},xcd, xc, xd, t, parms) = F_ml(xcd, xc, xd, t, parms)

immutable R_type; end
call(::Type{R_type},xc, xd, t, parms, sr) = R_ml(xc, xd, t, parms, sr)

immutable DX_type; end
call(::Type{DX_type},xc, xd, t, parms, ind_reaction) = Delta_ml(xc, xd, t, parms, ind_reaction)

xc0 = vec([p1["v(0)"]])
xd0 = vec([Int(p["N"]),    #Na closed
           0,         #Na opened
           Int(p["M"]),    #K closed
           0])         #K opened

nu = [[-1 1 0 0];[1 -1 0 1];[0 0 -1 1];[0 0 1 -1]]
parms = vec([p])
tf::Float64 = p1["t_end"]

reload("PDMP")
println(p)
result = chv(6,xc0,xd0, F_ml, R_ml,(x,y,t,p,id)->vec([0.]), nu , parms,0.0,0.01,false)
srand(123)
result = @time chv(2200,xc0,xd0, F_ml, R_ml,(x,y,t,p,id)->vec([0.]), nu , parms,0.0,tf,false)
dummy_t =  PDMP.chv_optim(2,xc0,xd0,F_type,R_type,DX_type,nu,parms,0.0,tf,false)
srand(123)
dummy_t =  @time PDMP.chv_optim(2200,xc0,xd0,F_type,R_type,DX_type,nu,parms,0.0,tf,false) #cpp= 100ms/2200 jumps
println("#jumps = ", length(dummy_t.time))
println(norm(dummy_t.time-result.time))
println("--> xc_f-xc_t = ",norm(dummy_t.xc-result.xc))
println("--> xd_f-xd_t = ",norm(dummy_t.xd-result.xd))
GR.plot(result.time,result.xc[1,:][:],"y",result.time, 0*result.xd[3,:][:] ,result.time,0*result.xd[1,:][:],title = string("#Jumps = ",length(result.time)))

using ProfileView
Profile.clear()
@profile for i=1:10 PDMP.chv_optim(3000,xc0,xd0,F_type,R_type,DX_type,nu,parms,0.0,tf,false) end
ProfileView.view()

