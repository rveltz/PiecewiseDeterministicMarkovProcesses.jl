push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia/")
cd("/Users/rveltz/work/prog_gd/julia/PDMP.jl/examples")
using PDMP, JSON, GR
GR.inline()
const p  = ( JSON.parsefile("ml.json"))["type II"]
const p1  = ( JSON.parsefile("ml.json"))

function F_ml(xcdot::Vector{Float64}, xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64})
  # vector field used for the continuous variable
  #compute the current, v = xc[1]
    fna::Float64 = p["g_Na"] * (p["v_Na"] - xc[1])
    fk::Float64  = p["g_K"]  * (p["v_K"]  - xc[1])
    fl::Float64  = p["g_L"]  * (p["v_L"]  - xc[1])
    r::Float64   = mod(xd[1],2)
  xcdot[1] = xd[2] / convert(Float64,p["N"]) * fna + xd[4]/convert(Float64,p["M"]) * fk + fl + p["I_app"]
  nothing
end

function R_ml(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64})
  return vec([ p["beta_na"] * exp(4.0 * p["gamma_na"] * xc[1] + 4.0 * p["k_na"]) * xd[1],
                 p["beta_na"] * xd[2],
                 p["beta_k"]  * exp(p["gamma_k"] * xc[1] + p["k_k"]) * xd[3],
                 p["beta_k"] * exp(-p["gamma_k"] * xc[1]  -p["k_k"]) * xd[4],
                 0.0]) # for printing
end

function Delta_ml(xc::Array{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64},ind_reaction::Int64)
  # this function return the jump in the continuous component
  return vec([0.])
end

xc0 = vec([p1["v(0)"]])
xd0 = vec([p["N"],    #Na closed
           0,         #Na opened
           p["M"],    #K closed
           0,         #K opened
           0])       #printing

nu = [[-1 1 0 0 0];[1 -1 0 1 0];[0 0 -1 1 0];[0 0 1 -1 0];[0 0 0 0 1]]
parms = [0.1,0.01]
tf = p1["t_end"]

reload("PDMP")
# Warm up
srand(1234)
result = chv(650,xc0,xd0, F_ml, R_ml,(x,y,t,p,id)->vec([0.]), nu , parms,0.0,0.01,false)
# Real run
srand(Int(floor(time()/10000)))
result = @time chv(2500,xc0,xd0, F_ml, R_ml,(x,y,t,p,id)->vec([0.]), nu , parms,0.0,tf,false)
GR.plot(result.time,[result.xc[1,:][:] 0*result.xd[1,:][:] 1*result.xd[2,:][:]],title = string("#Jumps = ",length(result.time)))



immutable F_type; end
call(::Type{F_type},xcd, xc, xd, t, parms) = F_ml(xcd, xc, xd, t, parms)

immutable R_type; end
call(::Type{R_type},xc, xd, t, parms) = R_ml(xc, xd, t, parms)

immutable DX_type; end
call(::Type{DX_type},xc, xd, t, parms, ind_reaction) = Delta_ml(xc, xd, t, parms, ind_reaction)

dummy_t =  PDMP.chv_optim(2,xc0,xd0,F_type,R_type,DX_type,nu,parms,0.0,tf,false)
dummy_t =  @time PDMP.chv_optim(3000,xc0,xd0,F_type,R_type,DX_type,nu,parms,0.0,tf,false)

@profile PDMP.chv(2500,xc0,xd0, F_ml, R_ml,(x,y,t,p,id)->vec([0.]), nu , parms,0.0,tf,false)
using ProfileView
ProfileView.view()

@code_llvm F_ml(xc0,xc0,xd0,t,parms)

