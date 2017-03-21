# push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia")
using JSON, PDMP

const p0  = convert(Dict{AbstractString,Float64}, JSON.parsefile("../examples/ml.json")["type II"])
const p1  = ( JSON.parsefile("../examples/ml.json"))
include("morris_lecar_variables.jl")
const p_ml = ml(p0)

function F_ml!(xcdot::Vector{Float64}, xc::Vector{Float64},xd::Array{Int64},t::Float64, parms::Vector)
  # vector field used for the continuous variable
  #compute the current, v = xc[1]
  xcdot[1] = xd[2] / p_ml.N * (p_ml.g_Na * (p_ml.v_Na - xc[1])) + xd[4] / p_ml.M * (p_ml.g_K  * (p_ml.v_K  - xc[1]))  + (p_ml.g_L  * (p_ml.v_L  - xc[1])) + p_ml.I_app
  nothing
end

function R_ml(xc::Vector{Float64},xd::Array{Int64},t::Float64, parms::Vector, sum_rate::Bool)
  if sum_rate==false
    return vec([p_ml.beta_na * exp(4.0 * p_ml.gamma_na * xc[1] + 4.0 * p_ml.k_na) * xd[1],
                p_ml.beta_na * xd[2],
                p_ml.beta_k * exp(p_ml.gamma_k * xc[1] + p_ml.k_k) * xd[3],
                p_ml.beta_k * exp(-p_ml.gamma_k * xc[1]  -p_ml.k_k) * xd[4]])
  else
    return (p_ml.beta_na * exp(4.0 * p_ml.gamma_na * xc[1] + 4.0 * p_ml.k_na) * xd[1] +
              p_ml.beta_na * xd[2] +
              p_ml.beta_k * exp( p_ml.gamma_k * xc[1] + p_ml.k_k) * xd[3] +
              p_ml.beta_k * exp(-p_ml.gamma_k * xc[1] - p_ml.k_k) * xd[4])
  end
end


xc0 = vec([p1["v(0)"]])
xd0 = vec([Int(p0["N"]),    #Na closed
           0,               #Na opened
           Int(p0["M"]),    #K closed
           0])              #K opened

nu_ml = [[-1 1 0 0];[1 -1 0 1];[0 0 -1 1];[0 0 1 -1]]
parms = vec([0.])

tf = p1["t_end"]
tf=350.

srand(123)
println("--> chv")
dummy_t =       PDMP.pdmp(6,   xc0,xd0, F_ml!, R_ml, nu_ml, parms,0.0,tf,false,ode=:cvode)
dummy_t = @time PDMP.pdmp(4500,xc0,xd0, F_ml!, R_ml, nu_ml, parms,0.0,tf,false,ode=:cvode)
dummy_t = @time PDMP.pdmp(4500,xc0,xd0, F_ml!, R_ml, nu_ml, parms,0.0,tf,false,ode=:lsoda)

srand(123)
println("--> chv_optim - call")
result =        PDMP.pdmp(2,   xc0,xd0,F_ml!,R_ml,nu_ml,parms,0.0,tf,false, algo=:chv_optim)
result =  @time PDMP.pdmp(4500,xc0,xd0,F_ml!,R_ml,nu_ml,parms,0.0,tf,false, algo=:chv_optim) #cpp = 100ms/2200 jumps
println("#jumps = (dummy / result) ", length(dummy_t.time),", ", length(result.time))

try
  println(norm(dummy_t.time-result.time))
  println("--> xc_f-xc_t = ",norm(dummy_t.xc-result.xc))
  println("--> xd_f-xd_t = ",norm(dummy_t.xd-result.xd))
end