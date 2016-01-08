using GR
GR.inline()
push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia/")
# import PDMP
reload("PDMP")

function F_tcp(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64})
  # vector field used for the continuous variable
  return vec([1.0])
end

function R_tcp(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64})
  # fonction de tau
  if xd[1] == 0
    return vec([1.0,0.0,100.0]) #transition 0->1
  else
    return vec([0.0,1.0, 100.0]) #transition 1->0
  end
end

function Delta_xc_tcp(xc::Array{Float64,1},xd::Array{Int64},t::Float64,parms::Vector{Float64},ind_reaction::Int64)
  # this function return the jump in the continuous component
  if ind_reaction==2
    return vec(-xc)
  end
  return vec(-xc*0)
end

immutable F_type; end
call(::Type{F_type},xc, xd, t, parms) = F_tcp(xc, xd, t, parms)

immutable R_type; end
call(::Type{R_type},xc, xd, t, parms) = R_tcp(xc, xd, t, parms)

immutable DX_type; end
call(::Type{DX_type},xc, xd, t, parms, ind_reaction) = Delta_xc_tcp(xc, xd, t, parms, ind_reaction)

xc0 = vec([0.0])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[-1 0];[0 1]]
parms = [0.1,0.01]
tf = 1000.


reload("PDMP")
println("Case with types:")
dummy_t =  PDMP.chv(2,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv(200000,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)

println("Case with types optimised:")
dummy_t =  PDMP.chv_optim(2,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv_optim(200000,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)

println("Case with functions:")
dummy_f =  PDMP.chv(2,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)
srand(1234)
dummy_f =  @time PDMP.chv(200000,xc0,xd0,F_tcp,R_tcp,Delta_xc_tcp,nu_tcp,parms,0.0,tf,false)

assert(norm(dummy_f.time-dummy_t.time)==0.0)
println("--> xc_f-xc_t = ",norm(dummy_f.xc-dummy_t.xc))
println("--> xd_f-xd_t = ",norm(dummy_f.xd-dummy_t.xd))

println("For simulations:")
srand(1234)
result = @time PDMP.chv(100000,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)

println(size(result.time))
ind = find(result.time.<140)
GR.plot(result.time[ind],[result.xc[1,:][ind]],colors=["b",".w"],title="Single neuron and sojourn time",legends=["T,d"])
type
