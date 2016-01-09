using GR
GR.inline()
push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia/")
# import PDMP
reload("PDMP")

function F_tcp(xcdot::Vector{Float64}, xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector)
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = 3. * xc[1]
  else
    xcdot[1] = -4. * xc[1]
  end
  nothing
end

function R_tcp(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector, sum_rate::Bool)
  # fonction de tau
  if sum_rate==false
    return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]]/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1)
  else
    return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1
  end
end

function Delta_xc_tcp(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector,ind_reaction::Int64)
  return vec([0.])
end

immutable F_type; end
call(::Type{F_type},xcd, xc, xd, t, parms) = F_tcp(xcd, xc, xd, t, parms)

immutable R_type; end
call(::Type{R_type},xc, xd, t, parms, sr) = R_tcp(xc, xd, t, parms, sr)

immutable DX_type; end
call(::Type{DX_type},xc, xd, t, parms, ind_reaction) = Delta_xc_tcp(xc, xd, t, parms, ind_reaction)

xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcp = [[1 0];[0 -1]]
parms = vec([100.])
tf = 10.

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

println(norm(dummy_f.time-dummy_t.time))
println("--> xc_f-xc_t = ",norm(dummy_f.xc-dummy_t.xc))
println("--> xd_f-xd_t = ",norm(dummy_f.xd-dummy_t.xd))

reload("PDMP")
println("For simulations:")
# srand(1234)
parms[1] = 1000.
xd0 = vec([0, 1001])
result = @time PDMP.chv_optim(1600,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)
println(size(result.time))
ind = find(result.time.<49)
GR.plot(result.time[ind],[result.xc[1,:][ind]],colors=["b"],title = string("#Jumps = ",length(result.time)))

using ProfileView
Profile.clear()
@profile PDMP.chv_optim(12000,xc0,xd0,F_type,R_type,DX_type,nu_tcp,parms,0.0,tf,false)
ProfileView.view()

