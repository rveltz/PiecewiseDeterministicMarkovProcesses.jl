using PDMP, Plots


function F_tcpf(xcdot::Vector{Float64}, xc::Vector{Float64},xd::Array{Int64},t::Float64, parms::Vector)
  # vector field used for the continuous variable
  if mod(xd[1],2)==0
    xcdot[1] = xc[1]
  else
    xcdot[1] = -xc[1]
  end
  nothing
end

function R_tcpf(xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector, sum_rate::Bool)
  # fonction de tau
  if sum_rate==false
    return vec([5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1, parms[1]])
  else
    return 5.0/(1.0 + exp(-xc[1]/1.0 + 5.0)) + 0.1 + parms[1]
  end
end

function Delta_xc_tcpf(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector,ind_reaction::Int64)
  return true
end

immutable F_type; end
call(::Type{F_type},xcd, xc, xd, t, parms) = F_tcpf(xcd, xc, xd, t, parms)

immutable R_type; end
call(::Type{R_type},xc, xd, t, parms, sr) = R_tcpf(xc, xd, t, parms, sr)

immutable DX_type; end
call(::Type{DX_type},xc, xd, t, parms, ind_reaction) = Delta_xc_tcpf(xc, xd, t, parms, ind_reaction)

xc0 = vec([0.05])
xd0 = vec([0, 1])

const nu_tcpf = [[1 0];[0 -1]]
parms = vec([1.]) # sampling rate
tf = 50.

println("Case with functions:")
dummy_f =  PDMP.chv(2,xc0,xd0,F_tcpf,R_tcpf,Delta_xc_tcpf,nu_tcpf,parms,0.0,tf,false)
srand(1234)
dummy_f =  @time PDMP.chv(2000,xc0,xd0,F_tcpf,R_tcpf,Delta_xc_tcpf,nu_tcpf,parms,0.0,tf,false)

println("Case with types:")
dummy_t =  PDMP.chv(2,xc0,xd0,F_type,R_type,DX_type,nu_tcpf,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv(2000,xc0,xd0,F_type,R_type,DX_type,nu_tcpf,parms,0.0,tf,false)

println("Case with optimised types:")
dummy_t =  PDMP.chv_optim(2,xc0,xd0,F_type,R_type,DX_type,nu_tcpf,parms,0.0,tf,false)
srand(1234)
dummy_t =  @time PDMP.chv_optim(2000,xc0,xd0,F_type,R_type,DX_type,nu_tcpf,parms,0.0,tf,false)

println("#jumps = ", length(dummy_f.time))
println(norm(dummy_f.time-dummy_t.time))
println("--> xc_f-xc_t = ",norm(dummy_f.xc-dummy_t.xc))
println("--> xd_f-xd_t = ",norm(dummy_f.xd-dummy_t.xd))

println("For simulations:")
tf = 10.
parms[1] = 1.0
result = @time PDMP.chv_optim(20000,xc0,xd0,F_type,R_type,DX_type,nu_tcpf,parms,0.0,tf,false)
println("--> stopping time == tf? (not more) ",maximum(result.time) == tf,maximum(result.time)," == ",tf)
println("#jumps = ", length(result.time))

ind = find(result.time.<2249)
Plots.plotlyjs()
Plots.plot(result.time[ind],result.xc[1,:][ind],color=:black)
Plots.plot!(result.time[ind],0*result.xd[1,:][ind],color=:red,title = string("TCP_fast #Jumps = ",length(result.time)))

