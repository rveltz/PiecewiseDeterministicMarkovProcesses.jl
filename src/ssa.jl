function F_dummy(xcdot::Vector{Float64}, xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  xcdot[1] = 0.
  nothing
end

function Delta_dummy(xc::Array{Float64,1}, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
  return true
end

function Phi_dummy(out::Array{Float64,2}, xc::Vector{Float64},xd::Array{Int64},t::Array{Float64},parms::Vector{Float64})
  # vector field used for the continuous variable
  # trivial dynamics
  out[1,:] .= xc
  out[2,:] .= xc
  nothing
end

"""
This function performs a simulation of a time-dependent jump process.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **R!** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. The returned vector has components. If sum_rate is `False`, one must return rate_vector, bound_ where bound_ is a bound on the total rate vector. In the case sum_rate is `True`, one must return total_rate,bound_ where total_rate is a `Float64` that is the sum of the rates. In any case, the function must return a couple (total_rates, bound) where bound is a bound for the total rate.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **ti** : the initial simulation time (`Float64`)
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
- **algo** : `[:chv,:rejection]` for selecting the algorithm
"""
function ssa{T}(n_max::Int64,xd0::Array{Int64,1},R!::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode = :lsoda,algo=:chv)
  @assert algo in [:chv,:chv_optim,:rejection,:rejection_exact]
  if algo==:chv
    return ssa_chv(n_max,xd0,R!,nu,parms,ti, tf,verbose,ode=ode)
  elseif algo==:rejection
    return ssa_rejection(n_max,xd0,R!,nu,parms,ti, tf,verbose,ode=ode)
  end
end


function ssa_chv{T}(n_max::Int64,xd0::Array{Int64,1},R!::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode = :lsoda)
  @assert ode in [:cvode,:lsoda]
  
  function rate_ssa(xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64}, sum_rate::Bool)
	rate = zeros(Float64,length(nu[:,1]))
	sr,bd = R!(rate,xd,t,parms,sum_rate)
	if sum_rate == false
	    return rate
	  else
	    return sr
	end
  end
  
  xc0 = [0.]
  return PDMP.chv(n_max,xc0,xd0,F_dummy,rate_ssa,Delta_dummy,nu,parms,ti,tf,verbose,ode=ode)
end

function ssa_rejection{T}(n_max::Int64,xd0::Array{Int64,1},R!::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode = :lsoda)
  @assert ode in [:cvode,:lsoda]
  xc0 = [0.]
  
  function Rate_ssa!(rate::Vector{Float64},xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64}, sum_rate::Bool)
	 return R!(rate,xd,t,parms,sum_rate) 
  end
  
  return PDMP.rejection_exact(n_max,xc0,xd0,Phi_dummy,Rate_ssa!,Delta_dummy,nu,parms,ti,tf,verbose)
end