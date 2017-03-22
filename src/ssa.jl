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
