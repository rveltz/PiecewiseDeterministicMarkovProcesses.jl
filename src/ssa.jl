function F_dummy(xcdot::Vector{Float64}, xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  xcdot[1] = 0.
  nothing
end

function Delta_dummy(xc, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
  return true
end