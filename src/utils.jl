@doc doc"""
  A type storing the status at the end of a call.
  """ ->
type pdmpStats
  termination_status::ASCIIString
  nsteps::Int64
end

@doc doc"""
  A type storing the call.
  """ ->
type pdmpArgs
  xc::Vector{Float64} # continuous variable
  xd::Vector{Int64}# discrete variable
  F::Any
  R::Any
  nu::Matrix{Int64}
  parms::Vector{Float64}
  tf::Float64
end

@doc doc"""
  This type stores the output, and comprises of:

      - **time** : a `Vector` of `Float64`, containing the times of simulated events.
      - **data** : a `Matrix` of `Float64`, containing the simulated states.
      - **stats** : an instance of `PDMPStats`.
      - **args** : arguments passed.
  """ ->
type pdmpResult
  time::Vector{Float64}
  data::Matrix{Float64}
  stats::pdmpStats
  args::pdmpArgs
end

@doc doc"""
  This takes a single argument of type `pdmpResult` and returns a `DataFrame`.
  """ ->
function pdmp_data(s::pdmpResult)
	println("--> Entry in pdmp_data_list")
  df = hcat(DataFrame(time=s.time),convert(DataFrame,s.data))
  df
end
