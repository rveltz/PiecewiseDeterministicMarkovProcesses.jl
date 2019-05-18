# Example of neural network
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, SparseArrays, Dates

const N	   = 100
const lambda_ = 0.24
const J	   = 0.98

function f(x)
  return  x.^8
end

function Phi(out::Array{Float64,2}, xc::Vector{Float64},xd::Array{Int64},t::Array{Float64},parms::Vector{Float64})
  # vector field used for the continuous variable
  # for this particular model, the empirical mean is constant between jumps
  xbar::Float64 = sum(xc) / N
  out[1,:] = copy(xc)
  out[2,:] = xbar .+ exp(-lambda_*(t[2]-t[1])) .* (xc .- xbar)
  nothing
end

function R_mf_rejet(rate::Vector{Float64},xc::Vector{Float64},xd::Array{Int64},t::Float64,parms::Vector{Float64}, sum_rate::Bool)
  bound = N * f(1.201)#1.5 works well
  # rate function
  if sum_rate == false
	for i=1:N
	  rate[i] = f(xc[i])
	end
	return 0., bound
  else
	res = 0.
	for i=1:N
	  res+=f(xc[i])
	end
	return res, bound
  end
end

function Delta_xc_mf(xc::Array{Float64,1},xd::Array{Int64},t::Float64,parms::Vector{Float64},ind_reaction::Int64)
  # this function return the jump in the continuous component
  for i=1:N
	xc[i] += J/N
  end
  xc[ind_reaction] = 0.0
  xd[ind_reaction] += 1
  return true
end

xc0 = rand(N)*0.2 .+ 0.5
xd0 = Vector{Int64}(zeros(N))

const nu_neur = SparseArrays.sparse(Array{Int64}(undef,N,N)*0)
parms = [0.1]
tf = 10_050.

println("--> Computing... (",string(Dates.now())[end-7:end],")")
result = @time PiecewiseDeterministicMarkovProcesses.rejection_exact(1,xc0,xd0,Phi,R_mf_rejet,Delta_xc_mf,nu_neur,parms,0.0,tf,false,false)
result = @time PiecewiseDeterministicMarkovProcesses.rejection_exact(40_000,xc0,xd0,Phi,R_mf_rejet,Delta_xc_mf,nu_neur,parms,0.0,tf,false,false,ind_save_d = 1:2,ind_save_c = 1:2)
