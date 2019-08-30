# Example of neural network
# using Revise
using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, SparseArrays
const PDMP = PiecewiseDeterministicMarkovProcesses

const N		 = 100
const lambda_= 0.24
const J		 = 0.98

function f(x)
	return	x.^8
end

function Phi(out::Array{Float64,2}, xc, xd, parms, t::Array{Float64})
	# vector field used for the continuous variable
	# for this particular model, the empirical mean is constant between jumps
	xbar::Float64 = sum(xc) / N
	out[1,:] = copy(xc)
	out[2,:] = xbar .+ exp(-lambda_*(t[2]-t[1])) .* (xc .- xbar)
	nothing
end

function R_mf_rejet(rate, xc, xd, parms, t::Float64, issum::Bool)
	bound = N * f(1.201)#1.5 works well
	# rate function
	if issum == false
		for i=1:N
			rate[i] = f(xc[i])
		end
		return -1., bound
	else
		res = 0.
		for i=1:N
			res += f(xc[i])
		end
		return res, bound
	end
end

function Delta_xc_mf(xc, xd, parms, t::Float64, ind_reaction::Int64)
	# this function return the jump in the continuous component
	for i=1:N
	xc[i] += J/N
	end
	xc[ind_reaction] = 0.0
	xd[ind_reaction] += 1
	return true
end

Random.seed!(1234)
xc0 = rand(N)*0.2 .+ 0.5
xd0 = Vector{Int64}(zeros(N))

const nu_neur = SparseArrays.sparse(Array{Int64}(undef,N,N)*0)
parms = [0.1]
tf = 10_050.

problem = PDMP.PDMPProblem(Phi,R_mf_rejet,Delta_xc_mf,nu_neur, xc0, xd0, parms, (0.0, tf))
	result = PDMP.solve(problem, PDMP.RejectionExact(); n_jumps = 10_000, ind_save_d = 1:2, ind_save_c = 1:2, verbose = false)

# result = @time PiecewiseDeterministicMarkovProcesses.rejection_exact(1,xc0,xd0,Phi,R_mf_rejet,Delta_xc_mf,nu_neur,parms,0.0,tf,false,false)
# result = @time PiecewiseDeterministicMarkovProcesses.rejection_exact(40_000,xc0,xd0,Phi,R_mf_rejet,Delta_xc_mf,nu_neur,parms,0.0,tf,false,false,ind_save_d = 1:2,ind_save_c = 1:2)
