"""
Finalising function. It is called at the end of each computed jump for the user to alter the saving, plotting... procedure.
"""
function finalize_dummy(rate, xc, xd, p, t)
	nothing
end

"""
Function to pre-allocate arrays contening the result.
"""
function allocate_arrays(ti	,xc0, xd0, n_max; rejection = false, ind_save_c=-1:1, ind_save_d=-1:1)
	if ind_save_c[1] == -1
		ind_save_c = 1:length(xc0)
	end

	if ind_save_d[1] == -1
		ind_save_d = 1:length(xd0)
	end

	if rejection
		X0  = copy(xc0)
		Xc  = copy(xc0)
	else
		# for the CVH method, needs to enlarge the state space
		X0 = copy(xc0); push!(X0,ti)
		Xc = copy(xc0)
	end
	Xd	 = copy(xd0)

	# arrays for storing history, pre-allocate storage
	t_hist  = [ti]
	xc_hist = VectorOfArray([copy(xc0)[ind_save_c]])
	xd_hist = VectorOfArray([copy(xd0)[ind_save_d]])
	res_ode = zeros(2, length(X0))

	return X0, Xc, Xd, t_hist, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c
end


"""
Function copied from Gillespie.jl and StatsBase

This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:
- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.
"""
function pfsample(w::vec, s::Tc, n::Int64) where {Tc, vec <: AbstractVector{Tc}}
	t = rand() * s
	i = 1
	cw = w[1]
	while cw < t && i < n
		i += 1
		@inbounds cw += w[i]
	end
	return i
end

pfsample(rate) = pfsample(rate, sum(rate), length(rate))

"""
This type stores the output composed of:
- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **xc** : containing the simulated states for the continuous variable.
- **xd** : containing the simulated states for the continuous variable.
- **rates** : containing the rates used during the simulation
"""
struct PDMPResult{Tc <: Real, vectype_xc, vectype_xd}
	time::Vector{Tc}
	xc::vectype_xc
	xd::vectype_xd
	rates::Vector{Tc}
	save_positions::Tuple{Bool, Bool}
	njumps::Int64
	nrejected::Int64
end

PDMPResult(time, xchist, xdhist) = PDMPResult(time, xchist, xdhist, eltype(xchist)[], (false, false), length(time), 0)
PDMPResult(time, xchist, xdhist, rates, savepos) = PDMPResult(time, xchist, xdhist, rates, savepos, length(time), 0)

PDMPResult(pb::PDMPProblem, savepos = (false, false)) = PDMPResult(copy(pb.time), copy(pb.Xc), copy(pb.Xd), copy(pb.rate_hist), savepos, pb.simjptimes.njumps, pb.simjptimes.fictitous_jumps)
