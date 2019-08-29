# Dummy functions to allow not specifying these characteristics
function F_dummy(ẋ, xc, xd, parms, t)
	fill!(ẋ, 0)
	nothing
end

mutable struct PDMPJumpTime{Tc <: Real, Td}
	tstop_extended::Tc
	lastjumptime::Tc
	njumps::Td

	# fields required for the rejection method
	lambda_star::Tc					# bound on the total rate
	ppf::Vector{Tc}
	reject::Bool					# boolean to know whether to reject or not the step
	fictitous_jumps::Td
end

struct PDMPCaracteristics{TF, TR, TJ, vecc, vecd, vecrate, Tparms}
	F::TF						# vector field for ODE between jumps
	R::TR			    		# rate function for jumps
	pdmpjump::TJ
	xc::vecc					# continuous variable
	xd::vecd					# discrete variable
	ratecache::vecrate			# to hold the rate vector for inplace computations. Also used to initialise rate as it can be an issue for StaticArrays
	parms::Tparms				# container to hold parameters to be passed to F,R,Delta

	function PDMPCaracteristics(F, R, Delta, nu::Tnu, xc0::vecc, xd0::vecd, parms::Tparms) where {Tc, Td, Tparms, Tnu <: AbstractMatrix{Td},
						vecc <: AbstractVector{Tc},
						vecd <: AbstractVector{Td}}
		jump = RateJump(nu, Delta)
		rate = zeros(Tc, size(nu, 1))
		return new{typeof(F), typeof(R), typeof(jump), vecc, vecd, typeof(rate), Tparms}(F, R, jump, copy(xc0), copy(xd0), rate, parms)
	end
end


function PDMPCaracteristics(F, R, nu::Tnu, xc0::vecc, xd0::vecd, parms::Tparms) where {Tc, Td, Tparms, Tnu <: AbstractMatrix{Td}, vecc <: AbstractVector{Tc}, vecd <: AbstractVector{Td}}
	return PDMPCaracteristics(F, R, Delta_dummy, nu, xc0, xd0, parms)
end

# """
# The function allows to reset the internal variables xc and xd. Call as follows:
# 	`reset!(pb, xc, xd)`
# """
# function reset!(pb::PDMPCaracteristics, xc, xd)
# 	pb.xc .= xc
# 	pb.xd .= xd
# end


struct PDMPProblem{Tc, Td, vectype_xc <: AbstractVector{Tc},
						vectype_xd <: AbstractVector{Td},
						vectype_rate,
						Tnu <: AbstractArray{Td},
						Tp, TF, TR, Tcar}
	tspan::Tuple{Tc, Tc}			    			# final simulation time interval
	simjptimes::PDMPJumpTime{Tc, Td}				# space to save result
	time::Vector{Float64}
	Xc::VectorOfArray{Tc, 2, Array{vectype_xc, 1}}		# continuous variable history
	Xd::VectorOfArray{Td, 2, Array{vectype_xd, 1}}		# discrete variable history
	# variables for debugging
	save_rate::Bool					# boolean for saving rates
	rate_hist::Vector{Tc}			# to save the rates for debugging purposes

	# structs for characteristics of the PDMP
	caract::Tcar
end

function PDMPProblem(xc0::vectype_xc,
		xd0::vectype_xd,
		rate::vectype_rate,
		F::TF, R::TR, DX::TD,
		nu::Tnu, parms::Tp,
		interval::Tuple{Tc, Tc}, savepre::Bool, verbose::Bool, alg::Talg, saverate = false) where {Tc, Td, vectype_xc <: AbstractVector{Tc}, vectype_xd <: AbstractVector{Td}, vectype_rate, Tnu <: AbstractMatrix{Td}, Tp, TF ,TR ,TD, Talg}
	ti, tf = interval
	ratecache = DiffCache(rate)
	caract = PDMPCaracteristics(F,R,nu,xc0,xd0,parms)
	return PDMPProblem{Tc, Td, vectype_xc, vectype_xd, typeof(ratecache), Tnu, Tp, TF, TR, typeof(caract)}(
			interval,
			PDMPJumpTime{Tc, Td}(-log(rand()), interval[1], 0, Tc(0), Vector{Tc}([0, 0]), false, 0),
			[ti],
			VectorOfArray([copy(xc0)]),
			VectorOfArray([copy(xd0)]),
			saverate, Tc[],
			caract)
end

# callable struct used in the iterator interface
function (prob::PDMPProblem)(u,t,integrator)
	t == prob.simjptimes.tstop_extended
end

# simplified constructors to PDMPProblem
function PDMPProblem(F::TF, R::TR, DX::TD, nu::Tnu,
				xc0::vecc, xd0::vecd, parms::Tp,
				interval;
				saverate = false) where {Tc, Td, Tnu <: AbstractMatrix{Td}, Tp, TF ,TR ,TD, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	ti, tf = interval
	rate = zeros(Tc, size(nu, 1))
	ratecache = copy(rate)#DiffCache(rate)
	caract = PDMPCaracteristics(F,R,DX,nu,xc0,xd0,parms)
	# custom type to collect all parameters in one structure
	return PDMPProblem{Tc, Td, vecc, vecd, typeof(ratecache), Tnu, Tp, TF, TR, typeof(caract)}(
			interval,
			PDMPJumpTime{Tc, Td}(-log(rand()), ti, 0, Tc(0), Vector{Tc}([0, 0]), false, 0),
			[ti],
			VectorOfArray([copy(xc0)]), VectorOfArray([copy(xd0)]),
			saverate, Tc[],
			caract)
end

function PDMPProblem(F, R, nu::Tnu, xc0::vecc, xd0::vecd, parms,
				interval; kwargs...) where {Tc, Td, Tnu <: AbstractMatrix{Td}, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	return PDMPProblem(F, R, Delta_dummy, nu, xc0, xd0, parms, interval; kwargs...)
end

function PDMPProblem(F, R, Delta, reaction_number::Int64, xc0::vecc, xd0::vecd, parms,
				interval; kwargs...) where {Tc, Td, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	return PDMPProblem(F, R, Delta, spzeros(Int64, reaction_number, length(xd0)), xc0, xd0, parms, interval; kwargs...)
end

# """
# The function allows to reset the internal variables xc and xd. Call as follows:
# 	`reset!(pb, xc, xd)`
# """
# function reset!(pb::PDMPProblem, xc, xd)
# 	pb.caract.xc .= xc
# 	pb.caract.xd .= xd
# end

"""
This type stores the output composed of:
- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **xc** : containing the simulated states for the continuous variable.
- **xd** : containing the simulated states for the continuous variable.
- **rates** : containing the rates used during the simulation
"""
struct PDMPResult{Tc <: Real,vectype_xc,vectype_xd}
	time::Vector{Tc}
	xc::vectype_xc
	xd::vectype_xd
	rates::Vector{Tc}
	save_positions::Tuple{Bool, Bool}
end

PDMPResult(time::Vector{Tc},xc::vectype_xc,xd::vectype_xd) where {Tc, vectype_xc, vectype_xd} = PDMPResult{Tc, vectype_xc, vectype_xd}(time, xc, xd, Tc[])


"""
Dummy flow to be used in rejection algorithm
"""
function Phi_dummy(out, xc, xd, parms, t)
	# vector field used for the continuous variable
	# trivial dynamics
	out[1,:] .= xc
	out[2,:] .= xc
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
	t_hist  = zeros(n_max)
	xc_hist = zeros(eltype(xc0), length(ind_save_c), n_max)
	xd_hist = zeros(eltype(xd0), length(ind_save_d), n_max)
	res_ode = zeros(2, length(X0))


	# initialise arrays
	t_hist[1] = ti
	xc_hist[:,1] .= xc0[ind_save_c]
	xd_hist[:,1] .= Xd[ind_save_d]
	return X0, Xc, Xd, t_hist, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c
end

"""
function to save data
"""
# function save_data(nsteps, X0, Xd, xc_hist, xd_hist, ind_save_d, ind_save_c)
# 	@inbounds for ii in eachindex(ind_save_c)
# 		xc_hist[ii,nsteps] = X0[ind_save_c[ii]]
# 	end
# 	@inbounds for ii in eachindex(ind_save_d)
# 		xd_hist[ii,nsteps] = Xd[ind_save_d[ii]]
# 	end
# end

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

# function filter_saveat!(save_at, ti, tf)
# 	filter!(x -> (x >= ti)*(x <= tf), save_at)
# 	if isempty(save_at) || save_at[end] < tf
# 		push!(save_at, tf)
# 	end
# end
