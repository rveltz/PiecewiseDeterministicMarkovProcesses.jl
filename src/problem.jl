using SparseArrays

# Dummy functions to allow not specifying these characteristics
function F_dummy(ẋ, xc, xd, parms, t)
	fill!(ẋ, 0)
	nothing
end

# Dummy flow to be used in rejection algorithm
function Phi_dummy(out, xc, xd, parms, t)
	# vector field used for the continuous variable
	# trivial dynamics
	out[1,:] .= xc
	out[2,:] .= xc
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

	xc::vecc					# current continuous variable
	xd::vecd					# current discrete variable

	xc0::vecc					# initial continuous variable
	xd0::vecd					# initial discrete variable
	ratecache::vecrate			# to hold the rate vector for inplace computations. Also used to initialise rate as this can be an issue for StaticArrays.jl
	parms::Tparms				# container to hold parameters to be passed to F, R, Delta

	function PDMPCaracteristics(F, R, Delta, nu::Tnu, xc0::vecc, xd0::vecd, parms::Tparms) where {Tc, Td, Tparms, Tnu <: AbstractMatrix{Td},
						vecc <: AbstractVector{Tc},
						vecd <: AbstractVector{Td}}
		jump = Jump(nu, Delta)
		rate = dualcache(get_rate_prototype(jump, Tc))
		ratefunction = VariableRate(R)
		return new{typeof(F), typeof(ratefunction), typeof(jump), vecc, vecd, typeof(rate), Tparms}(F, ratefunction, jump, copy(xc0), copy(xd0), copy(xc0), copy(xd0), rate, parms)
	end

	function PDMPCaracteristics(F, R::TR, Delta, nu::Tnu, xc0::vecc, xd0::vecd, parms::Tparms) where {Tc, Td, Tparms, Tnu <: AbstractMatrix{Td},
						vecc <: AbstractVector{Tc},
						vecd <: AbstractVector{Td},
						TR <: AbstractRate}
		jump = Jump(nu, Delta)
		rate = dualcache(get_rate_prototype(jump, Tc))
		return new{typeof(F), typeof(R), typeof(jump), vecc, vecd, typeof(rate), Tparms}(F, R, jump, copy(xc0), copy(xd0), copy(xc0), copy(xd0), rate, parms)
	end
end


function PDMPCaracteristics(F, R, nu::Tnu, xc0::vecc, xd0::vecd, parms::Tparms) where {Tc, Td, Tparms, Tnu <: AbstractMatrix{Td}, vecc <: AbstractVector{Tc}, vecd <: AbstractVector{Td}}
	return PDMPCaracteristics(F, R, Delta_dummy, nu, xc0, xd0, parms)
end

function init!(pb::PDMPCaracteristics)
	pb.xc .= pb.xc0
	pb.xd .= pb.xd0
	init!(pb.R)
end

struct PDMPProblem{Tc, Td, vectype_xc <: AbstractVector{Tc},
						vectype_xd <: AbstractVector{Td},
						Tcar}
	tspan::Vector{Tc}				    			# final simulation time interval, we use an array to be able to mutate it
	simjptimes::PDMPJumpTime{Tc, Td}				# space to save result
	time::Vector{Tc}
	Xc::VectorOfArray{Tc, 2, Array{vectype_xc, 1}}	# continuous variable history
	Xd::VectorOfArray{Td, 2, Array{vectype_xd, 1}}	# discrete variable history
	# variables for debugging
	rate_hist::Vector{Tc}							# to save the rates for debugging purposes
	caract::Tcar									# struct for characteristics of the PDMP
end

pushTime!(pb::PDMPProblem, t) = push!(pb.time, t)
pushXc!(pb::PDMPProblem, xc) = push!(pb.Xc, xc)
pushXd!(pb::PDMPProblem, xd) = push!(pb.Xd, xd)

function init!(pb::PDMPProblem)
	init!(pb.caract)
	pb.simjptimes.tstop_extended = -log(rand())
	pb.simjptimes.lastjumptime = pb.tspan[1]
	pb.simjptimes.njumps = 0
	pb.simjptimes.fictitous_jumps = 0
	resize!(pb.time, 1)
	resize!(pb.rate_hist, 1)
	resize!(pb.Xc.u, 1)
	resize!(pb.Xd.u, 1)
end

# callable struct used in the iterator interface
function (prob::PDMPProblem)(u, t, integrator)
	t == prob.simjptimes.tstop_extended
end

function PDMPProblem(F::TF, R::TR, DX::TD, nu::Tnu,
				xc0::vecc, xd0::vecd, parms::Tp,
				tspan) where {Tc, Td, Tnu <: AbstractMatrix{Td}, Tp, TF ,TR ,TD, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	ti, tf = tspan
	caract = PDMPCaracteristics(F, R, DX, nu, xc0, xd0, parms)
	return PDMPProblem{Tc, Td, vecc, vecd, typeof(caract)}(
			[ti, tf],
			PDMPJumpTime{Tc, Td}(Tc(0), ti, 0, Tc(0), Vector{Tc}([0, 0]), false, 0),
			[ti],
			VectorOfArray([copy(xc0)]), VectorOfArray([copy(xd0)]),
			Tc[],
			caract)
end

function PDMPProblem(F, R, nu::Tnu, xc0::vecc, xd0::vecd, parms,
				tspan; kwargs...) where {Tc, Td, Tnu <: AbstractMatrix{Td}, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	return PDMPProblem(F, R, Delta_dummy, nu, xc0, xd0, parms, tspan; kwargs...)
end

function PDMPProblem(F, R, Delta, reaction_number::Int64, xc0::vecc, xd0::vecd, parms,
				tspan; kwargs...) where {Tc, Td, vecc <: AbstractVector{Tc}, vecd <:  AbstractVector{Td}}
	return PDMPProblem(F, R, Delta, spzeros(Int64, reaction_number, length(xd0)), xc0, xd0, parms, tspan; kwargs...)
end
