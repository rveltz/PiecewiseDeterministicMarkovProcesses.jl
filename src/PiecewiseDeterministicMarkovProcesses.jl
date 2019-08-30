module PiecewiseDeterministicMarkovProcesses
	using Random, LinearAlgebra
	using LSODA, Sundials, DifferentialEquations, RecursiveArrayTools, DiffEqBase, SparseArrays
	using ForwardDiff
	import DiffEqBase: solve

	abstract type AbstractPDMPAlgorithm end
	abstract type AbstractCHV <: AbstractPDMPAlgorithm end
	abstract type AbstractCHVIterator <: AbstractCHV end
	abstract type AbstractRejection <: AbstractPDMPAlgorithm end
	abstract type AbstractRejectionExact <: AbstractRejection end
	abstract type AbstractRejectionIterator <: AbstractRejection end

	include("utilsforwarddiff.jl")
	include("problem.jl")
	include("utils.jl")
	include("chv.jl")
	include("chvdiffeq.jl")
	include("rejectiondiffeq.jl")
	include("rejection.jl")
	include("tau-leap.jl")

	export pdmp!,
		ssa,
		chv!,chv,
		rejection!,
		rejection_exact,
		chv_diffeq!,
		rejection_diffeq!,
		pdmpArgs,
		pdmpResult,
		pdmp_data,
		tauleap

	export PDMPProblem, CHV, Rejection, RejectionExact, solve
end # module
