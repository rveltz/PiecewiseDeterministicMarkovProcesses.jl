module PiecewiseDeterministicMarkovProcesses
	using Random, LinearAlgebra, SparseArrays#, Parameters
	import LSODA, Sundials
	import RecursiveArrayTools as RAT
	import JumpProcesses as JP
	import SciMLBase
	import SciMLBase: solve
	import PreallocationTools: dualcache, get_tmp

	abstract type AbstractPDMPAlgorithm end
	abstract type AbstractCHV <: AbstractPDMPAlgorithm end
	abstract type AbstractCHVIterator <: AbstractCHV end
	abstract type AbstractRejection <: AbstractPDMPAlgorithm end
	abstract type AbstractRejectionExact <: AbstractRejection end
	abstract type AbstractRejectionIterator <: AbstractRejection end

	include("jumps.jl")
	include("rate.jl")
	include("problem.jl")
	include("utils.jl")
	include("chvdiffeq.jl")
	include("utilsforwarddiff.jl")
	include("chv.jl")
	include("rejectiondiffeq.jl")
	include("rejection.jl")
	include("tau-leap.jl")
	include("diffeqwrap.jl")

	export chv_diffeq!,
			rejection_diffeq!,
			PDMPProblem,
			PDMPResult,
			ConstantRate, VariableRate, CompositeRate

	export PDMPProblem, CHV, Rejection, RejectionExact, solve
end # module
