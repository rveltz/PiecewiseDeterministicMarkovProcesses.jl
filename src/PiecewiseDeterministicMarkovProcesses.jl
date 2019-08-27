module PiecewiseDeterministicMarkovProcesses
	using Random, LinearAlgebra
	using LSODA, Sundials, DifferentialEquations, RecursiveArrayTools, DiffEqBase
	using ForwardDiff
	import DiffEqBase: solve

	abstract type AbstractPDMPAlgorithm end
	abstract type AbstractCHV <: AbstractPDMPAlgorithm end
	abstract type AbstractCHVIterator <: AbstractCHV end
	abstract type AbstractRejection <: AbstractPDMPAlgorithm end
	abstract type AbstractRejectionExact <: AbstractRejection end
	abstract type AbstractRejectionIterator <: AbstractRejection end

	include("utilsforwarddiff.jl")
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

	export PDMPProblem, CHV, Rejection, solve

"""
This function performs a pdmp simulation using the Change of Variable (CHV, see https://arxiv.org/abs/1504.06873) method or the rejection method.
	It takes the following arguments:

`pdmp!(xc0, xd0, F!, R!, DX, nu, parms, ti, tf; verbose::Bool = false, ode = :cvode, algo=:chv, n_jumps = 30_000, save_positions = (false, true), saverate = false, save_at = [])`

It takes the arguments:
- **xc0**: a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0**: a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F!**: an inplace `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system. `F!(xdot,xc,xd,t,parms)` returns `nothing`
- **R!**: an inplace `Function` or a callable type, which itself takes six arguments to represent the rate functions associated to the jumps;rate `Vector` of `Float64` holding the different reaction rates, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. `R!(rate,xc,xd,t,parms,sum_rate)` returns `Float64,Float64`
- **DX**: a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous/discrete variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump. `DX(xc,xd,t,parms,ind_rec)` returns `nothing`
- `nu`: a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- `parms` : data for the parameters of the system. It is passed to `F!`, `R!` and `DX`.
- `ti`: the initial simulation time (`Float64`)
- `tf`: the final simulation time (`Float64`)
- `verbose`: a `Bool` for printing verbose.
- `ode`: ode time stepper `:cvode`, `:lsoda` or any solver from `DifferentialEquations.jl`, like `CVODE_BDF()`.
- `n_jumps`: an `Int64` representing the maximum number of jumps to be computed. Can be put to `Inf64`.
- `ind_save_d`: a range to hold the indices of the discrete variable to be saved
- `ind_save_c`: a range to hold the indices of the continuous variable to be saved
- `save_positions = (false,true)` indicates whether to save the pre-jump (resp. post-jump) state
- `save_at`: an ordered list of times at which the user want to save the state. Work in progress: not available for every algorithm.
- `saverate::Bool` indicates whether the total rates at each jump is saved. Useful to study how to tidy the bound on the rates in the rejection algorithm.
- `return_pb::Bool` indicates if you want the solver to return the `struct` used to solve the problem
- `abstol = 1e-9`
- `restol = 1e-7`
"""
	function pdmp!(xc0::vecc, xd0::vecd,
					F, R, DX, nu::AbstractArray{Int64},
					parms, ti::Float64, tf::Float64;
					verbose::Bool = false,
					ode::Union{Symbol, DiffEqBase.AbstractODEAlgorithm} = :cvode, algo=:chv, n_jumps::Int64 = 30_000, ind_save_d=-1:1, ind_save_c=-1:1, dt=1.,
					save_at::vecc = Float64[], save_positions = (false,true), saverate = false, return_pb = false,reltol=1e-7, abstol=1e-9) where {vecc <: AbstractVector{Float64},
											  vecd <: AbstractVector{Int64}}

		@assert algo in [:chv, :rejection, :tauleap] "Call $algo() directly please, without passing by pdmp(). Indeed, the algo $algo() is specialized for speed and requires a particuliar interface."

		# hack to call DiffEq solver
		# if typeof(ode) != Symbol && algo==:chv
		# 	return chv_diffeq!(xc0, xd0, F, R, DX, nu, parms, ti, tf, verbose; ode = ode, save_positions = save_positions, n_jumps = n_jumps, saverate = saverate, reltol = reltol, abstol = abstol)
		# end

		if typeof(ode) != Symbol && algo==:rejection
			return rejection_diffeq!(xc0, xd0, F, R, DX, nu, parms, ti, tf, verbose; ode = ode, save_positions = save_positions, n_jumps = n_jumps, saverate = saverate)
		end

		# old solvers. There are kept here in order to be able to use LSODA()
		if algo==:chv
			# return PiecewiseDeterministicMarkovProcesses.chv!(xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c,n_max = n_jumps)
		elseif algo==:rejection
			return PiecewiseDeterministicMarkovProcesses.rejection!(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
		elseif algo==:rejection_exact
			return PiecewiseDeterministicMarkovProcesses.rejection_exact(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
		elseif algo==:tauleap
			return PiecewiseDeterministicMarkovProcesses.tauleap(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose=verbose,ode=ode,dt=dt)
		end
	end

	# simplified calls when not all options are passed
	pdmp!(xc0,xd0,F,R,nu,parms,ti,tf;kwargs...) = PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf;kwargs...)

	pdmp!(xc0,xd0,F,R,nu,parms,ti,tf;kwargs...) = PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf;kwargs...)

	pdmp!(xd0,R,nu,parms,ti,tf;kwargs...) = PiecewiseDeterministicMarkovProcesses.pdmp!([0.],xd0,F_dummy,R,Delta_dummy,nu,parms,ti, tf;kwargs...)

end # module
