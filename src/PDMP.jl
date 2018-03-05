module PDMP

using Sundials
using LSODA



export pdmp!,
	ssa,
	chv!,chv,
	rejection!,
	rejection_exact,
	chv_optim!,
	pdmpArgs,
	pdmpResult,
	pdmp_data,
	tauleap

include("utils.jl")
include("cvode.jl")
include("lsoda.jl")
include("chv.jl")
include("rejection.jl")
include("tau-leap.jl")


"""
This function performs a pdmp simulation using the Change of Variable (CHV, see https://arxiv.org/abs/1504.06873) method or the rejection method.
It takes the following arguments:

- **xc0**: a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0**: a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F!**: an inplace `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system. `F!(xdot,xc,xd,t,parms)` returns `nothing`
- **R!**: an inplace `Function` or a callable type, which itself takes six arguments to represent the rate functions associated to the jumps;rate `Vector` of `Float64` holding the different reaction rates, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise. `R!(rate,xc,xd,t,parms,sum_rate)` returns `Float64,Float64`
- **DX**: a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump. `DX(xc,xd,t,parms,ind_rec)` returns `nothing`
- **nu**: a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : data for the parameters of the system.
- **tf**: the final simulation time (`Float64`)
- **verbose**: a `Bool` for printing verbose.
- **ode**: ode time stepper :cvode or :lsoda.
- **n_jumps**: an `Int64` representing the maximum number of jumps to be computed.
- **ind_save_d**: a range to hold the indices of the discrete variable to be saved
- **ind_save_c**: a range to hold the indices of the continuous variable to be saved
"""
function pdmp!(xc0::AbstractVector{Float64},xd0::AbstractVector{Int64},F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::AbstractArray{Int64},parms,ti::Float64, tf::Float64;verbose::Bool = false,ode=:cvode,algo=:chv, n_jumps = 1000,ind_save_d=-1:1,ind_save_c=-1:1,dt=1.)
	@assert algo in [:chv,:chv_optim,:rejection,:tauleap] "Call $algo() directly please, without passing by pdmp(). Indded, the algo $algo() is specialized for speed and requires a particuliar interface."
	# determine if the Rate function is suited to rejection algorithms in which case we might want to take only
	# the first argument
	# res = R(xc0,xd0,0.,parms,false)
	# R_wrap = R

	# if length(res[1]) == size(nu)[1] # we have a rate function suited to the rejection algorithm, one could also have tested typeof(res[1]) == Vector
	# 	if algo==:rejection
	# 		@assert length(res) == 2 "You need the rate function to provide a global bound on the total rates to call a rejection algorithm, e.g. R(xc,xd,t,parms,sum_of_rates) must return [vector_of_rates,bound] or [sum(vector_of_rates),bound] depending on whether sum_of_rates == true.\n\n\n"
	# 	else
	# 		R_wrap(xc,xd,t,parms,verbose) = R(xc,xd,t,parms,verbose)[1]
	# 	end
	# else # we have a rate function suited to the CVH algorithm
	# 	@assert algo!=:rejection "You need the rate function to provide a global bound on the total rates to call a rejection algorithm, e.g. R(xc,xd,t,parms,sum_of_rates) must return [vector_of_rates,bound] or [sum(vector_of_rates),bound] depending on whether sum_of_rates == true.\n\n\n"
	# end
	if algo==:chv
		return PDMP.chv!(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	elseif algo==:chv_optim
		return PDMP.chv_optim!(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	elseif algo==:rejection
		return PDMP.rejection!(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	elseif algo==:rejection_exact
		return PDMP.rejection_exact(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
	elseif algo==:tauleap
		return PDMP.tauleap(n_jumps,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose=verbose,ode=ode,dt=dt)
	end
end

pdmp!(xc0,xd0,F,R,nu,parms,ti,tf;verbose = false,ode=:cvode,algo=:chv, n_jumps = 1000,ind_save_d=-1:1,ind_save_c=-1:1,dt=1.) = PDMP.pdmp!(xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf,verbose=verbose,ode=ode,algo=algo, n_jumps = n_jumps,ind_save_d=ind_save_d,ind_save_c=ind_save_c)

pdmp!(xc0,xd0,F,R,nu,parms,ti,tf;verbose = false,ode=:cvode,algo=:chv,n_jumps=1000,ind_save_d=-1:1,ind_save_c=-1:1,dt=1.) = PDMP.pdmp!(xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf,verbose=verbose,ode=ode,algo=algo,n_jumps = n_jumps,ind_save_d=ind_save_d,ind_save_c=ind_save_c)

pdmp!(xd0,R,nu,parms,ti,tf;verbose =  false,ode=:cvode,algo=:chv,n_jumps=1000,ind_save_d=-1:1,ind_save_c=-1:1,dt=1.) = PDMP.pdmp!([0.],xd0,F_dummy,R,Delta_dummy,nu,parms,ti, tf,verbose=verbose,ode=ode,algo=algo,n_jumps = n_jumps,ind_save_d=ind_save_d,ind_save_c=ind_save_c)
end # module
