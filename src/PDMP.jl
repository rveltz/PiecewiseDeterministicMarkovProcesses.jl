module PDMP

	using Distributions
	using StatsBase
	using DataFrames
	using DataArrays
	using Sundials
	using LSODA


	export pdmp,
		ssa,
		chv!,chv,
		rejection,
		rejection_exact,
		chv_optim!,
		pdmpArgs,
		pdmpResult,
		pdmp_data

	include("utils.jl")
	include("cvode.jl")
	include("lsoda.jl")
	include("chv.jl")
	include("rejection.jl")
	include("ssa.jl")

	function pdmp!{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv)
		@assert algo in [:chv,:chv_optim,:rejection] "Call $algo() directly please, without passing by pdmp(). Indded, the algo $algo() is specialized for speed and requires a particuliar interface."
		# determine if the Rate function is suited to rejection algorithms in which case we might want to take only
		# the first argument
		res = R(xc0,xd0,0.,parms,false)
		R_wrap = R

		if length(res[1]) == size(nu)[1] # we have a rate function suited to the rejection algorithm, one could also have tested typeof(res[1]) == Vector
			if algo==:rejection
				@assert (length(res) == 2) "You need the rate function to provide a global bound on the total rates to call a rejection algorithm, e.g. R(xc,xd,t,parms,sum_of_rates) must return [vector_of_rates,bound] or [sum(vector_of_rates),bound] depending on whether sum_of_rates == true. [algo = $algo]\n\n\n"
			else
				R_wrap(xc,xd,t,parms,verbose) = R(xc,xd,t,parms,verbose)[1]
			end
		else # we have a rate function suited to the CVH algorithm
			@assert (algo!=:rejection) "You need the rate function to provide a global bound on the total rates to call a rejection algorithm, e.g. R(xc,xd,t,parms,sum_of_rates) must return [vector_of_rates,bound] or [sum(vector_of_rates),bound] depending on whether sum_of_rates == true. [algo = $algo]\n\n\n"
		end
		if algo==:chv
			return PDMP.chv!(n_max,xc0,xd0,F,R_wrap,DX,nu,parms,ti, tf,verbose,ode=ode)
		elseif algo==:chv_optim
			return PDMP.chv_optim!(n_max,xc0,xd0,F,R_wrap,DX,nu,parms,ti, tf,verbose,ode=ode)
		elseif algo==:rejection
			PDMP.rejection(n_max,xc0,xd0,F,R_wrap,DX,nu,parms,ti, tf,verbose,ode=ode)
		# elseif algo==:rejection_exact
# 			return PDMP.rejection_exact(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode)
		end
	end

pdmp!{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv) = PDMP.pdmp!(n_max,xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf,verbose,ode=ode,algo=algo)

pdmp!{T}(n_max::Int64,xd0::Array{Int64,1},R::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv) =
PDMP.pdmp!(n_max,[0.],xd0,F_dummy,R,Delta_dummy,nu,parms,ti,tf,verbose,ode=ode,algo=algo)

function pdmp{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv)
	try
		# try to see if we have inplace vector field
		F(xc0,xd0,0.,parms)
	    function F_wrap(xcdot, xc, xd, t, parms )
	  	  xcdot .= F(xc,xd,t,parms)
	  	  nothing
	    end
		return PDMP.pdmp!(n_max,xc0,xd0,F_wrap,R,DX,nu,parms,ti, tf,verbose;ode=ode,algo=algo)
	catch e
		return PDMP.pdmp!(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose;ode=ode,algo=algo)
	end
end

pdmp{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv) = PDMP.pdmp(n_max,xc0,xd0,F,R,Delta_dummy,nu,parms,ti, tf,verbose,ode=ode,algo=algo)

function pdmp{T}(n_max::Int64,xd0::Array{Int64,1},R::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv)
	PDMP.pdmp(n_max,[0.],xd0,F_dummy,R,Delta_dummy,nu,parms,ti,tf,verbose,ode=ode,algo=algo)
end

end # module
