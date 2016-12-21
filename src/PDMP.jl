module PDMP

	using Distributions
	using StatsBase
	using DataFrames
	using DataArrays
	using Sundials
	using LSODA
	# using ProgressMeter

	export sample,
		chv,
		rejection,
		rejection_exact,
		chv_optim,
		pdmpArgs,
		pdmpResult,
		pdmp_data

	include("utils.jl")
	include("cvode.jl")
	include("lsoda.jl")
	include("chv.jl")
	include("rejection.jl")

	function sample{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Base.Callable,R::Base.Callable,DX::Base.Callable,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;ode=:cvode,algo=:chv)
		@assert algo in [:chv,:chv_optim,:rejection,:rejection_exact]
		if algo==:chv
			return chv(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode)
		elseif algo==:chv_optim
			return chv_optim(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode)
		elseif algo==:rejection
			rejection(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode)
		elseif algo==:rejection_exact
			return rejection_exact(n_max,xc0,xd0,F,R,DX,nu,parms,ti, tf,verbose,ode=ode)
		end
	end

end # module
