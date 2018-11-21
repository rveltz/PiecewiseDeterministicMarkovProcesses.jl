"""
A type storing the call.
"""

struct PDMPFunctions{TF,TR,TD}
	F::TF						# vector field for ODE between jumps
	R::TR			    		# rate function for jumps
	Delta::TD		    		# function to implement
end

mutable struct PDMPsimulation{Tc <: Real,Td}
	tstop_extended::Tc
	lastjumptime::Tc
	njumps::Td
end

struct PDMPProblem{Tc,Td,vectype_xc<:AbstractVector{Tc},vectype_xd<:AbstractVector{Td},Tnu<:AbstractArray{Td},Tp,TF,TR,TD}
	xc::vectype_xc					# continuous variable
	xd::vectype_xd					# discrete variable
	pdmpFunc::PDMPFunctions{TF,TR,TD}
	nu::Tnu
	parms::Tp			    		# container to hold parameters to be passed to F,R,Delta
	tf::Tc			    			# final simulation time
	rate::vectype_xc				# to hold the rate vector for inplace computations
	sim::PDMPsimulation{Tc,Td}
	# space to save result
	time::Vector{Float64}
	save_pre_jump::Bool				# save the pre jump?
	Xc::VectorOfArray{Tc,2,Array{vectype_xc,1}}		# continuous variable history
	Xd::VectorOfArray{Td,2,Array{vectype_xd,1}}		# discrete variable history
	verbose::Bool					# print message?
	PDMPProblem{Tc,Td,vectype_xc,vectype_xd,Tnu,Tp,TF,TR,TD}(
			xc0::vectype_xc,xd0::vectype_xd,
			F::TF,R::TR,DX::TD,
			nu::Tnu,parms::Tp,
			ti::Tc,tf::Tc,savepre::Bool,verbose::Bool) where {Tc,Td,vectype_xc<:AbstractVector{Tc},vectype_xd<:AbstractVector{Td},Tnu <: AbstractArray{Td},Tp,TF ,TR ,TD} = new(copy(xc0),copy(xd0),PDMPFunctions(F,R,DX),nu,
			parms,tf,zeros(Tc,size(nu,1)),PDMPsimulation{Tc,Td}(-log(rand()),ti,0),[ti],savepre,VectorOfArray([copy(xc0)]),VectorOfArray([copy(xd0)]),verbose)
end

# callable struct for the CHV method
function (prob::PDMPProblem{Tc,Td,vectype_xc,vectype_xd,Tnu,Tp,TF,TR,TD})(
	xdot,
	x,
	data,t::Tc) where {Tc,Td,vectype_xc,vectype_xd,Tnu<:AbstractArray{Td},Tp,TF,TR,TD}
	tau = x[end]
	sr = prob.pdmpFunc.R(prob.rate,x,prob.xd,tau,prob.parms,true)[1]
	prob.pdmpFunc.F(xdot,x,prob.xd,tau,prob.parms)
	xdot[end] = 1.0
	@inbounds for i in eachindex(xdot)
		xdot[i] = xdot[i] / sr
	end
	nothing
end

"""
This type stores the output, and comprises of:
- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **xc** : a `Matrix` of `Float64`, containing the simulated states for the continuous variable.
- **xd** : a `Matrix` of `Int	64`, containing the simulated states for the continuous variable.
- **stats** : an instance of `PDMPStats`.
- **args** : arguments passed.
"""
struct PDMPResult
	time::Vector{Float64}
	xc::Matrix{Float64}
	xd::Matrix{Int64}
end

"""
Dummy vector field to be used in gillespie algo
"""
function F_dummy(ẋ, xc, xd, t, parms::Ty) where Ty
	# vector field used for the continuous variable
	ẋ[1] = 0.
	nothing
end

"""
Dummy vector field to be used in gillespie algo
"""
function Delta_dummy(xc, xd, t, parms::Ty, ind_reaction::Int64) where Ty
	return true
end


"""
Dummy flow to be used in gillespie algo
"""
function Phi_dummy(out::Array{Float64,2}, xc::Vector{Float64},xd,t::Array{Float64},parms::Ty) where Ty
    # vector field used for the continuous variable
    # trivial dynamics
    out[1,:] .= xc
    out[2,:] .= xc
    nothing
end

"""
Function to pre-allocate arrays contening the result.
"""
function allocate_arrays(ti	,xc0,xd0,n_max,rejection = false;ind_save_c=-1:1,ind_save_d=-1:1)
	if ind_save_c[1] == -1
		ind_save_c = 1:length(xc0)
	end

	if ind_save_d[1] == -1
		ind_save_d = 1:length(xd0)
	end

	if rejection
		X0  = copy(xc0)
		Xc  = copy(X0)
	else
		# for the CVH method, needs to enlarge the state space
		X0 = copy(xc0); push!(X0,ti)
		Xc = copy(xc0)
	end
	Xd     = copy(xd0)

	# arrays for storing history, pre-allocate storage
	t_hist  = zeros(n_max)
	xc_hist = zeros(length(ind_save_c), n_max)
	xd_hist = zeros(length(ind_save_d), n_max)
	res_ode = zeros(2,length(X0))


	# initialise arrays
	t_hist[1] = ti
	xc_hist[:,1] .= copy(xc0)[ind_save_c]
	xd_hist[:,1] .= copy(Xd)[ind_save_d]
return X0, Xc, Xd, t_hist, xc_hist, xd_hist, res_ode, ind_save_d, ind_save_c
end

"""
function to save data
"""
function save_data(nsteps,X0,Xd,xc_hist,xd_hist,ind_save_d, ind_save_c)
	@inbounds for ii in eachindex(ind_save_c)
		xc_hist[ii,nsteps] = X0[ind_save_c[ii]]
    end
    @inbounds for ii in eachindex(ind_save_d)
		xd_hist[ii,nsteps] = Xd[ind_save_d[ii]]
    end
end

"
Function copied from Gillespie.jl and StatsBase

This function is a substitute for `StatsBase.sample(wv::WeightVec)`, which avoids recomputing the sum and size of the weight vector, as well as a type conversion of the propensity vector. It takes the following arguments:
- **w** : an `Array{Float64,1}`, representing propensity function weights.
- **s** : the sum of `w`.
- **n** : the length of `w`.
"
function pfsample(w::Array{Float64,1},s::Float64,n::Int64)
    t = rand() * s
    i = 1
    cw = w[1]
    while cw < t && i < n
        i += 1
        @inbounds cw += w[i]
    end
    return i
end
