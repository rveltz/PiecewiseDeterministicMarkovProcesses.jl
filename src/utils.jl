"""
A type storing the status at the end of a call.
"""
struct pdmpStats
	termination_status::String
	nsteps::Int64
end

"""
A type storing the call.
"""
struct pdmpArgs
	xc::Vector{Float64} # continuous variable
	xd::Vector{Int64}# discrete variable
	F::Any
	R::Any
	Delta::Any
	nu::Matrix{Int64}
	parms::Any
	tf::Float64
end

"""
This type stores the output, and comprises of:

- **time** : a `Vector` of `Float64`, containing the times of simulated events.
- **xc** : a `Matrix` of `Float64`, containing the simulated states for the continuous variable.
- **xd** : a `Matrix` of `Int	64`, containing the simulated states for the continuous variable.
- **stats** : an instance of `PDMPStats`.
- **args** : arguments passed.
"""
struct pdmpResult
	time::Vector{Float64}
	xc::Matrix{Float64}
	xd::Matrix{Int64}
	stats::pdmpStats
	args::pdmpArgs
end

"""
Dummy vector field to be used in gillespie algo
"""
function F_dummy(xcdot::Vector{Float64}, xc::Vector{Float64}, xd::Array{Int64}, t::Float64, parms::Vector{Float64})
	# vector field used for the continuous variable
	xcdot[1] = 0.
	nothing
end

"""
Dummy vector field to be used in gillespie algo
"""
function Delta_dummy(xc, xd::Array{Int64}, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
	return true
end

"""
Dummy flow to be used in gillespie algo
"""
function Phi_dummy(out::Array{Float64,2}, xc::Vector{Float64},xd,t::Array{Float64},parms)
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
		Xc = @view X0[1:end-1]
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