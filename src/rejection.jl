struct RejectionExact <: AbstractRejectionExact end

function solve(problem::PDMPProblem, algo::Rejection{Tode}; verbose::Bool = false, save_rejected=false, ind_save_d=-1:1, ind_save_c=-1:1, n_jumps = Inf64, reltol = 1e-7, abstol = 1e-9, save_positions = (false,true),) where {Tode <: Symbol}
	verbose && println("#"^30)
	ode = algo.ode
	@assert ode in [:cvode, :lsoda, :adams, :bdf]
	verbose && printstyled(color=:red,"--> Start rejection method\n")

	# initialise the problem. If I call twice this function, it should give the same result...
	init!(problem)

	# define the ODE flow
	if ode == :cvode || ode == :bdf
		Flow = (X0_,Xd_,tp_) -> Sundials.cvode(  (tt,x,xdot)->problem.caract.F(xdot,x,Xd,problem.caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :BDF)
	elseif	ode == :adams
		Flow = (X0_,Xd_,tp_) -> Sundials.cvode(  (tt,x,xdot)->problem.caract.F(xdot,x,Xd,problem.caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol, integrator = :Adams)
	elseif ode == :lsoda
		Flow = (X0_,Xd_,tp_) -> LSODA.lsoda((tt,x,xdot,data)->problem.caract.F(xdot,x,Xd,problem.caract.parms,tt), X0_, tp_, abstol = abstol, reltol = reltol)
	end

	ti, tf = problem.tspan
	# it is faster to pre-allocate arrays and fill it at run time
	n_jumps  += 1 #to hold initial vector
	nsteps  = 1
	npoints = 2 # number of points for ODE integration

	xc0 = problem.caract.xc
	xd0 = problem.caract.xd

	# Set up initial variables
	t = ti
	X0, _, Xd, t_hist, xc_hist, xd_hist, res_ode = allocate_arrays(ti, xc0, xd0, n_jumps; rejection = true)

	deltaxd = copy(problem.caract.pdmpjump.nu[1, :]) # declare this variable
	numpf   = size(problem.caract.pdmpjump.nu,1)	 # number of reactions
	rate	= zeros(numpf)  # vector of rates
	tp = [ti, tf]		   # vector to hold the time interval over which to integrate the flow

	#variables for rejection algorithm
	reject = true
	lambda_star = 0.0 # this is the bound for the rejection method
	ppf = problem.caract.R(rate, X0, Xd, problem.caract.parms, t, true)

	# @assert ppf[2] == R(rate,X0+0.1265987*cumsum(ones(length(X0))),Xd,t+0.124686489,parms,true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"
	# rate *= 0;ppf = R(rate,X0,Xd,t,parms,true)
	# @assert sum(rate) == 0 "You cannot modify the first argument of your rate function when sum_rate = true"

	δt = problem.simjptimes.tstop_extended

	while (t < tf) && (nsteps < n_jumps)
		verbose && println("--> step : ",nsteps," / ", n_jumps)
		reject = true
		while reject && (nsteps < n_jumps)
			tp .= [t, min(tf, t + δt / ppf[2]) ] #mettre un lambda_star?
			res_ode .= Flow(X0, Xd, tp)

			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end, ii]
			end
			verbose && println("----> δt = ", δt, ", t∈", tp, ", dt = ", tp[2]-tp[1], ", xc = ", X0)

			t = tp[end]
			ppf = problem.caract.R(rate, X0, Xd, problem.caract.parms, t, true)
			@assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf"
			if t == tf
				reject = false
			else
				reject = rand() < 1 - ppf[1] / ppf[2]
			end
			δt = -log(rand())
		end

		# there is a jump!
		ppf = problem.caract.R(rate, X0, Xd, problem.caract.parms, t, false)

		if (t < tf)
			verbose && println("----> Jump!, ratio = ", ppf[1] / ppf[2], ", xd = ", Xd)
			# make a jump
			ev = pfsample(rate, sum(rate), numpf)
			# deltaxd .= problem.caract.pdmpjump.nu[ev,:]
			#
			# # Xd = Xd .+ deltaxd
			# LinearAlgebra.BLAS.axpy!(1.0, deltaxd, Xd)
			#
			# # Xc = Xc .+ deltaxc
			# problem.caract.pdmpjump.Delta(X0, Xd, problem.caract.parms, t, ev)

			# we perform the jump
			affect!(problem.caract.pdmpjump, ev, X0, Xd, problem.caract.parms, t)

		end

		nsteps += 1
		t_hist[nsteps] = t
		@inbounds for ii in eachindex(X0)
			xc_hist[ii,nsteps] = X0[ii]
		end
		@inbounds for ii in eachindex(Xd)
			xd_hist[ii,nsteps] = Xd[ii]
		end
	end
	if verbose println("-->Done") end
	if verbose println("--> xd = ",xd_hist[:,1:nsteps]) end
	return PDMPResult(t_hist[1:nsteps], xc_hist[:,1:nsteps], xd_hist[:,1:nsteps], Float64[], save_positions)
end

function solve(problem::PDMPProblem, algo::Talgo; verbose::Bool = false, save_rejected = false, ind_save_d=-1:1, ind_save_c=-1:1, n_jumps = Inf64, xd_jump::Bool=true) where {Talgo <: AbstractRejectionExact}
	ti, tf = problem.tspan
	# it is faster to pre-allocate arrays and fill it at run time
	n_jumps += 1 #to hold initial vector
	nsteps = 1
	npoints = 2 # number of points for ODE integration
	njumps = 1

	# initialise the problem. If I call twice this function, it should give the same result...
	init!(problem)

	xc0 = problem.caract.xc
	xd0 = problem.caract.xd


	# Set up initial variables
	t = ti
	X0, _, Xd, t_hist, xc_hist, xd_hist, res_ode = allocate_arrays(ti, xc0, xd0, n_jumps, rejection = true, ind_save_d = ind_save_d, ind_save_c = ind_save_c )

	deltaxd = copy(problem.caract.pdmpjump.nu[1, :]) # declare this variable
	numpf   = size(problem.caract.pdmpjump.nu,1)	 # number of reactions
	rate_vector = zeros(numpf)#vector of rates
	tp = [0., 1.]

	reject = true
	nb_rejet = 0
	lambda_star = 0.0 # this is the bound for the rejection method
	tp = [0., 0.]
	lambda_star = problem.caract.R(rate_vector, X0, Xd, problem.caract.parms, t, true)[2]

	@assert lambda_star == problem.caract.R(rate_vector, X0, Xd, problem.caract.parms, t+rand(), true)[2] "Your rejection bound must be constant in between jumps, it cannot depend on time!!"

	δt = problem.simjptimes.tstop_extended

	while (t < tf) && (njumps < n_jumps)
		verbose && printstyled(color = :red, "--> step : $njumps, / $n_jumps, #reject = $nsteps / $nb_rejet\n" )
		reject = true
		nsteps = 1
		while (reject) && (nsteps < 10^6) && (t < tf)
			tp = [t, min(tf, t + δt / lambda_star)] 		# mettre un lambda_star?
			problem.caract.F(res_ode, X0, Xd, problem.caract.parms, tp) 	# we evolve the flow inplace

			@inbounds for ii in eachindex(X0)
				X0[ii] = res_ode[end, ii]
			end

			t = tp[end]
			ppf = problem.caract.R(rate_vector, X0, Xd, problem.caract.parms, t, true)
			@assert ppf[1] <= ppf[2] "(Rejection algorithm) Your bound on the total rate is wrong, $ppf"
			if t == tf
				reject = false
			else
				reject = rand() < 1 - ppf[1] / ppf[2]
			end
			δt = -log(rand())
			nsteps += 1
		end
		# keep track of nb of rejections
		nb_rejet += nsteps
		njumps += 1

		# there is a jump!
		ppf = problem.caract.R(rate_vector, X0, Xd, problem.caract.parms, t, false)

		if (t < tf)
			verbose && println("----> Jump!, ratio = ", ppf[1] / ppf[2], ", xd = ", Xd)
			# make a jump
			ev = pfsample(rate_vector, sum(rate_vector), numpf)

			# we perform the jump
			affect!(problem.caract.pdmpjump, ev, X0, Xd, problem.caract.parms, t)

		end

		nsteps += 1
		t_hist[njumps] = t
		xc_hist[:,njumps] = X0[ind_save_c]
		xd_hist[:,njumps] = Xd[ind_save_d]

	end
	println("njumps = ", njumps, " / rejections = ", nb_rejet)
	if verbose println("-->Done") end

	if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
	return PDMPResult(t_hist[1:njumps], xc_hist[:,1:njumps], xd_hist[:,1:njumps], Float64[], (false, true))
end
