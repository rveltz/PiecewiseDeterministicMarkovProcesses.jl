function rejection{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Function,R::Function,DX::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false;algo = :cvode)
  @assert algo in [:cvode,:lsoda]
  # it is faster to pre-allocate arrays and fill it at run time
  n_max += 1 #to hold initial vector
  nsteps = 1
  npoints = 2 # number of points for ODE integration

  # Args
  args = pdmpArgs(xc0,xd0,F,R,DX,nu,parms,tf)
  if verbose println("--> Args saved!") end

  # Set up initial variables
  t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec(xc0)
  xd0 = reshape(xd0,1,length(xd0))
  Xd  = deepcopy(xd0)
  deltaxc = copy(nu[1,:]) #declare this variable

  # arrays for storing history, pre-allocate storage
  t_hist  = Array(Float64, n_max)
  xc_hist = Array(Float64, length(xc0), n_max)
  xd_hist = Array(Int64,   length(xd0), n_max)
  res_ode = Array{Float64,2}

  # initialise arrays
  t_hist[nsteps] = t
  xc_hist[:,nsteps] = copy(xc0)
  xd_hist[:,nsteps] = copy(xd0)

  # Main loop
  termination_status = "finaltime"

  reject = true
  lambda_star = 0.0 # this is the bound for the rejection method
  tp = [0.,0.]
  ppf = R(X0,Xd,t,parms,true)
  while (t < tf) && (nsteps < n_max)
    if verbose println("--> step : ",nsteps," / ",n_max ) end
    reject = true
    while (reject) && (nsteps < n_max)
      tp = [t, min(tf, t - log(rand())/ppf[2]) ] #mettre un lambda_star?
      if algo==:cvode
		  res_ode = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd,t,parms), X0, tp, abstol = 1e-9, reltol = 1e-7)
      elseif algo==:lsoda
          _,res_ode = LSODA.lsoda((t,x,xdot,data)->F(xdot,x,Xd,t,parms), X0, tp, abstol = 1e-9, reltol = 1e-7)
	  end
      X0 = vec(res_ode[end,:])
      t = tp[end]
      ppf = R(X0,Xd,t,parms,true)
      if t == tf
        reject = false
      else
        reject = rand() < (1. - ppf[1] / ppf[2])
      end
      nsteps += 1
      t_hist[nsteps] = t
      xc_hist[:,nsteps] = copy(X0[1:end])
      xd_hist[:,nsteps] = copy(Xd)

    end
    # there is a jump!
    ppf = R(X0,Xd,t,parms,false)
    pf = WeightVec(convert(Array{Float64,1},ppf[1])) #this is to ease sampling

    if (t < tf)
      # make a jump
      ev = sample(pf)
      deltaxd = nu[ev,:]

      # Xd = Xd .+ deltaxd
      Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

      # Xc = Xc .+ deltaxc
      DX(X0,Xd,X0[end],parms,ev)
    end
  end
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
  result = pdmpResult(t_hist[1:nsteps],xc_hist[:,1:nsteps],xd_hist[:,1:nsteps],stats,args)
  return(result)
end

"""
This function performs a pdmp simulation using the rejection method when the flow is known analytically.
It takes the following arguments:

- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **Phi** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"""
function rejection_exact{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},Phi::Function,R::Function,DX::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false, xd_jump::Bool=true)
  # it is faster to pre-allocate arrays and fill it at run time
  n_max += 1 #to hold initial vector
  const nsteps = 1
  const npoints = 2 # number of points for ODE integration
  const njumps = 1

  # Args
  args = pdmpArgs(xc0,xd0,Phi,R,DX,nu,parms,tf)

  # Set up initial variables
  t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec(xc0)
  xd0 = reshape(xd0,1,length(xd0))
  Xd  = deepcopy(xd0)
  deltaxc = copy(nu[1,:]) #declare this variable

  # arrays for storing history, pre-allocate storage
  t_hist  = Array(Float64, n_max)
  xc_hist = Array(Float64, length(xc0), n_max)
  xd_hist = Array(Int64,   length(xd0), n_max)
  res_ode = Array(Float64,2,length(xc0))
  rate_vector = zeros(length(nu[1,:]))


  # initialise arrays
  t_hist[nsteps] = t
  xc_hist[:,nsteps] = copy(xc0)
  xd_hist[:,nsteps] = copy(xd0)

  # Main loop
  termination_status = "finaltime"

  const reject = true
  const nb_rejet::Int = 0
  const lambda_star = 0.0 # this is the bound for the rejection method
  tp = [0.,0.]
  lambda_star = R(rate_vector,X0,Xd,t,parms,true)[2]

  t_hist[njumps] = t
  xc_hist[:,njumps] = copy(X0[1:end])
  xd_hist[:,njumps] = copy(Xd)


  while (t < tf) && (njumps < n_max)
    if verbose println("--> step : ",njumps," / ",n_max ) end
    reject = true
    nsteps = 1
    while (reject) && (nsteps < 10^6) && (t < tf)

      tp = [t, t - log(rand())/lambda_star ] #mettre un lambda_star?
      Phi(res_ode,X0,Xd,tp,parms) # we evolve the flow
      X0 = vec(res_ode[end,:])
      t = tp[end]
      ppf = R(rate_vector,X0,Xd,t,parms,true) #we don't want the full rate vector, just the sum of rates
      @assert(ppf[1] <= ppf[2])
      reject = rand() <  (1. - ppf[1] / ppf[2])
      nsteps += 1
    end
    # keep track of nb of rejections
    nb_rejet += nsteps

    @assert(nsteps <= 10^6,"Error, too many rejections!!")
    njumps += 1
    t_hist[njumps] = t
    xc_hist[:,njumps] = copy(X0[1:end])
    xd_hist[:,njumps] = copy(Xd)
    # there is a jump!
    lambda_star = R(rate_vector,X0,Xd,t,parms,false)[2]
    pf = WeightVec(convert(Array{Float64,1},rate_vector)) #this is to ease sampling

    if (t < tf)
      # make a jump
      ev = sample(pf)

      if xd_jump
        deltaxd = nu[ev,:]
        # Xd = Xd .+ deltaxd
        Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)
      end
      # Xc = Xc .+ deltaxc
      DX(X0,Xd,t,parms,ev)
    end
  end
  println("njumps = ",njumps," / rejections = ", nb_rejet)
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  # if verbose println("--> xc = ",xd_hist[:,1:nsteps]) end
  result = pdmpResult(t_hist[1:njumps],xc_hist[:,1:njumps],xd_hist[:,1:njumps],stats,args)
  return(result)
end
