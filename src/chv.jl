function cvode_ode_wrapper(t::Float64, x, xdot, user_data)
  # Reminder: user_data = [F R Xd params]
  x = Sundials.asarray(x)
  xdot = Sundials.asarray(xdot)

  const r = 1.0 / user_data[2](x, user_data[3], t, user_data[4], true)::Float64
  user_data[1](xdot, x, user_data[3], t, user_data[4])
  const ly = length(xdot)
  for i in 1:ly
    xdot[i] = xdot[i] * r
  end
  xdot[end] = r
  return Int32(0)
end

function f_CHV{T}(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64}, xd::Array{Int64,2}, parms::Vector{T})
  # used for the exact method
  const r = 1.0 / R(x,xd,t,parms,true)::Float64
  F(xdot,x,xd,t,parms)
  xdot[end] = 1.0
  scale!(xdot, r)
  return Int32(0)
end


@doc doc"""
This function performs a pdmp simulation using the Change of Variable (CHV) method. It takes the following arguments:
- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"
""" ->
function chv{F,R,DX,T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},::Type{F},::Type{R},::Type{DX},nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false)
  # it is faster to pre-allocate arrays and fill it at run time
  n_max += 1 #to hold initial vector
  nsteps = 1
  npoints = 2 # number of points for ODE integration

  # Args
  args = pdmpArgs(xc0,xd0,F,R,nu,parms,tf)
  if verbose println("--> Args saved!") end

  # Set up initial variables
  t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec([xc0 t])
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
  nsteps += 1

  # Main loop
  termination_status = "finaltime"

  # prgs = Progress(int(tf), 1)

  while (t <= tf) && (nsteps<n_max)
    # update!(prgs, int(t))

    dt = -log(rand())
    if verbose println("--> t = ",t," - dt = ",dt) end

    res_ode = cvode(F,R,Xd,parms, X0, [0.0, dt], abstol = 1e-10, reltol = 1e-8)
    if verbose println("--> Sundials done!") end

    X0 = vec(res_ode[end,:])
    pf = R(X0[1:end-1],Xd,X0[end],parms, false)
    pf = WeightVec(convert(Array{Float64,1},pf)) #this is to ease sampling

    # jump time:
    t = res_ode[end,end]
    # Update event
    ev = sample(pf)
    deltaxd = nu[ev,:]

    # Xd = Xd .+ deltax
    Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

    # Xc = Xc .+ deltax
    DX(X0,Xd,X0[end],parms,ev)

    if verbose println("--> Which reaction? ",ev) end

    # save state
    t_hist[nsteps] = t
    xc_hist[:,nsteps] = copy(X0[1:end-1])
    xd_hist[:,nsteps] = copy(Xd)
    nsteps += 1
  end
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  if verbose println("--> xc = ",xd_hist[:,1:nsteps-1]) end
  result = pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
  return(result)
end


@doc doc"""
This function performs a pdmp simulation using the Change of Variable (CHV) method. It takes the following arguments:
- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"
""" ->
function chv{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Function,R::Function,DX::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false)
  # it is faster to pre-allocate arrays and fill it at run time
  n_max += 1 #to hold initial vector
  nsteps = 1
  npoints = 2 # number of points for ODE integration

  # Args
  args = pdmpArgs(xc0,xd0,F,R,nu,parms,tf)
  if verbose println("--> Args saved!") end

  # Set up initial variables
  t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec([xc0 t])
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
  nsteps += 1

  # Main loop
  termination_status = "finaltime"

  # prgs = Progress(int(tf), 1)

  while (t <= tf) && (nsteps<n_max)
    # update!(prgs, int(t))

    dt = -log(rand())
    if verbose println("--> t = ",t," - dt = ",dt) end

    res_ode = Sundials.cvode((t,x,xdot)->f_CHV(F,R,t,x,xdot,Xd,parms), X0, [0.0, dt], abstol = 1e-10, reltol = 1e-8)
    if verbose println("--> Sundials done!") end
    X0 = vec(res_ode[end,:])
    pf = R(X0[1:end-1],Xd,X0[end],parms, false)
    pf = WeightVec(convert(Array{Float64,1},pf)) #this is to ease sampling

    # jump time:
    t = res_ode[end,end]
    # Update event
    ev = sample(pf)
    deltaxd = nu[ev,:]

    # Xd = Xd .+ deltax
    Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

    # Xc = Xc .+ deltax
    DX(X0,Xd,X0[end],parms,ev)

    if verbose println("--> Which reaction? ",ev) end

    # save state
    t_hist[nsteps] = t
    xc_hist[:,nsteps] = copy(X0[1:end-1])
    xd_hist[:,nsteps] = copy(Xd)
    nsteps += 1
  end
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  if verbose println("--> xc = ",xd_hist[:,1:nsteps-1]) end
  result = pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
  return(result)
end


@doc doc"""
This function performs a pdmp simulation using the Change of Variable (CHV) method. It takes the following arguments:
- **n_max**: an `Int64` representing the maximum number of jumps to be computed.
- **xc0** : a `Vector` of `Float64`, representing the initial states of the continuous variable.
- **xd0** : a `Vector` of `Int64`, representing the initial states of the discrete variable.
- **F** : a `Function` or a callable type, which itself takes five arguments to represent the vector field; xdot a `Vector` of `Float64` representing the vector field associated to the continuous variable, xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time and parms, a `Vector` of `Float64` representing the parameters of the system.
- **R** : a `Function` or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and sum_rate a `Bool` being a flag asking to return a `Float64` if true and a `Vector` otherwise.
- **Delta** : a `Function` or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc `Vector` of `Float64` representing the current state of the continuous variable, xd `Vector` of `Int64` representing the current state of the discrete variable, t a `Float64` representing the current time, parms a `Vector` of `Float64` representing the parameters of the system and ind_rec an `Int64` representing the index of the discrete jump.
- **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
- **parms** : a `Vector` of `Float64` representing the parameters of the system.
- **tf** : the final simulation time (`Float64`)
- **verbose** : a `Bool` for printing verbose.
"
""" ->
function chv_optim{F,R,DX,T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},::Type{F},::Type{R},::Type{DX},nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false)
  # it is faster to pre-allocate arrays and fill it at run time
  n_max += 1 #to hold initial vector
  nsteps = 1
  npoints = 2 # number of points for ODE integration

  # Args
  args = pdmpArgs(xc0,xd0,F,R,nu,parms,tf)
  if verbose println("--> Args saved!") end

  # Set up initial variables
  t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec([xc0 t])
  xd0 = reshape(xd0,1,length(xd0))
  Xd  = deepcopy(xd0)
  deltaxc = copy(nu[1,:]) #declare this variable

  # arrays for storing history, pre-allocate storage
  t_hist  = Array(Float64, n_max)
  xc_hist = Array(Float64, length(xc0), n_max)
  xd_hist = Array(Int64,   length(xd0), n_max)
  res_ode = zeros(2, length(X0))

  # initialise arrays
  t_hist[nsteps] = t
  xc_hist[:,nsteps] = copy(xc0)
  xd_hist[:,nsteps] = copy(xd0)
  nsteps += 1

  # Main loop
  termination_status = "finaltime"

  # save Sundials solver
  mem = cvode_optim(F,R,Xd,parms, X0, [0.0, 1.0], abstol = 1e-10, reltol = 1e-8)

  while (t <= tf) && (nsteps<n_max)

    dt = -log(rand())
    if verbose println("--> t = ",t," - dt = ",dt) end

    evolve(res_ode, mem,F,R,Xd,parms, X0, [0.0, dt])
    if verbose println("--> Sundials done!") end

    X0 = vec(res_ode[end,:])
    pf = R(X0[1:end-1],Xd,X0[end],parms, false)
    pf = WeightVec(convert(Array{Float64,1},pf)) #this is to ease sampling

    # Update time
    t = res_ode[end,end]
    # Update event
    ev = sample(pf)
    deltaxd = nu[ev,:]

    # Xd = Xd .+ deltax
    Base.LinAlg.BLAS.axpy!(1.0, deltaxd, Xd)

    # Xc = Xc .+ deltax
    DX(X0,Xd,X0[end],parms,ev)

    if verbose println("--> Which reaction? ",ev) end

    # save state
    t_hist[nsteps] = t
    # copy cols: faster, cf. performance tips in JuliaLang
    xc_hist[:,nsteps] = X0[1:end-1]
    xd_hist[:,nsteps] = Xd
    nsteps += 1
  end

  Sundials.CVodeFree([mem])
  # collect the data
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  if verbose println("--> xc = ",xd_hist[:,1:nsteps-1]) end
  if verbose println("--> time = ",t_hist[1:nsteps-1]) end
  println("--> chv_optim, #jumps = ",length(t_hist[1:nsteps-1]))
  result = pdmpResult(t_hist[1:nsteps-1],xc_hist[:,1:nsteps-1],xd_hist[:,1:nsteps-1],stats,args)
  return(result)
end
