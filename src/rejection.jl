function rejection{T}(n_max::Int64,xc0::Vector{Float64},xd0::Array{Int64,1},F::Function,R::Function,DX::Function,nu::Matrix{Int64},parms::Vector{T},ti::Float64, tf::Float64,verbose::Bool = false)
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
    reject = true
    while (reject) && (nsteps < n_max)

      tp = [t, min(tf, t - log(rand())/ppf[2]) ] #mettre un lambda_star?
      res_ode = Sundials.cvode((t,x,xdot)->F(xdot,x,Xd,t,parms), X0, tp, abstol = 1e-9, reltol = 1e-7)
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
