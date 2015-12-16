immutable Feqn; end
immutable Reqn; end

function f_CHV(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64},xd::Array{Int64,2}, parms::Vector{Float64})
    # used for the exact method
    r::Float64 = sum(R(x,xd,t,parms));
    y = F(x,xd,t,parms)
    ly=length(y)
    for i in 1:ly
      xdot[i] = y[i]/r
    end
    xdot[ly] = 1.0/r
end

function f_CHV(F::Function,R::Function,t::Float64, x::Vector{Float64}, xdot::Vector{Float64},xd, parms::Vector{Float64})
    # used for the exact method
    r::Float64 = sum(R(x,xd,t,parms))
    y=F(x,xd,t,parms)
    ly=length(y)
    for i in 1:ly
      xdot[i] = y[i]/r
    end
    xdot[ly] = 1.0/r
end

function f_CHV(F::Type{Feqn},R::Type{Reqn},t::Float64, x::Vector{Float64}, xdot::Vector{Float64},xd, parms::Vector{Float64})
    # used for the exact method
    r::Float64 = sum(evaluate(R,x,xd,t,parms))
    y=evaluate(F,x,xd,t,parms)
    ly=length(y)
    for i in 1:ly
      xdot[i] = y[i]/r
    end
    xdot[ly] = 1.0/r
end

function f_CHV{f}(fr::Type{f},t::Float64, x::Vector{Float64}, xdot::Vector{Float64},xd, parms::Vector{Float64})
    # used for the exact method
    r::Float64 = sum(R(fr,x,xd,t,parms))
    y=F(fr,x,xd,t,parms)
    ly=length(y)
    for i in 1:ly
      xdot[i] = y[i]/r
    end
    xdot[ly] = 1.0/r
end


@doc doc"""
  This function performs a pdmp simulation using the Change of Variable (CHV) method.
  """ ->
function chv(xc0::Vector{Float64},xd0::Vector{Float64},F::Function,R::Function,nu::Matrix{Float64},parms::Vector{Float64},ti::Float64,tf::Float64,verbose::Bool = false)

	nsteps = 1
  npoints = 2 # number of points for ODE integration

	# Args
  args = pdmpArgs(xc0,xd0,F,R,nu,parms,tf)
  if verbose println("--> Args saved!") end

  # set up new ODE system for the simulation

  # Set up initial variables
	t::Float64 = ti
  xc0 = reshape(xc0,1,length(xc0))
  X0  = vec([xc0 t])
  xd0 = reshape(xd0,1,length(xd0))
  Xd  = vec(deepcopy(xd0))

  # arrays for storing history, pre-allocate storage
  xa = deepcopy(vec([X0[1:end-1], vec(Xd)]))
  ta = Array(Float64,0)
  push!(ta,t)
  f_chv = (t,x,xdot)->f_CHV(F,R,t,x,xdot,Xd,parms)
  # Main loop
  termination_status = "finaltime"
  while (t <= tf)
    dt = -log(rand())
    if verbose println("--> t = ",t," - dt = ",dt) end
    res_ode = cvode(f_chv, X0, [0.0, dt], abstol = 1e-10, reltol = 1e-8)
    if verbose println("--> Sundials done!") end
    X0 = vec(res_ode[end,:])
    pf = R(X0[1:end-1],Xd,X0[end],parms)
    pf = WeightVec(convert(Array{Float64,1},pf))
    # Update time
    if sum(pf) == 0.0
      termination_status = "zeroprop"
      break
    end
    # jumping time:
    t = res_ode[end,end]
    push!(ta,t)
    # Update event
    ev = sample(pf)
    deltax = nu[ev,:]
    Base.LinAlg.BLAS.axpy!(1.0, deltax, Xd)
    if verbose println("--> Which reaction? ",ev) end
    xa = vcat(xa,vec([X0[1:end-1], Xd]))
    # update nsteps
    nsteps += 1
  end
  if verbose println("-->Done") end
  stats = pdmpStats(termination_status,nsteps)
  result = pdmpResult(ta,reshape(xa,length(xc0)+length(xd0),nsteps),stats,args)
  return(result)
end
