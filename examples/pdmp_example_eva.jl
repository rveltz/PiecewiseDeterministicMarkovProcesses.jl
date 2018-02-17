using PDMP

function F_eva!(xcdot, xc, xd, t::Float64, parms::Vector{Float64})
  # vector field used for the continuous variable
  xcdot[1] = -xc[1]+1
  nothing
end

function R_eva(rate,xc, xd, t::Float64, parms, sum_rate::Bool)
  # rate function
  rate_print = 1.
  if sum_rate == false
    if xd[1] == 0
        rate[1] = 1.0
        rate[2] = 0.0
        rate[3] = rate_print
      return 0. #transition 0->1
    else
        rate[1] = 0.0
        rate[2] = 1.0
        rate[3] = rate_print
      return 0.0 #transition 1->0
    end
  else
    if xd[1] == 0
      return 1.0 + rate_print #transition 0->1
    else
      return 1.0 + rate_print #transition 1->0
    end
  end
end

function Delta_xc_eva(xc, xd, t::Float64, parms::Vector{Float64}, ind_reaction::Int64)
  # this function return the jump in the continuous component
  if ind_reaction==2
    xc[1] = 0.0
  end
  return true
end

xc0 = vec([0.0])
xd0 = vec([0, 1])

const nu_eva = [[1 0];[-1 0];[0 1]]
parms = [0.1,0.01]
tf = 100.

println("--> Case simple chv:")
dummy_t =  PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,n_jumps=1)
srand(1234)
dummy_t =  @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,n_jumps=200000)

println("--> Case chv optimised:")
dummy_t =  PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,algo=:chv_optim,n_jumps=20)
srand(1234)
dummy_t =  @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,algo=:chv_optim,n_jumps=200_000)

println("For simulations (lsoda):")
result = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=1)
srand(1234)
result = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=200000)

println("--> Case tauleap:")
resultt = PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=1,algo=:tauleap)
srand(1234)
resultt = @time PDMP.pdmp!(xc0,xd0,F_eva!,R_eva,Delta_xc_eva,nu_eva,parms,0.0,tf,ode=:lsoda,n_jumps=20000,algo=:tauleap,dt=0.01)


# Plots.plot(resultt.time,resultt.xc[1,:])
#
# Plots.plot(result.time,result.xc[1,:])
