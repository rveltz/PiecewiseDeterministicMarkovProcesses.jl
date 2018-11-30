using PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random

function R_sir!(rate,xc,xd,t::Float64,parms,sum_rate::Bool)
  (S,I,R,~) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  rate_display = 0.01
  if sum_rate == false
      rate[1] = infection
      rate[2] = recovery
      rate[3] = rate_display
      return 0.
  else
    return infection+recovery + rate_display
  end
end

function F_sir!(xdot,xc,xd,t::Float64,parms)
  # vector field used for the continuous variable
  xdot[1] = 0.0
  nothing
end

xc0 = [0.0]
xd0 = [99,10,0,0]
nu = [[-1 1 0 0];[0 -1 1 0];[0 0 0 1]]
parms = [0.1/100.0,0.01]
tf = 150.0



Random.seed!(1234)
dummy = PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_sir!,R_sir!,nu,parms,0.0,tf,n_jumps=1)
result = @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_sir!,R_sir!,nu,parms,0.0,tf,n_jumps=1000)
# this should throw an error:
# result = @time PiecewiseDeterministicMarkovProcesses.pdmp(1000,xc0,xd0,F_sir,R_sir,nu,parms,0.0,tf,algo=:rejection)

Random.seed!(1234)
# automatic determination of algorithm, here CHV for SSA
result = PiecewiseDeterministicMarkovProcesses.pdmp!(xd0,R_sir!,nu,parms,0.0,tf,n_jumps=1)
result = @time PiecewiseDeterministicMarkovProcesses.pdmp!(xd0,R_sir!,nu,parms,0.0,tf,n_jumps=1000)
