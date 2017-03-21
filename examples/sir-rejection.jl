push!(LOAD_PATH,"/Users/rveltz/work/prog_gd/julia")
using PDMP

function R_sir_rej(xc,xd,t,parms,sum_rate::Bool)
  (S,I,R,~) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  const rate_display = 0.000
  if sum_rate == false
    return [infection,recovery,rate_display], rate_display + 1.
  else
    return infection+recovery + rate_display, rate_display + 1.
  end
end

function F_sir_rej(xdot,xc,xd,t,parms)
  # vector field used for the continuous variable
  xdot[1] = 0.0
  nothing
end

xc0 = vec([0.0])
xd0 = vec([99,10,0,0])
nu = [[-1 1 0 0];[0 -1 1 0];[0 0 0 1]]
parms = [0.1/100.0,0.01]
tf = 150.0

reload("PDMP")

srand(1234)
dummy = PDMP.pdmp!(1,xc0,xd0,F_sir_rej,R_sir_rej,nu,parms,0.0,tf,false,algo=:rejection,ode=:cvode)
result = @time PDMP.pdmp!(1000,xc0,xd0,F_sir_rej,R_sir_rej,nu,parms,0.0,tf,false,algo=:rejection,ode=:cvode)
srand(1234)
dummy = PDMP.pdmp!(1,xc0,xd0,F_sir_rej,R_sir_rej,nu,parms,0.0,tf,false,algo=:rejection,ode=:lsoda)
result = @time PDMP.pdmp!(1000,xc0,xd0,F_sir_rej,R_sir_rej,nu,parms,0.0,tf,false,algo=:rejection,ode=:lsoda)

