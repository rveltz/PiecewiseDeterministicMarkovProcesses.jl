using PDMP

function R_sir(xc,xd,t,parms)
  (S,I,R) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  [infection,recovery]
end

function F_sir(xc,xd,t::Float64,parms::Vector{Float64})
  # vector field used for the continuous variable
  return vec([ 0.])
end

xc0 = vec([0.0])
xd0 = vec([999.0,1.0,0.0])
nu = [[-1.0 1.0 0.0];[0.0 -1.0 1.0]]
parms = [0.1/1000.0,0.01]
tf = 100.0

srand(1)
result = chv(xc0,xd0,F_sir,R_sir,nu,parms,0.0,0.01)
srand(1)
result = @time chv(xc0,xd0,F_sir,R_sir,nu,parms,0.0,tf)
println(result.stats)
