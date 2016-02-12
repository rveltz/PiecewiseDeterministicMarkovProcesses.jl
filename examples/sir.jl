using PDMP, GR
GR.inline()
reload("PDMP")

function R_sir(xc,xd,t::Float64,parms,sum_rate::Bool)
  (S,I,R) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  if sum_rate == false
    return [infection,recovery]
  elseusing PDMP, GR
GR.inline()
reload("PDMP")

function R_sir(xc,xd,t::Float64,parms,sum_rate::Bool)
  (S,I,R) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  if sum_rate == false
    return [infection,recovery]
  else
    return infection+recovery
  end
end

function F_sir(xdot,xc,xd,t::Float64,parms)
  # vector field used for the continuous variable
  xdot[1] = 0.0
  nothing
end

xc0 = vec([0.0])
xd0 = vec([99,10,0])
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/100.0,0.01]
tf = 1000.0

srand(1)
dummy = PDMP.chv(1,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
result = @time PDMP.chv(100000,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
    return infection+recovery
  end
end

function F_sir(xdot,xc,xd,t::Float64,parms)
  # vector field used for the continuous variable
  xdot[1] = 0.0
  nothing
end

xc0 = vec([0.0])
xd0 = vec([99,10,0])
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/100.0,0.01]
tf = 1000.0

srand(1)
dummy = PDMP.chv(1,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
result = @time PDMP.chv(100000,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
GR.plot(result.time,result.xd[1,:][:],"r",result.time, result.xd[2,:][:],"g",result.time, result.xd[3,:][:],"b",title = string("#Jumps = ",length(result.time)))




