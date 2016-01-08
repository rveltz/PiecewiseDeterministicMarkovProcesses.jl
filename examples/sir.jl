using GR
GR.inline()
push!(LOAD_PATH, "/Users/rveltz/work/prog_gd/julia/")
# import PDMP
reload("PDMP")

function R_sir(xc,xd,t::Float64,parms)
  (S,I,R) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  [infection,recovery]
end

function F_sir(xc,xd,t::Float64,parms)
  # vector field used for the continuous variable
  return vec([ 0.])
end

xc0 = vec([0.0])
xd0 = vec([999,1,0])
nu = [[-1 1 0];[0 -1 1]]
parms = [0.1/1000.0,0.01]
tf = 1000.0

srand(1)
dummy = PDMP.chv(1,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
srand(1)
result = @time PDMP.chv(10000,xc0,xd0,F_sir,R_sir,(x,y,t,p,id)->vec([0.]),nu,parms,0.0,tf,false)
GR.plot(result.time,[result.xd[1,:][:] result.xd[2,:][:] result.xd[3,:][:]],title = string("#Jumps = ",length(result.time)))


