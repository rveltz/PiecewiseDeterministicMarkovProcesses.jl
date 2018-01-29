push!(LOAD_PATH,"/Users/rveltz/work/prog_gd/julia")
using PDMP

function R_sir_rej!(rate,xc,xd,t,parms,sum_rate::Bool)
  (S,I,R,~) = xd
  (beta,mu) = parms
  infection = beta*S*I
  recovery = mu*I
  rate_display = parms[1]
  if sum_rate == false
    rate[1] = infection
    rate[2] = recovery
    rate[3] = rate_display
    return 0., rate_display + 3.5
  else
    return infection+recovery + rate_display, rate_display + 3.5
  end
end

xc0 = vec([0.0])
xd0 = vec([99,10,0,0])
nu = [[-1 1 0 0];[0 -1 1 0];[0 0 0 1]]
parms = [0.1/100.0,0.01]
tf = 150.0

srand(1234)
println("--> rejection algorithm for SSA")
dummy = PDMP.pdmp!(1,xd0,R_sir_rej!,nu,parms,0.0,tf,false,algo=:rejection)
result = @time PDMP.pdmp(1000,xd0,R_sir_rej!,nu,parms,0.0,tf,false,algo=:rejection)

srand(1234)
println("--> CHV algorithm for SSA")
dummy = PDMP.pdmp!(1,xd0,R_sir_rej,nu,parms,0.0,tf,false,algo=:chv)
result_chv = @time PDMP.pdmp(1000,xd0,R_sir_rej,nu,parms,0.0,tf,false,algo=:chv)

# using Plots
# gr()
# plot(result.time,result.xd[1,:])
#   plot!(result.time,result.xd[2,:])
#   plot!(result.time,result.xd[3,:])
#   plot!(result_chv.time,result_chv.xd[1,:],marker=:d,color=:blue)
#   plot!(result_chv.time,result_chv.xd[2,:],marker=:d,color=:red)
#   plot!(result_chv.time,result_chv.xd[3,:],marker=:d,color=:green,line=:step)
