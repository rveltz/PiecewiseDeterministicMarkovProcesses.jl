using LSODA
using PDMP, Plots, LinearAlgebra, Random

function n∞(v,v₃,v₄)
  ξ=(v-v₃)/v₄
  (1+tanh(ξ))/2
end

function τ(v,v₃,v₄,ϕ)
  ξ=(v-v₃)/v₄
  1/(ϕ*cosh(ξ/2))
end;

function m∞(v,v₁,v₂)
  (1+tanh((v-v₁)/v₂))/2
end

function α(v,v₃,v₄,ϕ)
    ξ=(v-v₃)/v₄
    ϕ*cosh(ξ/2)/(1+exp(-2*ξ))
end

function β(v,v₃,v₄,ϕ)
    ξ=(v-v₃)/v₄
    ϕ*cosh(ξ/2)/(1+exp(2*ξ))
end;

function f_ml_chv!(tt,x,xdot,data)
    (I,C,gL,gCa,gK,vL,vCa,vK,v₁,v₂,v₃,v₄,ϕ,Ntot,N) = data
    (v,t) = x
    tr = α(v,v₃,v₄,ϕ)*(Ntot-N)+β(v,v₃,v₄,ϕ)*N
    xdot[1] = (I-gL*(v-vL)-gCa*m∞(v,v₁,v₂)*(v-vCa)-gK*(N/Ntot)*(v-vK))/C/tr
    xdot[2] = 1.0/tr
    nothing
end;

function ml_chv(x0,parms,tf;n_jumps=500)
    (v,N) = x0
    (I,C,gL,gCa,gK,vL,vCa,vK,v₁,v₂,v₃,v₄,ϕ,Ntot) = parms
    eparms = [parms;N]
    t=0.0
    ta = Vector{Float64}()
    va = Vector{Float64}()
    Na = Vector{Float64}()
    push!(ta,t)
    push!(va,v)
    push!(Na,N)
    Flow(v_,t_,s_,eparms_) = LSODA.lsoda((tt,x,xdot,data)->f_ml_chv!(tt,x,xdot,eparms_),
        [v_;t_],
        [0.0,s_],
        abstol=1e-9,#Vector([1.e-10,1.e-6]),
        reltol=1e-7,
        nbsteps=10000)
    n = 1 #number of jumps, fairer to compare with PDMP
    while t<tf && n<n_jumps
        s = -log(rand())
        # res = LSODA.lsoda((tt,x,xdot,data)->f_ml_chv!(tt,x,xdot,eparms),
        #     [v;t],
        #     [0.0,s],
        #     abstol=1e-9,#Vector([1.e-10,1.e-6]),
        #     reltol=1e-7,
        #     nbsteps=10000)
        res = Flow(v,t,s,eparms)
        v = res[end,1]
        t = res[end,end]
        # Update N
        opn = α(v,v₃,v₄,ϕ)*(Ntot-N)
        cls = β(v,v₃,v₄,ϕ)*N
        p=opn/(opn+cls)
        if rand()<p
            N=N+1
            eparms[end]=N
        else
            N=N-1
            eparms[end]=N
        end
        push!(ta,t)
        push!(va,v)
        push!(Na,N)
        n += 1
    end
    return(ta,va,Na)
end

parms_chv=[100.0,20.0,2.0,4.4,8.0,-60.0,120.0,-84.0,-1.2,18.0,2.0,30,0.04,40]
x0_chv=[-50.0;20.0]
tf_chv=100000.;

srand(123)
sol_chv=ml_chv(x0_chv,parms_chv,tf_chv,n_jumps=660)
plot(sol_chv[1],sol_chv[2])


println("="^70)
srand(123)
@time begin
    for i in 1:100
        out=ml_chv(x0_chv,parms_chv,tf_chv,n_jumps=660)
    end
end

################################################################################
################################################################################
################################################################################

function f_ml_pdmp!(xcdot, xc, xd, t, parms)
  (I,C,gL,gCa,gK,vL,vCa,vK,v₁,v₂,v₃,v₄,ϕ,Ntot) = parms
  (v,) = xc
  (N,) = xd
  xcdot[1] = (I-gL*(v-vL)-gCa*m∞(v,v₁,v₂)*(v-vCa)-gK*(N/Ntot)*(v-vK))/C
  nothing
end

function r_ml_pdmp!(rate, xc, xd, t, parms, sum_rate)
  (I,C,gL,gCa,gK,vL,vCa,vK,v₁,v₂,v₃,v₄,ϕ,Ntot) = parms
  (v,) = xc
  (N,) = xd
  if sum_rate==false
      rate[1] = α(v,v₃,v₄,ϕ)*(Ntot-N)
      rate[2] = β(v,v₃,v₄,ϕ)*N
      return 0.
  else
      return α(v,v₃,v₄,ϕ)*(Ntot-N)+β(v,v₃,v₄,ϕ)*N
  end
end

rate_ =zeros(2)
xc0 = [-50.0]
xd0 = [20]
xd0 |> typeof |> println
nu_ml = reshape([[1];[-1]],2,1)
parms_chv_pdmp = [100.0,20.0,2.0,4.4,8.0,-60.0,120.0,-84.0,-1.2,18.0,2.0,30,0.04,40]
tf_pdmp = 100000.;

srand(123)
# sol_chv_pdmp=PDMP.pdmp(xc0, xd0, f_ml_pdmp!, r_ml_pdmp!, nu_ml, parms_chv_pdmp, 0.0, tf_pdmp, false, ode=:lsoda, n_jumps = 500);
sol_chv_pdmp=PDMP.chv!(660, xc0, xd0, f_ml_pdmp!, r_ml_pdmp!,PDMP.Delta_dummy, nu_ml, parms_chv_pdmp, 0.0, tf_pdmp, ode=:lsoda)
plot(sol_chv_pdmp.time,sol_chv_pdmp.xc[1,:]-sol_chv[2])

# 1.022143 seconds (8.14 M allocations: 443.023 MiB, 12.90% gc time)
# 1.072882 seconds (8.35 M allocations: 445.167 MiB, 12.10% gc time)
srand(123)
@time begin
    for i in 1:100
        # PDMP.pdmp(xc0, xd0, f_ml_pdmp!, r_ml_pdmp!, nu_ml, parms_chv_pdmp, 0.0, tf_pdmp, false, ode=:lsoda, n_jumps = 500);
        res_pdmp = PDMP.chv!(660, xc0, xd0, f_ml_pdmp!, r_ml_pdmp!,PDMP.Delta_dummy, nu_ml, parms_chv_pdmp, 0.0, tf_pdmp, ode=:lsoda);
    end
end

################################################################################
################################################################################
################################################################################
xd0chv=copy(x0_chv)
eparms = [parms_chv;20]
@time(f_ml_chv!(1.0,x0_chv,xd0chv,eparms))
    println("--> CHV result xdot = $xd0chv for x = $x0_chv")
    xc1 = [xc0;1.]
    xdotpdmp=copy(xc1)
    @time (PDMP.f_CHV!(f_ml_pdmp!, r_ml_pdmp!,0.,xc1,xdotpdmp,xd0,parms_chv_pdmp))
    println("-->PDMP result xdot = $xdotpdmp for x = $xc1, xd = $xd0")
    f_pdmp = @time  (tt,x,xdot)->PDMP.f_CHV!(f_ml_pdmp!, r_ml_pdmp!,tt,x,xdot,xd0,parms_chv_pdmp)
    @time f_pdmp(0.,xc1,xdotpdmp)
    f_chv = @time (tt,x,xdot,data)->f_ml_chv!(tt,x,xdot,eparms)
    @time f_chv(0.,xc1,xdotpdmp,eparms)
