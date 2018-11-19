- use [RecursiveArrayTools.jl](https://github.com/JuliaDiffEq/RecursiveArrayTools.jl) pour faire comme DiffEq package dans `PDMP.jl`


```julia
julia> integ.
EEst             dtpropose         isout             p                 t
accept_step      eigen_est         iter              q11               tdir
alg              erracc            just_hit_tstop    qold              tprev
cache            event_last_time   k                 reeval_fsal       u
dt               f                 kshortsize        saveiter          u_modified
dtacc            force_stepfail    last_event_error  saveiter_dense    uprev
dtcache          fsalfirst         last_stepfail     sol               uprev2
dtchangeable     fsallast          opts              success_iter
```

```
julia> include("/Users/rveltz/work/prog_gd/julia/dev/PDMP.jl/examples/tcp.jl")
--> inplace implementation,
 ----> cvode
[0.0, 1.77942]
  0.178954 seconds (1.27 M allocations: 42.461 MiB, 8.29% gc time)
 ----> lsoda
[0.0, 1.77942]
  0.039761 seconds (618.30 k allocations: 27.093 MiB, 18.10% gc time)
 ----> DiffEq
[0.0, 1.80003]
  0.008713 seconds (456 allocations: 32.188 KiB)
--> stopping ti
```

Ca n'alloue rien

```
probpdmp = PDMP.PDMPProblem{eltype(xc0),eltype(xd0),typeof(xc0),typeof(xd0),typeof(nu_tcp),typeof(parms),typeof(F_tcp!),typeof(R_tcp!),typeof(PDMP.Delta_dummy)}(xc0,xd0,F_tcp!, R_tcp!,PDMP.Delta_dummy,nu_tcp,parms,0.,tf,false,false)


using BenchmarkTools

xext = copy(xc0);push!(xext,0.)
xdext = copy(xext)
@benchmark probpdmp($xdext,$xext,1.,0.1)

```

```
using Profile, ProfileView
Profile.clear()
Random.seed!(1234)
result4 =  @profile PDMP.chv_diffeq!(xc0,xd0,
                F_tcp!,R_tcp!,PDMP.Delta_dummy,
                nu_tcp,parms,0.0,tf,false, n_jumps = 10000,ode = Tsit5(), save_positions = (false,false))
ProfileView.view()

```


julia> result4[1].time
10000-element Array{Float64,1}:
     0.0
     1.8000309514727386
     2.8604924438144486
     3.215069241279198
     9.268979277627087
    14.922571800250441
    16.143446475076303
    19.608371758286594
    19.82229503436417
    20.555376192292787
    21.910777385413176
    21.928158429061654
    
    
julia> result3.time
10000-element Array{Float64,1}:
     0.0
     1.7794207026823572
     2.826532479449366
     3.177121211156773
     9.192781711107573
    14.813270132833463
    16.017775265374073
    19.447264313819183
    19.660487514099163
    20.390544445409063
    21.741504121267866
    21.75885590162066
    
julia> result2.time
10000-element Array{Float64,1}:
     0.0
     1.779420704593002
     2.8265326138597446
     3.177123960582396
     9.192791598750683
    14.81328411466162
    16.017791949820154
    
```julia
using DifferentialEquations, OrdinaryDiffEq, Plots
f = @ode_def begin
    dGut = -k*Gut
end k
condition = function (u,t,integrator)
  t == 5
end
affect! = function (integrator)
  @show integrator.u[1] += 0.05
end

function evolve!(integrator)
	while integrator.t <5
		step!(integrator)
	end
	integrator.sol
end

cb = DiscreteCallback(condition, affect!)
u0 = [0.1]
prob = ODEProblem(f, u0, (0,50.), 2)
sol = solve(prob, Tsit5(), tstops = 5, callback=cb);
sol = @time solve(prob, Tsit5(), tstops=5, callback=cb);

sol = solve(prob, Tsit5(), tstops = 5);
sol = @time solve(prob, Tsit5(), tstops=5);

integ = init(prob, Tsit5(), tstops=5);
sol = @time evolve!(integ);

```