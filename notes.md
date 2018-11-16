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
  0.024538 seconds (149.75 k allocations: 5.009 MiB, 27.63% gc time)
 ----> lsoda
[0.0, 1.77942]
  0.004146 seconds (78.45 k allocations: 3.329 MiB)
 ----> DiffEq
[0.0, 1.77942]
  0.006381 seconds (137.60 k allocations: 2.185 MiB)
--> stopping time == tf? (not more) false
```