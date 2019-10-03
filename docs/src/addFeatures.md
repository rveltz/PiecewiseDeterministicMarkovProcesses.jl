## Specify a jump with a function
See `examples/pdmp_example_eva.jl` for an example.

## Rejection method stopped, recover data!

If you chose an upper bound for the rejection method that is too small and triggers an interruption like

```julia
ERROR: AssertionError: Error, your bound on the rates is not high enough!, [26.730756983739408, 20.0]
```

the `solve` does not return anything. However, in order to understand why your bound is too small, you would like to have a look at your trajectory up to the point where you bound failed. Don't worry, your computation is still in memory!

If your call is like this:

```
sol = solve(problem, Rejection(Tsit5()) )
```
then the trajectory is saved in the variables `problem.time`, `problem.Xc` and `problem.Xd`.