using PDMP
using Base.Test

# run cvode example
println("== start pdmp examples")

#include("../examples/morris_lecar.jl")
include("../examples/pdmp_example_eva.jl")
include("../examples/sir.jl")
include("../examples/tcp.jl")
include("../examples/tcp_fast.jl")

# run ida examples
#println("== start ida_Roberts example")
#include("../examples/ida_Roberts_simplified.jl")

#println("result at t=$(t[end]):")
#println(yout[end,:], "\n")

#println("== start ida_Heat2D example")
#include("../examples/ida_Heat2D.jl")

#println("result at t=$(t[end]):")
#println(yout[end,:], "\n")

# run kinsol example
#println("== start kinsol example")
#include("../examples/kinsol_mkin_simplified.jl")

#println("solution:")
#println(res)
#residual = ones(2)
#sysfn(res, residual)
#println("residual:")
#println(residual, "\n")

#@test abs(minimum(residual)) < 1e-5
