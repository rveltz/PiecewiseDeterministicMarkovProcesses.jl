using PDMP
using Base.Test

# cd(Pkg.dir("PDMP")*"/examples")

println("== start pdmp examples")

println("\n\n==== Example tcp ")
include("../examples/tcp.jl")

println("\n\n==== Example tcp fast with types ")
println("----> To make it interesting, this mathematical example can explode in finite time, hence the warning")
include("../examples/tcp_fast.jl")

println("\n\n==== Simple example of neuron model, would be better handled with rejection method")
include("../examples/pdmp_example_eva.jl")

println("\n\n==== Example sir ")
include("../examples/sir.jl")
include("../examples/sir-rejection.jl")


println("\n\n==== Simple example of neuron model, Morris-Leccar")
include("../examples/morris_lecar.jl")

