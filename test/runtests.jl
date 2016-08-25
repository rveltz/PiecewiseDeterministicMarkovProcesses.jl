using PDMP
using Base.Test

cd(Pkg.dir("PDMP")*"/examples")

println("== start pdmp examples")

println("\n\n==== Example tcp ")
include("../examples/tcp.jl")
@test isequal(result.time[end],200.)
@test isequal(result.xd[1,end],24)


println("\n\n==== Example tcp fast with types ")
println("----> To make it interesting, this mathematical example can explode in finite time, hence the warning")
include("../examples/tcp_fast.jl")
@test isequal(length(result.time),20)
@test isequal(result.xd[1,end],3)

println("\n\n==== Simple example of neuron model, would be better handled with rejection method")
include("../examples/pdmp_example_eva.jl")
@test isequal(result.time[end],100.)
@test isequal(result.xd[2,end],1025)


println("\n\n==== Example sir ")
include("../examples/sir.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],34)
@test isequal(result.xd[3,end],75)

include("../examples/sir-rejection.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],73)
@test isequal(result.xd[3,end],36)

println("\n\n==== Simple example of neuron model, Morris-Leccar")
include("../examples/morris_lecar.jl")
@test isequal(result.xd[1,end],25)
@test isequal(result.xd[2,end],0)
@test isequal(result.xd[3,end],53)
@test isequal(result.xd[4,end],9)
