using PDMP
using Base.Test

cd(Pkg.dir("PDMP")*"/examples")

println("== start pdmp examples")

println("\n\n==== Example tcp ")
include("../examples/tcp.jl")

println(result2.time[end],", ",result2.xd[1,end])
@test isequal(result2.time[end],200.)
@test isequal(result2.xd[1,end],30)

println("\n\n==== Simple example of neuron model")
include("../examples/pdmp_example_eva.jl")
@test isequal(result.time[end],100.)
@test isequal(result.xd[2,end],93)

println("\n\n==== Morris-Lecar model of neuron")
include("../examples/morris_lecar.jl")

println("\n\n==== Example sir ")
include("../examples/sir.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],36)
@test isequal(result.xd[3,end],73)

println("\n\n==== Example sir(rejection) ")
include("../examples/sir-rejection.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],32)
@test isequal(result.xd[3,end],77)

println("\n\n==== Example neural network ")
include("../examples/neuron_rejection_exact.jl")

println(result2.time[end],", ",result2.xd[1,end])
@test isequal(result2.time[end],200.)
@test isequal(result2.xd[1,end],30)
@test isequal(size(result2.xd)[1],2)

println("\n\n==== Example sir ")
include("../examples/sir.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],36)
@test isequal(result.xd[3,end],73)
