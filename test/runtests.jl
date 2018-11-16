using PDMP, Test, LinearAlgebra, Random, Pkg

cd(Pkg.dir("PDMP")*"/examples")

println("== start pdmp examples")

println("\n\n==== Example tcp ")
include("../examples/tcp.jl")

println(result2.time[end],", ",result2.xd[1,end])
@test isequal(result2.time[end],403.4606902370968)
@test isequal(result2.xd[1,end],61)
# @test isapprox(result2.time,result3.time)

println("\n\n==== Simple example of neuron model")
include("../examples/pdmp_example_eva.jl")
@test isequal(result.time[end],100.)
@test isequal(result.xd[2,end],91)

# println("\n\n==== Morris-Lecar model of neuron")
# include("../examples/morris_lecar.jl")

println("\n\n==== Example sir ")
include("../examples/sir.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],36)
@test isequal(result.xd[3,end],73)

println("\n\n==== Example sir(rejection) ")
include("../examples/sir-rejection.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],33)
@test isequal(result.xd[3,end],76)

println("\n\n==== Example neural network ")
include("../examples/neuron_rejection_exact.jl")

@test isequal(result.time[end],3176.5558333563877)
@test isequal(result.xd[1,end],408)
@test isequal(size(result.xd)[1],100)

println("\n\n==== Example sir ")
include("../examples/sir.jl")
@test isequal(result.xd[1,end],0)
@test isequal(result.xd[2,end],36)
@test isequal(result.xd[3,end],73)
