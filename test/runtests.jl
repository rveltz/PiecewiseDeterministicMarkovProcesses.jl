using PiecewiseDeterministicMarkovProcesses, Test, LinearAlgebra, Random, Pkg

cd(Pkg.dir("PiecewiseDeterministicMarkovProcesses")*"/examples")
# import PiecewiseDeterministicMarkovProcesses; joinpath(dirname(pathof(PiecewiseD
# eterministicMarkovProcesses)), "..", paths...)

println("== start pdmp examples")

println("\n\n==== Example tcp ")
include("../examples/tcp.jl")
@test prod(isless.(errors[6:end],1e-4))

println("\n\n==== Example pdmp explosion ")
include("../examples/pdmpExplosion.jl")
@test prod(isless.(errors[1:6],1e-4))

println("\n\n==== Example pdmp stiff ")
include("../examples/pdmpStiff.jl")
@test prod(isless.(errors,1e-3))

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
