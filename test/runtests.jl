using PiecewiseDeterministicMarkovProcesses, Test, LinearAlgebra, Random, DifferentialEquations

macro testS(label, args...)
	:(@testset $label begin @test $(args...); end)
end

include("simpleTests.jl")

@testset "Example TCP" begin
	include("../examples/tcp.jl")
	@test norm(errors[6:end], Inf64) < 1e-4
end

@testset "Example explosion of jump times" begin
	include("../examples/pdmpExplosion.jl")
	@test norm(errors[1:6], Inf64) < 1e-4
end

@testset "Example with stiff ODE part" begin
	include("pdmpStiff.jl")
	@test minimum(errors) < 1e-3
	@testS "Call many times the same problem" restime1 == res12.time
end

@testset "Controlling allocations" begin
	@test alloc1 == alloc2
end

@testset "Test Rate structures 1/2" begin
	include("testRatesCst.jl")
end

@testset "Test Rate structures 2/2" begin
	include("testRatesComposite.jl")
end

@testset "Example with 2d example, for autodiff" begin
	include("../examples/tcp2d.jl")
end

@testset "Rejection method" begin
	include("../examples/tcp_rejection.jl")
	@test norm(errors, Inf64) < 1e-5
end

@testset "Test number of fictitious jumps" begin
	@test res1.njumps == res2.njumps
	@test res1.nrejected == res2.nrejected
end

@testset "Neuron model" begin
	include("../examples/pdmp_example_eva.jl")
	@test result1.time[end] == 100.
	@test result1.xd[2,end] == 107
end

@testset "Example SIR" begin
	include("../examples/sir.jl")
	@test result.xd[1,end] == 0
	@test result.xd[2,end] == 29
	@test result.xd[3,end] == 80
end

@testset "Example SIR(rejection)" begin
	include("../examples/sir-rejection.jl")
	@test result.xd[1,end] == 0
	@test result.xd[2,end] == 26
	@test result.xd[3,end] == 83
end

@testset "Neural network" begin
	include("../examples/neuron_rejection_exact.jl")
	@test result.xd[1,end] == 100
	@test size(result.xd)[1] == 100
end

@testset "JumpProcesses Wrap" begin
	include("../examples/examplediffeqjumpwrapper.jl")
end
