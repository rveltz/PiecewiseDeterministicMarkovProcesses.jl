using PiecewiseDeterministicMarkovProcesses, Test, LinearAlgebra, Random, DifferentialEquations

@testset "Example TCP" begin
	include("../examples/tcp.jl")
	@test norm(errors[6:end], Inf64) < 1e-4
end

@testset "Example explosion of jump times" begin
	include("../examples/pdmpExplosion.jl")
	@test norm(errors[1:6], Inf64) < 1e-4
end

@testset "Example with stiff ODE part" begin
	include("../examples/pdmpStiff.jl")
	@test norm(errors, Inf64) < 1e-3
	@test restime1 == res12.time
end

@testset "Controlling allocations" begin
	# @test alloc1 == alloc2
end

@testset "Test Rate structures" begin
	include("testRatesCst.jl")
end

@testset "Example with 2d example, for autodiff" begin
	include("../examples/tcp2d.jl")
end

@testset "Rejection method" begin
	include("../examples/tcp_rejection.jl")
	@test norm(errors, Inf64) < 1e-5
end

@testset "Neuron model" begin
	include("../examples/pdmp_example_eva.jl")
	@test result1.time[end] == 100.
	@test result1.xd[2,end] == 113
end

@testset "Example SIR" begin
	include("../examples/sir.jl")
	@test result.xd[1,end] == 0
	@test result.xd[2,end] == 36
	@test result.xd[3,end] == 73
end

@testset "Example SIR(rejection)" begin
	include("../examples/sir-rejection.jl")
	@test result.xd[1,end] == 0
	@test result.xd[2,end] == 33
	@test result.xd[3,end] == 76
end

@testset "Neural network" begin
	include("../examples/neuron_rejection_exact.jl")
	@test result.xd[1,end] == 98
	@test size(result.xd)[1] == 2
end
