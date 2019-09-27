using PiecewiseDeterministicMarkovProcesses, Test, LinearAlgebra, Random, DifferentialEquations

@testset "Example TCP" begin
	include("../examples/tcp.jl")
	@test prod(isless.(errors[6:end],1e-4))
end

@testset "Example explosion of jump times" begin
	include("../examples/pdmpExplosion.jl")
	@test prod(isless.(errors[1:6],1e-4))
end

@testset "Example with stiff ODE part" begin
	include("../examples/pdmpStiff.jl")
	@test prod(isless.(errors,1e-3))
end

@testset "Controlling allocations" begin
	# @test alloc1 == alloc2
end

@testset "Example with 2d example, for autodiff" begin
	include("../examples/tcp2d.jl")
	# @test prod(isless.(errors,1e-3))
end

@testset "Rejection method" begin
	include("../examples/tcp_rejection.jl")
	@test prod(isless.(errors,1e-5))
end

@testset "Neuron model" begin
	include("../examples/pdmp_example_eva.jl")
	@test isequal(result1.time[end],100.)
	@test isequal(result1.xd[2,end],113)
end

@testset "Example SIR" begin
	include("../examples/sir.jl")
	@test isequal(result.xd[1,end],0)
	@test isequal(result.xd[2,end],36)
	@test isequal(result.xd[3,end],73)
end

@testset "Example SIR(rejection)" begin
	include("../examples/sir-rejection.jl")
	@test isequal(result.xd[1,end],0)
	@test isequal(result.xd[2,end],33)
	@test isequal(result.xd[3,end],76)
end

@testset "Neural network" begin
	include("../examples/neuron_rejection_exact.jl")
	@test isequal(result.time[end],857.6850502747997)
	@test isequal(result.xd[1,end],98)
	@test isequal(size(result.xd)[1],2)
end
