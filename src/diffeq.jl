using RecursiveArrayTools

include("problem.jl")

function (chv::CHV{Tode})(xdot::ExtendedJumpArray, x::ExtendedJumpArray, prob::Tpb, t) where {Tode, Tpb <: PDMPCaracteristics}
	tau = x.jump_u[1]
	sr = prob.R(prob.ratecache.rate, x.u, prob.xd, prob.parms, tau, true)[1]
	prob.F(xdot.u, x.u, prob.xd, prob.parms, tau)
	xdot.jump_u[1] = 1.0 / sr
	mul!(xdot.u, 1.0 / sr, xdot.u)
	return nothing
end

struct DiffeqJumpWrapper
	diffeqjumppb
	u
end

# encode the vector field
function (wrap::DiffeqJumpWrapper)(ẋ, xc, xd, p, t::Float64)
	# The next call requires ExtendedJumpArray
	wrap.diffeqjumppb.prob.f(ẋ, xc, p, t)
	nothing
end

# encode the rate function
function (wrap::DiffeqJumpWrapper)(rate, xc, xd, p, t::Float64, issum::Bool)
	for ii in eachindex(rate)
		rate[ii] = wrap.diffeqjumppb.variable_jumps[ii].rate(xc, p, t)
	end
	return sum(rate)
end

# encode the jump function
function (wrap::DiffeqJumpWrapper)(xc, xd, p, t::Float64, ind_reaction::Int64)
	# this is a hack to be able to call affect! from DiffEqJump which requires an integrator as an argument. But if the type of affect! is not enforced, it should work
	copyto!(wrap.u, xc.u)
	wrap.diffeqjumppb.variable_jumps[ind_reaction].affect!(wrap.u)
	copyto!(xc.u, wrap.u)
	nothing
end

function PDMPProblem(jpprob::JumpProblem)
	@assert isinplace(jpprob.prob) "The current interface requires the ODE to be written inplace"
	@assert jpprob.regular_jump == nothing
	@assert jpprob.massaction_jump == nothing

	pb_wrapper = DiffeqJumpWrapper(jpprob, copy(jpprob.prob.u0))

	# get PDMP characteristics
	F = (xdot,xc,xd,p,t::Float64) -> pb_wrapper(xdot,xc,xd,p,t)
	R = (rate,xc,xd,p,t::Float64,issum::Bool) -> pb_wrapper(rate,xc,xd,p,t,issum)
	Delta = (xc,xd,p,t::Float64,ind_reaction::Int) -> pb_wrapper(xc,xd,p,t,ind_reaction)
	xc0 = copy(jpprob.prob.u0)

	xd0 = [0]
	tspan = jpprob.prob.tspan
	p = jpprob.prob.p

	# determine the number of reactions
	nb_reactions = length(jpprob.variable_jumps)

	return PDMPProblem(F,R,Delta,nb_reactions,xc0,xd0,p,tspan)
end

function solve(jpprob::JumpProblem, algo::CHV{Tode};kwargs...) where {Tode <: DiffEqBase.DEAlgorithm}
	#issue with copy for eachindex(xc)...
	problem = PDMPProblem(jpprob)
	# X_extended = copy(jpprob.prob.u0)
	# resize!(X_extended.u, length(jpprob.prob.u0.u) + 1)
	X_extended = ArrayPartition(copy(jpprob.prob.u0), [0.])
	X_extended = ExtendedJumpArray(copy(jpprob.prob.u0), [1.0])
	return solve(problem, algo, X_extended; kwargs...)
end
