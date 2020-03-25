using DiffEqJump: AbstractAggregatorAlgorithm, NullAggregator

PDMPProblem(prob,jumps::ConstantRateJump;kwargs...) = PDMPProblem(prob,JumpSet(jumps);kwargs...)
PDMPProblem(prob,jumps::VariableRateJump;kwargs...) = PDMPProblem(prob,JumpSet(jumps);kwargs...)
PDMPProblem(prob,jumps::RegularJump;kwargs...) = PDMPProblem(prob,JumpSet(jumps);kwargs...)
PDMPProblem(prob,jumps::MassActionJump;kwargs...) = PDMPProblem(prob,JumpSet(jumps);kwargs...)
PDMPProblem(prob,jumps::DiffEqJump.AbstractJump...;kwargs...) = PDMPProblem(prob,JumpSet(jumps...);kwargs...)

PDMPProblem(prob,aggregator::AbstractAggregatorAlgorithm,jumps::ConstantRateJump;kwargs...) = PDMPProblem(prob,aggregator,JumpSet(jumps);kwargs...)
PDMPProblem(prob,aggregator::AbstractAggregatorAlgorithm,jumps::VariableRateJump;kwargs...) = PDMPProblem(prob,aggregator,JumpSet(jumps);kwargs...)
PDMPProblem(prob,aggregator::AbstractAggregatorAlgorithm,jumps::RegularJump;kwargs...) = PDMPProblem(prob,aggregator,JumpSet(jumps);kwargs...)
PDMPProblem(prob,aggregator::AbstractAggregatorAlgorithm,jumps::MassActionJump;kwargs...) = PDMPProblem(prob,aggregator,JumpSet(jumps);kwargs...)
PDMPProblem(prob,aggregator::AbstractAggregatorAlgorithm,jumps::DiffEqJump.AbstractJump...;kwargs...) = PDMPProblem(prob,aggregator,JumpSet(jumps...);kwargs...)
PDMPProblem(prob,jumps::JumpSet;kwargs...) = PDMPProblem(prob,NullAggregator(),jumps;kwargs...)

struct DiffeqJumpWrapper
	diffeqpb
	jumps
	u
end

# encode the vector field
function (wrap::DiffeqJumpWrapper)(ẋ, xc, xd, p, t::Float64)
	wrap.diffeqpb.f(ẋ, xc, p, t)
	nothing
end

# encode the rate function
function (wrap::DiffeqJumpWrapper)(rate, xc, xd, p, t::Float64, issum::Bool)
	for ii in eachindex(rate)
		rate[ii] = wrap.jumps.variable_jumps[ii].rate(xc, p, t)
	end
	return sum(rate)
end

# encode the jump function
function (wrap::DiffeqJumpWrapper)(xc, xd, p, t::Float64, ind_reaction::Int64)
	# this is a hack to be able to call affect! from DiffEqJump which requires an integrator as an argument. But if the type of affect! is not enforced, it should work
	@inbounds for ii=1:length(wrap.u)
		wrap.u[ii] = xc[ii]
	end
	wrap.jumps.variable_jumps[ind_reaction].affect!(wrap)
	@inbounds for ii=1:length(wrap.u)
		xc[ii] = wrap.u[ii]
	end
	nothing
end

function PDMPProblem(prob, aggregator::AbstractAggregatorAlgorithm, jumps::JumpSet;
	save_positions = typeof(prob) <: DiffEqBase.AbstractDiscreteProblem ? (false,true) : (true, true), kwargs...)

	@assert isinplace(prob) "The current interface requires the ODE to be written inplace"
	@assert jumps.regular_jump == nothing
	@assert jumps.massaction_jump == nothing

	pb_wrapper = DiffeqJumpWrapper(prob, jumps, copy(prob.u0))

	# get PDMP characteristics
	F = (xdot,xc,xd,p,t::Float64) -> pb_wrapper(xdot,xc,xd,p,t)
	R = (rate,xc,xd,p,t::Float64,issum::Bool) -> pb_wrapper(rate,xc,xd,p,t,issum)
	Delta = (xc,xd,p,t::Float64,ind_reaction::Int) -> pb_wrapper(xc,xd,p,t,ind_reaction)

	xc0 = copy(prob.u0)

	xd0 = [0]
	tspan = prob.tspan
	p = prob.p

	# determine the number of reactions
	nb_reactions = length(jumps.variable_jumps)

	return PDMPProblem(F,R,Delta,nb_reactions,xc0,xd0,p,tspan)
end
