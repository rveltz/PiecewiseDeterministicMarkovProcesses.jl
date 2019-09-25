# using RecursiveArrayTools
#
# include("problem.jl")
#
#
# struct DiffeqJumpWrapper
# 	diffeqjumppb
# 	u
# end
#
# # encode the vector field
# function (wrap::DiffeqJumpWrapper)(ẋ, xc, xd, p, t)
# 	# @assert 1==0 "WIP"
# 	# @show typeof(xc)
# 	# @show typeof(ẋ)
# 	# il faut des ExtendedJumpArray pour l'appel qui suit
# 	wrap.diffeqjumppb.prob.f(ẋ, xc, p, t)
# 	nothing
# end
#
# # encode the rate function
# function (wrap::DiffeqJumpWrapper)(rate, xc, xd, p, t, issum::Bool)
# 	for ii in eachindex(rate)
# 		rate[ii] = wrap.diffeqjumppb.variable_jumps[ii].rate(xc, p, t)
# 	end
# 	return sum(rate)
# end
#
# # encode the jump function
# function (wrap::DiffeqJumpWrapper)(xc, xd, p, t, ind_reaction::Int64)
# 	# this is a hack to be able to call affect! from DiffEqJump which require an integrator as an argument
# 	@show "jump function" xc xd wrap.u typeof(xc)
# 	wrap.u .= xc.u[1:end-1]
# 	wrap.diffeqjumppb.variable_jumps[ind_reaction].affect!(wrap)
# 	xc.u[1:end-1] .= wrap.u
# 	nothing
# end
#
# function PDMPProblem(jpprob::JumpProblem)
# 	@assert isinplace(jpprob.prob) "The current interface requires the ODE to be written inplace"
# 	@assert jpprob.regular_jump == nothing
# 	@assert jpprob.massaction_jump == nothing
#
# 	@show jpprob.prob.u0, typeof(jpprob.prob.u0)
#
# 	pb_wrapper = DiffeqJumpWrapper(jpprob, copy(jpprob.prob.u0.u))
#
# 	# get PDMP characteristics
# 	F = (xdot,xc,xd,p,t) -> pb_wrapper(xdot,xc,xd,p,t)
# 	R = (rate,xc,xd,p,t,issum) -> pb_wrapper(rate,xc,xd,p,t,issum)
# 	Delta = (xc,xd,p,t,ind_reaction) -> pb_wrapper(xc,xd,p,t,ind_reaction)
# 	xc0 = copy(jpprob.prob.u0)
#
# 	xd0 = [0]
# 	tspan = jpprob.prob.tspan
# 	p = jpprob.prob.p
#
# 	# determine the number of reactions
# 	nb_reactions = length(jpprob.variable_jumps)
#
# 	return PDMPProblem(F,R,Delta,nb_reactions,xc0,xd0,p,tspan)
# end
#
# function solve(jpprob::JumpProblem, algo::CHV{Tode};kwargs...) where {Tode <: DiffEqBase.DEAlgorithm}
# 	#issue with copy for eachindex(xc)...
# 	problem = PDMPProblem(jpprob)
# 	# X_extended = copy(jpprob.prob.u0)
# 	# resize!(X_extended.u, length(jpprob.prob.u0.u) + 1)
# 	X_extended = ArrayPartition(copy(jpprob.prob.u0), [0.])
# 	return solve(problem, algo, X_extended; kwargs...)
# end
