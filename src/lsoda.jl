"""
Allocate memory variable that contains the context to call LSODA.lsoda
"""
# function lsoda_ctx{T}(f::Base.Callable,r::Base.Callable,d::Array{Int64},p::Vector{T}, y0::Vector{Float64}, tspan::Vector{Float64}; reltol::Float64=1e-7, abstol::Float64=1e-8)
# 	neq = Int32(length(y0))
# 	userdata = nothing # no data passed by the function, we expected const type for now
# 	userfun = LSODA.UserFunctionAndData(f, userdata, neq)
#
# 	atol = ones{Float64}(neq)
# 	rtol = ones{Float64}(neq)
#
# 	yres = zeros(length(tspan), length(y0))
#
# 	if typeof(abstol) == Float64
# 		atol *= abstol
# 	else
# 		atol = copy(abstol)
# 	end
#
# 	if typeof(reltol) == Float64
# 		rtol *= reltol
# 	else
# 		rtol = copy(reltol)
# 	end
#
# 	opt = lsoda_opt_t()
# 	opt.ixpr = 0
# 	opt.rtol = pointer(rtol)
# 	opt.atol = pointer(atol)
# 	opt.itask = 1
#
# 	ctx = lsoda_context_t()
# 	ctx.function_ = LSODA.fex_c
# 	ctx.neq = neq
# 	ctx.state = 1
# 	ctx.data = pointer_from_objref(userfun)
# 	lsoda_prepare(ctx,opt)
# 	return ctx
# end
