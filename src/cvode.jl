"""
Allocate memory variable that contains the context to call Sundials.cvode
"""
function cvode_ctx{T}(f::Base.Callable,r::Base.Callable,d::Array{Int64},p::Vector{T}, y0::Vector{Float64}, t::Vector{Float64}; reltol::Float64=1e-7, abstol::Float64=1e-8)
    neq = length(y0)
    mem_ptr = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    (mem_ptr == C_NULL) && error("Failed to allocate CVODE solver object")
    mem = Sundials.Handle(mem_ptr)
    Sundials.@checkflag Sundials.CVodeInit(mem, cfunction(cvode_ode_wrapper, Cint, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Array{Any,1})), t[1], Sundials.nvector(y0))
    Sundials.@checkflag Sundials.CVodeSetUserData(mem, [f,r,d,p])
    Sundials.@checkflag Sundials.CVodeSStolerances(mem, reltol, abstol)
    # Sundials.@checkflag Sundials.CVDense(mem, neq)
    A = Sundials.SUNDenseMatrix(neq,neq)
    mat_handle = Sundials.MatrixHandle(A,Sundials.DenseMatrix())
    LS = Sundials.SUNDenseLinearSolver(convert(Sundials.N_Vector,y0),A)
    LS_handle = Sundials.LinSolHandle(LS,Sundials.Dense())
    Sundials.@checkflag Sundials.CVDlsSetLinearSolver(mem, LS, A)
    return (mem,mat_handle,LS_handle)
end

"""
This functions allows to save re-allocating internal variables to call Sundials.CVode unlike cvode() above.
"""
function cvode_evolve!{T}(yres::Array{Float64,2}, mem, f::Base.Callable,r::Base.Callable,d::Array{Int64},p::Vector{T},y0::Vector{Float64}, t::Vector{Float64})
    # How do I update the parameter d in mem??
    Sundials.CVodeReInit(mem,t[1],y0)
    Sundials.CVodeSetUserData(mem, [f,r,d,p])

    yres[1,:] = y0
    y = copy(y0)
    tout = [0.0]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        if flag != Sundials.CV_SUCCESS
            throw(KeyError("SUNDIALS_ERROR: CVODE failed with flag = ", flag))
        end
        # BAD!! should be copy cols, much faster
        yres[k,:] = y
    end
    nothing
end