function cvode(f::Function, y0::Vector{Float64}, t::Vector{Float64}; reltol::Float64=1e-4, abstol::Float64=1e-6)
    # f, Callable object to be optimized of the form f(y::Vector{Float64}, fy::Vector{Float64}, t::Float64)
    #    where `y` is the input vector, and `fy` is the
    # y0, Vector of initial values
    # t, Vector of time values at which to record integration results
    # reltol, Relative Tolerance to be used (default=1e-4)
    # abstol, Absolute Tolerance to be used (default=1e-6)
    # return: a solution matrix with time steps in `t` along rows and
    #         state variable `y` along columns
    neq = length(y0)
    mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Int32, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{Function})), t[1], Sundials.nvector(y0))
    flag = Sundials.CVodeSetUserData(mem, f)
    flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
    flag = Sundials.CVDense(mem, neq)
    yres = zeros(length(t), length(y0))
    yres[1,:] = y0
    y = copy(y0)
    tout = [0.0]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        yres[k,:] = y
    end
    Sundials.CVodeFree([mem])
    return yres
end

function cvode{f}(::Type{f}, y0::Vector{Float64}, t::Vector{Float64}; reltol::Float64=1e-4, abstol::Float64=1e-6)
    # f, Callable object to be optimized of the form f(y::Vector{Float64}, fy::Vector{Float64}, t::Float64)
    #    where `y` is the input vector, and `fy` is the
    # y0, Vector of initial values
    # t, Vector of time values at which to record integration results
    # reltol, Relative Tolerance to be used (default=1e-4)
    # abstol, Absolute Tolerance to be used (default=1e-6)
    # return: a solution matrix with time steps in `t` along rows and
    #         state variable `y` along columns
    neq = length(y0)
    mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Int32, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{Function})), t[1], Sundials.nvector(y0))
    flag = Sundials.CVodeSetUserData(mem, f)
    flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
    flag = Sundials.CVDense(mem, neq)
    yres = zeros(length(t), length(y0))
    yres[1,:] = y0
    y = copy(y0)
    tout = [0.0]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        yres[k,:] = y
    end
    Sundials.CVodeFree([mem])
    return yres
end

function cvode(f, y0::Vector{Float64}, t::Vector{Float64}; reltol::Float64=1e-4, abstol::Float64=1e-6)
    # f, Callable object/function to be optimized of the form f(y::Vector{Float64}, fy::Vector{Float64}, t::Float64)
    #    where `y` is the input vector, and `fy` is the
    # y0, Vector of initial values
    # t, Vector of time values at which to record integration results
    # reltol, Relative Tolerance to be used (default=1e-4)
    # abstol, Absolute Tolerance to be used (default=1e-6)
    # return: a solution matrix with time steps in `t` along rows and
    #         state variable `y` along columns
    neq = length(y0)
    mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Int32, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{Function})), t[1], Sundials.nvector(y0))
    flag = Sundials.CVodeSetUserData(mem, f)
    flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
    flag = Sundials.CVDense(mem, neq)
    yres = zeros(length(t), length(y0))
    yres[1,:] = y0
    y = copy(y0)
    tout = [0.0]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        yres[k,:] = y
    end
    Sundials.CVodeFree([mem])
    return yres
end
