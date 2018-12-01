using Revise, PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations

function F_fd!(ẋ, xc, xd, t, parms)
    # vector field used for the continuous variable
    if mod(xd[1],2)==0
        ẋ[1] = 1 + xd[1] + 0
    else
        ẋ[1] = -xc[1]
    end
    nothing
end

rate_tcp(x) = 1/x

function R_fd!(rate, xc, xd, t, parms, sum_rate::Bool)
    rate[1] = 1.0 + rate_tcp(xc[1]) * xc[1]
    if sum_rate==false
        return 0., 0.
    else
        return sum(rate), 0.
    end
end

xc0 = [ 1.0 ]
xd0 = [ 1.  ]

nu_fd = [[1 0];[0 -1]]
parms = [0.0]

# works:
Random.seed!(12)
    res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!, nu_fd, parms, 0.0, 10.0, n_jumps = 3,   ode = CVODE_BDF())

# fail because of autodiff
Random.seed!(12)
    res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!, nu_fd, parms, 0.0, 10.0, n_jumps = 3,   ode = Rosenbrock23())




# using StaticArrays
# sxc0 = @MVector [ 1.0 ]
# sxd0 = @MVector [1, 1]
# res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!, nu_fd, parms, 0.0, 1.0, n_jumps = 3,   ode = Tsit5())
