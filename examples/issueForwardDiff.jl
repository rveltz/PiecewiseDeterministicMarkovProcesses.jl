using Revise, PiecewiseDeterministicMarkovProcesses, LinearAlgebra, Random, DifferentialEquations, Sundials
const PDMP = PiecewiseDeterministicMarkovProcesses

function F_fd!(ẋ, xc, xd, t, parms)
	# vector field used for the continuous variable
	if mod(xd[1], 2)==0
		ẋ[1] = 1 + xd[1] + 0
	else
		ẋ[1] = -xc[1]
	end
	nothing
end

rate_tcp(x) = 1/x

function R_fd!(rate, xc, xd, t, parms, sum_rate::Bool)
	rate[1] = 1.0 + rate_tcp(xd[1]) * xc[1]
	if sum_rate == false
		return 0., 0.
	else
		return sum(rate), 0.
	end
end

Dummy! = PDMP.Delta_dummy

xc0 = [ 1.0 ]
xd0 = [ 1  ]

nu_fd = [[1 0];[0 -1]]
parms = [0.0]

# works:
Random.seed!(12)
	res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!,Dummy!, nu_fd, parms, 0.0, 10.0; algo = :chv, ode = CVODE_Adams())

Random.seed!(12)
	res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!,Dummy!, nu_fd, parms, 0.0, 100.0; algo = :chv, n_jumps = 3000,   ode = Tsit5(), save_positions=(false,false))

# fail because of autodiff
Random.seed!(12)
	res =  @time PDMP.pdmp!(xc0, xd0, F_fd!, R_fd!,Dummy!, nu_fd, parms, 0.0, 10000.0; algo = :chv, n_jumps = 3000,   ode = AutoTsit5(Rosenbrock23(autodiff=true)), save_positions=(false,false))




# using StaticArrays
# sxc0 = @MVector [ 1.0 ]
# sxd0 = @MVector [1]
# res =  @time PDMP.chv_diffeq!(sxc0, sxd0, F_fd!, R_fd!,Dummy!, nu_fd, parms, 0.0, 10.0,false; n_jumps = 30,   ode = Tsit5(),save_positions=(false, false))
#
#
#
# similar(sxc0, Size(4))
#
# pbfunc = PDMP.PDMPFunctions(F_fd!, R_fd!,Dummy!)
# pbsim = PDMP.PDMPsimulation{eltype(xc0), eltype(xd0)}(0.,0.,0, 0.,[0., 0.], false, 0)
#
# pb = PDMP.PDMPProblem{eltype(xc0),eltype(xd0),typeof(xc0),typeof(xd0),typeof(nu_fd),typeof(parms),typeof(F_fd!),typeof(R_fd!),typeof(Dummy!)}(copy(xc0), copy(xd0), F_fd!, R_fd!,Dummy!,nu_fd, parms, 0., 1.,false,false)
#
#
