var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.jl",
    "category": "page",
    "text": ""
},

{
    "location": "#PiecewiseDeterministicMarkovProcesses.jl-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.jl",
    "category": "section",
    "text": "PiecewiseDeterministicMarkovProcesses.jl is a Julia package that allows simulation of Piecewise Deterministic Markov Processes (PDMP); these encompass hybrid systems and jump processes, comprised of continuous and discrete components, as well as processes with time-varying rates. The aim of the package is to provide methods for the simulation of these processes that are \"exact\" up to the ODE integrator.We briefly recall facts about a simple class of PDMPs. They are described by a couple (x_cx_d) where x_c is solution of the differential equation fracdx_cdt = F(x_cx_dt). The second component x_d is a piecewise constant variable with type Int64. The jumps occurs at rates R(x_cx_dt). At each jump, x_d or x_c can be affected.We provide several methods for the simulation:a recent trick, called CHV, explained in paper-2015 which allows to implement the True Jump Method without the need to use event detection schemes for the ODE integrator. These event detections can be quite unstable as explained in paper-2015 and CHV provide a solution to this problem.\nrejection methods for which the user is asked to provide a bound on the reaction rates. These last methods are the most \"exact\" but not the fastest if the reaction rate bound is not tight. In case the flow is known analytically, a method is also provided.These methods require solving stiff ODEs (for CHV ) in an efficient manner. Sundials.jl and LSODA.jl are great, but other solvers can also be considered (see stiff ode solvers). Hence, the current package allows the use of all solvers in DifferentialEquations.jl thereby giving access to a wide range of solvers. In particular, we can test different solvers to see how precise they are. Here is an example from examples/pdmpStiff.jl for which an analytical expression is available which allows computation of the errorsComparison of solvers\n--> norm difference = 0.00019114008823351014  - solver = cvode\n--> norm difference = 0.00014770067837588385  - solver = lsoda\n--> norm difference = 0.00018404736432131585  - solver = CVODEBDF\n--> norm difference = 6.939603217404056e-5    - solver = CVODEAdams\n--> norm difference = 2.216652299580346e-5    - solver = tsit5\n--> norm difference = 2.2758951345736023e-6   - solver = rodas4P-noAutoDiff\n--> norm difference = 2.496987313804766e-6    - solver = rodas4P-AutoDiff\n--> norm difference = 0.0004373003700521849   - solver = RS23\n--> norm difference = 2.216652299580346e-5    - solver = AutoTsit5RS23note: ODE Solvers\nA lot of care have been taken to be sure that the algorithms do not allocate and hence are fast. This is based on an iterator interface of DifferentialEquations. If you chose save_positions = (false, false), the allocations should be independent from the requested jump number. However, the iterator solution is not yet available for LSODA in DifferentialEquations. Hence you can pass ode = :lsoda to access an old version of the algo (which allocates) or any other solver like ode = Tsit5() to access the new solver."
},

{
    "location": "#Installation-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Installation",
    "category": "section",
    "text": "To install this package, run the command add PiecewiseDeterministicMarkovProcesses"
},

{
    "location": "#Basic-example-with-CHV-method-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Basic example with CHV method",
    "category": "section",
    "text": "A strong requirement for the CHV method is that the total rate (i.e. sum(rate)) must be positive. This can be easily achieved by adding a dummy Poisson process with very low intensity (see next section).See also the examples directory for more involved examples. A simple example of jump process is given below. More precisely, we look at the following process of switching dynamics where X(t) = (x_c(t) x_d(t)) inmathbb Rtimeslbrace-11rbrace In between jumps, x_c evolves according to dot x_c(t) = x_d(t)x_c(t)  We first need to load the library.  using PiecewiseDeterministicMarkovProcessesWe then define a function that encodes the dynamics in between jumps. We need to provide the vector field of the ODE. Hence, we need to define a function that, given continuous state x_c and discrete state x_d at time t, returns the vector field. In addition some parameters can be passed with the variable parms.function F_tcp!(xcdot, xc, xd, t, parms)\n  # vector field used for the continuous variable\n  xcdot[1] = xd[1]*xc[1]\nend Let\'s consider a stochastic process with following transitions:Transition Rate Reaction number Jump\nx_dto x_d-2 if x_d0 1 1 [-2]\nx_dto x_d+2 if x_d0 1 2 [2]We implement these jumps using a 2x1 matrix nu of Integers, such that the jumps on each discrete component of xd is given by nu * xd. Hence, we have nu = reshape([[2];[-2]],2,1).	note: Implementing jumps\nThere are two ways to implement jumps. The first one is to provide a transition matrix nu to the solver but this will only implement jumps on the discrete variable xd and leaves xc unaffected. The more general way is to implement a function Delta!(xc, xd, t::Float64, parms, ind_reaction::Int64) in which you write the jump. See examples/pdmp_example_eva.jl for an example.These reactions with their rate are encoded in the following function.function R_tcp!(rate, xc, xd, t, parms, sum_rate::Bool)\n  # transition rates function for each transition\n  # in this case,  the transitions are xd->xd+2 or xd->xd-2\n  # sum_rate is a boolean which tells R_tcp if it needs to return the total reaction rates, this may \n  # i.e. the sum of the rates or the vector of the rates\n  if sum_rate == false\n      if xd[1] > 0\n          rate[1] = 0.\n          rate[2] = 1.\n      else\n      	  rate[1] = 1.\n          rate[2] = 0.\n      end\n      #we return 0. because nothing is supposed to be returned\n      return 0.\n  else\n  	# we return sum(rate) without altering rate as we are asked to do\n    return 1.\n  end\nend\n\n# initial conditions for the continuous/discrete variables\nxc0 = vec([0.05])\nxd0 = vec([1])\n\n# matrix of jumps for the discrete variables, analogous to chemical reactions\nnu = reshape([[2];[-2]],2,1)\n\n\n# parameters\nparms = [0.]\ntf = 25.\n\n# compile the program:\ndummy =  PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu,parms,0.0,tf,n_jumps=1)\n\n# compute a trajectory, in this case 100 jumps\nsrand(123)\nresult =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_tcp!,R_tcp!,nu,parms,0.0,tf,n_jumps=100)\n\n# plotting\nusing Plots\nPlots.plot(result.time, result.xd[1,:],line=:step,title = string(\"#Jumps = \",length(result.time)),label=\"Xd\")\nPlots.plot(result.time, result.xc\',title = string(\"#Jumps = \",length(result.time)),label=\"Xc\")This produces the following graph:(Image: TCP)"
},

{
    "location": "#Adding-more-sampling-points-in-between-jumps-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Adding more sampling points in between jumps",
    "category": "section",
    "text": "The current interface \"only\" returns the jump times. On may want to resolve the trajectory in between jumps. For example, in the previous example, in between two jumps, the trajectory should be exponential and not linear as shown. A simple trick to do this is to add a Poisson process to the reactions set with a given sampling rate. We have to modify nu, xcd0 and R_tcp! for this. The set of reactions is now the followingTransition Rate Jump\nx_d1to x_d1-2 if x_d10 1 [-2,0]\nx_d1to x_d1+2 if x_d10 1 [2,0]\nx_d2to x_d2+1 rate_save [0,1]Hence, we implement these jumps with the following matrix: nu2 = [[2 0];[-2 0];[0 1]].nu2 = [[2 0];[-2 0];[0 1]]\n# the second component is the Poisson process\nxd0 = vec([1, 0])\n\nfunction R_tcp2!(rate, xc, xd, t, parms, sum_rate::Bool)\n  # transition rates function for each transition\n  # in this case,  the transitions are xd->xd+2 or xd->xd-2\n  # sum_rate is a boolean which tells R_tcp if it needs to return the total reaction rates, this may \n  # i.e. the sum of the rates or the vector of the rates\n  rate_save = 10. #sampling rate in between true jumps\n  if sum_rate == false\n      if xd[1] > 0\n          rate[1] = 0.\n          rate[2] = 1.\n          rate[3] = rate_save #Poisson process used as sampling process\n      else\n          rate[1] = 1.\n          rate[2] = 0.\n          rate[3] = rate_save #Poisson process used as sampling process\n      end\n      #we return 0. because nothing is supposed to be returned\n      return 0.\n  else\n    # we see that we effectively return sum(rate) without altering rate because it is not asked to do so\n    return 1. + rate_save\n  end\nend\n\nsrand(123)  \nresult2 =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_tcp!,R_tcp2!,nu2,parms,0.0,tf,n_jumps=10000)\nPlots.plot(result2.time, result2.xc\',title = string(\"#Jumps = \",length(result2.time)),label=\"Xc2\")This gives the following result:(Image: TCP)"
},

{
    "location": "#Basic-example-with-the-rejection-method-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Basic example with the rejection method",
    "category": "section",
    "text": "The previous method is useful when the total rate function varies a lot. In the case where the total rate is mostly constant in between jumps, the rejection method is more appropriate. The rejection method assumes some a priori knowledge of the process one wants to simulate. In particular, the user must be able to provide a bound on the total rate. More precisely, the user must provide a constant bound in between jump. To use this method, one needs to return sum(rate), bound_rejection in the above function R_tcp!. Note that this means that in between jumps, one have:sum(rate)(t) <= bound_rejectionnu2 = [[2 0];[-2 0];[0 1]]\n# the second component is the Poisson process\nxd0 = vec([1, 0])\n\nfunction R_tcp2!(rate, xc, xd, t, parms, sum_rate::Bool)\n  # transition rates function for each transition\n  # in this case,  the transitions are xd->xd+2 or xd->xd-2\n  # sum_rate is a boolean which tells R_tcp if it needs to return the total reaction rates, this may \n  # i.e. the sum of the rates or the vector of the rates\n  rate_save       = 10.           # sampling rate in between true jumps\n  bound_rejection = 1.+rate_save  # bound on the total rate, here 0 + 1 + rate_save\n  if sum_rate == false\n      if xd[1] > 0\n          rate[1] = 0.\n          rate[2] = 1.\n          rate[3] = rate_save #Poisson process used as sampling process\n      else\n          rate[1] = 1.\n          rate[2] = 0.\n          rate[3] = rate_save #Poisson process used as sampling process\n      end\n      #we return 0. because nothing is supposed to be returned\n      return 0., bound_rejection\n  else\n    # we see that we effectively return sum(rate) without altering rate because it is not asked to do so\n    return 1. + rate_save, bound_rejection\n  end\nendWe can now simulate this process as followssrand(123)\nresult3 =  @time PiecewiseDeterministicMarkovProcesses.pdmp!(xc0,xd0,F_tcp!,R_tcp2!,nu2,parms,0.0,tf,n_jumps=10000,algo=:rejection)\nPlots.plot(result3.time, result3.xc\',title = string(\"#Jumps = \",length(result3.time)),label=\"rejection\")"
},

{
    "location": "#How-to-chose-a-simulation-method?-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "How to chose a simulation method?",
    "category": "section",
    "text": "The choice of the method CHV vs Rejection only depends on how much you know about the system. More precisely, if the total rate function does not vary much in between jumps, use the rejection method. For example, if the rate is R(x_c(t)) = 1+01cos(t),  then 1+01 will provide a tight bound to use for the rejection method and almost no (fictitious) jumps will be rejected. In all other cases, one should try the CHV method where no a priori knowledge of the rate function is requied."
},

{
    "location": "#Advanced-uses-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Advanced uses",
    "category": "section",
    "text": ""
},

{
    "location": "#Specify-a-jump-with-a-function-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Specify a jump with a function",
    "category": "section",
    "text": "See examples/pdmp_example_eva.jl for an example."
},

{
    "location": "#Application-programming-interface-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Application programming interface",
    "category": "section",
    "text": ""
},

{
    "location": "#PiecewiseDeterministicMarkovProcesses.pdmp!",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.pdmp!",
    "category": "function",
    "text": "This function performs a pdmp simulation using the Change of Variable (CHV, see https://arxiv.org/abs/1504.06873) method or the rejection method. 	It takes the following arguments:\n\npdmp!(xc0, xd0, F!, R!, DX, nu, parms, ti, tf; verbose::Bool = false, ode = :cvode, algo=:chv, n_jumps = 30_000, save_positions = (false, true), saverate = false, save_at = [])\n\nIt takes the arguments:\n\nxc0: a Vector of Float64, representing the initial states of the continuous variable.\nxd0: a Vector of Int64, representing the initial states of the discrete variable.\nF!: an inplace Function or a callable type, which itself takes five arguments to represent the vector field; xdot a Vector of Float64 representing the vector field associated to the continuous variable, xc Vector representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time and parms, a Vector of Float64 representing the parameters of the system. F!(xdot,xc,xd,t,parms) returns nothing\nR!: an inplace Function or a callable type, which itself takes six arguments to represent the rate functions associated to the jumps;rate Vector of Float64 holding the different reaction rates, xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and sumrate a Bool being a flag asking to return a Float64 if true and a Vector otherwise. `R!(rate,xc,xd,t,parms,sumrate)returnsFloat64,Float64`\nDX: a Function or a callable type, which itself takes five arguments to apply the jump to the continuous/discrete variable;xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and indrec an Int64 representing the index of the discrete jump. `DX(xc,xd,t,parms,indrec)returnsnothing`\nnu: a Matrix of Int64, representing the transitions of the system, organised by row.\nparms : data for the parameters of the system. It is passed to F!, R! and DX.\nti: the initial simulation time (Float64)\ntf: the final simulation time (Float64)\nverbose: a Bool for printing verbose.\node: ode time stepper :cvode, :lsoda or any solver from DifferentialEquations.jl, like CVODE_BDF().\nn_jumps: an Int64 representing the maximum number of jumps to be computed. Can be put to Inf64.\nind_save_d: a range to hold the indices of the discrete variable to be saved\nind_save_c: a range to hold the indices of the continuous variable to be saved\nsave_positions = (false,true) indicates whether to save the pre-jump (resp. post-jump) state\nsave_at: an ordered list of times at which the user want to save the state. Work in progress: not available for every algorithm.\nsaverate::Bool indicates whether the total rates at each jump is saved. Useful to study how to tidy the bound on the rates in the rejection algorithm.\n\n\n\n\n\n"
},

{
    "location": "#PiecewiseDeterministicMarkovProcesses.chv_diffeq!",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.chv_diffeq!",
    "category": "function",
    "text": "Implementation of the CHV method to sample a PDMP using the package DifferentialEquations. The advantage of doing so is to lower the number of calls to solve using an integrator method.\n\n\n\n\n\n"
},

{
    "location": "#PiecewiseDeterministicMarkovProcesses.chv!",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.chv!",
    "category": "function",
    "text": "chv!\n\nThis function performs a pdmp simulation using the Change of Variable (CHV) method see https://arxiv.org/abs/1504.06873. It takes the following arguments:\n\nn_max: an Int64 representing the maximum number of jumps to be computed.\nxc0 : a Vector of Float64, representing the initial states of the continuous variable.\nxd0 : a Vector of Int64, representing the initial states of the discrete variable.\nF! : an inplace Function or a callable type, which itself takes five arguments to represent the vector field; xdot a Vector of Float64 representing the vector field associated to the continuous variable, xc Vector representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time and parms, a Vector of Float64 representing the parameters of the system.\nR : an inplace Function or a callable type, which itself takes six arguments to represent the rate functions associated to the jumps;rate Vector of Float64 holding the different reaction rates, xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and sum_rate a Bool being a flag asking to return a Float64 if true and a Vector otherwise.\nDX : a Function or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and ind_rec an Int64 representing the index of the discrete jump.\nnu : a Matrix of Int64, representing the transitions of the system, organised by row.\nparms : data for the parameters of the system.\ntf : the final simulation time (Float64)\nverbose : a Bool for printing verbose.\node: ode time stepper, must be one of those: [:cvode,:lsoda,:Adams,:BDF]\nsave_at: array of ordered time at which the solution is required\n\n\n\n\n\n"
},

{
    "location": "#PiecewiseDeterministicMarkovProcesses.rejection_exact",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "PiecewiseDeterministicMarkovProcesses.rejection_exact",
    "category": "function",
    "text": "rejection_exact\n\nThis function performs a simulation using the rejection method when the flow is known analytically. It takes the following arguments:\n\nn_max: an Int64 representing the maximum number of jumps to be computed.\nxc0 : a Vector of Float64, representing the initial states of the continuous variable.\nxd0 : a Vector of Int64, representing the initial states of the discrete variable.\nPhi! : a Function or a callable type, which itself takes 6 arguments to represent the vector field; rate a Vector of Float64 representing the flow of the vector which needs to be filled with values of the rates, xdot a Vector of Float64 representing the vector field associated to the continuous variable, xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time and parms, a Vector of Float64 representing the parameters of the system, sumofrate a Bool stating if the function must return the total rate.\nR! : a Function or a callable type, which itself takes five arguments to represent the rate functions associated to the jumps;xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and sumrate a Bool being a flag asking to return a Float64 if true and a Vector otherwise. The returned vector has components. If sumrate is False, one must return ratevector, bound where bound_ is a bound on the total rate vector. In the case sumrate is True, one must return totalrate,bound_ where totalrate is a Float64 that is the sum of the rates. In any case, the function must return a couple (totalrates, bound) where bound is a bound for the total rate.\nDelta : a Function or a callable type, which itself takes five arguments to apply the jump to the continuous variable;xc Vector of Float64 representing the current state of the continuous variable, xd Vector of Int64 representing the current state of the discrete variable, t a Float64 representing the current time, parms a Vector of Float64 representing the parameters of the system and ind_rec an Int64 representing the index of the discrete jump.\nnu : a Matrix of Int64, representing the transitions of the system, organised by row.\nparms : data for the parameters of the system.\ntf : the final simulation time (Float64)\nverbose : a Bool for printing verbose.\n\n\n\n\n\n"
},

{
    "location": "#Functions-1",
    "page": "PiecewiseDeterministicMarkovProcesses.jl",
    "title": "Functions",
    "category": "section",
    "text": "pdmp!chv_diffeq!chv!rejection_exact"
},

]}
