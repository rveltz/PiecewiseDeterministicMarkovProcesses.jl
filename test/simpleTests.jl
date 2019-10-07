using PiecewiseDeterministicMarkovProcesses, Random
const PDMP = PiecewiseDeterministicMarkovProcesses
jp = PDMP.Jump(PDMP.Delta_dummy)
jp = PDMP.Jump(rand(Int64,2,2))

PDMP.finalize_dummy(1, 1, 1, 1, 1)
PDMP.PDMPResult(rand(2), rand(2), rand(2))
PDMP.PDMPResult(rand(2), rand(2), rand(2), rand(2), (false, false))

xc0 = rand(2)
xdot0 = rand(2,2)
PDMP.Phi_dummy(xdot0, xc0, 1, 1, 1)
