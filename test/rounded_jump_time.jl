module RoundedJumpTimeTests

using Test
using Random
using Sundials
import PiecewiseDeterministicMarkovProcesses as PDMP

@testset "CHV handles an isolated rounded physical jump time" begin
    start = 17530.0

    function simulate(first_rate, later_rate)
        drift!(du, u, xd, p, t) = (du[1] = 0.0)
        function rates!(rates, u, xd, p, t, as_sum)
            rates[1] = iszero(xd[1]) ? first_rate : later_rate
            return rates[1], 0.0
        end
        Random.seed!(71)
        problem = PDMP.PDMPProblem(drift!, rates!, ones(Int, 1, 1),
            [0.0], [0], nothing, (start, start + 0.1))
        return PDMP.solve(problem, PDMP.CHV(CVODE_BDF(linear_solver = :GMRES)); n_jumps = 3)
    end

    isolated = simulate(1e16, 1e4)
    @test isolated.time[1] == isolated.time[2] == start < isolated.time[3]
    @test vec(Array(isolated.xd)) == [0, 1, 2]

    ordinary = simulate(1e4, 1e4)
    @test all(diff(ordinary.time) .> 0)
    @test vec(Array(ordinary.xd)) == [0, 1, 2]

    @test_throws "Could not compute next jump time" simulate(1e16, 1e16)
end

end
