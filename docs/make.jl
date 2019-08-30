using Documenter, PiecewiseDeterministicMarkovProcesses

makedocs(doctest = false,
			sitename = "PiecewiseDeterministicMarkovProcesses.jl",
			pages = Any[
				"Home" => "index.md",
				"Tutorials" => "tutorials.md",
				"Problem Type" => "problem.md",
				"Solver Algorithms" => "solver.md",
				"Additional Features" => "addFeatures.md"
			]
)

deploydocs(
	repo   = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git",
)
