using Documenter, PiecewiseDeterministicMarkovProcesses

makedocs(doctest = false,
			sitename = "PiecewiseDeterministicMarkovProcesses.jl",
			pages = Any[
				"Home" => "index.md",
				"Tutorials" => "tutorials.md",
				"Problem Type" => "problem.md",
				"Solver Algorithms" => "solver.md",
				"FAQ" => "addFeatures.md",
				"Library" => "library.md"
			]
)

deploydocs(
	repo   = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git",
)
