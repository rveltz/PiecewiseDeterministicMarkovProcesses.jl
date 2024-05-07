using Pkg
cd(@__DIR__)
pkg" activate ."
pkg" dev LSODA Sundials Plots"

using Documenter, PiecewiseDeterministicMarkovProcesses

makedocs(modules = [PiecewiseDeterministicMarkovProcesses],
		authors = "Romain Veltz",
		doctest = false,
		pagesonly = true, # this is on Documenter#master, do not compile what is not in pages =
# 		draft = true,
		warnonly = true,
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
	push_preview=true, target="build", devbranch="master")
