using Documenter, PiecewiseDeterministicMarkovProcesses

makedocs(doctest = true,
			clean = true,
			format = :html,
			sitename = "PiecewiseDeterministicMarkovProcesses.jl"
			)

# ENV["DOCUMENTER_DEBUG"] = true
# ENV["TRAVIS_REPO_SLUG"] = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git"

deploydocs(
	# deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-cinder", "pygments"),
	repo   = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git",
	julia  = "1.0",
	osname = "linux",
    target = "build",
    deps   = nothing,
	make = nothing
)


