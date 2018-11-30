using Documenter, PiecewiseDeterministicMarkovProcesses

makedocs(
	# format = :html,
	sitename = "Piecewise Deterministic Markov Processes in Julia "
	)

# ENV["DOCUMENTER_DEBUG"] = true
# ENV["TRAVIS_REPO_SLUG"] = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git"

deploydocs(
	deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-cinder", "pygments"),
	repo   = "github.com/rveltz/PiecewiseDeterministicMarkovProcesses.jl.git",
	julia  = "1.0",
	osname = "linux",
)

