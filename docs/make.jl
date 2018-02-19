using Documenter, PDMP

makedocs(format = :html,sitename = "Piecewise Deterministic Markov Processes in Julia ")

ENV["DOCUMENTER_DEBUG"] = true

deploydocs(
	deps   = Deps.pip("mkdocs", "python-markdown-math", "mkdocs-cinder"),
	repo   = "github.com/rveltz/PDMP.jl.git",
	julia  = "0.6",
	osname = "linux",
	make = nothing
)

