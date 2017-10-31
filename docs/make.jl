using Documenter, PDMP

makedocs()

deploydocs(
    repo   = "github.com/rveltz/PDMP.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)

