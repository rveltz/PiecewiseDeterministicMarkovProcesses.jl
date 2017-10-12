using Documenter, PDMP

makedocs()

deploydocs(
    repo   = "github.com/sdwfrost/PDMP.jl.git",
    target = "build",
    deps   = nothing,
    make   = nothing
)

