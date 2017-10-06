using Documenter, LSODA

makedocs(
    modules=[LSODA],
    doctest = false
)

deploydocs(
    deps = Deps.pip("mkdocs","python-markdown-math"),
    repo   = "github.com/rveltz/LSODA.jl.git",
    julia  = "release"
)
