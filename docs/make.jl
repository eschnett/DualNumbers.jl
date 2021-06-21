# Generate documentation with this command:
# (cd docs && julia --color=yes make.jl)

push!(LOAD_PATH, "..")

using Documenter
using DualNumbers

makedocs(; sitename="DualNumbers", format=Documenter.HTML(), modules=[DualNumbers])

deploydocs(; repo="github.com/eschnett/DualNumbers.jl.git", devbranch="main", push_preview=true)
