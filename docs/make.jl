using Documenter
using SymbolicQuartetCF

makedocs(
    sitename = "SymbolicQuartetCF.jl",
    modules  = [SymbolicQuartetCF],
    pages = [
        "Home" => "index.md",
        "API"  => "api.md",
    ],
)

deploydocs(
    repo = "github.com/sungsik-kong/SymbolicQuartetCF.jl.git",
)