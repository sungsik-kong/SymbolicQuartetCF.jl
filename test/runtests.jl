#test written by Sungsik Kong 2025
using CSV
using DataFrames
using FileCmp
using PhyloNetworks
using SymbolicQuartetCF
using Test

@testset "SymbolicQuartetCF tests" begin
    include("test_symbolicQNC.jl")
end