module SymbolicQuartetCF

    using CSV
    using DataFrames
    using PhyloNetworks
    using StaticArrays
    import Random

    import PhyloNetworks: readnewick
    export readnewick

    const dpoints=10 # decimal points for all parameters when randomly generated
    const eLab="t_"
    const gLab="g_"
    const PN = PhyloNetworks

    export
    aloha,
    export_csv, 
    export_symbolic_format, 
    readTopologyrand,
    network_expectedCF_formulas, 
    makeEdgeLabel, 
    assignBinaryEdgeLengths

    include("misc.jl")
    include("symbolicCF.jl")
    include("formatting.jl")    

end # module
