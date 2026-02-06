module SymbolicQuartetCF

    using CSV
    using DataFrames
    using PhyloNetworks
    #using PhyloPlots
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
    read_topology_rand,
    network_expectedCF_formulas, 
    make_edge_label,
    clean_labels,
    reindex_edges


    include("misc.jl")
    include("symbolicCF.jl")
    include("formatting.jl")    

end # module
