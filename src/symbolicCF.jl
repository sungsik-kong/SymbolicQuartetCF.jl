"""
    network_expectedCF_formulas(net::HybridNetwork;
        showprogressbar=false,
        inheritancecorrelation=0,
        symbolic=false)

This is the main function of **symCF**. This function iterates over *all* quartets of taxa in `net` and
generates the symbolic and numerical formulas that compute **expected quartet concordance factors (CFs)** 
under the **network multispecies coalescent**.

# Mandatory argument
- `net::HybridNetwork`  
  A rooted/unrooted phylogenetic network (parsed in *HybridNetwork* object) with all parameters (edge lengths and inheritance probabilities).

# Keyword Arguments
- `showprogressbar::Bool=false`  
  Whether to print progress while computing CFs.

- `inheritancecorrelation::Float64=0`  
  The inheritance correlation parameter ρ. Must be in `[0,1]`.

- `symbolic::Bool=false`  
  Returns numerical expressions for each CF. If `true`, returns symbolic expressions.  
"""
function network_expectedCF_formulas(net::HybridNetwork; 
    showprogressbar=false, 
    inheritancecorrelation=0, 
    symbolic=false::Bool)
    
    # ESA -- possible fix to network edge number problem.  This requires the user to reset the edge numbers, and
    # therefore we do not need to add a ! to the function name.
    
    # Dec 17: commenting out line below for now, to test code
    # sort([e.number for e in net.edge])[end] == length(net.edge) || error("The edges are not numbered consecutively from 1 to $(length(net.edge)).  Please run PhyloNetworks.resetedgenumbers!() first.  Exiting.")
    #
    # This would need to added to make_edge_label and elsewhere too.

    # data frame for CFs and dictionary translation
    df = DataFrame(Split=String[], CF=String[]) 
    synth_e_dict = Dict()
    
    # Dictionary for inheritance probabilities (γ)
    hybridNodeNumbers = [n.number for n in net.node if n.hybrid]
    for (j, hybNode) in enumerate(hybridNodeNumbers)
        incoming = [e for e in net.edge if PhyloNetworks.getchild(e).number == hybNode]
        if length(incoming) != 2
            error("Hybrid node $hybNode has $(length(incoming)) incoming edges (expected 2).")
        end
        for (k, e) in enumerate(incoming)
            e.gamma = round(e.gamma, digits=dpoints)
            synth_e_dict[(PN.getparent(e).number, PN.getchild(e).number, e.gamma)] = k == 1 ? "$gLab{$j}" : "(1-$gLab{$j})"            
        end
    end
    
    # Inheritance correlation (ρ)
    synth_e_dict[(NaN,NaN,inheritancecorrelation)] = "rho"
    synth_e_dict[(NaN,NaN,round(1 - inheritancecorrelation, digits=dpoints))] = "1-rho"
    
    for e in net.edge
        synth_e_dict[(e.number,PN.getparent(e).number, PN.getchild(e).number)]="1" * repeat("0", e.number-1)
    end
    
    net.node[net.rooti].leaf && error("The root can't be a leaf.")
    PN.check_nonmissing_nonnegative_edgelengths(net,"Edge lengths are needed in coalescent units to calcualte expected CFs.")
    all(e.gamma >= 0.0 for e in net.edge) || error("some γ's are missing for hybrid edges: can't calculate expected CFs.")
    inheritancecorrelation >= 0 || error("the inheritance correlation should be non-negative")
    inheritancecorrelation <= 1 || error("the inheritance correlation should be <= 1")
    taxa = sort!(tiplabels(net))
    taxonnumber = Dict(taxa[i] => i for i in eachindex(taxa))
    ntax = length(taxa)
    qtype = MVector{3,Float64} # 3 floats: CF12_34, CF13_24, CF14_23; initialized at 0.0
    numq=binomial(ntax,4)
    quartet = Vector{PN.QuartetT{qtype}}(undef, numq)
    ts = [1,2,3,4]
    # get taxon sets
    for qi in 1:numq
        quartet[qi] = PN.QuartetT(qi, SVector{4}(ts), MVector(0.,0.,0.))
        # next: find the 4-taxon set with the next rank,
        #       faster than using the direct mapping function
        ind = findfirst(x -> x>1, diff(ts))
        if ind === nothing ind = 4; end
        ts[ind] += 1
        for j in 1:(ind-1)
            ts[j] = j
        end
    end
    if showprogressbar
        nstars = (numq < 50 ? numq : 50)
        nquarnets_perstar = (numq/nstars)
        println("Calculation quartet CFs for $numq quartets...")
        print("0+" * "-"^nstars * "+100%\n  ")
        stars = 0
        nextstar = Integer(ceil(nquarnets_perstar))
    end
    for qi in 1:numq
        
        network_expectedCF!(quartet[qi], net, taxa, taxonnumber, inheritancecorrelation, df, symbolic, synth_e_dict)
        if showprogressbar && qi >= nextstar
            print("*")
            stars += 1
            nextstar = Integer(ceil((stars+1) * nquarnets_perstar))
        end
    end
    showprogressbar && print("\n")
    
    return quartet, taxa, df
end

"""
    network_expectedCF!(quartet::QuartetT, net::HybridNetwork, taxa, taxonnumber,
            inheritancecorrelation)

Update `quartet.data` to contain the quartet concordance factors expected from
the multispecies coalescent along network `net` for the 4-taxon set `taxa[quartet.taxonnumber]`.
`taxa` should contain the tip labels in `net`. `quartet.taxonnumber` gives the
indices in `taxa` of the 4 taxa of interest. `taxonnumber` should be a dictionary
mapping taxon labels in to their indices in `taxa`, for easier lookup.

`net` is not modified.

For `inheritancecorrelation` see [`network_expectedCF`](@ref).
Its value should be between 0 and 1 (not checked by this internal function).
"""
function network_expectedCF!(quartet::PN.QuartetT{MVector{3,Float64}},
    net::HybridNetwork, 
    taxa, 
    taxonnumber, 
    inheritancecorrelation, 
    df,
    symbolic,
    synth_e_dict)
    # create an array that stores the CF formulas for ab|cd, ac|bd, ad|bc
    qCFp=String["","",""] 
    
    net = deepcopy(net)
    symCF_removedegree2nodes!(net, synth_e_dict)
    
    # delete all taxa except for the 4 in the quartet
    for taxon in taxa
        taxonnumber[taxon] in quartet.taxonnumber && continue
        Qdeleteleaf!(net, taxon, synth_e_dict, simplify=false, unroot=false) # would like unroot=true but deleteleaf! throws an error when the root is connected to 2 outgoing hybrid edges
    end
    
    q,qCFp=network_expectedCF_4taxa!(net, taxa[quartet.taxonnumber], inheritancecorrelation, qCFp, symbolic, synth_e_dict)
    quartet.data .= q
    
    # storing the equations to DataFrames
    qCFp = [filter(!isspace,x) for x in qCFp[1:3]]
    f = taxa[quartet.taxonnumber]
    pairs = [(1,2,3,4), (1,3,2,4), (1,4,2,3)]
    for (i,(a,b,c,d)) in enumerate(pairs)
        push!(df, ("$(f[a])$(f[b])|$(f[c])$(f[d])", string(qCFp[i])))
    end
    
    return quartet
end

"""
    network_expectedCF_4taxa!(net::HybridNetwork, fourtaxa, inheritancecorrelation)

Return the quartet concordance factors expected from the multispecies coalescent
along network `net`, where the 3 quartet topologies are ordered following the
ordering of taxon names in `fourtaxa`, that is: if `fourtaxa` is a,b,c,d,
then the concordance factors are listed in this order:

    (qCF(ab|cd), qCF(ac|bd), qCF(ad,bc))

Assumptions about `net`:
- has 4 taxa, and those are the same as `fourtaxa`
- no degree-2 nodes, except perhaps for the root
- edge lengths are non-missing
- hybrid edge γ's are non-missing

The network is modified as follows: what's above the LSA is removed,
the 2 edges incident to the root are fused (if the root is of degree 2),
and external degree-2 blobs are removed. `net` is then simplified recursively
by removing hybrid edges for the recursive calculation of qCFs.

For `inheritancecorrelation` see [`network_expectedCF`](@ref).
Its value should be between 0 and 1 (not checked by this internal function).
"""
function network_expectedCF_4taxa!(net::HybridNetwork, fourtaxa, inheritancecorrelation, qCFp, symbolic,synth_e_dict)
    #((((b:1.0,(a:1.0,e:1.0):1.0):1.0,(((c:1.0,d:1.0):1.0)#H5:1.0::0.5)#H4:1.0::0.5):1.0,#H4:1.0::0.5):1.0,#H5:1.0::0.5);
    #println(PhyloNetworks.writenewick(net))
    qCFp .*= "("         
    deleteaboveLSA!(net)
    # make sure the root is of degree 3+
    if length(net.node[net.rooti].edge) <= 2
        Qfuseedgesat!(net.rooti, net,synth_e_dict)
    end
    # find and delete degree-2 blobs along external edges
    bcc = biconnectedcomponents(net, true) # true: ignore trivial blobs
    entry = PN.biconnectedcomponent_entrynodes(net, bcc)
    # entryindex = indexin(entry, net.nodes_changed)
    entryindex = indexin(entry, net.vec_node)
    exitnodes = PN.biconnectedcomponent_exitnodes(net, bcc, false) # don't redo the preordering
    bloborder = sortperm(entryindex) # pre-ordering for blobs in their own blob tree
    
    function isexternal(ib) # is bcc[ib] of degree 2 and adjacent to an external edge?
        # yes if: 1 single exit adjacent to a leaf
        length(exitnodes[ib]) != 1 && return false
        ch = PN.getchildren(exitnodes[ib][1])
        return length(ch) == 1 && ch[1].leaf
    end
    
    for ib in reverse(bloborder)  
        isexternal(ib) || continue # keep bcc[ib] if not external of degree 2
        for he in bcc[ib]
            he.ismajor && continue
            # deletion of a hybrid can hide the deletion of another: check that he is still in net
            any(e -> e===he, net.edge) || continue
            # delete minor hybrid edge with options unroot=true: to make sure the
            # root remains of degree 3+, in case a degree-2 blob starts at the root
            # simplify=true: bc external blob
            net=Qdeletehybridedge!(net,he,synth_e_dict,false,true,false,true,false)
        end
    end
    ndes = 4 # number of taxa descendant from lowest hybrid node
    if net.numhybrids > 0
        preorder!(net)
        # find a lowest hybrid node and # of taxa below it
        # hyb = net.nodes_changed[findlast(n -> n.hybrid, net.nodes_changed)]
        hyb = net.vec_node[findlast(n -> n.hybrid, net.vec_node)]
        funneledge = [e for e in hyb.edge if PhyloNetworks.getparent(e) === hyb]
        ispolytomy = length(funneledge) > 1
        funneldescendants = union([PN.descendants(e) for e in funneledge]...)
        ndes = length(funneldescendants)
        n2 = (ispolytomy ? hyb : PhyloNetworks.getchild(funneledge[1]))
        ndes > 2 && n2.leaf && error("2+ descendants below the lowest hybrid, yet n2 is a leaf. taxa: $(fourtaxa)")
    end
    if ndes > 2 # simple formula for qCF: find cut edge and its length
        # inheritance correlation has no impact
        # pool of cut edges below. contains NO external edge, bc n2 not leaf (if reticulation), nice tree ow
        cutpool = (net.numhybrids == 0 ? net.edge :
        [e for e in n2.edge if PN.getparent(e) === n2])
        filter!(e -> !PhyloNetworks.getchild(e).leaf, cutpool)
        net.numhybrids > 0 || length(cutpool) <= 1 ||
        error("2+ cut edges, yet 4-taxon tree, degree-3 root and no degree-2 nodes. taxa: $(fourtaxa)")
        sistertofirst = 2    # arbitrarily correct if 3-way polytomy (no cut edge)
        internallength = 0.0 # correct if polytomy
        cut_edges_pa=0
        cut_edges_ch=0
        enumber=1
        for e in cutpool
            length(cutpool) < 3 || println("more than 2 edged merged")
            internallength += e.length
            enumber=e.number
            cut_edges_pa=PN.getparent(e).number
            cut_edges_ch=PN.getchild(e).number
            hwc = hardwiredcluster(e, fourtaxa)
            sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
        end
        minorcf = exp(-internallength)/3
        majorcf = 1.0 - 2 * minorcf
        qCF = (sistertofirst == 2 ? MVector{3,Float64}(majorcf,minorcf,minorcf) :
        (sistertofirst == 3 ? MVector{3,Float64}(minorcf,majorcf,minorcf) :
        MVector{3,Float64}(minorcf,minorcf,majorcf) ))
        #writing out the equations
        if symbolic
            minorcfp = "exp(-$(binary_to_tstring(synth_e_dict[(enumber,cut_edges_pa,cut_edges_ch)])))/3"
            majorcfp = "1-2*$minorcfp"
        else
            minorcfp = "exp(-$internallength)/3" 
            majorcfp = "1-2*$minorcfp"
        end
        #=to be added to symbolic=#
        (sistertofirst == 2 ? (qCFp[1]*="$majorcfp",qCFp[2]*="$minorcfp",qCFp[3]*="$minorcfp") :
        (sistertofirst == 3 ? (qCFp[1]*="$minorcfp",qCFp[2]*="$majorcfp",qCFp[3]*="$minorcfp") :
        (qCFp[1]*="$minorcfp",qCFp[2]*="$minorcfp",qCFp[3]*="$majorcfp") ))                      
        qCFp .*= ")" #kong: end qCF with an closing bracket      
        return qCF, qCFp
    end
    ndes > 0 || error("weird: hybrid node has no descendant taxa")
    # by now, there are 1 or 2 taxa below the lowest hybrid
    qCF = MVector{3,Float64}(0.0,0.0,0.0) # mutated later
    #parenthedge = [e for e in hyb.edge if getchild(e) === hyb]
        parenthedge = [e for e in hyb.edge if PhyloNetworks.getchild(e) === hyb]
        all(h.hybrid for h in parenthedge) || error("hybrid $(hyb.number) has a parent edge that's a tree edge")
        parenthnumber = [p.number for p in parenthedge]
        nhe = length(parenthedge)
        if ndes == 1 # weighted qCFs average of the nhe (often = 2) displayed networks
            # inheritance correlation has no impact
            for i in 1:nhe # keep parenthedge[i], remove all others
                # kong: add + if there are more than a single element
                if i>1 qCFp .*= "+" end
                gamma = parenthedge[i].gamma
                ei=parenthedge[i]
                simplernet = ( i < nhe ? deepcopy(net) : net ) # last case: to save memory allocation
                for j in 1:nhe
                    j == i && continue # don't delete hybrid edge i!
                    pe_index = findfirst(e -> e.number == parenthnumber[j], simplernet.edge)
                    simplernet=Qdeletehybridedge!(simplernet, simplernet.edge[pe_index],synth_e_dict, false,true,false,false,false) # ., unroot=true, ., simplify=false,.          
                end
                
                qCF0,qCFp = network_expectedCF_4taxa!(simplernet, fourtaxa, inheritancecorrelation,qCFp,symbolic,synth_e_dict)
                
                qCF .+= gamma .* qCF0
                
                if (symbolic) qCFp .*= "*$(synth_e_dict[(PN.getparent(ei).number, PN.getchild(ei).number, ei.gamma)])" 
                else qCFp .*= "*$gamma" end
                
            end
            qCFp .*= ")"
            return qCF, qCFp
        end
        # by now: 2 descendant below the lowest hybrid node: hardest case
        # weighted qCFs average of 3 networks: 2 displayed, 1 "parental" (unless same parents)
        sameparents = (inheritancecorrelation == 1)
        oneminusrho = 1 - inheritancecorrelation
        if symbolic 
            oneminusrhop="(1-$(synth_e_dict[(NaN,NaN,inheritancecorrelation)]))"
        else
            oneminusrhop="(1-$inheritancecorrelation)"
        end
        hwc = hardwiredcluster(parenthedge[1], fourtaxa)
        sistertofirst = findnext(x -> x == hwc[1], hwc, 2)
        internallength = ( ispolytomy ? 0.0 : funneledge[1].length)
        deepcoalprob = exp(-internallength)
        if symbolic deepcoalprobp="exp(-$internallength)" end
        # initialize qCF: when the 2 descendants coalesce before reaching the hybrid node
        qCF = (sistertofirst == 2 ? MVector{3,Float64}(1.0-deepcoalprob,0.0,0.0) :
        (sistertofirst == 3 ? MVector{3,Float64}(0.0,1.0-deepcoalprob,0.0) :
        MVector{3,Float64}(0.0,0.0,1.0-deepcoalprob) ))
        if symbolic 
            if iszero(internallength)
                deepcoalprobp = "exp(-0.0)" 
            else
                deepcoalprobp = "exp(-$(binary_to_tstring(synth_e_dict[(funneledge[1].number,PN.getparent(funneledge[1]).number,PN.getchild(funneledge[1]).number)])))" 
            end
        else 
            deepcoalprobp = "exp(-$internallength)" 
        end
        (sistertofirst == 2 ? (qCFp[1]*="(1-$deepcoalprobp)+",qCFp[2]*="",qCFp[3]*="") :
        (sistertofirst == 3 ? (qCFp[1]*="",qCFp[2]*="(1-$deepcoalprobp)+",qCFp[3]*="") :
        (qCFp[1]*="",qCFp[2]*="",qCFp[3]*="(1-$deepcoalprobp)+") ))                                     
        qCFp .*= ""
        # no coalescence on cut-edge: delete it and extract parental networks
        ispolytomy || PN.shrinkedge!(net, funneledge[1])
        # shrinkedge! requires PhyloNetworks v0.15.2
        childedge = [e for e in hyb.edge if PN.getparent(e) === hyb]
        length(childedge) == 2 || error("2-taxon subtree, but not 2 child edges after shrinking the cut edge.")
        all(PN.getchild(e).leaf for e in childedge) || error("2-taxon subtree, cut-edge shrunk, but the 2 edges aren't both external")
        childnumber = [e.number for e in childedge]
        for i in 1:nhe
            pgam=parenthedge[i].gamma
            ei=parenthedge[i]
            pgamtuple=(PN.getparent(ei).number,PN.getchild(ei).number,pgam)
            weighti = deepcoalprob * pgam      
            for j in (sameparents ? i : 1):i # if inheritancecorrelation=1 then i!=j has probability 0
                gammaj = parenthedge[j].gamma
                ei=parenthedge[j]
                gammajtuple=(PN.getparent(ei).number,PN.getchild(ei).number,gammaj)
                
                simplernet = ( i < nhe || j < nhe ? deepcopy(net) : net )
                # delete all hybedges other than i & j
                for k in 1:nhe
                    (k == i || k ==j) && continue # don't delete hybrid edges i or j
                    pe_index = findfirst(e -> e.number == parenthnumber[k], simplernet.edge)
                    simplernet=Qdeletehybridedge!(simplernet, simplernet.edge[pe_index],synth_e_dict,false,true,false,false,false) # ., unroot=true,., simplify=false,.
                end
                if i != j
                    # detach childedge[2] from hyb and attach it to hyb's parent j
                    pej_index = findfirst(e -> e.number == parenthnumber[j], simplernet.edge)
                    pej = simplernet.edge[pej_index]
                    pn = PN.getparent(pej)
                    hn = PN.getchild(pej) # hyb node, but in simplernet
                    ce2_index = findfirst(e -> e.number == childnumber[2], simplernet.edge)
                    ce2 = simplernet.edge[ce2_index]
                    PN.removeEdge!(hn,ce2)
                    hn_index = findfirst(x -> x === hn, ce2.node)
                    ce2.node[hn_index] = pn # ce2.ischild1 remains synchronized
                    push!(pn.edge, ce2)
                    # then delete hybedge j
                    simplernet=Qdeletehybridedge!(simplernet, pej,synth_e_dict, false,true,false,false,false) # ., unroot=true,., simplify=false,.)
                end
                #kong: add "+" if there are more than a single element in the equation
                if i>1 qCFp .*= "+" end
                #initialize qCFp for qCF_subnet
                qCFps=["","",""]
                qCF_subnet, qCFps = network_expectedCF_4taxa!(simplernet, fourtaxa, inheritancecorrelation, qCFps,symbolic,synth_e_dict)
                if i == j
                    prob = weighti * (gammaj * oneminusrho + inheritancecorrelation)
                    qCF .+= prob .* qCF_subnet
                    if symbolic
                        for i in 1:3
                            qCFp[i] *= "((($deepcoalprobp * $(synth_e_dict[pgamtuple])) * ($(synth_e_dict[gammajtuple]) * ($(synth_e_dict[(NaN,NaN,oneminusrho)])) + $(synth_e_dict[(NaN,NaN,inheritancecorrelation)]))) * ($(qCFps[i])))"
                        end
                    else
                        for i in 1:3
                            qCFp[i] *= "((($deepcoalprobp * $pgam) * ($gammaj * ($oneminusrhop) + $inheritancecorrelation)) * ($(qCFps[i])))"
                        end
                    end
                else # add subnetwork with flipped assignment of the 2 taxa to parents i & j
                    flipped_ij = (sistertofirst == 2 ? [1,3,2] :
                    (sistertofirst == 3 ? [3,2,1] : [2,1,3] ))
                    prob = weighti * gammaj * oneminusrho
                    qCF .+= prob .* (qCF_subnet .+ qCF_subnet[flipped_ij])
                    if symbolic
                        for i in 1:3
                            qCFp[i] *= "((($deepcoalprobp * $(synth_e_dict[pgamtuple])) * $(synth_e_dict[gammajtuple]) * ($(synth_e_dict[(NaN,NaN,oneminusrho)]))) * ($(qCFps[i])+$(qCFps[flipped_ij[i]])))"
                        end
                    else
                        for i in 1:3
                            qCFp[i] *= "((($deepcoalprobp * $pgam) * $gammaj * $oneminusrhop) * ($(qCFps[i])+$(qCFps[flipped_ij[i]])))"
                        end
                    end
                end
            end
        end
        qCFp .*= ")"
        return qCF, qCFp
    end
