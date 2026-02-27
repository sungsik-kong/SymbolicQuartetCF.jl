---
title: 'SymbolicQuartetCF.jl: A Julia package for symbolic computation of quartet concordance factors'
tags:
- Julia
- Phylogenetics
- Phylogenetic networks
- Network multispecies coalescent
- Identifiability
date: "28 February 2026"
output:
  html_document:
    df_print: paged
  pdf_document: default
authors:
- name: Sungsik Kong
  corresponding: yes
  orcid: "0000-0002-5882-2109"
  affiliation: 1, 2, 3
- name: Aviva K. Englander
  orcid: "0000-0002-2074-7684"
  affiliation: 2, 3
- name: Jennifer Garbett
  orcid: "0000-0003-0512-7003"
  affiliation: 4
- name: Ikenna Nometa
  orcid: "0000-0001-6461-8158"
  affiliation: 2, 5
- name: Elizabeth Allman
  orcid: "0000-0002-4866-9606"
  affiliation: 6
bibliography: paper.bib
aas-doi: #"10.3847/xxxxx <- update this with the DOI from AAS once you know it."
aas-journal: #"Astrophysical Journal <- The name of the AAS journal." 
affiliations:
- name: "The RIKEN Center for Interdisciplinary Theoretical and Mathematical Sciences, 2-1 Hirosawa, Wako, Saitama 351-0198, Japan"
  index: 1
- name: "The Institute for Computational and Experimental Research in Mathematics, 121 South Main Street, Box E, 11th Floor, Providence, RI 02903, USA"
  index: 2
- name: "Wisconsin Institute for Discovery, University of Wisconsin-Madison, 330 N. Orchard Street, Madison, WI 53715, USA"
  index: 3
- name: "Lenoir-Rhyne University, 625 7th Ave NE, Hickory, NC 28601, USA"
  index: 4
- name: "Georgia Institute of Technology, 225 North Ave, Atlanta, GA 30332, USA"
  index: 5
- name: "University of Alaska Fairbanks, 1792 Ambler Lane, Fairbanks, AK 99775, USA"
  index: 6
---

### Summary

Phylogenies are graphs depicting evolutionary relationships among biological entities, such as species or genes, typically called *taxa*.  Perhaps the most famous historical phylogeny is Darwin's sketch of a tree relating Galapagos finches [@Darwin1837Ithink]. Phylogenetic trees and networks (accommodating non-treelike biological processes [see @KongPonsKubatkoWicke2022]) are now routinely constructed from molecular data.

*Multilocus data*, thousands of gene trees (or locus trees) constructed from molecular data like genomic sequences, are often used to relate extant species with a tree or network.  These data typically show tremendous topological discordance, indicating that the evolutionary history of a stretch of the genome is different from the true species-level evolutionary history. A null model for describing this gene tree heterogeneity is the **Network Multispecies Coalescent Model** (NMSC).  By assuming the NMSC has generated a gene tree sample, many network inference methods and algorithms have been introduced [@PhyloNet_2009; @InferNetwork_MPL_software_2015; @solis-lemus2016; @PhyloNet_SystBiol_2018_network_inference;  @allman2019; @Ultrafast; @Allman2025_NANUQPlus].

To ease computation time, many NMSC-based inference methods are based on **quartets**, obtained from a large species network by focusing only on subsets of four taxa at a time. For each quartet $Q$, the frequencies at which particular 4-taxon tree topologies are observed on gene trees in a multilocus dataset are called empirical **quartet concordance factors** (qCFs), and are estimates of the expected concordance factors under the NMSC.  These data summaries are the input for several of the network inference methods, and require closed-form expressions of expected qCFs for algorithm development.

We introduce [`SymbolicQuartetCF.jl`](https://github.com/sungsik-kong/SymbolicQuartetCF.jl) ([`symCF`]() hereafter), a Julia package that computes symbolic expressions for expected quartet concordance factors under the NMSC for any model network $N$, building on some functionality of [`QuartetNetworkGoodnessFit.jl`](https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/) [@QuartetNetworkGoodnessFit.jl]. [`symCF`]() optionally saves symbolic formulas to a CSV file for downstream analyses, and can generate a plot of the model network $N$ with edges labeled by dynamically-assigned variable names used when generating the formulas.

Since formulas for qCFs have been used for proving numerous parameter identifiability results under the NMSC, including methods using computational algebraic geometry [@allman2024a; @allman2025b; @rhodes2025; @ExtractingDiamonds2025], output from [`symCF`]() can be formatted for use with `Singular`, `Macaulay2`, or `MATLAB` [@Singular2024; @M2; @matlab_2024b]. In particular, Singular and Macaulay2 files that define parameterized quartet concordance factor varieties corresponding to a model network $N$ under the NMSC can be written.  If the `MATLAB` option is chosen, a MATLAB script file that can compute the dimension of the variety $V$, either using a symbolic toolbox or numerically, is [written](https://github.com/sungsik-kong/SymbolicQuartetCF.jl/wiki/Usage) out.
Figure \autoref{fig:pipeline} illustrates the main pipeline:  model network $N \rightarrow$ quartet subgraphs $\rightarrow$ computation of $3 \binom{n}{4}$ symbolic qCFs.
 
The [`symCF`]() package can aid in proving novel parameter identifiability theorems under the NMSC for 
complex phylogenetic model networks, help in the development of new network inference methods
under the NMSC, and further strengthen the interdisciplinary bond 
between mathematical and computational phylogenetics.  

![Figure 1. Graphical workflow in [`symCF`](). The input species network topology on five taxa, $A,B,C,D,$ and $E$ (left), is decomposed into a set of quartets each of which contains four taxa (center), and [`symCF`]() generates symbolic formulas for the 3 expected qCFs for each of the $5\choose 4$ quartets (right). \label{fig:pipeline}](Figures/Figure-pipeline.png)
<!--{width=80% style="display:block; margin-left:auto; margin-right:auto;"}-->

### Statement of need

Inference of phylogenetic networks from genomic data is a timely and important problem in evolutionary biology. The enormous size of genome-scale datasets necessitates using data summaries like quartets and quartet methods to infer a network or tree.  Network inference is a nascent and growing field, and statistically sound methods are much needed. The NMSC is often taken as the null model for explaining the discordance seen in gene tree datasets, but what features of trees and networks are recoverable from gene tree data is not solved, though progress [@allman2024a] has been made for level-1 networks, those with cycle reticulations, under the NMSC, and some promising results for somewhat more complicated networks [@HoltgrefeLevel2025;@allman2025b].  Establishing parameter identifiability (and importantly, its failure) for more complex classes of networks from qCFs justifies that network methods from these data will be sound and consistent.  Such research investigations can not be performed without symbolic formulas for qCFs.

Quartet (and triplet) data derived from both sequence and gene tree data are used in composite-likelihood frameworks in existing NMSC-network inference packages [@InferNetwork_MPL_software_2015; @solis-lemus2016; @Kong2025phynest]. Many more inference frameworks will be developed in upcoming years, and tools like [`symCF`]() that aid in the testing of novel paradigms and proving identifiability theorems are essential.

Prior to [`symCF`](), whether needed for algorithm development or mathematical theorem proving, exact formulas for expected qCFs were derived by hand, a tedious and error-prone task.  The package [`symCF`]() provides several advances: First, in automatically generating symbolic qCF formulas, [`symCF`]() eliminates errors and accelerates the process. Second, such automation expands the complexity of network classes amenable to study under the NMSC since qCFs have complicated formulas that grow in complexity with that of $N$.  Moreover, provably sound inference under variants of the NMSC might become possible [*e.g.*, see @Fogg2023PhyloCoalSimulations]. Third, by outputting formulas in formats compatible with computer algebra packages Macaulay2, Singular, and MATLAB, [`symCF`]() provides researchers in algebraic geometry and algebraic statistics with a starting point for in-depth geometric investigations.

### Features

The following sections correspond to the implementations in [`symCF`]() version 0.1.7.
[`symCF`]() generates numeric and symbolic qCF formulas by extending an algorithm implemented in [`QuartetNetworkGoodnessFit.jl`](https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/) [@QuartetNetworkGoodnessFit.jl] that computes numerical expected qCFs from a species network under the NMSC. The main function of [`symCF`]() is `network_expectedCF_formulas` where the phylogenetic network $N$, a [`HybridNetwork`](https://juliaphylo.github.io/PhyloNetworks.jl/stable/lib/public/#PhyloNetworks.HybridNetwork) object, is the only mandatory argument. The input network should be a `HybridNetwork` object that has all numerical parameters like edge lengths in coalescent units and inheritance probabilities defined. Using a topological network (i.e., purely topological or with incomplete parameter assignment) encoded in Newick string format, [`symCF`]()’s function `read_topology_rand` can be used to read in the string and arbitrarily assign all numerical parameters, returning a `HybridNetwork` object for use with `network_expectedCF_formulas`.  See example below. For Newick string with numerical parameters, the `readnewick` function from the `Julia` package [`PhyloNetworks`](https://github.com/juliaphylo/PhyloNetworks.jl) ($\ge$ v1.2) creates a `HybridNetwork` object.

```
julia> using SymbolicQuartetCF

julia> net=read_topology_rand("((((b,(a,e)),(((c,d))#H5)#H4),#H4),#H5);")
PhyloNetworks.HybridNetwork, Rooted Network
14 edges
13 nodes: 5 tips, 2 hybrid nodes, 6 internal tree nodes.
tip labels: b, a, e, c, ...
((((b:1.253,(a:1.304,e:1.227):1.306):1.794,(((c:1.995,d:1.668):1.188)#H5:1.429::0.871)#H4:1.487::0.871):1.213,#H4:1.932::0.129):1.242,#H5:1.838::0.129);
```


The function `network_expectedCF_formulas` returns three values: 

1. the expected qCFs, 
2. a list of taxon names, and 
3. a `DataFrame` object with qCF formulas for all quartets. 

Note that values (1.) and (2.) are those returned by `network_expectedCF` from [`QuartetNetworkGoodnessFit.jl`](https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/). By default, the `DataFrame` returned (3.) contains a list of strings with unevaluated numerical formulas for the qCFs; that is, containing expressions like "exp(-1.2)".  By setting the argument `symbolic=true`, symbolic qCF formulas are returned.  These formulas can be computed with a correlated inheritance parameter $\rho$ [@Fogg2023PhyloCoalSimulations] by specifying a value for the optional argument `inheritancecorrelation` (`default=0`).

```
julia> cfs,taxa,formulas=network_expectedCF_formulas(net,symbolic=true)
(PhyloNetworks.QuartetT{StaticArraysCore.MVector{3, Float64}}[4-taxon set number 1; taxon numbers: 1,2,3,4
data: [0.9906628909602795, 0.004668554519860242, 0.004668554519860242], 4-taxon set number 2; taxon numbers: 1,2,3,5
data: [0.09026037037501333, 0.09026037037501333, 0.8194792592499733], 4-taxon set number 3; taxon numbers: 1,2,4,5
data: [0.09026037037501333, 0.09026037037501333, 0.8194792592499733], 4-taxon set number 4; taxon numbers: 1,3,4,5
data: [0.001264156380235584, 0.001264156380235584, 0.997471687239529], 4-taxon set number 5; taxon numbers: 2,3,4,5
data: [0.004668554519860242, 0.004668554519860242, 0.9906628909602795]], ["a", "b", "c", "d", "e"], 15×2 DataFrame
 Row │ Split   CF
     │ String  String
─────┼───────────────────────────────────────────
   1 │ ab|cd   ((1-exp(-t_{3}))+(((exp(-t_{3})*…
   2 │ ac|bd   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
   3 │ ad|bc   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
   4 │ ab|ce   ((exp(-t_{1})/3)*g_{1}+(exp(-t_{…
   5 │ ac|be   ((exp(-t_{1})/3)*g_{1}+(exp(-t_{…
   6 │ ae|bc   ((1-2*exp(-t_{1})/3)*g_{1}+(1-2*…
   7 │ ab|de   ((exp(-t_{1})/3)*g_{1}+(exp(-t_{…
   8 │ ad|be   ((exp(-t_{1})/3)*g_{1}+(exp(-t_{…
   9 │ ae|bd   ((1-2*exp(-t_{1})/3)*g_{1}+(1-2*…
  10 │ ac|de   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
  11 │ ad|ce   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
  12 │ ae|cd   ((1-exp(-t_{3}))+(((exp(-t_{3})*…
  13 │ bc|de   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
  14 │ bd|ce   ((((exp(-t_{3})*g_{1})*(g_{1}*(1…
  15 │ be|cd   ((1-exp(-t_{3}))+(((exp(-t_{3})*…)
  ```

The function `export_symbolic_format` converts symbolic formulas returned in `DataFrame` (named `formulas` above) for use in algebraic geometry applications using Macaulay2 (Gröbner basis computations, or using `MultigradedImplicitization`) or Singular. A third option creates a MATLAB script file that can compute (a lower bound for) the dimension of the parameterization map that defines the qCF variety.

Since symbolic edge length parameter names are generated internally, it can be difficult to determine which edge name (e.g., t4) corresponds to which particular edge in the model network. By combining the function `make_edge_label` and the Julia package [`PhyloPlots`](https://github.com/JuliaPhylo/PhyloPlots.jl), we visualize edge names directly on a network topology (Figure \autoref{fig:networkPlot}). Note that the terminal edges are not labeled here, as they do not play a role in qCF when a single lineage per taxon is sampled. To enhance readability, edge labels can be reindexed consecutively starting from one using the function `reindex_edges`. This function creates a copy of the original network with updated edge numbers, which can then be used for `make_edge_label` and plotting.

```
julia> reindexednet=reindex_edges(net)
julia> edgelabel=make_edge_label(reindexednet)
9×2 DataFrame
 Row │ number  label
     │ Int64   String
─────┼──────────────────────────
   1 │      1  t1
   2 │      2  t2
   3 │      3  t3
   4 │      9  t9, γ = g1
   5 │      8  t8, γ = g2
   6 │      4  t4
   7 │      7  t7, 1-γ = 1 - g2
   8 │      5  t5
   9 │      6  t6, 1-γ = 1 - g1

julia> using PhyloPlots
julia> plot(reindexednet,edgelabel=edgelabel)
```

![Figure 2. Visualizing symbolic parameter names on a topology using `make_edge_label` and `PhyloPlots`.](Figures/Figure-network-labelled.png)
<!--{width=50% style="display:block; margin-left:auto; margin-right:auto;" label="fig:networkPlot"}-->

### Author contributions

S.K. led the development of [`symCF`]() by extending `QuartetNetworkGoodnessFit.jl`. E.S.A. extended the functionalities of `symCF`, and A.K.E., J.G., and I.N. developed the computer algebra system applications.  All authors tested the code extensively and contributed to the preparation of this manuscript.  

### Acknowledgements

This material was supported in part by the National Science Foundation (NSF) under Grant No. DMS-1929284 to the Institute for Computational and Experimental Research in Mathematics in Providence, RI, where S.K., A.K.E., I.N., and E.S.A. were residents for the Theory, Methods, and Applications of Quantitative Phylogenomics semester program. This work was also supported by NSF DMS-2051760 to E.S.A., and DMS-194558 to Elizabeth Gross at the University of Hawai'i at Mānoa, where all authors participated in the workshop on Phylogenetic Networks. Additionally, throughout the course of this work, A.K.E. was supported by the Georgia Benkart Legacy Fund, the Wisconsin Alumni Research Foundation and NSF Award DMS-2023239. Finally, NSF DEB-2144367 to Claudia Solís-Lemus at the University of Wisconsin-Madison partially supported S.K.'s work as a post-doc.
We especially thank Claudia Solís-Lemus for suggesting to S.K. that he pursue this project.

### AI usage disclosure

No generative AI tools were used in the development of this software or in the writing of this manuscript, except for minor spell-checking.

### References
