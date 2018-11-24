========================================================================
    CRank: An approach for prioritizing network communities
========================================================================

This code implements a network community prioritization approach (CRank).
CRank is an automatic method for prioritizing network communities and
identifying the most promising ones for further experimentation. CRank
efficiently evaluates robustness and magnitude of structural features of
each community and then combines these features to obtain the community
prioritization. It can be used with any community detection method and
scales to large networks. It uses only information provided by the network
structure and does not requireany additional external metadata or labels.

Algorithm and the prioritization of model are described in the following
paper:
M. Zitnik, R. Sosic and J. Leskovec, Prioritizing Network Communities,
In review, 2017.

The code works under Windows with Visual Studio or Cygwin with GCC,
Mac OS X, Linux and other Unix variants with GCC. Make sure that a
C++ compiler is installed on the system. Visual Studio project files
and makefiles are provided. For makefiles, compile the code with
"make all".

Included are also binaries for Linux (crank_linux), Mac (crank_mac)
and Windows (crank_windows.exe). 

////////////////////////////////////////////////////////////////////////
Website:

Please check CRank's website for more information:

http://snap.stanford.edu/crank/index2.html


////////////////////////////////////////////////////////////////////////
Parameters:
   -c: Community affiliation data (default:'karate_communities.txt')
   -i: Input edgelist file name (default:'karate.txt')
   -o: Output file name (default:'prioritization.txt')
   -p: Fit an auxiliary network model (works with non-statistical community detection methods) (default: T)
   -in: Input file name (probabilities of nodes belonging to communities) (default: '')
   -ie: Input file name (probabilities of edges) (default: '')
   -ic: Input file name (probabilities of edges given communities) (default: '')
   -a: Network perturbation intensity alpha (default: 0.15)
   -pr: Relative size of temporary gold standard (p ~ prior probability) (default: 0.05)
   -b: Number of bins (-1: detect automatically, B = |C|/50) (default: -1)
   -mx: Maximum number of iterations for rank aggregation (default: 20)
   -s: Convergence criterion for rank aggregation (default: 1)


////////////////////////////////////////////////////////////////////////
Formatting of input files:
   -c: Each line represents a group and all members of the group are listed in
       a single line. The first entry in each line represents community name/id.
       If a node belongs to multiple groups, it is listed in multiple lines.
       For example, see: karate_communities.txt
   -i: Network edgelist. Each line represents one edge, given by its endpoints.
       For example, see: karate.txt
   -in: Each line represents a node-community affiliation and its
        probability. Each line has three entries: node id, community name/id,
        probability.
        For example, see: amazon_CProbaH.txt
   -ie: Each line represents an edge and its probability. Each
        line has three entries: node id, node id, probability.
        For example, see: amazon_EdgeProbaH.txt
   -ic: Each line represents an edge and its probability conditioned on
        nodes' joint affiliation with a community. Each line has four
        entries: node id, node id, community name/id, probability.
        For example, see: amazon_CEdgeProbaH.txt


////////////////////////////////////////////////////////////////////////
Format of output file:

Prioritization results are saved to a file specified by switch -o. Each
line represents a community, and has five entries: CRank prioritization
score, conductance score, modularity score, random score.

For example, see: karate_prioritization.txt


////////////////////////////////////////////////////////////////////////
Usage:

1) Prioritize 5 communities of Zachary's Karate club members (Cmt2 and Cmt4
represent two factions in the Karate club; Cmt1, Cmt3 and Cmt5 are
less meaningful groupings of club members):

./crank -i:karate.txt -c:karate_communities.txt -o:karate_prioritization.txt

2) Prioritize communities of the Amazon product co-purchasing network using
the probabilities returned by a statical community detection model:

./crank -i:amazon.txt -p:F -c:amazon_communities.txt -in:amazon_CProbaH.txt -ie:amazon_EdgeProbaH.txt -ic:amazon_CEdgeProbaH.txt -o:amazon_prioritization.txt

3) Prioritize communities of the PPI network using
the probabilities returned by a statical community detection model:

./crank -i:ppi2.txt -p:F -c:ppi2_communities.txt -in:ppi2_CProbaH.txt -ie:ppi2_EdgeProbaH.txt -ic:ppi2_CEdgeProbaH.txt -o:ppi2_prioritization.txt

4) Prioritize communities of the DBLP network using
the probabilities returned by a statical community detection model:

./crank -i:dblp.txt -p:F -c:dblp_communities.txt -in:dblp_CProbaH.txt -ie:dblp_EdgeProbaH.txt -ic:dblp_CEdgeProbaH.txt -o:dblp_prioritization.txt

5a) Prioritize communities of the medical drug network using
the probabilities returned by a statical community detection model:

./crank -i:drug_targets.txt -p:F -c:drug_targets_communities.txt -in:drug_targets_CProbaH.txt -ie:drug_targets_EdgeProbaH.txt -ic:drug_targets_CEdgeProbaH.txt -o:drug_targets_prioritization.txt -b:10

5b) Prioritize communities of the medical drug network using domain-specific information:

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_5.txt -o:drug_targets_prioritization_text_5.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_10.txt -o:drug_targets_prioritization_text_10.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_25.txt -o:drug_targets_prioritization_text_25.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_50.txt -o:drug_targets_prioritization_text_50.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_75.txt -o:drug_targets_prioritization_text_75.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_text_goldstandard_100.txt -o:drug_targets_prioritization_text_100.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_chemistry_goldstandard_5.txt -o:drug_targets_prioritization_chemistry_5.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_chemistry_goldstandard_10.txt -o:drug_targets_prioritization_chemistry_10.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_epistasis_goldstandard_5.txt -o:drug_targets_prioritization_epistasis_5.txt -b:10

./crank -i:drug_targets.txt -p:T -c:drug_targets_communities.txt -gs:drugs_epistasis_goldstandard_10.txt -o:drug_targets_prioritization_epistasis_10.txt -b:10

6) Prioritize communities of the synthetic planted community network using
the probabilities returned by a statical community detection model:

./crank -i:syn_30-30-30-30-30-30-30-30-30-30.txt -p:T -c:syn_30-30-30-30-30-30-30-30-30-30_communities.txt -in:syn_30-30-30-30-30-30-30-30-30-30_CProbaH.txt -ie:syn_30-30-30-30-30-30-30-30-30-30_EdgeProbaH.txt -ic:syn_30-30-30-30-30-30-30-30-30-30_CEdgeProbaH.txt -o:syn_30-30-30-30-30-30-30-30-30-30_prioritization.txt -b:3 -pr:0.2

./crank -i:syn_30-30-30-30-30-30-30-30-30-30--0.txt -p:T -c:syn_30-30-30-30-30-30-30-30-30-30--0_communities.txt -o:syn_30-30-30-30-30-30-30-30-30-30--0_prioritization.txt -b:3 -pr:0.2 -gs:syn_30-30-30-30-30-30-30-30-30-30--0_goldstandard.txt 

7) Prioritize communities of the synthetic planted community network using
the probabilities returned by a statical community detection model:

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K5.txt -gs:syn_groundtruth_goldstandard_5.txt

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K10.txt -gs:syn_groundtruth_goldstandard_10.txt -b:20

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K25.txt -gs:syn_groundtruth_goldstandard_25.txt -b:20 

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K50.txt -gs:syn_groundtruth_goldstandard_50.txt -b:20

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K75.txt -gs:syn_groundtruth_goldstandard_75.txt -b:20

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization_K100.txt -gs:syn_groundtruth_goldstandard_100.txt -b:20

./crank -i:syn_groundtruth.txt -p:T -c:syn_groundtruth_communities.txt -o:syn_groundtruth_prioritization.txt -gs:syn_groundtruth_goldstandard.txt -b:20

7) Prioritize communities of the synthetic planted community network using
the probabilities returned by a statical community detection model for annotated networks:

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A5_communities.txt -o:syn_groundtruth--0_prioritization_A5.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A10_communities.txt -o:syn_groundtruth--0_prioritization_A10.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A25_communities.txt -o:syn_groundtruth--0_prioritization_A25.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A50_communities.txt -o:syn_groundtruth--0_prioritization_A50.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A75_communities.txt -o:syn_groundtruth--0_prioritization_A75.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:T -c:syn_groundtruth--0A100_communities.txt -o:syn_groundtruth--0_prioritization_A100.txt -b:10

./crank -i:syn_groundtruth--0.txt -p:F -c:syn_groundtruth--0A100_communities.txt -o:syn_groundtruth--0_prioritization_A100.txt -in:syn_groundtruth--0A100_CProbaH.txt -ie:syn_groundtruth--0A100_EdgeProbaH.txt -ic:syn_groundtruth--0A100_CEdgeProbaH.txt -b:10