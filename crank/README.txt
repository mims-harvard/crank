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

The algorithm is described in:

M. Zitnik, R. Sosic and J. Leskovec, Prioritizing Network Communities,
Nature Communications, 9, 2544, 2018.

https://www.nature.com/articles/s41467-018-04948-5

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

http://snap.stanford.edu/crank


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

Prioritize 5 communities of Zachary's Karate club members (Cmt2 and Cmt4
represent two factions in the Karate club; Cmt1, Cmt3 and Cmt5 are
less meaningful groupings of club members):

./crank -i:karate.txt -c:karate_communities.txt -o:karate_prioritization.txt
