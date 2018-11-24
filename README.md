# CRank: Prioritizing network communities

#### Author: [Marinka Zitnik](http://stanford.edu/~marinka) (marinka@cs.stanford.edu)

#### [Project website](http://snap.stanford.edu/crank)

## Overview

This repository contains code necessary to run the CRank algorithm. CRank is an automatic unsupervised 
method for prioritizing network communities and identifying the most promising ones for further 
experimentation.
 
CRank can be used with any community detection method and scales to large networks. It uses only 
information provided by the network structure and does not require any additional external metadata or 
labels. However, when available, CRank can incorporate domain-specific metadata and labels to further boost 
performance. See our [paper](https://www.nature.com/articles/s41467-018-04948-5) for details on the algorithm.

<p align="center">
<img src="https://github.com/marinkaz/crank/blob/master/images/crank-overview.png" width="600" align="center">
</p>
  
## Usage: Zachary's Karate Club 

[Zachary's Karate Club](http://konect.uni-koblenz.de/networks/ucidata-zachary) network is a well known network, which consists of 34 nodes and 78 edges. 
It shows the relationships between Zachary Karate Club members. Each node represents a member of 
the karate club and each edge represents a ties between two members of the club.

A community detection method of user's choice takes as input the network and outputs a grouping of 
nodes into five communities, Cmt1, Cmt2, Cmt3, Cmt4, and Cmt5, as highlighted in the figure. Each 
community is given by a set of its member nodes. For example, community Cmt1 contains nodes 5, 6, 
7, 11, and 17. Notice that communities Cmt2 and Cmt4 split the friendship network among the Karate 
Club members into [two widely known factions](https://en.wikipedia.org/wiki/Zachary%27s_karate_club), 
whereas communities Cmt1, Cmt3 and Cmt5 represent less meaningful groups of the Karate Club members. 

<p align="center">
<img src="https://github.com/marinkaz/crank/blob/master/images/karate-input.png" width="600" align="center">
</p>

We now prioritize the detected communities using CRank. CRank takes as input the Zachary's Karate Club 
network given by its edge list (switch `-i:`) and node-community affiliations (switch `-c:`). As a 
result, CRank prioritizes the communities by ranking them by their aggregated prioritization score and 
saves the resulting prioritization to a file (switch `-o:`).
  
    $ ./crank -i:karate.txt -c:karate_communities.txt -o:karate_prioritization.txt
    
### Results

CRank assigns a score to each community and uses that score to determine the rank of a community in the 
final prioritization. The resulting prioritization of the five communities detected in the Zachary's network 
is shown in the figure below.

Communities Cmt4 and Cmt2 are placed at the top of the ranking, indicating that Cmt4 and Cmt2 are most 
promising communities for follow-up investigation. In contrast. Cmt5 is ranked last, 5 out of 5, 
suggesting that Cmt5 is the least promising community. Indeed, Cmt4 and Cmt2 correspond to two groups 
of people into which the karate club split after an argument between two teachers, a fact that is well 
known in the literature but that was not used for prioritization.

<p align="center">
<img src="https://github.com/marinkaz/crank/blob/master/images/karate-output.png" width="600" align="center">
</p> 
 
Notice that CRank allows the input to consist of only network and community affiliation data, given 
by switches `-i:` and `-c:`, respectively. As a result, CRank can be used with non-statistical community 
detection methods.  

See the [project website](http://snap.stanford.edu/crank) for more examples of usage. 

## Citing

If you find *CRank* useful for your research, please consider citing [this paper](https://www.nature.com/articles/s41467-018-04948-5):

    @article{Zitnik2018prioritizing,
      title     = {Prioritizing Network Communities},
      author    = {Zitnik, Marinka and Sosic, Rok and Leskovec, Jure},
      journal   = {Nature Communications},
      volume    = {9},
      number    = {1},
      pages     = {2544},
      year      = {2018}
    }

## Miscellaneous

Please send any questions you might have about the code and/or the 
algorithm to <marinka@cs.stanford.edu>.

This code implements several different variants, including the option to include user-specific community 
metrics and labels. Many prioritization variants are possible and what works 
best might depend on a concrete use case.  

## Requirements

CRank code is tested under Mac OS X, Linux and Windows systems. 

This is a C++ implementation. To compile the code, do:

    $ cd crank
    $ make all
    
CRank relies on [SNAP](https://github.com/snap-stanford/snap), a general-purpose network analysis 
and graph mining library. 

The code also includes `rra.py`, a simple Python implementation of CRank's rank aggregation algorithm.

## License

Decagon is licensed under the MIT License.