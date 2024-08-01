# Maximal Structural Balanced k-plex Enumeration in Signed Graphs

This project aims to enumerate all maximal (structural) balanced k-plexes in signed graphs.

MBPE and MBPE-opt are the executables, and are compiled on Ubuntu 18.04.5, with -O3 optimization.

All of our datasets are publically available, and most of them are downloaded from [SNAP](https://snap.stanford.edu/data/index.html) and [KONECT](http://konect.cc/networks/). 

## Running Format

./MBPE [1]input_graph [2]k [3]tau [4]Algorithm:A/B

**Running example for maximal balanced k-plex enumeration by our baseline algorithm**

./MBPE ./datasets/bitcoin.txt 2 3 B

**Running example for maximal balanced k-plex enumeration by our advanced algorithm**

./MBPE ./datasets/bitcoin.txt 2 3 A

**Running example for maximal balanced k-plex enumeration by our advanced algorithm with minimum-degree branching strategy**

./MBPE-opt ./datasets/bitcoin.txt 2 3 A

### Note

The value of k should be at least 2, and tau should not be smaller than k.

To test each technique, you can choose pre-defined parameters in the file "Utility.h" and re-compile the executables. Below are the functions of each parameter:

**\_CTprune\_**: Applies vertex and edge reductions on the original graph.

**\_ParVR\_** and **\_PP\_**: Applies the partition-based vertex reduction technique.

**\_VRinEnum\_**: Applies the subgraph reduction technique.

**\_PIVOT\_**: Applies the pivot technique on the baseline algorithm.

**\_mindegpivot\_**: Applies the minimum-degree branching technique.

## Graph Format

number of vertices \t number of edges \n

v0 \t v1 \t sign \n

v0 \t v2 \t sign \n

...
