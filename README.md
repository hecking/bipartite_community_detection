# bipartite_community_detection

This repository contains R scripts for clustering biparite networks. 
All scripts contain a method start() with example code. Please make sure that you have installed the necessary packages.
```R
install.packages("igraph", "dplyr", "foreach")
```
Nodes of the network must have an attribute 'type', indicating the mode of each node in the bipartite network (0 or 1).

## bipartite_cpm.R
This script implements the biclique percolation algorithm introduced by Lehman, Schwartz, and Hansen (2008)

Lehmann, S., Schwartz, M., & Hansen, L. K. (2008). Biclique communities. Physical review E, 78(1), 016108.

