# bipartite_community_detection

This repository contains R scripts for clustering biparite networks. 
All scripts contain a method start() with example code. Please make sure that you have installed the necessary packages.
```R
install.packages("igraph", "dplyr", "foreach")
```
**Nodes of the network must have an attribute 'type', indicating the mode of each node in the bipartite network (0 or 1).**
```R
 g <- read.csv("davis.csv") %>% graph_from_data_frame(directed=FALSE)
 V(g)[1:18]$type <- 1
 V(g)[1:18]$type <- 0
```

## bipartite_cpm.R
This script implements the biclique percolation algorithm introduced by Lehman, Schwartz, and Hansen (2008).

Example: Clusters based on 4,5 bicliques
`clusters <- cpm(g, 4, 5)`

Lehmann, S., Schwartz, M., & Hansen, L. K. (2008). Biclique communities. Physical review E, 78(1), 016108.

## bipartite_modularity_optimisation.R
This script implements an adaptation of the Louvain algorithm for bipartite networks.

```R
res <- bipartiteLouvaine(g)

subgroups <- lapply(unique(V(res)$cluster), function(c) {
   
  which(V(res)$cluster == c)
})
plot(res, mark.groups=subgroups)
```
- Hecking, T., Steinert, L., Göhnert, T., & Hoppe, H. U. (2014). Incremental clustering of dynamic bipartite networks. In Proceedings of the 1st European Network Intelligence Conference (pp. 9-16). IEEE.

- Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. Journal of statistical mechanics: theory and experiment, 2008(10), P10008.

