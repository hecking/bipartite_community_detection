require(igraph)
require(foreach)
require(dplyr)

computeBicliques <- function(graph, k, l) {
  
  graph <- connect.neighborhood(graph, 2)
  
  cliqueFrame <- data_frame(clique=max_cliques(graph))
  
  bicliques <- list()
  
  cliqueFrame %>%
    rowwise %>%
    do(biclique = {
      type0Members <- which(V(graph)[.$clique]$type == 0)
      type1Members <- which(V(graph)[.$clique]$type == 1)
      
      if ((length(type0Members) >= k) && (length(type1Members) >= l)) {
        
        list(m1=.$clique[type0Members], m2=.$clique[type1Members])
      }
    }) %>% 
    filter(!is.null(biclique))
}

cpm <- function(graph, k, l, runParallel=FALSE, numCores=4) {

  graph <- simplify(graph)
  
  clq <- computeBicliques(graph, k, l)

  if (runParallel) {
    
    require(doParallel)
    cl <- makeCluster(numCores)
    registerDoParallel(cl)
    edges <- foreach(i=1:(nrow(clq) - 1), .combine=rbind, .packages = "foreach") %dopar% {
        foreach(j=(i + 1):nrow(clq), .combine=rbind) %do% {
          if((length(intersect(clq$biclique[[i]]$m1, clq$biclique[[j]]$m1)) >= k-1) &&
             (length(intersect(clq$biclique[[i]]$m2, clq$biclique[[j]]$m2)) >= l-1)) {
            
            c(i,j)
          }
        }
      }
    stopCluster(cl)
  } else {
    
    edges <- foreach(i=1:(nrow(clq) - 1), .combine=rbind) %do% {
      
      foreach(j=(i + 1):nrow(clq), .combine=rbind) %do% {
        print(j)
        if((length(intersect(clq$biclique[[i]]$m1, clq$biclique[[j]]$m1)) >= k-1) &&
           (length(intersect(clq$biclique[[i]]$m2, clq$biclique[[j]]$m2)) >= l-1)) {

          c(i,j)
        } 
      }
    }
  }
  
    
  bicliqueCommunities <- list()
  if(!is.null(edges) && (nrow(edges) > 0)) {
    
    cliqueGraph <- graph_from_data_frame(edges, directed=FALSE, vertices = data_frame(name=1:nrow(clq)))
    
    V(cliqueGraph)$name <- 1:vcount(cliqueGraph)
    
    # Communities are the components of the biclique graph
    comps <- decompose(cliqueGraph)
    
    bicliqueCommunities <- lapply(comps, function(component) {
      
      clq$biclique[as.numeric(V(component)$name)] %>% unlist %>% unique
    })
  }
  
  bicliqueCommunities
}

start <- function() {
  inFile <- "your network as gml, graphml"
  fileFormat <- "gml" # or graphml
  k <- 3 # k parameter (see Lehmann paper)
  l <- 4 # l parameter
  
  # Nodes of the network must have an attribute 'type', indicating the mode of each node in the bipartite network (0 or 1).
  g <- read.graph(inFile, fileFormat)
 
  
  clusters <- cpm(g, k, l)
  # cpm(g, k, l, runParallel=TRUE, numCores=4): Parallel computation. (Might need modification depending on the OS)
}
