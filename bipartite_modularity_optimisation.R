require(igraph)
require(Matrix)
require(plyr)
require(doSNOW)
#source("../../../com_tracing/matching/matching.R")

numLabelsUsed <- 0

bipartiteModularityGain <- function(g, sOld, sNew, changedNode) {
  
  if (is.null(E(g)$weight)) {
    
    adj <- get.adjacency(g)
    # new
    l <- length(E(g))
  } else {
    
    adj <- get.adjacency(g, attr="weight")
    #new: considering edge weight
    l <- sum(E(g)$weight)
  }
  #l <- length(E(g))
  nodesInCold <- which(sOld == sOld[changedNode])
  nodesInColdType0 <- intersect(nodesInCold, which(V(g)$type == 0))
  nodesInColdType1 <- intersect(nodesInCold, which(V(g)$type == 1))
  
  nodesInCnew <- setdiff(which(sNew == sNew[changedNode]), changedNode)
  nodesInCnewType0 <- intersect(nodesInCnew, which(V(g)$type == 0))
  nodesInCnewType1 <- intersect(nodesInCnew, which(V(g)$type == 1))
  #nodesInC <- which(sNew == sNew[changedNode])
  deg <- rowSums(adj)#degree(g, weighted=TRUE)
  #l <- length(E(g))
  
  ki_in <- sum(adj[changedNode, nodesInCold])
  ki <- deg[changedNode]
  
  if (V(g)[changedNode]$type == 0) {
    
    sumKother <- sum(deg[nodesInColdType1])
  } else {
    
    sumKother <- sum(deg[nodesInColdType0])
  }
 
  if ((length(nodesInColdType0) == 0) || (length(nodesInColdType1) == 0)) {
    
    q0 <- 0
    #c_in <- 0
    #c_tot <- 0
    
  } else {
    c_in <- sum(adj[nodesInColdType0, nodesInColdType1])
    c_tot <- sum(sapply(nodesInColdType0, function(i) { #sum(deg[nodesInCold])
        
      sapply(nodesInColdType1, function(j) {
        
        deg[i] * deg[j]
      })
    }))
    
    
    q0 <- (((c_in - ki_in) / l) - ((c_tot - (ki * sumKother)) / l^2)) - ((c_in / l) - (c_tot / l^2))
  }
  
  #q0
  
  ki_in <- sum(adj[changedNode, nodesInCnew])
  ki <- deg[changedNode]
  
  if (V(g)[changedNode]$type == 0) {
    
    sumKother <- sum(deg[nodesInCnewType1])
  } else {
    
    sumKother <- sum(deg[nodesInCnewType0])
  }
  
  if ((length(nodesInCnewType0) == 0) || (length(nodesInCnewType1) == 0)) {
    
    #q1 <- 0
    c_in <- 0
    c_tot <- 0
  } else {
    
    c_in <- sum(adj[nodesInCnewType0, nodesInCnewType1])
    c_tot <- sum(sapply(nodesInCnewType0, function(i) { #sum(deg[nodesInCold])
      
      sapply(nodesInCnewType1, function(j) {
        
        deg[i] * deg[j]
      })
    }))
  }

  q1 <- (((c_in + ki_in) / l) - ((c_tot + (ki * sumKother)) / l^2)) - ((c_in / l) - (c_tot / l^2))
  #q1
  
  q0 + q1
 
}

bipartiteModularity <- function(g, s) {
  
  nodesType0 <- which(V(g)$type == 0)
  nodesType1 <- which(V(g)$type == 1)
  adj <- get.adjacency(g)
  
  deg <- rowSums(adj)
  
  sum(sapply(nodesType0, function(i) {
    
    sapply(nodesType1, function(j) {
      
      if (s[i] == s[j]) {
        
        adj[i,j] - ((deg[i] * deg[j]) / length(E(g)))
        
      } else {
        
        0
      } 
    })  
  })) / length(E(g))  
}

improvement <- TRUE

optimize <- function(g) {
  
  #print("optimize")
  cl <- 
  vType1 <- which(V(g)$type == 0)
  vType2 <- which(V(g)$type == 1)
  
  cl <- makeCluster(4, type = "SOCK", outfile="snow_log.log")
  
  clusterEvalQ(cl, require(igraph))
  clusterEvalQ(cl, require(Matrix))
  clusterExport(cl, "bipartiteModularityGain")
  clusterExport(cl, "g", envir = environment())
  registerDoSNOW(cl)
  
  update <- TRUE
  numUpdates <- 0
  while(update) {
    update <- FALSE

    ncPairs <- do.call(rbind, lapply(vType1, function(node) {
      adjClusters <- setdiff(V(g)[neighbors(g, node)]$cluster, V(g)[node]$cluster)
      if (length(adjClusters) == 0) {
        
        adjClusters <- NA
      }
      cbind(node, adjClusters)
    }))
    ncPairs <- as.data.frame(ncPairs)
    names(ncPairs) <- c("node", "cluster")
    
    bestFit <- 0
    fits <- ddply(ncPairs, "node", function(sdata) {
      print(paste("num nodes", length(V(g))))
      if (!is.na(sdata$cluster)) {
        
        gains <- sapply(sdata$cluster, function(c) {
          
          s1 <- V(g)$cluster
          s2 <- s1
          s2[sdata$node[1]] <- c
          bipartiteModularityGain(g, s1, s2, sdata$node[1])  
        }) 
        
        cbind(sdata, data.frame(gain=gains))  
      } else {
        
        cbind(sdata, data.frame(gain=0))
      }
      
    }, .parallel = TRUE)
    
    maxGain <- which(fits$gain == max(fits$gain))[1]
    if (maxGain > bestFit) {
      
      bestFit <- fits$gain[maxGain]
      V(g)$cluster[fits$node[maxGain]] <- fits$cluster[maxGain]
      
      update <- TRUE
      numUpdates <- numUpdates + 1
    }
    
    
    cPairs <- do.call(rbind, lapply(vType2, function(node) {
      adjClusters <- setdiff(V(g)[neighbors(g, node)]$cluster, V(g)[node]$cluster)
      if (length(adjClusters) == 0) {
        
        adjClusters <- NA
      }
      cbind(node, adjClusters)
    }))
    ncPairs <- as.data.frame(ncPairs)
    names(ncPairs) <- c("node", "cluster")
    
    bestFit <- 0
    fits <- ddply(ncPairs, "node", function(sdata) {
      
      if (!is.na(sdata$cluster)) {
        
        gains <- sapply(sdata$cluster, function(c) {
          
          s1 <- V(g)$cluster
          s2 <- s1
          s2[sdata$node[1]] <- c
          bipartiteModularityGain(g, s1, s2, sdata$node[1])  
        }) 
        
        cbind(sdata, data.frame(gain=gains))  
      } else {
        
        cbind(sdata, data.frame(gain=0))
      }
    }, .parallel = TRUE)
    
    maxGain <- which(fits$gain == max(fits$gain))[1]
    if (maxGain > bestFit) {
      
      bestFit <- fits$gain[maxGain]
      V(g)$cluster[fits$node[maxGain]] <- fits$cluster[maxGain]
      
      update <- TRUE
      numUpdates <- numUpdates + 1
    }
    
    print(paste("best fit:", bestFit))
  }
  stopCluster(cl)
  g
}

bipartiteLouvaine <- function(g, preclustering=FALSE) {
  
  if (!preclustering) {
    
    V(g)$cluster <- 1:length(V(g))  
  }
  print("first optimization")
  initGraph <- optimize(g)

  clusterNodeMapping <- list()
  scdLevelGraph <- graph.empty()
  i <- 1
  clustNames <- unique(V(initGraph)$cluster)
  #print("modularity: ")
  #print(bipartiteModularity(initGraph, V(initGraph)$cluster))
  scdLevelGraph <- add.vertices(scdLevelGraph, length(clustNames))
  sapply(1:length(clustNames), function(i) {
    
    nodesInCluster1 <- which(V(initGraph)$cluster == clustNames[i])
    clusterNodeMapping[[i]] <<- nodesInCluster1
    
    sapply(1:length(clustNames), function(j) {
  
      nodesInCluster2 <- which(V(initGraph)$cluster == clustNames[j])
      
      nb <- unlist(neighborhood(initGraph, 1, nodesInCluster2))
      #nb <- nb[which(!(nb %in% nodesInCluster2)] #setdiff(nb, nodesInCluster2)  WRONG!!!
      #linkWeight <- length(intersect(nb, nodesInCluster1))
      
      linkWeight <- length(which(V(initGraph)[nb]$cluster == clustNames[i]))
      
      if (i == j) {
        # because of the included 0 neighborhood.
        linkWeight <- (linkWeight - length(nodesInCluster2)) / 2
      } 
      
      if (linkWeight > 0) {
        
        scdLevelGraph <<- scdLevelGraph + edge(c(i, j), weight=linkWeight)  
      }
    })
  })
  
  scdLevelGraph <- as.undirected(scdLevelGraph)
  
  #plot(scdLevelGraph, edge.label=E(scdLevelGraph)$weight)
  mlRes <- multilevel.community(scdLevelGraph)

  sapply(unique(mlRes$membership), function(mem) {
    
    help <- which(mlRes$membership == mem)
    V(g)[unlist(clusterNodeMapping[help])]$cluster <<- mem
  })
  #print("modularity: ")
  #print(bipartiteModularity(g, V(g)$cluster))
  # workaround for the case that the initial clustering is already a local optimum.
  if (bipartiteModularity(initGraph, V(initGraph)$cluster) > bipartiteModularity(g, V(g)$cluster)) {
    
    V(g)$cluster <- V(initGraph)$cluster
  } else {
    
    print("second optimization")
    #print("modularity: ")
    g <- optimize(g)  
    #print(bipartiteModularity(g, V(g)$cluster))
  }
  
  g
}

start <- function(g) {
  
  improvement <<- TRUE
  res <- bipartiteLouvaine(g)

  subgroups <- lapply(unique(V(res)$cluster), function(c) {
    
    which(V(res)$cluster == c)
  })
  plot(res, mark.groups=subgroups)
  print(V(res)$cluster)
  print(bipartiteModularity(res, V(res)$cluster))
  mc <- multilevel.community(g)
  print(mc$membership)
  print(bipartiteModularity(g, mc$membership))
}


