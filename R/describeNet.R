.describeNet <- function(omega,graph_layout){
  
  #### NETWORK GRAPH
  base_graph <- qgraph(omega, DoNotPlot = TRUE)
  
  if (graph_layout == "spring" | graph_layout == "circle") {
    layout <- graph_layout
  } else if (graph_layout == "mds") {
    layout <- MDSnet(base_graph)
  } else if (graph_layout == "pca") {
    layout <- PCAnet(base_graph)
  } else if (graph_layout == "eigen") {
    layout <- EIGENnet(base_graph)
  } else {
    layout <- "spring"
    print(paste0("Graph layout '",graph_layout,"' is not supported. Use one of 'spring', 'circle', 'mds', 'pca' or 'eigen'.  Plotting the network using the default 'spring' layout"))
  }
      
  graph <- qgraph(base_graph, layout = layout, vsize = 7, threshold = 0, theme = "Borkulo", DoNotPlot = FALSE)
  
  
  #### LOCAL NETWORK METRICS
  # centrality coefficients
  centr <- centralityTable(omega, standardized = FALSE)
  centr <- as.data.frame(reshape2::dcast(centr, node ~ measure, value.var = "value"))
  
  # clustering coefficients
  clust <- clusteringTable(omega, signed = TRUE, standardized = FALSE)
  clust <- as.data.frame(reshape2::dcast(clust, node ~ measure, value.var = "value"))
  
  local_metrics <- list(centr,clust)
  names(local_metrics) <- c("centrality","clustering")
  
  
  #### GLOBAL NETWORK METRICS
  g_adj <- networktools::coerce_to_adjacency(graph)
  n_nodes <- ncol(g_adj)
  n_edges <- sum(g_adj[lower.tri(g_adj,diag=F)] != 0)
  
  ## clustering coefficient
  clust <- NetworkToolbox::clustcoeff(g_adj)$CC
  
  # mean clustering in 1000 lattice networks
  clust_latt <- sapply(1:1000, function(i) {
    set.seed(i)
    net <- NetworkToolbox::lattnet(nodes=n_nodes, edges=n_edges)
    clust_l <- NetworkToolbox::clustcoeff(net)$CC
    return(clust_l)
  })
  clust_mean <- mean(clust_latt)
  clust_low <- as.numeric(quantile(clust_latt, 0.025))
  clust_upp <- as.numeric(quantile(clust_latt, 0.975))
  
  ## average shortest path length
  path <- NetworkToolbox::pathlengths(g_adj)$ASPL
  
  # mean average path length in 1000 random networks
  path_rand <- sapply(1:1000, function(i) {
    set.seed(i)
    net <- NetworkToolbox::randnet(nodes=n_nodes, edges=n_edges)
    path_r <- NetworkToolbox::pathlengths(net)$ASPL
    return(path_r)
  })
  path_mean <- mean(path_rand)
  path_low <- as.numeric(quantile(path_rand, 0.025))
  path_upp <- as.numeric(quantile(path_rand, 0.975))
  
  ## smallworldness
  small_world <- (path_mean/path)-(clust/clust_mean)
  
  global_metrics <- setNames(c(clust,path,small_world,clust_mean,clust_low,clust_upp,path_mean,path_low,path_upp),
                             c("ClustCoef","AvgPathLength","SmallworldIndex","ClustCoef_sim_mean","ClustCoef_sim_low","ClustCoef_sim_upp","AvgPathLength_sim_mean","AvgPathLength_sim_low","AvgPathLength_sim_upp"))
  
  return(list(graph=graph,local_metrics=local_metrics,global_metrics=global_metrics))
}
