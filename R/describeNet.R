.describeNet <- function(omega,graph_layout){
  
  ### NETWORK GRAPH
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
  
  
  ### BASIC LOCAL AND GLOBAL NETWORK METRICS
  
  ## Local
  # centrality coefficients
  centr <- centralityTable(omega, standardized = FALSE)
  centr <- as.data.frame(reshape2::dcast(centr, node ~ measure, value.var = "value"))
  
  # clustering coefficients
  clust <- clusteringTable(omega, signed = TRUE, standardized = FALSE)
  clust <- as.data.frame(reshape2::dcast(clust, node ~ measure, value.var = "value"))
  
  local_metrics <- list(centr,clust)
  names(local_metrics) <- c("centrality","clustering")
  
  ## Global
  # smallwordness, average shortest path length, transitivity..
  global_metrics <- as.list(smallworldness(omega, B = 1000, up = 0.975, lo = 0.025))
  
  return(list(graph=graph,local_metrics=local_metrics,global_metrics=global_metrics))
}
