.describeNet <- function(omega,graph_layout){

  # network graph
  ifelse(graph_layout == "mds",
         layout <- mds(sim2diss(omega), type = "ordinal")$conf, #multidimensional scaling
         layout <- graph_layout)
  
  graph <- qgraph(omega, layout=layout, vsize = 7, threshold = 0, theme = "Borkulo")
  
  # centrality metrics
  centr <- centralityTable(graph)
  centr <- as.data.frame(reshape2::dcast(centr, node ~ measure, value.var = "value"))
  colnames(centr)[1] <- "Trait"
  
  return(list(graph=graph,centrality=centr))
}
