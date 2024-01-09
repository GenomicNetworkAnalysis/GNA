.runNet <- function(p_rg,traits,prune="bonf",alpha=0.05,threshold=10,graph_layout="spring"){
  
  # prune non significant edges
  p_rg$weight <- p_rg$Pcor_Estimate
  
  if(prune == "bayes"){
    if(bayes){
      if(anyNA(p_rg$BayesFactor10)){
        warning("Bayes factor estimation produced non-finite values. Consider using a different pruning method")
        }
    p_rg$weight[p_rg$BayesFactor10 < threshold] <- 0
      }else{print("You have specified pruning using bayes estimation but did not request bayes factors to be estimated.")
           }
  }
  
  if(prune == "bonf"){
    p_rg$weight[p_rg$Pcor_pvalue >= (0.05/nrow(p_rg))] <- 0
  }
  
  if(prune == "fdr"){
    p_rg$p_fdr <- p.adjust(p_rg$Pcor_pvalue, method='fdr')
    p_rg$weight[p_rg$p_fdr >= 0.05] <- 0
  }
  
  if(prune == "alpha"){
    p_rg$weight[p_rg$Pcor_pvalue >= alpha] <- 0
  }
  
  if(prune == "none"){
    p_rg$weight <- p_rg$weight
    print("You have selected not to prune the edges in your network. Please ensure this is correct.")
  }
  
  if(all(p_rg$weight==0)){
    warning("There are no significant edges using the current pruning threhsold. A network graph will not be created")
    out <- list(weights = NULL, graph = NULL, centrality = NULL)
}else{
  # create weights matrix
  val <- as.matrix(reshape2::dcast(p_rg, Trait1 ~ Trait2, value.var = "weight"))
  rownames(val) <- val[,1]
  mat <- matrix(nrow = length(traits), ncol = length(traits), dimnames = list(traits,traits))
  cols <- colnames(mat)[colnames(mat) %in% colnames(val)]
  rows <- rownames(mat)[rownames(mat) %in% rownames(val)]
  mat[rows, cols] <- val[rows, cols]
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- colnames(mat)
  diag(mat) <- 0
  
  # network graph
  ifelse(graph_layout == "mds",
         layout <- mds(sim2diss(mat), type = "ordinal")$conf, #multidimensional scaling
         layout <- graph_layout)
  
  network <- qgraph(mat, layout=layout, vsize = 7, threshold = 0, theme = "Borkulo")
  
  
  
  # centrality metrics
  centr <- centralityTable(network)
  centr <- as.data.frame(reshape2::dcast(centr, node ~ measure, value.var = "value"))
  centr <- cbind(centr[,"node"],apply(centr[,c("Betweenness","Closeness","Strength","ExpectedInfluence")], 2, as.numeric))
  colnames(centr)[1] <- "Trait"
  
  
  out <- list(mat,network,centr)
  names(out) <- c("weights","graph","centrality")
    }
  
  return(out)

}
