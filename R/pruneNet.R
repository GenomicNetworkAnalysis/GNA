.pruneNet <- function(model_out,prune="fdr",alpha=0.05,prunepower=FALSE){
  
  par <- model_out$parameters[model_out$parameters$matrix=="omega",] #pull parameters table
  par$p_adj <- p.adjust(par$p, method=prune) #adjusted pvals
  par$weight <- if_else(par$p_adj < alpha, par$est, 0) #prune non-significant edges
  
  if(prunepower){
   par$weight[par$power < .8] <- 0 
  }
  
  #update omega with pruned weights
  omega <- model_out$omega
  if(!all(par$weight == par$est)){
    for (i in 1:nrow(par)){
      omega[par$row[i],par$col[i]] <- par$weight[i]
      omega[par$col[i],par$row[i]] <- par$weight[i]
    }
  }
  return(omega)
}
