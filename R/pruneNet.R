.pruneNet <- function(model_out,traits,prune="bonf",alpha=0.05,threshold=10){
  
  par <- model_out$parameters[model_out$parameters$matrix=="omega"]
  omega <- model_out$omega
  
  # prune non significant edges
  
  if(prune == "bayes"){
    if(bayes){
      if(anyNA(par$BayesFactor10)){
        warning("Bayes factor estimation produced non-finite values. Consider using a different pruning method")
      }
      par$weight[par$BayesFactor10 < threshold] <- 0
    }else{print("You have specified pruning using bayes estimation but did not request bayes factors to be estimated.")
    }
  }
  
  if(prune == "bonf"){
    par$weight[par$p >= (0.05/nrow(par))] <- 0
  }
  
  if(prune == "fdr"){
    par$weight[p.adjust(par$p, method='fdr') >= 0.05] <- 0
  }
  
  if(prune == "alpha"){
    par$weight[par$p >= alpha] <- 0
  }
  
  #update omega with pruned weights
  omega[par$row,par$col] <- par$weight
  
  #reform parameter table
  params <- rbind(par, model_out$parameters[model_out$parameters$matrix != "omega"])
  
  return(list(parameters=params,omega_pruned=omega))
}
