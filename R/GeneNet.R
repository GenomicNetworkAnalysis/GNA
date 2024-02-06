GeneNet <- function(covstruc,traits=NULL,fix_omega="full",simruns=100,reestimate=TRUE,recursive=TRUE,prune="fdr",alpha=0.05,graph_layout="spring",toler=NULL,prunepower=FALSE){
  
  time<-proc.time()
  
  #use all traits in S matrix if nothing specified
  if(is.null(traits)){
    traits<-colnames(covstruc[[2]])
  }
  
  #restrict to variables listed in traits and smooth matrix
  covstruc<-.covQC(covstruc,traits)
  
  #estimate network paramterers
  print("Estimating network model.")
  model_out <- .runGGM(covstruc,fix_omega,toler)

  #simulations to estimate power for each partial rg
  if(is.numeric(simruns)){
  powerNet<-.simNet(covstruc,simruns,prune)
  model_out$parameters$power<-c(powerNet,rep(NA,ncol(covstruc$S_LD)))
  }

  #prune network
  if(prune == "none"){
    print("You have selected not to prune the edges in your network. Please ensure this is correct.")
    pruned_omega <- model_out$omega
  } else{
    if(prunepower){
       print(paste0("Pruning non-significant network edges ('",prune,"' adjusted p-value > ",alpha,") and edges with < 80% power in simulation"))
         pruned_omega <- .pruneNet(model_out,prune,alpha,prunepower)
      }
    else{print(paste0("Pruning non-significant network edges ('",prune,"' adjusted p-value > ",alpha,")"))
    pruned_omega <- .pruneNet(model_out,prune,alpha)
  }
      }
  model_results <- list(c(model_out, list(pruned_omega=pruned_omega)))

  #restimate the model (recursively)
  if(reestimate && (prune != "none")){
    iter <- 0
    repeat {
      model_out <- .runGGM(covstruc,fix_omega=pruned_omega,toler)
      if (recursive){
        iter <- iter+1
        print(paste0("Re-estimating the network model (iteration ",iter,")."))
        pruned_omega <- .pruneNet(model_out,prune,alpha)
        model_results <- c(model_results, list(c(model_out, list(pruned_omega=pruned_omega))))
        if (all(pruned_omega == model_out$omega)) break
      } else{
        print("Re-estimating the network model.")
        pruned_omega <- model_out$omega
        model_results <- c(model_results, list(c(model_out, list(pruned_omega=pruned_omega))))
        break
      }
    }
  } else{
    print("Network model not re-estimated")
  }
  
  
  #network description - plotting and centrality metrics
  if(all(pruned_omega == 0)){
    warning("There are no significant edges using the current pruning threhsold. A network graph will not be created")
    network <- NULL
  } else{
    network <- .describeNet(pruned_omega,graph_layout)
    }

  #function output
  output <- list(model_results=model_results, network=network)
  return(output)

  time_all<-proc.time()-time
  print(time_all[3])
}
