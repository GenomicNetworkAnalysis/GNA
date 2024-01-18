GeneNet <- function(covstruc,traits=NULL,fix_omega="full",simruns=100,prune="bonf",alpha=0.05,threshold=10,graph_layout="spring",bayes=TRUE,toler=NULL){
  
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

  #calculate Bayes Factor for each partial rg
  if(bayes){
  BF10 <- .BayesF(covstruc,p_rg,traits)
  p_rg <- merge(p_rg, BF10, by = c("Trait1","Trait2"))
  }

  #simulations to estimate power for each partial rg
  if(is.numeric(simruns)){
  powerNet<-.simNet(covstruc,simruns)
  model_out$parameters$power<-powerNet
  }

  #prune network
  if(prune == "none"){
    pruned_omega <- model_out$omega
    print("You have selected not to prune the edges in your network. Please ensure this is correct.")
  } else{
    pruned_omega <- .pruneNet(model_out,prune,alpha,threshold)
  }
  
  #network description - plotting and centrality metrics
  if(all(pruned_omega == 0)){
    network <- NULL
    warning("There are no significant edges using the current pruning threhsold. A network graph will not be created")
  } else{
    network <- .describeNet(pruned_omega,graph_layout)
    }

  #function output
  output <- list(model_out,pruned_omega,network)
  return(output)

  time_all<-proc.time()-time
  print(time_all[3])
}
