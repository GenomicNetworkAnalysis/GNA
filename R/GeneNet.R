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

  #prune and plot the network
  network <- .runNet(model_out,traits,prune,alpha,threshold,graph_layout)

  #function output
  output <- c(model_out, network)
  names(output) <- c("partial_rgs","weights","graph","centrality")
  
  return(output)

  time_all<-proc.time()-time
  print(time_all[3])
}
