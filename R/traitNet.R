traitNet <- function(covstruc,fix_omega="full",prune=TRUE,p.adjust="fdr",alpha=0.05,reestimate=TRUE,recursive=TRUE,estimation="ML",graph_layout="mds",traits=NULL,toler=NULL){
  
  time<-proc.time()
  
  #use all traits in S matrix if nothing specified
  if(is.null(traits)){
    traits<-colnames(covstruc[[2]])
  }
  
  #restrict to variables listed in traits and smooth matrix
  covstruc<-.covQC(covstruc,traits)
  
  #estimate saturated network parameters
  print("Estimating saturated network model.")
  model_out <- .runGGM(covstruc,fix_omega="full",saturated=NULL,estimation,toler)
  model_results <- list(saturated=model_out)
  
  #estimate base sparse network parameters (if fixed omega matrix provided)
  if(is.matrix(fix_omega)){
    print("Estimating base sparse network model.")
    model_out <- .runGGM(covstruc,fix_omega,saturated=model_results$saturated,estimation,toler)
    model_results <- c(model_results, list(base=model_out))
  }
  
  #network pruning and re-estimation
  #initial network pruning
  if(prune){
      print(paste0("Pruning non-significant network edges (alpha = ",alpha,", p-value adjust = '",p.adjust,"')"))
      pruned_omega <- .pruneNet(model_out,p.adjust,alpha)
    }
    
    # restimate the network (recursively)
    if(reestimate){
      model_iterations <- list()
      iter <- 0
      repeat {
        model_out <- .runGGM(covstruc,fix_omega=pruned_omega,saturated=model_results$saturated,estimation,toler)
        model_iterations <- c(model_iterations, list(model_out))
        
        if (recursive){
          iter <- iter+1
          print(paste0("Re-estimating the network model (iteration ",iter,")."))
          pruned_omega <- .pruneNet(model_out,p.adjust,alpha)
          if (all(pruned_omega == model_out$omega)){
            if(iter==1){
              names(model_iterations) <- c("sparse")
            } else{
              names(model_iterations) <- c(paste0("iteration",1:(iter-1)),"sparse")
            }
            break
          }
        } else{
          print("Re-estimating the network model.")
          pruned_omega <- model_out$omega
          names(model_iterations) <- "sparse"
          break
        }
      }
      model_results <- c(model_results, model_iterations)
    } else{
      print("Network model not re-estimated")
      model_results <- c(model_results, list(sparse=list(omega=pruned_omega)))
    }
  }else{
    print("You have selected not to prune the edges in your network. Please ensure this is correct.")
    pruned_omega <- model_out$omega
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
