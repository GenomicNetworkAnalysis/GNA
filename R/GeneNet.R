GeneNet <- function(covstruc,traits=NULL,simruns=100){
  
  time<-proc.time()
  
  #use all traits in S matrix if nothing specified
  if(is.null(traits)){
    traits<-colnames(covstruc$S_LD)
  }
  
  #restrict to variables listed in traits and smooth matrix
  covstruc<-.covQC(covstruc,traits)
  
  print("Estimating partial genetic correlations.")
  
  p_rg <- data.frame()
  trait_pairs <- as.matrix(t(combn(traits, 2)))
  for (i in 1:nrow(trait_pairs)){
    y <- trait_pairs[i,]
    x <- traits[!traits %in% y]
    model <- paste0(paste0(y, collapse = "+"), "~", paste0(x, collapse = "+"),"
                    F1=~",y[1],"
                    F2=~",y[2],"
                    F1~~F2
                    ",y[1],"~~0*",y[1],"
                    ",y[2],"~~0*",y[2])
    output <- .pcor(covstruc, model)
    p_rg <- rbind(p_rg,output)
    p_rg[i,"Trait1"] <- y[1]
    p_rg[i,"Trait2"] <- y[2]
  }
  
  p_rg$pvalue<-2*pnorm(abs(p_rg$Pcor_Estimate/p_rg$Pcor_SE),lower.tail=FALSE) 
  
  rownames(p_rg)<-NULL
  p_rg<-data.frame(p_rg)
  
  colnames(p_rg)<-c("Trait1","op","Trait2","Pcor_Estimate", "Pcor_SE", "Pcor_pvalue")

  #calculate Bayes Factor for each partial rg
  BayesF <- .BayesF(covstruc,p_rg,traits)
  p_rg <- merge(p_rg, BayesF, by = c("Trait1","Trait2"))

  #simulations to estimate power for each partial rg
  if(is.numeric(simruns)){
  powerNet<-.simNet(covstruc,simruns,trait_pairs)
  p_rg$power<-powerNet
  }
  
  time_all<-proc.time()-time
  print(time_all[3])
  
  return(p_rg)
}
