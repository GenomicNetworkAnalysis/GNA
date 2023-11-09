.pcor <-function(covstruc,model){ 

  ##read in the LD portion of the V (sampling covariance) matrix
  V_LD<-as.matrix(covstruc[[1]])
  
  ##read in the LD portion of the S (covariance) matrix
  S_LD<-as.matrix(covstruc[[2]])
  
  ##k = number of phenotypes in dataset (i.e., number of columns in LD portion of S matrix)
  k<-ncol(S_LD)
  
  ##size of V matrix used later in code to create diagonal V matrix
  z<-(k*(k+1))/2
  
  #take inverted diagonal of V for the weight matrix
  #note that this does not need reordering to match internal ordering of lavaan
  #as these partial correlation models always match the original ordering in V
  W<-diag(z)
  diag(W)<-diag(V_LD)
  W<-solve(W)
  
  ##run the model 
  Model1_Results <- sem(model, sample.cov = S_LD, estimator = "DWLS", std.lv=TRUE,WLS.V = W, sample.nobs = 2,optim.dx.tol = +Inf)

  #pull the delta matrix (this doesn't depend on N)
  ##note that while the delta matrix is reordered based on the ordering in the model specification
  ##that the lavaan output is also reordered so that this actually ensures that the results match up 
  S2.delt <- lavInspect(Model1_Results, "delta")
  
  ##weight matrix from stage 2. S2.W is not reordered by including something like model constraints
  S2.W <- lavInspect(Model1_Results, "WLS.V") 
  
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread <- solve(t(S2.delt)%*%S2.W%*%S2.delt)

  #create the "lettuce" part of the sandwich
  lettuce <- S2.W%*%S2.delt
    
  #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
  Ohtt <- bread %*% t(lettuce)%*%V_LD%*%lettuce%*%bread  
    
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.matrix(sqrt(diag(Ohtt)))

  results<-data.frame(inspect(Model1_Results, "list")[,c(2:4,8,14)])
  results<-subset(results, results$free != 0)                    
  results$free<-NULL
  results$SE<-SE
  colnames(results)=c("Trait1","op","Trait2","Pcor_Estimate", "Pcor_SE")
   
  results<-results[results$Trait1 == "F1" & results$Trait2 == "F2",]
  
  return(results)
    
}
