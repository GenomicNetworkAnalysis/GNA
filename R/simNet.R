.simNet<-function(covstruc,simruns,p.adjust,alpha){
  
  print("Beginning estimation of simulated partial correlations to produce power estimates for each edge weight.
        This step typically takes a few minutes.")
  
  #sample vector of genetic covariance matrix (S) estimates using:
  #observed S (genetic covariance matrix) as mean vector
  #observed V (sampling covariance matrix) as Sigma (the covariance among the estimaes)
  vectorS<- mvrnorm(simruns,covstruc$S_LD[lower.tri(covstruc$S_LD,diag = T)],Sigma = covstruc$V_LD)
  
  #create empty data.frame to store simulation results with rows = number of unique elements in S
  k<-ncol(covstruc$S_LD)
  simresults<-as.data.frame(matrix(ncol=simruns,nrow=(((k*(k+1))/2)-k)))
  
  #estimate partial correlations for each simnulation
  for(r in 1:simruns){
    
    #create empty matrix to store simulated results 
    simS<-matrix(0, ncol(covstruc$S_LD), ncol(covstruc$S_LD))
    
    #store vector for one simulation run in lower triangle
    simS[lower.tri(simS,diag=TRUE)] <-as.numeric(vectorS[r,])
    
    #reflect values above the diagonal
    simS <- simS + t(simS) - diag(diag(simS)) 
    
    #name columns
    colnames(simS)<-colnames(covstruc$S_LD)
    
    ##smooth to near positive definite if simulate Sis non-positive definite
    if(eigen(simS)$values[nrow(simS)] <= 0){
      simS<-as.matrix((nearPD(simS, corr = FALSE))$mat)
    }
    
    #pair simulated S matrix with observed V
    sim_covstruc<-list(V_LD=covstruc$V_LD,S_LD=simS)
    
    traits<-colnames(covstruc$S_LD)
    
    #estimate the partial correlations (network weights) for simulation
    p_rg_sim <- .runGGM(sim_covstruc,toler=NULL)
    
    #save only the partial correlation part of simulated output
    p_rg_sim<-p_rg_sim$parameters
    
    #remove the variance estimates
    p_rg_sim<-p_rg_sim %>% filter(row_number() <= n()-k)

    #save the trait pairs as a consistent first three columns across runs
    if(r == 1){
      rownames(p_rg_sim)<-paste0(p_rg_sim$trait1, p_rg_sim$op,p_rg_sim$trait2,sep="")
    }
  
    #store p-value for each partial correlation for given simulation 
    simresults[,r]<-p.adjust(p_rg_sim$p, method=p.adjust)
    
  }
  
  #calculate proportion of runs that were significant for each parameter (i.e., power)
  powerNet<-(apply(simresults, 1, function(row) sum(row < alpha)))/simruns
  
  return(powerNet)

}
