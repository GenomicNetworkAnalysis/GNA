gwasNET <- function(covstruc,SNPs,fix_omega="full",toler=NULL,TWAS=FALSE,parallel=TRUE,cores=NULL){
  
  time<-proc.time()

  #determine operating system; relevant for parallel runs
  Operating <- Sys.info()[['sysname']]
  
  #pull the betas and SEs from the sumstats output
  beta_SNP <- SNPs[,grep("beta.",fixed=TRUE,colnames(SNPs))]
  SE_SNP <- SNPs[,grep("se.",fixed=TRUE,colnames(SNPs))]
  
  #pull variance from output and save relevant columns
  if (TWAS) {
    SNPs$Gene <- as.character(SNPs$Gene)
    SNPs$Panel <- as.character(SNPs$Panel)
    
    #Gene variance = GREML h2 of gene from FUMA
    varSNP <- SNPs$HSQ
    SNPs <- SNPs[,1:3]
  } else {
    SNPs$A1 <- as.character(SNPs$A1)
    SNPs$A2 <- as.character(SNPs$A2)
    SNPs$SNP <- as.character(SNPs$SNP)
    
    #SNP variance = 2pq
    varSNP <- 2*SNPs$MAF*(1-SNPs$MAF)
    
    SNPs <- SNPs[,1:6]
    
  }
  
  #set sampling variance of SNP (or Gene) to small number
  #reflect near population fixed estiamte
  varSNPSE2 <- (.0005)^2
  
  #pull V, S, and I from LDSC output
  V_LD <- as.matrix(covstruc[[1]])
  S_LD <- as.matrix(covstruc[[2]])
  I_LD <- as.matrix(covstruc[[3]])
  
  #save the row/column positions in I; used for expanding S and V in later code
  coords <- which(I_LD != 'NA', arr.ind= T)
  
  #set univariate intercepts to 1 if estimated below 1
  diag(I_LD) <- ifelse(diag(I_LD)<= 1, 1, diag(I_LD))
    
  #number of phenotypes
  n_phenotypes <- ncol(beta_SNP)
  
  #create shell of full omega matrix
  fix_omegaFull <- diag(n_phenotypes+1)
  
  # If fixed omega matrix, recode - edges fixed to zero = 0; edge free = 1 (other integers encode equality constraints)
  #then expand to include SNP effects
  if (is.matrix(fix_omega)) {
    fix_omega[fix_omega != 0] <- 1
    diag(fix_omega) <- 0
    
    ##add the omega provided by user
    fix_omegaFull[(2:(n_phenotypes+1)),(2:(n_phenotypes+1))] <- fix_omega
    
    #create the vector for the SNP portion of omega
    omega_SNP<-c(0,rep(1,n_phenotypes))
    
    ##add in omega_SNP as first row/column of the expanded omega matrix
    fix_omegaFull[1:(n_phenotypes+1),1] <- omega_SNP
    fix_omegaFull[1,1:(n_phenotypes+1)] <- t(omega_SNP)
    print("Estimating a sparse network (some edges between traits are fixed to zero)")
  }else{
    fix_omegaFull = "full"
    print("Estimating a saturated network (all edges between traits are estimated)")
  }
  
  if(TWAS){
    print("Starting TWAS Network Estimation")
  } else {
    print("Starting GWAS Network Estimation")
  }
  
  #run in serial 
  if(!parallel){
 
    for (i in 1:nrow(beta_SNP)) {
      
      if(i == 1){
        cat(paste0("Running Network: ", i, "\n"))
      }else{
        if(i %% 1000==0) {
          cat(paste0("Running Network: ", i, "\n"))
        }
      }
      
     SNPrun   <- .gwasNET_main(i, cores=1, n_phenotypes, 1, I_LD, V_LD, S_LD, varSNPSE2, SNPs, beta_SNP, SE_SNP, varSNP, TWAS, toler,fix_omegaFull,coords)
     if(i == 1){
       #create empty data.frame to store results
       SNPnet <- as.data.frame(matrix(NA,ncol=ncol(SNPrun),nrow=nrow(SNPs)))
       colnames(SNPnet)<-colnames(SNPrun)
       }
      #store results from each run 
      SNPnet[i,]<-SNPrun[1,]
      
    }
}
 
  #run in parallel
  if(parallel){
    
  if(is.null(cores)){
    ##if no default provided use 1 less than the total number of cores available so your computer will still function
    int <- min(c(nrow(SNPs), detectCores() - 1))
  }else{
    if (cores > nrow(SNPs))
      warning(paste0("Provided number of cores was greater than number of SNPs, reverting to cores=",nrow(SNPs)))
    int <- min(c(cores, nrow(SNPs)))
  }
    
    ##specify the cores should have access to the local environment
    if (Operating != "Windows") {
      cl <- makeCluster(int, type="FORK")
    } else {
      cl <- makeCluster(int, type="PSOCK")
    }
    
    #register the cluster and terminate when code completes
    registerDoParallel(cl)
    on.exit(stopCluster(cl))
  
  #split the V_SNP and S_SNP matrices across cores to run these batches in parallel
  SNPs <- suppressWarnings(split(SNPs,1:int))
  beta_SNP <- suppressWarnings(split(beta_SNP,1:int))
  SE_SNP <- suppressWarnings(split(SE_SNP,1:int))
  varSNP <- suppressWarnings(split(varSNP,1:int))
  
  #separate parallel functions so that it will work on different OS
  if (Operating != "Windows") {
    SNPnet <- foreach(n = icount(int), .combine = 'rbind') %:%
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = "psychonetrics") %dopar% 
      .gwasNET_main(i, int, n_phenotypes, n, I_LD, V_LD, S_LD, varSNPSE2, SNPs[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], TWAS, toler,fix_omegaFull,coords)
  } else {
    #Util-functions have to be explicitly passed to the analysis function in PSOCK cluster
    utilfuncs <- list(); utilfuncs[[".get_V_SNP"]] <- .get_V_SNP; utilfuncs[[".get_Z_pre"]] <- .get_Z_pre; utilfuncs[[".get_V_full"]] <- .get_V_full; utilfuncs[[".sandwichSE"]] <- .sandwichSE
    
    SNPnet <- foreach(n = icount(int), .combine = 'rbind') %:%
      foreach (i=1:nrow(beta_SNP[[n]]), .combine='rbind', .packages = c("psychonetrics", "gdata"),
               .export=c(".gwasNET_main")) %dopar% {
       .gwasNET_main(i, int, n_phenotypes, n, I_LD, V_LD, S_LD,  varSNPSE2,SNPs[[n]], beta_SNP[[n]], SE_SNP[[n]], varSNP[[n]], TWAS, toler,fix_omegaFull,coords,utilfuncs)
               }
  }
  
  ##sort results so it is in order of the output lists provided for the function
  SNPnet <-  SNPnet[order(SNPnet$i),]

  }
  
  SNPnet$i <- NULL
  rownames(SNPnet)<-NULL
  
  time_all <- proc.time()-time
  print(time_all[3])
  
  return(SNPnet)
    
}
