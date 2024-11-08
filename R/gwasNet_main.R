.gwasNet_main<-function(i, cores, k, n, I_LD, V_LD, S_LD, varSNPSE2, SNPs, beta_SNP, SE_SNP, varSNP, TWAS, toler,fix_omega,coords,utilfuncs=NULL){
  
  # utilfuncs contains utility functions to enable this code to work on PSOC clusters (for Windows)
  if (!is.null(utilfuncs)) {
    for (j in names(utilfuncs)) {
      assign(j, utilfuncs[[j]], envir=environment())
    }
  }
  
  #determine the Z-statistics of each cell prior to any potential smoothing
  Z_pre <- .get_Z_pre(i, beta_SNP, SE_SNP, I_LD)
  
  #create the SNP portion of the V matrix
  V_SNP <- .get_V_SNP(SE_SNP, I_LD, varSNP, coords, k, i)

  #append the SNP-V matrix with the larger V matrix
  V_Full <- .get_V_full(k, V_LD, varSNPSE2, V_SNP)
  
  #smooth V if necessary
  if(eigen(V_Full)$values[nrow(V_Full)] <= 0){
    V_Full <- as.matrix((nearPD(V_Full, corr = FALSE))$mat)
    V_smooth <- 1
  }
  
  #create empty vector for S_SNP
  S_SNP <- vector(mode="numeric",length=k+1)
  
  #enter SNP variance from reference panel as first observation
  S_SNP[1] <- varSNP[i]
  
  #enter SNP covariances (standardized beta * SNP variance from refference panel)
  for (p in 1:k) {
    S_SNP[p+1] <- varSNP[i]*beta_SNP[i, p]
  }
  
  #create shell of the full S (observed covariance) matrix
  S_Full <- diag(k+1)
  
  ##add the LD portion of the S matrix
  S_Full[(2:(k+1)),(2:(k+1))] <- S_LD
  
  ##add in observed SNP variances as first row/column
  S_Full[1:(k+1),1] <- S_SNP
  S_Full[1,1:(k+1)] <- t(S_SNP)
  
  ##smooth to near positive definite if either V or S are non-positive definite
  ks <- nrow(S_Full)
  if(eigen(S_Full)$values[ks] <= 0){
    S_Full <- as.matrix((nearPD(S_Full, corr = FALSE))$mat)
    S_smooth <- 1
  }
  
  #save difference in Z-statis pre and post smoothing
  if(exists("S_smooth") | exists("V_smooth")){
    SE_smooth <- matrix(0, ks, ks)
    SE_smooth[lower.tri(SE_smooth,diag=TRUE)]  <- sqrt(diag(V_Full))
    Z_smooth <- (S_Fullrun/SE_smooth)[2:ks,1]
    Z_smooth <- max(abs(Z_smooth-Z_pre))
  }else{
    Z_smooth <- 0
  }
  
  #name the columns
  if(!(TWAS)){
    colnames(S_Full) <- c("SNP", colnames(S_LD))
  } else {
    colnames(S_Full) <- c("Gene", colnames(S_LD))
  }
  
  ##name rows like columns
  rownames(S_Full) <- colnames(S_Full)

  ### run the GGM 
  model <- varcov( type = "ggm", covs = S_Full, omega = fix_omega, nobs = 200, covtype = "ML", estimator = "ML", optimizer = "nlminb")
  Model_Results <- runmodel(model)
  
  #get sandwich corrected SEs
  SE<-.sandwichSE(Model_Results,V_Full,toler)
  
  ### extract results and calculate Z and p
   params <- data.frame(Model_Results@parameters)
  params$se <- SE
  params$Zstat <- params$est / params$se
  params$p <- 2 * pnorm(abs(params$est / params$se), lower.tail = FALSE)
  
  if(TWAS){
    params<-subset(params,(params$var1 == "Gene" | params$var2 == "Gene") & !(params$var1 == "Gene" & params$var2 == "Gene"))
    
    #save relevant columns
    params <- params[,c("var1","est","se","Zstat","p")]
    colnames(params) <- c("Trait","Genepcor_est","Genepcor_se","Zstat","p")
    
    #go from long to wide format
    params <- params  %>%
      pivot_wider(names_from = Trait, values_from = c("Genepcor_est", "Genepcor_se", "Zstat", "p"))
  }else{
  params<-subset(params,(params$var1 == "SNP" | params$var2 == "SNP") & !(params$var1 == "SNP" & params$var2 == "SNP"))

  #save relevant columns
  params <- params[,c("var1","est","se","Zstat","p")]
  colnames(params) <- c("Trait","SNPpcor_est","SNPpcor_se","Zstat","p")
  
  #go from long to wide format
  params <- params  %>%
    pivot_wider(names_from = Trait, values_from = c("SNPpcor_est", "SNPpcor_se", "Zstat", "p"))
  }

  # Split the column names by "_", get the last element of each, and order the column names based on these
  ordered_cols <- names(params)[order(sapply(strsplit(names(params), "_"), function(x) tail(x, n = 1)))]
  
  #Restructure parameter estimaets so the set of four estimates is in order for each trait
  params <- params %>% select(all_of(ordered_cols))
  
  #combine the SNP (or Gene) level information, Z_smooth, and the parameter estimates
  output<- cbind(n + (i-1) * cores, SNPs[i,], Z_smooth, params,row.names=NULL)
  
  if(TWAS){
  colnames(output)<-c("i", "Gene","Panel","HSQ","Z_smooth",colnames(params))
  }else{
    colnames(output)<-c("i", "SNP", "CHR", "BP", "MAF", "A1", "A2","Z_smooth",colnames(params))
  }
  
  return(output)
  
}
