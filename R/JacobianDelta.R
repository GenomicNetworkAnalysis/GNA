#### COMPUTE JACOBIAN DELTA MATRIX (for scaling of V matrix)
## adapted from Psychonetrics - needs to be cleaned

# model derivatives sigma in respect to omega
.d_sigma_omega <- function(L,delta_IminOinv,Dstar){
  res <- L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
  as.matrix(res)
}

# model derivatives sigma in respect to delta
.d_sigma_delta <- function(L,delta_IminOinv,In,A){
  res <- L %*% ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A
  as.matrix(res)
}


# Full Jacobian matrix ("Delta" in genomic sem code)
.d_jacobian_delta <- function(Model_Results){
  
  # extract model matrices / information
  omega <- getmatrix(Model_Results,"omega") #weight matrix
  delta <- getmatrix(Model_Results,"delta") #scaling matrix
  IminOinv <- solve(diag(ncol(omega)) - omega) #inverse of (identify matrix - omega)
  delta_IminOinv <- delta %*% IminOinv
  L <- as.matrix(Model_Results@extramatrices$L)
  Dstar <- as.matrix(Model_Results@extramatrices$Dstar)
  In <- as.matrix(Model_Results@extramatrices$In)
  A <- as.matrix(Model_Results@extramatrices$A)
  
  # !!can probably remove below as unlikely to be relevant in GeneNets???
  y <- Model_Results@submodel #"ggm"
  corinput <- Model_Results@sample@corinput #if input was correlation
  meanstructure <- Model_Results@meanstructure #logical (FALSE)
  mu <- getmatrix(Model_Results,"mu") #means 
  ifelse("tau" %in% names(Model_Results@modelmatrices$fullsample),
         tau <- getmatrix(Model_Results,"tau"), #tau (thresholds)
         tau <- matrix(NA,1,nrow(omega)))
  
  # Number of variables:
  nvar <- nrow(omega)
  
  # Number of means/thresholds:
  nMean_Thresh <- sum(!is.na(tau)) + sum(!is.na(mu))
  nThresh <- sum(!is.na(tau))
  
  # Number of observations:
  nobs <- nMean_Thresh + # Means
    (nvar * (nvar+1))/2 # Variances
  
  # Number of parameters is less if corinput is used or if meanstructure is ignored:
  npars <- nobs - corinput * nvar - (!meanstructure) * sum(!is.na(mu))
  
  # Mean part:
  meanPart <- seq_len(nMean_Thresh)
  
  # Variance part:
  varPart <- max(meanPart) + seq_len(nvar*(nvar+1)/2)    
  
  # Var part for parameters:
  varPartPars <- meanstructure * max(meanPart) + nThresh +  seq_len(nvar*(nvar+1)/2)    
  
  
  # Empty Jacobian:
  Jac <- matrix(0, nobs, npars)
  
  if (meanstructure || nThresh > 0){
    # Fill mean part with diagonal:
    Jac[meanPart,meanPart] <- as.matrix(Diagonal(nMean_Thresh))
  }
  
  
  # Now fill the sigma part:
  if (y == "ggm"){
    # Gaussian graphical model:
    netPart <- meanstructure*max(meanPart) + nThresh + seq_len(nvar*(nvar-1)/2)
    scalingPart <- max(netPart) + seq_len(nvar)
    
    Jac[varPart,netPart] <- .d_sigma_omega(L,delta_IminOinv,Dstar)
    Jac[varPart,scalingPart] <- .d_sigma_delta(L,delta_IminOinv,In,A)       
  }
  
  # Cut out the rows not needed
  # FIXME: Nicer to not have to compute these in the first place...
  if (corinput){
    keep <- c(rep(TRUE,nMean_Thresh),diag(nvar)[lower.tri(diag(nvar),diag=TRUE)]!=1)
    Jac <- Jac[keep,]
  }
  if (!meanstructure){
    
    if (all(is.na(tau))){
      Jac <- Jac[-(seq_len(nvar)), ] 
    } else if (any(colSums(is.na(tau)) == nrow(tau))) stop("Mix of continuous and ordinal variables is not yet supported.")
  }
  
  return(Jac)
}
