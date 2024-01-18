.runGGM <- function(covstruc,fix_omega="full",toler) {
  
  # read in V (sampling covariance) and S (covariance) matrices
  V_LD <- as.matrix(covstruc[[1]]) 
  S_LD <- as.matrix(covstruc[[2]])
  
  
  # If fixed omega matrix, recode - edges fixed to zero = 0; edge free = 1 (other integers encode equality constraints)
  if (is.matrix(fix_omega)) {
    fix_omega[fix_omega != 0] <- 1
    diag(fix_omega) <- 0
  }
  
  
  ### run the GGM
  model <- varcov( type = "ggm", covs = S_LD, omega = fix_omega, nobs = 100, covtype = "UB", estimator = "ML", optimizer = "nlminb",covtype="ML")
  Model_Results <- runmodel(model)
  
  
  ### get necessary model matrices / information
  omega <- getmatrix(Model_Results, "omega") #weight matrix
  delta <- getmatrix(Model_Results, "delta") #scaling matrix
  sigma <- getmatrix(Model_Results, "sigma") #mod-implied cov matrix
  # extra matrices
  mat <- lapply(Model_Results@extramatrices, as.matrix)
  L <- mat$L
  D <- mat$D #duplication matrix
  Dstar <- mat$Dstar
  In <- mat$In #k x k identity matrix
  A <- mat$A
  
  
  ### wls weight matrix from stage 2
  sigma_inv <- solve(sigma)
  S2.W <- 0.5 * t(D) %*% kronecker(sigma_inv, sigma_inv) %*% D
  
  
  ### Jacobian matrix of sigma in respect to theta (all model parameters) i.e "Delta" in genomic sem
  # Note: NEED TO CHECK ORDERING???
  IminOinv <- solve(diag(ncol(omega)) - omega) #inverse of (identify matrix - omega)
  delta_IminOinv <- delta %*% IminOinv
  
  # model derivatives sigma in respect to omega
  d_sigma_omega <- L %*% (delta_IminOinv %x% delta_IminOinv) %*% Dstar
  
  # model derivatives sigma in respect to delta
  d_sigma_delta <- L %*% ((delta_IminOinv %x% In) + (In %x% delta_IminOinv)) %*% A
  
  # full jacobian
  S2.delt <- cbind(d_sigma_omega, d_sigma_delta)
  
  
  ### sandwich correction
  #the "bread" part of the sandwich is the naive covariance matrix of parameter estimates that would only be correct if the fit function were correctly specified
  bread <- solve(t(S2.delt) %*% S2.W %*% S2.delt, tol = toler)
  
  #create the "lettuce" part of the sandwich
  lettuce <- S2.W %*% S2.delt
  
  #ohm-hat-theta-tilde is the corrected sampling covariance matrix of the model parameters
  Ohtt <- bread %*% t(lettuce) %*% V_LD %*% lettuce %*% bread
  
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.vector(sqrt(diag(Ohtt)))
  
  
  ### extract results
  results <- data.frame(Model_Results@parameters)[, c("var1","op","var2","est","se","p","matrix","par")]
  results$se <- SE
  results$se[results$par == 0] <- NA
  results$p <- 2 * pnorm(abs(results$est / results$se), lower.tail = FALSE)
  
  return(results)
}
