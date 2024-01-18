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
  model <- varcov( type = "ggm", covs = S_LD, omega = fix_omega, nobs = 200, covtype = "ML", estimator = "ML", optimizer = "nlminb")
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
  params <- data.frame(Model_Results@parameters)
  params$se <- SE
  params$se[params$par == 0] <- NA
  params$z <- params$est / params$se
  params$p <- 2 * pnorm(abs(params$est / params$se), lower.tail = FALSE)
  params <- params[,c("var1","op","var2","est","se","z","p","par","matrix","row","col")]
  colnames(params) <- c("trait1","op","trait2","est","se","z","p","free","matrix","row","col")
  
  #calculate model fit if there are pruned edges (otherwise fully saturated and not relevant)
  if(is.matrix(fix_omega)){

  #calculate model chi-square:
    
  #eigen values of the V matrix  
  Eig<-as.matrix(eigen(V_LD)$values)
  Eig2<-diag(ncol(V_LD))
  diag(Eig2)<-Eig
  
  #Pull P1 (the eigen vectors of V_eta)
  P1<-eigen(V_LD)$vectors
  
  #residual matrix: difference between model implied matrix and observed matrix
  resid<-S_LD-sigma
  
  #eta: the vector of unique elements of the residual matrix
  eta<-as.vector(lowerTriangle(resid,diag=TRUE))
  
  #matrix algebra weighting the vector of residuals by the precision of those residuals (i.e., P1 and Eig)
  model_chi<-t(eta)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta

  #degrees of freedom of model (how many unique edges are fixed to 0)
  df<-Model_Results@fitmeasures$df

  #calculate p-value for model chi-square
  model_chi_p<-pchisq(model_chi,df,lower.tail=FALSE)
  
  #calculate SRMR
  
  #unique values of genetic correlation  matrix
  obs <-  cov2cor(S_LD)[!lower.tri(cov2cor(S_LD))]
  
  #unique values of model implied correlation matrix
  imp <-  cov2cor(sigma)[!lower.tri(cov2cor(sigma))]
  
  #square root of the average squared residual 
  SRMR<-sqrt(mean((imp - obs)^2))
  
  #calculate CFI:
  
  #the difference between the model implied matrix of the independence model [only models variances]
  #and the observed matrix is the observed matrix with diagonals set to 0
  resid_CFI<-S_LD
  diag(resid_CFI)<-0
  
  #eta: the vector of unique elements of the residual matrix
  eta_CFI<-as.vector(lowerTriangle(resid_CFI,diag=TRUE))
  
  #matrix algebra weighting the vector of residuals by the precision of those residuals (i.e., P1 and Eig)
  CFI_chi<-t(eta_CFI)%*%P1%*%solve(Eig2)%*%t(P1)%*%eta_CFI
 
  ##df of independence Model
  k<-ncol(S_LD)
  dfCFI <- (((k * (k + 1))/2) - k)
  
  #calculate CFI
  CFI<-as.numeric(((CFI_chi-dfCFI)-(model_chi-df))/(CFI_chi-dfCFI))

  #calculate AIC
  AIC<-(model_chi + 2*Model_Results@fitmeasures$df)
  
  #combine model fit indices
  modelfit<-cbind(model_chi,df,model_chi_p,AIC, SRMR,CFI)
  colnames(modelfit)=c("model_chisquare","df","modelchi_pvalue","AIC", "SRMR", "CFI")
  
  }else{
    modelfit <- NULL
  }

  traitnames <- list(colnames(S_LD),colnames(S_LD))
  dimnames(omega) <- traitnames
  dimnames(delta) <- traitnames
  dimnames(sigma) <- traitnames
  
  return(list(parameters=params,modelfit=modelfit,omega=omega,delta=delta,sigma=sigma))
}
