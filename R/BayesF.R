.BayesF <- function(covstruc, model_out){
  
  # define function to calculate bayes factor .. (JZS Bayes factor - adapted from doi.org/10.3758/s13423-012-0295-x)
  jzs_bayes <- function(r2_0,r2_1,p0,p1,n){
    
    int <- function(r2,p,n,g) {
      f = exp((((n-1-p)/2) * log(1+g)) + ((-(n-1)/2) * log(1+(1-r2)*g)) + (-3/2 * log(g)) + (-n/(2*g)))
      return(ifelse(g == 0, 0, f))
    }
    
     #try catch so that if function returns non-finitive value it returns NA
      tryCatch({
   bf10 <- integrate(int,lower=0,upper=Inf,r2=r2_1,p=p1,n=n)$value / integrate(int,lower=0,upper=Inf,r2=r2_0,p=p0,n=n)$value
    return(bf10)
    }, error = function(e) {
      return(NA)
    })
}

  # apply function across all partial rgs
  p_rg <- na.omit(model_out$parameters[model_out$parameters$matrix == "omega",])
  cors <- cov2cor(covstruc[[2]])
  traits <- colnames(cors)
  out_bf10 <- data.frame()

  for (i in 1:nrow(p_rg)){
    pair <- p_rg[i,]
    y <- pair$trait1
    x <- pair$trait2
    tmp1 <- !traits %in% c(y, x)
    tmp2 <- !traits %in% y
    
    ryx_0 <- cors[tmp1, y]
    rxx_0 <- cors[tmp1, tmp1]
    ryx_1 <- cors[tmp2, y]
    rxx_1 <- cors[tmp2, tmp2]
    
    r2_0 <- as.numeric(t(ryx_0) %*% solve(rxx_0) %*% ryx_0)    #multiple R2 for model 0 (excluding x)
    r2_1 <- as.numeric(t(ryx_1) %*% solve(rxx_1) %*% ryx_1)    #multiple R2 for model 1 (full model)
    p0 <- as.numeric(length(traits)-2)
    p1 <- as.numeric(length(traits)-1)
    n <- as.numeric((((1-(pair$est^2))^2)/(pair$se^2))+length(traits)+1)
    
    BF10 <- jzs_bayes(r2_0, r2_1, p0, p1, n)
    
    out_bf10 <- rbind(out_bf10, c(y,x,r2_0,r2_1,BF10))
    }
  
  colnames(out_bf10) <- c("trait1","trait2","R2_mod0","R2_mod1","BayesFactor10")
  out_bf10[,c("R2_mod0","R2_mod1","BayesFactor10")] <- apply(out_bf10[,c("R2_mod0","R2_mod1","BayesFactor10")], 2, as.numeric)
  
  return(out_bf10)
}
