.sandwichSE<-function(Model_Results,V,toler){
  
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
  Ohtt <- bread %*% t(lettuce) %*% V %*% lettuce %*% bread
  
  #the lettuce plus inner "meat" (V) of the sandwich adjusts the naive covariance matrix by using the correct sampling covariance matrix of the observed covariance matrix in the computation
  SE <- as.vector(sqrt(diag(Ohtt)))
  
  return(SE)
}

.covQC <-function(covstruc,traits){ 
  
##read in the LD portion of the V (sampling covariance) matrix
V_LD<-as.matrix(covstruc[[1]])

##read in the LD portion of the S (covariance) matrix
S_LD<-as.matrix(covstruc[[2]])

##pull the column names specified in the munge function
S_names<-colnames(S_LD)
rownames(S_LD)<-colnames(S_LD)

##name columns of V to remove any variables not used in the current analysis
y<-expand.grid(S_names,S_names)
y<-y[!duplicated(apply(y,1,function(x) paste(sort(x),collapse=''))),]
V_Names<-paste(y$Var1,y$Var2,sep=" ")
colnames(V_LD)<-V_Names
rownames(V_LD)<-V_Names

##determine whether all variables in S are in the traits listed
##if not, remove them from S_LD and V_LD 
remove2<-c()
w<-1

##also for exact cases
for(i in 1:length(S_names)){
  S_names[[i]]<-paste0("\\b", S_names[[i]],"\\b",sep="")
}

for(i in 1:length(S_names)){
  b<-grepl(S_names[i], paste(traits,collapse= ","))
  if(b == FALSE){
    remove<-paste0("\\b", colnames(S_LD)[i],"\\b",sep="")
    remove2[w]<-i
    V_LD <- V_LD[-grep(pattern=remove[1],row.names(V_LD)),-grep(pattern=remove[1],colnames(V_LD))]
    w<-w+1
    if (!(is.matrix(V_LD))) {
      stop("None of the trait names in the LDSC output match names in the traits specified")
    }
  }
}

if(is.null(remove2) == FALSE){
  S_LD<-S_LD[-remove2,-remove2]
}

#define k and z after removing non-used variables
k<-ncol(S_LD)
z<-(k*(k+1))/2

##smooth to near positive definite if either V or S are non-positive definite
S_LDb<-S_LD
smooth1<-ifelse(eigen(S_LD)$values[nrow(S_LD)] <= 0, S_LD<-as.matrix((nearPD(S_LD, corr = FALSE))$mat), S_LD<-S_LD)
LD_sdiff<-max(abs(S_LD-S_LDb))

V_LDb<-V_LD
smooth2<-ifelse(eigen(V_LD)$values[nrow(V_LD)] <= 0, V_LD<-as.matrix((nearPD(V_LD, corr = FALSE))$mat), V_LD<-V_LD)
LD_sdiff2<-max(abs(V_LD-V_LDb))

SE_pre<-matrix(0, k, k)
SE_pre[lower.tri(SE_pre,diag=TRUE)] <-sqrt(diag(V_LDb))

SE_post<-matrix(0, k, k)
SE_post[lower.tri(SE_post,diag=TRUE)] <-sqrt(diag(V_LD))

Z_pre<-S_LDb/SE_pre
Z_post<-S_LD/SE_post
Z_diff<-(Z_pre-Z_post)
Z_diff[which(!is.finite(Z_diff))]<-0
Z_diff<-max(Z_diff)
rm(V_LDb,S_LDb,Z_pre,Z_post)


if(LD_sdiff > 0){
  print(paste("The S matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff, ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
}

if(LD_sdiff > .025){
  warning("A difference greater than .025 was observed pre- and post-smoothing in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
}

if(Z_diff > .025){
  warning("A difference greater than .025 was observed pre- and post-smoothing for Z-statistics in the genetic covariance matrix. This reflects a large difference and results should be interpreted with caution!! This can often result from including low powered traits, and you might consider removing those traits from the model. If you are going to run a multivariate GWAS we strongly recommend setting the smooth_check argument to true to check smoothing for each SNP.")
}

if(LD_sdiff2 > 0){
  print(paste("The V matrix was smoothed prior to model estimation due to a non-positive definite matrix. The largest absolute difference in a cell between the smoothed and non-smoothed matrix was ", LD_sdiff2, "As a result of the smoothing, the largest Z-statistic change for the genetic covariances was ", Z_diff,  ". We recommend setting the smooth_check argument to true if you are going to run a multivariate GWAS.", sep = " "))
}

return(list(V_LD=V_LD,S_LD=S_LD))

}

.get_V_full <- function(k, V_LD, varSNPSE2, V_SNP) {
    ##create shell of full sampling covariance matrix
    V_Full<-diag(((k+1)*(k+2))/2)

    ##input the ld-score regression region of sampling covariance from ld-score regression SEs
    V_Full[(k+2):nrow(V_Full),(k+2):nrow(V_Full)]<-V_LD

    ##add in SE of SNP variance as first observation in sampling covariance matrix
    V_Full[1,1]<-varSNPSE2

    ##add in SNP region of sampling covariance matrix
    V_Full[2:(k+1),2:(k+1)]<-V_SNP
    return(V_Full)
}

.get_V_SNP <- function(SE_SNP, I_LD, varSNP, coords, k, i) {
     V_SNP<-diag(k)
  
      for (p in 1:nrow(coords)) {
        x<-coords[p,1]
        y<-coords[p,2]
        if (x != y) {
          V_SNP[x,y]<-(SE_SNP[i,y]*SE_SNP[i,x]*I_LD[x,y]*sqrt(I_LD[x,x])*sqrt(I_LD[y,y])*varSNP[i]^2)}
        if (x == y) {
          V_SNP[x,x]<-(SE_SNP[i,x]*sqrt(I_LD[x,x])*varSNP[i])^2
        }
      }
  
    return(V_SNP)
}

.get_Z_pre <- function(i, beta_SNP, SE_SNP, I_LD) {

Z_pre<-beta_SNP[i,]/(SE_SNP[i,]*sqrt(diag(I_LD)))

    return(Z_pre)
}


.pruneNet <- function(model_out,p.adjust="fdr",alpha=0.05,prunepower=FALSE){
  
  par <- model_out$parameters[model_out$parameters$matrix=="omega",] #pull parameters table
  #adjusted pvals
  par$p_adj <- p.adjust(par$p, method=p.adjust) 
  #prune non-significant edges
  par$weight <- par$est
  par$weight[par$p_adj >= alpha] <- 0
   
  
  if(prunepower){
   par$weight[par$power < .8] <- 0 
  }
  
  #update omega with pruned weights
  omega <- model_out$omega
  if(!all(par$weight == par$est)){
    for (i in 1:nrow(par)){
      omega[par$row[i],par$col[i]] <- par$weight[i]
      omega[par$col[i],par$row[i]] <- par$weight[i]
    }
  }
  return(omega)
}
            
