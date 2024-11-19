### twas wrapper
twasNET <- function(covstruc,genes,fix_omega="full",toler=NULL,parallel=TRUE,cores=NULL){
  
  gwasNet(covstruc,SNPs=genes,fix_omega=fix_omega,toler=toler,TWAS=TRUE,parallel=parallel,cores=cores)
}
