### twas wrapper
twasNET <- function(covstruc,genes,fix_omega="full",toler=NULL,parallel=TRUE,cores=NULL){
  
  gwasNET(covstruc,SNPs=genes,fix_omega=fix_omega,toler=toler,TWAS=TRUE,parallel=parallel,cores=cores)
}
