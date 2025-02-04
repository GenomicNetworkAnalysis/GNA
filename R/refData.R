refData <- function(data, dir = getwd()){
  if (dir.exists(dir)){
    
    if (data=="eur" | data=="eas") {
      # ld scores
      file.copy(system.file("extdata", paste0(data,"_w_ld_chr"), package="GNA"), dir, recursive=TRUE)
      # hm3 snplist
      file.copy(system.file("extdata", "w_hm3.snplist", package="GNA"), dir, recursive=TRUE)
      # 1000g ref
      file.copy(list.files(system.file("extdata", paste0("reference_1000g_",data), package="GNA"), full.names = T), dir, recursive=TRUE)
      lapply(list.files(dir, full.names = T, pattern = ".gz"), R.utils::gunzip) 
      file.append(paste0(dir,"/reference.1000G.maf.0.005.",data,".txt"),paste0(dir,"/reference.1000G.chr",1:23,".maf.0.005.",data,".txt"))
      file.remove(paste0(dir,"/reference.1000G.chr",1:23,".maf.0.005.",data,".txt"))
      
    } else if (data=="example") {
      d <- file.path(dir,"example_data")
      dir.create(d, showWarnings = F)
      file.copy(list.files(system.file("extdata", "example_data", package="GNA"), full.names = T), d, recursive=TRUE)
      
    } else {
      print("data argument must be 'eur', 'eas' or 'example'")
    }
    
  } else {
    print(paste0("directory '",dir,"' does not exist"))
  }
}
