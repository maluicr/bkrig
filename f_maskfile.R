maskfile = function(blockobj){
  
  obj = blockobj[["outfile"]]
  na = blockobj[["gridpars"]]["NAs"]
  day = as.character(blockobj[["block.file"]]["day"])
  folder = as.character(blockobj[["block.file"]]["folder"])
    
  # write maskfile for dss
  
  # create mask vector
  mask = ifelse(obj == na, -1, 0)
  mask_zones = length(unique(mask))
  
  # prepare data to write file 
  nvars = 1
  namevars = "values"
  nval = length(mask)
  
  # create file path
  msk_name = "mask"
  msk_nameO = paste0(msk_name, ".out")
  fmsk = paste0(folder, "/", day, msk_nameO)
  
  if (file.exists(fmsk)){
    file.remove(fmsk)
  }
  
  # create mask file
  file.create(fmsk)
  
  # write file header
  cat(msk_name, file = fmsk, sep="\n")
  cat(nvars, file = fmsk, sep="\n", append = T)
  cat(namevars, file = fmsk, sep="\n", append = T)
  
  # write mask data
  write.table(mask, file = fmsk, append = T, row.names = F, col.names = F)
  
  
}


