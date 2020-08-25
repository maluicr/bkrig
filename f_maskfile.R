maskfile = function(blockobj){
  
  obj = unlist(blockobj[["outgrid"]]["values"], use.names = F)
  na = unlist(blockobj[["gridpars"]]["NAs"], use.names = F)
  day = as.character(blockobj[["file"]]["day"])
  folder = as.character(blockobj[["file"]]["folder"])
    
  # write maskfile for dss
  
  # create mask vector
  mask = ifelse(obj == na, -1, 0)
  val = unique(mask)
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
  
  listf = list(day = day, name = paste0(day, msk_nameO), folder = folder)
  listz = list(nzones = mask_zones, zoneval = val)
  
  return(list(file = listf, zones = listz))
}


