
blockfile = function (dfrate, gridimage){
  
  require("raster", "rgdal", "sp")
  
  library(raster)
  library(rgdal)
  library(sp)
  
  grd = gridimage
  
  # grid parameters
  ny = grd@nrows
  nx = grd@ncols
  resx = res(grd)[1]
  resy = res(grd)[2]
  ox = extent(grd)[1]
  oy = extent(grd)[3]
  
  # create matrix to store block coordinates
  matA = as.matrix( grd, nrow = ny, ncol = nx )
  matA2 = matrix(-999, nrow = ny, ncol = nx)
  
  # fill matrix with block coordinates
  for (i in 1:ny){
    matA2[i,] <- c(matA[ny-i+1,] ) 
  }
  
  # convert to vector str
  matA3 = as.data.frame(t(matA2))
  stac3 = stack(matA3)
  stacf = stac3$values
  
  # set NAs values to -999
  nas = -999
  stacf[is.na(stacf)] = nas
  
  # create array for blockfile
  grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = F)
  grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, by = resx), y = seq(oy, oy + (ny-1) * resy, resy))
  grx = grxy[,1]
  gry = grxy[,2]
  gridout = array (c (grdata, grx, gry), dim =c(ny, nx , 3))
  gridout = round(gridout, 4)
  
  # write blocksfile for dss
  
  # store vector grid id
  massid = unique(grd)
  
  # store length vector
  massn = length(massid)
  
  # create file path
  blk_name = "blockdata"
  blk_nameO = paste0(blk_name,".out")
  day = dfrate[[2]][1]
  folder = dfrate[[2]][3]
  fblk = paste0(folder, "/", day, blk_nameO)
  
  if (file.exists(fblk)){
    file.remove(fblk)  
  }
  
  # create notification file
  file.create(fblk)
  
  # write file header
  cat(blk_name, file = fblk, sep="\n")
  cat(massn, file = fblk, sep="\n", append = T)
  
  t0 = Sys.time()
  
  # write block id, rate, error and x, y, z
  for (i in 1 : massn){
    id = massid[i]
    # store nr blocks
    nc = sum(gridout[,,1] == id)
    
    # get vectors from irates function  
    oid = dfrate[[1]][, "oid"]
    rate = dfrate[[1]][, "rate"]
    err = dfrate[[1]][, "err"]
    tab = data.frame(oid, rate, err)
    
    # store rate and error
    t = tab[tab$oid == id, c("rate","err")]
    rt = as.numeric(t[1])
    er = as.numeric(t[2])
    
    # write block id
    cat(paste0("Block#", id), file = fblk, sep="\n", append = T)
    
    # write rate, error & nr blocks
    cat(rt, file = fblk, sep="\n", append = T)
    cat(er, file = fblk, sep="\n", append = T)
    cat(nc, file = fblk, sep="\n", append = T)
    
    # data.table::fwrite
    blk_list = list()
    require("data.table")
    library("data.table")
    for (k in 1:ny) {
      for(l in 1:nx){
        if(gridout[k, l, 1] == id) {
          d = paste(gridout[k, l, 2], "\t", gridout[k, l, 3], "\t",  0)
          blk_list[[length(blk_list) + 1]] = list (d)
        }
      }
    }
    data.table::fwrite(blk_list, file = fblk, append = T, sep="\n")
  }
  return(print(paste("File is stored at: ", fblk)))
}
