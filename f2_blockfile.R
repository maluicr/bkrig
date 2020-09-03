
# ----- description -----

# function for block data

# as input you should provide a georeferenced grid file
# with id region values at simulation locations, and 
# a list returned by funtion irates().

# ----- arguments ------

# rateobj, character, name of list, output of function irates() 
# gridimage, character,  name of georeferenced grid file (e.g. tif) 
# na.value, numeric, integer with grid value for "No data"  
 
# ------ function ------

blockfile = function (rateobj, gridimage, na.value = -999){
  
  if (!"raster" %in% installed.packages()) install.packages("raster")
  if (!"rgdal" %in% installed.packages()) install.packages("rgdal")
  if (!"sp" %in% installed.packages()) install.packages("sp")

  library(raster)
  library(rgdal)
  library(sp)
  
  # strings w/ path to store files
  day = as.character(rateobj[["file"]]["day"])
  folder = as.character(rateobj[["file"]]["folder"])
  
  grd = raster(gridimage)
  
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
  nas = na.value
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
  fblk = paste0(folder, "/", day, blk_nameO)
  
  if (file.exists(fblk)){
    file.remove(fblk)  
  }
  
  # create notification file
  file.create(fblk)
  
  # write file header
  cat(blk_name, file = fblk, sep="\n")
  cat(massn, file = fblk, sep="\n", append = T)
  
  # write block id, rate, error and x, y, z
  for (i in 1 : massn){
    id = massid[i]
    # store nr blocks
    nc = sum(gridout[,,1] == id)
    
    # get vectors from irates function  
    oid = rateobj[["rates"]][, "id"]
    rate = rateobj[["rates"]][, "rate"]
    err = rateobj[["rates"]][, "err"]
    dftab = data.frame(oid, rate, err)
    
    # store rate and error
    dft = dftab[dftab$oid == id, c("rate","err")]
    rt = as.numeric(dft[1])
    er = as.numeric(dft[2])
    
    # write block id
    cat(paste0("Block#", id), file = fblk, sep="\n", append = T)
    
    # write rate, error & nr blocks
    cat(rt, file = fblk, sep="\n", append = T)
    cat(er, file = fblk, sep="\n", append = T)
    cat(nc, file = fblk, sep="\n", append = T)
    
    # data.table::fwrite
    blk_list = list()
    if (!"data.table" %in% installed.packages()) install.packages("data.table")
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
  listgrid = grd
  listgridpars = list(nodes = c(nx, ny), resolution = c(resx, resy), origin = c(ox, oy), NAs = na.value)
  listfile = list(day = day, name = paste0(day, blk_nameO), folder = folder)
  listgridout = list(values = stacf, idblock = massid, nblock = massn)
  return(list(gridpars = listgridpars, outgrid = listgridout , file = listfile, ingrid = listgrid))
}
