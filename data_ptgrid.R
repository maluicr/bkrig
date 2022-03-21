

# author: m ribeiro; 
# date : 18-03-22; last revision : 18-03-2022
# email: manuel.ribeiro@tecnico.ulisboa.pt

# description:  
# data.frame to include in folder 'data' of blockdss package.
# data.frame is obtained from grid2k with 3 columns : x, y, grid2k 

# notes: 
# - grid2k values are oid_ (municipalities) values
# - spatial resolution is 2000 m

library(raster)
library(sp)

par(mfrow= c(1,1))
r <- raster("grid2k.tif")
plot(r)

xy <- sp::SpatialPixels(SpatialPoints(coordinates(r)))@coords
oid_ <- getValues(r)
ptgrid <- data.frame(x=xy[,1],y = xy[,2], oid_)
coordinates(ptgrid) = ~x+y
class(ptgrid)

# gridded data (pixels) 
gridded(ptgrid) = TRUE
class(ptgrid)
spplot(ptgrid, main = "municipality id's")


# code included in blockfile()
# hereafter

grd <- raster::raster(ptgrid)

# grid parameters
ny = grd@nrows
nx = grd@ncols
resx = raster::res(grd)[1]
resy = raster::res(grd)[2]
ox = raster::extent(grd)[1]
oy = raster::extent(grd)[3]

# create matrix to store block coordinates
matA = raster::as.matrix( grd, nrow = ny, ncol = nx )
matA2 = matrix(-999, nrow = ny, ncol = nx)

# fill matrix with block coordinates
for (i in 1:ny){
  matA2[i,] <- c(matA[ny-i+1,] )
}

# convert to vector str
matA3 = as.data.frame(t(matA2))
stac3 = raster::stack(matA3)
stacf = stac3$values

# create array for blockfile
grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = F)
grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, by = resx), y = seq(oy, oy + (ny-1) * resy, resy))
grx = grxy[,1]
gry = grxy[,2]
gridout = array (c (grdata, grx, gry), dim =c(ny, nx , 3))
gridout = round(gridout, 4)

# write blocksfile for dss

# store vector grid id
massid = raster::unique(grd)

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
  #if (!"data.table" %in% installed.packages()) install.packages("data.table")
  #library("data.table")
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

