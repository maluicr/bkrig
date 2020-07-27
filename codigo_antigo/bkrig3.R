library(gstat)
library(sp)
library(rgdal)
library(tidyverse)

# ---- folder to store results -----

day = 20200601
# create folder if exists is FALSE
sid = paste0(day)
if(!file.exists(sid)) dir.create(sid, recursive=T)  

# ----- shapefile cmassa concelhos pop 2018 -----

shp = readOGR(dsn = "gis","cc_pt_cmassa_fr_popres",encoding = "UTF-8", use_iconv = T)
names(shp@data) = c("x", "y", "objectid", "nome_cc", "cod",  "popres18")

# dss need z coordinate 3d
shp@data$z = 0

pop = as.data.frame(shp@data)
pop$popres18 = as.numeric(as.character(pop$popres18))
pop$objectid = as.numeric(as.character(pop$objectid))

# ----- grid pt tif to out -----

library(raster)
library(rgdal)
library(sp)

grd = raster("gis/grid2k.tif")

# grid parameters
ny = grd@nrows
nx = grd@ncols
resx = res(grd)[1]
resy = res(grd)[2]
ox = extent(grd)[1]
oy = extent(grd)[3]

# gslib format
matA = as.matrix( grd, nrow = ny, ncol = nx )
matA2 = matrix(-999, nrow = ny, ncol = nx)

for (i in 1:ny){
  matA2[i,] <- c(matA[ny-i+1,] ) 
}

matA3 = as.data.frame(t(matA2))
stac3 = stack(matA3)
stacf = stac3$values
stacf[is.na(stacf)] = -999

# for dss.exe format
grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = F)
grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, by = resx), y = seq(oy, oy + (ny-1) * resy, resy))
grx = grxy[,1]
gry = grxy[,2]

gridout = array (c (grdata, grx, gry), dim =c(ny, nx , 3))

# write.table(stacf, file = "grid2k.out", sep = "\t",row.names = FALSE, col.names=FALSE)

# ----- observed covid19 data -----

# dgs
dgs = read.table("dgsData20200601.txt", header = T, sep = ",", dec = ".", strip.white = T )
names(dgs) = c("nome_cc", "casos")

# merge with pop data
tab_all = merge(x = pop, y = dgs, by = c("nome_cc"), all.x = TRUE)

# crude rates
tab_all$rate = tab_all$casos / tab_all$popres18

# m = population-weighted mean of the n rates 
pop_tot = sum(tab_all$popres18, na.rm = T)
m = sum(tab_all$rate * tab_all$popres18, na.rm = T) / pop_tot

# error variance term (m/n_i)
tab_all$error = m / tab_all$popres18

# NA cases set to 2 cases
tab_all$casos = ifelse(is.na(tab_all$casos), 2, tab_all$casos)

# recalculate crude rates
tab_all$rate = tab_all$casos / tab_all$popres18

tabela = tab_all[, c("objectid", "x", "y", "z", "rate", "error")]

# ------ write file with counts, x,y,z ------

# prepare data to write file 
tabnotf = tabela[, c("x", "y", "z", "rate")]

nvars = ncol(tabnotf)
nobs = nrow(tabnotf)
namevars = names(tabnotf)

# create file
name = "notified"

if (file.exists(paste0(sid , "/", sid, name, ".out"))){
  file.remove(paste0(sid , "/", sid, name, ".out"))  
}

# start_time = Sys.time()
notf = file(paste0(sid , "/", sid, name, ".out"), open = "w")

# write filename
writeLines(name, con = notf, sep = "\n")

# write nr of variables
write(nvars, file = notf, append = T)

# write name of variables
for (i in 1:nvars) {
  write(namevars[i], file = notf, append = T)
}

# write data
for (i in 1:nobs){
  vec = as.numeric(tabnotf[i,])
  write(vec, file = notf, append = T)
}

close(notf)
closeAllConnections()

# ------ write file for dss ------

massid = unique(grd)
massn = length(massid)


# create file
name = "blockdata"

if (file.exists(paste0(sid , "/", sid, name, ".out"))){
  file.remove(paste0(sid , "/", sid, name, ".out"))  
}

# start_time = Sys.time()

logf = file(paste0(sid , "/", sid, name, ".out"), open = "w")

# write filename
writeLines(name, con = logf, sep = "\n")

# write nr of blocks
write(massn, file = logf, append = T)

# write block id, rate, error and x, y, z
for (i in 1 : massn){
  id = massid[i]
  # store nr blocks
  nc = sum(gridout[,,1] == id)
  
  # store rate and error
  t = tabela[tabela$objectid == id, c("rate", "error")]
  rt = as.numeric(t[1])
  er = as.numeric(t[2])
  
  # write block id
  writeLines(paste0("Block#", id), con = logf, sep = "\n")
  
  # write rate and error
  write(rt, file = logf, append = T)
  write(er, file = logf, append = T)
  
  # write nr blocks
  write(nc, file = logf, append = T)
  
  # write x, y, z (=0) coords block
  for (k in 1:ny) {
    for(l in 1:nx){
      if(gridout[k, l, 1] == id) {
        d = c(gridout[k,l,c(2,3)],0)
        write(d, file = logf, append = T)
      }
    }
  }
}

close(logf)
closeAllConnections()

# end_time = Sys.time()
# end_time - start_time
