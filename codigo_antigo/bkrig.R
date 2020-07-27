library(gstat)
library(sp)
library(rgdal)
library(tidyverse)

# ------------ folder to store results --------------------------------------------
day = 202000430
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
nx = grd@ncols
ny = grd@nrows
resx = res(grd)[1]
resy = res(grd)[2]
ox = extent(grd)[1]
oy = extent(grd)[3]

# gslib format
matA = as.matrix( grd, ncol = nx, nrow = ny)
matA2 = matrix(-999, ncol = nx, nrow = ny)

for (i in 1:ny){
  matA2[i,] <- c(matA[ny-i+1,] ) 
}

matA3 = as.data.frame(t(matA2))
stac3 = stack(matA3)
stacf = stac3$values
stacf[is.na(stacf)] = -999

# for dss.exe format
grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = T)
grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, resx), y = seq(oy, oy + (ny-1) * resy, resy))
grx = grxy[,1]
gry = grxy[,2]

gridout = array (c (grdata, grx, gry), dim =c(nx, ny , 3))

# write.table(stacf, file = "grid2k.out", sep = "\t",row.names = FALSE, col.names=FALSE)

# ----- observed data -----

# dgs
dgs = read.table("20200430_dgs.txt", header = T, sep = "\t", dec = ".", strip.white = T )
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

# ------ write file for dss ------

massid = unique(grd)
massn = length(massid)

# create file
name = "blockdata"
logf = file(paste0(sid , "/", sid, name, ".txt"),"w")
sink(logf, type = "output")

# write filename
writeLines(name, con = logf, sep = "\n")

# write nr of blocks
write(massn, file = logf, sep = "\n")

# create vector to store rate error by block id
t =c()

for (i in 1 : massn){
  
  # store nr blocks
  nc = sum(gridout[,,1] == massid[i])
  
  # write block id
  writeLines(paste0("Block#", massid[i]), con = logf, sep = "\n")
  
  # fill vector t
  t = tabela %>% dplyr::filter(objectid == massid[i]) %>%
    select("rate", "error") %>% as.numeric(t)
  
  # write vector t
  write(t[1], file = logf, sep = "\n")
  write(t[2], file = logf, sep = "\n")
  
  # write nr blocks
  write(nc, file = logf, sep = "\n")
  
}

close(logf)
closeAllConnections()



for (i in 1:massn){
  nc = sum(gridout[,,1] == massid[i])
  nbloco[i,] = c(massid[i], nc)
  }



for (k in 1:)
gridout[,,1] == massid[i] 

for (k in 1:nx) {
  for(l in 1:ny){
   if(gridout[k, l, 1] == massid[i]) {
     bloco = bloco
      c
      
      
    }
    
  }
}

gridout[1,1,1]
gridout[1,1,2]
gridout[1,1,3]

gridout[1,2,1]
gridout[1,2,2]
gridout[1,2,3]

dim(gridout)

gridout[1,,]==1
slice.index(gridout,2)
