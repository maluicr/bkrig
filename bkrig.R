library(gstat)
library(sp)
library(rgdal)
library(tidyverse)

# ---- folder to store results -----

day = 20200601

# create folder if exists is FALSE
sid = paste0(day)

# create folder for inputs
if(!file.exists(sid)) dir.create(sid, recursive=T)  

# store string with path and files prefix
wkdir = paste0(getwd(), "/", sid, "/", sid)

# create folder for outputs
sout = paste0("output", sid)
if(!file.exists(sout)) dir.create(sout, recursive=T)

# store string with path and files prefix   
wkdirout = paste0(getwd(), "/", sout, "/sim")

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

# ------ write maskfile for dss ------

# create mask vector
mask = ifelse(stacf==-999, -1, 0)
mask_zones = length(unique(mask))

# prepare data to write file 
nvars = 1
namevars = "values"
nval = length(mask)

# create file path
msk_name = "mask"
msk_nameO = paste0(msk_name, ".out")
fmsk = paste0(wkdir, msk_nameO)

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

# ----- observed covid19 data -----

# rate per phab habitants
phab = 10^4 

# dgs observed data
dgs = read.table("dgsData20200601.txt", header = T, sep = ",", dec = ".", strip.white = T )
names(dgs) = c("nome_cc", "casos")

# merge dgs with pop data
tab_all = merge(x = pop, y = dgs, by = c("nome_cc"), all.x = TRUE)

# crude rates
tab_all$rate = phab * tab_all$casos / tab_all$popres18

# compute overall mean rate
pop_tot = sum(tab_all$popres18, na.rm = T)
m = sum(tab_all$rate * tab_all$popres18, na.rm = T) / pop_tot

# error variance term (m/n_i)
tab_all$error = m / tab_all$popres18
tab_all$error = tab_all$error

# NA cases set to 2 cases
tab_all$casos = ifelse(is.na(tab_all$casos), 2, tab_all$casos)

# recalculate crude rates
tab_all$rate = phab * tab_all$casos / tab_all$popres18

tabela = tab_all[, c("objectid", "x", "y", "z", "rate", "error")]

# ------ write datafile for dss ------

# prepare data to write file  
tabnotf = tabela[, c("x", "y", "z", "rate")]

# store nr of variables
nvars = ncol(tabnotf)

# store nr of observations
nobs = nrow(tabnotf)

# store vector variable names
namevars = names(tabnotf)

# create file path
not_name = "notified"
not_nameO = paste0(not_name, ".out")
fnot = paste0(wkdir, not_nameO)

if (file.exists(fnot)){
  file.remove(fnot)  
}

# create notification file
file.create(fnot)

# write file header
cat(not_name, file = fnot, sep="\n")
cat(nvars, file = fnot, sep="\n", append = T)
cat(namevars, file = fnot, sep="\n", append = T)

# write notification data
write.table(format(tabnotf, digits = NULL, justify = "right"),
            file = fnot, quote = F, append = T, row.names = F, col.names = F)

# ------ experimental variograms -------

# get vectors

# incidence rates
rate = tab_all$rate
# coords x,y
ratexy = cbind(tab_all$x, tab_all$y)
# population
pop = tab_all$popres18 

# calculate distances

# lag distance
lagd = 10000
# lags up to
lagend = 100000
# nr lags
lag_n = floor(lagend / lagd) + 1
# vector of lags
lags = c()
for (i in 1:lag_n) {
  lags[i] = lagd * i - lagd  
}
# compute distance matrix
matd = as.matrix(dist(ratexy))

# calculate experimental variogram 

# create vectors to store results
nh = vector( mode = "integer", length = (lag_n-1)) # n pairs dist(h)
dh = vector( mode = "numeric", length = (lag_n-1)) # total dist(h)
mh = vector( mode = "numeric", length = (lag_n-1)) # mean dist(h)
num_gh = vector( mode = "numeric", length = (lag_n-1)) # numerator gamma(h)
den_gh = vector( mode = "numeric", length = (lag_n-1)) # denominator gamma(h)
gammah = vector( mode = "numeric", length = (lag_n-1)) # gamma(h)

# compute variogram
for (k in 1:(lag_n-1)) {
  for (i in 1:nobs) {
    for (j in 1:nobs) {
      if(matd[i,j] > lags[k] & matd[i,j] <= lags[k+1]){
        nh[k] = nh[k] + 1
        dh[k] = dh[k] + matd[i,j]
        num_gh[k] = num_gh[k] + ((pop[i] * pop[j]) * (rate[i]-rate[j])^2 - m) / (pop[i] + pop[j])
        den_gh[k] = den_gh[k] + pop[i] * pop[j] / (pop[i] + pop[j])
       }
    }
  }
  gammah[k] = num_gh[k]/(2*den_gh[k])
  mh[k] = dh[k] / nh[k]
}

plot(mh,gammah)
  

  

variance=0
poptotal=0

varw=pop(i)*pop(i)/(pop(i)+pop(i))
varp=varw*((racio(i)-media)**2)
poptotal=poptotal+varw
variance=variance+varp
do j=1, nobs
dist(i,j)=sqrt((coordx(i)-coordx(j))**2+(coordy(i)-coordy(j))**2)
w(i,j)=pop(i)*pop(j)/(pop(i)+pop(j))
gama1(i,j)=w(i,j)*((racio(i)-racio(j))**2)
write(30,*) i, j, dist(i,j), w(i,j), gama1(i,j)
end do
end do

variance=variance/poptotal

lag=5000
nv=20
do k=1,nv
nh(k)=0
soma1(k)=0
soma2(k)=0
distot(k)=0
end do

do i=1,nobs
do j=1,nobs
do k=1,nv
if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) nh(k)=nh(k)+1
if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) soma1(k)=soma1(k)+gama1(i,j)
if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) distot(k)=distot(k)+dist(i,j)
if ((dist(i,j)<=(lag*k)) .and. (dist(i,j)>lag*(k-1))) soma2(k)=soma2(k)+w(i,j)
end do          
end do
end do

write(20,*) variance
do k=1,nv
distot(k)=distot(k)/nh(k)
gama2(k)=(1/(2*soma2(k)))*soma1(k)
write(20,*) nh(k)/2, distot(k), gama2(k)
end do


# ------ create vector kriging pars -------

exampleFunction<-function(A=1,B=c(3,4)){
  return(unionvector=c(A,B))
} 
exampleFunction()

function (models = c("Exp", "Sph"), pepita = 0, range = 60000, psill= c()) {
  for (i in 1:models) {
    
    
  }
}

# only exponential and spherical model are implemented for now!
# functions are incomplete -> not defined for values h=0 and h > range
# to add.

spatcov = function(mod = c("Sph", "Exp"), nug = 0, C = 1, r = 1, mat){
  if(!is.matrix(mat)) stop("mat must be a matrix")
  dmat0=mat[-1,-1] # c[z(x)-z(x+h)]
  dmat1=mat[-1,1]  # c[z(x)-z(0)]
  n_0=nrow(dmat0)
  if (mod=="Sph"){
    c=matrix(nug+C*(1-1.5*(dmat0/r)-0.5*(dmat0/r)^3),nrow=n_0)
    c0=matrix(nug+C*(1-1.5*(dmat1/r)-0.5*(dmat1/r)^3),nrow=n_0)
  } else if (mod=="Exp"){
    c=matrix(nug+C*(exp(-3*dmat0/r)),nrow=n_0)
    c0=matrix(nug+C*(exp(-3*dmat1/r)),nrow=n_0)
  } else {
    stop("Covariance model must be 'Sph' or 'Exp'")
  }
  matK=rbind(c,1)
  matK=cbind(matK,1)
  ndim=dim(matK) 
  matK[ndim,ndim]=0
  matK=round(matK,3)
  matM=c(c0,1)
  matM=matrix(matM,ncol=1)
  matM=(round(matM,3))
  x=list(covK=matK,covM=matM)
  return(x)
}

spcov=spatcov(mod="Exp",C=10,r=10,mat=dmat)

# ------ write blocksfile for dss ------

# store vector grid id
massid = unique(grd)

# store length vector
massn = length(massid)

# create file path
blk_name = "blockdata"
blk_nameO = paste0(blk_name,".out")
fblk = paste0(wkdir, blk_nameO)

if (file.exists(fblk)){
  file.remove(fblk)  
}

# create notification file
file.create(fblk)

# write file header
cat(blk_name, file = fblk, sep="\n")
cat(massn, file = fblk, sep="\n", append = T)
cat(namevars, file = fnot, sep="\n", append = T)

t0 = Sys.time()

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
  cat(paste0("Block#", id), file = fblk, sep="\n", append = T)
  
  # write rate, error & nr blocks
  cat(rt, file = fblk, sep="\n", append = T)
  cat(er, file = fblk, sep="\n", append = T)
  cat(nc, file = fblk, sep="\n", append = T)
  
  # data.table::fwrite
  blk_list = list()
  library("data.table")
  for (k in 1:ny) {
    for(l in 1:nx){
      if(gridout[k, l, 1] == id) {
        d = paste0(gridout[k, l, 2], "\t", gridout[k, l, 3], "\t", 0)
        blk_list[[length(blk_list)+1]] = list (d)
      }
    }
  }
  data.table::fwrite(blk_list, file = fblk, append = T, sep="\n")
  # write x, y, z=0 coords block
  #for (k in 1:ny) {
  #  for(l in 1:nx){
  #    if(gridout[k, l, 1] == id) {
  #      d = paste0(gridout[k, l, 2], "\t", gridout[k, l, 3], "\t", 0)
  #       cat(d , file = fblk, append = T)
  #      cat("\n" , file = fblk, append = T)
  #    }
  #  }
  # }
}

t1 = Sys.time()
t1-t0

# ------ write ssdir.par for dss ------

# store number of simulations
nsim = 100

# store nr simulations for bias correction
nbias = 20

# store flag for (mean, var) correction (yes = 1, no = 0)
biascor = c(1,1)

# draw pseudo-random value 
pseudon = sample(10^8,1)

## Run external DSS exectuable and return realization
#
# [input]
#   simType (string) - SGEMS or GEOEAS type
#   avgCorr (0/1) - correct mean
#   varCorr (0/1) - correct variance 
#   usebihist (0/1) - use joint probability distributions
#   inputPath (string) - path for folder with input data and DSS executable
#   varIn (string) - name of the file with input experimental data
#   noSim - number of realizations, 
#   outputFilePath (string) - full path for output file name without extension
#   krigType (0/1/2/3/4/5) - kriging type
#   secVar (string) - fullpath for secondary variable
#   bihistFile (string) - fullpath for joint distribution file
#   bounds (noZones x 2) - min and max of variable to be simulated
#   XX (1 x 3) - Number of cells, origin, size in XX
#   YY (1 x 3) - Number of cells, origin, size in YY
#	  ZZ (1 x 3) - Number of cells, origin, size in ZZ
#   localCorr (string) - fullpath for collocated correlation coefficient file
#   globalCorr - global correlation coefficient between 2 variables
#   auxVar string) - fullpath for sec variable for joint distribution
#   krigANG (1 x 3) - azimuth, rake and dip
#   krigRANGE (1 x 3) - major, minor, vertical
#   varANG (1 x 3) - azimuth, rake and dip
#   varRANGE (noZones x 3) - major, minor, vertical
#   varType(noZones x 1) - variogram type
#   varNugget (noZones x 1) - nugget effect
#   zoneFileName(string) - fullpath for zone file 
#   noZones - number of zones in zonefile
#   noClasses - number of classes for joint distribution
#   rescale (0/1) - rescale secondary variable for co-simulation
#   lvmFile (string) - fullpath for local varying means
#   usePseudoHard (0/1) - for point distributions
#   pseudoHardFile (string) - fullpath for local pdfs
#   pseudoCorr 0/1) - correct local pdfs
#   covTab (1 x 3) - size of the covariance matrix to be stored in mem
# 
# [output]
#   var_out - realization
#
#
# Leonardo Azevedo - 2019, Manuel Ribeiro - 2020
# CERENA/Instituto Superior T?cnico (Portugal)
#

# create file path
ssd_name = "ssdir"
ssd_nameO = paste0(ssd_name, ".par")
fssd = paste0(wkdir, ssd_nameO)

if (file.exists(fssd)){
  file.remove(fssd)  
}

# create notification file
file.create(fssd)

# write file header
cat(not_name, file = fnot, sep="\n")
cat(nvars, file = fnot, sep="\n", append = T)
cat(namevars, file = fnot, sep="\n", append = T)

cat("#*************************************************************************************#", file = fssd, sep = "\n", append = T)
cat("#                                                                                     #", file = fssd, sep = "\n", append = T)
cat("#             PARALLEL DIRECT SEQUENCIAL SIMULATION PARAMETER FILE                    #", file = fssd, sep = "\n", append = T)
cat("#                                                                                     #", file = fssd, sep = "\n", append = T)                                                                                     #');
cat("#*************************************************************************************#", file = fssd, sep = "\n", append = T)
cat("#                                                                                     #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#         Remember, # - Comment; [GROUP] - parameter group; CAPS - parameter          #", file = fssd, sep = "\n", append = T)
cat("#         Also, no space allowed in paths/filenames                                   #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#                     here we define the hardata parameters                           #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)

cat("\n", file = fssd, append = T)
cat("[ZONES]", file = fssd, sep = "\n", append = T)
cat(paste0("ZONESFILE = ", msk_nameO), file = fssd, sep = "\n", append = T)
cat(paste0("NZONES = ", mask_zones), file = fssd, sep = "\n", append = T)
cat("", file = fssd, sep = "\n", append = T)

for (i in 1:mask_zones) {
  cat(paste0("[HARDDATA", i, "]" ), file = fssd, sep = "\n", append = T)
  cat(paste0("DATAFILE = ", fnot ), file = fssd, sep = "\n", append = T)
  cat(paste0("COLUMNS = ", nvars), file = fssd, sep = "\n", append = T)
  cat(paste0("XCOLUMN = ", grep("x", colnames(tabnotf))), file = fssd, sep = "\n", append = T)
  cat(paste0("YCOLUMN = ", grep("y", colnames(tabnotf))), file = fssd, sep = "\n", append = T)
  cat(paste0("ZCOLUMN = ", grep("z", colnames(tabnotf))), file = fssd, sep = "\n", append = T)
  cat(paste0("VARCOLUMN = ", grep("rate", colnames(tabnotf))), file = fssd, sep = "\n", append = T)
  cat("WTCOLUMN = 0", file = fssd, sep = "\n", append = T)
  cat(paste0("MINVAL = ", min(tabnotf$rate)), file = fssd, sep = "\n", append = T)
  cat(paste0("MAXVAL = ", max(tabnotf$rate)), file = fssd, sep = "\n", append = T)
  cat(paste0("USETRANS = 1"), file = fssd, sep = "\n", append = T)
  cat(paste0("TRANSFILE =  Cluster.trn"), file = fssd, sep = "\n", append = T)
}

cat("\n", file = fssd, append = T )

cat("[HARDDATA]", file = fssd, sep = "\n", append = T)
cat(paste0("ZMIN = ", min(tabnotf$rate)), file = fssd, sep = "\n", append = T)
cat(paste0("ZMAX = ", max(tabnotf$rate)), file = fssd, sep = "\n", append = T)
cat("LTAIL = 1", file = fssd, sep = "\n", append = T)
cat(paste0("LTPAR = ", min(tabnotf$rate)), file = fssd, sep = "\n", append = T)
cat("UTAIL = 1", file = fssd, sep = "\n", append = T)
cat(paste0("UTPAR = ", max(tabnotf$rate)), file = fssd, sep = "\n", append = T)

cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#                here we define parameters for the simulation                         #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)

cat("\n", file = fssd, append = T)
cat("[SIMULATION]", file = fssd, sep = "\n", append = T)

cat(paste0("OUTFILE = ", wkdirout), file = fssd, sep = "\n", append = T)
cat(paste0("NSIMS = ", nsim), file = fssd, sep = "\n", append = T)
cat(paste0("NTRY = ", nbias), file = fssd, sep = "\n", append = T)
cat(paste0("AVGCORR = ", biascor[1]), file = fssd, sep = "\n", append = T)
cat(paste0("VARCORR = ", biascor[2]), file = fssd, sep = "\n", append = T)

cat("\n", file = fssd, sep = "\n", append = T)

cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#        here we define the output grid (and secondary info grid)                     #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)

cat("\n", file = fssd, append = T)

cat("[GRID]", file = fssd, sep = "\n", append = T)
cat("# NX, NY and NZ are the number of blocks per direction", file = fssd, sep = "\n", append = T)
cat(paste0("NX = ", nx), file = fssd, sep = "\n", append = T)
cat(paste0("NY = ", ny), file = fssd, sep = "\n", append = T)
cat(paste0("NZ = ", 1), file = fssd, sep = "\n", append = T)

cat("# ORIGX, ORIGY and ORIGZ are the start coordinate for each direction", file = fssd, sep = "\n", append = T)

cat(paste0("ORIGX = ", ox), file = fssd, sep = "\n", append = T)
cat(paste0("ORIGY = ", oy), file = fssd, sep = "\n", append = T)
cat(paste0("ORIGZ = ", 0), file = fssd, sep = "\n", append = T)

cat("# SIZEX, SIZEY and SIZEZ is the size of blocks in each direction ", file = fssd, sep = "\n", append = T)

cat(paste0("SIZEX = ", resx ), file = fssd, sep = "\n", append = T)
cat(paste0("SIZEY = ", resy ), file = fssd, sep = "\n", append = T)
cat(paste0("SIZEZ = ", 1), file = fssd, sep = "\n", append = T)

cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#                 here we define some general parameters                              #", file = fssd, sep = "\n", append = T)
cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)

cat("[GENERAL]", file = fssd, sep = "\n", append = T)

cat(paste0("NULLVAL = ", nas), file = fssd, sep = "\n", append = T)
cat(paste0("SEED = ", pseudon), file = fssd, sep = "\n", append = T)
cat("USEHEADERS = 1", file = fssd, sep = "\n", append = T)
cat("FILETYPE = GEOEAS", file = fssd, sep = "\n", append = T)

cat("\n\n", file = fssd, append = T)

cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
cat("#                 here we define the parameters for search                            #", file = fssd, sep = "\n", append = T)
cat("#------------ -------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)

cat("[SEARCH]", file = fssd, sep = "\n", append = T)

# ------- data used to simulate a node -------

# min nr of observed samples
cat("NDMIN   = 1", file = fssd, sep = "\n", append = T)

# max nr of observed samples
cat("NDMAX   = 32", file = fssd, sep = "\n", append = T)

# max nr of previouly simulated nodes
cat("NODMAX  = 12", file = fssd, sep = "\n", append = T)
               
# Two-part search / data nodes flag
cat("SSTRAT  = 1", file = fssd, sep = "\n", append = T)

# Multiple grid simulation flag
cat("MULTS   = 0", file = fssd, sep = "\n", append = T)

# Nr of multiple grid refinements
cat("NMULTS   = 1", file = fssd, sep = "\n", append = T)  

# Nr of original data per octant
cat("NOCT    = 0", file = fssd, sep = "\n", append = T)

# Search radii in the major horizontal axe
cat(paste0("RADIUS1 =", poraquiqqcoisa), file = fssd, sep = "\n", append = T)
fprintf(fid, '#s \n',['RADIUS1 = ', num2str(krigRANGE(1,1)),'

# Search radii in the minimum horizontal direction']);
fprintf(fid, '#s \n',['RADIUS2 = ', num2str(krigRANGE(1,2)),'

# Search radii in the vertical direction']);
fprintf(fid, '#s \n',['RADIUS3 = ', num2str(krigRANGE(1,3)),'

# Orientation angle parameter of direction I (degrees)']);
fprintf(fid, '#s \n',['SANG1   = ', num2str(krigANG(1,1)),' 

# Orientation angle parameter of direction II (degrees)']);
fprintf(fid, '#s \n',['SANG2   = ', num2str(krigANG(1,2)),'  

# Orientation angle parameter of direction III (degrees)']);
fprintf(fid, '#s \n',['SANG3   = ', num2str(krigANG(1,3)),'  


fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define the kriging information, and secondary info when applicable');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[KRIGING]');
fprintf(fid, '#s \n',['KTYPE        = ', num2str(krigType),'                # Kriging type ;0=simple,1=ordinary,2=simple with locally varying mean,']);
fprintf(fid, '#s \n','# 3=external drif, 4=collo-cokrig global CC,5=local CC (KTYPE)');
fprintf(fid, '#s \n',['COLOCORR     = ', num2str(globalCorr),'                # Global CC to ktype=4 ']);
fprintf(fid, '#s \n',['SOFTFILE     = ', secVar,'                                  # Filename of the soft data']);
fprintf(fid, '#s \n',['LVMFILE      = ', lvmFile,'                                 # FOR KTYPE=2']);
fprintf(fid, '#s \n','NVARIL        = 1                                             # Number of columns in the secundary data file');
fprintf(fid, '#s \n','ICOLLVM       = 1                                             # Column number of secundary variable ');
fprintf(fid, '#s \n',['CCFILE       = ', localCorr,'                               # Filename of correlation file for local correlations (ktype=5)']);
fprintf(fid, '#s \n',['RESCALE      = ', num2str(rescale),'                         # rescale secondary variable (for ktype>=4)']);
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define the variogram to use. if more than 1, use [VARIOGRAM2]');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
for n=1:noZones
fprintf(fid, '#s \n',['[VARIOGRAMZ',num2str(n),']']);
fprintf(fid, '#s \n',['NSTRUCT  = ',num2str(size(varRANGE,2)/4),'                # Number of semivariograms structures (NST(1))']);
fprintf(fid, '#s \n',['NUGGET    =', num2str(varNugget(n,1)),'                                               # Nugget constant (C0(1))']);
for i=1:size(varRANGE,2)/4
fprintf(fid, '#s \n',['[VARIOGRAMZ',num2str(n),'S',num2str(i),']']);
fprintf(fid, '#s \n',['TYPE =', num2str(varType(n,1)), '                   # Struture type ;1=spherical,2=exponential,3=gaussian (IT(i))']);
fprintf(fid, '#s \n',['COV  =', num2str(varRANGE(n,4*i))  ,'               # C parameter "COV + NUGGET = 1.0" (CC(i))']);
fprintf(fid, '#s \n',['ANG1 =', num2str(varANG(n,1))  ,'                   # Geometric anisotropy angle I (ANG1(i))']);
fprintf(fid, '#s \n',['ANG2 =', num2str(varANG(n,2))  ,'                   # Geometric anisotropy angle II (ANG2(i))']);
fprintf(fid, '#s \n',['ANG3 =', num2str(varANG(n,3))  ,'                   # Geometric anisotropy angle III (ANG3(i))']);
fprintf(fid, '#s \n',['AA   =', num2str(varRANGE(n,4*i-3)),'               # Maximum horizontal range (AA(i))']);
fprintf(fid, '#s \n',['AA1  =', num2str(varRANGE(n,4*i-2)),'               # Minimum horizontal range (AA1) ']);
fprintf(fid, '#s \n',['AA2  =', num2str(varRANGE(n,4*i-1)),'               # Vertical range (AA2)']);
end
end
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define parameters for joint DSS ');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
for n=1:noZones
fprintf(fid, '#s \n',['[BIHIST',num2str(n),']']);
fprintf(fid, '#s \n',['USEBIHIST        = ', num2str(usebihist(1,1)),' #Use Bihist? 1-yes 0-no']);
fprintf(fid, '#s \n',['BIHISTFILE       = ', bihistFile(n,:),'   # bihistogram file']);
fprintf(fid, '#s \n',['NCLASSES         = ', num2str(noClasses),'                           # number of classes to use']);
fprintf(fid, '#s \n',['AUXILIARYFILE    = ', auxVar,'                   # auxiliary image']);            
end
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define the debug parameters - probably you wont need this');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[DEBUG]  ');
fprintf(fid, '#s \n','DBGLEVEL  = 1                 # 1 to 3, use higher than 1 only if REALLY needed');
fprintf(fid, '#s \n','DBGFILE   = debug.dbg         # File to write debug');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define parameters for COVARIANCE TABLE - reduce if memory is a problem ');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[COVTAB]');
fprintf(fid, '#s \n',['MAXCTX = ', num2str(covTab(1,1))]);
fprintf(fid, '#s \n',['MAXCTY = ', num2str(covTab(1,2))]);
fprintf(fid, '#s \n',['MAXCTZ = ', num2str(covTab(1,3))]);
fprintf(fid,'#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid,'#s \n','# here we define parameters for BLOCK KRIGING - if you are not block kriging useblocks should be 0 ');
fprintf(fid,'#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid,'#s \n','[BLOCKS]');
fprintf(fid,'#s \n',['USEBLOCKS  = ', num2str(useBlocks), '   # 1 use, 0 no']);
fprintf(fid,'#s \n',['BLOCKSFILE = ', blocksFile,'           # file']);
fprintf(fid,'#s \n',['MAXBLOCKS  = ', num2str(maxBlocks)]);
fprintf(fid, '#s \n','[PSEUDOHARD]');
fprintf(fid,'#s \n',['USEPSEUDO  = ', num2str(usePseudoHard), '  # 1 use, 0 no  pseudo hard data is point distributions that are simulated before all other nodes']);
fprintf(fid,'#s \n',['PSEUDOFILE = ', pseudoHardFile,'           # file']);
fprintf(fid,'#s \n',['PSEUDOCORR = ', num2str(pseudoCorr),'               # correct simulated value with point']);
fclose(fid);
##   RUN DSS
temp = [pwd,'/',inputPath,'/DSS.C.64.exe   ',inputPath,'/ssdir.par'];

# hack to silent mode in command window
[T,Status,results] = evalc('system(temp)');

##  IMPORT SIMULATIONS
if (strcmp(simType ,'SGEMS') == 1);

varStrut    = sgems_read([outputFilePath,'_1.sgems']);
var_out     = varStrut.data;
clear varStrut

else
  [~,var_out]=import_gslib([outputFilePath,'_1.out']);
end

end
