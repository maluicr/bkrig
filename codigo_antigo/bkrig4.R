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
stacf[is.na(stacf)] = -999

# create mask vector
mask = ifelse(stacf==-999, -1, 0)

# create array for blockfile
grdata = matrix (stacf, nrow = ny, ncol = nx, byrow = F)
grxy = expand.grid(x = seq(ox, ox + (nx-1) * resx, by = resx), y = seq(oy, oy + (ny-1) * resy, resy))
grx = grxy[,1]
gry = grxy[,2]
gridout = array (c (grdata, grx, gry), dim =c(ny, nx , 3))

# ------ write maskfile for dss ------

# prepare data to write file 
nvars = 1
namevars = "values"
nval = length(mask)

# create file
msk_name = "mask"
msk_nameO = paste0(msk_name, ".out")

if (file.exists(paste0(sid , "/", sid, msk_nameO))){
  file.remove(paste0(sid , "/", sid, msk_nameO))  
}

# start_time = Sys.time()
maskf = file(paste0(getwd(), "/", sid, "/", sid,  msk_nameO), open = "w")

# write file header
writeLines(msk_name, con = maskf, sep = "\n")
write(nvars, file = maskf, append = T)
write(namevars, file = maskf, append = T)

# write data
for (i in 1:nval){
  write(mask[i], file = maskf, append = T)
  }

msk_nzones = length(unique(mask))

close(maskf)
closeAllConnections()
# ----- observed covid19 data -----

phab = 10^4 

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
tab_all$error = round(tab_all$error, 3)

# NA cases set to 2 cases
tab_all$casos = ifelse(is.na(tab_all$casos), 2, tab_all$casos)

# recalculate crude rates
tab_all$rate = phab * tab_all$casos / tab_all$popres18

tabela = tab_all[, c("objectid", "x", "y", "z", "rate", "error")]

# ------ write datafile for dss ------

# prepare data to write file  
tabnotf = tabela[, c("x", "y", "z", "rate")]

nvars = ncol(tabnotf)
nobs = nrow(tabnotf)
namevars = names(tabnotf)

# create file
not_name = "notified"
not_nameO = paste0(not_name, ".out")

if (file.exists(paste0(sid , "/", sid, not_nameO))){
  file.remove(paste0(sid , "/", sid, not_nameO))  
}

# start_time = Sys.time()
notf = file(paste0(sid , "/", sid, not_nameO), open = "w")

# write filename
writeLines(not_name, con = notf, sep = "\n")

# write nr of variables
write(nvars, file = notf, append = T, sep = "\t")

# write name of variables
for (i in 1:nvars) {
  write(namevars[i], file = notf, ncolumns = nvars,  append = T, sep = "\t")
}

# write data
for (i in 1:nobs){
  vec = as.numeric(tabnotf[i,])
  write(vec, file = notf, ncolumns = nvars,  append = T, sep = "\t")
}

close(notf)
closeAllConnections()

# ------ write blocksfile for dss ------

massid = unique(grd)
massn = length(massid)


# create file
blk_name = "blockdata"
blk_nameO = paste0(blk_name,".out")

if (file.exists(paste0(sid , "/", sid, blk_nameO))){
  file.remove(paste0(sid , "/", sid, blk_nameO))  
}

# start_time = Sys.time()

logf = file(paste0(sid , "/", sid, blk_nameO), open = "w")

# write filename
writeLines(blk_name, con = logf, sep = "\n")

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

# ------ write ssdir.par for dss ------

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


ssd_name = "ssdir"
ssd_nameO = paste0(ssd_name, ".par")

if (file.exists(paste0(sid , "/", sid, ssd_nameO))){
  file.remove(paste0(sid , "/", sid, ssd_nameO))  
}

parsf = file(paste0(getwd(), "/", sid, "/", sid, ssd_nameO), open = "w")

writeLines("#*************************************************************************************#", con = parsf, sep = "\n")
writeLines("#                                                                                     #", con = parsf, sep = "\n")                                                                                     #');
writeLines("#             PARALLEL DIRECT SEQUENCIAL SIMULATION PARAMETER FILE                    #", con = parsf, sep = "\n")
writeLines("#                                                                                     #", con = parsf, sep = "\n")                                                                                     #');
writeLines("#*************************************************************************************#", con = parsf, sep = "\n")
writeLines("#                                                                                     #", con = parsf, sep = "\n")
writeLines("#-------------------------------------------------------------------------------------#", con = parsf, sep = "\n")
writeLines("#         Remember, # - Comment; [GROUP] - parameter group; CAPS - parameter          #", con = parsf, sep = "\n")
writeLines("#         Also, no space allowed in paths/filenames                                   #", con = parsf, sep = "\n")
writeLines("#-------------------------------------------------------------------------------------#", con = parsf, sep = "\n")
writeLines("#-------------------------------------------------------------------------------------#", con = parsf, sep = "\n")
writeLines("#                     here we define the hardata parameters                           #", con = parsf, sep = "\n")
writeLines("#-------------------------------------------------------------------------------------#", con = parsf, sep = "\n")
writeLines("", con = parsf, sep = "\n")
writeLines("[ZONES]", con = parsf, sep = "\n")
writeLines(paste0("ZONESFILE = ", msk_nameO), con = parsf, sep = "\n")
writeLines(paste0("NZONES = ", msk_nzones), con = parsf, sep = "\n")
writeLines("", con = parsf, sep = "\n")

for (i in 1:msk_nzones) {
  writeLines(paste0("[HARDDATA", i, "]" ), con = parsf, sep = "\n")
  writeLines(paste0("[DATAFILE = ", i , "]" ), con = parsf, sep = "\n")
  writeLines(paste0("[DATAFILE = ", i , "]" ), con = parsf, sep = "\n")
  
             fprintf(fid, '#s \n',['DATAFILE  =',varIn(n,:),'   # Hard Data file']);
fprintf(fid, '#s \n','COLUMNS   = 4                 # Number of columns in the data file (to eliminate)');
fprintf(fid, '#s \n','XCOLUMN   = 1                 # Column number of X coordinate');
fprintf(fid, '#s \n','YCOLUMN   = 2                 # Column number of Y coordinate');
fprintf(fid, '#s \n','ZCOLUMN   = 3                 # Column number of Z coordinate');
fprintf(fid, '#s \n','VARCOLUMN = 4                 # Column number for the variable');
fprintf(fid, '#s \n','WTCOLUMN  = 0                 # Column number for the declustering weight');
fprintf(fid, '#s \n',['MINVAL   = ', mat2str(bounds(n,1)),'          # Minimun threshold value']);
fprintf(fid, '#s \n',['MAXVAL   = ', mat2str(bounds(n,2)),'          # Maximum threshold value']);
fprintf(fid, '#s \n','USETRANS  = 1                 # Check if transformation file will be used');
fprintf(fid, '#s \n','TRANSFILE = Cluster.trn       # Transformation file');
}

close(parsf)

fprintf(fid, '\n');      
fprintf(fid, '#s \n','[HARDDATA]');
fprintf(fid, '#s \n',['ZMIN     = ',  mat2str(min(bounds(:,1))),'            # Minimum allowable data value ']);
fprintf(fid, '#s \n',['ZMAX     = ',  mat2str(max(bounds(:,2))),'            # Maximum allowable data value']);
fprintf(fid, '#s \n','LTAIL     = 1                   # Specify the back-transform implementation in the lower tail ');
fprintf(fid, '#s \n',['LTPAR    = ',  mat2str(min(bounds(:,1))),'            # Parameter for the ltail=2']);
fprintf(fid, '#s \n','UTAIL     = 1                   # Specify the back-transform implementation in the upper tail ');
fprintf(fid, '#s \n',['UTPAR    = ',  mat2str(max(bounds(:,2))),'            # Parameter for the utail=2']);
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define some parameters for the simulation');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[SIMULATION]');
fprintf(fid, '#s \n',['OUTFILE   = ', outputFilePath,'         # Filename of the resultant simulations']);
fprintf(fid, '#s \n',['NSIMS     = ',num2str(noSim),'   # Number of Simulations to generate']);
fprintf(fid, '#s \n','NTRY       = 10                   # Number of simulation for bias correction ');
fprintf(fid, '#s \n',['AVGCORR   = ', num2str(avgCorr), '       # Variance correction flag ']);
fprintf(fid, '#s \n',['VARCORR   = ', num2str(varCorr), '       # Average correction flag']);
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define the output grid (and secondary info grid)');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[GRID]');
fprintf(fid, '#s \n','# NX, NY and NZ are the number of blocks per direction ');
fprintf(fid, '#s \n',['NX        = ', num2str(XX(1,1))]);
fprintf(fid, '#s \n',['NY        = ',num2str(YY(1,1))]);
fprintf(fid, '#s \n',['NZ        = ',num2str(ZZ(1,1))]);
fprintf(fid, '#s \n','# ORIGX, ORIGY and ORIGZ are the start coordinate for each direction');
fprintf(fid, '#s \n',['ORIGX     = ', num2str(XX(1,2))]);
fprintf(fid, '#s \n',['ORIGY     = ', num2str(YY(1,2))]);
fprintf(fid, '#s \n',['ORIGZ     = ', num2str(ZZ(1,2))]);
fprintf(fid, '#s \n','# SIZEX, SIZEY and SIZEZ is the size of blocks in each direction ');
fprintf(fid, '#s \n',['SIZEX     = ', num2str(XX(1,3))]);
fprintf(fid, '#s \n',['SIZEY     = ', num2str(YY(1,3))]);
fprintf(fid, '#s \n',['SIZEZ     = ', num2str(ZZ(1,3))]);
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define some general parameters');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[GENERAL]');
fprintf(fid, '#s \n','NULLVAL       = -999.25                                     # Definition of the null value');
fprintf(fid, '#s \n',['SEED         = ', num2str(int32(rand*1000000)),'             # Seed for the pseudorandom number generator - if you want repeatable simulations']);
fprintf(fid, '#s \n','USEHEADERS    = 1                                             # Seed for the pseudorandom number generator - if you want repeatable simulations');
fprintf(fid, '#s \n',['FILETYPE     = ',simType, '                                  # accepts GEOEAS and SGEMS. default is sgems file type. used for output ']);
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','# here we define the parameters for search');
fprintf(fid, '#s \n','#-------------------------------------------------------------------------------------#');
fprintf(fid, '#s \n','[SEARCH]');
fprintf(fid, '#s \n','NDMIN   = 1                   # Minimum number of original data that should be used to simulate a node');
fprintf(fid, '#s \n','NDMAX   = 32                  # Maximum number of original data that should be used to simulate a node');
fprintf(fid, '#s \n','NODMAX  = 12                  # Maximum number of previouly simulated nodes to be used to simulate a new node');
fprintf(fid, '#s \n','SSTRAT  = 1                   # Two-part search / data nodes flag ');
fprintf(fid, '#s \n','MULTS   = 0                   # Multiple grid simulation flag ');
fprintf(fid, '#s \n','NMULTS  = 1                   # Number of multiple grid refinements ');
fprintf(fid, '#s \n','NOCT    = 0                   # Number of original data to use per octant');
fprintf(fid, '#s \n',['RADIUS1 = ', num2str(krigRANGE(1,1)),'                # Search radii in the maximum horizontal direction']);
fprintf(fid, '#s \n',['RADIUS2 = ', num2str(krigRANGE(1,2)),'                # Search radii in the minimum horizontal direction']);
fprintf(fid, '#s \n',['RADIUS3 = ', num2str(krigRANGE(1,3)),'                # Search radii in the vertical direction']);
fprintf(fid, '#s \n',['SANG1   = ', num2str(krigANG(1,1)),'                   # Orientation angle parameter of direction I (degrees)']); 
fprintf(fid, '#s \n',['SANG2   = ', num2str(krigANG(1,2)),'                   # Orientation angle parameter of direction II (degrees)']);
fprintf(fid, '#s \n',['SANG3   = ', num2str(krigANG(1,3)),'                   # Orientation angle parameter of direction III (degrees)']);
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
