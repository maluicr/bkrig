# ----- description -----

# function writes the parameters file for dss.exe
# and runs block kriging simulations

# as input you should provide objects returned by functions
# irates(), blockfile(), maskfile() and varmodel() and other 
# parameter values for  kriging and simulation processes 

# ----- arguments ------

# blockobj, string, name of list, output of function blockfile()
# maskobj, string, name of list, output of function maskfile()
# dfobj, string, name of list, output of function irates()
# varmobj, string, name of list, output of function varmodel()
# simulations, numeric, number of simulations
# nrbias, numeric, nr simulations for bias correction
# biascor, num vector, flag for (mean, variance) correction (yes = 1, no = 0)
# ndMin, numeric, min number of neighbour observations used in  kriging
# ndMax, numeric, max number of neighbour observations used in  kriging
# nodMax, numeric, max number of previously simulated nodes used in  kriging
# radius1, numeric, search radii in the major horizontal axe
# radius2, numeric, search radii in the axe orthogonal (horizontal) to radius1
# radius3, numeric, search radii in the vertical axe
# ktype, numeric, the kriging type to be used (available are: 0 = simple, 1 = ordinary)

# ------ function ------

ssdpars = function (blockobj, maskobj, dfobj, varmobj, simulations = 1, nrbias = 20, biascor = c(1,1),
                    ndMin = 1, ndMax = 32, nodMax = 12, radius1, radius2, radius3 = 1, ktype = 0) {

  day = as.character(blockobj[["file"]]["day"])
  folder = as.character(blockobj[["file"]]["folder"])
  
  # filenames
  rf = as.character(dfobj[["file"]]["name"])
  bf = as.character(blockobj[["file"]]["name"])
  mf = as.character(maskobj[["file"]]["name"])
  
  # hard data pars
  nvars = as.numeric(dfobj[["ssdirpars"]]["nvars"])
  xcolumn = as.numeric(dfobj[["ssdirpars"]]["xcolumn"])
  ycolumn = as.numeric(dfobj[["ssdirpars"]]["ycolumn"])
  zcolumn = as.numeric(dfobj[["ssdirpars"]]["zcolumn"])
  varcol = as.numeric(dfobj[["ssdirpars"]]["varcol"])
  minval = as.numeric(dfobj[["ssdirpars"]]["minval"])
  maxval = as.numeric(dfobj[["ssdirpars"]]["maxval"])
  
  # mask grid pars
  nzones = as.numeric(maskobj[["zones"]]["nzones"])
  
  # output grid pars
  nx = as.numeric(unlist(blockobj[["gridpars"]]["nodes"]))[1]
  ny = as.numeric(unlist(blockobj[["gridpars"]]["nodes"]))[2]
  ox = as.numeric(unlist(blockobj[["gridpars"]]["origin"]))[1]
  oy = as.numeric(unlist(blockobj[["gridpars"]]["origin"]))[2]
  rx = as.numeric(unlist(blockobj[["gridpars"]]["resolution"]))[1]
  ry = as.numeric(unlist(blockobj[["gridpars"]]["resolution"]))[2]
  
  # general pars
  nas = as.numeric(blockobj[["gridpars"]]["NAs"])
  
  # varigram pars
  nstruct = as.numeric(varmobj[["structures"]])
  nugget = as.numeric(varmobj[["parameters"]]["nugget"])
  range = as.numeric(varmobj[["parameters"]]["range"])
  psill = as.numeric(varmobj[["parameters"]]["psill"])
  nuggetp = nugget/(nugget + psill)
  psillp = psill/(nugget + psill)
  mtype = as.numeric(varmobj[["parameters"]]["modeltype"])
  
  # block kriging pars
  maxblocks = as.numeric(blockobj[["outgrid"]]["nblock"])
  
  
  # ------ write ssdir.par for dss ------
  
  # store number of simulations
  nsim = simulations
  
  # store nr simulations for bias correction
  nbias = nrbias
  
  # store flag for (mean, var) correction (yes = 1, no = 0)
  biascor = biascor
  
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
  # Leonardo Azevedo - 2019
  # CERENA/Instituto Superior Tecnico (Portugal)
  #
  
  # create file path
  ssd_name = "ssdir"
  ssd_nameO = paste0(ssd_name, ".par")
  fssd = paste0(folder, "/", day, ssd_nameO)
  
  if (file.exists(fssd)){
    file.remove(fssd)  
  }
  
  # create notification file
  file.create(fssd)
  
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
  cat(paste0("ZONESFILE = ", mf), file = fssd, sep = "\n", append = T)
  cat(paste0("NZONES = ", nzones), file = fssd, sep = "\n", append = T)
  cat("", file = fssd, sep = "\n", append = T)
  
  for (i in 1:nzones) {
    cat(paste0("[HARDDATA", i, "]" ), file = fssd, sep = "\n", append = T)
    cat(paste0("DATAFILE = ", rf), file = fssd, sep = "\n", append = T)
    cat(paste0("COLUMNS = ", nvars), file = fssd, sep = "\n", append = T)
    cat(paste0("XCOLUMN = ", xcolumn), file = fssd, sep = "\n", append = T)
    cat(paste0("YCOLUMN = ", ycolumn), file = fssd, sep = "\n", append = T)
    cat(paste0("ZCOLUMN = ", zcolumn), file = fssd, sep = "\n", append = T)
    cat(paste0("VARCOLUMN = ", varcol), file = fssd, sep = "\n", append = T)
    cat("WTCOLUMN = 0", file = fssd, sep = "\n", append = T)
    cat(paste0("MINVAL = ", minval), file = fssd, sep = "\n", append = T)
    cat(paste0("MAXVAL = ", maxval), file = fssd, sep = "\n", append = T)
    cat(paste0("USETRANS = 1"), file = fssd, sep = "\n", append = T)
    cat(paste0("TRANSFILE =  Cluster.trn"), file = fssd, sep = "\n", append = T)
  }
  
  cat("\n", file = fssd, append = T )
  
  cat("[HARDDATA]", file = fssd, sep = "\n", append = T)
  cat(paste0("ZMIN = ", minval), file = fssd, sep = "\n", append = T)
  cat(paste0("ZMAX = ", maxval), file = fssd, sep = "\n", append = T)
  cat("LTAIL = 1", file = fssd, sep = "\n", append = T)
  cat(paste0("LTPAR = ", minval), file = fssd, sep = "\n", append = T)
  cat("UTAIL = 1", file = fssd, sep = "\n", append = T)
  cat(paste0("UTPAR = ", maxval), file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#                here we define parameters for the simulation                         #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[SIMULATION]", file = fssd, sep = "\n", append = T)
  cat(paste0("OUTFILE = ", day, "sim"), file = fssd, sep = "\n", append = T)
  cat(paste0("NSIMS = ", nsim), file = fssd, sep = "\n", append = T)
  cat(paste0("NTRY = ", nbias), file = fssd, sep = "\n", append = T)
  cat(paste0("AVGCORR = ", biascor[1]), file = fssd, sep = "\n", append = T)
  cat(paste0("VARCORR = ", biascor[2]), file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define the output grid (and secondary info grid)                     #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
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
  
  cat(paste0("SIZEX = ", rx ), file = fssd, sep = "\n", append = T)
  cat(paste0("SIZEY = ", ry ), file = fssd, sep = "\n", append = T)
  cat(paste0("SIZEZ = ", 1), file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#                 here we define some general parameters                              #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[GENERAL]", file = fssd, sep = "\n", append = T)
  cat(paste0("NULLVAL = ", nas), file = fssd, sep = "\n", append = T)
  cat(paste0("SEED = ", pseudon), file = fssd, sep = "\n", append = T)
  cat("USEHEADERS = 1", file = fssd, sep = "\n", append = T)
  cat("FILETYPE = GEOEAS", file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#                 here we define the parameters for search                            #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[SEARCH]", file = fssd, sep = "\n", append = T)
  
  # ------- ssdir.par: search parameters -------
  
  # hard-coded values
  sstrat = 1 
  mults = 0
  nmults = 1
  noct = 0
  sang1 = 0
  sang2 = 0
  sang3 = 0
  
  # min nr of observed samples
  cat(paste0("NDMIN  = ", ndMin), file = fssd, sep = "\n", append = T)
  # max nr of observed samples
  cat(paste0("NDMAX  = ", ndMax), file = fssd, sep = "\n", append = T)
  # max nr of previouly simulated nodes
  cat(paste0("NODMAX = ", nodMax), file = fssd, sep = "\n", append = T)
  # Two-part search / data nodes flag
  cat(paste0("SSTRAT = ", sstrat), file = fssd, sep = "\n", append = T)
  # Multiple grid simulation flag
  cat(paste0("MULTS  = ", mults), file = fssd, sep = "\n", append = T)
  # Nr of multiple grid refinements
  cat(paste0("NMULTS = ", nmults), file = fssd, sep = "\n", append = T)  
  # Nr of original data per octant
  cat(paste0("NOCT = ", noct), file = fssd, sep = "\n", append = T)
  # Search radii in the major horizontal axe
  cat(paste0("RADIUS1 = ", radius1), file = fssd, sep = "\n", append = T)
  # Search radii in the ortogonal horizontal axe (to major)
  cat(paste0("RADIUS2 = ", radius2), file = fssd, sep = "\n", append = T)
  # Search radii in the vertical axe
  cat(paste0("RADIUS3 = ", radius3), file = fssd, sep = "\n", append = T)
  # Orientation angle parameter of direction I (degrees)
  cat(paste0("SANG1 = ", sang1), file = fssd, sep = "\n", append = T)
  # Orientation angle parameter of direction II (degrees)
  cat(paste0("SANG2 = ", sang2), file = fssd, sep = "\n", append = T)
  # Orientation angle parameter of direction III (degrees)
  cat(paste0("SANG3 = ", sang3), file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#   here we define the kriging information, and secondary info when applicable        #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  # ------- ssdir.par: krige info -------
  
  # hard-coded values
  colorcorr = 0
  softfile = "no file"
  lvmfile = "no file"
  nvaril = 1
  icollvm = 1
  ccfile = "no file"
  rescale = 0
  
  cat("[KRIGING]", file = fssd, sep = "\n", append = T)
  # Kriging type: 0 = simple, 1 = ordinary, 2 = simple with locally varying mean
  # 3 = external drift, 4 = collo-cokrig global CC, 5 = local CC
  cat(paste0("KTYPE = ", ktype), file = fssd, sep = "\n", append = T)
  # Global coef correlation (ktype = 4)
  cat(paste0("COLOCORR = ", colorcorr), file = fssd, sep = "\n", append = T)
  # Filename of the soft data (ktype = 2)
  cat(paste0("SOFTFILE = ", softfile), file = fssd, sep = "\n", append = T)
  # For ktype = 2
  cat(paste0("LVMFILE  = ", lvmfile), file = fssd, sep = "\n", append = T)
  # Number of columns in the secundary data file
  cat(paste0("NVARIL  = ", nvaril), file = fssd, sep = "\n", append = T)
  # Column number of secundary variable 
  cat(paste0("ICOLLVM  = ", icollvm), file = fssd, sep = "\n", append = T)
  # Filename of correlation file for local correlations (ktype = 5)
  cat(paste0("CCFILE  = ", ccfile), file = fssd, sep = "\n", append = T)
  # Rescale secondary variable (for ktype >= 4)
  cat(paste0("RESCALE  = ", rescale), file = fssd, sep = "\n", append = T)
  
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define the variogram to use. if more than 1, use [VARIOGRAM2]        #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  # ------- ssdir.par: variogram models -------
  
  for (j in 1 : nzones) {
    cat(paste0("[VARIOGRAMZ", j, "]"), file = fssd, sep = "\n", append = T)
    cat(paste0("NSTRUCT = ", nstruct), file = fssd, sep = "\n", append = T)
    cat(paste0("NUGGET = ", nuggetp), file = fssd, sep = "\n", append = T)
    for (i in 1 : nstruct){
      cat(paste0("[VARIOGRAMZ", j, "S", i, "]"), file = fssd, sep = "\n", append = T)
      # store struture type ; 1 = spherical, 2 = exponential
      cat(paste0("TYPE = ", mtype), file = fssd, sep = "\n", append = T)
      # C parameter "COV + NUGGET = 1.0" (CC(i))
      cat("COV = 1", file = fssd, sep = "\n", append = T)
      # Geometric anisotropy angle I (ANG1(i))
      cat("ANG1 = 0", file = fssd, sep = "\n", append = T)
      # Geometric anisotropy angle II (ANG2(i))
      cat("ANG2 = 0", file = fssd, sep = "\n", append = T)
      # Geometric anisotropy angle III (ANG3(i))
      cat("ANG3 = 0", file = fssd, sep = "\n", append = T)
      # Maximum horizontal range (AA(i))
      cat(paste0("AA = ", range), file = fssd, sep = "\n", append = T)
      # Minimum horizontal range (AA1)
      cat(paste0("AA1 = ", range), file = fssd, sep = "\n", append = T)
      # Vertical range (AA2)
      cat("AA2 = 1", file = fssd, sep = "\n", append = T)
    }
  }
  
  cat("#-------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define parameters for joint DSS        #", file = fssd, sep = "\n", append = T)
  cat("#-------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  # hard-coded values 
  usebihist = 0
  bihistfile = "no file"
  nclasses = 0
  auxfile = "no file"
  
  for (j in 1:nzones){
    cat(paste0("[BIHIST", j, "]"), file = fssd, sep = "\n", append = T)
    #Use Bihist? 1-yes 0-no'
    cat(paste0("USEBIHIST = ", usebihist), file = fssd, sep = "\n", append = T)
    # bihistogram file
    cat(paste0("BIHISTFILE = ", bihistfile), file = fssd, sep = "\n", append = T)
    # number of classes to use
    cat(paste0("NCLASSES  = ", nclasses), file = fssd, sep = "\n", append = T)
    # auxiliary image
    cat(paste0("AUXILIARYFILE  = ", auxfile), file = fssd, sep = "\n", append = T)
  }
  
  cat("#---------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define the debug parameters - probably you wont need this        #", file = fssd, sep = "\n", append = T)
  cat("#---------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[DEBUG]", file = fssd, sep = "\n", append = T)
  # 1 to 3, use higher than 1 only if REALLY needed
  cat("DBGLEVEL = 2", file = fssd, sep = "\n", append = T)
  # File to write debug
  cat("DBGFILE   = debug.dbg ", file = fssd, sep = "\n", append = T)
  
  cat("#---------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define parameters for COVARIANCE TABLE - reduce if memory is a problem       #", file = fssd, sep = "\n", append = T)
  cat("#---------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[COVTAB]", file = fssd, sep = "\n", append = T)
  cat(paste0("MAXCTX = ", range), file = fssd, sep = "\n", append = T)
  cat(paste0("MAXCTY = ", range), file = fssd, sep = "\n", append = T)
  cat("MAXCTZ = 1", file = fssd, sep = "\n", append = T)
  
  cat("#--------------------------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  cat("#        here we define parameters for BLOCK KRIGING - if you are not block kriging useblocks should be 0      #", file = fssd, sep = "\n", append = T)
  cat("#--------------------------------------------------------------------------------------------------------------#", file = fssd, sep = "\n", append = T)
  
  cat("[BLOCKS]", file = fssd, sep = "\n", append = T)
  cat("USEBLOCKS = 1", file = fssd, sep = "\n", append = T)
  cat(paste0("BLOCKSFILE = ", bf), file = fssd, sep = "\n", append = T)
  cat(paste0("MAXBLOCKS  = ", maxblocks), file = fssd, sep = "\n", append = T)
  cat("[PSEUDOHARD]", file = fssd, sep = "\n", append = T)
  # 1 use, 0 no  pseudo hard data is point distributions that are simulated before all other nodes
  cat("USEPSEUDO = 0", file = fssd, sep = "\n", append = T)
  # file
  cat(paste0("PSEUDOFILE = ", "no file"), file = fssd, sep = "\n", append = T)
  # correct simulated value with point
  cat("PSEUDOCORR = 0", file = fssd, sep = "\n", append = T)
  
  # ------ run dss ------
  
  wd = getwd()
  parfile = paste0(day,ssd_nameO)
  setwd("./input")
  system(paste0("./DSS.C.64.exe ", parfile))
  setwd(wd)

}

