

# ----- description -----

# function to transform simulations .out - .grd (grid, raster)

# as input you should provide the list returned by funtion blockfile().

# ----- arguments ------

# blockobj, string, name of list, output of function blockfile() 
# emaps, logical, if TRUE computes e-type and uncertainty maps, plots e-type map 

# ------ function ------

 outraster = function (blockobj, grids = F, emaps = T) {
   library(raster)
   
   day = blockobj[["file"]]["day"]
   folder = blockobj[["file"]]["folder"]
   
   # create prefix .out filenames
   simout = paste0(day, "sim_")
   
   # store list .out namefiles
   lf = list.files(paste0(folder,"/"), pattern ="\\.out$")
   dss_list = list(simnames = Filter(function(x) grepl(simout, x), lf))


   # store number of simulations
   nsims = length(dss_list[["simnames"]])
   ssims = stack()
   namesims = c(rep(0, nsims))

   # loop each simulation
   for (k in 1:nsims){
     print(dss_list[["simnames"]][k])
     out = read.table(file = paste0(folder, "/", dss_list[["simnames"]][k]), sep = " ", skip=3)
     out01 = as.data.frame(out)
     out01[out01 == -999] <- NA
     
     # set grid size
     # mc : nr columns, mr : nr rows
     mc = blockobj[["gridpars"]][["nodes"]][1]
     mr = blockobj[["gridpars"]][["nodes"]][2]
     
     xmin = blockobj[["ingrid"]]@extent[1]
     xmax = blockobj[["ingrid"]]@extent[2]
     ymin = blockobj[["ingrid"]]@extent[3]
     ymax = blockobj[["ingrid"]]@extent[4]
     
     crsname = as.character(blockobj[["ingrid"]]@crs)
     
     # create matrix
     out02 = matrix(0, nrow = mr, ncol = mc)
     
     count = 1
     for (i in 1:mr){
       for (j in 1:mc){
         out02[i, j]=out01[count, 1]
         count = count + 1
         }
     }
     
     # create matrix
     out03 = matrix(NA, nrow = mr, ncol = mc)
     
     # populate matrix
     for (i in 1:mr){
       out03[i, ]<- c(out02[mr - i + 1,] )
     }
     r = raster(out03)
     extent(r)=c(xmin, xmax, ymin, ymax)
     projection(r) = CRS(crsname)
     
     if(grids == T){
       writeRaster(r, filename = paste0(folder,"/", day, "sim_", k), overwrite = TRUE)
     }
     # add layer to stack
     ssims = stack(ssims, r)
     # change layer name
     names(ssims[[k]]) = paste0("sim_", k) 
   }
   etype = calc(ssims, fun = function(x) {quantile(x, probs = .5,na.rm=TRUE)})
   uncer = calc(ssims, fun = sd, na.rm = T)
   if (emaps == T){
     writeRaster(etype, filename = paste0(folder, "/", day, "etype"), overwrite = TRUE)
     writeRaster(uncer, filename = paste0(folder, "/", day, "uncertainty"), overwrite = TRUE)
   }
   listmaps = list(simulations = ssims, etype = etype, uncertainty = uncer )
   return(listmaps)
   }


