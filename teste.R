
# for runbkrig vignette

# create a data frame from covid19 data table
covid = read.table("covid19_data.txt", header = TRUE, sep = "\t", dec = ".")

# compute rates and error variance
rates = irates(dfobj = covid, oid = "id_region", xx = "xcoord", yy = "ycoord", zz = "zcoord", 
               cases = "ncases", pop = "poprisk", casesNA = 2, day = "20200601")

# create block file
block = blockfile(rates, "grid2k.tif")
# plot grid input
plot(block$ingrid)

# create mask file
mask = maskfile(block)

# compute experimental variogram
vexp = varexp(rates, lag = 7000, nlags = 15)

# compute theoretical variogram
vmod = varmodel(vexp, mod = "Exp", nug = 0, ran = 60000, sill = vexp[["weightsvar"]])

# plot experimental variogram
plot(vexp[["semivar"]], ylab = expression(paste(gamma, "(h)")), xlab = "h (in m)") 
# add sill
abline(h = vexp[["weightsvar"]], col ="red", lty = 2)
# add theoretical model
lines(vmod[["fittedval"]]) 

# run ssdir.par
ssdpars(blockobj = block, maskobj = mask, dfobj = rates, varmobj = vmod, 
        simulations = 5, radius1 = 60000, radius2 = 60000)

# export .out maps to raster
maps = outraster(block, emaps = T)

# plot simulations
spplot(maps[["simulations"]])

# plot etype and uncertainty
spplot(maps[["etype"]])
spplot(maps[["uncertainty"]])