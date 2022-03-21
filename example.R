
# author: m ribeiro; 
# date : 26-08-20; last revision : 21-03-2022
# email: manuel.ribeiro@tecnico.ulisboa.pt

# Please note: some packages may need to be installed!

# script example to compute block kriging covid-19 (Azevedo et al, 2020)


library(blockdss)
data(ptdata)
data(ptgrid)

# compute rates and error variance
rates = irates(dfobj = ptdata, oid = "oid_", xx = "x", yy = "y", zz = "t",
               cases = "ncases", pop = "pop19", casesNA = 1, day = "2021015")

coordinates(ptgrid) <- ~x+y
proj4string(ptgrid) <- CRS("+init=epsg:3763")
class(ptgrid)
gridded(ptgrid) <- T
class(ptgrid)

# create block file
block = blockfile(rates, ptgrid)

# plot grid input
plot(block$ingrid)

# create mask file
mask = maskfile(block)

# compute experimental variogram
vexp = varexp(rates, lag = 7000, nlags = 25)

# compute theoretical variogram
vmod = varmodel(vexp, mod = "sph", nug = 0, ran = 35000, sill = vexp[["weightsvar"]])

# plot experimental variogram
plot(vexp[["semivar"]][1:2], ylab = expression(paste(gamma, "(h)")), xlab = "h (in m)", main = "Semi-variogram") 

# add sill
abline(h = vexp[["weightsvar"]], col ="red", lty = 2)
# add theoretical model
lines(vmod[["fittedval"]]) 

# run ssdir.par
ssdpars(blockobj = block, maskobj = mask, dfobj = rates, varmobj = vmod, 
        simulations = 5, radius1 = 35000, radius2 = 35000)

# export .out maps to raster
maps = outraster(block, emaps = T)

# plot simulations
spplot(maps[["simulations"]])

# plot etype and uncertainty
spplot(maps[["etype"]])
spplot(maps[["uncertainty"]])
