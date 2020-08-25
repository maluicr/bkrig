
# script example to compute block kriging covid-19 (Azevedo et al, 2020)

# compute rates and error variance
tx = irates(df = "dados.txt", id = "objectid", x = "x", y = "y", z = "z", 
            cases = "casos", pop = "popres18", casesNA = 2, day = "20200302")

# create block file
rgrid = blockfile(tx, "gis/grid2k.tif")
plot(rgrid$ingrid)

# create mask file
m = maskfile(rgrid)

# compute experimental variogram
ve = varexp(tx, lag = 7000, nlags = 15)

# compute theoretical variogram
vm = varmodel(ve, mod = "Exp", nug = 0, ran = 60000, sill = ve[[1]])

# plot experimental variogram
plot(ve[["semivar"]], ylab = expression(paste(gamma, "(h)")), xlab = "h (in m)") 
# add sill
abline(h = ve[["dist"]], col ="red", lty = 2)
# add theoretical model
lines(vm[["fittedval"]]) 

# run ssdir.par
ssdpars(rgrid, m, tx, vm, simulations = 2, radius1 = 60000, radius2 = 60000 )