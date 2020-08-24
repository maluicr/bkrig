
# compute rates and error variance
tx = irates(df = "dados.txt", fid = "objectid", fx = "x", fy = "y", fz = "z", fcases = "casos", fpop = "popres18", fday = "20200302")

# create block file
rgrid = blockfile(tx, "gis/grid2k.tif")
plot(rgrid$ingrid)

# create mask file
mask = maskfile(rgrid)

# compute experimental variogram
ve = varexp(tx, lag = 7000, nlags = 15)

# compute theoretical variogram
vm = varmodel(ve, mod = "Exp", nug = 0, ran = 60000, sill = ve[[1]])

# plot experimental variogram
plot(ve[[2]], ylab = expression(paste(gamma, "(h)")), xlab = "h (in m)") 
# add sill
abline(h = ve[[1]], col ="red", lty = 2)
# add theoretical model
lines(vm[[3]]) 

