
# script example to compute block kriging covid-19 (Azevedo et al, 2020)

source("f_irates.R", echo = T)
source("f_blockfile.R", echo = T)
source("f_maskfile.R", echo = T)
source("f_varexp.R", echo = T)
source("f_varmodel.R", echo = T)
source("f_ssdir.R", echo = T)

# compute rates and error variance
covid = irates(df = "dados.txt", id = "objectid", x = "x", y = "y", z = "z", 
            cases = "casos", pop = "popres18", casesNA = 2, day = "20200415")

# create block file
block = blockfile(covid, "gis/grid2k.tif")
# plot grid input
plot(block$ingrid)

# create mask file
mask = maskfile(block)

# compute experimental variogram
vexp = varexp(covid, lag = 7000, nlags = 15)

# compute theoretical variogram
vmod = varmodel(vexp, mod = "Exp", nug = 0, ran = 60000, sill = vexp[["weightsvar"]])

# plot experimental variogram
plot(vexp[["semivar"]], ylab = expression(paste(gamma, "(h)")), xlab = "h (in m)") 
# add sill
abline(h = vexp[["weightsvar"]], col ="red", lty = 2)
# add theoretical model
lines(vmod[["fittedval"]]) 

# run ssdir.par
ssdpars(blockobj = block, maskobj = mask, dfobj = covid, varmobj = vmod, 
        simulations = 5, radius1 = 60000, radius2 = 60000)
