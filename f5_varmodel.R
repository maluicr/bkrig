
# ----- description -----

# calculates a variogram model

# as input you should provide:

# an object list returned by funtion varexp(), 
# the variogram model and the variogram parameters.

# ----- arguments ------

# varexp, character, name of object, output of function varexp() 
# mod, character, the variogram model type (available are: "Sph" or "Exp")
# nug, numeric, nugget-effect value of the variogram
# ran, numeric, range value of the variogram
# sill, numeric, sill (or partial sill) value of the variogram

# ------ function ------



varmodel = function (varexp, mod = c("Exp","Sph"), nug, ran , sill) {
  
  x = seq(1, ran , 1)
  xmax = seq(1, max(varexp[["semivar"]][,"dist"]), 1)
  if (mod=="Sph") {
    vtype = 1
    model = nug + sill * (1.5 * (x / ran) - 0.5 * (x / ran)^3)
    cutoff = max(xmax) - ran
    msill = rep(nug + sill, cutoff )
    model = c(0, model, msill)
  }
  if (mod=="Exp") {
    vtype = 2
    model = nug + sill * (1 - exp(- 3 * xmax / ran))
    model = c(0, model)
  }
 
  pars = data.frame (model = mod, modeltype = vtype, nugget = nug, range = ran, psill = sill)
  nstruc = 1 # for now only 1 (Spherical or Exponential)
  return(list(structures = nstruc, parameters = pars, fittedval = model))
}


