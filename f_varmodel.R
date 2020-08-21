
varmodel = function (varexp ,mod = c("Exp", "Sph"), nug, ran , sill) {
  x = seq(1, ran , 1)
  xmax = seq(1, max(varexp[[2]][,1]), 1)
  if (mod=="Sph") {
    model = nug + sill * (1.5 * (x / ran) - 0.5 * (x / ran)^3)
    cutoff = max(xmax) - ran
    msill = rep(nug + sill, cutoff )
    model = c(0, model, msill)
  }
  if (mod=="Exp") {
    model = nug + sill * (1 - exp(- 3 * xmax / ran))
    model = c(0, model)
  }
  pars = data.frame (model = mod, nugget = nug, range = ran, psill = sill)
  nstruc = 1 # for now only 1 (Spherical or Exponential)
  return(list(structures = nstruc, parameters = pars, fittedval = model))
}


