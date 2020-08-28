
# ----- description -----

# function for rates (/10^4) and variance-error terms

# as input you should provide a dataframe with covid19 data.
# columns should include the following data:

# id of region, x, y and z cartesian coordinates at region mass center,
# number of covid19 cases by region and population size by region

# ----- arguments ------

# dfobj, string, dataframe name with covid19 data
# id, character, field name for region id
# x, character, field name for x-coordinates 
# y, character, field name for y-coordinates 
# z, character, field name for z-coordinates 
# cases, character, field name for number of cases 
# pop, character, field name for population size
# casesNA, numeric, an integer used to replace rows with cases = NA,
# day, character, string indicating date (format "yyyymmdd") of covid19 data

# ------ function ------


irates = function(dfobj = NA, oid = NA, xx = NA, yy = NA, zz = NA, 
                  cases = NA, pop = NA, casesNA = 2, day = "20200301") {

  # rate per phab habitants 
  phab = 10^4 
  
  # index variables id, x, y, z, cases, risk pop
  ioid = grep(oid, colnames(dfobj))
  ix = grep(xx, colnames(dfobj))
  iy = grep(yy, colnames(dfobj))
  iz = grep(zz, colnames(dfobj))
  ic = grep(cases, colnames(dfobj))
  ip = grep(pop, colnames(dfobj))
  
  # compute overall mean rate (exclude rows w/ cases = NA)
  dfobjnas = subset(dfobj, dfobj[, ic]!= "NA")
  n = ncol(dfobjnas)
  dfobjnas$rate = phab * dfobjnas[, ic] / dfobjnas[, ip]
  poptnas = sum(dfobjnas[, ip])
  m = sum(dfobjnas[, "rate"] * dfobjnas[, ip]) / poptnas
  
  # error variance term (m/n_i)
  error = m / dfobj[, ip]
  
  # NA cases set to casesNA
  dfobj[, ic] = ifelse(is.na(dfobj[, ic]), casesNA, dfobj[, ic])
  
  # recalculate crude rates
  rate = phab * dfobj[, ic] / dfobj[, ip]
  
  tab = data.frame (dfobj[, ioid], x = dfobj[, ix], y = dfobj[, iy], z = dfobj[, iz], rate, error)
  
  # cases file for dss
  
  # set folder 
  foldin = "input"
  
  # create folder for inputs
  if(!file.exists(foldin)) dir.create(foldin, recursive = F)  
  
  # store string with path for input files
  wkin = paste0(getwd(), "/", foldin)
  
  # prepare data to write file  
  tabnotf = tab[, c("x", "y", "z", "rate")]
  
  # store nr of variables
  nvars = ncol(tabnotf)
  
  # store nr of observations
  nobs = nrow(tabnotf)
  
  # store vector variable names
  namevars = names(tabnotf)

  # create file path
  not_name = "notified"
  not_nameO = paste0(not_name, ".out")
  fnot = paste0(wkin, "/", day, not_nameO)
  
  if (file.exists(fnot)){
    file.remove(fnot)  
  }
  
  # create notification file
  file.create(fnot)
  
  # write file header
  cat(not_name, file = fnot, sep="\n")
  cat(nvars, file = fnot, sep="\n", append = T)
  cat(namevars, file = fnot, sep="\n", append = T)
  
  # write notification data
  write.table(format(tabnotf, digits = NULL, justify = "right"),
              file = fnot, quote = F, append = T, row.names = F, col.names = F)
  
  # pars for ssdir.par
  xcol = grep("x", colnames(tabnotf))
  ycol = grep("y", colnames(tabnotf))
  zcol = grep("z", colnames(tabnotf))
  varcol = grep("rate", colnames(tabnotf))
  minval = min(tabnotf$rate)
  maxval = max(tabnotf$rate)
  
  # return data.frame object for variogram calcs
  tabvgm = data.frame (id = dfobj[, ioid], x = dfobj[, ix], y = dfobj[, iy], 
                       z = dfobj[, iz], rate, err = error, pop = dfobj[, ip])
  listf = list(day = day, name = paste0(day, not_nameO), folder = wkin)
  listpars = list(nvars = nvars, xcolumn = xcol, ycolumn = ycol, zcolumn = zcol, 
                  varcol = varcol, minval = minval, maxval = maxval)
  return(list(rates = tabvgm, mrisk = m, file = listf, ssdirpars = listpars))
}

