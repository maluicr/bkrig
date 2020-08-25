

irates = function(df = NA, id = NA, x = NA, y = NA, z = NA, 
                  cases = NA, pop = NA, casesNA = 2, day = "20200301") {

  # rate per phab habitants
  phab = 10^4 
  
  # file ascii tab delimited
  df = read.table(df)

  # variables id, x, y, z, cases, risk pop
  id = df[, id]
  cx = df[, x]
  cy = df[, y]
  cz = df[, z]
  c = df[, cases]
  p = df[, pop]
  
  # compute overall mean rate (exclude rows w/ cases = NA)
  df_NAs = subset(df, df[, cases]!= "NA")
  df_NAs$rates = phab * df_NAs[, cases] / df_NAs[, pop]
  pop_tot_NAs = sum(df_NAs[, pop])
  m = sum(df_NAs[, rates] * df_NAs[, pop]) / pop_tot_NAs
  
  # error variance term (m/n_i)
  error = m / p
  
  # NA cases set to 2 cases
  df[, cases] = ifelse(is.na(df[, cases]), casesNA, df[, cases])
  c = df[, cases]
  
  # recalculate crude rates
  rate = phab * c / p
  
  tab = data.frame (id, x = x, y = y, z = z, rate, error)
  
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
  tabvgm = data.frame (oid = id, x = cx, y = cy, z = cz, rate, err = error, pop = p)
  listf = list(day = day, name = paste0(day, not_nameO), folder = wkin)
  listpars = list(xcolumn = xcol, ycolumn = ycol, zcolumn = zcol, varcol = varcol, minval = minval, maxval = maxval)
  return(list(rates = tabvgm, file = listf, ssdirpars = listpars))
}

