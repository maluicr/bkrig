
irates = function(df = NA, fid = NA, fx = NA, fy = NA, fz = NA, fcases = NA, fpop = NA, fday = "20200301"){

  # rate per phab habitants
  phab = 10^4 
  
  # file ascii tab delimited
  df = read.table(df)

  # variables id, x, y, z, cases, risk pop
  id = df[, fid]
  cx = df[, fx]
  cy = df[, fy]
  cz = df[, fz]
  c = df[, fcases]
  p = df[, fpop]
  
  # crude rates
  rate = phab * c / p
  
  # compute overall mean rate
  pt = sum(p, na.rm = T)
  m = sum(rate * p, na.rm = T) / pt
  
  # error variance term (m/n_i)
  error = m / p
  
  # NA cases set to 2 cases
  df[, fcases] = ifelse(is.na(df[, fcases]), 2, df[, fcases])
  c = df[, fcases]
  
  # recalculate crude rates
  rate = phab * c / p
  
  tab = data.frame (id, x = cx, y = cy, z = cz, rate, error)
  
  # cases file for dss
  
  # set folder 
  foldin = "input"
  
  # create folder for inputs
  if(!file.exists(foldin)) dir.create(foldin, recursive=T)  
  
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
  fnot = paste0(wkin, "/", fday, not_nameO)
  
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
  
  # return data.frame object for variogram calcs
  tabvgm = data.frame (x = cx, y = cy, rate, pop = p)
  return(list(date = fday, filename = fnot, rates = tabvgm))
}

