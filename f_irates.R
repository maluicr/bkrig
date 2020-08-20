
irates = function(df = NA, fid = NA, fx = NA, fy = NA, fz = NA, fcases = NA, fpop = NA){
  
  # rate per phab habitants
  phab = 10^4 
  
  # file
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
  c = df[, fcases]
  
  # recalculate crude rates
  rate = phab * c / p
  
  tab = data.frame (id, cx, cy, cz, rate, error)
  return(tab)
}