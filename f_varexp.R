varexp = function(dfobj, lag, nlags){

  # store nr of observations
  nobs = nrow(dfobj[["rates"]])
  
  # experimental variogram : get data
  rate = dfobj[["rates"]][,"rate"]
  ratexy = cbind(dfobj[["rates"]][c("x","y")])
  pop = dfobj[["rates"]][,"pop"]
  # store pop vector as double
  pop = as.double(pop)
  
  # weighted sample variance for sill estimation
  # no ref about this sill estimation
  # looked in goovaerts & cressie books, journel
  # at end, used fortran code from mjp
  
  # create integer num & den to calc weighted sample var
  nwsv = 0
  dwsv = 0
  
  # calc num and denominator 
  for (i in 1:nobs){
    nwsv = nwsv + pop[i]^2 / (2 * pop[i]) * (rate[i] - m)^2
    dwsv = dwsv + pop[i]^2 / (2 * pop[i])
  }
  
  # weighted sample variance
  wsvar = nwsv / dwsv
  
  # experimental variogram : calc distances
  # lag distance
  lagd = lag
  # cut distance
  lagend = lag * nlags
  # nr lags
  lagn = nlags + 1
  # vector of lags
  lags = c()
  for (i in 1:lagn) {
    lags[i] = lagd * i - lagd  
  }
  
  # compute distance matrix
  matd = as.matrix(dist(ratexy))
  
  # experimental variogram: create vectors to store results
  # n pairs dist(h)
  nh = vector( mode = "integer", length = (lagn-1)) 
  # total dist(h)
  dh = vector( mode = "numeric", length = (lagn-1)) 
  # mean dist(h)
  mh = vector( mode = "numeric", length = (lagn-1)) 
  # numerator gamma(h)
  num_gh = vector( mode = "double", length = (lagn-1)) 
  # denominator gamma(h)
  den_gh = vector( mode = "double", length = (lagn-1)) 
  # gamma(h)
  gammah = vector( mode = "double", length = (lagn-1)) 
  
  # experimental variogram: compute 
  for (k in 1:(lagn-1)) {
    for (i in 1:nobs) {
      for (j in 1:nobs) {
        if(matd[i,j] > lags[k] & matd[i,j] <= lags[k+1]){
          nh[k] = nh[k] + 1
          dh[k] = dh[k] + matd[i,j]
          num_gh[k] = num_gh[k] + ((pop[i] * pop[j]) * (rate[i]-rate[j])^2 - m) / (pop[i] + pop[j])
          den_gh[k] = den_gh[k] + pop[i] * pop[j] / (pop[i] + pop[j])
        }
      }
    }
    gammah[k] = num_gh[k]/(2*den_gh[k])
    mh[k] = dh[k] / nh[k]
  }
  v = data.frame(dist = mh,  semivariance = gammah)
  list(weightsvar = wsvar, semivar = v)
}  
 
