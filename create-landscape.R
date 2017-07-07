# Code to initiate a metacommunity with m patches on a lxl landscape


get.coords <- function(m, l){
  x = runif(n = m, min = 0, max = l)
  y = runif(n = m, min = 0, max = l)
  return(cbind(x,y))
}

set.random.environment <- function(m, l, e){
  E = matrix(NA, nrow = m, ncol = e)
  for(env in 1:ncol(E)) E[,env] = runif(n = m)
  return(E)
}

init.species <- function(m, s, e){
  # initialize metacom with prob occupancy of 0.5
  MC = matrix(rbinom(n = s * m, size = 1, prob = 0.5), nrow = m, ncol = s)
  
  # set species environmental prefs
  S.env = matrix(NA, nrow = s, ncol = e)
  for(env in 1:ncol(S.env)) S.env[,env] = runif(n = s)
  
  return(list(abunds = MC, envprefs = S.env))
  
}
