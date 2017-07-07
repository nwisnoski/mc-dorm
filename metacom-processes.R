# Functions that perform metacommunity processes

death <- function(MC, d){
  # randomly assign deaths with prob = d
  deaths = matrix(
    rbinom(n = prod(dim(MC$abunds)), size = 1, prob = 1-d),
    nrow = nrow(MC$abunds))
  
  # weight each death by deviation from env. optima (to do)
  
  MC$abunds = MC$abunds * deaths # element-wise multiplication by deaths (zeros)
  return(MC)
}

dispersal <- function(MC, disp){
  # randomly choose individuals to disperse
  
  dispersers = matrix(
    rbinom(n = prod(dim(MC$abunds)), size = 1, prob = disp),
    nrow = nrow(MC$abunds))
  
  MC$abunds = MC$abunds * (1-dispersers) # remove dispersers
  
  return(dispersers)
}

recruitment <- function(MC, dispersers){
  MC$abunds = MC$abunds + dispersers
  births = matrix(
    rbinom(n = prod(dim(MC$abunds)), size = 1, prob = b),
    nrow = nrow(MC$abunds))
  MC$abunds = MC$abunds + births
  return(MC)
}
