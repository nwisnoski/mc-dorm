library(gstat)
library(fields)
library(sp)

# Construct the landscape with X environmental vars
# set degree of spatial autocorrelation (higher range is courser correlation)
# set linear trends with model parameters
# http://santiago.begueria.es/2010/10/generating-spatially-correlated-random-fields-with-r/

stdize <- function(x, ...) {(x - min(x, ...)) / (max(x, ...) - min(x, ...))}

make.landscape <- function(xmax = 100, ymax = 100, coarse = 15, xslope = .005, yslope = .005, envvars = 1, sites){
  xy <- expand.grid(1:xmax, 1:ymax)
  names(xy) <- c('x','y')
  g.dummy <- gstat(formula=z~1+x+y, locations=~x+y, dummy=T, 
                   beta=c(1,xslope,yslope), model=vgm(psill=0.025, range=coarse, model='Exp'), nmax=20)
  yy <- predict(g.dummy, newdata=xy, nsim=envvars)
  gridded(yy) = ~x+y
  spplot(obj=yy)
  
  # Generate random points along the landscape
  geo.coords <- cbind(x = sample(c(1:xmax), size = sites, replace = F), 
                      y = sample(c(1:ymax), size = sites, replace = F))
  geo.dists <- dist(geo.coords)
  
  local.envs <- matrix(NA, nrow = M, ncol = envvars)
  colnames(local.envs) <- paste0("env",seq(1,envvars))
  for(i in 1:envvars){
    env.structure <- as.matrix(yy[i])
    local.envs[,i] <- stdize(env.structure[floor(geo.coords)])
  }
  
  out <- list(coords = geo.coords, env = local.envs)
  
  return(out)
}


create.edge.list <- function(L){
  edge.list = matrix(NA, nrow = 0, ncol = 2)
  for(i in 1:length(L)){
    for(j in 1:length(L[[i]])){
      edge.list = rbind(edge.list, cbind(i, L[[i]][j]))
    }
  }
  return(edge.list)
}
