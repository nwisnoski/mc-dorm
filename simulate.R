# This function runs through a metacommunity simulation
library(vegan)
source("parameters.R")
source("create-landscape.R")
source("metacom-processes.R")

XY <- get.coords(m, l)
XY.env <- set.random.environment(m, l, e)
MC <- init.species(m, s, e)
MC$envprefs


# step through one metacommunity run
timesteps <- 1000
out <- list()
for(i in 1:timesteps){
  MC = death(MC, d)
  dispersers <- dispersal(MC, disp)
  MC = recruitment(MC, dispersers)
  out[[i]] <- MC$abunds
}

rowSums(out[[1000]])

div <- lapply(X = out, FUN = function(x) rowSums(decostand(x, method = "pa")))
div <- matrix(unlist(div), nrow = m)
ts.div <- t(div)
plot(1:timesteps, ts.div[,9], type = "l")
