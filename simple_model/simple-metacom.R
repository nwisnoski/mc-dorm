library(tidyverse)
library(vegan)
library(progress)
library(vegetarian)
#library(viridis)
#library(adespatial)
#library(igraph)

###############################################################################
rm(list = ls())
set.seed(47405)

# Define model parameters
tsteps <- 50000      # Number of time steps in model
dt <- 1    # precision for model integration (step size)
M <- 20 # Number of sites
S <- 20 # Number of species
ext <- .01 # extinction thresh
disturb <- 0.0001


envs <- 1 # Number of environmental variables


# set up environment (E)
E <- matrix(seq(0, 1, length.out = M), nrow = M, ncol = envs)




###############################################################################
# Define metacommunity functions


# growth rate in a patch based on environmental match 
R.jx <- function(E, max.R, opt.envs, nbreadth){
  max.R * exp(-((E - opt.envs)^2/(2*nbreadth^2)))
}

competition <- function(N,a) {
  return( 1 / ( 1 + a*(sum(N))))
}

disperse <- function(d, N, M){
  
  dispersal_matrix <- matrix(1/(M-1), nrow = M, ncol = M)
  diag(dispersal_matrix) <- 0
  
  return((N * d) %*% dispersal_matrix)

}




###############################################################################
# Definte traits

# establish strength of species sorting/local control
opt.envs <- seq(0, 1, length.out = S) # species optimal environments
nbreadth <- .5 # as approach infinity, metacom approaches neutrality
max.R <- 1.2
a <- 4e-4 # strength of competition

d <- rep(.3, S) # Dispersal rates
decay <- rep(.00001, S) # Decay rate of dormant propagules
dorm <- rep(0.0, S) # Propensity to enter dormancy
activ <- rep(0.1, S) # Reactivation rate

###############################################################################

#d <- rep(0.001, S)
#disturb <- 0.001
disp.grad <- c(0, .001, .005, 0.01, 0.02, 0.05, 0.07, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
# disp.grad <- c(0, .005, 0.01, 0.05, 0.1, 0.5, 1)
# 
# # 1 - dispersal, 2-4 - alpha, beta, gamma, 
out.sum <- matrix(NA, nrow = length(disp.grad), ncol = 4)
i = 1
comms.out <- list(NA)
# loop over dipsersal rates
for(d in disp.grad){
  set.seed(47405)
  # initialize result array, Species X Sites X Time 
  out.N <- array(NA, c(S, M, tsteps), 
                 dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))
  out.D <- array(NA, c(S, M, tsteps),
                 dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))
  
  # add active species to sites, seed banks are empty
  out.N[,,1] <- 1
  out.D[,,1] <- 0
  
for(t in 1:(tsteps-1)){
  
  if(t == 1) pb <- progress_bar$new(total = tsteps, force = T)
  
  # update progress bar
  pb$update(ratio = t/tsteps)
  # get current abunds
  N.t <- out.N[,,t]
  D.t <- out.D[,,t]
  
  # calculates growth rates for each species (rows) in each site (cols)
  # dimensions = S x M
  R.t <- apply(X = E, MARGIN = 1, FUN = R.jx, max.R, opt.envs, nbreadth)
  
  # dormancy transitions
  # previous + dormancy - reactivation
  D.t1 <- D.t + (N.t * dorm) - (D.t * activ)
  # previous + reactivation - dormancy
  N.t1 <- N.t + (D.t * activ) - (N.t * dorm)
  
  # calculate competition at this time
  comp.t <- apply(N.t1, MARGIN = 2, FUN = competition, a = a)
  
  # growth and seed bank decay
  N.t2 <- R.t * N.t1 * comp.t
  D.t2 <- D.t1 - decay * D.t1
  
  # calculate immigration and then remove emigrants
  N.t3 <- N.t2 + disperse(d, N.t2, M) - (d * N.t2)
  D.t3 <- D.t2 + disperse(d, D.t2, M) - (d * D.t2)
  
  # patch disturbance
  N.t3[, which(rbernoulli(M, p = disturb) == 1)] <- 0
  
  out.N[,,t+1] <- ifelse(N.t3 > ext, N.t3, 0)
  out.D[,,t+1] <- ifelse(D.t3 > ext, D.t3, 0) 
  
}


# Creat path for sim output
sim.path <- file.path("figures", "sim_output", 
            paste0("dispersal", mean(d), "-dorm",mean(dorm),"-act",mean(activ),"-dist",mean(disturb)))
if(!dir.exists(sim.path)) dir.create(sim.path, recursive = T)
# extract SxS matrix
comm <- t(out.N[,,tsteps])
saveRDS(comm, file = file.path(sim.path, "model-output.rds"))
comms.out[[i]] <- comm
comm[is.na(comm)] <- 0
comm <- decostand(comm, method = "hellinger")
comm

#plot patch 6 just to show local dynamics example
autoplot.zoo(t(out.N[,6,]), facet = NULL) + 
  labs(x = "Time", y = "Abundance") + 
  theme_bw() +
  ggsave(file.path(sim.path, "local_dynamics.png"), width = 8, height = 6, units = "in", dpi = 500)

# specnumber(comm)
# specnumber(colSums(comm))
# ord <- rda(comm ~ E); plot(ord)
#percent.env <- as.numeric(eigenvals(ord)[1]/sum(eigenvals(ord)))
(alpha <- vegetarian::d(comm, lev = "alpha"))
(beta <- vegetarian::d(comm, lev = "beta"))
(gamma <- vegetarian::d(comm, lev = "gamma"))
out.sum[i,] <- c(d, alpha, beta, gamma)
i <- i + 1
}

#   comms.out[[i]] <- t(out.N[,,tsteps])
#   alpha <- mean(specnumber(t(out.N[,,tsteps])))
#   gamma <- specnumber(colSums(t(out.N[,,tsteps])))
#   beta <- round(gamma / alpha, 2)
#   div.part[i,] <- c(d, alpha, beta, gamma)
#   i <- i + 1
# }

out.sum
colnames(out.sum) <- c("dispersal", "alpha", "beta", "gamma")
as.data.frame(out.sum) %>% 
  gather(alpha, beta, gamma, key = scale, value = diversity) %>%
  ggplot(aes(dispersal, diversity, color = scale)) +
  geom_point(size = 2, alpha = 0.5) +
  geom_line() + 
  theme_bw() +
  ggsave(file.path("figures", "diversity-dispersal", 
paste0("dorm",mean(dorm),"-act",mean(activ),"-dist",mean(disturb),".png")), 
         width = 4, height = 3, units = "in", dpi = 500)
