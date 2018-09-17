library(tidyverse)
library(vegan)
library(progress)
library(viridis)
library(adespatial)
library(igraph)

###############################################################################
rm(list = ls())
set.seed(47405)

# Define model parameters
tsteps <- 5000      # Number of time steps in model
dt <- 1    # precision for model integration (step size)
M <- 20 # Number of sites
S <- 20 # Number of species
ext <- .01 # extinction thresh


envs <- 1 # Number of environmental variables

# Create Site by Species Matrix (N = active, D = seed bank)
N <- matrix(0, ncol = S, nrow = M)
D <- N * 0
N <- N + 0.1


# set up environment (E)
E <- matrix(seq(0, 1, length.out = M), nrow = M, ncol = envs)


# initialize result array, Species X Sites X Time 
out.N <- array(NA, c(S, M, tsteps), 
               dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))
out.D <- array(NA, c(S, M, tsteps),
               dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))

# add active species to sites, seed banks are empty
out.N[,,1] <- N
out.D[,,1] <- D

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
nbreadth <- .8 # as approach infinity, metacom approaches neutrality
max.R <- 1.2
a <- 4e-4 # strength of competition

d <- rep(.1, S) # Dispersal rates
decay <- rep(.0001, S) # Decay rate of dormant propagules

###############################################################################
# Run Simulation
# sim <- function(disp.grad){
#   div.part <- matrix(NA, nrow = length(disp.grad), ncol = 4)
#   i = 1
#   for(d in disp.grad){
#     
#     for(t in 1:(tsteps-1)){
#       
#       # get current abunds
#       N.t <- out.N[,,t]
#       D.t <- out.D[,,t]
#       
#       # calculates growth rates for each species (rows) in each site (cols)
#       # dimensions = S x M
#       R.t <- apply(X = E, MARGIN = 1, FUN = R.jx, max.R, opt.envs, nbreadth)
#       
#       # calculate competition at this time
#       comp.t <- apply(N.t, MARGIN = 2, FUN = competition, a = a)
#       
#       # growth and seed bank decay
#       N.t0 <- R.t * N.t * comp.t
#       D.t0 <- D.t - decay * D.t
#       
#       # calculate immigration and then remove emigrants
#       N.t1 <- N.t0 + disperse(d, N.t0, M) - (d * N.t0)
#       D.t1 <- D.t + disperse(d, D.t0, M) - (d * D.t0)
#       
#       out.N[,,t+1] <- ifelse(N.t1 > ext, N.t1, 0)
#       out.D[,,t+1] <- ifelse(D.t1 > ext, D.t1, 0) 
#       
#     }
#     
#     alpha <- mean(specnumber(t(out.N[,,tsteps])))
#     gamma <- specnumber(colSums(t(out.N[,,tsteps])))
#     beta <- round(gamma / alpha, 2)          
#     div.part[i,] <- c(d, alpha, beta, gamma)
#     i <- i + 1
#   }
#   return(div.part)
# }
# 
# 
# disp.grad <- seq(0.0001, 1, length.out = 10)
# divsum <- sim(disp.grad)
# colnames(divsum) <- c("dispersal", "alpha", "beta", "gamma")
# as.data.frame(divsum) %>%
#   gather(alpha, beta, gamma, key = scale, value = diversity) %>% 
#   ggplot(aes(dispersal, diversity, color = scale)) + 
#   geom_point() + 
#   geom_line()


d = rep(.15, S)
for(t in 1:(tsteps-1)){
  
  # get current abunds
  N.t <- out.N[,,t]
  D.t <- out.D[,,t]
  
  # calculates growth rates for each species (rows) in each site (cols)
  # dimensions = S x M
  R.t <- apply(X = E, MARGIN = 1, FUN = R.jx, max.R, opt.envs, nbreadth)
  
  # calculate competition at this time
  comp.t <- apply(N.t, MARGIN = 2, FUN = competition, a = a)
  
  # growth and seed bank decay
  N.t0 <- R.t * N.t * comp.t
  D.t0 <- D.t - decay * D.t
  
  # calculate immigration and then remove emigrants
  N.t1 <- N.t0 + disperse(d, N.t0, M) - (d * N.t0)
  D.t1 <- D.t + disperse(d, D.t0, M) - (d * D.t0)
  
  out.N[,,t+1] <- ifelse(N.t1 > ext, N.t1, 0)
  out.D[,,t+1] <- ifelse(D.t1 > ext, D.t1, 0) 
  
}

matplot(t(out.N[,,]), type = 'l')
t(out.N[,,tsteps])
