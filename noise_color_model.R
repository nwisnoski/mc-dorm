library(tidyverse)
library(vegan)
library(progress)
library(vegetarian)
library(zoo)
library(viridis)
#library(adespatial)
#library(igraph)

###############################################################################
rm(list = ls())
set.seed(47405)

# Define model parameters
tsteps <- 2500      # Number of time steps in model
dt <- 1    # precision for model integration (step size)
M <- 5 # Number of sites
S <- 20 # Number of species
ext <- .01 # extinction thresh
disturb <- 0.001
env.type <- "fluctuating" # "static", "fluctuating", "random" okay  
spatial.synchrony <- 0 #range from 0 to 1, what fraction of patches have the same environment

envs <- 1 # Number of environmental variables


# set up environment (E)
E <- matrix(seq(0, 1, length.out = M), nrow = M, ncol = envs)
env.ampl <- 1 # amplitude of env variability [0,1]
env.period <- 300 # bigger numbers, slower oscilations, more static env

E.j <- function(t, synch = spatial.synchrony){
  0.5*env.ampl*(
    sin((1:M)*2*pi/(M*(1-synch)) + (2*pi/env.period)*t) + 1)
  
}

env.dyn <- matrix(NA, nrow = tsteps, ncol = M)
for(t in 1:tsteps){env.dyn[t,] = E.j(t)}

autoplot.zoo(env.dyn, facet = NULL) + theme_minimal()

autoplot.zoo(arima.sim(model = list(ar = .7), n = (tsteps-1), sd = 0.1), facet = NULL) + theme_minimal()
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

d <- rep(.05, S) # Dispersal rates
decay <- rep(.000001, S) # Decay rate of dormant propagules
dorm <- rep(0.5, S) # Propensity to enter dormancy
activ <- rep(0.1, S) # Reactivation rate
ddcov <- 0 # 0 for negative, 1 for positive

###############################################################################
run.sim <- function(noise.col){
  set.seed(47405)
  # initialize result array, Species X Sites X Time 
  out.N <- array(NA, c(S, M, tsteps), 
                 dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))
  out.D <- array(NA, c(S, M, tsteps),
                 dimnames = list(c(paste0("sp",1:S)), c(paste0("site",1:M)), c(1:tsteps)))
  
  # add active species to sites, seed banks are empty
  out.N[,,1] <- 1
  out.D[,,1] <- 0
  
  env.noise <- arima.sim(model = list(ar = noise.col), n = tsteps, sd = 0.1)
  for(t in 1:(tsteps-1)){
    
    # if(t == 1) pb <- progress_bar$new(total = tsteps, force = T)
    
    # update progress bar
    # pb$update(ratio = t/(tsteps))
    # get current abunds
    N.t <- out.N[,,t]
    D.t <- out.D[,,t]
    
    # calculates growth rates for each species (rows) in each site (cols)
    # dimensions = S x M
    if(env.type == "fluctuating") E <- matrix(E.j(t), nrow = M, ncol = 1)
    if(env.type == "random") E <- matrix(runif(M), nrow = M, ncol = 1)
    R.t <- apply(X = (E + env.noise[t]), MARGIN = 1, FUN = R.jx, max.R, opt.envs, nbreadth)
    
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
    D.t3 <- D.t2 + ddcov*(disperse(d, D.t2, M) - ddcov*(d * D.t2))
    
    # patch disturbance
    N.t3[, which(rbernoulli(M, p = disturb) == 1)] <- 0
    
    out.N[,,t+1] <- ifelse(N.t3 > ext, N.t3, 0)
    out.D[,,t+1] <- ifelse(D.t3 > ext, D.t3, 0) 
  }
  return(out.N)
}


######
noise.grad <- seq(-.99,.99, length.out = 10)
dorm.grad <- c(0, .5)
disp.grad <- seq(0, .4, length.out = 50)
# # 1 - noisecol, 2 - dispersal, 3 - dorm, 4-6 - alpha, beta, gamma, 
out.sum <- matrix(NA, nrow = length(noise.grad)*length(disp.grad)*length(dorm.grad), ncol = 6)
i = 1
loops <- length(out.sum)

for(noise.col in noise.grad){
  # loop over dorm rates
  for(dorm in dorm.grad){
    
    # loop over dipsersal rates
    for(d in disp.grad){
      
      if(i == 1) pb <- progress_bar$new(total = loops, force = T)
      
      # update progress bar
      pb$update(ratio = i/loops)
      
      out.N <- run.sim(noise.col = noise.col)
      comm <- t(out.N[,,tsteps])
      comm[is.na(comm)] <- 0
      #comm <- decostand(comm, method = "hellinger")
      
      (alpha <- vegetarian::d(comm, lev = "alpha", q = 0))
      (beta <- vegetarian::d(comm, lev = "beta", q = 0))
      (gamma <- vegetarian::d(comm, lev = "gamma", q = 0))
      
      # write out dispersal, dormancy, and diversity
      out.sum[i,] <- c(noise.col, d, dorm, alpha, beta, gamma)
      i <- i + 1
    }
  }
}

out.sum
colnames(out.sum) <- c("AR", "Dispersal", "Dormancy", "Alpha", "Beta", "Gamma")
as.data.frame(out.sum) %>% 
  gather(Alpha, Beta, Gamma, key = Scale, value = Diversity) %>%
  mutate(AR_fact = factor(AR, levels = noise.grad, ordered = T)) %>% 
  mutate(Shift = ifelse(AR < 0, "Blue-shifted", "Red-shifted")) %>% 
  filter(Scale == "Gamma") %>% 
  mutate(Dormancy = factor(Dormancy, labels = c("No Dormancy", "Dormancy"))) %>% 
  ggplot(aes(Dispersal, Diversity, group = AR_fact)) +
  geom_point(aes(color = AR), alpha = 0.8) +
  geom_line(aes(color = AR), size = .75, alpha = 0.7) +
  facet_grid(Dormancy ~ Shift) +
  theme_minimal() + 
  scale_x_continuous(limits = c(0,.4)) +
  scale_color_gradientn("Noise Color", colours = colorspace::darken(rev(RColorBrewer::brewer.pal(name = "RdBu", n = 11)), 0)) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "grey80"),
        legend.position = c(.9, .25), strip.text = element_text(size = 16),
        axis.title = element_text(size = 18)) +
  ggsave("figures/env_noise_color.png", dpi = 600, width = 6, height = 3/4*6)

as.data.frame(out.sum) %>% select(-Beta, -Gamma) %>%
  rename(Diversity = Alpha) %>%
  mutate(Dormancy = factor(Dormancy, labels = c("No Dormancy", "Dormancy"))) %>% 
  ggplot(aes(x = Dispersal, y = AR, z = Diversity)) +
  geom_raster(aes(fill = Diversity)) +
  facet_wrap(~Dormancy) +
  scale_fill_viridis() +
  theme_minimal() +
  labs(y = "Noise color (AR parameter)") +
  ggsave("figures/env_noise_color_countour.png", dpi = 600, width = 4, height = 3)
  #geom_contour(aes(color = stat(level))) +