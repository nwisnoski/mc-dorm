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
tsteps <- 10000      # Number of time steps in model
dt <- 1    # precision for model integration (step size)
M <- 20 # Number of sites
S <- 20 # Number of species
ext <- .01 # extinction thresh
disturb <- 0.001
env.type <- "fluctuating" # "static", "fluctuating", "random" okay  
spatial.synchrony <- .999 #range from 0 to 1, what fraction of patches have the same environment

envs <- 1 # Number of environmental variables


# set up environment (E)
E <- matrix(seq(0, 1, length.out = M), nrow = M, ncol = envs)
env.ampl <- 1 # amplitude of env variability [0,1]
env.period <- 1000 # bigger numbers, slower oscilations, more static env

E.j <- function(t, synch = spatial.synchrony, stoch = F){
  0.5*env.ampl*(
    sin((1:M)*2*pi/(M*(1-synch)) + (2*pi/env.period)*t) + 1) + ifelse(stoch, abs(rnorm(1, sd = .1)), 0)
  
}

env.dyn <- matrix(NA, nrow = tsteps, ncol = M)
for(t in 1:tsteps){env.dyn[t,] = E.j(t)}
autoplot.zoo(env.dyn, facet = NULL) + theme_minimal()
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
decay <- rep(.00001, S) # Decay rate of dormant propagules
dorm <- rep(0.5, S) # Propensity to enter dormancy
activ <- rep(0.1, S) # Reactivation rate

ddcov <- 0 # 0 for negative, 1 for positive

###############################################################################

dorm.grad <- c(0, 0.01, 0.1, 0.5)
disp.grad <- c(0, 0.01, 0.1, 0.5, 1)
# dorm.grad <- c(0, .01, .1, 1)
# disp.grad <- c(0, 0.01, 0.1, 1)

# # 1 - dispersal, 2 - dorm, 3-5 - stability 
out.sum <- matrix(NA, nrow = length(disp.grad)*length(dorm.grad), ncol = 5)
i = 1
loops <- length(out.sum)
comms.out <- list(NA)

# loop over dorm rates
for(dorm in dorm.grad){
  
  # loop over dipsersal rates
  for(d in disp.grad){
    
    if(i == 1) pb <- progress_bar$new(total = loops, force = T)
    
    # update progress bar
    pb$update(ratio = i/loops)
    
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
      D.t3 <- D.t2 + ddcov*(disperse(d, D.t2, M) - (d * D.t2))
      
      # patch disturbance
      N.t3[, which(rbernoulli(M, p = disturb) == 1)] <- 0
      
      out.N[,,t+1] <- ifelse(N.t3 > ext, N.t3, 0)
      out.D[,,t+1] <- ifelse(D.t3 > ext, D.t3, 0) 
      
    }
    
    
    out.df <- as.data.frame.table(out.N, responseName = "abundance")
    colnames(out.df) <- c("species", "site", "time", "abundance")
    out.df$time <- as.numeric(out.df$time)
    
    if(sum(dorm)==0){
      stability.nodorm <- out.df %>% 
        community_stability(time.var = "time",
                            abundance.var = "abundance",
                            replicate.var = "site")
    }
    if(sum(dorm)!=0){
      stability.dorm <- out.df %>% 
        community_stability(time.var = "time",
                            abundance.var = "abundance",
                            replicate.var = "site")
    }
    
    # Creat path for sim output
    # sim.path <- file.path("figures", "sim_output", 
    #                       paste0("dispersal", mean(d), "-dorm",mean(dorm),"-act",mean(activ),"-dist",mean(disturb)))
    # if(!dir.exists(sim.path)) dir.create(sim.path, recursive = T)
    # extract SxS matrix
    # comm <- t(out.N[,,tsteps])
    # out.N
    
    # saveRDS(comm, file = file.path(sim.path, "model-output.rds"))
    # comms.out[[i]] <- comm
    # comm[is.na(comm)] <- 0
    # comm <- decostand(comm, method = "hellinger")
    # comm
    
    #plot patch 6 just to show local dynamics example
    autoplot.zoo(t(out.N[,10,]), facet = NULL) +
      labs(x = "Time", y = "Abundance") +
      theme_bw()
    #   ggsave(file.path(sim.path, "local_dynamics.png"), width = 8, height = 6, units = "in", dpi = 500)
    # 
    # specnumber(comm)
    # specnumber(colSums(comm))
    # ord <- rda(comm ~ E); plot(ord)
    #percent.env <- as.numeric(eigenvals(ord)[1]/sum(eigenvals(ord)))
    (alpha <- vegetarian::d(comm, lev = "alpha", q = 1))
    (beta <- vegetarian::d(comm, lev = "beta", q = 1))
    (gamma <- vegetarian::d(comm, lev = "gamma", q = 1))
    
    # write out dispersal, dormancy, and diversity
    out.sum[i,] <- c(d, dorm, alpha, beta, gamma)
    i <- i + 1
  }
}


data_frame(site = stability.dorm$site, 
           dorm = stability.dorm$stability, 
           no_dorm = stability.nodorm$stability) %>% 
  ggplot(aes(x = dorm, y = no_dorm, color = site)) + 
  geom_point() + 
  geom_abline(slope = 1, intercept = 0) + 
  coord_fixed()

out.sum
colnames(out.sum) <- c("Dispersal", "Dormancy", "Alpha", "Beta", "Gamma")
as.data.frame(out.sum) %>% 
  gather(Alpha, Beta, Gamma, key = Scale, value = Diversity) %>%
  mutate(Dormancy = factor(Dormancy, levels = dorm.grad, ordered = T)) %>% 
  ggplot(aes(Dispersal, Diversity, color = Dormancy)) +
  # geom_point(size = 2, alpha = 0.5) +
  geom_line(alpha = 0.8) + 
  facet_grid(Scale ~ ., scales = "free_y") +
  theme_minimal() + 
  scale_x_continuous(limits = c(0,1)) +
  scale_color_viridis(discrete = T, begin = .2, end = .8, direction = -1) + 
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "grey90"),
        legend.text = element_text(size = 8), legend.key.size = unit(0.3, "cm")) +
  ggsave(file.path("figures", "diversity-dispersal", paste0(ifelse(ddcov, "pos", "neg"),"-ddcov"),
                   paste0("diversity-dispersal_",env.type,"_dist",disturb,"_period",env.period,".png")), 
         width = 4, height = 4, units = "in", dpi = 500)


as.data.frame(out.sum) %>%
  gather(Alpha, Beta, Gamma, key = Scale, value = Diversity) %>%
  mutate(Dormancy = factor(Dormancy, levels = dorm.grad, ordered = T)) %>%
  filter(Scale == "Beta") %>%
  ggplot(aes(Dispersal, Diversity, color = Dormancy)) +
  # geom_point(size = 2, alpha = 0.5) +
  geom_line(alpha = 0.6) +
  theme_minimal() +
  scale_x_continuous(limits = c(0,0.2)) +
  scale_y_log10() +
  scale_color_viridis(discrete = T, begin = .2, end = .8, direction = -1) +
  theme(panel.grid = element_blank(), panel.background = element_rect(color = "grey90"))
