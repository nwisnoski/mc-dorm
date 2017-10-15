require(progress)
require(vegan)
rm(list=ls())

c.ij <- function(H.i, E.j){
  return(1.5 - abs(H.i - E.j))
}

E.j <- function(t){
  0.5*env.ampl*(
    sin((2*pi/env.period)*t + 1 + (1:M)*2*pi/M) + 1)
}

migrate <- function(N){
  (colSums(N) - N) * a / (M - 1)
}

div.part <- function(M){
  M.pa <- decostand(M, method = 'pa')
  alpha <- rowSums(M.pa)
  gamma <- sum(colSums(M.pa) >= 1)
  beta <- gamma / alpha
  return(list(alpha = alpha, beta = beta, gamma = gamma))
}

run.sim <- function(a){
  
  # set up metacom
  N <- matrix(n.ij, ncol = S, nrow = M)
  R <- rep(r.j, length = M) #Initial resources in each community
  
  out.dynamics <- list(NA)
  out.dynamics[[1]] <- N
  for(t.step in 1:length(timesteps)){
    # Init progress bar
    if(t.step == 1) pb <- progress_bar$new(total = stepnum, force = T)
    
    pb$update(ratio = t.step/length(timesteps))
    e.vec <- E.j(t.step)
    consumption <- sapply(env.optima, c.ij, E.j = e.vec)
    
    # Step through dynamics
    Nt1 <- N*(1 + dt*(R * e.ij * consumption - m.ij - a)) + dt*migrate(N)
    Rt1 <- R*(1 - dt*(rowSums(consumption*N) + l.j)) + dt*I.j
    
    # Update densities
    N <- Nt1 * (Nt1 > ext)
    R <- Rt1 * (Rt1 > ext)
    
    if(t.step > 10000){
      out.dynamics[[t.step-10000]] <- N
    }
  }
  return(out.dynamics)
}

### Set Parameters of Model
stepnum <- 10000
dt <- 0.1
timesteps <- seq(1,stepnum,dt)

M <- 10 # number of patches
S <- 20 # number of species
n.ij <- 10 # number of individuals in each patch at t0
r.j <- 30 # number of resources at t0

I.j <- 10 # resource inputs
l.j <- 3 # resource loss

e.ij <- 0.2 # conversion efficiencies
m.ij <- 0.3 # mortality rates

a <- .1 # dispersal rate
ext <- 0.1 # extinction threshold

env.period <- 50000 # period of env sinusoidal fluctuations
env.ampl <- 1 # amplitude of envrionment sinusoidal fluctuations
env.optima <- 1 - seq(0, env.ampl, by=env.ampl/(S-1))

## Run Simulation
mc.dynamics <- run.sim(a = a)

sp.dyn <- sapply(mc.dynamics, colSums)
sp.dyn <- (t(sp.dyn))
#head(sp.dyn)
matplot(sp.dyn, type = 'l')
div.part(mc.dynamics[[10000]])

