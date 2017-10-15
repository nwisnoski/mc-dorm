require(progress)
require(vegan)
rm(list=ls())

c.ij <- function(H.i, E.j){
  return((1.5 - abs(H.i - E.j))/10)
}

E.j <- function(t){
  0.5*env.ampl*(
    sin((1:M)*2*pi/M + (2*pi/env.period)*t) + 1)
}

migrate <- function(N){
  (colSums(N) - N) * a / (M - 1)
}

div.part <- function(M){
  M.pa <- decostand(M, method = 'pa')
  alpha <- rowSums(M.pa)
  gamma <- sum(colSums(M.pa) >= 1)
  beta <- gamma - alpha
  return(c(mean(alpha), mean(beta), gamma))
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

M <- 30 # number of patches
S <- 9 # number of species
n.ij <- 10 # number of individuals in each patch at t0
r.j <- 9 # number of resources at t0

I.j <- 150 # resource inputs
l.j <- 10 # resource loss

e.ij <- 0.2 # conversion efficiencies
m.ij <- 0.2 # mortality rates

a <- .01 # dispersal rate
ext <- 0.1 # extinction threshold

env.period <- 40000 # period of env sinusoidal fluctuations
env.ampl <- 1 # amplitude of envrionment sinusoidal fluctuations
env.optima <- 1 - seq(0, env.ampl, by=env.ampl/(S-1))

## Run Simulation
mc.dynamics <- run.sim(a = a)
#mc.dynamics <- mc.dynamics[seq(1, length(timesteps), length.out = stepnum)]

sp.dyn <- sapply(mc.dynamics, colSums)
sp.dyn <- (t(sp.dyn))
head(sp.dyn)
matplot(sp.dyn, type = 'l')
div.traj <- sapply(mc.dynamics, FUN = div.part)
div.traj <- t(div.traj)
colnames(div.traj) <- c("alpha", "beta", "gamma")
matplot(div.traj, type = 'l', col = c("red","green","blue"), lwd = 2)
div.traj[80000,]
