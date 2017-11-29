require(progress)
require(vegan)
require(vegetarian)
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
  alpha <- d(M, lev = "alpha")
  beta <- d(M, lev = "beta")
  gamma <- d(M, lev = "gamma")
  return(c(alpha, beta, gamma))
}

run.sim <- function(a){
  
  # set up metacom
  N <- matrix(n.ij, ncol = S, nrow = M)
  R <- rep(r.j, length = M) #Initial resources in each community
  
  out.dynamics <- list(NA)
  save.nums <- floor(seq(burnin/dt, length(timesteps), length.out = stepnum))
  for(t.step in 1:length(timesteps)){
    # Init progress bar
    if(t.step == 1) pb <- progress_bar$new(total = stepnum, force = T)
    
    pb$update(ratio = t.step/length(timesteps))
    e.vec <- E.j(t.step * dt)
    consumption <- sapply(env.optima, c.ij, E.j = e.vec)
    
    # Step through dynamics
    Nt1 <- N + dt*N*(e.ij * consumption * R - m.ij - a) + dt*migrate(N)
    Rt1 <- R + dt*I.j - dt*R*(l.j + rowSums(consumption*N)) 
    
    # Update densities
    N <- Nt1 * (Nt1 > ext)
    R <- Rt1 * (Rt1 > ext)
    
    if(t.step %in% save.nums){
      out.dynamics[[which(save.nums == t.step)]] <- N
    }
  }
  return(out.dynamics)
}

### Set Parameters of Model
stepnum <- 10000
dt <- 0.08
timesteps <- seq(0,stepnum,dt)
burnin <- 2000

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

env.period <- 50000 # period of env sinusoidal fluctuations
env.ampl <- 1 # amplitude of envrionment sinusoidal fluctuations
env.optima <- 1 - seq(0, env.ampl, by=env.ampl/(S-1))

## Run Simulation
mc.dynamics <- run.sim(a = a)
head(mc.dynamics)
tail(mc.dynamics)
sp.dyn <- sapply(mc.dynamics, colSums)
sp.dyn <- (t(sp.dyn))

png(filename = "figures/mc-dynamics.png", height = 6, width = 6, units = "in", res = 300)
matplot(sp.dyn, type = 'l', yaxt = "n", xaxt = "n", ylab = "", xlab = "", col = viridis::viridis(20))
axis(side = 1, labels = T, lwd.ticks = 2)
axis(side = 2, labels = T, las = 1, lwd.ticks = 2)
box(lwd = 2)
mtext(side = 1, line = 3, text = "Time", cex = 1.5)
mtext(side = 2, line = 3, text = "Abundance", cex = 1.5)
dev.off()

div.part(mc.dynamics[[10000]])

div.traj <- sapply(mc.dynamics, FUN = div.part)
div.traj <- t(div.traj)
colnames(div.traj) <- c("alpha", "beta", "gamma")
matplot(div.traj, type = 'l', col = c("red","green","blue"), lwd = 2)
div.traj[8000,]

disp.vec <- c(0.00001, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.75, 1)
disp.dynamics <- matrix(NA, nrow = length(disp.vec), ncol = 4)
i <- 1
for(disp in disp.vec){
  mc.dynamics <- run.sim(a = disp)
  mc.div <- div.part(mc.dynamics[[10000]])
  disp.dynamics[i,] <- c(disp, mc.div)
}
