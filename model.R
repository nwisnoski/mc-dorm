


sim.metacom <- function(species=3, dispersal=0.01, patches=3){
  
  N <- matrix(10, ncol = species, nrow = patches)
  R <- rep(10*species, patches) #Initial resources in each community
  D <- matrix(0, ncol = species, nrow = patches)
  
  Rin <- 100 # resource in
  Rout <- 10 # resource out
  eff <- 0.2 # conversion efficiency
  mort <- 0.2 # mortality
  Ext <- 0.1 # extinction Threshold
  d <- 0.01 # max dormancy rate
  s <- 2 # dormancy tracking; 0 = constant enter, higher = more specific
  r <- 0.01 # max resuscitation rate
  q <- 0 # How responsive reactivation is; 0 = random, higher = responsive
  u <- 100 # how much longer dormant propagules live than active
  env_period <- 40000 # period of env sinusoidal fluctuations
  env_ampl <- 1 # amplitude of envrionment sinusoidal fluctuations
  
  Tmax <- 140000 # number of time steps in Sim
  DT <- 0.08 # % size of discrete "time steps"
  
  
  # species environmental optima; partition 0:max amplitude by # of species
  env_optimum <- 1 - seq(0, env_ampl, by=env_ampl/(species-1))
  
  # dispersal conditions (global dispersal)
  dispersal_matrix <- matrix(1/(patches-1), nrow = patches, ncol = patches)
  diag(dispersal_matrix) <- 0
  calc.immigration <- function(N,a,dispersal_matrix){
    dispersal_matrix %*% N * rep(a,each=patches)
  }
  disp <- rep(dispersal,species)
  dispM <- matrix(rep(disp,each=patches),patches,species)
  
  Prod <- matrix(NA,species*patches,40000)
  Abund <- Prod
  Bank <- Prod
  
  for(TS in 1:Tmax){
    Immigrants <- calc.immigration(N,disp,dispersal_matrix)
    envt.v <- 0.5*env_ampl*(
      sin((2*pi/env_period)*TS + 1 + (1:patches)*2*pi/patches) + 1) #calculates current environment
    
    # Calculate consumption. sapply() function is subtracting eopt from envt
    consume <- 0.1 * (1.5-abs(sapply(env_optimum,'-',envt.v))) #calculates current comsumption
    dorm <- d * exp(-s*abs(sapply(env_optimum,'-',envt.v)))
    resus <- r * exp(-q*abs(sapply(env_optimum,'-',envt.v)))
    
    
    #abundance step
    Nt <- N*(1 + DT*(eff*R*consume - dispM - mort - dorm)) + 
      DT*Immigrants + D*DT*(resus)
    Dt <- D*(1 - DT*resus - mort/u) + N*DT*dorm
    
    #resource step
    Rt <- DT*Rin + R*(1 - DT*(Rout + rowSums(consume*N)))   
    
    N <- Nt * (Nt>Ext) # set to 0 if below extinction threshold
    D <- Dt * (Dt>Ext)
    R <- Rt
    
    if(TS>=100000){ #samples data after time step 100 000
      Prod[,(TS-100000)] <- c(t(eff*consume*R*N)) #Productivity
      Abund[,(TS-100000)] <- c(t(N)) #Abundance
      Bank[,(TS-100000)] <- c(t(D)) # Seed bank
    }
  } 
  
  Prod <- array(t(Prod),dim=c(40000,species,patches))
  Prod <- Prod[seq(1,40000,100),,] # take ever 100th value
  Abund <- array(t(Abund),dim=c(40000,species,patches))
  Abund <- Abund[seq(1,40000,100),,] # take every 100th value
  Bank <- array(t(Bank),dim=c(40000,species,patches))
  Bank <- Bank[seq(1,40000,100),,] # take every 100th value
  
  return(list(Prod=Prod, Abund=Abund, Bank=Bank))
}

#Metacommunity multifunctionality simulation####
runs <- 5 #number of replicates
species <- 2 #number of species
patches <- 3 #number of patches

#DispV <- c(0.0001,0.0005,0.001,0.0015,0.005,0.01,0.05,0.1,0.5,1) #dispersal rates

#calculate abundances
#SIH_data<-sapply(DispV,SIH,species=species,patches=patches)

mc_data <- sim.metacom(species = 3, dispersal = .01, patches = 3)
matplot(mc_data$Abund[,3,], type = "l", lwd = 3)
matplot(mc_data$Bank[,3,], type = "l", lwd = 3)
