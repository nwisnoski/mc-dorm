env.dynamics = matrix(NA, nrow = 0, ncol = length(env[,1]))
for (t in 1:tsteps){
  env.dynamics <- rbind(env.dynamics, asynchronize.env(env = env[,1], t, period = 1000, moran = F))
}
matplot(env.dynamics, type = 'l')
