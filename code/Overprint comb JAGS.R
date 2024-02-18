# JAGS model file
model {
  
  # evaluate R1.mix against LA-ICP-MS En10 data
  for (i in 1:n.En10){
    
    R.En10[i] ~ dnorm(R1.En10[En10.bins[i]], 1/R.sd.En10[i]^2)
    
  }
  
  # evaluate R1.mix against LA-ICP-MS En9 data
  for (i in 1:n.En9){
    
    R.En9[i] ~ dnorm(R1.En9[En9.bins[i]], 1/R.sd.En9[i]^2)
    
  }
  
  # evaluate R1.mix against micromill data
  for (i in 1:n.Rm3.5b){
    
    R.Rm3.5b[i] ~ dnorm(R1.Rm3.5b[Rm3.5b.bins[i]], 1/R.sd.Rm3.5b[i]^2)
    
  }
  
  # evaluate R1.mix against hand drill data
  for (i in 1:n.drill){
    
    R.drill[i] ~ dnorm(R1.drill[drill.bins[i]], 1/R.sd.drill[i]^2)
    
  }
  
  # evaluate R1.m against LA-ICP-MS En1 data
  for (i in 1:n.r2){
    
    R.r2[i] ~ dnorm(R1.m[r2.bins[i]], 1/R.sd.r2[i]^2)
    
  }
  
  # evaluate R1.m against dentine data
  for (i in 1:n.r1){
    
    R.r1[i] ~ dnorm(R1.m[r1.bins[i]], 1/R.sd.r1[i]^2)
    
  }

  # establish evaluation bins
  # measurement record
  drill.bins = trunc(drill.tl/bin.thin + switch)
  Rm3.5b.bins = trunc(Rm3.5b.tl/bin.thin + switch)
  En10.bins = trunc(En10.tl/bin.thin + switch)
  En9.bins = trunc(En9.tl/bin.thin + switch)
  
  # En1 record
  r2.bins = trunc(r2.tl/bin.thin + switch)
  
  # dentine record
  r1.bins = trunc(r1.tl/bin.thin + switch)
  
  # mass balance equation 
  R1.drill = R1.m * (1 - pr.drill) + Raft.mod * pr.drill
  R1.Rm3.5b = R1.m * (1 - pr.Rm3.5b) + Raft.mod * pr.Rm3.5b
  R1.En9 = R1.m * (1 - pr.En9) + Raft.mod * pr.En9
  R1.En10 = R1.m * (1 - pr.En10) + Raft.mod * pr.En10
  
  # proportion of post switch over print
  pr.drill ~ dbeta(1, 1) #uninformative prior
  pr.Rm3.5b ~ dbeta(1, 1) #uninformative prior
  pr.En10 ~ dbeta(1, 1) #uninformative prior
  pr.En9 ~ dbeta(1, 1) #uninformative prior
  
  for (i in 2:t){
    #serum ratio
    R1.m[i] <- R1.m[i - 1] + b * (R2.m[i - 1] - R1.m[i - 1]) + a * (Rin[i - 1] - R1.m[i - 1])
    
    #bone ratios
    R2.m[i] <- R2.m[i - 1] + c * (R1.m[i - 1] - R2.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  R2.m[1] ~ dnorm(Rin[1], Sr.pre.2) 
  R1.m[1] ~ dnorm(Rin[1], Sr.pre.1)
  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.2 ~ dgamma(Sr.pre.shape, Sr.pre.rate.2)  
  Sr.pre.rate.2 <- 5e-7
  
  Sr.pre.1 ~ dgamma(Sr.pre.shape, Sr.pre.rate.1) 
  
  Sr.pre.shape <- 100
  Sr.pre.rate.1 <- 5e-6

  
  for(i in 1:t){
    
    Rin[i] ~ dnorm(ifelse(i > switch, Raft.mod, Rpri.mod), ifelse(i > switch, Raft.pre, Rpri.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Raft.mean, Rin.m.pre)
  }

  #suspected date of the switch with some uncertainty
  switch ~ dnorm(400/bin.thin, 1 )
  
  # priors
  Rpri.mod ~ dnorm(Rpri, Rpri.pre)
  
  Raft.mod ~ dnorm(Raft, Raft.pre)
  
  Raft.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Raft)
  Rpri.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Rpri)
  Sr.pre.rate.Raft <- 2e-5
  Sr.pre.rate.Rpri <- 2e-6
  
  Rpri = CA.Sr  # California 87Sr/86Sr
  
  Raft = UT.Sr  # Utah 87Sr/86Sr
  
  #sampling from parameter posteriors from the calibration
  # scale these parameters with bin.thinning parameter
  a = exp(params[1]) * bin.thin
  
  b = exp(params[2]) * bin.thin
  
  c = exp(params[3]) * bin.thin
  
  #supply the model with estimated parameters
  params ~ dmnorm.vcov(params.mu, params.vcov)
  
}