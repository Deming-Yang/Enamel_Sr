# JAGS model file
model {
  
  # evaluate R1.mix against enamel data
  
  for (i in 1:n.mea){
    
    R.mea[i] ~ dnorm(R1.eva[i], 1/(R.sd.mea[i])^2)
  }
  
  # generate indeces for data referencing
  scwax.bins = t - trunc(age.sc.d/age.res)
  
  
  # evaluate R1.m against dentine data
  
  for (i in 1:n.sc.wax){
    #short chain wax
    scwax.d2H.dat[i] ~ dnorm(scwax.d2H[scwax.bins[i]], 1/scwax.d2H.sd[i]^2)
    
  }
  

  
  #Data model priors for tissue growth
  for (i in 2:t){
    Ivo.rate[i] ~ dnorm(Ivo.rate.mean, 1/Ivo.rate.sd^2) #ivory growth rate, micron/day
    dist[i] <- dist[i - 1] - Ivo.rate[i] #simulate daily distance increment
  }
  #
  
  # mass balance equation 
  
  R1.mix <- R1.m * (1 - r) + UT.Sr * r
  
  # proportion of post switch over print
  r ~ dbeta(1, 1) #uninformative prior
  
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
    
    Rin[i] ~ dnorm(ifelse(i > switch, Raft.mean, Rpri.mean), ifelse(i > switch, Raft.pre, Rpri.pre))
    #When i is greater than the switch point, Rin ~ dnorm(Raft.mean, Rin.m.pre)
  }
  
  # model movement in Misha's intake
  switch ~ dcat(pi)
  
  pi <- c(pi.int, pi.switch, pi.end)
  
  pi.end <- rep(0, t - date - err.date)
  
  pi.switch <- rep(1, 1 + 2 * err.date)
  
  pi.int <- rep(0, date - err.date - 1)
  
  #uncertainty of the date of switch = +- the number of days
  err.date <- 2 
  
  #suspected date of the switch
  date <- 85
  
  
  # priors
  Rpri.mean ~ dnorm(Rpri, Rpri.pre)
  
  Raft.mean ~ dnorm(Raft, Raft.pre)
  
  Raft.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Raft)
  Rpri.pre ~ dgamma(Sr.pre.shape, Sr.pre.rate.Rpri)
  Sr.pre.rate.Raft <- 2e-5
  Sr.pre.rate.Rpri <- 2e-6
  
  Rpri ~ dnorm(CA.Sr, 1 / CA.Sr.sd^2)  # California 87Sr/86Sr
  
  
  Raft ~ dnorm(UT.Sr, 1 / UT.Sr.sd^2)  # Utah 87Sr/86Sr


  #sampling from parameter posteriors from the calibration
  a <- a.post[indx]
  b <- b.post[indx]
  c <- c.post[indx]
  
  indx ~ dcat(rep(1, post.leng))


}