model {
  ######inversion######
  for (i in 1:n.mea){
    #averaging data
    R1.eva[i] = inprod(R1.m, ((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2) < dist))/sum(((dist.mea[i] + s.intv/2) < dist) - ((dist.mea[i] - s.intv/2)  < dist))
    #evaluate its ratio
    R.mea[i] ~ dnorm(R1.eva[i], 1/(R.sd.mea[i])^2)
  }
  
  #Data model priors for molar growth
  #grow by the day
  # establish a time-length relationship
  for (i in 2:t){
    
    dist[i] = dist[i - 1] - rate[i - 1] * bin.thin #simulate distance increment ber time step
    
    rate[i] = interp.lin(dist[i], length.v, rate.v)
    
  }
  
  rate[1] = interp.lin(dist[1], length.v, rate.v)

  dist[1] = max.dist.mea #maximum distance from the cervix in mm

  for (i in 2:t){
    #serum ratio
    R1.m[i] = R1.m[i - 1] + b * (R2.m[i - 1] - R1.m[i - 1]) + a * (Rin.m[i - 1] - R1.m[i - 1])
    
    #bone ratios
    R2.m[i] = R2.m[i - 1] + c * (R1.m[i - 1] - R2.m[i - 1])
  }
  
  # assuming the starting value of the two pools are close to the starting Rin value
  #e.g., close to equilibrium with Rin
  R2.m[1] ~ dnorm(Rin.m[1], Sr.pre.2) #use a different error term for pool 2
  R1.m[1] ~ dnorm(Rin.m[1], Sr.pre.1)
  
  #precision for average Sr measurements in bone should be much smaller
  Sr.pre.2 ~ dgamma(Sr.pre.shape, Sr.pre.rate.2)  
  Sr.pre.rate.2 = 5e-7
  
  Sr.pre.1 ~ dgamma(Sr.pre.shape, Sr.pre.rate.1) 
  
  Sr.pre.shape = 100
  Sr.pre.rate.1 = 5e-6
  
  #generate null input series
  for (i in 2:t){
    
    Rin.m[i] = Rin.m[i - 1] + Rin.m.cps[i]
    
    Rin.m.cps[i] ~ dt(0, Rin.m.pre, 1) T(-6e-3, 6e-3) #Brownian motion, cauchy error term

  }
  # initiate the series with an reasonable prior
  Rin.m[1] ~ dnorm(Rin.int, Rin.m.pre) #allowed some variation
  
  Rin.int ~ dnorm(0.706, 1/0.01^2)  #a weakly informative initial value
  
  #initial change per step
  Rin.m.cps[1] ~ dt(0, Rin.m.pre, 1) T(-6e-3, 6e-3)
  
  Rin.m.pre ~ dgamma(Rin.m.pre.shp, Rin.m.pre.rate)
  Rin.m.pre.shp = 100
  Rin.m.pre.rate = 2.5e-8 * bin.thin
  
  #sampling from parameter posteriors from the calibration
  a = exp(params[1]) * bin.thin
  
  b = exp(params[2]) * bin.thin
  
  c = exp(params[3]) * bin.thin
  
  #supply the model with estimated parameters
  params ~ dmnorm.vcov(params.mu, params.vcov)
  
}