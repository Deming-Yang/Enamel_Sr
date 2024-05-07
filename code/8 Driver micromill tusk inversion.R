library(coda)
library(rjags)
library(R2jags)
library(mcmcplots)
library(bayestestR)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)

source("code/1 Helper functions.R")

############################# inversion of Misha's tusk micromill data using parameters #######

# have to adjust the sd, If used as is,
# the model won't converge on the a reasonable set of reaction rate parameters
# It seems that 20 time the uncertainty would converge on reasonable reaction rates
R.sd.mea <- tusk.mill.tl$sd * 10

# add 10mm distance to all data to avoid 0 dist towards the end of the sequence
dist.mea <- misha.tusk.micromill.dist/1e3 +10

R.mea <- tusk.mill.tl$Sr

n.mea = length(R.mea)

# micromill interval
s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading history

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
# for the tusk, the growth was considered constant at 14.7 microns/day
length.v <- seq(0.1, 50, by = 0.1)
rate.v <- rep(14.7e-3, length(length.v))

bin.thin <- 0.2 * s.intv/min(rate.v)

M640.bt <- bin.thin

d.offset.M640 <- 2* s.intv/mean(rate.v)

# rate.sd <- 0.1 * mean(rate.v)  # mm per day

t = round(1.35 * max(M640.micromill.tl)/bin.thin)

t.M640b <-t


parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, bin.thin = bin.thin,
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e5
n.burnin = 2e4
n.thin = 30

#Run it
post.misha.M640b = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 40 min

save(post.misha.M640b, file = "out/post.misha.M640b.RData")

load("out/post.misha.M640b.RData")

denplot(as.mcmc(post.misha.M640b), parms = c("a","b","c"))


post.misha.M640b.R1.m.89 <- MCMC.CI.bound(post.misha.M640b$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.M640b.Rin.m.89 <- MCMC.CI.bound(post.misha.M640b$BUGSoutput$sims.list$Rin.m, 0.89)


# calculate RXN equivalent parameters

#flux ratio
M640b.flux.ratio <- post.misha.M640b$BUGSoutput$sims.list$a/post.misha.M640b$BUGSoutput$sims.list$b

M640b.f1 <- M640b.flux.ratio/(1+M640b.flux.ratio)

M640b.f2 <- 1- M640b.f1

#pool ratio
M640b.lambda1 <- post.misha.M640b$BUGSoutput$sims.list$a/M640.bt /M640b.f1

M640b.lambda2 <- post.misha.M640b$BUGSoutput$sims.list$c/M640.bt * M640b.f1

# Extract MAPE estimates
map_estimate(M640b.f1)[[1]]
# 0.7126433

M640b.f1.ci1 <- hdi(M640b.f1, .89)[[2]]
# 0.5911979

M640b.f1.ci2 <- hdi(M640b.f1, .89)[[3]]
# 0.7837885

map_estimate(M640b.f2)[[1]]
# 0.2873567

M640b.f2.ci1 <- hdi(M640b.f2, .89)[[2]]
# 0.2162115

M640b.f2.ci2 <- hdi(M640b.f2, .89)[[3]]
# 0.4088021

map_estimate(M640b.lambda1)[[1]]
# 0.03842481

M640b.lambda1.ci1 <- hdi(M640b.lambda1, .89)[[2]]
# 0.02972423

M640b.lambda1.ci2 <- hdi(M640b.lambda1, .89)[[3]]
# 0.04848269

map_estimate(M640b.lambda2)[[1]]
# 0.002703624

M640b.lambda2.ci1 <- hdi(M640b.lambda2, .89)[[2]]
# 0.001583127

M640b.lambda2.ci2 <- hdi(M640b.lambda2, .89)[[3]]
# 0.004168115

############## Model sensitivity demonstrations ##############

# here if the original range of uncertainty is used, the model will overfit the data
R.sd.mea <- tusk.mill.tl$sd 

# add 10mm distance to all data to avoid 0 dist towards the end of the sequence
dist.mea <- misha.tusk.micromill.dist/1e3 +10

R.mea <- tusk.mill.tl$Sr

n.mea = length(R.mea)

# micromill interval
s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading history

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
# for the tusk, the growth was considered constant at 14.7 microns/day
length.v <- seq(0.1, 50, by = 0.1)
rate.v <- rep(14.7e-3, length(length.v))

bin.thin <- 0.2 * s.intv/min(rate.v)

M640.bt <- bin.thin

d.offset.M640 <- 2* s.intv/mean(rate.v)

# rate.sd <- 0.1 * mean(rate.v)  # mm per day

t = round(1.35 * max(M640.micromill.tl)/bin.thin)

t.M640b <-t


parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, bin.thin = bin.thin,
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e5
n.burnin = 2e4
n.thin = 30

#Run it
post.misha.M640b.sens = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                              parameters.to.save = parameters, 
                                              data = dat, n.chains=5, n.iter = n.iter, 
                                              n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 40 min

save(post.misha.M640b.sens, file = "out/post.misha.M640b.sens.RData")

load("out/post.misha.M640b.sens.RData")

denplot(as.mcmc(post.misha.M640b.sens), parms = c("a","b","c"))


post.misha.M640b.sens.R1.m.89 <- MCMC.CI.bound(post.misha.M640b.sens$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.M640b.sens.Rin.m.89 <- MCMC.CI.bound(post.misha.M640b.sens$BUGSoutput$sims.list$Rin.m, 0.89)

# check parameter a
# plot(density(exp(log.a)), col = "blue", ylim = c(0,110),
#      main="Posteriors of parameter a")
# lines(density(post.misha.M640b.sens$BUGSoutput$sims.list$a/M640.bt),col="red")
# legend(0.03,110,c("LA-ICP-MS tusk","Micromill tusk"),
#        lwd = c(1,1), col=c("blue","red"))