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

############################# inversion of Misha's molar LA-ICP-MS data using parameters #######

Enamel1.Sr <- x.pt.avg(proc.Enamel1.rm.f$avg, 0, 100)

Enamel1.length <- x.pt.avg(proc.Enamel1.rm.new.x, 0, 100)

n.mea = length(Enamel1.length$avg.sr)

dist.mea <- Enamel1.length$avg.sr
R.mea <- Enamel1.Sr$avg.sr
n.mea = length(R.mea)
R.sd.mea <- Enamel1.Sr$sd.sr

s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
length.v <- rev(ref.length.v)
rate.v <- rev(Ref.growth.v)

rate.sd <- 0.005 #mm per day

bin.thin <- 0.15 * s.intv/mean(rate.v)

# use timeline reconstruction results to estimate t

# t = round(1.2 * max(days.cumm.en1)/bin.thin)

t = round(1.2 * max(days.cumm.en1))

t.En1 <-t

#posterior samples of parameters from Misha's calibration


parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, 
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.En1 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param molar.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 15 hours

save(post.misha.En1, file = "out/post.misha.En1.RData")

load("out/post.misha.En1.RData")

post.misha.En1$BUGSoutput$summary

#check posterior density of parameter a:

plot(density(post.misha.En1$BUGSoutput$sims.list$a),col="red")

post.misha.En1.R1.m.89 <- MCMC.CI.bound(post.misha.En1$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En1.Rin.m.89 <- MCMC.CI.bound(post.misha.En1$BUGSoutput$sims.list$Rin.m, 0.89)

plot(en1.tl$tl + 410, en1.tl$Sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(0,1500),ylim=c(0.705,0.714),
     xlab="Days from Misha's move",
     main = "Tusk micromill vs molar E micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.Rin.m.89[[1]], lwd = 2)
lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.Rin.m.89[[2]], lty = 2)
lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.Rin.m.89[[3]], lty = 2)

####### test with bin thinning#############
Enamel1.Sr <- x.pt.avg(proc.Enamel1.rm.f$avg, 0, 100)

Enamel1.length <- x.pt.avg(proc.Enamel1.rm.new.x, 0, 100)

n.mea = length(Enamel1.length$avg.sr)

dist.mea <- Enamel1.length$avg.sr
R.mea <- Enamel1.Sr$avg.sr
n.mea = length(R.mea)
R.sd.mea <- Enamel1.Sr$sd.sr

s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
length.v <- rev(ref.length.v)
rate.v <- rev(Ref.growth.v)

rate.sd <- 0.005 #mm per day

bin.thin <- 0.2 * s.intv/mean(rate.v)

t = round(1.2 * max(days.cumm.en1)/bin.thin)

t.En1.bt <-t

parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, bin.thin = bin.thin,
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 4e3
n.thin = 6

#Run it
post.misha.En1.bt = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param bt.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 51 hours

save(post.misha.En1.bt, file = "out/post.misha.En1.bt.RData")

load("out/post.misha.En1.bt.RData")

post.misha.En1.bt$BUGSoutput$summary

#check posterior density of parameter a:

plot(density(post.misha.En1.bt$BUGSoutput$sims.list$a/bin.thin),col="red")

post.misha.En1.bt.R1.m.89 <- MCMC.CI.bound(post.misha.En1.bt$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En1.bt.Rin.m.89 <- MCMC.CI.bound(post.misha.En1.bt$BUGSoutput$sims.list$Rin.m, 0.89)

plot(en1.tl$tl + 410, en1.tl$Sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(0,1500),ylim=c(0.705,0.714),
     xlab="Days from Misha's move",
     main = "Tusk micromill vs molar E micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.bt.Rin.m.89[[1]], lwd = 2)
lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.bt.Rin.m.89[[2]], lty = 2)
lines(1:t -(2* s.intv/mean(rate.v)), post.misha.En1.bt.Rin.m.89[[3]], lty = 2)