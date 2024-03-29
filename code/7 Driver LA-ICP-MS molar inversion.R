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

####### molar LA-ICP-MS inversion with data reduction and bin thinning#############

####### data set 1, enamel transect 1#############

En1.50avg.tl <- x.pt.avg(proc.Enamel1.rm.f$avg, en1.tl$tl, 50)

Enamel1.Sr <- x.pt.avg(proc.Enamel1.rm.f$avg, 0, 50)

Enamel1.length <- x.pt.avg(proc.Enamel1.rm.new.x, 0, 50)

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

rate.sd <- 0.05 * mean(rate.v) # mm per day

bin.thin <- 0.25 * s.intv/min(rate.v)

d.offset.En1 <- 2* s.intv/mean(rate.v)

En1.bt <- bin.thin

t = round(1.2 * max(days.cumm.en1)/bin.thin)

t.En1 <-t

parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, bin.thin = bin.thin,
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e4
n.burnin = 2e4
n.thin = 30

#Run it
post.misha.En1 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                               parameters.to.save = parameters, 
                                               data = dat, n.chains=5, n.iter = n.iter, 
                                               n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 1.2 hours

save(post.misha.En1, file = "out/post.misha.En1.RData")

load("out/post.misha.En1.RData")

denplot(as.mcmc(post.misha.En1), parms = c("a","b","c"))

post.misha.En1.R1.m.89 <- MCMC.CI.bound(post.misha.En1$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En1.Rin.m.89 <- MCMC.CI.bound(post.misha.En1$BUGSoutput$sims.list$Rin.m, 0.89)

####### End of data set 1 #############

