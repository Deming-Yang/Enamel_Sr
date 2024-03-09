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

n.mea = length(tusk.mill.tl$Sr)

# remove the last 2 data point to avoid 0 dist

R.sd.mea <- tusk.mill.tl$sd[1:(n.mea-2)]
dist.mea <- misha.tusk.micromill.dist[1:(n.mea-2)]/1e3
R.mea <- tusk.mill.tl$Sr[1:(n.mea-2)]

n.mea = length(R.mea)

s.intv <- 0.547 # mm micromill interval

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
# for the tusk, the growth was considered constant at 14.7 microns/day
length.v <- seq(0, 25, by = 1)
rate.v <- rep(14.7e-3, length(length.v))

rate.sd <- 0.15 * 14.7e-3 # mm per day

# make sure that bin.thin is larger than 2
bin.thin <- 0.25 * s.intv/mean(rate.v)

# use timeline reconstruction results to estimate t

t = round(1.2 * (max(tusk.mill.tl$tl) - min(tusk.mill.tl$tl))/bin.thin)

#posterior samples of parameters from Misha's calibration


parameters <- c("rate", "dist","Rin.m.pre",
                "R1.m","Rin.m", "R2.m","a","b", "c")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            s.intv = s.intv, max.dist.mea = max.dist.mea, bin.thin = bin.thin,
            length.v = length.v, rate.v = rate.v, rate.sd = rate.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 5e3
n.burnin = 2e3
n.thin = 1

#Run it
post.misha.M640b = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param bt.R", 
                                                     parameters.to.save = parameters, 
                                                     data = dat, n.chains=5, n.iter = n.iter, 
                                                     n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 51 hours

save(post.misha.M640b, file = "out/post.misha.M640b.RData")

load("out/post.misha.M640b.RData")

post.misha.M640b$BUGSoutput$summary

#check posterior density of parameter a:

plot(density(post.misha.M640b$BUGSoutput$sims.list$a),col="red")
