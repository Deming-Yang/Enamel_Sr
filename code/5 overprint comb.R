library(readxl)
library(dplyr)
library(coda)
library(R2jags)
library(ggplot2)
library(mcmcplots)
library(bayesplot)
library(bayestestR)

source("./code/1 Helper functions.R")

# modeling the amount of maturation overprint in samples. 
# this was done on four data series: LA-ICP-MS Enamel 9, Enamel 10, Drill, and micromill series
# the tusk dentine micromill data (solution method) was used to generate a reference timeline

# adding CA and UT Sr posterior from Yang et al. 2023
CA.Sr <- 0.70597
UT.Sr <- 0.71115

# uncertainty CA & UT Sr 
CA.Sr.sd <- (0.706243 - CA.Sr)/2

UT.Sr.sd <- (0.7120174 - UT.Sr)/2

# uncertainty CA & UT Sr 
# CA.Sr.sd <- 0.0005
# 
# UT.Sr.sd <- 0.001

# data reduction for LA-ICP-MS

days.cumm.en9.al <- days.cumm.en9 - 205 
days.cumm.en10.al <- days.cumm.en10 - 195 

#visualize timeline adjustment
par(mfrow=c(1,1))
plot(days.cumm.en9.al, proc.Enamel9.rm.f$avg, col= alpha("lightcyan4", 0.9),
     type = "l", lwd=2, xlim=c(-400,500),ylim=c(0.704,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")

#CA and UT Sr measurements from Yang et al. 2023
abline(h = CA.Sr)

abline(h = UT.Sr)

lines(days.cumm.en10.al, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.9),lwd=2)

# data thinning, to reduce computational demand
En1.50avg.tl <- x.pt.avg(proc.Enamel1.rm.f$avg, days.cumm.en1.al, 50)

En9.50avg.tl <- x.pt.avg(proc.Enamel9.rm.f$avg, days.cumm.en9.al, 50)

En10.50avg.tl <- x.pt.avg(proc.Enamel10.rm.f$avg, days.cumm.en10.al, 50)

#get the 1/2 average time interval as undertainty for timeline matching
En1.tl.sd <- min(base::diff(En1.50avg.tl$avg.tl))/2

En9.tl.sd <- min(base::diff(En9.50avg.tl$avg.tl))/2
En10.tl.sd <- min(base::diff(En10.50avg.tl$avg.tl))/2


# identify the max and min number of days to model
# use the tusk micromill series as the reference
# trim off Sr history beyond that of the tusk series
max(tusk.mill.tl$tl) - min(tusk.mill.tl$tl) #>1378 days

# the drill series is longer than the tusk series
drill.tl.f <- filter(drill.tl,
                          drill.tl$tl > min(tusk.mill.tl$tl) &
                          drill.tl$tl < max(tusk.mill.tl$tl))

# the En1 series is also longer than the tusk series
En1.50avg.tl.f <- filter(En1.50avg.tl,
                     En1.50avg.tl$avg.tl > min(tusk.mill.tl$tl) &
                       En1.50avg.tl$avg.tl < max(tusk.mill.tl$tl))


drill.tl.sd <- min(base::diff(drill.tl.f$tl))/2

Rm3.5b.mill.tl.sd <- min(base::diff(Rm3.5b.mill.tl$tl))/2 

# parameter estimates from Misha's turnover model, generated from Yang et al., 2023
# all parameters have been log-transformed
turnover.params <- read.csv("./data/turnover_param_posterior.csv")

log.a <- turnover.params$log.a
log.b <- turnover.params$log.b
log.c <- turnover.params$log.c

# use a multi-normal distribution as a prior
# prior parameter estimates using means and a variance-covariance matrix
# means
turnover.params.mu <- c(mean(log.a), mean(log.b), mean(log.c))

turnover.params <- data.frame(log.a, log.b, log.c)

# calculate v.cov matrix
turnover.params.vcov <- var(turnover.params)

####### calculate post-move Sr overprint for drilled samples #######

####### start compiling parameters #######
# used to reduce computational demand

Rm3.5b.tl.sd <- Rm3.5b.mill.tl.sd #uncertainty in the timeline

bin.thin <- min(drill.tl.sd, Rm3.5b.tl.sd, En1.tl.sd, En9.tl.sd, En10.tl.sd)

t <- round(1600/bin.thin) # set a higher number of days to accommodate 

# use a universal uncertainty among the substrates
# use a large combined error to ensure parameter convergence on the tusk data,
# also estimate the common behavior between the data series, without
# overfitting the data

r.sd <- 0.0005

# 1 reference data set 1: dentine micromill
n.r1 <- length(tusk.mill.tl$Sr)

# reference data, timeline, Sr, and sd
R.r1 <- tusk.mill.tl$Sr

r1.tl <- tusk.mill.tl$tl

R.sd.r1 <- rep(r.sd, n.r1)

# 2 reference data set 2: Enamel 9 transect
n.En1 <- length(En1.50avg.tl.f$avg.sr)

# measurement data, timeline, Sr, and sd
R.En1 <- En1.50avg.tl.f$avg.sr

En1.tl <- En1.50avg.tl.f$avg.tl

# 3 conventional drill data set
n.drill <- length(drill.tl.f$Sr)

# measurement data, timeline, Sr, and sd
R.drill <- drill.tl.f$Sr

R.sd.drill <-rep(r.sd, n.drill)

drill.tl <- drill.tl.f$tl

# 4 micromill enamel data set
n.Rm3.5b <- length(Rm3.5b.mill.tl$Sr)

# measurement data, timeline, Sr, and sd
R.Rm3.5b <- Rm3.5b.mill.tl$Sr

R.sd.Rm3.5b <-rep(r.sd, n.Rm3.5b)

Rm3.5b.tl<- Rm3.5b.mill.tl$tl

# 5 Enamel 9 transect
n.En9 <- length(En9.50avg.tl$avg.sr)

# measurement data, timeline, Sr, and sd
R.En9 <- En9.50avg.tl$avg.sr

R.sd.En9 <-rep(r.sd, n.En9)

En9.tl <- En9.50avg.tl$avg.tl

# 6 Enamel 10 transect
n.En10 <- length(En10.50avg.tl$avg.sr)

# measurement data, timeline, Sr, and sd
R.En10 <- En10.50avg.tl$avg.sr

# R.sd.En10 <- En10.50avg.tl$sd.sr * 5

R.sd.En10 <-rep(r.sd, n.En10)

En10.tl <- En10.50avg.tl$avg.tl

####### end compiling parapeters #######

####### Preparing the model run #######
parameters <- c("switch", "Rin", "R1.m", "R2.m", "a","b","c",
                "R1.drill","R1.Rm3.5b","R1.En9", "R1.En10",
                "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                "Rpri.mod", "Raft.mod")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            tl.sd = tl.sd, bin.thin = bin.thin, t = t,
            n.r1 = n.r1, R.r1 = R.r1, R.sd.r1 = rep(r.sd, n.r1), r1.tl = r1.tl,
            n.r2 = n.En1, R.r2 = R.En1, R.sd.r2 = rep(r.sd, n.r2), r2.tl = En1.tl,
            n.drill = n.drill, R.drill = R.drill, R.sd.drill = rep(r.sd, n.drill), drill.tl = drill.tl,
            n.Rm3.5b = n.Rm3.5b, R.Rm3.5b = R.Rm3.5b, R.sd.Rm3.5b = rep(r.sd, n.Rm3.5b), Rm3.5b.tl = Rm3.5b.tl,
            n.En9 = n.En9, R.En9 = R.En9, R.sd.En9 = rep(r.sd, n.En9), En9.tl = En9.tl,
            n.En10 = n.En10, R.En10 = R.En10, R.sd.En10 = rep(r.sd, n.En10), En10.tl = En10.tl,
            UT.Sr = UT.Sr, CA.Sr = CA.Sr)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e5 #100 k iterations
n.burnin = 4e4 #40 k burnin
n.thin = 60

#Run it
post.comb = do.call(jags.parallel,list(model.file = "./code/Overprint comb JAGS.R", 
                                        parameters.to.save = parameters, 
                                        data = dat, n.chains=5, n.iter = n.iter, 
                                        n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~3 hours

save(post.comb, file = "./out/post.comb.RData")

#check convergence
denplot(as.mcmc(post.comb), parms = c("a","b","c",
                                     "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                                     "Rpri.mod", "Raft.mod"))

# prepare timeline of R1 and mixed R1 for plotting
post.comb.R1m.89 <- MCMC.CI.bound(post.comb$BUGSoutput$sims.list$R1.m, 0.89)

post.comb.R1.drill.89 <- MCMC.CI.bound(post.comb$BUGSoutput$sims.list$R1.drill, 0.89)

post.comb.R1.Rm3.5b.89 <- MCMC.CI.bound(post.comb$BUGSoutput$sims.list$R1.Rm3.5b, 0.89)

post.comb.R1.En9.89 <- MCMC.CI.bound(post.comb$BUGSoutput$sims.list$R1.En9, 0.89)

post.comb.R1.En10.89 <- MCMC.CI.bound(post.comb$BUGSoutput$sims.list$R1.En10, 0.89)

# record MAP and CI of estimated proportion
pr.drill.map <- map_estimate(post.comb$BUGSoutput$sims.list$pr.drill)[[1]]
# 0.2296887
pr.drill.ci1 <- hdi(post.comb$BUGSoutput$sims.list$pr.drill, .89)[[2]]
# 0.1776457
pr.drill.ci2 <- hdi(post.comb$BUGSoutput$sims.list$pr.drill, .89)[[3]]
# 0.2868468

pr.Rm3.5b.map <- map_estimate(post.comb$BUGSoutput$sims.list$pr.Rm3.5b)[[1]]
# 0.07984524
pr.Rm3.5b.ci1 <- hdi(post.comb$BUGSoutput$sims.list$pr.Rm3.5b, .89)[[2]]
# 0.03396436
pr.Rm3.5b.ci2 <- hdi(post.comb$BUGSoutput$sims.list$pr.Rm3.5b, .89)[[3]]
# 0.1290464

pr.En9.map <- map_estimate(post.comb$BUGSoutput$sims.list$pr.En9)[[1]]
# 0.06070172
pr.En9.ci1 <- hdi(post.comb$BUGSoutput$sims.list$pr.En9, .89)[[2]]
# 0.02131612
pr.En9.ci2 <- hdi(post.comb$BUGSoutput$sims.list$pr.En9, .89)[[3]]
# 0.09912549

pr.En10.map <- map_estimate(post.comb$BUGSoutput$sims.list$pr.En10)[[1]]
# 0.1313421
pr.En10.ci1 <- hdi(post.comb$BUGSoutput$sims.list$pr.En10, .89)[[2]]
# 0.09397739
pr.En10.ci2 <- hdi(post.comb$BUGSoutput$sims.list$pr.En10, .89)[[3]]
# 0.1748904

# record Sr ratios prior to and after the movement
Sr.pri.map <- map_estimate(post.comb$BUGSoutput$sims.list$Rpri.mod)[[1]]

Sr.aft.map <- map_estimate(post.comb$BUGSoutput$sims.list$Raft.mod)[[1]]

# record reaction rate parameters for the forward model
# they need to be scaled back with bin.thin for the forward model
MAP.a <- map_estimate(post.comb$BUGSoutput$sims.list$a)[[1]]/bin.thin
MAP.b <- map_estimate(post.comb$BUGSoutput$sims.list$b)[[1]]/bin.thin
MAP.c <- map_estimate(post.comb$BUGSoutput$sims.list$c)[[1]]/bin.thin

####### end of model run #######


###################################
####### Sensitivity test #######
####### Preparing the model run #######
####### start compiling parameters #######

# here we focus on using measurement uncertainties to compile the priors
R.sd.r1 <- tusk.mill.tl$sd

R.sd.En1 <- En1.50avg.tl.f$sd.sr

R.sd.drill <- drill.tl.f$sd

R.sd.Rm3.5b <- Rm3.5b.mill.tl$sd 

R.sd.En9 <- En9.50avg.tl$sd.sr 

R.sd.En10 <- En10.50avg.tl$sd.sr

# note that the solution method is associated with much lower uncertainties
# than LA-ICP-MS

# other input data are the same compared to the default model

####### end compiling parapeters #######

####### Preparing the model run #######
parameters <- c("switch", "Rin", "R1.m", "R2.m", "a","b","c",
                "R1.drill","R1.Rm3.5b","R1.En9", "R1.En10",
                "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                "Rpri.mod", "Raft.mod")

dat = list( params.mu = turnover.params.mu, params.vcov = turnover.params.vcov, 
            tl.sd = tl.sd, bin.thin = bin.thin, t = t,
            n.r1 = n.r1, R.r1 = R.r1, R.sd.r1 = R.sd.r1, r1.tl = r1.tl,
            n.r2 = n.En1, R.r2 = R.En1, R.sd.r2 = R.sd.En1, r2.tl = En1.tl,
            n.drill = n.drill, R.drill = R.drill, R.sd.drill = R.sd.drill, drill.tl = drill.tl,
            n.Rm3.5b = n.Rm3.5b, R.Rm3.5b = R.Rm3.5b, R.sd.Rm3.5b = R.sd.Rm3.5b, Rm3.5b.tl = Rm3.5b.tl,
            n.En9 = n.En9, R.En9 = R.En9, R.sd.En9 = R.sd.En9, En9.tl = En9.tl,
            n.En10 = n.En10, R.En10 = R.En10, R.sd.En10 = R.sd.En10, En10.tl = En10.tl,
            UT.Sr = UT.Sr, CA.Sr = CA.Sr)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e5 #100 k iterations
n.burnin = 4e4 #40 k burnin
n.thin = 60

#Run it
post.sens = do.call(jags.parallel,list(model.file = "./code/Overprint comb JAGS.R", 
                                       parameters.to.save = parameters, 
                                       data = dat, n.chains=5, n.iter = n.iter, 
                                       n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~3 hours

save(post.sens, file = "./out/post.sens.RData")

#check convergence
denplot(as.mcmc(post.sens), parms = c("a","b","c",
                                      "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                                      "Rpri.mod", "Raft.mod"))

# results show issues with convergence (denplot), and
# distorted posterior distributions in a, b and c parameters 

####### end of sensitivity test #######
