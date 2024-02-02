library(readxl)
library(dplyr)
library(coda)
library(R2jags)
library(ggplot2)
library(mcmcplots)
library(bayesplot)

# modeling the amount of maturation overprint in samples. 
# this was done on four data series: LA-ICP-MS Enamel 9, Enamel 10, Drill, and micromill series
# the tusk dentine micromill data (solution method) was used to generate a reference timeline

CA.Sr <- 0.70597
UT.Sr <- 0.71115

# posterior CA Sr CI 
CA.Sr.sd <- (0.706243 - CA.Sr)/2

UT.Sr.sd <- (0.7120174 - UT.Sr)/2

# data reduction for LA-ICP-MS

days.cumm.en9.al <- days.cumm.en9 - 205 
days.cumm.en10.al <- days.cumm.en10 - 195 

#visualize timeline adjustment
par(mfrow=c(1,1))
plot(-10, -10, col= alpha("lightcyan4", 0),
     pch=16, cex=1, xlim=c(-500,1000),ylim=c(0.704,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")

#CA and UT Sr measurements from Yang et al. 2023
abline(h = CA.Sr)

abline(h = UT.Sr)

lines(days.cumm.en9.al, proc.Enamel9.rm.f$avg, col= alpha("orange4", 0.9),lwd=2)

lines(days.cumm.en10.al, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.9),lwd=2)


En9.25avg.tl <- x.pt.avg(proc.Enamel9.rm.f$avg, days.cumm.en9.al, 25)

En10.25avg.tl <- x.pt.avg(proc.Enamel10.rm.f$avg, days.cumm.en10.al, 25)

#get the 1/2 average time interval as undertainty for timeline matching
En9.tl.sd <- mean(base::diff(En9.25avg.tl$avg.tl))/2
En10.tl.sd <- mean(base::diff(En10.25avg.tl$avg.tl))/2

# identify the max and min number of days to model
# hand drill samples represent the longest time span of samples
max(days.cumm.drill.al) - min(days.cumm.drill.al) #>2200 days

t <- 2500 # set a higher number of days to accomodate 

post.misha.pc2p3$BUGSoutput$sims.list$R1.m[1:30000,1:750]

# parameter estimates from Misha's turnover model, generated from Yang et al., 2023
turnover.params <- read.csv("./data/turnover_param_posterior")

a.post <- turnover.params$a
b.post <- turnover.params$b
c.post <- turnover.params$c
post.leng <- length(a.post)

n.mea <- length(En9.25avg.tl$avg.sr)

R.mea <- En9.25avg.tl$avg.sr

R.sd.mea <- En9.25avg.tl$sd.sr

tl.sd <- En9.tl.sd

parameters <- c("Ivo.rate","dist", "R1.m","Rin.m", "R2.m","a","b","c","exp.ab",
                "Rin.m.pre","a.m","b.m","c.m","Body.mass.m", "Body.mass", "Rin.m.cps.ac")

dat = list( post.leng=post.leng, 
            a.post=a.post, b.post=b.post, c.post=c.post,
            tl.sd = tl.sd,
            R.mea = R.mea, dist.mea = dist.mea, R.sd.mea = R.sd.mea, t = t, n.mea = n.mea)

#Start time
t1 = proc.time()

set.seed(t1[3])
n.iter = 1e4
n.burnin = 5e3
n.thin = 1

#Run it
post.misha.invmamm.i = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param mammi.R", 
                                                  parameters.to.save = parameters, 
                                                  data = dat, n.chains=5, n.iter = n.iter, 
                                                  n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 28 hours

save(post.misha.invmamm.i, file = "out/post.misha.invmamm.i.RData")



par(mfrow=c(1,4))
# Panel 1 molar dentine vs micromill tusk dentine 
plot(-1000, -1, col= "gray24",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from moves",
     main = "Tusk micromill vs molar D LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

polygon(c(days.cumm.den.al, rev(days.cumm.den.al)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, 
          rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(days.cumm.den.al,dent.rm.f2$avg,col= "gray24", lwd=2)
points(misha.tusk.micromill.tl.al, misha.tusk.micromill$corr..87Sr.86Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

# Panel 2 molar enamel1 vs micromill tusk dentine 
plot(-1000, -1, col= "gray24",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from moves",
     main = "Tusk micromill vs molar E1 LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

polygon(c(days.cumm.en1.al, rev(days.cumm.en1.al)), 
        c(proc.Enamel1.rm.f$avg + proc.Enamel1.rm.f$sd, 
          rev(proc.Enamel1.rm.f$avg - proc.Enamel1.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(days.cumm.en1.al, proc.Enamel1.rm.f$avg,col = alpha("orange", 0.9), lwd=2)
points(misha.tusk.micromill.tl.al, misha.tusk.micromill$corr..87Sr.86Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

# Panel 3 molar enamel vs micromill tusk dentine 
plot(days.cumm.drill.al, rev(Drill.no$corr..87Sr.86Sr), col= "red4",
     pch=16, cex = 2,
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from moves",
     main = "Tusk micromill vs molar E drill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)
points(misha.tusk.micromill.tl.al, misha.tusk.micromill$corr..87Sr.86Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

# Panel 4 molar enamel vs micromill molar enamel 
plot(Rm3.5b.mill.tl.al, rev(Rm3.5b.mill.no$corr..87Sr.86Sr), col= "cyan4",
     pch=16, cex = 2,
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from moves",
     main = "Tusk micromill vs molar E micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)
points(misha.tusk.micromill.tl.al, misha.tusk.micromill$corr..87Sr.86Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))