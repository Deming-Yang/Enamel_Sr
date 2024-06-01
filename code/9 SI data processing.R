
############# Time line reconstructions with individual data series in Fig 3D #############

# referenced enamel length at 52.7 mm
ref.en.length <- 85.7 # in mm

ref.tan <- approx(x = ref.length.v, y = Pred.tan.Rm3.5, xout = ref.en.length)$y

# average enamel growth rate 55.3e-3 # mm/day, which is true for most of the crown
ref.en.growth <- 55.3e-3 # mm/day

# estimate a vector of enamel growth rates at the established length interval
Ref.growth.v <- ref.en.growth * ref.tan / Pred.tan.Rm3.5 # mm/day

# linear approximation of measured growth rates from the reference vector
# for selected data based on the molar measurements

# 1. for the dentine trasect
rate.denx <- approx(x = ref.length.v, y = Ref.growth.v, xout = dent.rm.new.x) # mm per day

# 2. for selected enamel LA trasect 1
rate.en1x <- approx(x = ref.length.v, y = Ref.growth.v, xout = proc.Enamel1.rm.new.x) 

# for the hand drill samples, x need to be reversed
Drill.newx <- rev(Drill.no$Dist..From.cervix + 0.05) # to avoid 0 in log calculations

rate.Drill <- approx(x = ref.length.v, y = Ref.growth.v, xout = Drill.newx) # mm per day

# calculate the differential of the sample lengths
# for the dentine trasect
interv.denx <- c(base::diff(-1*dent.rm.new.x)[1], base::diff(-1*dent.rm.new.x))

# for selected enamel LA trasects: 1, 5, and 10
interv.en1x <- c(base::diff(-1*proc.Enamel1.rm.new.x)[1], base::diff(-1*proc.Enamel1.rm.new.x))

# for the hand drill samples, the order has to be reversed
interv.Drillx <- c(-1*base::diff(Drill.newx)[1], -1*base::diff(Drill.newx))

# calculate time intervals
time.interv.den <- interv.denx / rate.denx$y
time.interv.en1 <- interv.en1x / rate.en1x$y
time.interv.drill <- interv.Drillx / rate.Drill$y

# calculate cumulative days
days.cumm.den <- base::cumsum(time.interv.den)
days.cumm.en1 <- base::cumsum(time.interv.en1)
days.cumm.drill <- base::cumsum(time.interv.drill)

# align with the switch, this step was done manually
days.cumm.den.al <- days.cumm.den - 910
days.cumm.en1.al <- days.cumm.en1 - 410 
days.cumm.drill.al <- days.cumm.drill - 530

### 3 enamel micromill translation ###
# transform this into time, not vertical distance, use trigonometry
# micromill angle is approx. 3.2 degrees (Uno et al., 2020)
# enamel daily secretion rate is ~3.5 microns/day (Uno et al., 2020)

Rm3.5b.mill.tl <- rev(Rm3.5b.mill.no$Dist..From.EDJ/3.5) # in days
Rm3.5b.mill.tl.al <- Rm3.5b.mill.tl - 120

# use micromill tusk dentine data in the comparison
# in mm, micromill was done at 0.5 mm interval
# from the pulp cavity, which is the newest dentine
maxl.tusk <- max(misha.tusk.micromill$Position)

#this is the total length of dentine milled, 500 micron interval
misha.tusk.micromill.dist <- (maxl.tusk - misha.tusk.micromill$Position) * 500 

# assuming a constant tusk growth rate at 14.7 microns/day (Uno et al., 2020)
misha.tusk.micromill.tl <- misha.tusk.micromill.dist/14.7 #mm/day

# align with the switch
misha.tusk.micromill.tl.al <- misha.tusk.micromill.tl - 335

#################################################
# prep for Fig S4
# examine the relationships between the residuals between 
# hand drill 87Sr/86Sr and the mixed serum as enamel input

# filter the data and focus on the timeline between 0 and 600

drill.tl.f2 <- dplyr::filter(Edrill.tl, Edrill.tl$tl > 0 & Edrill.tl$tl < 600)

pred.drill <- approx(x = bin.thin.oc*(1:t.oc) - 400, y = post.comb.R1.drill.89[[1]]
                     , xout = drill.tl.f2$tl)$y
# prediction: larger depths, lower 87Sr/86Sr, more deviation from predicted value
res.drill <- drill.tl.f2$Sr - pred.drill

res.tib <- tibble(x = drill.tl.f2$depth, y = res.drill)

lm.res <- lm(data = res.tib, y ~ x)
summary(lm.res) # results are significant

# calculate 95% confidence intervals for plotting 
newx = seq(min(res.tib$x), max(res.tib$x), by = 0.01)
conf_interval <- predict(lm.res, newdata=data.frame(x=newx), interval="confidence",
                         level = 0.95)

# comparing LA-ICP-MS tusk dentine and molar dentine series

misha.raw <- read.csv("data/Misha M640b raw.csv")

misha.raw.Sr <- misha.raw$Corr..87Sr.86Sr
misha.raw.dist <- misha.raw$dist

##calculate 25 point average to reduce data amount
n.avg <- 25 #now many data points to average
n <- floor(length(misha.raw$Corr..87Sr.86Sr)/n.avg) #99 data points, discard the last <50 data points

#initiate vectors
n.avg.misha.25.sr <- rep(0, n)
n.sd.misha.25.sr <- rep(0, n)

n.avg.misha.25.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.25.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.25.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.25.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}

#timeline reconstruction with simple growth rate
n.avg.misha.25.tl <- rev(n.avg.misha.25.dist/14.7) #mm/day

#mannually align data
n.avg.misha.25.tl.al <- n.avg.misha.25.tl - 670

########################################################################################
############## Code below is not used in the article ################################### 
########################################################################################

####### Inversion of more data transects #######
####### data set 2, enamel transect 3#############

En3.50avg.tl <- x.pt.avg(proc.Enamel3.rm.f$mov.avg, en3.tl$tl, 50)

Enamel3.Sr <- x.pt.avg(proc.Enamel3.rm.f$mov.avg, 0, 50)

Enamel3.length <- x.pt.avg(proc.Enamel3.rm.new.x, 0, 50)

n.mea = length(Enamel3.length$avg.sr)

dist.mea <- Enamel3.length$avg.sr
R.mea <- Enamel3.Sr$avg.sr
n.mea = length(R.mea)
R.sd.mea <- Enamel3.Sr$sd.sr

s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
length.v <- rev(ref.length.v)
rate.v <- rev(Ref.growth.v)

rate.sd <- 0.05 * mean(rate.v) # mm per day

bin.thin <- 0.25 * s.intv/min(rate.v)

d.offset.En3 <- 2* s.intv/mean(rate.v)

En3.bt <- bin.thin

t = round(1.2 * max(days.cumm.en3)/bin.thin)

t.En3 <-t

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
post.misha.En3 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 mins

save(post.misha.En3, file = "out/post.misha.En3.RData")

load("out/post.misha.En3.RData")

denplot(as.mcmc(post.misha.En3), parms = c("a","b","c"))

post.misha.En3.R1.m.89 <- MCMC.CI.bound(post.misha.En3$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En3.Rin.m.89 <- MCMC.CI.bound(post.misha.En3$BUGSoutput$sims.list$Rin.m, 0.89)


####### data set 3, enamel transect 5#############

En5.50avg.tl <- x.pt.avg(proc.Enamel5.rm.f$mov.avg, en5.tl$tl, 50)

Enamel5.Sr <- x.pt.avg(proc.Enamel5.rm.f$mov.avg, 0, 50)

Enamel5.length <- x.pt.avg(proc.Enamel5.rm.new.x, 0, 50)

n.mea = length(Enamel5.length$avg.sr)

dist.mea <- Enamel5.length$avg.sr
R.mea <- Enamel5.Sr$avg.sr
n.mea = length(R.mea)
R.sd.mea <- Enamel5.Sr$sd.sr

s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
length.v <- rev(ref.length.v)
rate.v <- rev(Ref.growth.v)

rate.sd <- 0.05 * mean(rate.v) # mm per day

bin.thin <- 0.25 * s.intv/min(rate.v)

d.offset.En5 <- 2* s.intv/mean(rate.v)

En5.bt <- bin.thin

t = round(1.2 * max(days.cumm.en5)/bin.thin)

t.En5 <-t

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
post.misha.En5 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 mins

save(post.misha.En5, file = "out/post.misha.En5.RData")

load("out/post.misha.En5.RData")

denplot(as.mcmc(post.misha.En5), parms = c("a","b","c"))

post.misha.En5.R1.m.89 <- MCMC.CI.bound(post.misha.En5$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En5.Rin.m.89 <- MCMC.CI.bound(post.misha.En5$BUGSoutput$sims.list$Rin.m, 0.89)


####### data set 4, enamel transect 7#############
En7.50avg.tl <- x.pt.avg(proc.Enamel7.rm.f$mov.avg, en7.tl$tl, 50)

Enamel7.Sr <- x.pt.avg(proc.Enamel7.rm.f$mov.avg, 0, 50)

Enamel7.length <- x.pt.avg(proc.Enamel7.rm.new.x, 0, 50)

n.mea = length(Enamel7.length$avg.sr)

dist.mea <- Enamel7.length$avg.sr
R.mea <- Enamel7.Sr$avg.sr
n.mea = length(R.mea)
R.sd.mea <- Enamel7.Sr$sd.sr

s.intv <- - 1 * mean(diff(dist.mea))

max.dist.mea <- max(dist.mea) + 2 * s.intv # create a bit of leading enamel

# use a growth length to rate relationship for growth simulations
# the xgrid has to be in ascending order
length.v <- rev(ref.length.v)
rate.v <- rev(Ref.growth.v)

rate.sd <- 0.05 * mean(rate.v) # mm per day

bin.thin <- 0.25 * s.intv/min(rate.v)

d.offset.En7 <- 2* s.intv/mean(rate.v)

En7.bt <- bin.thin

t = round(1.2 * max(days.cumm.en7)/bin.thin)

t.En7 <-t

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
post.misha.En7 = do.call(jags.parallel,list(model.file = "code/Sr inversion JAGS param.R", 
                                            parameters.to.save = parameters, 
                                            data = dat, n.chains=5, n.iter = n.iter, 
                                            n.burnin = n.burnin, n.thin = n.thin))

#Time taken
proc.time() - t1 #~ 20 mins

save(post.misha.En7, file = "out/post.misha.En7.RData")

load("out/post.misha.En7.RData")

denplot(as.mcmc(post.misha.En7), parms = c("a","b","c"))

post.misha.En7.R1.m.89 <- MCMC.CI.bound(post.misha.En7$BUGSoutput$sims.list$R1.m, 0.89)

post.misha.En7.Rin.m.89 <- MCMC.CI.bound(post.misha.En7$BUGSoutput$sims.list$Rin.m, 0.89)

############# Export all BITS data output to .csv #############

# Misha's Sr intake from tusk micromill (M640b)

M640b.Rin.tl <- 1:t.M640b * M640.bt -d.offset.M640-365

M640b.Rin.MAPE <- post.misha.M640b.Rin.m.89[[1]]

M640b.Rin.CIL <- post.misha.M640b.Rin.m.89[[2]]

M640b.Rin.CIH <- post.misha.M640b.Rin.m.89[[3]]

M640b.Rin.est <- data.frame(timeline = M640b.Rin.tl, Sr.est = M640b.Rin.MAPE, 
           Sr.est.5pct = M640b.Rin.CIL, Sr.est.95pct = M640b.Rin.CIH)

write.csv(M640b.Rin.est, file = "out/Sr intake est M640b.csv")

#############
# Misha's Sr intake from molar LA-ICP-MS (transect 1)

ICPMS.en1.Rin.tl <- 1:t.En1 * En1.bt - d.offset.En1 - offset.en1

ICPMS.en1.Rin.MAPE <- post.misha.En1.Rin.m.89[[1]]

ICPMS.en1.Rin.CIL <- post.misha.En1.Rin.m.89[[2]]

ICPMS.en1.Rin.CIH <- post.misha.En1.Rin.m.89[[3]]

ICPMS.en1.Rin.est <- data.frame(timeline = ICPMS.en1.Rin.tl, Sr.est = ICPMS.en1.Rin.MAPE, 
                            Sr.est.5pct = ICPMS.en1.Rin.CIL, Sr.est.95pct = ICPMS.en1.Rin.CIH)

write.csv(ICPMS.en1.Rin.est, file = "out/Sr intake est LA-ICP-MS en1.csv")


