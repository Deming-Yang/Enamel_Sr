library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(segmented)

source("./code/1 Helper functions.R")

######################## change point analysis ########################

############ Dentine combined ###################

dent.rm.f <- filter(dent.rm.proj, dent.rm.proj$avg > 0.703)

fit_lm.D.rm = lm(avg ~ 1 + new.x, data = dent.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.D.rm = segmented(fit_lm.D.rm, seg.Z = ~new.x, npsi = 3)  # Three change points along new.x

summary(fit_segmented.D.rm)

cp1.D.rm <- fit_segmented.D.rm$psi[5] #change point estimates, can use this to get x for plotting
cp1.D.rm.err <- fit_segmented.D.rm$psi[8]

cp2.D.rm <- fit_segmented.D.rm$psi[6] #change point estimates, can use this to get x for plotting
cp2.D.rm.err <- fit_segmented.D.rm$psi[9]

D.rm.sl1 <- fit_segmented.D.rm$coefficients[2] + fit_segmented.D.rm$coefficients[3]+ fit_segmented.D.rm$coefficients[4] #second slope

D.rm.sl2 <- fit_segmented.D.rm$coefficients[2]+ fit_segmented.D.rm$coefficients[3]+ fit_segmented.D.rm$coefficients[4] + 
  fit_segmented.D.rm$coefficients[5]#third slope

# D.rm.sl3 <- fit_segmented.D.rm$coefficients[2]+ fit_segmented.D.rm$coefficients[3]+ 
#   fit_segmented.D.rm$coefficients[4] + fit_segmented.D.rm$coefficients[5]#third slope
############Enamel 1###################
proc.Enamel1.rm.f <- filter(proc.Enamel1.rm.proj, proc.Enamel1.rm.proj$avg > 0.703)
fit_lm.E1.rm = lm(avg ~ 1 + new.x, data = proc.Enamel1.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E1.rm = segmented(fit_lm.E1.rm, seg.Z = ~new.x, npsi = 2)  # Two change points along y

summary(fit_segmented.E1.rm)

cp1.E1.rm <- fit_segmented.E1.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E1.rm.err <- fit_segmented.E1.rm$psi[5]

cp2.E1.rm <- fit_segmented.E1.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E1.rm.err <- fit_segmented.E1.rm$psi[6]

# #extract slopes and intercepts
E1.rm.sl0 <- fit_segmented.E1.rm$coefficients[2]

E1.rm.sl1 <- fit_segmented.E1.rm$coefficients[2] + fit_segmented.E1.rm$coefficients[3] #second slope

E1.rm.sl2 <- fit_segmented.E1.rm$coefficients[2]+ fit_segmented.E1.rm$coefficients[3]+ fit_segmented.E1.rm$coefficients[4] #third slope

############Enamel 2###################
proc.Enamel2.rm.f <- filter(proc.Enamel2.rm.proj, proc.Enamel2.rm.proj$avg > 0.703)
fit_lm.E2.rm = lm(avg ~ 1 + new.x, data = proc.Enamel2.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E2.rm = segmented(fit_lm.E2.rm, seg.Z = ~new.x, npsi = 2)  # Two change points along y

summary(fit_segmented.E2.rm)

cp1.E2.rm <- fit_segmented.E2.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E2.rm.err <- fit_segmented.E2.rm$psi[5]

cp2.E2.rm <- fit_segmented.E2.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E2.rm.err <- fit_segmented.E2.rm$psi[6]

# #extract slopes and intercepts
E2.rm.sl0 <- fit_segmented.E2.rm$coefficients[2]

E2.rm.sl1 <- fit_segmented.E2.rm$coefficients[2] + fit_segmented.E2.rm$coefficients[3] #second slope

E2.rm.sl2 <- fit_segmented.E2.rm$coefficients[2]+ fit_segmented.E2.rm$coefficients[3]+ fit_segmented.E2.rm$coefficients[4] #third slope

############Enamel 3###################
proc.Enamel3.rm.f <- filter(proc.Enamel3.rm.proj, proc.Enamel3.rm.proj$avg > 0.703)
fit_lm.E3.rm = lm(avg ~ 1 + new.x, data = proc.Enamel3.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E3.rm = segmented(fit_lm.E3.rm, seg.Z = ~new.x, npsi = 2)  # Two change points along y

summary(fit_segmented.E3.rm)

cp1.E3.rm <- fit_segmented.E3.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E3.rm.err <- fit_segmented.E3.rm$psi[5]

cp2.E3.rm <- fit_segmented.E3.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E3.rm.err <- fit_segmented.E3.rm$psi[6]

# #extract slopes and intercepts
E3.rm.sl0 <- fit_segmented.E3.rm$coefficients[2]

E3.rm.sl1 <- fit_segmented.E3.rm$coefficients[2] + fit_segmented.E3.rm$coefficients[3] #second slope

E3.rm.sl2 <- fit_segmented.E3.rm$coefficients[2]+ fit_segmented.E3.rm$coefficients[3]+ fit_segmented.E3.rm$coefficients[4] #third slope

############Enamel 4###################
proc.Enamel4.rm.f <- filter(proc.Enamel4.rm.proj, proc.Enamel4.rm.proj$avg > 0.703)
fit_lm.E4.rm = lm(avg ~ 1 + new.x, data = proc.Enamel4.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E4.rm = segmented(fit_lm.E4.rm, seg.Z = ~new.x, npsi = 4)  # four change points along y to make sense

summary(fit_segmented.E4.rm)

fit_segmented.E4.rm$psi

cp1.E4.rm <- fit_segmented.E4.rm$psi[5] #change point estimates, can use this to get x for plotting
cp1.E4.rm.err <- fit_segmented.E4.rm$psi[9]

cp2.E4.rm <- fit_segmented.E4.rm$psi[6] #change point estimates, can use this to get x for plotting
cp2.E4.rm.err <- fit_segmented.E4.rm$psi[10]

#extract slopes and intercepts
E4.rm.sl0 <- fit_segmented.E4.rm$coefficients[2]

E4.rm.sl1 <- fit_segmented.E4.rm$coefficients[2] + fit_segmented.E4.rm$coefficients[3] #second slope

E4.rm.sl2 <- fit_segmented.E4.rm$coefficients[2]+ fit_segmented.E4.rm$coefficients[3]+ fit_segmented.E4.rm$coefficients[4]+
  fit_segmented.E4.rm$coefficients[5]#third slope

############Enamel 5###################
proc.Enamel5.rm.f <- filter(proc.Enamel5.rm.proj, proc.Enamel5.rm.proj$avg > 0.703)
fit_lm.E5.rm = lm(avg ~ 1 + new.x, data = proc.Enamel5.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E5.rm = segmented(fit_lm.E5.rm, seg.Z = ~new.x, npsi = 2)  # Two change points along y

summary(fit_segmented.E5.rm)

cp1.E5.rm <- fit_segmented.E5.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E5.rm.err <- fit_segmented.E5.rm$psi[5]

cp2.E5.rm <- fit_segmented.E5.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E5.rm.err <- fit_segmented.E5.rm$psi[6]

#extract slopes and intercepts
E5.rm.sl0 <- fit_segmented.E5.rm$coefficients[2]

E5.rm.sl1 <- fit_segmented.E5.rm$coefficients[2] + fit_segmented.E5.rm$coefficients[3] #second slope

E5.rm.sl2 <- fit_segmented.E5.rm$coefficients[2]+ fit_segmented.E5.rm$coefficients[3]+ fit_segmented.E5.rm$coefficients[4] #third slope

############Enamel 6###################
proc.Enamel6.rm.f <- filter(proc.Enamel6.rm.proj, proc.Enamel6.rm.proj$avg > 0.703)

fit_lm.E6.rm = lm(avg ~ 1 + new.x, data = proc.Enamel6.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E6.rm = segmented(fit_lm.E6.rm, seg.Z = ~new.x, npsi = 3)  # one change points along y

summary(fit_segmented.E6.rm)

cp1.E6.rm <- fit_segmented.E6.rm$psi[5] #change point estimates, can use this to get x for plotting
cp1.E6.rm.err <- fit_segmented.E6.rm$psi[8]

cp2.E6.rm <- fit_segmented.E6.rm$psi[6] #change point estimates, can use this to get x for plotting
cp2.E6.rm.err <- fit_segmented.E6.rm$psi[9]

#extract slopes and intercepts

E6.rm.sl0 <- fit_segmented.E6.rm$coefficients[2] + fit_segmented.E6.rm$coefficients[3] #second slope

E6.rm.sl1 <- fit_segmented.E6.rm$coefficients[2]+ fit_segmented.E6.rm$coefficients[3]+ fit_segmented.E6.rm$coefficients[4] #third slope

E6.rm.sl2 <- fit_segmented.E6.rm$coefficients[2]+ fit_segmented.E6.rm$coefficients[3]+ fit_segmented.E6.rm$coefficients[4] +
  fit_segmented.E6.rm$coefficients[5]#third slope
############Enamel 7###################
proc.Enamel7.rm.f <- filter(proc.Enamel7.rm.proj, proc.Enamel7.rm.proj$avg > 0.703)
fit_lm.E7.rm = lm(avg ~ 1 + new.x, data = proc.Enamel7.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E7.rm = segmented(fit_lm.E7.rm, seg.Z = ~new.x, npsi = 3)  # one change point to make it work

summary(fit_segmented.E7.rm)

cp1.E7.rm <- fit_segmented.E7.rm$psi[4] #change point estimates, can use this to get x for plotting
cp1.E7.rm.err <- fit_segmented.E7.rm$psi[7] 

cp2.E7.rm <- fit_segmented.E7.rm$psi[5] #change point estimates, can use this to get x for plotting
cp2.E7.rm.err <- fit_segmented.E7.rm$psi[8] 

#extract slopes and intercepts
E7.rm.sl0 <- fit_segmented.E7.rm$coefficients[2]

E7.rm.sl1 <- fit_segmented.E7.rm$coefficients[2] + fit_segmented.E7.rm$coefficients[3] #second slope

E7.rm.sl2 <- fit_segmented.E7.rm$coefficients[2]+ fit_segmented.E7.rm$coefficients[3]+ fit_segmented.E7.rm$coefficients[4]+
 fit_segmented.E7.rm$coefficients[5]#third slope

############Enamel 8###################
proc.Enamel8.rm.f <- filter(proc.Enamel8.rm.proj, proc.Enamel8.rm.proj$avg > 0.703)
fit_lm.E8.rm = lm(avg ~ 1 + new.x, data = proc.Enamel8.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E8.rm = segmented(fit_lm.E8.rm, seg.Z = ~new.x, npsi = 2)  # one change point to make it work

summary(fit_segmented.E8.rm)

cp1.E8.rm <- fit_segmented.E8.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E8.rm.err <- fit_segmented.E8.rm$psi[5]

cp2.E8.rm <- fit_segmented.E8.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E8.rm.err <- fit_segmented.E8.rm$psi[6]

#extract slopes and intercepts
E8.rm.sl0 <- fit_segmented.E8.rm$coefficients[2]

E8.rm.sl1 <- fit_segmented.E8.rm$coefficients[2] + fit_segmented.E8.rm$coefficients[3] #second slope

E8.rm.sl2 <- fit_segmented.E8.rm$coefficients[2]+ fit_segmented.E8.rm$coefficients[3]+ fit_segmented.E8.rm$coefficients[4] #third slope


############Enamel 9###################
proc.Enamel9.rm.f <- filter(proc.Enamel9.rm.proj, proc.Enamel9.rm.proj$avg > 0.703)
fit_lm.E9.rm = lm(avg ~ 1 + new.x, data = proc.Enamel9.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E9.rm = segmented(fit_lm.E9.rm, seg.Z = ~new.x, npsi = 2)  # one change point to make it work

summary(fit_segmented.E9.rm)

cp1.E9.rm <- fit_segmented.E9.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E9.rm.err <- fit_segmented.E9.rm$psi[5]

cp2.E9.rm <- fit_segmented.E9.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E9.rm.err <- fit_segmented.E9.rm$psi[6]

#extract slopes and intercepts
E9.rm.sl0 <- fit_segmented.E9.rm$coefficients[2]

E9.rm.sl1 <- fit_segmented.E9.rm$coefficients[2] + fit_segmented.E9.rm$coefficients[3] #second slope

E9.rm.sl2 <- fit_segmented.E9.rm$coefficients[2]+ fit_segmented.E9.rm$coefficients[3]+ fit_segmented.E9.rm$coefficients[4] #third slope

############Enamel 10###################
proc.Enamel10.rm.f <- filter(proc.Enamel10.rm.proj, proc.Enamel10.rm.proj$avg > 0.703)
fit_lm.E10.rm = lm(avg ~ 1 + new.x, data = proc.Enamel10.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.E10.rm = segmented(fit_lm.E10.rm, seg.Z = ~new.x, npsi = 2)  # one change point to make it work

summary(fit_segmented.E10.rm)



cp1.E10.rm <- fit_segmented.E10.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E10.rm.err <- fit_segmented.E10.rm$psi[5]

cp2.E10.rm <- fit_segmented.E10.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E10.rm.err <- fit_segmented.E10.rm$psi[6]

#extract slopes and intercepts
E10.rm.sl0 <- fit_segmented.E10.rm$coefficients[2]

E10.rm.sl1 <- fit_segmented.E10.rm$coefficients[2] + fit_segmented.E10.rm$coefficients[3] #second slope

E10.rm.sl2 <- fit_segmented.E10.rm$coefficients[2]+ fit_segmented.E10.rm$coefficients[3]+ fit_segmented.E10.rm$coefficients[4] #third slope


#compile slopes of the first transition
Sr.trans.sl1 <- c(D.rm.sl1, E1.rm.sl1, E2.rm.sl1, E3.rm.sl1, E4.rm.sl1, E5.rm.sl1, E6.rm.sl1, E7.rm.sl1,
                  E8.rm.sl1, E9.rm.sl1, E10.rm.sl1)

#compile slopes of the second transition
Sr.trans.sl2 <- c(D.rm.sl2, E1.rm.sl2, E2.rm.sl2, E3.rm.sl2, E4.rm.sl2, E5.rm.sl2, E6.rm.sl2, E7.rm.sl2,
                  E8.rm.sl2, E9.rm.sl2, E10.rm.sl2)


#################compiling first change points within enamel#################
cp1.E.rm.x <- c(cp1.E1.rm,
      cp1.E2.rm,
      cp1.E3.rm,
      cp1.E4.rm,
      cp1.E5.rm,
      cp1.E6.rm,
      cp1.E7.rm,
      cp1.E8.rm,
      cp1.E9.rm,
      cp1.E10.rm)

cp1.E.rm.x.err <- c(cp1.E1.rm.err,
             cp1.E2.rm.err,
             cp1.E3.rm.err,
             cp1.E4.rm.err,
             cp1.E5.rm.err,
             cp1.E6.rm.err,
             cp1.E7.rm.err,
             cp1.E8.rm.err,
             cp1.E9.rm.err,
             cp1.E10.rm.err)

#second change point within enamel
cp2.E.rm.x <- c(cp2.E1.rm,
                cp2.E2.rm,
                cp2.E3.rm,
                cp2.E4.rm,
                cp2.E5.rm,
                cp2.E6.rm,
                cp2.E7.rm,
                cp2.E8.rm,
                cp2.E9.rm,
                cp2.E10.rm)

cp2.E.rm.x.err <- c(cp2.E1.rm.err,
                    cp2.E2.rm.err,
                    cp2.E3.rm.err,
                    cp2.E4.rm.err,
                    cp2.E5.rm.err,
                    cp2.E6.rm.err,
                    cp2.E7.rm.err,
                    cp2.E8.rm.err,
                    cp2.E9.rm.err,
                    cp2.E10.rm.err)

#get y values from enamel transects
cp.D.rm.x <- approx(x = dent.rm.f$y, y = dent.rm.f$x, xout = cp.D.rm[1])$y

cp1.E.rm.y <- rep(0,10) #initialize vector
cp1.E.rm.y[1] <- approx(x = proc.Enamel1.rm.f$new.x, y = proc.Enamel1.rm.f$new.y, xout = cp1.E.rm.x[1])$y
cp1.E.rm.y[2] <- approx(x = proc.Enamel2.rm.f$new.x, y = proc.Enamel2.rm.f$new.y, xout = cp1.E.rm.x[2])$y
cp1.E.rm.y[3] <- approx(x = proc.Enamel3.rm.f$new.x, y = proc.Enamel3.rm.f$new.y, xout = cp1.E.rm.x[3])$y
cp1.E.rm.y[4] <- approx(x = proc.Enamel4.rm.f$new.x, y = proc.Enamel4.rm.f$new.y, xout = cp1.E.rm.x[4])$y
cp1.E.rm.y[5] <- approx(x = proc.Enamel5.rm.f$new.x, y = proc.Enamel5.rm.f$new.y, xout = cp1.E.rm.x[5])$y
cp1.E.rm.y[6] <- approx(x = proc.Enamel6.rm.f$new.x, y = proc.Enamel6.rm.f$new.y, xout = cp1.E.rm.x[6])$y
cp1.E.rm.y[7] <- approx(x = proc.Enamel7.rm.f$new.x, y = proc.Enamel7.rm.f$new.y, xout = cp1.E.rm.x[7])$y
cp1.E.rm.y[8] <- approx(x = proc.Enamel8.rm.f$new.x, y = proc.Enamel8.rm.f$new.y, xout = cp1.E.rm.x[8])$y
cp1.E.rm.y[9] <- approx(x = proc.Enamel9.rm.f$new.x, y = proc.Enamel9.rm.f$new.y, xout = cp1.E.rm.x[9])$y
cp1.E.rm.y[10] <- approx(x = proc.Enamel10.rm.f$new.x, y = proc.Enamel10.rm.f$new.y, xout = cp1.E.rm.x[10])$y

cp2.E.rm.y <- rep(0,10) #initialize vector
cp2.E.rm.y[1] <- approx(x = proc.Enamel1.rm.f$new.x, y = proc.Enamel1.rm.f$new.y, xout = cp2.E.rm.x[1])$y
cp2.E.rm.y[2] <- approx(x = proc.Enamel2.rm.f$new.x, y = proc.Enamel2.rm.f$new.y, xout = cp2.E.rm.x[2])$y
cp2.E.rm.y[3] <- approx(x = proc.Enamel3.rm.f$new.x, y = proc.Enamel3.rm.f$new.y, xout = cp2.E.rm.x[3])$y
cp2.E.rm.y[4] <- approx(x = proc.Enamel4.rm.f$new.x, y = proc.Enamel4.rm.f$new.y, xout = cp2.E.rm.x[4])$y
cp2.E.rm.y[5] <- approx(x = proc.Enamel5.rm.f$new.x, y = proc.Enamel5.rm.f$new.y, xout = cp2.E.rm.x[5])$y
cp2.E.rm.y[6] <- approx(x = proc.Enamel6.rm.f$new.x, y = proc.Enamel6.rm.f$new.y, xout = cp2.E.rm.x[6])$y
cp2.E.rm.y[7] <- approx(x = proc.Enamel7.rm.f$new.x, y = proc.Enamel7.rm.f$new.y, xout = cp2.E.rm.x[7])$y
cp2.E.rm.y[8] <- approx(x = proc.Enamel8.rm.f$new.x, y = proc.Enamel8.rm.f$new.y, xout = cp2.E.rm.x[8])$y
cp2.E.rm.y[9] <- approx(x = proc.Enamel9.rm.f$new.x, y = proc.Enamel9.rm.f$new.y, xout = cp2.E.rm.x[9])$y
cp2.E.rm.y[10] <- approx(x = proc.Enamel10.rm.f$new.x, y = proc.Enamel10.rm.f$new.y, xout = cp2.E.rm.x[10])$y

cp1.E.proj <- data.frame(x = cp1.E.rm.x, y = cp1.E.rm.y)
cp2.E.proj <- data.frame(x = cp2.E.rm.x, y = cp2.E.rm.y)

#it seems that only the last 2 data points are not conforming to a linear trend
#remove the points for E8-10
cp1.E.proj.inner <- cp1.E.proj[1:7,]

cp2.E.proj.inner <- cp2.E.proj[1:7,]

#calculate slope using lm
lm.cp1.E.proj.inner <- lm(y~x, data=cp1.E.proj.inner)
summary(lm.cp1.E.proj.inner)
#slope = -0.05838
#Std. Error = 0.003673

#calculate slope using lm
lm.cp2.E.proj.inner <- lm(y~x, data=cp2.E.proj.inner)
summary(lm.cp2.E.proj.inner)
#slope = -0.05824
#Std. Error = 0.001748

#the slopes are within the std error range from each other

#calculate appositional angle
atan(abs(lm.cp1.E.proj.inner$coefficients[2]))/pi*180
atan(abs(lm.cp2.E.proj.inner$coefficients[2]))/pi*180
#~3.37 degrees, which is identical to the appositional angle measured by Uno 2012 (3.3 +- 0.5)

#use mean of the two angles in the forward model
fwd.appo.sl <- mean(c(lm.cp1.E.proj.inner$coefficients[2],lm.cp2.E.proj.inner$coefficients[2]))

# compile change points for plots
#####use the sf package to plot these values
all.dat.rm <- rbind(dent.rm.f, proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                    proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                    proc.Enamel9.rm.f, proc.Enamel10.rm.f)

all.dat.rm.proj <- rbind(dent.rm.f, proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                         proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                         proc.Enamel9.rm.f, proc.Enamel10.rm.f)

all.enamel.rm <- rbind(proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                       proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                       proc.Enamel9.rm.f, proc.Enamel10.rm.f)

all.enamel.rm.proj <- rbind(proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                            proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                            proc.Enamel9.rm.f, proc.Enamel10.rm.f)

require(sf)
###########real world tooth geometry###########
sf.all.dat.rm <- st_as_sf(all.dat.rm,  agr = NA_agr_,
                          coords = c("y","x"),
                          dim = "XYZ",
                          remove = TRUE,
                          na.fail = TRUE,
                          sf_column_name = NULL)

###########EDJ flattened geometry########### in supplementary
sf.all.dat.rm.proj <- st_as_sf(all.dat.rm.proj,  agr = NA_agr_,
                               coords = c("new.x","new.y"),
                               dim = "XY",
                               remove = TRUE,
                               na.fail = TRUE,
                               sf_column_name = NULL)

sf.all.enamel.rm <- st_as_sf(all.enamel.rm,  agr = NA_agr_,
                             coords = c("y","x"),
                             dim = "XY",
                             remove = TRUE,
                             na.fail = TRUE,
                             sf_column_name = NULL)

sf.all.enamel.rm.proj <- st_as_sf(all.enamel.rm.proj,  agr = NA_agr_,
                                  coords = c("new.x","new.y"),
                                  dim = "XY",
                                  remove = TRUE,
                                  na.fail = TRUE,
                                  sf_column_name = NULL)




#record change positions within enamel:
cp1.loc <- data.frame(avg = rep(0.7,length(cp1.E.rm.x)), x = cp1.E.rm.x, y = cp1.E.rm.y)
sf.cp1.loc <- st_as_sf(cp1.loc,  agr = NA_agr_,
                       coords = c("x","y"),
                       dim = "XYZ",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)

cp2.loc <- data.frame(avg = rep(0.7,length(cp2.E.rm.x)), x = cp2.E.rm.x, y = cp2.E.rm.y)
sf.cp2.loc <- st_as_sf(cp2.loc,  agr = NA_agr_,
                       coords = c("x","y"),
                       dim = "XYZ",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)


#####################plot slopes of the first segment#####
# data.frame(slope1 = Sr.trans.sl1, slope2 = Sr.trans.sl2)
# 
# plot(1:11,Sr.trans.sl1) #no significant difference in the first slope
# plot(1:11,Sr.trans.sl2) #lower slopes for the outer two transects
# 
# newdataD <- data.frame(new.x = c(cp1.D.rm, cp2.D.rm))
# segmented.D.Sr <- predict.segmented(fit_segmented.D.rm, newdata = newdataD)
# Tr1.D <- segmented.D.Sr[2] - segmented.D.Sr[1]
# 
# Tr1.E <- rep(0,10)
# 
# newdataE1 <- data.frame(new.x = c(cp1.E1.rm, cp2.E1.rm))
# segmented.E1.Sr <- predict.segmented(fit_segmented.E1.rm, newdata = newdataE1)
# Tr1.E[1] <- segmented.E1.Sr[2] - segmented.E1.Sr[1]
# 
# newdataE2 <- data.frame(new.x = c(cp1.E2.rm, cp2.E2.rm))
# segmented.E2.Sr <- predict.segmented(fit_segmented.E2.rm, newdata = newdataE2)
# Tr1.E[2] <- segmented.E2.Sr[2] - segmented.E2.Sr[1]
# 
# newdataE3 <- data.frame(new.x = c(cp1.E3.rm, cp2.E3.rm))
# segmented.E3.Sr <- predict.segmented(fit_segmented.E3.rm, newdata = newdataE3)
# Tr1.E[3] <- segmented.E3.Sr[2] - segmented.E3.Sr[1]
# 
# newdataE4 <- data.frame(new.x = c(cp1.E4.rm, cp2.E4.rm))
# segmented.E4.Sr <- predict.segmented(fit_segmented.E4.rm, newdata = newdataE4)
# Tr1.E[4] <- segmented.E4.Sr[2] - segmented.E4.Sr[1]
# 
# newdataE5 <- data.frame(new.x = c(cp1.E5.rm, cp2.E5.rm))
# segmented.E5.Sr <- predict.segmented(fit_segmented.E5.rm, newdata = newdataE5)
# Tr1.E[5] <- segmented.E5.Sr[2] - segmented.E5.Sr[1]
# 
# newdataE6 <- data.frame(new.x = c(cp1.E6.rm, cp2.E6.rm))
# segmented.E6.Sr <- predict.segmented(fit_segmented.E6.rm, newdata = newdataE6)
# Tr1.E[6] <- segmented.E6.Sr[2] - segmented.E6.Sr[1]
# 
# newdataE7 <- data.frame(new.x = c(cp1.E7.rm, cp2.E7.rm))
# segmented.E7.Sr <- predict.segmented(fit_segmented.E7.rm, newdata = newdataE7)
# Tr1.E[7] <- segmented.E7.Sr[2] - segmented.E7.Sr[1]
# 
# newdataE8 <- data.frame(new.x = c(cp1.E8.rm, cp2.E8.rm))
# segmented.E8.Sr <- predict.segmented(fit_segmented.E8.rm, newdata = newdataE8)
# Tr1.E[8] <- segmented.E8.Sr[2] - segmented.E8.Sr[1]
# 
# newdataE9 <- data.frame(new.x = c(cp1.E9.rm, cp2.E9.rm))
# segmented.E9.Sr <- predict.segmented(fit_segmented.E9.rm, newdata = newdataE9)
# Tr1.E[9] <- segmented.E9.Sr[2] - segmented.E9.Sr[1]
# 
# 
# newdataE10 <- data.frame(new.x = c(cp1.E10.rm, cp2.E10.rm))
# segmented.E10.Sr <- predict.segmented(fit_segmented.E10.rm, newdata = newdataE10)
# Tr1.E[10] <- segmented.E10.Sr[2] - segmented.E10.Sr[1]
# 
# Tr1.E

# there is a higher fraction of UT Sr in E10 than in E9. Negligable in E8 