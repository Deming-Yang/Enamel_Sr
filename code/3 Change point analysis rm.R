library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(segmented)

source("code/1 Helper functions.R")

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant molar/Enamel_Sr")

############Dentine combined###################
dent.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)
dent.rm.f <- filter(dent.rm, dent.rm$avg > 0.703)

fit_lm.D.rm = lm(avg ~ 1 + y, data = dent.rm.f)  # intercept-only model
fit_segmented.D.rm = segmented(fit_lm.D.rm, seg.Z = ~y, npsi = 3)  # Three change points along y

summary(fit_segmented.D.rm)

cp.D.rm <- fit_segmented.D.rm$psi[5] #change point estimates, can use this to get x for plotting
cp.D.rm.err <- fit_segmented.D.rm$psi[8]

D.rm.sl1 <- fit_segmented.D.rm$coefficients[2] + fit_segmented.D.rm$coefficients[3] #second slope

D.rm.sl2 <- fit_segmented.D.rm$coefficients[2]+ fit_segmented.D.rm$coefficients[3]+ fit_segmented.D.rm$coefficients[4] #third slope

D.rm.sl3 <- fit_segmented.D.rm$coefficients[2]+ fit_segmented.D.rm$coefficients[3]+ 
  fit_segmented.D.rm$coefficients[4] + fit_segmented.D.rm$coefficients[5]#third slope
############Enamel 1###################
proc.Enamel1.rm.f <- filter(proc.Enamel1.rm, proc.Enamel1.rm$avg > 0.703)
fit_lm.E1.rm = lm(avg ~ 1 + y, data = proc.Enamel1.rm.f)  # intercept-only model
fit_segmented.E1.rm = segmented(fit_lm.E1.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E1.rm)

cp1.E1.rm <- fit_segmented.E1.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E1.rm.err <- fit_segmented.E1.rm$psi[5]

cp2.E1.rm <- fit_segmented.E1.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E1.rm.err <- fit_segmented.E1.rm$psi[6]

#extract slopes and intercepts
E1.rm.sl1 <- fit_segmented.E1.rm$coefficients[2]

E1.rm.sl2 <- fit_segmented.E1.rm$coefficients[2] + fit_segmented.E1.rm$coefficients[3] #second slope

E1.rm.sl3 <- fit_segmented.E1.rm$coefficients[2]+ fit_segmented.E1.rm$coefficients[3]+ fit_segmented.E1.rm$coefficients[4] #third slope

############Enamel 2###################
proc.Enamel2.rm.f <- filter(proc.Enamel2.rm, proc.Enamel2.rm$avg > 0.703)
fit_lm.E2.rm = lm(avg ~ 1 + y, data = proc.Enamel2.rm.f)  # intercept-only model
fit_segmented.E2.rm = segmented(fit_lm.E2.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E2.rm)

cp1.E2.rm <- fit_segmented.E2.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E2.rm.err <- fit_segmented.E2.rm$psi[5]

cp2.E2.rm <- fit_segmented.E2.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E2.rm.err <- fit_segmented.E2.rm$psi[6]

#extract slopes and intercepts
E2.rm.sl1 <- fit_segmented.E2.rm$coefficients[2]

E2.rm.sl2 <- fit_segmented.E2.rm$coefficients[2] + fit_segmented.E2.rm$coefficients[3] #second slope

E2.rm.sl3 <- fit_segmented.E2.rm$coefficients[2]+ fit_segmented.E2.rm$coefficients[3]+ fit_segmented.E2.rm$coefficients[4] #third slope

############Enamel 3###################
proc.Enamel3.rm.f <- filter(proc.Enamel3.rm, proc.Enamel3.rm$avg > 0.703)
fit_lm.E3.rm = lm(avg ~ 1 + y, data = proc.Enamel3.rm.f)  # intercept-only model
fit_segmented.E3.rm = segmented(fit_lm.E3.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E3.rm)

cp1.E3.rm <- fit_segmented.E3.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E3.rm.err <- fit_segmented.E3.rm$psi[5]

cp2.E3.rm <- fit_segmented.E3.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E3.rm.err <- fit_segmented.E3.rm$psi[6]

#extract slopes and intercepts
E3.rm.sl1 <- fit_segmented.E3.rm$coefficients[2]

E3.rm.sl2 <- fit_segmented.E3.rm$coefficients[2] + fit_segmented.E3.rm$coefficients[3] #second slope

E3.rm.sl3 <- fit_segmented.E3.rm$coefficients[2]+ fit_segmented.E3.rm$coefficients[3]+ fit_segmented.E3.rm$coefficients[4] #third slope

############Enamel 4###################
proc.Enamel4.rm.f <- filter(proc.Enamel4.rm, proc.Enamel4.rm$avg > 0.703)
fit_lm.E4.rm = lm(avg ~ 1 + y, data = proc.Enamel4.rm.f)  # intercept-only model
fit_segmented.E4.rm = segmented(fit_lm.E4.rm, seg.Z = ~y, npsi = 4)  # four change points along y to make sense

summary(fit_segmented.E4.rm)

fit_segmented.E4.rm$psi

cp1.E4.rm <- fit_segmented.E4.rm$psi[5] #change point estimates, can use this to get x for plotting
cp1.E4.rm.err <- fit_segmented.E4.rm$psi[9]

cp2.E4.rm <- fit_segmented.E4.rm$psi[6] #change point estimates, can use this to get x for plotting
cp2.E4.rm.err <- fit_segmented.E4.rm$psi[10]

#extract slopes and intercepts
E4.rm.sl1 <- fit_segmented.E4.rm$coefficients[2]

E4.rm.sl2 <- fit_segmented.E4.rm$coefficients[2] + fit_segmented.E4.rm$coefficients[3] #second slope

E4.rm.sl3 <- fit_segmented.E4.rm$coefficients[2]+ fit_segmented.E4.rm$coefficients[3]+ fit_segmented.E4.rm$coefficients[4] #third slope

############Enamel 5###################
proc.Enamel5.rm.f <- filter(proc.Enamel5.rm, proc.Enamel5.rm$avg > 0.703)
fit_lm.E5.rm = lm(avg ~ 1 + y, data = proc.Enamel5.rm.f)  # intercept-only model
fit_segmented.E5.rm = segmented(fit_lm.E5.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E5.rm)

cp1.E5.rm <- fit_segmented.E5.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E5.rm.err <- fit_segmented.E5.rm$psi[5]

cp2.E5.rm <- fit_segmented.E5.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E5.rm.err <- fit_segmented.E5.rm$psi[6]

#extract slopes and intercepts
E5.rm.sl1 <- fit_segmented.E5.rm$coefficients[2]

E5.rm.sl2 <- fit_segmented.E5.rm$coefficients[2] + fit_segmented.E5.rm$coefficients[3] #second slope

E5.rm.sl3 <- fit_segmented.E5.rm$coefficients[2]+ fit_segmented.E5.rm$coefficients[3]+ fit_segmented.E5.rm$coefficients[4] #third slope

############Enamel 6###################
proc.Enamel6.rm.c <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel6ext2.rm)
proc.Enamel6.rm.f <- filter(proc.Enamel6.rm.c, proc.Enamel6.rm.c$avg > 0.703)

fit_lm.E6.rm = lm(avg ~ 1 + y, data = proc.Enamel6.rm.f)  # intercept-only model
fit_segmented.E6.rm = segmented(fit_lm.E6.rm, seg.Z = ~y, npsi = 2)  # one change points along y

summary(fit_segmented.E6.rm)

cp1.E6.rm <- fit_segmented.E6.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E6.rm.err <- fit_segmented.E6.rm$psi[5]

cp2.E6.rm <- fit_segmented.E6.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E6.rm.err <- fit_segmented.E6.rm$psi[6]

#extract slopes and intercepts
E6.rm.sl1 <- fit_segmented.E6.rm$coefficients[2]

E6.rm.sl2 <- fit_segmented.E6.rm$coefficients[2] + fit_segmented.E6.rm$coefficients[3] #second slope

############Enamel 7###################
proc.Enamel7.rm.c <- rbind(proc.Enamel7.rm, proc.Enamel7ext.rm)
proc.Enamel7.rm.f <- filter(proc.Enamel7.rm.c, proc.Enamel7.rm.c$avg > 0.703)
fit_lm.E7.rm = lm(avg ~ 1 + y, data = proc.Enamel7.rm.f)  # intercept-only model
fit_segmented.E7.rm = segmented(fit_lm.E7.rm, seg.Z = ~y, npsi = 2)  # one change point to make it work

summary(fit_segmented.E7.rm)

cp1.E7.rm <- fit_segmented.E7.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E7.rm.err <- fit_segmented.E7.rm$psi[5] 

cp2.E7.rm <- fit_segmented.E7.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E7.rm.err <- fit_segmented.E7.rm$psi[6] 

#extract slopes and intercepts
E7.rm.sl1 <- fit_segmented.E7.rm$coefficients[2]

E7.rm.sl2 <- fit_segmented.E7.rm$coefficients[2] + fit_segmented.E7.rm$coefficients[3] #second slope

############Enamel 8###################
proc.Enamel8.rm.c <- rbind(proc.Enamel8.rm, proc.Enamel8ext.rm)
proc.Enamel8.rm.f <- filter(proc.Enamel8.rm.c, proc.Enamel8.rm.c$avg > 0.703)
fit_lm.E8.rm = lm(avg ~ 1 + y, data = proc.Enamel8.rm.f)  # intercept-only model
fit_segmented.E8.rm = segmented(fit_lm.E8.rm, seg.Z = ~y, npsi = 2)  # one change point to make it work

summary(fit_segmented.E8.rm)

cp1.E8.rm <- fit_segmented.E8.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E8.rm.err <- fit_segmented.E8.rm$psi[5]

cp2.E8.rm <- fit_segmented.E8.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E8.rm.err <- fit_segmented.E8.rm$psi[6]

#extract slopes and intercepts
E8.rm.sl1 <- fit_segmented.E8.rm$coefficients[2]

E8.rm.sl2 <- fit_segmented.E8.rm$coefficients[2] + fit_segmented.E8.rm$coefficients[3] #second slope

############Enamel 9###################
proc.Enamel9.rm.c <- rbind(proc.Enamel9.rm, proc.Enamel9ext2.rm)
proc.Enamel9.rm.f <- filter(proc.Enamel9.rm.c, proc.Enamel9.rm.c$avg > 0.703)
fit_lm.E9.rm = lm(avg ~ 1 + y, data = proc.Enamel9.rm.f)  # intercept-only model
fit_segmented.E9.rm = segmented(fit_lm.E9.rm, seg.Z = ~y, npsi = 2)  # one change point to make it work

summary(fit_segmented.E9.rm)

cp1.E9.rm <- fit_segmented.E9.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E9.rm.err <- fit_segmented.E9.rm$psi[5]

cp2.E9.rm <- fit_segmented.E9.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E9.rm.err <- fit_segmented.E9.rm$psi[6]

#extract slopes and intercepts
E9.rm.sl1 <- fit_segmented.E9.rm$coefficients[2]

E9.rm.sl2 <- fit_segmented.E9.rm$coefficients[2] + fit_segmented.E9.rm$coefficients[3] #second slope

############Enamel 10###################
proc.Enamel10.rm.f <- filter(proc.Enamel10.rm, proc.Enamel10.rm$avg > 0.703)
fit_lm.E10.rm = lm(avg ~ 1 + y, data = proc.Enamel10.rm.f)  # intercept-only model
fit_segmented.E10.rm = segmented(fit_lm.E10.rm, seg.Z = ~y, npsi = 2)  # one change point to make it work

summary(fit_segmented.E10.rm)

cp1.E10.rm <- fit_segmented.E10.rm$psi[3] #change point estimates, can use this to get x for plotting
cp1.E10.rm.err <- fit_segmented.E10.rm$psi[5]

cp2.E10.rm <- fit_segmented.E10.rm$psi[4] #change point estimates, can use this to get x for plotting
cp2.E10.rm.err <- fit_segmented.E10.rm$psi[6]

#extract slopes and intercepts
E10.rm.sl1 <- fit_segmented.E10.rm$coefficients[2]

E10.rm.sl2 <- fit_segmented.E10.rm$coefficients[2] + fit_segmented.E10.rm$coefficients[3] #second slope


D.rm.sl2
E1.rm.sl2 #enamel is is about the same as dentine for slope 2
E2.rm.sl2
E3.rm.sl2
E4.rm.sl2
E5.rm.sl2
E6.rm.sl2
E7.rm.sl2

#compile slopes of the first transition
Sr.trans.sl1 <- c(D.rm.sl2, E1.rm.sl2, E2.rm.sl2, E3.rm.sl2, E4.rm.sl2, E5.rm.sl2, E6.rm.sl2, E7.rm.sl2)

#compile slopes of the second transition
D.rm.sl3 #dentine is steeper than enamel for slope 3
E1.rm.sl3
E2.rm.sl3
E3.rm.sl3
E4.rm.sl3
E5.rm.sl3

Sr.trans.sl2 <- c(D.rm.sl3, E1.rm.sl3, E2.rm.sl3, E3.rm.sl3, E4.rm.sl3, E5.rm.sl3)

#################compiling first change points within enamel#################
cp1.E.rm.y <- c(cp1.E1.rm,
      cp1.E2.rm,
      cp1.E3.rm,
      cp1.E4.rm,
      cp1.E5.rm,
      cp1.E6.rm,
      cp1.E7.rm,
      cp1.E8.rm,
      cp1.E9.rm,
      cp1.E10.rm)

cp1.E.rm.y.err <- c(cp1.E1.rm.err,
             cp1.E2.rm.err,
             cp1.E3.rm.err,
             cp1.E4.rm.err,
             cp1.E5.rm.err,
             cp1.E6.rm.err,
             cp1.E7.rm.err,
             cp1.E8.rm.err,
             cp1.E9.rm.err,
             cp1.E10.rm.err)


cp2.E.rm.y <- c(cp2.E1.rm,
                cp2.E2.rm,
                cp2.E3.rm,
                cp2.E4.rm,
                cp2.E5.rm,
                cp2.E6.rm,
                cp2.E7.rm,
                cp2.E8.rm,
                cp2.E9.rm,
                cp2.E10.rm)

cp2.E.rm.y.err <- c(cp2.E1.rm.err,
                    cp2.E2.rm.err,
                    cp2.E3.rm.err,
                    cp2.E4.rm.err,
                    cp2.E5.rm.err,
                    cp2.E6.rm.err,
                    cp2.E7.rm.err,
                    cp2.E8.rm.err,
                    cp2.E9.rm.err,
                    cp2.E10.rm.err)


#get x values from enamel transects
cp.D.rm.x <- approx(x = dent.rm.f$y, y = dent.rm.f$x, xout = cp.D.rm[1])$y

cp1.E.rm.x <- rep(0,10) #initialize vector
cp1.E.rm.x[1] <- approx(x = proc.Enamel1.rm.f$y, y = proc.Enamel1.rm.f$x, xout = cp1.E.rm.y[1])$y
cp1.E.rm.x[2] <- approx(x = proc.Enamel2.rm.f$y, y = proc.Enamel2.rm.f$x, xout = cp1.E.rm.y[2])$y
cp1.E.rm.x[3] <- approx(x = proc.Enamel3.rm.f$y, y = proc.Enamel3.rm.f$x, xout = cp1.E.rm.y[3])$y
cp1.E.rm.x[4] <- approx(x = proc.Enamel4.rm.f$y, y = proc.Enamel4.rm.f$x, xout = cp1.E.rm.y[4])$y
cp1.E.rm.x[5] <- approx(x = proc.Enamel5.rm.f$y, y = proc.Enamel5.rm.f$x, xout = cp1.E.rm.y[5])$y
cp1.E.rm.x[6] <- approx(x = proc.Enamel6.rm.f$y, y = proc.Enamel6.rm.f$x, xout = cp1.E.rm.y[6])$y
cp1.E.rm.x[7] <- approx(x = proc.Enamel7.rm.f$y, y = proc.Enamel7.rm.f$x, xout = cp1.E.rm.y[7])$y
cp1.E.rm.x[8] <- approx(x = proc.Enamel8.rm.f$y, y = proc.Enamel8.rm.f$x, xout = cp1.E.rm.y[8])$y
cp1.E.rm.x[9] <- approx(x = proc.Enamel9.rm.f$y, y = proc.Enamel9.rm.f$x, xout = cp1.E.rm.y[9])$y
cp1.E.rm.x[10] <- approx(x = proc.Enamel10.rm.f$y, y = proc.Enamel10.rm.f$x, xout = cp1.E.rm.y[10])$y

cp2.E.rm.x <- rep(0,10) #initialize vector
cp2.E.rm.x[1] <- approx(x = proc.Enamel1.rm.f$y, y = proc.Enamel1.rm.f$x, xout = cp2.E.rm.y[1])$y
cp2.E.rm.x[2] <- approx(x = proc.Enamel2.rm.f$y, y = proc.Enamel2.rm.f$x, xout = cp2.E.rm.y[2])$y
cp2.E.rm.x[3] <- approx(x = proc.Enamel3.rm.f$y, y = proc.Enamel3.rm.f$x, xout = cp2.E.rm.y[3])$y
cp2.E.rm.x[4] <- approx(x = proc.Enamel4.rm.f$y, y = proc.Enamel4.rm.f$x, xout = cp2.E.rm.y[4])$y
cp2.E.rm.x[5] <- approx(x = proc.Enamel5.rm.f$y, y = proc.Enamel5.rm.f$x, xout = cp2.E.rm.y[5])$y
cp2.E.rm.x[6] <- approx(x = proc.Enamel6.rm.f$y, y = proc.Enamel6.rm.f$x, xout = cp2.E.rm.y[6])$y
cp2.E.rm.x[7] <- approx(x = proc.Enamel7.rm.f$y, y = proc.Enamel7.rm.f$x, xout = cp2.E.rm.y[7])$y
cp2.E.rm.x[8] <- approx(x = proc.Enamel8.rm.f$y, y = proc.Enamel8.rm.f$x, xout = cp2.E.rm.y[8])$y
cp2.E.rm.x[9] <- approx(x = proc.Enamel9.rm.f$y, y = proc.Enamel9.rm.f$x, xout = cp2.E.rm.y[9])$y
cp2.E.rm.x[10] <- approx(x = proc.Enamel10.rm.f$y, y = proc.Enamel10.rm.f$x, xout = cp2.E.rm.y[10])$y


###########calculate shifts in x between enamel 2 3 and 4
head(proc.Enamel2.rm)
head(proc.Enamel3.rm)
head(proc.Enamel4.rm)
#300 micron shifts in x, ~5000 micron shifts in y
atan(3/50)/pi*180 #~3.43 degrees appositional angle, exactly the same as Uno 2020
