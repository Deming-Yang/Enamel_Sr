library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(mcp)
library(EnvCpt)
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
fit_lm.E1.rm = lm(avg ~ 1 + y, data = proc.Enamel1.rm)  # intercept-only model
fit_segmented.E1.rm = segmented(fit_lm.E1.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E1.rm)

cp.E1.rm <- fit_segmented.E1.rm$psi[3] #change point estimates, can use this to get x for plotting
cp.E1.rm.err <- fit_segmented.E1.rm$psi[5]

#extract slopes and intercepts
E1.rm.sl1 <- fit_segmented.E1.rm$coefficients[2]

E1.rm.sl2 <- fit_segmented.E1.rm$coefficients[2] + fit_segmented.E1.rm$coefficients[3] #second slope

E1.rm.sl3 <- fit_segmented.E1.rm$coefficients[2]+ fit_segmented.E1.rm$coefficients[3]+ fit_segmented.E1.rm$coefficients[4] #third slope

############Enamel 2###################
fit_lm.E2.rm = lm(avg ~ 1 + y, data = proc.Enamel2.rm)  # intercept-only model
fit_segmented.E2.rm = segmented(fit_lm.E2.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E2.rm)

cp.E2.rm <- fit_segmented.E2.rm$psi[3] #change point estimates, can use this to get x for plotting
cp.E2.rm.err <- fit_segmented.E2.rm$psi[5]

#extract slopes and intercepts
E2.rm.sl1 <- fit_segmented.E2.rm$coefficients[2]

E2.rm.sl2 <- fit_segmented.E2.rm$coefficients[2] + fit_segmented.E2.rm$coefficients[3] #second slope

E2.rm.sl3 <- fit_segmented.E2.rm$coefficients[2]+ fit_segmented.E2.rm$coefficients[3]+ fit_segmented.E2.rm$coefficients[4] #third slope

############Enamel 3###################
fit_lm.E3.rm = lm(avg ~ 1 + y, data = proc.Enamel3.rm)  # intercept-only model
fit_segmented.E3.rm = segmented(fit_lm.E3.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E3.rm)

cp.E3.rm <- fit_segmented.E3.rm$psi[3] #change point estimates, can use this to get x for plotting
cp.E3.rm.err <- fit_segmented.E3.rm$psi[5]

#extract slopes and intercepts
E3.rm.sl1 <- fit_segmented.E3.rm$coefficients[2]

E3.rm.sl2 <- fit_segmented.E3.rm$coefficients[2] + fit_segmented.E3.rm$coefficients[3] #second slope

E3.rm.sl3 <- fit_segmented.E3.rm$coefficients[2]+ fit_segmented.E3.rm$coefficients[3]+ fit_segmented.E3.rm$coefficients[4] #third slope

############Enamel 4###################
fit_lm.E4.rm = lm(avg ~ 1 + y, data = proc.Enamel4.rm)  # intercept-only model
fit_segmented.E4.rm = segmented(fit_lm.E4.rm, seg.Z = ~y, npsi = 4)  # four change points along y to make sense

summary(fit_segmented.E4.rm)

fit_segmented.E4.rm$psi

cp.E4.rm <- fit_segmented.E4.rm$psi[5] #change point estimates, can use this to get x for plotting
cp.E4.rm.err <- fit_segmented.E4.rm$psi[9]

#extract slopes and intercepts
E4.rm.sl1 <- fit_segmented.E4.rm$coefficients[2]

E4.rm.sl2 <- fit_segmented.E4.rm$coefficients[2] + fit_segmented.E4.rm$coefficients[3] #second slope

E4.rm.sl3 <- fit_segmented.E4.rm$coefficients[2]+ fit_segmented.E4.rm$coefficients[3]+ fit_segmented.E4.rm$coefficients[4] #third slope

############Enamel 5###################
fit_lm.E5.rm = lm(avg ~ 1 + y, data = proc.Enamel5.rm)  # intercept-only model
fit_segmented.E5.rm = segmented(fit_lm.E5.rm, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E5.rm)

cp.E5.rm <- fit_segmented.E5.rm$psi[3] #change point estimates, can use this to get x for plotting
cp.E5.rm.err <- fit_segmented.E5.rm$psi[5]

#extract slopes and intercepts
E5.rm.sl1 <- fit_segmented.E5.rm$coefficients[2]

E5.rm.sl2 <- fit_segmented.E5.rm$coefficients[2] + fit_segmented.E5.rm$coefficients[3] #second slope

E5.rm.sl3 <- fit_segmented.E5.rm$coefficients[2]+ fit_segmented.E5.rm$coefficients[3]+ fit_segmented.E5.rm$coefficients[4] #third slope

############Enamel 6###################
proc.Enamel6.rm.c <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm)
proc.Enamel6.rm.f <- filter(proc.Enamel6.rm.c, proc.Enamel6.rm.c$avg > 0.703)

fit_lm.E6.rm = lm(avg ~ 1 + y, data = proc.Enamel6.rm.f)  # intercept-only model
fit_segmented.E6.rm = segmented(fit_lm.E6.rm, seg.Z = ~y, npsi = 1)  # one change points along y

summary(fit_segmented.E6.rm)

cp.E6.rm <- fit_segmented.E6.rm$psi[2] #change point estimates, can use this to get x for plotting
cp.E6.rm.err <- fit_segmented.E6.rm$psi[3]

#extract slopes and intercepts
E6.rm.sl1 <- fit_segmented.E6.rm$coefficients[2]

E6.rm.sl2 <- fit_segmented.E6.rm$coefficients[2] + fit_segmented.E6.rm$coefficients[3] #second slope

############Enamel 7###################
proc.Enamel7.rm.f <- filter(proc.Enamel7.rm, proc.Enamel7.rm$avg > 0.703)
fit_lm.E7.rm = lm(avg ~ 1 + y, data = proc.Enamel7.rm.f)  # intercept-only model
fit_segmented.E7.rm = segmented(fit_lm.E7.rm, seg.Z = ~y, npsi = 1)  # one change point to make it work

summary(fit_segmented.E7.rm)

cp.E7.rm <- fit_segmented.E7.rm$psi[2] #change point estimates, can use this to get x for plotting
cp.E7.rm.err <- fit_segmented.E7.rm$psi[3]

#extract slopes and intercepts
E7.rm.sl1 <- fit_segmented.E7.rm$coefficients[2]

E7.rm.sl2 <- fit_segmented.E7.rm$coefficients[2] + fit_segmented.E7.rm$coefficients[3] #second slope


D.rm.sl2
E1.rm.sl2 #enamel is is about the same as dentine for slope 2
E2.rm.sl2
E3.rm.sl2
E4.rm.sl2
E5.rm.sl2
E6.rm.sl2
E7.rm.sl2

D.rm.sl3 #dentine is steeper than enamel for slope 3
E1.rm.sl3
E2.rm.sl3
E3.rm.sl3
E4.rm.sl3
E5.rm.sl3

#################compiling first change points within enamel#################
cp.E.rm.y <- c(cp.E1.rm,
      cp.E2.rm,
      cp.E3.rm,
      cp.E4.rm,
      cp.E5.rm,
      cp.E6.rm,
      cp.E7.rm)

cp.E.rm.y.err <- c(cp.E1.rm.err,
             cp.E2.rm.err,
             cp.E3.rm.err,
             cp.E4.rm.err,
             cp.E5.rm.err,
             cp.E6.rm.err,
             cp.E7.rm.err)
#get x values from enamel transects
cp.D.rm.x <- approx(x = dent.rm.f$y, y = dent.rm.f$x, xout = cp.D.rm[1])$y

cp.E.rm.x <- rep(0,7) #initialize vector
cp.E.rm.x[1] <- approx(x = proc.Enamel1.rm$y, y = proc.Enamel1.rm$x, xout = cp.E.rm.y[1])$y
cp.E.rm.x[2] <- approx(x = proc.Enamel2.rm$y, y = proc.Enamel2.rm$x, xout = cp.E.rm.y[2])$y
cp.E.rm.x[3] <- approx(x = proc.Enamel3.rm$y, y = proc.Enamel3.rm$x, xout = cp.E.rm.y[3])$y
cp.E.rm.x[4] <- approx(x = proc.Enamel4.rm$y, y = proc.Enamel4.rm$x, xout = cp.E.rm.y[4])$y
cp.E.rm.x[5] <- approx(x = proc.Enamel5.rm$y, y = proc.Enamel5.rm$x, xout = cp.E.rm.y[5])$y
cp.E.rm.x[6] <- approx(x = proc.Enamel6.rm$y, y = proc.Enamel6.rm$x, xout = cp.E.rm.y[6])$y
cp.E.rm.x[7] <- approx(x = proc.Enamel7.rm$y, y = proc.Enamel7.rm$x, xout = cp.E.rm.y[7])$y


###########calculate shifts in x between enamel 2 3 and 4
head(proc.Enamel2.rm)
head(proc.Enamel3.rm)
head(proc.Enamel4.rm)
#300 micron shifts in x, ~5000 micron shifts in y
atan(3/50)/pi*180 #~3.43 degrees appositional angle, exactly the same as Uno 2020
