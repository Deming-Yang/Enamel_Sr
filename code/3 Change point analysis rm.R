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

#####################transform all coordinates to evaluate Sr slopes###########
dent.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)
#first transform dentine transect
dent.rm.s <- dent.rm[order(dent.rm$y),]

dent.rm.s.diff.x <- diff(dent.rm.s$x)
dent.rm.s.diff.y <- diff(dent.rm.s$y)

#calculate distance between adjacent points
dent.rm.s.dist <- sqrt(dent.rm.s.diff.x^2 + dent.rm.s.diff.y^2)

dent.rm.new.x <- c(0, cumsum(dent.rm.s.dist))

dent.rm.new.y <- rep(-150, length(dent.rm.new.x)) #dentine is -150 microns away from EDJ

dent.rm.proj <- tibble(dent.rm, data.frame(new.x = dent.rm.new.x, new.y = dent.rm.new.y))

#######Enamel 1######
proc.Enamel1.rm.s <- proc.Enamel1.rm[order(proc.Enamel1.rm$y),]

E1.proj <- proj.transect(proc.Enamel1.rm.s, dent.rm.s)

proc.Enamel1.rm.proj <- tibble(proc.Enamel1.rm.s, E1.proj)

#######Enamel 2######
proc.Enamel2.rm.s <- proc.Enamel2.rm[order(proc.Enamel2.rm$y),]

E2.proj <- proj.transect(proc.Enamel2.rm.s, dent.rm.s)

proc.Enamel2.rm.proj <- tibble(proc.Enamel2.rm.s, E2.proj)

#######Enamel 3######
proc.Enamel3.rm.s <- proc.Enamel3.rm[order(proc.Enamel3.rm$y),]

E3.proj <- proj.transect(proc.Enamel3.rm.s, dent.rm.s)

proc.Enamel3.rm.proj <- tibble(proc.Enamel3.rm.s, E3.proj)

#######Enamel 4######
proc.Enamel4.rm.s <- proc.Enamel4.rm[order(proc.Enamel4.rm$y),]

E4.proj <- proj.transect(proc.Enamel4.rm.s, dent.rm.s)

proc.Enamel4.rm.proj <- tibble(proc.Enamel4.rm.s, E4.proj)

#######Enamel 5######
proc.Enamel5.rm.s <- proc.Enamel5.rm[order(proc.Enamel5.rm$y),]

E5.proj <- proj.transect(proc.Enamel5.rm.s, dent.rm.s)

proc.Enamel5.rm.proj <- tibble(proc.Enamel5.rm.s, E5.proj)

#######Enamel 6######
proc.Enamel6.rm <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel6ext2.rm)

proc.Enamel6.rm.s <- proc.Enamel6.rm[order(proc.Enamel6.rm$y),]

E6.proj <- proj.transect(proc.Enamel6.rm.s, dent.rm.s)

proc.Enamel6.rm.proj <- tibble(proc.Enamel6.rm.s, E6.proj)

#######Enamel 7######
proc.Enamel7.rm <- rbind(proc.Enamel7.rm, proc.Enamel7ext.rm)

proc.Enamel7.rm.s <- proc.Enamel7.rm[order(proc.Enamel7.rm$y),]

E7.proj <- proj.transect(proc.Enamel7.rm.s, dent.rm.s)

proc.Enamel7.rm.proj <- tibble(proc.Enamel7.rm.s, E7.proj)

#######Enamel 8######
proc.Enamel8.rm <- rbind(proc.Enamel8.rm, proc.Enamel8ext.rm)

proc.Enamel8.rm.s <- proc.Enamel8.rm[order(proc.Enamel8.rm$y),]

E8.proj <- proj.transect(proc.Enamel8.rm.s, dent.rm.s)

proc.Enamel8.rm.proj <- tibble(proc.Enamel8.rm.s, E8.proj)

#######Enamel 9######
proc.Enamel9.rm <- rbind(proc.Enamel9.rm, proc.Enamel9ext2.rm)

proc.Enamel9.rm.s <- proc.Enamel9.rm[order(proc.Enamel9.rm$y),]

E9.proj <- proj.transect(proc.Enamel9.rm.s, dent.rm.s)

proc.Enamel9.rm.proj <- tibble(proc.Enamel9.rm.s, E9.proj)

#######Enamel 10######
proc.Enamel10.rm.s <- proc.Enamel10.rm[order(proc.Enamel10.rm$y),]

E10.proj <- proj.transect(proc.Enamel10.rm.s, dent.rm.s)

proc.Enamel10.rm.proj <- tibble(proc.Enamel10.rm.s, E10.proj)

# #######Enamel 1######
# E1.proj <- proj.transect(proc.Enamel1.rm$x, proc.Enamel1.rm$y, dent.rm)
# 
# proc.Enamel1.rm.proj <- tibble(proc.Enamel1.rm, E1.proj)
# 
# #######Enamel 2######
# E2.proj <- proj.transect(proc.Enamel2.rm$x, proc.Enamel2.rm$y, dent.rm)
# 
# proc.Enamel2.rm.proj <- tibble(proc.Enamel2.rm, E2.proj)
# 
# #######Enamel 3######
# E3.proj <- proj.transect(proc.Enamel3.rm$x, proc.Enamel3.rm$y, dent.rm)
# 
# proc.Enamel3.rm.proj <- tibble(proc.Enamel3.rm, E3.proj)
# 
# #######Enamel 4######
# E4.proj <- proj.transect(proc.Enamel4.rm$x, proc.Enamel4.rm$y, dent.rm)
# 
# proc.Enamel4.rm.proj <- tibble(proc.Enamel4.rm, E4.proj)
# 
# #######Enamel 5######
# E5.proj <- proj.transect(proc.Enamel5.rm$x, proc.Enamel5.rm$y, dent.rm)
# 
# proc.Enamel5.rm.proj <- tibble(proc.Enamel5.rm, E5.proj)
# 
# #######Enamel 6######
# proc.Enamel6.rm <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel6ext2.rm)
# 
# proc.Enamel6.rm.s <- proc.Enamel6.rm[order(proc.Enamel6.rm$y),]
# 
# E6.proj <- proj.transect(proc.Enamel6.rm.s$x, proc.Enamel6.rm.s$y, dent.rm)
# 
# proc.Enamel6.rm.proj <- tibble(proc.Enamel6.rm.s, E6.proj)
# 
# #######Enamel 7######
# proc.Enamel7.rm <- rbind(proc.Enamel7.rm, proc.Enamel7ext.rm)
# 
# proc.Enamel7.rm.s <- proc.Enamel7.rm[order(proc.Enamel7.rm$y),]
# 
# E7.proj <- proj.transect(proc.Enamel7.rm.s$x, proc.Enamel7.rm.s$y, dent.rm)
# 
# proc.Enamel7.rm.proj <- tibble(proc.Enamel7.rm.s, E7.proj)
# 
# #######Enamel 8######
# proc.Enamel8.rm <- rbind(proc.Enamel8.rm, proc.Enamel8ext.rm)
# 
# proc.Enamel8.rm.s <- proc.Enamel8.rm[order(proc.Enamel8.rm$y),]
# 
# E8.proj <- proj.transect(proc.Enamel8.rm.s$x, proc.Enamel8.rm.s$y, dent.rm)
# 
# proc.Enamel8.rm.proj <- tibble(proc.Enamel8.rm.s, E8.proj)
# 
# #######Enamel 9######
# proc.Enamel9.rm <- rbind(proc.Enamel9.rm, proc.Enamel9ext2.rm)
# 
# proc.Enamel9.rm.s <- proc.Enamel9.rm[order(proc.Enamel9.rm$y),]
# 
# E9.proj <- proj.transect(proc.Enamel9.rm.s$x, proc.Enamel9.rm.s$y, dent.rm)
# 
# proc.Enamel9.rm.proj <- tibble(proc.Enamel9.rm.s, E9.proj)
# 
# #######Enamel 10######
# E10.proj <- proj.transect(proc.Enamel10.rm$x, proc.Enamel10.rm$y, dent.rm)
# 
# proc.Enamel10.rm.proj <- tibble(proc.Enamel10.rm, E10.proj)


######################## change point analysis ########################

############Dentine combined###################

dent.rm.f <- filter(dent.rm.proj, dent.rm.proj$avg > 0.703)

fit_lm.D.rm = lm(avg ~ 1 + new.x, data = dent.rm.f)  # intercept-only model
set.seed(1234)
fit_segmented.D.rm = segmented(fit_lm.D.rm, seg.Z = ~new.x, npsi = 3)  # Three change points along new.x

summary(fit_segmented.D.rm)

cp.D.rm <- fit_segmented.D.rm$psi[5] #change point estimates, can use this to get x for plotting
cp.D.rm.err <- fit_segmented.D.rm$psi[8]

D.rm.sl1 <- fit_segmented.D.rm$coefficients[2] + fit_segmented.D.rm$coefficients[3] #second slope

D.rm.sl2 <- fit_segmented.D.rm$coefficients[2]+ fit_segmented.D.rm$coefficients[3]+ fit_segmented.D.rm$coefficients[4] #third slope

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

E4.rm.sl2 <- fit_segmented.E4.rm$coefficients[2]+ fit_segmented.E4.rm$coefficients[3]+ fit_segmented.E4.rm$coefficients[4] #third slope

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

D.rm.sl2
E1.rm.sl2 #enamel is is about the same as dentine for slope 2
E2.rm.sl2
E3.rm.sl2
E4.rm.sl2
E5.rm.sl2
E6.rm.sl2
E7.rm.sl2

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

par(mfrow=c(1,2))
plot(cp1.E.rm.x, cp1.E.rm.y, main = "Change point 1")
plot(cp2.E.rm.x, cp2.E.rm.y, main = "Change point 2")

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
#~3.3 degrees, which is identical to the appositional angle measured by Uno 2012

#use mean of the two angles in the forward model
fwd.appo.sl <- mean(c(lm.cp1.E.proj.inner$coefficients[2],lm.cp2.E.proj.inner$coefficients[2]))


