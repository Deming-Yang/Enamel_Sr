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

##############try mcp package################
#model structure really matters!
model = list(
  corr.87Sr.86Sr ~ 1,  # Plateau in the first segment (int_1)
   ~ y,    # Joined slope (time_2) at cp_1
   ~ y      # new intercept (int_3) at cp_2
)

proc.Enamel1_mcp = mcp(model, data = proc.Enamel1, par_x = "y")

summary(proc.Enamel1_mcp)

plot(proc.Enamel1_mcp)

plot_pars(proc.Enamel1_mcp)


############try segmented package################
############Dentine combined###################
dent <- rbind(proc.EDJ.dentin1, proc.EDJ.dentin2)
dent.f <- filter(dent, dent$corr.87Sr.86Sr > 0.703)

fit_lm.D = lm(corr.87Sr.86Sr ~ 1 + y, data = dent.f)  # intercept-only model
fit_segmented.D = segmented(fit_lm.D, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.D)

#preliminary plot
plot(fit_segmented.D)
lines.segmented(fit_segmented.D)
points.segmented(fit_segmented.D)
points(dent.f$y, dent.f$corr.87Sr.86Sr)

#extract slopes and intercepts
D.sl1 <- fit_segmented.D$coefficients[2]

D.sl2 <- fit_segmented.D$coefficients[2] + fit_segmented.D$coefficients[3] #second slope

D.sl3 <- fit_segmented.D$coefficients[2]+ fit_segmented.D$coefficients[3]+ fit_segmented.D$coefficients[4] #third slope


############Enamel 1###################
fit_lm.E1 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel1)  # intercept-only model
fit_segmented.E1 = segmented(fit_lm.E1, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E1)

#preliminary plot
plot(fit_segmented.E1)
lines.segmented(fit_segmented.E1)
points.segmented(fit_segmented.E1)
points(proc.Enamel1$y, proc.Enamel1$corr.87Sr.86Sr)

#extract slopes and intercepts
E1.sl1 <- fit_segmented.E1$coefficients[2]

E1.sl2 <- fit_segmented.E1$coefficients[2] + fit_segmented.E1$coefficients[3] #second slope

E1.sl3 <- fit_segmented.E1$coefficients[2]+ fit_segmented.E1$coefficients[3]+ fit_segmented.E1$coefficients[4] #third slope

############Enamel 2###################
fit_lm.E2 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel2)  # intercept-only model
fit_segmented.E2 = segmented(fit_lm.E2, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E2)

#preliminary plot
plot(fit_segmented.E2)
lines.segmented(fit_segmented.E2)
points.segmented(fit_segmented.E2)
points(proc.Enamel2$y, proc.Enamel2$corr.87Sr.86Sr)

#extract slopes and intercepts
E2.sl1 <- fit_segmented.E2$coefficients[2]

E2.sl2 <- fit_segmented.E2$coefficients[2] + fit_segmented.E2$coefficients[3] #second slope

E2.sl3 <- fit_segmented.E2$coefficients[2]+ fit_segmented.E2$coefficients[3]+ fit_segmented.E2$coefficients[4] #third slope

############Enamel 3###################
fit_lm.E3 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel3)  # intercept-only model
fit_segmented.E3 = segmented(fit_lm.E3, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E3)

#preliminary plot
plot(fit_segmented.E3)
lines.segmented(fit_segmented.E3)
points.segmented(fit_segmented.E3)
points(proc.Enamel3$y, proc.Enamel3$corr.87Sr.86Sr)

#extract slopes and intercepts
E3.sl1 <- fit_segmented.E3$coefficients[2]

E3.sl2 <- fit_segmented.E3$coefficients[2] + fit_segmented.E3$coefficients[3] #second slope

E3.sl3 <- fit_segmented.E3$coefficients[2]+ fit_segmented.E3$coefficients[3]+ fit_segmented.E3$coefficients[4] #third slope

############Enamel 4###################
fit_lm.E4 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel4)  # intercept-only model
fit_segmented.E4 = segmented(fit_lm.E4, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E4)

#preliminary plot
plot(fit_segmented.E4)
lines.segmented(fit_segmented.E4)
points.segmented(fit_segmented.E4)
points(proc.Enamel4$y, proc.Enamel4$corr.87Sr.86Sr)

#extract slopes and intercepts
E4.sl1 <- fit_segmented.E4$coefficients[2]

E4.sl2 <- fit_segmented.E4$coefficients[2] + fit_segmented.E4$coefficients[3] #second slope

E4.sl3 <- fit_segmented.E4$coefficients[2]+ fit_segmented.E4$coefficients[3]+ fit_segmented.E4$coefficients[4] #third slope

############Enamel 5###################
fit_lm.E5 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel5)  # intercept-only model
fit_segmented.E5 = segmented(fit_lm.E5, seg.Z = ~y, npsi = 2)  # Two change points along y

summary(fit_segmented.E5)

#preliminary plot
plot(fit_segmented.E5)
lines.segmented(fit_segmented.E5)
points.segmented(fit_segmented.E5)
points(proc.Enamel5$y, proc.Enamel5$corr.87Sr.86Sr)

#extract slopes and intercepts
E5.sl1 <- fit_segmented.E5$coefficients[2]

E5.sl2 <- fit_segmented.E5$coefficients[2] + fit_segmented.E5$coefficients[3] #second slope

E5.sl3 <- fit_segmented.E5$coefficients[2]+ fit_segmented.E5$coefficients[3]+ fit_segmented.E5$coefficients[4] #third slope

############Enamel 6###################
proc.Enamel6.c <- rbind(proc.Enamel6, proc.Enamel6ext)
proc.Enamel6.f <- filter(proc.Enamel6.c, proc.Enamel6.c$corr.87Sr.86Sr > 0.703)

fit_lm.E6 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel6.f)  # intercept-only model
fit_segmented.E6 = segmented(fit_lm.E6, seg.Z = ~y, npsi = 1)  # one change points along y

summary(fit_segmented.E6)

#preliminary plot
plot(fit_segmented.E6)
lines.segmented(fit_segmented.E6)
points.segmented(fit_segmented.E6)
points(proc.Enamel6.f$y, proc.Enamel6.f$corr.87Sr.86Sr)

#extract slopes and intercepts
E6.sl1 <- fit_segmented.E6$coefficients[2]

E6.sl2 <- fit_segmented.E6$coefficients[2] + fit_segmented.E6$coefficients[3] #second slope

############Enamel 7###################
proc.Enamel7.f <- filter(proc.Enamel7, proc.Enamel7$corr.87Sr.86Sr > 0.703)
fit_lm.E7 = lm(corr.87Sr.86Sr ~ 1 + y, data = proc.Enamel7.f)  # intercept-only model
fit_segmented.E7 = segmented(fit_lm.E7, seg.Z = ~y, npsi = 1)  # one change point to make it work

summary(fit_segmented.E7)

#preliminary plot
plot(fit_segmented.E7)
lines.segmented(fit_segmented.E7)
points.segmented(fit_segmented.E7)
points(proc.Enamel7.f$y, proc.Enamel7.f$corr.87Sr.86Sr)

#extract slopes and intercepts
E7.sl1 <- fit_segmented.E7$coefficients[2]

E7.sl2 <- fit_segmented.E7$coefficients[2] + fit_segmented.E7$coefficients[3] #second slope


D.sl2
E1.sl2 #enamel is is about the same as dentine for slope 2
E2.sl2
E3.sl2
E5.sl2
E6.sl2
E7.sl2

D.sl3 #dentine is steeper than enamel for slope 3
E1.sl3
E2.sl3
E3.sl3
E5.sl3

fit_segmented.D$psi #change point estimates, can use this to get x for plotting
fit_segmented.E1$psi #change point estimates, can use this to get x for plotting
fit_segmented.E2$psi #change point estimates, can use this to get x for plotting
fit_segmented.E3$psi #change point estimates, can use this to get x for plotting
fit_segmented.E4$psi #change point estimates, can use this to get x for plotting
fit_segmented.E5$psi
fit_segmented.E6$psi
fit_segmented.E7$psi
