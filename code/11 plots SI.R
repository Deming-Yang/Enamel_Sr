library(scales)
library(lattice)
library(zoo)
library(changepoint)
library(segmented)
library(viridisLite)

#################  Fig S2 #################  
################# LA-ICP-MS enamel-dentine comparison #################
################# overlay with shifted x values in change points ###############
par(mfrow=c(2,5))
#enamel 1, no x- shift is needed
plot(dent.rm.f2$new.x, dent.rm.f2$mov.avg, main = "Enamel 1",col= "gray24",
     type="l", lwd=2, xlim=c(4e4,8e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x, rev(dent.rm.f2$new.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En1.plot <- na.omit(proc.Enamel1.rm.f)
polygon(c(proc.En1.plot$new.x, rev(proc.En1.plot$new.x)), 
        c(proc.En1.plot$mov.avg + proc.En1.plot$sd, 
          rev(proc.En1.plot$mov.avg - proc.En1.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel1.rm.f$new.x, proc.Enamel1.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E1.rm, add = T)
lines.segmented(fit_segmented.E1.rm)
points.segmented(fit_segmented.E1.rm)

#enamel 2, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[2]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 2",col= "gray24",
     type="l", lwd=2, xlim=c(3.5e4,7.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En2.plot <- na.omit(proc.Enamel2.rm.f)
polygon(c(proc.En2.plot$new.x, rev(proc.En2.plot$new.x)), 
        c(proc.En2.plot$mov.avg + proc.En2.plot$sd, 
          rev(proc.En2.plot$mov.avg - proc.En2.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel2.rm.f$new.x, proc.Enamel2.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E2.rm, add = T)
lines.segmented(fit_segmented.E2.rm)
points.segmented(fit_segmented.E2.rm)

#enamel 3, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[3] 

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 3",col= "gray24",
     type="l", lwd=2, xlim=c(3e4,7e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En3.plot <- na.omit(proc.Enamel3.rm.f)
polygon(c(proc.En3.plot$new.x, rev(proc.En3.plot$new.x)), 
        c(proc.En3.plot$mov.avg + proc.En3.plot$sd, 
          rev(proc.En3.plot$mov.avg - proc.En3.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel3.rm.f$new.x, proc.Enamel3.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E3.rm, add = T)
lines.segmented(fit_segmented.E3.rm)
points.segmented(fit_segmented.E3.rm)

#enamel 4, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[4]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 4",col= "gray24",
     type="l", lwd=2, xlim=c(2.5e4,6.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En4.plot <- na.omit(proc.Enamel4.rm.f)
polygon(c(proc.En4.plot$new.x, rev(proc.En4.plot$new.x)), 
        c(proc.En4.plot$mov.avg + proc.En4.plot$sd, 
          rev(proc.En4.plot$mov.avg - proc.En4.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel4.rm.f$new.x, proc.Enamel4.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E4.rm, add = T)
lines.segmented(fit_segmented.E4.rm)
points.segmented(fit_segmented.E4.rm)

#enamel 5, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[5]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 5",col= "gray24",
     type="l", lwd=2, xlim=c(1.5e4,5.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En5.plot <- na.omit(proc.Enamel5.rm.f)
polygon(c(proc.En5.plot$new.x, rev(proc.En5.plot$new.x)), 
        c(proc.En5.plot$mov.avg + proc.En5.plot$sd, 
          rev(proc.En5.plot$mov.avg - proc.En5.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel5.rm.f$new.x, proc.Enamel5.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E5.rm, add = T)
lines.segmented(fit_segmented.E5.rm)
points.segmented(fit_segmented.E5.rm)

#enamel 6, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[6]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 6",col= "gray24",
     type="l", lwd=2, xlim=c(1.2e4,5.2e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En6.plot <- na.omit(proc.Enamel6.rm.f)
polygon(c(proc.En6.plot$new.x, rev(proc.En6.plot$new.x)), 
        c(proc.En6.plot$mov.avg + proc.En6.plot$sd, 
          rev(proc.En6.plot$mov.avg - proc.En6.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel6.rm.f$new.x, proc.Enamel6.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E6.rm, add = T)
lines.segmented(fit_segmented.E6.rm)
points.segmented(fit_segmented.E6.rm)

#enamel 7, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[7]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 7",col= "gray24",
     type="l", lwd=2, xlim=c(1e4,5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En7.plot <- na.omit(proc.Enamel7.rm.f)
polygon(c(proc.En7.plot$new.x, rev(proc.En7.plot$new.x)), 
        c(proc.En7.plot$mov.avg + proc.En7.plot$sd, 
          rev(proc.En7.plot$mov.avg - proc.En7.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel7.rm.f$new.x, proc.Enamel7.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E7.rm, add = T)
lines.segmented(fit_segmented.E7.rm)
points.segmented(fit_segmented.E7.rm)

#enamel 8, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[8]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 8",col= "gray24",
     type="l", lwd=2, xlim=c(0.5e4,4.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En8.plot <- na.omit(proc.Enamel8.rm.f)
polygon(c(proc.En8.plot$new.x, rev(proc.En8.plot$new.x)), 
        c(proc.En8.plot$mov.avg + proc.En8.plot$sd, 
          rev(proc.En8.plot$mov.avg - proc.En8.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel8.rm.f$new.x, proc.Enamel8.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E8.rm, add = T)
lines.segmented(fit_segmented.E8.rm)
points.segmented(fit_segmented.E8.rm)

#enamel 9, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[9]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 9",col= "gray24",
     type="l", lwd=2, xlim=c(0.5e4,4.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En9.plot <- na.omit(proc.Enamel9.rm.f)
polygon(c(proc.En9.plot$new.x, rev(proc.En9.plot$new.x)), 
        c(proc.En9.plot$mov.avg + proc.En9.plot$sd, 
          rev(proc.En9.plot$mov.avg - proc.En9.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel9.rm.f$new.x, proc.Enamel9.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E9.rm, add = T)
lines.segmented(fit_segmented.E9.rm)
points.segmented(fit_segmented.E9.rm)

#enamel 10, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[10]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg, main = "Enamel 10",col= "gray24",
     type="l", lwd=2, xlim=c(0e4,4e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$mov.avg + dent.rm.f2$sd, rev(dent.rm.f2$mov.avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$mov.avg,col= "gray24", lwd=2)

proc.En10.plot <- na.omit(proc.Enamel10.rm.f)
polygon(c(proc.En10.plot$new.x, rev(proc.En10.plot$new.x)), 
        c(proc.En10.plot$mov.avg + proc.En10.plot$sd, 
          rev(proc.En10.plot$mov.avg - proc.En10.plot$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel10.rm.f$new.x, proc.Enamel10.rm.f$mov.avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr.1, lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr.2, lty = 2, lwd = 1.5)
plot(fit_segmented.E10.rm, add = T)
lines.segmented(fit_segmented.E10.rm)
points.segmented(fit_segmented.E10.rm)

#################  Fig S3 ################# 

################# Multi-substrate comparison with x axis transformed to time axis###############

# Fig S3 A compare laser, hand drill and micromill transects using the molar plate geometry
par(mfrow=c(1,1)) #700*390
plot(-10, -10, col= alpha("lightcyan4", 0),
     pch=16, cex=1, xlim=c(100,0),ylim=c(0.704,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")

#CA and UT Sr measurements from Yang et al. 2023
abline(h = CA.Sr)

abline(h = UT.Sr)
# dentine transect
# polygon(c(dent.rm.new.x, rev(dent.rm.new.x)), 
#         c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
#         col = "gray60", border = NA)
# lines(dent.rm.new.x, dent.rm.f2$avg,col= "gray24", lwd=2)

# Enamel transect 1
polygon(c(proc.Enamel1.rm.new.x, rev(proc.Enamel1.rm.new.x)), 
        c(proc.Enamel1.rm.f$avg + proc.Enamel1.rm.f$sd, 
          rev(proc.Enamel1.rm.f$avg - proc.Enamel1.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel1.rm.new.x, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.9),lwd=2)

# Enamel transect 7
polygon(c(proc.Enamel7.rm.new.x, rev(proc.Enamel7.rm.new.x)), 
        c(proc.Enamel7.rm.f$avg + proc.Enamel7.rm.f$sd, 
          rev(proc.Enamel7.rm.f$avg - proc.Enamel7.rm.f$sd)), 
        col = alpha("orange3", 0.3), border = NA)
lines(proc.Enamel7.rm.new.x, proc.Enamel7.rm.f$avg, col= alpha("orange3", 0.9),lwd=2)

# Enamel transect 10
polygon(c(proc.Enamel10.rm.new.x, rev(proc.Enamel10.rm.new.x)), 
        c(proc.Enamel10.rm.f$avg + proc.Enamel10.rm.f$sd, 
          rev(proc.Enamel10.rm.f$avg - proc.Enamel10.rm.f$sd)), 
        col = alpha("orange4", 0.3), border = NA)
lines(proc.Enamel10.rm.new.x, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.9),lwd=2)

points(Drill.no$Dist..From.cervix, Drill.no$corr..87Sr.86Sr, pch=16, cex = 1.2, col ="red4")

# Fig S3 B
# micromill enamel on its own 
par(mfrow=c(1,1))
plot(Rm3.5b.mill.no$Dist..From.EDJ, Rm3.5b.mill.no$corr..87Sr.86Sr, col= "cyan4",
     pch=16, cex = 1.2,
     xlim=c(0,2200),ylim=c(0.704,0.712),
     xlab="Distance from EDJ (microns)",
     main = "Enamel micromill data of Rm3.5",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)


#################  Fig S4 #################  
############# Figure modeling growth rate along the length of the crown #############
# visualize the model fit
par(mfrow=c(1, 2))
plot(Rm3.5.angle$dist, Rm3.5.angle$angle, ylim = c(2,5), xlim = c(0,100),
     ylab = "Alpha (degrees)", xlab ="Distance from cervix (mm)",
     main = "Angle of enamel apposition", pch =16)
lines(ref.length.v, Pred.ang.Rm3.5)

plot(Rm3.5.angle$dist, Rm3.5.angle$tan, ylim = c(0.04,0.09), xlim = c(0,100),
     ylab = "Tangent(alpha)", xlab ="Distance from cervix (mm)", pch =16)
lines(ref.length.v, Pred.tan.Rm3.5)
############# end of growth rate plot ###########

#################  Fig S5 ################# 
# Fig S5 was made with InkScape

#################  Fig S6 #################
# comparing reconstructed timelines of three data series
# 1. tusk micromill, 2. molar dentine ICPMS, 3. Tusk dentine ICPMS
par(mfrow=c(1,1))
plot(-1000, -1, col= "gray24",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Tusk micromill vs LA-ICP-MS molar D and tusk D",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

polygon(c(dent.tl$tl, rev(dent.tl$tl)), 
        c(dent.tl$Sr + dent.tl$sd, 
          rev(dent.tl$Sr - dent.tl$sd)), 
        col = "gray60", border = NA)
lines(dent.tl$tl, dent.tl$Sr, col= "gray24", lwd=2)


# misha's tusk dentine LA-ICP-MS with the same interference correction
polygon(c(n.avg.misha.25.tl.al, rev(n.avg.misha.25.tl.al)), 
        c(n.avg.misha.25.sr + n.sd.misha.25.sr, 
          rev(n.avg.misha.25.sr - n.sd.misha.25.sr)), 
        col = alpha("orange", 0.3), border = NA)
lines(n.avg.misha.25.tl.al, n.avg.misha.25.sr, col = alpha("orange", 0.9), lwd=2)

points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

#################  Fig S7 #################

par(mfrow=c(3,3))
plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.2, window = 300",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha , misha.sim.en.1, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.2, window = 600",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.2, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.2, window = 900",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.3, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.4, window = 300",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.4, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.4, window = 600",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.5, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.4, window = 900",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.6, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.6, window = 300",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.7, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.6, window = 600",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.8, col = "red")

plot(1:t.fwd.misha, misha.sim.serum[[1]], type = "l",
     main = "f.ma = 0.6, window = 900",
     xlab = "Timeline (days)", ylab = "87Sr/86Sr")
lines(1:t.fwd.misha, misha.sim.en.9, col = "red")

#################  Fig S8 #################
# density plots showing posterior
denplot(as.mcmc(post.comb), parms = c("a","b","c",
                                      "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                                      "Rpri.mod", "Raft.mod"))

#################  Fig S8 #################
# compare estimated fraction of post-movement overprint,
# in each of the following selected data series
# part I: model-data comparison
par(mfrow=c(2,5)) # 1200 * 550

# 2 micromill tusk dentine and modeled serum, assumes no overprint
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2,  
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "Ref. data 1: micromill tusk dentine")
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[2]], lty = 2)
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[3]], lty = 2)
points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
       pch = 18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

# 3 LA-ICP MS EN 9 and mixed R
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2, col = "gray56", 
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "LA-ICP-MS En 9")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[2]], lty = 2, col = "gray48")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[3]], lty = 2, col = "gray48")
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.En9.89[[1]],lwd = 2, lty = 2)
lines(En9.tl, R.En9,
      col = "orange3", lwd = 2)

# 4 LA-ICP MS EN 10 and mixed R
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2, col = "gray56", 
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "LA-ICP-MS En 10")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[2]], lty = 2, col = "gray48")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[3]], lty = 2, col = "gray48")
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.En10.89[[1]],lwd = 2, lty = 2)
lines(En10.tl, R.En10,
      col = "orange4", lwd = 2)

# 5 hand drill and mixed R
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2, col = "gray56", 
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "Enamel conventional drill")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[2]], lty = 2, col = "gray48")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[3]], lty = 2, col = "gray48")
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.drill.89[[1]],lwd = 2, lty = 2)
points(Edrill.tl.f$tl, Edrill.tl.f$Sr, 
       pch = 16, cex = 2, col = alpha("red4", 0.7))

# 6 Rm3.5b micromill and mixed R
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2, col = "gray56", 
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "Enamel micromill")
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[2]], lty = 2, col = "gray48",)
# lines(bin.thin*1:t - 400, post.comb.R1m.89[[3]], lty = 2, col = "gray48",)
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.Rm3.5b.89[[1]],lwd = 2, lty = 2)
points(Rm3.5b.tl$tl, Rm3.5b.tl$Sr, 
       pch = 16, cex = 2, col = alpha("cyan4", 0.5))

# 1 LA-ICP MS and modeled serum, assumes no overprint
plot(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.706, 0.712),
     lwd = 2 , 
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     main = "Ref. data 2: LA-ICP-MS En 1")
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[2]], lty = 2)
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[3]], lty = 2)
lines(En1.50avg.tl.f$avg.tl, En1.50avg.tl.f$avg.sr,
      col = "orange", lwd = 2)

# Part II: estimated fraction of post-movement overprint
# record densities
d.pr.En9 <- density(post.comb$BUGSoutput$sims.list$pr.En9)
d.pr.En10 <- density(post.comb$BUGSoutput$sims.list$pr.En10)
d.pr.drill <- density(post.comb$BUGSoutput$sims.list$pr.drill)
d.pr.Rm3.5b <- density(post.comb$BUGSoutput$sims.list$pr.Rm3.5b)

plot(d.pr.En9$y, d.pr.En9$x, type = "l",
     ylim = c(0,0.35), col = "orange3", lwd = 2,
     xlab = "Density", ylab = "Fraction",
     main = "Est. f. post-relocation: En 9")
abline(h = pr.En9.map, lwd = 2)
abline(h = pr.En9.ci1, lty = 2)
abline(h = pr.En9.ci2, lty = 2)

plot(d.pr.En10$y, d.pr.En10$x, type = "l",
     ylim = c(0,0.35), col = "orange4", lwd = 2,
     xlab = "Density", ylab = "Fraction",
     main = "Est. f. post-relocation: En10")
abline(h = pr.En10.map, lwd = 2)
abline(h = pr.En10.ci1, lty = 2)
abline(h = pr.En10.ci2, lty = 2)

plot(d.pr.drill$y, d.pr.drill$x, type = "l",
     ylim = c(0,0.35), col = "red4", lwd = 2,
     xlab = "Density", ylab = "Fraction",
     main = "Est. f. post-relocation: Drill")
abline(h = pr.drill.map, lwd = 2)
abline(h = pr.drill.ci1, lty = 2)
abline(h = pr.drill.ci2, lty = 2)

plot(d.pr.Rm3.5b$y, d.pr.Rm3.5b$x, type = "l",
     ylim = c(0,0.35), col = "cyan4", lwd = 2,
     xlab = "Density", ylab = "Fraction",
     main = "Est. f. post-relocation: Micromill")
abline(h = pr.Rm3.5b.map, lwd = 2)
abline(h = pr.Rm3.5b.ci1, lty = 2)
abline(h = pr.Rm3.5b.ci2, lty = 2)

#################  Fig S9 ################# 
# linear regression showing the negative relationship between
# model residual and sampling depth
# supporting the notion that sampling depth affects sample 87Sr/86Sr
# by incorporating heterogeneous 87Sr/86Sr within the thickness of enamel 
par(mfrow=c(1, 1))
plot(data = res.tib, y ~ x, ylim = c(-5e-4,2e-4),
     xlab = "Drill depth (mm)", ylab = "Residual: data - model",
     main = "Sampling depth vs 87Sr/86Sr residuals",
     pch = 16, col = "red4")

abline(lm.res, lwd = 2, col="blue3")
lines(newx, conf_interval[,2], col="blue3", lty=2)
lines(newx, conf_interval[,3], col="blue3", lty=2)
# draw reference lines
abline(h = 0)
# this is the mean depth of the sampling groove
abline(v = mean(drill.tl.f2$depth)) 

# most of the data between timeline days 0 and 600 have lower 87Sr/86Sr,
# because they are drilled deeper

# no apparent correlation between sample depth and residual of 87Sr/86Sr

#################  Fig S10 #################
# density plots showing posterior of the sensitivity test
denplot(as.mcmc(post.sens), parms = c("a","b","c",
                                      "pr.drill","pr.Rm3.5b","pr.En9", "pr.En10",
                                      "Rpri.mod", "Raft.mod"))

#################  Fig S11 ################# 
# examine modeled serum ratios
post.sens.R1m.89 <- MCMC.CI.bound(post.sens$BUGSoutput$sims.list$R1.m, 0.89)

# preliminary plot
# 2 micromill tusk dentine and modeled serum, assumes no overprint
plot(bin.thin.oc*1:t.oc - 400, post.sens.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.705, 0.712),
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     lwd = 2, main = "Modeled serum 87Sr/86Sr")
lines(bin.thin.oc*1:t.oc - 400, post.sens.R1m.89[[2]], lty = 2)
lines(bin.thin.oc*1:t.oc - 400, post.sens.R1m.89[[3]], lty = 2)
abline(h = Sr.pri.map)
abline(h = Sr.aft.map)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch = 18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
#################  end Fig S11 #################

#################  Fig S12 ################# 
# simulated enamel Sr variation with the forwar model

par(mfrow=c(4, 1))
par(mar = c(5, 4.5, 3, 2))
plot(Misha.sim[[3]], 
     col = viridis(100),
     main="Elephant annual migration, res = 10 microns", 
     xlab = "Length", 
     ylab = "Thickness")

plot(Misha.sim2[[3]], 
     col = viridis(100),
     main="Elephant annual migration, innermost 0.2mm of enamel, res = 10 microns", 
     xlab = "Length", 
     ylab = "Thickness")

plot(sheep.sim[[3]], 
     col = viridis(100),
     main="Sheep annual migration, res = 10 microns", 
     xlab = "Length", 
     ylab = "Thickness")

plot(sheep.sim2[[3]], 
     col = viridis(100),
     main="Sheep migration: twice a year, res = 10 microns", 
     xlab = "Length", 
     ylab = "Thickness")


#################  Fig S13 ################# 
# simulated Sr variation of samples, with the forwar model

par(mfrow=c(4, 1))
par(mar = c(4, 4, 3, 2))
plot(Misha.sim.dist, Misha.sim[[1]], pch = 16, col ="red4",
     ylim = c(syn.low, syn.high), main = "Elephant annual migration, sampling depth = 0.8mm",
     xlab ="Crown length (mm)", ylab = "87Sr/86Sr")
# simulated EDJ: sheep.sim[[2]]
lines(Misha.sim[[2]]$x, Misha.sim[[2]]$Sr, lwd=2, col = "#00b4ffff")
lines(Misha.sim[[2]]$x, Misha.sim[[2]]$Intake)
legend(91,0.711,c("Intake","Serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

# calculate the fraction of intake recovered using conventional drilling

(max(Misha.sim[[1]]) - min(Misha.sim[[1]]))/(syn.high - syn.low)
#31%

plot(Misha.sim.dist, Misha.sim2[[1]], pch = 16, col ="red4",
     ylim = c(syn.low, syn.high), main = "Elephant annual migration, innermost 0.2mm of enamel",
     xlab ="Crown length (mm)", ylab = "87Sr/86Sr")
# simulated EDJ: sheep.sim[[2]]
lines(Misha.sim2[[2]]$x, Misha.sim2[[2]]$Sr, lwd=2, col = "#00b4ffff")
lines(Misha.sim2[[2]]$x, Misha.sim2[[2]]$Intake)

(max(Misha.sim2[[1]]) - min(Misha.sim2[[1]]))/(syn.high - syn.low)
# much higher amplitude at 55.5%

plot(sheep.sim.dist, sheep.sim[[1]], pch = 16, col ="red4",
     ylim = c(syn.low, syn.high), main = "Sheep annual migration",
     xlab ="Crown length (mm)", ylab = "87Sr/86Sr")
# simulated EDJ: sheep.sim[[2]]
lines(sheep.sim[[2]]$x, sheep.sim[[2]]$Sr, lwd=2, col = "#00b4ffff")
lines(sheep.sim[[2]]$x, sheep.sim[[2]]$Intake)


(max(sheep.sim[[1]]) - min(sheep.sim[[1]]))/(syn.high - syn.low)
# 70% of the variation from an annual variation 5-month winter place, 5-months summer place 

plot(sheep.sim.dist, sheep.sim2[[1]], pch = 16, col ="red4",
     ylim = c(syn.low, syn.high), main = "Sheep migration: twice a year",
     xlab ="Crown length (mm)", ylab = "87Sr/86Sr")
# simulated EDJ: sheep.sim[[2]]
lines(sheep.sim2[[2]]$x, sheep.sim2[[2]]$Sr, lwd=2, col = "#00b4ffff")
lines(sheep.sim2[[2]]$x, sheep.sim2[[2]]$Intake)

(max(sheep.sim2[[1]]) - min(sheep.sim2[[1]]))/(syn.high - syn.low)
# 46% of the variation from a seasonal migration 2-month one place, 2-months the other place


#################  Fig S14 ################# 
# Comparing three isotope tracers in tusk and molar enamel##########
# adding c13C and d18O to the curve
par(mfrow=c(3,2))
# tusk dentine micromill Sr
plot(tusk.mill.tl$tl, tusk.mill.tl$Sr,
     xlim=c(-400,700), ylim=c(0.706,0.712),
     pch = 18, cex = 2.2, col ="#00b4ffff",
     xlab="Days from Misha's move",
     main = "Tusk micromill 87Sr/86Sr",
     ylab = "87Sr/86Sr")
abline(v = 0, lty = 2)

plot(Edrill.tl$tl, Edrill.tl$Sr, col= "red4",
     pch=16, cex = 2,
     xlim=c(-400,700), ylim=c(0.706,0.712),
     xlab="Days from Misha's move",
     main = "Molar conventional drill 87Sr/86Sr",
     ylab = "87Sr/86Sr") 
abline(v = 0, lty = 2)

plot(M640.micromill.tl.al, rev(M640.micromill$d13C), 
     xlim=c(-400,700),ylim=c(-15,-4),
     col = "orange",
     pch=18, cex = 2,
     xlab="Days from Misha's move",
     main = "Tusk micromill d13C",
     ylab = "d13C")
abline(v = 0, lty = 2)


plot(Edrill.tl$tl, rev(Drill.no$d13C), 
     xlim=c(-400,700),ylim=c(-15,-4),
     col = "orange",
     pch=16, cex = 2,
     xlab="Days from Misha's move",
     main = "Molar conventional drill d13C",
     ylab = "d13C")
abline(v = 0, lty = 2)

plot(M640.micromill.tl.al, rev(M640.micromill$d18O), 
     xlim=c(-400,700), ylim=c(-15,-7),
     col = "lightblue4",
     pch=18, cex = 2,
     xlab="Days from Misha's move",
     main = "Tusk micromill d18O",
     ylab = "d18O")
abline(v = 0, lty = 2)

plot(Edrill.tl$tl, rev(Drill.no$d18O), 
     xlim=c(-400,700), ylim=c(-15,-7),
     col = "lightblue4",
     pch=16, cex = 2,
     xlab="Days from Misha's move",
     main = "Molar conventional drill d18O",
     ylab = "d18O")
abline(v = 0, lty = 2)

#################  Fig S15 #################
# posteriors of parameter a among BITS runs

par(mfrow=c(1, 2)) #900 * 450

plot(density(exp(log.a)), col = "blue", ylim = c(0,110),
     main="Posteriors of parameter a", xlab = "Parameter a")
lines(density(post.misha.En1$BUGSoutput$sims.list$a/En1.bt),col="red")
legend(0.03,110,c("LA-ICP-MS tusk","LA-ICP-MS enamel"),
       lwd = c(1,1), col=c("blue","red"))


plot(density(exp(log.a)), col = "blue", ylim = c(0,110),
     main="Posteriors of parameter a", xlab = "Parameter a")
lines(density(post.misha.M640b$BUGSoutput$sims.list$a/M640.bt),col="red")
legend(0.03,110,c("LA-ICP-MS tusk","Micromill tusk"),
       lwd = c(1,1), col=c("blue","red"))

#################  Fig S16 #################
# compare estimated Sr intake using BITS vs RnxProg approaches

Rxn.misha <- read.csv("data/RxnProg_Misha ivory.csv")

par(mfrow=c(1, 1))
plot(tusk.mill.tl$tl, tusk.mill.tl$Sr, pch=18, col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on tusk dentine micromill",
     ylab = "87Sr/86Sr") 

lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[1]], lwd = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[2]], lty = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[3]], lty = 2)

lines(Rxn.misha$calculate.date, Rxn.misha$calculated.87Sr.86Sr.input, 
      lty = 1, col = "red", lwd = 2)

legend(600,0.707,c("BITS est. intake","RxnProg est. intake"),lwd = c(2,2), col=c("black","red"))


# plot(density(exp(log.a)), col = "blue", ylim = c(0,110),
#      main="Posteriors of parameter a", xlab = "Parameter a")
# lines(density(post.misha.En3$BUGSoutput$sims.list$a/En3.bt),col="red")
# legend(0.03,110,c("LA-ICP-MS tusk","LA-ICP-MS enamel"),
#        lwd = c(1,1), col=c("blue","red"))

#check posterior density of parameter a:
# this tends to produce higher posterior of rate parameter a


# plot(En3.50avg.tl$avg.tl, En3.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
#      lwd = 2,
#      xlim=c(-400,1000),ylim=c(0.705,0.713),
#      xlab="Days from Misha's move",
#      main = "Estimated intake based on molar enamel LA-ICP-MS",
#      ylab = "87Sr/86Sr") 
# abline(h = CA.Sr)
# abline(h = UT.Sr)
# 
# lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
#       post.misha.En3.Rin.m.89[[1]], lwd = 2)
# lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
#       post.misha.En3.Rin.m.89[[2]], lty = 2)
# lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
#       post.misha.En3.Rin.m.89[[3]], lty = 2)
# legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))
# 
# plot(En5.50avg.tl$avg.tl, En5.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
#      lwd = 2,
#      xlim=c(-400,1000),ylim=c(0.705,0.713),
#      xlab="Days from Misha's move",
#      main = "Estimated intake based on molar enamel LA-ICP-MS",
#      ylab = "87Sr/86Sr") 
# abline(h = CA.Sr)
# abline(h = UT.Sr)
# 
# lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
#       post.misha.En5.Rin.m.89[[1]], lwd = 2)
# lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
#       post.misha.En5.Rin.m.89[[2]], lty = 2)
# lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
#       post.misha.En5.Rin.m.89[[3]], lty = 2)
# legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))
# 
# 
# plot(En7.50avg.tl$avg.tl, En7.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
#      lwd = 2,
#      xlim=c(-400,1000),ylim=c(0.705,0.713),
#      xlab="Days from Misha's move",
#      main = "Estimated intake based on molar enamel LA-ICP-MS",
#      ylab = "87Sr/86Sr") 
# abline(h = CA.Sr)
# abline(h = UT.Sr)
# 
# lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
#       post.misha.En7.Rin.m.89[[1]], lwd = 2)
# lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
#       post.misha.En7.Rin.m.89[[2]], lty = 2)
# lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
#       post.misha.En7.Rin.m.89[[3]], lty = 2)
# legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))

