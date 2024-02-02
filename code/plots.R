library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(changepoint)
library(ggplot2)
library(gstat)
library(stars)
library(segmented)

###########real world tooth geometry###########
#############plot 1 all dentine and enamel data with transition points#######
par(mfrow=c(1,1))
ggplot()+
  geom_sf(data = sf.all.dat.rm , aes(color = avg))+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")

###########EDJ flattened geometry###########
##############both transitions#################
par(mfrow=c(1,1))
ggplot()+
  geom_sf(data = sf.all.dat.rm.proj , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp2.loc, color = "black",pch=18,cex=3)+ 
  geom_sf(data = sf.cp1.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")

################# LA-ICP-MS enamel-dentine comparison #################
################# overlay with shifted x values in change points ###############
par(mfrow=c(2,5))
#enamel 1, no x- shift is needed
plot(dent.rm.f2$new.x, dent.rm.f2$avg, main = "Enamel 1",col= "gray24",
     type="l", lwd=2, xlim=c(4e4,8e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x, rev(dent.rm.f2$new.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x, dent.rm.f2$avg,col= "gray24", lwd=2)


polygon(c(proc.Enamel1.rm.f$new.x, rev(proc.Enamel1.rm.f$new.x)), 
        c(proc.Enamel1.rm.f$avg + proc.Enamel1.rm.f$sd, 
          rev(proc.Enamel1.rm.f$avg - proc.Enamel1.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel1.rm.f$new.x, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E1.rm, add = T)
lines.segmented(fit_segmented.E1.rm)
points.segmented(fit_segmented.E1.rm)

#enamel 2, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[2]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 2",col= "gray24",
     type="l", lwd=2, xlim=c(3.5e4,7.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel2.rm.f$new.x, rev(proc.Enamel2.rm.f$new.x)), 
        c(proc.Enamel2.rm.f$avg + proc.Enamel2.rm.f$sd, 
          rev(proc.Enamel2.rm.f$avg - proc.Enamel2.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel2.rm.f$new.x, proc.Enamel2.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E2.rm, add = T)
lines.segmented(fit_segmented.E2.rm)
points.segmented(fit_segmented.E2.rm)

#enamel 3, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[3] 

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 3",col= "gray24",
     type="l", lwd=2, xlim=c(3e4,7e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel3.rm.f$new.x, rev(proc.Enamel3.rm.f$new.x)), 
        c(proc.Enamel3.rm.f$avg + proc.Enamel3.rm.f$sd, 
          rev(proc.Enamel3.rm.f$avg - proc.Enamel3.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel3.rm.f$new.x, proc.Enamel3.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E3.rm, add = T)
lines.segmented(fit_segmented.E3.rm)
points.segmented(fit_segmented.E3.rm)

#enamel 4, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[4]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 4",col= "gray24",
     type="l", lwd=2, xlim=c(2.5e4,6.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel4.rm.f$new.x, rev(proc.Enamel4.rm.f$new.x)), 
        c(proc.Enamel4.rm.f$avg + proc.Enamel4.rm.f$sd, 
          rev(proc.Enamel4.rm.f$avg - proc.Enamel4.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel4.rm.f$new.x, proc.Enamel4.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E4.rm, add = T)
lines.segmented(fit_segmented.E4.rm)
points.segmented(fit_segmented.E4.rm)

#enamel 5, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[5]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 5",col= "gray24",
     type="l", lwd=2, xlim=c(1.5e4,5.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel5.rm.f$new.x, rev(proc.Enamel5.rm.f$new.x)), 
        c(proc.Enamel5.rm.f$avg + proc.Enamel5.rm.f$sd, 
          rev(proc.Enamel5.rm.f$avg - proc.Enamel5.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel5.rm.f$new.x, proc.Enamel5.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E5.rm, add = T)
lines.segmented(fit_segmented.E5.rm)
points.segmented(fit_segmented.E5.rm)

#enamel 6, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[6]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 6",col= "gray24",
     type="l", lwd=2, xlim=c(1.2e4,5.2e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel6.rm.f$new.x, rev(proc.Enamel6.rm.f$new.x)), 
        c(proc.Enamel6.rm.f$avg + proc.Enamel6.rm.f$sd, 
          rev(proc.Enamel6.rm.f$avg - proc.Enamel6.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel6.rm.f$new.x, proc.Enamel6.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E6.rm, add = T)
lines.segmented(fit_segmented.E6.rm)
points.segmented(fit_segmented.E6.rm)

#enamel 7, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[7]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 7",col= "gray24",
     type="l", lwd=2, xlim=c(1e4,5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel7.rm.f$new.x, rev(proc.Enamel7.rm.f$new.x)), 
        c(proc.Enamel7.rm.f$avg + proc.Enamel7.rm.f$sd, 
          rev(proc.Enamel7.rm.f$avg - proc.Enamel7.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel7.rm.f$new.x, proc.Enamel7.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E7.rm, add = T)
lines.segmented(fit_segmented.E7.rm)
points.segmented(fit_segmented.E7.rm)

#enamel 8, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[8]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 8",col= "gray24",
     type="l", lwd=2, xlim=c(0.5e4,4.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel8.rm.f$new.x, rev(proc.Enamel8.rm.f$new.x)), 
        c(proc.Enamel8.rm.f$avg + proc.Enamel8.rm.f$sd, 
          rev(proc.Enamel8.rm.f$avg - proc.Enamel8.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel8.rm.f$new.x, proc.Enamel8.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E8.rm, add = T)
lines.segmented(fit_segmented.E8.rm)
points.segmented(fit_segmented.E8.rm)

#enamel 9, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[9]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 9",col= "gray24",
     type="l", lwd=2, xlim=c(0.5e4,4.5e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel9.rm.f$new.x, rev(proc.Enamel9.rm.f$new.x)), 
        c(proc.Enamel9.rm.f$avg + proc.Enamel9.rm.f$sd, 
          rev(proc.Enamel9.rm.f$avg - proc.Enamel9.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel9.rm.f$new.x, proc.Enamel9.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E9.rm, add = T)
lines.segmented(fit_segmented.E9.rm)
points.segmented(fit_segmented.E9.rm)

#enamel 10, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[10]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 10",col= "gray24",
     type="l", lwd=2, xlim=c(0e4,4e4),ylim=c(0.704,0.712)) #thick line
polygon(c(dent.rm.f2$new.x - shift.x, rev(dent.rm.f2$new.x - shift.x)), 
        c(dent.rm.f2$avg + dent.rm.f2$sd, rev(dent.rm.f2$avg - dent.rm.f2$sd)), 
        col = "gray60", border = NA)
lines(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg,col= "gray24", lwd=2)

polygon(c(proc.Enamel10.rm.f$new.x, rev(proc.Enamel10.rm.f$new.x)), 
        c(proc.Enamel10.rm.f$avg + proc.Enamel10.rm.f$sd, 
          rev(proc.Enamel10.rm.f$avg - proc.Enamel10.rm.f$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(proc.Enamel10.rm.f$new.x, proc.Enamel10.rm.f$avg, col= alpha("orange", 0.8),lwd=2)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E10.rm, add = T)
lines.segmented(fit_segmented.E10.rm)
points.segmented(fit_segmented.E10.rm)

################# Multi-substrate comparison with x axis transformed to time axis###############

# Fig 3 A compare laser, hand drill and micromill transects using the molar plate geometry
par(mfrow=c(1,1))
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

# Enamel transect 5
polygon(c(proc.Enamel5.rm.new.x, rev(proc.Enamel5.rm.new.x)), 
        c(proc.Enamel5.rm.f$avg + proc.Enamel5.rm.f$sd, 
          rev(proc.Enamel5.rm.f$avg - proc.Enamel5.rm.f$sd)), 
        col = alpha("orange3", 0.3), border = NA)
lines(proc.Enamel5.rm.new.x, proc.Enamel5.rm.f$avg, col= alpha("orange3", 0.9),lwd=2)

# Enamel transect 10
polygon(c(proc.Enamel10.rm.new.x, rev(proc.Enamel10.rm.new.x)), 
        c(proc.Enamel10.rm.f$avg + proc.Enamel10.rm.f$sd, 
          rev(proc.Enamel10.rm.f$avg - proc.Enamel10.rm.f$sd)), 
        col = alpha("orange4", 0.3), border = NA)
lines(proc.Enamel10.rm.new.x, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.9),lwd=2)

points(Drill.no$Dist..From.cervix, Drill.no$corr..87Sr.86Sr, pch=16, cex = 1.2, col ="red4")

# Fig 3 B
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

# Fig 3 C after converting to time, 
# compare laser, hand drill and micromill transects using the molar plate geometry

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