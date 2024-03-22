library(scales)
library(lattice)

#################  Fig S2 #################  
############# compare position of the change points to appositional angle ###########
par(mfrow=c(1, 2))
plot(cp1.E.rm.x, cp1.E.rm.y, main = "Change point 1", xlim = c(50000,10000),
     xlab = "Distance from the top of EDJ (micron)",
     ylab = "Distance from the EDJ (micron)")
#this is based on measurements by Uno et al 2020, with an appositional angle at 3.3 +- 0.5 degrees
abline(a = 3.02e+03, b = -tan(3.3/180*pi), lwd = 6, col ="pink3")
# abline(a = 3.02e+03, b = -tan(3.25/180*pi), lty = 2, col = "red3")
# abline(a = 3.02e+03, b = -tan(3.35/180*pi), lty = 2, col = "red3")
abline(lm.cp1.E.proj.inner,lwd = 2, col = "blue3")
points(cp1.E.rm.x[1:8], cp1.E.rm.y[1:8], pch = 16) 
points(cp1.E.rm.x[9:10], cp1.E.rm.y[9:10], col = "red2", pch = 16) #mark the ones that don't conform to the angle

plot(cp2.E.rm.x, cp2.E.rm.y, main = "Change point 2",xlim = c(55000,15000),
     xlab = "Distance from the top of EDJ (micron)",
     ylab = "Distance from the EDJ (micron)")
#this is based on measurements by Uno et al 2020, with an appositional angle at 3.3 degrees
abline(a = 3.317e+03, b = -tan(3.3/180*pi), lwd = 6, col ="pink3")
# abline(a = 3.317e+03, b = -tan(3.25/180*pi), lty = 2, col = "red3")
# abline(a = 3.317e+03, b = -tan(3.35/180*pi), lty = 2, col = "red3")
abline(lm.cp2.E.proj.inner,lwd = 2, col = "blue3")
points(cp2.E.rm.x[1:8], cp2.E.rm.y[1:8], pch = 16) #mark the ones that don't conform to the angle
points(cp2.E.rm.x[9:10], cp2.E.rm.y[9:10], col = "red2", pch = 16) #mark the ones that don't conform to the angle

############# end of change point comparisons ###########

#################  Fig S3 #################  
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


#################  Fig S4 ################# 
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

#################  Fig S7 ################# 
# examine modeled serum ratios
post.sens.R1m.89 <- MCMC.CI.bound(post.sens$BUGSoutput$sims.list$R1.m, 0.89)

# preliminary plot
# 2 micromill tusk dentine and modeled serum, assumes no overprint
plot(bin.thin*1:t - 400, post.sens.R1m.89[[1]], type = "l",
     xlim = c(-400, 1100), ylim = c(0.705, 0.712),
     xlab = "Timeline (days)", ylab = "87Sr/86Sr",
     lwd = 2, main = "Modeled serum 87Sr/86Sr")
lines(bin.thin*1:t - 400, post.sens.R1m.89[[2]], lty = 2)
lines(bin.thin*1:t - 400, post.sens.R1m.89[[3]], lty = 2)
abline(h = Sr.pri.map)
abline(h = Sr.aft.map)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch = 18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
#################  end Fig S7 #################

#################  Fig S8 ################# 
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


#################  Fig S9 ################# 
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

#################  Fig S10 #################

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

########## Fig S11 Comparing three isotope tracers in tusk and molar enamel##########
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

plot(drill.tl$tl, drill.tl$Sr, col= "red4",
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


plot(drill.tl$tl, rev(Drill.no$d13C), 
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

plot(drill.tl$tl, rev(Drill.no$d18O), 
     xlim=c(-400,700), ylim=c(-15,-7),
     col = "lightblue4",
     pch=16, cex = 2,
     xlab="Days from Misha's move",
     main = "Molar conventional drill d18O",
     ylab = "d18O")
abline(v = 0, lty = 2)


#################  Fig S12 #################
par(mfrow=c(1, 2)) #1100 * 450
# reconstructed Sr input signal from LA-ICP-MS molar enamel
plot(En1.50avg.tl$avg.tl, En1.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on molar enamel LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[1]], lwd = 2)
lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[2]], lty = 2)
lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[3]], lty = 2)
legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))


# reconstructed Sr input signal from micromilled tusk dentine
plot(tusk.mill.tl$tl, tusk.mill.tl$Sr, pch=18, col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on tusk dentine micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[1]], lwd = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[2]], lty = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[3]], lty = 2)

plot(En3.50avg.tl$avg.tl, En3.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on molar enamel LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
      post.misha.En3.Rin.m.89[[1]], lwd = 2)
lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
      post.misha.En3.Rin.m.89[[2]], lty = 2)
lines(1:t.En3 * En3.bt - d.offset.En3 - offset.en3, 
      post.misha.En3.Rin.m.89[[3]], lty = 2)
legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))

plot(En5.50avg.tl$avg.tl, En5.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on molar enamel LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
      post.misha.En5.Rin.m.89[[1]], lwd = 2)
lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
      post.misha.En5.Rin.m.89[[2]], lty = 2)
lines(1:t.En5 * En5.bt - d.offset.En5 - offset.en5, 
      post.misha.En5.Rin.m.89[[3]], lty = 2)
legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))


plot(En7.50avg.tl$avg.tl, En7.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on molar enamel LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
      post.misha.En7.Rin.m.89[[1]], lwd = 2)
lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
      post.misha.En7.Rin.m.89[[2]], lty = 2)
lines(1:t.En7 * En7.bt - d.offset.En7 - offset.en7, 
      post.misha.En7.Rin.m.89[[3]], lty = 2)
legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))



#################  Fig S13 #################
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

plot(density(exp(log.a)), col = "blue", ylim = c(0,110),
     main="Posteriors of parameter a", xlab = "Parameter a")
lines(density(post.misha.En3$BUGSoutput$sims.list$a/En3.bt),col="red")
legend(0.03,110,c("LA-ICP-MS tusk","LA-ICP-MS enamel"),
       lwd = c(1,1), col=c("blue","red"))

#check posterior density of parameter a:
# this tends to produce higher posterior of rate parameter a
