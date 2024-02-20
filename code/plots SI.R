library(scales)

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
# linear regression showing the posivie relationship between
# model resitual and sampling depth
# supporting the notion that sampling depth affects sample 87Sr/86Sr
# by incorporating heterogeneous 87Sr/86Sr within the thickness of enamel 
par(mfrow=c(1, 1))
plot(data = res.tib, y ~ x, ylim = c(-2e-4,5e-4),
     xlab = "Drill depth (mm)", ylab = "Residual: model - data",
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
     lwd = 2, main = "")
lines(bin.thin*1:t - 400, post.sens.R1m.89[[2]], lty = 2)
lines(bin.thin*1:t - 400, post.sens.R1m.89[[3]], lty = 2)
abline(h = CA.Sr)
abline(h = UT.Sr)
points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
       pch = 18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
#################  end Fig S7 #################

# preliminary plots shows substantial over-fitting of the data (solution method),
# which is expected due to the small sds. 
# this suggests that sds should be adjusted to reflect real-world data and model uncertainties

############# Time line reconstructions with individual data series in Fig 3D #############
# Fig SI after converting to time, 
# add +- sd in growth rates to create 67% confidence intervals around the reconstructed timeline
# okay, there is no easy way of doing this...
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












################# contour map of the enamel block#######################
#method 1 for each transet, create a precited data series

#enamel 1
xgrid <- seq(en.bbox[1], en.bbox[3], by = 300) #consider necessary data density

#simulate uncertainty 
sd.Sr.en <- 0.00005

#y axis values using a mean
sim.Enamel1.y <- mean(proc.Enamel1.rm.f$new.y)

sim.Enamel1.x1 <- xgrid[which(xgrid < min(proc.Enamel1.rm.f$new.x))]
sim.Enamel1.x2 <-xgrid[which(xgrid > max(proc.Enamel1.rm.f$new.x))]

pred.en1.segm1 <- predict.segmented(fit_segmented.E1.rm, newdata = data.frame(new.x = sim.Enamel1.x1))
#add uncertainty:
sim.en1.segm1.Sr <- pred.en1.segm1 + rnorm(length(pred.en1.segm1), 0, sd.Sr.en)

#combine segments:
sim.en1 <- data.frame(avg = c(sim.en1.segm1.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel1.x1),
                      new.y = c(rep(sim.Enamel1.y, length(sim.Enamel1.x1))))


sim.en1.ext <- rbind(sim.en1, proc.Enamel1.rm.f)

#enamel 2
#y axis values using a mean
sim.Enamel2.y <- mean(proc.Enamel2.rm.f$new.y)

sim.Enamel2.x1 <- xgrid[which(xgrid < min(proc.Enamel2.rm.f$new.x))]
sim.Enamel2.x2 <-xgrid[which(xgrid > max(proc.Enamel2.rm.f$new.x))]

pred.en2.segm1 <- predict.segmented(fit_segmented.E2.rm, newdata = data.frame(new.x = sim.Enamel2.x1))
#check unrealistic values
pred.en2.segm1[which(pred.en2.segm1 < 0.706)] <- 0.706
#add uncertainty:
sim.en2.segm1.Sr <- pred.en2.segm1 + rnorm(length(pred.en2.segm1), 0, sd.Sr.en)

pred.en2.segm2 <- predict.segmented(fit_segmented.E2.rm, newdata = data.frame(new.x = sim.Enamel2.x2))
#add uncertainty:
sim.en2.segm2.Sr <- pred.en2.segm2 + rnorm(length(pred.en2.segm2), 0, sd.Sr.en)

#combine segments:
sim.en2 <- data.frame(avg = c(sim.en2.segm1.Sr, sim.en2.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel2.x1, sim.Enamel2.x2),
                      new.y = c(rep(sim.Enamel2.y, length(sim.Enamel2.x1)),rep(sim.Enamel2.y, length(sim.Enamel2.x2))))


sim.en2.ext <- rbind(sim.en2, proc.Enamel2.rm.f)

#enamel 3
#y axis values using a mean
sim.Enamel3.y <- mean(proc.Enamel3.rm.f$new.y)

sim.Enamel3.x1 <- xgrid[which(xgrid < min(proc.Enamel3.rm.f$new.x))]
sim.Enamel3.x2 <-xgrid[which(xgrid > max(proc.Enamel3.rm.f$new.x))]

pred.en3.segm1 <- predict.segmented(fit_segmented.E3.rm, newdata = data.frame(new.x = sim.Enamel3.x1))
#check unrealistic values
pred.en3.segm1[which(pred.en3.segm1 < 0.706)] <- 0.706
#add uncertainty:
sim.en3.segm1.Sr <- pred.en3.segm1 + rnorm(length(pred.en3.segm1), 0, sd.Sr.en)

pred.en3.segm2 <- predict.segmented(fit_segmented.E3.rm, newdata = data.frame(new.x = sim.Enamel3.x2))
#add uncertainty:
sim.en3.segm2.Sr <- pred.en3.segm2 + rnorm(length(pred.en3.segm2), 0, sd.Sr.en)

#combine segments:
sim.en3 <- data.frame(avg = c(sim.en3.segm1.Sr, sim.en3.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel3.x1, sim.Enamel3.x2),
                      new.y = c(rep(sim.Enamel3.y, length(sim.Enamel3.x1)),rep(sim.Enamel3.y, length(sim.Enamel3.x2))))


sim.en3.ext <- rbind(sim.en3, proc.Enamel3.rm.f)

#enamel 4
#y axis values using a mean
sim.Enamel4.y <- mean(proc.Enamel4.rm.f$new.y)

sim.Enamel4.x1 <- xgrid[which(xgrid < min(proc.Enamel4.rm.f$new.x))]
sim.Enamel4.x2 <-xgrid[which(xgrid > max(proc.Enamel4.rm.f$new.x))]

pred.en4.segm1 <- predict.segmented(fit_segmented.E4.rm, newdata = data.frame(new.x = sim.Enamel4.x1))
#add uncertainty:
sim.en4.segm1.Sr <- pred.en4.segm1 + rnorm(length(pred.en4.segm1), 0, sd.Sr.en)

pred.en4.segm2 <- predict.segmented(fit_segmented.E4.rm, newdata = data.frame(new.x = sim.Enamel4.x2))
#check unrealistic values
pred.en4.segm2[which(pred.en4.segm2 > 0.7118)] <- 0.7118
#add uncertainty:
sim.en4.segm2.Sr <- pred.en4.segm2 + rnorm(length(pred.en4.segm2), 0, sd.Sr.en)

#combine segments:
sim.en4 <- data.frame(avg = c(sim.en4.segm1.Sr, sim.en4.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel4.x1, sim.Enamel4.x2),
                      new.y = c(rep(sim.Enamel4.y, length(sim.Enamel4.x1)),rep(sim.Enamel4.y, length(sim.Enamel4.x2))))


sim.en4.ext <- rbind(sim.en4, proc.Enamel4.rm.f)

#enamel 5
#y axis values using a mean
sim.Enamel5.y <- mean(proc.Enamel5.rm.f$new.y)

sim.Enamel5.x1 <- xgrid[which(xgrid < min(proc.Enamel5.rm.f$new.x))]
sim.Enamel5.x2 <-xgrid[which(xgrid > max(proc.Enamel5.rm.f$new.x))]

pred.en5.segm1 <- predict.segmented(fit_segmented.E5.rm, newdata = data.frame(new.x = sim.Enamel5.x1))
#add uncertainty:
sim.en5.segm1.Sr <- pred.en5.segm1 + rnorm(length(pred.en5.segm1), 0, sd.Sr.en)

pred.en5.segm2 <- predict.segmented(fit_segmented.E5.rm, newdata = data.frame(new.x = sim.Enamel5.x2))
#check unrealistic values
pred.en5.segm2[which(pred.en5.segm2 > 0.7118)] <- 0.7118
#add uncertainty:
sim.en5.segm2.Sr <- pred.en5.segm2 + rnorm(length(pred.en5.segm2), 0, sd.Sr.en)

#combine segments:
sim.en5 <- data.frame(avg = c(sim.en5.segm1.Sr, sim.en5.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel5.x1, sim.Enamel5.x2),
                      new.y = c(rep(sim.Enamel5.y, length(sim.Enamel5.x1)),rep(sim.Enamel5.y, length(sim.Enamel5.x2))))


sim.en5.ext <- rbind(sim.en5, proc.Enamel5.rm.f)

#enamel 6
#y axis values using a mean
sim.Enamel6.y <- mean(proc.Enamel6.rm.f$new.y)

sim.Enamel6.x1 <- xgrid[which(xgrid < min(proc.Enamel6.rm.f$new.x))]
sim.Enamel6.x2 <-xgrid[which(xgrid > max(proc.Enamel6.rm.f$new.x))]

pred.en6.segm1 <- predict.segmented(fit_segmented.E6.rm, newdata = data.frame(new.x = sim.Enamel6.x1))
#add uncertainty:
sim.en6.segm1.Sr <- pred.en6.segm1 + rnorm(length(pred.en6.segm1), 0, sd.Sr.en)

pred.en6.segm2 <- predict.segmented(fit_segmented.E6.rm, newdata = data.frame(new.x = sim.Enamel6.x2))
#add uncertainty:
sim.en6.segm2.Sr <- pred.en6.segm2 + rnorm(length(pred.en6.segm2), 0, sd.Sr.en)

#combine segments:
sim.en6 <- data.frame(avg = c(sim.en6.segm1.Sr, sim.en6.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel6.x1, sim.Enamel6.x2),
                      new.y = c(rep(sim.Enamel6.y, length(sim.Enamel6.x1)),
                                rep(sim.Enamel6.y, length(sim.Enamel6.x2))))


sim.en6.ext <- rbind(sim.en6, proc.Enamel6.rm.f)

#enamel 7
#y axis values using a mean
sim.Enamel7.y <- mean(proc.Enamel7.rm.f$new.y)

sim.Enamel7.x1 <- xgrid[which(xgrid < min(proc.Enamel7.rm.f$new.x))]
sim.Enamel7.x2 <-xgrid[which(xgrid > max(proc.Enamel7.rm.f$new.x))]

pred.en7.segm1 <- predict.segmented(fit_segmented.E7.rm, newdata = data.frame(new.x = sim.Enamel7.x1))
#add uncertainty:
sim.en7.segm1.Sr <- pred.en7.segm1 + rnorm(length(pred.en7.segm1), 0, sd.Sr.en)

pred.en7.segm2 <- predict.segmented(fit_segmented.E7.rm, newdata = data.frame(new.x = sim.Enamel7.x2))
#add uncertainty:
sim.en7.segm2.Sr <- pred.en7.segm2 + rnorm(length(pred.en7.segm2), 0, sd.Sr.en)

#combine segments:
sim.en7 <- data.frame(avg = c(sim.en7.segm1.Sr, sim.en7.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel7.x1, sim.Enamel7.x2),
                      new.y = c(rep(sim.Enamel7.y, length(sim.Enamel7.x1)),
                                rep(sim.Enamel7.y, length(sim.Enamel7.x2))))


sim.en7.ext <- rbind(sim.en7, proc.Enamel7.rm.f)

#enamel 8
#y axis values using a mean
sim.Enamel8.y <- mean(proc.Enamel8.rm.f$new.y)

sim.Enamel8.x2 <-xgrid[which(xgrid > max(proc.Enamel8.rm.f$new.x))]

pred.en8.segm2 <- predict.segmented(fit_segmented.E8.rm, newdata = data.frame(new.x = sim.Enamel8.x2))
#check unrealistic values
pred.en8.segm2[which(pred.en8.segm2 > 0.7118)] <- 0.7118
#add uncertainty:
sim.en8.segm2.Sr <- pred.en8.segm2 + rnorm(length(pred.en8.segm2), 0, sd.Sr.en)

#combine segments:
sim.en8 <- data.frame(avg = c(sim.en8.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c( sim.Enamel8.x2),
                      new.y = c(rep(sim.Enamel8.y, length(sim.Enamel8.x2))))


sim.en8.ext <- rbind(sim.en8, proc.Enamel8.rm.f)

#enamel 9
#y axis values using a mean
sim.Enamel9.y <- mean(proc.Enamel9.rm.f$new.y)

sim.Enamel9.x1 <- xgrid[which(xgrid < min(proc.Enamel9.rm.f$new.x))]
sim.Enamel9.x2 <-xgrid[which(xgrid > max(proc.Enamel9.rm.f$new.x))]

pred.en9.segm1 <- predict.segmented(fit_segmented.E9.rm, newdata = data.frame(new.x = sim.Enamel9.x1))
#add uncertainty:
sim.en9.segm1.Sr <- pred.en9.segm1 + rnorm(length(pred.en9.segm1), 0, sd.Sr.en)

pred.en9.segm2 <- predict.segmented(fit_segmented.E9.rm, newdata = data.frame(new.x = sim.Enamel9.x2))
#add uncertainty:
sim.en9.segm2.Sr <- pred.en9.segm2 + rnorm(length(pred.en9.segm2), 0, sd.Sr.en)

#combine segments:
sim.en9 <- data.frame(avg = c(sim.en9.segm1.Sr, sim.en9.segm2.Sr),
                      x = NA, y = NA,
                      new.x = c(sim.Enamel9.x1, sim.Enamel9.x2),
                      new.y = c(rep(sim.Enamel9.y, length(sim.Enamel9.x1)),
                                rep(sim.Enamel9.y, length(sim.Enamel9.x2))))


sim.en9.ext <- rbind(sim.en9, proc.Enamel9.rm.f)

#enamel 10
#y axis values using a mean
sim.Enamel10.y <- mean(proc.Enamel10.rm.f$new.y)

sim.Enamel10.x1 <- xgrid[which(xgrid < min(proc.Enamel10.rm.f$new.x))]
sim.Enamel10.x2 <-xgrid[which(xgrid > max(proc.Enamel10.rm.f$new.x))]

pred.en10.segm1 <- predict.segmented(fit_segmented.E10.rm, newdata = data.frame(new.x = sim.Enamel10.x1))
#add uncertainty:
sim.en10.segm1.Sr <- pred.en10.segm1 + rnorm(length(pred.en10.segm1), 0, sd.Sr.en)

pred.en10.segm2 <- predict.segmented(fit_segmented.E10.rm, newdata = data.frame(new.x = sim.Enamel10.x2))
#add uncertainty:
sim.en10.segm2.Sr <- pred.en10.segm2 + rnorm(length(pred.en10.segm2), 0, sd.Sr.en)

#combine segments:
sim.en10 <- data.frame(avg = c(sim.en10.segm1.Sr, sim.en10.segm2.Sr),
                       x = NA, y = NA,
                       new.x = c(sim.Enamel10.x1, sim.Enamel10.x2),
                       new.y = c(rep(sim.Enamel10.y, length(sim.Enamel10.x1)),
                                 rep(sim.Enamel10.y, length(sim.Enamel10.x2))))


sim.en10.ext <- rbind(sim.en10, proc.Enamel10.rm.f)

#convert to sf object
all.enamel.sim <- rbind(sim.en1.ext, sim.en2.ext, sim.en3.ext, sim.en4.ext, 
                        sim.en5.ext, sim.en6.ext, sim.en7.ext, sim.en8.ext,
                        sim.en9.ext, sim.en10.ext)

all.enamel.sim.xy <- data.frame(avg = all.enamel.sim$avg, x = all.enamel.sim$new.x, y = all.enamel.sim$new.y,
                                new.x = all.enamel.sim$new.x, new.y = all.enamel.sim$new.y)

sf.all.enamel.sim.xy <- st_as_sf(all.enamel.sim.xy,  agr = NA_agr_,
                                 coords = c("new.x","new.y"),
                                 dim = "XY",
                                 remove = TRUE,
                                 na.fail = TRUE,
                                 sf_column_name = NULL)

#check data for unrealistic result
ggplot()+
  geom_sf(data = sf.all.enamel.sim , aes(color = avg))+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")


##################solution 2: use geostats to predict missing values#################

#solution 2: use geostats to predict missing values and Enamel only
en.bbox <- st_bbox(sf.all.enamel.sim.xy)

en.bbox[1] <- 800

en.bbox[2] <- 200

en.bbox[3] <- 77000

en.bbox[4] <- 3200

en.grid = st_as_stars(en.bbox, dx = 100, dy = 100)

variog.Sr <- variogram(avg ~ x * y , data = sf.all.enamel.sim.xy)

plot(variog.Sr, numbers = TRUE)

modSill = 6e-7
modRange = 5000
modNugget = 4e-7
Sr.vgm1 = vgm(psill=modSill, "Sph", range=modRange, nugget=modNugget)

Sr.vgm2 = fit.variogram(variog.Sr, Sr.vgm1) 
plot(variog.Sr, Sr.vgm2)

en.Sr.pred.uk <- krige(formula = avg ~ x * y, sf.all.enamel.sim.xy, en.grid, Sr.vgm2, 
                       nmin=40, nmax=100, maxdist=5e3)

# cont.en.Sr.pred.uk <- contour(en.Sr.pred.uk, levels = seq(0.705,0.712,0.001), plot = F)

cont.en.Sr.pred.uk <- st_contour(en.Sr.pred.uk, contour_lines =T, breaks = seq(0.705,0.712,0.001))
par(mfrow=c(2,1))


# plot(en.Sr.pred.uk, axes = T, breaks = seq(0.705,0.712,0.001), nbreaks = 8, col = viridis) #this works now!

plot(en.Sr.pred.uk, axes = T, breaks = seq(0.705,0.712,0.0005), nbreaks = 15, col = viridis) #this works now!

plot(cont.en.Sr.pred.uk, col = "black",add = T)#????