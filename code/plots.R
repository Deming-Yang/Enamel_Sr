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

##################plot moving averages###################
# par(mfrow=c(3,3))
# 
# plot(nrow(proc.EDJ.dentin2.rm):1,proc.EDJ.dentin2.rm$avg, ylim=c(0.702,0.712),type="l")
# plot(nrow(proc.EDJ.dentin1.rm):1,proc.EDJ.dentin1.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel1.rm),proc.Enamel1.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel2.rm),proc.Enamel2.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel3.rm),proc.Enamel3.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel4.rm),proc.Enamel4.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel5.rm),proc.Enamel5.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel6ext.rm),proc.Enamel6ext.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel7.rm),proc.Enamel7.rm$avg, ylim=c(0.702,0.712),type="l")

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

###########EDJ flattened geometry###########
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




#plot change positions within enamel:
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

#################one transition################
ggplot()+
  geom_sf(data = sf.all.enamel.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp1.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")


#################plot 2 Sr ratio transitions using moving averages###############
par(mfrow=c(2,6))
# preliminary plot
plot(dent.rm.f$new.x, dent.rm.f$avg, main = "Dentine",col="darkgray",
     pch=16,ylim=c(0.705,0.712))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.D.rm, add = T)
lines.segmented(fit_segmented.D.rm)
points.segmented(fit_segmented.D.rm)

#preliminary plot
plot(proc.Enamel1.rm.f$new.x, proc.Enamel1.rm.f$avg, main = "Enamel 1",col="darkgray", 
     pch=16,xlim=c(4e4,8e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E1.rm, add = T)
lines.segmented(fit_segmented.E1.rm)
points.segmented(fit_segmented.E1.rm)

#preliminary plot
plot(proc.Enamel2.rm.f$new.x, proc.Enamel2.rm.f$avg, main = "Enamel 2",col="darkgray", 
     pch=16,xlim=c(2e4,6e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E2.rm, add=T)
lines.segmented(fit_segmented.E2.rm)
points.segmented(fit_segmented.E2.rm)

#preliminary plot
plot(proc.Enamel3.rm.f$new.x, proc.Enamel3.rm.f$avg, main = "Enamel 3",col="darkgray"
     , pch=16,xlim=c(2e4,6e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E3.rm, add=T)
lines.segmented(fit_segmented.E3.rm)
points.segmented(fit_segmented.E3.rm)

#preliminary plot
plot(proc.Enamel4.rm.f$new.x, proc.Enamel4.rm.f$avg, main = "Enamel 4",col="darkgray", 
     pch=16,xlim=c(2e4,6e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E4.rm, add=T)
lines.segmented(fit_segmented.E4.rm)
points.segmented(fit_segmented.E4.rm)

#preliminary plot
plot(proc.Enamel5.rm.f$new.x, proc.Enamel5.rm.f$avg, main = "Enamel 5",col="darkgray"
     , pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E5.rm, add =T)
lines.segmented(fit_segmented.E5.rm)
points.segmented(fit_segmented.E5.rm)

#preliminary plot
plot(proc.Enamel6.rm.f$new.x, proc.Enamel6.rm.f$avg, main = "Enamel 6",col="darkgray"
     , pch=16,xlim=c(0.5e4,4.5e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E6.rm, add =T)
lines.segmented(fit_segmented.E6.rm)
points.segmented(fit_segmented.E6.rm)

#preliminary plot
plot(proc.Enamel7.rm.f$new.x, proc.Enamel7.rm.f$avg, main = "Enamel 7",col="darkgray", 
     pch=16,xlim=c(0.5e4,4.5e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E7.rm, add =T)
lines.segmented(fit_segmented.E7.rm)
points.segmented(fit_segmented.E7.rm)

plot(proc.Enamel8.rm.f$new.x, proc.Enamel8.rm.f$avg, main = "Enamel 8",col="darkgray", 
     pch=16,xlim=c(0.5e4,4.5e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E8.rm, add =T)
lines.segmented(fit_segmented.E8.rm)
points.segmented(fit_segmented.E8.rm)

plot(proc.Enamel9.rm.f$new.x, proc.Enamel9.rm.f$avg, main = "Enamel 9",col="darkgray", 
     pch=16,xlim=c(0,4e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E9.rm, add =T)
lines.segmented(fit_segmented.E9.rm)
points.segmented(fit_segmented.E9.rm)

plot(proc.Enamel10.rm.f$new.x, proc.Enamel10.rm.f$avg, main = "Enamel 10",col="darkgray", 
     pch=16,xlim=c(0,4e4),ylim=c(0.705,0.711))
abline(h = 0.706, lty = 2, lwd = 1.5)
abline(h = 0.709, lty = 2, lwd = 1.5)
plot(fit_segmented.E10.rm, add =T)
lines.segmented(fit_segmented.E10.rm)
points.segmented(fit_segmented.E10.rm)

#################plot 3 Sr value overlay by shifting x values transitions using moving averages###############
#filter out 
dent.rm.fa <- filter(dent.rm.f, dent.rm.f$new.x < cp2.D.rm & dent.rm.f$avg < 0.7085)
dent.rm.fb <- filter(dent.rm.f, dent.rm.f$new.x > cp2.D.rm & dent.rm.f$avg > 0.7085)

dent.rm.f2 <- rbind(dent.rm.fa, dent.rm.fb)

par(mfrow=c(2,5))
#enamel 1, no x- shift is needed
plot(dent.rm.f2$new.x, dent.rm.f2$avg, main = "Enamel 1",col= alpha("lightcyan4", 0.2),
     pch=16, cex=2, xlim=c(4e4,8e4),ylim=c(0.705,0.711))
points(proc.Enamel1.rm.f$new.x, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E1.rm, add = T)
lines.segmented(fit_segmented.E1.rm)
points.segmented(fit_segmented.E1.rm)

#enamel 2, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[2]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 2",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(3.5e4,7.5e4),ylim=c(0.705,0.711))
points(proc.Enamel2.rm.f$new.x, proc.Enamel2.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E2.rm, add = T)
lines.segmented(fit_segmented.E2.rm)
points.segmented(fit_segmented.E2.rm)

#enamel 3, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[3] 

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 3",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(3e4,7e4),ylim=c(0.705,0.711))
points(proc.Enamel3.rm.f$new.x, proc.Enamel3.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E3.rm, add = T)
lines.segmented(fit_segmented.E3.rm)
points.segmented(fit_segmented.E3.rm)

#enamel 4, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[4]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 4",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(2.5e4,6.5e4),ylim=c(0.705,0.711))
points(proc.Enamel4.rm.f$new.x, proc.Enamel4.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E4.rm, add = T)
lines.segmented(fit_segmented.E4.rm)
points.segmented(fit_segmented.E4.rm)

#enamel 5, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[5]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 5",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(2e4,6e4),ylim=c(0.705,0.711))
points(proc.Enamel5.rm.f$new.x, proc.Enamel5.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E5.rm, add = T)
lines.segmented(fit_segmented.E5.rm)
points.segmented(fit_segmented.E5.rm)

#enamel 6, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[6]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 6",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(1.5e4,5.5e4),ylim=c(0.705,0.711))
points(proc.Enamel6.rm.f$new.x, proc.Enamel6.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E6.rm, add = T)
lines.segmented(fit_segmented.E6.rm)
points.segmented(fit_segmented.E6.rm)

#enamel 7, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[7]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 7",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(1e4,5e4),ylim=c(0.705,0.711))
points(proc.Enamel7.rm.f$new.x, proc.Enamel7.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E7.rm, add = T)
lines.segmented(fit_segmented.E7.rm)
points.segmented(fit_segmented.E7.rm)

#enamel 8, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[8]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 8",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(0.5e4,4.5e4),ylim=c(0.705,0.711))
points(proc.Enamel8.rm.f$new.x, proc.Enamel8.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E8.rm, add = T)
lines.segmented(fit_segmented.E8.rm)
points.segmented(fit_segmented.E8.rm)

#enamel 9, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[9]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 9",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(0.5e4,4.5e4),ylim=c(0.705,0.711))
points(proc.Enamel9.rm.f$new.x, proc.Enamel9.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E9.rm, add = T)
lines.segmented(fit_segmented.E9.rm)
points.segmented(fit_segmented.E9.rm)

#enamel 10, x- shift is needed, need to shift dentine values
shift.x <- cp1.D.rm - cp1.E.rm.x[10]

plot(dent.rm.f2$new.x - shift.x, dent.rm.f2$avg, main = "Enamel 10",col= alpha("lightcyan4", 0.2), 
     pch=16, cex=2, xlim=c(0,4e4),ylim=c(0.705,0.711))
points(proc.Enamel10.rm.f$new.x, proc.Enamel10.rm.f$avg, col= alpha("orange", 0.1),pch=16)

abline(h = segmented.D.Sr[1], lty = 2, lwd = 1.5)
abline(h = segmented.D.Sr[2], lty = 2, lwd = 1.5)
plot(fit_segmented.E10.rm, add = T)
lines.segmented(fit_segmented.E10.rm)
points.segmented(fit_segmented.E10.rm)


#################plot 4 contour map of the enamel block#######################
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

###########plot loess model#########################
#create a raster for the enamel block
# bbox <- st_bbox(sf.all.enamel.sim)
# 
# cell_size <- 100
# 
# xgrid <- seq(0, round(bbox$xmax), by=cell_size)
# ygrid <- seq(0, round(bbox$ymax), by=cell_size)
# 
# #use a two-dimensional loess model to predict missing values in each transect!!!!!! this doesn't work
# all.enamel.sim.loess <- loess(avg ~ new.x + new.y, data = all.enamel.sim, 
#                                degree = 2, span=0.18,
#                                normalize = T, family = "gaussian")
# #this takes ~20 sec
# 
# data.fit <-  expand.grid(new.y = ygrid, new.x = xgrid)
# 
# pred.en.sim.loess <- predict(all.enamel.sim.loess, newdata =data.fit)
# 
# # colnames(pred.en.sim.loess) <- xgrid
# # 
# # rownames(pred.en.sim.loess) <- xgrid
# 
# #rast object
# r.en.sim.Sr <- rast(pred.en.sim.loess,resolution=100) #x-axis is flipped!
# 
# image(r.en.sim.Sr, axes = TRUE)
# 
# #stars object
# st.en.sim.Sr <- st_as_stars(pred.en.sim.loess)
# 
# str(attr(st_dimensions(st.en.sim.Sr), "raster"))
# 
# attr(st_dimensions(st.en.sim.Sr), "raster")[[1]] = c(0,0)
# 
# attr(st_dimensions(st.en.sim.Sr), "raster")[[2]] = c("new.x","new.y")
# 
# image(st.en.sim.Sr, axes = TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# st.en.sim.Sr
# 
# as.data.frame(pred.en.sim.loess, colnames())
# 
# rast(as.data.frame(pred.en.sim.loess))
# 
# st.en.sim.Sr <- st_as_stars(pred.en.sim.loess)
# 
# st.en.sim.Sr$A1[[1]] <-pred.en.sim.loess
# 
# st.en.sim.Sr <- st_as_stars(st.en.sim.Sr)
# 
# 
# #convert predictions to sf object for plotting
# 
# #need to convert to a raster
# 
# plot(st.en.sim.Sr)
# contour(x = xgrid, y = ygrid, z = pred.en.sim.loess, xlab = "microns", ylab = "microns", add = T)
# 
# ggplot()+geom_stars()
# 
# ggplot()+
#   geom_stars(data = st.en.sim.Sr , aes(color = avg))+
#   scale_color_viridis_c()+
#   theme(legend.position = "bottom")
# 
# contour(x = xgrid, y = ygrid, z = pred.en.sim.loess, xlab = "microns", ylab = "microns")
# 
# r.en.Sr <- rast(pred.en.sim.loess, xmn=0, xmx=max(xgrid), ymn=0, ymx=max(ygrid))
# 
# plot(r.en.Sr)

##################solution 2: use geostats to predict missing values#################
# all.bbox <- st_bbox(sf.all.dat.rm.proj)
# 
# all.grid = st_as_stars(all.bbox, dx = 100, dy = 100)
# 
# variog.Sr <- variogram(avg ~ x + y , sf.all.dat.rm.proj)
# 
# plot(variog.Sr, numbers = TRUE)
# 
# modSill = 1.3e-6
# modRange = 2500
# modNugget = 5e-7
# Sr.vgm1 = vgm(psill=modSill, "Sph", range=modRange, nugget=modNugget)
# 
# Sr.vgm2 = fit.variogram(variog.Sr, Sr.vgm1) 
# plot(variog.Sr, Sr.vgm2)
# 
# en.Sr.pred.ok <- krige(avg ~ x + y, sf.all.dat.rm.proj, all.grid, Sr.vgm2, 
#                       nmin=20, nmax=200, maxdist=5e3)
# 
# plot(en.Sr.pred.ok["var1.pred"])#doesn't seem to work!

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

# st_bbox(en.Sr.pred.uk)
# 
# contour(en.Sr.pred.uk$var1.pred, levels = seq(0.705,0.712,0.001))
# 
# cont.en.Sr.pred.uk <- st_contour(en.Sr.pred.uk, contour_lines =F, breaks = seq(0.706,0.711,0.001),
#                                  as_points=FALSE, merge=TRUE)

# xgrid <- seq(0, round(bbox$xmax), by=cell_size)
# ygrid <- seq(0, round(bbox$ymax), by=cell_size)
# 
# #use a two-dimensional loess model to smooth the values
# all.enamel.sim.loess <- loess(avg ~ new.x + new.y, data = all.enamel.sim, 
#                               degree = 2, span=0.18,
#                               normalize = T, family = "gaussian")
# #this takes ~20 sec
# 
# data.fit <-  expand.grid(new.y = ygrid, new.x = xgrid)
# 
# pred.en.sim.loess <- predict(all.enamel.sim.loess, newdata =data.fit)
# 
# 
# bbox <- st_bbox(sf.all.dat.rm.proj)
# 
# cell_size <- 100
# 
# xgrid <- seq(bbox$xmin, bbox$xmax, by=cell_size)
# ygrid <- seq(0, bbox$ymax, by=cell_size)
# 
# r.en.Sr <- raster(ncols = length(xgrid), nrows = length(ygrid),
#                   xmn=0, xmx=max(xgrid), ymn=0, ymx=max(ygrid))
# 
# 
# gs.enamel.rm.proj <- gstat(formula=avg~1, locations=~new.x + new.y, data=sf.all.enamel.rm)
# 
# 
# 
# 
# 

#replace values that are too high or too low

#create cutoff points using 100 point average
