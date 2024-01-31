library(scales)

############# Validation, compare micromilled tusk dentine data with corrected LA-ICP-MS data #############
#the micromill was done at 500-micron increment, but a geometric conversion is needed (explained in SI)
#here the conversion is done by multiplying 547 (microns)
misha.micromill.dist <- misha.micromill$Position * 1.6 - 5 #best effort

misha.micromill.sr <- misha.micromill$corr..87Sr.86Sr
misha.micromill.sr.err <- misha.micromill$comb..Err

plot(misha.micromill.dist, misha.micromill.sr,type = "l",col="#00b4ffff",lwd=1.5,
     xlim=c(100,0),ylim=c(0.705,0.712),xlab="Distance from pulp cavity (microns)",ylab="87Sr/86Sr",
     main="Comparing data from the LA-ICP-MS and the solution methods")

#preliminary plot
plot(Drill.no$Dist..From.cervix, Drill.no$corr..87Sr.86Sr,ylim = c(0.706,0.712),pch=16,col ="red")

lines(84-(1:n), avg.Sr.samp, ylim = c(0.706,0.712),xlim = c(0,100),col = "blue",lwd =2)

#preliminary plot comparing to dentine transect:
#convert dist from crown top to dist from cervix
dent.rm.new.x<- 92-dent.rm.f2$new.x/1e3
proc.Enamel10.rm.new.x<- 92-proc.Enamel10.rm.f$new.x/1e3
proc.Enamel5.rm.new.x<- 92-proc.Enamel5.rm.f$new.x/1e3
proc.Enamel1.rm.new.x<- 92-proc.Enamel1.rm.f$new.x/1e3

# transform this into distance along EDJ, not vertical distance, use trigonometry
# top of the Rm3.5b is 81mm from the cervix
# angle is approx. 3.2 degrees (Uno 2012)

Rm3.5b.mill.dist <- 81 - Rm3.5b.mill.no$Dist..From.EDJ/tan(3.2/180*pi)/1e3 # in mm

plot(Drill.no$Dist..From.cervix, Drill.no$corr..87Sr.86Sr,ylim = c(0.706,0.712),pch=16,col ="red")
points(Rm3.5b.mill.dist, Rm3.5b.mill.no$corr..87Sr.86Sr,col = "blue" )

plot(dent.rm.new.x, dent.rm.f2$avg, col= alpha("lightcyan4", 0.2),
     pch=16, cex=1, xlim=c(100,0),ylim=c(0.705,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")
points(proc.Enamel1.rm.new.x, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.1),pch=16)
points(proc.Enamel5.rm.new.x, proc.Enamel5.rm.f$avg, col= alpha("orange3", 0.1),pch=16)
points(proc.Enamel10.rm.new.x, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.1),pch=16)
points(Drill.no$Dist..From.cervix, Drill.no$corr..87Sr.86Sr, pch=18, cex = 2, col ="red")




########## use growth rates to match tusk series and enamel LA-ICP-MS series
# to estimate time line, linear interpolations have to be made based on the growth curve
# assume that the top of the crown is at 100 mm to the cervix
# assume that the enamel growth rate at the top of the crown is 55.3 micron/day
# note that all the calculations depend on the 55.3 micron/day extension rate estimate

Ref.proc.en1x <- tan((-0.346*log(100) + 4.478)/180*pi)*55.3e-3 # mm/day

# create reference grid: 100mm to 0.1 mm 
ref.dist <- seq(100,0.1, by = -0.05) # 0.05mm interval

ref.rate <- Ref.proc.en1x/tan((-0.346*log(ref.dist) + 4.478)/180*pi) # mm per day

rate.en1x <- approx(x = ref.dist, y = ref.rate, xout = proc.Enamel1.rm.new.x) # mm per day 

interv.en1x <- c(base::diff(-1*proc.Enamel1.rm.new.x)[1], base::diff(-1*proc.Enamel1.rm.new.x))

# calculate time interval
time.interv <- interv.en1x / rate.en1x$y

days.cumm.en1 <- cumsum(time.interv) #cumulative days

# align with the switch
days.cumm.en1.al <- days.cumm.en1 - 400

# use micromill tusk dentine data in the comparison
misha.micromill.dist <- (42 - misha.micromill$Position) *0.5 # in mm, micromill was done at 0.5 mm interval
misha.micromill.tl <- misha.micromill.dist*1000/14.7 # assuming a constant growth rate

# align with the switch
misha.micromill.tl.al <- misha.micromill.tl - 330

# plot time line of Sr isotope change
plot(days.cumm.en1.al, proc.Enamel1.rm.f$avg, col= "gray24",
     type="l", lwd=2, xlim=c(-400,800),ylim=c(0.705,0.712),
     xlab="Days from moves",
     main = "Conventional vs LA-ICP-MS",
     ylab = "87Sr/86Sr") #thick line
polygon(c(days.cumm.en1.al, rev(days.cumm.en1.al)), 
        c(proc.Enamel1.rm.f$avg + proc.Enamel1.rm.f$sd, 
          rev(proc.Enamel1.rm.f$avg - proc.Enamel1.rm.f$sd)), 
        col = "gray60", border = NA)
lines(days.cumm.en1.al, proc.Enamel1.rm.f$avg,col= "gray24", lwd=2)

points(misha.micromill.tl.al, misha.micromill.sr,pch=18,cex = 2, col ="#00b4ffff")


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