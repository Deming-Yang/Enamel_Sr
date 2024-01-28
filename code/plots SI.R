library(scales)

############# Supplementary plot showing the relative position of the change points ###########
par(mfrow=c(1,2))
plot(cp1.E.rm.x, cp1.E.rm.y, main = "Change point 1", xlim = c(50000,10000),
     xlab = "Distance to reference point (micron)", ylab = "Enamel thickness (micron)")
#this is based on measurements by Uno et al 2020, with an appositional angle at 3.3 +- 0.5 degrees
abline(a = 3.02e+03, b = -tan(3.3/180*pi), lwd = 6, col =alpha("red4", 0.5))
# abline(a = 3.02e+03, b = -tan(3.25/180*pi), lty = 2, col = "red3")
# abline(a = 3.02e+03, b = -tan(3.35/180*pi), lty = 2, col = "red3")
abline(lm.cp1.E.proj.inner,lwd = 2, col = "blue3")
points(cp1.E.rm.x[1:8], cp1.E.rm.y[1:8], pch = 16) 
points(cp1.E.rm.x[9:10], cp1.E.rm.y[9:10], col = "red2", pch = 16) #mark the ones that don't conform to the angle

plot(cp2.E.rm.x, cp2.E.rm.y, main = "Change point 2",xlim = c(55000,15000),
      xlab = "Distance to reference point (micron)", ylab = "Enamel thickness (micron)")
#this is based on measurements by Uno et al 2020, with an appositional angle at 3.3 degrees
abline(a = 3.317e+03, b = -tan(3.3/180*pi), lwd = 6, col = alpha("red4", 0.5))
# abline(a = 3.317e+03, b = -tan(3.25/180*pi), lty = 2, col = "red3")
# abline(a = 3.317e+03, b = -tan(3.35/180*pi), lty = 2, col = "red3")
abline(lm.cp2.E.proj.inner,lwd = 2, col = "blue3")
points(cp2.E.rm.x[1:8], cp2.E.rm.y[1:8], pch = 16) #mark the ones that don't conform to the angle
points(cp2.E.rm.x[9:10], cp2.E.rm.y[9:10], col = "red2", pch = 16) #mark the ones that don't conform to the angle

############# end of change point comparisons ###########

############# Validation, compare micromilled tusk dentine data with corrected LA-ICP-MS data #############
misha.micromill <- read.csv("data/Misha dentin micromill.csv")

misha.micromill<-na.omit(misha.micromill)

#the micromill was done at 500-micron increment, but a geometric conversion is needed (explained in SI)
#here the conversion is done by multiplying 547 (microns)
misha.micromill.dist <- misha.micromill$Position *1.6-5

misha.micromill.sr <- misha.micromill$corr..87Sr.86Sr
misha.micromill.sr.err <- misha.micromill$comb..Err

plot(misha.micromill.dist, misha.micromill.sr,type = "l",col="#00b4ffff",lwd=1.5,
     xlim=c(100,0),ylim=c(0.705,0.712),xlab="Distance from pulp cavity (microns)",ylab="87Sr/86Sr",
     main="Comparing data from the LA-ICP-MS and the solution methods")

#preliminary plot
plot(Drill.no$Dist..From.cervix, Drill.no$X87Sr.86Sr,ylim = c(0.706,0.712),pch=16,col ="red")

lines(84-(1:n), avg.Sr.samp, ylim = c(0.706,0.712),xlim = c(0,100),col = "blue",lwd =2)

#preliminary plot comparing to dentine transect:
#convert dist from crown top to dist from cervix
dent.rm.new.x<- 92-dent.rm.f2$new.x/1e3
proc.Enamel10.rm.new.x<- 92-proc.Enamel10.rm.f$new.x/1e3
proc.Enamel5.rm.new.x<- 92-proc.Enamel5.rm.f$new.x/1e3
proc.Enamel1.rm.new.x<- 92-proc.Enamel1.rm.f$new.x/1e3


plot(dent.rm.new.x, dent.rm.f2$avg, col= alpha("lightcyan4", 0.2),
     pch=16, cex=1, xlim=c(100,0),ylim=c(0.705,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")
points(proc.Enamel1.rm.new.x, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.1),pch=16)
points(proc.Enamel5.rm.new.x, proc.Enamel5.rm.f$avg, col= alpha("orange3", 0.1),pch=16)
points(proc.Enamel10.rm.new.x, proc.Enamel10.rm.f$avg, col= alpha("orange4", 0.1),pch=16)
points(Drill.no$Dist..From.cervix, Drill.no$X87Sr.86Sr, pch=18, cex = 2, col ="red")

#### matching tusk micromill and molar LA-ICP-MS
# molar scaling is a model!
# y = -0.346*log(dist) + 4.478 #alpha or appositional angle, in degrees
# so the enamel extension rate should be cotan()
# top of the crown is 55.3 microns per day
## aling all curves to the switch event as time 0 
dent.rm.new.xt<- 92-dent.rm.f2$new.x/1e3 - 44
proc.Enamel1.rm.new.xt<- 92-proc.Enamel1.rm.f$new.x/1e3 - 44

tan((-0.346*log(proc.Enamel1.rm.new.x) + 4.478)/180*pi) # x in mm

plot(dent.rm.new.xt, dent.rm.f2$avg, col= alpha("lightcyan4", 0.2),
     pch=16, cex=1, xlim=c(45,-65),ylim=c(0.705,0.712),
     xlab="Distance from cervix (mm)",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")
points(proc.Enamel1.rm.new.xt, proc.Enamel1.rm.f$avg, col= alpha("orange", 0.1),pch=16)
misha.micromill.dist <- misha.micromill$Position *2.2-26
points(misha.micromill.dist, misha.micromill.sr,pch=18, cex = 2, col ="#00b4ffff")
# note that the molar sequence is not completely matched up with the tusk


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
misha.micromill.tl.al <- misha.micromill.tl - 325

# plot time line of Sr isotope change
plot(days.cumm.en1.al, proc.Enamel1.rm.f$avg, col= alpha("lightcyan4", 0.2),
     pch=16, cex=1, xlim=c(-400,800),ylim=c(0.705,0.712),
     xlab="Days from moves",
     main = "Conventional vs LAICP-MS",
     ylab = "87Sr/86Sr")
points(misha.micromill.tl.al, misha.micromill.sr,pch=18,cex = 2, col ="#00b4ffff")

