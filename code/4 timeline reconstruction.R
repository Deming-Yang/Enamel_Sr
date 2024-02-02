
#filter out 
dent.rm.fa <- filter(dent.rm.f, dent.rm.f$new.x < cp2.D.rm & dent.rm.f$avg < 0.7085)
dent.rm.fb <- filter(dent.rm.f, dent.rm.f$new.x > cp2.D.rm & dent.rm.f$avg > 0.7085)

dent.rm.f2 <- rbind(dent.rm.fa, dent.rm.fb)


### 0 hand drill samples don't need translation: already measured in mm from the cervix ###

### 1 LA-ICP-MS translation ###
# convert dist from crown top to dist from cervix, also from microns to mm
# choose the dentine transet and 3 enamel transects as references
# 88 mm is the top of the crown

dent.rm.new.x <- 88 - dent.rm.f2$new.x/1e3
proc.Enamel1.rm.new.x <- 88 - proc.Enamel1.rm.f$new.x/1e3
proc.Enamel5.rm.new.x <- 88 - proc.Enamel5.rm.f$new.x/1e3
proc.Enamel9.rm.new.x <- 88 - proc.Enamel9.rm.f$new.x/1e3
proc.Enamel10.rm.new.x <- 88 - proc.Enamel10.rm.f$new.x/1e3

# adding CA and UT Sr posterior from Yang et al. 2023
CA.Sr <- 0.70597
UT.Sr <- 0.71115

# posterior CA Sr CI 
CA.Sr.sd <- (0.706243 - CA.Sr)/2

UT.Sr.sd <- (0.7120174 - UT.Sr)/2

########## use growth rates to match tusk series and enamel LA-ICP-MS series
# to estimate time line, linear interpolations have to be made based on the growth curve
# assume that the top of the crown is at 100 mm to the cervix
# assume that the enamel growth rate at the top of the crown is 55.3 micron/day
# note that all the calculations depend on the 55.3 micron/day extension rate estimate

# create reference grid: 101mm to 0 mm 
ref.length.v <- seq(101,0, by = -0.1) # 0.1mm interval

# use appositional angle to get enamel growth rate
# Uno et al. 2020 used a logarithmic relationship
# chech regression of the logarithmic model
Rm3.5.angle$log.dist <- log(Rm3.5.angle$dist) # create log term

logM.tan <- lm(tan ~ log.dist , data = Rm3.5.angle)
summary(logM.tan)
#         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0781921  0.0022340  35.001 0.000815 ***
#   log.dist    -0.0060427  0.0006361  -9.499 0.010902 * 

Pred.tan.Rm3.5 <- predict(logM.tan, list(log.dist = log(ref.length.v)))

logM.tan$coefficients[[1]]
logM.tan$coefficients[[2]]

# visualize the model fit
plot(Rm3.5.angle$dist, Rm3.5.angle$tan, ylim = c(0.04,0.09), xlim = c(0,100),
     ylab = "Tangent(alpha)", xlab ="Distance from cervix (mm)")
lines(ref.length.v, Pred.tan.Rm3.5)

# referenced enamel length at 52.7 mm
ref.en.length <- 85.7 # in mm

ref.tan <- approx(x = ref.length.v, y = Pred.tan.Rm3.5, xout = ref.en.length)$y

# average enamel growth rate 55.3e-3 # mm/day, which is true for most of the crown
ref.en.growth <- 55.3e-3 # mm/day

# estimate a vector of enamel growth rates at the established length interval
Ref.growth.v <- ref.en.growth * ref.tan / Pred.tan.Rm3.5 # mm/day

# linear approximation of measured growth rates from the reference vector
# for selected data based on the molar measurements

# 1. for the dentine trasect
rate.denx <- approx(x = ref.length.v, y = Ref.growth.v, xout = dent.rm.new.x) # mm per day

# 2. for selected enamel LA trasect 1
rate.en1x <- approx(x = ref.length.v, y = Ref.growth.v, xout = proc.Enamel1.rm.new.x) 

# for the hand drill samples, x need to be reversed
Drill.newx <- rev(Drill.no$Dist..From.cervix + 0.05) # to avoid 0 in log calculations

rate.Drill <- approx(x = ref.length.v, y = Ref.growth.v, xout = Drill.newx) # mm per day

# calculate the differential of the sample lengths
# for the dentine trasect
interv.denx <- c(base::diff(-1*dent.rm.new.x)[1], base::diff(-1*dent.rm.new.x))

# for selected enamel LA trasects: 1, 5, and 10
interv.en1x <- c(base::diff(-1*proc.Enamel1.rm.new.x)[1], base::diff(-1*proc.Enamel1.rm.new.x))

# for the hand drill samples, the order has to be reversed
interv.Drillx <- c(-1*base::diff(Drill.newx)[1], -1*base::diff(Drill.newx))

# calculate time intervals
time.interv.den <- interv.denx / rate.denx$y
time.interv.en1 <- interv.en1x / rate.en1x$y
time.interv.drill <- interv.Drillx / rate.Drill$y

# calculate cumulative days
days.cumm.den <- base::cumsum(time.interv.den)
days.cumm.en1 <- base::cumsum(time.interv.en1)
days.cumm.drill <- base::cumsum(time.interv.drill)

# align with the switch, this step was done manually
days.cumm.den.al <- days.cumm.den - 910
days.cumm.en1.al <- days.cumm.en1 - 410 
days.cumm.drill.al <- days.cumm.drill - 530

### 3 enamel micromill translation ###
# transform this into time, not vertical distance, use trigonometry
# micromill angle is approx. 3.2 degrees (Uno et al., 2020)
# enamel daily secretion rate is ~3.5 microns/day (Uno et al., 2020)

Rm3.5b.mill.tl <- rev(Rm3.5b.mill.no$Dist..From.EDJ/3.5) # in days
Rm3.5b.mill.tl.al <- Rm3.5b.mill.tl - 120

# use micromill tusk dentine data in the comparison
# in mm, micromill was done at 0.5 mm interval
# from the pulp cavity, which is the newest dentine
maxl.tusk <- max(misha.tusk.micromill$Position)

#this is the total length of dentine milled, 500 micron interval
misha.tusk.micromill.dist <- (maxl.tusk - misha.tusk.micromill$Position) * 500 

# assuming a constant tusk growth rate at 14.7 microns/day (Uno et al., 2020)
misha.tusk.micromill.tl <- misha.tusk.micromill.dist/14.7 #mm/day

# align with the switch
misha.tusk.micromill.tl.al <- misha.tusk.micromill.tl - 335

# additional timeline reconstructions for En9 and En10

# calculate rates
rate.en9x <- approx(x = ref.length.v, y = Ref.growth.v, xout = proc.Enamel9.rm.new.x) 
rate.en10x <- approx(x = ref.length.v, y = Ref.growth.v, xout = proc.Enamel10.rm.new.x) 

# calculate sample intervals 
interv.en9x <- c(base::diff(-1*proc.Enamel9.rm.new.x)[1], base::diff(-1*proc.Enamel9.rm.new.x))
interv.en10x <- c(base::diff(-1*proc.Enamel10.rm.new.x)[1], base::diff(-1*proc.Enamel10.rm.new.x))

# calculate time intervals
time.interv.en9 <- interv.en9x / rate.en9x$y
time.interv.en10 <- interv.en10x / rate.en10x$y

# calculate cumulative days
days.cumm.en9 <- base::cumsum(time.interv.en9)
days.cumm.en10 <- base::cumsum(time.interv.en10)

days.cumm.en9.al <- days.cumm.en9 - 175 
days.cumm.en10.al <- days.cumm.en10 - 175 