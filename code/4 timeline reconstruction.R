
##### compile data from difference sources ####
#filter out 
dent.rm.fa <- filter(dent.rm.f, dent.rm.f$new.x < cp2.D.rm & dent.rm.f$avg < 0.7085)
dent.rm.fb <- filter(dent.rm.f, dent.rm.f$new.x > cp2.D.rm & dent.rm.f$avg > 0.7085)

dent.rm.f2 <- rbind(dent.rm.fa, dent.rm.fb)


### 0 hand drill samples don't need translation: already measured in mm from the cervix ###

### 1 LA-ICP-MS translation ###
# convert dist from crown top to dist from cervix, also from microns to mm
# choose the dentine transet and 3 enamel transects as references
# 88 mm is the top of the crown
# do this for selected transects

dent.rm.new.x <- 88 - dent.rm.f2$new.x/1e3
proc.Enamel1.rm.new.x <- 88 - proc.Enamel1.rm.f$new.x/1e3
proc.Enamel3.rm.new.x <- 88 - proc.Enamel3.rm.f$new.x/1e3
proc.Enamel5.rm.new.x <- 88 - proc.Enamel5.rm.f$new.x/1e3
proc.Enamel7.rm.new.x <- 88 - proc.Enamel7.rm.f$new.x/1e3
proc.Enamel9.rm.new.x <- 88 - proc.Enamel9.rm.f$new.x/1e3
proc.Enamel10.rm.new.x <- 88 - proc.Enamel10.rm.f$new.x/1e3

########## use growth rates to match tusk series and enamel LA-ICP-MS series
# to estimate time line, linear interpolations have to be made based on the growth curve
# assume that the top of the crown is at 100 mm to the cervix
# assume that the enamel growth rate at the top of the crown is 55.3 micron/day
# note that all the calculations depend on the 55.3 micron/day extension rate estimate

# create enamel growth reference grid: 101mm to 0.2 mm 
ref.length.v <- seq(101, 0.2, by = -0.2) # 0.1mm interval

# use appositional angle to get enamel growth rate
# Uno et al. 2020 used a logarithmic relationship
# chech regression of the logarithmic model
Rm3.5.angle$log.dist <- log(Rm3.5.angle$dist) # create log term

logM.ang <- lm(angle ~ log.dist , data = Rm3.5.angle)
summary(logM.ang) # check model fit

# generate pridiction grid
Pred.ang.Rm3.5 <- predict(logM.ang, list(log.dist = log(ref.length.v)))

# also check model fit with log-transformed data
logM.tan <- lm(tan ~ log.dist , data = Rm3.5.angle)
summary(logM.tan)
#         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  0.0781921  0.0022340  35.001 0.000815 ***
#   log.dist    -0.0060427  0.0006361  -9.499 0.010902 * 

# generate prediction grid
Pred.tan.Rm3.5 <- predict(logM.tan, list(log.dist = log(ref.length.v)))

# chekc coefficients
logM.tan$coefficients[[1]] 
logM.tan$coefficients[[2]]

# referenced enamel length at 85.7 mm, 
# for the average enamel extention rate, Uno et al 2020.
ref.en.length <- 85.7 # in mm

# use linear interpolation to approximate growth rate from the reference grid
ref.tan <- approx(x = ref.length.v, y = Pred.tan.Rm3.5, 
                  xout = ref.en.length)$y

# average enamel growth rate 55.3e-3 # mm/day, which is true for most of the crown,
# From Uno et al., 2020
ref.en.growth <- 55.3e-3 # mm/day

# estimate a vector of enamel growth rates at the established length interval
Ref.growth.v <- ref.en.growth * ref.tan / Pred.tan.Rm3.5 # mm/day
# this will be used for estimating the growth rates at sampling locations
# for all enamel LA-ICP-MS transects and drill locations

# linear approximation of measured growth rates from the reference vector
# for selected data based on the molar measurements
# 1. for the dentine trasect
rate.denx <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = dent.rm.new.x) # mm per day

# 2. for selected enamel LA trasect 1, 3, 5, 7, and 10
rate.en1x <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = proc.Enamel1.rm.new.x) 

rate.en3x <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = proc.Enamel3.rm.new.x) 

rate.en5x <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = proc.Enamel5.rm.new.x) 

rate.en7x <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = proc.Enamel7.rm.new.x) 

# for the hand drill samples, x need to be reversed for consistency
Drill.newx <- rev(Drill.no$Dist..From.cervix + 0.05) # to avoid 0 in log calculations

rate.Drill <- approx(x = ref.length.v, y = Ref.growth.v, 
                     xout = Drill.newx) # mm per day

# calculate the differential of the sample lengths
# this will be used to calculate time intervals
# for the dentine trasect
interv.denx <- c(base::diff(-1*dent.rm.new.x)[1], 
                 base::diff(-1*dent.rm.new.x))

# for selected enamel LA trasects: 1, 3, 5, 7, and 10
interv.en1x <- c(base::diff(-1*proc.Enamel1.rm.new.x)[1], 
                 base::diff(-1*proc.Enamel1.rm.new.x))

interv.en3x <- c(base::diff(-1*proc.Enamel3.rm.new.x)[1], 
                 base::diff(-1*proc.Enamel3.rm.new.x))

interv.en5x <- c(base::diff(-1*proc.Enamel5.rm.new.x)[1], 
                 base::diff(-1*proc.Enamel5.rm.new.x))

interv.en7x <- c(base::diff(-1*proc.Enamel7.rm.new.x)[1], 
                 base::diff(-1*proc.Enamel7.rm.new.x))

# for the hand drill samples, the order has to be reversed
interv.Drillx <- c(-1*base::diff(Drill.newx)[1], -1*base::diff(Drill.newx))

# calculate time intervals for all transects and drill samples
time.interv.den <- interv.denx / rate.denx$y
time.interv.en1 <- interv.en1x / rate.en1x$y
time.interv.en3 <- interv.en3x / rate.en3x$y
time.interv.en5 <- interv.en5x / rate.en5x$y
time.interv.en7 <- interv.en7x / rate.en7x$y
time.interv.drill <- interv.Drillx / rate.Drill$y

# calculate cumulative days to build the timeline
days.cumm.den <- base::cumsum(time.interv.den)
days.cumm.en1 <- base::cumsum(time.interv.en1)
days.cumm.en3 <- base::cumsum(time.interv.en3)
days.cumm.en5 <- base::cumsum(time.interv.en5)
days.cumm.en7 <- base::cumsum(time.interv.en7)
days.cumm.drill <- base::cumsum(time.interv.drill)

# align the timeline with the switch, this step was done manually by
# checking that the first abrupt change in 87Sr/86Sr in each series
# is around 0, which is the "reference date" among the timelines

# record the offset numbers for plotting
offset.den <- 890
offset.en1 <- 410
offset.en3 <- 185
offset.en5 <- 155
offset.en7 <- 175
offset.drill <- 535

# apply the offset numbers
days.cumm.den.al <- days.cumm.den - offset.den
days.cumm.en1.al <- days.cumm.en1 - offset.en1 
days.cumm.en3.al <- days.cumm.en3 - offset.en3
days.cumm.en5.al <- days.cumm.en5 - offset.en5
days.cumm.en7.al <- days.cumm.en7 - offset.en7
days.cumm.drill.al <- days.cumm.drill - offset.drill


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

# this is the total length of dentine milled, originally 500 micron interval
# following Yang et al 2023, because of a geometric correction,
# the equivalent sampling interval is 547 microns per day

misha.tusk.micromill.dist <- (maxl.tusk - misha.tusk.micromill$Position) * 547 

# assuming a constant tusk radial growth rate at 14.7 microns/day (Uno et al., 2020)
misha.tusk.micromill.tl <- misha.tusk.micromill.dist/14.7 #mm/day

# align with the switch, this was also done manually
misha.tusk.micromill.tl.al <- misha.tusk.micromill.tl - 365

# a second tusk section was micromilled, which is the M640 slab (Uno et al., 2020)
# This slab has tusk dentine d13C and d18O records
# use a similar framework to reconstruct timeline
M640.tusk <- max(M640.micromill$dist)

# this is the total length of dentine milled, originally measured in microns
M640.micromill.dist <- rev(M640.tusk - M640.micromill$dist) * 1.0934

# assuming a constant tusk radial growth rate at 14.7 microns/day (Uno et al., 2020)
M640.micromill.tl <- M640.micromill.dist/19 #mm/day, Uno 2012

# align with the switch. This was done manually
M640.micromill.tl.al <- M640.micromill.tl - 680

# additional timeline reconstructions for En9 and En10
# calculate rates
rate.en9x <- approx(x = ref.length.v, y = Ref.growth.v, 
                    xout = proc.Enamel9.rm.new.x) 
rate.en10x <- approx(x = ref.length.v, y = Ref.growth.v, 
                     xout = proc.Enamel10.rm.new.x) 

# calculate sample intervals 
interv.en9x <- c(base::diff(-1*proc.Enamel9.rm.new.x)[1], 
                 base::diff(-1*proc.Enamel9.rm.new.x))
interv.en10x <- c(base::diff(-1*proc.Enamel10.rm.new.x)[1], 
                  base::diff(-1*proc.Enamel10.rm.new.x))

# calculate time intervals
time.interv.en9 <- interv.en9x / rate.en9x$y
time.interv.en10 <- interv.en10x / rate.en10x$y

# calculate cumulative days
days.cumm.en9 <- base::cumsum(time.interv.en9)
days.cumm.en10 <- base::cumsum(time.interv.en10)

days.cumm.en9.al <- days.cumm.en9 - 175 
days.cumm.en10.al <- days.cumm.en10 - 175 

## compile data frames for plotting later
dent.tl <- tibble(tl = days.cumm.den.al, Sr = dent.rm.f2$avg, 
                  sd = dent.rm.f2$sd)
en1.tl <- tibble(tl = days.cumm.en1.al, Sr = proc.Enamel1.rm.f$avg, 
                 sd = proc.Enamel1.rm.f$sd)

en3.tl <- tibble(tl = days.cumm.en3.al, Sr = proc.Enamel3.rm.f$avg, 
                 sd = proc.Enamel3.rm.f$sd)

en5.tl <- tibble(tl = days.cumm.en5.al, Sr = proc.Enamel5.rm.f$avg, 
                 sd = proc.Enamel5.rm.f$sd)

en7.tl <- tibble(tl = days.cumm.en7.al, Sr = proc.Enamel7.rm.f$avg, 
                 sd = proc.Enamel7.rm.f$sd)

en9.tl <- tibble(tl = days.cumm.en9.al, Sr = proc.Enamel9.rm.f$avg, 
                 sd = proc.Enamel9.rm.f$sd)
en10.tl <- tibble(tl = days.cumm.en10.al, Sr = proc.Enamel10.rm.f$avg, 
                  sd = proc.Enamel10.rm.f$sd)

Edrill.tl <- tibble(tl = days.cumm.drill.al, 
                   Sr = rev(Drill.no$corr..87Sr.86Sr),
                   sd = rev(Drill.no$comb..Err),
                   depth = rev(Drill.no$Sample.depth)) # sample depth is useful later

Rm3.5b.mill.tl <- tibble(tl = Rm3.5b.mill.tl.al, Sr = rev(Rm3.5b.mill.no$corr..87Sr.86Sr),
                         sd = rev(Rm3.5b.mill.no$comb..Err))

tusk.mill.tl <- tibble(tl = rev(misha.tusk.micromill.tl.al), 
                       Sr = rev(misha.tusk.micromill$corr..87Sr.86Sr),
                       sd = rev(misha.tusk.micromill$comb..Err))
