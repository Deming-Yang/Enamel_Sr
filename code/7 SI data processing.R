
############# Time line reconstructions with individual data series in Fig 3D #############

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

#################################################
# prep for Fig S4
# examine the relationships between the resituals between 
# hand drill 87Sr/86Sr and the mixed serum as enamel input

# filter the data and focus on the timeline between 0 and 600

drill.tl.f2 <- filter(drill.tl, drill.tl$tl > 0 & drill.tl$tl < 600)

pred.drill <- approx(x = bin.thin*1:t - 400, y = post.comb.R1.drill.89[[1]]
                     , xout = drill.tl.f2$tl)$y
# prediction: larger depths, lower 87Sr/86Sr, more deviation from predicted value
res.drill <- drill.tl.f2$Sr - pred.drill

res.tib <- tibble(x = drill.tl.f2$depth, y = res.drill)

lm.res <- lm(data = res.tib, y ~ x)
summary(lm.res) # results are significant

# Coefficients:
#               Estimate  Std. Error t value  Pr(>|t|)    
# (Intercept) -0.0007788  0.0002157  -3.611   0.002838 ** 
# x            0.0012311  0.0002577   4.777   0.000295 ***

# calculate 95% confidence intervals for plotting 
newx = seq(min(res.tib$x), max(res.tib$x), by = 0.01)
conf_interval <- predict(lm.res, newdata=data.frame(x=newx), interval="confidence",
                         level = 0.95)

# comparing LA-ICP-MS tusk dentine and molar dentine series

misha.raw <- read.csv("data/Misha M640b raw.csv")

misha.raw.Sr <- misha.raw$Corr..87Sr.86Sr
misha.raw.dist <- misha.raw$dist

##calculate 25 point average to reduce data amount
n.avg <- 25 #now many data points to average
n <- floor(length(misha.raw$Corr..87Sr.86Sr)/n.avg) #99 data points, discard the last <50 data points

#initiate vectors
n.avg.misha.25.sr <- rep(0, n)
n.sd.misha.25.sr <- rep(0, n)

n.avg.misha.25.dist <- rep(0, n)

for(i in 1:n){
  x <- ((i-1)*n.avg + 1):(i*n.avg)
  n.avg.misha.25.sr[i] <- mean(misha.raw.Sr[x])
  n.sd.misha.25.sr[i] <- sd(misha.raw.Sr[x])
  
  n.avg.misha.25.dist[i] <- mean(misha.raw.dist[x])#mean and median are the same
}

#timeline reconstruction with simple growth rate
n.avg.misha.25.tl <- rev(n.avg.misha.25.dist/14.7) #mm/day

#mannually align data
n.avg.misha.25.tl.al <- n.avg.misha.25.tl - 670

