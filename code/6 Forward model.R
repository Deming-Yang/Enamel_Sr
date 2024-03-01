library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(segmented)
library(raster)

source("./code/1 Helper functions.R")

#########Forward model to demonstrate signal attenuation, using synthetic inputs#######

# experiment 1: Misha's enamel using misha's molar parameters

#########Step 1: create synthetic intake and serum histories #########

####begin forward model####
#1 extract turnover parameters, use MAP in the forward model
a <- 0.0169 #reaction rate estimates from Yang et al. 2023
b <- 0.0141
c <- 0.0041

#2 number of days in the simulation
#t <- 700 #omitted here due to synthetic input data 

#3 generate input series with fixed duration#
#Mid -> low -> mid -> High -> Mid
syn.mid <- 0.709
syn.high <- 0.711
syn.low <- 0.707

syn.input.150 <- c(rep(syn.mid,30), rep(syn.high,150),rep(syn.mid,30),rep(syn.low,150)  )
syn.input.150.7y <- rep(syn.input.150,7)
#4 generate serum and bone series based on input series and turnover parameters#
res150.7 <- forw.m(t = length(syn.input.150.7y), input = syn.input.150.7y, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)

# plot input and serum turnover history
plot(0,0, xlim = c(1,length(syn.input.150.7y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="150-day ends & 30-day intermediate switches")
lines(1:length(syn.input.150.7y), syn.input.150.7y)
lines(1:length(syn.input.150.7y), res150.7[[1]],lwd=2, col = "#00b4ffff")
legend(200,0.711,c("Intake","Serum"),lwd = c(1,2), col=c("black","#00b4ffff"))


###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

leng <- 100 #mm

thick <- 3 #mm

ap.angle <- 3.5 

# assume that growth rate is constant, at 55.3 microns/day
rate <- 55.3 #100 microns per day

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

history <- res150.7[[1]] #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.8 #0.8 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

Misha.sim <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)

# visualization
# simulated enamel samples: sheep.sim[[1]]
plot(sample.dist, Misha.sim[[1]],
     ylim = c(syn.low, syn.high))
# simulated EDJ: sheep.sim[[2]]
lines(Misha.sim[[2]]$x, Misha.sim[[2]]$Sr)
abline(h = syn.low, lty = 2)
abline(h = syn.high, lty = 2)

# # simulated enamel block
# plot(Misha.sim[[3]]) 

############ end of experiment 1############


############ experiment 2 Simulate turnover and sampling within sheep enamel ############
Body.mass.misha <- 3000 #kg

Body.mass.sheep <- 30 #kg

exp.scl <- -0.25 #reaction rate scales negatively with body mass

a.m <- 0.0169 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)
b.m <- 0.0141 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)
c.m <- 0.0041 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)

####begin forward model with a winter-summer migration pattern####

syn.input.120 <- c(rep(syn.mid,60), rep(syn.high,120),rep(syn.mid,60),rep(syn.low,120)  )
syn.input.120.3y <- rep(syn.input.120,3)
#4 generate serum and bone series based on input series and turnover parameters#
res120.3 <- forw.m(t = length(syn.input.120.3y), input = syn.input.120.3y, 
                   a = a.m, b = b.m, c = c.m, R1.int = NULL, R2.int = NULL)

# plot input and serum turnover history
plot(0,0, xlim = c(1,length(syn.input.120.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="150-day ends & 30-day intermediate switches")
lines(1:length(syn.input.120.3y), syn.input.120.3y)
lines(1:length(syn.input.120.3y), res120.3[[1]],lwd=2, col = "#00b4ffff")
legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

#build an enamel length at 80 mm with redundancy
leng <- 40 #mm

# sheep enamel is ~0.7 mm thick, Zazzo et al. 2012
# build an enamel at 0.4 mm thickness
thick <- 0.8 #mm

ap.angle <- 3 #Zazzo et al. 2012

# assume that growth rate is constant, at 55.3 microns/day
# growth rate at 35e3 microns in 350 days, Zazzo et al. 2012
rate <- 35e3/ 350 #100 microns per day

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

history <- res120.3[[1]] #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.4 #0.4 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

sheep.sim <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)


# visualization

# simulated enamel samples: sheep.sim[[1]]
plot(sample.dist, sheep.sim[[1]],
     ylim = c(syn.low, syn.high))
# simulated EDJ: sheep.sim[[2]]
lines(sheep.sim[[2]]$x, sheep.sim[[2]]$Sr)
abline(h = syn.low, lty = 2)
abline(h = syn.high, lty = 2)

############ end of experiment 2############

############ experiment 3 Simulate turnover and sampling within sheep enamel ############

####begin forward model with more frequent migration pattern####
# one movement every 3 months, with two months at ends, and one month intermediate
syn.input.60 <- c(rep(syn.mid,30), rep(syn.high,60),rep(syn.mid, 30),rep(syn.low,60) )
syn.input.60.3y <- rep(syn.input.60,6)
#4 generate serum and bone series based on input series and turnover parameters#
res60.3 <- forw.m(t = length(syn.input.60.3y), input = syn.input.60.3y, 
                   a = a.m, b = b.m, c = c.m, R1.int = NULL, R2.int = NULL)

# plot input and serum turnover history
plot(0,0, xlim = c(1,length(syn.input.60.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="60-day ends & 30-day intermediate switches")
lines(1:length(syn.input.60.3y), syn.input.60.3y)
lines(1:length(syn.input.60.3y), res60.3[[1]],lwd=2, col = "#00b4ffff")
legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

#build an enamel length at 80 mm with redundancy
leng <- 40 #mm

# sheep enamel is ~0.7 mm thick, Zazzo et al. 2012
# build an enamel at 0.4 mm thickness
thick <- 0.8 #mm

ap.angle <- 3 #Zazzo et al. 2012

# assume that growth rate is constant, at 55.3 microns/day
# growth rate at 35e3 microns in 350 days, Zazzo et al. 2012
rate <- 35e3/ 350 #100 microns per day

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

history <- res60.3[[1]] #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.4 #0.4 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

sheep.sim2 <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)


# visualization

# simulated enamel samples: sheep.sim[[1]]
plot(sample.dist, sheep.sim2[[1]],
     ylim = c(syn.low, syn.high))
# simulated EDJ: sheep.sim[[2]]
lines(sheep.sim2[[2]]$x, sheep.sim2[[2]]$Sr)
abline(h = syn.low, lty = 2)
abline(h = syn.high, lty = 2)

############ end of experiment 3############

############ experiment 4 Simulate turnover and sampling within sheep enamel ############

####begin forward model with more frequent migration pattern####
# one movement every 1 months, with two months at ends, and one month intermediate
syn.input.30 <- c(rep(syn.mid,30), rep(syn.high,30),rep(syn.mid, 30),rep(syn.low,30) )
syn.input.30.3y <- rep(syn.input.30,9)
#4 generate serum and bone series based on input series and turnover parameters#
res30.3 <- forw.m(t = length(syn.input.30.3y), input = syn.input.30.3y, 
                  a = a.m, b = b.m, c = c.m, R1.int = NULL, R2.int = NULL)

# plot input and serum turnover history
plot(0,0, xlim = c(1,length(syn.input.30.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="60-day ends & 30-day intermediate switches")
lines(1:length(syn.input.30.3y), syn.input.30.3y)
lines(1:length(syn.input.30.3y), res30.3[[1]],lwd=2, col = "#00b4ffff")
legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

#build an enamel length at 80 mm with redundancy
leng <- 40 #mm

# sheep enamel is ~0.7 mm thick, Zazzo et al. 2012
# build an enamel at 0.4 mm thickness
thick <- 0.8 #mm

ap.angle <- 3 #Zazzo et al. 2012

# assume that growth rate is constant, at 55.3 microns/day
# growth rate at 35e3 microns in 350 days, Zazzo et al. 2012
rate <- 35e3/ 350 #100 microns per day

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

history <- res30.3[[1]] #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.4 #0.4 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

sheep.sim3 <- sampling.sim(leng = leng, thick = thick, 
                           angle = ap.angle, rate = rate,
                           history = history, 
                           sample.wd = sample.wd, sample.dp = sample.dp,
                           sample.dist = sample.dist,
                           res = pix.res)


# visualization

# simulated enamel samples: sheep.sim[[1]]
plot(sample.dist, sheep.sim3[[1]],
     ylim = c(syn.low, syn.high))
# simulated EDJ: sheep.sim[[2]]
lines(sheep.sim3[[2]]$x, sheep.sim3[[2]]$Sr)
abline(h = syn.low, lty = 2)
abline(h = syn.high, lty = 2)


# ggplot() + 
#   geom_tile(data=r.Sr, aes(x=x, y=y, fill=value))+
#   scale_fill_viridis() +
#   theme(legend.position="bottom") +
#   theme(legend.key.width=unit(2, "cm"))
# 
# 
# 
# ##plot the Sr grid###
# sf.all.dat <- st_as_sf(all.dat,  agr = NA_agr_,
#                        coords = c("y","x"),
#                        dim = "XYZ",
#                        remove = TRUE,
#                        na.fail = TRUE,
#                        sf_column_name = NULL)


###############forward model for Sr ratio transition at different enamel depths
# 
# sample.grid #sampling grid
# 
# n <- floor(ncol(sim.en.Sr)/ncol(sample.grid))
# 
# avg.Sr.samp <- rep(0,n) #initiate vector
# 
# for(i in 1:n.samp){
#   #testing the upper right corner cell of the grid
#   if(is.na(En.grid[1,i*ncol(sample.grid)])==T){
#     #when there is NA at the upper right corner cell
#     j <- sum(!is.na(En.grid[,i*ncol(sample.grid)]))##record the number of valid cells in this column
#     if(j > nrow(sample.grid)){#more cells to sample
#       #increase depth of the sampling cell
#       l <- nrow(En.grid)-j+1
#       temp <- En.grid[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid 
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#     }
#     else{#not enough cells to sample, sample to the bottom
#       
#       temp <- En.grid[(nrow(En.grid)-nrow(sample.grid)+1):nrow(En.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
#     }
#   }else{#y,x
#     temp <- En.grid[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#     temp <- replace(temp, is.na(temp), 0)
#     temp.prod <- temp * sample.grid 
#     avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#   }
# }
