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

########### costom functions associated with the forward modeling ###########

########### fn 10 model P1 and P2 values ###########
forw.m <- function(t, input, a, b, c, R1.int, R2.int){
  if(length(input != t)){
    #warning("Length of input has to be same as t. Vector input is being recycled")
    Rin <- rep(0,t) #Initiate Rin
    Rin <- Rin + input #recycle input vector
  }
  else{
    Rin <- input
  }
  R2 <- rep(0,t) #initiate 
  R1 <- rep(0,t) #initiate 
  if(is.null(R1.int)){
    R1[1] <- input[1]
  }
  if(is.null(R2.int)){
    R2[1] <- input[1]
  }
  else{
    R2[1] <- R2.int
    R1[1] <- R1.int
  }
  for (i in 2:t){ #generate serum and bone series
    R1[i] <- R1[i - 1] + b * (R2[i - 1] - R1[i - 1]) + a * (Rin[i - 1] - R1[i - 1])
    
    R2[i] <- R2[i - 1] + c * (R1[i - 1] - R2[i - 1])
  }
  return(list(R1, R2, input))
}

########### fn 11 simulate sampling grid using width, depth and grid resolution ###########
# all nubers in microns
sim.samp.grid <- function(sample.wd, sample.dp, res){
  # initiate the sampling matrix
  nrow <- round(sample.dp/res)
  
  ncol <- round(sample.wd/res)
  
  sample.grid <- matrix(1, ncol = ncol, nrow = nrow)
  
  ncol.bot <- round(sample.wd/res/2)
  
  # set cutoffs for diagonals
  cutoff <- round(ncol.bot * 0.2)
  
  for(i in nrow : (nrow - ncol.bot)){
    for(j in 1 : ncol){
      x <- nrow - i
      y <- ncol.bot - (cutoff + x)
      
      if(j < y | j > ncol - y){
        sample.grid[i,j] <- 0 
      }
    }
    
  }
  
  return(sample.grid)
}

########### fn 12 simulate enamel block with given parameters###########
# leng = length of enamel block
# thick = thickness of enamel block
# angle = appositional angle
# rate = enamel extension rate
# history = simulated serum Sr ratio history
# sample.wd = width of the sampling groove
# sample.dp = depth of the sampling groove
# sample.dist = dist measurement from the TOP of the crown

sampling.sim <- function(leng, thick, angle, rate,
                         history, sample.wd, sample.dp,
                         sample.dist, res = 10){
  
  total.leng <- leng*2e3 #micron 
  
  total.thick <- (thick/2) * 1e3 #micron
  
  tan.angle <- tan(pi/180*angle)
  
  serum <- history[[1]]
  
  intake <- history[[3]]
  
  n.days <- length(serum)
  
  local.ref <- 0:n.days * rate/res
  # this is at daily resolution
  # so horizontal pixel resolution is higher than daily Sr history
  
  Serum.ref <- c(serum[1], serum)
  
  intake.ref <- c(intake[1], intake)
  
  #for each row of cells, the reference data set will be the serum value at 1:length
  
  cell.loc <- 1:(total.leng/res)
  
  Sr.ref <- approx(x = local.ref, y = Serum.ref, xout = cell.loc)$y
  
  intake.ref <- approx(x = local.ref, y = intake.ref, xout = cell.loc)$y
  
  #simulate appositional angle
  n.shift <- 1/tan.angle #so every vertical cell, shift x towards the back
  
  n.h.en <- total.leng/res
  
  n.v.en <- total.thick/res
  
  Sr.grid <- matrix(0, ncol = n.h.en, nrow = n.v.en)
  
  Sr.grid[n.v.en,] <- Sr.ref
  
  for(i in (n.v.en-1):1){
    x.shift <- round((n.v.en - i) * n.shift)
    
    Sr.aft.temp <- rep(NA, x.shift)
    
    Sr.int.temp <- Sr.ref[(x.shift + 1):n.h.en]
    
    Sr.grid[i,] <- c(Sr.int.temp, Sr.aft.temp)
    
  }
  
  # subset the first half of the matrix for sampling simulations
  en.subset <- 1:(total.leng/res/2)
  
  En.grid <- Sr.grid[,en.subset]
  
  r.Sr <- raster(ncol = ncol(En.grid), nrow = nrow(En.grid),
                 xmn=0, xmx=ncol(En.grid), ymn=0, ymx=nrow(En.grid))
  
  values(r.Sr) <- En.grid
  
  ############ Step 3: build the sampling matrix ######
  #####simulate sample averaging######
  
  # see custom function 11
  sample.grid <- sim.samp.grid(sample.wd*1e3, sample.dp*1e3, res)
  
  #simulate sampling, with gaps between sampling grooves
  
  # mm between sampling grooves
  
  n.samp <- length(sample.dist)
  
  sample.loc <- round(sample.dist * 1e3 / res)
  
  avg.Sr.samp <- rep(0, n.samp) #initiate vector
  
  ############ Step 4: aggregate samples ######
  
  for(i in 1:n.samp){
    
    temp <- En.grid[1:nrow(sample.grid),
                    round(sample.loc[i] - ncol(sample.grid)/2 ):
                      round(sample.loc[i] + ncol(sample.grid)/2 - 1)]
    temp <- replace(temp, is.na(temp), 0)
    temp.prod <- temp * sample.grid 
    avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
    
  }
  
  sim.EDJ <- data.frame(x = en.subset*res/1e3,
                        Sr = Sr.ref[en.subset], Intake = intake.ref[en.subset])
  
  return(list(sim.Sr = avg.Sr.samp,  sim.EDJ = sim.EDJ, sim.En = r.Sr))
  
}

#########Forward model to demonstrate signal attenuation, using synthetic inputs#######

######### experiment 1: Misha's enamel using misha's molar parameters#########

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
# plot(0,0, xlim = c(1,length(syn.input.150.7y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
#      main="150-day ends & 30-day intermediate switches")
# lines(1:length(syn.input.150.7y), syn.input.150.7y)
# lines(1:length(syn.input.150.7y), res150.7[[1]],lwd=2, col = "#00b4ffff")
# legend(200,0.711,c("Intake","Serum"),lwd = c(1,2), col=c("black","#00b4ffff"))


###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

leng <- 100 #mm

thick <- 3 #mm

ap.angle <- 3.5 

# assume that growth rate is constant, at 55.3 microns/day
rate <- 55.3 #100 microns per day

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

history <- res150.7 #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.8 #0.8 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

Misha.sim.dist <- sample.dist

Misha.sim <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)

############ end of experiment 1############

######### experiment 2: Misha's enamel but with shallower grooves #########

###### Step 2: build simulated enamel block using assumptions of enamel maturation#######

history <- res150.7 #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.4 #0.4 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

Misha.sim.dist <- sample.dist

Misha.sim2 <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)


############ experiment 2 Simulate turnover and sampling within sheep enamel ############
Body.mass.misha <- 3000 #kg

Body.mass.sheep <- 30 #kg

exp.scl <- -0.25 #reaction rate scales negatively with body mass

a.m <- 0.0169 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)
b.m <- 0.0141 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)
c.m <- 0.0041 * ((Body.mass.sheep/Body.mass.misha)^exp.scl)

####begin forward model with a winter-summer migration pattern####

syn.input.150 <- c(rep(syn.mid,30), rep(syn.high,150),rep(syn.mid,30),rep(syn.low,150)  )
syn.input.150.3y <- rep(syn.input.150,3)
#4 generate serum and bone series based on input series and turnover parameters#
res150.3 <- forw.m(t = length(syn.input.150.3y), input = syn.input.150.3y, 
                   a = a.m, b = b.m, c = c.m, R1.int = NULL, R2.int = NULL)

# plot input and serum turnover history
# plot(0,0, xlim = c(1,length(syn.input.120.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
#      main="150-day ends & 30-day intermediate switches")
# lines(1:length(syn.input.150.3y), syn.input.150.3y)
# lines(1:length(syn.input.150.3y), res150.3[[1]],lwd=2, col = "#00b4ffff")
# legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

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

history <- res150.3 #simulated serum 87Sr/86Sr history

sample.wd <- 0.8 #0.8 mm wide sampling groove
sample.dp <- 0.4 #0.4 mm sampling depth

sampling.intv <- 2

n.samp <- leng/sampling.intv

sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2

sheep.sim.dist <- sample.dist

sheep.sim <- sampling.sim(leng = leng, thick = thick, 
                          angle = ap.angle, rate = rate,
                          history = history, 
                          sample.wd = sample.wd, sample.dp = sample.dp,
                          sample.dist = sample.dist,
                          res = pix.res)

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
# plot(0,0, xlim = c(1,length(syn.input.60.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
#      main="60-day ends & 30-day intermediate switches")
# lines(1:length(syn.input.60.3y), syn.input.60.3y)
# lines(1:length(syn.input.60.3y), res60.3[[1]],lwd=2, col = "#00b4ffff")
# legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))

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

history <- res60.3 #simulated serum 87Sr/86Sr history

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

############ end of experiment 3############

# ############ experiment 4 Simulate turnover and sampling within sheep enamel ############
# 
# ####begin forward model with more frequent migration pattern####
# # one movement every 1 months, with two months at ends, and one month intermediate
# syn.input.30 <- c(rep(syn.mid,30), rep(syn.high,30),rep(syn.mid, 30),rep(syn.low,30) )
# syn.input.30.3y <- rep(syn.input.30,9)
# #4 generate serum and bone series based on input series and turnover parameters#
# res30.3 <- forw.m(t = length(syn.input.30.3y), input = syn.input.30.3y, 
#                   a = a.m, b = b.m, c = c.m, R1.int = NULL, R2.int = NULL)
# 
# # plot input and serum turnover history
# plot(0,0, xlim = c(1,length(syn.input.30.3y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
#      main="60-day ends & 30-day intermediate switches")
# lines(1:length(syn.input.30.3y), syn.input.30.3y)
# lines(1:length(syn.input.30.3y), res30.3[[1]],lwd=2, col = "#00b4ffff")
# legend(200,0.711,c("Intake","serum"),lwd = c(1,2), col=c("black","#00b4ffff"))
# 
# ###### Step 2: build simulated enamel block using assumptions of enamel maturation#######
# 
# #build an enamel length at 80 mm with redundancy
# leng <- 40 #mm
# 
# # sheep enamel is ~0.7 mm thick, Zazzo et al. 2012
# # build an enamel at 0.4 mm thickness
# thick <- 0.8 #mm
# 
# ap.angle <- 3 #Zazzo et al. 2012
# 
# # assume that growth rate is constant, at 55.3 microns/day
# # growth rate at 35e3 microns in 350 days, Zazzo et al. 2012
# rate <- 35e3/ 350 #100 microns per day
# 
# pix.res <- 10 #micron
# #each vertical 10 microns averages 3 days worth of serum history
# 
# history <- res30.3 #simulated serum 87Sr/86Sr history
# 
# sample.wd <- 0.8 #0.8 mm wide sampling groove
# sample.dp <- 0.4 #0.4 mm sampling depth
# 
# sampling.intv <- 2
# 
# n.samp <- leng/sampling.intv
# 
# sample.dist <-  sampling.intv * 1:n.samp - sampling.intv/2
# 
# sheep.sim3 <- sampling.sim(leng = leng, thick = thick, 
#                            angle = ap.angle, rate = rate,
#                            history = history, 
#                            sample.wd = sample.wd, sample.dp = sample.dp,
#                            sample.dist = sample.dist,
#                            res = pix.res)
# 
# 
# # visualization
# 
# # simulated enamel samples: sheep.sim[[1]]
# plot(sample.dist, sheep.sim3[[1]],
#      ylim = c(syn.low, syn.high))
# # simulated EDJ: sheep.sim[[2]]
# lines(sheep.sim3[[2]]$x, sheep.sim3[[2]]$Sr, lwd=2, col = "#00b4ffff")
# lines(sheep.sim3[[2]]$x, sheep.sim3[[2]]$Intake)


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
