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
syn.input.150.5y <- rep(syn.input.150,7)
#4 generate serum and bone series based on input series and turnover parameters#
res150 <- forw.m(t = length(syn.input.150.5y), input = syn.input.150.5y, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)

#generate Sr intake change of 1 year
syn.input.120 <- c(rep(syn.mid,60), rep(syn.high,120),rep(syn.mid,60),rep(syn.low,120))

# repeat the pattern for multiple years
syn.input.120.5y <- rep(syn.input.120, 7)


#4 generate serum and bone series based on input series and turnover parameters#
res120 <- forw.m(t = length(syn.input.120.5y), input = syn.input.120.5y, a = a, b = b, c = c, R1.int = NULL, R2.int = NULL)


# plot input and serum turnover history
plot(0,0, xlim = c(1,length(syn.input.150.5y)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
     main="150-day ends & 30-day intermediate switches")
lines(1:length(syn.input.150.5y), syn.input.150.5y)
lines(1:length(syn.input.150.5y), res150[[1]],lwd=2, col = "#00b4ffff")
legend(200,0.711,c("Intake","Ivory"),lwd = c(1,2), col=c("black","#00b4ffff"))

# plot(0,0, xlim = c(1,length(syn.input.120)), ylim = c(syn.low, syn.high), xlab = "days", ylab ="87Sr/86Sr",
#      main="120-day ends & 60-day intermediate switches")
# lines(1:length(syn.input.120), syn.input.120)
# lines(1:length(syn.input.120), res120[[1]],lwd=2, col = "#00b4ffff")
# legend(200,0.711,c("Intake","Ivory"),lwd = c(1,2), col=c("black","#00b4ffff"))


###### Step 2: build simulated enamel block using assumptions of enamel maturation#######
# res120[[1]] is the input

# smallest unit of the simulated enamel block: 100 microns/cell
# 100 mm long enamel -> 1000 cell long, 
# 30 mm thick enamel with outer half removed (exclude overprint) -> 15 cell tall

#build an enamel length at 200 mm
total.leng <- 200e3 #micron 

# build an enamel at 1.5 mm thickness
total.thick <- 1.5e3 #micron 

ap.angle <- 3.5

tan.angle <- tan(pi/180*ap.angle)

# misha's enamel geometry: appositional angle: 3.5 degrees, tan = 0.06116262
# so the y-axis shift rate is 55.3 * 0.06116262 = 3.382293 microns
# use a spatial resolution that is 3 times that of 3.382293
# about 10 microns

pix.res <- 10 #micron
#each vertical 10 microns averages 3 days worth of serum history

# assume that growth rate is constant, at 55.3 microns/day
rate <- 55.3

#each cell horizontally represents 3 days 

# then the crown locations of the bottom row of enamel cells (at EDJ)
# would be length * rate

n.days <- length(res150[[1]])

local.ref <- 0:n.days * rate/pix.res
# this is at daily resolution
# so horizontal pixel resolution is higher than daily Sr history

Serum.ref <- c(res150[[1]][1], res150[[1]])

#for each row of cells, the reference data set will be the serum value at 1:length

cell.loc <- 1:(total.leng/pix.res)

Sr.ref <- approx(x = local.ref, y = Serum.ref, xout = cell.loc)$y
#20k cells horizontally

#simulate appositional angle
n.shift <- 1/tan.angle #so every vertical cell, shift x towards the back for 16.35 cells

n.h.en <- total.leng/pix.res

n.v.en <- total.thick/pix.res

# define enamel matrix

Sr.grid <- matrix(0, ncol = n.h.en, nrow = n.v.en)

Sr.grid[n.v.en,] <- Sr.ref

for(i in (n.v.en-1):1){
  x.shift <- round((n.v.en - i) * n.shift)
  
  Sr.aft.temp <- rep(NA, x.shift)
  
  Sr.int.temp <- Sr.ref[(x.shift + 1):n.h.en]
  
  Sr.grid[i,] <- c(Sr.int.temp, Sr.aft.temp)

}

# subset the first 100 mm of the matrix for sampling simulations
en.subset <- 1:10000

En.grid <- Sr.grid[,en.subset]

r.Sr <- raster(ncol = ncol(En.grid), nrow = nrow(En.grid),
               xmn=0, xmx=ncol(En.grid), ymn=0, ymx=nrow(En.grid))

values(r.Sr) <- En.grid

# visualize simulated enamel
plot(r.Sr, ylim = c(0,150))

############ Step 3: build the sampling matrix ######
#####simulate sample averaging######
sample.wd <- 800 #0.8 mm wide sampling groove
sample.dp <- 800 #0.8 mm sampling depth

# see custom function 11
sample.grid <- sim.samp.grid(sample.wd, sample.dp, pix.res)

#simulate complete sampling, no gaps between sampling grooves
n.samp <- round(ncol(En.grid)/ncol(sample.grid)) #simulate complete sampling

avg.Sr.samp <- rep(0, n.samp) #initiate vector

############ Step 4: aggregate samples ######

for(i in 1:n.samp){
  #testing the upper right corner cell of the grid
  if(is.na(En.grid[1,i*ncol(sample.grid)])==T){
    #when there is NA at the upper right corner cell
    j <- sum(!is.na(En.grid[,i*ncol(sample.grid)]))##record the number of valid cells in this column
    if(j > nrow(sample.grid)){#more cells to sample
      #increase depth of the sampling cell
      l <- nrow(En.grid)-j+1
      temp <- En.grid[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid 
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
    }
    else{#not enough cells to sample, sample to the bottom
      
      temp <- En.grid[(nrow(En.grid)-nrow(sample.grid)+1):nrow(En.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
    }
  }else{#y,x
    temp <- En.grid[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
    temp <- replace(temp, is.na(temp), 0)
    temp.prod <- temp * sample.grid 
    avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
  }
}

# visualization
plot(1:length(avg.Sr.samp)*sample.wd/pix.res*1e-3, avg.Sr.samp,
     ylim = c(syn.low, syn.high))
lines(en.subset*1e-3, Sr.ref[en.subset])



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

# ######approach 1 use Kriging silulated enamel block######
# dim(en.Sr.pred.uk) #753 pix long, 28 pix thick
# 
# en.Sr.pred.uk$var1.pred
# 
# sim.en.Sr <- t(en.Sr.pred.uk$var1.pred)
# 
# sample.res <- 100
# 
# ####simulate sample averaging######
# sample.wd <- 1e3 #1mm wide sampling groove
# sample.dp <- 1e3 #1mm wide sampling depth
# 
# sample.grid <- matrix(1,ncol = sample.wd/sample.res, nrow = sample.dp/sample.res)
# 
# sample.grid[7:10,1] <- 0
# sample.grid[9:10,2] <- 0
# sample.grid[10,3] <- 0
# sample.grid[10,4] <- 0
# sample.grid[7:10,10] <- 0
# sample.grid[9:10,9] <- 0
# sample.grid[10,8] <- 0
# sample.grid[10,7] <- 0
# 
# sample.grid #sampling grid
# 
# n <- floor(ncol(sim.en.Sr)/ncol(sample.grid))
# 
# avg.Sr.samp <- rep(0,n) #initiate vector
# 
# for(i in 1:n){#100 samples
#   #testing the upper right corner cell of the grid
#   if(is.na(sim.en.Sr[1,i*ncol(sample.grid)])==T){
#     #when there is NA at the upper right corner cell
#     j <- sum(!is.na(sim.en.Sr[,i*ncol(sample.grid)]))##record the number of valid cells in this column
#     if(j > nrow(sample.grid)){#more cells to sample
#       #increase depth of the sampling cell
#       l <- nrow(sim.en.Sr)-j+1
#       temp <- sim.en.Sr[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid 
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#     }
#     else{#not enough cells to sample, sample to the bottom
#       
#       temp <- sim.en.Sr[(nrow(sim.en.Sr)-nrow(sample.grid)+1):nrow(sim.en.Sr),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
#     }
#   }else{#y,x
#     temp <- sim.en.Sr[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#     temp <- replace(temp, is.na(temp), 0)
#     temp.prod <- temp * sample.grid 
#     avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#   }
# }



# ####simulate sample averaging, with higher depths######
# sample.wd <- 1e3 #1mm wide sampling groove
# sample.dp <- 0.6e3 #1mm wide sampling depth
# 
# make.avg.matrix <- function(sample.wd, sample.dp, sample.res){
#   
#   ncol <- sample.wd/sample.res
#   nrow <- sample.dp/sample.res
#   
#   sample.grid <- matrix(1, ncol, nrow)
#   
#   sample.grid[7:10,1] <- 0
#   sample.grid[7:10,10] <- 0
#   
#   sample.grid[9:10,2] <- 0
#   sample.grid[9:10,9] <- 0
#   
#   sample.grid[10,3] <- 0
#   sample.grid[10,4] <- 0
#   sample.grid[10,8] <- 0
#   sample.grid[10,7] <- 0
#   
# }
# 
# sample.grid <- matrix(1,ncol = sample.wd/sample.res, nrow = sample.dp/sample.res)
# 
# sample.grid #sampling grid
# 
# n <- floor(ncol(sim.en.Sr)/ncol(sample.grid))
# 
# avg.Sr.samp <- rep(0,n) #initiate vector
# 
# for(i in 1:n){#100 samples
#   #testing the upper right corner cell of the grid
#   if(is.na(sim.en.Sr[1,i*ncol(sample.grid)])==T){
#     #when there is NA at the upper right corner cell
#     j <- sum(!is.na(sim.en.Sr[,i*ncol(sample.grid)]))##record the number of valid cells in this column
#     if(j > nrow(sample.grid)){#more cells to sample
#       #increase depth of the sampling cell
#       l <- nrow(sim.en.Sr)-j+1
#       temp <- sim.en.Sr[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid 
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#     }
#     else{#not enough cells to sample, sample to the bottom
#       
#       temp <- sim.en.Sr[(nrow(sim.en.Sr)-nrow(sample.grid)+1):nrow(sim.en.Sr),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#       temp <- replace(temp, is.na(temp), 0)
#       temp.prod <- temp * sample.grid
#       avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
#     }
#   }else{#y,x
#     temp <- sim.en.Sr[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
#     temp <- replace(temp, is.na(temp), 0)
#     temp.prod <- temp * sample.grid 
#     avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
#   }
# }
# 
# plot(Drill.no$Dist..From.cervix, Drill.no$X87Sr.86Sr,ylim = c(0.706,0.712),pch=16,col ="red")
# 
# lines(85-(1:n), avg.Sr.samp, ylim = c(0.706,0.712),xlim = c(0,100),col = "blue",lwd =2)
