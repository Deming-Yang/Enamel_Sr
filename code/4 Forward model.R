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

source("code/1 Helper functions.R")

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant molar/Enamel_Sr")

###############forward model for Sr ratio transition at different enamel depths

######build simulated enamel block#######
#enamel geometry: appositional angle: 3.5 degrees, atan = 3/50

#smallest unit of the block: 100 microns -> 1000 pix long, 28 pix thick

pix.res <- 100

total.leng <- 100e3 #micron
total.thick <- 2.8e3 #micron 

####build transitional block with length 25mm
# tr1.leng <- 5e3
# tr2.leng <- 30e3

tr2.leng <- 5e3
tr1.leng <- 30e3

########################part I apposition######################

#####identify starting and end value of Sr ratio

Sr.int <- 0.7117#same as dentine and tusk data

Sr.end <- 0.706  #from ICPMS data, Misha's tusk, Yang et al 2023

#####identify slopes
sl.2 <- -6e-7 #estimated slope 1

sl.1 <- -9e-8 #estimated slope 2

####build the transitional array####
tr1.x <- 1:(tr1.leng/pix.res)
Sr.tr1 <- Sr.int + (tr1.x*pix.res)*sl.1

tr2.x <- 1:(tr2.leng/pix.res)
Sr.tr2 <- Sr.tr1[length(Sr.tr1)] + (tr2.x*pix.res)*sl.2

Sr.tr.all <- c(Sr.tr1, Sr.tr2)
plot(1:length(Sr.tr.all), Sr.tr.all)

#####set location of the initial transition at EDJ####
loc.tr <- 20e3
Sr.int.EDJ <- rep(Sr.int, loc.tr/pix.res)

n.end.EDJ <- (total.leng - (tr1.leng+tr2.leng+loc.tr))/pix.res
Sr.end.EDJ <- rep(Sr.end,n.end.EDJ)

Sr.EDJ.all <- c(Sr.int.EDJ, Sr.tr.all, Sr.end.EDJ)
#this is the reference value, all others will be shifted based on this
plot(1:length(Sr.EDJ.all), Sr.EDJ.all)

#simulate appositional angle
n.shift <- 1/3*50 #so every one pix, shift x towards the back for 16.67 pix

Sr.grid <- matrix(0,ncol = total.leng/pix.res, nrow = total.thick/pix.res)

#mineral density grid
Md.grid <- matrix(0,ncol = total.leng/pix.res, nrow = total.thick/pix.res)

#appositional fraction of mineral density: 0.65
f.appo <- 0.65

Sr.grid[total.thick/pix.res,] <- Sr.EDJ.all #28

Md.grid[total.thick/pix.res,] <- rep(f.appo, ncol(Md.grid))#28

for(i in (total.thick/pix.res -1):1){#28 shifts
  
  Sr.int.temp1 <- rep(NA, (round((total.thick/pix.res-i)*n.shift)))
  
  Md.int.temp1 <- rep(NA, (round((total.thick/pix.res-i)*n.shift)))
  
  Sr.int.temp2 <- rep(Sr.int, loc.tr/pix.res)
  
  Md.int.temp2 <- rep(f.appo, loc.tr/pix.res)
  
  
  n.end.temp <- (total.leng - (tr1.leng+tr2.leng+loc.tr))/pix.res- round((total.thick/pix.res-i)*n.shift)
  Sr.end.temp <- rep(Sr.end,n.end.temp)
  Md.end.temp <- rep(f.appo,n.end.temp)
  
  Sr.grid[i,] <- c(Sr.int.temp1, Sr.int.temp2, Sr.tr.all, Sr.end.temp)
  
  Md.grid[i,] <- c(Md.int.temp1, Md.int.temp2, rep(f.appo, length(Sr.tr.all)), Md.end.temp)
  
}

r.Sr <- raster(ncol = total.leng/pix.res, nrow = total.thick/pix.res,
               xmn=0, xmx=total.leng/pix.res, ymn=0, ymx=total.thick/pix.res)

r.Md <- raster(ncol = total.leng/pix.res, nrow = total.thick/pix.res,
               xmn=0, xmx=total.leng/pix.res, ymn=0, ymx=total.thick/pix.res)

values(r.Sr) <- Sr.grid

plot(r.Sr,ylim=c(0,28))

values(r.Md) <- Md.grid

plot(r.Md,ylim=c(0,28))

########################part II maturation######################
#can use linear interpolation, or loess with predictions
syn.dat <- data.frame(index = c(1:ncol(Sr.grid)), input = rev(Sr.grid[total.thick/pix.res,]))
#synthetic hair data, in distance from root

#use loess function to smooth and downsample the longer vector
plot(1:ncol(Sr.grid), syn.dat$input, type="l")
loess.Sr <- loess(input ~ index, data=syn.dat, span=0.02)
lines(syn.dat$index, predict(loess.Sr),col="red")

#use loess.Sr to rescale x axis to time
#enamel growth: 55.3 microns/day, resolution is 100 microns, so one pixel is around 2 days

n.days <- round(ncol(Sr.grid)/0.553) #n days represented in enamel 1808 days


input.Sr <- predict(loess.Sr, newdata = 1:n.days*ncol(Sr.grid)/n.days)
input.Sr[1] <- input.Sr[2]

#obtain timeline of Sr input into enamel
tl <- data.frame(index = 1:n.days, input = input.Sr)

plot(tl$index, tl$input, type="l")

#then model enamel formation
#apposition is done

#maturation: angle, length, and weight distribution
#iteration 1: angle is the same as apposition, immediately starting
#duration is evaluating using lm =700 pxiels or 70mm/ 70e3 microns
#Duration to weight distribution: 70e3/55.3

days.to.add <- round(70e3/55.3) #1266 days before start of the experiment!

#need mineral density map!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
for(i in (total.thick/pix.res -1):1){#28 shifts

Sr.int.temp1 <- rep(NA, (round((total.thick/pix.res-i)*n.shift)))

Sr.int.temp2 <- rep(Sr.int, loc.tr/pix.res)

n.end.temp <- (total.leng - (tr1.leng+tr2.leng+loc.tr))/pix.res- round((total.thick/pix.res-i)*n.shift)
Sr.end.temp <- rep(Sr.end,n.end.temp)

Sr.grid[i,] <- c(Sr.int.temp1, Sr.int.temp2, Sr.tr.all, Sr.end.temp)

}
input.Sr.ext <- c(rep(input.Sr[1],days.to.add),input.Sr)

#extended timeline to accommodate maturation averaging
tl.ext <- data.frame(index = 1:(n.days+days.to.add), input = input.Sr.ext)

plot(tl.ext$index, tl.ext$input, type="l")

n.mat <- n.days+days.to.add


#also assume evenly distributed weights? Yes


#####simulate sample averaging######
sample.wd <- 1e3 #1mm wide sampling groove
sample.dp <- 1e3 #1mm wide sampling depth

sample.grid <- matrix(1,ncol = sample.wd/pix.res, nrow = sample.dp/pix.res)

sample.grid[7:10,1] <- 0
sample.grid[9:10,2] <- 0
sample.grid[10,3] <- 0
sample.grid[10,4] <- 0
sample.grid[7:10,10] <- 0
sample.grid[9:10,9] <- 0
sample.grid[10,8] <- 0
sample.grid[10,7] <- 0

sample.grid #sampling grid

n <- (ncol(Sr.grid)/ncol(sample.grid))

avg.Sr.samp <- rep(0,n) #initiate vector

for(i in 1:n){#100 samples
  #testing the upper right corner cell of the grid
  if(is.na(Sr.grid[1,i*ncol(sample.grid)])==T){
    #when there is NA at the upper right corner cell
    j <- sum(!is.na(Sr.grid[,i*ncol(sample.grid)]))##record the number of valid cells in this column
    if(j > nrow(sample.grid)){#more cells to sample
      #increase depth of the sampling cell
      l <- nrow(Sr.grid)-j+1
      temp <- Sr.grid[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid 
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
    }
    else{#not enough cells to sample, sample to the bottom
      
      temp <- Sr.grid[(nrow(Sr.grid)-nrow(sample.grid)+1):nrow(Sr.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
    }
  }else{#y,x
    temp <- Sr.grid[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
    temp <- replace(temp, is.na(temp), 0)
    temp.prod <- temp * sample.grid 
    avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
  }
}

plot((1:n)*10-5, avg.Sr.samp, ylim = c(0.706,0.712))
lines(1:ncol(Sr.grid), Sr.grid[1,])
lines(1:ncol(Sr.grid), Sr.grid[nrow(Sr.grid),])

##############try a different sampling matrix##########
#####simulate sample averaging######
sample.wd <- 1e3 #1mm wide sampling groove
sample.dp <- 1.5e3 #1mm wide sampling depth

sample.grid <- matrix(1,ncol = sample.wd/pix.res, nrow = sample.dp/pix.res)

#thik about ways to generate this matrix
sample.grid[(nrow(sample.grid)-3):nrow(sample.grid),1] <- 0
sample.grid[(nrow(sample.grid)-1):nrow(sample.grid),2] <- 0
sample.grid[nrow(sample.grid),3] <- 0
sample.grid[nrow(sample.grid),4] <- 0
sample.grid[(nrow(sample.grid)-3):nrow(sample.grid),10] <- 0
sample.grid[(nrow(sample.grid)-1):nrow(sample.grid),9] <- 0
sample.grid[nrow(sample.grid),8] <- 0
sample.grid[nrow(sample.grid),7] <- 0

sample.grid #sampling grid

n <- (ncol(Sr.grid)/ncol(sample.grid))

avg.Sr.samp <- rep(0,n) #initiate vector

for(i in 1:n){#100 samples
  #testing the upper right corner cell of the grid
  if(is.na(Sr.grid[1,i*ncol(sample.grid)])==T){
    #when there is NA at the upper right corner cell
    j <- sum(!is.na(Sr.grid[,i*ncol(sample.grid)]))##record the number of valid cells in this column
    if(j > nrow(sample.grid)){#more cells to sample
      #increase depth of the sampling cell
      l <- nrow(Sr.grid)-j+1
      temp <- Sr.grid[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid 
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
    }
    else{#not enough cells to sample, sample to the bottom

      temp <- Sr.grid[(nrow(Sr.grid)-nrow(sample.grid)+1):nrow(Sr.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
    }
  }else{#y,x
    temp <- Sr.grid[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
    temp <- replace(temp, is.na(temp), 0)
    temp.prod <- temp * sample.grid 
    avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
  }
}

plot((1:n)*10-5, avg.Sr.samp, ylim = c(0.706,0.712))
lines(1:ncol(Sr.grid), Sr.grid[1,])
lines(1:ncol(Sr.grid), Sr.grid[nrow(Sr.grid),])

##############try a different sampling matrix##########
#####simulate sample averaging######
sample.wd <- 1e3 #1mm wide sampling groove
sample.dp <- 2.8e3 #2.8mm sampling depth

sample.grid <- matrix(1,ncol = sample.wd/pix.res, nrow = sample.dp/pix.res)

#thik about ways to generate this matrix
sample.grid[(nrow(sample.grid)-3):nrow(sample.grid),1] <- 0
sample.grid[(nrow(sample.grid)-1):nrow(sample.grid),2] <- 0
sample.grid[nrow(sample.grid),3] <- 0
sample.grid[nrow(sample.grid),4] <- 0
sample.grid[(nrow(sample.grid)-3):nrow(sample.grid),10] <- 0
sample.grid[(nrow(sample.grid)-1):nrow(sample.grid),9] <- 0
sample.grid[nrow(sample.grid),8] <- 0
sample.grid[nrow(sample.grid),7] <- 0

sample.grid #sampling grid

n <- (ncol(Sr.grid)/ncol(sample.grid))

avg.Sr.samp <- rep(0,n) #initiate vector

for(i in 1:n){#100 samples
  #testing the upper right corner cell of the grid
  if(is.na(Sr.grid[1,i*ncol(sample.grid)])==T){
    #when there is NA at the upper right corner cell
    j <- sum(!is.na(Sr.grid[,i*ncol(sample.grid)]))##record the number of valid cells in this column
    if(j > nrow(sample.grid)){#more cells to sample
      #increase depth of the sampling cell
      l <- nrow(Sr.grid)-j+1
      temp <- Sr.grid[l:(l+nrow(sample.grid)-1),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid 
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
    }
    else{#not enough cells to sample, sample to the bottom
      
      temp <- Sr.grid[(nrow(Sr.grid)-nrow(sample.grid)+1):nrow(Sr.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
      temp <- replace(temp, is.na(temp), 0)
      temp.prod <- temp * sample.grid
      avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))#here the sample grid is incorrect!
    }
  }else{#y,x
    temp <- Sr.grid[1:nrow(sample.grid),((i-1)*ncol(sample.grid)+1):(i*ncol(sample.grid))]
    temp <- replace(temp, is.na(temp), 0)
    temp.prod <- temp * sample.grid 
    avg.Sr.samp[i]<- sum(temp.prod)/sum(as.integer(temp.prod>0))
  }
}

plot((1:n)*10-5, avg.Sr.samp, ylim = c(0.706,0.712))
lines(1:ncol(Sr.grid), Sr.grid[1,])
lines(1:ncol(Sr.grid), Sr.grid[nrow(Sr.grid),])

#####add concentration matrix####



###think about rounding of corners

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


