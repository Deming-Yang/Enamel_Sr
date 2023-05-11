library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)

source("code/1 Helper functions.R")

setwd("C:/Users/ydmag/Google Drive/U of U/Elephant molar/Enamel_Sr")

EDJ.dentin1 <- read.csv("data/EDJ_dentin1.csv")
EDJ.dentin2 <- read.csv("data/EDJ_dentin2.csv")
Enamel1 <- read.csv("data/Enamel1.csv")
Enamel2 <- read.csv("data/Enamel2.csv")
Enamel3 <- read.csv("data/Enamel3.csv")
Enamel4 <- read.csv("data/Enamel4.csv")
Enamel5 <- read.csv("data/Enamel5.csv")
Enamel6 <- read.csv("data/Enamel6.csv")
Enamel6ext <- read.csv("data/Enamel6ext.csv")
Enamel7 <- read.csv("data/Enamel7.csv")

pl.EDJ.dentin1 <- read.csv("data/pl_EDJ_dentin1.csv")
pl.EDJ.dentin2 <- read.csv("data/pl_EDJ_dentin2.csv")
pl.Enamel1 <- read.csv("data/pl_Enamel1.csv")
pl.Enamel2 <- read.csv("data/pl_Enamel2.csv")
pl.Enamel3 <- read.csv("data/pl_Enamel3.csv")
pl.Enamel4 <- read.csv("data/pl_Enamel4.csv")
pl.Enamel5 <- read.csv("data/pl_Enamel5.csv")
pl.Enamel6 <- read.csv("data/pl_Enamel6.csv")
pl.Enamel6ext <- read.csv("data/pl_Enamel6ext.csv")
pl.Enamel7 <- read.csv("data/pl_Enamel7.csv")

############data processing using custom function###########

#dentine at EDJ
proc.EDJ.dentin1 <- data.process(EDJ.dentin1, pl.EDJ.dentin1)
proc.EDJ.dentin1
proc.EDJ.dentin1.rm <- moving.avg(proc.EDJ.dentin1, 25)

proc.EDJ.dentin2 <- data.process(EDJ.dentin2, pl.EDJ.dentin2)
proc.EDJ.dentin2
proc.EDJ.dentin2.rm <- moving.avg(proc.EDJ.dentin2, 25)

#enamel at EDJ
proc.Enamel1 <- data.process(Enamel1, pl.Enamel1)
proc.Enamel1
proc.Enamel1.rm <- moving.avg(proc.Enamel1, 25)

#enamel transect 2
proc.Enamel2 <- data.process(Enamel2, pl.Enamel2)
proc.Enamel2
proc.Enamel2.rm <- moving.avg(proc.Enamel2, 25)

#enamel transect 3
proc.Enamel3 <- data.process(Enamel3, pl.Enamel3)
proc.Enamel3
proc.Enamel3.rm <- moving.avg(proc.Enamel3, 25)

#enamel transect 4
proc.Enamel4 <- data.process(Enamel4, pl.Enamel4)
proc.Enamel4
proc.Enamel4.rm <- moving.avg(proc.Enamel4, 25)

#enamel transect 5
proc.Enamel5 <- data.process(Enamel5, pl.Enamel5)
proc.Enamel5
proc.Enamel5.rm <- moving.avg(proc.Enamel5, 25)

#enamel transect 6
proc.Enamel6 <- data.process(Enamel6, pl.Enamel6)
proc.Enamel6
proc.Enamel6.rm <- moving.avg(proc.Enamel6, 25)


#enamel transect 6
proc.Enamel6ext <- data.process(Enamel6ext, pl.Enamel6ext)
proc.Enamel6ext
proc.Enamel6ext.rm <- moving.avg(proc.Enamel6ext, 25)

proc.Enamel7 <- data.process(Enamel7, pl.Enamel7)
proc.Enamel7
proc.Enamel7.rm <- moving.avg(proc.Enamel7, 25)
#####use the sf package to plot those values

all.dat <- rbind(proc.EDJ.dentin1, proc.EDJ.dentin2, proc.Enamel1, proc.Enamel2, proc.Enamel3,
                 proc.Enamel4, proc.Enamel5, proc.Enamel6, proc.Enamel6ext, proc.Enamel7)

all.dat.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm, proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
      proc.Enamel4.rm, proc.Enamel5.rm, proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel7.rm)

all.enamel.rm <- rbind(proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
                    proc.Enamel4.rm, proc.Enamel5.rm, proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel7.rm)
all.dentine.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)

require(sf)
sf.all.dat <- st_as_sf(all.dat,  agr = NA_agr_,
                       coords = c("x","y"),
                       dim = "XYZ",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)

sf.all.dat.rm <- st_as_sf(all.dat.rm,  agr = NA_agr_,
                       coords = c("x","y"),
                       dim = "XYZ",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)

sf.all.enamel.rm <- st_as_sf(all.enamel.rm,  agr = NA_agr_,
                          coords = c("x","y"),
                          dim = "XYZ",
                          remove = TRUE,
                          na.fail = TRUE,
                          sf_column_name = NULL)

sf.all.dentine.rm <- st_as_sf(all.dentine.rm,  agr = NA_agr_,
                             coords = c("x","y"),
                             dim = "XYZ",
                             remove = TRUE,
                             na.fail = TRUE,
                             sf_column_name = NULL)


proc.Enamel6.rm

plot(sf.all.dat, breaks = seq(0.703, 0.715,by = 0.001), pch=16)

plot(sf.all.dat.rm, breaks = seq(0.700, 0.716,by = 0.001), pch=16)

plot(sf.all.enamel.rm, breaks = seq(0.700, 0.716,by = 0.001), pch=16)

plot(sf.all.dentine.rm,add = T)

# plot extent
plot(st_geometry(bb_sol))



#convert time stamp to time in seconds

#also need to correct data for Sr 88! This should be symple! Use equation from the other paper
EDJ.dentin1.ef$X87Sr.86Sr..5.
#plot uncorrected data
plot(1:length(EDJ.dentin1.ef$X87Sr.86Sr..5.),EDJ.dentin1.ef$X87Sr.86Sr..5.)
points(1:length(EDJ.dentin1.ef$X87Sr.86Sr..5.),
       Corr.87.86.to.88(EDJ.dentin1.ef$X88Sr, EDJ.dentin1.ef$X87Sr.86Sr..5.),col = "red")

tibble(EDJ.dentin1.ef, Corr.87.86.to.88(EDJ.dentin1.ef$X88Sr, EDJ.dentin1.ef$X87Sr.86Sr..5.))

EDJ.dentin1.ef$X88Sr




