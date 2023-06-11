library(scales)
library(MASS)
library(viridisLite)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(changepoint)

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
Enamel6ext2 <- read.csv("data/Enamel6ext2.csv")
Enamel7 <- read.csv("data/Enamel7.csv")
Enamel7ext <- read.csv("data/Enamel7ext.csv")
Enamel8 <- read.csv("data/Enamel8.csv")
Enamel8ext <- read.csv("data/Enamel8ext.csv")
Enamel9 <- read.csv("data/Enamel9.csv")
Enamel9ext2 <- read.csv("data/Enamel9ext2.csv")
Enamel10 <- read.csv("data/Enamel10.csv")

pl.EDJ.dentin1 <- read.csv("data/pl_EDJ_dentin1.csv")
pl.EDJ.dentin2 <- read.csv("data/pl_EDJ_dentin2.csv")
pl.Enamel1 <- read.csv("data/pl_Enamel1.csv")
pl.Enamel2 <- read.csv("data/pl_Enamel2.csv")
pl.Enamel3 <- read.csv("data/pl_Enamel3.csv")
pl.Enamel4 <- read.csv("data/pl_Enamel4.csv")
pl.Enamel5 <- read.csv("data/pl_Enamel5.csv")
pl.Enamel6 <- read.csv("data/pl_Enamel6.csv")
pl.Enamel6ext <- read.csv("data/pl_Enamel6ext.csv")
pl.Enamel6ext2 <- read.csv("data/pl_Enamel6ext2.csv")
pl.Enamel7 <- read.csv("data/pl_Enamel7.csv")
pl.Enamel7ext <- read.csv("data/pl_Enamel7ext.csv")
pl.Enamel8 <- read.csv("data/pl_Enamel8.csv")
pl.Enamel8ext <- read.csv("data/pl_Enamel8ext.csv")
pl.Enamel9 <- read.csv("data/pl_Enamel9.csv")
pl.Enamel9ext2 <- read.csv("data/pl_Enamel9ext2.csv")
pl.Enamel10 <- read.csv("data/pl_Enamel10.csv")

############data processing using custom function###########

#dentine at EDJ
proc.EDJ.dentin1 <- data.process(EDJ.dentin1, pl.EDJ.dentin1, "d")
proc.EDJ.dentin1
proc.EDJ.dentin1.rm <- moving.avg(proc.EDJ.dentin1, 25)

proc.EDJ.dentin2 <- data.process(EDJ.dentin2, pl.EDJ.dentin2, "d")
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


#enamel transect 6 ext
proc.Enamel6ext <- data.process(Enamel6ext, pl.Enamel6ext)
proc.Enamel6ext
proc.Enamel6ext.rm <- moving.avg(proc.Enamel6ext, 25)

proc.Enamel7 <- data.process(Enamel7, pl.Enamel7)
proc.Enamel7
proc.Enamel7.rm <- moving.avg(proc.Enamel7, 25)

###################################### transform corrdinates ###############
x.transf <- 100
y.transf <- -890

proc.Enamel6ext2 <- data.process(Enamel6ext2, pl.Enamel6ext2)
proc.Enamel6ext2.x.t <- proc.Enamel6ext2$x + x.transf
proc.Enamel6ext2.y.t <- proc.Enamel6ext2$y + y.transf

proc.Enamel6ext2.t <- tibble(corr.87Sr.86Sr = proc.Enamel6ext2$corr.87Sr.86Sr, 
                             x = proc.Enamel6ext2.x.t, y = proc.Enamel6ext2.y.t)

proc.Enamel6ext2.rm <- moving.avg(proc.Enamel6ext2.t, 25)

proc.Enamel7ext <- data.process(Enamel7ext, pl.Enamel7ext)
proc.Enamel7ext.x.t <- proc.Enamel7ext$x + x.transf
proc.Enamel7ext.y.t <- proc.Enamel7ext$y + y.transf

proc.Enamel7ext.t <- tibble(corr.87Sr.86Sr = proc.Enamel7ext$corr.87Sr.86Sr, 
                             x = proc.Enamel7ext.x.t, y = proc.Enamel7ext.y.t)

proc.Enamel7ext.rm <- moving.avg(proc.Enamel7ext.t, 25)

proc.Enamel8 <- data.process(Enamel8, pl.Enamel8)
proc.Enamel8.x.t <- proc.Enamel8$x + x.transf
proc.Enamel8.y.t <- proc.Enamel8$y + y.transf

proc.Enamel8.t <- tibble(corr.87Sr.86Sr = proc.Enamel8$corr.87Sr.86Sr, 
                            x = proc.Enamel8.x.t, y = proc.Enamel8.y.t)

proc.Enamel8.rm <- moving.avg(proc.Enamel8.t, 25)

proc.Enamel8ext <- data.process(Enamel8ext, pl.Enamel8ext)
proc.Enamel8ext.x.t <- proc.Enamel8ext$x + x.transf
proc.Enamel8ext.y.t <- proc.Enamel8ext$y + y.transf

proc.Enamel8ext.t <- tibble(corr.87Sr.86Sr = proc.Enamel8ext$corr.87Sr.86Sr, 
                         x = proc.Enamel8ext.x.t, y = proc.Enamel8ext.y.t)

proc.Enamel8ext.rm <- moving.avg(proc.Enamel8ext.t, 25)


proc.Enamel9 <- data.process(Enamel9, pl.Enamel9)
proc.Enamel9.x.t <- proc.Enamel9$x + x.transf
proc.Enamel9.y.t <- proc.Enamel9$y + y.transf

proc.Enamel9.t <- tibble(corr.87Sr.86Sr = proc.Enamel9$corr.87Sr.86Sr, 
                         x = proc.Enamel9.x.t, y = proc.Enamel9.y.t)

proc.Enamel9.rm <- moving.avg(proc.Enamel9.t, 25)


proc.Enamel9ext2 <- data.process(Enamel9ext2, pl.Enamel9ext2)
proc.Enamel9ext2.x.t <- proc.Enamel9ext2$x + x.transf
proc.Enamel9ext2.y.t <- proc.Enamel9ext2$y + y.transf

proc.Enamel9ext2.t <- tibble(corr.87Sr.86Sr = proc.Enamel9ext2$corr.87Sr.86Sr, 
                             x = proc.Enamel9ext2.x.t, y = proc.Enamel9ext2.y.t)

proc.Enamel9ext2.rm <- moving.avg(proc.Enamel9ext2.t, 25)


proc.Enamel10 <- data.process(Enamel10, pl.Enamel10)
proc.Enamel10.x.t <- proc.Enamel10$x + x.transf
proc.Enamel10.y.t <- proc.Enamel10$y + y.transf

proc.Enamel10.t <- tibble(corr.87Sr.86Sr = proc.Enamel10$corr.87Sr.86Sr, 
                         x = proc.Enamel10.x.t, y = proc.Enamel10.y.t)

proc.Enamel10.rm <- moving.avg(proc.Enamel10.t, 25)