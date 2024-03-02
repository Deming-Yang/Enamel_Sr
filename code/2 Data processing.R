library(scales)
library(MASS)
library(viridisLite)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(changepoint)

source("./code/1 Helper functions.R")

EDJ.dentin1 <- read.csv("./data/EDJ_dentin1.csv")
EDJ.dentin2 <- read.csv("./data/EDJ_dentin2.csv")
Enamel1 <- read.csv("./data/Enamel1.csv")
Enamel2 <- read.csv("./data/Enamel2.csv")
Enamel3 <- read.csv("./data/Enamel3.csv")
Enamel4 <- read.csv("./data/Enamel4.csv")
Enamel5 <- read.csv("./data/Enamel5.csv")
Enamel6 <- read.csv("./data/Enamel6.csv")
Enamel6ext <- read.csv("./data/Enamel6ext.csv")
Enamel6ext2 <- read.csv("./data/Enamel6ext2.csv")
Enamel7 <- read.csv("./data/Enamel7.csv")
Enamel7ext <- read.csv("./data/Enamel7ext.csv")
Enamel8 <- read.csv("./data/Enamel8.csv")
Enamel8ext <- read.csv("./data/Enamel8ext.csv")
Enamel9 <- read.csv("./data/Enamel9.csv")
Enamel9ext2 <- read.csv("./data/Enamel9ext2.csv")
Enamel10 <- read.csv("./data/Enamel10.csv")

pl.EDJ.dentin1 <- read.csv("./data/pl_EDJ_dentin1.csv")
pl.EDJ.dentin2 <- read.csv("./data/pl_EDJ_dentin2.csv")
pl.Enamel1 <- read.csv("./data/pl_Enamel1.csv")
pl.Enamel2 <- read.csv("./data/pl_Enamel2.csv")
pl.Enamel3 <- read.csv("./data/pl_Enamel3.csv")
pl.Enamel4 <- read.csv("./data/pl_Enamel4.csv")
pl.Enamel5 <- read.csv("./data/pl_Enamel5.csv")
pl.Enamel6 <- read.csv("./data/pl_Enamel6.csv")
pl.Enamel6ext <- read.csv("./data/pl_Enamel6ext.csv")
pl.Enamel6ext2 <- read.csv("./data/pl_Enamel6ext2.csv")
pl.Enamel7 <- read.csv("./data/pl_Enamel7.csv")
pl.Enamel7ext <- read.csv("./data/pl_Enamel7ext.csv")
pl.Enamel8 <- read.csv("./data/pl_Enamel8.csv")
pl.Enamel8ext <- read.csv("./data/pl_Enamel8ext.csv")
pl.Enamel9 <- read.csv("./data/pl_Enamel9.csv")
pl.Enamel9ext2 <- read.csv("./data/pl_Enamel9ext2.csv")
pl.Enamel10 <- read.csv("./data/pl_Enamel10.csv")

#load Rm3.5 hand drill data
Drill <- read.csv("./data/Rm3.5 hand drill.csv")

Drill.no <- na.omit(Drill)

#load Rm3.5 micromill data
Rm3.5b.mill <- read.csv("./data/Rm3.5b micromill.csv")
Rm3.5b.mill.no <- na.omit(Rm3.5b.mill)
# sort from largest dist to smallest
desc(Rm3.5b.mill.no)

# load Misha tusk dentine micromill data (Yang et al. 2023)
misha.tusk.micromill <- read.csv("data/Misha dentin micromill.csv")

misha.tusk.micromill<-na.omit(misha.tusk.micromill)

Rm3.5.angle <- read.csv("data/Rm3.5 appo angle.csv")

#load M640 tusk dentine micromill C and O isotope data (Uno 2012)
M640.micromill <- read.csv("data/M640 tusk C and O.csv")

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

#####################transform all coordinates to evaluate Sr slopes###########
dent.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)
#first transform dentine transect
dent.rm.s <- dent.rm[order(dent.rm$y),]

dent.rm.s.diff.x <- diff(dent.rm.s$x)
dent.rm.s.diff.y <- diff(dent.rm.s$y)

#calculate distance between adjacent points
dent.rm.s.dist <- sqrt(dent.rm.s.diff.x^2 + dent.rm.s.diff.y^2)

dent.rm.new.x <- c(0, cumsum(dent.rm.s.dist))

dent.rm.new.y <- rep(-150, length(dent.rm.new.x)) #dentine is -150 microns away from EDJ

dent.rm.proj <- tibble(dent.rm, data.frame(new.x = dent.rm.new.x, new.y = dent.rm.new.y))

# dent.rm.proj <- tibble(dent.rm, data.frame(new.x = dent.rm.new.x, new.y = dent.rm.new.y))
#######Enamel 1######
proc.Enamel1.rm.s <- proc.Enamel1.rm[order(proc.Enamel1.rm$y),]

E1.proj <- proj.transect(proc.Enamel1.rm.s, dent.rm.s)

proc.Enamel1.rm.proj <- tibble(proc.Enamel1.rm.s, E1.proj)

#######Enamel 2######
proc.Enamel2.rm.s <- proc.Enamel2.rm[order(proc.Enamel2.rm$y),]

E2.proj <- proj.transect(proc.Enamel2.rm.s, dent.rm.s)

proc.Enamel2.rm.proj <- tibble(proc.Enamel2.rm.s, E2.proj)

#######Enamel 3######
proc.Enamel3.rm.s <- proc.Enamel3.rm[order(proc.Enamel3.rm$y),]

E3.proj <- proj.transect(proc.Enamel3.rm.s, dent.rm.s)

proc.Enamel3.rm.proj <- tibble(proc.Enamel3.rm.s, E3.proj)

#######Enamel 4######
proc.Enamel4.rm.s <- proc.Enamel4.rm[order(proc.Enamel4.rm$y),]

E4.proj <- proj.transect(proc.Enamel4.rm.s, dent.rm.s)

proc.Enamel4.rm.proj <- tibble(proc.Enamel4.rm.s, E4.proj)

#######Enamel 5######
proc.Enamel5.rm.s <- proc.Enamel5.rm[order(proc.Enamel5.rm$y),]

E5.proj <- proj.transect(proc.Enamel5.rm.s, dent.rm.s)

proc.Enamel5.rm.proj <- tibble(proc.Enamel5.rm.s, E5.proj)

#######Enamel 6######
proc.Enamel6.rm <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel6ext2.rm)

proc.Enamel6.rm.s <- proc.Enamel6.rm[order(proc.Enamel6.rm$y),]

E6.proj <- proj.transect(proc.Enamel6.rm.s, dent.rm.s)

proc.Enamel6.rm.proj <- tibble(proc.Enamel6.rm.s, E6.proj)

#######Enamel 7######
proc.Enamel7.rm <- rbind(proc.Enamel7.rm, proc.Enamel7ext.rm)

proc.Enamel7.rm.s <- proc.Enamel7.rm[order(proc.Enamel7.rm$y),]

E7.proj <- proj.transect(proc.Enamel7.rm.s, dent.rm.s)

proc.Enamel7.rm.proj <- tibble(proc.Enamel7.rm.s, E7.proj)

#######Enamel 8######
proc.Enamel8.rm <- rbind(proc.Enamel8.rm, proc.Enamel8ext.rm)

proc.Enamel8.rm.s <- proc.Enamel8.rm[order(proc.Enamel8.rm$y),]

E8.proj <- proj.transect(proc.Enamel8.rm.s, dent.rm.s)

proc.Enamel8.rm.proj <- tibble(proc.Enamel8.rm.s, E8.proj)

#######Enamel 9######
proc.Enamel9.rm <- rbind(proc.Enamel9.rm, proc.Enamel9ext2.rm)

proc.Enamel9.rm.s <- proc.Enamel9.rm[order(proc.Enamel9.rm$y),]

E9.proj <- proj.transect(proc.Enamel9.rm.s, dent.rm.s)

proc.Enamel9.rm.proj <- tibble(proc.Enamel9.rm.s, E9.proj)

#######Enamel 10######
proc.Enamel10.rm.s <- proc.Enamel10.rm[order(proc.Enamel10.rm$y),]

E10.proj <- proj.transect(proc.Enamel10.rm.s, dent.rm.s)

proc.Enamel10.rm.proj <- tibble(proc.Enamel10.rm.s, E10.proj)