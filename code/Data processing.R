library(coda)
library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(lubridate)
library(signal)

############helper functions##############
###fn 1: convert time stamp to seconds####
time.total <- function(df){
  if(is.na(df$milisecond[1])==T){
    df.h.s <- as.numeric(df$hour) * 3600
    df.m.s <- as.numeric(df$minute) * 60
    df.s.s <- as.numeric(df$second)
    df.s.total <- df.h.s + df.m.s + df.s.s
  }else{
    df.h.s <- as.numeric(df$hour) * 3600
    df.m.s <- as.numeric(df$minute) * 60
    df.s.s <- as.numeric(df$second)
    df.s.ms <- as.numeric(df$milisecond) * 1e-3
    
    df.s.total <- df.h.s + df.m.s + df.s.s + df.s.ms
  }
  return(df.s.total)
}

##fn 2: linear interpolate x and y from time stamp after conversion####
li.interp <- function (df, ref){
  #initiate vectors
  x.li.int <- rep(0, nrow(df))
  y.li.int <- rep(0, nrow(df))
  
  require(signal)
  for(i in 1:nrow(df)){
    
    x.li.int[i] <- interp1(ref$tt, ref$x, xi=df$tt[i], method = "linear", extrap=T)
    y.li.int[i] <- interp1(ref$tt, ref$y, xi=df$tt[i], method = "linear", extrap=T)
    
  }
  
  return(tibble(df, x.li.int, y.li.int))
  
}

##fn 3: correct 87/86 ratio by Sr88 voltage

Corr.87.86.to.88 <- function(Sr88, Sr87.86){
  accuracy <- 0.002333/Sr88 + 0.999462
  corr.Sr87.86 <- Sr87.86/accuracy
  return(corr.Sr87.86)
}


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


EDJ.dentin1.ef <- EDJ.dentin1 %>% filter(X88Sr > 0.2)
EDJ.dentin1.ef <- select(EDJ.dentin1.ef, -Column1) #remove last column

EDJ.dentin1.ef <- na.omit(EDJ.dentin1.ef)

EDJ.dentin1.ef$X87Sr.86Sr..5.<- as.numeric(EDJ.dentin1.ef$X87Sr.86Sr..5.)
EDJ.dentin1.ef$X87Sr.86Sr..8.<- as.numeric(EDJ.dentin1.ef$X87Sr.86Sr..8.)
EDJ.dentin1.ef$X88Sr.86Sr..9.<- as.numeric(EDJ.dentin1.ef$X88Sr.86Sr..9.)

#get time stamp from the Time coloumn
EDJ.dentin1.ef <- EDJ.dentin1.ef %>%
  separate_wider_delim(Time, delim = ":", names = c("hour", "minute", "second", "milisecond"),
                       too_few = "align_start")


#convert time stamp to time in seconds

#also need to correct data for Sr 88! This should be symple! Use equation from the other paper
EDJ.dentin1.ef$X87Sr.86Sr..5.
#plot uncorrected data
plot(1:length(EDJ.dentin1.ef$X87Sr.86Sr..5.),EDJ.dentin1.ef$X87Sr.86Sr..5.)
points(1:length(EDJ.dentin1.ef$X87Sr.86Sr..5.),
       Corr.87.86.to.88(EDJ.dentin1.ef$X88Sr, EDJ.dentin1.ef$X87Sr.86Sr..5.),col = "red")

tibble(EDJ.dentin1.ef, Corr.87.86.to.88(EDJ.dentin1.ef$X88Sr, EDJ.dentin1.ef$X87Sr.86Sr..5.))

EDJ.dentin1.ef$X88Sr

Corr.87.86.to.88 <- function(Sr88, Sr87.86){
  accuracy <- 0.002333/Sr88 + 0.999462
  corr.Sr87.86 <- Sr87.86/accuracy
  return(corr.Sr87.86)
}


tt.EDJ.dentin1.ef <- time.total(EDJ.dentin1.ef)
#calculate time interval between time stamps
diff(time.total(EDJ.dentin1.ef)) #~0.527 sec per acquisition

#position log has fewer data points than data points
#so data points have to be interpolated to the position log

####intake of position log####
#pl stands for "position log"
pl.EDJ.dentin1 <- read.csv("data/pl_EDJ_dentin1.csv")
pl.EDJ.dentin1.ef <- pl.EDJ.dentin1 %>%
  separate_wider_delim(Time, delim = ":", names = c("hour", "minute", "second", "milisecond"),
                       too_few = "align_start")

tt.pl.EDJ.dentin1.ef <- time.total(pl.EDJ.dentin1.ef)

tt.EDJ.dentin1.ef[1]
tt.pl.EDJ.dentin1.ef[1] #about 6 seconds of difference.
tt.EDJ.dentin1.ef[1] - tt.pl.EDJ.dentin1.ef[1]


tt.EDJ.dentin1.ef[length(tt.EDJ.dentin1.ef)]
tt.pl.EDJ.dentin1.ef[length(tt.pl.EDJ.dentin1.ef)] #about 7.5 seconds of difference.

tt.EDJ.dentin1.ef[length(tt.EDJ.dentin1.ef)] - tt.pl.EDJ.dentin1.ef[length(tt.pl.EDJ.dentin1.ef)]

#Q: is pl wider?
tt.pl.EDJ.dentin1.ef[length(tt.pl.EDJ.dentin1.ef)] - tt.pl.EDJ.dentin1.ef[1]

tt.EDJ.dentin1.ef[length(tt.EDJ.dentin1.ef)] - tt.EDJ.dentin1.ef[1]
#pl is wider! So some choice of time should work!

#try a 7 second delay

EDJ.dentin1.ef.tt <- tibble(EDJ.dentin1.ef, tt = tt.EDJ.dentin1.ef -7)

EDJ.dentin1.ef.tt$tt

pl.EDJ.dentin1.ef.tt <- tibble(pl.EDJ.dentin1.ef, tt = tt.pl.EDJ.dentin1.ef)

pl.EDJ.dentin1.ef.tt$tt

EDJ.dentin1.ef.tt.li <- li.interp(df = EDJ.dentin1.ef.tt, ref = pl.EDJ.dentin1.ef.tt)

EDJ.dentin1.ef.tt.li$x.li.int
EDJ.dentin1.ef.tt.li$y.li.int

plot(EDJ.dentin1.ef.tt.li$x.li.int, EDJ.dentin1.ef.tt.li$y.li.int)

####next step: packaging this into a function for data proccessing#####



