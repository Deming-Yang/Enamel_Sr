library(MASS)
library(dplyr)
library(tidyr)
library(sf)

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
  
  for(i in 1:nrow(df)){
    
    x.li.int[i] <- approx(x = ref$tt, y = ref$x, xout = df$tt[i])$y
    y.li.int[i] <- approx(x = ref$tt, y = ref$y, xout = df$tt[i])$y
    
  }
  
  return(data.frame(x = x.li.int, y = y.li.int))
  
}

##fn 3: correct 87/86 ratio by Sr88 voltage
#Horstwood et al. (2008)
Corr1.87.86.to.88 <- function(Sr88, Sr87.86){
  accuracy <- 0.002333/Sr88 + 0.999462
  corr.Sr87.86 <- Sr87.86/accuracy
  return(corr.Sr87.86)
}

#Lugli et al. (2017), rhino enamel
Corr2.87.86.to.88 <- function(Sr88, Sr87.86){
  accuracy <- 0.0005954/Sr88 + 1.0000077
  corr.Sr87.86 <- Sr87.86/accuracy
  return(corr.Sr87.86)
}

################data processing function###########
data.process <- function(df, log, type = "e"){
  require(tidyr)
  require(dplyr)
  
  #convert character to number
  df$X88Sr<- as.numeric(df$X88Sr)
  df$X87Sr.86Sr..5.<- as.numeric(df$X87Sr.86Sr..5.)

  df.ef <- filter(.data = df, df$X88Sr > 0.2)
  
  df.ef <- na.omit(df.ef)
  
  #correction by Sr 88 voltage
  if(type == "d"){#dentine
    X87Sr.86Sr.corr <- Corr1.87.86.to.88(df.ef$X88Sr, df.ef$X87Sr.86Sr..5.)
  }
  else{
    X87Sr.86Sr.corr <- Corr2.87.86.to.88(df.ef$X88Sr, df.ef$X87Sr.86Sr..5.)
  }
  
  #get time stamp from the Time coloumn
  df.ef <-  separate_wider_delim(df.ef, Time, delim = ":", names = c("hour", "minute", "second", "milisecond"),
                                 too_few = "align_start")
  
  log.ef <- separate_wider_delim(log,Time, delim = ":", names = c("hour", "minute", "second", "milisecond"),
                                 too_few = "align_start")
  
  #calculate total time in both data frames
  tt.df.ef <- time.total(df.ef)
  
  tt.log.ef <- time.total(log.ef)
  
  t.gap <- mean(c(tt.df.ef[1]-tt.log.ef[1], tt.df.ef[length(tt.df.ef)]-tt.log.ef[length(tt.log.ef)]))
  
  #built in 7 second offset
  df.ef.tt <- tibble(df.ef, tt = tt.df.ef - t.gap)
  
  log.ef.tt <- tibble(log.ef, tt = tt.log.ef)
  
  #linear interpolate
  
  df.ef.tt.li <- li.interp(df = df.ef.tt, ref = log.ef.tt)
  
  #create a new data frame to store the corrected and interpolated data:
  corr.df <- tibble(corr.87Sr.86Sr = X87Sr.86Sr.corr, x =df.ef.tt.li$x, y = df.ef.tt.li$y)
  
  return(na.omit(corr.df))
}

###########fn 4 moving average and sd, n points########
moving.avg <- function(df, n){
  require(zoo)
  rm.dat <- rollmean(df$corr.87Sr.86Sr, n)
  rsd.dat <- rollapply(df$corr.87Sr.86Sr, n, sd)
  n.remove <- floor(n/2)
  n.x <- df$x[(n.remove + 1): (nrow(df)-n.remove)]
  n.y <- df$y[(n.remove + 1): (nrow(df)-n.remove)]
  tibble(avg = rm.dat, sd = rsd.dat, x = n.x, y = n.y)
}


###########fn 5 projecting a point to a curve and transform its x and y ########
proj.pt <- function(x, y, ref){
  
  ref.s <- ref[order(ref$y),]
  
  ref.s.diff.x <- diff(ref.s$x)
  ref.s.diff.y <- diff(ref.s$y)
  
  #calculate distance between adjacent points
  ref.s.dist <- sqrt(ref.s.diff.x^2 + ref.s.diff.y^2)
  
  new.x.axis <- c(0, cumsum(ref.s.dist)) #"faltten" the EDJ with new "x" axis values
  
  dist.to.xy <- sqrt((x - ref.s$x)^2 + (y - ref.s$y)^2)
  
  new.y <- min(dist.to.xy) #smallest element between xy and the curve
  
  indx <- which(dist.to.xy == new.y) #extract index of the point
  
  new.x <- new.x.axis[indx]
  
  return(c(new.x, new.y))
  
}

###########fn 6 projecting a transect to a curve and transform its x and y ########

proj.transect <- function(trans, ref){
  
  n <- nrow(trans)
  
  diff.x <- diff(ref$x)
  diff.y <- diff(ref$y)
  
  dist <- sqrt(diff.x^2 + diff.y^2)
  
  ref.x <- c(0, cumsum(dist))
  
  new.x <- rep(0, n)
  
  dist.to.xy <- sqrt((trans$x[1] - ref$x)^2 + (trans$y[1] - ref$y)^2)
  
  indx <- which(dist.to.xy == min(dist.to.xy)) #extract index of the point
  
  #aligning the first data point to new x axis
  
  trans.diff.x <- diff(trans$x)
  trans.diff.y <- diff(trans$y)
  
  trans.dist <- sqrt(trans.diff.x^2 + trans.diff.y^2)
  
  #initiate vectors
  new.x <- c(ref.x[indx], ref.x[indx] + cumsum(trans.dist))
  
  new.y <- rep(0, n)
  
  #calculate distance between adjacent points
  for(i in 1:n){
    
    dist.to.ref <- sqrt((trans$x[i] - ref$x)^2 + (trans$y[i] - ref$y)^2)
    
    new.y[i] <- min(dist.to.ref) #smallest element between xy and the curve
    
  }
  
  return(data.frame(new.x, new.y))
 
}

###########fn 7 calculate 50 point average to reduce data amount for LA-ICP-MS data ########
x.pt.avg <- function(df.Sr, df.tl, n.avg){
  # three objects
  # one is the Sr ratio
  
  n <- floor(length(df.Sr)/n.avg) #discard some data towards the end of the df
  
  #initiate vectors

  n.avg.sr <- rep(0, n)
  n.sd.sr <- rep(0, n)
  
  n.avg.tl <- rep(0, n)
  
  for(i in 1:n){
    x <- ((i-1)*n.avg + 1):(i*n.avg)
    
    n.avg.sr[i] <- mean(df.Sr[x])
    
    n.sd.sr[i] <- sd(df.Sr[x])
    
    n.avg.tl[i] <- mean(df.tl[x])
  }
  return(data.frame(avg.sr = n.avg.sr, sd.sr = n.sd.sr, avg.tl = n.avg.tl))
  
}

########### fn 8 MCMC credible interval calculations ###########
MCMC.CI.bound <- function (MCMC.res, CI){
  require(KernSmooth)
  require(bayestestR)
  dim.MCMC <- dim(MCMC.res)
  #typically, the first element is # of iterations
  #the second element is the time series
  map.res <- rep(0, dim.MCMC[2])
  hdi.high <- rep(0, dim.MCMC[2])
  hdi.low <- rep(0, dim.MCMC[2])
  
  for(i in 1:dim.MCMC[2]){
    map.res[i] <- map_estimate(MCMC.res[,i], method = "KernSmooth")
    hdi <- hdi(MCMC.res[,i], ci = CI)
    hdi.low[i] <- hdi$CI_low
    hdi.high[i] <- hdi$CI_high
  }
  
  return (list(map.res, hdi.low, hdi.high, CI))
}

########### fn 9 simulate enamel with maturation overprint###########
# the model is after Passey et al 2002.
# with 1-D values and two parameters
# f.ma = fraction of maturation overprint
# ma.wind = maturation window (days)
# the function returns a 1-D "history" of sim enamel
# equivalent to data from the innermost enamel

forw.m.en <- function(serum, f.ma, ma.wind){
  if (f.ma >1 | f.ma < 0){
    warning("f must be between 0 and 1")
    break
  } else{
    
    serum.rm <- rollmean(serum, ma.wind)
    
    serum.rm.sft <- c(serum.rm, rep(NA, ma.wind - 1))
    
    sim.en <- serum * (1 - f.ma) + serum.rm.sft * f.ma
    
    return(sim.en)
  }
}

