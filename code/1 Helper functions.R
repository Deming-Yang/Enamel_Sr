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

########### fn 9 initiate switches ###########
initiate.switch <- function(t, n.switch, day.switch, a, gap, duration){
  
  #initial vector
  v.init <- rep(a, (day.switch -1))
  v.switch <- rep(a+gap, duration)
  v.back <- rep(a, duration)
  sw.d <- n.switch*duration
  
  v.switches <- c(v.switch, v.back)
  
  rep.switches <- rep(v.switches,ceiling(n.switch/2))
  
  if((day.switch -1 + sw.d) < t){#switch is completed before the total number of days
    t.end <- t-(day.switch -1 + sw.d)
    
    if((n.switch %% 2) == 0){#even number of switches
      r.end <- a
    }
    else{#odd number of switches
      r.end <- a + gap
    }
    v.end <- rep(r.end, t.end)
    all.input <- c(v.init, rep.switches[1:sw.d], v.end)
  }
  else{#switch is not completed before the total number of days
    #no need to worry about end values
    all.input <- c(v.init, rep.switches[1:sw.d])#take advantage of the vector recycling feature
  }
  Sr.input <- all.input[1:t]
  return(Sr.input)
}

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
