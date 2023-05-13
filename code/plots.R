library(scales)
library(MASS)
library(viridisLite)
library(EnvStats)
library(dplyr)
library(tidyr)
library(sf)
library(zoo)
library(changepoint)
library(ggplot2)

par(mfrow=c(3,3))

plot(nrow(proc.EDJ.dentin2.rm):1,proc.EDJ.dentin2.rm$avg, ylim=c(0.702,0.712),type="l")
plot(nrow(proc.EDJ.dentin1.rm):1,proc.EDJ.dentin1.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel1.rm),proc.Enamel1.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel2.rm),proc.Enamel2.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel3.rm),proc.Enamel3.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel4.rm),proc.Enamel4.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel5.rm),proc.Enamel5.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel6ext.rm),proc.Enamel6ext.rm$avg, ylim=c(0.702,0.712),type="l")

plot(1:nrow(proc.Enamel7.rm),proc.Enamel7.rm$avg, ylim=c(0.702,0.712),type="l")


#####use the sf package to plot those values

all.dat <- rbind(proc.EDJ.dentin1, proc.EDJ.dentin2, proc.Enamel1, proc.Enamel2, proc.Enamel3,
                 proc.Enamel4, proc.Enamel5, proc.Enamel6, proc.Enamel6ext, proc.Enamel7)

all.dat.rm <- rbind(dent.rm.f, proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm, proc.Enamel4.rm, 
                    proc.Enamel5.rm, proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel7.rm.f)

all.enamel.rm <- rbind(proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
                       proc.Enamel4.rm, proc.Enamel5.rm, proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel7.rm.f)
all.dentine.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)


all.enamel4.rm <- rbind(proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
                        proc.Enamel4.rm)
require(sf)
sf.all.dat <- st_as_sf(all.dat,  agr = NA_agr_,
                       coords = c("y","x"),
                       dim = "XYZ",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)

sf.all.dat.rm <- st_as_sf(all.dat.rm,  agr = NA_agr_,
                          coords = c("y","x"),
                          dim = "XYZ",
                          remove = TRUE,
                          na.fail = TRUE,
                          sf_column_name = NULL)

sf.all.enamel.rm <- st_as_sf(all.enamel.rm,  agr = NA_agr_,
                             coords = c("y","x"),
                             dim = "XYZ",
                             remove = TRUE,
                             na.fail = TRUE,
                             sf_column_name = NULL)

sf.all.dentine.rm <- st_as_sf(all.dentine.rm,  agr = NA_agr_,
                              coords = c("y","x"),
                              dim = "XYZ",
                              remove = TRUE,
                              na.fail = TRUE,
                              sf_column_name = NULL)

sf.all.enamel4.rm <- st_as_sf(all.enamel4.rm,  agr = NA_agr_,
                              coords = c("y","x"),
                              dim = "XYZ",
                              remove = TRUE,
                              na.fail = TRUE,
                              sf_column_name = NULL)


proc.Enamel6.rm

plot(sf.all.dat, breaks = seq(0.704, 0.715,by = 0.001), pch=16,main ="all")

plot(sf.all.dat.rm, breaks = seq(0.704, 0.712,by = 0.001), pch=16,main ="avg")

plot(sf.all.enamel.rm, breaks = seq(0.704, 0.712,by = 0.001), pch=16,main ="all enamel",cex=0.4)

plot(sf.all.enamel4.rm, breaks = seq(0.704, 0.712,by = 0.001), pch=16,main ="enamel 1-4",cex=0.5)


plot(sf.all.enamel.rm, breaks = c(0.700,seq(0.704, 0.712,by = 0.001)), pch=16,main ="all enamel",cex=0.4)
plot(sf.cp.loc,add=T,pch=16,cex=2)

plot(sf.cp.loc, col ="black",pch=17)

#plot change positions within enamel:
cp.loc <- data.frame(avg = rep(0.7,length(cp.E.rm.x)), x = cp.E.rm.x, y = cp.E.rm.y)
sf.cp.loc <- st_as_sf(cp.loc,  agr = NA_agr_,
                      coords = c("y","x"),
                      dim = "XYZ",
                      remove = TRUE,
                      na.fail = TRUE,
                      sf_column_name = NULL)

cp.D.loc <- data.frame(avg = rep(0.7,length(cp.D.rm.x)), x = cp.D.rm.x, y = cp.D.rm)
sf.cp.D.loc <- st_as_sf(cp.D.loc,  agr = NA_agr_,
                      coords = c("y","x"),
                      dim = "XYZ",
                      remove = TRUE,
                      na.fail = TRUE,
                      sf_column_name = NULL)
#############plot 1 all enamel data with transition points#######
ggplot()+
  geom_sf(data = sf.all.enamel.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp.loc, color = "black",pch=18,cex=3)+ 
  theme(legend.position = "bottom")

#############plot 2 all dentine and enamel data with transition points#######
ggplot()+
  geom_sf(data = sf.all.dat.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp.loc, color = "black",pch=18,cex=3)+ 
  geom_sf(data = sf.cp.D.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")