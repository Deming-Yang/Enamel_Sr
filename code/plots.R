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

##################plot moving averages###################
# par(mfrow=c(3,3))
# 
# plot(nrow(proc.EDJ.dentin2.rm):1,proc.EDJ.dentin2.rm$avg, ylim=c(0.702,0.712),type="l")
# plot(nrow(proc.EDJ.dentin1.rm):1,proc.EDJ.dentin1.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel1.rm),proc.Enamel1.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel2.rm),proc.Enamel2.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel3.rm),proc.Enamel3.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel4.rm),proc.Enamel4.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel5.rm),proc.Enamel5.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel6ext.rm),proc.Enamel6ext.rm$avg, ylim=c(0.702,0.712),type="l")
# 
# plot(1:nrow(proc.Enamel7.rm),proc.Enamel7.rm$avg, ylim=c(0.702,0.712),type="l")

#####use the sf package to plot these values

all.dat <- rbind(proc.EDJ.dentin1, proc.EDJ.dentin2, proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                       proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                       proc.Enamel9.rm.f, proc.Enamel10.rm.f)

all.dat.rm <- rbind(dent.rm.f, proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                    proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                    proc.Enamel9.rm.f, proc.Enamel10.rm.f)

# all.enamel.rm <- rbind(proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
#                        proc.Enamel4.rm, proc.Enamel5.rm, proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel7.rm.f)

all.enamel.rm <- rbind(proc.Enamel1.rm.f, proc.Enamel2.rm.f, proc.Enamel3.rm.f, proc.Enamel4.rm.f, 
                       proc.Enamel5.rm.f, proc.Enamel6.rm.f, proc.Enamel7.rm.f, proc.Enamel8.rm.f,
                       proc.Enamel9.rm.f, proc.Enamel10.rm.f)


all.dentine.rm <- rbind(proc.EDJ.dentin1.rm, proc.EDJ.dentin2.rm)

all.enamel4.rm <- rbind(proc.Enamel1.rm, proc.Enamel2.rm, proc.Enamel3.rm,
                        proc.Enamel4.rm)

#try seeing how the two scans can be put together
all.enamel6.rm <- rbind(proc.Enamel6.rm, proc.Enamel6ext.rm, proc.Enamel6ext2.rm)

all.enamel7.rm <- rbind(proc.Enamel7.rm, proc.Enamel7ext.rm)

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

sf.all.enamel6.rm <- st_as_sf(all.enamel6.rm,  agr = NA_agr_,
                              coords = c("y","x"),
                              dim = "XYZ",
                              remove = TRUE,
                              na.fail = TRUE,
                              sf_column_name = NULL)

sf.all.enamel7.rm.f <- st_as_sf(all.enamel7.rm,  agr = NA_agr_,
                              coords = c("y","x"),
                              dim = "XYZ",
                              remove = TRUE,
                              na.fail = TRUE,
                              sf_column_name = NULL)

#plot change positions within enamel:
cp1.loc <- data.frame(avg = rep(0.7,length(cp1.E.rm.x)), x = cp1.E.rm.x, y = cp1.E.rm.y)
sf.cp1.loc <- st_as_sf(cp1.loc,  agr = NA_agr_,
                      coords = c("y","x"),
                      dim = "XYZ",
                      remove = TRUE,
                      na.fail = TRUE,
                      sf_column_name = NULL)

cp2.loc <- data.frame(avg = rep(0.7,length(cp2.E.rm.x)), x = cp2.E.rm.x, y = cp2.E.rm.y)
sf.cp2.loc <- st_as_sf(cp2.loc,  agr = NA_agr_,
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


#############plot 1 all dentine and enamel data with transition points#######
ggplot()+
  geom_sf(data = sf.all.dat.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp2.loc, color = "black",pch=18,cex=3)+ 
  geom_sf(data = sf.cp1.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")

#################one transition################
ggplot()+
  geom_sf(data = sf.all.enamel.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp1.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")

##############both transitions#################
ggplot()+
  geom_sf(data = sf.all.enamel.rm , aes(color = avg))+
  scale_color_viridis_c()+
  geom_sf(data = sf.cp2.loc, color = "black",pch=18,cex=3)+ 
  geom_sf(data = sf.cp1.loc, color = "red",pch=18,cex=3)+ 
  theme(legend.position = "bottom")

ggplot()+
  geom_sf(data = sf.all.enamel6.rm , aes(color = avg))+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")

ggplot()+
  geom_sf(data = sf.all.enamel7.rm.f , aes(color = avg))+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")


#################plot 2 Sr ratio transitions using moving averages###############
par(mfrow=c(2,5))
#preliminary plot
# plot(dent.rm.f$y, dent.rm.f$avg, main = "Dentine",col="darkgray", 
#      pch=16,xlim=c(4e4,8e4),ylim=c(0.705,0.711))
# plot(fit_segmented.D.rm, add = T)
# lines.segmented(fit_segmented.D.rm)
# points.segmented(fit_segmented.D.rm)

#preliminary plot
plot(proc.Enamel1.rm.f$y, proc.Enamel1.rm.f$avg, main = "Enamel 1",col="darkgray", 
     pch=16,xlim=c(4e4,8e4),ylim=c(0.705,0.711))
plot(fit_segmented.E1.rm, add = T)
lines.segmented(fit_segmented.E1.rm)
points.segmented(fit_segmented.E1.rm)

#preliminary plot
plot(proc.Enamel2.rm.f$y, proc.Enamel2.rm.f$avg, main = "Enamel 2",col="darkgray", 
     pch=16,xlim=c(4e4,8e4),ylim=c(0.705,0.711))
plot(fit_segmented.E2.rm, add=T)
lines.segmented(fit_segmented.E2.rm)
points.segmented(fit_segmented.E2.rm)

#preliminary plot
plot(proc.Enamel3.rm.f$y, proc.Enamel3.rm.f$avg, main = "Enamel 3",col="darkgray"
     , pch=16,xlim=c(3e4,7e4),ylim=c(0.705,0.711))
plot(fit_segmented.E3.rm, add=T)
lines.segmented(fit_segmented.E3.rm)
points.segmented(fit_segmented.E3.rm)

#preliminary plot
plot(proc.Enamel4.rm.f$y, proc.Enamel4.rm.f$avg, main = "Enamel 4",col="darkgray", 
     pch=16,xlim=c(3e4,7e4),ylim=c(0.705,0.711))
plot(fit_segmented.E4.rm, add=T)
lines.segmented(fit_segmented.E4.rm)
points.segmented(fit_segmented.E4.rm)

#preliminary plot
plot(proc.Enamel5.rm.f$y, proc.Enamel5.rm.f$avg, main = "Enamel 5",col="darkgray"
     , pch=16,xlim=c(2e4,6e4),ylim=c(0.705,0.711))
plot(fit_segmented.E5.rm, add =T)
lines.segmented(fit_segmented.E5.rm)
points.segmented(fit_segmented.E5.rm)

#preliminary plot
plot(proc.Enamel6.rm.f$y, proc.Enamel6.rm.f$avg, main = "Enamel 6",col="darkgray"
     , pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
plot(fit_segmented.E6.rm, add =T)
lines.segmented(fit_segmented.E6.rm)
points.segmented(fit_segmented.E6.rm)

#preliminary plot
plot(proc.Enamel7.rm.f$y, proc.Enamel7.rm.f$avg, main = "Enamel 7",col="darkgray", 
     pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
plot(fit_segmented.E7.rm, add =T)
lines.segmented(fit_segmented.E7.rm)
points.segmented(fit_segmented.E7.rm)

# plot(dent.rm.f$y, dent.rm.f$avg, main = "Dentine full length",col="darkgray", 
#      pch=16,ylim=c(0.705,0.711))
# plot(fit_segmented.D.rm, add = T)
# lines.segmented(fit_segmented.D.rm)
# points.segmented(fit_segmented.D.rm)

# plot(proc.Enamel1.rm$y, proc.Enamel1.rm$avg, main = "Enamel 1 full length",col="darkgray", 
#      pch=16,ylim=c(0.705,0.711))
# plot(fit_segmented.E1.rm, add = T)
# lines.segmented(fit_segmented.E1.rm)
# points.segmented(fit_segmented.E1.rm)

plot(proc.Enamel8.rm.f$y, proc.Enamel8.rm.f$avg, main = "Enamel 8",col="darkgray", 
     pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
plot(fit_segmented.E8.rm, add =T)
lines.segmented(fit_segmented.E8.rm)
points.segmented(fit_segmented.E8.rm)

plot(proc.Enamel9.rm.f$y, proc.Enamel9.rm.f$avg, main = "Enamel 9",col="darkgray", 
     pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
plot(fit_segmented.E9.rm, add =T)
lines.segmented(fit_segmented.E9.rm)
points.segmented(fit_segmented.E9.rm)

plot(proc.Enamel10.rm.f$y, proc.Enamel10.rm.f$avg, main = "Enamel 10",col="darkgray", 
     pch=16,xlim=c(1e4,5e4),ylim=c(0.705,0.711))
plot(fit_segmented.E10.rm, add =T)
lines.segmented(fit_segmented.E10.rm)
points.segmented(fit_segmented.E10.rm)