library(scales)
library(viridisLite)
library(dplyr)
library(sf)
library(ggplot2)

###########real world tooth geometry###########
#############plot 1 all dentine and enamel data with transition points#######
par(mfrow=c(1,1))
ggplot()+
  geom_sf(data = sf.all.dat.rm , aes(color = mov.avg))+
  scale_color_viridis_c()+
  theme(legend.position = "bottom")

###########EDJ flattened geometry###########
##############both transitions#################
n.y <- 4

all.dat.narm.proj <- mutate(all.dat.narm.proj, Y2 = n.y * new.y)

cp1.loc <- mutate(cp1.loc, Y2 = n.y * y)

cp2.loc <- mutate(cp2.loc, Y2 = n.y * y)

sf.all.dat.y2 <- st_as_sf(all.dat.narm.proj,  agr = NA_agr_,
                               coords = c("new.x","Y2"),
                               dim = "XY",
                               remove = TRUE,
                               na.fail = TRUE,
                               sf_column_name = NULL)

sf.cp1.y2 <- st_as_sf(cp1.loc,  agr = NA_agr_,
                       coords = c("x","Y2"),
                       dim = "XY",
                       remove = TRUE,
                       na.fail = TRUE,
                       sf_column_name = NULL)

sf.cp2.y2 <- st_as_sf(cp2.loc,  agr = NA_agr_,
                      coords = c("x","Y2"),
                      dim = "XY",
                      remove = TRUE,
                      na.fail = TRUE,
                      sf_column_name = NULL)

ggplot()+
  geom_sf(data = sf.all.dat.y2 , aes(color = mov.avg),cex=3.5)+
  scale_color_viridis_c()+

  geom_abline(intercept = n.y * 3.02e+03, slope = n.y * -tan(3.3/180*pi), color="gray64",
              linewidth=2.5)+
  geom_abline(intercept = lm.cp1.E.proj.inner$coefficients[1] * n.y,
              slope = lm.cp1.E.proj.inner$coefficients[2] * n.y,
              color="black",
              linetype="dashed", linewidth=1)+
  geom_abline(intercept = n.y * 3.317e+03, slope = n.y * -tan(3.3/180*pi), color="gray64",
              linewidth=2.5)+
  geom_abline(intercept = lm.cp2.E.proj.inner$coefficients[1] * n.y,
              slope = lm.cp2.E.proj.inner$coefficients[2] * n.y,
              color="black",
              linetype="dashed", linewidth=1)+
  
  geom_sf(data = sf.cp2.y2, color = "black",pch=18,cex=4)+ 
  geom_sf(data = sf.cp1.y2, color = "red",pch=18,cex=4)+ 
  theme(legend.position = "bottom")

############ end of change point comparisons ###########

############# Fig 3 ############# 
# Fig 3 after converting to time, 
# compare laser, hand drill and micromill transects using the molar plate geometry

par(mfrow=c(1,5)) #1200 * 400

# Panel 1 molar enamel1 vs micromill tusk dentine 
plot(-1000, -1, col= "gray56",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Reference turnover",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr, lty = 2)
abline(h = UT.Sr, lty = 2)

polygon(c(en1.tl$tl, rev(en1.tl$tl)), 
        c(en1.tl$Sr + en1.tl$sd, 
          rev(en1.tl$Sr - en1.tl$sd)), 
        col = alpha("orange", 0.3), border = NA)
lines(days.cumm.en1.al, proc.Enamel1.rm.f$mov.avg,col = alpha("orange", 0.9), lwd=2)

# reference fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]],lwd = 2)
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[2]], lty = 5)
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[3]], lty = 5)

graphics::text(400, 0.707, "f = 0.00", cex = 1.3)

points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
       pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))

# Panel 2 molar enamel9 vs micromill tusk dentine 
plot(-1000, -1, col= "gray56",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Overprint: En9",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr, lty = 2)
abline(h = UT.Sr, lty = 2)

polygon(c(en9.tl$tl, rev(en9.tl$tl)), 
        c(en9.tl$Sr + en9.tl$sd, 
          rev(en9.tl$Sr - en9.tl$sd)), 
        col = alpha("orange3", 0.3), border = NA)
lines(days.cumm.en9.al, proc.Enamel9.rm.f$mov.avg,col = alpha("orange3", 0.9), lwd=2)

# reference fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]],lwd = 2, col = "gray56")
# sample fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.En9.89[[1]],lwd = 2, lty = 2)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
graphics::text(400, 0.707, "f = 0.06", cex = 1.3)

# Panel 3 molar enamel10 vs micromill tusk dentine 
plot(-1000, -1, col= "gray56",
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Overprint: En10",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr, lty = 2)
abline(h = UT.Sr, lty = 2)

polygon(c(en10.tl$tl, rev(en10.tl$tl)), 
        c(en10.tl$Sr + en10.tl$sd, 
          rev(en10.tl$Sr - en10.tl$sd)), 
        col = alpha("orange4", 0.3), border = NA)
lines(days.cumm.en10.al, proc.Enamel10.rm.f$mov.avg,col = alpha("orange4", 0.9), lwd=2)

# reference fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]],lwd = 2, col = "gray56")
# sample fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.En10.89[[1]],lwd = 2, lty = 2)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
graphics::text(400, 0.707, "f = 0.13", cex = 1.3)

# Panel 4 molar enamel vs micromill tusk dentine 
plot(Edrill.tl$tl, Edrill.tl$Sr, col= "red4",
     pch=16, cex = 2,
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Overprint: drill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr, lty = 2)
abline(h = UT.Sr, lty = 2)

# reference fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]],lwd = 2, col = "gray56")
# sample fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.drill.89[[1]],lwd = 2, lty = 2)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
graphics::text(400, 0.707, "f = 0.23", cex = 1.3)

# Panel 5 molar enamel vs micromill molar enamel 
plot(Rm3.5b.mill.tl$tl, Rm3.5b.mill.tl$Sr, col= "cyan4",
     pch=16, cex = 2,
     xlim=c(-400,700),ylim=c(0.705,0.7115),
     xlab="Days from Misha's move",
     main = "Overprint: micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr, lty = 2)
abline(h = UT.Sr, lty = 2)

# reference fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1m.89[[1]],lwd = 2, col = "gray56")
# sample fit
lines(bin.thin.oc*1:t.oc - 400, post.comb.R1.Rm3.5b.89[[1]],lwd = 2, lty = 2)
# points(tusk.mill.tl$tl, tusk.mill.tl$Sr,
#        pch=18, cex = 2.2, col = alpha("#00b4ffff", 0.8))
graphics::text(400, 0.707, "f = 0.08", cex = 1.3)

#################  Fig 4 #################
############# estimated intake ###########

par(mfrow=c(1, 2)) #1100 * 450
# reconstructed Sr input signal from LA-ICP-MS molar enamel
plot(En1.50avg.tl$avg.tl, En1.50avg.tl$avg.sr, type= "l", col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on molar enamel LA-ICP-MS",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[1]], lwd = 2)
lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[2]], lty = 2)
lines(1:t.En1 * En1.bt - d.offset.En1 - offset.en1, post.misha.En1.Rin.m.89[[3]], lty = 2)
legend(400,0.709,c("Est. intake","Measured"),lwd = c(2,2), col=c("black","#00b4ffff"))


# reconstructed Sr input signal from micromilled tusk dentine
plot(tusk.mill.tl$tl, tusk.mill.tl$Sr, pch=18, col= "#00b4ffff",
     lwd = 2,
     xlim=c(-400,1000),ylim=c(0.705,0.713),
     xlab="Days from Misha's move",
     main = "Estimated intake based on tusk dentine micromill",
     ylab = "87Sr/86Sr") 
abline(h = CA.Sr)
abline(h = UT.Sr)

lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[1]], lwd = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[2]], lty = 2)
lines(1:t.M640b * M640.bt -d.offset.M640-365, post.misha.M640b.Rin.m.89[[3]], lty = 2)
