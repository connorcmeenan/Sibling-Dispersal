### building on Connor Meenan + Katrina Catalano's work with full siblings

library(tidyverse)
library(geosphere)

### read list of half-sib pairs

halfsibs<-read_csv("data/Clownfish_HalfSib_data.csv")
head(halfsibs)

### import fish metadata

metadata <-readRDS("data/fish_meta.rds")

### filter half-sub pairs to include only pair with >95% probability

halfsibs_filtered <- halfsibs %>%
  filter(Probability>=0.95)

nrow(halfsibs_filtered)
## 556 pairs

## add metadata
halfsibs_metadata <- halfsibs_filtered %>%
  left_join(metadata,by=c("OffspringID1"="ligation_id")) %>%
  left_join(metadata,by=c("OffspringID2"="ligation_id"),suffix=c("_sib1","_sib2"))
  
head(halfsibs_metadata)

## remove any pairs with missing metadata
halfsibs_noNA=na.omit(halfsibs_metadata)
nrow(halfsibs_noNA)
## 487 pairs

## Make two vectors with points for sibling 1 and sibling 2
sib1_points<-as.matrix(halfsibs_noNA[,c('lon_sib1','lat_sib1')])

sib2_points<-as.matrix(halfsibs_noNA[,c('lon_sib2','lat_sib2')])

## Calculate distances!
halfsibs_noNA$distance <- distHaversine(sib1_points,sib2_points,r=6378.137)

head(halfsibs_noNA)

#### read in full-sibling and parent-offspring pair data 

sib_dist <- readRDS("data/sib_dist.rds")

po_dist <- read.csv("data/parentage_results_allyears.csv")

## simple density plot of pairwise distance ditributions 
plot(density(na.omit(sib_dist$dist_km)))
plot(density(po_dist$dist_par_km))
plot(density(halfsibs_noNA$distance))

## grab pairwise distances for each category (full-sib,1000 unrelated pairs,parent-offspring,half-sib) 
trimmed_sibs <- sib_dist %>%
  dplyr::select(sib1_fish_indiv, sib2_fish_indiv, dist_km)
sibling_distances <- trimmed_sibs %>%
  dplyr::select(dist_km)
sibling_distances <- sibling_distances %>%
  rename(distance_km = "dist_km")

pairs_no_rep <- read.csv(file = "data/PairsWithNoRep.csv")
no_reps <- pairs_no_rep %>%
  dplyr::select(distance_km)

po_distances <-  po_dist %>%
  dplyr::select(dist_par_km)
po_distances <- po_distances %>%
  rename(distance_km = "dist_par_km")

halfsibs_distances <- halfsibs_noNA %>%
  dplyr::select(distance)
halfsibs_distances <- halfsibs_distances %>%
  rename(distance_km = "distance")

halfsibs_distances$distances <- 'half_sibling_pairs'
no_reps$distances <- 'unrelated_fish_pairs'
sibling_distances$distances <- 'sibling_fish_pairs'
po_distances$distances <- 'parent_offspring_pairs'

### plot histogram of distances - difficult to see all categories...
Distances <- rbind(no_reps,po_distances,halfsibs_distances,sibling_distances)
dist_hist_1 <- ggplot(Distances, aes(distance_km, fill = distances)) + geom_histogram(alpha = 0.5, position = 'identity') + scale_fill_manual(values=c("grey20", "grey60","grey90","black")) + theme_bw()
print(dist_hist_1)
dev.off()

### violin plot of distances
violin_1 <- ggplot(data = Distances, aes(x=distances, y=distance_km)) + geom_violin(aes(fill=distances)) + scale_fill_manual(values=c("grey50", "grey80","grey20", "grey60"), labels=c('Half-Sibling Pairs','Parent-Offspring Pairs','Sibling fish pairs', 'Unrelated fish pairs'), name="") + theme_bw() + ylim(0,30) + ylab("Distance between individuals (km)") + xlab('Fish pairs')
print(violin_1)

#### simulating a distribution of halfsib distances under random dispersal
### randomly select 2 PO distances
### add or subtract (equal probability) smaller from larger
### repeat for the number of observed half-sib pairs

halfsib_sim_dist=vector()

for (i in 1:length(halfsibs_noNA$distance)){
  dists=sample(po_dist$dist_par_km,2)
  dists_ord=dists[order(dists,decreasing=TRUE)]
  dists_ord
  dists_comp=ifelse(runif(1)<0.5,dists_ord[1]-dists_ord[2],dists_ord[1]+dists_ord[2])
  dists_comp
  halfsib_sim_dist[i]=dists_comp
}

### compare simulated distribution to observed distribution

par(mfrow=c(1,2))
plot(density(halfsibs_noNA$distance))
plot(density(halfsib_sim_dist))

ks.test(halfsibs_noNA$distance,halfsib_sim_dist)

### some of the simulated values are larger than the extent of the study area / maximum observed distance
## for now just remove these, but there might be a better solution

halfsib_sim_dist_trunc=halfsib_sim_dist[which(halfsib_sim_dist<max(halfsibs_noNA$distance))]

par(mfrow=c(1,2))
plot(density(halfsibs_noNA$distance))
plot(density(halfsib_sim_dist_trunc))

ks.test(halfsibs_noNA$distance,halfsib_sim_dist_trunc)

### generate a random distribution of full-sibs using the same methodology

fullsib_sim_dist=vector()

for (i in 1:length(na.omit(sib_dist$dist_km))){
  dists=sample(po_dist$dist_par_km,2)
  dists_ord=dists[order(dists,decreasing=TRUE)]
  dists_ord
  dists_comp=ifelse(runif(1)<0.5,dists_ord[1]-dists_ord[2],dists_ord[1]+dists_ord[2])
  dists_comp
  fullsib_sim_dist[i]=dists_comp
}

par(mfrow=c(1,2))
plot(density(na.omit(sib_dist$dist_km)))
plot(density(fullsib_sim_dist))

ks.test(na.omit(sib_dist$dist_km),fullsib_sim_dist)


#### just generate 500 random sibling distances and remove any that are > maximum PO distance

n500_sim_dist=vector()

for (i in 1:500){
  dists=sample(po_dist$dist_par_km,2)
  dists_ord=dists[order(dists,decreasing=TRUE)]
  dists_ord
  dists_comp=ifelse(runif(1)<0.5,dists_ord[1]-dists_ord[2],dists_ord[1]+dists_ord[2])
  dists_comp
  n500_sim_dist[i]=dists_comp
}

n500_sim_dist_trunc=n500_sim_dist[which(n500_sim_dist<max(po_dist$dist_par_km))]


### plot empirical cumulative distribution functions

par(mfrow=c(1,1))

plot(ecdf(pairs_no_rep$distance_km),main = "",xlab="Distance",cex=0,lwd=2,col="black")
plot(ecdf(po_dist$dist_par_km),cex=0,lwd=2,col="purple",add=T)
legend(x=15,y=0.4,legend=c("Unrelated Pairs","Parent-Offspring Pairs"),col=c("black","purple"),lty=c(1,1),lwd=2)

plot(ecdf(pairs_no_rep$distance_km),main = "",xlab="Distance",cex=0,lwd=2,col="black")
plot(ecdf(po_dist$dist_par_km),cex=0,lwd=2,col="purple",add=T)
plot(ecdf(sib_dist$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_noNA$distance),cex=0,lwd=2,col="blue",add=T)

plot(ecdf(n500_sim_dist_trunc),cex=0,lwd=2,col="gold",add=T)

legend(x=15,y=0.4,legend=c("Unrelated Pairs","Parent-Offspring Pairs","Full Sibling Pairs","Half Sibling Pairs","Simulated Sibling Pairs"),col=c("black","purple","red","blue","gold"),lty=c(1,1,1,1,1),lwd=2)


ks.test(na.omit(sib_dist$dist_km),n500_sim_dist_trunc)
ks.test(halfsibs_noNA$distance,n500_sim_dist_trunc)
ks.test(po_dist$dist_par_km,n500_sim_dist_trunc)

ks.test(halfsibs_noNA$distance,pairs_no_rep$distance_km)


##### looking at difference in capture dates and differences in size

as.Date(halfsibs_noNA$date_sib1)-as.Date(halfsibs_noNA$date_sib2)
min(na.omit(as.vector(as.Date(halfsibs_noNA$date_sib1)-as.Date(halfsibs_noNA$date_sib2))))
max(na.omit(as.vector(as.Date(halfsibs_noNA$date_sib1)-as.Date(halfsibs_noNA$date_sib2))))
plot(as.Date(halfsibs_noNA$date_sib1)-as.Date(halfsibs_noNA$date_sib2),halfsibs_noNA$size_sib1 - halfsibs_noNA$size_sib2, ylab="Size Difference",xlab="Capture Date Difference")

as.Date(sib_dist$sib1_date)-as.Date(sib_dist$sib2_date)
min(na.omit(as.Date(sib_dist$sib1_date)-as.Date(sib_dist$sib2_date)))
max(na.omit(as.Date(sib_dist$sib1_date)-as.Date(sib_dist$sib2_date)))
plot(as.Date(sib_dist$sib1_date)-as.Date(sib_dist$sib2_date),sib_dist$sib1_size - sib_dist$sib2_size, ylab="Size Difference",xlab="Capture Date Difference")

### when were pairs first captured?

plot(sib_dist$sib1_year,sib_dist$sib2_year)

plot(halfsibs_noNA$year_sib1,halfsibs_noNA$year_sib2)

halfsibs_yearmatch <- halfsibs_noNA %>% filter(year_sib1 == year_sib2)

as.Date(halfsibs_yearmatch$date_sib1)-as.Date(halfsibs_yearmatch$date_sib2)
min(na.omit(as.vector(as.Date(halfsibs_yearmatch$date_sib1)-as.Date(halfsibs_yearmatch$date_sib2))))
max(na.omit(as.vector(as.Date(halfsibs_yearmatch$date_sib1)-as.Date(halfsibs_yearmatch$date_sib2))))
plot(as.Date(halfsibs_yearmatch$date_sib1)-as.Date(halfsibs_yearmatch$date_sib2),halfsibs_yearmatch$size_sib1 - halfsibs_yearmatch$size_sib2, ylab="Size Difference",xlab="Capture Date Difference")

### try to filter by year and plot ECDF

sibs_2013 <- sib_dist %>% filter(sib1_year==2013)
halfsibs_2013 <- halfsibs_yearmatch %>% filter(year_sib1==2013)
po_2013 <- po_dist %>% filter(offs_year==2013)

plot(ecdf(po_2013$dist_par_km),cex=0,lwd=2,col="purple",main = "",xlab="Distance")
plot(ecdf(sibs_2013$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_2013$distance),cex=0,lwd=2,col="blue",add=T)

## plot other years

sibs_2012 <- sib_dist %>% filter(sib1_year==2012)
halfsibs_2012 <- halfsibs_yearmatch %>% filter(year_sib1==2012)
po_2012 <- po_dist %>% filter(offs_year==2012)

plot(ecdf(po_2012$dist_par_km),cex=0,lwd=2,xlim=c(0,30), col="purple",main = "",xlab="Distance")
plot(ecdf(sibs_2012$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_2012$distance),cex=0,lwd=2,col="blue",add=T)

sibs_2014 <- sib_dist %>% filter(sib1_year==2014)
po_2014 <- po_dist %>% filter(offs_year==2014)

plot(ecdf(po_2014$dist_par_km),cex=0,lwd=2,xlim=c(0,30),col="purple",main = "",xlab="Distance")
plot(ecdf(sibs_2014$dist_km),cex=0,lwd=2,col="red",add=T)

sibs_2015 <- sib_dist %>% filter(sib1_year==2015)
halfsibs_2015 <- halfsibs_yearmatch %>% filter(year_sib1==2015)
po_2015 <- po_dist %>% filter(offs_year==2015)

plot(ecdf(po_2015$dist_par_km),cex=0,lwd=2,xlim=c(0,30),col="purple",main = "",xlab="Distance")
plot(ecdf(sibs_2015$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_2015$distance),cex=0,lwd=2,col="blue",add=T)

halfsibs_2016 <- halfsibs_yearmatch %>% filter(year_sib1==2016)
po_2016 <- po_dist %>% filter(offs_year==2016)

plot(ecdf(po_2016$dist_par_km),cex=0,lwd=2,xlim=c(0,30),col="purple",main = "",xlab="Distance")
plot(ecdf(halfsibs_2016$distance),cex=0,lwd=2,col="blue",add=T)

sibs_2017 <- sib_dist %>% filter(sib1_year==2017)
halfsibs_2017 <- halfsibs_yearmatch %>% filter(year_sib1==2017)
po_2017 <- po_dist %>% filter(offs_year==2017)

plot(ecdf(po_2017$dist_par_km),cex=0,lwd=2,xlim=c(0,30),col="purple",main = "",xlab="Distance")
plot(ecdf(sibs_2017$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_2017$distance),cex=0,lwd=2,col="blue",add=T)

halfsibs_2018 <- halfsibs_yearmatch %>% filter(year_sib1==2018)
po_2018 <- po_dist %>% filter(offs_year==2018)

plot(ecdf(po_2018$dist_par_km),cex=0,lwd=2,xlim=c(0,30),col="purple",main = "",xlab="Distance")
plot(ecdf(halfsibs_2018$distance),cex=0,lwd=2,col="blue",add=T)

## plot size differences by year

## think about what to do with half-sibs from different years

## simulations using PO distances from same year as siblings

## did we have any parents of the full/half siblings in the PO dataset (for identifying dispersal origin)

names(sib_dist)

names(po_dist)

po_dist$par_sample_id
sib_dist$sib1_sample_id

intersect(po_dist$offs_sample_id,sib_dist$sib1_sample_id)
intersect(po_dist$offs_sample_id,sib_dist$sib2_sample_id)

sib_po<-intersect(po_dist$offs_sample_id,c(sib_dist$sib1_sample_id,sib_dist$sib2_sample_id))

po_dist$offs_sample_id %in% sib_po

po_dist$par_sample_id[po_dist$offs_sample_id %in% sib_po]

metadata$lat[which(metadata$sample_id=="APCL13_650")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL13_650")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL13_650")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL13_650")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.85,10.95),xlim=c(124.7,124.75))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## plot of APCL12_117

metadata$lat[which(metadata$sample_id=="APCL12_117")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL12_117")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL12_177")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL12_117")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.85,10.95),xlim=c(124.70,124.75))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## the following only has one identifiable match (APCL12_093)

metadata$lat[which(metadata$sample_id=="APCL12_093")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL12_093")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL12_093")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL12_093")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.80,10.93),xlim=c(124.70,124.74))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## plot of APCL13_560

metadata$lat[which(metadata$sample_id=="APCL13_560")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL13_560")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL13_560")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL13_560")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.59,10.91),xlim=c(124.71,124.78))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]


## plot of APCL13_295

metadata$lat[which(metadata$sample_id=="APCL13_295")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL13_295")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL13_295")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL13_295")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.77,10.93),xlim=c(124.71,124.79))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## plot of APCL12_154

metadata$lat[which(metadata$sample_id=="APCL12_154")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL12_154")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL12_154")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL12_154")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.79,10.91),xlim=c(124.715,124.735))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## plot of APCL15_356001

metadata$lat[which(metadata$sample_id=="APCL15_356001")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL15_356001")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL15_356001")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL15_356001")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.850,10.865),xlim=c(124.72,124.727))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]

## plot of APCL15_363819

metadata$lat[which(metadata$sample_id=="APCL15_363819")]

parent1.y= metadata$lat[which(metadata$sample_id=="APCL15_363819")]
parent1.x= metadata$lon[which(metadata$sample_id=="APCL15_363819")]

offsp1=po_dist$offs_sample_id[which(po_dist$par_sample_id=="APCL15_363819")]

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

plot(x=NULL,y=NULL,ylim=c(10.855,10.88),xlim=c(124.708,124.725))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")

po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[1])]
po_dist$dist_par_km[which(po_dist$offs_sample_id==offsp1[2])]
