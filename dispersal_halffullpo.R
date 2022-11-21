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
library(wesanderson)
ecdf_palette=wes_palette("Darjeeling1")

par(mfrow=c(1,1))

plot(ecdf(pairs_no_rep$distance_km),main = "",xlab="Distance",cex=0,lwd=2,col="black")
plot(ecdf(po_dist$dist_par_km),cex=0,lwd=2,col="purple",add=T)
legend(x=15,y=0.4,legend=c("Unrelated Pairs","Parent-Offspring Pairs"),col=c("black","purple"),lty=c(1,1),lwd=2)

plot(ecdf(pairs_no_rep$distance_km),main = "",xlab="Distance",cex=0,lwd=2,col=ecdf_palette[1])
plot(ecdf(po_dist$dist_par_km),cex=0,lwd=2,col=ecdf_palette[2],add=T)
plot(ecdf(sib_dist$dist_km),cex=0,lwd=2,col=ecdf_palette[3],add=T)
plot(ecdf(halfsibs_noNA$distance),cex=0,lwd=2,col=ecdf_palette[4],add=T)

plot(ecdf(n500_sim_dist_trunc),cex=0,lwd=2,col=ecdf_palette[5],add=T)

legend(x=15,y=0.4,legend=c("Unrelated Pairs","Parent-Offspring Pairs","Full Sibling Pairs","Half Sibling Pairs","Simulated Sibling Pairs"),col=c("black","purple","red","blue","gold"),lty=c(1,1,1,1,1),lwd=2)


ks.test(na.omit(sib_dist$dist_km),n500_sim_dist_trunc)
ks.test(halfsibs_noNA$distance,n500_sim_dist_trunc)
ks.test(po_dist$dist_par_km,n500_sim_dist_trunc)

ks.test(halfsibs_noNA$distance,pairs_no_rep$distance_km)

### calculating mean, median, sd, CV for Parent Offspring
po_dist
po_dist$dist_par_km
##mean
mean(po_dist$dist_par_km)
##median
median(po_dist$dist_par_km)
##standard deviation
sd(po_dist$dist_par_km)
#CV
100*sd(po_dist$dist_par_km)/mean(po_dist$dist_par_km)

### calculating mean, median, sd, CV for Parent Offspring
po_dist
po_dist$dist_par_km
##mean
mean(po_dist$dist_par_km)
##median
median(po_dist$dist_par_km)
##standard deviation
sd(po_dist$dist_par_km)
#CV
100*sd(po_dist$dist_par_km)/mean(po_dist$dist_par_km)

### calculating mean, median, sd, CV for Half Siblings
halfsibs_noNA
halfsibs_noNA$distance
##mean
mean(halfsibs_noNA$distance)
##median
median(halfsibs_noNA$distance)
##standard deviation
sd(halfsibs_noNA$distance)
#CV
100*sd(halfsibs_noNA$distance)/mean(halfsibs_noNA$distance)

### calculating mean, median, sd, CV for Full Siblings
na.omit(sib_dist)
na.omit(sib_dist$dist_km)
##mean
mean(na.omit(sib_dist$dist_km))
##median
median(na.omit(sib_dist$dist_km))
##standard deviation
sd(na.omit(sib_dist$dist_km))
#CV
100*sd(na.omit(sib_dist$dist_km))/mean(na.omit(sib_dist$dist_km))

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

## Plotting PO data

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

## "APCL13_295" missing from metadata, but there is a third sibling?

offsp1_1.y= metadata$lat[which(metadata$sample_id==offsp1[1])]
offsp1_1.x= metadata$lon[which(metadata$sample_id==offsp1[1])]

offsp1_2.y= metadata$lat[which(metadata$sample_id==offsp1[2])]
offsp1_2.x= metadata$lon[which(metadata$sample_id==offsp1[2])]

offsp1_3.y= metadata$lat[which(metadata$sample_id==offsp1[3])]
offsp1_3.x= metadata$lon[which(metadata$sample_id==offsp1[3])]

plot(x=NULL,y=NULL,ylim=c(10.78,10.93),xlim=c(124.70,124.8))
points(parent1.x,parent1.y,col="red")
points(offsp1_1.x,offsp1_1.y,col="blue")
points(offsp1_2.x,offsp1_2.y,col="blue")
points(offsp1_3.x,offsp1_3.y,col="blue")

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

### I think this might be included in both the PO and sibling lists?

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

## offfspring  are the same as APCL12_154 (both parents sampled!)

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

## putting all of the PO coordinates in one script

po_dist$par_sample_id
sib_dist$sib1_sample_id

#sib1s in parent-offspring list
intersect(po_dist$offs_sample_id,sib_dist$sib1_sample_id)
#sib2s in parent-offspring list
intersect(po_dist$offs_sample_id,sib_dist$sib2_sample_id)

#parents of all of these siblings
po_dist$par_sample_id[po_dist$offs_sample_id %in% sib_po]

#either sib1 or sib2 in parent-offspring list
sibs_to_map<-intersect(po_dist$offs_sample_id,c(sib_dist$sib1_sample_id,sib_dist$sib2_sample_id))

sib_parent_df=data.frame(matrix(ncol=12,nrow=length(sibs_to_map)))
colnames(sib_parent_df)=c("parent1.id","parent1.site","parent1.y","parent1.x","parent2.id","parent2.site","parent2.y","parent2.x","offsp.id","offsp.site","offsp.y","offsp.x")

#which offspring have two parents?
which(po_dist$offs_sample_id==sibs_to_map[2])
which(po_dist$offs_sample_id==sibs_to_map[6])
which(po_dist$offs_sample_id==sibs_to_map[9])
which(po_dist$offs_sample_id==sibs_to_map[10])


for(i in 1:length(sibs_to_map)){
  sib_parent_df$parent1.id[i]=po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[1]]
  sib_parent_df$parent1.site[i]=metadata$site[which(metadata$sample_id==sib_parent_df$parent1.id[i])]
  sib_parent_df$parent1.y[i]=metadata$lat[which(metadata$sample_id==sib_parent_df$parent1.id[i])]
  sib_parent_df$parent1.x[i]=metadata$lon[which(metadata$sample_id==sib_parent_df$parent1.id[i])]
  sib_parent_df$parent2.id[i]=ifelse(is.na(po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[2]]),"NA",po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[2]])
  sib_parent_df$parent2.site[i]=ifelse(is.na(po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[2]]),"NA",metadata$site[which(metadata$sample_id==sib_parent_df$parent1.id[i])])
  sib_parent_df$parent2.y[i]=ifelse(is.na(po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[2]]),"NA",metadata$lat[which(metadata$sample_id==sib_parent_df$parent1.id[i])])
  sib_parent_df$parent2.x[i]=ifelse(is.na(po_dist$par_sample_id[which(po_dist$offs_sample_id==sibs_to_map[i])[2]]),"NA",metadata$lon[which(metadata$sample_id==sib_parent_df$parent1.id[i])])
  sib_parent_df$offsp.id[i]=sibs_to_map[i]
  sib_parent_df$offsp.site[i]=metadata$site[which(metadata$sample_id==sib_parent_df$offsp.id[i])]
  sib_parent_df$offsp.y[i]=metadata$lat[which(metadata$sample_id==sib_parent_df$offsp.id[i])]
  sib_parent_df$offsp.x[i]=metadata$lon[which(metadata$sample_id==sib_parent_df$offsp.id[i])]
}

## making maps!

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("ggplot2")
theme_set(theme_bw())
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")

# plotting world map
ggplot(data = world) +
  geom_sf()

#plotting just philippines
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(110, 130), ylim = c(5, 20), expand = FALSE)
  
#plotting just the study site plus one parent-offspring pair 
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.5, 125), ylim = c(10.7, 11), expand = FALSE) +
  geom_point(aes(x=parent1.x,y=parent1.y),colour="blue") +
  geom_point(aes(x=offsp1_1.x,y=offsp1_1.y),colour="red") 

#plotting just the study site plus all parent-offspring pairs
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.6, 124.9), ylim = c(10.6, 10.9), expand = FALSE) +
  geom_curve(data=sib_parent_df,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc")))
  
#plotting just the study site plus all parent-offspring pairs
ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.7, 124.74), ylim = c(10.84, 10.88), expand = FALSE) +
  geom_curve(data=sib_parent_df,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc"))) 

# select offspring of APCL13_560 and plot
APCL13_560_offspring=sib_parent_df[which(sib_parent_df$parent1.id=="APCL13_560"),]

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.7, 124.95), ylim = c(10.63, 10.88), expand = FALSE) +
  geom_segment(data=APCL13_560_offspring,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc"))) 

# select offspring of APCL12_093 and plot
APCL12_093_offspring=sib_parent_df[which(sib_parent_df$parent1.id=="APCL12_093"),]

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.7, 124.8), ylim = c(10.78, 10.88), expand = FALSE) +
  geom_segment(data=APCL12_093_offspring,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc"))) 

# select offspring from short-distance dispersals and plot
sib_parent_short=sib_parent_df[which(sib_parent_df$parent1.id!="APCL12_093" & sib_parent_df$parent1.id!="APCL13_560"),]

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.705, 124.73), ylim = c(10.85, 10.88), expand = FALSE) +
  geom_segment(data=sib_parent_short,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc"))) 

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.705, 124.73), ylim = c(10.85, 10.88), expand = FALSE) +
  geom_curve(data=sib_parent_short,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y),arrow = arrow(length = unit(0.01, "npc"))) 

#coloured by parent ID

ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(124.705, 124.73), ylim = c(10.853, 10.876), expand = FALSE) +
  geom_curve(data=sib_parent_short,aes(x=parent1.x,y=parent1.y,xend=offsp.x,yend=offsp.y,colour=parent1.id),arrow = arrow(length = unit(0.05, "npc"))) +
  geom_point(data=sib_parent_short,aes(x=parent1.x,y=parent1.y,colour=parent1.id),size=0.1)



## Github Push test

## plot size differences by year

## think about what to do with half-sibs from different years

## simulations using PO distances from same year as siblings

## did we have any parents of the full/half siblings in the PO dataset (for identifying dispersal origin)

######## find pairs that have the same or very similar lat + long + dates (= sampled on the same anemone?)

sib_dist

##### find # of decimal points
decimalplaces <- function(x) {
  ifelse(abs(x - round(x)) > .Machine$double.eps^0.5,
         nchar(sub('^\\d+\\.', '', sub('0+$', '', as.character(x)))),
         0)
}

decimalplaces(c(234.1, 3.7500, 1.345, 3e-15))

decimalplaces(sib_dist$sib1_lat)
decimalplaces(sib_dist$sib1_lon)
decimalplaces(sib_dist$sib2_lat)
decimalplaces(sib_dist$sib2_lon)

### 5-12 places!!

### finding sibs with same time of capture
which(sib_dist$sib1_time_date==sib_dist$sib2_time_date)
### finding sibs with same lat
which(sib_dist$sib1_lat==sib_dist$sib2_lat)
### finding sibs with same lat
which(sib_dist$sib1_lon==sib_dist$sib2_lon)

### finding half-sibs with same time of capture
which(halfsibs_noNA$time_date_sib1==halfsibs_noNA$time_date_sib2)
### finding half-sibs with same lat
which(halfsibs_noNA$lat_sib1==halfsibs_noNA$lat_sib2)
### finding half-sibs with same lat
which(halfsibs_noNA$lon_sib1==halfsibs_noNA$lon_sib2)

### finding PO with same time of capture
which(po_dist$par_time_date==po_dist$offs_time_date)
### finding PO with same lat
which(po_dist$par_lat==po_dist$offs_lat)
### finding PO with same lat
which(po_dist$par_lon==po_dist$offs_lon)




#### 2014 sibling group

sib_dist %>% 
  filter(sib1_year==2014) %>% 
  filter(sib1_site=="Palanas") -> sib_2014Palanas

sib_2014Palanas$dist_km

unique(c(sib_2014Palanas$sib1_sample_id,sib_2014Palanas$sib2_sample_id)) -> sib_2014Palanas_IDs

sib_dist %>% 
  filter(!sib1_sample_id %in% sib_2014Palanas_IDs) %>% 
  filter(!sib2_sample_id %in% sib_2014Palanas_IDs) -> sib_no2014Palanas

sib_no2014Palanas$dist_km


##### plotting 


### plot empirical cumulative distribution functions

par(mfrow=c(1,1))

plot(ecdf(pairs_no_rep$distance_km),main = "",xlab="Distance",cex=0,lwd=2,col="black")
plot(ecdf(po_dist$dist_par_km),cex=0,lwd=2,col="purple",add=T)
plot(ecdf(sib_2014Palanas$dist_km),cex=0,lwd=2,col="red",add=T)
plot(ecdf(halfsibs_noNA$distance),cex=0,lwd=2,col="blue",add=T)
plot(ecdf(sib_no2014Palanas$dist_km),cex=0,lwd=2,col="gold",add=T)

legend(x=15,y=0.5,legend=c("Unrelated Pairs","Parent-Offspring Pairs","Full Sibling Pairs (2014 Palanas)","Half Sibling Pairs","Full Sibling Pairs (Others)"),col=c("black","purple","red","blue","gold"),lty=c(1,1,1,1,1),lwd=2)

### load in anemone data

load("data/anem_db.RData")
load("data/fish_db.RData")


sib_dist_noNA=na.omit(sib_dist)

sib_dist_noNA$sib1_anemID=NA

for (i in 1:nrow(sib_dist_noNA)){ 
  sib_dist_noNA$sib1_anemID[i]=anem_db$anem_id[which(anem_db$anem_table_id==fish_db$anem_table_id[which(fish_db$sample_id==sib_dist_noNA$sib1_sample_id[i])])]
}

sib_dist_noNA$sib2_anemID=NA

for (i in 1:nrow(sib_dist_noNA)){ 
  sib_dist_noNA$sib2_anemID[i]=anem_db$anem_id[which(anem_db$anem_table_id==fish_db$anem_table_id[which(fish_db$sample_id==sib_dist_noNA$sib2_sample_id[i])])]
}

#find siblings captured on the same anemone!
which(sib_dist_noNA$sib1_anemID==sib_dist_noNA$sib2_anemID)