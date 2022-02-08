### building on Connor Meenan + Katrina Catalano's work with full siblings

library(tidyverse)
library(geosphere)

### read list of half-sib pairs

halfsibs<-read.csv("~/Documents/GitHub/Sibling-Dispersal/data/Clownfish_HalfSib_data.csv")
head(halfsibs)

### import fish metadata

metadata <-readRDS("~/Documents/GitHub/Sibling-Dispersal/data/fish_meta.rds")

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

sib_dist <- readRDS("~/Documents/GitHub/Sibling-Dispersal/data/sib_dist.rds")

po_dist <- read.csv("~/Documents/GitHub/Sibling-Dispersal/data/parentage_results_allyears.csv")

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
