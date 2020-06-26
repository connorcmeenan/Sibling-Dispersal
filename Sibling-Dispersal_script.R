getwd()
[1] "/Users/connormeenan"
setwd("/Users/connormeenan/desktop/Sibling-Dispersal")
getwd()
[1] "/Users/connormeenan/Desktop/Sibling-Dispersal"
fish_meta <- readRDS("data/fish_meta.rds")
getwd()
[1] "/Users/connormeenan/Desktop/Sibling-Dispersal"
head(fish_meta)

#first attempt at table and disatnce calculations
trimmed_fish_meta <- fish_meta %>%
select(fish_indiv, lat, lon)
head(trimmed_fish_meta)

indiv_fish_dist <- trimmed_fish_meta %>%
mutate(fish_indiv_1=fish_indiv) %>%
mutate(fish_indiv1_lat=lat) %>%
mutate(fish_indiv1_lon=lon)
head(indiv_fish_dist)

indiv_fish_dist <- indiv_fish_dist %>%
rename(fish_indiv_2 = fish_indiv) %>%
rename(fish_indiv2_lat = lat) %>%
rename(fish_indiv2_lon = lon)
head(indiv_fish_dist)

indiv_fish_dist_reordered <- indiv_fish_dist[c(4,1,5,6,2,3)]
head(indiv_fish_dist_reordered)

indiv_fish_dist <- indiv_fish_dist_reordered
head(indiv_fish_dist)

#calculate the distance 
alldists <- rdist.earth(as.matrix(indiv_fish_dist[,c('fish_indiv1_lon', 'fish_indiv1_lat')]), as.matrix(indiv_fish_dist[,c('fish_indiv2_lon', 'fish_indiv2_lat')]), miles=FALSE, R=6371)
head(alldists)

indiv_fish_dist$fish_km_between_fish1_fish2 <- diag(alldists)
head(indiv_fish_dist)

hist(indiv_fish_dist$fish_km_between_fish1_fish2)

#backwards way of selecting data with only unique combinations of distances from above alldists matrix
install.packages('matrixcalc')
library(matrixcalc)
upper_tri_alldists <- upper.triangle(alldists)
head(upper_tri_alldists)
hist(upper_tri_alldists)