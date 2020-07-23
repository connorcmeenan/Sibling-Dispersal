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

#redid distances
trimmed_fish <- readRDS("trimmed_fish_meta_copy.rds")
fish_1 <- trimmed_fish %>%
rename(fish_indiv_1 = fish_indiv) %>%
rename(fish_indiv1_lat = lat) %>%
rename(fish_indiv1_lon = lon)
head(fish_1)

fish_2 <- trimmed_fish %>%
rename(fish_indiv_2 = fish_indiv) %>%
rename(fish_indiv2_lat = lat) %>%
rename(fish_indiv2_lon = lon)
head(fish_2)

library(spectralGP)
alldists_redo <- rdist.earth(as.matrix(fish_1[,c('fish_indiv1_lon', 'fish_indiv1_lat')]), as.matrix(fish_2[,c('fish_indiv2_lon', 'fish_indiv2_lat')]), miles=FALSE, R=6371)

#installed reshape2
install.packages('reshape2')
library('reshape2')
melted_alldists_redo <- melt(alldists_redo)
fish_pairs_df1 <- melted_alldists_redo %>%
rename(fish_1 = Var1) %>%
rename(fish_2 = Var2) %>%
rename(dist_km = value)
head(fish_pairs_df1)

#turn matrix into df
alldists <- readRDS('alldists_copy.rds')
fish_pairs_df <- as.data.frame(t(alldists))
head(fish_pairs_df)

# name columns and add column with row names
trimmed <- readRDS('trimmed_fish_meta_copy.rds')
id_vector <- pull(trimmed, fish_indiv)
named_columns <- fish_pairs_df
names(named_columns) <- id_vector
head(named_columns)
named_rows_columns <- named_columns
named_rows_columns <- add_column(named_columns, id_vector, .before = "985153000371280")

# make all possible combos of indiv fish
install.packages('gtools')
library(gtools)
permutated_fish_no_reps <- permutations(n=2903, r=2, v=id_vector)
no_reps_fish_pairs <- as.data.frame(t(permutated_fish_no_reps))

anti_join(fish_dists, sib_dist, by=c(id_vector="sib1_fish_indiv", fish_2id="sib2_fish_indiv "))
#make correct df
gathered_2 <- gather(fish_dists, 2:2904, key="fish_2id", value="id")
named_fish_dists <- gathered_2 %>%
rename(fish_1id = id_vector) %>%
rename(distance_km = id)
head(named_fish_dists)

#remove fish compared to self
filtered_fish_dist <- named_fish_dist %>%
filter(fish_1id != fish_2id)
head(filtered_fish_dist)

#remove "duplicated" fish


#remove sibling pairs
sib_dist <- readRDS('sib_dist.rds')
df_wo_sibs <- anti_join(filtered_fish_dist, sib_dist, by=c(fish_1id="sib1_fish_indiv", fish_2id="sib2_fish_indiv"))

#remove duplicates using the matrix with deleted upper triangle
upper_triag <- readRDS(upper_tri_alldists.rds)
upper_triag <- readRDS('upper_tri_alldists.rds')
upper_triag_df <- as.data.frame(t(upper_triag))
trimmed <- readRDS('trimmed_fish_meta_copy.rds')
id_vector <- pull(trimmed, fish_indiv)
named_upper_triag <- upper_triag_df
names(named_upper_triag) <- id_vector
head(named_upper_triag)
named_rows_columns_upper_triag <- named_upper_triag
upper_triag_named <- add_column(named_rows_columns_upper_triag, id_vector, .before = "985153000371280")
gathered_upper <- gather(upper_triag_named, 2:2904, key="fish2_id", value="distance_km")
fish_no_duplicates <- subset(filtered_upper, distance_km!=0)
final_df_missing_13k_rows <- anti_join(fish_no_duplicates, sib_dist, by=c(fish1_id="sib1_fish_indiv", fish2_id="sib2_fish_indiv"))
