#script for masking the tracks using the land layer on the cluster
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R
#Elham Nourani. Feb. 6. 2020. Radolfzell, Germany


args <- (commandArgs(trailingOnly = TRUE)) # provides access to a copy of the command line supplied when this R session was invoked
eval(parse(text = args)) #eval evalueates an R expression in a specified environment
n <- as.numeric(as.character(Line)) # the object “Line” comes from the .slrm file. This is the index

#install packages that I couldnt install using conda install
#none in the current script

#open libraries
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
#library(parallel)
#library(lutz)
#library(move)
#library(mapview)
#library(rWind)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

#load files. make sure they are all stored in the working directory
load("ocean_0_60.RData") #ocean
load("twz_sf.RData") #twz_sf
load("tmpz_sf.RData") #tmpz_sf
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset
#load("lines.RData") #lines. prepped in data_prep_track_based
load("land_0_60_1km_buffer.RData") #land_1km
load("land_0_60.RData") #land_0_60


# #####STEP 1: remove tracks with no points over water #####
 #convert dataset to sf
 coordinates(dataset) <- ~ location.long + location.lat
 proj4string(dataset) <- wgs
 dataset_sf <- st_as_sf(dataset)
 
 dataset_sea <- dataset_sf %>% 
  st_intersects(ocean)

#keep only tracks with some points over the sea
tracks_sea <- dataset[dataset$track %in% dataset_sea$track,]

coordinates(tracks_sea)<-~location.long+location.lat
proj4string(tracks_sea)<-wgs
tracks_sea_sf <- st_as_sf(tracks_sea) #convert to sf object

# ##### STEP 2: convert tracks to spatial lines #####

#convert tracks to lines
#only keep tracks that have at least three points
less_than_three_point <- tracks_sea %>% #find out which tracks have only one or two points. the two point tracks are only in East Asia
  group_by(track) %>%
  summarise(n = length(track)) %>%
  filter(n < 3)

lines <- tracks_sea %>%
  filter(!(track %in% less_than_three_point$track)) %>%
  group_by(track) %>%
  arrange(date_time) %>%
  summarize(species = head(species,1),do_union = F) %>%
  st_cast("LINESTRING")


##### STEP 5: filter for sea-crossing segments #####

#only 1 km buffer
sea_lines_no_buffer <- lines %>% 
  st_difference(land_1km)

save(sea_lines_no_buffer, file = "R_files/sea_lines_1km.RData")

sea_seg_no_buffer <- sea_lines_no_buffer %>% 
  st_cast("MULTILINESTRING") %>% 
  st_cast("LINESTRING")

save(sea_seg_no_buffer, file = "R_files/sea_segs_1km.RData")

#number of points per segment
less_than_two<- sea_seg_no_buffer %>% #find out which tracks have only one or two points. the two point tracks are only in East Asia
  group_by(track) %>% 
  summarise(n = length(track)) %>% 
  filter(n < 3) 



###test to see why the squiggliness
#tr <- lines[lines$track == "337_2015",]
tr <- dataset[dataset$track %in% c("337_2015_autumn","337_2015_spring"),]
tr <- dataset[dataset$track %in% "337_2015_spring",]
tr <- dataset[dataset$track %in% "337_2015_autumn",]

coordinates(tr) <-~ location.long+location.lat
proj4string(tr) <- wgs

tr_sf <- tr %>%  
  st_as_sf(tr) %>%  
  arrange(date_time) %>%
  summarize(species = head(species,1),do_union = F) %>%
  st_cast("LINESTRING")

tr_sp <- as(tr_sf, "Spatial")


# b <- Sys.time()
# tr_land <- tr_sf %>% 
#   st_difference(land_1km)
# Sys.time() - b
# 
# b <- Sys.time()
# tr_land2 <- tr_sf %>% 
#   st_intersection(ocean)
# Sys.time() - b
# 
# b <- Sys.time()
# tr_land2 <- tr_sf %>% 
#   st_intersects(ocean)
# Sys.time() - b

#tr_sp <- as(tr, "Spatial")
land_sp <- as(land_1km,"Spatial")

b <- Sys.time()
r <- erase(tr_sp,ocean_sp)
Sys.time() - b

ocean_sp <- as(ocean,"Spatial")

b <- Sys.time()
r <- raster::intersect(as(tr_sp,"SpatialLines"),as(ocean_sp,"SpatialPolygons"))
Sys.time() - b


b <- Sys.time()
r <- gIntersection(as(tr_sp,"SpatialLines"),as(ocean_sp,"SpatialPolygons"))
Sys.time() - b



###try creating the lines using the sp package

tr_l <- as(tr,"SpatialLines")

#go over the list of tracks, for each track, filter out land


