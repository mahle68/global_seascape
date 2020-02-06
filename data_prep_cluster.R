#script for masking the tracks using the land layer on the cluster
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R
#Elham Nourani. Feb. 6. 2020. Radolfzell, Germany


args <- (commandArgs(trailingOnly = TRUE)) # provides access to a copy of the command line supplied when this R session was invoked
eval(parse(text = args)) #eval evalueates an R expression in a specified environment
n <- as.numeric(as.character(Line)) # the object “Line” comes from the .slrm file. This is the index

#install packages that I couldnt install using conda install

#open libraries
library(tidyverse)
library(lubridate)
#library(move)
#library(mapview)
#library(rWind)
library(sf)
#library(parallel)
#library(lutz)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

#load files. make sure they are all stored in the working directory
load("ocean_0_60.RData") #ocean
load("twz_sf.RData") #twz_sf
load("tmpz_sf.RData") #tmpz_sf
load("R_files/all_spp_unfiltered_updated_lc_0_removed.RData") #dataset
load("R_files/lines.RData") #lines


#####STEP 1: remove tracks with no points over water #####
#convert dataset to sf
coordinates(dataset) <- ~ location.long + location.lat
proj4string(dataset) <- wgs
dataset_sf <- st_as_sf(dataset)

dataset_sea <- dataset_sf %>%  #consider a one km buffer to remove tracks following the coast.
  st_intersection(ocean)

#keep only tracks with some points over the sea
tracks_sea <- dataset[dataset$track %in% dataset_sea$track,]

coordinates(tracks_sea)<-~location.long+location.lat
proj4string(tracks_sea)<-wgs
tracks_sea_sf <- st_as_sf(tracks_sea) #convert to sf object

##### STEP 2: convert tracks to spatial lines #####

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
tr <- lines337_2015
