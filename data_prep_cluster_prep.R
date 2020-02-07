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
library(sf)
library(raster)
library(parallel)
#library(parallel)
#library(lutz)
#library(move)
#library(mapview)
#library(rWind)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

#load files. make sure they are all stored in the working directory
load("twz_sf.RData") #twz_sf
load("tmpz_sf.RData") #tmpz_sf
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset
load("land_0_60.RData") #land_0_60
load("ocean_0_60.RData") #ocean

##### STEP 1: find tracks with at least some points over the sea #####
dataset_sp <- dataset
coordinates(dataset_sp) <- ~location.long + location.lat
proj4string(dataset_sp) <- wgs


intersect(rr,ocean_sp)
#dataset_sf <- dataset_sp %>% 
#  st_as_sf()

#pts_sea <- erase(dataset_sp,as(ocean,"Spatial"))



# pts_sea <- dataset_sf %>%  
#   st_difference(ocean)
# 
# save(pts_sea,file = "all_tracks_pts_sea.RData")

##### STEP 1: convert tracks to spatial lines and remove portions over land #####
ocean_sp <- as(ocean,"Spatial")
land_sp <- as(land_0_60,"Spatial")

track_ls<-split(dataset,dataset$track)
track_ls<-track_ls[lapply(track_ls,nrow)>0]


mycl <- makeCluster(9) #total number of tracks is 369, so 41 will be sent to each core

clusterExport(mycl, c("track_ls", "land_sp","ocean_sp","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(tidyverse)
  library(lubridate)
  library(sf)
  library(raster)
})

Lines_ls<-lapply(track_ls,function(x){
  #find out if the track has any points over water
  over_sea <- intersect(x,ocean_sp) #track_ls needs to be spatial for this to work
  #if the track has any point over water, convert to spatial line and subset for sea
  if(nrow(over_sea) != 0){
  line<-SpatialLines(list(Lines(Line(matrix(x@coords,ncol=2)), ID=x$track[1])),proj4string = wgs)
  line_sea <- erase(line,land_sp)
  line_sea$track <- x$track[1]
  } else {
    line_sea <- NA
  }
  
  line_sea
})

stopCluster(mycl)

save(Lines_ls,file = "Lines_ls_no_land.RData")




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

rr <- dataset_sp[dataset_sp$track == "84430_2009_autumn",]
rr2 <- dataset_sp[dataset_sp$track == "90757_2009_autumn",]

int <- intersect(rr,ocean_sp)
int2 <- intersect(rr2,ocean_sp)
