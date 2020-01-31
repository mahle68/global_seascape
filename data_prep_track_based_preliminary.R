#script for preparing all input data simultaneously. using water-crossing segments as sampling units
#Elham Nourani,
#Jan. 27. 2020. Radolfzell, Germany.
#Jan 31. 2020. Rstablish methodology with a small portion of the data. then apply to all.

library(multidplyr)
library(tidyverse)
library(lubridate)
library(move)
library(mapview)
library(rWind)
library(sf)
library(parallel)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")


load("R_files/land_15km.RData") #called land_15 km
load("R_files/land_0_60.RData") #called land_0_60... no buffer


land <- shapefile("/home/enourani/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp")

#land_0_60_prec <- st_set_precision(land_0_60,0.001)

alt_pts_temporal <- function(date_time,n_days) {
  #same year, same hour, only day changes
  alt_pts_before <- vector()
  alt_pts_after <- vector()
  
  for (i in 1:as.integer(n_days/2)) {
    alt_pts_before[i] <- as.character(date_time - days(i)) 
    
    alt_pts_after[i] <- as.character(date_time + days(i)) 
  }
  
  rbind( data.frame(dt = alt_pts_before, period =  "before", stringsAsFactors = F),
         data.frame(dt = alt_pts_after, period =  "after", stringsAsFactors = F))
}
land_1km <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp") %>%
  st_crop(y = c(xmin = -180, xmax = 180, ymin = 0, ymax = 60)) %>%
  st_transform(meters_proj) %>%
  st_buffer(dist = units::set_units(1000, 'm')) %>%
  st_transform(wgs) %>%
  st_union()

##### STEP 1: open data and take a sample #####
load("R_files/all_spp_unfiltered_updated_lc_0_removed.RData") #dataset; data prepared in all_data_prep_analyze

smpl <- dataset %>%
  group_by(species) %>%
  dplyr::select(track) %>% 
  sample_n(3)
smpl <- dataset %>% 
  filter(track %in% smpl$track)

coordinates(smpl)<-~location.long+location.lat
proj4string(smpl)<-wgs
smpl_sf <- st_as_sf(smpl) #convert to sf object



##### STEP 2: create an ocean layer #####
load("R_files/land_0_60.RData") #called land_0_60... no buffer

pts <-rbind(c(-180,0),c(180,0),c(180,60),c(-180,60),c(-180,0))
pol <- coords2Polygons(pts, ID = 1)
proj4string(pol) <- wgs

ocean <- mask(pol, )
pol <- st_polygon(list(rbind(c(-180,0),c(180,0),c(180,60),c(-180,60),c(-180,0))))  #create a polygon
  st_set_crs(pol) <- 4326

  pol <- st_polygon(list(rbind(c(-180,0),c(180,0),c(180,60),c(-180,60),c(-180,0)))) %>%  #create a polygon
    st_set_crs(wgs)


land_1km <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp") %>%
  st_crop(y = c(xmin = -180, xmax = 180, ymin = 0, ymax = 60)) %>%
  st_transform(meters_proj) %>%
  st_buffer(dist = units::set_units(1000, 'm')) %>%
  st_transform(wgs) %>%
  st_union()

##### STEP 3: remove tracks with no points over the sea #####
#open points over the sea
load("R_files/all_spp_spatial_filtered_updated.RData") # dataset_sea; created in all_data_prep_analyze... includes tracks that are now removed because of lc, etc.

#consider replacing the above with this
# dataset_sea_1km <- dataset %>% 
#   drop_na(c("location.long", "location.lat")) %>% 
#   st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
#   st_difference(land_1km) #change the buffer from 15 km to 1 km. maybe a track has one point near shore, but the sea-crossing is long
# 
# save(dataset_sea_1km, file = "R_files/")

dataset_tracks_sea <- dataset[dataset$track %in% dataset_sea$track,]



##### STEP 4: convert tracks to spatial lines #####
#convert tracks to lines

#only keep tracks that have at least three points
less_than_three_point <- dataset_sf %>% #find out which tracks have only one or two points. the two point tracks are only in East Asia
  group_by(track) %>% 
  summarise(n = length(track)) %>% 
  filter(n < 3) 

lines <- dataset_sf %>% 
  filter(!(track %in% less_than_three_point$track)) %>% 
  group_by(track) %>%
  arrange(date_time) %>% 
  summarize(species = head(species,1),do_union = F) %>% 
  st_cast("LINESTRING")

save(lines, file = "R_files/lines.RData")

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




#using multidplyr produces error. just use parallel
mycl <- makeCluster(detectCores() - 2)

clusterExport(mycl, c("lines", "land_1km")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(dplyr)
  library(sf)
  library(raster)
})

b <-Sys.time()
sea_lines_no_buffer_ls <- parLapply(cl = mycl, X = split(lines,lines$track), fun = st_difference, land_1km) #this didn't finish after 21 hours. somthing wrong with the code
Sys.time()-b

stopCluster(mycl)

save(sea_lines_no_buffer_ls, file = "R_files/sea_lines_no_bufffer_ls.RData")

#convert multilinestrings to linestrings
sea_seg <- sea_lines %>% 
  st_cast("MULTILINESTRING") %>% 
  st_cast("LINESTRING")

##### STEP 6: filter sea-crossing segments: length, dist to coast #####





#with the 15 km buffers
# sea_lines <- lines %>% #this takes forever, do it on the cluster using multidplyr ... error: no applicable method for 'st_difference' applied to an object of class "multidplyr_party_df"
#   st_difference(land_15km)
# 
# save(sea_lines, file = "R_files/sea_lines_15_km.RData")


  
save(sea_seg, file = "R_files/sea_segments_15_km.RData")

#without the 15 km buffer
sea_lines_no_buffer <- lines %>% 
  st_difference(land_1km)

save(sea_lines_no_buffer, file = "R_files/sea_lines_1km.RData")

sea_seg_no_buffer <- sea_lines_no_buffer %>% 
  st_cast("MULTILINESTRING") %>% 
  st_cast("LINESTRING")

save(sea_seg_no_buffer, file = "R_files/sea_segments_1km.RData")

# #using multidplyr produces error. just use parallel
# mycl <- makeCluster(detectCores() - 2)
# 
# clusterExport(mycl, c("lines", "land_0_60")) #define the variable that will be used within the function
# 
# clusterEvalQ(mycl, {
#   library(dplyr)
#   library(sf)
# })
# 
# sea_lines_no_buffer_ls <- parLapply(cl = mycl, X = split(lines,lines$track), fun = st_difference, land_0_60)
# 
# stopCluster(mycl)
# 
# save(sea_lines_no_buffer_ls, file = "R_files/sea_lines_no_bufffer_ls.RData")
# 

  
##### STEP 4: filter sea-crossing segments based on length and distance from coast #####



##### STEP 4: temporal filter for autocorrelation #####

load("R_files/all_spp_spatial_filtered_updated.RData") #named dataset_sea
dataset_sea <- dataset_sea %>% 
  dplyr::select(2,3,4,14,15,16,33,38,39,40)

#check for duplicated time-stamps
rows_to_delete <- sapply(getDuplicatedTimestamps(x = as.factor(dataset_sea$track),timestamps = dataset_sea$date_time,sensorType = "gps"),"[",2) #get the second row of each pair of duplicate rows
dataset_sea <- dataset_sea[-rows_to_delete,]

#thin the data to have one-hour difference?
#create a move object
dataset_sea <- dataset_sea[order(dataset_sea$track, dataset_sea$date_time),]
mv <- move(x = st_coordinates(dataset_sea)[,"X"],y = st_coordinates(dataset_sea)[,"Y"],time = dataset_sea$date_time,
           data = as.data.frame(dataset_sea,row.names = c(1:nrow(dataset_sea))),animal = dataset_sea$track,proj = wgs)

#start the cluster....dataset is not big enough to need cluster computing

#clusterExport(mycl, "mv") #define the variable that will be used within the function

#clusterEvalQ(mycl, {
#  library(move)
#  library(lubridate)
#  library(dplyr)
#  library(raster)
#})

sp_obj_ls <- lapply(split(mv),function(one_track){ #for each track within the group
  
  thinned_track <- one_track %>%
    thinTrackTime(interval = as.difftime(1, units = 'hours'), 
                  tolerance = as.difftime(15, units = 'mins'))
  
  #convert back to a move object (from move burst)
  thinned_track <- as(thinned_track,"Move")
  thinned_track$track <- one_track@idData$track #reassign the track
  thinned_track$species <- one_track@idData$species

  thinned_track
})

#stopCluster(mycl)

sp <- do.call(rbind,sp_obj_ls)
sf <- st_as_sf(sp)

##### STEP 5: produce alternative points in time #####

#timestamps at midnight cause a problem becuase hour becomes NA. add 5 minutes to all timestamps at midnight
sf <- sf %>% 
  mutate(date_time = ifelse(hour(date_time) == 0 & minute(date_time) == 0, as.character(date_time + minutes(5)), as.character(date_time))) %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))


#for each point, create alternative points a week before and a week after the observed point. year and hour dont change.
alts <- sf %>%
  mutate(#date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
    obs_id = row_number()) %>%
  slice(rep(row_number(),15)) %>% #copy each row 60 times. 1 used, 60 alternative
  arrange(date_time) %>%
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 15 == 0, 1, 0)),  #assign used and available values
         lon = st_coordinates(.)[,"X"],
         lat = st_coordinates(.)[,"Y"]) %>% 
  st_drop_geometry()

mycl <- makeCluster(detectCores() - 2)

clusterExport(mycl, c("alts", "alt_pts_temporal")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(lubridate)
  library(dplyr)
  library(sf)
})


alt_ls <- parLapply(cl = mycl, X = split(alts,alts$obs_id),fun = function(x){ #didnt manage to write this part using dplyr and purrr
  alt_times <- alt_pts_temporal(x$date_time[1],14)
  x %>%
    mutate(timestamp = as.POSIXct(strptime(c(as.character(x$date_time[1]),alt_times$dt),format = "%Y-%m-%d %H:%M:%S"),tz = "UTC",),
           period = c("now",alt_times$period))
})

stopCluster(mycl)

alt_cmpl <- do.call(rbind,alt_ls)

alt_cmpl <- alt_cmpl %>% 
  dplyr::select(date_time,month, season, timestamp, zone, track, species, obs_id, used, lon, lat, period)

save(alt_cmpl,file = "R_files/all_spp_temp_sp_filtered_15km_alt_14days.RData") #with 5 minutes added to 00:00; extra columns removed

##### STEP 6: annotate all points #####

load("R_files/all_spp_temp_sp_filtered_15km_alt_14days.RData") #called alt_cmpl

#prep for track annotation on movebank
alt_cmpl_mb <- alt_cmpl %>%
  mutate(timestamp = paste(as.character(timestamp),"000",sep = "."))

#rename columns
colnames(alt_cmpl_mb)[c(10,11)] <- c("location-long","location-lat")

write.csv(alt_cmpl_mb,"R_files/all_spp_temp_sp_filtered_15km_alt_14days.csv") #with 5 minutes added to 00:00

#downloaded from movebank
ann <- read.csv("movebank_annotation/all_spp_temp_sp_filtered_15km_alt_14days.csv-6653387681147029371.csv", stringsAsFactors = F) %>%
  drop_na() %>%# NA values are for the 2019 tracks. with a transition to ERA5, I should be able to use this data as well
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>%
  mutate(wspd = uv2ds(u925,v925)[,2],
         wdir = uv2ds(u925,v925)[,1],
         year = year(timestamp)) %>% 
  mutate(delta_t = sst - t2m)

#identify observed id's with less than 15 observations (the rest were removed because they produced NAs in the annotation step)
NA_obs_ids <- ann %>% 
  group_by(obs_id) %>% 
  summarise(count = n()) %>% 
  filter(count < 15) %>% 
  .$obs_id

ann <- ann %>% 
  filter(!(obs_id %in% NA_obs_ids))


save(ann, file = "R_files/all_spp_temp_sp_filtered_15km_alt_ann_14days.RData")


##### MAP all data #####
load("R_files/all_spp_unfiltered.RData") #called dataset

dataset <- dataset %>% 
  mutate(color = ifelse(species == "OHB", "cornflowerblue",
                        ifelse(species == "GFB","darksalmon",
                               ifelse(species == "PF", "firebrick1",
                                      "darkseagreen3"))))

X11(width = 15, height = 8)
tiff("/home/enourani/ownCloud/Work/safi_lab_meeting/presentation_jan17/all_tracks.tiff", width = 15, height = 8, units = "in", res = 500)
maps::map("world",fill = TRUE, col = "grey30", border = F)
points(dataset$location.long,dataset$location.lat, col= alpha(dataset$color,0.3), pch = 16, cex = 0.5)
legend(x = -170, y = -40, legend = c("Grey-faced buzzard","Osprey","Oriental honey buzzard","Peregrine falcon"),
       col = c("darksalmon","darkseagreen3","cornflowerblue","firebrick1"), pch = 16, bty = "n", cex = 0.9)
abline(h = 0, lty = 2,lwd = 0.2, col = "grey50")
abline(h = 30, lty = 2, lwd = 0.2, col = "grey50")
abline(h = 60, lty = 2, lwd = 0.2, col = "grey50")
text(x = -175, y = 32, "30° N", col = "grey50", cex = 0.7)
text(x = -175, y = 62, "60° N", col = "grey50", cex = 0.7)

dev.off()
