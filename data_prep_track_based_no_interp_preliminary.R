#script for preparing all input data simultaneously. using water-crossing segments as sampling units
#Elham Nourani,
#Jan. 27. 2020. Radolfzell, Germany.
#Jan 31. 2020. Establish methodology with a small portion of the data. then apply to all.
#Feb. 5. 2020. skip the interpolation step. it is not a good idea to assign arbitrary timestamps to points... the study is temporal. keep the temporal aspect as empirical as possible. 

library(multidplyr)
library(tidyverse)
library(lubridate)
library(move)
library(mapview)
library(rWind)
library(sf)
library(parallel)
library(lutz)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

alt_pts_temporal_hr <- function(date_time,n_alts, n_days, n_hours) { #n_alts is the number of alternative points. n_days is the number of alternative days
  #same year, hour changes
  alt_pts_before <- vector()
  alt_pts_after <- vector()
  #n_other_days <- n_alts - 23
  
  for (i in 1:as.integer((n_days/2)*n_hours)) {
    alt_pts_before[i] <- as.character(date_time - hours(i)) 
    
    alt_pts_after[i] <- as.character(date_time + days(i)) 
  }
  
  rbind( data.frame(dt = alt_pts_before, period =  "before", stringsAsFactors = F),
         data.frame(dt = alt_pts_after, period =  "after", stringsAsFactors = F))
}


##### STEP 1: create an ocean layer #####
load("R_files/land_0_60.RData") #called land_0_60... no buffer
pts <-rbind(c(-180,0),c(180,0),c(180,60),c(-180,60),c(-180,0))
pol <- coords2Polygons(pts, ID = 1)
proj4string(pol) <- wgs

ocean <- st_as_sf(pol) %>% 
  st_difference(land_0_60)
save(ocean,file = "R_files/ocean_0_60.RData")
load("R_files/ocean_0_60.RData")

##### STEP1.2: create two zone polygons #####
pts_twz <-rbind(c(-180,0),c(180,0),c(180,30),c(-180,30),c(-180,0))
pol_twz <- coords2Polygons(pts_twz, ID = 1)
proj4string(pol_twz) <- wgs
twz_sf <- st_as_sf(pol_twz)

pts_tmpz <-rbind(c(-180,30),c(180,30),c(180,60),c(-180,60),c(-180,30))
pol_tmpz <- coords2Polygons(pts_tmpz, ID = 1)
proj4string(pol_tmpz) <- wgs
tmpz_sf <- st_as_sf(pol_tmpz)

##### STEP 2: open data and take a sample #####
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
save(smpl_sf, file = "R_files/sample_tracks_pts.RData")

##### STEP 3: remove tracks with no points over the sea #####
load("R_files/sample_tracks_pts.RData")

smpl_sea <- smpl_sf %>%  #consider a one km buffer to remove tracks following the coast.
  st_intersection(ocean)

save(smpl_sea,file = "R_files/sample_tracks_pts_sea.RData")

#extract the tracks from the general filter of tracks over water
load ("R_files/sea_lines_1km.RData") # called sea_lines_no_buffer. created in data_prep_track_based
load("R_files/sample_tracks_pts_sea.RData")


lines <- sea_lines_no_buffer %>% 
  filter(track %in% smpl_sea$track)

##### STEP 4: filter segments #####

#remove segments shorter than 30 km
segs_30_km <- lines %>% 
  st_cast("MULTILINESTRING") %>% 
  st_cast("LINESTRING") %>% #convert from track to segment
  #st_transform(meters_proj) %>% 
  mutate(length = as.numeric(st_length(.))) %>% 
  filter(length >= 30000)

#remove straight lines with only two points
segs_filtered <- segs_30_km %>% 
  mutate(n = npts(.,by_feature = T)) %>% 
  #rowwise() %>% 
  #mutate(n = dim(st_coordinates(geometry))[1]) %>%  #count the number of points per line
  #ungroup() %>% 
  filter(!(length >= 30000 & n <=2))  #filter segments longer than 100 km and with only two points
  
#assign zone to each segment
segs_filtered$twz <- as.numeric(st_within(segs_filtered,twz_sf))
segs_filtered$tmpz <- as.numeric(st_within(segs_filtered,tmpz_sf))

segs_filtered <- segs_filtered %>% 
  mutate(zone = ifelse(is.na(twz) == TRUE & is.na(tmpz) ==T, "both",
                       ifelse(is.na(twz) == TRUE & tmpz == 1, "tmpz",
                              "twz"))) %>% 
  dplyr::select(-c("twz","tmpz"))

save(segs_filtered, file = "R_files/sample_filtered_segs.RData")

##### STEP 5: annotate with date-time #####

#convert segments to points
segs_pts <- segs_filtered %>% 
  mutate(seg_id = seq(1:nrow(.))) %>% 
  st_cast("POINT")

segs_pts_dt <- lapply(split(segs_pts,segs_pts$track), function(x){
  track <- x$track[1]
  x_df <- as(x,"Spatial") %>% 
    as.data.frame() %>% 
    rename(location.long = coords.x1,
           location.lat = coords.x2) %>% 
    mutate(season = NA,
           date_time = NA)
  
  data_df <- dataset[dataset$track == track,]
  
  for (i in 1:nrow(x_df)){
    pt <- x_df[i,]
    intersect <- data_df[round(data_df$location.long,1) == round(pt$location.long,1) & round(data_df$location.lat,1) == round(pt$location.lat,1),]
    if (nrow(intersect) == 1){ #if the point overlaps with an observed point over the track, assign the date_time and season
      x_df[i,"season"] <- intersect$season
      x_df[i,"date_time"] <- as.character(intersect$date_time+ min(1))
      
    } else { #if the point does not overlap with an observed point over the track, find the nearest neighbor and assign the date_time and season
      nearest <- data_df %>% 
        rowwise() %>% 
        mutate(lon_diff = abs(location.long-pt$location.long),
               lat_diff = abs(location.lat-pt$location.lat)) %>% 
        mutate(closeness = lon_diff + lat_diff) %>% 
        arrange(closeness) %>% 
        slice(1) #extract the first row
      
      x_df[i,"season"] <- nearest$season
      x_df[i,"date_time"] <- as.character(nearest$date_time + min(1)) #to make sure even 00:00 is retained.
    }
  }
  x_df
}) %>% 
  reduce(rbind)

save(segs_pts_dt, file = "R_files/sample_segs_pts_dt.RData")

##### STEP 6: create alternative tracks in time #####

hours_to_add <- c(0,cumsum(rep(1,672/2)),cumsum(rep(-1,(672/2))))

pts_alt <- segs_pts_dt %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         tz = tz_lookup_coords(location.lat,location.long)) %>% 
  as.data.frame() %>% 
  mutate(obs_id = row_number()) %>% 
  slice(rep(row_number(),673)) %>%  #paste each row 695 times for alternative points: 28days *24 hours + 23 hours in observed day
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 673 == 0, 1, 0))) %>% 
  group_by(obs_id) %>% 
  arrange(obs_id) %>% 
  mutate(hours_to_add = hours_to_add) %>% 
  mutate(alt_date_time = date_time + hours(hours_to_add)) %>%  #use hours to add as an id for alternative segments
  ungroup()
    
save(pts_alt, file = "R_files/sample_alt_pts_alt_time.RData")


##### STEP 7: annotate all points #####

load("R_files/sample_alt_pts_alt_time.RData") #called pts_alt

#prep for track annotation on movebank
pts_alt_mb <- pts_alt %>%
  mutate(timestamp = paste(as.character(alt_date_time),"000",sep = "."))

#rename columns
colnames(pts_alt_mb)[c(7,8)] <- c("location-long","location-lat")

write.csv(pts_alt_mb,"R_files/sample_pts_mb.csv") 


##########################################
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
