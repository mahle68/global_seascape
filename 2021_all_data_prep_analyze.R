#script for preparing all input data simultaneously. previously done separately for each species/region.
#Elham Nourani,
#Feb. 9 2021. Radolfzell, Germany.
#no limit to northern hemisphere. in the current data, I will only have one point for the eleonora's falcon from E. Africa to Madagascar. 
#only redo data for 2019 onwards and for the latitudinal zones that were not looked before.
#update Mar 1: in the first attempt, couldn't find the coords2Lines function, so used sf to create lines and segmentate the seacrossing tracks. that somehow didn't return the
#correct lines.... so, using Orcs::coords2Lines, with the original code from track_based_prep_analyze_daily.R

library(tidyverse)
library(readxl) #read_excel()
library(lubridate)
library(move)
library(mapview)
library(rWind)
library(sf)
library(parallel)
library(Orcs) #new home of coords2lines


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
wgs_sf <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

load("R_files/2021/ocean.RData")
load("R_files/2021/land_no_buffer.RData")

##### STEP 0: prep spatial layers #####
land_15km <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp") %>% 
  st_crop(y = c(xmin = -180, xmax = 180, ymin = -38, ymax = 70)) %>% 
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(15000, 'm')) %>% 
  st_transform(wgs) %>% 
  st_union()

save(land_15km, file = "R_files/2021/land_15km.RData")

land_no_buffer <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp") %>% #this should be used in place of land_0_60
  st_crop(y = c(xmin = -180, xmax = 180, ymin = -38, ymax = 70)) %>% 
  st_transform(wgs) %>% 
  st_union()

save(land_no_buffer, file = "R_files/2021/land_no_buffer.RData") #called land_no_buffer

#create an ocean layer
ocean <- st_polygon(list(rbind(c(-180,-38), c(180,-38), c(180,70), c(-180,70),c(-180,-38)))) %>% 
  st_sfc(crs = wgs) %>% 
  st_difference(land_no_buffer)

save(ocean, file = "R_files/2021/ocean.RData") #called ocean

#ocean_sp <- as(ocean,"Spatial")
#land_sp <- as(land_no_buffer,"Spatial")
land_b <- buffer(land_sp,width = 0.001)


##### STEP 1: read in the data and filter for adult birds and season and databse-specific filters #####

#read in meta-data for peregrine falcon and american osprey. extract IDs for adult birds
pf_ad <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey_Americas/Peregrines Ivan.csv", stringsAsFactors = F) %>% 
  filter(Age == "ad")
ao_ad <-read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Osprey_Americas/ROB mig data190411.csv", stringsAsFactors = F) %>% 
  filter(Age2 == "a")

#also assign date_time, year, month, track, species, location.lat, location.long

OHB_files <- list.files("data/Oriental_honey_buzzard",pattern = ".xls",full.names = T)
OHB <- lapply(OHB_files,read_excel,1,col_types = c("numeric","date","numeric","numeric","numeric","skip","text",rep("numeric",8))) %>%
  reduce(full_join) %>%
  rename(date_time = 'date(gmt)',lon = longitud,lat = latitude) %>%
  mutate(yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other"))) %>%  #11 Sep-20 Oct; spring between 1-5 May
  mutate(track = paste(ptt,year,season,sep = "_"),
         species = "OHB",
         ind = as.character(ptt),
         sensor.type = "ptt") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season == "autumn" & #no sea-crossing in spring
           class %in% c("1","2","3")) #filter for location classes

OHB_GPS_aut <- read.csv("data/Tracking_of_the_migration_of_Oriental_Honey_Buzzards.csv") %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "OHB",
         yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other")),  #11 Sep-20 Oct; spring between 1-5 May
         track = paste(tag.local.identifier, year,season,sep = "_"),
         ind = as.character(tag.local.identifier)) %>% 
  filter(season == "autumn")

###
GFB_files <- list.files("data/Grey_faced_buzzard/",pattern = ".csv",full.names = T)
GFB <- lapply(GFB_files,read.csv,stringsAsFactors = F) %>%
  reduce(full_join) %>% #is locdate is in UTC
  mutate(locdate,date_time = as.POSIXct(strptime(locdate,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate(month = month(date_time),
         year = year(date_time),
         season = ifelse(month %in% c(3,4),"spring",ifelse(month %in% c(10),"autumn","other"))) %>% 
  mutate(track = paste(platform,year,season,sep = "_"),
         species = "GFB",
         ind = as.character(platform),
         sensor.type = "ptt") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season %in% c("spring","autumn") &
           class %in% c("1","2","3")) #filter for location classes

###make sure migration season is correctly defined
PF <- read.csv("data/LifeTrack Peregrine falcon_2021.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,36,38,39) %>% #remove columns that are not needed
  filter(individual.local.identifier %in% pf_ad$animal.id) %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "PF",
         season = ifelse (month %in% c(9,10), "autumn", "other")) %>% 
  mutate(track = paste(tag.local.identifier, year, season,sep = "_"),
         ind = as.character(tag.local.identifier)) %>% 
  filter(season != "other")

OE <- read.csv("data/Osprey in Mediterranean (Corsica, Italy, Balearics)_2021.csv", stringsAsFactors = F) %>% 
  dplyr::filter(grepl("ad",individual.local.identifier,, ignore.case = T) & !grepl("juv",individual.local.identifier,, ignore.case = T)) %>% #extract adult data
  dplyr::select(1,3:5,16,35:37) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "O",
         season = ifelse(month %in% c(2:4),"spring",ifelse (month %in% c(8:10), "autumn","other"))) %>% 
  mutate(track = paste(tag.local.identifier, year,season, sep = "_"),
         ind = as.character(tag.local.identifier)) %>% 
  filter(season != "other")

OA <- read.csv("data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,13,48:52) %>% #remove columns that are not needed
  filter(sensor.type == "gps" | sensor.type == "argos-doppler-shift" & argos.lc %in% c("1","2","3"),
         individual.local.identifier %in% ao_ad$Bird) %>%
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "O",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  mutate(track = paste(tag.local.identifier, year,season,sep = "_"),
         ind = as.character(tag.local.identifier)) %>% 
  filter(season != "other")

EF_aut <- read_excel("data/eleonoras_falcon.xlsx", sheet = 2) %>% 
  mutate(date_time = as.POSIXct(strptime(`timestamp (UTC)`,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(timestamp = `timestamp (UTC)`) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "EF",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  mutate(track = paste(ID.individual, year,season,sep = "_"),
         ind = as.character(ID.individual),
         sensor.type = "gps") %>% 
  filter(season == "autumn")

##### STEP 2: merge all data, assign zone #####

dataset <- list(OHB,GFB, PF, OE, OA, EF_aut, OHB_GPS_aut) %>% 
  reduce(full_join, by = c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species", "ind", "sensor.type")) %>% 
  mutate(zone = ifelse(between(location.lat, 0, 30) | between(location.lat, 0, -30), "tradewind",
                       ifelse(between(location.lat, 30,60) | between(location.lat, -30,-60), "temperate",
                              ifelse(between(location.lat, -30,30), "tropical",
                              "arctic")))) %>% 
  dplyr::select(c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species", "zone", "ind", "sensor.type")) %>% 
  as.data.frame()

save(dataset, file = "R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData")

##### STEP 3: convert tracks to lines ##### 

load("R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #called dataset

# #only keep gps data (ie.get rid of GFB)
# dataset <- dataset %>% 
#   filter(sensor.type == "gps")
# 
# ocean_sp <- as(ocean,"Spatial")
# land_sp <- as(land_no_buffer,"Spatial")
# land_b <-buffer(land_sp, width=0.001)
# 
# dataset <- dataset[complete.cases(dataset$location.long),]
# 
# coordinates(dataset) <- ~ location.long + location.lat
# proj4string(dataset) <- wgs
# 
# track_ls <- split(dataset,dataset$track)
# track_ls <- track_ls[lapply(track_ls,nrow) > 1] #remove tracks with one point
# 
# (b <- Sys.time())
# mycl <- makeCluster(10) #total number of tracks is 369, so 41 will be sent to each core
# 
# clusterExport(mycl, c("track_ls", "land_sp","ocean_sp","wgs")) #define the variable that will be used within the function
# 
# clusterEvalQ(mycl, {
#   library(raster)
#   library(mapview)
#   library(Orcs)
# })
# 
# Lines_ls <- parLapply(mycl, track_ls, function(x){
#   #find out if the track has any points over water
#   over_sea <- intersect(x,ocean_sp) #track_ls needs to be spatial for this to work
#   #if the track has any point over water, convert to spatial line and subset for sea
#   if(nrow(over_sea) != 0){
#     x <- x[order(x$date_time),]
#     line <- coords2Lines(x@coords, ID = x$track[1], proj4string = wgs)
#     line_sea <- erase(line,land_sp)
#     line_sea$track <- x$track[1]
#   } else {
#     line_sea <- NA
#   }
#   
#   line_sea
# })
# 
# 
# Sys.time() - b #takes 45 min
# 
# stopCluster(mycl)
# 
# save(Lines_ls,file = "R_files/2021/Lines_ls_no_land.RData") 
# 
# ##### STEP 2: break up tracks into sea-crossing segments and filter #####
# 
# #remove elements with 0 elements (tracks with no sea-crossing)
# Lines_ls_no_na <- Lines_ls[lapply(Lines_ls,is.na) == FALSE] 
# 
# #only keep the track column (some objects have an ID column)
# Lines_ls_no_na <- lapply(Lines_ls_no_na,"[",,"track")
# 
# #convert to one object
# lines <- do.call(rbind,Lines_ls_no_na)
# 
# #filter segments
# segs_filtered <- st_as_sf(lines) %>% #convert to sf object
#   st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
#   mutate(length = as.numeric(st_length(.)),
#          n = npts(.,by_feature = T)) %>% 
#   filter(n > 2 & length >= 30000) #remove sea-crossing shorter than 30 km and segment with less than 2 points 
# 
# segs_filtered$track <- as.character(segs_filtered$track)
# 
# save(segs_filtered,file = "Segs_no_land_filtered.RData") 


#############################################3

 dataset_sea <- dataset %>%
   filter(year >= 2009 & season == "autumn" & sensor.type == "gps") %>% 
   drop_na(c("location.long", "location.lat")) %>% 
   st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
   group_by(track) %>% 
   filter(n() > 1) %>%  #remove tracks with only one point
   arrange(date_time) %>% 
   summarise(track = head(track,1),
             species = head(species,1),
             zone = head(zone, 1), do_union = F) %>% 
   st_cast("LINESTRING")
 
 save(dataset_sea, file = "R_files/2021/all_spp_2009_2020_complete_lines.RData")


##### STEP 4: subset for tracks over the sea and assign segments ##### 

load("R_files/2021/all_spp_2009_2020_complete_lines.RData") #dataset_sea

segs_filtered <- dataset_sea %>% 
  st_difference(land_no_buffer)

save(segs_filtered, file = "R_files/2021/all_spp_2009_2020_all_segments.RData")



###################### for some reason what is below doesnt manage to filter out short sea crossings
segs <- segs_filtered %>% 
  as("Spatial") %>% 
  st_as_sf() %>% 
  st_cast("LINESTRING") %>% #create one line object for each segment, instead of line multilines that include multiple segments for each track
  mutate(length = as.numeric(st_length(.)),
         n = npts(., by_feature = T)) %>% 
  filter(n > 2 & length >= 30000)  #remove sea-crossing shorter than 30 km and segment with less than 2 points 

save(segs, file = "R_files/2021/all_spp_2009_2020_filtered_segments.RData")


##### STEP 5: annotate with date-time #####

#open data points and segments
load("R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #called dataset
load("R_files/2021/all_spp_2009_2020_filtered_segments.RData") #segs

#only keep data points from tracks that are represented in the over-sea segments (this will also get rid of spring tracks)
points_with_segs <- dataset %>% 
  drop_na(location.long) %>% 
  filter(track %in% unique(segs$track)) 


mycl <- makeCluster(4)
clusterExport(mycl, c("points_with_segs", "segs", "wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(raster)
  library(tidyverse)
})

(b <- Sys.time())

#points_oversea <- parLapply(mycl, split(points_with_segs, points_with_segs$track), function(x){

points_oversea <- lapply(split(points_with_segs, points_with_segs$track), function(x){

  seg <- segs[segs$track == x$track[1],]

  oversea <- x %>%
    st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>%
    st_intersection(seg, tolerance = 0.0001) %>%
    dplyr::select(-c(track.1, zone.1, species.1))

  oversea
}) %>%
  reduce(rbind)

Sys.time() - b #50 secs

stopCluster(mycl)


save(points_oversea, file = "R_files/2021/all_spp_2009_2020_overwater_points.RData")


# ##### STEP 5: annotate with date-time #####
# 
# #convert segments to points
# segs_pts <- segs %>% 
#   mutate(seg_id = seq(1:nrow(.))) %>% 
#   st_cast("POINT")
# 
# #create a buffer around the dataset points to make polygons. then overlay
# dataset_buff <- points_with_segs %>% 
#   drop_na() %>% 
#   st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
#   st_transform(meters_proj) %>% 
#   st_buffer(dist = units::set_units(10, 'm')) %>% 
#   st_transform(wgs) 
# 
# save(dataset_buff,file = "R_files/2021/dataset_10m_buffer.RData")
# 
# load("R_files/2021/dataset_10m_buffer.RData")
# 
# #for each segs_pts point, find the index of the dataset_buff polygon that it intersects, then extract that row from dataset and add to segs_pts
# mycl <- makeCluster(9) 
# 
# clusterExport(mycl, c("segs_pts", "dataset_buff")) #define the variable that will be used within the function
# 
# clusterEvalQ(mycl, {
#   library(sf)
#   library(raster)
#   library(tidyverse)
# })
# 
# b <- Sys.time()
# 
# segs_ann <- parLapply(mycl,split(segs_pts,segs_pts$track), function(x){ #separate by track first to break up the job into smaller chunks
#   
#   data <- dataset_buff[dataset_buff$track == x$track[1],]
#   #track_ann <- apply(x,1,function(y){ #for each point on the track
#   #x2 <- list()
#   track_ann <- lapply(split(x,rownames(x)), function(y){ #for each point on the track
#     #for (i in 1:nrow(x)){
#     #   y <- x[i,]
#     inter <- st_intersection(y,data)
#     
#     if(nrow(inter) == 0){ #if there are no intersections, find the nearest neighbor
#       nearest <- data[st_nearest_feature(y,data),]
#       # x$date_time[i] <- as.character(nearest$date_time)
#       # x$season[i] <- nearest$season
#       # x$species[i] <- nearest$species
#       
#       y <- y %>% 
#         full_join(st_drop_geometry(nearest))
#       y
#     } else { #if there is an intersection, just return the intersection result
#       # x$date_time[i] <- as.character(inter$date_time)
#       # x$season[i] <- inter$season
#       # x$species[i] <- inter$species
#       inter %>% 
#         dplyr::select(-track.1)
#     }
#     #}
#   }) %>% 
#     reduce(rbind)
#   
#   track_ann
#   
# }) %>% 
#   reduce(rbind)
# Sys.time() - b
# 
# stopCluster(mycl)
# 
# save(segs_ann, file = "segs_dt.RData")



##### STEP 6: append eleonora's falcon from spain #####

EF_spain <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Eleonoras_falcon/EF_autumn-sea-crossings_RESAMPLED-half-hourly.csv", stringsAsFactors = F) %>% 
  mutate(date_time = as.POSIXct(strptime(dt,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         species = "EF") %>% 
  rename(location.long = long,
         location.lat = lat,
         year = yr,
         track = tripID,
         ind = dev) %>% 
  mutate(month = month(date_time),
         zone = "tropical",
         season = "autumn",
         n = NA,
         length = NA,
         sensor.type = "gps") %>% 
  dplyr::select(-c("X", "alt", "dt")) %>% 
  filter(year >= 2013) %>% #should I remove this filter?!!!!!!!!!!!!!!
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_difference(land_no_buffer)


all_oversea <- rbind(points_oversea, EF_spain)

save(all_oversea, file = "R_files/2021/all_2009_2020_overwater_points.RData")


### some summary stats

all_oversea %>%  #this includes ptts
  group_by(species) %>% 
  summarize( n_ind = n_distinct(ind),
             n_tracks = n_distinct(track),
            start_yr = min(year),
            end_yr = max(year))
