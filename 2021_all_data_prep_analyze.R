#script for preparing all input data simultaneously. previously done separately for each species/region.
#Elham Nourani,
#Feb. 9 2021. Radolfzell, Germany.
#no limit to northern hemisphere. in the current data, I will only have one point for the eleonora's falcon from E. Africa to Madagascar. 
#only redo data for 2019 onwards and for the latitudinal zones that were not looked before.

library(tidyverse)
library(readxl) #read_excel()
library(lubridate)
library(move)
library(mapview)
library(rWind)
library(sf)
library(parallel)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")


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
land_b<-buffer(land_sp,width=0.001)


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
         species = "OHB") %>% 
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
         track = paste(tag.local.identifier, year,season,sep = "_")) %>% 
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
         species = "GFB") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season %in% c("spring","autumn") &
           class %in% c("1","2","3")) #filter for location classes

###make sure migration season is correctly defined
PF <- read.csv("data/LifeTrack Peregrine falcon_2021.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,38,39) %>% #remove columns that are not needed
  filter(individual.local.identifier %in% pf_ad$animal.id) %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "PF",
         season = ifelse (month %in% c(9,10), "autumn", "other")) %>% 
  mutate(track = paste(tag.local.identifier, year, season,sep = "_")) %>% 
  filter(season != "other")

OE <- read.csv("data/Osprey in Mediterranean (Corsica, Italy, Balearics)_2021.csv", stringsAsFactors = F) %>% 
  dplyr::filter(grepl("ad",individual.local.identifier,, ignore.case = T) & !grepl("juv",individual.local.identifier,, ignore.case = T)) %>% #extract adult data
  dplyr::select(1,3:5,16,35:37) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "O",
         season = ifelse(month %in% c(2:4),"spring",ifelse (month %in% c(8:10), "autumn","other"))) %>% 
  mutate(track = paste(tag.local.identifier, year,season, sep = "_")) %>% 
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
  mutate(track = paste(tag.local.identifier, year,season,sep = "_")) %>% 
  filter(season != "other")

EF_aut <- read_excel("data/eleonoras_falcon.xlsx", sheet = 2) %>% 
  mutate(date_time = as.POSIXct(strptime(`timestamp (UTC)`,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(timestamp = `timestamp (UTC)`) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "EF",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  mutate(track = paste(ID.individual, year,season,sep = "_")) %>% 
  filter(season == "autumn")

##### STEP 2: merge all data, assign zone #####

dataset <- list(OHB,GFB, PF, OE, OA, EF_aut, OHB_GPS_aut) %>% 
  reduce(full_join, by = c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
  mutate(zone = ifelse(between(location.lat, 0, 30) | between(location.lat, 0, -30), "tradewind",
                       ifelse(between(location.lat, 30,60) | between(location.lat, -30,-60), "temperate",
                              ifelse(between(location.lat, -30,30), "tropical",
                              "arctic")))) %>% 
  dplyr::select(c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species", "zone")) %>% 
  as.data.frame()

save(dataset, file = "R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData")

##### STEP 3: convert tracks to lines ##### 

load("R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #called dataset

dataset_sea <- dataset %>%  
  filter(year >= 2013 & season == "autumn") %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  group_by(track) %>% 
  filter(n() > 1) %>%  #remove tracks with only one point
  arrange(date_time) %>% 
  summarise(track = head(track,1),
            species = head(species,1),
            zone = head(zone, 1), do_union = F) %>% 
  st_cast("LINESTRING")

save(dataset_sea, file = "R_files/2021/all_spp_2013_2020_complete_lines.RData")


##### STEP 4: subset for tracks over the sea and assign segments ##### 


segs_filtered <- dataset_sea %>% 
  st_difference(land_no_buffer) %>% 
  #st_as_sf(lines) %>% #convert to sf object
  #st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(n > 2 & length >= 30000) #remove sea-crossing shorter than 30 km and segment with less than 2 points 

save(segs_filtered, file = "R_files/2021/all_spp_2013_2020_all_segments.RData")

segs <- segs_filtered %>% 
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(n > 2 & length >= 30000)

save(segs, file = "R_files/2021/all_spp_2013_2020_filtered_segments.RData")

##### STEP 3: filter out points over land ##### 

load("R_files/2021/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData")

dataset_sea <- dataset %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_difference(land_no_buffer)

save(dataset_sea, file = "R_files/2021/all_spp_spatial_filtered_updated_no_buffer2.RData")

##### STEP 4: remove duplicated points #####

load("R_files/2021/all_spp_spatial_filtered_updated_no_buffer2.RData") #named dataset_sea

#check for duplicated time-stamps
rows_to_delete <- sapply(getDuplicatedTimestamps(x = as.factor(dataset_sea$track),timestamps = dataset_sea$date_time,sensorType = "gps"),"[",2) #get the second row of each pair of duplicate rows
dataset_sea <- dataset_sea[-rows_to_delete,]

##### STEP 5: convert tracks to spatial lines (to have sea-crossing segments) and remove portions over land #####
          
#redo all 2013 onwards. keep only autumn
new_data <- dataset_sea %>%
  #filter(year >= 2019 & zone %in% unique(dataset_sea$zone) | year %in% unique(dataset_sea$year) & zone %in% c("arctic", "tropical")) #this is not a good idea. I need the track that has points in the arctic and the tropics, not only the points
  filter(year >= 2013 & season == "autumn")

track_ls <- split(new_data,new_data$track)
track_ls <- track_ls[lapply(track_ls,nrow) > 1] 


load("R_files/2021/ocean.RData")
load("R_files/2021/land_no_buffer.RData")


#new way.
lines <- new_data %>% 
  group_by(track) %>% 
  filter(n() > 1) %>%  #remove tracks with only one point
  arrange(date_time) %>% 
  summarise(track = head(track,1),
            species = head(species,1),
            zone = head(zone, 1), do_union = F) %>% 
  st_cast("LINESTRING") %>% 
  st_difference(land_no_buffer)
  


mycl <- makeCluster(10) 

clusterExport(mycl, c("track_ls", "land_no_buffer","ocean","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(raster)
  library(sf)
  library(tidyverse)
})

(b <- Sys.time())
Lines_ls <- parLapply(mycl,track_ls,function(x){
  #find out if the track has any points over water
  #over_sea <- st_intersection(x,ocean) #track_ls needs to be spatial for this to work. update: track_ls is already subsetted for sea. this step seems redundant
  #if the track has any point over water, convert to spatial line and subset for sea
  #if(nrow(over_sea) != 0){
    line_sea <- x %>%
      arrange(date_time) %>% 
      summarise(track = head(track,1),
                species = head(species,1),
                zone = head(zone, 1), do_union = F) %>% 
      #summarise_at(c("track", "species", "zone"), ~ head(.,1), do_union = FALSE) %>% #do_union is important to prevent squiggly lines. but with summarize_at, still order of points is messed up
      st_cast("LINESTRING") %>% 
      st_difference(land_no_buffer)
    
  #} else {
   # line_sea <- NA
  #}
  
  line_sea
})


stopCluster(mycl)

Sys.time() - b # 5.767824


save(Lines_ls, file = "R_files/2021/Lines_ls_2019_20_arctic_tropical.RData")

##### STEP 2: break up tracks into sea-crossing segments and filter #####

load("R_files/2021/Lines_ls_2019_20_arctic_tropical.RData")

#remove elements with 0 elements (tracks with no sea-crossing)
Lines_ls_no_na <- Lines_ls[lapply(Lines_ls,length) > 0]

#only keep the track column (some objects have an ID column)
#Lines_ls_no_na <- lapply(Lines_ls_no_na,"[",,"track")

#convert to one object
lines <- do.call(rbind,Lines_ls_no_na)

#filter segments
segs_filtered <- st_as_sf(lines) %>% #convert to sf object
  st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(n > 2 & length >= 30000) #remove sea-crossing shorter than 30 km and segment with less than 2 points 

segs_filtered$season <- sapply(strsplit(segs_filtered$track, "_"), "[",3)

segs_filtered$track <- as.character(segs_filtered$track)

aut <- segs_filtered[segs_filtered$season == "autumn",]

##### OHB GPS DATA PREP #####

OHB_GPS_aut <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/data/Tracking_of_the_migration_of_Oriental_Honey_Buzzards.csv") %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "OHB",
         yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other")),  #11 Sep-20 Oct; spring between 1-5 May
         track = paste(tag.local.identifier, year,season,sep = "_")) %>% 
  filter(season == "autumn" & location.lat >= 0)

#from track_based_prep_analyze_daily.R
load("R_files/land_0_60.RData") #land_0_60
load("R_files/ocean_0_60.RData") #ocean


##### STEP 1: convert tracks to spatial lines and remove portions over land #####
ocean_sp <- as(ocean,"Spatial")
land_sp <- as(land_0_60,"Spatial")
land_b<-buffer(land_sp,width=0.001)

coordinates(OHB_GPS_aut) <- ~location.long + location.lat
proj4string(OHB_GPS_aut) <- wgs

track_ls <- split(OHB_GPS_aut,OHB_GPS_aut$track)
track_ls <- track_ls[lapply(track_ls,nrow)>1] 

b <- Sys.time()
mycl <- makeCluster(6) 

clusterExport(mycl, c("track_ls", "land_sp","ocean_sp","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(raster)
  library(mapview)
})

Lines_ls<-parLapply(mycl,track_ls,function(x){
  #find out if the track has any points over water
  over_sea <- intersect(x,ocean_sp) #track_ls needs to be spatial for this to work
  #if the track has any point over water, convert to spatial line and subset for sea
  if(nrow(over_sea) != 0){
    x <- x[order(x$date_time),]
    line<-coords2Lines(x@coords, ID=x$track[1],proj4string = wgs)
    line_sea <- erase(line,land_sp)
    line_sea$track <- x$track[1]
  } else {
    line_sea <- NA
  }
  
  line_sea
})

stopCluster(mycl)

Sys.time() - b # 2 min

##### STEP 2: break up tracks into sea-crossing segments and filter #####

#remove elements with 0 elements (tracks with no sea-crossing)
Lines_ls_no_na <- Lines_ls[lapply(Lines_ls,is.na) == FALSE] 

#only keep the track column (some objects have an ID column)
Lines_ls_no_na <- lapply(Lines_ls_no_na,"[",,"track")

#convert to one object
lines <- do.call(rbind,Lines_ls_no_na)

#filter segments
segs_filtered<- st_as_sf(lines) %>% #convert to sf object
  st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(n > 2 & length >= 30000) #remove sea-crossing shorter than 30 km and segment with less than 2 points 

segs_filtered$track <- as.character(segs_filtered$track)

##### STEP 3: annotate with date-time #####
#convert segments to points
segs_pts <- segs_filtered %>% 
  mutate(seg_id = seq(1:nrow(.))) %>% 
  st_cast("POINT")

#create a buffer around the dataset points to make polygons. then overlay
dataset_buff <- OHB_GPS_aut %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(10, 'm')) %>% 
  st_transform(wgs) 

#for each segs_pts point, find the index of the dataset_buff polygon that it intersects, then extract that row from dataset and add to segs_pts
mycl <- makeCluster(6) 
clusterExport(mycl, c("segs_pts", "dataset_buff")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(raster)
  library(tidyverse)
})

b <- Sys.time()

segs_ann_OHB <- parLapply(mycl,split(segs_pts,segs_pts$track), function(x){ #separate by track first to break up the job into smaller chunks
  data <- dataset_buff[dataset_buff$track == x$track[1],]
  #track_ann <- apply(x,1,function(y){ #for each point on the track
  #x2 <- list()
  track_ann <- lapply(split(x,rownames(x)), function(y){ #for each point on the track
    #for (i in 1:nrow(x)){
    #   y <- x[i,]
    inter <- st_intersection(y,data)
    
    if(nrow(inter) == 0){ #if there are no intersections, find the nearest neighbor
      nearest <- data[st_nearest_feature(y,data),]
      # x$date_time[i] <- as.character(nearest$date_time)
      # x$season[i] <- nearest$season
      # x$species[i] <- nearest$species
      
      y <- y %>% 
        full_join(st_drop_geometry(nearest))
      y
    } else { #if there is an intersection, just return the intersection result
      # x$date_time[i] <- as.character(inter$date_time)
      # x$season[i] <- inter$season
      # x$species[i] <- inter$species
      inter %>% 
        dplyr::select(-track.1)
    }
    #}
  }) %>% 
    reduce(rbind)
  
  track_ann
  
}) %>% 
  reduce(rbind)
Sys.time() - b

stopCluster(mycl)


save(segs_ann_OHB, file = "R_files/segs_OHB_dt.RData")


##### EF GPS DATA PREP #####

EF_aut <- read_excel("/home/enourani/ownCloud/Work/Projects/delta_t/data/eleonoras_falcon.xlsx", sheet = 2) %>% 
  mutate(date_time = as.POSIXct(strptime(`timestamp (UTC)`,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(timestamp = `timestamp (UTC)`) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "EF",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  mutate(track = paste(ID.individual, year,season,sep = "_")) %>% 
  filter(season == "autumn" & location.lat >= 0) #few tracks between african continent and madagascar.... maybe not worth including?... hmmm

#from track_based_prep_analyze_daily.R
load("R_files/land_0_60.RData") #land_0_60
load("R_files/ocean_0_60.RData") #ocean


##### STEP 1: convert tracks to spatial lines and remove portions over land #####
ocean_sp <- as(ocean,"Spatial")
land_sp <- as(land_0_60,"Spatial")
land_b<-buffer(land_sp,width=0.001)

coordinates(EF_aut) <- ~location.long + location.lat
proj4string(EF_aut) <- wgs

track_ls <- split(EF_aut,EF_aut$track)
track_ls <- track_ls[lapply(track_ls,nrow) > 1] 


Lines_ls <- lapply(track_ls, function(x){  
  #find out if the track has any points over water
  over_sea <- intersect(x,ocean_sp) #track_ls needs to be spatial for this to work
  #if the track has any point over water, convert to spatial line and subset for sea
  if(nrow(over_sea) != 0){
    x <- x[order(x$date_time),]
    line <- coords2Lines(x@coords, ID = x$track[1])
    proj4string(line) <- wgs
    line_sea <- erase(line,land_sp)
    line_sea$track <- x$track[1]
  } else {
    line_sea <- NA
  }
  
  line_sea
})



##### STEP 2: break up tracks into sea-crossing segments and filter #####

#remove elements with 0 elements (tracks with no sea-crossing)
Lines_ls_no_na <- Lines_ls[lapply(Lines_ls,is.na) == FALSE] 

#only keep the track column (some objects have an ID column)
Lines_ls_no_na <- lapply(Lines_ls_no_na,"[",,"track")

#convert to one object
lines <- do.call(rbind,Lines_ls_no_na)

#filter segments
segs_filtered<- st_as_sf(lines) %>% #convert to sf object
  st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(n > 2 & length >= 30000) #remove sea-crossing shorter than 30 km and segment with less than 2 points 

segs_filtered$track <- as.character(segs_filtered$track)

##### STEP 3: annotate with date-time #####
#convert segments to points
segs_pts <- segs_filtered %>% 
  mutate(seg_id = seq(1:nrow(.))) %>% 
  st_cast("POINT")

#create a buffer around the dataset points to make polygons. then overlay
dataset_buff <- EF_aut %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(10, 'm')) %>% 
  st_transform(wgs) 

#for each segs_pts point, find the index of the dataset_buff polygon that it intersects, then extract that row from dataset and add to segs_pts

segs_ann_EF <- lapply(split(segs_pts,segs_pts$track), function(x){ #separate by track first to break up the job into smaller chunks
  
  data <- dataset_buff[dataset_buff$track == x$track[1],]
  
  track_ann <- lapply(split(x,rownames(x)), function(y){ #for each point on the track
    inter <- st_intersection(y,data)
    if(nrow(inter) == 0){ #if there are no intersections, find the nearest neighbor
      nearest <- data[st_nearest_feature(y,data),]
      y <- y %>% 
        full_join(st_drop_geometry(nearest))
      y
    } else { #if there is an intersection, just return the intersection result
      inter %>% 
        dplyr::select(-track.1)
    }
  }) %>% 
    reduce(rbind)
  
  
  
  track_ann
  
}) %>% 
  reduce(rbind)


save(segs_ann_EF, file = "R_files/segs_EF_dt.RData")

maps::map("world", xlim = c(-2,50), ylim = c(-30,50))
