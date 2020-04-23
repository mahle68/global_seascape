#script for preparing all input data simultaneously. previously done separately for each species/region
#Elham Nourani,
#Dec. 31. 2019. Radolfzell, Germany.
#update Jan 29, 2020: remove location class 0 from all ptt data
#update Apr 7, 2020: prep OHB GPS data (end of script)
#update Apr 23, 2020: prep Greek EF data (end of script)

library(tidyverse)
library(readxl) #read_excel()
library(lubridate)
library(move)
library(mapview)
library(rWind)
library(sf)
library(parallel)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")


load("R_files/land_15km.RData") #called land_15 km
load("R_files/land_0_60.RData") #called land_0_60; this doesnt have the 15 km buffer

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
PF <- read.csv("data/LifeTrack Peregrine falcon.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,38,39) %>% #remove columns that are not needed
  filter(individual.local.identifier %in% pf_ad$animal.id) %>% 
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         species = "PF",
         season = ifelse (month %in% c(9,10), "autumn", "other")) %>% 
  mutate(track = paste(tag.local.identifier, year, season,sep = "_")) %>% 
  filter(season != "other")

OE <- read.csv("data/Osprey in Mediterranean (Corsica, Italy, Balearics).csv", stringsAsFactors = F) %>% 
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



##### STEP 2: merge all data, assign zone (only northern hemisphere) #####

dataset <- list(OHB,GFB, PF, OE, OA) %>% 
  reduce(full_join, by = c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
  mutate(zone = ifelse(between(location.lat, 0, 30), "tradewind",
                       ifelse(between(location.lat, 30,60), "temperate",
                              "other"))) %>% 
  filter(zone != "other") %>% 
  dplyr::select(c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
  as.data.frame()

save(dataset, file = "R_files/all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData")

##### STEP 3: filter out points over land ##### 
#from here on is not updated with the new version of dataset (new track id) Feb.6

dataset_sea <- dataset %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_difference(land_0_60)

save(dataset_sea, file = "R_files/all_spp_spatial_filtered_updated_no_buffer.RData")

##### STEP 4: temporal filter for autocorrelation #####

load("R_files/all_spp_spatial_filtered_updated_no_buffer.RData") #named dataset_sea
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
