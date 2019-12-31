#script for preparing all input data simultaneously. previously done separately for each species/region
#Elham Nourani,
#Dec. 31. 2019. Radolfzell, Germany.

library(readxl) #read_excel()
library(lubridate)
library(move)
library(mapview)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")
mycl <- makeCluster(detectCores() - 2)

load("R_files/land_15km.RData") #called land_15 km

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

#also assign date_time, year, month, track, species, location.lat, location.long

OHB_files <- list.files("data/Oriental_honey_buzzard",pattern = ".xls",full.names = T)
OHB <- lapply(OHB_files,read_excel,1,col_types = c("numeric","date","numeric","numeric","numeric","skip","text",rep("numeric",8))) %>%
  reduce(full_join) %>%
  rename(date_time = 'date(gmt)',lon = longitud,lat = latitude) %>%
  mutate(yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other")), #11 Sep-20 Oct; spring between 1-5 May
         track = paste(ptt,year,sep = "_"),
         species = "OHB") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season == "autumn" & #no sea-crossing in spring
           class %in% c("0","1","2","3")) #filter for location classes


GFB_files <- list.files("data/Grey_faced_buzzard/",pattern = ".csv",full.names = T)
GFB <- lapply(GFB_files,read.csv,stringsAsFactors = F) %>%
  reduce(full_join) %>% #is locdate is in UTC
  mutate(locdate,date_time = as.POSIXct(strptime(locdate,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  mutate(month = month(date_time),
         year = year(date_time),
         season = ifelse(month %in% c(3,4),"spring",ifelse(month %in% c(10),"autumn","other")),
         track = paste(platform,year,sep = "_"),
         species = "GFB") %>% 
  rename(location.long = lon,
         location.lat = lat) %>% 
  filter(season %in% c("spring","autumn") &
           class %in% c("0","1","2","3")) #filter for location classes

###pf have unknown age! ask Ivan.... also make sure migration season is correctly defined
PF <- read.csv("data/LifeTrack Peregrine falcon.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,16,38,39) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "PF",
         season = ifelse (month %in% c(9,10), "autumn", "other")) %>% 
  filter(season != "other")

OE <- read.csv("data/Osprey in Mediterranean (Corsica, Italy, Balearics).csv", stringsAsFactors = F) %>% 
  filter(grepl("ad",individual.local.identifier,, ignore.case = T) & !grepl("juv",individual.local.identifier,, ignore.case = T)) %>% #extract adult data
  dplyr::select(1,3:5,16,35:37) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "O",
         season = ifelse(month %in% c(2:4),"spring",ifelse (month %in% c(8:10), "autumn","other"))) %>% 
  filter(season != "other")

OA <- read.csv("data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,48:52) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         year = year(date_time),
         track = paste(tag.local.identifier, year,sep = "_"),
         species = "O",
         season = ifelse(month %in% c(3,4),"spring",ifelse (month %in% c(9,10), "autumn", "other"))) %>% 
  filter(season != "other") %>% 
  filter(sensor.type == "gps",#keep only the gps tag data
         individual.local.identifier %in% c("Holly","Hackett","Daphne","Shanawdithit","Gundersen","Wausau",
                                            "Crabby","Charlie","Roger Tory")) #keep only adults (based on what i found on the website)



##### STEP 2: merge all data, assign zone (only northern hemisphere) #####

dataset <- list(OHB,GFB, PF, OE, OA) %>% 
  reduce(full_join, by = c("location.long", "location.lat", "date_time", "track", "month", "year" , "season", "species")) %>% 
    mutate(zone = ifelse(between(location.lat, 0, 30), "tradewind",
                       ifelse(between(location.lat, 30,60), "temperate",
                              "other"))) %>% 
  filter(zone != "other")

##### STEP 3: filter out points over land #####

dataset_sea <- dataset %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  st_difference(land_15km)

save(dataset_sea, file = "R_files/all_spp_spatial_filtered.RData")

##### STEP 4: temporal filter for autocorrelation #####

#check for duplicated time-stamps
rows_to_delete <- sapply(getDuplicatedTimestamps(x = as.factor(dataset_sea$track),timestamps = dataset_sea$date_time,sensorType = "gps"),"[",2) #get the second row of each pair of duplicate rows
dataset_sea <- dataset_sea[-rows_to_delete,]

#thin the data to have one-hour difference?
#create a move object
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

##### STEP 6: annotate alternative points #####

load("R_files/all_spp_temp_sp_filtered_15km_alt_14days.RData") #called alt_cmpl

#prep for track annotation on movebank
alt_cmpl_mb <- alt_cmpl %>%
  mutate(timestamp = paste(as.character(timestamp),"000",sep = "."))

#rename columns
colnames(alt_cmpl_mb)[c(10,11)] <- c("location-long","location-lat")

write.csv(alt_cmpl_mb,"R_files/all_spp_temp_sp_filtered_15km_alt_14days.csv") #with 5 minutes added to 00:00

#downloaded from movebank
ann <- read.csv("movebank_annotation/PF_temp_sp_filtered_15km_alt_14days.csv-5813773760251615680.csv", stringsAsFactors = F) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m) %>%
  drop_na() #remove NAs...for european osprey, the NA values are for the 2019 track. with a transition to ERA5, I should be able to use this data as well

save(ann, file = "R_files/PF_temp_sp_filtered_15km_alt_ann_14days.RData")
