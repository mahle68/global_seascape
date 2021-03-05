# --------EMPIRICAL STUDY------- ospreys in the american flyway #### 

#Elham Nourani; Dec. 11. 2019; Radolfzell, Germany.

#description:
#preparing the osprey dataset for meta-analysis of relationship between delta_T and wind and sea-crossing behavior in soaring birds

#### ---Libraries, ftns, and misc ####
library(dplyr)
library(purrr)
library(lubridate)
library(readxl)
library(raster)
library(maptools)
library(rgdal)
library(sf)
library(mapview)
library(move)
library(readr) #read_csv()
library(tidyverse) #nest()
library(lme4)
library(survival)
library(rstanarm)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")
mycl <- makeCluster(detectCores() - 1)


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


#### ---PREP geographic layers ####

#creating the land_america layer
land_am_15km <- st_read("/home/mahle68/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp") %>% 
  st_crop(y = c(xmin = -132, xmax = -19.5, ymin = -17, ymax = 53)) %>% 
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(15000, 'm')) %>% 
  st_transform(wgs) %>% 
  st_union()

save(land_am_15km,file = "R_files/land_americas_15km.RData")

#### ---STEP 1: filter for adults #####
#adults were identified based on the individual bios on http://www.ospreytrax.com/html_files/Shanawdithit.html#


OA <- read.csv("data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,48:52) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         track = paste(tag.local.identifier, year(date_time),sep = "_"),
         species = "Osprey",
         zone = "tradewind") %>%  #set up a season variable
  filter(sensor.type == "gps",#keep only the gps tag data
         individual.local.identifier %in% c("Holly","Hackett","Daphne","Shanawdithit","Gundersen","Wausau",
                                            "Crabby","Charlie","Roger Tory")) #keep only adults (based on what i found on the website)


#### ---STEP 2: spatial filter for points over the sea ####

OA_sp <- OA %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)


save(OA_sp,file = "R_files/OA_GPS_pre_filter.RData")

load("R_files/land_americas_15km.RData") #named land_am_15km
load("R_files/OA_GPS_pre_filter.RData")

sea_data <- st_difference(OA_sp,land_am_15km)
save(sea_data, file = "R_files/OA_GPS_sea.RData")

#### ---STEP 3: temporal filter for migration season #####
#assume migration is: autumn: Sep-Oct; spring: March-April...based on amount of data over water
sea_data_migr <-sea_data %>% 
  mutate(season = ifelse(month %in% c(3,4),"spring",
                         ifelse (month %in% c(9,10), "autumn",
                                 "other"))) %>% 
  filter(season != "other")

save(sea_data_migr, file = "R_files/OA_GPS_sea_migration.RData")

#### ---STEP 4: temporal filter for autocorrelation #####

load("R_files/OA_GPS_sea_migration.RData") #called sea_data_migr

#check for duplicated time-stamps
getDuplicatedTimestamps(x = as.factor(sea_data$track),timestamps = sea_data$date_time,sensorType = "gps") #three cases of duplicated timestamps. the points are very close to each other though.

#thin the data to have one-hour difference?
#create a move object
mv <- move(x = st_coordinates(sea_data)[,"X"],y = st_coordinates(sea_data)[,"Y"],time = sea_data$date_time,
           data = as.data.frame(sea_data,row.names = c(1:nrow(sea_data))),animal = sea_data$track,proj = wgs)

#start the cluster
clusterExport(mycl, "mv") #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(move)
  library(lubridate)
  library(dplyr)
  library(raster)
})

sp_obj_ls <- parLapply(cl = mycl,split(mv),function(one_track){ #for each track within the group
  
  thinned_track <- one_track %>%
    thinTrackTime(interval = as.difftime(1, units = 'hours'), 
                  tolerance = as.difftime(15, units = 'mins'))
  
  #convert back to a move object (from move burst)
  thinned_track <- as(thinned_track,"Move")
  thinned_track$track <- one_track@idData$track #reassign the track
  thinned_track
})

stopCluster(mycl)

sp <- do.call(rbind,sp_obj_ls)
sf <- st_as_sf(sp)



#### ---STEP 4: separate temperate and tradewind zone data #####

sf <-sf %>% 
  mutate(zone = ifelse(st_coordinates(.)[,"Y"] <= 30, "tradewind",
                       ifelse(between(st_coordinates(.)[,"Y"], 30,60), "temperate",
                              "other")))

save(sf,file = "R_files/AO_sea_mgr_15_1hr.RData")

#### ---STEP 5: produce alternative points in time #####

#open data
load("R_files/AO_sea_mgr_15_1hr.RData") #called sf. with original times

#timestamps at midnight cause a problem becuase hour becomes NA. add 5 minutes to all timestamps at midnight
sf <- sf %>% 
  mutate(date_time = ifelse(hour(date_time) == 0 & minute(date_time) == 0, as.character(date_time + minutes(5)), as.character(date_time))) %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"))


#for each point, create a alternative points a week before and a week after the observed point. year and hour dont change.
alts <- sf %>%
  mutate(#date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
    obs_id = row_number()) %>%
  slice(rep(row_number(),15)) %>% #copy each row 60 times. 1 used, 60 alternative
  arrange(date_time) %>%
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 15 == 0, 1, 0)),  #assign used and available values
         lon = st_coordinates(.)[,"X"],
         lat = st_coordinates(.)[,"Y"]) %>% 
  st_drop_geometry()%>% 
  filter(individual.local.identifier %in% c("Holly","Hackett","Daphne","Shanawdithit","Gundersen","Wausau",
                                            "Crabby","Charlie","Roger Tory"))


alt_ls <- lapply( split(alts,alts$obs_id),function(x){ #didnt manage to write this part using dplyr and purrr
  alt_times <- alt_pts_temporal(x$date_time[1],14)
  x %>%
    mutate(timestamp = as.POSIXct(strptime(c(as.character(x$date_time[1]),alt_times$dt),format = "%Y-%m-%d %H:%M:%S"),tz = "UTC",),
           period = c("now",alt_times$period))
})

alt_cmpl <- do.call(rbind,alt_ls)


save(alt_cmpl,file = "R_files/AO_temp_sp_filtered_15km_alt_14days.RData") #with 5 minutes added to 00:00

#### ---STEP 6: annotate alternative points with delta T #####

load("R_files/AO_temp_sp_filtered_15km_alt_14days.RData") #called alt_cmpl

#prep for track annotation on movebank
alt_cmpl_mb <- alt_cmpl %>%
  mutate(timestamp = paste(as.character(timestamp),"000",sep = "."))

#rename columns
colnames(alt_cmpl_mb)[c(11,12)] <- c("location-long","location-lat")

write.csv(alt_cmpl_mb,"R_files/AO_temp_sp_filtered_15km_alt_14days.csv") #with 5 minutes added to 00:00

#downloaded from movebank
ann <- read.csv("movebank_annotation/AO_temp_sp_filtered_15km_alt_14days.csv-2568627422738596992.csv", stringsAsFactors = F) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m) %>%
  drop_na() #remove NAs

save(ann, file = "R_files/AO_temp_sp_filtered_15km_alt_ann_14days.RData")


#### ---STEP 7: visualizations #####

load("R_files/AO_temp_sp_filtered_15km_alt_ann_14days.RData") #called ann

spring <- ann %>% 
  filter(month %in% c(3,4))

autumn <- ann %>% 
  filter(month %in% c(9,10))

for(i in list(spring,autumn)){
  
  #are used and available density plots different for each season
  X11()
  par(mfrow = c(1,3))
  plot(density(ann[ann$period == "now","delta_t"]),col = "red", main = "delta t")
  lines(density(ann[ann$period != "now","delta_t"]),col = "green")
  legend("topleft",legend = c("used","available"), col = c("red","green"),lty = 1, bty = "n", cex = 0.9)
  
  plot(density(ann[ann$period == "now","u925"]),col = "red", main = "u-wind")
  lines(density(ann[ann$period != "now","u925"]),col = "green")
  
  plot(density(ann[ann$period == "now","v925"]),col = "red", main = "v-wind")
  lines(density(ann[ann$period != "now","v925"]),col = "green")
  
} #not much seasonal difference. there seems to be selection for u-wind


#### ---STEP 8: analysis #####
#---------------------------------------make sure to model only the tradewind zone. and perhaps seasonal
#standardize the variables to have comparable effect sizes
ann_sc <- ann %>%
  mutate_at(c("delta_t","u925","v925"), scale) %>% 
  as.data.frame

#one glm for both seasons together... for zone-specific models, make sure to scale the data only based on values in that zone :p
model <- glm(used ~ delta_t + u925 + v925 , family = binomial, data = ann_sc[ann_sc$zone == "tradewind",] ) #higher selection on u-wind

model <- glmer(used ~ delta_t + u925 + v925 + (1 | obs_id), family = binomial, data = ann_sc[ann_sc$zone == "tradewind",]  ) #singularity error

#clogit for both seasons together
model_cl <- clogit(used ~ delta_t + u925 + v925 + strata(obs_id), data = ann_sc[ann_sc$zone == "tradewind",] )
model_cl

model_cl2 <- clogit(used ~ delta_t + u925 + v925 + strata(obs_id), data = ann_sc[ann_sc$zone == "temperate",] )
model_cl2

#ann_sc <- ann_sc[order(ann_sc,ann_sc$obs_id),]
model_cl_b <- stan_clogit(used ~ delta_t + u925 + v925, strata = obs_id, data = ann_sc)


#### ---STEP 9: modeling #####
#model only for autumn
#autumn_one_week_before <- ann %>%
#  filter(period %in% c("now","before") & week %in% c("obs_day", "week_one") & season == "autumn" ) %>%
#select(-delta_t) %>%
#  map_if(is.factor, as.character) %>%
#  as.data.frame()


#model <- glmer(used ~ delta_T + (1 | obs_id) + (1 | species), family = binomial, data = autumn_one_week_before)

#model <- glm(used ~ delta_T  , family = binomial, data = autumn_one_week_before)

#autumn four weeks before
##autumn_four_week_before <- ann %>%
#  filter(period %in% c("now","before") & season == "autumn" ) %>%
#select(-delta_t) %>%
#  map_if(is.factor, as.character) %>%
#  as.data.frame()


model <- glmer(used ~ delta_T + (1 | obs_id) + (1 | species), family = binomial, data = autumn_four_week_before)

model <- glm(used ~ delta_T  , family = binomial, data = autumn_one_week_before)

#create training and testing set?
model <- glmer(used ~ delta_t + (1 | obs_id) + (1 | species), family = binomial, data = ann)

model <- glmer(used ~ delta_t + (1 | obs_id), family = binomial, data = ann)

model2 <- glm(used ~ scale(delta_t) , family = binomial, data = ann)

#use conditional logistic regressoin
library(survival)
form1a <- (used ~ scale(delta_T) + strata(obs_id))

#build the model using all the data
model <- clogit(form1a, data = autumn_one_week_before)
model


#evaluate the model



####---PLOTTINTG #####
windows()
plot(land_asia)
points(data_sp[[1]],col = data_sp[[1]]$month,pch = 16,cex = 0.5)
points(GFB$lon,GFB$lat,col = as.factor(GFB$season))
points(GFB[GFB$season == "spring","lon"],GFB[GFB$season == "spring","lat"],pch = 16,cex = 0.4,col = "green")
points(GFB[GFB$season == "autumn","lon"],GFB[GFB$season == "autumn","lat"],pch = 16,cex = 0.4,col = "blue")
points(GFB[GFB$month == 10,"lon"],GFB[GFB$month == 10,"lat"],pch = 16,cex = 0.4,col = "blue")
points(GFB[GFB$month == 9,"lon"],GFB[GFB$month == 9,"lat"],pch = 16,cex = 0.4,col = "orange")
points(GFB[GFB$month == 3,"lon"],GFB[GFB$month == 3,"lat"],pch = 16,cex = 0.4,col = "pink")
points(GFB[GFB$month == 4,"lon"],GFB[GFB$month == 4,"lat"],pch = 16,cex = 0.4,col = "purple")

lapply(data_sp,points,pch = 16,cex = 0.3,col = "blue")
lapply(sea_data,points,pch = 16,cex = 0.5,col = "orange")

library(mapview)
mapview(x, col.regions = sf.colors(10))
mapview(pts_over_water, col.regions = "red")

points(ann[is.na(ann$delta_T),c("location.long","location.lat")])


lapply(sea_data_15,mapview)


plot(land_am_15km)
points(sp,cex = 0.3, col = factor(sp$month))
points(alt_cmpl[alt_cmpl$used == 1,c("lon","lat")], pch = 16, cex= 0.3, col = factor(alt_cmpl[alt_cmpl$used == 1,"zone"]))
#### end ####
