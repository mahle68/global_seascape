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

setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

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

#### ---STEP 1: temporal filter for migratory seasons and filter for location classes (for PTTs) ####

OA <- read.csv("data/Osprey_Americas/Osprey Bierregaard North and South America.csv", stringsAsFactors = F) %>% 
  dplyr::select(1,3:5,48:52) %>% #remove columns that are not needed
  mutate(date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  mutate(month = month(date_time),
         track = paste(tag.local.identifier, year(date_time),sep = "_"),
         species = "Osprey",
         zone = "tradewind") %>%  #set up a season variable
  filter(sensor.type == "gps") #keep only the gps tag data


OA_sp <- OA %>% 
  drop_na(c("location.long", "location.lat")) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs)


save(OA_sp,file = "R_files/OA_GPS_pre_filter.RData")

#filter temporally.... for migration season.....................................

#### ---STEP 3: temporal filter for autocorrelation #####

#check for duplicated time-stamps
getDuplicatedTimestamps(x = as.factor(OA_sp$track),timestamps = OA_sp$date_time,sensorType = "gps") #three cases of duplicated timestamps. the points are very close to each other though.


mapview(OA_sp[c(164559:164560),])
mapview(OA_sp[c(164561:164562),])
mapview(OA_sp[c(164563:164564),])

OA_sp<- OA_sp[-c(164559,164561,164563),] #remove the first records of duplicate timestamps;


#thin the data to have one-hour difference?
#create a move object
mv <- move(x = st_coordinates(OA_sp)[,"X"],y = st_coordinates(OA_sp)[,"Y"],time = OA_sp$date_time,data = as.data.frame(OA_sp,row.names = c(1:nrow(OA_sp))),animal = OA_sp$track,proj = wgs)
save(mv,file = "R_files/OA_GPS_pre_filter_mv.RData")
  
load("R_files/OA_GPS_pre_filter_mv.RData")
  
  #start the cluster
mycl <- makeCluster(detectCores() - 1)
  
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


###################################
sea_data_15_1hr[[1]]$species <- "GFB"
sea_data_15_1hr[[2]]$species <- "OHB"

sea_data_15_1hr <- sea_data_15_1hr[[1]] %>%
  full_join(sea_data_15_1hr[[2]], by = c("coords.x1","coords.x2","date_time","month","class","species","track")) %>% 
  dplyr::select(c(coords.x1,coords.x2,date_time,month,class,species,track)) %>% 
  rename(long = coords.x1, lat = coords.x2)

save(sea_data_15_1hr,file = "R_files/GFB_HB_sea_data_15_1hr.RData")



#### ---STEP 2: spatial filter for points over the sea ##### 
load("R_files/land_americas_15km.RData") #named land_am_15km
load("R_files/OA_GPS_pre_filter.RData")

sea_data <- st_difference(OA_sp,land_am_15km)
save(sea_data, file = "R_files/OA_GPS_sea.RData")

#sea_data <- lapply(OA_sp,function(x){
#  pts_over_land <- st_interover(x,land_am_15km) #make sure the spatialpolygon layer is not a dataframe
#  pts_over_sea <- x[is.na(pts_over_land),]
#  pts_over_sea
#})

#save the points
save(sea_data,file = "R_files/GFB_HB_temp_sp_filtered.RData")


#filter spatially with 15 km buffer
data_sf <- lapply(data_sp,st_as_sf,crs = wgs) #convert points to sf object

sea_data_15 <- lapply(data_sf,function(x){
  pts_over_water <- st_difference(x,land_asia_15km)
  pts_over_water
})

#save the points
save(sea_data_15,file = "R_files/GFB_HB_temp_sp_filtered_15km.RData")

#temporal and spatial filter: points over water
load("R_files/GFB_HB_temp_sp_filtered_15km.RData") #called sea_data_15






#temporally fitlered data for GBF and OHB
load("R_files/GFB_HB_temp_filtered.RData") #called data_sp



#### ---STEP 3: temporal filter for autocorrelation #####

#check for duplicated time-stamps
getDuplicatedTimestamps(x = as.factor(sea_data_15[[1]]$track),timestamps = sea_data_15[[1]]$date_time,sensorType = "ptt")
getDuplicatedTimestamps(x = as.factor(sea_data_15[[2]]$track),timestamps = sea_data_15[[2]]$date_time,sensorType = "ptt") #two cases of duplicated timestamps. the points are very close to each other though.
mapview(sea_data_15[[2]][c(111:114),])
sea_data_15[[2]] <- sea_data_15[[2]][-c(17,112,275),] #remove the first records of duplicate timestamps; also remove the very unlikely record near india (row 275)


#create a move list
move_ls <- lapply(sea_data_15,function(x){
  mv <- move(x = st_coordinates(x)[,"X"],y = st_coordinates(x)[,"Y"],time = x$date_time,data = as.data.frame(x,row.names = c(1:nrow(x))),animal = x$track,proj = wgs)
  mv
})

#thin the data to have one-hour difference?
sea_data_15_1hr <- lapply(move_ls,function(group){ #for each of the groups (each species) 
  
  sp_obj_ls <- lapply(split(group),function(one_track){ #for each track within the group
    
    thinned_track <- one_track %>%
      thinTrackTime(interval = as.difftime(1, units = 'hours'), 
                    tolerance = as.difftime(15, units = 'mins'))
    
    #convert back to a move object (from move burst)
    thinned_track <- as(thinned_track,"Move")
    thinned_track$track <- one_track@idData$track #reassign the track
    thinned_track
  })
  
  sp <- do.call(rbind,sp_obj_ls)
  #st_as_sf(sp) #conver the spatialdataframe to an sf object
  as.data.frame(sp)
}) 

sea_data_15_1hr[[1]]$species <- "GFB"
sea_data_15_1hr[[2]]$species <- "OHB"

sea_data_15_1hr <- sea_data_15_1hr[[1]] %>%
  full_join(sea_data_15_1hr[[2]], by = c("coords.x1","coords.x2","date_time","month","class","species","track")) %>% 
  dplyr::select(c(coords.x1,coords.x2,date_time,month,class,species,track)) %>% 
  rename(long = coords.x1, lat = coords.x2)

save(sea_data_15_1hr,file = "R_files/GFB_HB_sea_data_15_1hr.RData")


#### ---STEP 4: produce alternative points in time #####

#open data
load("R_files/GFB_HB_sea_data_15_1hr.RData") #called sea_data_15_1hr
#load("R_files/GFB_HB_temp_sp_filtered_15km_ann.RData") #called ann_df ...this is from STEP 4 of the previous version. just combine the two datasets without annotation

#for each point, create a alternative points a week before and a week after the observed point. year and hour dont change.
alts <- sea_data_15_1hr %>%
  mutate(#date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         obs_id = row_number()) %>%
  slice(rep(row_number(),15)) %>% #copy each row 60 times. 1 used, 60 alternative
  arrange(date_time) %>%
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 15 == 0, 1, 0))) #assign used and available values
  

  alt_ls <- lapply( split(alts,alts$obs_id),function(x){ #didnt manage to write this part using dplyr and purrr
    alt_times <- alt_pts_temporal(x$date_time[1],14)
    x %>%
      mutate(timestamp = as.POSIXct(strptime(c(as.character(x$date_time[1]),alt_times$dt),format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
      period = c("now",alt_times$period))
  })

  alt_cmpl <- do.call(rbind,alt_ls)
  
save(alt_cmpl,file = "R_files/GFB_HB_temp_sp_filtered_15km_alt_14days.RData")
  
#### ---STEP 5: annotate alternative points with delta T #####

load("R_files/GFB_HB_temp_sp_filtered_15km_alt_14days.RData") #called ann_df_alt_cmpl

#prep for track annotation on movebank
alt_cmpl_mb <- alt_cmpl %>%
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) 
  
#rename columns
colnames(alt_cmpl_mb)[c(1,2)] <- c("location-long","location-lat")

write.csv(alt_cmpl_mb,"R_files/GFB_HB_temp_sp_filtered_15km_alt_14days.csv")
  
#downloaded from movebank
ann <- read.csv("movebank_annotation/GFB_HB_temp_sp_filtered_15km_alt_14days.csv-799772115122349617.csv", stringsAsFactors = F) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m) %>%
  drop_na() #remove NAs
  
save(ann, file = "R_files/GFB_HB_temp_sp_filtered_15km_alt_ann_14days.RData")


#### ---STEP 6: visualizations #####
  
load("R_files/GFB_HB_temp_sp_filtered_15km_alt_ann_14days.RData") #called ann

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

  
#### ---STEP 7: analysis #####

#standardize the variables to have comparable effect sizes
ann_sc <- ann %>%
  mutate_at(c("delta_t","u925","v925"), scale) %>% 
  as.data.frame

#one glm for both seasons together
model <- glm(used ~ delta_t + u925 + v925 , family = binomial, data = ann_sc ) #higher selection on u-wind

model <- glmer(used ~ delta_t + u925 + v925 + (1 | obs_id) + (1 | species), family = binomial, data = ann ) #singularity error

#clogit for both seasons together
model_cl <- clogit(used ~ delta_t + u925 + v925 + strata(obs_id), data = ann_sc)
model_cl

#ann_sc <- ann_sc[order(ann_sc,ann_sc$obs_id),]
        model_cl_b <- stan_clogit(used ~ delta_t + u925 + v925, strata = obs_id, data = ann_sc)


#### ---STEP 8: modeling #####
#model only for autumn
autumn_one_week_before <- ann %>%
  filter(period %in% c("now","before") & week %in% c("obs_day", "week_one") & season == "autumn" ) %>%
  #select(-delta_t) %>%
  map_if(is.factor, as.character) %>%
  as.data.frame()
  

model <- glmer(used ~ delta_T + (1 | obs_id) + (1 | species), family = binomial, data = autumn_one_week_before)
  
model <- glm(used ~ delta_T  , family = binomial, data = autumn_one_week_before)

#autumn four weeks before
autumn_four_week_before <- ann %>%
  filter(period %in% c("now","before") & season == "autumn" ) %>%
  #select(-delta_t) %>%
  map_if(is.factor, as.character) %>%
  as.data.frame()


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
#### end ####
