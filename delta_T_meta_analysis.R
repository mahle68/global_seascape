# --------PILOT STUDY------- Oriental honey buzzard and Grey-faced Buzzard #### 

#Elham Nourani; Oct. 30. 2019; Radolfzell, Germany.

#description:
#Scripts for meta-analysis of relationship between delta_T and sea-crossing behavior in soaring birds


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

setwd("C:/Users/mahle/ownCloud/Work/Projects/delta_t")
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

#creating landmass layer
land <- shapefile("C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files/countries.shp")
land$landmass <- 1
landmasses <- unionSpatialPolygons(land,IDs = land$landmass)
landmasses$land <- 1
writeOGR(landmasses, "C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files", "world_landmass", driver = "ESRI Shapefile")

#creating the land_asia layer
land <- shapefile("C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files/world_landmass.shp") #open world landmass shapefile (codes for making the layer at the end of file)
land_asia <- crop(land,extent(90,143,-9,42.5))
land_asia <- as(land_asia,"SpatialPolygons")
save(land_asia,file = "R_files/land_east_asia.RData")

land_asia_sf <- st_as_sf(land_asia,crs = wgs) #convert to sf object
land_asia_sf_m <- st_transform(land_asia_sf,meters_proj)

med_sf1 <- st_buffer(land_asia_sf_m, dist = units::set_units(15000, 'm')) #put a buffer around the land polygon, to eliminate short sea-crossings. buffer is set to 15 km to match the 30 km res of env data
land_asia_15km <- st_transform(med_sf1,wgs)

save(land_asia_15km,file = "R_files/land_east_asia_15km.RData")

#geographic files
load("R_files/land_east_asia.RData") #called land_asia
load("R_files/land_east_asia_15km.RData") #called land_asia_15km



#### ---STEP 1: temporal filter for migratory seasons and filter for location classes (for PTTs) ####

OHB_files <- list.files("data/Oriental_honey_buzzard",pattern = ".xls",full.names = T)
OHB <- lapply(OHB_files,read_excel,1,col_types = c("numeric","date","numeric","numeric","numeric","skip","text",rep("numeric",8))) %>%
  reduce(full_join) %>%
  rename(date_time = 'date(gmt)',lon = longitud,lat = latitude) %>%
  mutate(yday = yday(date_time)) %>%
  mutate(season = ifelse(between(yday,253,294),"autumn",ifelse(month == 5,"spring","other")),
         track = paste(ptt,year,sep = "_"),
         species = "OHB") %>% #11 Sep-20 Oct; spring between 1-5 May
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
         species = "GFB") %>% #(Nourani et al 2017 for autumn; no sea-crossing in spring based on visually exploring the data.)
  filter(season %in% c("spring","autumn") &
           class %in% c("0","1","2","3")) #filter for location classes

#put the two dataframes together as a list
data <- list(GFB, as.data.frame(OHB))

data_sp <- lapply(data,function(x){
  coordinates(x) <- ~lon+lat
  proj4string(x) <- wgs
  x
})

save(data_sp,file = "R_files/GFB_HB_temp_filtered.RData")

#temporally fitlered data for GBF and OHB
load("R_files/GFB_HB_temp_filtered.RData") #called data_sp



#### ---STEP 2: spatial filter for points over the sea ##### 
sea_data <- lapply(data_sp,function(x){
  pts_over_land <- over(x,land_asia) #make sure the spatialpolygon layer is not a dataframe
  pts_over_sea <- x[is.na(pts_over_land),]
  pts_over_sea
})

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



#### ---STEP 3: temporal filter for autocorrelation #####

#check for duplicated time-stamps
getDuplicatedTimestamps(x = as.factor(sea_data_15[[1]]$track),timestamps = sea_data_15[[1]]$date_time,sensorType = "ptt")
getDuplicatedTimestamps(x = as.factor(sea_data_15[[2]]$track),timestamps = sea_data_15[[2]]$date_time,sensorType = "ptt") #two cases of duplicated timestamps. the points are very close to each other though.
mapview(sea_data_15[[2]][c(111:114),])
sea_data_15[[2]] <- sea_data_15[[2]][-c(17,112,273),] #remove the first records of duplicate timestamps; also remove the very unlikely record near india (row 273)


#create a move list
move_ls <- lapply(sea_data_15,function(x){
  mv <- move(x = st_coordinates(x)[,"X"],y = st_coordinates(x)[,"Y"],time = x$date_time,data = as.data.frame(x,row.names = c(1:nrow(x))),animal = x$track,proj = wgs)
  mv
})

#thin the data to have one-hour difference?
sea_data_15_1hr <- lapply(move_ls,function(group){ #for each of the groups (each species) 
  
  sp_obj_ls <- lapply(split(group),function(track){ #for each track within the group
    
    thinned_track <- track %>%
      thinTrackTime(interval = as.difftime(1, units = 'hours'), 
                    tolerance = as.difftime(15, units = 'mins'))
    
    #convert back to a move object (from move burst)
    thinned_track <- as(thinned_track,"Move")
    thinned_track
  })
  
  sp <- do.call(rbind,sp_obj_ls)
  sp$track <- track@idData$track #reassign the track
  st_as_sf(sp) #conver the spatialdataframe to an sf object
})

save(sea_data_15_1hr,file = "R_files/GFB_HB_sea_data_15_1hr.RData")

#temporal filter for autocorrelation (1-hourly thinned data)
load("R_files/GFB_HB_sea_data_15_1hr.RData") #called sea_data_15_1hr



#### ---STEP 4: produce alternative points in time #####
#.... within the migration period or outside? or both.... or, separately for before and after.

#open data
load("R_files/GFB_HB_temp_sp_filtered_15km_ann.RData") #called ann_df ...this is from STEP 4 of the previous version. just combine the two datasets without annotation

#for each point, create a alternative points a week before and a week after the observed point. year and hour dont change.
ann_df_alt <- ann_df %>%
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         obs_id = row_number()) %>%
  slice(rep(row_number(),61)) %>% #copy each row 60 times. 1 used, 60 alternative
  arrange(date_time) %>%
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 61 == 0, 1, 0))) #assign used and available values
  
  ann_alt_ls <- lapply( split(ann_df_alt,ann_df_alt$obs_id),function(x){ #didnt manage to write this part using dplyr and purrr
    alt_times <- alt_pts_temporal(x$date_time[1],60)
    x %>%
      mutate(timestamp = as.POSIXct(strptime(c(as.character(x$date_time[1]),alt_times$dt),format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
      period = c("now",alt_times$period))
  })

  ann_df_alt_cmpl <- do.call(rbind,ann_alt_ls)
  
save(ann_df_alt_cmpl,file = "R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_60days.RData")
  
#### ---STEP 5: annotate alternative points with delta T #####
load("R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_60days.RData") #called ann_df_alt_cmpl

#prep for track annotation on movebank
ann_df_alt_cmpl_mb <- ann_df_alt_cmpl %>%
  dplyr::select(-contains("ECMWF")) %>% #remove the already existing movebank columns
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) 
  
#rename columns
colnames(ann_df_alt_cmpl_mb)[c(8,9)] <- c("location-long","location-lat")

write.csv(ann_df_alt_cmpl_mb,"R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_60days.csv")
  
#downloaded from movebank
ann <- read.csv("movebank_annotation/GFB_HB_temp_sp_filtered_15km_ann_alt_60days.csv-5317423982503990720.csv", stringsAsFactors = F) %>%
  mutate(delta_T = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature - ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.) %>%
  drop_na() #remove NAs
  
save(ann, file = "R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_ann_60days.RData")


#### ---STEP 6: analysis #####
  
load("R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_ann_60days.RData") #called ann

#are used and available density plots different
plot(density(ann[ann$period == "now","delta_T"]),col = "red")
lines(density(ann[ann$period == "before","delta_T"]),col = "green")
lines(density(ann[ann$period == "after","delta_T"]),col = "blue")

legend("topleft",legend = c("used","available-before","available-after"), col = c("red","green","blue"),pch = 20, bty = "n", cex = 0.9)

#are used and available density plots different... season-specific
windows()
par(mfrow = c(1,2))
plot(density(ann[ann$period == "now" & ann$month %in% c(3,4),"delta_T"]),col = "red", main = "spring")
lines(density(ann[ann$period == "before" & ann$month %in% c(3,4),"delta_T"]),col = "green")
lines(density(ann[ann$period == "after" & ann$month %in% c(3,4),"delta_T"]),col = "blue")

legend("topleft",legend = c("used","available-before","available-after"), col = c("red","green","blue"),pch = 20, bty = "n", cex = 0.9)

plot(density(ann[ann$period == "now" & ann$month %in% c(9,10),"delta_T"]),col = "red", main = "autumn")
lines(density(ann[ann$period == "before" & ann$month %in% c(9,10),"delta_T"]),col = "green")
lines(density(ann[ann$period == "after" & ann$month %in% c(9,10),"delta_T"]),col = "blue")

#are used and available density plots different.. weekly plots
#assign a week variable
ann <- ann %>%
  mutate(timestamp,date_time = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  group_by(obs_id,period) %>%
  mutate(position = 1:n()) %>%
  mutate(week = ifelse( used == 1, "obs_day",
                        ifelse(position <= 8, "week_one",
                               ifelse(between(position, 9, 16), "week_two",
                                      ifelse(between(position, 17, 24), "week_three",
                                             "week_four")))),
         season = ifelse(month %in% c(3,4),"spring","autumn")) %>%
  ungroup() %>%
  as.data.frame()


save(ann, file = "R_files/GFB_HB_temp_sp_filtered_15km_ann_alt_ann_60days_weeks.RData")

#create weekly plots
par(mfrow = c(2,4))
ann %>%
  filter(used == 0) %>%
  group_by(season, period, week)
 
#library(sm)
#sm.density.compare(ann$delta_T,ann$used,xlab = "delta T")
#legend(locator(1),levels(ann$used))

#create training and testing set?
model <- glmer(used ~ delta_T + (1 | obs_id) + (1 | species), family = binomial, data = ann)

model <- glmer(used ~ delta_T + (1 | obs_id), family = binomial, data = ann)

model2 <- glm(used ~ scale(delta_T) , family = binomial, data = ann)

#use conditional logistic regressoin
library(survival)
form1a <- (used ~ scale(delta_T) + strata(obs_id))

#build the model using all the data
model <- clogit(form1a, data = ann)
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