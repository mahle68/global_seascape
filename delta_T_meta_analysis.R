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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
alt_pts_week <- function(date_time) {
  #same year, same hour, only day changes
  alt_pts_before <- vector()
  alt_pts_after <- vector()
  
  for (i in 1:7) {
    alt_pts_before[i] <- as.character(date_time - days(i))  
    alt_pts_after[i] <- as.character(date_time + days(i)) 
  }
  c(alt_pts_before,alt_pts_after)
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



#### ---STEP 4: annotae with delta T #####
#change column names to match movebank generic file standards and save to file

sea_data_15_1hr_movebank <- lapply(sea_data_15_1hr,function(x){
  x_df <- data.frame(x)
  x_df <- x_df[,-which(names(x_df) == "geometry")]
  #add miliseconds
  x_df$timestamp <- paste(as.character(x_df$date_time),"000",sep = ".")
  #coordinate columns
  x_df$X <- st_coordinates(x)[,"X"]
  x_df$Y <- st_coordinates(x)[,"Y"]
  x_df
})

colnames(sea_data_15_1hr_movebank[[1]])[c(9,10)] <- c("location-long","location-lat")
colnames(sea_data_15_1hr_movebank[[2]])[c(14,15)] <- c("location-long","location-lat")

#files to submit to movebank
write.csv(sea_data_15_1hr_movebank[[1]],"R_files/movebank_sea_data_15_1hr_GFB.csv")
write.csv(sea_data_15_1hr_movebank[[2]],"R_files/movebank_sea_data_15_1hr_OHB.csv")

#open annotated file
file_ls <- list.files("movebank_annotation",pattern = ".csv",full.names = T)

ann_df <- file_ls %>%
  map(read.csv, stringsAsFactors = F) %>%
  map(dplyr::select, com_col) %>% #define com_col as a vector of common columns between the two files: com_col <- intersect(colnames(ann_ls[[1]]),colnames(ann_ls[[2]]))
  reduce(rbind) %>%
  mutate(delta_t = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature - ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.)

save(ann_df,file = "R_files/GFB_HB_temp_sp_filtered_15km_ann.RData") #save



#### ---STEP 5: explore delta T values ####
lapply(ann_ls,summary)


#### ---STEP 6: produce alternative points in time #####
#.... within the migration period or outside? or both.... or, separately for before and after.

#open data
load("R_files/GFB_HB_temp_sp_filtered_15km_ann.RData") #called ann_df 

#for each point, create a alternative points a week before and a week after the observed point. year and hour dont change.
ann_df_alt <- ann_df %>%
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC"),
         obs_id = row_number()) %>%
  slice(rep(row_number(),15)) %>% #copy each row 15 times. 1 used, 14 alternative
  arrange(date_time) %>%
  mutate(used = ifelse(row_number() == 1,1,
                       ifelse((row_number() - 1) %% 15 == 0, 1, 0))) #assign used and available values
  
  ann_alt_ls <- lapply( split(ann_df_alt,ann_df_alt$obs_id),function(x){ #didnt manage to write this part using dplyr and purrr
    alt_times <- alt_pts_week(x$date_time[1])
    x$date_time[-1 ] <- as.POSIXct(strptime(alt_times,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")
    x
  })

  ann_df_alt_cmpl <- do.call(rbind,ann_alt_ls)
save(ann_df_alt_cmpl,file = "R_files/GFB_HB_temp_sp_filtered_15km_ann_alt.RData")
  
#prep for track annotation on movebank
ann_df_alt_cmpl_mb <- ann_df_alt_cmpl %>%
  dplyr::select(-contains("ECMWF")) %>% #remove the already existing movebank columns
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) 
  
#rename columns
colnames(ann_df_alt_cmpl_mb)[c(8,9)] <- c("location-long","location-lat")

write.csv(ann_df_alt_cmpl_mb,"R_files/GFB_HB_temp_sp_filtered_15km_ann_alt.csv")

  

#### ---STEP 7: annotate alternative points with delta T #####
  
  
#### ---STEP 8: analysis #####
  
load("R_files/GFB_HB_temp_sp_filtered_15km_ann_alt.RData") #called ann_df_alt_cmpl

model <- glmer(used ~ delta_t + (1 | obs_id), family = binomial, data = ann_df_alt_cmpl)
  
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

lapply(sea_data_15,mapview)
#### end ####