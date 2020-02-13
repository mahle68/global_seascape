#script for masking the tracks using the land layer on the cluster
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R
#Elham Nourani. Feb. 6. 2020. Radolfzell, Germany


args <- (commandArgs(trailingOnly = TRUE)) # provides access to a copy of the command line supplied when this R session was invoked
eval(parse(text = args)) #eval evalueates an R expression in a specified environment
n <- as.numeric(as.character(Line)) # the object “Line” comes from the .slrm file. This is the index

#install packages that I couldnt install using conda install


#open libraries
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(parallel)
library(mapview)
library(lutz)
library(RNCEP)

#options(digits = 10) #this only affects what is printed. not what is in the data. use round

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

source("wind_support_Kami.R")
NCEP.loxodrome.mod <- function (lat1, lat2, lon1, lon2) {
  deg2rad <- pi/180
  acot <- function(x) {
    return(atan(1/x))
  }
  lat1 <- deg2rad * lat1
  lat2 <- deg2rad * lat2
  lon1 <- deg2rad * lon1
  lon2 <- deg2rad * lon2
  deltaLon <- lon2 - lon1
  pi4 <- pi/4
  Sig1 <- log(tan(pi4 + lat1/2))
  Sig2 <- log(tan(pi4 + lat2/2))
  deltaSig <- Sig2 - Sig1
  if (deltaLon == 0 && deltaSig > 0) {
    head <- 0
  }
  else if (deltaLon == 0 && deltaSig < 0) {
    head <- 180
  }
  else if (deltaSig == 0 && deltaLon > 0) {
    head <- 90
  }
  else if (deltaSig == 0 && deltaLon < 0) {
    head <- 270
  }
  else if (deltaSig < 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig < 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 180
  }
  else if (deltaSig > 0 && deltaLon > 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi
  }
  else if (deltaSig > 0 && deltaLon < 0) {
    head <- acot(deltaSig/deltaLon) * 180/pi + 360
  }
  else head <- NA
  return(head)
}


#load files. make sure they are all stored in the working directory
load("twz_sf.RData") #twz_sf
load("tmpz_sf.RData") #tmpz_sf
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset
load("land_0_60.RData") #land_0_60
load("ocean_0_60.RData") #ocean

##### STEP 1: convert tracks to spatial lines and remove portions over land #####
ocean_sp <- as(ocean,"Spatial")
land_sp <- as(land_0_60,"Spatial")
land_b<-buffer(land_sp,width=0.001)

coordinates(dataset) <- ~location.long + location.lat
proj4string(dataset) <- wgs

track_ls <- split(dataset,dataset$track)
track_ls <- track_ls[lapply(track_ls,nrow)>1] #remove tracks with one point

b <- Sys.time()
mycl <- makeCluster(9) #total number of tracks is 369, so 41 will be sent to each core

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

Sys.time() - b #takes 45 min

save(Lines_ls,file = "Lines_ls_no_land.RData") 

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

save(segs_filtered,file = "Segs_no_land_filtered.RData") 

#assign zone to each segment
segs_filtered$twz <- as.numeric(st_within(segs_filtered,twz_sf))
segs_filtered$tmpz <- as.numeric(st_within(segs_filtered,tmpz_sf))

segs_filtered <- segs_filtered %>% 
  mutate(zone = ifelse(is.na(twz) == TRUE & is.na(tmpz) ==T, "both",
                       ifelse(is.na(twz) == TRUE & tmpz == 1, "tmpz",
                              "twz"))) %>% 
  dplyr::select(-c("twz","tmpz"))

save(segs_filtered, file = "filtered_segs.RData")

##### STEP 3: annotate with date-time #####

load("filtered_segs.RData") #segs_filtered. make sure track is character and not factor. so that there are no empty tracks

#convert segments to points
segs_pts <- segs_filtered %>% 
  mutate(seg_id = seq(1:nrow(.))) %>% 
  st_cast("POINT")

#create a buffer around the dataset points to make polygons. then overlay
dataset_buff <- dataset %>% 
  st_as_sf(coords = c("location.long","location.lat"), crs = wgs) %>% 
  st_transform(meters_proj) %>% 
  st_buffer(dist = units::set_units(10, 'm')) %>% 
  st_transform(wgs) 

save(dataset_buff,file = "dataset_10m_buffer.RData")
load("dataset_10m_buffer.RData")

#for each segs_pts point, find the index of the dataset_buff polygon that it intersects, then extract that row from dataset and add to segs_pts
mycl <- makeCluster(9) 

clusterExport(mycl, c("segs_pts", "dataset_buff")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(raster)
  library(tidyverse)
})

b <- Sys.time()

segs_ann <- parLapply(mycl,split(segs_pts,segs_pts$track), function(x){ #separate by track first to break up the job into smaller chunks
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

save(segs_ann, file = "segs_dt.RData")

##### STEP 4: create alternative tracks in time #####

load("segs_dt.RData") #segs_ann
segs_df <- segs_ann %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) 

#add daily alternatives for two weeks before and two weeks after the point
days_to_add <- c(0,cumsum(rep(1,14)),cumsum(rep(-1,14)))

pts_alt <- segs_df %>% 
  mutate(obs_id = row_number()) %>% 
  slice(rep(row_number(),29)) %>%  #paste each row 29 time for 29 days
  #mutate(used = ifelse(row_number() == 1,1,
  #                     ifelse((row_number() - 1) %% 29 == 0, 1, 0))) %>% 
  arrange(obs_id) %>% 
  group_by(obs_id) %>% 
  #arrange(obs_id) %>% 
  mutate(days_to_add = days_to_add) %>% 
  mutate(alt_date_time = date_time + days(days_to_add)) %>%  #use days to add as an id for alternative segments
  ungroup()

  save(pts_alt, file = "alt_pts_alt_time.RData")


#adding hourly alternatives takes too long. try daily.
# hours_to_add <- c(0,cumsum(rep(1,672/2)),cumsum(rep(-1,(672/2))))
# 
# pts_alt <- segs_ann %>% 
#   mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% #,
#          #tz = tz_lookup_coords(st_coordinates(.)[,2],st_coordinates(.)[,1])) #no need for local time, because i decided not to limit it to daytime.
#   as.data.frame() %>% 
#   mutate(obs_id = row_number()) %>% 
#   slice(rep(row_number(),673)) %>%  #paste each row 695 times for alternative points: 28days *24 hours + 23 hours in observed day
#   mutate(used = ifelse(row_number() == 1,1,
#                        ifelse((row_number() - 1) %% 673 == 0, 1, 0))) %>% 
#   group_by(obs_id) %>% 
#   arrange(obs_id) %>% 
#   mutate(hours_to_add = hours_to_add) %>% 
#   mutate(alt_date_time = date_time + hours(hours_to_add)) %>%  #use hours to add as an id for alternative segments
#   ungroup()
# 
# save(pts_alt, file = "alt_pts_alt_time.RData")


##### STEP 5: annotate all points #####

load("alt_pts_alt_time.RData") #called pts_alt

#prep for track annotation on movebank
pts_alt_mb <- pts_alt %>%
  mutate(timestamp = paste(as.character(alt_date_time),"000",sep = ".")) %>% 
  as.data.frame()

#rename columns
colnames(pts_alt_mb)[c(11,12)] <- c("location-long","location-lat")

write.csv(pts_alt_mb,"alt_pts_mb.csv") 

# this is over 9 million rows. break up into 10 files to upload to movebank
pts_alt_mb$chuncks <-c(rep(1,1e6),rep(2,1e6),rep(3,1e6),rep(4,1e6),rep(5,1e6),rep(6,1e6),rep(7,1e6),
                       rep(8,1e6),rep(9,1e6),rep(10,nrow(pts_alt_mb)-9e6))

lapply(split(pts_alt_mb,pts_alt_mb$chuncks),function(x){
  write.csv(x,paste("alt_pts_mb_chunk_",x$chuncks[1],".csv",sep = ""))
})

#downloaded from movebank
file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/alt_segs/", full.names = T, pattern = ".csv$") # $ means end of string
pts_ann <- lapply(file_ls,read.csv, stringsAsFactors = F) %>% 
  reduce(rbind) %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u950 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v950 = ECMWF.Interim.Full.Daily.PL.V.Wind,
         u10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.U.Component.,
         v10m = ECMWF.Interim.Full.Daily.SFC.Wind..10.m.above.Ground.V.Component.) %>%
  mutate(delta_t = sst - t2m)


save(pts_ann, file = "alt_pts_mb_ann.RData")

#investigate NAs for sst
X11();maps::map("world")
points(pts_ann[is.na(pts_ann$sst),c("location.long","location.lat")],pch = 16,cex= 0.3,col = "red") #either 2019 (sad face) or data from points over land

pts_ann_no_NA <- pts_ann %>% 
  drop_na() #the sst NANs and also some all NA rows...

#add wind support and crosswind

segs_w <- lapply(split(pts_ann,pts_ann$seg_id), function(x){ #for each segment
  x_w <- lapply(split(x, x$days_to_add), function(y){ #for each alternative version of the segment
    #calculate heading from each point to the endpoint
    if(nrow(y) < 2){
      y_w <- mutate(heading = NA,
                    wind_support_950 = NA,
                    cross_wind_950 = NA,
                    wind_support_10m = NA,
                    cross_wind_10m = NA)
      y_w
    } else {
      y_w <- y %>% 
        mutate(heading = NCEP.loxodrome.mod(lat1=location.lat,lat2=tail(location.lat,1),lon1=location.long,lon2=tail(location.long,1)),
               wind_support_950 = wind_support(u=u950,v=v950,heading= heading),
               cross_wind_950 = cross_wind(u=u950,v=v950,heading= heading),
               wind_support_10m = wind_support(u=u10m,v=v10m,heading= heading),
               cross_wind_10m = cross_wind(u=u10m,v=v10m,heading= heading))
      y_w
    }
    
  }) %>% 
    reduce(rbind)
}) 

save(segs_w, file = "alt_pts_ann_w.RData") # a list

##### STEP 6: compare observed to best condition #####

load("alt_pts_ann_w.RData") #segs_w
#make sure all strata have 29 versions? or doesn't matter at this stage?

#are wind 950 and wind10m correlated?
#calculate variables
segs_avg <- lapply(segs_w,function(x){
  x_avg <- x %>% 
    group_by(days_to_add) %>% 
    summarise(avg_ws_950 = mean(wind_support_950, na.rm = T), 
              avg_abs_cs_950 = mean(abs(cross_wind_950), na.rm = T),
              avg_ws_10 = mean(wind_support_10m, na.rm = T),
              avg_abs_cs_10 = mean(abs(cross_wind_10m), na.rm = T),
              avg_delta_t = mean(delta_t, na.rm = T),
              cu_ws_950 = sum(abs(cross_wind_950), na.rm = T),
              cu_abs_cs_950 = sum(wind_support_950, na.rm = T),
              cu_ws_10 = sum(wind_support_10m, na.rm = T),
              cu_abs_cs_10 = sum(abs(cross_wind_10m), na.rm = T),
              cu_delta_t = sum(delta_t, na.rm = T),
              length = head(length,1),
              zone = head(zone,1),
              track = head(track,1),
              seg_id = head(seg_id,1),
              species = head(species,1),
              obs_id = head(obs_id,1),
              season = head(season,1)) %>% 
    ungroup()
  
  x_avg
})

#calc observed statistics

obs_st <- lapply(segs_avg,function(x){ #for each segment
  obs_d_avg_delta_t <- x[x$days_to_add == 0, "avg_delta_t"] - colMeans(x[x$days_to_add != 0, "avg_delta_t"])
  
  obs_d_avg_ws_950 <- x[x$days_to_add == 0, "avg_ws_950"] - colMeans(x[x$days_to_add != 0, "avg_ws_950"])
  obs_d_avg_cw_950 <- x[x$days_to_add == 0, "avg_abs_cs_950"] - colMeans(x[x$days_to_add != 0, "avg_abs_cs_950"])
  
  obs_d_avg_ws_10 <- x[x$days_to_add == 0, "avg_ws_10"] - colMeans(x[x$days_to_add != 0, "avg_ws_10"])
  obs_d_avg_cw_10 <- x[x$days_to_add == 0, "avg_abs_cs_10"] - colMeans(x[x$days_to_add != 0, "avg_abs_cs_10"])
  
  
  
})



#identify observed id's with less than 15 observations (the rest were removed because they produced NAs in the annotation step)
NA_obs_ids <- ann %>% 
  group_by(obs_id) %>% 
  summarise(count = n()) %>% 
  filter(count < 673) %>% 
  .$obs_id #none! interesting..... no 2019 data was included in the sample ;)

ann <- ann %>% 
  filter(!(obs_id %in% NA_obs_ids))


save(ann, file = "R_files/sample_alt_ann.RData")

