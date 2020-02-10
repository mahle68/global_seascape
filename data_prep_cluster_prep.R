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
#library(parallel)
#library(lutz)getwd()
#library(move)
#library(mapview)
#library(rWind)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

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

save(segs_filtered,file = "Segs_no_land_filtered.RData") 

#assign zone to each segment
segs_filtered$twz <- as.numeric(st_within(segs_filtered,twz_sf))
segs_filtered$tmpz <- as.numeric(st_within(segs_filtered,tmpz_sf))

segs_filtered <- segs_filtered %>% 
  mutate(zone = ifelse(is.na(twz) == TRUE & is.na(tmpz) ==T, "both",
                       ifelse(is.na(twz) == TRUE & tmpz == 1, "tmpz",
                              "twz"))) %>% 
  dplyr::select(-c("twz","tmpz"))

save(segs_filtered, file = "R_files/sample_filtered_segs.RData")

##### STEP 3: annotate with date-time #####

#convert segments to points
segs_pts <- segs_filtered %>% 
  mutate(seg_id = seq(1:nrow(.))) %>% 
  st_cast("POINT") 

# lapply(split(segs_pts,segs_pts$track), function(x){ #for each set of points, find overlapping observed points and extract the datetime
#   pts<-interp_pts_ls[[i]]
#   name<-names(interp_pts_ls)[[i]]
#   observed_pts<-data_sea_region[data_sea_region$track == name,] #extract observation points from the track
#   intersections<-intersect(observed_pts,pts) #extract only points over the sea
#   
  if(length(intersections) == 0){ #if there are no intersections between the interpolated and observed points, find the nearest neighbor and extract the time
    pts_dt_ls<-apply(pts@coords,1,function(t){
      coord<-data.frame(x=t[1],y=t[2])
      t<-SpatialPoints(coord, proj4string=wgs)
      proj4string(t)<-wgs
      nearest_neighbor_indx <- pointDistance(t,observed_pts, lonlat=TRUE) %>%
        order()%>%
        head(1) #extract indices for the two points closest to the start point
      
      t$dt<-observed_pts@data[nearest_neighbor_indx,"dt"] #extract datetime of the nearest observed point and assign it to the interpolated point
      t$track<-name
      t
    })
  } else{
    #for each interpolated point, find the nearest observed point and assing the datetime 
    #assign the time of the closest point to the interpolated point
    pts_dt_ls<-apply(pts@coords,1,function(t){
      coord<-data.frame(x=t[1],y=t[2])
      t<-SpatialPoints(coord, proj4string=wgs)
      proj4string(t)<-wgs
      nearest_neighbor_indx <- pointDistance(t,intersections, lonlat=TRUE) %>%
        order()%>%
        head(1) #extract indices for the two points closest to the start point
      
      t$dt<-intersections@data[nearest_neighbor_indx,"dt"] #extract datetime of the nearest observed point and assign it to the interpolated point
      t$track<-name
      t
    })
  }
  pts_dt<-do.call(rbind, pts_dt_ls)
  interp_pts_dt_ls[[i]]<-pts_dt
  
}



#####################################
segs_pts_dt <- lapply(split(segs_pts,segs_pts$track), function(x){
  track <- as.character(x$track[1])
  x_sp <- as(x,"Spatial") 
  #x_sp$date_time <- NA
  #x_sp$season <- NA
  
  data_sp <- dataset[dataset$track == track,] #make sure dataset is spatial
  
  for (i in 1:nrow(x_sp)){
    pt <- x_sp[i,]
    #intersect <- data_df[round(data_df$location.long,2) == round(pt$location.long,2) & round(data_df$location.lat,2) == round(pt$location.lat,2),]
    #intersection <- raster::intersect(pt,data_sp)
    
    if(length(intersect) == 0){ #if there are no intersections between the interpolated and observed points, find the nearest neighbor and extract the time
      #pts_dt_ls<-apply(pt@coords,1,function(t){
      #  coord<-data.frame(x=t[1],y=t[2])
      #  t<-SpatialPoints(coord, proj4string=wgs)
      #  proj4string(t)<-wgs
        nearest_neighbor_indx <- pointDistance(pt,data_sp, lonlat=TRUE) %>%
          order()%>%
          head(1) #extract indices for the two points closest to the start point
        
        pt$date_time<-data_sp@data[nearest_neighbor_indx,"date_time"] #extract datetime of the nearest observed point and assign it to the interpolated point
        pt$track<-as.character(pt$track)
        pt$season <- data_sp@data[nearest_neighbor_indx,"season"] 
        pt
      } 
    #else{
      # #for each interpolated point, find the nearest observed point and assing the datetime 
      # #assign the time of the closest point to the interpolated point
      # pts_dt_ls<-apply(pt@coords,1,function(t){
      #   coord<-data.frame(x=t[1],y=t[2])
      #   t<-SpatialPoints(coord, proj4string=wgs)
      #   proj4string(t)<-wgs
      #   nearest_neighbor_indx <- pointDistance(t,intersections, lonlat=TRUE) %>%
      #     order()%>%
      #     head(1) #extract indices for the two points closest to the start point
      #   
      #   t$dt<-intersections@data[nearest_neighbor_indx,"dt"] #extract datetime of the nearest observed point and assign it to the interpolated point
      #   t$track<-name
      #   t
      # })
    
    if (nrow(intersect) == 1){ #if the point overlaps with an observed point over the track, assign the date_time and season
      x_df[i,"season"] <- intersect$season
      x_df[i,"date_time"] <- as.character(intersect$date_time+ min(1))
    }
    
    
    if(length(intersect) == 0){ #if there are no intersections between the interpolated and observed points, find the nearest neighbor and extract the time
      pts_dt_ls<-apply(pt@coords,1,function(t){
        coord<-data.frame(x=t[1],y=t[2])
        t<-SpatialPoints(coord, proj4string=wgs)
        proj4string(t)<-wgs
        nearest_neighbor_indx <- pointDistance(t,observed_pts, lonlat=TRUE) %>%
          order()%>%
          head(1) #extract indices for the two points closest to the start point
        
        t$dt<-observed_pts@data[nearest_neighbor_indx,"dt"] #extract datetime of the nearest observed point and assign it to the interpolated point
        t$track<-name
        t
      })
    } else{
      #for each interpolated point, find the nearest observed point and assing the datetime 
      #assign the time of the closest point to the interpolated point
      pts_dt_ls<-apply(pts@coords,1,function(t){
        coord<-data.frame(x=t[1],y=t[2])
        t<-SpatialPoints(coord, proj4string=wgs)
        proj4string(t)<-wgs
        nearest_neighbor_indx <- pointDistance(t,intersections, lonlat=TRUE) %>%
          order()%>%
          head(1) #extract indices for the two points closest to the start point
        
        t$dt<-intersections@data[nearest_neighbor_indx,"dt"] #extract datetime of the nearest observed point and assign it to the interpolated point
        t$track<-name
        t
      })
    }
    
    
 
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

##### STEP 7: create alternative tracks in time #####

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


##### STEP 8: annotate all points #####

load("R_files/sample_alt_pts_alt_time.RData") #called pts_alt

#prep for track annotation on movebank
pts_alt_mb <- pts_alt %>%
  mutate(timestamp = paste(as.character(alt_date_time),"000",sep = "."))

#rename columns
colnames(pts_alt_mb)[c(7,8)] <- c("location-long","location-lat")

write.csv(pts_alt_mb,"R_files/sample_pts_mb.csv") 

#downloaded from movebank
ann <- read.csv("movebank_annotation/sample_pts_mb.csv-3131529835871517968/sample_pts_mb.csv-3131529835871517968.csv", stringsAsFactors = F) %>%
  #drop_na() %>%# NA values are for the 2019 tracks. with a transition to ERA5, I should be able to use this data as well
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.FC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.FC.Temperature..2.m.above.Ground.,
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
  filter(count < 673) %>% 
  .$obs_id #none! interesting..... no 2019 data was included in the sample ;)

ann <- ann %>% 
  filter(!(obs_id %in% NA_obs_ids))


save(ann, file = "R_files/sample_alt_ann.RData")



