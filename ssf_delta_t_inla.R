#script to estimate the step selection function for water-crossing raptors.
#each segment is analyzed separately, so first I tried not burstifying the segments (1 hour continuous) because I'd lose a lot of points of already short segments
#but that made the distribution of turning angles and step lengths problematic. some step lenghts are too large. so, back to burstifying (Apr. 6)
#April 2. 2020. Radolfzell, Germany. Elham Nourani, PhD
#update APril 7. The ptt data (OHB and GFB) are too coarse and after thinning and burstification, no data point remains. So, processed OHB GPS data
#separately and will add to the analysis instead of OHB ptt. for GFB, use data from Open Science paper.
#update April 28. I tried the analysis with 2 hourly intervals, but much of the data was lost, and the boxplots did not show much difference. so, back to 
#one-hour intervals, but with a tolerance of 30 min (instead of 15)
#update April 28. tried 10-year avg and variances. very similar to 40 yr values. continue with 40 yr.
#also distance to coast was added, but boxplots show littel variance

library(tidyverse)
library(move)
library(sf)
library(circular)
library(CircStats)
library(fitdistrplus)
library(RNCEP)
library(lubridate)
library(mapview)
library(parallel)
library(tidyr)
library(corrr)
library(lme4)
library(MuMIn)
library(mgcv)
library(survival)
library(INLA)
library(ggregplot) #devtools::install_github("gfalbery/ggregplot")
library(maptools)


setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

#meters_proj <- CRS("+proj=moll +ellps=WGS84")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

par(mar(c(0,0,0,0)), oma = c(0,0,0,0))
maps::map("world",fil = TRUE,col = "grey85", border=NA) 
maps::map("world", ylim = c(5,35), xlim= c(120,140), fil = TRUE,col = "grey85", border=NA) # east asia

NCEP.loxodrome.na <- function (lat1, lat2, lon1, lon2) {
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
  else {
    head <-NA}
  return(head)
}

source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")


# STEP 1: prepare the input data#####
#open segments (not tracks, because tracks may be intersected by land)
load("segs_dt.RData") #segs_ann; prepared in track_based_prep_analyze_daily.R; filtered to be over 30 km and have more than 2 points
load("segs_OHB_dt.RData") #segs_ann_OHB; annotated OHB GPS data for autumn. prepared in all_data_prep_analyze.R
load("segs_EF_dt.RData") #segs_ann_EF; annotated Greek EF data for autumn (0-60 lat.) prepared in all_data_prep_analyze.R

segs_OHB_df <- segs_ann_OHB %>% 
  dplyr::select(intersect(colnames(segs_ann),colnames(segs_ann_OHB))) %>% 
  mutate(group = "OHB") %>% 
  as("Spatial") %>% 
  as.data.frame()

segs_EF_df <- segs_ann_EF %>% 
  dplyr::select(intersect(colnames(segs_ann),colnames(segs_ann_EF))) %>% 
  mutate(group = "EF") %>% 
  as("Spatial") %>%
  as.data.frame()

#remove spring and give different values to Osprey and Peregrine in each flyway
segs <- segs_ann %>% 
  dplyr::filter(season == "autumn") %>% 
  mutate(group = ifelse(species == "OHB", "OHB",
                        ifelse(species == "GFB", "GFB",
                               ifelse(species == "O" & st_coordinates(.)[,1] < -30, "O_A",
                                      ifelse(species == "O" & st_coordinates(.)[,1] > -30, "O_E",
                                             ifelse(species == "PF" & st_coordinates(.)[,1] < -30, "PF_A",
                                                    "PF_E")))))) %>% 
  as("Spatial") %>%
  as.data.frame() %>% 
  full_join(segs_OHB_df) %>% 
  full_join(segs_EF_df) %>% 
  mutate(unique_seg_id = paste(track, seg_id, sep = "_")) %>% 
  dplyr::arrange(group,seg_id,date_time)

#remove duplicate points
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(segs$seg_id),timestamps = segs$date_time,sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows
segs <- segs[-rows_to_delete,]

# STEP 2: interpolation#####

##decided not to do it. GFB points are days apart
#convert to lines
coordinates(segs) <- ~ x + y 
proj4string(segs) <- wgs


#for each consecutive points, calculate the distance between each set of points
#interpolate only for GFB

segs_GFB <- segs[segs$species == "GFB",]
segs_GFB <- spTransform(segs_GFB, meters_proj) %>% 
  st_as_sf()%>% 
  group_by(seg_id) %>% 
  mutate(elapsed_time = (lead(date_time) - date_time)/3600,
    distance_to_next = sf::st_distance(geometry,lead(geometry, default = empty), by_element = TRUE)) %>% 
  as.data.frame()
  
empty <- st_as_sfc("POINT(EMPTY)")
lapply(split(segs_GFB, segs_GFB$seg_id), function(x){
    x %>% 
    mutate(distance_to_next = sf::st_distance(
      geometry, 
      lead(geometry, default = empty), 
      by_element = TRUE))
})


segs_lines <- lapply(split(segs,segs$unique_seg_id), function(x){  #use segment ID :p
  line <- coords2Lines(x@coords, ID = x$unique_seg_id[1])
  proj4string(line) <- wgs
  line
})

#interpolate points along sea-crossing sections of trajectories (codes taken from track_analysis12.R from honey buzzard sea-crossing)

interp_pts_ls<-lapply(segs_lines,function(x){ # add a condition to only interpolate if there is a gap of certain amount....
  n <- as.integer(SpatialLinesLengths(x)/100) #determine the number of points to be interpolated. every 10 km, so the hypotenuse of the 7 km cells
  interp_pts <- spsample(x, n = n, type = "regular") #find the location of the new points
  
  interp_pts
})

lapply(interp_pts_ls, points,cex=0.051,pch=16,col="orange")

#add time data
coordinates(data_f) <-~long+lat
sea_region<-extent(-7,52,29,46) #only the northern shore
data_sea_region<-crop(data_f,sea_region) #extract points within the sea area, to save time

interp_pts_dt_ls<-rep(list(NA),length(interp_pts_ls)) #create empty list
names(interp_pts_dt_ls)<-names(interp_pts_ls) #assign names

for(i in 1:length(interp_pts_ls)){ #for each set of interpolated points, find overlapping observed points and extract the datetime
  pts<-interp_pts_ls[[i]]
  name<-names(interp_pts_ls)[[i]]
  observed_pts<-data_sea_region[data_sea_region$track == name,] #extract observation points from the track
  intersections<-intersect(observed_pts,pts) #extract only points over the sea
  
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

# STEP 3: prepare alternative steps#####


#for each species/flyway, thin the data, burstify, and produce alternative steps
#create a move list

move_ls<-lapply(split(segs,segs$group),function(x){
  x<-as.data.frame(x)
  mv<-move(x = x$x,y = x$y,time = x$date_time,data = x,animal = x$unique_seg_id,proj = wgs)
  mv
})
move_ls <- move_ls[-2] #remove GFB

start_time <- Sys.time()

used_av_ls_1hr <- lapply(move_ls,function(group){ #each group is a species/flyway combo
  #group <- mv_OHB #use this as group
  sp_obj_ls<-lapply(split(group),function(seg){ #sp_obj_ls will have the filtered and bursted segments
    
    #--STEP 1: thin the data to 1-hourly intervals
    seg_th<-seg%>%
      thinTrackTime(interval = as.difftime(1, units='hours'),
                    tolerance = as.difftime(30, units='mins')) #the unselected bursts are the large gaps between the selected ones
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    seg_th$selected <- c(as.character(seg_th@burstId),NA) #assign selected as a variable
    seg_th$burst_id <-c(1,rep(NA,nrow(seg_th)-1)) #define value for first row
    
    if(nrow(seg_th@data) == 1){
      seg_th@data$burst_id <- seg_th$burst_id
    } else {for(i in 2:nrow(seg_th@data)){
      
      if(i== nrow(seg_th@data)){
        seg_th@data$burst_id[i]<-NA
      } else
        if(seg_th@data[i-1,"selected"] == "selected"){
          seg_th@data$burst_id[i]<-seg_th@data[i-1,"burst_id"]
        } else {
          seg_th@data$burst_id[i]<-seg_th@data[i-1,"burst_id"]+1
        }
    }
    }
    #convert back to a move object (from move burst)
    seg_th <- as(seg_th,"Move")
    
    #--STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
    burst_ls<-split(seg_th,seg_th$burst_id)
    burst_ls<-Filter(function(x) length(x)>=3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls<-lapply(burst_ls,function(burst){
      burst$step_length<-c(distance(burst),NA) #
      burst$turning_angle<-c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp<-do.call(rbind,burst_ls)
    
    #reassign values
    
    if(length(bursted_sp) >= 1){
      bursted_sp$track<-seg@idData$track
      bursted_sp$group<-seg@idData$group
    }
    
    bursted_sp$seg_id<-seg@idData$seg_id 
    bursted_sp
  }) %>% 
    Filter(function(x) length(x) > 1, .) #remove segments with no observation (these have only one obs due to the assignment of segment id)
  
  #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  bursted_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>% 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl<-bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot
  # X11();par(mfrow=c(1,2))
  # hist(sl,freq=F,main="",xlab = "Step length (km)")
  # plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
  #                         rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  # 
  # hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  # plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  # 
  
  #--STEP 5: produce alternative steps
  used_av_seg <- lapply(sp_obj_ls, function(seg){ #for each segment
    
    used_av_burst <- lapply(split(seg,seg$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
        
      used_av_step <- lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point<-burst[this_point-1,] #this is the previous point, for calculating turning angle.
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #randomly generate 69 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 69, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl<-rgamma(n= 69, shape=fit.gamma1$estimate[[1]], rate= fit.gamma1$estimate[[2]])*1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing<-NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                     previous_point@coords[,1], current_point@coords[,1])
        
        #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
        current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
        rnd <- data.frame(lon = current_point_m@coords[,1] + rsl*cos(rta),lat = current_point_m@coords[,2] + rsl*sin(rta)) #for this to work, lat and lon should be in meters as well. boo. coordinates in meters?
        
        #covnert back to lat-lon proj
        rnd_sp<-rnd
        coordinates(rnd_sp)<-~lon+lat
        proj4string(rnd_sp)<-meters_proj
        rnd_sp<-spTransform(rnd_sp,wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(),70)) %>% #paste each row 69 times for the used and alternative steps
          mutate(x = c(head(x,1),rnd_sp@coords[,1]),
                 y = c(head(y,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,69)))  %>% #one hour after the start point of the step
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1=current_point$y,lat2=y,lon1=current_point$x,lon2= x)) %>% 
          as.data.frame()
        
        df
        
      }) %>% 
        reduce(rbind)
      used_av_step
    }) %>% 
      reduce(rbind)
    used_av_burst
  }) %>% 
    reduce(rbind)
  used_av_seg
})
  Sys.time() - start_time

save(used_av_ls_1hr, file = "ssf_input_all_plus_EF_1hr.RData")


# X11()
# maps::map("world",xlim = c(-75,-70), ylim = c(15,25),fil = TRUE,col = "ivory") #flyway
# points(burst,col = "grey", pch = 16, cex = 0.5)
# points(previous_point,col = "green", pch = 16, cex = 1)
# points(current_point,col = "red", pch = 16, cex = 1)
# points(rnd_sp, col = "orange", pch = 16, cex = 0.5)
# points(used_point, col = "purple", pch = 16, cex = 1)
# 
# plot(current_point)
# text(y~x, labels=row.names(df),data=df, cex=0.5, font=2)
# points(y~x, data = df)

#plotting
#r <- mapview(burst)
#r + mapview(current_point,color = "red") + mapview(previous_point, color = "green")

## investigate whether to use 1 hourly or 2 hourly data: based on spread, amount of data lost

# used_av_all_1hr <- lapply(used_av_ls_1hr, function(x){
#   x %>% 
#     dplyr::select(c("date_time", "x", "y", "burst_id", "track", "group", "seg_id", "step_id", "used", "heading")) %>% #later, add a unique step id: paste track, seg_id, burst_id and step_id. lol
#     mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
#     as.data.frame()
# }) %>% 
#   reduce(rbind)

used_av_all_1hr <- lapply(used_av_ls_1hr, function(x){
  x %>% 
    dplyr::select(c("date_time", "x", "y", "burst_id", "track", "group", "seg_id", "step_id", "used", "heading")) %>% #later, add a unique step id: paste track, seg_id, burst_id and step_id. lol
    mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
    as.data.frame()
}) %>% 
  reduce(rbind)
# 
# #have a look
# X11();par(mfrow= c(2,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_all_1hr[used_av_all_1hr$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_all_1hr[used_av_all_1hr$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_all_1hr2[used_av_all_1hr2$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_all_1hr2[used_av_all_1hr2$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_all_2hr[used_av_all_2hr$used == 0,c("x","y")], pch = 16, cex = 0.4, col = "gray55")
# points(used_av_all_2hr[used_av_all_2hr$used == 1,c("x","y")], pch = 16, cex = 0.4, col = "orange")


# STEP 4: annotate data (movebank)#####
#rename columns
colnames(used_av_all_1hr)[c(2,3)] <- c("location-long","location-lat")

write.csv(used_av_all_1hr, "ssf_input_all_1hr.csv")

#open annotated data and add wind support and crosswind
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/ssf_input_all_1hr.csv-5721307433824845494/ssf_input_all_1hr.csv-5721307433824845494.csv",
                stringsAsFactors = F) %>%
  drop_na() %>%
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>%
  mutate(row_id = row_number(),
         delta_t = sst - t2m,
         wind_support= wind_support(u=u925,v=v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading))
  
#add distance to coast (annotate separately)
ann2 <- ann %>% mutate(timestamp = paste(as.character(date_time),"000",sep = "."))
colnames(ann2)[c(3,4)] <- c("location-long","location-lat")
write.csv(ann2, "annotate_only_dist_coast.csv")
dist <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/annotate_only_dist_coast.csv-4138218756109615066/annotate_only_dist_coast.csv-4138218756109615066.csv",
                        stringsAsFactors = F) %>% 
  rename(dist_coast = NASA.Distance.to.Coast..Signed.) %>% 
  dplyr::select(c(row_id,dist_coast))

# STEP 5: annotate data (prep variance layer)#####
#prep a dataframe with 41 rows corresponding to 41 years (1979-2019), for each point. then i can calculate variance of delta t over 41 years for each point
df_40 <- ann %>% 
  dplyr::select(-c(v925,u925,t2m,sst,delta_t)) %>% 
  slice(rep(row_number(),41)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1979:2019)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
  as.data.frame()

str_sub(df_40$timestamp,1,4) <- df_40$year #replace original year with years from 1979-2019
colnames(df_40)[c(3,4)] <- c("location-long","location-lat") #rename columns to match movebank format

#break up into two parts. over 1 million rows
df_40_1 <- df_40 %>% 
  slice(1:900000)
write.csv(df_40_1, "ssf_40_all_spp_1hr_1.csv")

df_40_2 <- df_40 %>% 
  slice(900001:1800000)
write.csv(df_40_2, "ssf_40_all_spp_1hr_2.csv")

df_40_3 <- df_40 %>% 
  slice(1800001:2700000)
write.csv(df_40_3, "ssf_40_all_spp_1hr_3.csv")

df_40_4 <- df_40 %>% 
  slice(2700001:3600000)
write.csv(df_40_4, "ssf_40_all_spp_1hr_4.csv")

df_40_5 <- df_40 %>% 
  slice(3600001:nrow(.))
write.csv(df_40_5, "ssf_40_all_spp_1hr_5.csv")


#calculate variance delta-t for each point and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/all_ssf_40yrs_1hr/",pattern = ".csv", full.names = T) 

ann_40 <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u=u925,v=v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading)) %>% 
  group_by(row_id) %>%   
  summarise(avg_delta_t_40 = mean(delta_t,na.rm = T), 
            avg_ws_40 = mean(wind_support, na.rm = T),
            avg_cw_40 = mean(abs(cross_wind), na.rm = T),
            avg_u925_40 = mean(u925,na.rm = T),
            avg_v925_40 = mean(v925,na.rm = T),
            var_delta_t_40 = var(delta_t,na.rm = T),
            var_u925_40 = var(u925,na.rm = T),
            var_v925_40 = var(v925,na.rm = T),
            var_ws_40 = var(wind_support, na.rm = T),
            var_cw_40 = var(abs(cross_wind),na.rm = T)) #%>% 
  #full_join(ann, by = "row_id")
  
#calculate 10 year averages
ann_cmpl <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u=u925,v=v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading)) %>% 
  filter(between(year,2009,2019)) %>% 
  group_by(row_id) %>% 
  summarise(avg_delta_t_10 = mean(delta_t,na.rm = T), 
            avg_ws_10 = mean(wind_support, na.rm = T),
            avg_cw_10 = mean(abs(cross_wind), na.rm = T),
            avg_u925_10 = mean(u925,na.rm = T),
            avg_v925_10 = mean(v925,na.rm = T),
            var_delta_t_10 = var(delta_t,na.rm = T),
            var_u925_10 = var(u925,na.rm = T),
            var_v925_10 = var(v925,na.rm = T),
            var_ws_10 = var(wind_support, na.rm = T),
            var_cw_10 = var(abs(cross_wind),na.rm = T)) %>% 
  full_join(ann, by = "row_id") %>% 
  full_join(ann_40, by = "row_id") %>% 
  full_join(dist, by = "row_id")
  

#assign unique step-ids and species
ann_cmpl <- ann_cmpl %>% 
  mutate(stratum = paste(track, seg_id, burst_id, step_id, sep = "_")) %>% 
  rowwise() %>% 
  mutate(species = strsplit(group, "_")[[1]][1],
         lat_zone = ifelse(location.lat > 30, "tmpz","twz")) %>% 
  as.data.frame()

ann_cmpl$species <- factor(ann_cmpl$species) #for the purpose of plotting

save(ann_cmpl, file = "ssf_input_ann_1hr_10yrly_dist_coast.RData")

# STEP 6: data exploration#####


#plot

X11(width = 15, height = 10);par(mfrow= c(3,2), oma = c(0,0,3,0))
for(i in c("avg_ws_40", "avg_cw_40","avg_delta_t_40")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "avg_ws" & j == "tmpz"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 
}
mtext("40-yr averages at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(3,2), oma = c(0,0,3,0))
for(i in c("var_ws", "var_cw","var_delta_t")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "var_ws" & j == "tmpz"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 
}
mtext("40-yr variances at each point", side = 3, outer = T, cex = 1.3)

## 10 yr values

X11(width = 15, height = 10);par(mfrow= c(3,2), oma = c(0,0,3,0))
for(i in c("avg_ws_10", "avg_cw_10","avg_delta_t_10")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "avg_ws" & j == "tmpz"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 
}
mtext("10-yr averages at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(3,2), oma = c(0,0,3,0))
for(i in c("var_ws_10", "var_cw_10","var_delta_t_10")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "var_ws" & j == "tmpz"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,"species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 
}
mtext("10-yr variances at each point", side = 3, outer = T, cex = 1.3)
##

X11(width = 15, height = 10);par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in c("wind_support", "cross_wind","delta_t", "dist_coast")){
  for(j in c("tmpz", "twz")){ 
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
  if(i == "wind_support" & j == "tmpz"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
    boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$lat_zone == j,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$lat_zone == j,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 
}
mtext("values at timestamp of each point", side = 3, outer = T, cex = 1.3)

###conclusion: 
#correlation
ann_cmpl %>% 
  dplyr::select(c("var_ws_40","var_cw_40","var_delta_t_40","avg_ws_40","avg_cw_40","avg_delta_t_40", 
                  "wind_support","cross_wind","delta_t","location.lat")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

#z-transform
all_data <- ann_cmpl %>% 
  #group_by(species) # z_transform for each species separately. or not? ... huh!
  mutate_at(c("var_ws_40","var_cw_40","var_delta_t_40","wind_support","cross_wind","delta_t"),
            list(z = ~scale(.))) %>%
  as.data.frame()

save(all_data, file = "ssf_input_ann_1hr_z.RData")

# STEP 7: modeling#####
load("ssf_input_ann_1hr_z.RData") #all_data

#############using glm
#instantaneous model
form_inst <- formula(used ~ lat_zone * delta_t_z + lat_zone * wind_support_z + lat_zone * cross_wind_z)
form_inst_m <- formula(used ~ lat_zone * delta_t_z + lat_zone * wind_support_z + lat_zone * cross_wind_z + (1 | stratum) + (1 |species))

m1_inst <- glmer(form_inst_m, family = binomial(link = "cloglog"), data = all_data) #did not converge

m2_inst <- gamm(used ~ s(wind_support) + s(cross_wind) + delta_t,
           random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

m3_inst <- gamm(used ~ s(wind_support_z) + s(cross_wind_z) + delta_t_z,
                random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #lower AIC than without scaling


#contemporaneous model
form_cnt <- formula(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * var_cw_z)
form_cnt_m <- formula(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * var_cw_z + (1 | stratum) + (1 |species))

m1_cnt <- glmer(form_cnt_m, family = binomial(link = "cloglog"), data = all_data)

m2_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t,
                random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

all_data$lat_zone_f <- factor(all_data$lat_zone)
all_data$species_f <- factor(all_data$species)
all_data$stratum_f <- factor(all_data$stratum)

m3_cnt <- gamm(used ~ s(var_ws, by = lat_zone_f) + s(var_cw, by = lat_zone_f) + var_delta_t * lat_zone_f,
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #singluar

m4_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t * lat_zone_f,
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #interaction not significant

m5_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t + f(lat_zone_f),
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

#inspo from https://drmowinckels.io/blog/gamm-random-effects/
m6_cnt <- gamm(used ~ s(var_ws, bs = "cr") + s(var_cw, bs = "cr") + s(var_delta_t, bs = "cr") +
                 s(stratum_f, bs = "re") + s(species_f, bs = "re") + s(lat_zone_f, bs = "re"), #random effects can be assigned this way, instead of using the random argument
               family = binomial(link = "cloglog"), data = all_data , correlation=corAR1()) #what's the deal with the correlation?? sooo slow



form <- formula(used ~ lat_zone * delta_t_z + lat_zone * var_delta_t_z + lat_zone * wind_support_z + lat_zone * var_ws_z + 
                  lat_zone * cross_wind_z + lat_zone * var_cw_z )

 


############# gamm + poisson distribution (as suggested by Muff et al 2018)

#the weight implements different variances per species. the by commamd
# in the smoother ensures that we have one smoother for each bird species (p. 368)... but data format needs to be changed
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000)

form <- formula(used ~ delta_t_z + wind_support_z + (1|stratum_f) + (0 + delta_t_z | species) + (0 + wind_support_z | species))

TMBstr <- glmmTMB(form, family = poisson, data = all_data, doFit = F)
#fix the standard deviation of the first random term, which is the (1|stratum) component in the model equation
TMBstr$parameters$theta[1] <- log(1e3)
TMBstr$mapArg <- list(theta = factor(c(NA,1:2))) #I have no idea what is happening here. refer to Muff et al

mTMB <- glmmTMB:::fitTMB(TMBstr)

gamm(form, control = lmc, method = "REML", weights = varIdent(form = ~1 | species), family = poisson, data = all_data) #muff et al say dont use weigths for ssf

############# inla
##codes from muff et al supp material https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.r?sequence=22&isAllowed=y
##and https://ourcodingclub.github.io/tutorials/inla/

#my stuff
#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = factor(species),
         species2 = factor(species),
         species3 = factor(species),
         species4 = factor(species),
         species5 = factor(species),
         species6 = factor(species),
         stratum = factor(stratum),
         lat_zone = factor(lat_zone),
         lat_zone1 = factor(lat_zone),
         lat_zone2 = factor(lat_zone),
         lat_zone3 = factor(lat_zone),
         lat_zone4 = factor(lat_zone),
         lat_zone5 = factor(lat_zone))

#priors are set using hyper. theta is name assigned to the hyperparametr; initial: initial values of the hyperparameter in the interla scale; fixed: whether to keep the hyperparameter fixed
#prior: name of prior distribution; param: parameter values of prior distribution (mean, precision)

####### random effect for species and latitudinal zone
formula3 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  #f(species4, var_delta_t_z, model = "iid",
  #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(lat_zone1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(lat_zone2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(lat_zone3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
 # f(lat_zone4, var_delta_t_z, model = "iid",
 #  hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(lat_zone5, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

m1 <- inla(formula3, family ="Poisson",  #random effect for species
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
             data = all_data) #,
           #control.compute = list(mlik = T,dic = T)) #model evaluation criteria

#view fixed effects
m1$summary.fixed

#summary of posterior distrbution
m1$summary.hyperpar

#plot coefficients
Efxplot(list(m1,mf)) + theme_bw()

#model selection
resp <- "used" # Response variable

covar <- c("delta_t_z", 
           "wind_support_z", 
           "cross_wind_z", 
           "var_delta_t_z",
           "var_ws_z")

HostModelSel <- INLAModelSel(resp, covar, "species", "iid", "poisson", all_data)

Finalcovar <- HostModelSel$Removed[[length(HostModelSel$Removed)]] #only var_delta_t gets thrown out

#write the final formula using the final covars
f_final <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

#final model
mf <- inla(f_final, family ="Poisson",  #random effect for species
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           control.compute = list(mlik = T,dic = T))

summary(mf)
Efxplot(mf) + theme_bw()


###### include a nested random effect for individuals within species
#create a model matrix for the nested effect
#create a var for individual id
all_data <- all_data %>% 
  rowwise() %>% 
  mutate(ind = strsplit(track, "_")[[1]][1]) %>% 
  as.data.frame()

Z_ind <- as(model.matrix(~ 0 + species:ind, data = all_data), "Matrix") #this notation means that ind is nested in species
all_data$indinsp <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]})) #create an index variable
all_data$indinsp2 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))
all_data$indinsp3 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))
all_data$indinsp4 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))

formula4 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4 <- inla(formula4, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,
           #control.compute = list(mlik = T,dic = T))
#high estimates of precision could indicate that these parameters are poorly identified and may require a more informative prior
summary(m4)

Efxplot(list(m1,mf,m4)) + theme_bw()


#extracting point estimates of the random effects
m4$summary.random$species1 #species-specific for delta_t
m4$summary.random$species2 #species-specific for wind support
m4$summary.random$species3 #species-specific for cross_wind
m4$summary.random$species4 #species-specific for variance of wind support


#### non-scaled vars... leads to many warnings
#*** WARNING *** Eigenvalue 0 of the Hessian is -19306.1 < 0
#*** WARNING *** Set this eigenvalue to 16840.4
#*** WARNING *** This have consequence for the accurancy of
#*** WARNING *** the approximations; please check!!!
#*** WARNING *** R-inla: Use option inla(..., control.inla = list(h = h.value), ...) 
#*** WARNING *** R-inla: to chose a different  `h.value'.
formula5 <- used ~ -1 + delta_t + wind_support + cross_wind + var_ws +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m5 <- inla(formula5, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) 

summary(m5)

Efxplot(list(m1,mf,m4,m5)) + theme_bw()


#extracting point estimates of the random effects
m5$summary.random$species1 #species-specific for delta_t
m5$summary.random$species2 #species-specific for wind support
m5$summary.random$species3 #species-specific for cross_wind
m5$summary.random$species4 #species-specific for variance of wind support


####model with smooth terms for latitude. alteratively, use the spde method. will produce narrower CI. also many warnings!!!!!!! also, latitude is correlated with some atm vars
#without the inla.group, it will be assumed that measurements are regular
all_data$ws_z_grp <- inla.group(all_data$wind_support_z, n = 20, method = "quantile")
summary(all_data$ws_z_grp )

all_data$lat_grp <- inla.group(all_data$location.lat, n = 20, method = "quantile")
summary(all_data$lat_grp )

formula6 <-  used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z +
  f(lat_grp, model = "rw1", constr = F) + #use rw1 or rw2 for smooth terms
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m6 <- inla(formula6, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,

Efxplot(list(m1,mf,m4,m6)) + theme_bw()


#############
#including species as a fixed effect. does it make sense though? lol. NO!!
formula7 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z + species +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m7 <- inla(formula7, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,
#control.compute = list(mlik = T,dic = T))
#high estimates of precision could indicate that these parameters are poorly identified and may require a more informative prior
summary(m7)

Efxplot(list(m1,mf,m4,m7)) + theme_bw()


######
# include random effect of species, nested for individual, and latitudinal zone
formula8 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))+
  f(lat_zone1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(lat_zone2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(lat_zone3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(lat_zone5, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m8 <- inla(formula8, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,

summary(m8)

Efxplot(list(m1,mf,m4,m8)) + theme_bw()


######
# include random effect of species, nested for individual, and interaction term for wind support and delta t
formula9 <- used ~ -1 + delta_t_z * wind_support_z + cross_wind_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m9 <- inla(formula9, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,

summary(m9)

Efxplot(list(m1,mf,m4,m9)) + theme_bw()


##### remove cross wind
formula4b <- used ~ -1 + delta_t_z + wind_support_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4b <- inla(formula4b, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) 

summary(m4b) #improved compared to m4 :D

Efxplot(list(m1,mf,m4,m4b)) + theme_bw()

##### include var of cross wind. not important. coef close to zero
formula4c <- used ~ -1 + delta_t_z + wind_support_z + var_cw_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, var_cw_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, var_cw_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4c <- inla(formula4c, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) 

Efxplot(list(m1,mf,m4,m4b,m4c)) + theme_bw()


##### remove delta t :( .... marginal log-likelihood is not too different from the model with delta t (m4b)
formula4d <- used ~ -1 + wind_support_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4d <- inla(formula4d, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data) 

summary(m4b) #improved compared to m4 :D

Efxplot(list(m1,mf,m4,m4b, m4d)) + theme_bw()


### spde in inla #####
#create a mesh
#consider ocean_sp as mesh
load("ocean_0_60.RData") #ocean
load("ssf_input_ann_z.RData") #all_data

ocean_sp <- as(ocean, "Spatial")

pts <- all_data[all_data$used == 1,]
coordinates(pts) <-~ location.long + location.lat 
proj4string(pts) <- wgs
#pts_m <- spTransform(pts, meters_proj)

#define boundary
bdy <- unionSpatialPolygons(
  as(ocean_sp, "SpatialPolygons "), rep (1, length(ocean_sp)) #this method produces a uniform mesh
)
bdy2 <- bdy@polygons[[1]]@Polygons[[2]]@coords

mesh <- inla.mesh.2d(loc.domain = bdy, max.edge = c(15,50), offset = c(10,25))
par(mar = c(0, 0, 0, 0))
plot(mesh, asp = 1, main = "")
lines(bdy)

mesh_a <- inla.mesh.2d(loc.domain = pts, max.edge = c(5, 10)) #regular grid
mesh_b <- inla.mesh.2d(pts, max.edge = c(5, 10))

X11();par(mfrow = c(2,1), mar = c(0, 0, 0, 0))
plot(mesh_a, main = "")
points(pts, pch = 16, cex = 0.4, col = "orange")
plot(mesh_b, main = "")
plot(mesh_b, main = "")
plot(ocean_sp, add = T
     )

