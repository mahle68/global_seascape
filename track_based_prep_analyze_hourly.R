#script for preparing sea-crossing segments for analysis and analyzing them... producing alternative tracks hours apart
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R and track_based_prep_analyze_daily.R
#Elham Nourani. Feb. 18. 2020. Radolfzell, Germany


#open libraries
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(parallel)
library(mapview)
library(lutz)
library(RNCEP)
library(corrr)


wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

source("wind_support_Kami.R")

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

##### STEP 1: create alternative tracks in time #####

load("segs_dt.RData") #segs_ann

segs_df <- segs_ann %>% 
  as("Spatial") %>% 
  as.data.frame()

#add 6-hourly alternatives for 3 days before the observed point
hours_to_add <- c(0,cumsum(rep(-6,18))) #six hourly data addition for three days prior to the observed segment

pts_alt <- segs_df %>% 
  mutate(obs_id = row_number()) %>% 
  slice(rep(row_number(),19)) %>%  #paste each row 18 time for 18 additional hours
  arrange(obs_id) %>% 
  group_by(obs_id) %>% 
  mutate(hours_to_add = hours_to_add) %>% 
  mutate(alt_date_time = date_time + hours(hours_to_add)) %>%  #use hours to add as an id for alternative segments
  ungroup() %>% 
  mutate(used = ifelse(hours_to_add == 0,1,0)) %>% 
  as.data.frame()
  
  save(pts_alt, file = "alt_pts_alt_6_hr.RData")


##### STEP 2: annotate all points #####

load("alt_pts_alt_6_hr.RData") #called pts_alt

#prep for track annotation on movebank
pts_alt_mb <- pts_alt %>%
  mutate(timestamp = paste(as.character(alt_date_time),"000",sep = ".")) %>% 
  as.data.frame()

#rename columns
colnames(pts_alt_mb)[c(11,12)] <- c("location-long","location-lat")

# this is over 6 million rows. break up into 7 files to upload to movebank
pts_alt_mb$chuncks <-c(rep(1,1e6),rep(2,1e6),rep(3,1e6),rep(4,1e6),rep(5,1e6),rep(6,1e6),rep(7,nrow(pts_alt_mb)-6e6))

lapply(split(pts_alt_mb,pts_alt_mb$chuncks),function(x){
  write.csv(x,paste("alt_pts_6hr_mb_chunk_",x$chuncks[1],".csv",sep = ""))
})

#downloaded from movebank

file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/alt_segs_6hr/", full.names = T, pattern = ".csv$") # $ means end of string
pts_ann <- lapply(file_ls,read.csv, stringsAsFactors = F) %>% 
  reduce(rbind) %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u950 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v950 = ECMWF.Interim.Full.Daily.PL.V.Wind,
         w950 = ECMWF.Interim.Full.Daily.PL.Pressure.Vertical.Velocity) %>%
  mutate(delta_t = sst - t2m)


save(pts_ann, file = "alt_pts_mb_ann_6hr.RData")
load ("alt_pts_mb_ann_6hr.RData")

#investigate NAs for sst
X11();maps::map("world")
points(pts_ann[is.na(pts_ann$sst),c("location.long","location.lat")],pch = 16,cex= 0.3,col = "red") #either 2019 (sad face) or data from points over land

pts_ann <- pts_ann %>% 
  drop_na() #the sst NANs and also some all NA rows...

#add wind support and crosswind

segs_w <- lapply(split(pts_ann,pts_ann$seg_id), function(x){ #for each segment
  x_w <- lapply(split(x, x$hours_to_add), function(y){ #for each alternative version of the segment
    #calculate heading from each point to the endpoint
    if(nrow(y) < 2){
      y_w <- y %>%   
        mutate(heading = NA,
               wind_support_950 = NA,
               cross_wind_950 = NA)
      y_w
    } else {
      y_w <- y %>% 
        mutate(heading = NCEP.loxodrome.na(lat1=location.lat,lat2=tail(location.lat,1),lon1=location.long,lon2=tail(location.long,1)),
               wind_support_950 = wind_support(u=u950,v=v950,heading= heading),
               cross_wind_950 = cross_wind(u=u950,v=v950,heading= heading))
      y_w
    }
    

    
  }) %>% 
    reduce(rbind)
}) 

save(segs_w, file = "alt_pts_ann_w_6hr.RData") # a list


##### STEP 3: calculate average conditions along segments #####

load("alt_pts_ann_w_6hr.RData") #segs_w

#calculate variables
segs_avg <- lapply(segs_w,function(x){ #for each seg_id
  x_avg <- x %>% 
    group_by(hours_to_add) %>% 
    summarise(avg_ws_950 = mean(wind_support_950, na.rm = T), 
              avg_abs_cw_950 = mean(abs(cross_wind_950), na.rm = T),
              avg_delta_t = mean(delta_t, na.rm = T),
              min_omega_950 = min(w950,na.rm = T),
              cu_ws_950 = sum(abs(cross_wind_950), na.rm = T),
              cu_abs_cw_950 = sum(wind_support_950, na.rm = T),
              cu_delta_t = sum(delta_t, na.rm = T),
              length = head(length,1),
              zone = head(zone,1),
              track = head(track,1),
              seg_id = head(seg_id,1),
              species = head(species,1),
              obs_id = head(obs_id,1),
              season = head(season,1),
              used = head(used,1)) %>% 
    ungroup()
  
  x_avg
})

save(segs_avg,file = "alt_pts_ann_w_avg_6hr.RData")



##### STEP 4: clogit for observed vs available #####


segs_avg_df <- segs_avg %>% 
  reduce(rbind) %>% 
  group_by(season,zone) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  ungroup() %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()

#correlations
segs_avg_df[,c(2:5)] %>%
  correlate() %>% 
  stretch() %>% 
  filter(abs(r)>0.5)

formula <- used ~ avg_delta_t_z + avg_ws_950_z + min_omega_950_z + avg_abs_cw_950_z + strata(seg_id)
formula <- used ~ avg_delta_t_z + avg_ws_950_z  + strata(seg_id)

full_models <- lapply(split(segs_avg_df,list(segs_avg_df$season,segs_avg_df$zone)), function(x){
  model <- clogit(formula, data = x)
  model
})

ws_delta_t_models <- lapply(split(segs_avg_df,list(segs_avg_df$season,segs_avg_df$zone)), function(x){
  model <- clogit(used ~ avg_delta_t_z + avg_ws_950_z + strata(seg_id), 
                  data = x)
  model
})

delta_t_models <- lapply(split(segs_avg_df,list(segs_avg_df$season,segs_avg_df$zone)), function(x){
  model <- clogit(used ~ avg_delta_t_z + strata(seg_id), 
                  data = x)
  model
})

#----------------------------------
#zone specific, not season
df_zone <- segs_avg %>% 
  reduce(rbind) %>% 
  group_by(zone) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  ungroup() %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()


full_models <- lapply(split(df_zone,df_zone$zone), function(x){
  model <- clogit(formula, data = x)
  summary(model)
})

#----------------------------------
#season specific, not 
df_season <- segs_avg %>% 
  reduce(rbind) %>% 
  group_by(season) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  ungroup() %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()


full_models <- lapply(split(df_season,df_season$season), function(x){
  model <- clogit(formula, data = x)
  model
})

#----------------------------------
#species specific, not season
df_sp <- segs_avg %>% 
  reduce(rbind) %>% 
  group_by(species) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  ungroup() %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()


full_models <- lapply(split(df_sp,df_sp$species), function(x){
  model <- clogit(formula, data = x)
  summary(model)
})

#---------------------------------
#only species as predictor
df <- segs_avg %>% 
  reduce(rbind) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()


glm(used ~ species, family= "binomial", data = df)

glm(used ~ avg_delta_t_z + min_omega_950_z + avg_abs_cw_950_z + avg_ws_950_z + species , family = binomial, data = ann_sc[ann_sc$zone == "tradewind",] ) #higher selection on u-wind
r.squaredLR(model)

model <- glm(used ~ avg_delta_t_z + min_omega_950_z + avg_abs_cw_950_z + avg_ws_950_z + zone + season + species , family = binomial, data = df ) #higher selection on u-wind
r.squaredLR(model)


model <- glmer(used ~ avg_delta_t_z + min_omega_950_z + avg_abs_cw_950_z + avg_ws_950_z + (1 | seg_id) + (1|species), family = binomial, data = df ) 
model <- glmer(used ~ avg_delta_t_z + avg_ws_950_z + (1 | seg_id) + (1|species), family = binomial, data = df ) 

model <- glmer(used ~ avg_delta_t_z + min_omega_950_z + avg_abs_cw_950_z + avg_ws_950_z + (1 | seg_id) + (1|species), family = binomial, data = df )

#---------------------------------
#not using scaled variables, interactions, second order, etc.
model <- clogit(used ~ poly(avg_delta_t,2) + avg_ws_950 + strata(seg_id), data = df)

model <- clogit(used ~ avg_delta_t_z*avg_ws_950_z + strata(seg_id), data = df)

#------------------------------------
#only one day before. season specific... still negative.
df_season <- segs_avg %>% 
  reduce(rbind) %>% 
  filter(hours_to_add %in% c(-21:0)) %>% 
  group_by(season) %>% 
  mutate_at(c(2:5),
            list(z=~scale(.))) %>% 
  ungroup() %>% 
  drop_na() %>%  #why do i have NAs anyway?
  as.data.frame()


full_models <- lapply(split(df_season,df_season$season), function(x){
  model <- clogit(formula, data = x)
  summary(model)
})
