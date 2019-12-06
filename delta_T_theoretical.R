# This script investigates theoretically the hypothesis that there is higher vatiation in wind than delta T in the trade wind zone and therefore route selection (time/space)
# should depend more on the wind, while the opposite is true in the temperate zone.
# by: Elham Nourani. Dec, 4, 2019. Radolfzell, Germany.

library(sf)
library(tidyverse)
library(raster)
library(mapview)
library(maptools)
library(lwgeom)
library(lubridate)
library(lutz) #local time zone assignment

setwd("C:/Users/mahle/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")

source("R_files/alt_pts_temporal.R")

# STEP 1: create a dataset for autumn #####

#select points in space
ocean <- st_read("C:/Users/mahle/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp")

set.seed(555)
twz <- st_crop(ocean,xmin = -180, ymin = 0, xmax = 180, ymax = 30)  %>% #trade-wind zone N
  st_sample(100) %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(zone = "tradewind")

  tmz <- st_crop(ocean,xmin = -180, ymin = 30, xmax = 180, ymax = 60)  %>% #temperate zone N
  st_sample(100) %>%  #get rid of any points over the caspian sea.. if any points fall on it ;)
  st_coordinates() %>% 
  as.data.frame() %>%
  mutate(zone = "temperate")

#assingn time values to the points and create alternative points.
dataset <- list(twz,tmz) %>%
    reduce(rbind) %>%
    mutate(tz = tz_lookup_coords(lat = Y, lon = X, method = "accurate")) %>%  #find the time zone
    rowwise %>%
    mutate(local_date_time = as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(8:10, 1),sample(1:30, 1),sep = "-"),
                                                               paste(sample(11:15, 1),"00","00",sep = ":"), sep = " "), tz = tz))) %>%
    mutate(date_time = as.POSIXct(local_date_time, format = "%Y-%m-%d %H:%M:%OS",tz = tz)) %>% 
    ungroup() %>% #stop applying functions rowwise
    mutate(obs_id = row_number()) %>%
    slice(rep(row_number(),15)) %>% #copy each row 15 times. 1 used, 14 alternative
    arrange(date_time) %>%
    #group_by(obs_id) %>%
    mutate(used = ifelse(row_number() == 1,1,
                         ifelse((row_number() - 1) %% 15 == 0, 1, 0))) %>% #assign used and available values
    as.data.frame()

data_ls <- lapply( split(dataset,dataset$obs_id),function(x){ #didnt manage to write this part using dplyr and purrr
    alt_times <- alt_pts_temporal(x$date_time[1],15)
    alt_x <- cbind(x,alt_times)
    alt_x
  })

dataset_alt <- do.call(rbind,data_ls)

dataset_alt[which(is.na(dataset_alt$timestamp)),] #make sure there are no NA values for date_time

save(dataset_alt,file = "R_files/thr_dataset_14_alt_days.RData")

 ########### CONSIDER A 15 KM BUFFER FOR CHOOSING THE POINTS
# STEP 2: annotate each point with delta T and wind #####

load("R_files/thr_dataset_14_alt_days.RData")

#prep for track annotation on movebank
dataset_mb <- dataset_alt %>%
  mutate(timestamp = paste(as.character(timestamp),"000",sep = ".")) 

#rename columns
colnames(dataset_mb)[c(1,2)] <- c("location-long","location-lat")

write.csv(dataset_mb,"R_files/thr_dataset_14_alt.csv") #request track annotation with sst and t2m (nearest neighbour),u, v and omega at 925 (bilinear)
  
#downloaded from movebank
dataset_env <- read.csv("movebank_annotation/thr_dataset_14_alt.csv-818519971622438358.csv", stringsAsFactors = F) %>%
    drop_na() %>% #remove NAs
    rename( t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
            sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
            u_925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
            v_925 = ECMWF.Interim.Full.Daily.PL.V.Wind,
            w_925 = ECMWF.Interim.Full.Daily.PL.Pressure.Vertical.Velocity) %>%
  mutate(delta_t = sst - t2m) 

save(dataset_env, file = "R_files/thr_dataset_14_alt_env.RData")
  
# STEP 3: compare variances between the two zones ##### 

dataset_env_delta <- lapply(split(dataset_env, dataset_env$obs_id), function(x){
  obs <- x[1,]
  #  x$delta_delta_t <- obs$delta_t - x$delta_t
  #  x$delta_u <- obs$u_925 - x$u_925
  #  x$delta_v <- obs$v_925 - x$v_925
  #  x$delta_w <- obs$w_925 - x$w_925
  
  x <- x %>% 
    mutate(delta_delta_t = obs$delta_t - delta_t,
           delta_u = obs$u_925 - u_925,
           delta_v = obs$v_925 - v_925,
           delta_w = obs$w_925 - w_925)
}) %>%
  reduce(rbind)

par(mfrow = c(2,2))
boxplot(delta_delta_t ~ zone, data = dataset_env_delta)
boxplot(delta_w ~ zone, data = dataset_env_delta)
boxplot(delta_u ~ zone, data = dataset_env_delta)
boxplot(delta_v ~ zone, data = dataset_env_delta)

#look at the variance of alternative values, regardless of the choice (week)
dataset_env_alt_var <- dataset_env %>%
  group_by(obs_id) %>%
  slice(-1) %>%  #remove the first row of each obs_id (i.e. remove the used value)
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            av_w_var = var(w_925),
            zone = head(zone,1))

windows()
par(mfrow = c(4,2))
boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_var)
title("one week before and after")
boxplot(log(av_w_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
  
#look at the variance of alternative values only BEFORE, regardless of the choice (week)
dataset_env_alt_before_var <- dataset_env %>%
  filter(period != "after") %>% 
  group_by(obs_id) %>%
  slice(-1) %>%  #remove the first row of each obs_id (i.e. remove the used value)
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            av_w_var = var(w_925),
            zone = head(zone,1))

#par(mfrow = c(2,2))
boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
title("one week before")
boxplot(log(av_w_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
  
###
#look at the variance of alternative values, regardless of the choice (month)
dataset_env_alt_var <- dataset_env %>%
  group_by(obs_id) %>%
  slice(-1) %>%  #remove the first row of each obs_id (i.e. remove the used value)
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            av_w_var = var(w_925),
            zone = head(zone,1))

windows()
par(mfrow = c(4,2))
boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_var)
title("one week before and after")
boxplot(log(av_w_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_var)#, log = "y")

#look at the variance of alternative values only BEFORE, regardless of the choice (month)
dataset_env_alt_before_var <- dataset_env %>%
  filter(period != "after") %>% 
  group_by(obs_id) %>%
  slice(-1) %>%  #remove the first row of each obs_id (i.e. remove the used value)
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            av_w_var = var(w_925),
            zone = head(zone,1))

#par(mfrow = c(2,2))
boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
title("one week before")
boxplot(log(av_w_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_before_var)#, log = "y")

# MAPPING #####
mapview(list(twz,tmz))
