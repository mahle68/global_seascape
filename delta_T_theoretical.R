# This script investigates theoretically the hypothesis that there is higher vatiation in wind than delta T in the trade wind zone and therefore route selection (time/space)
# should depend more on the wind, while the opposite is true in the temperate zone.
# first version: points were categorized as used vs. available. second version: points are considered as available.
# also: add a 15 km buffer within the ocean layer before selecting random points.
# by: Elham Nourani. Dec, 4, 2019. Radolfzell, Germany.

library(sf)
library(tidyverse)
library(raster)
library(mapview)
library(maptools)
library(lwgeom)
library(lubridate)
library(lutz) #local time zone assignment

setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")

source("R_files/alt_pts_temporal.R")

# STEP 1: create a dataset for autumn #####

#select points in space
ocean <- st_read("C:/Users/mahle/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp") %>% 
  slice(2) #remove the caspian sea
  
#ocean_2 <- st_read("C:/Users/mahle/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp") %>% 
#  slice(2) %>% #remove the caspian sea
#  st_buffer(dist = -0.2) #create a buffer of 0.2 arcdegrees within the layer #selecting sample points is considerably faster if this step is skipped

twz <- st_crop(ocean,xmin = -180, ymin = 0, xmax = 180, ymax = 30)  %>% #trade-wind zone N
  st_sample(500) %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(zone = "tradewind")

tmz <- st_crop(ocean,xmin = -180, ymin = 30, xmax = 180, ymax = 60)  %>% #temperate zone N
  st_sample(500) %>%  
  st_coordinates() %>% 
  as.data.frame() %>%
  mutate(zone = "temperate")


#assingn time values to the points and create alternative points.
dataset <- list(twz,tmz) %>%
  reduce(rbind) %>%
  mutate(obs_id = row_number()) %>%
  mutate(tz = tz_lookup_coords(lat = Y, lon = X, method = "accurate")) %>%  #find the time zone
  rowwise() %>%
  mutate(local_date_time = as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(8:10, 1),sample(1:30, 1),sep = "-"),
                                                             paste(sample(11:15, 1),"00","00",sep = ":"), sep = " "), tz = tz))) %>%
  mutate(date_time = as.POSIXct(local_date_time, format = "%Y-%m-%d %H:%M:%OS",tz = tz)) %>%  #date_time in UTC
  ungroup() %>% #stop applying functions rowwise
  slice(rep(row_number(),14)) %>% #copy each row 14 times
  group_by(obs_id) %>% 
  mutate(date_time = ifelse(row_number() == 1, as.character(date_time),
                           as.character(date_time - lubridate::days(row_number() - 1)))) %>% 
  as.data.frame()

save(dataset,file = "R_files/thr_dataset_14_days.RData")

# STEP 2: annotate each point with delta T and wind #####

load("R_files/thr_dataset_14_days.RData")

#prep for track annotation on movebank
dataset_mb <- dataset %>%
  mutate(timestamp = paste(date_time,"000",sep = ".")) 

#rename columns
colnames(dataset_mb)[c(1,2)] <- c("location-long","location-lat")

write.csv(dataset_mb,"R_files/thr_dataset_14_days.csv") #request track annotation with sst and t2m (nearest neighbour),u, v and omega at 925 (bilinear)
  
#downloaded from movebank
dataset_env <- read.csv("movebank_annotation/thr_dataset_14_day.csv-6470248576325268580.csv", stringsAsFactors = F) %>%
    drop_na() %>% #remove NAs
    rename( t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
            sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
            u_925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
            v_925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>%
  mutate(delta_t = sst - t2m) 

save(dataset_env, file = "R_files/thr_dataset_14_alt_env.RData")
  
# STEP 3: compare variances between the two zones ##### 

load("R_files/thr_dataset_14_alt_env.RData")

#density plots for the two zones
X11()
par(mfrow = c(1,3))
plot(density(dataset_env[dataset_env$zone == "tradewind","delta_t"]),col = "red", main = "delta t")
lines(density(dataset_env[dataset_env$zone == "temperate","delta_t"]),col = "blue")
legend("topleft",legend = c("tradewind","temperate"), col = c("red","blue"),lty = 1, bty = "n", cex = 0.9)
plot(density(dataset_env[dataset_env$zone == "tradewind","u_925"]),col = "red", main = "u wind")
lines(density(dataset_env[dataset_env$zone == "temperate","u_925"]),col = "blue")
plot(density(dataset_env[dataset_env$zone == "tradewind","v_925"]),col = "red", main = "v wind")
lines(density(dataset_env[dataset_env$zone == "temperate","v_925"]),col = "blue")


#variance of alternative values
dataset_env_alt_var <- dataset_env %>%
  group_by(obs_id) %>%
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            zone = head(zone,1))

X11()
par(mfrow = c(1,3))
boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_var)
title("one week before and after")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_var)#, log = "y")

###print out the plots
pdf("theoretical_results.pdf", height = 7, width = 9)
par(mfrow = c(2,3))
plot(density(dataset_env[dataset_env$zone == "tradewind","delta_t"]),col = "red", main = "delta t")
lines(density(dataset_env[dataset_env$zone == "temperate","delta_t"]),col = "blue")
legend("topleft",legend = c("tradewind","temperate"), col = c("red","blue"),lty = 1, bty = "n", cex = 0.9)
plot(density(dataset_env[dataset_env$zone == "tradewind","u_925"]),col = "red", main = "u wind")
lines(density(dataset_env[dataset_env$zone == "temperate","u_925"]),col = "blue")
plot(density(dataset_env[dataset_env$zone == "tradewind","v_925"]),col = "red", main = "v wind")
lines(density(dataset_env[dataset_env$zone == "temperate","v_925"]),col = "blue")

boxplot(log(av_delta_t_var) ~ zone, data = dataset_env_alt_var)
title("one week before and after")
boxplot(log(av_u_var) ~ zone, data = dataset_env_alt_var)#, log = "y")
boxplot(log(av_v_var) ~ zone, data = dataset_env_alt_var)#, log = "y")


dev.off()

# MAPPING #####
mapview(list(twz,tmz))
mapview(twz$X,twz$Y)
