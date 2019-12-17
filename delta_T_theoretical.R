# This script investigates theoretically the hypothesis that there is higher vatiation in wind than delta T in the trade wind zone and therefore route selection (time/space)
# should depend more on the wind, while the opposite is true in the temperate zone.
# first version: points were categorized as used vs. available. second version: points are considered as available.
# also: add a 15 km buffer within the ocean layer before selecting random points.
# by: Elham Nourani. Dec, 4, 2019. Radolfzell, Germany..

library(sf)
library(tidyverse)
library(raster)
library(mapview)
library(maptools)
library(lwgeom)
library(lubridate)
library(lutz) #local time zone assignment

#raincloud plot
library(cowplot)
library(dplyr)
library(readr)

source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/R_rainclouds.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/summarySE.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/simulateData.R")
 
setwd("/home/mahle68/ownCloud/Work/Projects/delta_t")
setwd("/home/enourani/ownCloud/Work/Projects/delta_t")

wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")

source("R_files/alt_pts_temporal.R")

# STEP 1: create a dataset for both autumn and spring #####

#select points in space
ocean <- st_read("/home/mahle68/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp") %>% 
  slice(2) #remove the caspian sea


twz <- st_crop(ocean,xmin = -180, ymin = 0, xmax = 180, ymax = 30)  %>% #trade-wind zone N
  st_sample(1000) %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(zone = "tradewind")

tmz <- st_crop(ocean,xmin = -180, ymin = 30, xmax = 180, ymax = 60)  %>% #temperate zone N
  st_sample(1000) %>%  
  st_coordinates() %>% 
  as.data.frame() %>%
  mutate(zone = "temperate")


#assingn time values to the points and create alternative points.
dataset <- list(twz,tmz) %>%
  reduce(rbind) %>%
  mutate(obs_id = row_number()) %>%
  mutate(tz = tz_lookup_coords(lat = Y, lon = X, method = "accurate")) %>%  #find the time zone
  rowwise() %>%
  mutate(local_date_time = ifelse(between(obs_id,1,1000), as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(8:10, 1),sample(1:30, 1),sep = "-"),
                                                             paste(sample(11:15, 1),"05","00",sep = ":"), sep = " "), tz = tz)),
                                  as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(3:4, 1),sample(1:30, 1),sep = "-"),
                                                                    paste(sample(11:15, 1),"05","00",sep = ":"), sep = " "), tz = tz)))) %>%
  mutate(date_time = as.POSIXct(local_date_time, format = "%Y-%m-%d %H:%M:%OS",tz = tz),
         season = ifelse(between(obs_id,1,1000),"autumn", "spring")) %>%  #date_time in UTC
  ungroup() %>% #stop applying functions rowwise
  slice(rep(row_number(),14)) %>% #copy each row 14 times
  group_by(obs_id) %>% 
  mutate(date_time = ifelse(row_number() == 1, as.character(date_time),
                           as.character(date_time - lubridate::days(row_number() - 1)))) %>% 
  as.data.frame()

save(dataset,file = "R_files/thr_dataset_14_days_spr_aut.RData")

# STEP 2: annotate each point with delta T and wind #####

load("R_files/thr_dataset_14_days_spr_aut.RData")

#prep for track annotation on movebank
dataset_mb <- dataset %>%
  mutate(timestamp = paste(date_time,"000",sep = ".")) 

#rename columns
colnames(dataset_mb)[c(1,2)] <- c("location-long","location-lat")

write.csv(dataset_mb,"R_files/thr_dataset_14_days_spr_aut.csv") #request track annotation with sst and t2m (nearest neighbour; forecast database),u, v and omega at 925 (bilinear)
  
#downloaded from movebank
dataset_env <- read.csv("movebank_annotation/thr_dataset_14_day.csv-6470248576325268580.csv", stringsAsFactors = F) %>%
    drop_na() %>% #remove NAs
    rename( t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
            sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
            u_925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
            v_925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>%
  mutate(delta_t = sst - t2m) 

save(dataset_env, file = "R_files/thr_dataset_14_alt_env.RData")
  
# STEP 3: compare distribution and variances between the two zones ##### 

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

##### raincloud plots
p2 <- ggplot(dataset_env,aes(x=zone,y=delta_t))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0),adjust =2)+
  geom_point(position = position_jitter(width = .15), size = .25)+
  ylab('delta t')+xlab('zone')+theme_cowplot()+
  ggtitle('Basic plot for delta t')

ap2 <- ggplot(dataset_env,aes(x=zone,y=delta_t,fill=zone,col=zone))+
  geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .4,adjust =4)+
  geom_point(position=position_jitter(width = .15),size = 1, alpha = 0.4)+ylab('Score')+
  scale_fill_brewer(palette = "Dark2")+
  scale_colour_brewer(palette = "Dark2")+
  guides(fill = FALSE, col = FALSE)+
  ggtitle('Added jitter')
all_plot <- plot_grid(ap1, ap2, labels="AUTO")


ggsave('../figs/tutorial_R/2basic.png', width = w, height = h)


#restructure the dataframe
delta_t <- dataset_env %>% 
  dplyr::select(-c(u_925,v_925)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = delta_t)
wind_u <- dataset_env %>% 
  dplyr::select(-c(delta_t,v_925)) %>% 
  mutate(variable = "wind_u_925") %>% 
  dplyr::rename(score = u_925)
wind_v <- dataset_env %>% 
  dplyr::select(-c(delta_t,u_925)) %>% 
  mutate(variable = "wind_v_925") %>% 
  dplyr::rename(score = v_925)

new_data <- rbind(delta_t,wind_u,wind_v)


p10 <- ggplot(new_data, aes(x = variable, y = score, fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = zone),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(x = variable, y = score, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  ggtitle("Figure 10: Repeated Measures Factorial Rainclouds")

ggsave('../figs/tutorial_R/10repanvplot.png', width = w, height = h)




# MAPPING #####
mapview(list(twz,tmz))
mapview(twz$X,twz$Y)
