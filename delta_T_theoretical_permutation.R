# This script investigates theoretically the hypothesis that there is higher vatiation in wind than delta T in the trade wind zone and therefore route selection (time/space)
# should depend more on the wind, while the opposite is true in the temperate zone.
# first version: points were categorized as used vs. available. second version: points are considered as available.
# also: add a 15 km buffer within the ocean layer before selecting random points.
# by: Elham Nourani. Dec, 4, 2019. Radolfzell, Germany..
#update: using permutation tests. Jan. 12. 2020

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

cv <- function(x){
  cv <- sd(x, na.rm = T) / mean(x, na.rm = T)
  return(cv)
}

rel_sd <- function(x){
  rel_sd <- sd(x,na.rm = T) / abs(mean(x, na.rm = T))
  return(rel_sd)
}

source("R_files/alt_pts_temporal.R")

# STEP 1: create a dataset for both autumn and spring #####

#select points in space
ocean <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp") %>% 
  slice(2) #remove the caspian sea


twz <- st_crop(ocean,xmin = -180, ymin = 0, xmax = 180, ymax = 30)  %>% #trade-wind zone N
  st_sample(5000) %>%
  st_coordinates() %>% 
  as.data.frame() %>% 
  mutate(zone = "tradewind")

tmz <- st_crop(ocean,xmin = -180, ymin = 30, xmax = 180, ymax = 60)  %>% #temperate zone N
  st_sample(5000) %>%  
  st_coordinates() %>% 
  as.data.frame() %>%
  mutate(zone = "temperate")


#assingn time values to the points and create alternative points.
dataset <- list(twz,tmz) %>%
  reduce(rbind) %>%
  group_by(zone) %>% 
  mutate(obs_id = row_number()) %>%
  ungroup() %>% 
  mutate(tz = tz_lookup_coords(lat = Y, lon = X, method = "accurate")) %>%  #find the time zone
  mutate(season = ifelse(between(obs_id,1,2500), "spring","autumn")) %>% #assign half the rows of each zone to one season and the other half to the other season
  rowwise() %>%
  #mutate(local_date_time = ifelse(between(obs_id,1,2500) | between(obs_id,7500,10000), as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(8:10, 1),sample(1:30, 1),sep = "-"),
  mutate(local_date_time = ifelse(season == "autumn", as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(8:10, 1),sample(1:30, 1),sep = "-"),
                                                                                        paste(sample(11:15, 1),"05","00",sep = ":"), sep = " "), tz = tz)),
                                  as.character(as.POSIXlt(x = paste(paste(sample(2007:2018, 1),month = sample(3:4, 1),sample(1:30, 1),sep = "-"),
                                                                    paste(sample(11:15, 1),"05","00",sep = ":"), sep = " "), tz = tz)))) %>%
  mutate(date_time = as.POSIXct(local_date_time, format = "%Y-%m-%d %H:%M:%OS",tz = tz)) %>% 
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

write.csv(dataset_mb,"R_files/thr_dataset_14_days_spr_aut.csv") #request track annotation with sst and t2m (nearest neighbour; non-forecast database),u, v and omega at 925 (bilinear)

#downloaded from movebank
dataset_env <- read.csv("movebank_annotation/thr_dataset_14_days_spr_aut.csv-481867360243216892/thr_dataset_14_days_spr_aut.csv-481867360243216892.csv",
                        stringsAsFactors = F) %>%
  drop_na() %>% #remove NAs
  rename( t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
          sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
          u_925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
          v_925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>%
  mutate(delta_t = sst - t2m) 

save(dataset_env, file = "R_files/thr_dataset_14_alt_env_spr_aut.RData")

# STEP 3: compare distribution and variances between the two zones ##### 

load("R_files/thr_dataset_14_alt_env_spr_aut.RData") #named dataset_env

#variance of alternative values
dataset_env_alt_var <- dataset_env %>%
  group_by(zone,obs_id) %>%
  summarise(av_delta_t_var = var(delta_t),
            av_u_var = var(u_925),
            av_v_var = var(v_925),
            season = head(season,1))

##### raincloud plots for variances #####

#restructure the dataframe
delta_t_var <- dataset_env_alt_var %>% 
  dplyr::select(-c(av_u_var,av_v_var)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = av_delta_t_var)
wind_u_var <- dataset_env_alt_var %>% 
  dplyr::select(-c(av_delta_t_var,av_v_var)) %>% 
  mutate(variable = "wind_u_925") %>% 
  dplyr::rename(score = av_u_var)
wind_v_var <- dataset_env_alt_var %>% 
  dplyr::select(-c(av_delta_t_var,av_u_var)) %>% 
  mutate(variable = "wind_v_925") %>% 
  dplyr::rename(score = av_v_var)

new_data_var <- rbind(delta_t_var,wind_u_var,wind_v_var)

X11()
ggplot(new_data_var, aes(x = variable, y = score, fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = zone),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~.) +
  ggtitle("variance of atmospheric conditions")

###print out the plots

for(i in unique(dataset_env$season)){
  data <- dataset_env[dataset_env$season == i,]
  pdf(paste("theoretical_results_",i, ".pdf",sep = ""), height = 7, width = 9)
  par(mfrow = c(2,3))
  plot(density(data[data$zone == "tradewind","delta_t"]),col = "red", main = "delta t")
  lines(density(data[data$zone == "temperate","delta_t"]),col = "blue")
  legend("topleft",legend = c("tradewind","temperate"), col = c("red","blue"),lty = 1, bty = "n", cex = 0.9)
  plot(density(data[data$zone == "tradewind","u_925"]),col = "red", main = "u wind")
  lines(density(data[data$zone == "temperate","u_925"]),col = "blue")
  plot(density(data[data$zone == "tradewind","v_925"]),col = "red", main = "v wind")
  lines(density(data[data$zone == "temperate","v_925"]),col = "blue")
  
  vars<- dataset_env_alt_var[dataset_env_alt_var$season == i,]
  boxplot(log(av_delta_t_var) ~ zone, data = vars)
  title("one week before and after_variances")
  boxplot(log(av_u_var) ~ zone, data = vars)#, log = "y")
  boxplot(log(av_v_var) ~ zone, data = vars)#, log = "y")
  
  dev.off()
}

##### raincloud plots for values #####

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

X11()
ggplot(new_data[new_data$season == "spring",], aes(x = variable, y = score, fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = zone),position = position_jitter(width = .05), size = 0.5, shape = 19, alpha=0.05)+
  geom_boxplot(aes(x = variable, y = score, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  ylim(-30,30) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic() +
  ggtitle("Atmospheric conditions in Spring")

X11()
ggplot(new_data[new_data$season == "autumn",], aes(x = variable, y = score, fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = zone),position = position_jitter(width = .05), size = 0.5, shape = 19, alpha = 0.05)+
  geom_boxplot(aes(x = variable, y = score, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  ylim(-30,30) + 
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  ggtitle("Atmospheric conditions in Autumn")

#plot both seasons in one
X11()
ggplot(new_data, aes(x = variable, y = score, fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = zone),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~.) +
  ggtitle("Atmospheric conditions")

#plot absolute values of the wind components. i dont care about directionality, just the strength
X11()
ggplot(new_data[new_data$variable != "delta_t",], aes(x = variable, y = abs(score), fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = abs(score), colour = zone),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = abs(score), fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~.)
ggtitle("Atmospheric conditions-absolute values for wind")

#plot scaled values, so that all variables are comparable
X11()
ggplot(new_data, aes(x = variable, y = scale(score), fill = zone)) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = scale(score), colour = zone),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = scale(score), fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~.) +
  ggtitle("Atmospheric conditions-scaled")





# MAPPING #####
mapview(list(twz,tmz))
mapview(twz$X,twz$Y)

# STEP 4: permutation test: is, within each zone, the seasonal variation in each variable higher than expected by chance? ##### 
load("R_files/thr_dataset_14_alt_env_spr_aut.RData") #named dataset_env

#Q: is variation in spread of values in spring and autumn higher than expected by chance?
#calculate observed cv... cv scales the measure of spread by the mean... so , sd proportional to the mean
obs_d_rsd <- dataset_env %>% 
  group_by(zone, season) %>% 
  summarise(rsd_delta_t = rel_sd(delta_t),
            rsd_u_925 = rel_sd(u_925),
            rsd_v_925 = rel_sd(v_925)) %>% 
  group_by(zone) %>% 
  summarise(delta_t_d_rsd_obs = abs(diff(rsd_delta_t)),
            u_d_rsd_obs = abs(diff(rsd_u_925)),
            v_d_rsd_obs = abs(diff(rsd_v_925)))

#produce random datasets, randomizing wihtin zone. shuffle season because that's the pattern that I want to remove
permutations <- 1000

#create alternative datasets  
#new_data_tmpz_ls <- lapply(1:permutations, function(x){ #for each permutation
#  new_obs_id_ls <-lapply(split(ann_z[ann_z$zone == "temperate",], ann_z[ann_z$zone == "temperate","obs_id"]), function(y){ #for each obs_id
#    rnd_obs <- sample(1:15,1) #draw a random number to be the new index for used row
#    y[rnd_obs,"used"] <- 1
#    y[-rnd_obs,"used"] <- 0
#    y
#  })
#  new_obs_id <- do.call(rbind, new_obs_id_ls)
#  new_obs_id
#}) 





