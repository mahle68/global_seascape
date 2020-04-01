#script for preparing sea-crossing segments for analysis and analyzing them..
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R and track_based_prep_analyze_daily.R
#difference with previous version: here, the entire track is one sampling unit, not each water-crossing segment.
#Elham Nourani. Mar. 30. 2020. Radolfzell, Germany


#open libraries
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(parallel)
library(mapview)
library(lutz)
library(RNCEP)


wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

##### STEP 1: calculate average conditions along each track (put all segments within each track together! :p ) #####

load("alt_pts_ann_w.RData") #segs_w

#create a dataframe
df <- segs_w %>% 
  reduce(rbind) %>% 
  as.data.frame()

#calculate variables
#detach(package:plyr) #if plyr is installed after dplyr, summarise will give the overall summary and ignore grouped structure
tracks_avg <- df %>% 
  group_by(track, days_to_add) %>% 
  summarise(avg_ws_950 = mean(wind_support_950, na.rm = T), 
            avg_abs_cw_950 = mean(abs(cross_wind_950), na.rm = T),
            avg_ws_10 = mean(wind_support_10m, na.rm = T),
            avg_abs_cw_10 = mean(abs(cross_wind_10m), na.rm = T),
            avg_delta_t = mean(delta_t, na.rm = T),
            cu_ws_950 = sum(abs(cross_wind_950), na.rm = T),
            cu_abs_cw_950 = sum(wind_support_950, na.rm = T),
            cu_ws_10 = sum(wind_support_10m, na.rm = T),
            cu_abs_cw_10 = sum(abs(cross_wind_10m), na.rm = T),
            cu_delta_t = sum(delta_t, na.rm = T),
            avg_u_950 = mean(u950, na.rm = T),
            avg_u_10 = mean(u10m, na.rm = T),
            avg_v_950 = mean(v950, na.rm = T),
            avg_v_10 = mean(v10m, na.rm = T),
            var_delta_t = var(delta_t, na.rm = T),
            var_ws_950 = var(wind_support_950, na.rm = T),
            var_abs_cw_950 = var(cross_wind_950, na.rm = T),
            length = sum(length,na.rm = T)/1000,
            #zone = head(zone,1),
            #track = head(track,1),
            #seg_id = head(seg_id,1),
            species = head(species,1),
            #obs_id = head(obs_id,1),
            season = head(season,1)) %>% 
  mutate(avg_ws_950_l = cu_ws_950/length,
         avg_abs_cw_950_l = cu_abs_cw_950/length,
         avg_delta_t_l = cu_delta_t/length) %>% 
  ungroup() %>% 
  mutate(used = ifelse(days_to_add == 0, 1,0))

save(tracks_avg,file = "alt_tracks_ann_w_avg.RData")


##### STEP 7: plot distribution of used vs available #####

#### all two weeks before and after
#for average values (sum over length of track)
#restructure the dataframe
delta_t_l <- tracks_avg %>% 
  dplyr::select(c(avg_delta_t_l,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = avg_delta_t_l)
wind_s_l <- tracks_avg %>% 
  dplyr::select(c(avg_ws_950_l,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = avg_ws_950_l)
c_wind_l <- tracks_avg %>% 
  dplyr::select(c(avg_abs_cw_950_l,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = avg_abs_cw_950_l)

new_data_l <- rbind(delta_t_l,wind_s_l,c_wind_l)


#plot both seasons in one
X11()
ggplot(new_data_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

X11() #only delta t
ggplot(delta_t_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~.) +
  ggtitle("accumulated values over length of track")


#restructure the dataframe. variance
delta_t <- tracks_avg %>% 
  dplyr::select(c(var_delta_t,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = var_delta_t)
wind_s <- tracks_avg %>% 
  dplyr::select(c(var_ws_950,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = var_ws_950)
c_wind <- tracks_avg %>% 
  dplyr::select(c(var_abs_cw_950,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = var_abs_cw_950)

new_data_v <- rbind(delta_t,wind_s,c_wind)


#plot both seasons in one
X11()
ggplot(new_data_v, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("variances")

###################################
#### one week before and after

#for average values (sum over length of track)
#restructure the dataframe
delta_t_l <- tracks_avg %>%
  filter(between(days_to_add,-7,7)) %>% 
  dplyr::select(c(avg_delta_t_l,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = avg_delta_t_l)
wind_s_l <- tracks_avg %>% 
  filter(between(days_to_add,-7,7)) %>% 
  dplyr::select(c(avg_ws_950_l,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = avg_ws_950_l)
c_wind_l <- tracks_avg %>% 
  filter(between(days_to_add,-7,7)) %>% 
  dplyr::select(c(avg_abs_cw_950_l,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = avg_abs_cw_950_l)

new_data_l <- rbind(delta_t_l,wind_s_l,c_wind_l)


X11()
ggplot(new_data_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

X11() #only delta t
ggplot(delta_t_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")


#restructure the dataframe. variance
delta_t <- tracks_avg %>% 
  filter(between(days_to_add,-7,7)) %>%
  dplyr::select(c(var_delta_t,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = var_delta_t)
wind_s <- tracks_avg %>% 
  filter(between(days_to_add,-7,7)) %>%
  dplyr::select(c(var_ws_950,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = var_ws_950)
c_wind <- tracks_avg %>% 
  filter(between(days_to_add,-7,7)) %>%
  dplyr::select(c(var_abs_cw_950,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = var_abs_cw_950)

new_data_v <- rbind(delta_t,wind_s,c_wind)

#variances
X11()
ggplot(new_data_v, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("variances")

X11() #only delta t
ggplot(delta_t, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

###################################
#### two weeks before

#for average values (sum over length of track)
#restructure the dataframe
delta_t_l <- tracks_avg %>%
  filter(between(days_to_add,0,14)) %>% 
  dplyr::select(c(avg_delta_t_l,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = avg_delta_t_l)
wind_s_l <- tracks_avg %>% 
  filter(between(days_to_add,0,14)) %>% 
  dplyr::select(c(avg_ws_950_l,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = avg_ws_950_l)
c_wind_l <- tracks_avg %>% 
  filter(between(days_to_add,0,14)) %>% 
  dplyr::select(c(avg_abs_cw_950_l,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = avg_abs_cw_950_l)

new_data_l <- rbind(delta_t_l,wind_s_l,c_wind_l)


X11()
ggplot(new_data_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

X11() #only delta t
ggplot(delta_t_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")


#restructure the dataframe. variance
delta_t <- tracks_avg %>% 
  filter(between(days_to_add,0,14)) %>%
  dplyr::select(c(var_delta_t,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = var_delta_t)
wind_s <- tracks_avg %>% 
  filter(between(days_to_add,0,14)) %>%
  dplyr::select(c(var_ws_950,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = var_ws_950)
c_wind <- tracks_avg %>% 
  filter(between(days_to_add,0,14)) %>%
  dplyr::select(c(var_abs_cw_950,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = var_abs_cw_950)

new_data_v <- rbind(delta_t,wind_s,c_wind)

#variances
X11()
ggplot(new_data_v, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("variances")

X11() #only delta t
ggplot(delta_t, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")


###################################
#### one week before

#for average values (sum over length of track)
#restructure the dataframe
delta_t_l <- tracks_avg %>%
  filter(between(days_to_add,0,7)) %>% 
  dplyr::select(c(avg_delta_t_l,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = avg_delta_t_l)
wind_s_l <- tracks_avg %>% 
  filter(between(days_to_add,0,7)) %>% 
  dplyr::select(c(avg_ws_950_l,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = avg_ws_950_l)
c_wind_l <- tracks_avg %>% 
  filter(between(days_to_add,0,7)) %>% 
  dplyr::select(c(avg_abs_cw_950_l,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = avg_abs_cw_950_l)

new_data_l <- rbind(delta_t_l,wind_s_l,c_wind_l)


X11()
ggplot(new_data_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

X11() #only delta t
ggplot(delta_t_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")


#restructure the dataframe. variance
delta_t <- tracks_avg %>% 
  filter(between(days_to_add,0,7)) %>%
  dplyr::select(c(var_delta_t,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = var_delta_t)
wind_s <- tracks_avg %>% 
  filter(between(days_to_add,0,7)) %>%
  dplyr::select(c(var_ws_950,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = var_ws_950)
c_wind <- tracks_avg %>% 
  filter(between(days_to_add,0,7)) %>%
  dplyr::select(c(var_abs_cw_950,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = var_abs_cw_950)

new_data_v <- rbind(delta_t,wind_s,c_wind)

#variances
X11()
ggplot(new_data_v, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("variances")

X11() #only delta t
ggplot(delta_t, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

###################################
#### 3 days before.. not very good

#for average values (sum over length of track)
#restructure the dataframe
delta_t_l <- tracks_avg %>%
  filter(between(days_to_add,0,3)) %>% 
  dplyr::select(c(avg_delta_t_l,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = avg_delta_t_l)
wind_s_l <- tracks_avg %>% 
  filter(between(days_to_add,0,3)) %>% 
  dplyr::select(c(avg_ws_950_l,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = avg_ws_950_l)
c_wind_l <- tracks_avg %>% 
  filter(between(days_to_add,0,3)) %>% 
  dplyr::select(c(avg_abs_cw_950_l,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = avg_abs_cw_950_l)

new_data_l <- rbind(delta_t_l,wind_s_l,c_wind_l)


X11()
ggplot(new_data_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

X11() #only delta t
ggplot(delta_t_l, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")


#restructure the dataframe. variance
delta_t <- tracks_avg %>% 
  filter(between(days_to_add,0,3)) %>%
  dplyr::select(c(var_delta_t,species, season,used)) %>% 
  mutate(variable = "delta_t") %>% 
  dplyr::rename(score = var_delta_t)
wind_s <- tracks_avg %>% 
  filter(between(days_to_add,0,3)) %>%
  dplyr::select(c(var_ws_950,species, season,used)) %>% 
  mutate(variable = "wind_support") %>% 
  dplyr::rename(score = var_ws_950)
c_wind <- tracks_avg %>% 
  filter(between(days_to_add,0,3)) %>%
  dplyr::select(c(var_abs_cw_950,species, season,used)) %>% 
  mutate(variable = "abs_cross_wind") %>% 
  dplyr::rename(score = var_abs_cw_950)

new_data_v <- rbind(delta_t,wind_s,c_wind)

#variances
X11()
ggplot(new_data_v, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("variances")

X11() #only delta t
ggplot(delta_t, aes(x = variable, y = score, fill = as.factor(used))) +
  geom_flat_violin(aes(fill = as.factor(used)),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = score, colour = as.factor(used)),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = score, fill = as.factor(used)),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic()+
  facet_grid(season~species) +
  ggtitle("accumulated values over length of track")

##### STEP 8: compare observed to alternatives #####
#calc observed statistics

load("alt_tracks_ann_w_avg.RData") #tracks_avg

avg_7 <- tracks_avg %>% 
  filter(between(days_to_add, -7, 0)) %>% 
  as.data.frame()

permutations <- 100

obs_7 <- lapply(split(avg_7,avg_7$species), function(species){
  lapply(split(species, species$track), function(track){ #because the length of track is the same in all alternative ones, use the sum of delta t
    obs <- data.frame(stat = track[track$used == 1, "cu_delta_t"] - mean(track[track$used == 0, "cu_delta_t"]),
               coefvar = sd(track[, "cu_delta_t"])/mean(track[, "cu_delta_t"], na.rm = T)*100) #including both used and alternatives
    
    rnd <- lapply(1:permutations, function(p){
      rnd_obs <- sample(c(1:7),1) #draw a random number to be the new index for used row.... dont let row 8 be chosen. it is the actual observed row
      track[rnd_obs, "cu_delta_t"] - mean(track[-rnd_obs, "cu_delta_t"])
    }) %>% 
      reduce(rbind)
    
    p_value <- sum(obs$stat <= rnd) / permutations 
    
    data.frame(obs_stat = obs$stat,coefvar = obs$coefvar, p_value,track = track$track[1],species = species$species[1], season = track$season[1])
  }) %>% 
    reduce(rbind)
}) %>% 
  reduce(rbind)

plot(density(obs_7$p_value))
boxplot(p_value~species, data = obs_7)
plot(p_value~coefvar, data = obs_7)
