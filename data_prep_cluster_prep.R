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


wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

source("wind_support_Kami.R")



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
segs_avg <- lapply(segs_w,function(x){ #for each seg_id
  x_avg <- x %>% 
    group_by(days_to_add) %>% 
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

save(segs_avg,file = "alt_pts_ann_w_avg.RData")


#calc observed statistics
load("alt_pts_ann_w_avg.RData")

obs_st <- lapply(segs_avg,function(x){ #for each segment
  obs_stats <- x[1,c(12:18)]
  obs_stats$obs_d_avg_delta_t <- x[x$days_to_add == 0, "avg_delta_t"] - mean(x[x$days_to_add != 0, "avg_delta_t"])
  
  obs_stats$obs_d_avg_ws_950 <- x[x$days_to_add == 0, "avg_ws_950"] - mean(x[x$days_to_add != 0, "avg_ws_950"])
  obs_stats$obs_d_avg_cw_950 <- x[x$days_to_add == 0, "avg_abs_cw_950"] - mean(x[x$days_to_add != 0, "avg_abs_cw_950"])
  
  obs_stats$obs_d_avg_ws_10 <- x[x$days_to_add == 0, "avg_ws_10"] - mean(x[x$days_to_add != 0, "avg_ws_10"])
  obs_stats$obs_d_avg_cw_10 <- x[x$days_to_add == 0, "avg_abs_cw_10"] - mean(x[x$days_to_add != 0, "avg_abs_cw_10"])
  
  obs_stats
})

#create alternative datasets, calc random statistics, calculate p-values #####
permutations <- 1000

#convert tibbles to dfs
segs_avg <- lapply(segs_avg,as.data.frame)

#create alternative datasets... first do a general one. with no distinction for season or zone.  
 #for each permutation
rnd_st <-lapply(segs_avg, function(x){ #for each seg_id

  new_data <- lapply(1:permutations, function(p){
    
    rnd_stats <- x[1,c(12:18)]
    rnd_obs <- sample(c(1:14,16:29),1) #draw a random number to be the new index for used row.... dont let row 15 be chosen. it is the actual observed row
    
    #rnd_stats$obs_d_avg_delta_t <- x[rownames(x) == as.character(rnd_obs), "avg_delta_t"] - colMeans(x[rownames(x) != as.character(rnd_obs), "avg_delta_t"]) #now the row with rnd_obs index is considered used.
    
    rnd_stats$obs_d_avg_delta_t <- x[rnd_obs, "avg_delta_t"] - mean(x[-rnd_obs, "avg_delta_t"], na.rm = T)
    
    rnd_stats$obs_d_avg_ws_950 <- x[rnd_obs, "avg_ws_950"] - mean(x[-rnd_obs, "avg_ws_950"], na.rm = T)
    rnd_stats$obs_d_avg_cw_950 <- x[rnd_obs, "avg_abs_cw_950"] - mean(x[-rnd_obs, "avg_abs_cw_950"], na.rm = T)
    
    rnd_stats$obs_d_avg_ws_10 <- x[rnd_obs, "avg_ws_10"] - mean(x[-rnd_obs, "avg_ws_10"], na.rm = T)
    rnd_stats$obs_d_avg_cw_10 <- x[rnd_obs, "avg_abs_cw_10"] - mean(x[-rnd_obs, "avg_abs_cw_10"], na.rm = T)
    
    rnd_stats
  }) %>% 
    reduce(rbind)
  
  new_data
}) 

#names(rnd_st) <- paste(lapply(rnd_st, "[",1,"season"), lapply(rnd_st, "[",1,"zone"), sep = "_")

#names(rnd_st) <- as.character(lapply(rnd_st, "[",1,"seg_id"))

## calculate p-values 

p_vals <-lapply(obs_st, function(x){ #for each seg_id
 #extract the corresponding random statistics
  rnds <- rnd_st[names(rnd_st) == as.character(x$seg_id)] %>% 
    reduce(rbind)
   
  P_values <- x[1,c(1:7)]
  P_values$p_delta_t <- sum(as.numeric(x$obs_d_avg_delta_t) <= rnds$obs_d_avg_delta_t) / permutations 
  P_values$p_ws_950 <- sum(as.numeric(x$obs_d_avg_ws_950) <= rnds$obs_d_avg_ws_950) / permutations 
  P_values$p_cw_950 <- sum(as.numeric(x$obs_d_avg_cw_950) <= rnds$obs_d_avg_cw_950) / permutations 
  P_values$p_ws_10 <- sum(as.numeric(x$obs_d_avg_ws_10) <= rnds$obs_d_avg_ws_10) / permutations 
  P_values$p_cw_10 <- sum(as.numeric(x$obs_d_avg_cw_10) <= rnds$obs_d_avg_cw_10) / permutations 
  
  P_values
}) %>% 
  reduce(rbind)

#plot all p_values to see
#restructure the data
dt <- p_vals %>% 
  dplyr::select(-c(p_cw_10,p_ws_10,p_cw_950,p_ws_950)) %>% 
  dplyr::rename(value = p_delta_t) %>% 
  mutate(variable = "delta_t")
ws10 <- p_vals %>% 
  dplyr::select(-c(p_cw_10,p_delta_t,p_cw_950,p_ws_950)) %>% 
  dplyr::rename(value = p_ws_10) %>% 
  mutate(variable = "wind_support_10m")
ws950 <- p_vals %>% 
  dplyr::select(-c(p_cw_10,p_ws_10,p_cw_950,p_delta_t)) %>% 
  dplyr::rename(value = p_ws_950) %>% 
  mutate(variable = "wind_support_950")
cw10 <- p_vals %>% 
  dplyr::select(-c(p_delta_t,p_ws_10,p_cw_950,p_ws_950)) %>% 
  dplyr::rename(value = p_cw_10) %>% 
  mutate(variable = "cross_wind_10")
cw950 <- p_vals %>% 
  dplyr::select(-c(p_cw_10,p_ws_10,p_delta_t,p_ws_950)) %>% 
  dplyr::rename(value = p_cw_950) %>% 
  mutate(variable = "cross_wind_950")

plotting_data <- rbind(dt,ws10,ws950,cw10,cw950)

#plot
par(mfrow = c(1,2))
boxplot(value ~ variable, data = plotting_data[plotting_data$season == "spring",])


#raincloud plots
library(ggplot2)
library(cowplot)

source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/R_rainclouds.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/summarySE.R")
source("/home/enourani/ownCloud/Work/R_source_codes/RainCloudPlots-master/tutorial_R/simulateData.R")

X11()
ggplot(plotting_data[plotting_data$zone != "both",], aes(x = variable, y = value, fill = zone)) +
  #ylim(0,150) +
  geom_flat_violin(aes(fill = zone),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA)+
  geom_point(aes(x = as.numeric(factor(variable))-.15, y = value, colour = zone),position = position_jitter(width = .05), size = 1, shape = 19, alpha = 0.1)+
  geom_boxplot(aes(x = variable, y = value, fill = zone),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  theme_classic(base_size = 20) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  facet_grid(season~.)



#plot to check
hist(rnds$obs_d_avg_delta_t, breaks = 50, col = "lightgrey")
abline(v = x$obs_d_avg_delta_t, col = "red") #random network is homogenous, that is why cv goes down the more randomization we do


####### clogit
#test: autumn in tmpz
names(segs_avg)<- paste(names(segs_avg),lapply(segs_avg, "[",1,"season"), lapply(segs_avg, "[",1,"zone"), sep = "_")


dat <- segs_avg[grep("autumn_tmpz",names(segs_avg))] %>% 
  reduce(rbind) %>% 
  mutate(used = ifelse(days_to_add == 0, 1,0))
  

#calc correlations
library(corrr)
dat[,c(2:11)] %>%
  correlate() %>% 
  stretch() %>% 
  filter(abs(r)>0.5) #wind at 10 and 950 are correlated. use just 950

dat[,c(2:3,6)] %>% #only 950 level. only averages
  correlate() %>% 
  stretch() %>% 
  filter(abs(r)>0.5) 

#scale the variables, zone-specific
dat_z<-dat%>%
  group_by(zone,season) %>% #cosider season later 
  mutate_at(c(2:3,6),
            list(z=~scale(.)))%>%
  as.data.frame()

formula <- used ~ avg_delta_t_z + avg_ws_950_z + avg_abs_cw_950_z + strata(seg_id)
#formula2 <- used ~ delta_t_z + u925_z + strata(obs_id)

################################################################ 
#all data clogit

formula <- used ~ avg_delta_t_z + avg_ws_950_z + avg_abs_cw_950_z + strata(seg_id)

for (season in c("spring","autumn")){
  for(zone in c("twz","tmpz","both")){
    
    data <-  segs_avg[grep(paste(season,zone,sep = "_"),names(segs_avg))] %>% 
      reduce(rbind) %>% 
      mutate(used = ifelse(days_to_add == 0, 1,0))
    
    data_z<-data%>% 
      mutate_at(c(2,3,6),
                list(z=~scale(.)))%>%
      as.data.frame()
    
    model1 <- clogit(formula, data = data_z)
    summary(model1)
  }
}



