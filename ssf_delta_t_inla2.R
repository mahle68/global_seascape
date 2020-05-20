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
#update May 4. tried to include a prior for delta t to show the expectation of postive coefficient. but failed. lol. seek advice!
#update May 14. include deviation from long-term average as a variable. also, include coef of variation instead of variance... also calculate wind speed

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
library(brinla)


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
rsd <- function(x){
  cv <- sd(x, na.rm = T)/abs(mean(x, na.rm = T))
  rsd <- cv*100
  return(rsd)
}

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

# STEP 2: prepare alternative steps#####


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


# STEP 3: annotate data (observed time)#####
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
         cross_wind= cross_wind(u=u925,v=v925,heading=heading),
         wind_speed = sqrt(u925^2 + v925^2), # Pythagorean Theorem
         stratum = paste(track, seg_id, burst_id, step_id, sep = "_"),
         lat_zone = ifelse(location.lat > 30, "tmpz","twz")) %>% 
  rowwise() %>% 
  mutate(species = strsplit(group, "_")[[1]][1]) %>% 
  ungroup() %>% 
  as.data.frame()
          

# STEP 4: annotate data (40 year data)#####
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


#calculate long-term metrics and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/all_ssf_40yrs_1hr/",pattern = ".csv", full.names = T) 

ann_cmpl <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u=u925,v=v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading),
         abs_cross_wind = abs(cross_wind(u=u925,v=v925,heading=heading)),
         wind_speed = sqrt(u925^2 + v925^2)) %>% 
  group_by(row_id) %>% 
  summarise_at(c("delta_t", "wind_speed", "wind_support", "abs_cross_wind", "u925", "v925"), #before calculating these, investigate why we have NAs??
               list(avg = ~mean(., na.rm = T), var = ~var(., na.rm = T), rsd = ~rsd(.))) %>% 
  full_join(ann, by = "row_id") %>% 
  as.data.frame()

save(ann_cmpl, file = "ssf_input_ann_cmpl_1hr.RData")

# STEP 5: data exploration#####

load("ssf_input_ann_cmpl_1hr.RData") #ann_cmpl

#for plotting
ann_cmpl$species <- factor(ann_cmpl$species)
ann_cmpl$abs_cross_wind <- abs(ann_cmpl$cross_wind)

#plot
X11(width = 15, height = 10);par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in c("wind_speed_avg", "abs_cross_wind_avg","delta_t_avg", "wind_support_avg")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "wind_speed_avg" & j == "tmpz"){
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

X11(width = 15, height = 10);par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in c("wind_speed_var", "abs_cross_wind_var","delta_t_var", "wind_support_var")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "wind_speed_var" & j == "tmpz"){
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

X11(width = 15, height = 10);par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in c("wind_speed_rsd", "abs_cross_wind_rsd","delta_t_rsd", "wind_support_rsd")){
  for(j in c("tmpz", "twz")){ 
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
    if(i == "wind_speed_rsd" & j == "tmpz"){
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
mtext("40-yr relative standard deviation (%) at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(4,2), oma = c(0,0,3,0))
for(i in c("wind_speed", "abs_cross_wind","delta_t", "wind_support")){
  for(j in c("tmpz", "twz")){ 
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
  if(i == "wind_speed" & j == "tmpz"){
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

#plot relationship between wind and delta t
X11(); par(mfrow=c(3,2))
plot(delta_t ~ wind_support, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)

plot(delta_t_var ~ wind_support_var, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t_var ~ wind_speed_var, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t_var ~ abs_cross_wind_var, data = all_data[all_data$used ==1,], col= all_data[all_data$used ==1,"species"], pch = 16, cex = 0.7)

X11(); par(mfrow=c(3,1))
plot(delta_t ~ wind_support, data = all_data, col= factor(all_data$used), pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed, data = all_data, col= factor(all_data$used), pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind, data = all_data, col = factor(all_data$used), pch = 16, cex = 0.7)

X11(); par(mfrow=c(3,1))
plot(delta_t ~ wind_support_var, data = all_data, col=factor(all_data$used), pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed_var, data = all_data, col= factor(all_data$used), pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind_var, data = all_data, col= factor(all_data$used), pch = 16, cex = 0.7)



###conclusion: 
#correlation
ann_cmpl %>% 
  dplyr::select(c(2:5,8:11,36:39)) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat. avg delta_t and delta_t. avg_ws and var_delta_t
#correlated: wind support var & wind speed var and cross wind var, crosswind var and wind speed var.  

#z-transform
all_data <- ann_cmpl %>% 
  #group_by(species) 
  mutate_at(c(2:5,8:11,36:39),
            list(z = ~scale(.))) %>%
  as.data.frame() 

save(all_data, file = "ssf_input_ann_1hr_z.RData")

# STEP 6: modeling#####
load("ssf_input_ann_1hr_z.RData") #all_data

############# inla
##inspiration from muff et al supp material https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.r?sequence=22&isAllowed=y
##and https://ourcodingclub.github.io/tutorials/inla/

#priors are set using hyper. theta is name assigned to the hyperparametr; initial: initial values of the hyperparameter in the interla scale; fixed: whether to keep the hyperparameter fixed
#prior: name of prior distribution; param: parameter values of prior distribution (mean, precision)
#lat zone is correlated with other variables. dont inlcude

#quick and dirty clogit
m1 <- clogit(used ~ delta_t_z * wind_support_z + cross_wind_z + delta_t_z * wind_speed_z + 
                        strata(stratum), data = all_data )
m1

m2 <- clogit(used ~ delta_t_z * wind_support_z + delta_t_z * cross_wind_z + delta_t_var_z * wind_support_var_z + 
               strata(stratum), data = all_data )
m2 #lower AIC

m3 <- clogit( used ~ delta_t_z * wind_support_z + delta_t_var_z * wind_support_var_z +
               strata(stratum), data = all_data )
m3 #higher AIC than m2

#choose the osprey as the reference level and compare other vars to that
all_data$species <- relevel(all_data$species, "O")

#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = species,
         species2 = species,
         species3 = species,
         species4 = species,
         species5 = species,
         species6 = species,
         stratum = factor(stratum))

# include a nested random effect for individuals within species
#create a var for individual id
all_data <- all_data %>% 
  rowwise() %>% 
  mutate(ind = strsplit(track, "_")[[1]][1]) %>% 
  as.data.frame()

#create a model matrix for the nested effect
Z_ind <- as(model.matrix(~ 0 + species:ind, data = all_data), "Matrix") #this notation means that ind is nested in species
all_data$indinsp <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]})) #create an index variable
all_data$indinsp2 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))
all_data$indinsp3 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))
all_data$indinsp4 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))
all_data$indinsp5 <- as.factor(apply(Z_ind, 1, function(x){names(x)[x == 1]}))

# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

### species as fixed effect
formula0 <- used ~ -1 + delta_t_z * wind_support_z +  delta_t_z * cross_wind_z + delta_t_var_z * wind_support_var_z + species +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp5, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m0 <- inla(formula0, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,
#control.compute = list(mlik = T,dic = T))
#high estimates of precision could indicate that these parameters are poorly identified and may require a more informative prior

save(m0, file = "inla_models/m_0.RData")

summary(m0)

###full model
formula1 <- used ~ -1 + delta_t_z * wind_support_z +  delta_t_z * cross_wind_z + delta_t_var_z * wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp5, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m1 <- inla(formula1, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,
#control.compute = list(mlik = T,dic = T))
#high estimates of precision could indicate that these parameters are poorly identified and may require a more informative prior

save(m1, file = "inla_models/full_m.RData")

summary(m1)

#plot the posterior densities
bri.hyperpar.plot(m1) #summary of hyperparameters in SD scale (converts precision to sd)
Efxplot(m1) + theme_bw()


  ##remove crosswind. add interaction with lat zone.... results look weird
formula2b <- used ~ -1 + delta_t_z * lat_zone + wind_support_z + var_delta_t_40_z * lat_zone + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m2b <- inla(formula2b, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m2b, file = "inla_models/m2b.RData")

Efxplot(list(m1,m2,m2b)) + theme_bw()

##remove crosswind. improves compared to m1
formula2b <- used ~ -1 + delta_t_z *  wind_support_z + var_delta_t_40_z * var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m2 <- inla(formula2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m2, file = "inla_models/m2.RData")

Efxplot(list(m1,m2)) + theme_bw()

##remove crosswind and interaction of delta t and wind support... model does not improve compared to m2
formula3 <- used ~ -1 + delta_t_z + wind_support_z + var_delta_t_40_z * var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m3 <- inla(formula3, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m3, file = "inla_models/m3.RData")

Efxplot(list(m1,m2,m3)) + theme_bw()

##remove crosswind and interaction of variances of delta t and wind support... higher mlik than m3
formula3b <- used ~ -1 + delta_t_z * wind_support_z + var_delta_t_40_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m3b <- inla(formula3b, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m3b, file = "inla_models/m3b.RData")

Efxplot(list(m1,m2,m3,m3b)) + theme_bw()

##remove crosswind and interaction of delta t and wind support and delta t... slight improvement compared to m2 and m3
formula4 <- used ~ -1 + wind_support_z + var_delta_t_40_z * var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4 <- inla(formula4, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m4, file = "inla_models/m4.RData")

Efxplot(list(m1,m2,m3,m4)) + theme_bw()

##remove crosswind and delta t, but keep interaction between wind support and delta t...
formula4b <- used ~ -1 + delta_t_z : wind_support_z + var_delta_t_40_z * var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4b <- inla(formula4b, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m4b, file = "inla_models/m4.RData")

Efxplot(list(m1,m2,m3,m4, m4b)) + theme_bw()

##remove crosswind and interaction of delta t and wind support and delta t .. improves... but maybe too extreme
formula5 <- used ~ -1 + wind_support_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m5 <- inla(formula5, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m5, file = "inla_models/m5.RData")

Efxplot(list(m1,m2,m3,m4, m5)) + theme_bw()


load("inla_models/full_m.RData")
load("inla_models/m2.RData")
load("inla_models/m3.RData")
load("inla_models/m4.RData")
load("inla_models/m5.RData")


##only include variances
formula5 <- used ~ -1 + wind_support_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m5 <- inla(formula5, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data)

save(m5, file = "inla_models/m5.RData")

Efxplot(list(m1,m2,m3,m4, m5)) + theme_bw()


#######################################################
#model with positive prior for delta t... remove some vars for faster convergence
formula2 <- used ~ -1 + wind_support_z + var_delta_t_40_z + var_ws_40_z +
  f(delta_t_z, model = "iid",
    hyper = list(theta = list(initial = log(1), fixed = F, prior = "pc.prec", param = c(3,0.05)))) +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) +
  f(species1, delta_t_z, model = "iid",  
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec",param = c(3,0.05)))) +
  f(species2, wind_support_z,  model = "iid",
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec",param = c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec", param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec",param = c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec",param = c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper = list(theta = list(initial = log(1),fixed = F,prior = "pc.prec",param = c(3,0.05))))

m2 <- inla(formula2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           #control.predictor=list(compute=TRUE),
           data = all_data)


##change prior in teh model call
m3 <- inla(formula1, family ="Poisson", 
           control.fixed = list(
             mean = list(delta_t_z = 0.3, default = mean.beta),
             prec = list(default = prec.beta)),
           #control.predictor=list(compute=TRUE),
           data = all_data)


summary(m3)
bri.hyperpar.plot(m3) 
bri.fixed.plot(m3)

Efxplot(list(m1,m2,m3)) + theme_bw()

#posterior mode of variances
marginals_mode <- sapply(m2$marginals.hyperpar,
                         function(x)
                           inla.mmarginal(inla.tmarginal(function(x) 1/x, x)))

names(marginals_mode) <- sapply(as.vector(as.character(names(marginals_mode))),
                                function(y) gsub("Precision", x=y, "Mode of variance"))
marginals_mode

#posterior mean of variances
marginals_mean <- sapply(m2$marginals.hyperpar,
                         function(x)
                           inla.emarginal(function(x) x, inla.tmarginal(function(x) 1/x, x)))
names(marginals_mean) <- sapply(as.vector(as.character(names(marginals_mean))),
                                function(y) gsub("Precision", x=y, "Mean of variance"))
marginals_mean

####### random effect for species
formula3 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_delta_t_40_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_delta_t_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, var_ws_40_z, model = "iid",
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

formula4 <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m4 <- inla(formula4, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,
           #control.compute = list(mlik = T,dic = T))
#high estimates of precision could indicate that these parameters are poorly identified and may require a more informative prior
summary(m4) #smalle rmarginal log-likelihood than m1

Efxplot(list(m1,m4)) + theme_bw()


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
# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 


# include random effect of species, nested for individual, and interaction term for wind support and delta t
formula9 <- used ~ -1 + delta_t_z * wind_support_z + cross_wind_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m9 <- inla(formula9, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,

summary(m9)

Efxplot(list(m1,m4,m9)) + theme_bw()

# remove cross_wind
formula9a <- used ~ -1 + delta_t_z * wind_support_z + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m9a <- inla(formula9a, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data) #,

summary(m9a)

Efxplot(list(m9,m9a)) + theme_bw()
Efxplot(m9a) + theme_bw()

# add quad for ws..... not much diff than m9
formula9b <- used ~ -1 + delta_t_z * wind_support_z + I(wind_support_z^2) + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m9b <- inla(formula9b, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data) #,

summary(m9b)

Efxplot(list(m9,m9a,m9b)) + theme_bw()

# different priors
formula9c <- used ~ -1 + delta_t_z * wind_support_z + I(wind_support_z^2) + var_ws_40_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid",  
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(indinsp2, wind_support_z,  model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(indinsp4, var_ws_40_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

m9c <- inla(formula9c, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data) #,

summary(m9c )

Efxplot(list(m9,m9a,m9b)) + theme_bw()

##### remove cross wind
 

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

