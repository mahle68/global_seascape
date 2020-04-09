#script to estimate the step selection function for water-crossing raptors.
#each segment is analyzed separately, so first I tried not burstifying the segments (1 hour continuous) because I'd lose a lot of points of already short segments
#but that made the distribution of turning angles and step lengths problematic. some step lenghts are too large. so, back to burstifying (Apr. 6)
#April 2. 2020. Radolfzell, Germany. Elham Nourani, PhD
#update APril 7. The ptt data (OHB and GFB) are too coarse and after thinning and burstification, no data point remains. So, processed OHB GPS data
#separately and will add to the analysis instead of OHB ptt. for GFB, use data from Open Science paper.


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

#meters_proj <- CRS("+proj=moll +ellps=WGS84")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

maps::map("world",fil=TRUE,col="ivory") #flyway

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


setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")
source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/wind_support_Kami.R")


# STEP 1: prepare the input data#####
#open segments (not tracks, because tracks may be intersected by land)
load("segs_dt.RData") #segs_ann; prepared in track_based_prep_analyze_daily.R; filtered to be over 30 km and have more than 2 points
load("segs_OHB_dt.RData") #segs_ann_OHB; annotated OHB GPS data for autumn. prepared in all_data_prep_analyze.R

segs_OHB_df <- segs_ann_OHB %>% 
  dplyr::select(intersect(colnames(segs_ann),colnames(segs_ann_OHB))) %>% 
  mutate(group = "OHB_GPS") %>% 
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
                                                    "PF_E"))))))  %>% 
  as("Spatial") %>%
  as.data.frame() %>% 
  full_join(segs_OHB_df)%>% 
  dplyr::arrange(group,seg_id,date_time)


#create a move list
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(segs$seg_id),timestamps = segs$date_time,sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows
segs <- segs[-rows_to_delete,]

#maybe, unique seg_id needs to be track_id*seg_id? but it seems that seg_id within each species are unique

move_ls<-lapply(split(segs,segs$group),function(x){
  x<-as.data.frame(x)
  mv<-move(x = x$x,y = x$y,time = x$date_time,data = x,animal = x$seg_id,proj = wgs)
  mv
})


#for each species/flyway, thin the data, burstify, and produce alternative steps

# STEP 2: prepare alternative steps#####
start_time <- Sys.time()
used_av_ls <- lapply(move_ls[-c(1,4)],function(group){ #each group is a species/flyway combo
  #group <- mv_OHB #use this as group
  sp_obj_ls<-lapply(split(group),function(seg){ #sp_obj_ls will have the filtered and bursted segments
    
    #--STEP 1: thin the data to 1-hourly intervals
    seg_th<-seg%>%
      thinTrackTime(interval = as.difftime(1, units='hours'),
                    tolerance = as.difftime(15, units='mins')) #the unselected bursts are the large gaps between the selected ones
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
  mu<-mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl<-bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot
  #X11();par(mfrow=c(1,2))
  #hist(sl,freq=F,main="",xlab = "Step length (m)")
  #plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
  #                        rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  #hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  #plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  
  #--STEP 5: produce alternative steps
  used_av_seg <- lapply(sp_obj_ls, function(seg){ #for each segment
    
    used_av_burst <- lapply(split(seg,seg$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
        
      used_av_step <- lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point<-burst[this_point-1,] #this is the previous point, for calculating turning angle.
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #randomly generate 49 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 49, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl<-rgamma(n= 49, shape=fit.gamma1$estimate[[1]], rate= fit.gamma1$estimate[[2]])*1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
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
          slice(rep(row_number(),50)) %>% #paste each row 49 times for the used and alternative steps
          mutate(x = c(head(x,1),rnd_sp@coords[,1]),
                 y = c(head(y,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,49)))  %>% #one hour after the start point of the step
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

save(used_av_ls, file = "ssf_input_all.RData")


#get rid of alt. points that fall over land (do this later for all alternative poinst together :p)... actually, dont do this now. after track annotation, those with NA for
#sst can be easily removed :p

X11()
maps::map("world",xlim = c(-75,-70), ylim = c(15,25),fil = TRUE,col = "ivory") #flyway
points(burst,col = "grey", pch = 16, cex = 0.5)
points(previous_point,col = "green", pch = 16, cex = 1)
points(current_point,col = "red", pch = 16, cex = 1)
points(rnd_sp, col = "orange", pch = 16, cex = 0.5)
points(used_point, col = "purple", pch = 16, cex = 1)

plot(current_point)
text(y~x, labels=row.names(df),data=df, cex=0.5, font=2)
points(y~x, data = df)

#plotting
#r <- mapview(burst)
#r + mapview(current_point,color = "red") + mapview(previous_point, color = "green")

# STEP 3: annotate data (movebank)#####

used_av_all <- lapply(used_av_ls, function(x){
  x %>% 
    dplyr::select(c("date_time", "x", "y", "burst_id", "track", "group", "seg_id", "step_id", "used", "heading")) %>% #later, add a unique step id: paste track, seg_id, burst_id and step_id. lol
    mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
    as.data.frame()
}) %>% 
  reduce(rbind)

#rename columns
colnames(used_av_all)[c(2,3)] <- c("location-long","location-lat")

write.csv(used_av_all, "ssf_input_all.csv")

#open annotated data and add wind support and crosswind
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/ssf_input_all.csv-6333944159911063870/ssf_input_all.csv-6333944159911063870.csv",
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
  

# STEP 3: annotate data (prep variance layer)#####
#prep a dataframe with 41 rows corresponding to 41 years (1979-2019), for each point. then i can calculate variance of delta t over 41 years for each point
df_40 <- ann %>% 
  dplyr::select(-c(v925,u925,t2m,sst,delta_t)) %>% 
  slice(rep(row_number(),41)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1979:2019)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = "."))

str_sub(df_40$timestamp,1,4) <- df_40$year #replace original year with years from 1979-2019
colnames(df_40)[c(3,4)] <- c("location-long","location-lat") #rename columns to match movebank format

#break up into two parts. over 1 million rows
df_40_1 <- df_40 %>% 
  slice(1:700000)
write.csv(df_40_1, "ssf_40_all_spp_1.csv")

df_40_2 <- df_40 %>% 
  slice(700001:1400000)
write.csv(df_40_2, "ssf_40_all_spp_2.csv")

df_40_3 <- df_40 %>% 
  slice(1400001:nrow(.))
write.csv(df_40_3, "ssf_40_all_spp_3.csv")

#calculate variance delta-t for each point and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/movebank_annotation/all_ssf_40yrs",pattern = ".csv", full.names = T) 

ann_cmpl <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.Interim.Full.Daily.SFC.Sea.Surface.Temperature ,
         t2m = ECMWF.Interim.Full.Daily.SFC.Temperature..2.m.above.Ground.,
         u925 = ECMWF.Interim.Full.Daily.PL.U.Wind,
         v925 = ECMWF.Interim.Full.Daily.PL.V.Wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u=u925,v=v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading)) %>% 
  group_by(row_id) %>%   
  summarise(avg_delta_t = mean(delta_t,na.rm = T), 
            avg_ws = mean(wind_support, na.rm = T),
            avg_cw = mean(abs(cross_wind), na.rm = T),
            avg_u925 = mean(u925,na.rm = T),
            avg_v925 = mean(v925,na.rm = T),
            var_delta_t = var(delta_t,na.rm = T),
            var_u925 = var(u925,na.rm = T),
            var_v925 = var(v925,na.rm = T),
            var_ws = var(wind_support, na.rm = T),
            var_cw = var(abs(cross_wind),na.rm = T)) %>% 
  full_join(ann, by = "row_id")
  
#assign unique step-ids

# STEP 4: glmm#####

