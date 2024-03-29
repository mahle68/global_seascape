# Scripts for generating alternative steps for sea-crossing tracks
# This is script 1 of 6 for reproducing the results of Nourani et al 2021, ProcB.
# Elham Nourani, PhD. Jun.10. 2021
#-----------------------------------------------------------------

#to do: add the wind support functions to functions.R

library(tidyverse)
library(lubridate)
library(move)
library(sp)
library(sf)
library(parallel)
library(maptools)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

# ---------- STEP 1: load data #####

#load sea-crossing segments of migratory tracks
load("2021/public/move_ls.RData") #move_ls; this list contains sea-crossing tracks. There is one move object for each unique species-flyway

wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84") #Mollweide projection (in meters) for accurate calculation of length

source("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/global_seascape/functions.R")

# ---------- STEP 2: generate alternative steps #####

hr <- 60 #minutes; determine the sub-sampling interval
tolerance <- 15 #minutes; tolerance for sub-samplling
n_alt <- 150 #number of alternative steps.

#prepare cluster for parallel computation
mycl <- makeCluster(10) #the number of CPUs to use (adjust this based on your machine)

clusterExport(mycl, c("move_ls", "hr", "tolerance", "n_alt","wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the ParLapply call

clusterEvalQ(mycl, { #the packages that will be used within the ParLapply call
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
})

(b <- Sys.time()) 

used_av_ls <- parLapply(mycl, move_ls, function(group){ # for each group (ie. unique species-flyway combination)
  
  sp_obj_ls <- lapply(split(group), function(track){ #for each track,
    
    #--STEP 1: thin the data 
    track_th <- track %>%
      thinTrackTime(interval = as.difftime(hr, units = 'mins'),
                    tolerance = as.difftime(tolerance, units = 'mins')) #the unselected bursts are the large gaps between the selected ones
    
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
    
    if(nrow(track_th@data) == 1){
      track_th@data$burst_id <- track_th$burst_id
    } else {for(i in 2:nrow(track_th@data)){
      #if(i == nrow(track_th@data)){
      #  track_th@data$burst_id[i] <- NA #why?!
      #} else
        if(track_th@data[i-1,"selected"] == "selected"){
          track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"]
        } else {
          track_th@data$burst_id[i] <- track_th@data[i-1,"burst_id"] + 1
        }
    }
    }
    
    #convert back to a move object (from move burst)
    track_th <- as(track_th,"Move")
    
    #--STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst.
    burst_ls <- split(track_th, track_th$burst_id)
    burst_ls <- Filter(function(x) length(x) >= 3, burst_ls) #remove bursts with less than 3 observations
    
    burst_ls <- lapply(burst_ls, function(burst){
      burst$step_length <- c(distance(burst),NA) 
      burst$turning_angle <- c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp <- do.call(rbind, burst_ls)
    
    #reassign values
    if(length(bursted_sp) >= 1){
      bursted_sp$track<-track@idData$track
      bursted_sp$species<-track@idData$species
    }
    
    #bursted_sp$track<-track@idData$seg_id 
    
    bursted_sp
    
  }) %>% 
    Filter(function(x) length(x) > 1, .) #remove segments with no observation
  
  #--STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  bursted_df <- sp_obj_ls %>%  
    reduce(rbind) %>% 
    as.data.frame() %>% 
    dplyr::select(-c("coords.x1","coords.x2"))
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan. OR use circular::mean.circular
  mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot turning angle and step length distributions
  pdf(paste0("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/ssf_plots/", hr, "_", tolerance, "/",group@idData$group[1], ".pdf"))
  par(mfrow=c(1,2))
  hist(sl, freq=F, main="", xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  dev.off()
  
  #diagnostic plots for step length distribution
  pdf(paste0("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/ssf_plots/", hr, "_", tolerance, "/",group@idData$group[1], "_diag.pdf"))
  plot(fit.gamma1)
  dev.off()
  
  #--STEP 5: produce alternative steps
  used_av_track <- lapply(sp_obj_ls, function(track){ #for each trackment
    
    lapply(split(track,track$burst_id),function(burst){ #for each burst,
      
      #assign unique step id
      burst$step_id <- 1:nrow(burst)
      
      
      lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point <- burst[this_point-1,] #this is the previous point, for calculating turning angle.
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                          previous_point@coords[,1], current_point@coords[,1])
        
        current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
        
        #randomly generate n alternative points
        rnd <- data.frame(turning_angle = as.vector(rvonmises(n = n_alt, mu = mu, kappa = kappa)), #randomly generate n step lengths and turning angles
                          step_length = rgamma(n = n_alt, shape = fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000) %>% 
          #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
          mutate(lon = current_point_m@coords[,1] + step_length*cos(turning_angle),
                 lat = current_point_m@coords[,2] + step_length*sin(turning_angle))
        
        
        #covnert back to lat-lon proj
        rnd_sp <- rnd
        coordinates(rnd_sp) <- ~lon+lat
        proj4string(rnd_sp) <- meters_proj
        rnd_sp <- spTransform(rnd_sp, wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(), n_alt+1)) %>% #paste each row n_alt times for the used and alternative steps
          mutate(location.long = c(head(x,1),rnd_sp@coords[,1]), #the coordinates were called x and y in the previous version
                 location.lat = c(head(y,1),rnd_sp@coords[,2]),
                 turning_angle = c(head(turning_angle,1),deg(rnd_sp$turning_angle)),
                 step_length = c(head(step_length,1),rnd_sp$step_length),
                 used = c(1,rep(0,n_alt)))  %>%
          dplyr::select(-c("x","y")) %>% 
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1 = current_point@coords[,2], lat2 = location.lat, lon1 = current_point@coords[,1], lon2 = location.long)) %>% 
          as.data.frame()
        
        df
        
      }) %>% 
        reduce(rbind)
      
    }) %>% 
      reduce(rbind)
    
  }) %>% 
    reduce(rbind)
  used_av_track
})

Sys.time() - b 

stopCluster(mycl) 


save(used_av_ls, file = "2021/public/ssf_input_all_60_15_150.RData")

# ---------- STEP 3: annotation#####

#to submit the data for annotation on Movebank, we need to create csv files with nrow() < 1e6. The timestamps should have miliseconds, and lat and lon columns should be names according to the Movebank requirements.

#create one dataframe with movebank specifications
used_av_df <- lapply(c(1:length(used_av_ls)), function(i){
  
  data <- used_av_ls[[i]] %>% 
    dplyr::select( c("date_time", "location.lat", "location.long", "selected", "species",  "burst_id", "step_length", "turning_angle", "track", "step_id", "used", "heading")) %>% 
    mutate(timestamp = paste(as.character(date_time),"000",sep = "."), #the movebank
           group = names(used_av_ls)[[i]]) %>% 
    rowwise() %>% 
    mutate(ind = strsplit(track, "_")[[1]][1], #create variable for individua id
           stratum = paste(track, burst_id, step_id, sep = "_")) %>% #create unique stratum id
    as.data.frame()
}) %>% 
  reduce(rbind)

save(used_av_df, file = "2021/public/ssf_input_all_df_60_15_150.RData")


# #have a look
# X11();par(mfrow= c(2,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_ls[used_av_ls$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_ls[used_av_ls$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_df_1hr2[used_av_df_1hr2$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_df_1hr2[used_av_df_1hr2$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_all_2hr[used_av_all_2hr$used == 0,c("x","y")], pch = 16, cex = 0.4, col = "gray55")
# points(used_av_all_2hr[used_av_all_2hr$used == 1,c("x","y")], pch = 16, cex = 0.4, col = "orange")

#rename lat and lon columns
colnames(used_av_df)[c(2,3)] <- c("location-lat", "location-long")

write.csv(used_av_df, "2021/public/ssf_input_df_60_15_150.csv") 

#submit this file to movebank and annotate with the following variables:
#   Name: ECMWF ERA5 PL U Wind
# Description: Velocity of the east-west (zonal) component of wind. Positive values indicate west to east flow.
# Unit: m/s
# No data values: NaN (interpolated)
# Interpolation: bilinear
# Z-Dimension: 925.0
# 
# Name: ECMWF ERA5 PL V wind
# Description: Velocity of the north-south (meridional) component of wind. Positive values indicate south to north flow.
# Unit: m/s
# No data values: NaN (interpolated)
# Interpolation: bilinear
# Z-Dimension: 925.0
# 
# Name: ECMWF ERA5 SL Sea Surface Temperature
# Description: The temperature of sea water near the surface (SST). In this product (ECMWF ERA5), this parameter is a foundation SST, which means there are no variations due to the daily cycle of the sun (diurnal variations). Before September 2007, SST from the HadISST2 dataset is used and from September 2007 onwards, the OSTIA dataset is used.
# Unit: K
# No data values: NaN (interpolated)
# Interpolation: nearest-neighbour
# 
# Name: ECMWF ERA5 SL Temperature (2 m above Ground)
# Description: Air temperature 2 m above the ground or water surface. Calculated by interpolating between the lowest model level and the earth's surface, taking account of the atmospheric conditions.
# Unit: K
# No data values: NaN (interpolated)
# Interpolation: nearest-neighbour

#remove for public
# summary stats
used_av_df %>% 
  #group_by(species) %>% 
  group_by(group) %>% 
  summarise(yrs_min = min(year(date_time)),
            yrs_max = max(year(date_time)),
            n_ind = n_distinct(ind),
            n_tracks = n_distinct(track))

#visual inspection
data_sf <- used_av_df %>% 
  filter(used == 1) %>% 
  st_as_sf(coords = c(3,2), crs = wgs)

mapview(data_sf, zcol = "species")


#---- after annotation in movebank, download the annotated file and read into R
ann <- read.csv("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/public/annotation/ssf_input_df_60_15_150.csv-2201090471065357181/ssf_input_df_60_15_150.csv-2201090471065357181.csv",
                stringsAsFactors = F) %>% 
  drop_na() 

#extract startum IDs for those that have less than 50 alternative points over the sea
less_than_50 <- ann %>% 
  filter(used == 0) %>% 
  group_by(stratum) %>% 
  summarise(n = n()) %>% 
  filter(n < 50) #all are EF from Spain. It doesnt hurt to have less of that 

# retain 50 alternative steps per stratum
used <- ann %>% 
  filter(!(stratum %in% less_than_50$stratum)) %>% 
  filter(used == 1)

used_avail_50 <- ann %>% 
  filter(!(stratum %in% less_than_50$stratum)) %>% 
  filter(used == 0) %>% 
  group_by(stratum) %>% 
  sample_n(50, replace = F) %>% 
  ungroup() %>% 
  full_join(used) #append the used levels

#make sure all strata have 51 points
no_used <- used_avail_50 %>% 
  summarise(n = n()) %>% 
  filter(n < 51) # should be zero, and is! ;)

ann_50 <- used_avail_50 %>%
  filter(!(stratum %in% no_used$stratum)) %>% 
  mutate(timestamp,timestamp = as.POSIXct(strptime(timestamp,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>%
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.wind) %>%
  mutate(row_id = row_number(),
         delta_t = sst - t2m,
         wind_support= wind_support(u = u925,v = v925,heading=heading),
         cross_wind= cross_wind(u=u925,v=v925,heading=heading),
         wind_speed = sqrt(u925^2 + v925^2),
         abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading))) %>% 
  st_as_sf(coords = c("location.long", "location.lat"), crs = wgs) %>% 
  mutate(s_elev_angle = solarpos(st_coordinates(.), timestamp, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                           ifelse(s_elev_angle > 40, "high", "low"))) %>% 
  as("Spatial") %>% 
  as.data.frame()


save(ann_50, file = "2021/public/ssf_input_annotated_60_15_50.RData")

# long-term annotation (40 year variances)

#prep a dataframe with 40 rows corresponding to 40 years (1981,2020), for each point
df_40 <- ann_50 %>%# make sure this is not grouped!
  dplyr::select(-c(v925,u925,t2m,sst,delta_t)) %>% 
  slice(rep(row_number(),40)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1981:2020)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% #add miliseconds per Movebank's requirements
  as.data.frame()

str_sub(df_40$timestamp,1,4) <- df_40$year #replace original year with years from 1981-2020
colnames(df_40)[c(23,24)] <- c("location-long","location-lat") #rename columns to match movebank format


#Movebank track annotation service has a limit of one million rows per file.
#break up the dataframe into chunks with max 1 million rows
n_chunks <- ceiling(nrow(df_40)/1e6)
nrow_chunks <- round(nrow(df_40)/n_chunks)

r <- rep(1: n_chunks,  each = nrow_chunks)

chunks <- split(df_40, r)

#save as csv
lapply(c(1:length(chunks)), function(i){
  
  write.csv(chunks[[i]], paste0("2021/public/ssf_input_40yr_60_15_50_", i, "_.csv"))
  
})


#---- after annotation in movebank, download the annotated file and append to the file with instantaneous annotations

load("2021/public/ssf_input_annotated_60_15_50.RData") #ann_50

#calculate long-term metrics and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/public/annotation/40yr/",pattern = ".csv", recursive = T, full.names = T) 

ann_cmpl <- lapply(ann_40_ls, read.csv, stringsAsFactors = F) %>% 
  reduce(full_join) %>% 
  rename(sst = ECMWF.ERA5.SL.Sea.Surface.Temperature,
         t2m = ECMWF.ERA5.SL.Temperature..2.m.above.Ground.,
         u925 = ECMWF.ERA5.PL.U.Wind,
         v925 = ECMWF.ERA5.PL.V.wind) %>% 
  mutate(delta_t = sst - t2m,
         wind_support= wind_support(u = u925, v = v925, heading = heading),
         cross_wind= cross_wind(u = u925, v = v925, heading = heading),
         abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading)),
         wind_speed = sqrt(u925^2 + v925^2)) %>% 
  group_by(row_id) %>% 
  summarise_at(c("delta_t", "wind_speed", "wind_support", "abs_cross_wind", "u925", "v925"), #before calculating these, investigate why/if we have NAs??
               list(var = ~var(., na.rm = T))) %>% 
  ungroup() %>% 
  full_join(ann_50, by = "row_id") %>% 
  rename(location.long = coords.x1,
         location.lat = coords.x2) %>% 
  rowwise() %>% 
  mutate(species = strsplit(group, "_")[[1]][1],
         zone = ifelse(between(location.lat, 0, 30) | between(location.lat, 0, -30), "tradewind",
                       ifelse(between(location.lat, 30,60) | between(location.lat, -30,-60), "temperate",
                              ifelse(between(location.lat, -30,30), "tropical",
                                     "arctic")))) %>% 
  ungroup() %>% 
  as.data.frame()

save(ann_cmpl, file = "2021/public/ssf_input_ann_cmpl_60_15.RData")





