#script to estimate the step selection function for water-crossing raptors.
#mostly copied from ssf_delta_t_inla2.R. For previous steps, see 2021_all_data_prep_analyze.R
#Feb 22. 2021. Radolfzell, Germany. Elham Nourani, PhD.

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
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

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

# ---------- STEP 1: prepare the input data#####
#open over water points prepared in 2021_all_data_prep_analyze.R
load("2021/all_2009_2020_overwater_points.RData") #all_oversea
#load("2021/ocean.RData") #ocean (prepared in 2021_all_data_prep_analyze.R)
#load("2021/land_no_buffer.RData") #land_no_buffer

# ---------- STEP 2: create move object#####

#remove duplicated points
oversea_df <- all_oversea %>%
  mutate(group = ifelse(species == "OHB", "OHB",
                        ifelse(species == "O" & st_coordinates(.)[,1] < -30, "O_A",
                               ifelse(species == "O" & st_coordinates(.)[,1] > -30, "O_E",
                                      ifelse(species == "EF" & st_coordinates(.)[,2] > 0, "EF_G",
                                             ifelse(species == "EF" & st_coordinates(.)[,2] < 0, "EF_S",
                                                     ifelse(species == "PF" & st_coordinates(.)[,1] < -30, "PF_A",
                                                            "PF_E")))))))  %>% 
  arrange(track,date_time) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(x = coords.x1,
         y = coords.x2)
  
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(oversea_df$track),timestamps = oversea_df$date_time),"[",-1)) #get all but the first row of each set of duplicate rows
oversea_df <- oversea_df[-rows_to_delete,]

#create a mov object. one per group. this will lead to more accurrate estimation of step lengths and turning angles
move_ls <- lapply(split(oversea_df,oversea_df$group),function(x){
  x <- as.data.frame(x)
  mv <- move(x = x$x, y = x$y, time = x$date_time, data = x, animal = x$track, proj = wgs)
  mv
})

# ---------- STEP 3: generate alternative steps#####

mycl <- makeCluster(7) 
clusterExport(mycl, c("move_ls", "wgs", "meters_proj", "NCEP.loxodrome.na")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(sf)
  library(sp)
  library(tidyverse)
  library(move)
  library(CircStats)
  library(circular)
  library(fitdistrplus)
})

(b <- Sys.time())

used_av_ls_90_60 <- parLapply(mycl, move_ls, function(group){ #each group
  
  sp_obj_ls <- lapply(split(group), function(track){
  
    #--STEP 1: thin the data to 1-hourly intervals
    track_th <- track %>%
      thinTrackTime(interval = as.difftime(60, units='mins'),
                    tolerance = as.difftime(30, units='mins')) #the unselected bursts are the large gaps between the selected ones
    #--STEP 2: assign burst IDs (each chunk of track with 1 hour intervals is one burst... longer gaps will divide the brusts) 
    track_th$selected <- c(as.character(track_th@burstId),NA) #assign selected as a variable
    track_th$burst_id <-c(1,rep(NA,nrow(track_th)-1)) #define value for first row
    
    if(nrow(track_th@data) == 1){
      track_th@data$burst_id <- track_th$burst_id
    } else {for(i in 2:nrow(track_th@data)){
      
      if(i== nrow(track_th@data)){
        track_th@data$burst_id[i] <- NA
      } else
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
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
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
  #calculate the averages (mu).steps: 1) convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu <- mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl <- bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot
  pdf(paste0("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/ssf_plots/",group@idData$group[1], ".pdf"))
  par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
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
        
        #randomly generate 50 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 150, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl <- rgamma(n = 150, shape=fit.gamma1$estimate[[1]], rate = fit.gamma1$estimate[[2]]) * 1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing <- NCEP.loxodrome.na(previous_point@coords[,2], current_point@coords[,2],
                                          previous_point@coords[,1], current_point@coords[,1])
        
        #find the gepgraphic location of each alternative point; calculate bearing to the next point: add ta to the bearing of the previous point
        current_point_m <- spTransform(current_point, meters_proj) #convert to meters proj
        rnd <- data.frame(lon = current_point_m@coords[,1] + rsl*cos(rta),lat = current_point_m@coords[,2] + rsl*sin(rta)) #for this to work, lat and lon should be in meters as well. boo. coordinates in meters?
        
        #covnert back to lat-lon proj
        rnd_sp <- rnd
        coordinates(rnd_sp) <- ~lon+lat
        proj4string(rnd_sp) <- meters_proj
        rnd_sp <- spTransform(rnd_sp, wgs)
        
        #put used and available points together
        df <- used_point@data %>%  
          slice(rep(row_number(),151)) %>% #paste each row 150 times for the used and alternative steps
          mutate(x = c(head(x,1),rnd_sp@coords[,1]),
                 y = c(head(y,1),rnd_sp@coords[,2]),
                 used = c(1,rep(0,150)))  %>%
          rowwise() %>% 
          mutate(heading = NCEP.loxodrome.na(lat1 = current_point$y, lat2 = y, lon1 = current_point$x, lon2 = x)) %>% 
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

Sys.time() - b #1.368372 mins

stopCluster(mycl) 
  
  
save(used_av_ls_90_60, file = "2021/ssf_input_all__90_60_150.RData")

#remove alternative points over land, and randomly select 50 points from those that remain

# points_water <- lapply(used_av_ls_1hr, function(x){
#   
#   pts <- x %>% 
#     st_as_sf(coords = c("x","y"), crs = wgs) %>% 
#     #st_intersection(ocean)
#     st_difference(land_no_buffer)
#   
# })


# ---------- STEP 2: annotate#####

load("2021/ssf_input_all_90_60_150.RData") #used_av_ls_90_60

#create one dataframe with movebank specs
used_av_df_90_60 <- lapply(c(1:length(used_av_ls_90_60)), function(i){
  
  data <- used_av_ls_90_60[[i]] %>% 
    dplyr::select( c("date_time", "x", "y", "selected",  "burst_id", "step_length", "turning_angle", "track", "step_id", "used", "heading")) %>% #later, add a unique step id: paste track, burst_id and step_id. lol
    mutate(timestamp = paste(as.character(date_time),"000",sep = "."),
           group = names(used_av_ls_90_60)[[i]]) %>% 
    rowwise() %>% 
    mutate(ind = strsplit(track, "_")[[1]][1],
           stratum = paste(track, burst_id, step_id, sep = "_")) %>% 
    as.data.frame()
}) %>% 
  reduce(rbind)

save(used_av_df_90_60, file = "2021/ssf_input_all_df_90_60_150.RData")


# #have a look
# X11();par(mfrow= c(2,1), mar = c(0,0,0,0), oma = c(0,0,0,0))
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_df_1hr[used_av_df_1hr$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_df_1hr[used_av_df_1hr$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_df_1hr2[used_av_df_1hr2$used == 0,c("x","y")], pch = 16, cex = 0.2, col = "gray55")
# points(used_av_df_1hr2[used_av_df_1hr2$used == 1,c("x","y")], pch = 16, cex = 0.2, col = "orange")
# 
# 
# maps::map("world",fil = TRUE,col = "grey85", border=NA) 
# points(used_av_all_2hr[used_av_all_2hr$used == 0,c("x","y")], pch = 16, cex = 0.4, col = "gray55")
# points(used_av_all_2hr[used_av_all_2hr$used == 1,c("x","y")], pch = 16, cex = 0.4, col = "orange")

#rename columns
colnames(used_av_df_90_60)[c(2,3)] <- c("location-long","location-lat")

write.csv(used_av_df_90_60, "2021/ssf_input_df_90_60_150.csv")

# summary stats
used_av_df_1hr %>% 
  group_by(species) %>% 
  #group_by(group) %>% 
  summarise(yrs_min = min(year(date_time)),
            yrs_max = max(year(date_time)),
            n_ind = n_distinct(ind),
            n_tracks = n_distinct(track))

#visual inspection
data_sf <- used_av_df_1hr %>% 
  filter(used == 1) %>% 
  st_as_sf(coords = c(3,4), crs = wgs)

mapview(data_sf, zcol = "species", viewer.suppress = TRUE)

# --- after movebank
#open annotated data and add wind support and crosswind
ann <- read.csv("2021/annotations/ssf_input_df_90_60_150.csv-5816495146751889418/ssf_input_df_90_60_150.csv-5816495146751889418.csv",
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


save(ann_50, file = "2021/ssf_input_annotated_90_60_50.RData")

#annotate with heat flux and blh, to calculate w*... could have been done with the previous annotation round.
df_for_ann <- ann_50 %>% 
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) 

colnames(df_for_ann)[c(26,27)] <- c("location-long","location-lat") #rename columns to match movebank format

write.csv(df_for_ann, "2021/df_for_ann_all_spp_90_60.csv")

# long-term annotation data (40 year data)
#prep a dataframe with 40 rows corresponding to 40 years (1981,2020), for each point. then i can calculate variance of delta t over 41 years for each point
df_40 <- ann_50 %>% #make sure this is not grouped!
  dplyr::select(-c(v925,u925,t2m,sst,delta_t)) %>% 
  slice(rep(row_number(),40)) %>% 
  group_by(row_id) %>% 
  mutate(year = c(1981:2020)) %>%
  ungroup() %>%
  mutate(timestamp = paste(as.character(date_time),"000",sep = ".")) %>% 
  as.data.frame()

str_sub(df_40$timestamp,1,4) <- df_40$year #replace original year with years from 1979-2019
colnames(df_40)[c(21,22)] <- c("location-long","location-lat") #rename columns to match movebank format

#break up into two parts. over 1 million rows
df_40_1 <- df_40 %>% 
  slice(1:800000)
write.csv(df_40_1, "2021/ssf_40_all_spp_90_60_1.csv")

df_40_2 <- df_40 %>% 
  slice(800001:1600000)
write.csv(df_40_2, "2021/ssf_40_all_spp_90_60_2.csv")

df_40_3 <- df_40 %>% 
  slice(1600001 : nrow(df_40))
write.csv(df_40_3, "2021/ssf_40_all_spp_90_60_3.csv")

#---- after movebank

load("2021/ssf_input_annotated_90_60_50.RData") #ann_50
ann_50 <- ann_50 %>% 
  mutate(abs_cross_wind = abs(cross_wind(u = u925, v = v925, heading = heading)))

#calculate long-term metrics and merge with previously annotated data
ann_40_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/annotations/40_yrs_90_60/",pattern = ".csv", full.names = T) 

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
  summarise_at(c("delta_t", "wind_speed", "wind_support", "abs_cross_wind", "u925", "v925"), #before calculating these, investigate why we have NAs??
               list(avg = ~mean(., na.rm = T), var = ~var(., na.rm = T), rsd = ~rsd(.))) %>% 
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

save(ann_cmpl, file = "2021/ssf_input_ann_cmpl_90_60.RData")


# ---------- STEP 4: data exploration#####

# --- 4.1: plot variances against latitude

#for each stratum, calculate variances and means across all points, not only the available points. then create a scatter plot of 
var_data <- ann_cmpl %>% 
  group_by(stratum) %>% 
  arrange(desc(used), .by_group = T) %>% 
  summarise(m_wspt = mean(wind_support),
            var_wspt = var(wind_support),
            m_delta_t = mean(delta_t),
            var_delta_t = var(delta_t),
            m_wspd = mean(wind_speed),
            var_wspd = var(wind_speed),
            used_lat = head(location.lat,1),
            species = head(species,1),
            zone = head(zone,1)) #the first row is the used point
  
#plot
par(mfrow = c(3,1))
plot(var_data$used_lat, var_data$var_wspt, main = "variation in wind support", col = as.factor(var_data$species))
plot(var_data$used_lat, var_data$var_wspd, main = "variation in wind speed", col = as.factor(var_data$species))
plot(var_data$used_lat, var_data$var_delta_t, main = "variation in delta t", col = as.factor(var_data$species))

par(mfrow = c(3,1))
plot(var_data$used_lat, var_data$m_wspt, main = "average wind support")
plot(var_data$used_lat, var_data$m_wspd, main = "average wind speed")
plot(var_data$used_lat, var_data$m_delta_t, main = "average delta t")


# --- 4.2

load("2021/ssf_input_ann_cmpl_90_60.RData") #ann_cmpl

#for plotting
ann_cmpl$species <- factor(ann_cmpl$species)


#plot
X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
for(i in c("wind_speed_avg", "abs_cross_wind_avg","delta_t_avg", "wind_support_avg")){
    
    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
    if(i == "wind_speed_avg"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  } 

mtext("40-yr averages at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
for(i in c("wind_speed_var", "abs_cross_wind_var","delta_t_var", "wind_support_var")){

  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  if(i == "wind_speed_var"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
} 

mtext("40-yr variances at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
for(i in c("wind_speed_rsd", "abs_cross_wind_rsd","delta_t_rsd", "wind_support_rsd")){
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  if(i == "wind_speed_rsd"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
} 

mtext("40-yr relative standard deviation (%) at each point", side = 3, outer = T, cex = 1.3)

X11(width = 15, height = 10);par(mfrow= c(3,3), oma = c(0,0,3,0))
for(i in c("wind_speed", "abs_cross_wind","delta_t", "wind_support")){
  
  for(j in unique(ann_cmpl$zone)){ 
    
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = paste(i,"(",j,")",sep = " "), xlab = "", ylab = "")
  if(i == "wind_speed"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$zone == j, i] ~ ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$zone == j,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$zone == j, "species"])) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$zone == j, i] ~ ann_cmpl[ann_cmpl$used == 0 & ann_cmpl$zone == j, "species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 & ann_cmpl$zone == j, "species"])) + 0.15)
} 
}
mtext("values at timestamp of each point", side = 3, outer = T, cex = 1.3)


X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
for(i in c("wind_speed","delta_t", "wind_support", "abs_cross_wind")){

    boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
    if(i == "wind_speed"){
      legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
    }
    boxplot(ann_cmpl[ann_cmpl$used == 1, i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], 
            xaxt = "n", add = T, boxfill = "orange",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1, "species"])) - 0.15)
    boxplot(ann_cmpl[ann_cmpl$used == 0, i] ~ ann_cmpl[ann_cmpl$used == 0, "species"], 
            xaxt = "n", add = T, boxfill = "grey",
            boxwex = 0.25, at = 1:length(unique(ann_cmpl[ann_cmpl$used == 1 , "species"])) + 0.15)

}
mtext("values at timestamp of each point", side = 3, outer = T, cex = 1.3)


######## just looking at the temperate and the tradewind zones

old_zones <- ann_cmpl[ann_cmpl$zone %in% c("tradewind","temperate"),]


X11(width = 15, height = 10);par(mfrow= c(2,2), oma = c(0,0,3,0))
for(i in c("wind_speed", "abs_cross_wind","delta_t", "wind_support")){
  
  boxplot(old_zones[,i] ~ old_zones[,"species"], data = old_zones, boxfill = NA, border = NA, main = i, xlab = "", ylab = "")
  if(i == "wind_speed_rsd"){
    legend("bottomleft", legend = c("used","available"), fill = c("orange","gray"), bty = "n")
  }
  boxplot(old_zones[old_zones$used == 1,i] ~ old_zones[old_zones$used == 1,"species"], 
          xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(old_zones$species)) - 0.15)
  boxplot(old_zones[old_zones$used == 0,i] ~ old_zones[old_zones$used == 0,"species"], 
          xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(old_zones$species)) + 0.15)
} 

mtext("values at timestamp of each point", side = 3, outer = T, cex = 1.3)

##############################

#plot relationship between wind and delta t
X11(); par(mfrow=c(3,2))
plot(delta_t ~ wind_support, data = ann_cmpl[ann_cmpl$used ==1,], col= factor(ann_cmpl[ann_cmpl$used ==1,"species"]), pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed, data = ann_cmpl[ann_cmpl$used ==1,], col= factor(ann_cmpl[ann_cmpl$used ==1,"species"]), pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind, data = ann_cmpl[ann_cmpl$used ==1,], col= factor(ann_cmpl[ann_cmpl$used ==1,"species"]), pch = 16, cex = 0.7)

plot(delta_t_var ~ wind_support_var, data = ann_cmpl[ann_cmpl$used ==1,], col= ann_cmpl[ann_cmpl$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t_var ~ wind_speed_var, data = ann_cmpl[ann_cmpl$used ==1,], col= ann_cmpl[ann_cmpl$used ==1,"species"], pch = 16, cex = 0.7)
plot(delta_t_var ~ abs_cross_wind_var, data = ann_cmpl[ann_cmpl$used ==1,], col= ann_cmpl[ann_cmpl$used ==1,"species"], pch = 16, cex = 0.7)

X11(); par(mfrow=c(3,1))
plot(delta_t ~ wind_support, data = ann_cmpl, col= factor(ann_cmpl$used), pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed, data = ann_cmpl, col= factor(ann_cmpl$used), pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind, data = ann_cmpl, col = factor(ann_cmpl$used), pch = 16, cex = 0.7)

X11(); par(mfrow=c(3,1))
plot(delta_t ~ wind_support_var, data = ann_cmpl, col=factor(ann_cmpl$used), pch = 16, cex = 0.7)
plot(delta_t ~ wind_speed_var, data = ann_cmpl, col= factor(ann_cmpl$used), pch = 16, cex = 0.7)
plot(delta_t ~ abs_cross_wind_var, data = ann_cmpl, col= factor(ann_cmpl$used), pch = 16, cex = 0.7)


#### exploration. run a quick clogit to see if the results are what I expect

library(survival)
form_original <- formula(used ~  delta_t_z * wind_speed_z + wind_support_z +
                    strata(stratum))

form_1 <- formula(used ~ zone * delta_t_z * wind_speed_z + zone * delta_t_z *  wind_support_z +
                      strata(stratum))

form_2 <- formula(used ~ zone * delta_t_z +  zone * wind_support_z +
                    strata(stratum))


m1 <- clogit(form_1, data = all_data) #when zone is added, delta t coeff becomes positive. still not sig, but positive. but wind support is negative. but interaction of wind support with zones is positive. so, only negative in the arctic



all_old <- ann_old %>% 
  #group_by(species) 
  mutate_at(c(2:5,8:11,36:39),
            list(z = ~scale(.))) %>%
  as.data.frame() 

m2 <- clogit(form_1, data = all_old)


#look at only the two old zones in new data
all_2_zones <- old_zones %>% 
  #group_by(species) 
  mutate_at(c(2:5,8:11,40:43),
            list(z = ~scale(.))) %>%
  as.data.frame() 

m3 <- clogit(form_1, data = all_2_zones)



#-------
#extract number of strata per zone and consider removing the arctic

all_data %>% 
  group_by (zone) %>% 
  summarise(n_s = n_distinct(stratum),
            n_t = n_distinct(track),
            n_i = n_distinct(ind)) #the arctic only has 26 strata from 4 tracks of 3 individuals.... let's remove it!

data_no_arctic <- all_data %>% 
  filter(zone != "arctic")

m4 <- clogit(form_original, data = data_no_arctic) #coefficients are not too different from using all the data. 


#----------------------------------------------------------------------
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
  mutate_at(c(2:5,8:11,40:43),
            list(z = ~scale(.))) %>%
  as.data.frame() 

save(all_data, file = "ssf_input_ann_1hr_z.RData")

# STEP 6: modeling#####
load("ssf_input_ann_1hr_z.RData") #all_data

#create a var for individual id
all_data <- all_data %>% 
  rowwise() %>% 
  mutate(ind = strsplit(track, "_")[[1]][1]) %>% 
  as.data.frame()

#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = species,
         species2 = species,
         species3 = species,
         species4 = species,
         species5 = species,
         species6 = species,
         ind1 = factor(ind),
         ind2 = factor(ind),
         ind3 = factor(ind),
         ind4 = factor(ind),
         ind5 = factor(ind),
         ind6 = factor(ind),
         lat_zone1 = factor(lat_zone),
         lat_zone2 = factor(lat_zone),
         lat_zone3 = factor(lat_zone),
         lat_zone4 = factor(lat_zone),
         lat_zone5 = factor(lat_zone),
         lat_zone6 = factor(lat_zone),
         stratum = factor(stratum)) #%>% 
  dplyr::select(c("used","delta_t_z","wind_speed_z","wind_support_z","wind_support_var_z","delta_t_var_z","species1","species2", "species3", "ind1", "ind2", "ind3", "stratum"))

#save file for cluster computing
save(all_data, file = "/home/enourani/ownCloud/Work/cluster_computing/global_seascape_proj/inla_ssf/inla_ssf.RData")

### trying out predicting with INLA (use model m1d)
# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

#create a sample for testing model formulas in short time
str_to_keep <- sample(unique(all_data$stratum),150)
sample <- all_data[all_data$stratum %in% str_to_keep,] %>% 
  dplyr::select(c("used","delta_t_z","wind_speed_z","wind_support_z", "species1","species2", "species3", "ind1", "ind2", "ind3", "stratum"))

#add missing data to the sample for prediction (purpose is to make interaction plots. also see inla_plots.R)

n <- 200
new <- data.frame(used = rep(NA,n),
                       delta_t_z = as.numeric(sample(c(-3:3), n, replace = T)),
                       wind_speed_z = as.numeric(sample(c(-3:3), n, replace = T)), #keep wind speed values at -2,0 and 2
                       wind_support_z = rep(0,n), #keep wind support values at 0
                       stratum = factor(sample(c(1:3), n, replace = T)),
                       species1 = sample(all_data$species1, n, replace =T),
                       ind1 = sample(all_data$ind1, n, replace = T)) %>% 
  mutate(species2 = species1,
         species3 = species1,
         ind2 = ind1,
         ind3 = ind1) 
new_data <- new %>% 
  full_join(sample)

formula1d <- used ~ -1 + wind_speed_z * delta_t_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1d_sample <- inla(formula1d, family = "Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = new_data, #use the sample dataset 
            num.threads = 10,
            #control.family = list(link='test1'), 
            control.predictor = list(link = 1,compute=TRUE), #this means that NA values will be predicted. link can also be set. but will make the predictions Inf (response is binary but family is poisson.. what is the correct link!!??) # “apply the first link function to everything”.
            control.compute = list(openmp.strategy="huge", config = TRUE))#, mlik = T, waic = T)) 
Sys.time() - b
#setting the link to 1 produces NaNs and not setting it leads to numbers of a strange range (negatives)..hmmm

#extract predicted values
used_na <- which(is.na(new_data$used))
m1d_sample$summary.fitted.values[used_na,]

#plot interaction effects
ws_low <- which(is.na(new_data$used) & new_data$wind_speed_z == -1)
ws_high <- which(is.na(new_data$used) & new_data$wind_speed_z == 1)
ws_zero <- which(is.na(new_data$used) & new_data$wind_speed_z == 0)

X11();par(mfrow = c(1,3))
plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = -1" ,
      xlab="delta_t" )
points(new_data[ws_low,"delta_t_z"], m1d_sample$summary.fitted.values[ws_low,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 0" ,
      xlab="delta_t" )
points(new_data[ws_zero,"delta_t_z"], m1d_sample$summary.fitted.values[ws_zero,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_sample$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 1" ,
      xlab="delta_t" )
points(new_data[ws_high,"delta_t_z"], m1d_sample$summary.fitted.values[ws_high,"mean"])


#create a raster of predictions
preds <- data.frame(delta_t = new_data[is.na(new_data$used) ,"delta_t_z"],
                    wind_speed = new_data[is.na(new_data$used) ,"wind_speed_z"],
                    pres = m1d_sample$summary.fitted.values[used_na,"mean"])

avg_preds <- preds %>% 
  group_by(delta_t,wind_speed) %>% 
  summarise(avg_pres = mean(pres)) %>% 
  as.data.frame()
  
coordinates(avg_preds) <-~ delta_t + wind_speed
gridded(avg_preds) <- TRUE
plot(raster(avg_preds),ylab = "wind speed (z)", xlab = "delta_t (z)")

#try a 2-D plot with plotly
fig <-plot_ly(x = new_data[is.na(new_data$used) ,"delta_t_z"], 
              y = new_data[is.na(new_data$used) ,"wind_speed_z"])


dt.seq <- -1:1
for ( wsp in -1:1 ) {
  dt <- all_data[all_data$wind_speed_z == wsp,]
  plot( used ~ delta_t_z , data=dt , col=rangi2 ,
        main=paste("wspeed =",w) , xaxp=c(-1,1,2) , ylim=c(0,362) ,
        xlab="delta_t (centered)" )
  mu <- link( mc1 , data=data.frame(water.c=w,shade.c=shade.seq) ) #dont know the equivalent of this in inla
  mu.mean <- apply( mu , 2 , mean )
  mu.PI <- apply( mu , 2 , PI , prob=0.97 )
  lines( -3:3 , m1d_sample$summary.fitted.values[used_na,"mean"] )
  lines( shade.seq , mu.PI[1,] , lty=2 )
  lines( shade.seq , mu.PI[2,] , lty=2 )
}

#instead of making predictions, we can also sample from the posterior distribution. (following Virgilio's book)
#make sure config = TRUE in the inla call
samp1 <- inla.posterior.sample(100, m1d_sample, selection = list(wind_speed_z= 1, delta_t_z = 1)) #selection can be used to determine which variables we want
#1 in the above call means that we want to keep the first element in effect 
samp2 <- inla.posterior.sample.eval(function(...) {wind_speed_z * delta_t_z}, samp1)
summary(as.vector(samp2))


##try on complete dataset
########################
#model with wind speed instead of ws and cw...... delta_t var * wspeed var produces loads of warnings!
formula1c <- used ~ -1 + delta_t_z * wind_speed_z + delta_t_var_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species6, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind5, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind6, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1c <- inla(formula1c, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.predictor = list(link = 1), #required to set the right link (i.e., the logit function) 
                                                #to have the fitted values in the appropriate scale (i.e., the expit of the linear predictor).
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T, cpo = T))
Sys.time() - b

save(m1c, file = "inla_models/m1c.RData")

load("inla_models/m1c.RData")
summary(m1c)

#plot the posterior densities
bri.hyperpar.plot(m1c) #summary of hyperparameters in SD scale (converts precision to sd)
Efxplot(m1c) + theme_bw()

#remove var of delta t
  formula1e <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species6, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind6, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1e <- inla(formula1e, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T, cpo = T))
Sys.time() - b

save(m1e, file = "inla_models/m1e.RData")

load("inla_models/m1e.RData")
summary(m1c)

#plot the posterior densities
bri.hyperpar.plot(m1e) #summary of hyperparameters in SD scale (converts precision to sd)
Efxplot(list(m1c,m1e,m1d)) + theme_bw()


#remove both vars
formula1d <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1d <- inla(formula1d, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = F, cpo = F))
Sys.time() - b

save(m1d, file = "inla_models/m1d.RData")

load("inla_models/m1d.RData")
summary(m1d)

#plot the posterior densities
bri.hyperpar.plot(m1d) #summary of hyperparameters in SD scale (converts precision to sd)
Efxplot(list(m1c,m1e,m1d)) + theme_bw()



#######
#model with wind speed and wind support
#also predict with this model
new_data <- all_data %>% 
  full_join(new)

formula1d <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1d_pred <- inla(formula1d, family ="Poisson", 
                 control.fixed = list(
                   mean = mean.beta,
                   prec = list(default = prec.beta)),
                 data = all_data,
                 num.threads = 10,
                 control.predictor = list(compute=TRUE),
                 control.compute = list(openmp.strategy="huge", config = TRUE, mlik = T, waic = T, cpo = T)) 
Sys.time() - b

save(m1d_pred, file = "inla_models/m1d_pred_waic.RData")

save(m1d_pred, file = "inla_models/m1d_pred.RData")

load("inla_models/m1d.RData")
summary(m1d)

#plot the posterior densities 
bri.hyperpar.plot(m1d) #summary of hyperparameters in SD scale (converts precision to sd)
Efxplot(m1d) + theme_bw()
Efxplot(list(m1c,m1d)) + theme_bw()

#extract predicted values
used_na <- which(is.na(new_data$used))
m1d_pred$summary.fitted.values[used_na,]


#remove species. plots showed that species are not very different
formula1f <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
m1f <- inla(formula1f, family ="Poisson", 
            control.fixed = list(
              mean = mean.beta,
              prec = list(default = prec.beta)),
            data = all_data,
            num.threads = 10,
            control.predictor = list(compute=TRUE),
            control.compute = list(openmp.strategy="huge", config = TRUE))#, mlik = T, waic = T, cpo = T)) 
Sys.time() - b


Efxplot(list(m1c,m1e, m1d, m1f)) + theme_bw()  
  
  
################################
#plot interaction effects
ws_low <- which(is.na(new_data$used) & new_data$wind_speed_z == -1)
ws_high <- which(is.na(new_data$used) & new_data$wind_speed_z == 1)
ws_zero <- which(is.na(new_data$used) & new_data$wind_speed_z == 0)

X11();par(mfrow = c(1,3))
plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_pred$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = -1" ,
      xlab="delta_t" )
points(new_data[ws_low,"delta_t_z"], m1d_pred$summary.fitted.values[ws_low,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_pred$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 0" ,
      xlab="delta_t" )
points(new_data[ws_zero,"delta_t_z"], m1d_pred$summary.fitted.values[ws_zero,"mean"])

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], m1d_pred$summary.fitted.values[used_na,"mean"],
      type="n",
      main="wind speed = 1" ,
      xlab="delta_t" )
points(new_data[ws_high,"delta_t_z"], m1d_pred$summary.fitted.values[ws_high,"mean"])



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

