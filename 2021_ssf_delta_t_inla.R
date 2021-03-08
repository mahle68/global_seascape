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


# --- 4.2: boxplots

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


# --- 4.3:plot relationship between wind and delta t

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


# ---------- STEP 5: clogit to get a feel for things #####

load("2021/ssf_input_ann_90_60_z.RData") #all_data

#### exploration. run a quick clogit to see if the results are what I expect

form_original <- formula(used ~  delta_t_z * wind_speed_z + wind_support_z +
                    strata(stratum))

form_6 <- formula(used ~  delta_t_z * wind_support_z + 
                           strata(stratum))

form_1 <- formula(used ~ lat_at_used * delta_t_z + lat_at_used * wind_support_z +
                      strata(stratum))

form_2 <- formula(used ~ lat_at_used * delta_t_z +  lat_at_used * wind_support_z +
                    strata(stratum))

form_3 <- formula(used ~  delta_t_z * wind_speed_z + delta_t_z * wind_support_z +
                           strata(stratum))

form_4 <- formula(used ~  delta_t_z * wind_support_z +
                    strata(stratum))

form_5 <- formula(used ~  delta_t_z * lat_at_used + wind_support_z * lat_at_used +
                    strata(stratum))

m1 <- clogit(form_1, data = all_data) #when zone is added, delta t coeff becomes positive. still not sig, but positive. but wind support is negative. but interaction of wind support with zones is positive. so, only negative in the arctic
m1a <- clogit(form_original, data = all_data)
m1b <- clogit(form_2, data = all_data)
m1c <- clogit(form_3, data= all_data)
m1d <- clogit(form_4, data= all_data)
m1e <- clogit(form_5, data = all_data)
m1f <- clogit(form_6, data = all_data)


###zone_specific... conclusion: the order of importance is pretty much the same. only the direction of delta_t:windspeed changes... but still small. delta t is negative all through (not sig. in tw)

twtmp <- clogit(form_original, data = all_data[all_data$lat_at_used %in% c("temperate", "tradewind"),]) #similar to all_data
tmp <- clogit(form_original, data = all_data[all_data$lat_at_used == "temperate",])
tw <- clogit(form_original, data = all_data[all_data$lat_at_used == "tradewind",])
arctr <- clogit(form_original, data = all_data[all_data$lat_at_used %in% c("tropical", "arctic"),])

#---
#extract number of strata per zone and consider removing the arctic

all_data %>% 
  group_by (zone) %>% 
  summarise(n_s = n_distinct(stratum),
            n_t = n_distinct(track),
            n_i = n_distinct(ind)) #the arctic only has 26 strata from 4 tracks of 3 individuals.... let's remove it!


# ---------- STEP 6: INLA #####

load("2021/ssf_input_ann_cmpl_90_60.RData") #ann_cmpl

#correlation
ann_cmpl %>% 
  dplyr::select(c(2:5,8:11,38:41,45,46)) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat. avg delta_t and delta_t. avg_ws and var_delta_t
#correlated: wind support var & wind speed var and cross wind var, crosswind var and wind speed var.  

#z-transform
all_data <- ann_cmpl %>% 
  #group_by(species) 
  mutate_at(c(2:5,8:11,38:41,45,46),
            list(z = ~scale(.))) %>%
  arrange(stratum, desc(used)) %>% 
  group_by(stratum) %>%  
  mutate(lat_at_used = head(zone,1)) %>%  #add a variable for latitudinal zone. This will assign the lat zone of the used point to the entire stratum
  ungroup() %>% 
  as.data.frame() 

save(all_data, file = "2021/ssf_input_ann_90_60_z.RData")


#check to make sure each stratum has one lat zone value
all_data %>% 
  group_by(stratum) %>% 
  summarize(n = n_distinct(lat_at_used),
            n_z = n_distinct(zone)) %>% 
  filter(n > 1)


load("2021/ssf_input_ann_90_60_z.RData") #all_data

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
         zone1 = factor(lat_at_used),
         zone2 = factor(lat_at_used),
         zone3 = factor(lat_at_used),
         zone4 = factor(lat_at_used),
         zone5 = factor(lat_at_used),
         zone6 = factor(lat_at_used),
         stratum = factor(stratum)) %>% 
  dplyr::select(c("used","stratum","delta_t_z","wind_speed_z","wind_support_z","wind_support_var_z", "abs_cross_wind_z","delta_t_var_z",
                  "species1","species2", "species3", "species4","species5","ind1", "ind2", "ind3", "ind4", "ind5", "zone1", "zone2",
                  "zone3","zone4","zone5","zone6","location.lat"))


# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

formulaM1 <- used ~ -1 + delta_t_z * wind_speed_z + delta_t_var_z + wind_support_z + wind_support_var_z +
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
  f(species5, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind5, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

(b <- Sys.time())
M1 <- inla(formulaM1, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10, #This depends on your computer
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T))

Sys.time() - b #3.3 hours

save(M1, file = "2021/inla_models/m1.RData")

load("2021/inla_models/m1.RData")
summary(M1)



#remove variance of delta t
formulaM2 <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(species1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(species2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind1, delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ind2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(ind4, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))


M2 <- inla(formulaM2, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T, cpo = T))

#summary(M2)
save(M2, file = "2021/inla_models/m2.RData")

load("2021/inla_models/m2.RData")

#remove variance of wind support
formulaM3 <- used ~ -1 + delta_t_z * wind_speed_z + wind_support_z +
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


M3 <- inla(formulaM3, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10,
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = F, cpo = F),
           control.predictor = list(compute=TRUE)) #to be able to calculate linear combo to then plot the interaction term)

save(M3, file = "2021/inla_models/m3.RData")
#summary(M3)

load("2021/inla_models/m3.RData")
  

# latitudinal zone as a random effect (on the full model. instead of species and individual for now)
formulaM0 <- used ~ -1 + delta_t_z * wind_speed_z + delta_t_var_z + wind_support_z + wind_support_var_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + 
  f(zone1, delta_t_z, model = "iid", 
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(zone2, wind_speed_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(zone3, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(zone4, delta_t_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(zone5, wind_support_var_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))
  
(b <- Sys.time())
M0 <- inla(formulaM0, family ="Poisson", 
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
           data = all_data,
           num.threads = 10, #This depends on your computer
           control.compute = list(openmp.strategy = "huge", config = TRUE, mlik = T, waic = T))

Sys.time() - b #3.3 hours

save(M0, file = "2021/inla_models/m0.RData")


# ---------- STEP 7: plots #####
#FIGURE 2: posterior means of fixed effects ####

ModelList <- list(M1,M2,M3)
graphlist<-list()
for(i in 1:length(ModelList)){
  model<-ModelList[[i]]
  
  graph<-as.data.frame(summary(model)$fixed)
  colnames(graph)[which(colnames(graph)%in%c("0.025quant","0.975quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("0.05quant","0.95quant"))]<-c("Lower","Upper")
  colnames(graph)[which(colnames(graph)%in%c("mean"))]<-c("Estimate")
  
  graph$Model<-i
  graph$Factor<-rownames(graph)
  
  graphlist[[i]]<-graph
}

graph <- bind_rows(graphlist)

graph$Sig <- with(graph, ifelse(Lower*Upper>0, "*", ""))

graph$Model <- as.factor(graph$Model)

position <- ifelse(length(unique(graph$Model))  ==  1, "none", "right")

VarOrder <- rev(unique(graph$Factor))
VarNames <- VarOrder

graph$Factor <- factor(graph$Factor, levels = VarOrder)
levels(graph$Factor) <- VarNames

min<-min(graph$Lower,na.rm = T)
max<-max(graph$Upper,na.rm = T)

graph$Factor_n <- as.numeric(graph$Factor)

X11(width = 3.5, height = 3)

par(mfrow=c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)

#plot(0, type = "n", bty = "l",labels = FALSE, tck = 0, ann = F, xlim = c(-6,8), ylim = c(0,6.3))
plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-6,8), ylim = c(0,6.3), xlab = "Estimate", ylab = "")

#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)
#add points and error bars
points(graph[graph$Model == 1, "Estimate"], graph[graph$Model == 1,"Factor_n"] - 0.2, col = "salmon2", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 1, "Lower"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       graph[graph$Model == 1, "Upper"], graph[graph$Model == 1,"Factor_n"] - 0.2,
       col = "salmon2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(graph[graph$Model == 2, c("Estimate","Factor")], col = "palegreen3", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 2, "Lower"], graph[graph$Model == 2,"Factor_n"],
       graph[graph$Model == 2, "Upper"], graph[graph$Model == 2,"Factor_n"],
       col = "palegreen3", code = 3, length = 0.03, angle = 90)

points(graph[graph$Model == 3, "Estimate"], graph[graph$Model == 3,"Factor_n"] + 0.2, col = "steelblue1", pch = 19, cex = 1.3)
arrows(graph[graph$Model == 3, "Lower"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       graph[graph$Model == 3, "Upper"], graph[graph$Model == 3,"Factor_n"] + 0.2,
       col = "steelblue1", code = 3, length = 0.03, angle = 90)
#add axes
axis(side= 1, at= c(-5,0,5), labels= c("-5", "0", "5"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:6), #line=-4.8, 
     labels= c( expression(paste(Delta,"t"," : wind speed")),
                "wind support var","wind support",#expression(paste(Delta,"t"," var")),
                "wind speed", expression(paste(Delta,"t")),NA),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2 ) # text perpendicular to axis label 

#add legend
legend(x = 5.3, y = 0.8, legend=c("model 3", "model 2", "model 1"), col = c("steelblue1","palegreen3","salmon2"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.75)







#SUPPLEMENTARY FIGURE 1: species-specific coefficients ####
#for the best model (M3); original code by Virgilio Gomez-Rubio (Bayesian inference with INLA, 2020)
species_names <- unique(all_data$species)

tab_dt <- data.frame(ID = as.factor(M3$summary.random$species1$ID),
                     mean = M3$summary.random$species1$mean,
                     IClower = M3$summary.random$species1[, 4],
                     ICupper = M3$summary.random$species1[, 6])

tab_wspd <- data.frame(ID = as.factor(M3$summary.random$species2$ID),
                       mean = M3$summary.random$species2$mean,
                       IClower = M3$summary.random$species2[, 4],
                       ICupper = M3$summary.random$species2[, 6])

tab_wspt <- data.frame(ID = as.factor(M3$summary.random$species3$ID),
                       mean = M3$summary.random$species3$mean,
                       IClower = M3$summary.random$species3[, 4],
                       ICupper = M3$summary.random$species3[, 6])

X11(width = 3.5, height = 3)

par(mfrow = c(1,1), bty="n", #no box around the plot
    #cex.axis= 0.75, #x and y labels have 0.75% of the default size
    #font.axis= 0.75, #3: axis labels are in italics
    #cex.lab = 0.75,
    cex = 0.7,
    oma = c(0,3.5,0,0),
    mar = c(3, 3.5, 0.5, 1),
    bty = "l"
)


plot(0, type = "n", labels = FALSE, tck = 0, xlim = c(-4,4), ylim = c(0,4.3), xlab = "", ylab = "")
#add vertical line for zero
abline(v = 0, col = "grey30",lty = 2)

points(tab_dt$mean, as.numeric(tab_dt$ID) - 0.2, col = "lightcoral", pch = 19, cex = 1.3)
arrows(tab_dt$IClower, as.numeric(tab_dt$ID) - 0.2,
       tab_dt$ICupper, as.numeric(tab_dt$ID) - 0.2,
       col = "lightcoral", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspd$mean, as.numeric(tab_wspd$ID), col = "yellowgreen", pch = 19, cex = 1.3)
arrows(tab_wspd$IClower, as.numeric(tab_wspd$ID),
       tab_wspd$ICupper, as.numeric(tab_wspd$ID),
       col = "yellowgreen", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

points(tab_wspt$mean, as.numeric(tab_wspt$ID) + 0.2, col = "paleturquoise2", pch = 19, cex = 1.3)
arrows(tab_wspt$IClower, as.numeric(tab_wspt$ID) + 0.2,
       tab_wspt$ICupper, as.numeric(tab_wspt$ID) + 0.2,
       col = "paleturquoise2", code = 3, length = 0.03, angle = 90) #angle of 90 to make the arrow head as straight as a line

axis(side= 1, at= c(-2,0,2), labels= c("-2", "0", "2"), 
     tick=T ,col = NA, col.ticks = 1, tck=-.015)

axis(side= 2, at= c(1:4), #line=-4.8, 
     #labels= c("Falco eleonorae", "Pandion haliaetus", "Pernis ptilorhynchus", "Falco peregrinus"),
     labels= c(expression(italic("F. eleonorae")), expression(italic("P. haliaetus")),
               expression(italic("P. ptilorhynchus")),
               expression(italic("F. peregrinus"))),
     tick=T ,col = NA, col.ticks = 1, # NULL would mean to use the defult color specified by "fg" in par
     tck=-.015 , #tick marks smaller than default by this proportion
     las=2 ) # text perpendicular to axis label 
  
#add legend
legend(x = 1.8, y = 0.6, legend=c("wind support", "wind speed", expression(paste(Delta,"t"))), 
       col = c("paleturquoise2","yellowgreen","lightcoral"), #coords indicate top-left
       pch = 19, bg="white",bty="n", cex = 0.75)

##############plot the interaction term
#https://groups.google.com/g/r-inla-discussion-group/c/DwPyw_fOUGI
dt = inla.rmarginal(10000, M3$marginals.fixed$delta_t_z)
wspd = inla.rmarginal(10000, M3$marginals.fixed$wind_speed_z)
dt_wspd = inla.rmarginal(10000, M3$marginals.fixed$delta_t_z) + inla.rmarginal(10000, M3$marginals.fixed$`delta_t_z:wind_speed_z`)

X11()
plot(density(dt), lty=2, col = "blue")
lines(density(wspd), lty=2, col = "red")
lines(density(dt_wspd), lty=1)



#from inla_plots.R
#use a sample
str_to_keep <- sample(unique(all_data$stratum),150)
sample <- all_data[all_data$stratum %in% str_to_keep,]

#add missing data to the sample for prediction (purpose is to make interaction plots. also see inla_plots.R)

n <- 200
new <- data.frame(used = rep(NA,n),
                  delta_t_z = as.numeric(sample(c(-3:3), n, replace = T)),
                  wind_speed_z = as.numeric(sample(c(-3:3), n, replace = T)), #keep wind speed values at -2,0 and 2
                  wind_support_z = rep(0,n), #keep wind support values at 0
                  stratum = factor(sample(c(1:3), n, replace = T)),
                  species1 = sample(sample$species1, n, replace =T),
                  ind1 = sample(sample$ind1, n, replace = T)) %>% 
  mutate(species2 = species1,
         species3 = species1,
         ind2 = ind1,
         ind3 = ind1) 

new_data <- new %>% 
  full_join(sample)

#model 
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
                   control.predictor = list(compute=TRUE), #this means that NA values will be predicted. link can also be set. but will make the predictions Inf (response is binary but family is poisson.. what is the correct link!!??) # “link = 1: apply the first link function to everything”.
                   control.compute = list(openmp.strategy="huge", config = TRUE))#, mlik = T, waic = T)) 
Sys.time() - b

#extract predicted values
used_na <- which(is.na(new_data$used))
m1d_sample$summary.fitted.values[used_na,]

#plot interaction effects. this works. but the range of y values doesnt make sense, even when I take the exp()
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



X11();par(mfrow = c(1,3))
plot( new_data[which(is.na(new_data$used)),"delta_t_z"], exp(m1d_sample$summary.fitted.values[used_na,"mean"]),
      type="n",
      main="wind speed = -1" ,
      xlab="delta_t" )
points(new_data[ws_low,"delta_t_z"], exp(m1d_sample$summary.fitted.values[ws_low,"mean"]))

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], exp(m1d_sample$summary.fitted.values[used_na,"mean"]),
      type="n",
      main="wind speed = 0" ,
      xlab="delta_t" )
points(new_data[ws_zero,"delta_t_z"], exp(m1d_sample$summary.fitted.values[ws_zero,"mean"]))

plot( new_data[which(is.na(new_data$used)),"delta_t_z"], exp(m1d_sample$summary.fitted.values[used_na,"mean"]),
      type="n",
      main="wind speed = 1" ,
      xlab="delta_t" )
points(new_data[ws_high,"delta_t_z"], exp(m1d_sample$summary.fitted.values[ws_high,"mean"]))


#instead of three plots, try to do one




## try the linear combo


lc = inla.make.lincombs(
  dt = new_data$delta_t_z,
  "delta_t_z:wind_speed_z" = new_data$wind_speed_z * new_data$delta_t_z)

(b <- Sys.time())
m1d_sample <- inla(formula1d, family = "Poisson", 
                   control.fixed = list(
                     mean = mean.beta,
                     prec = list(default = prec.beta)),
                   data = new_data, #use the sample dataset 
                   num.threads = 10,
                   #control.family = list(link='test1'), 
                   control.predictor = list(compute=TRUE), lincomb = lc, #this means that NA values will be predicted. link can also be set. but will make the predictions Inf (response is binary but family is poisson.. what is the correct link!!??) # “link = 1: apply the first link function to everything”.
                   control.compute = list(openmp.strategy="huge", config = TRUE)) 
Sys.time() - b
#i get an error regarding the lc

