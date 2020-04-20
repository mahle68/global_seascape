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
library(corrr)
library(lme4)
library(MuMIn)
library(mgcv)
library(survival)
library(INLA)
library(ggregplot) #devtools::install_github("gfalbery/ggregplot")


setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")

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
  
save(ann_cmpl, file = "ssf_input_ann.RData")

#assign unique step-ids and species
load("ssf_input_ann.RData")

ann_cmpl <- ann_cmpl %>% 
  mutate(stratum = paste(track, seg_id, burst_id, step_id, sep = "_")) %>% 
  rowwise() %>% 
  mutate(species = strsplit(group, "_")[[1]][1],
         lat_zone = ifelse(location.lat > 30, "tmpz","twz")) %>% 
  as.data.frame()

# STEP 4: data exploration#####

#plot
X11();par(mfrow= c(3,1))
for(i in c("avg_ws","avg_cw","avg_delta_t")){
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i)
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  
}

X11();par(mfrow= c(3,1))
for(i in c("var_ws", "var_cw","var_delta_t")){
  
  boxplot(ann_cmpl[,i] ~ ann_cmpl[,"species"], data = ann_cmpl, boxfill = NA, border = NA, main = i)
  boxplot(ann_cmpl[ann_cmpl$used == 1,i] ~ ann_cmpl[ann_cmpl$used == 1,"species"], xaxt = "n", add = T, boxfill = "orange",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) - 0.15)
  boxplot(ann_cmpl[ann_cmpl$used == 0,i] ~ ann_cmpl[ann_cmpl$used == 0,"species"], xaxt = "n", add = T, boxfill = "grey",
          boxwex = 0.25, at = 1:length(unique(ann_cmpl$species)) + 0.15)
  
}

#correlation
ann_cmpl %>% 
  dplyr::select(c("var_ws","var_cw","var_delta_t","wind_support","cross_wind","delta_t","location.lat")) %>% 
  correlate() %>% 
  stretch() %>% 
  filter(abs(r) > 0.6) #correlated: var_cw with location.lat and var_delta_t with location.lat

#z-transform
all_data <- ann_cmpl %>% 
  #group_by(species) # z_transform for each species separately. or not? ... huh!
  mutate_at(c("var_ws","var_cw","var_delta_t","wind_support","cross_wind","delta_t"),
            list(z = ~scale(.))) %>%
  as.data.frame()

save(all_data, file = "ssf_input_ann_z.RData")

# STEP 5: modeling#####
load("ssf_input_ann_z.RData") #all_data

#############using glm
#instantaneous model
form_inst <- formula(used ~ lat_zone * delta_t_z + lat_zone * wind_support_z + lat_zone * cross_wind_z)
form_inst_m <- formula(used ~ lat_zone * delta_t_z + lat_zone * wind_support_z + lat_zone * cross_wind_z + (1 | stratum) + (1 |species))

m1_inst <- glmer(form_inst_m, family = binomial(link = "cloglog"), data = all_data) #did not converge

m2_inst <- gamm(used ~ s(wind_support) + s(cross_wind) + delta_t,
           random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

m3_inst <- gamm(used ~ s(wind_support_z) + s(cross_wind_z) + delta_t_z,
                random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #lower AIC than without scaling


#contemporaneous model
form_cnt <- formula(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * var_cw_z)
form_cnt_m <- formula(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * var_cw_z + (1 | stratum) + (1 |species))

m1_cnt <- glmer(form_cnt_m, family = binomial(link = "cloglog"), data = all_data)

m2_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t,
                random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

all_data$lat_zone_f <- factor(all_data$lat_zone)
all_data$species_f <- factor(all_data$species)
all_data$stratum_f <- factor(all_data$stratum)

m3_cnt <- gamm(used ~ s(var_ws, by = lat_zone_f) + s(var_cw, by = lat_zone_f) + var_delta_t * lat_zone_f,
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #singluar

m4_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t * lat_zone_f,
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data ) #interaction not significant

m5_cnt <- gamm(used ~ s(var_ws) + s(var_cw) + var_delta_t + f(lat_zone_f),
               random = list(stratum = ~1, species = ~1), family = binomial(link = "cloglog"), data = all_data )

#inspo from https://drmowinckels.io/blog/gamm-random-effects/
m6_cnt <- gamm(used ~ s(var_ws, bs = "cr") + s(var_cw, bs = "cr") + s(var_delta_t, bs = "cr") +
                 s(stratum_f, bs = "re") + s(species_f, bs = "re") + s(lat_zone_f, bs = "re"), #random effects can be assigned this way, instead of using the random argument
               family = binomial(link = "cloglog"), data = all_data , correlation=corAR1()) #what's the deal with the correlation?? sooo slow



form <- formula(used ~ lat_zone * delta_t_z + lat_zone * var_delta_t_z + lat_zone * wind_support_z + lat_zone * var_ws_z + 
                  lat_zone * cross_wind_z + lat_zone * var_cw_z )

 


############# gamm + poisson distribution (as suggested by Muff et al 2018)

#the weight implements different variances per species. the by commamd
# in the smoother ensures that we have one smoother for each bird species (p. 368)... but data format needs to be changed
lmc <- lmeControl(niterEM = 5000, msMaxIter = 1000)

form <- formula(used ~ delta_t_z + wind_support_z + (1|stratum_f) + (0 + delta_t_z | species) + (0 + wind_support_z | species))

TMBstr <- glmmTMB(form, family = poisson, data = all_data, doFit = F)
#fix the standard deviation of the first random term, which is the (1|stratum) component in the model equation
TMBstr$parameters$theta[1] <- log(1e3)
TMBstr$mapArg <- list(theta = factor(c(NA,1:2))) #I have no idea what is happening here. refer to Muff et al

mTMB <- glmmTMB:::fitTMB(TMBstr)

gamm(form, control = lmc, method = "REML", weights = varIdent(form = ~1 | species), family = poisson, data = all_data) #muff et al say dont use weigths for ssf

############# inla
##codes from muff et al supp material https://conservancy.umn.edu/bitstream/handle/11299/204737/Otters_SSF.r?sequence=22&isAllowed=y
##and https://ourcodingclub.github.io/tutorials/inla/

#my stuff
#repeat variabels that will be used as random slopes
all_data <- all_data %>% 
  mutate(species1 = factor(species),
         species2 = factor(species),
         species3 = factor(species),
         species4 = factor(species),
         species5 = factor(species),
         species6 = factor(species),
         stratum = factor(stratum),
         lat_zone = factor(lat_zone))

formula2 <- used ~ -1 + lat_zone * delta_t_z + lat_zone * wind_support_z + lat_zone * cross_wind_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

formula <- used ~ -1 + delta_t_z + wind_support_z + cross_wind_z + var_delta_t_z + var_ws_z +
  f(stratum, model = "iid", 
    hyper = list(theta = list(initial = log(1e-6),fixed = T))) + # stratum-specific  intercepts  are  implicitly estimated by modelling them as a random intercept with a fixed variance log(1e-6)(why is it a log?)
  f(species1, delta_t_z, model = "iid",  # what are values? i thought they correspond to the number of alternative steps, but I get an error when setting it to 1:49
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + #param corresponds to the precision priors assinged to the random slopes. see text
  f(species2, wind_support_z,  model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species3, cross_wind_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species4, var_delta_t_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) +
  f(species5, var_ws_z, model = "iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05))))

# Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 

m1 <- inla(formula, family ="Poisson",  #random effect for species
           control.fixed = list(
             mean = mean.beta,
             prec = list(default = prec.beta)),
             data = all_data)

#view fixed effects
m1$summary.fixed

#summary of posterior distrbution
m1$summary.hyperpar

#plot coefficients
Efxplot(m1) + theme_tufte() # theme_bw()

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
           data = all_data)

summary(mf)
Efxplot(mf) + theme_bw()



#from muff:
#' To fit the model with random slopes in INLA, we need to generate new (but identical) variables of individual ID (ID cannot be used multiple times in the model formula):
dat$ANIMAL_ID1 <- dat$ANIMAL_ID
dat$ANIMAL_ID2 <- dat$ANIMAL_ID
dat$ANIMAL_ID3 <- dat$ANIMAL_ID

#' Set the model formula as for the fixed-effects model, but now add three random slope terms, namely for river width and for the two levels of the categorical variable (STAU1 and REST1) which are not the reference categories. The priors for precision of the three random slopes are PC(3,0.05), while the intercept variance is again fixed:
formula.random <- Loc ~  -1 + STAU1 + REST1 + Sohlenbrei +  
  Breaks_Dis +
  f(str_ID,model="iid",hyper=list(theta = list(initial=log(1e-6),fixed=T))) +
  f(ANIMAL_ID1,Sohlenbrei,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID2,STAU1,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) + 
  f(ANIMAL_ID3,REST1,values=1:9,model="iid",
    hyper=list(theta=list(initial=log(1),fixed=F,prior="pc.prec",param=c(3,0.05)))) 

#' Fit the model
#' #' Set mean and precision for the priors of slope coefficients
mean.beta <- 0
prec.beta <- 1e-4 
r.inla.random <- inla(formula.random, family ="Poisson", data=dat, 
                      control.fixed = list(
                        mean = mean.beta,
                        prec = list(default = prec.beta)
                      )
)

#' The summary for the posterior distribution of the fixed effects:
r.inla.random$summary.fixed 

#' Since variances are parameterized and treated as precisions, the summary of the respective posterior distributions is given for the precisions:
r.inla.random$summary.hyperpar


#' Source R functions for calculating posterior means 
#' and medians of the precisions.
source("inla_emarginal.R")
source("inla_mmarginal.R")
inla_emarginal(r.inla.random)
inla_mmarginal(r.inla.random)



############################################
B1 <- glm(form, 
          family = binomial, data = all_data)

B2 <- glmer(used ~ lat_zone * var_delta_t_z + lat_zone * avg_delta_t_z + lat_zone * avg_cw_z + lat_zone * avg_ws_z + 
              (1 | unique_step_id) + (1 |species), 
            family = binomial, data = all_data) #singular


#adjust fixed structure
B2 <- glmer(used ~ lat_zone * avg_delta_t_z + lat_zone * avg_cw_z + lat_zone * avg_ws_z + 
              (1 | unique_step_id) + (1 |species), 
            family = binomial, data = all_data) #did not converge

#move lat from fixed to random
B2 <- glmer(used ~ var_delta_t_z + avg_delta_t_z + avg_cw_z + avg_ws_z + 
              (1 | unique_step_id) + (1 |species) + (1|lat_zone), 
            family = binomial, data = all_data) #



B3 <- glmer(used ~ var_delta_t + avg_delta_t + avg_cw + avg_ws + lat_zone + 
              (1 + var_delta_t + avg_delta_t + avg_cw + avg_ws| unique_step_id) + 
              (1 + var_delta_t + avg_delta_t + avg_cw + avg_ws|species), family = binomial, data = all_data)



#step 1: beyond optimal model with all variables, comparing different interaction structures


#all species together first! lol


var_model <- glmer(used ~ var_ws_z + var_cw_z + var_delta_t_z + var_ws_z:lat_zone + var_cw_z:lat_zone + var_delta_t_z:lat_zone +
                     (1 | unique_step_id) + (1|species), 
                   family = binomial, data = all_data )
summary(var_model)
r.squaredLR(var_model)


library(interactions)
interact_plot(var_model, pred = var_delta_t_z, modx = lat_zone, plot.points = TRUE)


avg_model <- glmer(used ~ avg_ws_z + var_cw_z + avg_delta_t_z + avg_ws_z:lat_zone + avg_cw_z:lat_zone + avg_delta_t_z:lat_zone +
                     (1 | unique_step_id) + (1|species), 
                   family = binomial, data = all_data )
summary(avg_model)
r.squaredLR(avg_model)

#model with only delta t
delta_t_model <- glmer(used ~ avg_delta_t_z + var_delta_t_z + (avg_delta_t_z|lat_zone) + (var_delta_t_z|lat_zone) +
                     (1 | unique_step_id) + (1|species), 
                   family = binomial, data = all_data )
summary(var_model)
r.squaredLR(var_model)

#gam... do I need to scale the variables?
all_data$fzone <- factor(all_data$lat_zone)
m1 <- gamm(used ~ s(scale(avg_u925),scale(avg_v925)) + avg_delta_t, 
           random = list(species = ~1), family = binomial, data = all_data )




lapply(data_z, function(x){

  
})


############## try geeglm (p. 319)
load("ssf_input_ann_z.RData") #all_data

library(geepack)

#reminder: avg and var cw were correlated.

geem1 <- geeglm(used ~ var_delta_t_z , 
                family = binomial(link = "cloglog"), data = all_data,
                id = stratum, corstr = "exchangeable") #crashes rstudio

############## try clogit I guess. lol

#make sure to have equal number of available steps for each stratum
strata_l <- all_data %>% 
  group_by(stratum) %>% 
  summarise(n= n()) %>% 
  dplyr::select(n) %>% 
  as.data.frame() #28 strata have lengths less than 50, 11 have lengths less than 49. keep the 49s...



form_all <- formula(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * avg_delta_t_z + 
                      lat_zone * avg_ws_z + lat_zone * avg_cw_z +
                      strata(stratum))

form1 <- formula(used ~ var_delta_t_z + var_ws_z + var_cw_z +  
             strata(stratum))
form2 <- formula(used ~ avg_delta_t_z + avg_ws_z + avg_cw_z +  
                   strata(stratum))
m1 <- clogit(form1, data = all_data)
m2 <- clogit(form2, data = all_data)

m_all <- clogit(form_all, data = all_data)

models <- lapply(split(all_data, all_data$species), function(x){
  m <- clogit(used ~ scale(avg_ws) + scale(avg_delta_t) + scale(delta_t) + scale(wind_support) + scale(cross_wind) + scale(avg_cw) +
                strata(stratum), data = x)
  summary(m)
})

for (i in c("O","OHB","PF")){
  x <- all_data[]
}

library(TwoStepCLogit)
m3 <- Ts.estim(used ~ lat_zone * var_delta_t_z + lat_zone * var_ws_z + lat_zone * avg_delta_t_z + 
                 lat_zone * avg_ws_z + lat_zone * avg_cw_z +  strata(stratum) + cluster(species),
           random = ~ var_delta_t_z + var_ws_z + avg_delta_t_z + avg_ws_z, 
           data = all_data, 
           D="UN(1)")

lapply(data_ls_z,function(x){
  
  #model with wind and average updraft conditions
  Ts.estim(y ~ cos(ta_) + log(sl_) + wind_support_z + cross_wind_z + dist_target_z + clim_cross_wind_z + 
             clim_wind_support_z +  strata(step_id_) + cluster(id),
           random = ~ wind_support_z + cross_wind_z + dist_target_z + clim_cross_wind_z + 
             clim_wind_support_z, 
           data = x, 
           D="UN(1)") #D is the structure of the output matrix
})



################################
#following Zuur et al protocol (p. 130).. 
# strata as a random effect. species as a randomm effect. latitude as an interaction term, because I actually care about it.

#step 1 : no random effects at first
form <- formula(used ~ lat_zone * delta_t_z + lat_zone * var_delta_t_z + lat_zone * wind_support_z + lat_zone * var_ws_z + 
                  lat_zone * cross_wind_z + lat_zone * var_cw_z )
M.1m <- lm(form, 
           data = all_data)

plot(M.1m, select = c(1))
E <- rstandard(M.1m)
boxplot(E ~ species, data = all_data, axes = F) #if all the data clouds around zero, no need to include species as a random effect. but that is not really the case here
abline(0,0); axis(2)
X11(); boxplot(E ~ unique_step_id, data = all_data, axes = F) 
abline(0,0); axis(2)

#step 2: make the same model with gls, to allow for comparisons with the mixed effect models
M.gls <- gls(form, data = all_data)

#step 3 and 4: choose a variance structure and build model (later try to add latitude as a varident structure. ch. 4 of book)
#try only random effect on intercept
M1.lme <- lme(form, random = ~ 1| species, method = "REML", data = all_data)

#step 5: compare models
anova(M.gls,M1.lme)
