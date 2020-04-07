#script to estimate the step selection function for water-crossing raptors.
#each segment is analyzed separately, so first I tried not burstifying the segments (1 hour continuous) because I'd lose a lot of points of already short segments
#but that made the distribution of turning angles and step lengths problematic. some step lenghts are too large. so, back to burstifying (Apr. 6)
#April 2. 2020. Radolfzell, Germany.
#Elham Nourani

library(dplyr)
library(purrr)
library(move)
library(sf)
library(circular)
library(CircStats)
library(fitdistrplus)
library(RNCEP)
library(lubridate)
library(mapview)

#meters_proj <- CRS("+proj=moll +ellps=WGS84")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

maps::map("world",fil=TRUE,col="ivory") #flyway

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/")


# STEP 1: prepare the input data#####
#open segments (not tracks, because tracks may be intersected by land)
load("segs_dt.RData") #segs_ann; prepared in track_based_prep_analyze_daily.R; filtered to be over 30 km and have more than 2 points

#remove spring and give different values to Osprey and Peregrine in each flyway
segs <- segs_ann %>% 
  dplyr::filter(season == "autumn") %>% 
  mutate(group = ifelse(species == "OHB", "OHB",
                        ifelse(species == "GFB", "GFB",
                               ifelse(species == "O" & st_coordinates(.)[,1] < -30, "O_A",
                                      ifelse(species == "O" & st_coordinates(.)[,1] > -30, "O_E",
                                             ifelse(species == "PF" & st_coordinates(.)[,1] < -30, "PF_A",
                                                    "PF_E")))))) %>% 
  dplyr::arrange(group,seg_id,date_time) %>% 
  as("Spatial") %>%
  as.data.frame() 

#create a move list
rows_to_delete <- unlist(sapply(getDuplicatedTimestamps(x = as.factor(segs$seg_id),timestamps = segs$date_time,sensorType = "gps"),"[",-1)) #get all but the first row of each set of duplicate rows
segs <- segs[-rows_to_delete,]

move_ls<-lapply(split(segs,segs$group),function(x){
  x<-as.data.frame(x)
  mv<-move(x=x$x,y=x$y,time=x$date_time,data=x,animal=x$seg_id,proj=wgs)
  mv
})

#for each species/flyway, thin the data, burstify, and produce alternative steps
lapply(move_ls,function(group){ #each group is a species/flyway combo
  
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
  X11();par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (m)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  
  #--STEP 5: produce alternative steps
  used_av_seg <- lapply(sp_obj_ls, function(seg){ #for each segment
    
    used_av_burst <- lapply(split(seg,seg$burst_id),function(burst){ #for each burst,
      
      used_av_step <- lapply(c(2:(length(burst)-1)), function(this_point){ #first point has no bearing to calc turning angle, last point has no used endpoint.
        
        current_point<- burst[this_point,]
        previous_point<-burst[this_point-1,]
        used_point <- burst[this_point+1,] #this is the next point. the observed end-point of the step starting from the current_point
        
        #randomly generate 49 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 49, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl<-rgamma(n= 49, shape=fit.gamma1$estimate[[1]], rate= fit.gamma1$estimate[[2]])*1000  #generate random step lengths from the gamma distribution. make sure unit is meters
        
        #calculate bearing of previous point
        #prev_bearing<-bearing(previous_point,current_point) #am I allowing negatives?... no, right? then use NCEP.loxodrome
        prev_bearing<-NCEP.loxodrome(previous_point@coords[,2], current_point@coords[,2],
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
                 used = c(1,rep(0,49))) #one hour after the start point of the step
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

#get rid of alt. points that fall over land (do this later for all alternative poinst together :p)



maps::map("world",xlim = c(-75,-70), ylim = c(15,25),fil = TRUE,col = "ivory") #flyway
points(burst,col = "grey", pch = 16, cex = 0.5)
points(previous_point,col = "green", pch = 16, cex = 1)
points(current_point,col = "red", pch = 16, cex = 1)
points(rnd_sp, col = "orange", pch = 16, cex = 0.5)
points(used_point, col = "purple", pch = 16, cex = 1)

#plotting
#r <- mapview(burst)
#r + mapview(current_point,color = "red") + mapview(previous_point, color = "green")



#check for duplicated time-stamps
rows_to_delete <- sapply(getDuplicatedTimestamps(x = as.factor(segs$seg_id),timestamps = segs$date_time,sensorType = "gps"),"[",2) #get the second row of each pair of duplicate rows
segs <- segs[-rows_to_delete,]

mv <- move(x = segs$x,y = segs$y,time = segs$date_time,
           data = segs,animal = segs$seg_id,proj = wgs)

segs_th_an <- lapply(split(mv),function(one_seg){ #for each segment within the group
  seg_th <- one_seg %>%
    thinTrackTime(interval = as.difftime(1, units = 'hours'),
                  tolerance = as.difftime(15, units = 'mins'))
  
  #convert back to a move object (from move burst)
  seg_th <- as(seg_th,"Move")
  
  #calculate step lengths and turning angles
  if(nrow(seg_th) == 1){
    seg_th$step_length <- NA
    seg_th$turning_angle <- NA
  } else {
    seg_th$step_length <- c(distance(seg_th),NA)
    seg_th$turning_angle <- c(NA,turnAngleGc(seg_th),NA) #if the segment has less than three points, all will be NA
  }
  
  seg_th %>% 
    as.data.frame() %>% 
    dplyr::select(c(x, y, track, seg_id, date_time, length, step_length, turning_angle, species))
  
}) %>%
  reduce(rbind)


# STEP 2: summarise distribution of step length and turning angle #####
   #per species
distr <- lapply(split(segs_th_an, segs_th_an$species), function(species){
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu<-mean.circular(rad(species$turning_angle[complete.cases(species$turning_angle)])) #unit is considered in radians
  kappa <- est.kappa(rad(species$turning_angle[complete.cases(species$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl<-species$step_length[complete.cases(species$step_length) & species$step_length > 0]/1000 #remove 0s and NAs 
  
  #some step lengths are over 1200 km!!?? what?
  
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  
  #plot
  X11();par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 1200, col = "blue")
  
  hist(rad(species$turning_angle[complete.cases(species$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
})

  #----------------STEP 4: estimate step length and turning angle distributions
  #put everything in one df
  sp_obj_ls<-Filter(function(x) length(x)>1, sp_obj_ls) #remove tracks with no observation (these have only one obs due to the assignment of track id)
  bursted_df<-do.call(rbind,lapply(sp_obj_ls,as.data.frame))
  bursted_df<-bursted_df[,-which(names(bursted_df) %in% c("coords.x1","coords.x2","coords.x1.1","coords.x2.1"))]#remove repeated variables for lat and long
  
  #estimate von Mises parameters for turning angles
  #calculate the averages (mu).steps: 1)convert to radians. step 2) calc mean of the cosines and sines. step 3) take the arctan.OR use circular::mean.circular
  mu<-mean.circular(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  kappa <- est.kappa(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]))
  
  #estimate gamma distribution for step lengths and CONVERT TO KM!!! :p
  sl<-bursted_df$step_length[complete.cases(bursted_df$step_length) & bursted_df$step_length > 0]/1000 #remove 0s and NAs
  fit.gamma1 <- fitdist(sl, distr = "gamma", method = "mle")
  
  #plot
  X11();par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(species$turning_angle[complete.cases(species$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
  plot(function(x) dvonmises(x, mu = mu, kappa = kappa), add = TRUE, from = -3.5, to = 3.5, col = "red")
  
  #------------------------------------------------------------------------------
  #----------------STEP 5: produce alternative steps
  alt_steps_ls<-lapply(sp_obj_ls, function(track){ #within each track
    alt_steps_burst<-lapply(split(track,track$burst_id),function(burst){ #within each burst,
      
      for(this_point in 2:length(burst)){ #for each point. start from the second point. the first point has no bearing
        current_point<- burst[this_point,]
        previous_point<-burst[this_point-1,]
        
        #randomly generate 49 step lengths and turning angles
        rta <- as.vector(rvonmises(n = 49, mu = mu, kappa = kappa)) #generate random turning angles with von mises distribution (in radians)
        rsl<-rgamma(n= 49, shape=fit.gamma1$estimate[[1]], rate= fit.gamma1$estimate[[2]])  #generate random step lengths from the gamma distribution
        
        #calculate bearing of previous point
        prev_bearing<-bearing(previous_point,current_point)
        #find the gepgraphic location of each alternative point
        #calculate bearing for the current point. this is the bearing of previous step plus turning angle. then, subtract 360 if it is larger than 360.
        
        #BE CAREFUL! the angle here should be the bearing, not ta_
        #ALSO, find a good way to convert lat and lon to meters in a consistent way.
        
        #e turning angle. so, calculate bearing:add ta to the bearing of the previous point?
        rnd<- data.frame(long=start_point$long+slr*cos(tar),lat=start_point$lat+slr*sin(tar)) #for this to work, lat and lon should be in meters as well. boo. coordinates in meters?
        
        #assign date_time....two hours after the start point of the step
        rnd$date_time<- start_point$date_time+ hours(2)
        
        #check the random points. covnert back to latlon for plotting
        rnd_sp<-rnd
        coordinates(rnd_sp)<-~long+lat
        proj4string(rnd_sp)<-meters_proj
        rnd_sp<-spTransform(rnd_sp,wgs)
        
        #also convert  starting point back to lat-long
        coordinates(start_point)<-~long+lat
        proj4string(start_point)<-meters_proj
        start_point_wgs<-spTransform(start_point,wgs)
      }
    })
  })

  
# STEP 3: generate random steps. remove random steps over land #####

# STEP 4: annotate ####

# STEP 5: glmm! #####