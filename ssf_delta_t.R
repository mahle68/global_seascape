#script to estimate the step selection function for water-crossing raptors.
#each segment is analyzed separately, so I'm not burstifying the segments (1 hour continuous) because I will lose a lot of points of already short segments
#April 2. 2020. Radolfzell, Germany.
#Elham Nourani

library(dplyr)
library(purrr)
library(move)
library(sf)

#meters_proj <- CRS("+proj=moll +ellps=WGS84")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")

# STEP 1: prepare the input data#####
#open segments (not tracks, because tracks may be intersected by land)
load("segs_dt.RData") #segs_ann; prepared in track_based_prep_analyze_daily.R; filtered to be over 30 km and have more than 2 points

#remove spring
segs <- segs_ann %>% 
  dplyr::filter(season == "autumn") %>% 
  dplyr::arrange(species,seg_id,date_time) %>% 
  as("Spatial") %>%
  as.data.frame() 

#thin hourly
#create a move object

#check for duplicated time-stamps
rows_to_delete <- sapply(getDuplicatedTimestamps(x = as.factor(segs$seg_id),timestamps = segs$date_time,sensorType = "gps"),"[",2) #get the second row of each pair of duplicate rows
segs <- segs[-rows_to_delete,]

#create a move list
# move_ls<-lapply(split(segs,segs$species),function(x){
#   x<-as.data.frame(x)
#   mv<-move(x=x$x,y=x$y,time=x$date_time,data=x,animal=x$seg_id,proj=wgs)
#   mv
# })

mv <- move(x = segs$x,y = segs$y,time = segs$date_time,
           data = segs,animal = segs$seg_id,proj = wgs)

segs_th <- lapply(split(mv),function(one_seg){ #for each track within the group

  seg_th <- one_seg %>%
    thinTrackTime(interval = as.difftime(1, units = 'hours'),
                  tolerance = as.difftime(15, units = 'mins'))

  #convert back to a move object (from move burst)
  seg_th <- as(seg_th,"Move")
  seg_th$track <- one_seg@idData$track #reassign the track
  seg_th$seg_id <- one_seg@idData$seg_id
  seg_th$length <- one_seg@idData$length
  seg_th$species <- one_seg@idData$species

  seg_th
}) %>%
  reduce(rbind)


# STEP 2: summarise distribution of step length and turning angle #####

#calculate step lengths and turning angles

segs_th_an <- lapply(split(segs_th, segs_th$species),)
burst_ls<-split(thinned_track,thinned_track$burst_id)
burst_ls<-Filter(function(x) length(x)>=3, burst_ls) #remove bursts with less than 3 observations


burst_ls<-lapply(burst_ls,function(burst){
  burst$step_length<-c(distance(burst),NA) #
  burst$turning_angle<-c(NA,turnAngleGc(burst),NA)
  burst
})

#put burst_ls into one dataframe
bursted_sp<-do.call(rbind,burst_ls)
bursted_sp$track<-track@idData$track #reassign the track id
bursted_sp











lapply(move_ls,function(species){ #for each of the groups (country-continent specific)
  
  sp_obj_ls<-lapply(split(species),function(seg){ 
  
  # thinned_track<-seg %>%
  #   thinTrackTime(interval = as.difftime(2, units='hours'),
  #                 tolerance = as.difftime(15, units='mins'))#the unselected bursts are the large gaps between the selected ones
  # 
  #   #---STEP 2: assign burst IDs (each chunk of track with 2 hour intervals is one burst... longer gaps will divide the brusts) 
  #   thinned_track$selected <- c(as.character(thinned_track@burstId),NA) #assign selected as a variable
  #   
  #   thinned_track$burst_id<-c(1,rep(NA,nrow(thinned_track)-1)) #define value for first row
  #   
  #   for(i in 2:nrow(thinned_track@data)){
  #     
  #     if(i== nrow(thinned_track@data)){
  #       thinned_track@data$burst_id[i]<-NA
  #     } else
  #       if(thinned_track@data[i-1,"selected"] == "selected"){
  #         thinned_track@data$burst_id[i]<-thinned_track@data[i-1,"burst_id"]
  #       } else {
  #         thinned_track@data$burst_id[i]<-thinned_track@data[i-1,"burst_id"]+1
  #       }
  #   }
  #   
  #   #convert back to a move object (from move burst)
  #   thinned_track <- as(thinned_track,"Move")
    
    #----------------STEP 3: calculate step lengths and turning angles 
    #sl_ and ta_ calculations should be done for each burst. converting to a move burst doesnt make this automatic. so just split manually
    burst_ls<-split(thinned_track,thinned_track$burst_id)
    burst_ls<-Filter(function(x) length(x)>=3, burst_ls) #remove bursts with less than 3 observations
    
    
    burst_ls<-lapply(burst_ls,function(burst){
      burst$step_length<-c(distance(burst),NA) #
      burst$turning_angle<-c(NA,turnAngleGc(burst),NA)
      burst
    })
    
    #put burst_ls into one dataframe
    bursted_sp<-do.call(rbind,burst_ls)
    bursted_sp$track<-track@idData$track #reassign the track id
    bursted_sp
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
  windows();par(mfrow=c(1,2))
  hist(sl,freq=F,main="",xlab = "Step length (km)")
  plot(function(x) dgamma(x, shape = fit.gamma1$estimate[[1]],
                          rate = fit.gamma1$estimate[[2]]), add = TRUE, from = 0.1, to = 150, col = "blue")
  
  hist(rad(bursted_df$turning_angle[complete.cases(bursted_df$turning_angle)]),freq=F,main="",xlab="Turning angles (radians)")
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