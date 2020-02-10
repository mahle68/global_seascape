#script for converting tracking points to spatial lines and removing parts of the tracks that fall over land
#follows up on data_prep_track_based.R and data_prep_track_based_no_interp_preliminary.R
#Elham Nourani. Feb. 6. 2020. Radolfzell, Germany

#not an array job

#open libraries
library(tidyverse)
library(lubridate)
library(sf)
library(raster)
library(parallel)


#define variables
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")

#load files. make sure they are all stored in the working directory
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset
load("land_0_60.RData") #land_0_60

##### STEP 1: convert land polygon from sf to sp #####
land_sp <- as(land_0_60,"Spatial")

##### STEP 2:split the data into a list of tracks and remove tracks with no points #####
track_ls<-split(tr,tr$track)
track_ls<-track_ls[lapply(track_ls,nrow)>0]

##### STEP 3: for each track, create a spatial line and remove portions over land #####

mycl <- makeCluster(9) #total number of tracks is 369, so 41 will be sent to each core

clusterExport(mycl, c("track_ls", "land_sp", "wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(raster)
})

Lines_ls<-lapply(track_ls,function(x){
  line<-SpatialLines(list(Lines(Line(matrix(x@coords,ncol=2)), ID=x$track[1])),proj4string = wgs)
  line_sea <- erase(line,land_sp)
  line_sea$track <- x$track[1]
  line_sea
})

stopCluster(mycl)

##### STEP 4: save the results #####
save(Lines_ls,file = "Lines_ls_no_land.RData")
