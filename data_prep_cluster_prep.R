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
#library(parallel)
#library(lutz)getwd()
#library(move)
#library(mapview)
#library(rWind)

wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/") #remove this before submitting to cluster

#load files. make sure they are all stored in the working directory
load("twz_sf.RData") #twz_sf
load("tmpz_sf.RData") #tmpz_sf
load("all_spp_unfiltered_updated_lc_0_removed_new_track_id.RData") #dataset
load("land_0_60.RData") #land_0_60
load("ocean_0_60.RData") #ocean

##### STEP 1: convert tracks to spatial lines and remove portions over land #####
ocean_sp <- as(ocean,"Spatial")
land_sp <- as(land_0_60,"Spatial")
coordinates(dataset) <- ~location.long + location.lat
proj4string(dataset) <- "wgs"

track_ls <- split(dataset,dataset$track)
track_ls <- track_ls[lapply(track_ls,nrow)>1] #remove tracks with one point


mycl <- makeCluster(9) #total number of tracks is 369, so 41 will be sent to each core

clusterExport(mycl, c("track_ls", "land_sp","ocean_sp","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  library(raster)
})

Lines_ls<-parLapply(mycl,track_ls,function(x){
  #find out if the track has any points over water
  over_sea <- intersect(x,ocean_sp) #track_ls needs to be spatial for this to work
  #if the track has any point over water, convert to spatial line and subset for sea
  if(nrow(over_sea) != 0){
  line<-SpatialLines(list(Lines(Line(matrix(x@coords,ncol=2)), ID=x$track[1])),proj4string = wgs)
  line_sea <- erase(line,land_sp)
  line_sea$track <- x$track[1]
  } else {
    line_sea <- NA
  }
  
  line_sea
})

stopCluster(mycl)

save(Lines_ls,file = "Lines_ls_no_land.RData") #each track has a 

##### STEP 5: break up tracks into sea-crossing segments and filter #####

#remove elements with 0 elements (tracks with no sea-crossing)
Lines_ls_no_na <- Lines_ls[lapply(Lines_ls,is.na) == FALSE] 

#only keep the track column
Lines_ls_no_na <- lapply(Lines_ls_no_na,"[",,"track")

#convert to one object
lines <- do.call(rbind,Lines_ls_no_na)

#filter segments
segs_filtered<- st_as_sf(lines) %>% #convert to sf object
  st_cast("LINESTRING") %>% #convert to linestring (separate the segments)
  mutate(length = as.numeric(st_length(.)),
         n = npts(.,by_feature = T)) %>% 
  filter(length >= 30000 & n <= 2) #remove sea-crossing shorter than 30 km and segment with less than 2 points 
#there are still problematic tracks in the above... :/

