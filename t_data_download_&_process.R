#script for downloading and processing Era-inerim temperature data for constructing global energy seascapes
#Nov. 2019. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)


##### STEP1: connect to python #####

path_to_python <- "C:/Users/mahle/AppData/Local/Programs/Python/Python37"
use_python(path_to_python)

#import python CDS-API
py_install("ecmwf-api-client") 

#import the python library ecmwfapi
ecmwf <- import_from_path("ecmwfapi",path = "C:/Users/mahle/Anaconda3/envs/r-reticulate/Lib/site-packages")

server = ecmwf$ECMWFDataServer() #start the connection... make sure cds key file is in Documents

##### STEP2: download temp data #####

#non_parallel
for (yr in as.character(c(1979:2018))) {
  
  dates <- paste(yr,"-01-01/to/",yr,"-12-31",sep = "") 
  target <- paste(yr,"_sst_t2m.nc",sep = "")
  
  yr_query <- r_to_py(list(
    area = "60/-180/-60/180", #N/W/S/E
    class = 'ei',
    dataset = "interim",
    date = dates,
    expver = "1",
    grid = "0.75/0.75",
    levtype = "sfc",
    param = "34.128/167.128",
    step = "0",
    stream = "oper",
    time = "00:00:00/06:00:00/12:00:00/18:00:00", 
    type = "an",
    format = "netcdf",
    target = paste("D:/ERA_Interim_temp_all/",target,sep = "")
  ))
  server$retrieve(yr_query)
}


##### STEP3: process  temp data #####

setwd("D:/ERA_Interim_temp_all/")

vname <- c("sst","t2m")
file_list <- list.files(pattern = ".nc",full.names = TRUE)

#start the cluster
mycl <- makeCluster(detectCores() - 1)

clusterExport(mycl, list("vname","file_list")) #define the variable that will be used within the function

clusterEvalQ(mycl, {

  library(ncdf4)
  library(lubridate)
  library(dplyr)
})


data_list <- parLapply(cl = mycl,file_list,function(x){
  
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'latitude')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'longitude')
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the hours into date + hour
  timestamp <- as_datetime(c(t*60*60),origin = "1900-01-01")
  
  #put everything in a large df
  row_names <- expand.grid(lon,lat,timestamp)
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,vname[1])), nrow = nlon * nlat * nt, ncol = 1), #array to vector to matrix
    matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1)))
  
  colnames(var_df) <- c("lon","lat","date_time",vname)   #set column names
  
  #remove points over land (NAs)
  
  land_df <- var_df %>%
    na.omit() %>%
    mutate(delta_t = sst - t2m,
           yday = yday(date_time),
           hour = hour(date_time)) %>%
    data.frame()
  
  sample <- land_df %>% 
    sample_n(50000)
  
  save(sample,file = paste0("C:/Users/mahle/ownCloud/Work/Projects/delta_t/ERA_INTERIM_data_global_sampled/",
                          paste(year(timestamp)[1],min(sample$yday),max(sample$yday),sep = "_"),".RData"))
  gc()
  gc()
})

stopCluster(mycl)


##### Quick plotting #####
windows()
map("world",fill = TRUE,col = "beige")
points(sample$lon,sample$lat,cex = 0.3,pch = 16,col = "blue")
#####


