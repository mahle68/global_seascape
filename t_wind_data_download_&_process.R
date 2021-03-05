#script for downloading and processing Era-inerim temperature data for constructing global energy seascapes
#Nov. 2019. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

##### STEP1: connect to python #####

#emcwf <- import("ecmwfapi")
#path_to_python <- "/usr/bin/python3.7"
#use_python(path_to_python)


#import python CDS-API
#py_install("ecmwf-api-client") 
#py_install


#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/ecmwfapi"
ecmwf <- import_from_path("ecmwfapi", path = path)
#ecmwf <- import_from_path("ecmwfapi",path = "C:/Users/mahle/Anaconda3/envs/r-reticulate/Lib/site-packages")

server = ecmwf$ECMWFDataServer() #start the connection... make sure cds key file is in Documents

##### STEP2: download temp data #####

#non_parallel
for (yr in as.character(c(1979:2018))) {
  
  dates <- paste(yr,"-02-01/to/",yr,"-10-31",sep = "") 
  target <- paste(yr,"_sst_t2m.nc",sep = "")
  
  yr_query <- r_to_py(list(
    area = "60/-180/0/180", #N/W/S/E
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
    target = paste("ERA_INTERIM_data_0_60/",target,sep = "")
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


##### STEP4: download wind data #####

#non_parallel
for (yr in as.character(c(1979:2018))) {
  
  dates <- paste(yr,"-02-01/to/",yr,"-10-31",sep = "") 
  target <- paste(yr,"_wind.nc",sep = "")
  
  yr_query <- r_to_py(list(
    area = "60/-180/0/180", #N/W/S/E
    class = 'ei',
    dataset = "interim",
    date = dates,
    expver = "1",
    grid = "0.75/0.75",
    levelist = "925",
    levtype = "pl",
    param = "131.128/132.128/135.128",
    step = "0",
    stream = "oper",
    time = "00:00:00/06:00:00/12:00:00/18:00:00", 
    type = "an",
    resol = "av",
    format = "netcdf",
    target = paste("ERA_INTERIM_data_0_60/",target,sep = "")
  ))
  server$retrieve(yr_query)
}

#######################################
#download ecmwf data using ecmwfr
library(ecmwfr)
#https://cran.r-project.org/web/packages/ecmwfr/vignettes/webapi_vignette.html

wf_set_key(user = "mahle68@gmail.com",
           key = "ff8a52b16561c01bf792ca70e4a806b0",
           service = "webapi")


wf_get_key(user = "mahle68@gmail.com",
           service= "webapi")

wf_set_key(service = "webapi")

years<-c(1979:2019) 
for (yr in 1:length(years)){
  this_year<-years[yr]
  dates<-paste(this_year,"-02-1/to/",this_year,"-10-31",sep = "") 
  target<-paste(this_year,"_wind.nc",sep="")
  wf_request(user="mahle68@gmail.com",
             request = list(
               area="60/-180/0/180", #N/W/S/E
               class='ei',
               dataset= "interim",
               date= dates,
               expver= "1",
               grid= "0.75/0.75",
               levelist= "925",
               levtype="pl", #pressure level
               param= "131.128/132.128/135.128",
               resol="av",
               step= "0",
               stream="oper",
               time="06:00:00/12:00:00/18:00:00", 
               type="an",
               format= "netcdf",
               target=target),
             transfer = T,
             path="ERA_INTERIM_data_0_60/")
}




##### Quick plotting #####
windows()
map("world",fill = TRUE,col = "beige")
points(sample$lon,sample$lat,cex = 0.3,pch = 16,col = "blue", )
points(ann_z[ann_z$month == 2,c("location.long","location.lat")], pch = 16, cex = 0.4)
plot(st_geometry(land_eur_15km))
#####


setwd("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_data_0_60")
