#download ERA-interim data for construction of the regional GAMs
#Elham Nourani. Feb 25. 2021. Radolfzell, DE
#update: first try was to dowanload Era-5 data, but it ate up all my memory, so going back to era-interim (see 2021_Era_5_download_prep.R)

library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)
#test...

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

load("R_files/2021/extent_ls_regional_gam.RData") #the regional extents. extent_ls (from 2021_regional_gams.R)

##### STEP 1: connect to cdsapi server

#import the python library ecmwfapi
path<-"/home/enourani/.local/lib/python2.7/site-packages/"

use_python(path_to_python)

#py_install("ecmwf-api-client") 

#import the python library ecmwfapi
ecmwf <- import_from_path("ecmwfapi", path = path) #replace path_to_library with the relevant path on your system

server <- ecmwf$ECMWFDataServer() #start the connection

##### STEP 2: request and download data #####

output_path <- "/home/enourani/Documents/ERA_interim_zones/"


lapply(c(2:length(extent_ls)), function(zone){
  x <- extent_ls[[zone]]
  
for (yr in as.character(c(1979:2018))) {
  
  dates <- paste(yr,"-01-01/to/",yr,"-12-31",sep = "") 
  target <- paste(yr, "_",  names(extent_ls)[[zone]],"_sst_t2m.nc",sep = "")
  
  yr_query <- r_to_py(list(
    area = paste(x[4], x[1], x[2], x[3], sep = "/"), #N/W/S/E
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
    target = paste0(output_path, target)
  ))
  server$retrieve(yr_query)
}
})

##### STEP 3: process the data #####
setwd("/home/enourani/Documents/ERA_interim_zones/")

vname <- c("sst","t2m")
file_list <- list.files(pattern = ".nc",full.names = TRUE)

#start the cluster
mycl <- makeCluster(detectCores() - 2)

clusterExport(mycl, list("vname","file_list")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  
  library(ncdf4)
  library(lubridate)
  library(tidyverse)
})


data_list <- parLapply(cl = mycl, file_list,function(x){
  
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
           hour = hour(date_time),
           year = year(date_time),
           region = str_sub(strsplit(x, "/")[[1]][2],start = 6,end = -12)) %>%
    data.frame()
  
  sample <- land_df %>% 
    sample_n(5000, replace = F)
  
  save(sample,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_interim_regional_sampled/",
                            str_sub(strsplit(x, "/")[[1]][2],end = -4),".RData"))
  gc()
  gc()
})

stopCluster(mycl)


##### STEP 4: put everything together #####

file_ls <- list.files("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_interim_regional_sampled/",".RData",full.names = TRUE)

data_df <- sapply(file_ls, function(x) mget(load(x)), simplify = TRUE) %>%
  reduce(rbind)

save(data_df, file = "/home/enourani/ownCloud/Work/Projects/delta_t/R_files/2021/ecmwf_regions_5ksample.RData")
