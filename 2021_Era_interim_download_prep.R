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
