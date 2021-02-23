#download ERA-5 data for construction of the regional GAMs
#Elham Nourani. Feb 22. 2021. Radolfzell, DE
#update: first try was to dowanload all the data for a large extent. local system ran out of space. so, now, just download data for the regions
# of interest (in 2021_regional_gams.R) and only 6 levels for each day (4-hourly)


library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

load("R_files/2021/extent_ls_regional_gam.RData") #the regional extents. extent_ls

##### STEP 1: connect to cdsapi server

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

##### STEP 2: request and download data #####

path <- "/home/enourani/Documents/ERA5_zones/"

yrs <- as.character(c(seq(1981,2020)))
mns <- c(str_pad(seq(1:9),2,"left","0"),"10","11","12")


lapply(c(1:length(extent_ls)), function(zone){
  x <- extent_ls[[zone]]
for(yr in yrs){
#for(mn in mns){ #for each month
  
  query<- r_to_py(list(
    product_type = "reanalysis",
    area =  paste(x[4], x[1], x[2], x[3], sep = "/"),     #"60/-180/0/180", #N/W/S/E
    format = "netcdf",
    variable = c("2m_temperature", "sea_surface_temperature"),
    year = yr,
    month = c(str_pad(seq(1:9),2,"left","0"),"10","11","12"),
    day = str_pad(1:31,2,"left","0"),
    time = str_c(seq(0,23,4),"00",sep=":") %>% str_pad(5,"left","0"), #four-hourly data. 6 levels
    dataset = "reanalysis-era5-single-levels"
  ))
  
 server$retrieve("reanalysis-era5-single-levels",
                query,
                target = paste0(path,"sst_t2m_", names(extent_ls)[[zone]], "_", yr, ".nc")) 
}#}

})

