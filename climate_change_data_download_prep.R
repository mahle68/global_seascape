#script for downloading and processing CMIP5 data from climate data store
# Jan 21. 2020. Elham Nourani. Radolfzell, Germany
#update Jan 23: only include HadGEM2-CC and ACCESS1.0. they both have decent resolution and the variables that I want. BUT, ACCESS1.0 needs regridding, as it has rotated pole grids.
#when one var, use nc format, when multiple, use zip.

library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)
library(raster)
library(PCICt) #to deal with 360 day calendar

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

##### STEP1: connect to python #####

#import the python library cdsapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client() #start the connection... make sure cds key file is in Documents

##### STEP2: download wind data #####

#HADGEM2_CC
for (period in c("200512-203011", "203012-205511", "205512-208011", "208012-209912", "210001-210012")) { #RCP 8.5 not available for 2100
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    period <- "210001-210012"
    rcp <- "rcp_8_5"
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/hadgem2_cc/", paste(period, rcp,"wind.zip", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = c("u_component_of_wind", "v_component_of_wind"),
      model = "hadgem2_cc",
      ensemble_member = "r1i1p1",
      period = period,
      format = "zip",
      dataset = "projections-cmip5-monthly-pressure-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-pressure-levels",
                    query,
                    target = target)
  }
}

#ACCESS1.0
for (period in c("200601-205512", "205601-210012")) {
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/access_1_0/", paste(period, rcp,"wind.zip", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = c("u_component_of_wind", "v_component_of_wind"),
      model = "access1_0",
      ensemble_member = "r1i1p1",
      period = period,
      format = "zip",
      dataset = "projections-cmip5-monthly-pressure-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-pressure-levels",
                    query,
                    target = target)
  }
}

##### STEP3: download temperature data #####

#HADGEM2_CC
for (period in c("200512-209911", "209912-210012")) {
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/hadgem2_cc/", paste(period, rcp,"ssf.nc", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = "sea_surface_temperature",
      model = "hadgem2_cc",
      ensemble_member = "r1i1p1",
      period = period,
      format = "netcdf",
      dataset = "projections-cmip5-monthly-single-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-single-levels",
                    query,
                    target = target)
  }
}

for (period in c("200512-203011", "203012-205511", "205512-208011", "208012-209912", "210001-210012")) {
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/hadgem2_cc/", paste(period, rcp,"t2m.nc", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = "2m_temperature",
      model = "hadgem2_cc",
      ensemble_member = "r1i1p1",
      period = period,
      format = "netcdf",
      dataset = "projections-cmip5-monthly-single-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-single-levels",
                    query,
                    target = target)
  }
}

#ACCESS1.0
period <- "200601-210012"
for (rcp in c("rcp_4_5","rcp_8_5")){
  
  target <- paste("/home/enourani/Documents/climate_data_store_downloads/access_1_0/", paste(period, rcp,"temp.zip", sep = "_"),sep = "")
  
  query <- r_to_py (list(
    experiment = rcp,
    variable = c("2m_temperature", "sea_surface_temperature"),
    model = "access1_0",
    ensemble_member = "r1i1p1",
    period = period,
    format = "zip",
    dataset = "projections-cmip5-monthly-single-levels"
  ))
  
  server$retrieve("projections-cmip5-monthly-single-levels",
                  query,
                  target = target)
}

##### STEP4: process temperature data #####

#HadGEM2-CC_ temperature
var <- "tas"
file_list <- list.files("/home/enourani/Documents/climate_data_store_downloads/hadgem2_cc",pattern = "t2m",full.names = TRUE)


lapply(file_list,function(x){

    nc <- nc_open(x)
  
  rcp <- ncatt_get(nc,0,"experiment_id")$value
  model <-ncatt_get(nc,0,"model_id")$value
    
  #extract lon and lat
  lat <- ncvar_get(nc,'lat')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'lon') #- 180 #longitude is from 0-360, subtract 180 to change the range to -180-180
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the "days since" into date 
  timestamp <- as.PCICt(t*(60*60*24), origin = "1859-12-01", cal = "360_day") %>%  #x has to be in seconds, my t is in days
    as.character() %>% 
    as_date()
  
  
  
  #put everything in a large df
  row_names <- expand.grid(lon,lat,timestamp)
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,var)), nrow = nlon * nlat * nt, ncol = 1) #array to vector to matrix
    #matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1)
  ))
  
  nc_close(nc)
  
  colnames(var_df) <- c("lon","lat","date",var)   #set column names
  
  #only keep data for latitude 0-60 and months 
  var_df_filter <- var_df %>%
    mutate(month = month(var_df$date),
           year = year(var_df$date)) %>% 
    filter(month %in% c(2:5,8:11) & between(lat,0,60))
  
  save(var_df_filter,file = paste0("/home/enourani/Documents/climate_data_store_downloads/all_processed/",
                                   paste(year(timestamp)[1],year(timestamp)[length(timestamp)], var, rcp, model, sep = "_"),".RData"))
  gc()
  gc()
})
#ACCESS1.0_ temperature at 2m (sst has roated pole grid system. use NCL to regrid)

var <- "tas"
file_list <- list.files("/home/enourani/Documents/climate_data_store_downloads/access_1_0/",pattern = var,full.names = TRUE)


lapply(file_list,function(x){
  
  rcp <- strsplit(x,"_")[[1]][9]
  model <- strsplit(x,"_")[[1]][8]
  
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'lat')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'lon') #- 180 #longitude is from 0-360, subtract 180 to change the range to -180-180
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the "days since" into date 
  timestamp <- as_date(t,origin = "0001-01-01") #monthly from 2006 to 2100
  
  #put everything in a large df
  row_names <- expand.grid(lon,lat,timestamp)
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,var)), nrow = nlon * nlat * nt, ncol = 1) #array to vector to matrix
    #matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1)
  ))
  
  colnames(var_df) <- c("lon","lat","date",var)   #set column names
  
  #only keep data for latitude 0-60 and months 
  var_df_filter <- var_df %>%
    mutate(month = month(var_df$date),
           year = year(var_df$date)) %>% 
    filter(month %in% c(2:5,8:11) & between(lat,0,60))
  
  save(var_df_filter,file = paste0("/home/enourani/Documents/climate_data_store_downloads/all_processed/",
                                   paste(year(timestamp)[1],year(timestamp)[length(timestamp)], var, rcp, model, sep = "_"),".RData"))
  gc()
  gc()
})
  



## visualize just to make sure
smpl <- var_df_filter[var_df_filter$year == 2015,]
coordinates(smpl) <- ~ lon + lat
gridded(smpl) <- TRUE
smpl_r <- raster(smpl,layer = 2)
maps::map("world")
plot(smpl_r, add= T)
