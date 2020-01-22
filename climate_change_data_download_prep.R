#script for downloading and processing CMIP5 data from climate data store
# Jan 21. 2020. Elham Nourani. Radolfzell, Germany

library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

##### STEP1: connect to python #####

#import the python library cdsapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client() #start the connection... make sure cds key file is in Documents

##### STEP2: download wind data #####

#MPI_ESM_MR
for (period in c("200601-200912", "201001-201912", "204001-204912", "205001-205912", "209001-210012")) {
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/mpi_esm_mr/", paste(period, rcp,"wind.nc", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = c("u_component_of_wind", "v_component_of_wind"),
      model = "mpi_esm_mr",
      ensemble_member = "r1i1p1",
      period = period,
      format = "netcdf",
      dataset = "projections-cmip5-monthly-pressure-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-pressure-levels",
                    query,
                    target = target)
  }
}

#HADGEM2_CC
for (period in c("200512-203011", "203012-205511", "205512-208011", "208012-209912", "210001-210012")) {
  for (rcp in c("rcp_4_5","rcp_8_5")){
    
    target <- paste("/home/enourani/Documents/climate_data_store_downloads/hadgem2_cc/", paste(period, rcp,"wind.nc", sep = "_"),sep = "")
    
    query <- r_to_py (list(
      experiment = rcp,
      variable = c("u_component_of_wind", "v_component_of_wind"),
      model = "hadgem2_cc",
      ensemble_member = "r1i1p1",
      period = period,
      format = "netcdf",
      dataset = "projections-cmip5-monthly-pressure-levels"
    ))
    
    server$retrieve("projections-cmip5-monthly-pressure-levels",
                    query,
                    target = target)
  }
}
##### STEP3: download temperature data #####

#MPI_ESM_MR
for (rcp in c("rcp_4_5","rcp_8_5")){
  period <- "200601-210012"
  target <- paste("/home/enourani/Documents/climate_data_store_downloads/mpi_esm_mr/", paste(period, rcp,"temp.nc", sep = "_"),sep = "")
  
  query <- list(
    experiment = rcp,
    variable = c("2m_temperature", "sea_surface_temperature"),
    model = "mpi_esm_mr",
    ensemble_member = "r1i1p1",
    period = period,
    format = "netcdf",
    dataset = "projections-cmip5-monthly-single-levels"
  )
  
  server$retrieve("projections-cmip5-monthly-single-levels",
                  query,
                  target = target)
}

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

##### STEP4: process wind data #####


vname <- c("sst","t2m")
file_list <- list.files("/home/enourani/Documents/climate_data_store_downloads/mpi_esm_mr/",pattern = ".nc",full.names = TRUE)

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

