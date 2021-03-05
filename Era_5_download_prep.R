#download ERA-5 data for manual interpolation of global sea-crossing data
#Elham Nourani. Feb 18. 2020. Radolfzell, DE


library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

##### STEP 1: connect to cdsapi server

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

##### STEP 2: extract years and months for which I need data #####

#create alternative hourly segments (in track_based_prep_analyze_hourly.R, i did 6-hourly)

load("segs_dt.RData") #segs_ann .....do this on the data that has the alternatives. the number of months could increase!!!

segs_df <- segs_ann %>% 
  as("Spatial") %>% 
  as.data.frame()

hours_to_add <- c(0,cumsum(rep(-1,72))) #six hourly data addition for three days prior to the observed segment

pts_alt <- segs_df %>% 
  mutate(obs_id = row_number()) %>% 
  slice(rep(row_number(),73)) %>%  #paste each row 18 time for 18 additional hours
  arrange(obs_id) %>% 
  group_by(obs_id) %>% 
  mutate(hours_to_add = hours_to_add) %>% 
  mutate(alt_date_time = date_time + hours(hours_to_add)) %>%  #use hours to add as an id for alternative segments
  ungroup() %>% 
  mutate(used = ifelse(hours_to_add == 0,1,0)) %>% 
  as.data.frame()


ym <- pts_alt %>% 
  group_by(year,month) %>% 
  summarise(freq = n())  %>% 
  as.data.frame()
  
##### STEP 3: request and download data #####

yms <- lapply(split(ym,ym$year),function(x){ #create a vector of months for each year
  mths <- str_pad(x$month,2,"left","0")
  list(mths = mths, yr = x$year[1])
})


lapply(yms[-c(1:9)],function(x){
  yr<- x$yr
  mnths <- x$mths
  
  lapply(mnths,function(mn){
    
    query<- r_to_py(list(
      product_type = "reanalysis",
      area="60/-180/0/180", #N/W/S/E #limit to the tropical and temperate zones (exclude the arctic and antarctic circles)
      format = "netcdf",
      variable = c("2m_temperature", "sea_surface_temperature"),
      year = yr,
      month = mn,
      day = str_pad(1:31,2,"left","0"),
      time = str_c(0:23,"00",sep=":") %>% str_pad(5,"left","0"),
      dataset = "reanalysis-era5-single-levels"
    ))
    
    server$retrieve("reanalysis-era5-single-levels",
                    query,
                    target=paste0("/home/enourani/Documents/ERA_5_for_track_analysis/",paste(yr,mn,sep="_"),"_t2m_sst.nc"))

    
  })
  
})
