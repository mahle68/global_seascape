#download ERA-5 data for construction of the regional GAMs
#Elham Nourani. Feb 22. 2021. Radolfzell, DE


library(tidyverse)
library(reticulate)
library(ncdf4)
library(sf)
library(lubridate)
library(parallel)

setwd("/home/enourani/ownCloud/Work/Projects/delta_t/")

load("2021/extent_ls_regional_gam.RData") #the regional extents. extent_ls

##### STEP 1: connect to cdsapi server

#import the python library ecmwfapi
path <- "/home/enourani/.local/lib/python2.7/site-packages/"
cdsapi <- import_from_path("cdsapi", path = path)

server = cdsapi$Client()

##### STEP 2: request and download data #####

path <- "/home/enourani/Documents/ERA5_zones/"

yrs <- as.character(c(seq(1981,2020)))
mns <- c(str_pad(seq(1:9),2,"left","0"),"10","11","12")

for(yr in yrs){
for(mn in mns){ #for each month
  
  
  
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
                target = paste0(path,"sst_t2m_", yr, "_",  mn,".nc")) 
}}


# request <- r_to_py(list(
#   product_type = "reanalysis",
#   format = "netcdf",
#   variable = c("2m_temperature", "sea_surface_temperature"),
#   year = c("1980", "1981", "1982", "1983", "1984", "1985", "1986", "1987", "1988", "1989", "1990", "1991", "1992", "1993", "1994", "1995", "1996", "1997", "1998", "1999", "2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008", "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020"),
#   month = i,
#   day = c("01", "02", "03", "04", "05", "06", "07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", "31"),
#   time = c("00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", "10:00", "11:00", "12:00", "13:00", "14:00", "15:00", "16:00", "17:00", "18:00", "19:00", "20:00", "21:00", "22:00", "23:00"),
#   area = c(70, -180, -22, 180),
#   dataset = "reanalysis-era5-single-levels",
# ))






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
