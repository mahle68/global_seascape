#Global energy seascape regional GAMs

#This is my attempt to put all the steps in one clean script. The script is not fully reproducable as you need to download the temperature data and shapefiles
#separately. It is mainly aimed at clarifying the methodology further and perhaps providing some inspiration for your own work. 
#Jul 2. 2020. Radolfzell, Germany.
#by Elham Nourani, PhD. enourani@ab.mpg.de

#load libraries
library(tidyverse) #for data wrangling- mostly dplyr
library(reticulate) #to enable using python in R (for downloading data from ECMWF)
library(ncdf4) #to manipulate netcdf data from ECMWF
library(sf) #spatial data wrangling (vectors)
library(raster) #spatial data wrangling (rasters)
library(lubridate) #to deal with time
library(parallel) #for parallel computing
library(mgcv) #for Generalized Additive Modeling
library(itsadug) #for effect plots
library(berryFunctions) #to add small plots to the main map or add.image from fields package

#define variables
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")

#define working directory
setwd("mywd")

## STEP 1: data donwload ####

#1.1-------------------------------
# connect to the ecmwf server

#to donwload data from ECMWF directly from R, follow this tutorial to make sure you have python and the ECMWF API: 
#https://dominicroye.github.io/en/2018/access-to-climate-reanalysis-data-from-r/

path_to_python <- "path_to_python/Python37" #replace path-to-python with path to python on your system
use_python(path_to_python)

#import python CDS-API
py_install("ecmwf-api-client") 

#import the python library ecmwfapi
ecmwf <- import_from_path("ecmwfapi",path = "path_to_library") #replace path_to_library with the relevant path on your system

server <- ecmwf$ECMWFDataServer() #start the connection

#1.2----------------------
# download the data

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
    target = paste("directory_of_choice/",target,sep = "")
  ))
  server$retrieve(yr_query)
}

#1.3---------------------
# process the data

setwd("directory_of_choice/")

vname <- c("sst","t2m")
file_list <- list.files(pattern = ".nc",full.names = TRUE)

#start the cluster (if you are doing this on multiple cores)
mycl <- makeCluster(detectCores() - 1)

clusterExport(mycl, list("vname","file_list")) #import variables that will be used by all cores

clusterEvalQ(mycl, { #import packages that will be used by all cores
  
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
  
  #take a sample of the data to make further analysis computationally less expensive
  
  sample <- land_df %>% 
    sample_n(50000)
  
  save(sample,file = paste0("new_directory_of_choice/",
                            paste(year(timestamp)[1],min(sample$yday),max(sample$yday),sep = "_"),".RData"))
})

stopCluster(mycl) #stop the cluster






## STEP 2: data preparation ####

#2.1------------
# merge the data
file_ls <- list.files("new_directory_of_choice",".RData",full.names = TRUE)

data_df <- file_ls %>%
  map(read_csv) %>% 
  reduce(rbind) %>%
  mutate(year = year(date_time))

#2.2--------------
# spatial filter and add solar elevation angle 

#download waterbodies layer from: https://www.naturalearthdata.com/downloads/10m-physical-vectors/10m-ocean/ and unzip
ocean <- st_read("ne_110m_ocean.shp") %>% 
  slice(2) %>% #remove the caspian sea
  st_crop(xmin = -180, ymin = 0, xmax = 180, ymax = 60) #crop to geographic extents of the study

data_sf <- data_df %>% 
  filter(between(lat,0,60)) %>% #filter for lat zone 0-60
  as.data.frame() %>% 
  st_as_sf(coords = c("lon","lat"), crs = wgs) %>% 
  st_intersection(ocean) %>% #filter out lakes
  mutate(s_elev_angle = solarpos(st_coordinates(.), date_time, proj4string=CRS("+proj=longlat +datum=WGS84"))[,2]) %>% #calculate solar elevation angle
  mutate(sun_elev = ifelse(s_elev_angle < -6, "night", #create a categorical variable for teh position of the sun
                           ifelse(s_elev_angle > 40, "high", "low")),
         month = month(date_time))



#2.3----------------------
# create regional datasets

#dowload the shapefile for world waterbodies (IHO Sea Areas): https://www.marineregions.org/downloads.php
#This will come in handy later to remove unwanted seas from our regional maps

world_seas <- st_read("World_Seas_IHO_v3.shp")

East_asia <- data %>% 
  st_crop(xmin = 99, xmax = 130, ymin = 2, ymax = 35) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#for indian ocean region, we need to remove the Persian gulf
pg <- world_seas %>%
  filter(NAME == "Persian Gulf") #prepare the persian gulf layer to use for masking

Indian_ocean <- data %>% 
  st_crop(xmin = 43, xmax = 78, ymin = 9, ymax = 26) %>% 
  st_difference(pg) %>%  #remove points over the persian gulf
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#for the Americas, we need to remove the pacific ocean
np_ocean <- world_seas %>% 
  filter(NAME == "North Pacific Ocean") %>% 
  st_crop(xmin = -118, xmax = -60, ymin = 0, ymax = 29) 

Americas <- data %>% 
  st_crop(xmin = -98, xmax = -54, ymin = 5.7, ymax = 48) %>% 
  st_difference(np_ocean) %>% 
  st_difference(outlier) %>%  #some points remain on the Pacific ocean. If this also happens to you, identify these outliers and remove them.
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#for Europe, include the Mediterranean and the Baltic seas ###############FIX THE FOLLOWING USING WOLRD_SEAS

eur_sea <- world_seas %>% 
  filter(NAME %in% c("Balearic (Iberian Sea)", "Alboran Sea", "Aegean Sea", "Adriatic Sea","Mediterranean Sea - Eastern Basin", 
                     "Mediterranean Sea - Western Basin","Tyrrhenian Sea", "Ionian Sea", "Ligurian Sea", "Baltic Sea",
                     "Gulf of Riga", "Gulf of Bothnia", "Gulf of Finland")) %>% 
  st_union() 
         
Europe <- data %>% 
  st_crop(xmin = -11, xmax = 37, ymin = 30, ymax = 62) %>% 
  st_intersection(eur_sea) %>% 
  as("Spatial") %>% 
  as.data.frame() %>% 
  rename(lon = coords.x1,
         lat = coords.x2) %>% 
  mutate(sun_elev_f = factor(sun_elev)) %>% 
  as.data.frame()

#put all regional data in a list
data_ls <- list(East_asia = East_asia, Americas = Americas, Indian_ocean = Indian_ocean, Europe = Europe)






## STEP 3: modeling ####

#3.1-------------------------
# modeling and model checking

#start the cluster (if parallel computing)
mycl <- makeCluster(detectCores() - 1) 

clusterExport(mycl, "data_ls") 

clusterEvalQ(mycl, {
  library(mgcv)
})

models_ls <- lapply(data_ls, function(x){
  #bam is just gam, but for large datasets
  bam(delta_t ~ s(lat,lon, by = sun_elev_f, k = 100) +
        s(yday, by = sun_elev_f, bs = "cc") +
        s(year, bs = "re") + #add year as a random effect
        sun_elev_f , method = "REML", data = x, cluster = mycl)
  
})

stopCluster(mycl)


#model checking
lapply(models_ls, function(x){
  X11()
  par(mfrow= c(2,2), oma = c(0,0,3,0))
  gam.check(x)
})

#3.2--------------------------------
#predict with the model (for the migratory season)

#timing of migration points over the open sea (from tracking data; for Amur falcon, extracted from available sources)
timing <- list(OHB = c(260:294), #Oriental honey buzzard
               GFB = c(277:299), #Grey-faced buzzard
               AF = c(319:323),  #Amur falcon
               EF = c(288:301), #Eleonora's falcon
               PF_EU = c(261:289), #Peregrine falcon-Europe
               O_EU = c(222:277), #Osprey-Europe
               PF_A = c(279:305), #Peregrine falcon-Americas
               O_A = c(244:299)) #Osprey-Americas

#create a list of timings for each region
timing_areas <- list(East_asia = c(min(c(timing$GFB,timing$OHB)):max(c(timing$GFB,timing$OHB))),
                     Americas = c(min(c(timing$O_A,timing$PF_A)):max(c(timing$O_A,timing$PF_A))),
                     Indian_ocean = timing$AF,
                     Europe = c(min(c(timing$O_EU,timing$PF_EU,timing$EF)):max(c(timing$O_EU,timing$PF_EU,timing$EF)))) #consider separating this into regions


#make predictions
mycl <- makeCluster(detectCores() - 1) 

clusterExport(mycl, list("timing_areas","models_ls", "data_ls", "wgs"))

clusterEvalQ(mycl, {
  library(mgcv)
  library(raster)
  library(dplyr)
  library(sf)
  library(fields) #for Thin-plate spline interpolation
})

preds <- parLapply(cl = mycl, c(names(models_ls)),function(x){
  d <- data_ls[[x]] %>% 
    filter(yday %in% timing_areas[[x]])
  m <- models_ls[[x]]
  
  pred <- data.frame(pred = as.numeric(predict(m,d)), lon = d$lon, lat = d$lat)
  
  coordinates(pred) <- ~lon+lat
  gridded(pred) <- T
  r <- raster(pred)
  proj4string(r) <- wgs
  
  #interpolate. for visualization purposes
  surf.1 <- Tps(as.matrix(as.data.frame(r,xy = T)[,c(1,2)],col = 2),as.data.frame(r,xy = T)[,3])
  
  grd <- expand.grid(x = seq(from = extent(r)[1],to = extent(r)[2],by = 0.1),
                     y = seq(from = extent(r)[3],to = extent(r)[4],by = 0.1))
  
  grd$coords <- matrix(c(grd$x,grd$y),ncol=2)
  
  surf.1.pred <- predict.Krig(surf.1,grd$coords)
  interpdf <- data.frame(grd$coords,surf.1.pred)
  
  colnames(interpdf)<-c("lon","lat","delta_t")
  
  coordinates(interpdf) <- ~lon+lat
  gridded(interpdf) <- TRUE
  interpr <- raster(interpdf)
  proj4string(interpr) <- wgs
  
  return(interpr)
  
})

stopCluster(mycl)

names(preds) <- names(models_ls)

#thin plate splines have also interpolated values of delta-t over land! let's remove them

#open a costline shapefile from https://www.arcgis.com/home/item.html?id=5cf4f223c4a642eb9aa7ae1216a04372 and unzip

coastline <- st_read("/continent_shapefile/continent.shp") %>% 
  st_crop(xmin = -135, xmax = 157, ymin = -61, ymax = 73) %>%
  st_union()

#mask all regions with the coastline layer
preds_filt <- lapply(preds, function(x){
  x_f <- raster::mask(x,as(coastline,"Spatial"), inverse = T) 
  x_f
})

#now for each region, further crop as necessary

#Americas
preds_filt$Americas <- mask(preds_filt$Americas, np_ocean, inverse = T)

#Europe
preds_filt$Europe <- mask(preds_filt$Europe, as(eur_sea,"Spatial"), inverse = F)

#Indian ocean
pg_rs <- world_seas %>% 
  filter(NAME %in% c("Persian Gulf","Red Sea"))

preds_filt$Indian_ocean <- mask(preds_filt$Indian_ocean, pg_rs, inverse = T)

#East Asia
ea <- world_seas %>% 
  filter(NAME %in% c("Andaman or Burma Sea", "Malacca Strait")) %>% 
  st_buffer(0.1) #some points on the edges remain after filtering, so make a buffer to make sure they all go away

preds_filt$East_asia <- mask(preds_filt$East_asia, ea, inverse = T)




## STEP 4: plotting ####

#4.1---------
#effect plots

#prepare a plotting function for summed effects (incl. intercept) mgcv default plots are partial.



#4.2--------------
#regional GAM maps

#to make the legend later, I need a raster layer that includes all possible delta-t values in my prediction layers, so from -1 to 5.
#use one of the prected rasters as the base and assign the desired min and max values to it

