#global seascapes 
#using wind and delta t downloaded in t_wind_data_download_&_process.R
#script follows global_seascape.R
#Elham Nourani. Radolfzell, Germany. Feb. 17. 2020

#tutorials
#https://www.r-bloggers.com/advanced-sab-r-metrics-parallelization-with-the-mgcv-package/
  
library(parallel)
library(dplyr)
library(ncdf4)
library(raster)
library(maps)
library(lubridate)
library(maptools)
library(purrr)
library(miceadds)
library(sf)
library(dismo)

wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")

#input points for the analyses will be as follows. just to keep in mind and match months, etc.

load("R_files/all_spp_temp_sp_filtered_15km_alt_14days.RData") #alt_cmpl....make sure this is updated! in all_data_prep_analyze.R
pts <- alt_cmpl %>% 
  filter(period == "now")

table(pts$season,pts$month) #spring: c(2:4); autumn: c(8:9)
table(year(pts$date_time)) #2003:2019

##### STEP 1: process the data-create yearly averges (season and zone-specific) #####
## temp data #####
setwd("/home/enourani/Documents/ERA_INTERIM_data_0_60")

vname <- c("sst","t2m")
file_list <- list.files(pattern = "sst_t2m.nc",full.names = TRUE)

load("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/day_indices.RData") #days

#start the cluster
mycl <- makeCluster(detectCores() - 2)
mycl <- makeCluster(6)
clusterExport(mycl, list("vname","file_list","days")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  
  library(ncdf4)
  library(lubridate)
  library(dplyr)
})


parLapply(cl = mycl,file_list,function(x){
  
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
  row_names <- expand.grid(lon,lat,as.character(timestamp))
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,vname[1])), nrow = nlon * nlat * nt, ncol = 1), #array to vector to matrix
    matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1)))
  
  colnames(var_df) <- c("lon","lat","date_time",vname)   #set column names
  
  #remove points over land (NAs)
  sea_df <- var_df %>%
    na.omit() %>% #this omits all points over land
    mutate(row_n = row_number()) %>% 
    filter(row_n %in% days$row_n) %>% #keep only rows that are daytime according to "days"
    #mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
    mutate(delta_t = sst - t2m,
           yday = yday(date_time),
           hour = hour(date_time),
           month = month(date_time),
           year = year(date_time),
           zone = ifelse(lat <= 30,"twz","tmpz"),
           season = ifelse(month %in% c(2:4),"spring",
                           ifelse(month %in% c(8:10), "autumn", "other")))
  
  #take averages
  sea_avg <- sea_df %>% 
    filter(season != "other") %>% 
    group_by(zone, season,lon,lat) %>% 
    summarise(avg_delta_t = mean(delta_t,na.rm = T)) %>% 
    as.data.frame()
  
  
  #identify rows at daytime. only do this for the first file and use for others to subset days
  # sea_df_crds <- as.matrix(sea_df[,c("lon","lat")])
  # dn <- sea_df %>%
  #   mutate(daynight = ifelse(solarpos(sea_df_crds,date_time)[,2] < -6, "night","day")) %>% #second vlaue returned by solarpos is the solar elevation angle. first one is azimuth
  #   mutate(row_n = row_number()) 
  # days <- dn %>%
  #   filter(daynight == "day")
  #   
  # save(days, file = "/home/enourani/ownCloud/Work/Projects/delta_t/R_files/day_indices.RData")
  
  save(sea_avg,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/",
                            paste(year(timestamp)[1],"avg_delta_t",sep = "_"),".RData"))

  
})

stopCluster(mycl)

#-------------------------------------------------------------------------------
## wind data #####
setwd("/home/enourani/Documents/ERA_INTERIM_data_0_60")

vname <- c("u","v")
file_list <- list.files(pattern = "wind.nc",full.names = TRUE)

load("/home/enourani/ownCloud/Work/Projects/delta_t/R_files/day_indices_wind.RData") #days

#start the cluster

mycl <- makeCluster(8)
clusterExport(mycl, list("vname","file_list","days")) #define the variable that will be used within the function

clusterEvalQ(mycl, {
  
  library(ncdf4)
  library(lubridate)
  library(dplyr)
})


parLapply(cl = mycl,file_list,function(x){
  
  nc <- nc_open(x)
  
  #extract lon and lat
  lat <- ncvar_get(nc,'latitude')
  nlat <- dim(lat) 
  lon <- ncvar_get(nc,'longitude') #this is from 0 - 360. convert to -180 to 180
  nlon <- dim(lon) 
  
  #extract the time
  t <- ncvar_get(nc, "time")
  nt <- dim(t)
  
  #convert the hours into date + hour
  timestamp <- as_datetime(c(t*60*60),origin = "1900-01-01")
  
  #put everything in a large df
  row_names <- expand.grid(lon,lat,as.character(timestamp))
  
  var_df <- data.frame(cbind(
    row_names,
    matrix(as.vector(ncvar_get(nc,vname[1])), nrow = nlon * nlat * nt, ncol = 1), #array to vector to matrix
    matrix(as.vector(ncvar_get(nc,vname[2])), nrow = nlon * nlat * nt, ncol = 1)))
  
  colnames(var_df) <- c("lon","lat","date_time",vname)   #set column names
  
  #  #identify rows at daytime. only do this for the first file and use for others to subset days
  # var_df_crds <- as.matrix(var_df[,c("lon","lat")])
  # dn <- var_df %>%
  #   mutate(date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
  #   mutate(daynight = ifelse(solarpos(var_df_crds,date_time)[,2] < -6, "night","day")) %>% #second vlaue returned by solarpos is the solar elevation angle. first one is azimuth
  #   mutate(row_n = row_number())
  # days <- dn %>%
  #   filter(daynight == "day")
  # 
  # save(days, file = "/home/enourani/ownCloud/Work/Projects/delta_t/R_files/day_indices_wind.RData")
  #  
  
  var_avg <- var_df %>%
    mutate(row_n = row_number(),
           date_time = as.POSIXct(strptime(date_time,format = "%Y-%m-%d %H:%M:%S"),tz = "UTC")) %>% 
    filter(row_n %in% days$row_n) %>% #keep only rows that are daytime according to "days"
    mutate(yday = yday(date_time),
           hour = hour(date_time),
           month = month(date_time),
           year = year(date_time),
           zone = ifelse(lat <= 30,"twz","tmpz"),
           season = ifelse(month %in% c(2:4),"spring",
                           ifelse(month %in% c(8:10), "autumn", "other"))) %>% 
    filter(season != "other") %>% 
    group_by(zone, season,lon,lat) %>% 
    summarise(avg_u = mean(u,na.rm = T),
              avg_v = mean(v,na.rm = T)) %>% 
    as.data.frame()

  save(var_avg,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/",
                             paste(year(timestamp)[1],"avg_wind",sep = "_"),".RData"))
  
  
})

stopCluster(mycl)



##### STEP 2: create averages for the study period: 2003-2018 #####
#open files, filter land and lakes, take average 

ocean <- st_read("/home/enourani/ownCloud/Work/GIS_files/ne_110m_ocean/ne_110m_ocean.shp")[2,] #remove caspian sea
setwd("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg")

#temperature
t_files <- list.files(pattern = "delta_t",full.names = TRUE)[c(25:40)] #filter for years that correspond to empirical data

t_avg_3_18 <- lapply(t_files,load.Rdata2) %>%
  reduce(rbind) %>% 
  group_by(zone, season, lon,lat) %>% 
  summarise(avg_delta_t = mean(avg_delta_t,na.rm = T)) %>% 
  ungroup()

#wind
w_files <- list.files(pattern = "wind",full.names = TRUE)[c(25:40)] #filter for years that correspond to empirical data

w_avg_3_18 <- lapply(w_files,load.Rdata2) %>%
  reduce(rbind) %>% 
  group_by(zone, season, lon,lat) %>% 
  summarise(avg_u = mean(avg_u,na.rm = T),
            avg_v = mean(avg_v,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(lon = lon -180) #longitude is from 0 - 360

#create raster layers
for (s in (c("spring","autumn"))){
  for (z in c("twz","tmpz")){
    #delta_t
    t <- t_avg_3_18 %>% 
      filter(season == s & zone == z) %>% 
      as.data.frame()
    
    coordinates(t) <-~ lon + lat
    gridded(t) <- TRUE
    dr <- raster(t,layer = 3)
    
    #subset lakes and land
    dr_masked <- mask(dr,ocean) #remove lakes
    
    save(dr_masked,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_delta_t.RData", sep = "_")))
    
    #-------------------
    #wind
    w <- w_avg_3_18 %>% 
      filter(season == s & zone == z) %>% 
      as.data.frame()
    coordinates(w) <-~ lon + lat
    gridded(w) <- TRUE
    ur <- raster(w,layer = 3)
    ur_masked <- mask(ur, ocean) #remove lakes and land
    vr <- raster(w,layer = 4)
    vr_masked <- mask(vr,ocean)
    
    save(ur_masked,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_u.RData", sep = "_")))
    save(vr_masked,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_v.RData", sep = "_")))

  }
}


##### STEP 3: build MaxEnt model #####

setwd("/home/enourani/ownCloud/Work/Projects/delta_t")

#open presence points. makes sure to discard 2019... no env data on Era-Interim
load("R_files/all_spp_temp_sp_filtered_15km_alt_14days.RData") #alt_cmpl....make sure this is updated! in all_data_prep_analyze.R
presence_ls <- alt_cmpl %>% 
  filter(period == "now"& year(timestamp) != 2019) %>% 
  group_split(season,zone) 

names(presence_ls) <- c("autumn_tmpz","autumn_twz","spring_tmpz","spring_twz")
presence_ls <- lapply(presence_ls, "[",,c("lon","lat"))  #only keep lon and lat columns

#open raster data and create one raster stack per model
raster_ls <- lapply(names(presence_ls), function(x){
  files <- list.files("ERA_INTERIM_yearly_avg/rasters", pattern = x, full.names = TRUE)
  rs <- lapply(files,load.Rdata2) 
  stack(rs)
})
names(raster_ls) <- names(presence_ls)

#generate background points

bg_ls <- lapply(presence_ls, function(x){
  #create a buffer around the points to select background points from
  buff <- circles(x, d=1000000, lonlat=T) #create buffer around the points
  buff <- intersect(buff@polygons,as(ocean, "Spatial")) #make sure the buffer is not over land
  bg <- spsample(buff, 10000, type='random', iter=1000)
  bg <- as.data.frame(bg)
  colnames(bg) <- c("lon","lat")
  bg
})
names(bg_ls) <- names(presence_ls)


#mask_tmpz <- subset(raster_ls[[1]],1) #use one raster layer as a masking layer
#mask_twz <- subset(raster_ls[[2]],1)
#points(randomPoints(mask_tmpz,1000))

#build the model
#make sure to set removeDuplicates to TRUE so >1 observations withn the same grid cell are removed

maxentmodel<-maxent(raster_ls[[4]],p=as.data.frame(presence_ls[[4]]), a=bg_ls[[4]], removeDuplicates=T,args=c("responsecurves","jackknife","nothreshold","nohinge","product","noautofeature","maximumiterations=1000" ),
                    path="R_files/maxent")


