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
  #mutate(lon = lon -180) #longitude is from 0 - 360

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
    dr_masked <- mask(dr,ocean) 
    save(dr_masked,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_delta_t.RData", sep = "_")))
    
    #-------------------
    #wind
    w <- w_avg_3_18 %>% 
      filter(season == s & zone == z) %>% 
      as.data.frame()
    coordinates(w) <-~ lon + lat
    gridded(w) <- TRUE
    proj4string(w) <- wgs
    ur <- raster(w,layer = 3)
    ur_masked <- mask(ur, ocean)
    vr <- raster(w,layer = 4)
    vr_masked <- mask(vr,ocean)
    
    save(ur,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_u.RData", sep = "_")))
    save(vr,file = paste0("/home/enourani/ownCloud/Work/Projects/delta_t/ERA_INTERIM_yearly_avg/rasters/",paste(s,z,"2003_2018_v.RData", sep = "_")))

  }
}

#remove lakes

###step 3: remove lakes! #####

load("R_files/processed_era_interim_data/samples_with_local_time_hour.RData") #called df_lt
df_lt$group <- rep(1:3, length.out = nrow(df_lt), each = ceiling(nrow(df_lt)/3))

land <- shapefile("C:/Users/mahle/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp")

ocean <- 

land$land <- 1
land <- unionSpatialPolygons(land,IDs = land2$land)
land_sf <- st_as_sf(land,crs = wgs)

#start the cluster
clusterExport(mycl, list("land_sf","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, {#load packages for each node
  library(sf)
  #library(raster)
  library(dplyr)
})



ls_lt_lakes <- parLapply(cl = mycl,split(df_lt,df_lt$group),function(x){
  
  df_lt_lakes <- x %>%
    st_as_sf(coords = c("lon", "lat"), crs = wgs) %>%
    mutate(lon = as.numeric(st_coordinates(.)[,"X"]),
           lat = as.numeric(st_coordinates(.)[,"Y"])) %>%
    st_intersection(y = land_sf) %>%
    st_drop_geometry() %>%
    as.data.frame()
  
  df_lt_no_lake
  
})

stopCluster(mycl)
save(ls_lt_lakes,file = "R_files/processed_era_interim_data/samples_with_local_time_hour_over_lakes_ls.RData")


#errors following the cluster computing attempt. trying it in QGIS



###step 4: visualizations #####
load("R_files/processed_era_interim_data/sample_of_samples_no_lakes.RData")

#plot delta_t agains yday, separtely for each latitudinal zone
windows(12,13)
par(mfrow = c(2,1),
    par(oma = c(5.1, 0,0,0), xpd = NA))

#noon
plot(delta_t ~ yday, data = sample2[sample2$local_hour == 12,],
     type = "n",
     xlab = "day of the year",
     ylab = "delta T",
     bty = "n",
     ylim = c(-2.2,5), xlim = c(0,365), main = "global delta_t at local noon")
axis(1,seq(0,350,50), c(seq(0,350,50)))
lines(x = c(0,350),y = c(0,0),lty = 2, col = "gray")
lines(lowess(sample2[between(sample2$lat,0,30) & sample2$local_hour == 12,c("yday","delta_t")]), col = "deepskyblue",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-30,0) & sample2$local_hour == 12,c("yday","delta_t")]), col = "forestgreen",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,30,60) & sample2$local_hour == 12,c("yday","delta_t")]), col = "darkviolet",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-60,-30) & sample2$local_hour == 12,c("yday","delta_t")]), col = "firebrick",lwd = 2) 

#mid-day
plot(delta_t ~ yday, data = sample2[sample2$local_hour %in% c(11:15),],
     type = "n",
     xlab = "day of the year",
     ylab = "delta T",
     bty = "n",
     ylim = c(-1,2), xlim = c(0,365), main = "global delta_t at local mid-day (11:00 - 15:00)")
axis(1,seq(0,350,50), c(seq(0,350,50)))
lines(x = c(0,350),y = c(0,0),lty = 2, col = "gray")
lines(lowess(sample2[between(sample2$lat,0,30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","delta_t")]), col = "deepskyblue",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-30,0) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","delta_t")]), col = "forestgreen",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,30,60) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","delta_t")]), col = "darkviolet",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-60,-30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","delta_t")]), col = "firebrick",lwd = 2) 


legend("bottom",legend = c("0°- 30° N","0°- 30° S","30°- 60° N","30°- 60° S"),horiz = T, 
       col = c("deepskyblue", "forestgreen", "darkviolet", "firebrick"), bty = "n", cex = 0.9,lty = 1, lwd = 2,inset = c(0,-0.6))



#plot temperature
#sst
plot(sst ~ yday, data = sample2[sample2$local_hour %in% c(11:15),],
     type = "n",
     xlab = "day of the year",
     ylab = "delta T",
     bty = "n",
     ylim =  c(270,310), xlim = c(0,365), main = "global sst at local mid-day (11:00 - 15:00)")
axis(1,seq(0,350,50), c(seq(0,350,50)))
lines(lowess(sample2[between(sample2$lat,0,30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","sst")]), col = "deepskyblue",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-30,0) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","sst")]), col = "forestgreen",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,30,60) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","sst")]), col = "darkviolet",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-60,-30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","sst")]), col = "firebrick",lwd = 2) 

#t2m
plot(t2m ~ yday, data = sample2[sample2$local_hour %in% c(11:15),],
     type = "n",
     xlab = "day of the year",
     ylab = "delta T",
     bty = "n",
     ylim = c(270,310), xlim = c(0,365), main = "global t2m at local mid-day (11:00 - 15:00)")
axis(1,seq(0,350,50), c(seq(0,350,50)))
lines(lowess(sample2[between(sample2$lat,0,30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","t2m")]), col = "deepskyblue",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-30,0) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","t2m")]), col = "forestgreen",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,30,60) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","t2m")]), col = "darkviolet",lwd = 2) 
lines(lowess(sample2[between(sample2$lat,-60,-30) & between(as.numeric(sample2$local_hour), 11, 15),c("yday","t2m")]), col = "firebrick",lwd = 2) 

legend("bottom",legend = c("0°- 30° N","0°- 30° S","30°- 60° N","30°- 60° S"),horiz = T, 
       col = c("deepskyblue", "forestgreen", "darkviolet", "firebrick"), bty = "n", cex = 0.9,lty = 1, lwd = 2,inset = c(0,-0.6))


###step 5: modeling-testing-subset testing and training set #####
#include in the model: delta T ~ s(yday)+ s(lon+lat) + local hour of day (as a random effect or just a factor?)
#perhaps an interaction term for local hour and latitude to take into account seasonality?

#for modleing, use bam() instead of gam(), it breaks the data into chunks and has a lower memory footprint. also has the option to parallelize
#work on a sample first


sample <- df_lt %>% 
  sample_n(3000)

#open data
load("R_files/processed_era_interim_data/sample_of_samples_no_lakes.RData")

#break up into testing and testing set
train_set <- sample2 %>%
  sample_frac(0.6, replace = F)

test_set <- sample2 %>%
  anti_join(train_set)
  

#model on a cluster
clusterExport(mycl, "train_set") 

model <- bam(delta_t ~ s(lon,lat) + s(yday,bs="cc") + factor(local_hour) , data = train_set, cluster = mycl) #in the example: method="GCV.Cp",

stopCluster(mycl)

#model without a cluster
model <- bam(delta_t ~ s(lon,lat) + s(yday,bs="cc") + factor(local_hour) , data = train_set)
AIC(model) 
plot(model)


#separate smooths for hour
model2 <- bam(delta_t ~ s(lon,lat) + s(yday,bs="cc",by = factor(hour)) + factor(hour) , data = train_set)
AIC(model2) 
plot(model2)

#compare models with Wald test
anova(model, model2)

#evaluate the best model
gam.check(model)

#test the model
test_pred <- predict.bam(model2,test_set)
plot(test_pred,test_set$delta_t)
cor(test_pred,test_set$delta_t)

##### from previous codes

#-----------------------------------------------------------
#effect plots 
#-----------------------------------------------------------
#previous attempts are in the previous version ;)

pdf("delta_T_gam/yday_smoother_fig_model1_3.pdf", height = 4, width = 4)
#x11()
#dev.new(width=4, height=4)

par(mfrow = c(1,1), bty="n", #no box around the plot
    cex.axis = 0.75, #x and y labels have 0.75% of the default size
    font.axis = 3, #axis labels are in italics
    cex.lab = 1
)


plot(model2,select = 2, xlab = "Julian date",ylab = "Delta T", ylim = c(-2,2), shade = T, shade.col = "grey75",
     scheme= 0,se = 12,bty = "l",labels = FALSE, tck = 0) #plot only the second smooth term

#rect(xleft=268,ybottom=-2.3,xright=296,ytop=2.3, col="#96CDCD30",border=NA) #for juv...#last two digits are for transparency
#rect(xleft=242,ybottom=-2.3,xright=277,ytop=2.3, col="#FA807230",border=NA) #for adults

rect(xleft = 268,ybottom=-2.3,xright = 296,ytop = 1.14, col = "#96CDCD30",border = NA) #for juv...#last two digits are for transparency
rect(xleft = 213,ybottom=0.35,xright = 268,ytop = 1.14, col = "#96CDCD30",border = NA)

rect(xleft=242,ybottom=-2.3,xright=277,ytop=0.6, col="#FA807230",border=NA) #for adults
rect(xleft=213,ybottom=-0.5,xright=242,ytop=0.6, col="#FA807230",border=NA)

axis(side= 1, at= c(213,240,260,280,305), line=-0, labels= c(213,240,260,280,305), 
     tick=T , col.ticks = 1, col=NA, tck=-.015)
axis(side= 2, at= c(-1.5,0,1.5), line=0, labels= c(-1.5,0,1.5),
     tick=T , col.ticks = 1,col=NA, tck=-.015, 
     las=2) # text perpendicular to axis label 
legend(x=212, y=2, legend=c("juveniles", "adults"), fill=c("#96CDCD40","#FA807240"), #coords indicate top-left
       cex=0.8, bg="white",bty="n")

dev.off()

#----------------------------------------------------------------------------------------
# Predict with the model!
#----------------------------------------------------------------------------------------
#make predictions with the model for the periods of adult and juvenile migration
#juvenile Sep. 25-Oct. 23 (i.e. ydays 268-296); adults Aug 30-Oct. 4 (i.e. ydays 242-277)
#create a set for each..... conclusion: the two prediction maps look very similar. so, just do one map for the entire autumn
seasons <- list(
  juv = data_aut_sea[data_aut_sea$yday %in% c(268:296),],
  adlt = data_aut_sea[data_aut_sea$yday %in% c(242:277),]
)

#predict
preds_ls<-lapply(seasons, function(x){
  pred<-data.frame(pred=as.numeric(predict(model,x)),lon=x$lon,lat=x$lat)
  
  coordinates(pred)<-~lon+lat
  gridded(pred)<-T
  r<-raster(pred)
  
  return(r)
})


save(preds_ls,file="delta_T_gam/predictions_with_model1.RData")
load("delta_T_gam/predictions_with_model1.RData")


###step 6: modleing-full model #####

###step 7: visualization #####
#predict with the full model for autumn migration season

# interpolate to 1 km for visualization purposes



#####
windows()
maps::map("world",fill = TRUE,col = "gray",border = "gray")
points(sample$lon,sample$lat,cex = 0.5,pch = 16,col = "red")
points(sample2$lon,sample2$lat,cex = 0.3,pch = 16,col = "green")
points(train_set$lon,train_set$lat,cex = 0.5,pch = 16,col = "red")
points(test_set$lon,test_set$lat,cex = 0.5,pch = 16,col = "blue")

plot(st_geometry(land_pol))
lines(land)
plot(lakes,col = "red",add = T)