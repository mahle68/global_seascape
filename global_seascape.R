#global seascapes 
#local time conversion: tried googleways, the API key doesnt work because my bank account hasn't been verified. tried using a shapefile of time zones and
#extracting the time zone id for each spatial point. then found a package lutz that gives the time zone based on lat and lon and seems to be more accurate 
#than the shapefile
#using RStudio verison control to keep track of updates.

#tutorials
#https://www.r-bloggers.com/advanced-sab-r-metrics-parallelization-with-the-mgcv-package/
  
library(mgcv)
library(parallel)
library(dplyr)
#library(googleway) #for local time calculations.. needs google API
#library(XML)
library(sf)
library(raster)
library(maps)
library(lubridate)
library(lutz)
library(maptools)
library(purrr)
#library(multidplyr)



setwd("C:/Users/mahle/ownCloud/Work/Projects/delta_t/")


#define variables
wgs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
mycl <- makeCluster(detectCores() - 1) 


###step 1: open data ####
file_ls <- list.files("ERA_INTERIM_data_global_sampled",".RData",full.names = TRUE) #data processed and sampled in ecmwf_in_R.R, actually, in t_data_download_&_process

data_df <- file_ls %>%
  map(read_csv) %>% 
  reduce(rbind) %>%
  mutate(year = year(date_time))

save(data_df,file = "R_files/processed_era_interim_data/samples_df.RData")



###step 2: convert hour of day to local time #####
load("R_files/processed_era_interim_data/samples_df.RData") #called data_df

df_lt <- data_df %>%
  mutate(tz = tz_lookup_coords(lat,lon)) %>%
  rowwise() %>%
  mutate(local_date_time = as.character(as.POSIXlt(x = date_time, tz = tz))) %>% #has to be character. otherwise I get an error for POSIXlt format not being supported
  mutate(local_hour = strsplit(strsplit(local_date_time,split = " ")[[1]][2],split = ":")[[1]][1]) %>%
  as.data.frame()#local hour is NA when it should be 00

#conver NA local_hours to 00
df_lt[is.na(df_lt$local_hour),"local_hour"] <- "00"

save(df_lt,file = "R_files/processed_era_interim_data/samples_with_local_time_hour.RData")


###step 3: remove lakes! #####

load("R_files/processed_era_interim_data/samples_with_local_time_hour.RData") #called df_lt
df_lt$group <- rep(1:3, length.out = nrow(df_lt), each = ceiling(nrow(df_lt)/3))

land <- shapefile("C:/Users/mahle/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp")
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