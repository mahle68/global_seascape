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
library(multidplyr)



setwd("C:/Users/mahle/ownCloud/Work/Projects/delta_t/")
wgs<-CRS("+proj=longlat +datum=WGS84 +no_defs")



###step 1: open data ####
file_ls<-list.files("ERA_INTERIM_data_global_sampled",".RData",full.names= TRUE) #data processed and sampled in ecmwf_in_R.R
data_df<-data.frame()
for(i in file_ls){
  load(i)
  data_df<-rbind(data_df,sample)
  data_df
}

data_df$year<-year(data_df$date_time)

save(data_df,file="R_files/processed_era_interim_data/samples_df.RData")


###step 2: convert hour of day to local time #####
load("R_files/processed_era_interim_data/samples_df.RData") #called data_df

df_lt<-data_df%>%
  mutate(tz=tz_lookup_coords(lat,lon))%>%
  rowwise()%>%
  mutate(local_date_time= as.character(as.POSIXlt(x=date_time, tz = tz)))%>% #has to be character. otherwise I get an error for POSIXlt format not being supported
  mutate(local_hour=strsplit(strsplit(local_date_time,split=" ")[[1]][2],split=":")[[1]][1])%>%
  as.data.frame()#local hour is NA when it should be 00

#conver NA local_hours to 00
df_lt[is.na(df_lt$local_hour),"local_hour"]<-"00"

save(df_lt,file="R_files/processed_era_interim_data/samples_with_local_time_hour.RData")


#time zone calc using shapefile. not used #####
tzs <- st_read("C:/Users/mahle/ownCloud/Work/GIS_files/time_zones/combined-shapefile-with-oceans.shp", quiet = TRUE)

data_df_lt<-data_df%>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
  mutate(lon=as.numeric(st_coordinates(.)[,"X"]),
         lat=as.numeric(st_coordinates(.)[,"Y"]))%>%
  st_join(tzs)%>%
  st_drop_geometry()%>%
  rowwise()%>%
  mutate(local_date_time= as.character(as.POSIXlt(x=date_time, tz = as.character(tzid))))%>% #has to be character. otherwise I get an error for POSIXlt format not being supported
  as.data.frame()

save(data_df_lt,file="R_files/processed_era_interim_data/samples_with_local_time.RData")
#####



###step 3: remove lakes! #####

load("R_files/processed_era_interim_data/samples_with_local_time_hour.RData") #called df_lt
df_lt$group<-rep(1:3, length.out = nrow(df_lt), each = ceiling(nrow(df_lt)/3))

land<-shapefile("C:/Users/mahle/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp")
land$land<-1
land<-unionSpatialPolygons(land,IDs=land2$land)
land_sf<-st_as_sf(land,crs=wgs)

#start the cluster
mycl<-makeCluster(detectCores()-1) #define workers :p

clusterExport(mycl, list("land_sf","wgs")) #define the variable that will be used within the function

clusterEvalQ(mycl, { #load packages for each node
  library(sf)
  #library(raster)
  library(dplyr)
})



ls_lt_no_lake<-parLapply(cl=mycl,split(df_lt,df_lt$group),function(x){
  
  df_lt_no_lake<-x%>%
    st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
    mutate(lon=as.numeric(st_coordinates(.)[,"X"]),
           lat=as.numeric(st_coordinates(.)[,"Y"]))%>%
    st_difference(y=land_sf)%>%
    st_drop_geometry()%>%
    as.data.frame()
  
  df_lt_no_lake
  
})


save(ls_lt_no_lake,file="R_files/processed_era_interim_data/samples_with_local_time_hour_no_lakes_ls.RData")



#remove points that overlap with the lakes
df_lt_no_lake<-df_lt%>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
  mutate(lon=as.numeric(st_coordinates(.)[,"X"]),
         lat=as.numeric(st_coordinates(.)[,"Y"]))%>%
  st_difference(y=land_sf)%>%
  st_drop_geometry()%>%
  as.data.frame()


#run in parallel
cluster <- new_cluster(2)
cluster_assign_each(cluster, filename = c("land", "df_lt"))
cluster_send(cluster, my_data <- vroom::vroom(filename))


sample2<-sample%>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
  mutate(lon=as.numeric(st_coordinates(.)[,"X"]),
         lat=as.numeric(st_coordinates(.)[,"Y"]))%>%
  st_difference(y=land_sf)%>%
  st_drop_geometry()%>%
  as.data.frame()




#sample_sf<-sample%>%
#  st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
#  st_difference(y=land_sf)

#land3<-st_read("C:/Users/mahle/ownCloud/Work/GIS_files/ne_10m_land/ne_10m_land.shp")%>%
#st_combine()

sample_sf<-st_as_sf(coords = c("lon", "lat"), crs = wgs)

sample_years<-sample%>%
  group_by(year)%>%
  partition(cluster)

d<-Sys.time()
sample2<-sample_years%>%
  st_as_sf(coords = c("lon", "lat"), crs = wgs)%>%
  mutate(lon=as.numeric(st_coordinates(.)[,"X"]),
         lat=as.numeric(st_coordinates(.)[,"Y"]))%>%
  st_difference(y=land_sf)%>%
  st_drop_geometry()%>%
  as.data.frame()%>%
  collect()
Sys.time()-d







save(sample2,file="R_files/processed_era_interim_data/sample_of_samples_no_lakes.RData")



land$land<-1
land<-unionSpatialPolygons(land,IDs=land$land)
land_sf<-st_as_sf(land,crs=wgs)




load("R_files/processed_era_interim_data/samples_with_local_time_hour.RData") #called df_lt
lakes<-shapefile("C:/Users/mahle/ownCloud/Work/GIS_files/world_lakes/ne_10m_lakes.shp")
lakes2<-st_buffer(lakes, dist = 2)

lakes$lake<-1
lakes<-unionSpatialPolygons(lakes,IDs=lakes$lake)
lakes_sf<-st_as_sf(lakes,crs=wgs)



save(df_lt_no_lake,file="R_files/processed_era_interim_data/samples_with_local_time_hour_no_lakes.RData") #caspian sea is included :(


land_l<-shapefile("C:/Users/mahle/ownCloud/Work/GIS_files/world_continents/continent_ln.shp") #this is a spatial lines file
land_p<-SpatialLines2PolySet(land_l)%>%
PolySet2SpatialPolygons()


land_pol <- st_polygonize(land)
land_pol <- as(land_pol, "Spatial") 

land_pol$land<-1

land<-unionSpatialPolygons(land_pol,IDs=land_pol$land)

class(shp_airports)
###step 4: modeling-testing-subset testing and training set #####

sample<-df_lt%>%
  sample_n(3000)



#include in the model: delta T ~ s(yday)+ s(lon+lat) + local hour of day (as a random effect or just a factor?)
#perhaps an interaction term for local hour and latitude to take into account seasonality?

#for modleing, use bam() instead of gam(), it breaks the data into chunks and has a lower memory footprint. also has the option to parallelize
#work on a sample first


mycl<-makeCluster(detectCores()-1) #define workers :p

clusterExport(mycl, "sample") #define the variable that will be used within the function


model<-bam(delta_t ~ s(lon,lat) + s(yday) + factor(local_hour) , data=sample, cluster=mycl) #in the example: method="GCV.Cp",


stopCluster(mycl)


model<-bam(delta_t ~ s(lon,lat) + s(yday) + factor(local_hour) , data=sample)
AIC(model) 
plot(model)

#separate smooths for hour
model2<-bam(delta_t ~ s(lon,lat) + s(yday,by=factor(hour)) + factor(hour) , data=sample, cluster=mycl)
AIC(model2) 
plot(model2)

#compare models with Wald test
anova(model, model2)

#save the best model

##### from previous codes

#-----------------------------------------------------------
#effect plots 
#-----------------------------------------------------------
#previous attempts are in the previous version ;)

pdf("delta_T_gam/yday_smoother_fig_model1_3.pdf", height=4, width=4)
#x11()
#dev.new(width=4, height=4)

par(mfrow=c(1,1), bty="n", #no box around the plot
    cex.axis= 0.75, #x and y labels have 0.75% of the default size
    font.axis= 3, #axis labels are in italics
    cex.lab=1
)


plot(model2,select=2, xlab="Julian date",ylab="Delta T", ylim= c(-2,2), shade=T, shade.col= "grey75",
     scheme=0,se=12,bty="l",labels = FALSE, tck=0) #plot only the second smooth term

#rect(xleft=268,ybottom=-2.3,xright=296,ytop=2.3, col="#96CDCD30",border=NA) #for juv...#last two digits are for transparency
#rect(xleft=242,ybottom=-2.3,xright=277,ytop=2.3, col="#FA807230",border=NA) #for adults

rect(xleft=268,ybottom=-2.3,xright=296,ytop=1.14, col="#96CDCD30",border=NA) #for juv...#last two digits are for transparency
rect(xleft=213,ybottom=0.35,xright=268,ytop=1.14, col="#96CDCD30",border=NA)

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
seasons<-list(
  juv=data_aut_sea[data_aut_sea$yday %in% c(268:296),],
  adlt=data_aut_sea[data_aut_sea$yday %in% c(242:277),]
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
map("world",fill=TRUE,col="gray",border="gray")
points(sample$lon,sample$lat,cex=0.5,pch=16,col="red")
points(sample2$lon,sample2$lat,cex=0.3,pch=16,col="green")
plot(st_geometry(land_pol))
lines(land)
plot(lakes,col="red",add=T)