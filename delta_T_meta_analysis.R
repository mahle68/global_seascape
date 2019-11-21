#Scripts for meta-analysis of relationship between delta_T and sea-crossing behavior in soaring birds
#Elham Nourani
#Oct. 30. 2019
#Radolfzell, Germany


library(dplyr)
library(purrr)
library(lubridate)
library(readxl)
library(raster)
library(maptools)
library(rgdal)
library(sf)

setwd("C:/Users/mahle/ownCloud/Work/Projects/delta_t")
wgs<-"+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
meters_proj <- CRS("+proj=moll +ellps=WGS84")



#--------PILOT STUDY---------------- Oriental honey buzzard

#geographic files
load("R_files/land_east_asia.RData") #called land_asia
load("R_files/land_east_asia_15km.RData") #called land_asia_15km

#temporally fitlered data for GBF and OHB
load("R_files/GFB_HB_temp_filtered.RData") #called data_sp



#open files (raw. not filtered.)
#---STEP 1: temporal filter for migratory seasons and filter for location classes (for PTTs).

OHB_files<-list.files("data/Oriental_honey_buzzard",pattern=".xls",full.names=T)
OHB<-lapply(OHB_files,read_excel,1,col_types=c("numeric","date","numeric","numeric","numeric","skip","text",rep("numeric",8)))%>%
  reduce(full_join)%>%
  rename(date_time='date(gmt)',lon=longitud,lat=latitude)%>%
  mutate(yday=yday(date_time))%>%
  mutate(season=ifelse(between(yday,253,294),"autumn",ifelse(month==5,"spring","other")))%>% #11 Sep-20 Oct; spring between 1-5 May
  filter(season == "autumn" & #no sea-crossing in spring
           class %in% c("0","1","2","3")) #filter for location classes


GFB_files<-list.files("data/Grey_faced_buzzard/",pattern=".csv",full.names=T)
GFB<-lapply(GFB_files,read.csv,stringsAsFactors=F)%>%
  reduce(full_join)%>% #is locdate is in UTC
  mutate(locdate,date_time= as.POSIXct(strptime(locdate,format="%Y-%m-%d %H:%M:%S"),tz="UTC"))%>%
  mutate(month=month(date_time),
         year=year(date_time),
         season= ifelse(month %in% c(3,4),"spring",ifelse(month %in% c(10),"autumn","other")))%>% #(Nourani et al 2017 for autumn; no sea-crossing in spring based on visually exploring the data.)
  filter(season %in% c("spring","autumn") &
         class %in% c("0","1","2","3")) #filter for location classes

#put the two dataframes together as a list
data<-list(GFB, as.data.frame(OHB))

data_sp<-lapply(data,function(x){
  coordinates(x)<-~lon+lat
  proj4string(x)<-wgs
  x
  })

save(data_sp,file="R_files/GFB_HB_temp_filtered.RData")

#---STEP 2: spatial filter for points over the sea. 
sea_data<-lapply(data_sp,function(x){
  pts_over_land<-over(x,land_asia) #make sure the spatialpolygon layer is not a dataframe
  pts_over_sea<-x[is.na(pts_over_land),]
  pts_over_sea
})

#save the points
save(sea_data,file="R_files/GFB_HB_temp_sp_filtered.RData")

#------------------------------------
#filter spatially with 15 km buffer
data_sf<-lapply(data_sp,st_as_sf,crs=wgs) #convert points to sf object

st_join(pts, poly, join = st_intersects)
as.data.frame(st_join(pts, poly, join = st_intersects))[2] %>% setNames("ID")

sea_data_15<-lapply(data_sp,function(x){
  pts_over_water<-st_difference(x,land_asia_15km)
  pts_over_water
})



sea_data_15<-lapply(data_sp,function(x){
  pts_over_land<-over(x,land_asia_15km) #make sure the spatialpolygon layer is not a dataframe
  pts_over_sea<-x[is.na(pts_over_land),]
  pts_over_sea
})

#save the points
save(sea_data_15,file="R_files/GFB_HB_temp_sp_filtered_15km.RData")


#---STEP 3: temporal filter for autocorrelation


#---STEP 4: annotae with delta T




#------------------misc

#creating landmass layer
land<-shapefile("C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files/countries.shp")
land$landmass<-1
landmasses<-unionSpatialPolygons(land,IDs=land$landmass)
landmasses$land<-1
writeOGR(landmasses, "C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files", "world_landmass", driver="ESRI Shapefile")

#creating the land_asia layer
land<-shapefile("C:/Users/mahle/Desktop/work/Max_Planck_postdoc/GIS files/world_landmass.shp") #open world landmass shapefile (codes for making the layer at the end of file)
land_asia<-crop(land,extent(90,143,-9,42.5))
land_asia<-as(land_asia,"SpatialPolygons")
save(land_asia,file="R_files/land_east_asia.RData")

land_asia_sf<-st_as_sf(land_asia,crs=wgs) #convert to sf object
land_asia_sf_m<-st_transform(land_asia_sf,meters_proj)

med_sf1<-st_buffer(land_asia_sf_m, dist = units::set_units(15000, 'm')) #put a buffer around the land polygon, to eliminate short sea-crossings. buffer is set to 15 km to match the 30 km res of env data
land_asia_15km<-st_transform(med_sf1,wgs)

save(land_asia_15km,file="R_files/land_east_asia_15km.RData")

#-----------------plotting
windows()
plot(land_asia)
points(data_sp[[1]],col=data_sp[[1]]$month,pch=16,cex=0.5)
points(GFB$lon,GFB$lat,col=as.factor(GFB$season))
points(GFB[GFB$season == "spring","lon"],GFB[GFB$season == "spring","lat"],pch=16,cex=0.4,col="green")
points(GFB[GFB$season == "autumn","lon"],GFB[GFB$season == "autumn","lat"],pch=16,cex=0.4,col="blue")
points(GFB[GFB$month == 10,"lon"],GFB[GFB$month == 10,"lat"],pch=16,cex=0.4,col="blue")
points(GFB[GFB$month == 9,"lon"],GFB[GFB$month == 9,"lat"],pch=16,cex=0.4,col="orange")
points(GFB[GFB$month == 3,"lon"],GFB[GFB$month == 3,"lat"],pch=16,cex=0.4,col="pink")
points(GFB[GFB$month == 4,"lon"],GFB[GFB$month == 4,"lat"],pch=16,cex=0.4,col="purple")

lapply(data_sp,points,pch=16,cex=0.3,col="blue")
lapply(sea_data,points,pch=16,cex=0.5,col="orange")

library(mapview)
mapview(x, col.regions = sf.colors(10))
mapview(pts_over_water, col.regions = "red")
